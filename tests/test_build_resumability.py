"""The network build must survive being interrupted.

A full build is hours of API calls, so it has to be runnable across several
sittings. These tests use fake layers so nothing touches the network.
"""

import os
import pickle
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "maintenance"))

import build_networks as bn                            # noqa: E402


class Layer:
    """Stand-in for a populated omics layer."""

    def __init__(self, name, size=0):
        self.name = name
        self.size = size

    def __eq__(self, other):
        return (self.name, self.size) == (other.name, other.size)


class TestCheckpointRoundTrip:
    def test_a_saved_layer_comes_back(self, tmp_path):
        layer = Layer("proteome", 3112)
        bn.save_checkpoint(str(tmp_path), "mouse", "proteome_base", layer)
        assert bn.load_checkpoint(str(tmp_path), "mouse", "proteome_base") == layer

    def test_a_missing_checkpoint_is_none(self, tmp_path):
        assert bn.load_checkpoint(str(tmp_path), "mouse", "nothing") is None

    def test_checkpoints_are_per_organism(self, tmp_path):
        bn.save_checkpoint(str(tmp_path), "mouse", "pathways", Layer("m"))
        assert bn.load_checkpoint(str(tmp_path), "human", "pathways") is None

    def test_a_corrupt_checkpoint_is_rebuilt_not_fatal(self, tmp_path):
        path = os.path.join(bn.checkpoint_dir(str(tmp_path), "mouse"), "x.pkl")
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "wb") as handle:
            handle.write(b"not a pickle")
        assert bn.load_checkpoint(str(tmp_path), "mouse", "x") is None

    def test_writing_is_atomic(self, tmp_path):
        """No .tmp file is left behind for a later run to trip over."""
        bn.save_checkpoint(str(tmp_path), "mouse", "pathways", Layer("m"))
        directory = bn.checkpoint_dir(str(tmp_path), "mouse")
        assert not [f for f in os.listdir(directory) if f.endswith(".tmp")]


class TestStep:
    def test_a_step_runs_once_and_is_reused(self, tmp_path):
        calls = []

        def build():
            calls.append(1)
            return Layer("proteome", len(calls))

        first = bn.step(str(tmp_path), "mouse", "proteome_base", build)
        second = bn.step(str(tmp_path), "mouse", "proteome_base", build)

        assert len(calls) == 1, "the step ran again instead of resuming"
        assert first == second

    def test_resume_false_forces_a_rebuild(self, tmp_path):
        calls = []

        def build():
            calls.append(1)
            return Layer("proteome", len(calls))

        bn.step(str(tmp_path), "mouse", "proteome_base", build)
        bn.step(str(tmp_path), "mouse", "proteome_base", build, resume=False)
        assert len(calls) == 2

    def test_a_failing_step_leaves_no_checkpoint(self, tmp_path):
        """A crash must not be remembered as a completed step."""
        def build():
            raise RuntimeError("API down")

        with pytest.raises(RuntimeError):
            bn.step(str(tmp_path), "mouse", "proteome_string", build)

        assert bn.load_checkpoint(str(tmp_path), "mouse", "proteome_string") is None

    def test_later_steps_resume_independently(self, tmp_path):
        """Killing during STRING must not cost the proteome populate."""
        done = []

        bn.step(str(tmp_path), "mouse", "proteome_base",
                lambda: done.append("base") or Layer("base"))

        with pytest.raises(RuntimeError):
            bn.step(str(tmp_path), "mouse", "proteome_string",
                    lambda: (_ for _ in ()).throw(RuntimeError("killed")))

        # Rerun: base is reused, string is retried.
        bn.step(str(tmp_path), "mouse", "proteome_base",
                lambda: done.append("base") or Layer("base"))
        bn.step(str(tmp_path), "mouse", "proteome_string",
                lambda: done.append("string") or Layer("string"))

        assert done == ["base", "string"], f"unexpected work: {done}"


class TestClearCheckpoints:
    def test_clearing_removes_them(self, tmp_path):
        bn.save_checkpoint(str(tmp_path), "mouse", "pathways", Layer("m"))
        bn.clear_checkpoints(str(tmp_path), "mouse")
        assert bn.load_checkpoint(str(tmp_path), "mouse", "pathways") is None

    def test_clearing_an_absent_organism_is_harmless(self, tmp_path):
        bn.clear_checkpoints(str(tmp_path), "nothing_here")


class TestReactionTableReuse:
    def test_an_existing_table_is_read_not_refetched(self):
        """The KEGG reaction table is ~12,000 REST calls; it must be reusable."""
        table = os.path.join(ROOT, "data", "master_reactions.csv")
        if not os.path.exists(table):
            pytest.skip("no master_reactions.csv in this checkout")

        frame = bn.load_reactions_from_file(table)
        assert len(frame) > 1000
        assert isinstance(frame.iloc[0]["substrates"], list)
        assert "reversible" in frame.columns


class TestStaleCheckpointsInAChain:
    """Adding --brenda to an organism whose ChIP-Atlas step was checkpointed
    ran BRENDA, then loaded the older ChIP-Atlas checkpoint over the result."""

    def test_a_rebuilt_step_invalidates_later_steps_in_its_chain(self, tmp_path):
        import build_networks as bn

        calls = []
        chain = {"rebuilt": False}
        bn.save_checkpoint(str(tmp_path), "yeast", "base", {"v": "old base"})
        bn.save_checkpoint(str(tmp_path), "yeast", "chip", {"v": "old chip, no brenda"})

        base = bn.step(str(tmp_path), "yeast", "base", lambda: calls.append("base") or {"v": "new"},
                       chain=chain)
        brenda = bn.step(str(tmp_path), "yeast", "brenda",
                         lambda: calls.append("brenda") or {"v": "with brenda"}, chain=chain)
        chip = bn.step(str(tmp_path), "yeast", "chip",
                       lambda: calls.append("chip") or {"v": "chip on top of brenda"}, chain=chain)

        assert base == {"v": "old base"}                 # still loaded: nothing before it ran
        assert calls == ["brenda", "chip"]               # chip rebuilt, not loaded
        assert chip == {"v": "chip on top of brenda"}

    def test_a_fully_checkpointed_chain_still_resumes(self, tmp_path):
        import build_networks as bn

        chain = {"rebuilt": False}
        for name in ("base", "chip"):
            bn.save_checkpoint(str(tmp_path), "mouse", name, {"v": name})
        for name in ("base", "chip"):
            assert bn.step(str(tmp_path), "mouse", name, lambda: 1 / 0, chain=chain) == {"v": name}
