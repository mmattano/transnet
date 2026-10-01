"""Identifier handling in the database-enrichment steps.

Each of these encodes a bug that produced an *empty result with no error*: the
build finished, reported success, and shipped a network missing a whole class
of edge. They are cheap unit tests precisely because the failures were silent.
"""

import pytest

from transnet.biology.elements import Protein
from transnet.biology.layers import Proteome


def _proteome(*proteins):
    layer = Proteome()
    layer.proteins = list(proteins)
    layer.ncbi_organism = "10090"
    return layer


class TestChipAtlasSymbolExtraction:
    """`protein.gene` is a list; `str()` on it yields "['Actb']"."""

    def test_symbols_are_flattened_not_stringified(self, monkeypatch):
        captured = {}

        def fake_get_chip_data(genome, distance, proteins=None):
            captured["proteins"] = proteins
            raise RuntimeError("stop here; we only care about the argument")

        monkeypatch.setattr(
            "transnet.biology.layers.get_ChIP_data", fake_get_chip_data
        )
        layer = _proteome(
            Protein(uniprot_id="P1", gene=["Actb"]),
            Protein(uniprot_id="P2", gene=["Myc", "Bhlhe39"]),
        )
        layer.get_transcription_factor_targets(genome_ChIP="mm10")

        assert captured["proteins"] == ["Actb", "Bhlhe39", "Myc"], (
            "gene symbols must be flattened, not str()'d as lists"
        )
        assert not any("[" in symbol for symbol in captured["proteins"])

    def test_a_missing_gene_field_is_skipped(self, monkeypatch):
        captured = {}

        def fake_get_chip_data(genome, distance, proteins=None):
            captured["proteins"] = proteins
            raise RuntimeError("stop")

        monkeypatch.setattr(
            "transnet.biology.layers.get_ChIP_data", fake_get_chip_data
        )
        layer = _proteome(
            Protein(uniprot_id="P1", gene=["Actb"]),
            Protein(uniprot_id="P2", gene=None),
            Protein(uniprot_id="P3", gene=float("nan")),
        )
        layer.get_transcription_factor_targets(genome_ChIP="mm10")
        assert captured["proteins"] == ["Actb"]


class TestChipAtlasExperimentColumn:
    """`get_ChIP_exps` returns named columns; `exp_info[0]` raised KeyError: 0.

    The experiment list is only consulted when a cell-type restriction is
    asked for -- without one it selects every experiment for the genome,
    which is what `get_ChIP_data` already returned -- so these tests pass a
    cell type to reach that code.
    """

    def test_named_experiment_column_is_used(self, monkeypatch):
        import pandas as pd

        # Columns are bare experiment accessions, as ChIP-Atlas returns them.
        scores = pd.DataFrame({"SRX1726750": [500.0]}, index=["TargetGene"])

        monkeypatch.setattr(
            "transnet.biology.layers.get_ChIP_data",
            lambda genome, distance, proteins=None: (
                scores, {"Actb": ["SRX1726750"]}
            ),
        )
        monkeypatch.setattr(
            "transnet.biology.layers.get_ChIP_exps",
            lambda **kwargs: pd.DataFrame({
                "exp_id": ["SRX1726750"], "genome": ["mm10"],
            }),
        )

        layer = _proteome(Protein(uniprot_id="P1", gene=["Actb"]))
        layer.get_transcription_factor_targets(
            genome_ChIP="mm10", cell_type_ChIP="Liver"
        )

        assert layer.proteins[0].transcription_factor_targets, (
            "an integer column index silently discarded every TF target"
        )

    def test_the_experiment_list_is_not_fetched_without_a_cell_filter(
        self, monkeypatch
    ):
        """Fetching it downloads the full cross-genome experimentList.tab
        (100k+ rows) to perform a no-op filter."""
        import pandas as pd

        scores = pd.DataFrame({"SRX1726750": [500.0]}, index=["TargetGene"])
        monkeypatch.setattr(
            "transnet.biology.layers.get_ChIP_data",
            lambda genome, distance, proteins=None: (
                scores, {"Actb": ["SRX1726750"]}
            ),
        )

        def should_not_run(**kwargs):
            raise AssertionError(
                "get_ChIP_exps was called without a cell-type restriction"
            )

        monkeypatch.setattr(
            "transnet.biology.layers.get_ChIP_exps", should_not_run
        )

        layer = _proteome(Protein(uniprot_id="P1", gene=["Actb"]))
        assert layer.get_transcription_factor_targets(genome_ChIP="mm10")
        assert layer.proteins[0].transcription_factor_targets == ["TargetGene"]


class TestChipAtlasScoreThreshold:
    """ChIP-Atlas lists every gene with any detectable binding.

    Taking the union of everything non-zero across a factor's experiments
    therefore keeps essentially the whole file -- ~15,000 targets per factor,
    which made transcriptional regulation 91% of the mouse network. Targets
    are scored by their mean across experiments instead, which is exactly
    ChIP-Atlas's own "{TF}|Average" consensus score.
    """

    @staticmethod
    def _layer(monkeypatch, scores, gene_to_experiment):
        monkeypatch.setattr(
            "transnet.biology.layers.get_ChIP_data",
            lambda genome, distance, proteins=None: (
                scores, gene_to_experiment
            ),
        )
        return _proteome(Protein(uniprot_id="P1", gene=["Actb"]))

    def test_a_single_strong_peak_does_not_make_a_target(self, monkeypatch):
        import pandas as pd

        # Bound at 900 in one of ten experiments: mean 90, below the default.
        scores = pd.DataFrame(
            [[900.0] + [0.0] * 9, [200.0] * 10],
            index=["OneOffPeak", "Consistent"],
            columns=[f"SRX{i}" for i in range(10)],
        )
        layer = self._layer(
            monkeypatch, scores, {"Actb": list(scores.columns)}
        )
        layer.get_transcription_factor_targets(genome_ChIP="mm10")

        assert layer.proteins[0].transcription_factor_targets == ["Consistent"]

    def test_threshold_zero_restores_the_permissive_behaviour(self, monkeypatch):
        import pandas as pd

        scores = pd.DataFrame(
            [[900.0] + [0.0] * 9, [200.0] * 10],
            index=["OneOffPeak", "Consistent"],
            columns=[f"SRX{i}" for i in range(10)],
        )
        layer = self._layer(
            monkeypatch, scores, {"Actb": list(scores.columns)}
        )
        layer.get_transcription_factor_targets(
            genome_ChIP="mm10", score_threshold=0
        )

        assert set(layer.proteins[0].transcription_factor_targets) == {
            "OneOffPeak", "Consistent"
        }

    def test_scores_are_recorded_and_targets_ranked(self, monkeypatch):
        import pandas as pd

        scores = pd.DataFrame(
            [[300.0, 300.0], [800.0, 800.0]],
            index=["Weaker", "Stronger"],
            columns=["SRX0", "SRX1"],
        )
        layer = self._layer(
            monkeypatch, scores, {"Actb": ["SRX0", "SRX1"]}
        )
        layer.get_transcription_factor_targets(genome_ChIP="mm10")

        protein = layer.proteins[0]
        assert protein.transcription_factor_targets == ["Stronger", "Weaker"]
        assert protein.transcription_factor_target_scores == {
            "Stronger": 800.0, "Weaker": 300.0
        }

    def test_max_targets_per_tf_caps_promiscuous_factors(self, monkeypatch):
        import pandas as pd

        scores = pd.DataFrame(
            {"SRX0": [500.0, 400.0, 300.0, 200.0]},
            index=["A", "B", "C", "D"],
        )
        layer = self._layer(monkeypatch, scores, {"Actb": ["SRX0"]})
        layer.get_transcription_factor_targets(
            genome_ChIP="mm10", max_targets_per_tf=2
        )

        assert layer.proteins[0].transcription_factor_targets == ["A", "B"]


class TestStringFailureIsReported:
    """A timeout must be distinguishable from 'ran, found nothing'."""

    def test_a_failed_lookup_returns_false(self, monkeypatch):
        def boom(*args, **kwargs):
            raise TimeoutError("read timed out")

        monkeypatch.setattr(
            "transnet.biology.layers.string_map_identifiers", boom
        )
        layer = _proteome(Protein(uniprot_id="P1", gene=["Actb"]))
        assert layer.get_interaction_partners() is False

    def test_a_successful_lookup_returns_true(self, monkeypatch):
        monkeypatch.setattr(
            "transnet.biology.layers.string_map_identifiers",
            lambda protein_list, species: {"P1": "10090.ENSMUSP1"},
        )
        monkeypatch.setattr(
            "transnet.biology.layers.string_get_interactions",
            lambda protein_list, species, cutoff_score: {
                "10090.ENSMUSP1": ["10090.ENSMUSP2"]
            },
        )
        monkeypatch.setattr(
            "transnet.biology.layers.translate_string_dict",
            lambda mapping_dict, interactions_dict: {"P1": ["P2"]},
        )
        layer = _proteome(Protein(uniprot_id="P1", gene=["Actb"]))
        assert layer.get_interaction_partners() is True
        assert layer.proteins[0].interaction_partners == ["P2"]


class TestStringSymbolFallback:
    """STRING resolves UniProt for mouse but not yeast; fall back to symbols."""

    def test_symbols_are_tried_when_accessions_fail(self, monkeypatch):
        seen = []

        def fake_map(protein_list, species):
            seen.append(list(protein_list))
            # First call (UniProt accessions) resolves nothing.
            if seen[0] == protein_list and len(seen) == 1:
                return {}
            return {"PGI1": "4932.YBR196C"}

        monkeypatch.setattr("transnet.biology.layers.string_map_identifiers", fake_map)
        monkeypatch.setattr(
            "transnet.biology.layers.string_get_interactions",
            lambda protein_list, species, cutoff_score: {"4932.YBR196C": []},
        )
        monkeypatch.setattr(
            "transnet.biology.layers.translate_string_dict",
            lambda mapping_dict, interactions_dict: {},
        )

        layer = _proteome(Protein(uniprot_id="P12345", gene=["PGI1"]))
        layer.ncbi_organism = "4932"
        layer.get_interaction_partners()

        assert len(seen) == 2, "no gene-symbol retry after UniProt found nothing"
        assert seen[1] == ["PGI1"]

    def test_a_nan_gene_field_does_not_crash_the_fallback(self, monkeypatch):
        monkeypatch.setattr(
            "transnet.biology.layers.string_map_identifiers",
            lambda protein_list, species: {},
        )
        layer = _proteome(Protein(uniprot_id="P1", gene=float("nan")))
        layer.ncbi_organism = "4932"
        assert layer.get_interaction_partners() is None or True   # must not raise


class TestChipAtlasRetries:
    """A timeout used to drop a factor silently.

    `get_ChIP_data` logged a warning and moved on, so a factor lost to a
    transient network error was indistinguishable in the finished network
    from one that genuinely has no targets -- and the build step was still
    checkpointed as complete.
    """

    def test_a_timeout_is_retried_and_then_succeeds(self, monkeypatch):
        import requests
        from transnet.api import chip_atlas

        calls = {"n": 0}

        class Response:
            status_code = 200
            text = "ok"

        def flaky(url, timeout=None):
            calls["n"] += 1
            if calls["n"] < 3:
                raise requests.Timeout("read timed out")
            return Response()

        monkeypatch.setattr(chip_atlas.requests, "get", flaky)
        # don't actually sleep between retries
        import time as _t
        monkeypatch.setattr(_t, "sleep", lambda *_: None)

        response = chip_atlas._get_with_retry("http://x", label="Sox2")
        assert response.status_code == 200
        assert calls["n"] == 3, "should have retried twice before succeeding"

    def test_exhausted_retries_raise_rather_than_return_nothing(self, monkeypatch):
        import requests
        import time as _t
        from transnet.api import chip_atlas

        monkeypatch.setattr(
            chip_atlas.requests, "get",
            lambda url, timeout=None: (_ for _ in ()).throw(
                requests.Timeout("read timed out")
            ),
        )
        monkeypatch.setattr(_t, "sleep", lambda *_: None)

        with pytest.raises(requests.RequestException, match="after 3 attempts"):
            chip_atlas._get_with_retry("http://x", label="Sox2")

    def test_a_404_is_not_retried(self, monkeypatch):
        from transnet.api import chip_atlas

        calls = {"n": 0}

        class NotFound:
            status_code = 404
            text = ""

        def counting(url, timeout=None):
            calls["n"] += 1
            return NotFound()

        monkeypatch.setattr(chip_atlas.requests, "get", counting)
        response = chip_atlas._get_with_retry("http://x", label="Nope")
        assert response.status_code == 404
        assert calls["n"] == 1, "404 means no such factor; retrying is pointless"

    def test_failed_factors_are_recorded_on_the_layer(self, monkeypatch):
        import pandas as pd

        scores = pd.DataFrame({"SRX0": [500.0]}, index=["TargetGene"])
        monkeypatch.setattr(
            "transnet.biology.layers.get_ChIP_data",
            lambda genome, distance, proteins=None: (
                scores, {"Actb": ["SRX0"]}, ["Sox2", "Myc"]
            ),
        )
        layer = _proteome(Protein(uniprot_id="P1", gene=["Actb"]))
        assert layer.get_transcription_factor_targets(genome_ChIP="mm10")
        assert layer.chip_atlas_failed_factors == ["Sox2", "Myc"]

    def test_the_two_value_return_is_still_accepted(self, monkeypatch):
        """Older callers and test doubles return (scores, gene_to_experiment)."""
        import pandas as pd

        scores = pd.DataFrame({"SRX0": [500.0]}, index=["TargetGene"])
        monkeypatch.setattr(
            "transnet.biology.layers.get_ChIP_data",
            lambda genome, distance, proteins=None: (scores, {"Actb": ["SRX0"]}),
        )
        layer = _proteome(Protein(uniprot_id="P1", gene=["Actb"]))
        assert layer.get_transcription_factor_targets(genome_ChIP="mm10")
        assert layer.chip_atlas_failed_factors == []
        assert layer.proteins[0].transcription_factor_targets == ["TargetGene"]


class TestPubChemToKeggUsesEveryChebiId:
    """A compound's neutral form and zwitterion are separate ChEBI entries.

    KEGG cross-references only some of them, so keeping just the first id
    mychem returns silently dropped ~20% of a metabolite panel: 106 of 171
    mouse metabolites resolved instead of 144.
    """

    def test_a_later_chebi_id_still_resolves(self, monkeypatch):
        import pandas as pd
        from transnet.biology.layers import Metabolome

        # CID 5951: mychem returns the zwitterion first; only the neutral
        # form (the second entry) is cross-referenced by KEGG.
        monkeypatch.setattr(
            "transnet.biology.layers.chemical_info_converter",
            lambda ids: pd.DataFrame({
                "chebi_ids": ["CHEBI:57305"],
                "chebi_ids_all": [["CHEBI:57305", "CHEBI:16811"]],
            }),
        )
        monkeypatch.setattr(
            "transnet.biology.layers.chebi_to_kegg",
            lambda chebis: pd.DataFrame({
                "chebi_compounds": ["CHEBI:57305", "CHEBI:16811"],
                "kegg_compounds": [None, "C00065"],
            }),
        )
        monkeypatch.setattr(
            "transnet.biology.layers.kegg_list_compounds",
            lambda: pd.DataFrame({
                "kegg_compounds": ["C00065"], "name": ["L-Serine"],
            }),
        )

        layer = Metabolome()
        layer.add_experimental_data(
            pd.DataFrame({"metabolite_id": ["5951"]}),
            metabolite_column_name="metabolite_id", id_type="pubchem",
        )
        layer.populate()

        resolved = [m for m in layer.metabolites if m.pubchem_id]
        assert len(resolved) == 1, (
            "only the first ChEBI id was tried, so the compound was lost"
        )
        assert resolved[0].kegg_compound_id == "C00065"
        assert resolved[0].pubchem_id == "5951"

    def test_the_old_single_id_frame_still_works(self, monkeypatch):
        """chemical_info_converter without the chebi_ids_all column."""
        import pandas as pd
        from transnet.biology.layers import Metabolome

        monkeypatch.setattr(
            "transnet.biology.layers.chemical_info_converter",
            lambda ids: pd.DataFrame({"chebi_ids": ["CHEBI:16811"]}),
        )
        monkeypatch.setattr(
            "transnet.biology.layers.chebi_to_kegg",
            lambda chebis: pd.DataFrame({
                "chebi_compounds": ["CHEBI:16811"], "kegg_compounds": ["C00065"],
            }),
        )
        monkeypatch.setattr(
            "transnet.biology.layers.kegg_list_compounds",
            lambda: pd.DataFrame({
                "kegg_compounds": ["C00065"], "name": ["L-Serine"],
            }),
        )

        layer = Metabolome()
        layer.add_experimental_data(
            pd.DataFrame({"metabolite_id": ["5951"]}),
            metabolite_column_name="metabolite_id", id_type="pubchem",
        )
        layer.populate()
        assert [m.pubchem_id for m in layer.metabolites if m.pubchem_id] == ["5951"]


class TestKeggTranscriptomeCarriesNcbiIds:
    """A gene with no NCBI id cannot be joined to its protein.

    `Transcriptome.populate(kegg_api=True)` built genes from the KEGG gene
    list without resolving NCBI ids, and `fill_gene_info()` -- which would
    have backfilled them -- is not called by populate. The Transcriptome
    layer was then present and wired downward by transcriptional_regulation,
    but had no `translation` edges up to the Proteome at all: on rat, 0 of
    26,431 genes carried an ncbi_id.
    """

    @staticmethod
    def _patch(monkeypatch, idtable):
        import pandas as pd

        monkeypatch.setattr(
            "transnet.biology.layers.kegg_list_genes",
            lambda org: pd.DataFrame({
                "gene_id": ["24152", "134485287"],
                "name(s)": ["Acly", "Zfp431"],
                "description": ["ATP citrate lyase", "zinc finger protein"],
            }),
        )
        monkeypatch.setattr(
            "transnet.biology.layers.kegg_conv_ncbi_idtable", lambda org: idtable
        )

    def test_ncbi_ids_are_resolved_from_the_kegg_id_table(self, monkeypatch):
        import pandas as pd
        from transnet.biology.layers import Transcriptome

        self._patch(monkeypatch, pd.DataFrame({
            "ncbi_id": ["24152", "134485287"],
            "kegg_id": ["24152", "134485287"],
        }))

        layer = Transcriptome()
        layer.populate(kegg_organism="rno", kegg_api=True)

        assert len(layer.genes) == 2
        assert [g.ncbi_id for g in layer.genes] == ["24152", "134485287"], (
            "genes built from the KEGG list carried no NCBI id, so no "
            "translation edge could ever be built"
        )

    def test_a_failed_conversion_leaves_genes_but_warns(self, monkeypatch):
        """The layer is still usable; only the join is lost."""
        from transnet.biology.layers import Transcriptome

        def boom(org):
            raise RuntimeError("KEGG unreachable")

        self._patch(monkeypatch, None)
        monkeypatch.setattr(
            "transnet.biology.layers.kegg_conv_ncbi_idtable", boom
        )

        layer = Transcriptome()
        layer.populate(kegg_organism="rno", kegg_api=True)
        assert len(layer.genes) == 2
        assert all(g.ncbi_id is None for g in layer.genes)

    def test_translation_edges_follow_from_the_ncbi_id(self, monkeypatch):
        """The end the bug was actually felt at."""
        import pandas as pd
        from transnet.biology.layers import Transcriptome, Proteome
        from transnet.biology.elements import Protein
        from transnet import Transnet

        self._patch(monkeypatch, pd.DataFrame({
            "ncbi_id": ["24152", "134485287"],
            "kegg_id": ["24152", "134485287"],
        }))

        transcriptome = Transcriptome()
        transcriptome.populate(kegg_organism="rno", kegg_api=True)

        proteome = Proteome()
        proteome.proteins = [Protein(uniprot_id="P16638", entrez_id=["24152"])]

        network = Transnet(transcriptome=transcriptome, proteome=proteome)
        edges = network.gene_protein_interaction()

        assert [(e["source"], e["target"]) for e in edges] == [("24152", "P16638")]


def test_the_chip_atlas_experiment_list_is_downloaded_once(tmp_path, monkeypatch):
    """It is 200 MB; three helpers read it, and each used to fetch it again.
    The columns follow the real file: id, genome, antigen class, antigen, cell
    type class, cell type."""
    from transnet.api import chip_atlas

    calls = []

    class Response:
        text = "SRX1\tmm10\tTFs and others\tFoxo1\tLiver\tHepatocytes\n"

        def raise_for_status(self):
            pass

    def fake_get(*args, **kwargs):
        calls.append(args)
        return Response()

    monkeypatch.setenv("TRANSNET_CHIP_CACHE", str(tmp_path))
    monkeypatch.setattr(chip_atlas, "_get_with_retry", fake_get)
    assert chip_atlas.list_ChIP_cell_type_classes("mm10") == ["Liver"]
    assert chip_atlas.list_ChIP_cell_types("mm10", "Liver") == ["Hepatocytes"]
    assert len(calls) == 1
