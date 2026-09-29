"""The organism registry: what defines a network that can be built.

The registry used to live inside ``maintenance/build_networks.py``, which meant
adding an organism required editing a script in the repository and nothing
checked that an entry was complete. These tests cover the contract the builder
relies on.
"""

import importlib
import importlib.util
from pathlib import Path

import pytest

from transnet.organisms import (
    ORGANISMS,
    REQUIRED_KEYS,
    available_organisms,
    organism_config,
    register_organism,
)

ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture(autouse=True)
def _restore_registry():
    """Tests that register an organism must not leak it into the others."""
    original = {name: dict(config) for name, config in ORGANISMS.items()}
    yield
    ORGANISMS.clear()
    ORGANISMS.update(original)


class TestTheBuiltOrganisms:

    def test_all_five_are_present(self):
        assert available_organisms() == ["ecoli", "human", "mouse", "rat", "yeast"]

    @pytest.mark.parametrize("name", ["human", "mouse", "rat", "yeast", "ecoli"])
    def test_every_entry_carries_the_required_keys(self, name):
        missing = [key for key in REQUIRED_KEYS if key not in ORGANISMS[name]]
        assert not missing, f"{name} is missing {missing}"

    def test_yeast_overrides_the_uniprot_taxon(self):
        """UniProt files reviewed yeast entries under the strain, not the
        species: 6,733 proteins against 43. Losing this override builds a
        network that looks complete and is not."""
        yeast = ORGANISMS["yeast"]
        assert yeast["uniprot_org"] == "559292"
        assert yeast["uniprot_org"] != yeast["ncbi_org"]

    def test_ecoli_declares_what_cannot_exist(self):
        """ChIP-Atlas carries no bacterial genome, so the absence of those edges
        is a fact about the databases rather than a failed download."""
        ecoli = ORGANISMS["ecoli"]
        assert ecoli["transcriptome_source"] == "kegg"
        assert "transcriptional_regulation" in ecoli["expected_absent"]

    def test_rat_is_pinned_to_the_assembly_chip_atlas_has(self):
        assert ORGANISMS["rat"]["genome_chip"] == "rn6"


class TestRegistering:

    def test_a_new_organism_can_be_added(self):
        register_organism(
            "zebrafish", kegg_org="dre", organism_full="Danio rerio",
            ncbi_org="7955", ensembl_org="danio_rerio", ensembl_release=109,
            genome_chip="danRer11",
        )
        assert organism_config("zebrafish")["kegg_org"] == "dre"
        assert "zebrafish" in available_organisms()

    def test_it_will_not_silently_replace_a_built_organism(self):
        """A typo colliding with a built organism would redefine the network a
        study depends on."""
        with pytest.raises(ValueError, match="already registered"):
            register_organism(
                "rat", kegg_org="xxx", organism_full="Not a rat", ncbi_org="1",
                ensembl_org="x", ensembl_release=109, genome_chip="x",
            )

    def test_replacing_is_possible_when_asked_for(self):
        register_organism(
            "rat", replace=True, kegg_org="rno", organism_full="Rattus norvegicus",
            ncbi_org="10116", ensembl_org="rattus_norvegicus",
            ensembl_release=110, genome_chip="rn6",
        )
        assert ORGANISMS["rat"]["ensembl_release"] == 110

    def test_a_missing_key_is_named(self):
        with pytest.raises(ValueError, match="genome_chip"):
            register_organism(
                "incomplete", kegg_org="xxx", organism_full="X", ncbi_org="1",
                ensembl_org="x", ensembl_release=109,
            )

    def test_an_unknown_organism_lists_the_alternatives(self):
        with pytest.raises(KeyError, match="Configured:"):
            organism_config("platypus")


def test_the_builder_reads_the_same_registry():
    """``build_networks.py`` must resolve every organism the registry holds,
    since that is what ``--organisms`` accepts."""
    spec = importlib.util.spec_from_file_location(
        "build_networks", ROOT / "maintenance" / "build_networks.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.ORGANISM_CONFIG is ORGANISMS
