"""The MoTrPAC post-translational modification loaders.

The raw ``.rda`` objects are a gitignored cache, so the tests that need them are
marked ``network`` and skip when the cache is empty. The rest check the parts
that hold regardless: which tissue has which assay, what object names a study
needs, and the site/protein split that mapping depends on.
"""

import pandas as pd
import pytest

from transnet.datasets import (
    MOTRPAC_PTM_ASSAYS,
    MOTRPAC_RAW,
    MOTRPAC_TISSUES,
    load_motrpac_ptm,
    motrpac_objects,
)


class TestAvailability:

    def test_every_ptm_tissue_is_a_known_tissue(self):
        for assay, tissues in MOTRPAC_PTM_ASSAYS.items():
            unknown = set(tissues) - set(MOTRPAC_TISSUES)
            assert not unknown, f"{assay} names tissues we cannot load: {unknown}"

    def test_acetylation_and_ubiquitination_are_heart_and_liver_only(self):
        """Stated here because a study asking for them elsewhere gets an error,
        and the error should be the documented behaviour rather than a surprise."""
        assert MOTRPAC_PTM_ASSAYS["ACETYL"] == ("HEART", "LIVER")
        assert MOTRPAC_PTM_ASSAYS["UBIQ"] == ("HEART", "LIVER")

    def test_asking_for_an_assay_a_tissue_lacks_is_an_error(self):
        with pytest.raises(ValueError, match="not measured in CORTEX"):
            load_motrpac_ptm("CORTEX", "ACETYL")

    def test_an_unknown_assay_is_an_error(self):
        with pytest.raises(ValueError, match="assay must be one of"):
            load_motrpac_ptm("HEART", "METHYL")


class TestObjectNames:

    def test_the_design_and_metabolome_are_always_included(self):
        names = motrpac_objects(["HEART"], ["PROT"])
        assert "PHENO" in names and "METAB_NORM_DATA_FLAT" in names

    def test_a_tissue_without_an_assay_is_skipped_rather_than_named(self):
        names = motrpac_objects(["CORTEX", "HEART"], ["ACETYL"])
        assert "ACETYL_HEART_NORM_DATA" in names
        assert not any("CORTEX" in name for name in names)

    def test_the_file_label_differs_from_the_tissue_name(self):
        """SKM_GN is SKMGN in the filenames and SKM-GN inside the data; getting
        this wrong silently loads nothing."""
        names = motrpac_objects(["SKM_GN"], ["PHOSPHO"])
        assert "PHOSPHO_SKMGN_NORM_DATA" in names


@pytest.mark.network
class TestAgainstTheRealData:
    """Needs the .rda cache; fetch it with fetch_motrpac_raw()."""

    @pytest.fixture(autouse=True)
    def _require_cache(self):
        if not (MOTRPAC_RAW / "PHOSPHO_HEART_NORM_DATA.rda").exists():
            pytest.skip("MoTrPAC phospho cache absent")

    def test_the_contrast_has_the_documented_columns(self):
        table = load_motrpac_ptm("HEART", "PHOSPHO")
        assert list(table.columns) == [
            "feature", "protein", "site", "log2FC", "pvalue", "se", "padj",
        ]

    def test_the_site_is_split_off_the_refseq_accession(self):
        table = load_motrpac_ptm("HEART", "PHOSPHO")
        row = table.iloc[0]
        assert row["feature"] == f"{row['protein']}_{row['site']}"
        assert row["protein"].startswith(("NP_", "XP_", "YP_"))

    def test_a_protein_can_carry_several_sites(self):
        table = load_motrpac_ptm("HEART", "PHOSPHO")
        assert table["protein"].nunique() < len(table)

    def test_padj_is_a_probability(self):
        table = load_motrpac_ptm("HEART", "ACETYL")
        assert table["padj"].between(0, 1).all()


class TestAttachingSitesToTheNetwork:
    """``map_modification_sites`` turns a site table into Signaling nodes.

    A site is a regulator of its protein's activity, not a measure of its
    abundance, so each one becomes its own node with an edge into the protein,
    and two sites on one protein stay separate because they can move opposite
    ways.
    """

    @staticmethod
    def _graph():
        import networkx as nx

        graph = nx.MultiDiGraph()
        graph.add_node("P1", layer="Proteome", node_type="Protein",
                       name="Enzyme one", symbol="Eno1")
        graph.add_node("P2", layer="Proteome", node_type="Protein",
                       name="Enzyme two", symbol="Hk2")
        return graph

    @staticmethod
    def _table():
        return pd.DataFrame({
            "protein": ["NP_1", "NP_1", "NP_2", "NP_absent"],
            "site": ["S10s", "T20t", "S5s", "S1s"],
            "log2FC": [1.5, -1.2, 0.1, 3.0],
            "padj": [0.001, 0.01, 0.9, 0.001],
        })

    def test_sites_become_signaling_nodes_wired_to_their_protein(self):
        from transnet import map_modification_sites

        graph = self._graph()
        report = map_modification_sites(
            graph, self._table(), id_map={"NP_1": "P1", "NP_2": "P2"})

        assert report["n_mapped"] == 3 and report["n_proteins"] == 2
        assert graph.nodes["P1_S10s"]["layer"] == "Signaling"
        assert graph.has_edge("P1_S10s", "P1")

    def test_two_sites_on_one_protein_stay_separate_and_can_disagree(self):
        from transnet import map_modification_sites

        graph = self._graph()
        map_modification_sites(graph, self._table(),
                               id_map={"NP_1": "P1", "NP_2": "P2"})
        assert graph.nodes["P1_S10s"]["regulated"] == 1
        assert graph.nodes["P1_T20t"]["regulated"] == -1

    def test_an_unchanged_site_is_attached_but_not_regulated(self):
        from transnet import map_modification_sites

        graph = self._graph()
        map_modification_sites(graph, self._table(),
                               id_map={"NP_1": "P1", "NP_2": "P2"})
        assert graph.nodes["P2_S5s"]["regulated"] == 0

    def test_a_protein_absent_from_the_network_is_skipped(self):
        from transnet import map_modification_sites

        graph = self._graph()
        report = map_modification_sites(graph, self._table(),
                                        id_map={"NP_1": "P1", "NP_2": "P2"})
        assert report["n_sites"] == 4 and report["n_mapped"] == 3
        assert not any("absent" in str(node) for node in graph)

    def test_the_edges_are_unsigned_by_default(self):
        """The effect of a site on activity is not generally known, and an
        invented sign would propagate into every downstream direction."""
        from transnet import map_modification_sites

        graph = self._graph()
        map_modification_sites(graph, self._table(), id_map={"NP_1": "P1"})
        signs = {d["sign"] for _, _, d in graph.edges(data=True)}
        assert signs == {0}

    def test_the_sites_feed_the_phospho_axis(self):
        """The point of attaching them: a reaction whose enzyme amount held
        steady while its modification state moved becomes visible."""
        import networkx as nx

        from transnet import map_modification_sites, reaction_regulation_table

        graph = self._graph()
        graph.add_node("R1", layer="Reactions", node_type="Reaction",
                       name="a reaction", reversible=False)
        graph.add_edge("P1", "R1", key="catalysis", edge_type="catalysis", sign=1)
        graph.nodes["P1"].update(measured=True, regulated=0, log2fc=0.0)

        map_modification_sites(graph, self._table(), id_map={"NP_1": "P1"})
        row = reaction_regulation_table(graph).iloc[0]
        assert row["gene_axis"] == 0
        # Two sites on P1 moved opposite ways, so no single direction is honest:
        # the count is what says the enzyme is phospho-regulated at all.
        assert row["phospho_axis"] == 0
        assert row["n_phosphosites_changed"] == 2
        assert "P1_S10s" in row["phospho_axis_via"]

    def test_no_proteome_is_reported_rather_than_crashing(self):
        import networkx as nx

        from transnet import map_modification_sites

        report = map_modification_sites(nx.MultiDiGraph(), self._table())
        assert report == {"n_sites": 4, "n_mapped": 0, "n_changed": 0,
                          "n_proteins": 0}
