"""reaction regulation axes/per-pathway balance: reaction regulation-axis attribution and per-pathway summary."""

import networkx as nx
import pandas as pd
import pytest

from transnet.analysis.transomics import (
    map_omics_to_network,
    reaction_regulation_table,
    regulation_axis_summary,
)


@pytest.fixture
def mapped_network(standard_network, omics_tables):
    map_omics_to_network(
        standard_network,
        {k: v for k, v in omics_tables.items() if k != "Signaling"},
        id_column="id", log2fc_column="log2FC", qvalue_column="padj",
    )
    return standard_network


class TestAxisAttribution:
    def test_hexokinase_axes_disagree(self, mapped_network):
        # HK1 protein is up (gene axis activating) while its product G6P is up
        # and inhibits it (metabolite axis inhibiting).  That is the textbook
        # controversial reaction.
        table = reaction_regulation_table(mapped_network)
        row = table.set_index("reaction").loc["R00299"]

        assert row["gene_axis"] == 1
        assert row["metabolite_axis"] == -1
        assert bool(row["controversial"])
        assert row["consensus"] == 0

    def test_protein_evidence_is_preferred_over_transcript(self, mapped_network):
        table = reaction_regulation_table(mapped_network)
        row = table.set_index("reaction").loc["R00299"]
        assert row["gene_axis_evidence"] == "protein"
        assert "P19367" in row["gene_axis_via"]

    def test_allosteric_regulator_is_named_with_its_effect(self, mapped_network):
        table = reaction_regulation_table(mapped_network)
        row = table.set_index("reaction").loc["R00299"]
        assert "C00668" in row["allosteric_regulators"]
        assert "inhibited" in row["allosteric_regulators"]

    def test_activating_allostery_gives_a_positive_metabolite_axis(self, mapped_network):
        # F1,6BP is up and activates PFK.
        table = reaction_regulation_table(mapped_network)
        row = table.set_index("reaction").loc["R00756"]
        assert "C00354" in row["allosteric_regulators"]
        assert row["n_activators"] >= 1

    def test_agreeing_axes_produce_a_consensus_and_no_controversy(self, standard_network):
        # HK1 protein down, and its substrate glucose down too: less enzyme and
        # less substrate both push the reaction the same way.
        map_omics_to_network(
            standard_network,
            {
                "Proteome": pd.DataFrame({
                    "id": ["P19367"], "log2FC": [-1.0], "padj": [0.01]}),
                "Metabolome": pd.DataFrame({
                    "id": ["C00031"], "log2FC": [-1.5], "padj": [0.01]}),
            },
            id_column="id", log2fc_column="log2FC", qvalue_column="padj",
        )
        row = reaction_regulation_table(standard_network).set_index("reaction").loc["R00299"]
        assert row["gene_axis"] == -1
        assert row["metabolite_axis"] == -1
        assert not bool(row["controversial"])
        assert row["consensus"] == -1

    def test_a_metabolite_with_two_opposing_roles_cancels_out(self, standard_network):
        # F1,6BP is both the product of PFK and its allosteric activator. When
        # it falls, "less product" pushes the reaction forward while "less
        # activator" pushes it back. The metabolite axis correctly reports no
        # net direction rather than picking one arbitrarily -- a distinction the
        # old flat, unsigned graph could not even represent.
        map_omics_to_network(
            standard_network,
            {"Metabolome": pd.DataFrame({
                "id": ["C00354"], "log2FC": [-1.5], "padj": [0.01]})},
            id_column="id", log2fc_column="log2FC", qvalue_column="padj",
        )
        row = reaction_regulation_table(standard_network).set_index("reaction").loc["R00756"]
        assert row["n_activators"] == 1
        assert row["n_inhibitors"] == 1
        assert row["metabolite_axis"] == 0

    def test_unmapped_network_yields_no_regulated_reactions(self, standard_network):
        table = reaction_regulation_table(standard_network)
        assert len(table) == 2
        assert (table["gene_axis"] == 0).all()
        assert (table["metabolite_axis"] == 0).all()
        assert not table["controversial"].any()

    def test_allosteric_only_mode_ignores_substrate_changes(self, mapped_network):
        full = reaction_regulation_table(mapped_network)
        allosteric = reaction_regulation_table(
            mapped_network, include_substrate_product=False
        )
        row = allosteric.set_index("reaction").loc["R00299"]
        assert row["substrates_changed"] == ""
        # glucose (a changed substrate) counted in the full run but not here
        assert (
            full.set_index("reaction").loc["R00299"]["n_inhibitors"]
            >= row["n_inhibitors"]
        )


class TestGeneAxisDegradation:
    """The gene axis must report the evidence it actually had."""

    def test_transcript_only_evidence_is_labelled_gene_protein(self, standard_network):
        map_omics_to_network(
            standard_network,
            {"Transcriptome": pd.DataFrame({
                "id": ["3098"], "log2FC": [1.5], "padj": [0.01]})},
            id_column="id", log2fc_column="log2FC", qvalue_column="padj",
        )
        row = reaction_regulation_table(standard_network).set_index("reaction").loc["R00299"]
        assert row["gene_axis"] == 1
        assert row["gene_axis_evidence"] == "gene_protein"

    def test_no_proteome_falls_back_to_gene_catalysis(self, network_builder):
        graph = network_builder("minimal").generate_graph()
        map_omics_to_network(
            graph,
            {"Transcriptome": pd.DataFrame({
                "id": ["3098"], "log2FC": [1.5], "padj": [0.01]})},
            id_column="id", log2fc_column="log2FC", qvalue_column="padj",
        )
        row = reaction_regulation_table(graph).set_index("reaction").loc["R00299"]
        assert row["gene_axis"] == 1
        assert row["gene_axis_evidence"] == "gene"

    def test_metabolite_axis_works_with_no_gene_evidence_at_all(self, network_builder):
        graph = network_builder("minimal").generate_graph()
        map_omics_to_network(
            graph,
            {"Metabolome": pd.DataFrame({
                "id": ["C00031"], "log2FC": [-1.5], "padj": [0.01]})},
            id_column="id", log2fc_column="log2FC", qvalue_column="padj",
        )
        row = reaction_regulation_table(graph).set_index("reaction").loc["R00299"]
        assert row["gene_axis"] == 0
        assert row["gene_axis_evidence"] is None
        # glucose down -> less substrate -> reaction pushed down
        assert row["metabolite_axis"] == -1

    def test_no_reactions_layer_returns_an_empty_typed_table(self, network_builder):
        graph = network_builder("no_reactions").generate_graph()
        table = reaction_regulation_table(graph)
        assert table.empty
        assert "controversial" in table.columns


class TestAxisSummary:
    def test_summary_counts_controversial_reactions(self, mapped_network):
        table = reaction_regulation_table(mapped_network)
        summary = regulation_axis_summary(table)
        assert len(summary) == 1
        assert summary.iloc[0]["n_controversial"] >= 1
        assert 0.0 <= summary.iloc[0]["fraction_controversial"] <= 1.0

    def test_pathway_map_splits_the_summary(self, mapped_network):
        table = reaction_regulation_table(mapped_network)
        summary = regulation_axis_summary(
            table, pathway_map={"R00299": "glycolysis", "R00756": "glycolysis"}
        )
        assert summary.iloc[0]["pathway"] == "glycolysis"

    def test_a_reaction_in_several_pathways_counts_in_each(self, mapped_network):
        table = reaction_regulation_table(mapped_network)
        regulated = table[(table["gene_axis"] != 0) | (table["metabolite_axis"] != 0)]
        reaction = regulated["reaction"].iloc[0]
        summary = regulation_axis_summary(
            table, pathway_map={reaction: ["glycolysis", "gluconeogenesis"]})
        counts = summary.set_index("pathway")["n_reactions"]
        assert counts["glycolysis"] == counts["gluconeogenesis"] == 1
        assert counts["unassigned"] == len(regulated) - 1

    def test_the_example_pathway_map_covers_every_reaction(self):
        from transnet import load_example_network, load_example_pathways

        reactions = {n for n, d in load_example_network().nodes(data=True)
                     if d.get("layer") == "Reactions"}
        assert reactions <= set(load_example_pathways())

    def test_empty_table_gives_empty_summary_with_columns(self):
        summary = regulation_axis_summary(pd.DataFrame())
        assert summary.empty
        assert "fraction_controversial" in summary.columns


class TestCurrencyMetabolites:
    """Cofactors carry no specificity as substrates, but can be real effectors."""

    def test_a_currency_substrate_is_ignored_by_default(self, standard_network):
        # ATP is a substrate of hexokinase and a currency metabolite.
        map_omics_to_network(
            standard_network,
            {"Metabolome": pd.DataFrame({
                "id": ["C00002"], "log2FC": [2.0], "padj": [0.001]})},
            id_column="id", log2fc_column="log2FC", qvalue_column="padj",
        )
        row = reaction_regulation_table(standard_network).set_index("reaction").loc["R00299"]
        assert row["metabolite_axis"] == 0
        assert row["substrates_changed"] == ""

    def test_the_same_substrate_counts_when_exclusion_is_off(self, standard_network):
        map_omics_to_network(
            standard_network,
            {"Metabolome": pd.DataFrame({
                "id": ["C00002"], "log2FC": [2.0], "padj": [0.001]})},
            id_column="id", log2fc_column="log2FC", qvalue_column="padj",
        )
        row = reaction_regulation_table(
            standard_network, exclude_currency_metabolites=False
        ).set_index("reaction").loc["R00299"]
        assert row["metabolite_axis"] != 0
        assert "C00002" in row["substrates_changed"]

    def test_allosteric_roles_survive_exclusion(self, standard_network):
        # ATP is currency *and* an allosteric inhibitor of PFK. The allosteric
        # role is the finding; suppressing it would defeat the analysis.
        map_omics_to_network(
            standard_network,
            {"Metabolome": pd.DataFrame({
                "id": ["C00002"], "log2FC": [2.0], "padj": [0.001]})},
            id_column="id", log2fc_column="log2FC", qvalue_column="padj",
        )
        row = reaction_regulation_table(standard_network).set_index("reaction").loc["R00756"]
        assert "C00002" in row["allosteric_regulators"]
        assert row["metabolite_axis"] == -1

    def test_a_non_currency_substrate_still_counts(self, standard_network):
        map_omics_to_network(
            standard_network,
            {"Metabolome": pd.DataFrame({
                "id": ["C00031"], "log2FC": [-1.5], "padj": [0.001]})},
            id_column="id", log2fc_column="log2FC", qvalue_column="padj",
        )
        row = reaction_regulation_table(standard_network).set_index("reaction").loc["R00299"]
        assert "C00031" in row["substrates_changed"]
        assert row["metabolite_axis"] == -1

    def test_the_set_can_be_overridden(self, standard_network):
        map_omics_to_network(
            standard_network,
            {"Metabolome": pd.DataFrame({
                "id": ["C00031"], "log2FC": [-1.5], "padj": [0.001]})},
            id_column="id", log2fc_column="log2FC", qvalue_column="padj",
        )
        row = reaction_regulation_table(
            standard_network, currency_metabolites=["C00031"]
        ).set_index("reaction").loc["R00299"]
        assert row["substrates_changed"] == ""


class TestTranscriptSupport:
    """A changed protein is always the strongest evidence, so
    ``gene_axis_evidence`` reads "protein" whether or not the transcript also
    moved. Pklr (protein +2.0, transcript +1.74) was reported as if its
    regulation were post-transcriptional."""

    @staticmethod
    def _network(gene_state, protein_state, gene_measured=True):
        import networkx as nx
        graph = nx.MultiDiGraph()
        graph.add_node("g", layer="Transcriptome")
        if gene_measured:
            graph.nodes["g"].update(measured=True, regulated=gene_state,
                                    log2fc=float(gene_state))
        graph.add_node("p", layer="Proteome", measured=True,
                       regulated=protein_state, log2fc=float(protein_state))
        graph.add_node("R1", layer="Reactions")
        graph.add_edge("g", "p", edge_type="translation", sign=1)
        graph.add_edge("p", "R1", edge_type="catalysis", sign=1)
        return reaction_regulation_table(graph).set_index("reaction").loc["R1"]

    def test_transcript_and_protein_together_are_transcriptional(self):
        row = self._network(1, 1)
        assert row["gene_axis_evidence"] == "protein"
        assert row["gene_axis_transcript_support"] == True  # noqa: E712

    def test_protein_without_its_transcript_is_post_transcriptional(self):
        row = self._network(0, 1)
        support = row["gene_axis_transcript_support"]
        assert support is not None and support == False  # noqa: E712

    def test_no_measured_transcript_is_unknown(self):
        row = self._network(0, 1, gene_measured=False)
        assert pd.isna(row["gene_axis_transcript_support"])


class TestPhosphorylationAxis:
    """Enzyme phosphorylation, reported beside the gene axis rather than in it.

    A phosphosite's effect on catalysis is site-specific and usually unrecorded,
    so the direction the site moved and its effect on the reaction are separate
    columns, and the effect is 0 unless the edge is signed. Merging the two would
    assert a direction the data does not carry.
    """

    @staticmethod
    def _network(site_state=1, sign=0, protein_state=0):
        graph = nx.MultiDiGraph()
        graph.add_node("R1", layer="Reactions", node_type="Reaction",
                       name="a reaction", reversible=False)
        graph.add_node("P1", layer="Proteome", node_type="Protein", name="Enzyme",
                       measured=True, regulated=protein_state,
                       log2fc=float(protein_state))
        graph.add_edge("P1", "R1", key="catalysis", edge_type="catalysis", sign=1)
        if site_state is not None:
            graph.add_node("K1_S100", layer="Signaling", node_type="SignalingProtein",
                           name="Kinase S100", measured=True, regulated=site_state,
                           log2fc=float(site_state))
            graph.add_edge("K1_S100", "P1", key="phosphorylation",
                           edge_type="phosphorylation", sign=sign)
        return graph

    def test_an_unsigned_site_reports_its_direction_but_no_effect(self):
        row = reaction_regulation_table(self._network(site_state=1, sign=0)).iloc[0]
        assert row["phospho_axis"] == 1
        assert row["phospho_axis_effect"] == 0
        assert row["phospho_axis_via"] == "K1_S100"

    def test_a_signed_activating_site_gets_an_effect(self):
        row = reaction_regulation_table(self._network(site_state=1, sign=1)).iloc[0]
        assert row["phospho_axis"] == 1 and row["phospho_axis_effect"] == 1

    def test_a_signed_inhibitory_site_reverses_the_effect(self):
        row = reaction_regulation_table(self._network(site_state=1, sign=-1)).iloc[0]
        assert row["phospho_axis"] == 1 and row["phospho_axis_effect"] == -1

    def test_phosphorylation_never_moves_the_gene_axis(self):
        """The whole point of separate columns: an unchanged enzyme stays 0."""
        row = reaction_regulation_table(self._network(site_state=-1, sign=-1)).iloc[0]
        assert row["gene_axis"] == 0
        assert row["phospho_axis"] == -1

    def test_enzyme_amount_steady_with_a_moved_site_is_visible(self):
        table = reaction_regulation_table(self._network(site_state=1, sign=0))
        row = table.iloc[0]
        assert row["gene_axis"] == 0 and row["phospho_axis"] != 0

    def test_a_network_without_a_signaling_layer_is_unchanged(self):
        row = reaction_regulation_table(self._network(site_state=None)).iloc[0]
        assert row["phospho_axis"] == 0
        assert row["phospho_axis_effect"] == 0
        assert row["phospho_axis_via"] == ""

    def test_an_unchanged_site_contributes_nothing(self):
        row = reaction_regulation_table(self._network(site_state=0)).iloc[0]
        assert row["phospho_axis"] == 0 and row["phospho_axis_via"] == ""

    def test_disagreeing_sites_give_no_direction_but_a_count(self):
        graph = self._network(site_state=1, sign=0)
        graph.add_node("K2_T50", layer="Signaling", node_type="SignalingProtein",
                       name="Kinase T50", measured=True, regulated=-1, log2fc=-1.0)
        graph.add_edge("K2_T50", "P1", key="phosphorylation",
                       edge_type="phosphorylation", sign=0)
        row = reaction_regulation_table(graph).iloc[0]
        assert row["phospho_axis"] == 0
        assert row["n_phosphosites_changed"] == 2

    def test_no_changed_site_counts_zero(self):
        row = reaction_regulation_table(self._network(site_state=0)).iloc[0]
        assert row["n_phosphosites_changed"] == 0
