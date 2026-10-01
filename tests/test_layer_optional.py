"""No analysis may require a particular layer to be present.

A study with only transcriptomics and metabolomics is a normal study, not a
degraded input.  Every function in the trans-omics API must run on every layer
configuration, return a correctly typed result, and say in its output which
evidence it actually had.

The ``configured_network`` fixture is parametrised over:

``full``
    Every layer, Signaling included.
``no_signaling``
    The common case: transcriptome, proteome, metabolome, reactions.
``minimal``
    Transcriptome + metabolome + reactions, no Proteome at all -- so no
    catalysis, no PPI, no allosteric edges exist.
"""

import pandas as pd
import pytest

from transnet.analysis.transomics import (
    assign_temporal_parameters,
    compare_transomic_networks,
    cross_layer_connectivity,
    downstream_influence,
    expression_concordance,
    hierarchical_propagation,
    layer_coverage,
    map_omics_to_network,
    path_consistency_summary,
    reaction_regulation_table,
    regulated_nodes,
    regulation_axis_summary,
    responsive_subnetwork,
    split_by_response_class,
    temporal_network_structure,
    trace_regulatory_paths,
    transcription_factor_activity,
    transomic_hubs,
)


@pytest.fixture
def mapped(configured_network, omics_tables):
    """Each configuration, with every omics table offered to it.

    Tables for layers the network lacks are deliberately included: supplying
    phosphoproteomics to a network with no Signaling layer must be reported,
    not raised.
    """
    map_omics_to_network(
        configured_network, omics_tables,
        id_column="id", log2fc_column="log2FC", qvalue_column="padj",
    )
    return configured_network


class TestEveryAnalysisRunsOnEveryConfiguration:
    def test_mapping(self, configured_network, omics_tables):
        report = map_omics_to_network(
            configured_network, omics_tables,
            id_column="id", log2fc_column="log2FC", qvalue_column="padj",
        )
        assert not report.per_layer.empty
        assert report.per_layer["n_matched"].sum() > 0

    def test_cross_layer_connectivity(self, mapped):
        result = cross_layer_connectivity(mapped)
        assert 0.0 <= result["cross_layer_fraction"] <= 1.0
        assert not result["matrix"].empty

    def test_layer_coverage(self, mapped):
        coverage = layer_coverage(mapped)
        assert not coverage.empty
        assert coverage["n_regulated"].sum() > 0

    def test_responsive_subnetwork(self, mapped):
        sub = responsive_subnetwork(mapped)
        assert sub.number_of_nodes() > 0

    def test_reaction_regulation_table(self, mapped):
        table = reaction_regulation_table(mapped)
        assert "controversial" in table.columns
        assert len(table) == 2   # both fixture reactions, in every configuration

    def test_regulation_axis_summary(self, mapped):
        summary = regulation_axis_summary(reaction_regulation_table(mapped))
        assert "fraction_controversial" in summary.columns

    def test_trace_regulatory_paths(self, mapped):
        paths = trace_regulatory_paths(mapped)
        assert "consistent" in paths.columns
        assert not paths.empty
        assert path_consistency_summary(paths) is not None

    def test_transomic_hubs(self, mapped):
        hubs = transomic_hubs(mapped, top_percent=50)
        assert not hubs.empty
        assert "versatility" in hubs.columns

    def test_temporal(self, mapped):
        table = assign_temporal_parameters(
            mapped,
            {"Metabolome": pd.DataFrame({
                "id": ["C00668", "C00031"],
                "0": [0.0, 0.0], "10": [0.8, -0.1],
                "30": [1.0, -0.5], "60": [1.0, -1.0]})},
            time_columns=[0, 10, 30, 60], id_column="id",
        )
        assert len(table) == 2
        assert "degree_vs_thalf" in temporal_network_structure(mapped)
        assert "thresholds" in split_by_response_class(mapped)

    def test_propagation(self, mapped):
        seed = next(
            n for n, d in mapped.nodes(data=True) if d.get("layer") == "Metabolome"
        )
        scores = hierarchical_propagation(mapped, {seed: 1.0})
        assert "predicted_direction" in scores.columns
        assert "agrees" in downstream_influence(mapped, {seed: 1.0}).columns

    def test_comparison(self, mapped):
        result = compare_transomic_networks(mapped, mapped)
        assert result["summary"]["edge_jaccard"] == pytest.approx(1.0)

    def test_regulated_nodes(self, mapped):
        assert len(regulated_nodes(mapped)) > 0

    def test_transcription_factor_activity(self, mapped):
        # typed and non-raising on every configuration, with or without TFs
        table = transcription_factor_activity(mapped, min_targets=1)
        assert "inferred_activity" in table.columns

    def test_expression_concordance(self, mapped):
        result = expression_concordance(mapped)
        assert set(result["counts"]) == {
            "concordant", "protein_only", "transcript_only", "discordant", "unchanged",
            "tested_difference", "protein_beyond_transcript",
        }


class TestEvidenceIsReportedHonestly:
    """A thinner network must produce visibly weaker claims, not silent ones."""

    def test_gene_axis_evidence_matches_the_layers_available(self, mapped):
        configuration = mapped.graph["configuration"]
        row = reaction_regulation_table(mapped).set_index("reaction").loc["R00299"]

        if configuration == "minimal":
            # No Proteome: the transcript is the only enzyme-level evidence.
            assert row["gene_axis_evidence"] == "gene"
        else:
            assert row["gene_axis_evidence"] == "protein"

    def test_metabolite_axis_needs_only_metabolome_and_reactions(self, mapped):
        table = reaction_regulation_table(mapped)
        # G6P is measured in every configuration and is a product of R00299,
        # so the metabolite axis is informative regardless of the upper layers.
        assert (table["metabolite_axis"] != 0).any()

    def test_allosteric_regulators_appear_only_where_brenda_edges_exist(self, mapped):
        configuration = mapped.graph["configuration"]
        table = reaction_regulation_table(mapped)
        has_allostery = table["allosteric_regulators"].str.len().gt(0).any()
        assert has_allostery == (configuration != "minimal")

    def test_path_tracing_reports_the_hierarchy_it_used(self, mapped):
        configuration = mapped.graph["configuration"]
        paths = trace_regulatory_paths(mapped)
        first_layer = paths.iloc[0]["layers"].split(" -> ")[0]

        expected = {
            "full": "Signaling",
            "no_signaling": "Proteome",
            "minimal": "Transcriptome",
        }[configuration]
        assert first_layer == expected


class TestAbsentLayersAreWarningsNotErrors:
    def test_data_for_an_absent_layer_is_reported(self, configured_network, omics_tables):
        report = map_omics_to_network(
            configured_network, omics_tables,
            id_column="id", log2fc_column="log2FC", qvalue_column="padj",
        )
        configuration = configured_network.graph["configuration"]
        by_layer = report.per_layer.set_index("layer")

        if configuration != "full":
            assert by_layer.loc["Signaling", "n_matched"] == 0
            assert "layer absent" in by_layer.loc["Signaling", "note"]
        if configuration == "minimal":
            assert by_layer.loc["Proteome", "n_matched"] == 0

    def test_targeting_an_absent_layer_returns_an_empty_typed_frame(self, mapped):
        if mapped.graph["configuration"] == "full":
            pytest.skip("this configuration has a Signaling layer")
        paths = trace_regulatory_paths(mapped, target_layer="Signaling")
        assert paths.empty
        assert list(paths.columns)

    def test_propagating_to_an_absent_layer_is_empty_not_fatal(self, mapped):
        if mapped.graph["configuration"] == "full":
            pytest.skip("this configuration has a Signaling layer")
        result = downstream_influence(
            mapped, {"C00668": 1.0}, target_layer="Signaling"
        )
        assert result.empty
        assert "agrees" in result.columns

    def test_network_with_no_reactions_layer_degrades_cleanly(self, network_builder):
        graph = network_builder("no_reactions").generate_graph()
        map_omics_to_network(
            graph,
            {"Proteome": pd.DataFrame({
                "id": ["P19367"], "log2FC": [1.0], "padj": [0.01]})},
            id_column="id", log2fc_column="log2FC", qvalue_column="padj",
        )
        table = reaction_regulation_table(graph)
        assert table.empty and "controversial" in table.columns
        assert not layer_coverage(graph).empty
        assert cross_layer_connectivity(graph)["cross_layer_fraction"] >= 0.0
