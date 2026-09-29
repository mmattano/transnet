"""Mapping, subnetworks, signed paths, topology, temporal structure, propagation."""

import numpy as np
import pandas as pd
import pytest

from transnet.analysis.transomics import (
    assign_temporal_parameters,
    compare_transomic_networks,
    cross_layer_connectivity,
    downstream_influence,
    hierarchical_propagation,
    layer_coverage,
    map_omics_to_network,
    path_consistency_summary,
    regulated_nodes,
    responsive_subnetwork,
    split_by_response_class,
    temporal_network_structure,
    trace_regulatory_paths,
    transomic_hubs,
)


@pytest.fixture
def mapped_network(standard_network, omics_tables):
    map_omics_to_network(
        standard_network,
        {k: v for k, v in omics_tables.items() if k != "Signaling"},
        id_column="id", log2fc_column="log2FC", qvalue_column="padj",
    )
    return standard_network


class TestMapping:
    def test_values_and_direction_land_on_nodes(self, mapped_network):
        node = mapped_network.nodes["C00668"]
        assert node["log2fc"] == 2.0
        assert node["qvalue"] == 0.001
        assert node["regulated"] == 1

    def test_non_significant_features_are_not_called_regulated(self, mapped_network):
        # ATP has padj 0.6
        assert mapped_network.nodes["C00002"]["regulated"] == 0
        assert mapped_network.nodes["C00002"]["log2fc"] == 0.2

    def test_report_counts_matches_per_layer(self, standard_network, omics_tables):
        report = map_omics_to_network(
            standard_network,
            {k: v for k, v in omics_tables.items() if k != "Signaling"},
            id_column="id", log2fc_column="log2FC", qvalue_column="padj",
        )
        row = report.per_layer.set_index("layer").loc["Metabolome"]
        assert row["n_supplied"] == 4
        assert row["n_matched"] == 4
        assert row["n_up"] == 2 and row["n_down"] == 1
        assert report.match_rate == pytest.approx(1.0)

    def test_a_table_for_a_missing_layer_is_reported_not_raised(
        self, standard_network, omics_tables
    ):
        report = map_omics_to_network(
            standard_network, {"Signaling": omics_tables["Signaling"]},
            id_column="id", log2fc_column="log2FC", qvalue_column="padj",
        )
        row = report.per_layer.set_index("layer").loc["Signaling"]
        assert row["n_matched"] == 0
        assert "layer absent" in row["note"]

    def test_unmatched_identifiers_are_surfaced(self, standard_network):
        report = map_omics_to_network(
            standard_network,
            {"Metabolome": pd.DataFrame({
                "id": ["glucose", "G6P"], "log2FC": [1.0, 1.0], "padj": [0.01, 0.01]})},
            id_column="id", log2fc_column="log2FC", qvalue_column="padj",
        )
        assert report.per_layer.iloc[0]["n_matched"] == 0
        assert "glucose" in report.unmatched["Metabolome"]
        assert "identifier-type mismatch" in str(report)
        assert "build_name_id_map" in str(report)

    def test_a_smaller_network_is_not_reported_as_an_id_mismatch(
        self, standard_network
    ):
        """A low match rate has two causes; the report must tell them apart."""
        metabolites = [
            n for n, d in standard_network.nodes(data=True)
            if d.get("layer") == "Metabolome"
        ]
        # Every network node is covered, plus features the network lacks.
        table = pd.DataFrame({
            "id": metabolites + [f"absent_{i}" for i in range(len(metabolites) * 3)],
        })
        table["log2FC"] = 1.0
        table["padj"] = 0.01

        report = map_omics_to_network(
            standard_network, {"Metabolome": table},
            id_column="id", log2fc_column="log2FC", qvalue_column="padj",
        )
        row = report.per_layer.iloc[0]
        assert row["match_fraction"] < 0.5
        assert row["layer_covered"] == pytest.approx(1.0)

        text = str(report)
        assert "the network is simply smaller than the assay" in text
        assert "identifier-type mismatch" not in text

    def test_layer_coverage_counts_distinct_nodes(self, standard_network):
        """Two aliases resolving to one node must not push coverage past 100%."""
        report = map_omics_to_network(
            standard_network,
            {"Metabolome": pd.DataFrame({
                "id": ["a", "b"], "log2FC": [1.0, 1.0], "padj": [0.01, 0.01]})},
            id_column="id", log2fc_column="log2FC", qvalue_column="padj",
            id_map={"Metabolome": {"a": "C00668", "b": "C00668"}},
        )
        row = report.per_layer.iloc[0]
        assert row["n_matched"] == 2
        assert row["layer_covered"] <= 1.0

    def test_id_map_translates_foreign_identifiers(self, standard_network):
        report = map_omics_to_network(
            standard_network,
            {"Metabolome": pd.DataFrame({
                "id": ["glucose"], "log2FC": [1.0], "padj": [0.01]})},
            id_column="id", log2fc_column="log2FC", qvalue_column="padj",
            id_map={"Metabolome": {"glucose": "C00031"}},
        )
        assert report.per_layer.iloc[0]["n_matched"] == 1
        assert standard_network.nodes["C00031"]["regulated"] == 1

    def test_regulated_nodes_filters_by_layer_and_direction(self, mapped_network):
        up = regulated_nodes(mapped_network, layers=["Metabolome"], direction=1)
        assert set(up) == {"C00668", "C00354"}
        assert "C00031" in regulated_nodes(mapped_network, direction=-1)


class TestResponsiveSubnetwork:
    def test_only_regulated_molecules_and_their_connectors_survive(self, mapped_network):
        sub = responsive_subnetwork(mapped_network)
        for node in sub.nodes:
            layer = sub.nodes[node].get("layer")
            regulated = sub.nodes[node].get("regulated", 0)
            assert regulated != 0 or layer == "Reactions"

    def test_reactions_are_kept_as_connectors(self, mapped_network):
        sub = responsive_subnetwork(mapped_network)
        assert "R00299" in sub.nodes

    def test_connectors_can_be_switched_off(self, mapped_network):
        sub = responsive_subnetwork(mapped_network, include_connectors=False)
        assert all(sub.nodes[n].get("regulated", 0) != 0 for n in sub.nodes)

    def test_no_regulated_nodes_returns_an_empty_graph(self, standard_network):
        sub = responsive_subnetwork(standard_network)
        assert sub.number_of_nodes() == 0

    def test_edge_attributes_are_preserved(self, mapped_network):
        sub = responsive_subnetwork(mapped_network)
        for _, _, data in sub.edges(data=True):
            assert "edge_type" in data and "sign" in data


class TestComparison:
    def test_identical_networks_have_jaccard_one(self, mapped_network):
        result = compare_transomic_networks(mapped_network, mapped_network)
        assert result["summary"]["edge_jaccard"] == pytest.approx(1.0)
        assert result["summary"]["n_regulation_shifts"] == 0

    def test_regulation_shifts_are_detected(self, network_builder, omics_tables):
        first = network_builder("no_signaling").generate_graph()
        second = network_builder("no_signaling").generate_graph()
        map_omics_to_network(
            first, {"Metabolome": pd.DataFrame({
                "id": ["C00668"], "log2FC": [2.0], "padj": [0.01]})},
            id_column="id", log2fc_column="log2FC", qvalue_column="padj")
        map_omics_to_network(
            second, {"Metabolome": pd.DataFrame({
                "id": ["C00668"], "log2FC": [-2.0], "padj": [0.01]})},
            id_column="id", log2fc_column="log2FC", qvalue_column="padj")

        result = compare_transomic_networks(first, second, "wt", "ko")
        shifts = result["regulation_shifts"]
        assert len(shifts) == 1
        assert shifts.iloc[0]["node"] == "C00668"
        assert shifts.iloc[0]["wt"] == 1 and shifts.iloc[0]["ko"] == -1

    def test_edge_types_are_broken_out(self, mapped_network):
        result = compare_transomic_networks(mapped_network, mapped_network)
        assert "edge_type" in result["edges_by_type"].columns
        assert set(result["edges_by_type"]["edge_type"]) >= {"catalysis", "product"}


class TestPaths:
    def test_sign_is_the_product_of_edge_signs(self, mapped_network):
        paths = trace_regulatory_paths(
            mapped_network, source_layer="Proteome", target_layer="Metabolome"
        )
        assert not paths.empty
        # HK1 -> R00299 -> G6P: catalysis (+1) then product (+1)
        row = paths[paths["path"] == "P19367 -> R00299 -> C00668"]
        assert len(row) == 1
        assert row.iloc[0]["sign"] == 1
        assert row.iloc[0]["unsigned_steps"] == 0

    def test_inhibitory_step_flips_the_predicted_sign(self, mapped_network):
        paths = trace_regulatory_paths(
            mapped_network, sources=["C00668"], targets=["C00008"], max_length=3
        )
        # G6P -| R00299 -> ADP passes through an inhibition
        via_inhibition = paths[paths["edge_types"].str.contains("allosteric_inhibition")]
        assert len(via_inhibition) >= 1
        assert (via_inhibition["sign"] == -1).all()

    def test_consistency_compares_prediction_with_measurement(self, mapped_network):
        paths = trace_regulatory_paths(
            mapped_network, source_layer="Proteome", target_layer="Metabolome"
        )
        for _, row in paths.iterrows():
            assert row["predicted"] == row["sign"] * row["source_regulated"]
            assert row["consistent"] == bool(row["observed"] and row["predicted"] == row["observed"])

    def test_a_decreased_source_predicts_a_decrease_along_a_positive_path(self):
        """Fbp1 down -> FBPase -> fructose-6-phosphate was reported as
        predicting an increase, because the source's direction was ignored."""
        import networkx as nx
        graph = nx.MultiDiGraph()
        graph.add_node("E", layer="Proteome", regulated=-1, responsive=True)
        graph.add_node("R", layer="Reactions")
        graph.add_node("M_down", layer="Metabolome", regulated=-1, responsive=True)
        graph.add_node("M_up", layer="Metabolome", regulated=1, responsive=True)
        graph.add_edge("E", "R", edge_type="catalysis", sign=1)
        graph.add_edge("R", "M_down", edge_type="product", sign=1)
        graph.add_edge("R", "M_up", edge_type="product", sign=1)

        paths = trace_regulatory_paths(graph, source_layer="Proteome",
                                       target_layer="Metabolome").set_index("target")
        assert paths.loc["M_down", "predicted"] == -1
        assert bool(paths.loc["M_down", "consistent"]) is True
        assert bool(paths.loc["M_up", "consistent"]) is False

    def test_a_measured_unchanged_intermediate_blocks_the_path(self):
        """In MoTrPAC muscle, 247 of the 359 "consistent" paths ran through
        glutathione, which was measured and had not changed. A molecule that
        did not move passed nothing on, so those paths were not evidence."""
        import networkx as nx
        graph = nx.MultiDiGraph()
        graph.add_node("E", layer="Proteome", regulated=1, responsive=True, measured=True)
        graph.add_node("R1", layer="Reactions")
        graph.add_node("R2", layer="Reactions")
        graph.add_node("Mid", layer="Metabolome", regulated=0, responsive=False, measured=True)
        graph.add_node("Unmeasured", layer="Metabolome")
        graph.add_node("End", layer="Metabolome", regulated=1, responsive=True, measured=True)
        graph.add_edge("E", "R1", edge_type="catalysis", sign=1)
        graph.add_edge("R1", "Mid", edge_type="product", sign=1)
        graph.add_edge("R1", "Unmeasured", edge_type="product", sign=1)
        graph.add_edge("Mid", "R2", edge_type="allosteric_activation", sign=1)
        graph.add_edge("Unmeasured", "R2", edge_type="allosteric_activation", sign=1)
        graph.add_edge("R2", "End", edge_type="product", sign=1)

        paths = trace_regulatory_paths(graph, source_layer="Proteome",
                                       target_layer="Metabolome", max_length=4)
        through = set(paths["path"])
        assert "E -> R1 -> Unmeasured -> R2 -> End" in through   # untested, so kept
        assert "E -> R1 -> Mid -> R2 -> End" not in through      # contradicted
        assert (paths["unchanged_intermediates"] == 0).all()

        kept = trace_regulatory_paths(graph, source_layer="Proteome",
                                      target_layer="Metabolome", max_length=4,
                                      allow_unchanged_intermediates=True)
        blocked = kept[kept["path"] == "E -> R1 -> Mid -> R2 -> End"]
        assert len(blocked) == 1 and blocked.iloc[0]["unchanged_intermediates"] == 1

    def test_the_endpoints_own_state_does_not_block_its_paths(self):
        """Only intermediates are checked: the target's measurement is the
        thing being predicted, and the source's direction drives it."""
        import networkx as nx
        graph = nx.MultiDiGraph()
        graph.add_node("E", layer="Proteome", regulated=1, responsive=True, measured=True)
        graph.add_node("R", layer="Reactions")
        graph.add_node("M", layer="Metabolome", regulated=1, responsive=True, measured=True)
        graph.add_edge("E", "R", edge_type="catalysis", sign=1)
        graph.add_edge("R", "M", edge_type="product", sign=1)
        assert len(trace_regulatory_paths(graph, source_layer="Proteome",
                                          target_layer="Metabolome")) == 1

    def test_a_target_is_scored_once_however_many_paths_reach_it(self):
        """Paths sharing a hub are not independent observations."""
        import networkx as nx
        graph = nx.MultiDiGraph()
        graph.add_node("E", layer="Proteome", regulated=1, responsive=True, measured=True)
        graph.add_node("M", layer="Metabolome", regulated=-1, responsive=True, measured=True)
        for i in range(5):
            graph.add_node(f"R{i}", layer="Reactions")
            graph.add_edge("E", f"R{i}", edge_type="catalysis", sign=1)
            graph.add_edge(f"R{i}", "M", edge_type="product", sign=1 if i else -1)
        paths = trace_regulatory_paths(graph, source_layer="Proteome",
                                       target_layer="Metabolome")
        assert len(paths) == 5 and int(paths["consistent"].sum()) == 1

        summary = path_consistency_summary(paths)
        assert len(summary) == 1                    # one molecule, one verdict
        assert summary.iloc[0]["predicted"] == 1    # 4 paths against 1
        assert bool(summary.iloc[0]["agrees"]) is False

    def test_source_layer_is_inferred_when_not_given(self, mapped_network):
        paths = trace_regulatory_paths(mapped_network)
        assert not paths.empty
        assert paths["layers"].str.startswith("Proteome").all()

    def test_absent_target_layer_returns_an_empty_typed_frame(self, mapped_network):
        paths = trace_regulatory_paths(mapped_network, target_layer="Signaling")
        assert paths.empty
        assert "consistent" in paths.columns

    def test_summary_reports_the_shortest_consistent_path(self, mapped_network):
        paths = trace_regulatory_paths(
            mapped_network, source_layer="Proteome", target_layer="Metabolome"
        )
        summary = path_consistency_summary(paths)
        assert "shortest_consistent_path" in summary.columns
        assert (summary["n_consistent"] <= summary["n_paths"]).all()


class TestTopology:
    def test_cross_layer_fraction_is_reported(self, standard_network):
        result = cross_layer_connectivity(standard_network)
        assert 0.0 <= result["cross_layer_fraction"] <= 1.0
        assert result["matrix"].loc["Metabolome", "Reactions"] > 0

    def test_edge_type_table_names_the_layer_pairs(self, standard_network):
        result = cross_layer_connectivity(standard_network)
        types = result["edge_types"].set_index("edge_type")
        assert "Metabolome->Reactions" in types.loc["substrate", "layer_pairs"]

    def test_coverage_counts_measured_and_regulated(self, mapped_network):
        coverage = layer_coverage(mapped_network).set_index("layer")
        assert coverage.loc["Metabolome", "n_measured"] == 4
        assert coverage.loc["Metabolome", "n_regulated"] == 3

    def test_hubs_report_cross_layer_reach(self, mapped_network):
        hubs = transomic_hubs(mapped_network, top_percent=50)
        assert {"degree", "cross_layer_degree", "n_layers_touched",
                "versatility", "is_hub"} <= set(hubs.columns)
        assert hubs["is_hub"].any()

    def test_versatility_rewards_spanning_more_layers(self, standard_network):
        hubs = transomic_hubs(standard_network, responsive_only=False)
        multi = hubs[hubs["n_layers_touched"] > 1]
        single = hubs[hubs["n_layers_touched"] == 1]
        if len(multi) and len(single):
            comparable = single[single["degree"] >= multi["degree"].min()]
            if len(comparable):
                assert multi["versatility"].max() > comparable["versatility"].min()

    def test_hubs_fall_back_to_all_nodes_without_data(self, standard_network):
        hubs = transomic_hubs(standard_network)
        assert not hubs.empty


class TestTemporal:
    @pytest.fixture
    def timecourse(self):
        # G6P rises fast, glucose falls slowly.
        return {"Metabolome": pd.DataFrame({
            "id": ["C00668", "C00031"],
            "0": [0.0, 0.0],
            "10": [0.8, -0.1],
            "30": [1.0, -0.5],
            "60": [1.0, -1.0],
        })}

    def test_t_half_lands_on_nodes(self, standard_network, timecourse):
        table = assign_temporal_parameters(
            standard_network, timecourse,
            time_columns=[0, 10, 30, 60], id_column="id",
        )
        assert len(table) == 2
        assert standard_network.nodes["C00668"]["t_half"] is not None
        # fast responder reaches half-max before the slow one
        fast = standard_network.nodes["C00668"]["t_half"]
        slow = standard_network.nodes["C00031"]["t_half"]
        assert fast < slow

    def test_structure_reports_degree_vs_thalf(self, standard_network, timecourse):
        assign_temporal_parameters(
            standard_network, timecourse,
            time_columns=[0, 10, 30, 60], id_column="id")
        result = temporal_network_structure(standard_network)
        assert "degree_vs_thalf" in result
        assert "per_layer_thalf" in result

    def test_neighbour_correlation_when_profiles_are_supplied(
        self, standard_network, timecourse
    ):
        assign_temporal_parameters(
            standard_network, timecourse,
            time_columns=[0, 10, 30, 60], id_column="id")
        profiles = pd.DataFrame(
            np.random.RandomState(0).normal(size=(len(standard_network), 4)),
            index=list(standard_network.nodes),
        )
        result = temporal_network_structure(standard_network, profiles)
        assert result["neighbour_correlation"]["n_pairs"] > 0
        assert 0.0 <= result["neighbour_correlation"]["fraction_strong"] <= 1.0

    def test_response_classes_split_the_network(self, standard_network, timecourse):
        assign_temporal_parameters(
            standard_network, timecourse,
            time_columns=[0, 10, 30, 60], id_column="id")
        result = split_by_response_class(standard_network)
        assert result["thresholds"]["t_half"] is not None
        assert set(result["classes"]["response_class"]) <= {"fast", "slow"}

    def test_no_temporal_data_is_handled(self, standard_network):
        result = split_by_response_class(standard_network)
        assert result["subnetworks"] == {}


class TestPropagation:
    def test_scores_flow_downstream_from_a_seed(self, mapped_network):
        scores = hierarchical_propagation(mapped_network, {"P19367": 1.0})
        assert not scores.empty
        reached = set(scores["node"])
        assert "R00299" in reached and "C00668" in reached

    def test_inhibitory_edges_flip_the_predicted_direction(self, mapped_network):
        scores = hierarchical_propagation(
            mapped_network, {"C00668": 1.0}
        ).set_index("node")
        # G6P inhibits R00299, so a rise in G6P predicts the reaction going down
        assert scores.loc["R00299", "predicted_direction"] == -1

    def test_sign_can_be_ignored(self, mapped_network):
        signed = hierarchical_propagation(
            mapped_network, {"C00668": 1.0}).set_index("node")
        unsigned = hierarchical_propagation(
            mapped_network, {"C00668": 1.0}, respect_sign=False).set_index("node")
        assert signed.loc["R00299", "score"] < 0
        assert unsigned.loc["R00299", "score"] > 0

    def test_missing_seeds_are_dropped_not_fatal(self, mapped_network):
        scores = hierarchical_propagation(
            mapped_network, {"P19367": 1.0, "NOT_A_NODE": 1.0})
        assert not scores.empty

    def test_no_valid_seed_returns_an_empty_typed_frame(self, mapped_network):
        scores = hierarchical_propagation(mapped_network, {"NOT_A_NODE": 1.0})
        assert scores.empty
        assert "predicted_direction" in scores.columns

    def test_downstream_influence_compares_with_measurement(self, mapped_network):
        result = downstream_influence(mapped_network, {"P19367": 1.0})
        assert (result["layer"] == "Metabolome").all()
        assert "agrees" in result.columns
        measured = result[result["observed"] != 0]
        assert measured["agrees"].notna().all()

    def test_unmeasured_targets_are_na_not_false(self, mapped_network):
        result = downstream_influence(mapped_network, {"P19367": 1.0})
        unmeasured = result[result["observed"] == 0]
        assert unmeasured["agrees"].isna().all()


class TestPathEdgeFiltering:
    """Protein associations are not regulatory steps."""

    def test_ppi_is_excluded_by_default(self, mapped_network):
        paths = trace_regulatory_paths(
            mapped_network, source_layer="Proteome", target_layer="Metabolome"
        )
        assert not paths.empty
        assert not paths["edge_types"].str.contains("protein_interaction").any()

    def test_ppi_can_be_kept(self, mapped_network):
        with_ppi = trace_regulatory_paths(
            mapped_network, source_layer="Proteome", target_layer="Metabolome",
            exclude_edge_types=(),
        )
        without = trace_regulatory_paths(
            mapped_network, source_layer="Proteome", target_layer="Metabolome",
        )
        assert len(with_ppi) >= len(without)

    def test_an_arbitrary_type_can_be_excluded(self, mapped_network):
        paths = trace_regulatory_paths(
            mapped_network, source_layer="Proteome", target_layer="Metabolome",
            exclude_edge_types=("product",),
        )
        assert not paths["edge_types"].str.contains("product").any()

    def test_excluding_everything_yields_no_paths(self, mapped_network):
        paths = trace_regulatory_paths(
            mapped_network, source_layer="Proteome", target_layer="Metabolome",
            exclude_edge_types=(
                "catalysis", "product", "substrate", "translation",
                "allosteric_activation", "allosteric_inhibition",
                "transcriptional_regulation", "protein_interaction",
            ),
        )
        assert paths.empty
        assert "consistent" in paths.columns


class TestPathNodeFiltering:
    """A path through water joins reactions that have nothing to do with each other."""

    def test_currency_metabolites_are_routed_around_by_default(self, mapped_network):
        paths = trace_regulatory_paths(
            mapped_network, source_layer="Proteome", target_layer="Metabolome"
        )
        # ADP and ATP are currency; no path may pass through them.
        for currency in ("C00008", "C00002"):
            assert not paths["path"].str.contains(currency).any()

    def test_keeping_them_finds_more_paths(self):
        """On a network with shared cofactors, hopping creates spurious paths."""
        from transnet import load_example_network, load_example_omics

        graph = load_example_network()
        map_omics_to_network(
            graph, load_example_omics(),
            id_column="id", log2fc_column="log2FC", qvalue_column="padj",
        )
        without = trace_regulatory_paths(
            graph, source_layer="Proteome", target_layer="Metabolome")
        with_currency = trace_regulatory_paths(
            graph, source_layer="Proteome", target_layer="Metabolome",
            exclude_nodes=())
        assert len(with_currency) > len(without), (
            "cofactor hopping should create extra, spurious paths"
        )

    def test_an_explicit_endpoint_overrides_the_exclusion(self, mapped_network):
        paths = trace_regulatory_paths(
            mapped_network, sources=["C00002"], target_layer="Metabolome",
            max_length=3, regulated_only=False,
        )
        assert not paths.empty, (
            "naming a currency metabolite as a source must keep it"
        )
        assert (paths["source"] == "C00002").all()

    def test_the_exclusion_set_can_be_overridden(self, mapped_network):
        paths = trace_regulatory_paths(
            mapped_network, source_layer="Proteome", target_layer="Metabolome",
            exclude_nodes=("C00668",),
        )
        assert not paths["path"].str.contains("C00668").any()


class TestPathVerdictsInput:
    """``path_verdicts`` takes traced paths, and must not reject them.

    A guard that keyed on ``target`` looked right and was wrong: the traced paths
    and the per-target summary both have that column, so the guard rejected every
    legitimate call. ``n_paths`` exists only in the summary.
    """

    def test_traced_paths_are_accepted(self, mapped_network):
        from transnet import trace_regulatory_paths
        from transnet.analysis import path_verdicts

        paths = trace_regulatory_paths(mapped_network, target_layer="Metabolome",
                                       max_length=4)
        summary, agree, tested = path_verdicts(paths)
        assert tested >= 0 and agree <= tested

    def test_passing_the_summary_instead_says_so(self, mapped_network):
        import pytest

        from transnet import trace_regulatory_paths
        from transnet.analysis import path_verdicts

        paths = trace_regulatory_paths(mapped_network, target_layer="Metabolome",
                                       max_length=4)
        summary, _, _ = path_verdicts(paths)
        with pytest.raises(ValueError, match="not the per-target summary"):
            path_verdicts(summary)
