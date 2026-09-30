"""The general graph analyses, on the typed network they now run against.

These were exported but never called or tested: `detect_communities` imported
an undeclared package and returned `{}` on any machine without it, and the
rest raised `NetworkXNotImplemented` on a directed multigraph -- which is the
only kind of graph this package builds.
"""

import networkx as nx
import pandas as pd
import pytest

from transnet.analysis.network_analysis import (
    compute_centrality_measures,
    compute_network_statistics,
    detect_communities,
    enrichment_analysis,
    find_active_modules,
    identify_hubs,
)


@pytest.fixture
def mapped(configured_network, omics_tables):
    from transnet import map_omics_to_network

    map_omics_to_network(configured_network, omics_tables, id_column="id",
                         log2fc_column="log2FC", qvalue_column="padj")
    return configured_network


class TestOnADirectedMultigraph:
    """Every one of these must work on the network the builder produces."""

    def test_statistics(self, mapped):
        stats = compute_network_statistics(mapped)
        assert stats["nodes"] == mapped.number_of_nodes()
        assert stats["connected_components"] >= 1
        assert 0 <= stats["largest_component_ratio"] <= 1

    def test_centrality(self, mapped):
        result = compute_centrality_measures(mapped, top_n=3)
        assert set(result) >= {"degree_centrality", "betweenness_centrality",
                               "closeness_centrality"}
        assert all(len(v) <= 3 for v in result.values())

    def test_hubs(self, mapped):
        hubs = identify_hubs(mapped, top_n=5)
        assert hubs and hubs == sorted(hubs, key=lambda pair: -pair[1])

    def test_communities(self, mapped):
        communities = detect_communities(mapped, method="louvain")
        assert communities, "louvain returned nothing -- the old import fell back to {}"
        assert set(communities) <= set(mapped.nodes())

    def test_every_community_method_agrees_on_the_node_set(self, mapped):
        for method in ("louvain", "label_propagation", "greedy_modularity"):
            communities = detect_communities(mapped, method=method)
            assert set(communities) == set(mapped.nodes()), method

    def test_active_modules(self, mapped):
        modules = find_active_modules(mapped, score_attr="log2fc", n_modules=2)
        assert modules
        for module in modules:
            assert set(module) <= set(mapped.nodes())

    def test_active_modules_are_connected_and_cross_layers(self):
        """Reactions carry no measurement, so growth has to step through them
        to join an enzyme to its metabolites."""
        from transnet import load_example_network, load_example_omics, map_omics_to_network

        graph = load_example_network()
        map_omics_to_network(graph, load_example_omics(), id_column="id",
                             log2fc_column="log2FC", qvalue_column="padj")
        largest = find_active_modules(graph, n_modules=1)[0]
        assert nx.is_connected(graph.subgraph(largest).to_undirected())
        assert {"Proteome", "Reactions", "Metabolome"} <= {
            graph.nodes[n]["layer"] for n in largest}

    def test_active_modules_leave_the_graph_alone(self, mapped):
        before = {n: dict(d) for n, d in mapped.nodes(data=True)}
        find_active_modules(mapped, score_attr="log2fc")
        assert {n: dict(d) for n, d in mapped.nodes(data=True)} == before


class TestPathLengthGuard:
    """Diameter is quadratic; an interactome would take hours."""

    def test_skipped_above_the_threshold(self):
        graph = nx.path_graph(40)
        stats = compute_network_statistics(graph, exact_paths_below=10)
        assert stats["diameter"] is None and stats["avg_path_length"] is None

    def test_computed_below_it(self):
        graph = nx.path_graph(6)
        stats = compute_network_statistics(graph, exact_paths_below=100)
        assert stats["diameter"] == 5


class TestEnrichment:

    def test_benjamini_hochberg_uses_the_rank_after_sorting(self):
        """The denominator once came from the pre-sort insertion order, so
        every q-value this function returned was wrong."""
        background = [f"g{i}" for i in range(100)]
        annotations = {"hit": background[:10], "miss": background[50:60]}
        result = enrichment_analysis(background[:10], background, annotations)

        ranked = result.sort_values("P_value").reset_index(drop=True)
        assert (ranked["FDR"].diff().dropna() >= -1e-12).all(), "q-values not monotone"
        assert (ranked["FDR"] >= ranked["P_value"] - 1e-12).all()
