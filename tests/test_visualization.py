"""The trans-omics figures build on every layer configuration."""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt              # noqa: E402
import pandas as pd                          # noqa: E402
import pytest                                # noqa: E402

from transnet.analysis.transomics import (   # noqa: E402
    cross_layer_connectivity,
    map_omics_to_network,
    reaction_regulation_table,
    regulation_axis_summary,
)
from transnet.visualization import (         # noqa: E402
    EDGE_STYLES,
    plot_layer_connectivity,
    plot_regulation_axes,
    plot_transomic_network,
)
from transnet.biology.schema import EDGE_TYPES   # noqa: E402


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


@pytest.fixture
def mapped(configured_network, omics_tables):
    map_omics_to_network(
        configured_network, omics_tables,
        id_column="id", log2fc_column="log2FC", qvalue_column="padj",
    )
    return configured_network


class TestTransomicNetworkFigure:
    def test_builds_for_every_layer_configuration(self, mapped):
        figure = plot_transomic_network(mapped)
        assert figure.axes

    def test_layer_order_follows_the_hierarchy(self, full_network):
        figure = plot_transomic_network(full_network)
        labels = [t.get_text() for t in figure.axes[0].texts]
        assert labels.index("Signaling") < labels.index("Metabolome")

    def test_absent_layers_are_simply_not_drawn(self, standard_network):
        figure = plot_transomic_network(standard_network)
        labels = [t.get_text() for t in figure.axes[0].texts]
        assert "Signaling" not in labels

    def test_colour_by_layer_is_available(self, mapped):
        assert plot_transomic_network(mapped, node_color_by="layer").axes

    def test_graph_without_layers_raises_a_clear_error(self):
        import networkx as nx
        graph = nx.MultiDiGraph()
        graph.add_edge("a", "b")
        with pytest.raises(ValueError, match="layer"):
            plot_transomic_network(graph)


class TestBackboneSelection:
    """BRENDA lists every effector an enzyme was tested with. Laboratory
    reagents such as p-chloromercuribenzoate were deciding which reactions the
    MoTrPAC muscle figure showed."""

    @staticmethod
    def _graph():
        import networkx as nx
        graph = nx.MultiDiGraph()
        graph.add_node("enzyme", layer="Proteome", regulated=1, responsive=True, measured=True)
        for name in ("endogenous", "reagent"):
            graph.add_node(name, layer="Metabolome", regulated=1, responsive=True, measured=True)
        graph.add_node("real", layer="Reactions")
        graph.add_node("tested_only", layer="Reactions")
        graph.add_node("product", layer="Metabolome", regulated=1, responsive=True, measured=True)
        for reaction in ("real", "tested_only"):
            graph.add_edge("enzyme", reaction, edge_type="catalysis", sign=1)
            graph.add_edge(reaction, "product", edge_type="product", sign=1)
        # only "endogenous" is made or consumed by a reaction here
        graph.add_edge("endogenous", "real", edge_type="substrate", sign=1)
        graph.add_edge("endogenous", "real", edge_type="allosteric_inhibition", sign=-1)
        graph.add_edge("reagent", "tested_only", edge_type="allosteric_inhibition", sign=-1)
        return graph

    def test_effectors_outside_the_organisms_metabolism_are_not_drawn(self):
        from transnet.visualization import transomic_backbone
        backbone = transomic_backbone(self._graph())
        assert "endogenous" in backbone
        assert "reagent" not in backbone

    def test_they_can_be_drawn_on_request(self):
        from transnet.visualization import transomic_backbone
        backbone = transomic_backbone(self._graph(), metabolic_effectors_only=False)
        assert "reagent" in backbone

    def test_substrates_and_products_are_never_filtered(self):
        """The pool is about allosteric annotation, not about being unusual."""
        from transnet.visualization import transomic_backbone
        graph = self._graph()
        assert "product" in transomic_backbone(graph)


class TestPromiscuousEnzymes:
    """Fifteen glutathione-transferase reactions, one per xenobiotic, filled
    the MoTrPAC muscle figure and hid every other regulated reaction."""

    @staticmethod
    def _graph(n_reactions=10):
        import networkx as nx
        graph = nx.MultiDiGraph()
        graph.add_node("promiscuous", layer="Proteome", regulated=1, responsive=True, measured=True)
        graph.add_node("specific", layer="Proteome", regulated=1, responsive=True, measured=True)
        graph.add_node("cofactor", layer="Metabolome", regulated=1, responsive=True, measured=True)
        for i in range(n_reactions):
            graph.add_node(f"R{i}", layer="Reactions")
            graph.add_edge("promiscuous", f"R{i}", edge_type="catalysis", sign=1)
            graph.add_edge("cofactor", f"R{i}", edge_type="substrate", sign=1)
        graph.add_node("Rspecific", layer="Reactions")
        graph.add_edge("specific", "Rspecific", edge_type="catalysis", sign=1)
        graph.add_edge("cofactor", "Rspecific", edge_type="substrate", sign=1)
        return graph

    def test_one_enzyme_cannot_fill_the_figure(self):
        from transnet.visualization import transomic_backbone
        backbone = transomic_backbone(self._graph(), max_reactions=6,
                                      max_reactions_per_enzyme=3)
        reactions = [n for n, d in backbone.nodes(data=True) if d.get("layer") == "Reactions"]
        assert len(reactions) == 4                 # 3 promiscuous + the other enzyme's
        assert "Rspecific" in reactions

    def test_the_cap_can_be_lifted(self):
        from transnet.visualization import transomic_backbone
        backbone = transomic_backbone(self._graph(), max_reactions=6,
                                      max_reactions_per_enzyme=0)
        reactions = [n for n, d in backbone.nodes(data=True) if d.get("layer") == "Reactions"]
        assert len(reactions) == 6


class TestRegulationAxesFigure:
    def test_builds_from_a_summary(self, mapped):
        summary = regulation_axis_summary(reaction_regulation_table(mapped))
        if summary.empty:
            pytest.skip("no regulated reactions in this configuration")
        assert plot_regulation_axes(summary).axes

    def test_empty_summary_raises_a_helpful_error(self):
        with pytest.raises(ValueError, match="map_omics_to_network"):
            plot_regulation_axes(pd.DataFrame())


class TestLayerConnectivityFigure:
    def test_builds_from_connectivity(self, standard_network):
        assert plot_layer_connectivity(
            cross_layer_connectivity(standard_network)
        ).axes

    def test_empty_matrix_raises(self):
        import networkx as nx
        with pytest.raises(ValueError, match="empty"):
            plot_layer_connectivity(cross_layer_connectivity(nx.MultiDiGraph()))


def test_every_edge_type_has_a_line_style():
    """A new edge type must not silently fall back to an untyped default."""
    missing = sorted(set(EDGE_TYPES) - set(EDGE_STYLES))
    assert not missing, f"edge types with no line style: {missing}"


class TestLayerOrderInFigures:
    """Figures read top to bottom as Transcriptome, Proteome, Reactions,
    Metabolome. They showed Proteome on top while the regulatory order
    (transcription factors above genes) was reused for display."""

    def test_stacked_network_is_drawn_in_display_order(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from transnet import load_example_network
        from transnet.visualization import plot_transomic_network

        figure = plot_transomic_network(load_example_network())
        labels = {
            text.get_text(): text.get_position()[1]
            for text in figure.axes[0].texts
            if text.get_text() in ("Transcriptome", "Proteome", "Reactions", "Metabolome")
        }
        plt.close(figure)
        top_to_bottom = sorted(labels, key=labels.get, reverse=True)
        assert top_to_bottom == ["Transcriptome", "Proteome", "Reactions", "Metabolome"]

    def test_connectivity_matrix_is_in_display_order(self):
        from transnet import cross_layer_connectivity, load_example_network
        matrix = cross_layer_connectivity(load_example_network())["matrix"]
        assert list(matrix.index) == ["Transcriptome", "Proteome", "Reactions", "Metabolome"]


class TestBackboneKeepsTheSignalingLayer:
    """``transomic_backbone`` walks up from a selected enzyme to its sites.

    It selects around reactions and keeps each reaction's enzymes, their genes
    and their metabolites. It never walked *up* past the Proteome, so a network
    with a Signaling layer produced a figure with no Signaling layer in it, and
    a study that measured a phosphoproteome could not see it in any figure.
    """

    @staticmethod
    def _signalling_network():
        from transnet import (load_example_network, load_example_omics,
                              load_example_phosphoproteomics, map_omics_to_network)

        graph = load_example_network(signaling=True)
        tables = dict(load_example_omics())
        tables["Signaling"] = load_example_phosphoproteomics()
        map_omics_to_network(graph, tables, id_column="id",
                             log2fc_column="log2FC", qvalue_column="padj")
        return graph

    def _layers(self, graph):
        return {data.get("layer") for _, data in graph.nodes(data=True)}

    def test_the_signaling_layer_reaches_the_figure(self):
        from transnet.visualization import transomic_backbone

        graph = self._signalling_network()
        assert "Signaling" in self._layers(graph), "fixture has no Signaling layer"
        assert "Signaling" in self._layers(transomic_backbone(graph, max_reactions=8))

    def test_it_can_be_switched_off(self):
        from transnet.visualization import transomic_backbone

        selected = transomic_backbone(self._signalling_network(), max_reactions=8,
                                      max_signalling_per_enzyme=0)
        assert "Signaling" not in self._layers(selected)

    def test_a_network_without_signalling_is_unaffected(self):
        from transnet import (load_example_network, load_example_omics,
                              map_omics_to_network)
        from transnet.visualization import transomic_backbone

        graph = load_example_network()
        map_omics_to_network(graph, load_example_omics(), id_column="id",
                             log2fc_column="log2FC", qvalue_column="padj")
        selected = transomic_backbone(graph, max_reactions=8)
        assert selected.number_of_nodes() > 0
        assert "Signaling" not in self._layers(selected)


class TestCommunityFigure:
    def _figure(self, mapped, **kwargs):
        from transnet.analysis import detect_communities
        from transnet.visualization import plot_community_network

        communities = detect_communities(mapped, method="louvain")
        return plot_community_network(mapped, communities, **kwargs)

    def test_the_drawn_community_is_labelled(self, mapped):
        drawing = self._figure(mapped, label_top=5).axes[1]
        assert 0 < len(drawing.texts) <= 5

    def test_it_has_a_legend(self, mapped):
        drawing = self._figure(mapped).axes[1]
        entries = [t.get_text() for t in drawing.get_legend().get_texts()]
        assert "increased" in entries and "decreased" in entries

    def test_labels_can_be_switched_off(self, mapped):
        assert not self._figure(mapped, label_top=0).axes[1].texts


class TestControversialFigure:
    """It used to call a helper that did not exist, which went unnoticed
    because the walkthrough only ever passed it a table with no controversy."""

    @staticmethod
    def _mapped():
        from transnet import load_example_network, load_example_omics

        graph = load_example_network()
        map_omics_to_network(graph, load_example_omics("insulin_sensitive"), id_column="id",
                             log2fc_column="log2FC", qvalue_column="padj")
        return graph

    def test_draws_the_controversial_reactions(self):
        from transnet.visualization import plot_controversial_reactions

        graph = self._mapped()
        table = reaction_regulation_table(graph)
        assert table["controversial"].any()
        figure = plot_controversial_reactions(graph, table)
        assert figure.axes[0].patches

    def test_contributions_are_signed_by_effect(self):
        from transnet.visualization.findings import reaction_contributions

        graph = self._mapped()
        rows = reaction_contributions(graph, "R00299").set_index("role")
        product_or_inhibitor = rows.loc[rows.index.isin(["product", "inhibitor"])]
        assert (product_or_inhibitor["contribution"]
                == -product_or_inhibitor["log2fc"]).all()
        assert set(rows["axis"]) == {"enzyme", "metabolite"}
