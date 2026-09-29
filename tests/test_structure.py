"""Motifs, bottlenecks and the convergence null.

These are the analyses that read the network itself rather than data through
it, so the tests build networks whose answer is known by construction.
"""

import networkx as nx
import numpy as np
import pandas as pd
import pytest

from transnet import (
    convergence_significance,
    regulatory_motifs,
    structural_vulnerability,
)


def _reaction(graph, reaction, enzyme, substrate, product):
    graph.add_node(reaction, layer="Reactions", name=reaction)
    graph.add_node(enzyme, layer="Proteome", name=enzyme, symbol=enzyme)
    graph.add_node(substrate, layer="Metabolome", name=substrate)
    graph.add_node(product, layer="Metabolome", name=product)
    graph.add_edge(enzyme, reaction, edge_type="catalysis", sign=1)
    graph.add_edge(substrate, reaction, edge_type="substrate", sign=1)
    graph.add_edge(reaction, product, edge_type="product", sign=1)


class TestRegulatoryMotifs:

    def test_product_inhibition_is_found_without_being_named(self):
        """Hexokinase makes G6P, and G6P inhibits hexokinase. The motif is
        defined by direction and sign, so it is found wherever it occurs."""
        graph = nx.MultiDiGraph()
        _reaction(graph, "R1", "Hk1", "glucose", "G6P")
        graph.add_edge("G6P", "R1", edge_type="allosteric_inhibition", sign=-1)

        motifs = regulatory_motifs(graph)["motifs"]
        inhibition = motifs[motifs["motif"] == "product_inhibition"]
        assert len(inhibition) == 1
        assert inhibition.iloc[0]["metabolite"] == "G6P"
        assert inhibition.iloc[0]["enzyme"] == "Hk1"
        assert inhibition.iloc[0]["sign_product"] == -1

    def test_product_activation_is_not_called_inhibition(self):
        graph = nx.MultiDiGraph()
        _reaction(graph, "R1", "E1", "A", "B")
        graph.add_edge("B", "R1", edge_type="allosteric_activation", sign=1)
        counts = regulatory_motifs(graph)["counts"]
        assert counts.get("product_activation") == 1
        assert "product_inhibition" not in counts

    def test_currency_metabolites_do_not_manufacture_feedback(self):
        """ATP is a product of half the network; counting it would make every
        reaction look self-regulating."""
        graph = nx.MultiDiGraph()
        _reaction(graph, "R1", "E1", "A", "C00002")          # ATP
        graph.add_edge("C00002", "R1", edge_type="allosteric_inhibition", sign=-1)

        assert regulatory_motifs(graph)["counts"] == {}
        kept = regulatory_motifs(graph, exclude_currency=False)["counts"]
        assert kept.get("product_inhibition") == 1

    def test_feed_forward_needs_two_targets_on_one_reaction(self):
        graph = nx.MultiDiGraph()
        _reaction(graph, "R1", "P1", "A", "B")
        graph.add_node("TF", layer="Proteome", name="TF", symbol="TF")
        for gene, protein in (("g1", "P1"), ("g2", "P2")):
            graph.add_node(gene, layer="Transcriptome", name=gene)
            graph.add_edge("TF", gene, edge_type="transcriptional_regulation", sign=0)
            graph.add_node(protein, layer="Proteome", name=protein, symbol=protein)
            graph.add_edge(gene, protein, edge_type="translation", sign=1)
        graph.add_edge("P2", "R1", edge_type="catalysis", sign=1)

        motifs = regulatory_motifs(graph)["motifs"]
        assert (motifs["motif"] == "feed_forward").sum() == 1

    def test_responsive_only_narrows_a_census_to_this_response(self):
        graph = nx.MultiDiGraph()
        _reaction(graph, "R1", "E1", "A", "B")
        graph.add_edge("B", "R1", edge_type="allosteric_inhibition", sign=-1)

        assert regulatory_motifs(graph)["counts"].get("product_inhibition") == 1
        assert regulatory_motifs(graph, responsive_only=True)["counts"] == {}

        graph.nodes["B"]["regulated"] = 1
        graph.nodes["E1"]["regulated"] = 1
        graph.nodes["A"]["regulated"] = -1
        assert regulatory_motifs(graph, responsive_only=True)["counts"] == {
            "product_inhibition": 1}


class TestStructuralVulnerability:

    def test_the_molecule_joining_two_halves_is_found(self):
        graph = nx.MultiDiGraph()
        for node in "abcde":
            graph.add_node(node, layer="Proteome", name=node, regulated=1)
        # a - b - BRIDGE - d - e
        graph.add_node("BRIDGE", layer="Proteome", name="BRIDGE", regulated=1)
        for u, v in (("a", "b"), ("b", "BRIDGE"), ("BRIDGE", "d"), ("d", "e")):
            graph.add_edge(u, v, edge_type="protein_interaction", sign=0)

        table = structural_vulnerability(graph)
        assert table.iloc[0]["node"] == "BRIDGE"
        assert table.iloc[0]["fragments"] == 2
        assert 0 < table.iloc[0]["largest_loss"] < 1

    def test_a_ring_has_no_single_point_of_failure(self):
        graph = nx.MultiDiGraph()
        nodes = list("abcdef")
        for node in nodes:
            graph.add_node(node, layer="Proteome", regulated=1)
        for u, v in zip(nodes, nodes[1:] + nodes[:1]):
            graph.add_edge(u, v, edge_type="protein_interaction", sign=0)
        assert structural_vulnerability(graph).empty

    def test_unresponsive_molecules_are_excluded_by_default(self):
        graph = nx.MultiDiGraph()
        for node in "abc":
            graph.add_node(node, layer="Proteome", regulated=0)
        graph.add_edge("a", "b", edge_type="protein_interaction")
        graph.add_edge("b", "c", edge_type="protein_interaction")
        assert structural_vulnerability(graph).empty
        assert not structural_vulnerability(graph, responsive_only=False).empty


class TestConvergenceSignificance:

    @staticmethod
    def _network(convergent_reactions, total_reactions=20):
        """Reactions whose enzyme and metabolite both changed, among others."""
        graph = nx.MultiDiGraph()
        for i in range(total_reactions):
            changed = i < convergent_reactions
            _reaction(graph, f"R{i}", f"E{i}", f"S{i}", f"P{i}")
            graph.nodes[f"E{i}"]["regulated"] = 1 if changed else 0
            graph.nodes[f"S{i}"]["regulated"] = -1 if changed else 0
            graph.nodes[f"P{i}"]["regulated"] = 0
        return graph

    def test_planted_convergence_beats_the_null(self):
        result = convergence_significance(self._network(14), n_randomisations=200)
        assert result["observed"] == 14
        assert result["observed"] > result["null_mean"]
        assert result["p_value"] < 0.05

    def test_scattered_changes_do_not(self):
        graph = self._network(0)
        rng = np.random.default_rng(0)
        molecules = [n for n, d in graph.nodes(data=True) if d.get("layer") != "Reactions"]
        for node in rng.choice(molecules, size=len(molecules) // 3, replace=False):
            graph.nodes[node]["regulated"] = 1
        result = convergence_significance(graph, n_randomisations=200)
        assert result["p_value"] > 0.05

    def test_the_null_distribution_is_returned_for_drawing(self):
        result = convergence_significance(self._network(5), n_randomisations=50)
        assert len(result["null"]) == 50
