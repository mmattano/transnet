"""transcription-factor activity transcription factor activity and transcript-protein concordance expression concordance.

Hand-built networks with a known answer, so the tests check the inference and
not just that a table comes back.
"""

import networkx as nx
import pandas as pd
import pytest

from transnet import expression_concordance, transcription_factor_activity


def _gene(graph, node, regulated=0, log2fc=None, measured=True):
    graph.add_node(node, layer="Transcriptome", name=node)
    if measured:
        graph.nodes[node].update(
            measured=True, regulated=regulated,
            responsive=bool(regulated),
            log2fc=log2fc if log2fc is not None else float(regulated),
        )


def _protein(graph, node, regulated=0, log2fc=None, measured=True, symbol=None):
    graph.add_node(node, layer="Proteome", name=node, symbol=symbol or node)
    if measured:
        graph.nodes[node].update(
            measured=True, regulated=regulated,
            responsive=bool(regulated),
            log2fc=log2fc if log2fc is not None else float(regulated),
        )


def _tf_network(up_targets=10, down_targets=0, background=60):
    """One factor whose targets are all responsive, in a quiet background."""
    graph = nx.MultiDiGraph()
    _protein(graph, "TF1", measured=False, symbol="Srebf1")
    _protein(graph, "TF2", measured=False, symbol="Quiet")
    for i in range(up_targets):
        _gene(graph, f"up{i}", regulated=1)
        graph.add_edge("TF1", f"up{i}", edge_type="transcriptional_regulation", sign=0)
    for i in range(down_targets):
        _gene(graph, f"down{i}", regulated=-1)
        graph.add_edge("TF1", f"down{i}", edge_type="transcriptional_regulation", sign=0)
    for i in range(background):
        _gene(graph, f"bg{i}")
        if i < 10:
            graph.add_edge("TF2", f"bg{i}", edge_type="transcriptional_regulation", sign=0)
    return graph


class TestTranscriptionFactorActivity:

    def test_a_factor_whose_targets_all_rise_is_inferred_active(self):
        table = transcription_factor_activity(_tf_network())
        top = table.iloc[0]
        assert top["factor"] == "TF1" and top["name"] == "Srebf1"
        assert top["n_responsive_targets"] == 10 and top["n_up"] == 10
        assert top["q_value"] < 0.05
        assert top["inferred_activity"] == 1

    def test_a_factor_whose_targets_are_quiet_is_not_implicated(self):
        table = transcription_factor_activity(_tf_network()).set_index("factor")
        assert table.loc["TF2", "n_responsive_targets"] == 0
        assert table.loc["TF2", "inferred_activity"] == 0

    def test_mixed_directions_implicate_without_a_direction(self):
        """Targets split up and down: the factor matters, its sign is unclear."""
        table = transcription_factor_activity(
            _tf_network(up_targets=5, down_targets=5)
        ).set_index("factor")
        assert table.loc["TF1", "q_value"] < 0.05
        assert table.loc["TF1", "inferred_activity"] == 0

    def test_the_factors_own_level_is_reported_separately(self):
        graph = _tf_network()
        graph.nodes["TF1"].update(measured=True, regulated=0, log2fc=0.0)
        row = transcription_factor_activity(graph).set_index("factor").loc["TF1"]
        # active through its targets, unchanged itself: post-translational
        assert row["inferred_activity"] == 1 and row["factor_regulated"] == 0

    def test_unmeasured_factor_level_is_none_not_zero(self):
        row = transcription_factor_activity(_tf_network()).set_index("factor").loc["TF1"]
        assert row["factor_regulated"] is None

    def test_min_confidence_drops_weak_binding(self):
        graph = _tf_network()
        for _, _, data in graph.edges(data=True):
            data["confidence"] = 120.0
        assert transcription_factor_activity(graph, min_confidence=500).empty

    def test_min_confidence_drops_unscored_edges(self):
        """A network with no binding scores must not pass every threshold.

        Reading a network back from CSV gives missing scores as NaN, and every
        comparison against NaN is False, so a plain ``score < threshold`` test
        kept unscored edges at any threshold and quietly reported a filtered
        result that had not been filtered.
        """
        import math

        for missing in (None, float("nan"), math.nan, ""):
            graph = _tf_network()
            for _, _, data in graph.edges(data=True):
                data["confidence"] = missing
            assert transcription_factor_activity(graph, min_confidence=1).empty, \
                f"confidence={missing!r} passed a threshold it cannot meet"

    def test_unscored_edges_are_kept_when_no_threshold_is_asked_for(self):
        graph = _tf_network()
        for _, _, data in graph.edges(data=True):
            data["confidence"] = float("nan")
        assert not transcription_factor_activity(graph).empty

    def test_no_transcriptome_returns_an_empty_typed_frame(self):
        graph = nx.MultiDiGraph()
        _protein(graph, "P1", regulated=1)
        table = transcription_factor_activity(graph)
        assert table.empty and "inferred_activity" in table.columns

    def test_no_tf_edges_returns_an_empty_typed_frame(self):
        graph = nx.MultiDiGraph()
        for i in range(10):
            _gene(graph, f"g{i}", regulated=1)
        assert transcription_factor_activity(graph).empty


class TestExpressionConcordance:

    @staticmethod
    def _pairs():
        graph = nx.MultiDiGraph()
        cases = {
            "same": (1, 1),            # concordant
            "prot": (0, 1),            # protein only
            "tran": (-1, 0),           # transcript only
            "flip": (1, -1),           # discordant
            "flat": (0, 0),            # unchanged
        }
        for name, (g, p) in cases.items():
            _gene(graph, f"g_{name}", regulated=g)
            _protein(graph, f"p_{name}", regulated=p)
            graph.add_edge(f"g_{name}", f"p_{name}", edge_type="translation", sign=1)
        return graph

    def test_every_category_is_classified(self):
        result = expression_concordance(self._pairs())
        by_protein = result["table"].set_index("protein")["category"]
        assert by_protein["p_same"] == "concordant"
        assert by_protein["p_prot"] == "protein_only"
        assert by_protein["p_tran"] == "transcript_only"
        assert by_protein["p_flip"] == "discordant"
        assert by_protein["p_flat"] == "unchanged"
        assert result["counts"] == {
            "concordant": 1, "protein_only": 1, "transcript_only": 1,
            "discordant": 1, "unchanged": 1,
            "tested_difference": 0, "protein_beyond_transcript": 0,
        }

    def test_pairs_missing_a_measurement_are_skipped(self):
        graph = self._pairs()
        _gene(graph, "g_lone", regulated=1)
        _protein(graph, "p_lone", measured=False)
        graph.add_edge("g_lone", "p_lone", edge_type="translation", sign=1)
        result = expression_concordance(graph)
        assert "p_lone" not in set(result["table"]["protein"])

    def test_correlation_is_reported(self):
        result = expression_concordance(self._pairs())
        assert -1.0 <= result["correlation"] <= 1.0

    def test_no_proteome_gives_an_empty_result(self):
        graph = nx.MultiDiGraph()
        _gene(graph, "g", regulated=1)
        result = expression_concordance(graph)
        assert result["table"].empty
        assert sum(result["counts"].values()) == 0


class TestProteinBeyondTranscript:
    """Category counts depend on the transcript test's power; the difference
    test does not."""

    @staticmethod
    def _pair(gene_fc, gene_se, protein_fc, protein_se, gene_regulated=0, protein_regulated=1):
        graph = nx.MultiDiGraph()
        _gene(graph, "g", regulated=gene_regulated, log2fc=gene_fc)
        _protein(graph, "p", regulated=protein_regulated, log2fc=protein_fc)
        graph.nodes["g"]["se"] = gene_se
        graph.nodes["p"]["se"] = protein_se
        graph.add_edge("g", "p", edge_type="translation", sign=1)
        return expression_concordance(graph)

    def test_a_protein_far_beyond_a_flat_transcript_is_called(self):
        result = self._pair(0.0, 0.1, 1.5, 0.1)
        assert result["counts"]["protein_beyond_transcript"] == 1

    def test_an_underpowered_transcript_is_not_evidence(self):
        """Transcript not significant (large se) but moved as far as the protein."""
        result = self._pair(1.4, 0.9, 1.5, 0.1)
        row = result["table"].iloc[0]
        assert row["category"] == "protein_only"
        assert not row["protein_beyond_transcript"]

    def test_without_standard_errors_nothing_is_tested(self):
        graph = nx.MultiDiGraph()
        _gene(graph, "g", regulated=0, log2fc=0.0)
        _protein(graph, "p", regulated=1, log2fc=1.5)
        graph.add_edge("g", "p", edge_type="translation", sign=1)
        assert expression_concordance(graph)["counts"]["tested_difference"] == 0
