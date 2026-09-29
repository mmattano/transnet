"""Summaries that keep a result honest.

Three things every study in this package needs and none of the analysis
functions should own: a chance baseline for directional predictions, a way to
score traced paths once per molecule rather than once per path, and a hub
ranking that shows what ChIP-Atlas binding does to it.
"""

from typing import Tuple

import pandas as pd

__all__ = ["versus_chance", "path_verdicts", "hub_rankings"]


def versus_chance(successes: int, trials: int) -> str:
    """'62%, better than the 50% of a coin flip (p = ...)' -- or not.

    A share of correct directional predictions means nothing without the
    baseline: with no information, half of all up/down calls are right.
    """
    from scipy.stats import binomtest

    if not trials:
        return "no predictions to test"
    share = successes / trials
    p_value = binomtest(successes, trials, 0.5).pvalue
    if p_value > 0.05:
        verdict = "no better than the 50% expected by chance"
    elif share > 0.5:
        verdict = "better than the 50% expected by chance"
    else:
        verdict = "worse than the 50% expected by chance"
    return f"{share:.0%}, {verdict} (binomial p = {p_value:.2g})"


def path_verdicts(paths):
    """Score traced paths once per target molecule.

    Paths are not independent: one metabolite can be reached by dozens of
    paths that share most of their steps, so a rate over paths counts the same
    claim many times and its binomial p-value is far too small. This scores
    each molecule once, on the direction most of its paths predict.

    Returns
    -------
    tuple
        ``(summary, agree, tested)`` -- the per-target table from
        :func:`~transnet.path_consistency_summary`, how many molecules the
        network's paths call correctly, and how many it calls at all.
    """
    from transnet import path_consistency_summary

    # ``n_paths`` exists only in the per-target summary; the traced paths and the
    # summary both carry ``target``, so that column cannot tell them apart.
    if "n_paths" in getattr(paths, "columns", []):
        raise ValueError(
            "path_verdicts takes the traced paths, not the per-target summary; "
            "pass the output of trace_regulatory_paths()"
        )
    summary = path_consistency_summary(paths)
    if summary.empty:
        return summary, 0, 0
    tested = summary[(summary["observed"] != 0) & (summary["predicted"] != 0)]
    return summary, int(tested["agrees"].sum()), len(tested)


def hub_rankings(graph, top_percent: float = 2.0):
    """Trans-omic hubs two ways, and how much ChIP-Atlas drives the first.

    Hubs are ranked by degree within the responsive network (Morita et al.).
    ChIP-Atlas lists thousands of targets for any protein that has been
    ChIP-sequenced, chromatin factors included, so the top of that ranking is
    often decided by binding breadth. The second ranking drops
    ``transcriptional_regulation`` edges to show the hubs of the metabolic and
    protein relationships underneath. Both are reported; neither is hidden.

    Returns ``(hubs, hubs_without_binding, top_hub_binding_share)``.
    """
    from transnet import responsive_subnetwork, transomic_hubs

    responsive = responsive_subnetwork(graph)
    hubs = transomic_hubs(responsive, top_percent=top_percent)
    without = responsive.copy()
    without.remove_edges_from([(u, v, k) for u, v, k, d in responsive.edges(keys=True, data=True)
                               if d.get("edge_type") == "transcriptional_regulation"])
    hubs_without = transomic_hubs(without, top_percent=top_percent)
    share = float("nan")
    if not hubs.empty:
        top = hubs.sort_values("cross_layer_degree", ascending=False).iloc[0]["node"]
        edges = [d.get("edge_type") for _, _, d in responsive.edges(top, data=True)] + \
                [d.get("edge_type") for _, _, d in responsive.in_edges(top, data=True)]
        if edges:
            share = sum(e == "transcriptional_regulation" for e in edges) / len(edges)
    return hubs, hubs_without, share
