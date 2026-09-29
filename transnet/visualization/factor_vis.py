"""Figures for factor analysis: every factor, and what it is associated with.

A factor's number is a label, not a rank -- NMF factors come in no order -- so
these figures never select factors by number. They show all of them side by
side with the three things that decide which matter: which part of the design
they follow, how much of the data they explain, and whether their strongest
features are connected on the trans-omic network.
"""

from typing import List, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from transnet.visualization.palette import (
    BASELINE, GRID, INK, LAYER_COLORS, MUTED, SECONDARY, SURFACE, style_axes,
)

__all__ = ["plot_factor_overview", "plot_factor_scores"]


def _stars(q):
    if pd.isna(q):
        return ""
    return "***" if q <= 0.001 else "**" if q <= 0.01 else "*" if q <= 0.05 else ""


def plot_factor_overview(association: pd.DataFrame, variance: pd.DataFrame,
                         coherence: Optional[pd.DataFrame] = None,
                         title: str = "Every factor: what it follows, explains and connects"):
    """One row per factor, three aligned panels.

    1. Partial eta-squared for each design term from
       :func:`~transnet.analysis.factor.factor_design_association`, with
       Benjamini-Hochberg stars.
    2. Variance accounted for across all layers, from
       :func:`~transnet.analysis.factor.factor_variance_explained`.
    3. Network coherence: how many more cross-layer links join the factor's
       top features than random features, from
       :func:`~transnet.analysis.factor.factor_network_coherence`.
    """
    factors = list(dict.fromkeys(association["factor"]))
    terms = list(dict.fromkeys(association["term"]))
    eta = association.pivot(index="factor", columns="term", values="partial_eta_squared").loc[factors, terms]
    q = association.pivot(index="factor", columns="term", values="q_value").loc[factors, terms]

    panels = 3 if coherence is not None else 2
    widths = [1.1 * len(terms), 1.6] + ([1.6] if coherence is not None else [])
    figure, axes = plt.subplots(1, panels, figsize=(2.0 + sum(widths) * 1.35, 0.55 * len(factors) + 2.1),
                                gridspec_kw={"width_ratios": widths}, sharey=True)
    heat = axes[0]
    heat.imshow(eta.to_numpy(dtype=float), cmap="Purples", vmin=0, vmax=max(0.6, float(np.nanmax(eta.to_numpy()))),
                aspect="auto")
    for i, factor in enumerate(factors):
        for j, term in enumerate(terms):
            value = eta.loc[factor, term]
            dark = value > 0.35
            heat.text(j, i, f"{value:.2f}{_stars(q.loc[factor, term])}", ha="center", va="center",
                      fontsize=9, color=SURFACE if dark else INK)
    heat.set_xticks(range(len(terms)))
    heat.set_xticklabels(terms, fontsize=9, color=SECONDARY)
    heat.set_yticks(range(len(factors)))
    heat.set_yticklabels(factors, fontsize=9.5, color=INK)
    heat.tick_params(length=0)
    for side in heat.spines.values():
        side.set_visible(False)
    heat.set_title("partial eta-squared\n(* q<=0.05, ** q<=0.01, *** q<=0.001)",
                   fontsize=9, color=SECONDARY, loc="left")

    share = variance[variance["layer"] == "all layers"].set_index("factor").reindex(factors)["variance_accounted"]
    bars = axes[1]
    style_axes(bars, grid_axis="x")
    bars.barh(range(len(factors)), share.to_numpy(dtype=float), color=MUTED, height=0.6,
              edgecolor=SURFACE, linewidth=1)
    for i, value in enumerate(share):
        bars.text(value + 0.004, i, f"{value:.1%}", va="center", fontsize=8.5, color=INK)
    bars.set_xlim(0, max(float(share.max()) * 1.35, 0.05))
    bars.set_title("variance accounted\n(all layers)", fontsize=9, color=SECONDARY, loc="left")
    bars.tick_params(axis="y", length=0)

    if coherence is not None:
        panel = axes[2]
        style_axes(panel, grid_axis="x")
        table = coherence.set_index("factor").reindex(factors)
        fold = table["fold_enrichment"].fillna(0).to_numpy(dtype=float)
        significant = (table["q_value"] <= 0.05).to_numpy()
        panel.barh(range(len(factors)), fold, color=[LAYER_COLORS["Reactions"] if s else MUTED
                                                     for s in significant],
                   height=0.6, edgecolor=SURFACE, linewidth=1)
        panel.axvline(1, color=BASELINE, linewidth=1)
        for i, (value, links, qv) in enumerate(zip(fold, table["n_links"], table["q_value"])):
            panel.text(value + 0.03, i, f"{value:.1f}x ({int(links)} links){_stars(qv)}",
                       va="center", fontsize=8.5, color=INK)
        panel.set_xlim(0, max(fold.max() * 1.6, 2))
        panel.set_title("network coherence\n(links vs random features)", fontsize=9,
                        color=SECONDARY, loc="left")
        panel.tick_params(axis="y", length=0)

    figure.suptitle(title, x=0.01, ha="left", fontsize=12.5, fontweight="bold", color=INK)
    figure.tight_layout()
    return figure


def plot_factor_scores(factors: pd.DataFrame, design: pd.DataFrame, x: str,
                       x_order: Sequence[str], hue: Optional[str] = None,
                       association: Optional[pd.DataFrame] = None, columns: int = 3,
                       title: str = "Factor scores across the design"):
    """Small multiples, one per factor: scores by ``x``, split by ``hue``.

    The two ``hue`` groups are told apart by marker fill -- solid and hollow --
    and a mean line each, so no colour is needed for them. Each panel is
    headed with its strongest design term from ``association``.
    """
    common = factors.index.intersection(design.index)
    names = list(factors.columns)
    rows = int(np.ceil(len(names) / columns))
    figure, axes = plt.subplots(rows, columns, figsize=(3.6 * columns, 2.9 * rows), squeeze=False)
    groups = [None] if hue is None else sorted(design.loc[common, hue].astype(str).unique())
    styles = [dict(facecolor=INK, edgecolor=INK), dict(facecolor=SURFACE, edgecolor=SECONDARY)]
    positions = {level: i for i, level in enumerate(x_order)}
    rng = np.random.default_rng(0)
    for k, factor in enumerate(names):
        axes_k = axes[k // columns][k % columns]
        style_axes(axes_k, grid_axis="y")
        for g, group in enumerate(groups):
            index = common if group is None else common[design.loc[common, hue].astype(str) == group]
            levels = design.loc[index, x].astype(str)
            xs = levels.map(positions).to_numpy(dtype=float)
            offset = (g - (len(groups) - 1) / 2) * 0.22
            ys = factors.loc[index, factor].to_numpy(dtype=float)
            axes_k.scatter(xs + offset + rng.uniform(-0.05, 0.05, len(xs)), ys, s=16, linewidths=0.9,
                           label=group, zorder=3, **styles[g % 2])
            means = pd.Series(ys, index=levels.to_numpy()).groupby(level=0).mean().reindex(x_order)
            axes_k.plot(np.arange(len(x_order)) + offset, means.to_numpy(), color=INK if g % 2 == 0 else SECONDARY,
                        linewidth=1.4, linestyle="-" if g % 2 == 0 else "--", zorder=2)
        axes_k.set_xticks(range(len(x_order)))
        axes_k.set_xticklabels(x_order, fontsize=8)
        heading = factor
        if association is not None:
            rows_f = association[association["factor"] == factor].sort_values("partial_eta_squared",
                                                                              ascending=False)
            if not rows_f.empty:
                top = rows_f.iloc[0]
                heading += (f"  -  {top['term']} eta2={top['partial_eta_squared']:.2f}"
                            f"{_stars(top['q_value'])}")
        axes_k.set_title(heading, loc="left", fontsize=9.5, color=INK)
        if k == 0 and hue is not None:
            axes_k.legend(frameon=False, fontsize=8, labelcolor=SECONDARY, loc="best")
    for k in range(len(names), rows * columns):
        axes[k // columns][k % columns].axis("off")
    figure.suptitle(title, x=0.01, ha="left", fontsize=12.5, fontweight="bold", color=INK)
    figure.tight_layout()
    return figure


def plot_factor_network(graph, loadings, factor: str,
                        id_maps: Optional[dict] = None, top_n: int = 15,
                        nodes: Optional[Sequence[str]] = None, max_nodes: int = 90,
                        title: Optional[str] = None):
    """A factor's features, drawn where they connect on the network.

    The overview figures say whether a factor is coherent; this says what it
    is *about*. Pass ``nodes`` -- the linked top features that
    :func:`~transnet.analysis.factors.factor_network_coherence` returns for the
    factor -- to draw exactly the part of the factor that the network joins.
    Without it, each layer's ``top_n`` strongest features are drawn.

    The reactions joining the selected molecules are added, or the figure is a
    row of disconnected dots, and the whole is capped at ``max_nodes``
    (strongest loadings first) so a large factor stays legible.
    """
    from transnet.visualization.transomics_vis import plot_transomic_network

    weight: dict = {}
    for layer, frame in loadings.items():
        if factor not in frame.columns:
            continue
        mapping = (id_maps or {}).get(layer, {})
        for feature, value in frame[factor].items():
            node = mapping.get(str(feature), str(feature))
            weight[node] = max(weight.get(node, 0.0), float(value))

    if nodes is not None:
        selected = [node for node in nodes if node in graph]
    else:
        selected = []
        for layer, frame in loadings.items():
            if factor not in frame.columns:
                continue
            mapping = (id_maps or {}).get(layer, {})
            ranked = frame[factor].sort_values(ascending=False).head(top_n)
            selected += [mapping.get(str(f), str(f)) for f in ranked.index]
        selected = [node for node in dict.fromkeys(selected) if node in graph]

    if not selected:
        raise ValueError(f"none of {factor}'s features are in the network")

    molecules = [n for n in selected if graph.nodes[n].get("layer") != "Reactions"]
    molecules = sorted(molecules, key=lambda n: -weight.get(n, 0.0))[:max_nodes]
    chosen = set(molecules)
    connectors = {
        node for node, data in graph.nodes(data=True)
        if data.get("layer") == "Reactions"
        and sum(1 for neighbour in set(graph.predecessors(node)) | set(graph.successors(node))
                if neighbour in chosen) >= 2
    }
    # keep only the most-joined reactions, so the reaction row stays readable
    ranked_reactions = sorted(
        connectors,
        key=lambda r: -sum(1 for n in set(graph.predecessors(r)) | set(graph.successors(r))
                           if n in chosen))[:max(20, max_nodes // 3)]
    subgraph = graph.subgraph(chosen | set(ranked_reactions))
    return plot_transomic_network(
        subgraph, title=title or f"{factor}: its top features, connected on the network")


__all__.append("plot_factor_network")
