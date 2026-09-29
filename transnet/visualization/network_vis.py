"""Ordinary graph figures, drawn on the typed trans-omic network.

These are the counterparts to :mod:`transnet.analysis.network_analysis`: the
census, the communities and a values heatmap. They say nothing about
regulation -- :func:`~transnet.visualization.plot_transomic_network` is the
figure for that -- but they answer the questions a reader asks first. Does
this network have one component or twenty? Do its communities span layers, or
does each sit inside one?

Everything here reads the current schema (``layer``, ``edge_type``, ``sign``,
``regulated``) and the shared palette, so a figure from this module sits
beside the trans-omic ones without a colour clashing.
"""

from typing import Dict, Optional, Sequence

import logging

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

from transnet.biology.schema import order_layers
from transnet.visualization.palette import (
    BASELINE, DOWN, GRID, INK, LAYER_COLORS, MUTED, SECONDARY, SURFACE, UP,
    style_axes,
)

logger = logging.getLogger(__name__)

__all__ = [
    "plot_network_metrics",
    "plot_community_network",
    "plot_values_heatmap",
]

_UNKNOWN = MUTED


def _layer(graph, node) -> str:
    return str(graph.nodes[node].get("layer") or "Unknown")


def plot_network_metrics(graph, title: str = "What this network is made of"):
    """Degree distribution, component sizes and nodes per layer.

    The census that belongs beside
    :func:`~transnet.analysis.network_analysis.compute_network_statistics`: a
    trans-omic network with one giant component and a long-tailed degree
    distribution behaves very differently from one that fragmented because a
    layer failed to map.
    """
    undirected = graph.to_undirected() if graph.is_directed() else graph
    degrees = np.array([d for _, d in graph.degree()], dtype=float)
    components = sorted((len(c) for c in nx.connected_components(undirected)), reverse=True)
    per_layer = pd.Series([_layer(graph, n) for n in graph.nodes()]).value_counts()
    per_layer = per_layer.reindex(order_layers(per_layer.index)).dropna()

    figure, axes = plt.subplots(1, 3, figsize=(13, 3.6))

    style_axes(axes[0], grid_axis="y")
    positive = degrees[degrees > 0]
    if positive.size:
        bins = np.logspace(0, np.log10(positive.max() + 1), 30)
        axes[0].hist(positive, bins=bins, color=MUTED, edgecolor=SURFACE, linewidth=0.6)
        axes[0].set_xscale("log")
    axes[0].set_xlabel("degree")
    axes[0].set_ylabel("molecules")
    axes[0].set_title(f"median {np.median(degrees):.0f}, max {degrees.max():.0f}",
                      fontsize=9.5, color=SECONDARY, loc="left")

    style_axes(axes[1], grid_axis="y")
    shown = components[:12]
    axes[1].bar(range(len(shown)), shown, color=MUTED, edgecolor=SURFACE, linewidth=0.6)
    axes[1].set_yscale("log")
    axes[1].set_xlabel("component, largest first")
    axes[1].set_ylabel("molecules")
    largest = shown[0] / graph.number_of_nodes() if shown else 0
    axes[1].set_title(f"{len(components)} components; largest holds {largest:.0%}",
                      fontsize=9.5, color=SECONDARY, loc="left")

    style_axes(axes[2], grid_axis="x")
    axes[2].barh(range(len(per_layer)), per_layer.to_numpy(),
                 color=[LAYER_COLORS.get(l, _UNKNOWN) for l in per_layer.index],
                 edgecolor=SURFACE, linewidth=0.8)
    axes[2].set_yticks(range(len(per_layer)))
    axes[2].set_yticklabels(per_layer.index, fontsize=9)
    axes[2].invert_yaxis()
    axes[2].set_xlabel("molecules")
    axes[2].set_title("per layer", fontsize=9.5, color=SECONDARY, loc="left")

    figure.suptitle(title, x=0.01, ha="left", fontsize=12.5, fontweight="bold", color=INK)
    figure.tight_layout()
    return figure


def plot_community_network(graph, communities: Dict[str, int], max_communities: int = 8,
                           title: str = "Do the communities span layers?"):
    """Community composition by layer, and the largest communities drawn.

    A community that lives inside one layer is a statement about how well that
    layer is annotated. A community holding genes, proteins and metabolites is
    a candidate module, and the question worth asking of
    :func:`~transnet.analysis.network_analysis.detect_communities`.
    """
    frame = pd.DataFrame({
        "node": list(communities),
        "community": [communities[n] for n in communities],
        "layer": [_layer(graph, n) for n in communities if n in graph],
    })
    if frame.empty:
        raise ValueError("no communities to draw")

    sizes = frame["community"].value_counts().head(max_communities)
    shown = frame[frame["community"].isin(sizes.index)]
    composition = (shown.groupby(["community", "layer"]).size()
                   .unstack(fill_value=0).loc[sizes.index])
    composition = composition[order_layers(composition.columns)]

    figure, axes = plt.subplots(1, 2, figsize=(12, 4.2),
                                gridspec_kw={"width_ratios": [1.15, 1]})

    style_axes(axes[0], grid_axis="x")
    left = np.zeros(len(composition))
    for layer in composition.columns:
        values = composition[layer].to_numpy(dtype=float)
        axes[0].barh(range(len(composition)), values, left=left, label=layer,
                     color=LAYER_COLORS.get(layer, _UNKNOWN),
                     edgecolor=SURFACE, linewidth=1.2)
        left += values
    axes[0].set_yticks(range(len(composition)))
    axes[0].set_yticklabels([f"community {c}" for c in composition.index], fontsize=9)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("molecules")
    spanning = int((composition > 0).sum(axis=1).gt(1).sum())
    axes[0].set_title(f"{spanning} of {len(composition)} shown communities span "
                      f"more than one layer", fontsize=9.5, color=SECONDARY, loc="left")
    axes[0].legend(frameon=False, fontsize=8, labelcolor=SECONDARY, ncol=2)

    biggest = sizes.index[0]
    members = [n for n in shown[shown["community"] == biggest]["node"] if n in graph][:150]
    sub = graph.subgraph(members)
    position = nx.spring_layout(sub.to_undirected() if sub.is_directed() else sub, seed=0)
    axes[1].set_axis_off()
    for u, v in (sub.to_undirected() if sub.is_directed() else sub).edges():
        axes[1].plot(*zip(position[u], position[v]), color=GRID, linewidth=0.6, zorder=1)
    for node in sub.nodes():
        axes[1].scatter(*position[node], s=42, zorder=2,
                        color=LAYER_COLORS.get(_layer(graph, node), _UNKNOWN),
                        edgecolors=SURFACE, linewidths=0.8)
    axes[1].set_title(f"community {biggest}, the largest ({sizes.iloc[0]} molecules)",
                      fontsize=9.5, color=SECONDARY, loc="left")

    figure.suptitle(title, x=0.01, ha="left", fontsize=12.5, fontweight="bold", color=INK)
    figure.tight_layout()
    return figure


def plot_values_heatmap(values: pd.DataFrame, graph=None, max_rows: int = 40,
                        label: str = "log2 fold change",
                        title: str = "The same molecules across conditions"):
    """Molecules x conditions, ordered by layer.

    For the comparisons a single contrast cannot show: six MoTrPAC tissues
    side by side, or several clinical groupings of one cohort. Rows are
    ordered by layer when a ``graph`` is given, so the block structure is
    readable rather than alphabetical.
    """
    frame = values.dropna(how="all")
    if frame.empty:
        raise ValueError("nothing to draw")
    frame = frame.loc[frame.abs().max(axis=1).sort_values(ascending=False).index[:max_rows]]

    layers = None
    if graph is not None:
        layers = pd.Series({n: _layer(graph, n) for n in frame.index if n in graph})
        if not layers.empty:
            order = sorted(frame.index, key=lambda n: (
                order_layers(layers.unique()).index(layers.get(n, "Unknown"))
                if layers.get(n, "Unknown") in list(layers.unique()) else 99, str(n)))
            frame = frame.loc[order]

    limit = float(np.nanmax(np.abs(frame.to_numpy()))) or 1.0
    height = max(3.0, 0.24 * len(frame) + 1.6)
    figure, axes = plt.subplots(figsize=(1.35 * frame.shape[1] + 5.0, height))
    image = axes.imshow(frame.to_numpy(dtype=float), cmap="RdBu_r", vmin=-limit, vmax=limit,
                        aspect="auto")

    names = []
    for node in frame.index:
        name = str(graph.nodes[node].get("symbol") or graph.nodes[node].get("name") or node) \
            if graph is not None and node in graph else str(node)
        names.append(name.split(";")[0][:26])
    axes.set_yticks(range(len(frame)))
    axes.set_yticklabels(names, fontsize=8)
    axes.set_xticks(range(frame.shape[1]))
    axes.set_xticklabels(frame.columns, fontsize=9, rotation=30, ha="right")
    axes.tick_params(length=0)
    for side in axes.spines.values():
        side.set_visible(False)
    if layers is not None and not layers.empty:
        for i, node in enumerate(frame.index):
            axes.text(-0.9, i, "", va="center")
            axes.scatter(-0.75, i, s=28, clip_on=False,
                         color=LAYER_COLORS.get(layers.get(node, "Unknown"), _UNKNOWN))
    bar = figure.colorbar(image, ax=axes, fraction=0.025, pad=0.02)
    bar.set_label(label, fontsize=9, color=SECONDARY)
    bar.outline.set_visible(False)
    axes.set_title(title, loc="left", fontsize=12.5, fontweight="bold", color=INK, pad=12)
    figure.tight_layout()
    return figure
