"""An interactive counterpart to the static trans-omic figure.

:func:`~transnet.visualization.plot_transomic_network` is for print: fixed
layout, labels in ink, a legend. This is for looking: hover a node to see what
it is and what it did, zoom into a dense region, and hand the result to
someone who does not have Python -- ``figure.write_html("network.html")``
writes a single self-contained file that opens in any browser.

Requires plotly (``pip install transnet[viz]``).
"""

from typing import Dict, Optional

import logging

import networkx as nx
import numpy as np

from transnet.biology.schema import order_layers
from transnet.visualization.palette import (
    DOWN, FLOW, LAYER_COLORS, SURFACE, UNCHANGED, UNMEASURED, UP,
)

logger = logging.getLogger(__name__)

__all__ = ["plot_transomic_network_interactive"]

_LAYOUTS = {
    "spring": lambda g: nx.spring_layout(g, seed=42, k=1.5 / max(1, np.sqrt(g.number_of_nodes()))),
    "kamada_kawai": nx.kamada_kawai_layout,
    "circular": nx.circular_layout,
    "spectral": nx.spectral_layout,
    "layered": None,      # handled below: layers stacked, as in the 2.5D figure
}

#: Regulation is drawn in the sign colours; mass flow stays grey, so the
#: coloured edges are the handful that carry a regulatory claim.
_SIGNED_EDGE_TYPES = frozenset({
    "allosteric_activation", "allosteric_inhibition",
    "transcriptional_regulation", "phosphorylation", "kinase_tf",
})


def _node_fill(data) -> str:
    regulated = int(data.get("regulated", 0) or 0)
    if regulated > 0:
        return UP
    if regulated < 0:
        return DOWN
    measured = data.get("measured") or data.get("responsive") or any(
        data.get(key) is not None for key in ("value", "log2fc", "qvalue"))
    return UNCHANGED if measured else UNMEASURED


def _hover(graph, node) -> str:
    data = graph.nodes[node]
    name = data.get("symbol") or data.get("name") or node
    lines = [f"<b>{str(name)[:60]}</b>", f"{node}",
             f"layer: {data.get('layer', 'Unknown')}",
             f"degree: {graph.degree(node)}"]
    if data.get("log2fc") is not None:
        lines.append(f"log2FC: {float(data['log2fc']):+.2f}")
    if data.get("qvalue") is not None:
        lines.append(f"q: {float(data['qvalue']):.2g}")
    return "<br>".join(lines)


def _layered_positions(graph):
    """Layers stacked top to bottom, nodes spread within each."""
    layers = order_layers({str(d.get("layer") or "Unknown") for _, d in graph.nodes(data=True)})
    position, spacing = {}, 1.0 / max(1, len(layers) - 1 if len(layers) > 1 else 1)
    for row, layer in enumerate(layers):
        members = [n for n, d in graph.nodes(data=True)
                   if str(d.get("layer") or "Unknown") == layer]
        for i, node in enumerate(sorted(members, key=str)):
            x = (i + 0.5) / max(1, len(members))
            position[node] = np.array([x, 1.0 - row * spacing])
    return position


def plot_transomic_network_interactive(graph, layout: str = "spring",
                                       max_nodes: int = 1500,
                                       title: str = "Trans-omic network"):
    """Draw the network as an interactive plotly figure.

    Parameters
    ----------
    graph : networkx.MultiDiGraph
        A trans-omic network, ideally with data mapped onto it.
    layout : str
        ``spring``, ``kamada_kawai``, ``circular``, ``spectral`` or
        ``layered`` -- the last stacks the layers the way the static 2.5D
        figure does, which is the one to use when the point is cross-layer
        structure rather than clustering.
    max_nodes : int
        Refuse to draw more than this; a browser handles a few thousand nodes,
        not a whole interactome. Pass a subnetwork instead -- see
        :func:`~transnet.responsive_subnetwork`.
    title : str

    Returns
    -------
    plotly.graph_objects.Figure
        Call ``.write_html(path)`` for a standalone file, or ``.show()``.
    """
    try:
        import plotly.graph_objects as go
    except ImportError as error:              # pragma: no cover - optional dependency
        raise ImportError(
            "The interactive view needs plotly: pip install 'transnet[viz]'"
        ) from error

    if graph.number_of_nodes() == 0:
        raise ValueError("nothing to draw: the network is empty")
    if graph.number_of_nodes() > max_nodes:
        raise ValueError(
            f"{graph.number_of_nodes():,} nodes is more than a browser draws well "
            f"(max_nodes={max_nodes:,}). Pass a subnetwork -- responsive_subnetwork(graph) "
            f"or transomic_backbone(graph) -- or raise max_nodes deliberately."
        )

    if layout == "layered":
        position = _layered_positions(graph)
    else:
        if layout not in _LAYOUTS:
            raise ValueError(f"unknown layout '{layout}'; choose from {sorted(_LAYOUTS)}")
        undirected = graph.to_undirected() if graph.is_directed() else graph
        position = _LAYOUTS[layout](undirected)

    # Edges: one trace per colour, so regulation and mass flow are separable
    # in the legend and can be switched off independently.
    groups: Dict[str, list] = {"activating": [], "inhibiting": [], "flow": []}
    for u, v, data in graph.edges(data=True):
        if u not in position or v not in position:
            continue
        edge_type = str(data.get("edge_type") or "")
        sign = int(data.get("sign", 0) or 0)
        if edge_type in _SIGNED_EDGE_TYPES and sign > 0:
            key = "activating"
        elif edge_type in _SIGNED_EDGE_TYPES and sign < 0:
            key = "inhibiting"
        else:
            key = "flow"
        x0, y0 = position[u]
        x1, y1 = position[v]
        groups[key] += [(x0, y0), (x1, y1), (None, None)]

    traces = []
    for key, colour, width, name in (
        ("flow", FLOW, 0.5, "mass flow / catalysis"),
        ("activating", UP, 1.0, "activating"),
        ("inhibiting", DOWN, 1.0, "inhibiting"),
    ):
        points = groups[key]
        if not points:
            continue
        traces.append(go.Scatter(
            x=[p[0] for p in points], y=[p[1] for p in points],
            mode="lines", line=dict(width=width, color=colour),
            hoverinfo="none", name=name, legendgroup="edges",
        ))

    # Nodes: one trace per layer, coloured by measured change, ringed in the
    # layer colour -- the same encoding as the static figure.
    by_layer: Dict[str, list] = {}
    for node, data in graph.nodes(data=True):
        if node in position:
            by_layer.setdefault(str(data.get("layer") or "Unknown"), []).append(node)

    for layer in order_layers(by_layer):
        members = by_layer[layer]
        traces.append(go.Scatter(
            x=[position[n][0] for n in members],
            y=[position[n][1] for n in members],
            mode="markers", name=layer,
            marker=dict(
                size=[6 + min(26, graph.degree(n) * 1.5) for n in members],
                color=[_node_fill(graph.nodes[n]) for n in members],
                line=dict(width=1.6, color=LAYER_COLORS.get(layer, "#898781")),
            ),
            text=[_hover(graph, n) for n in members],
            hoverinfo="text",
        ))

    figure = go.Figure(data=traces, layout=go.Layout(
        title=dict(text=title, x=0.01, xanchor="left"),
        showlegend=True, hovermode="closest",
        margin=dict(b=20, l=10, r=10, t=48),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        paper_bgcolor=SURFACE, plot_bgcolor=SURFACE,
    ))
    logger.info(f"Interactive view: {graph.number_of_nodes():,} nodes, "
                f"{graph.number_of_edges():,} edges, layout '{layout}'")
    return figure
