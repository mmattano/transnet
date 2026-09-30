"""Figures of the trans-omic network.

* :func:`plot_transomic_network`: the layers as stacked planes, with edges
  between them; line style shows the edge type and colour its sign.
* :func:`plot_regulation_axes`: per pathway, the reactions each regulation
  axis activates or inhibits, and the fraction that are controversial.
* :func:`plot_layer_connectivity`: a heatmap of edge counts between layers.
* :func:`transomic_backbone`: selects a readable part of a large network to
  draw.
"""

import logging
from typing import Dict, List, Optional, Sequence

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

from transnet.analysis.transomics.mapping import is_measured as _measured
from transnet.biology.schema import (
    DISPLAY_ORDER, available_layers, metabolic_pool,
)

logger = logging.getLogger(__name__)

__all__ = [
    "plot_transomic_network",
    "transomic_backbone",
    "plot_regulation_axes",
    "plot_layer_connectivity",
    "LAYER_COLORS",
    "EDGE_STYLES",
]

from transnet.visualization.palette import (  # noqa: E402
    DOWN, FLOW as FLOW_COLOR, INK, LAYER_COLORS, SECONDARY, SIGN_COLORS,
    UNCHANGED as UNCHANGED_COLOR, UNDIRECTED as UNDIRECTED_COLOR,
    UNMEASURED as UNMEASURED_COLOR, UP, style_axes,
)

#: Line style per relationship, so an edge's kind is readable without a legend
#: lookup: solid for mass flow, dashed for regulation, dotted for association.
EDGE_STYLES = {
    "substrate": "-",
    "product": "-",
    "catalysis": "--",
    "gene_catalysis": "--",
    "translation": "-",
    "transcriptional_regulation": "--",
    "phosphorylation": "--",
    "kinase_tf": "--",
    "allosteric_activation": ":",
    "allosteric_inhibition": ":",
    "protein_interaction": ":",
    "enzymatic": "--",
    "unknown": "-",
}

#: Edge types whose sign is a *regulatory* direction rather than a structural
#: default.  A substrate edge carries sign +1 because consuming more substrate
#: drives the reaction forward, not because anything was measured -- colouring
#: those by sign would drown the figure in one colour and hide the handful of
#: edges where the sign is the finding.
REGULATORY_EDGE_TYPES = frozenset({
    "allosteric_activation", "allosteric_inhibition",
    "transcriptional_regulation", "phosphorylation", "kinase_tf",
})

def _label(name, node, limit: int = 18, symbol=None) -> str:
    """A short label from a database description.

    KEGG stores compound synonyms as ``"NADPH; TPNH; Reduced ..."`` and UniProt
    stores full protein descriptions; printing either verbatim turns the figure
    into a wall of text. A gene symbol wins when there is one; otherwise take
    the first clause.
    """
    if symbol and str(symbol) not in ("", "nan", "None"):
        return str(symbol)[:limit]
    text = str(name or node).split(";")[0].strip()
    if " (" in text:
        text = text.split(" (")[0].strip()
    return text[:limit]


#: Edges drawn between planes in the stacked view. Protein-protein
#: interactions are left out: they run *within* the Proteome plane, and at
#: interactome scale they are thousands of edges that hide the cross-layer
#: flow the figure exists to show.
_BACKBONE_EDGE_TYPES = frozenset({
    "translation", "catalysis", "substrate", "product",
    "allosteric_activation", "allosteric_inhibition",
    "transcriptional_regulation", "phosphorylation", "kinase_tf",
})


def _incident(graph, reaction, metabolite, edge_types) -> List[dict]:
    """Edge records joining a reaction and one of its metabolites, either way."""
    records = []
    for data in (graph.get_edge_data(metabolite, reaction) or {}).values():
        records.append(data)
    for data in (graph.get_edge_data(reaction, metabolite) or {}).values():
        records.append(data)
    return [d for d in records if d.get("edge_type") in edge_types]


def _regulated(graph, node) -> int:
    return int(graph.nodes[node].get("regulated", 0) or 0)


def transomic_backbone(
    graph,
    max_reactions: int = 20,
    max_enzymes_per_reaction: int = 3,
    max_metabolites_per_reaction: int = 4,
    include_currency: bool = False,
    metabolic_effectors_only: bool = True,
    max_reactions_per_enzyme: int = 3,
    max_signalling_per_enzyme: int = 2,
):
    """Select a readable part of a large network, built around reactions.

    Keeping the highest-degree nodes of an organism-wide network keeps
    interactome hubs and drops every reaction. This instead ranks reactions by
    how much regulation meets there (a changed enzyme and a changed metabolite
    count most) and keeps each reaction with its enzymes, their genes, and its
    metabolites.

    Parameters
    ----------
    graph : networkx.MultiDiGraph
        A mapped trans-omic network.
    max_reactions : int
        Reactions to keep.
    max_enzymes_per_reaction, max_metabolites_per_reaction : int
        Limits per reaction.
    include_currency : bool
        Keep currency metabolites such as ATP and water (off by default).
    metabolic_effectors_only : bool
        Draw only allosteric regulators that occur in the organism's
        metabolism, not laboratory compounds recorded by BRENDA. Affects the
        figure only.
    max_reactions_per_enzyme : int
        Reactions to keep per enzyme, so one enzyme family with many similar
        reactions cannot fill the figure.
    max_signalling_per_enzyme : int
        Signaling nodes to draw above each enzyme; 0 leaves the layer out.

    Returns
    -------
    networkx.MultiDiGraph
        The selected nodes, with the edges between layers.
    """
    from transnet.biology.schema import CURRENCY_METABOLITES

    currency = set() if include_currency else set(CURRENCY_METABOLITES)
    pool = metabolic_pool(graph) if metabolic_effectors_only else None
    layer = {n: d.get("layer") for n, d in graph.nodes(data=True)}

    def drawable(metabolite, edge_type) -> bool:
        if metabolite in currency:
            return False
        if pool is None or not str(edge_type).startswith("allosteric"):
            return True
        return metabolite in pool

    def neighbours(node, edge_types, incoming=True, outgoing=True):
        found = []
        if incoming:
            found += [u for u, _, d in graph.in_edges(node, data=True)
                      if d.get("edge_type") in edge_types]
        if outgoing:
            found += [v for _, v, d in graph.out_edges(node, data=True)
                      if d.get("edge_type") in edge_types]
        return list(dict.fromkeys(found))

    def rank(nodes):
        # regulated first, then measured, then by degree for stability
        return sorted(nodes, key=lambda n: (
            -abs(_regulated(graph, n)),
            -int(_measured(graph.nodes[n])),
            -graph.degree(n), str(n)))

    metabolite_types = {"substrate", "product",
                        "allosteric_activation", "allosteric_inhibition"}
    scored = []
    for reaction in (n for n, l in layer.items() if l == "Reactions"):
        enzymes = neighbours(reaction, {"catalysis"}, outgoing=False)
        metabolites = [
            m for m in neighbours(reaction, metabolite_types)
            if any(drawable(m, d.get("edge_type"))
                   for d in _incident(graph, reaction, m, metabolite_types))
        ]
        reg_enzymes = sum(1 for e in enzymes if _regulated(graph, e))
        reg_metabolites = sum(1 for m in metabolites if _regulated(graph, m))
        allosteric = sum(
            1 for u, _, d in graph.in_edges(reaction, data=True)
            if str(d.get("edge_type", "")).startswith("allosteric")
            and _regulated(graph, u) and drawable(u, d.get("edge_type"))
        )
        score = (4 * (reg_enzymes > 0 and reg_metabolites > 0)
                 + reg_enzymes + reg_metabolites + allosteric)
        if score:
            scored.append((score, reaction, enzymes, metabolites))

    scored.sort(key=lambda row: (-row[0], str(row[1])))
    seen_enzymes: Dict[frozenset, int] = {}
    selected = []
    for row in scored:
        if len(selected) >= max_reactions:
            break
        catalysts = frozenset(row[2])
        if catalysts and max_reactions_per_enzyme:
            if seen_enzymes.get(catalysts, 0) >= max_reactions_per_enzyme:
                continue
            seen_enzymes[catalysts] = seen_enzymes.get(catalysts, 0) + 1
        selected.append(row)

    keep = []
    for _, reaction, enzymes, metabolites in selected:
        keep.append(reaction)
        chosen_enzymes = rank(enzymes)[:max_enzymes_per_reaction]
        keep += chosen_enzymes
        for enzyme in chosen_enzymes:
            keep += [g for g in neighbours(enzyme, {"translation"}, outgoing=False)
                     if layer.get(g) == "Transcriptome"][:1]
            # The Signaling layer sits above the Proteome, so walking up from a
            # selected enzyme is the only way a phosphosite or kinase reaches the
            # figure. Without this the layer is silently absent from every study
            # figure even when the study measured it.
            signalling = [s for s in neighbours(enzyme, {"phosphorylation"},
                                                outgoing=False)
                          if layer.get(s) == "Signaling"]
            keep += rank(signalling)[:max_signalling_per_enzyme]
        keep += rank(metabolites)[:max_metabolites_per_reaction]

    keep = list(dict.fromkeys(keep))
    backbone = graph.subgraph(keep).copy()
    drop = [(u, v, k) for u, v, k, d in backbone.edges(keys=True, data=True)
            if d.get("edge_type") not in _BACKBONE_EDGE_TYPES]
    backbone.remove_edges_from(drop)
    logger.info(
        f"Backbone: {len(scored)} reactions carry regulatory evidence; "
        f"kept {len(selected)} with "
        f"{backbone.number_of_nodes()} nodes, {backbone.number_of_edges()} edges"
    )
    return backbone


def _aligned_positions(graph, layers, spread: float = 5.0, sweeps: int = 4):
    """Place each node over or under the neighbours it connects to.

    A layered-graph (Sugiyama-style) barycentre heuristic: start from the
    Reactions plane, place each adjacent plane at the mean position of its
    already-placed neighbours, then sweep down and up a few times so every
    plane is refined against both of its neighbours. That is what removes most
    edge crossings; a single outward pass left metabolites ordered against
    reactions only. Nodes are spaced evenly within a plane in their barycentre
    order. Returns ``{node: x}`` in ``[-spread, spread]``.
    """
    by_layer = {l: [n for n, d in graph.nodes(data=True) if d.get("layer") == l]
                for l in layers}
    present = [l for l in layers if by_layer.get(l)]
    if not present:
        return {}
    undirected = graph.to_undirected(as_view=True)
    x: Dict[str, float] = {}

    def place(layer, order):
        nodes = list(order)
        if len(nodes) == 1:
            x[nodes[0]] = 0.0
            return
        for i, node in enumerate(nodes):
            x[node] = -spread + 2 * spread * i / (len(nodes) - 1)

    def barycentre_order(layer, against):
        wanted = {}
        for position, node in enumerate(by_layer[layer]):
            placed = [x[m] for m in undirected.neighbors(node)
                      if m in x and graph.nodes[m].get("layer") in against]
            # unanchored nodes keep their current slot
            wanted[node] = (float(np.mean(placed)) if placed
                            else x.get(node, float(position)))
        return sorted(by_layer[layer], key=lambda n: (wanted[n], str(n)))

    anchor = "Reactions" if "Reactions" in present else present[0]
    place(anchor, sorted(by_layer[anchor], key=str))
    index = present.index(anchor)
    for layer in [present[i] for i in range(index - 1, -1, -1)] + present[index + 1:]:
        i = present.index(layer)
        neighbours = {present[j] for j in (i - 1, i + 1) if 0 <= j < len(present)}
        place(layer, barycentre_order(layer, neighbours))

    for _ in range(sweeps):
        for direction in (present, list(reversed(present))):
            for layer in direction:
                i = present.index(layer)
                neighbours = {present[j] for j in (i - 1, i + 1)
                              if 0 <= j < len(present)}
                place(layer, barycentre_order(layer, neighbours))
    return x


def _ordered_layers(graph, layer_order: Optional[Sequence[str]]) -> List[str]:
    present = available_layers(graph)
    if layer_order is not None:
        return [layer for layer in layer_order if layer in present]
    ordered = [layer for layer in DISPLAY_ORDER if layer in present]
    return ordered + [layer for layer in present if layer not in ordered]


def plot_transomic_network(
    graph,
    layer_order: Optional[Sequence[str]] = None,
    node_color_by: str = "regulated",
    label_top_n: int = 8,
    layer_spacing: float = 3.0,
    skew: float = 0.35,
    figsize=(14, 10),
    title: str = "Trans-omic network",
    seed: int = 0,
    layout: str = "auto",
):
    """Draw the network as stacked layer planes.

    Each layer is laid out in its own plane, and the planes are stacked in the
    order regulation flows, so edges between layers run vertically. Reduce large
    networks first with :func:`transomic_backbone`.

    Parameters
    ----------
    graph : networkx.Graph
        A trans-omic network. Absent layers are skipped.
    layer_order : sequence of str, optional
        Top-to-bottom order. Default: ``DISPLAY_ORDER``.
    node_color_by : {"regulated", "layer"}
        Colour by measured direction (red up, blue down, grey unchanged) or by
        layer.
    label_top_n : int
        Label this many of the best-connected nodes per layer.
    layer_spacing : float
        Vertical distance between planes.
    skew : float
        Tilt of the planes; 0 draws flat rows.
    figsize : tuple
    title : str
    seed : int
        Layout seed.
    layout : {"auto", "aligned", "spring"}
        ``"aligned"`` places nodes above or below the reactions they connect
        to; ``"spring"`` lays out each plane independently; ``"auto"`` aligns
        when there is a Reactions layer.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import networkx as nx

    layers = _ordered_layers(graph, layer_order)
    if not layers:
        raise ValueError("Graph has no nodes carrying a 'layer' attribute")

    positions: Dict[str, tuple] = {}
    if layout == "auto":
        layout = "aligned" if "Reactions" in layers else "spring"
    busiest = max(
        (sum(1 for _, d in graph.nodes(data=True) if d.get("layer") == l)
         for l in layers), default=1)
    # about 0.9 units per node keeps rotated labels clear of each other
    spread = max(5.0, 0.45 * busiest)
    if layout == "aligned":
        aligned_x = _aligned_positions(graph, layers, spread=spread)
        for index, layer in enumerate(layers):
            nodes = sorted(
                (n for n, d in graph.nodes(data=True) if d.get("layer") == layer),
                key=lambda n: aligned_x.get(n, 0.0),
            )
            depth = -index * layer_spacing
            for position, node in enumerate(nodes):
                # a small alternating stagger keeps neighbouring labels apart
                y = 0.35 if position % 2 else -0.35
                positions[node] = (aligned_x.get(node, 0.0) + y * skew,
                                   depth + y * skew)
    for index, layer in enumerate(layers if layout != "aligned" else []):
        nodes = [n for n, d in graph.nodes(data=True) if d.get("layer") == layer]
        if not nodes:
            continue

        subgraph = graph.subgraph(nodes)
        try:
            within = nx.spring_layout(
                nx.Graph(subgraph), seed=seed, k=1.2, iterations=60
            )
        except Exception:                                  # pragma: no cover
            within = {node: (i, 0.0) for i, node in enumerate(nodes)}

        xs = np.array([within[n][0] for n in nodes], dtype=float)
        ys = np.array([within[n][1] for n in nodes], dtype=float)
        x_range, y_range = np.ptp(xs), np.ptp(ys)
        if x_range > 0:
            xs = (xs - xs.min()) / x_range * 10 - 5
        if y_range > 0:
            ys = (ys - ys.min()) / y_range * 2 - 1

        depth = -index * layer_spacing
        for node, x, y in zip(nodes, xs, ys):
            positions[node] = (x + y * skew * 2, depth + y * skew)

    figure, axes = plt.subplots(figsize=figsize)

    # Layer planes, drawn first so edges and nodes sit on top.
    for index, layer in enumerate(layers):
        depth = -index * layer_spacing
        axes.add_patch(mpatches.FancyBboxPatch(
            (-spread - 2.5, depth - 1.1), 2 * spread + 5.0, 2.2,
            boxstyle="round,pad=0.15",
            facecolor=LAYER_COLORS.get(layer, SECONDARY),
            alpha=0.07, edgecolor=LAYER_COLORS.get(layer, SECONDARY),
            linewidth=1.0, zorder=0,
        ))
        axes.text(
            -spread - 3.0, depth, layer, ha="right", va="center",
            fontsize=11, fontweight="bold",
            color=INK, zorder=5,
        )

    drawn_edges = set()
    for u, v, data in graph.edges(data=True):
        if u not in positions or v not in positions:
            continue
        edge_type = data.get("edge_type", "unknown")
        sign = int(data.get("sign", 0) or 0)
        same_layer = graph.nodes[u].get("layer") == graph.nodes[v].get("layer")
        regulatory = edge_type in REGULATORY_EDGE_TYPES
        drawn_edges.add(("regulatory", sign) if regulatory else ("flow", 0))
        axes.annotate(
            "", xy=positions[v], xytext=positions[u],
            arrowprops=dict(
                arrowstyle="-|>" if graph.is_directed() else "-",
                linestyle=EDGE_STYLES.get(edge_type, "-"),
                color=SIGN_COLORS[sign] if regulatory else FLOW_COLOR,
                alpha=(0.85 if regulatory else 0.35) * (0.6 if same_layer else 1.0),
                linewidth=(1.4 if regulatory else 0.8),
                shrinkA=4, shrinkB=4,
                connectionstyle="arc3,rad=0.08",
            ),
            zorder=2 if regulatory else 1,
        )

    drawn_states = set()
    for layer in layers:
        nodes = [
            n for n in positions if graph.nodes[n].get("layer") == layer
        ]
        if not nodes:
            continue

        colors, sizes = [], []
        for node in nodes:
            attributes = graph.nodes[node]
            state = attributes.get("regulated", 0) or 0
            measured = _measured(attributes)
            if node_color_by == "regulated":
                if state:
                    colors.append(SIGN_COLORS[state])
                    drawn_states.add(state)
                elif attributes.get("responsive"):
                    colors.append(UNDIRECTED_COLOR)
                    drawn_states.add("undirected")
                elif measured:
                    colors.append(UNCHANGED_COLOR)
                    drawn_states.add("unchanged")
                else:
                    colors.append(UNMEASURED_COLOR)
                    drawn_states.add("unmeasured")
            else:
                colors.append(LAYER_COLORS.get(layer, SECONDARY))
            # Size carries the measurement, not the degree: sizing by degree
            # made interactome hubs enormous and shrank metabolites to dots
            # whose colour could not be read.
            sizes.append(260 if state else (150 if measured else 90))

        axes.scatter(
            [positions[n][0] for n in nodes],
            [positions[n][1] for n in nodes],
            s=sizes, c=colors,
            marker="s" if layer == "Reactions" else "o",
            edgecolors=LAYER_COLORS.get(layer, SECONDARY),
            linewidths=1.6, zorder=3,
        )

        # Label regulated nodes first, then the best-connected. Alternate the
        # offset so neighbouring labels in a dense layer do not overlap.
        limit = len(nodes) if len(nodes) <= max(label_top_n, 24) else label_top_n
        labelled = sorted(
            nodes,
            key=lambda n: (-abs(graph.nodes[n].get("regulated", 0) or 0),
                           -graph.degree(n)),
        )[:limit]
        crowded = len(nodes) > 10
        for position, node in enumerate(sorted(labelled, key=lambda n: positions[n][0])):
            attributes = graph.nodes[node]
            text = _label(attributes.get("name"), node,
                          limit=26 if layer in ("Reactions", "Metabolome") else 18,
                          symbol=attributes.get("symbol"))
            above = position % 2 == 0
            axes.annotate(
                text, positions[node], fontsize=7,
                xytext=(0, 12 if above else -14),
                textcoords="offset points",
                # rotate in crowded planes so neighbouring labels cannot
                # overlap; anchor at the node so the rotation fans outward
                rotation=30 if crowded else 0,
                ha=("left" if above else "right") if crowded else "center",
                va="bottom" if above else "top",
                rotation_mode="anchor", zorder=4,
            )

    def dot(color, label):
        return plt.Line2D([], [], color=color, marker="o", linestyle="",
                          markersize=9, markeredgecolor=SECONDARY, label=label)

    handles = []
    if 1 in drawn_states:
        handles.append(dot(SIGN_COLORS[1], "increased"))
    if -1 in drawn_states:
        handles.append(dot(SIGN_COLORS[-1], "decreased"))
    if "undirected" in drawn_states:
        handles.append(dot(UNDIRECTED_COLOR, "changed, direction unknown"))
    if "unchanged" in drawn_states:
        handles.append(dot(UNCHANGED_COLOR, "measured, not changed"))
    if "unmeasured" in drawn_states:
        handles.append(dot(UNMEASURED_COLOR, "not measured"))
    if ("flow", 0) in drawn_edges:
        handles.append(plt.Line2D([], [], color=FLOW_COLOR, linestyle="-",
                                  label="mass flow / catalysis"))
    if ("regulatory", 1) in drawn_edges:
        handles.append(plt.Line2D([], [], color=SIGN_COLORS[1], linestyle=":",
                                  linewidth=1.6, label="activating (allosteric)"))
    if ("regulatory", -1) in drawn_edges:
        handles.append(plt.Line2D([], [], color=SIGN_COLORS[-1], linestyle=":",
                                  linewidth=1.6, label="inhibiting (allosteric)"))
    if ("regulatory", 0) in drawn_edges:
        handles.append(plt.Line2D([], [], color=SIGN_COLORS[0], linestyle="--",
                                  label="regulation, sign unknown"))
    axes.legend(
        handles=handles, loc="lower right", fontsize=8, framealpha=0.95,
        title="node fill = measured change, ring = layer\n"
              "coloured edges = regulation, grey = mass flow",
        title_fontsize=7.5,
    )

    axes.set_title(title, fontsize=13, fontweight="bold", color=INK)
    axes.set_xlim(-spread - 5.5, spread + 3.0)
    axes.set_ylim(-len(layers) * layer_spacing - 1.0, 1.8)
    axes.axis("off")
    figure.tight_layout()
    return figure


def plot_regulation_axes(
    summary: pd.DataFrame,
    figsize=(11, 6),
    title: str = "Regulation of metabolic reactions, by axis",
):
    """Per-pathway bars of activation and inhibition, split by regulation axis.

    Left panel: the gene-expression axis (how much enzyme).  Right panel: the
    metabolite axis (how hard the enzyme works).  The marked fraction is where
    the two disagree -- reactions no single-axis analysis would flag.

    Parameters
    ----------
    summary : pandas.DataFrame
        Output of
        :func:`~transnet.analysis.transomics.regulation_axis_summary`.
    figsize : tuple
    title : str

    Returns
    -------
    matplotlib.figure.Figure
    """
    if summary.empty:
        raise ValueError(
            "Nothing to plot: the regulation summary is empty. Map omics data "
            "with map_omics_to_network() before building the table."
        )

    data = summary.sort_values("n_reactions")
    y = np.arange(len(data))
    figure, (left, right, extra) = plt.subplots(
        1, 3, figsize=figsize, sharey=True,
        gridspec_kw={"width_ratios": [3, 3, 1.4]},
    )

    for axis, (up, down, label) in [
        (left, ("gene_activated", "gene_inhibited", "Gene-expression axis")),
        (right, ("metabolite_activated", "metabolite_inhibited", "Metabolite axis")),
    ]:
        axis.barh(y, data[up], color=SIGN_COLORS[1], label="activated")
        axis.barh(y, -data[down], color=SIGN_COLORS[-1], label="inhibited")
        axis.axvline(0, color=SECONDARY, linewidth=0.8)
        axis.set_title(label, fontsize=11)
        axis.set_xlabel("reactions")
        axis.grid(axis="x", alpha=0.25)

    left.set_yticks(y)
    left.set_yticklabels(data["pathway"])
    left.legend(fontsize=8, loc="lower left")

    extra.barh(y, data["fraction_controversial"], color=LAYER_COLORS["Reactions"])
    extra.set_title("axes\ndisagree", fontsize=10)
    extra.set_xlabel("fraction")
    extra.set_xlim(0, 1)
    extra.grid(axis="x", alpha=0.25)

    figure.suptitle(title, fontsize=13, fontweight="bold")
    figure.tight_layout()
    return figure


def plot_layer_connectivity(
    connectivity: Dict[str, object],
    figsize=(8, 6.5),
    title: str = "Cross-layer connectivity",
):
    """Heatmap of the layer x layer edge-count matrix.

    Parameters
    ----------
    connectivity : dict
        Output of
        :func:`~transnet.analysis.transomics.cross_layer_connectivity`.
    figsize : tuple
    title : str

    Returns
    -------
    matplotlib.figure.Figure
    """
    matrix = connectivity["matrix"]
    if matrix.empty:
        raise ValueError("Nothing to plot: the connectivity matrix is empty")

    figure, axes = plt.subplots(figsize=figsize)
    values = matrix.to_numpy(dtype=float)

    image = axes.imshow(
        np.where(values > 0, values, np.nan),
        cmap="viridis", aspect="auto",
    )

    axes.set_xticks(range(len(matrix.columns)))
    axes.set_xticklabels(matrix.columns, rotation=35, ha="right")
    axes.set_yticks(range(len(matrix.index)))
    axes.set_yticklabels(matrix.index)
    axes.set_xlabel("target layer")
    axes.set_ylabel("source layer")

    threshold = np.nanmax(values) * 0.6 if np.nanmax(values) else 0
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            if values[i, j] > 0:
                axes.text(
                    j, i, int(values[i, j]), ha="center", va="center",
                    fontsize=9,
                    color="white" if values[i, j] < threshold else "black",
                )

    figure.colorbar(image, ax=axes, label="edges", shrink=0.8)
    fraction = connectivity.get("cross_layer_fraction", 0.0)
    axes.set_title(
        f"{title}\n{fraction:.0%} of edges cross between layers",
        fontsize=12, fontweight="bold",
    )
    figure.tight_layout()
    return figure
