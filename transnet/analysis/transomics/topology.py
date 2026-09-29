"""Layer-aware topology: cross-layer connectivity and trans-omic hubs.

Two measurements that only mean something on a multi-layer network:

* **cross-layer connectivity** -- how much of the network actually crosses
  between omic layers, and how well each layer is covered by measurement.  A
  network that is 95% within-layer edges is a set of stacked single-omics
  networks, not a trans-omic one, and this is how you find that out.
* **trans-omic hubs** -- molecules that are hubs *across* layers.  Morita et al.
  found ATP and AMP dominating the starvation-responsive liver network at
  degrees above 200, coordinating reactions throughout metabolism; a
  within-layer degree ranking would never surface them as such.

References
----------
Morita K, et al. Structural robustness and temporal vulnerability of the
starvation-responsive metabolic network in healthy and obese mouse liver.
*Science Signaling* 18, 2025.

De Domenico M, Sole-Ribalta A, Omodei E, Gomez S, Arenas A. Ranking in
interconnected multilayer networks reveals versatile nodes. *Nature
Communications* 6:6868, 2015.

Sugimoto H, Morita K, Li D, et al. iTraNet: a web-based platform for integrated
trans-omics network visualization and analysis. *Bioinformatics Advances*
4(1):vbae141, 2024.
"""

from typing import Dict, List, Optional, Sequence

import logging

import numpy as np
import pandas as pd

from transnet.analysis.transomics.mapping import is_responsive
from transnet.biology.schema import available_edge_types, available_layers

logger = logging.getLogger(__name__)

__all__ = [
    "cross_layer_connectivity",
    "layer_coverage",
    "transomic_hubs",
]


def cross_layer_connectivity(graph) -> Dict[str, object]:
    """Measure how much of the network crosses between omic layers.

    Parameters
    ----------
    graph : networkx.Graph

    Returns
    -------
    dict
        ``matrix`` : pandas.DataFrame
            Layer x layer edge counts, source layer on rows.
        ``edge_types`` : pandas.DataFrame
            Count of each relationship, with the layers it connects.
        ``coverage`` : pandas.DataFrame
            Per layer: node count, how many carry measured data, and how many
            are differentially regulated.
        ``cross_layer_fraction`` : float
            Share of edges that connect two different layers.  The headline
            number: this is what makes the network trans-omic.
    """
    layers = available_layers(graph)
    matrix = pd.DataFrame(0, index=layers, columns=layers, dtype=int)

    cross = total = 0
    for u, v, data in graph.edges(data=True):
        source_layer = graph.nodes[u].get("layer")
        target_layer = graph.nodes[v].get("layer")
        total += 1
        if source_layer in matrix.index and target_layer in matrix.columns:
            matrix.loc[source_layer, target_layer] += 1
        if source_layer != target_layer:
            cross += 1

    type_rows = []
    for etype, count in available_edge_types(graph).items():
        pairs = {
            (graph.nodes[u].get("layer"), graph.nodes[v].get("layer"))
            for u, v, d in graph.edges(data=True)
            if (d.get("edge_type") or "unknown") == etype
        }
        type_rows.append({
            "edge_type": etype,
            "n_edges": count,
            "layer_pairs": "; ".join(
                f"{a}->{b}" for a, b in sorted(pairs, key=lambda p: (str(p[0]), str(p[1])))
            ),
            "cross_layer": all(a != b for a, b in pairs) if pairs else False,
        })

    result = {
        "matrix": matrix,
        "edge_types": pd.DataFrame(type_rows),
        "coverage": layer_coverage(graph),
        "cross_layer_fraction": (cross / total) if total else 0.0,
    }

    logger.info(
        f"{cross}/{total} edges cross layers ({result['cross_layer_fraction']:.1%})"
    )
    return result


def layer_coverage(graph) -> pd.DataFrame:
    """Per-layer node counts, measurement coverage and regulated counts.

    Returns
    -------
    pandas.DataFrame
        Columns ``layer``, ``n_nodes``, ``n_measured``, ``measured_fraction``,
        ``n_regulated``, ``n_up``, ``n_down``, ``mean_degree``.
    """
    rows = []
    for layer in available_layers(graph):
        nodes = [n for n, d in graph.nodes(data=True) if d.get("layer") == layer]
        # q-value-only tables (an omnibus test) are measurements too
        measured = [
            n for n in nodes
            if graph.nodes[n].get("measured")
            or any(graph.nodes[n].get(key) is not None
                   for key in ("log2fc", "value", "qvalue"))
        ]
        states = [graph.nodes[n].get("regulated", 0) or 0 for n in nodes]
        n_responsive = sum(1 for n in nodes if is_responsive(graph.nodes[n]))
        degrees = [graph.degree(n) for n in nodes]
        rows.append({
            "layer": layer,
            "n_nodes": len(nodes),
            "n_measured": len(measured),
            "measured_fraction": len(measured) / len(nodes) if nodes else 0.0,
            "n_responsive": n_responsive,
            "n_regulated": sum(1 for s in states if s),
            "n_up": sum(1 for s in states if s > 0),
            "n_down": sum(1 for s in states if s < 0),
            "mean_degree": float(np.mean(degrees)) if degrees else 0.0,
        })
    return pd.DataFrame(rows, columns=[
        "layer", "n_nodes", "n_measured", "measured_fraction",
        "n_responsive", "n_regulated", "n_up", "n_down", "mean_degree",
    ])


def transomic_hubs(
    graph,
    top_percent: float = 2.0,
    responsive_only: bool = True,
    layers: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """Rank molecules by how much they connect *across* layers.

    Follows Morita et al., who define hubs as the top 2% by degree among the
    nodes that responded to the perturbation, then reads off which of those
    coordinate multiple metabolic pathways.  Alongside plain degree this
    reports two multi-layer quantities:

    * ``n_layers_touched`` -- how many distinct layers a node's neighbours span;
    * ``versatility`` -- the node's share of total degree, weighted by the
      number of layers it reaches, in the spirit of De Domenico's multilayer
      versatility.  A node with 50 neighbours in one layer scores below a node
      with 50 neighbours spread across four.

    Parameters
    ----------
    graph : networkx.Graph
    top_percent : float
        Degree percentile defining a hub, as a percentage.
    responsive_only : bool
        Rank only among differentially regulated nodes (default), matching the
        published definition.  Falls back to all nodes, with a log message,
        when no data has been mapped.
    layers : sequence of str, optional
        Restrict the ranking to these layers.

    Returns
    -------
    pandas.DataFrame
        Columns ``node``, ``layer``, ``name``, ``degree``, ``cross_layer_degree``,
        ``n_layers_touched``, ``layers_touched``, ``versatility``, ``regulated``,
        ``is_hub``, sorted with the strongest hubs first.
    """
    wanted = set(layers) if layers is not None else set(available_layers(graph))
    candidates = [
        n for n, d in graph.nodes(data=True) if d.get("layer") in wanted
    ]

    if responsive_only:
        responsive = [n for n in candidates if is_responsive(graph.nodes[n])]
        if responsive:
            candidates = responsive
        else:
            logger.info(
                "No regulated nodes found; ranking hubs across all nodes instead. "
                "Map omics data first for the published definition."
            )

    if not candidates:
        return pd.DataFrame(columns=[
            "node", "layer", "name", "degree", "cross_layer_degree",
            "n_layers_touched", "layers_touched", "versatility",
            "regulated", "is_hub",
        ])

    rows = []
    for node in candidates:
        own_layer = graph.nodes[node].get("layer")
        neighbours = set(graph.predecessors(node)) | set(graph.successors(node)) \
            if graph.is_directed() else set(graph.neighbors(node))

        neighbour_layers = {
            graph.nodes[n].get("layer") for n in neighbours
        }
        neighbour_layers.discard(None)
        cross_layer_degree = sum(
            1 for n in neighbours if graph.nodes[n].get("layer") != own_layer
        )
        rows.append({
            "node": node,
            "layer": own_layer,
            "name": graph.nodes[node].get("name"),
            "degree": len(neighbours),
            "cross_layer_degree": cross_layer_degree,
            "n_layers_touched": len(neighbour_layers),
            "layers_touched": ";".join(sorted(str(x) for x in neighbour_layers)),
            "regulated": graph.nodes[node].get("regulated", 0) or 0,
        })

    table = pd.DataFrame(rows)
    total_degree = table["degree"].sum()
    n_layers = max(len(wanted), 1)
    table["versatility"] = (
        (table["degree"] / total_degree if total_degree else 0.0)
        * (table["n_layers_touched"] / n_layers)
    )

    if len(table):
        cutoff = np.percentile(table["degree"], 100.0 - top_percent)
        table["is_hub"] = table["degree"] >= cutoff
    else:
        table["is_hub"] = False

    table = table.sort_values(
        ["is_hub", "versatility", "degree"], ascending=False
    ).reset_index(drop=True)

    n_hubs = int(table["is_hub"].sum())
    logger.info(
        f"{n_hubs} hub molecules (top {top_percent}% by degree among "
        f"{len(table)} candidates); highest degree {int(table['degree'].max())}"
    )
    return table[[
        "node", "layer", "name", "degree", "cross_layer_degree",
        "n_layers_touched", "layers_touched", "versatility",
        "regulated", "is_hub",
    ]]
