"""Extract and compare responsive trans-omic networks.

The core construct of a trans-omic study is not the whole reference network but
the *responsive* part of it: the molecules that changed, in every layer, plus
the regulatory edges that connect them.  That subnetwork is the object papers
draw, count and compare between conditions.

References
----------
Kawata K, et al. Trans-omic Analysis Reveals Selective Responses to Induced and
Basal Insulin across Signaling, Transcriptional, and Metabolic Networks.
*iScience* 7:212-229, 2018.

Egami R, et al. Trans-omic analysis reveals obesity-associated dysregulation of
inter-organ metabolic cycles between the liver and skeletal muscle. *iScience*
24(3):102217, 2021.
"""

from typing import Dict, List, Optional, Sequence

import logging

import pandas as pd

from transnet.analysis.transomics.mapping import is_responsive
from transnet.biology.schema import available_edge_types, available_layers

logger = logging.getLogger(__name__)

__all__ = [
    "responsive_subnetwork",
    "compare_transomic_networks",
]


def responsive_subnetwork(
    graph,
    layers: Optional[Sequence[str]] = None,
    include_connectors: bool = True,
    connector_layers: Sequence[str] = ("Reactions",),
    keep_isolated: bool = False,
):
    """Extract the differentially regulated trans-omic network.

    Parameters
    ----------
    graph : networkx.Graph
        A network with omics data mapped on.
    layers : sequence of str, optional
        Restrict the responsive set to these layers.  Defaults to all present.
    include_connectors : bool
        Keep unmeasured nodes from ``connector_layers`` that sit *between* two
        responsive nodes.  Reactions are rarely measured directly, but dropping
        them would sever the enzyme-metabolite link that makes the network
        trans-omic, so they are kept by default.
    connector_layers : sequence of str
        Layers eligible to be kept as connectors.
    keep_isolated : bool
        Keep responsive nodes with no surviving edges.  Off by default, so the
        result is the connected regulatory picture.

    Returns
    -------
    networkx.Graph
        A subgraph of the same type as ``graph``, with node and edge attributes
        preserved.  Empty if no node is regulated.
    """
    wanted = set(layers) if layers is not None else set(available_layers(graph))

    responsive = {
        n for n, d in graph.nodes(data=True)
        if d.get("layer") in wanted and is_responsive(d)
    }

    if not responsive:
        logger.warning(
            "No differentially regulated nodes; map omics data with "
            "map_omics_to_network() before extracting a responsive subnetwork."
        )
        return graph.__class__()

    selected = set(responsive)

    if include_connectors:
        connector_set = set(connector_layers)
        for node, data in graph.nodes(data=True):
            if node in selected or data.get("layer") not in connector_set:
                continue
            neighbours = (
                set(graph.predecessors(node)) | set(graph.successors(node))
                if graph.is_directed() else set(graph.neighbors(node))
            )
            # Keep the connector only if it actually joins two responsive
            # molecules -- otherwise it adds an unmeasured dead end.
            if len(neighbours & responsive) >= 2:
                selected.add(node)

    sub = graph.subgraph(selected).copy()

    if not keep_isolated:
        isolated = [n for n in sub.nodes if sub.degree(n) == 0]
        sub.remove_nodes_from(isolated)

    logger.info(
        f"Responsive trans-omic network: {sub.number_of_nodes()} nodes "
        f"({len(responsive)} regulated, "
        f"{sub.number_of_nodes() - len(responsive & set(sub.nodes))} connectors), "
        f"{sub.number_of_edges()} edges across layers {available_layers(sub)}"
    )
    return sub


def compare_transomic_networks(
    graph1,
    graph2,
    name1: str = "condition_1",
    name2: str = "condition_2",
) -> Dict[str, object]:
    """Compare two trans-omic networks layer by layer and relationship by relationship.

    Where a generic differential-network comparison reports shared and unique
    edges, this reports *which kinds of regulation* were gained and lost -- a
    condition that loses its allosteric edges but keeps its transcriptional ones
    is a different biological story from the reverse, and the edge-type
    breakdown is what tells them apart.

    Parameters
    ----------
    graph1, graph2 : networkx.Graph
    name1, name2 : str
        Labels used in the returned tables.

    Returns
    -------
    dict
        ``nodes_by_layer`` : pandas.DataFrame
            Per layer: node counts in each network, shared, and unique to each.
        ``edges_by_type`` : pandas.DataFrame
            Per relationship: edge counts, shared, and unique to each.
        ``regulation_shifts`` : pandas.DataFrame
            Nodes whose regulated direction differs between the two networks.
        ``summary`` : dict
            Node and edge totals plus Jaccard similarity.
    """
    def _edge_key(graph):
        return {
            (str(u), str(v), (d.get("edge_type") or "unknown"))
            for u, v, d in graph.edges(data=True)
        }

    nodes1, nodes2 = set(graph1.nodes), set(graph2.nodes)
    edges1, edges2 = _edge_key(graph1), _edge_key(graph2)

    layers = sorted(set(available_layers(graph1)) | set(available_layers(graph2)))
    layer_rows = []
    for layer in layers:
        in1 = {n for n in nodes1 if graph1.nodes[n].get("layer") == layer}
        in2 = {n for n in nodes2 if graph2.nodes[n].get("layer") == layer}
        layer_rows.append({
            "layer": layer,
            f"n_{name1}": len(in1),
            f"n_{name2}": len(in2),
            "n_shared": len(in1 & in2),
            f"unique_to_{name1}": len(in1 - in2),
            f"unique_to_{name2}": len(in2 - in1),
        })

    types = sorted(
        set(available_edge_types(graph1)) | set(available_edge_types(graph2))
    )
    type_rows = []
    for etype in types:
        in1 = {e for e in edges1 if e[2] == etype}
        in2 = {e for e in edges2 if e[2] == etype}
        type_rows.append({
            "edge_type": etype,
            f"n_{name1}": len(in1),
            f"n_{name2}": len(in2),
            "n_shared": len(in1 & in2),
            f"gained_in_{name2}": len(in2 - in1),
            f"lost_in_{name2}": len(in1 - in2),
        })

    shift_rows = []
    for node in nodes1 & nodes2:
        state1 = graph1.nodes[node].get("regulated", 0) or 0
        state2 = graph2.nodes[node].get("regulated", 0) or 0
        if state1 != state2:
            shift_rows.append({
                "node": node,
                "layer": graph1.nodes[node].get("layer"),
                "name": graph1.nodes[node].get("name"),
                name1: state1,
                name2: state2,
                "log2fc_" + name1: graph1.nodes[node].get("log2fc"),
                "log2fc_" + name2: graph2.nodes[node].get("log2fc"),
            })

    union_nodes = len(nodes1 | nodes2)
    union_edges = len(edges1 | edges2)

    result = {
        "nodes_by_layer": pd.DataFrame(layer_rows),
        "edges_by_type": pd.DataFrame(type_rows),
        "regulation_shifts": pd.DataFrame(shift_rows),
        "summary": {
            f"n_nodes_{name1}": len(nodes1),
            f"n_nodes_{name2}": len(nodes2),
            f"n_edges_{name1}": len(edges1),
            f"n_edges_{name2}": len(edges2),
            "node_jaccard": len(nodes1 & nodes2) / union_nodes if union_nodes else 0.0,
            "edge_jaccard": len(edges1 & edges2) / union_edges if union_edges else 0.0,
            "n_regulation_shifts": len(shift_rows),
        },
    }

    logger.info(
        f"Compared {name1} vs {name2}: edge Jaccard "
        f"{result['summary']['edge_jaccard']:.2f}, "
        f"{len(shift_rows)} nodes changed regulation direction"
    )
    return result
