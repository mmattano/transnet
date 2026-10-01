"""Read and write typed trans-omic networks as CSV.

A saved network is two files -- ``interactions.csv`` carrying the typed edges
and ``nodes.csv`` carrying node identity -- written by
:meth:`transnet.Transnet.save_network` and read back here into a
:class:`networkx.MultiDiGraph` with the full schema intact.
"""

import logging
import os
from typing import Optional

import networkx as nx
import pandas as pd

from transnet.biology.schema import INTERACTION_COLUMNS, available_edge_types

logger = logging.getLogger(__name__)

__all__ = ["read_network", "write_network"]


#: Legacy column names, from networks written before the typed schema.
_LEGACY_COLUMNS = {
    "SourceNode": "source",
    "TargetNode": "target",
    "SourceLayer": "source_layer",
    "TargetLayer": "target_layer",
    "InteractionType": "edge_type",
    "Weight": "weight",
    "Comment": "edge_type",
}


def read_network(path: str, nodes_file: Optional[str] = None) -> nx.MultiDiGraph:
    """Load a saved trans-omic network into a typed directed multigraph.

    Parameters
    ----------
    path : str
        Directory containing ``interactions.csv`` (and optionally
        ``nodes.csv``), or the path to an interactions CSV directly.
    nodes_file : str, optional
        Explicit path to the node table, when it is not beside the edges.

    Returns
    -------
    networkx.MultiDiGraph
        Nodes carry ``layer``, ``node_type``, ``name`` and, for reactions,
        ``reversible``.  Edges carry the full edge schema.

    Notes
    -----
    Networks written before the typed schema are read too: their columns are
    renamed, and the missing ``sign`` / ``role`` / ``stoichiometry`` fields are
    filled with neutral defaults.  Such a network will report every edge as
    sign 0, which analyses treat as "unknown", not "no effect" -- rebuild it
    with :meth:`transnet.Transnet.save_network` to recover the signs.
    """
    if os.path.isdir(path):
        edges_path = os.path.join(path, "interactions.csv")
        default_nodes = os.path.join(path, "nodes.csv")
    else:
        edges_path = path
        default_nodes = os.path.join(os.path.dirname(path), "nodes.csv")

    if not os.path.exists(edges_path):
        raise FileNotFoundError(f"No interactions file at {edges_path}")

    edges = pd.read_csv(edges_path, low_memory=False)
    edges = edges.rename(columns={
        old: new for old, new in _LEGACY_COLUMNS.items() if old in edges.columns
    })

    missing = [c for c in ("source", "target") if c not in edges.columns]
    if missing:
        raise ValueError(
            f"{edges_path} is missing required column(s) {missing}; "
            f"found {list(edges.columns)}"
        )

    legacy = [c for c in INTERACTION_COLUMNS if c not in edges.columns]
    if legacy:
        logger.warning(
            f"{edges_path} predates the typed edge schema (missing {legacy}). "
            f"Rebuild it with Transnet.save_network() to recover edge signs "
            f"and roles."
        )
    for column, default in (
        ("edge_type", "unknown"), ("role", ""), ("sign", 0),
        ("directed", True), ("weight", 1.0), ("stoichiometry", None),
        ("ec", ""), ("source_db", ""), ("evidence", ""), ("confidence", None),
        ("source_layer", "Unknown"), ("target_layer", "Unknown"),
    ):
        if column not in edges.columns:
            edges[column] = default

    graph = nx.MultiDiGraph()

    nodes_path = nodes_file or default_nodes
    if os.path.exists(nodes_path):
        nodes = pd.read_csv(nodes_path, low_memory=False)
        for _, row in nodes.iterrows():
            attrs = {
                "layer": row.get("Layer"),
                "node_type": row.get("Type"),
                "name": row.get("Name"),
            }
            # The gene symbol is what labels a protein in a figure; keep it
            # when the table has one (older exports do not).
            symbol = row.get("Symbol")
            if symbol is not None and pd.notna(symbol):
                attrs["symbol"] = str(symbol)
            reversible = row.get("Reversible")
            if pd.notna(reversible):
                attrs["reversible"] = bool(reversible)
            graph.add_node(str(row["ID"]), **attrs)
    else:
        logger.info(f"No node table at {nodes_path}; layers taken from the edges")

    for _, row in edges.iterrows():
        for node, layer in (
            (row["source"], row["source_layer"]),
            (row["target"], row["target_layer"]),
        ):
            node = str(node)
            if node not in graph:
                graph.add_node(node, layer=layer, node_type=layer, name=node)

        attributes = {
            column: row[column] for column in INTERACTION_COLUMNS
            if column not in ("source", "target", "source_layer", "target_layer")
        }
        attributes["sign"] = int(attributes.get("sign", 0) or 0)
        graph.add_edge(
            str(row["source"]), str(row["target"]),
            key=row["edge_type"], **attributes,
        )

    logger.info(
        f"Loaded network from {path}: {graph.number_of_nodes()} nodes, "
        f"{graph.number_of_edges()} edges, "
        f"edge types {list(available_edge_types(graph))}"
    )
    return graph


def write_network(graph: nx.Graph, path: str) -> None:
    """Write a trans-omic graph to ``interactions.csv`` + ``nodes.csv``.

    The inverse of :func:`read_network`, for graphs that were modified after
    loading (a responsive subnetwork, say) rather than rebuilt from layers.

    Parameters
    ----------
    graph : networkx.Graph
    path : str
        Directory to write into; created if missing.
    """
    os.makedirs(path, exist_ok=True)

    edge_rows = []
    for u, v, data in graph.edges(data=True):
        row = {"source": u, "target": v,
               "source_layer": graph.nodes[u].get("layer"),
               "target_layer": graph.nodes[v].get("layer")}
        for column in INTERACTION_COLUMNS:
            if column not in row:
                row[column] = data.get(column)
        edge_rows.append(row)

    pd.DataFrame(edge_rows, columns=INTERACTION_COLUMNS).to_csv(
        os.path.join(path, "interactions.csv"), index=False
    )

    node_rows = [
        {
            "ID": node,
            "Name": data.get("name"),
            "Symbol": data.get("symbol"),
            "Type": data.get("node_type"),
            "Layer": data.get("layer"),
            "FC": data.get("log2fc"),
            "P_value": data.get("qvalue"),
            "Reversible": data.get("reversible"),
        }
        for node, data in graph.nodes(data=True)
    ]
    pd.DataFrame(node_rows).to_csv(
        os.path.join(path, "nodes.csv"), index=False
    )

    logger.info(
        f"Wrote {graph.number_of_nodes()} nodes and "
        f"{graph.number_of_edges()} edges to {path}"
    )
