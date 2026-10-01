"""Export a network for other tools.

* :func:`to_cytoscape_json`: Cytoscape Desktop and cytoscape.js.
* :func:`to_arena3d`: Arena3D Web, a 3-D multilayer network viewer.
* :func:`to_transomics2cytoscape`: tables for the Bioconductor package
  transomics2cytoscape, which draws stacked layers in Cytoscape.

Each takes a ``Transnet`` object or a graph. Plain CSV files are written by
:func:`transnet.io.write_network`.
"""

from __future__ import annotations

import json
import textwrap
import zipfile
from io import BytesIO
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import networkx as nx
import pandas as pd


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _get_graph(network) -> nx.Graph:
    """Accept either a NetworkX graph or a Transnet object."""
    if isinstance(network, nx.Graph):
        return network
    if hasattr(network, "graph") and isinstance(network.graph, nx.Graph):
        return network.graph
    raise TypeError(
        f"Expected a NetworkX graph or Transnet object, got {type(network)}"
    )


def _get_name(network, default: str = "transnet") -> str:
    if hasattr(network, "name") and network.name:
        return str(network.name)
    return default


def _layer_of(G: nx.Graph, node: str) -> str:
    return G.nodes[node].get("layer", "Unknown") if node in G else "Unknown"


def _unique_layers(G: nx.Graph) -> List[str]:
    return sorted({d.get("layer", "Unknown") for _, d in G.nodes(data=True)})


def _save_or_return(data: bytes, path: Optional[Union[str, Path]]) -> bytes:
    """Write *data* to *path* if given, always return the bytes."""
    if path is not None:
        Path(path).write_bytes(data)
    return data


# ---------------------------------------------------------------------------
# Arena3D Web export
# ---------------------------------------------------------------------------

# Arena3D Web JSON schema (v2):
#
#  {
#    "version": "2.0",
#    "name":    "<network name>",
#    "layers": [
#      {
#        "id":    0,
#        "name":  "Transcriptome",
#        "color": "#4C72B0",
#        "nodes": [{"id": "BRCA1", "label": "BRCA1", "size": 5, "color": "#4C72B0"}],
#        "edges": [{"source": "BRCA1", "target": "TP53", "weight": 0.9}]
#      },
#      ...
#    ],
#    "inter_layer_edges": [
#      {"source_layer": 0, "target_layer": 1,
#       "source_node": "BRCA1", "target_node": "BRCA1_PROTEIN",
#       "weight": 1.0}
#    ]
#  }

_ARENA3D_LAYER_COLORS = {
    "Transcriptome": "#4C72B0",
    "Proteome":      "#DD8452",
    "Metabolome":    "#55A868",
    "Reactions":     "#C44E52",
    "Pathways":      "#8172B2",
    "Unknown":       "#888888",
}

_ARENA3D_LAYER_ORDER = [
    "Transcriptome", "Proteome", "Metabolome", "Reactions", "Pathways"
]


def to_arena3d(
    network,
    path: Optional[Union[str, Path]] = None,
    *,
    include_inter_layer: bool = True,
    max_nodes_per_layer: int = 2000,
    score_attr: Optional[str] = None,
) -> bytes:
    """Export a TransNet network to Arena3D Web JSON format.

    Parameters
    ----------
    network : Transnet or nx.Graph
        The network to export.
    path : str or Path, optional
        If given, write the JSON to this file.
    include_inter_layer : bool
        Whether to include cross-layer edges in the ``inter_layer_edges`` list.
        Intra-layer edges are always included in their respective layer objects.
    max_nodes_per_layer : int
        Cap the number of nodes per layer (highest-degree nodes are kept).
        Keeps the JSON manageable for the web viewer.
    score_attr : str, optional
        Node attribute to use as ``size`` (e.g. ``"propagated_score"``).

    Returns
    -------
    bytes
        UTF-8 encoded JSON.
    """
    G = _get_graph(network)
    name = _get_name(network)

    layers_in_graph = _unique_layers(G)
    # Preserve preferred order, append any extra layers at the end
    layer_order = [lyr for lyr in _ARENA3D_LAYER_ORDER if lyr in layers_in_graph]
    layer_order += [lyr for lyr in layers_in_graph if lyr not in layer_order]

    layer_id_map: Dict[str, int] = {lyr: i for i, lyr in enumerate(layer_order)}

    # Collect nodes per layer, pruned to max_nodes_per_layer by degree
    nodes_per_layer: Dict[str, List[str]] = {lyr: [] for lyr in layer_order}
    for n, d in G.nodes(data=True):
        layer = d.get("layer", "Unknown")
        if layer not in nodes_per_layer:
            nodes_per_layer[layer] = []
        nodes_per_layer[layer].append(n)

    # Trim to top-degree nodes per layer
    retained: set = set()
    for layer, node_list in nodes_per_layer.items():
        if len(node_list) > max_nodes_per_layer:
            node_list = sorted(node_list, key=lambda n: G.degree(n), reverse=True)[
                :max_nodes_per_layer
            ]
        nodes_per_layer[layer] = node_list
        retained.update(node_list)

    # Build layer objects
    arena_layers: List[Dict[str, Any]] = []
    for layer in layer_order:
        color = _ARENA3D_LAYER_COLORS.get(layer, "#888888")
        node_list = nodes_per_layer.get(layer, [])

        arena_nodes = []
        for n in node_list:
            attrs = G.nodes[n]
            size = 5
            if score_attr and score_attr in attrs:
                size = max(3, min(20, int(attrs[score_attr] * 20)))
            arena_nodes.append({
                "id":    str(n),
                "label": str(attrs.get("label", n)),
                "size":  size,
                "color": color,
            })

        # Intra-layer edges only
        arena_edges = []
        for u, v, d in G.edges(data=True):
            if u in retained and v in retained:
                if _layer_of(G, u) == layer and _layer_of(G, v) == layer:
                    arena_edges.append({
                        "source": str(u),
                        "target": str(v),
                        "weight": float(d.get("weight", 1.0)),
                    })

        arena_layers.append({
            "id":    layer_id_map[layer],
            "name":  layer,
            "color": color,
            "nodes": arena_nodes,
            "edges": arena_edges,
        })

    # Build inter-layer edges
    inter_layer_edges: List[Dict[str, Any]] = []
    if include_inter_layer:
        for u, v, d in G.edges(data=True):
            if u in retained and v in retained:
                lu, lv = _layer_of(G, u), _layer_of(G, v)
                if lu != lv:
                    inter_layer_edges.append({
                        "source_layer": layer_id_map.get(lu, -1),
                        "target_layer": layer_id_map.get(lv, -1),
                        "source_node":  str(u),
                        "target_node":  str(v),
                        "weight":       float(d.get("weight", 1.0)),
                    })

    doc = {
        "version":            "2.0",
        "name":               name,
        "layers":             arena_layers,
        "inter_layer_edges":  inter_layer_edges,
    }

    data = json.dumps(doc, ensure_ascii=False, indent=2).encode("utf-8")
    return _save_or_return(data, path)


# ---------------------------------------------------------------------------
# Transomics2cytoscape export
# ---------------------------------------------------------------------------

# The Bioconductor package transomics2cytoscape
# expects the following inputs to createTransomicsNetwork():
#
#   networkDataList  — a named list of igraph objects OR data frames
#   transomicsEdge   — a data frame with columns:
#                        layer1, nodeID1, layer2, nodeID2, [weight]
#
# We export a ZIP bundle containing:
#   layers/
#     <LayerName>_nodes.tsv   — columns: id, label, [degree, ...]
#     <LayerName>_edges.tsv   — columns: source, target, weight, interaction
#   inter_layer_edges.tsv     — columns: layer1, nodeID1, layer2, nodeID2, weight
#   README.md                 — R code snippet to load with transomics2cytoscape


_T2C_README_TEMPLATE = textwrap.dedent("""\
    # Transomics2cytoscape Bundle — {network_name}

    ## Loading in R

    ```r
    library(transomics2cytoscape)
    library(RCy3)

    # Read layer networks (as igraph objects via igraph::read_graph or
    # as data frames via read.delim, then convert with igraph::graph_from_data_frame)
    layer_files <- list.files("layers", pattern="_edges\\\\.tsv$", full.names=TRUE)

    networkDataList <- lapply(layer_files, function(f) {{
      edges_df <- read.delim(f, stringsAsFactors=FALSE)
      igraph::graph_from_data_frame(edges_df, directed=FALSE)
    }})
    names(networkDataList) <- sub("_edges\\\\.tsv$", "", basename(layer_files))

    # Read inter-layer edges
    transomicsEdge <- read.delim("inter_layer_edges.tsv", stringsAsFactors=FALSE)

    # Connect to Cytoscape (must be running with CyREST on port 1234)
    cytoscapePing()

    net <- createTransomicsNetwork(
      networkDataList = networkDataList,
      transomicsEdge  = transomicsEdge
    )
    ```

    ## File descriptions

    * `layers/<LayerName>_nodes.tsv` — node attributes for each biological layer
    * `layers/<LayerName>_edges.tsv` — intra-layer edges (source, target, weight)
    * `inter_layer_edges.tsv`        — cross-layer connections
      * Columns: `layer1`, `nodeID1`, `layer2`, `nodeID2`, `weight`

    ## Layers in this bundle

    {layer_list}

    Generated by TransNet (https://github.com/mmattano/transnet)
""")


def to_transomics2cytoscape(
    network,
    output_dir: Optional[Union[str, Path]] = None,
    *,
    as_zip: bool = True,
    zip_path: Optional[Union[str, Path]] = None,
) -> Optional[bytes]:
    """Export to a Transomics2cytoscape-compatible bundle.

    Parameters
    ----------
    network : Transnet or nx.Graph
        The network to export.
    output_dir : str or Path, optional
        Directory to write files into.  Created if it does not exist.
        Ignored when *as_zip* is True and *zip_path* / return value is used.
    as_zip : bool
        If True (default) return a ZIP archive as ``bytes`` (and optionally
        write it to *zip_path*).  If False, write individual files to
        *output_dir*.
    zip_path : str or Path, optional
        When *as_zip* is True, also save the ZIP to this path.

    Returns
    -------
    bytes or None
        ZIP bytes when *as_zip* is True; None when writing to *output_dir*.
    """
    G = _get_graph(network)
    name = _get_name(network)
    layers = _unique_layers(G)

    # Partition nodes by layer
    nodes_by_layer: Dict[str, List[str]] = {lyr: [] for lyr in layers}
    for n, d in G.nodes(data=True):
        layer = d.get("layer", "Unknown")
        nodes_by_layer.setdefault(layer, []).append(n)

    # Build per-layer node and edge tables
    layer_nodes: Dict[str, pd.DataFrame] = {}
    layer_edges: Dict[str, pd.DataFrame] = {}
    for layer, node_list in nodes_by_layer.items():
        node_set = set(node_list)
        n_df = pd.DataFrame([
            {
                "id":     n,
                "label":  G.nodes[n].get("label", n),
                "degree": G.degree(n),
                "layer":  layer,
            }
            for n in node_list
        ])
        layer_nodes[layer] = n_df

        edge_rows = []
        for u, v, d in G.edges(data=True):
            if u in node_set and v in node_set:
                edge_rows.append({
                    "source":      str(u),
                    "target":      str(v),
                    "weight":      float(d.get("weight", 1.0)),
                    # The builder writes "edge_type"; "interaction" was
                    # never set by anything, so every edge exported as
                    # the placeholder "interacts".
                    "interaction": str(
                        d.get("edge_type") or d.get("interaction") or "interacts"
                    ),
                })
        layer_edges[layer] = pd.DataFrame(
            edge_rows,
            columns=["source", "target", "weight", "interaction"]
        )

    # Build inter-layer edge table
    inter_rows = []
    for u, v, d in G.edges(data=True):
        lu, lv = _layer_of(G, u), _layer_of(G, v)
        if lu != lv:
            inter_rows.append({
                "layer1":  lu,
                "nodeID1": str(u),
                "layer2":  lv,
                "nodeID2": str(v),
                "weight":  float(d.get("weight", 1.0)),
            })
    inter_df = pd.DataFrame(inter_rows, columns=["layer1", "nodeID1", "layer2", "nodeID2", "weight"])

    # README
    layer_bullet = "\n".join(
        f"* {lyr} ({len(nodes_by_layer[lyr])} nodes)" for lyr in layers
    )
    readme = _T2C_README_TEMPLATE.format(
        network_name=name,
        layer_list=layer_bullet,
    )

    def _tsv(df: pd.DataFrame) -> bytes:
        return df.to_csv(sep="\t", index=False).encode("utf-8")

    if as_zip:
        buf = BytesIO()
        with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as zf:
            for layer in layers:
                safe = layer.replace(" ", "_")
                zf.writestr(f"layers/{safe}_nodes.tsv", _tsv(layer_nodes[layer]))
                zf.writestr(f"layers/{safe}_edges.tsv", _tsv(layer_edges[layer]))
            zf.writestr("inter_layer_edges.tsv", _tsv(inter_df))
            zf.writestr("README.md", readme.encode("utf-8"))
        data = buf.getvalue()
        if zip_path is not None:
            Path(zip_path).write_bytes(data)
        return data

    else:
        out = Path(output_dir) if output_dir else Path("transomics2cytoscape_bundle")
        (out / "layers").mkdir(parents=True, exist_ok=True)
        for layer in layers:
            safe = layer.replace(" ", "_")
            layer_nodes[layer].to_csv(out / "layers" / f"{safe}_nodes.tsv", sep="\t", index=False)
            layer_edges[layer].to_csv(out / "layers" / f"{safe}_edges.tsv", sep="\t", index=False)
        inter_df.to_csv(out / "inter_layer_edges.tsv", sep="\t", index=False)
        (out / "README.md").write_text(readme, encoding="utf-8")
        return None


# ---------------------------------------------------------------------------
# Cytoscape JSON export (cytoscape.js / File → Import → Network from JSON)
# ---------------------------------------------------------------------------

def to_cytoscape_json(
    network,
    path: Optional[Union[str, Path]] = None,
    *,
    include_propagated_scores: Optional[pd.DataFrame] = None,
) -> bytes:
    """Export to Cytoscape.js JSON format.

    Compatible with Cytoscape Desktop's *File → Import → Network from File*
    (select JSON format) and with the ``cytoscape.js`` JavaScript library.

    Parameters
    ----------
    network : Transnet or nx.Graph
        The network to export.
    path : str or Path, optional
        Write the JSON to this file.
    include_propagated_scores : pd.DataFrame, optional
        If provided (nodes × factors DataFrame from network propagation),
        factor scores are embedded as node attributes for use in Cytoscape
        styles.

    Returns
    -------
    bytes
        UTF-8 encoded JSON.
    """
    G = _get_graph(network)

    elements: Dict[str, List[Dict]] = {"nodes": [], "edges": []}

    score_cols: List[str] = (
        list(include_propagated_scores.columns)
        if include_propagated_scores is not None
        else []
    )

    for n, d in G.nodes(data=True):
        node_data: Dict[str, Any] = {
            "id":     str(n),
            "label":  str(d.get("label", n)),
            "layer":  d.get("layer", "Unknown"),
            "degree": G.degree(n),
        }
        # Embed propagated scores if available
        if include_propagated_scores is not None and n in include_propagated_scores.index:
            for col in score_cols:
                node_data[f"score_{col}"] = float(include_propagated_scores.loc[n, col])
        elements["nodes"].append({"data": node_data})

    for i, (u, v, d) in enumerate(G.edges(data=True)):
        lu, lv = _layer_of(G, u), _layer_of(G, v)
        elements["edges"].append({
            "data": {
                "id":           f"e{i}",
                "source":       str(u),
                "target":       str(v),
                "weight":       float(d.get("weight", 1.0)),
                "interaction":  d.get("edge_type") or d.get("interaction") or "interacts",
                "sign":         int(d.get("sign", 0) or 0),
                "ec":           d.get("ec", ""),
                "source_layer": lu,
                "target_layer": lv,
                "cross_layer":  lu != lv,
            }
        })

    doc = {"elements": elements}
    data = json.dumps(doc, ensure_ascii=False, indent=2).encode("utf-8")
    return _save_or_return(data, path)
