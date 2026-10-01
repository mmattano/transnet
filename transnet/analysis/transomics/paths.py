"""Signed regulatory paths: predict a molecule's direction and check it.

The sign of a path is the product of its edge signs. Multiplied by the
measured direction of the starting molecule, it predicts the direction of the
end molecule, which is then compared with the measurement.
"""

from typing import Dict, List, Optional, Sequence

import logging

import pandas as pd

from transnet.analysis.transomics.mapping import is_measured, is_responsive
from transnet.biology.schema import (
    CURRENCY_METABOLITES,
    LAYER_HIERARCHY,
    available_layers,
    top_layer_present,
)

logger = logging.getLogger(__name__)

__all__ = [
    "trace_regulatory_paths",
    "path_consistency_summary",
]


def _edge_records(graph, u, v) -> List[dict]:
    """All edge attribute dicts between ``u`` and ``v``."""
    data = graph.get_edge_data(u, v)
    if data is None:
        return []
    if graph.is_multigraph():
        return list(data.values())
    return [data]


def trace_regulatory_paths(
    graph,
    source_layer: Optional[str] = None,
    target_layer: str = "Metabolome",
    sources: Optional[Sequence[str]] = None,
    targets: Optional[Sequence[str]] = None,
    max_length: int = 6,
    regulated_only: bool = True,
    max_paths: int = 10000,
    exclude_edge_types: Optional[Sequence[str]] = None,
    exclude_nodes: Optional[Sequence[str]] = None,
    allow_unchanged_intermediates: bool = False,
) -> pd.DataFrame:
    """Trace signed paths between layers and score them against the data.

    Parameters
    ----------
    graph : networkx.Graph
        A trans-omic network, ideally directed, with omics data mapped onto it.
    source_layer : str, optional
        Layer the paths start from. Default: the highest layer present.
    target_layer : str
        Layer the paths end in.
    sources, targets : sequence of str, optional
        Explicit start and end nodes, instead of whole layers.
    max_length : int
        Maximum number of edges in a path.
    regulated_only : bool
        Start and end only at molecules that changed. False traces paths on the
        structure alone.
    max_paths : int
        Stop after this many paths.
    exclude_edge_types : sequence of str, optional
        Edge types to skip. Default ``("protein_interaction",)``, which has
        neither direction nor sign. Pass ``()`` to keep every edge.
    exclude_nodes : sequence of str, optional
        Molecules to route around. Default: the currency metabolites (ATP,
        water, NAD and similar), which connect almost every reaction to every
        other. Excluded molecules can still be named in ``sources`` or
        ``targets``. Pass ``()`` to keep them.
    allow_unchanged_intermediates : bool
        Keep paths through molecules that were measured and did not change.
        Off by default, because the data contradict such paths. Unmeasured
        intermediates are always kept.

    Returns
    -------
    pandas.DataFrame
        One row per path: ``source``, ``target``, ``path``, ``length``,
        ``layers``, ``edge_types``; ``sign`` (product of the edge signs);
        ``unsigned_steps`` (edges with unknown sign); ``unchanged_intermediates``;
        ``source_regulated``; ``predicted`` (``sign * source_regulated``);
        ``observed``; and ``consistent`` (``predicted == observed``).
    """
    columns = [
        "source", "target", "path", "length", "layers", "edge_types",
        "sign", "unsigned_steps", "unchanged_intermediates",
        "source_regulated", "predicted", "observed", "consistent",
    ]

    layers_present = available_layers(graph)

    if source_layer is None:
        source_layer = top_layer_present(graph)
        if source_layer is None:
            logger.warning("Network has no recognised layers; cannot trace paths")
            return pd.DataFrame(columns=columns)
        if source_layer == target_layer:
            # Nothing above the target layer to trace from.
            higher = [
                layer for layer in LAYER_HIERARCHY
                if layer in layers_present and layer != target_layer
            ]
            if not higher:
                logger.warning(
                    f"Only layer present is '{target_layer}'; no regulatory "
                    f"hierarchy to trace through"
                )
                return pd.DataFrame(columns=columns)
            source_layer = higher[0]
        logger.info(
            f"No source layer given; tracing from '{source_layer}' "
            f"(layers present: {layers_present})"
        )

    for layer, role in ((source_layer, "source"), (target_layer, "target")):
        if layer not in layers_present:
            logger.warning(
                f"{role.capitalize()} layer '{layer}' absent from this network "
                f"(present: {layers_present}); no paths to trace"
            )
            return pd.DataFrame(columns=columns)

    def _endpoints(explicit, layer):
        if explicit is not None:
            return [n for n in explicit if n in graph]
        nodes = [n for n, d in graph.nodes(data=True) if d.get("layer") == layer]
        if regulated_only:
            regulated = [n for n in nodes if is_responsive(graph.nodes[n])]
            if regulated:
                return regulated
            logger.info(
                f"No regulated nodes in '{layer}'; tracing structural paths instead"
            )
        return nodes

    source_nodes = _endpoints(sources, source_layer)
    target_nodes = _endpoints(targets, target_layer)

    if not source_nodes or not target_nodes:
        logger.warning("No endpoint nodes available for path tracing")
        return pd.DataFrame(columns=columns)

    if not graph.is_directed():
        logger.warning(
            "Graph is undirected; paths cannot respect regulatory direction. "
            "Build the network with generate_graph(directed=True)."
        )

    target_set = set(target_nodes)
    rows = []
    blocked = 0

    import networkx as nx

    if exclude_edge_types is None:
        exclude_edge_types = ("protein_interaction",)
    excluded = set(exclude_edge_types)

    if exclude_nodes is None:
        exclude_nodes = CURRENCY_METABOLITES
    # An endpoint the caller asked for by name is never routed around.
    dropped_nodes = (
        set(map(str, exclude_nodes))
        - set(map(str, sources or ()))
        - set(map(str, targets or ()))
    )

    if excluded:
        keep = [
            (u, v, key) for u, v, key, data in graph.edges(keys=True, data=True)
            if (data.get("edge_type") or "unknown") not in excluded
        ] if graph.is_multigraph() else [
            (u, v) for u, v, data in graph.edges(data=True)
            if (data.get("edge_type") or "unknown") not in excluded
        ]
        search = graph.edge_subgraph(keep).copy() if keep else graph.__class__()
        removed = graph.number_of_edges() - search.number_of_edges()
        if removed:
            logger.info(
                f"Excluding {removed:,} {sorted(excluded)} edge(s) from path "
                f"tracing; they are associations rather than regulatory steps"
            )
    else:
        search = graph

    present = dropped_nodes & set(search.nodes)
    if present:
        search = search.copy()
        search.remove_nodes_from(present)
        logger.info(
            f"Routing around {len(present)} currency metabolite(s); a path "
            f"through water or ATP joins reactions that are unrelated"
        )
        source_nodes = [n for n in source_nodes if n in search]
        target_nodes = [n for n in target_nodes if n in search]
        target_set = set(target_nodes)
        if not source_nodes or not target_nodes:
            logger.warning(
                "Every endpoint was a currency metabolite; pass exclude_nodes=() "
                "to keep them"
            )
            return pd.DataFrame(columns=columns)
    for source in source_nodes:
        if len(rows) >= max_paths:
            break
        if source not in search:
            continue
        try:
            walker = nx.all_simple_paths(
                search, source, target_set, cutoff=max_length
            )
        except (nx.NodeNotFound, nx.NetworkXNoPath):
            continue

        for path in walker:
            if len(rows) >= max_paths:
                break

            sign = 1
            unsigned = 0
            edge_types = []
            ok = True

            for u, v in zip(path[:-1], path[1:]):
                records = _edge_records(search, u, v)
                if not records:
                    ok = False
                    break
                # Prefer the signed relationship where several connect the pair,
                # since that is the one carrying regulatory meaning.
                record = max(records, key=lambda r: abs(int(r.get("sign", 0) or 0)))
                edge_sign = int(record.get("sign", 0) or 0)
                edge_types.append(record.get("edge_type", "unknown"))
                if edge_sign == 0:
                    unsigned += 1
                else:
                    sign *= edge_sign

            if not ok:
                continue

            # A molecule between the endpoints that was measured and did not
            # move cannot have carried the change: the path is contradicted by
            # its own data, not evidence for the endpoint.
            unchanged = sum(
                1 for node in path[1:-1]
                if graph.nodes[node].get("layer") != "Reactions"
                and is_measured(graph.nodes[node])
                and not is_responsive(graph.nodes[node])
            )
            if unchanged and not allow_unchanged_intermediates:
                blocked += 1
                continue

            target = path[-1]
            observed = graph.nodes[target].get("regulated", 0) or 0
            source_state = int(graph.nodes[source].get("regulated", 0) or 0)
            predicted = sign * source_state
            rows.append({
                "source": source,
                "target": target,
                "path": " -> ".join(str(p) for p in path),
                "length": len(path) - 1,
                "layers": " -> ".join(
                    str(graph.nodes[p].get("layer", "?")) for p in path
                ),
                "edge_types": " -> ".join(edge_types),
                "sign": sign,
                "unsigned_steps": unsigned,
                "unchanged_intermediates": unchanged,
                "source_regulated": source_state,
                "predicted": predicted,
                "observed": observed,
                "consistent": bool(observed and predicted == observed),
            })

    table = pd.DataFrame(rows, columns=columns)

    if blocked:
        logger.info(
            f"Dropped {blocked:,} path(s) running through a molecule that was "
            f"measured and did not change; pass "
            f"allow_unchanged_intermediates=True to keep them"
        )

    if not table.empty:
        confident = table[table["unsigned_steps"] == 0]
        logger.info(
            f"Traced {len(table)} paths from {source_layer} to {target_layer}; "
            f"{len(confident)} fully signed, "
            f"{int(table['consistent'].sum())} consistent with the measured direction"
        )
    else:
        logger.info(f"No paths found from {source_layer} to {target_layer}")

    return table


def path_consistency_summary(paths: pd.DataFrame) -> pd.DataFrame:
    """One verdict per target molecule from a set of traced paths.

    Paths that share most of their steps are not independent, so consistency
    should be counted per molecule, not per path.

    Parameters
    ----------
    paths : pandas.DataFrame
        Output of :func:`trace_regulatory_paths`.

    Returns
    -------
    pandas.DataFrame
        One row per target: the number of paths, fully signed paths and
        consistent paths; ``predicted``, the direction most paths predict (0 if
        they are split); ``agrees``, whether it matches the measurement; and the
        shortest consistent path.
    """
    columns = [
        "target", "observed", "n_paths", "n_fully_signed",
        "n_consistent", "fraction_consistent", "predicted", "agrees",
        "shortest_consistent_path",
    ]
    if paths.empty:
        return pd.DataFrame(columns=columns)

    rows = []
    for target, group in paths.groupby("target"):
        consistent = group[group["consistent"]]
        shortest = (
            consistent.sort_values("length").iloc[0]["path"]
            if not consistent.empty else ""
        )
        votes = group["predicted"]
        up, down = int((votes > 0).sum()), int((votes < 0).sum())
        predicted = 1 if up > down else -1 if down > up else 0
        observed = group["observed"].iloc[0]
        rows.append({
            "target": target,
            "observed": observed,
            "n_paths": len(group),
            "n_fully_signed": int((group["unsigned_steps"] == 0).sum()),
            "n_consistent": len(consistent),
            "fraction_consistent": len(consistent) / len(group),
            "predicted": predicted,
            "agrees": bool(observed and predicted == observed),
            "shortest_consistent_path": shortest,
        })

    return (
        pd.DataFrame(rows, columns=columns)
        .sort_values("n_consistent", ascending=False)
        .reset_index(drop=True)
    )
