"""Directed, sign-aware propagation down the trans-omic hierarchy.

Ordinary network propagation treats every edge as an undirected similarity. In a
trans-omic network that throws away the two things the edges were built to
carry: a transcription factor regulates its target and not the reverse, and an
allosteric inhibitor pushes its reaction the *other* way.

:func:`hierarchical_propagation` walks the regulatory hierarchy in its own
direction and multiplies by the edge sign, so a score arriving at a reaction
through an inhibitor arrives negative.  The undirected random-walk-with-restart
in :mod:`transnet.analysis.network_propagation` remains the right tool when you
want diffusion-based similarity rather than regulatory flow.
"""

from typing import Dict

import logging

import numpy as np
import pandas as pd

from transnet.biology.schema import available_layers

logger = logging.getLogger(__name__)

__all__ = [
    "hierarchical_propagation",
    "downstream_influence",
]


def hierarchical_propagation(
    graph,
    seeds: Dict[str, float],
    alpha: float = 0.85,
    max_iter: int = 100,
    tol: float = 1e-8,
    respect_sign: bool = True,
) -> pd.DataFrame:
    """Propagate seed scores forward along directed, signed regulatory edges.

    Parameters
    ----------
    graph : networkx.Graph
        A trans-omic network.  A directed graph is strongly preferred: on an
        undirected graph this degenerates to signed diffusion and says so.
    seeds : dict
        ``{node: score}``.  Signed scores are meaningful -- seed a
        down-regulated kinase with a negative value and its inhibitory targets
        come out positive.
    alpha : float
        Share of the score that propagates onward at each step; ``1 - alpha``
        stays at the seed.
    max_iter, tol : int, float
        Power-iteration controls.
    respect_sign : bool
        Multiply by each edge's ``sign`` while propagating.  Unsigned edges
        (sign 0) pass the score through unchanged rather than zeroing it, since
        sign 0 means "unknown", not "no effect".

    Returns
    -------
    pandas.DataFrame
        Columns ``node``, ``layer``, ``name``, ``score``, ``abs_score``,
        ``predicted_direction``, sorted by absolute score.  Empty if no seed is
        in the graph.

    Examples
    --------
    >>> scores = hierarchical_propagation(G, {"P31749": 1.0})   # active AKT1
    >>> scores.head()
    """
    present = {str(n): float(v) for n, v in seeds.items() if str(n) in graph}
    missing = len(seeds) - len(present)
    if missing:
        logger.warning(
            f"{missing}/{len(seeds)} seed nodes are not in the network and were "
            f"dropped; check that the identifiers match the graph's."
        )
    if not present:
        logger.warning("No seed node found in the network")
        return pd.DataFrame(columns=[
            "node", "layer", "name", "score", "abs_score", "predicted_direction",
        ])

    if not graph.is_directed():
        logger.warning(
            "Graph is undirected; propagation cannot follow regulatory direction. "
            "Build with generate_graph(directed=True) for hierarchical flow."
        )

    nodes = list(graph.nodes)
    index = {node: i for i, node in enumerate(nodes)}
    n = len(nodes)

    # Signed transition matrix, column-normalised by out-degree so that a hub
    # does not simply dominate by fan-out.
    from scipy import sparse

    rows, cols, values = [], [], []
    out_degree = np.zeros(n)
    edge_iter = graph.out_edges(data=True) if graph.is_directed() else graph.edges(data=True)
    for u, v, data in edge_iter:
        sign = int(data.get("sign", 0) or 0) if respect_sign else 1
        if sign == 0:
            sign = 1  # unknown sign: pass the score through undisturbed
        weight = float(data.get("weight", 1.0) or 1.0)
        i, j = index[u], index[v]
        rows.append(j)
        cols.append(i)
        values.append(sign * weight)
        out_degree[i] += abs(weight)
        if not graph.is_directed():
            rows.append(i)
            cols.append(j)
            values.append(sign * weight)
            out_degree[j] += abs(weight)

    if not values:
        logger.warning("Network has no edges; returning the seeds unchanged")
        matrix = sparse.csr_matrix((n, n))
    else:
        matrix = sparse.csr_matrix((values, (rows, cols)), shape=(n, n))
        scale = np.zeros(n)
        np.divide(1.0, out_degree, out=scale, where=out_degree > 0)
        matrix = matrix @ sparse.diags(scale)

    seed_vector = np.zeros(n)
    for node, value in present.items():
        seed_vector[index[node]] = value

    scores = seed_vector.copy()
    for _ in range(max_iter):
        updated = alpha * (matrix @ scores) + (1 - alpha) * seed_vector
        if np.abs(updated - scores).sum() < tol:
            scores = updated
            break
        scores = updated

    table = pd.DataFrame({
        "node": nodes,
        "layer": [graph.nodes[node].get("layer") for node in nodes],
        "name": [graph.nodes[node].get("name") for node in nodes],
        "score": scores,
    })
    table["abs_score"] = table["score"].abs()
    table["predicted_direction"] = np.sign(table["score"]).astype(int)
    table = table[table["abs_score"] > 0].sort_values(
        "abs_score", ascending=False
    ).reset_index(drop=True)

    logger.info(
        f"Propagated {len(present)} seeds to {len(table)} nodes across "
        f"layers {available_layers(graph)}"
    )
    return table


def downstream_influence(
    graph,
    seeds: Dict[str, float],
    target_layer: str = "Metabolome",
    **kwargs,
) -> pd.DataFrame:
    """Propagate from seeds and report only what reaches one layer.

    A convenience wrapper for the usual trans-omic question: given these
    signaling or transcriptional changes, which metabolites does the network
    predict will move, and in which direction?

    Parameters
    ----------
    graph : networkx.Graph
    seeds : dict
        ``{node: score}``.
    target_layer : str
        Layer to report.  Absent layers produce an empty table and a warning.
    **kwargs
        Passed to :func:`hierarchical_propagation`.

    Returns
    -------
    pandas.DataFrame
        As :func:`hierarchical_propagation`, restricted to ``target_layer`` and
        with an ``observed`` column where measured data is present, so predicted
        and measured directions can be compared directly.
    """
    if target_layer not in available_layers(graph):
        logger.warning(
            f"Target layer '{target_layer}' absent from this network "
            f"(present: {available_layers(graph)})"
        )
        return pd.DataFrame(columns=[
            "node", "layer", "name", "score", "abs_score",
            "predicted_direction", "observed", "agrees",
        ])

    scores = hierarchical_propagation(graph, seeds, **kwargs)
    if scores.empty:
        return scores.assign(observed=[], agrees=[])

    result = scores[scores["layer"] == target_layer].copy()
    result["observed"] = [
        graph.nodes[node].get("regulated", 0) or 0 for node in result["node"]
    ]
    # Unmeasured targets are pd.NA rather than False: the network predicts a
    # direction for them, but there is nothing to agree or disagree with.
    result["agrees"] = [
        pd.NA if observed == 0 else bool(observed == predicted)
        for observed, predicted in zip(
            result["observed"], result["predicted_direction"]
        )
    ]

    measured = result[result["observed"] != 0]
    if len(measured):
        logger.info(
            f"{int(measured['agrees'].sum())}/{len(measured)} measured "
            f"{target_layer} nodes move in the predicted direction"
        )
    return result.reset_index(drop=True)
