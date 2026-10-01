"""Diffusion on the network, and a permutation test for connectedness.

random_walk_with_restart
    Spread seed scores over the network until they settle.
network_enrichment_permutation
    Test whether a set of molecules is more connected than random sets of
    the same size.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import networkx as nx
from scipy import sparse
from scipy.stats import rankdata

logger = logging.getLogger(__name__)

__all__ = [
    "random_walk_with_restart",
    "network_enrichment_permutation",
]


# ---------------------------------------------------------------------------
# Random Walk with Restart
# ---------------------------------------------------------------------------

def random_walk_with_restart(
    node_scores: Dict[str, float],
    graph: nx.Graph,
    alpha: float = 0.85,
    max_iter: int = 200,
    tol: float = 1e-8,
    as_ranks: bool = False,
) -> Dict[str, float]:
    """Spread seed scores over the network by random walk with restart.

    Iterates ``F = alpha * W * F + (1 - alpha) * F0`` until it converges, where
    ``W`` is the column-normalised adjacency matrix and ``F0`` the normalised
    seed scores. Nodes close to many seeds end with high scores. Direction and
    sign are ignored.

    Parameters
    ----------
    node_scores : dict
        Seed scores by node id. Nodes not in the graph are ignored.
    graph : networkx.Graph
        Directed graphs are treated as undirected.
    alpha : float
        Probability of continuing the walk rather than restarting (0-1).
    max_iter : int
        Maximum number of iterations.
    tol : float
        Convergence threshold.
    as_ranks : bool
        Return scores as ranks scaled to 0-1.

    Returns
    -------
    dict
        Score for every node in the graph.
    """
    if not (0 < alpha < 1):
        raise ValueError(f"alpha must be in (0, 1), got {alpha}")

    # Work on undirected graph without self-loops
    G = nx.Graph(graph)
    G.remove_edges_from(nx.selfloop_edges(G))

    if G.number_of_nodes() == 0:
        return {}

    nodes = list(G.nodes())
    n = len(nodes)
    node_idx = {node: i for i, node in enumerate(nodes)}

    # Build column-normalised adjacency (W)
    adj = nx.to_scipy_sparse_array(G, nodelist=nodes, format="csr", dtype=float)
    col_sums = np.asarray(adj.sum(axis=0)).flatten()
    # Avoid division by zero for isolated nodes
    col_sums[col_sums == 0] = 1.0
    # Column-normalise: D^{-1} A
    inv_d = sparse.diags(1.0 / col_sums)
    W = adj.dot(inv_d)

    # Build F0 from node_scores — only seeds present in the graph
    f0 = np.zeros(n)
    for node, score in node_scores.items():
        if node in node_idx:
            f0[node_idx[node]] = max(0.0, float(score))

    total = f0.sum()
    if total == 0:
        logger.warning(
            "random_walk_with_restart: no seed nodes found in graph. "
            "Returning uniform scores."
        )
        return {node: 1.0 / n for node in nodes}

    f0 /= total  # L1-normalise

    # Iterate
    f = f0.copy()
    for iteration in range(max_iter):
        f_prev = f
        f = alpha * W.dot(f) + (1.0 - alpha) * f0
        if np.max(np.abs(f - f_prev)) < tol:
            logger.debug(f"RWR converged in {iteration + 1} iterations")
            break
    else:
        logger.debug(f"RWR reached max_iter={max_iter} without converging")

    if as_ranks:
        # Convert to [0,1] ranks (1 = highest score)
        ranks = rankdata(f) / n
        return dict(zip(nodes, ranks))

    return dict(zip(nodes, f.tolist()))


# ---------------------------------------------------------------------------
# Network enrichment via permutation
# ---------------------------------------------------------------------------

def network_enrichment_permutation(
    top_features: List[str],
    graph: nx.Graph,
    n_permutations: int = 1000,
    seed: int = 42,
) -> Dict:
    """Test whether a set of molecules is more connected than random sets.

    Counts the edges among ``top_features`` and compares with random sets of
    the same size drawn from the graph.

    Parameters
    ----------
    top_features : list of str
        The molecules to test.
    graph : networkx.Graph
    n_permutations : int
        Number of random sets.
    seed : int
        Random seed.

    Returns
    -------
    dict
        ``observed_edges``, ``mean_null_edges``, ``std_null_edges``,
        ``p_value`` (one-sided), ``connectivity_density`` and
        ``n_features_in_graph``.
    """
    G = nx.Graph(graph)
    all_nodes = list(G.nodes())

    # Restrict to features present in graph
    in_graph = [f for f in top_features if f in G]
    n_features = len(in_graph)

    if n_features < 2:
        logger.warning(
            f"network_enrichment_permutation: only {n_features} feature(s) "
            f"found in graph (need ≥ 2). Returning p_value=1.0."
        )
        return {
            "observed_edges": 0,
            "mean_null_edges": 0.0,
            "std_null_edges": 0.0,
            "p_value": 1.0,
            "connectivity_density": 0.0,
            "n_features_in_graph": n_features,
        }

    # Observed: edges in induced subgraph
    observed_edges = G.subgraph(in_graph).number_of_edges()
    max_possible = n_features * (n_features - 1) / 2
    connectivity_density = observed_edges / max_possible if max_possible > 0 else 0.0

    # Null distribution via permutation
    rng = np.random.default_rng(seed)
    null_edges = np.empty(n_permutations, dtype=int)

    for i in range(n_permutations):
        null_set = rng.choice(all_nodes, size=n_features, replace=False).tolist()
        null_edges[i] = G.subgraph(null_set).number_of_edges()

    # Empirical p-value (one-sided: how often does null ≥ observed?)
    p_value = (np.sum(null_edges >= observed_edges) + 1) / (n_permutations + 1)

    return {
        "observed_edges": int(observed_edges),
        "mean_null_edges": float(np.mean(null_edges)),
        "std_null_edges": float(np.std(null_edges)),
        "p_value": float(p_value),
        "connectivity_density": float(connectivity_density),
        "n_features_in_graph": n_features,
    }


# ---------------------------------------------------------------------------
# Propagate factor loadings
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# End-to-end summary
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _apply_bh_fdr(df: pd.DataFrame, p_col: str = "p_value") -> pd.DataFrame:
    """Add a BH-FDR corrected column to a DataFrame.

    The Benjamini-Hochberg procedure controls the false discovery rate at the
    specified level by adjusting each p-value as::

        fdr_i = p_i * n / rank_i

    where rank is ascending by p-value and the result is capped at 1.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing a p-value column.
    p_col : str
        Name of the p-value column.

    Returns
    -------
    pd.DataFrame
        Copy of *df* with an additional ``fdr`` column.
    """
    df = df.copy()
    n = len(df)
    if n == 0:
        df["fdr"] = pd.Series(dtype=float)
        return df

    sorted_idx = df[p_col].argsort()
    ranks = np.empty(n)
    ranks[sorted_idx] = np.arange(1, n + 1)

    fdr_values = np.minimum(df[p_col].values * n / ranks, 1.0)

    # Ensure monotonicity (cumulative minimum from largest p)
    fdr_sorted = fdr_values[sorted_idx]
    for i in range(n - 2, -1, -1):
        fdr_sorted[i] = min(fdr_sorted[i], fdr_sorted[i + 1])
    fdr_values[sorted_idx] = fdr_sorted

    df["fdr"] = fdr_values
    return df
