"""Network propagation and factor-enrichment for TransNet.

This module is the statistical-network bridge — the core unique feature that
separates TransNet from pure statistical tools (MOFA, DIABLO) and from pure
network tools (PathwayCommons, etc.).

The key insight: multi-omics integration (NMF/PCA/FA) identifies factors that
capture variance across omics layers.  Network propagation then asks: which
network regions do these factors activate?  The result is a score for every
network node that reflects both its loading in a factor *and* its topological
proximity to other highly-loaded nodes.

Workflow
--------
1. Fit ``MultiOmicsIntegrator`` on your omics data.
2. Build (or load) a ``Transnet`` network.
3. Call ``network_factor_summary(integrator, network)`` to get per-factor
   propagated scores and connectivity enrichment p-values.

Functions
---------
random_walk_with_restart
    Propagate an initial score vector through a network using RWR.
network_enrichment_permutation
    Permutation test: are the top-loading features of a factor significantly
    more connected in the network than expected by chance?
propagate_factor_loadings
    Propagate each factor's loadings through the full network, returning
    network-wide scores for every node.
network_factor_summary
    End-to-end convenience function.
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
    "propagate_factor_loadings",
    "network_factor_summary",
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
    """Propagate an initial score vector through a network using RWR.

    Random Walk with Restart (RWR) is a standard diffusion method in network
    biology.  Starting from a set of seed nodes with initial scores, it
    iteratively spreads those scores to neighbours, letting each node "absorb"
    the influence of its topological context.

    The update rule is::

        F_{t+1} = alpha * W * F_t + (1 - alpha) * F_0

    where ``W`` is the column-normalised adjacency matrix and ``F_0`` is the
    initial (restart) score vector, until convergence.

    Parameters
    ----------
    node_scores : Dict[str, float]
        Seed scores keyed by node ID.  Nodes not in the graph are silently
        ignored.  Scores are L1-normalised before propagation.
    graph : nx.Graph
        The network to propagate through.  Directed graphs are converted to
        undirected; self-loops are removed.
    alpha : float
        Restart probability (0 < alpha < 1).  Higher values weight topology
        more strongly; lower values stay closer to the seed nodes.
    max_iter : int
        Maximum number of iterations before convergence is declared by force.
    tol : float
        L-inf convergence threshold.
    as_ranks : bool
        If True, return rank-normalised scores in [0, 1] instead of raw
        propagated probabilities.

    Returns
    -------
    Dict[str, float]
        Propagated scores for every node in the graph.

    Notes
    -----
    Isolated nodes (no edges) receive only their restart probability and
    are not influenced by propagation.
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
    """Test whether a set of features is significantly connected in the network.

    Null hypothesis: the observed number of edges among ``top_features`` is no
    greater than expected for a random same-size node set drawn from the graph.

    Parameters
    ----------
    top_features : List[str]
        The feature set to test (e.g. top-loading genes/proteins for a factor).
    graph : nx.Graph
        The network.  Only nodes present in the graph are used.
    n_permutations : int
        Number of random permutations for the null distribution.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    Dict with keys:

    ``observed_edges`` : int
        Number of edges in the subgraph induced by *top_features*.
    ``mean_null_edges`` : float
        Mean edges in the permuted null sets.
    ``std_null_edges`` : float
        Standard deviation of null edges.
    ``p_value`` : float
        Empirical one-sided p-value (probability that null ≥ observed).
    ``connectivity_density`` : float
        Observed edges / maximum possible edges for the given set size.
    ``n_features_in_graph`` : int
        Number of top_features found in the graph.

    Notes
    -----
    For very small sets (< 2 nodes in graph), the test is undefined and
    ``p_value=1.0`` is returned.
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
