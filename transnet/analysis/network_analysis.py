"""
Network analysis tools for transnet networks.
"""

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from typing import Dict, List, Optional, Tuple, Union
import logging

logger = logging.getLogger(__name__)


def _simple_undirected(G: nx.Graph) -> nx.Graph:
    """A plain undirected graph, for algorithms defined only on one.

    The trans-omic network is a directed multigraph: connectivity, community
    detection and closeness are not defined on it. Collapsing direction and
    parallel edges is the standard reading for those questions, and doing it
    here means no caller has to remember.
    """
    if G.is_multigraph() or G.is_directed():
        simple = nx.Graph()
        simple.add_nodes_from(G.nodes(data=True))
        simple.add_edges_from((u, v) for u, v in G.edges())
        return simple
    return G

def compute_network_statistics(G: nx.Graph, exact_paths_below: int = 5000) -> Dict:
    """Basic statistics for a network: size, density, degrees, connectivity.

    Parameters
    ----------
    G : networkx.Graph
        Any graph; a directed multigraph is collapsed for the measures that
        need an undirected simple graph.
    exact_paths_below : int
        Average shortest path length and diameter are quadratic in the number
        of nodes, which is hours on an interactome. They are computed only
        when the largest component is smaller than this, and reported as
        ``None`` otherwise rather than silently skipped.

    Returns
    -------
    dict
    """
    simple = _simple_undirected(G)
    stats: Dict = {
        'nodes': G.number_of_nodes(),
        'edges': G.number_of_edges(),
        'density': nx.density(simple),
    }
    if stats['nodes'] == 0:
        return stats

    degrees = np.array([d for _, d in G.degree()], dtype=float)
    stats['avg_degree'] = float(degrees.mean())
    stats['median_degree'] = float(np.median(degrees))
    stats['min_degree'] = float(degrees.min())
    stats['max_degree'] = float(degrees.max())

    components = list(nx.connected_components(simple))
    stats['connected_components'] = len(components)
    largest = max(components, key=len)
    stats['largest_component_size'] = len(largest)
    stats['largest_component_ratio'] = len(largest) / stats['nodes']

    if 1 < len(largest) < exact_paths_below:
        component = simple.subgraph(largest)
        stats['avg_path_length'] = nx.average_shortest_path_length(component)
        stats['diameter'] = nx.diameter(component)
    else:
        stats['avg_path_length'] = None
        stats['diameter'] = None
        if len(largest) >= exact_paths_below:
            logger.info(
                f"Largest component has {len(largest):,} nodes; skipping path "
                f"length and diameter (raise exact_paths_below to force them)"
            )

    stats['avg_clustering'] = nx.average_clustering(simple)
    return stats

def identify_hubs(G: nx.Graph, top_n: int = 10) -> List[Tuple[str, int]]:
    """
    Identify hub nodes based on degree.
    
    Parameters
    ----------
    G : nx.Graph
        NetworkX graph
    top_n : int
        Number of top hubs to return
        
    Returns
    -------
    List[Tuple[str, int]]
        List of (node, degree) tuples for the top hubs
    """
    degrees = dict(G.degree())
    sorted_nodes = sorted(degrees.items(), key=lambda x: x[1], reverse=True)
    return sorted_nodes[:top_n]

def compute_centrality_measures(G: nx.Graph, top_n: int = 10) -> Dict[str, List[Tuple[str, float]]]:
    """
    Compute various centrality measures for nodes.
    
    Parameters
    ----------
    G : nx.Graph
        NetworkX graph
    top_n : int
        Number of top nodes to return for each measure
        
    Returns
    -------
    Dict[str, List[Tuple[str, float]]]
        Dictionary mapping centrality measures to lists of (node, value) tuples
    """
    results = {}
    
    # Degree centrality
    simple = _simple_undirected(G)

    degree_cent = nx.degree_centrality(simple)
    results['degree_centrality'] = sorted(degree_cent.items(), key=lambda x: x[1], reverse=True)[:top_n]
    
    # Betweenness centrality
    betweenness_cent = nx.betweenness_centrality(simple)
    results['betweenness_centrality'] = sorted(betweenness_cent.items(), key=lambda x: x[1], reverse=True)[:top_n]
    
    # Closeness centrality
    closeness_cent = nx.closeness_centrality(simple)
    results['closeness_centrality'] = sorted(closeness_cent.items(), key=lambda x: x[1], reverse=True)[:top_n]
    
    # Eigenvector centrality
    try:
        eigenvector_cent = nx.eigenvector_centrality(simple, max_iter=1000)
        results['eigenvector_centrality'] = sorted(eigenvector_cent.items(), key=lambda x: x[1], reverse=True)[:top_n]
    except nx.PowerIterationFailedConvergence:
        logger.warning("Eigenvector centrality calculation did not converge")
        results['eigenvector_centrality'] = []
    
    return results

def detect_communities(G: nx.Graph, method: str = 'louvain') -> Dict[str, int]:
    """
    Detect communities in a network.
    
    Parameters
    ----------
    G : nx.Graph
        NetworkX graph
    method : str
        Community detection method ('louvain', 'label_propagation', or 'greedy_modularity')
        
    Returns
    -------
    Dict[str, int]
        Dictionary mapping nodes to community IDs
    """
    G = _simple_undirected(G)

    try:
        if method == 'louvain':
            # networkx ships Louvain since 3.0; the external ``python-louvain``
            # package this used to import is not a dependency, so the whole
            # function used to return {} on any machine without it.
            from networkx.algorithms import community
            communities = community.louvain_communities(G, seed=0)
            return {node: i for i, comm in enumerate(communities) for node in comm}

        elif method == 'label_propagation':
            from networkx.algorithms import community
            communities = community.label_propagation_communities(G)
            partition = {}
            for i, comm in enumerate(communities):
                for node in comm:
                    partition[node] = i
            return partition
        
        elif method == 'greedy_modularity':
            from networkx.algorithms import community
            communities = community.greedy_modularity_communities(G)
            partition = {}
            for i, comm in enumerate(communities):
                for node in comm:
                    partition[node] = i
            return partition
        
        else:
            logger.error(f"Unknown community detection method: {method}")
            return {}
            
    except ImportError as e:
        logger.error(f"Required package not found: {e}")
        return {}


def enrichment_analysis(node_set: List[str], 
                      background_set: List[str], 
                      annotations: Dict[str, List[str]],
                      method: str = 'hypergeometric') -> pd.DataFrame:
    """
    Perform enrichment analysis on a set of nodes.
    
    Parameters
    ----------
    node_set : List[str]
        Set of nodes to analyze
    background_set : List[str]
        Background set of nodes
    annotations : Dict[str, List[str]]
        Dictionary mapping annotation terms to lists of nodes
    method : str
        Statistical method ('hypergeometric' or 'fisher')
        
    Returns
    -------
    pd.DataFrame
        Dataframe of enrichment results
    """
    from scipy import stats
    
    results = []
    node_set = set(node_set)
    background_set = set(background_set)
    
    N = len(background_set)  # Total population
    n = len(node_set)        # Selected population
    
    for term, annotated_nodes in annotations.items():
        annotated_nodes = set(annotated_nodes)
        
        # Count statistics
        K = len(annotated_nodes.intersection(background_set))  # Total annotated
        k = len(annotated_nodes.intersection(node_set))        # Selected annotated
        
        # Skip if no overlap
        if k == 0:
            continue
        
        # Calculate p-value
        if method == 'hypergeometric':
            p_value = stats.hypergeom.sf(k-1, N, K, n)
        elif method == 'fisher':
            # Fisher's exact test
            contingency_table = [
                [k, K-k],
                [n-k, N-n-K+k]
            ]
            _, p_value = stats.fisher_exact(contingency_table)
        else:
            logger.error(f"Unknown method: {method}")
            continue
        
        # Calculate enrichment ratio
        enrichment_ratio = (k / n) / (K / N)
        
        results.append({
            'Term': term,
            'Annotated_in_selection': k,
            'Annotated_total': K,
            'Enrichment_ratio': enrichment_ratio,
            'P_value': p_value
        })
    
    # Create dataframe and sort by p-value
    result_df = pd.DataFrame(results)
    if not result_df.empty:
        # Benjamini-Hochberg. The rank must be the position after sorting, not
        # the original insertion index that sort_values carries along, and the
        # result must be made monotone so a test never gets a smaller FDR than
        # a more significant one.
        result_df = result_df.sort_values('P_value').reset_index(drop=True)
        n = len(result_df)
        ranks = np.arange(1, n + 1)
        raw = result_df['P_value'].to_numpy() * n / ranks
        result_df['FDR'] = np.minimum.accumulate(raw[::-1])[::-1].clip(max=1.0)

    return result_df


def find_active_modules(G: nx.Graph,
                        score_attr: str = 'log2fc',
                        n_modules: int = 5,
                        absolute: bool = True,
                        min_size: int = 2) -> List[List[str]]:
    """Connected groups of strongly changed molecules.

    Each measured node is scored by how far its value lies above the average
    (a z-score of ``|value|`` when ``absolute``). Nodes without a value, such as
    reactions, score 0 and act as connectors. A module starts at the highest
    scoring node and grows one neighbour at a time while the total score
    rises. A connector is added together with the best node behind it, so a
    module can run from an enzyme through its reaction to a metabolite. The
    module's nodes are then removed and the next module is grown.

    Parameters
    ----------
    G : nx.Graph
        Network with a numeric node attribute ``score_attr``.
    score_attr : str
        Node attribute holding the measurement, usually ``"log2fc"``.
    n_modules : int
        Maximum number of modules to return.
    absolute : bool
        Score increases and decreases alike. False scores only increases.
    min_size : int
        Modules smaller than this are discarded.

    Returns
    -------
    list of list of str
        Node ids per module, largest total score first.
    """
    values = {n: G.nodes[n].get(score_attr) for n in G.nodes()}
    values = {n: (abs(v) if absolute else v) for n, v in values.items()
              if v is not None and not (isinstance(v, float) and np.isnan(v))}
    if not values:
        logger.error("No scores found in the network")
        return []
    std_score = np.std(list(values.values()))
    if std_score == 0:
        logger.error("All scores are identical, can't compute Z-scores")
        return []
    mean_score = np.mean(list(values.values()))
    z_scores = {n: (values.get(n, mean_score) - mean_score) / std_score
                if n in values else 0.0 for n in G.nodes()}

    remaining = _simple_undirected(G)
    modules = []
    while len(modules) < n_modules and remaining.number_of_nodes():
        seed = max(remaining.nodes(), key=lambda n: z_scores[n])
        if z_scores[seed] <= 0:
            break
        module = _grow_module(remaining, seed, z_scores)
        remaining.remove_nodes_from(module)
        if len(module) >= min_size:
            modules.append(module)
    return modules


def _grow_module(G: nx.Graph, seed: str, z_scores: Dict[str, float]) -> List[str]:
    """Greedy growth with a one-step look-ahead through connectors."""
    module = [seed]
    members = {seed}
    while True:
        frontier = {m for n in module for m in G.neighbors(n)} - members
        best_gain, best_nodes = 0.0, None
        for candidate in frontier:
            gain = z_scores[candidate]
            if gain > best_gain:
                best_gain, best_nodes = gain, [candidate]
            if gain <= 0:
                behind = [m for m in G.neighbors(candidate) if m not in members]
                if behind:
                    after = max(behind, key=lambda m: z_scores[m])
                    two_step = gain + z_scores[after]
                    if two_step > best_gain:
                        best_gain, best_nodes = two_step, [candidate, after]
        if best_nodes is None:
            return module
        module.extend(best_nodes)
        members.update(best_nodes)


