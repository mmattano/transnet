"""
Network analysis tools for transnet networks.
"""

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from typing import Dict, List, Optional, Tuple, Union
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def compute_network_statistics(G: nx.Graph) -> Dict:
    """
    Compute basic statistics for a network.
    
    Parameters:
    -----------
    G : nx.Graph
        NetworkX graph
        
    Returns:
    --------
    Dict
        Dictionary of network statistics
    """
    stats = {}
    
    # Basic statistics
    stats['nodes'] = G.number_of_nodes()
    stats['edges'] = G.number_of_edges()
    stats['density'] = nx.density(G)
    
    # Node degree statistics
    degrees = [d for _, d in G.degree()]
    stats['avg_degree'] = np.mean(degrees)
    stats['median_degree'] = np.median(degrees)
    stats['min_degree'] = np.min(degrees)
    stats['max_degree'] = np.max(degrees)
    
    # Connected components
    components = list(nx.connected_components(G))
    stats['connected_components'] = len(components)
    largest_cc = max(components, key=len)
    stats['largest_component_size'] = len(largest_cc)
    stats['largest_component_ratio'] = len(largest_cc) / stats['nodes']
    
    # Average path length and diameter (only for largest component)
    if len(largest_cc) > 1:
        largest_cc_graph = G.subgraph(largest_cc)
        stats['avg_path_length'] = nx.average_shortest_path_length(largest_cc_graph)
        stats['diameter'] = nx.diameter(largest_cc_graph)
    else:
        stats['avg_path_length'] = 0
        stats['diameter'] = 0
    
    # Clustering coefficient
    stats['avg_clustering'] = nx.average_clustering(G)
    
    return stats

def identify_hubs(G: nx.Graph, top_n: int = 10) -> List[Tuple[str, int]]:
    """
    Identify hub nodes based on degree.
    
    Parameters:
    -----------
    G : nx.Graph
        NetworkX graph
    top_n : int
        Number of top hubs to return
        
    Returns:
    --------
    List[Tuple[str, int]]
        List of (node, degree) tuples for the top hubs
    """
    degrees = dict(G.degree())
    sorted_nodes = sorted(degrees.items(), key=lambda x: x[1], reverse=True)
    return sorted_nodes[:top_n]

def compute_centrality_measures(G: nx.Graph, top_n: int = 10) -> Dict[str, List[Tuple[str, float]]]:
    """
    Compute various centrality measures for nodes.
    
    Parameters:
    -----------
    G : nx.Graph
        NetworkX graph
    top_n : int
        Number of top nodes to return for each measure
        
    Returns:
    --------
    Dict[str, List[Tuple[str, float]]]
        Dictionary mapping centrality measures to lists of (node, value) tuples
    """
    results = {}
    
    # Degree centrality
    degree_cent = nx.degree_centrality(G)
    results['degree_centrality'] = sorted(degree_cent.items(), key=lambda x: x[1], reverse=True)[:top_n]
    
    # Betweenness centrality
    betweenness_cent = nx.betweenness_centrality(G)
    results['betweenness_centrality'] = sorted(betweenness_cent.items(), key=lambda x: x[1], reverse=True)[:top_n]
    
    # Closeness centrality
    closeness_cent = nx.closeness_centrality(G)
    results['closeness_centrality'] = sorted(closeness_cent.items(), key=lambda x: x[1], reverse=True)[:top_n]
    
    # Eigenvector centrality
    try:
        eigenvector_cent = nx.eigenvector_centrality(G, max_iter=1000)
        results['eigenvector_centrality'] = sorted(eigenvector_cent.items(), key=lambda x: x[1], reverse=True)[:top_n]
    except nx.PowerIterationFailedConvergence:
        logger.warning("Eigenvector centrality calculation did not converge")
        results['eigenvector_centrality'] = []
    
    return results

def detect_communities(G: nx.Graph, method: str = 'louvain') -> Dict[str, int]:
    """
    Detect communities in a network.
    
    Parameters:
    -----------
    G : nx.Graph
        NetworkX graph
    method : str
        Community detection method ('louvain', 'label_propagation', or 'greedy_modularity')
        
    Returns:
    --------
    Dict[str, int]
        Dictionary mapping nodes to community IDs
    """
    try:
        if method == 'louvain':
            import community as community_louvain
            partition = community_louvain.best_partition(G)
            return partition
        
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

def plot_network(G: nx.Graph, 
                node_colors: Optional[Dict[str, str]] = None, 
                node_sizes: Optional[Dict[str, float]] = None,
                edge_colors: Optional[Dict[Tuple[str, str], str]] = None,
                title: str = "Network Visualization",
                layout: str = "spring",
                figsize: Tuple[int, int] = (12, 10)) -> plt.Figure:
    """
    Plot a network with customizable node and edge colors/sizes.
    
    Parameters:
    -----------
    G : nx.Graph
        NetworkX graph
    node_colors : Dict[str, str], optional
        Dictionary mapping nodes to colors
    node_sizes : Dict[str, float], optional
        Dictionary mapping nodes to sizes
    edge_colors : Dict[Tuple[str, str], str], optional
        Dictionary mapping edges to colors
    title : str
        Plot title
    layout : str
        Layout algorithm ('spring', 'circular', 'kamada_kawai', 'spectral', 'random', 'shell')
    figsize : Tuple[int, int]
        Figure size
        
    Returns:
    --------
    plt.Figure
        Matplotlib figure
    """
    plt.figure(figsize=figsize)
    
    # Choose layout
    if layout == 'spring':
        pos = nx.spring_layout(G)
    elif layout == 'circular':
        pos = nx.circular_layout(G)
    elif layout == 'kamada_kawai':
        pos = nx.kamada_kawai_layout(G)
    elif layout == 'spectral':
        pos = nx.spectral_layout(G)
    elif layout == 'random':
        pos = nx.random_layout(G)
    elif layout == 'shell':
        pos = nx.shell_layout(G)
    else:
        logger.warning(f"Unknown layout: {layout}, using spring layout")
        pos = nx.spring_layout(G)
    
    # Set default node colors (by layer)
    if node_colors is None:
        node_colors = {}
        layer_colors = {
            'Transcriptome': 'blue',
            'Proteome': 'red',
            'Metabolome': 'green',
            'Reactions': 'orange'
        }
        
        for node in G.nodes():
            layer = G.nodes[node].get('layer', 'Unknown')
            node_colors[node] = layer_colors.get(layer, 'gray')
    
    # Set default node sizes
    if node_sizes is None:
        node_sizes = {node: 100 for node in G.nodes()}
    
    # Set default edge colors
    if edge_colors is None:
        edge_colors = {edge: 'black' for edge in G.edges()}
    
    # Draw the network
    nx.draw_networkx_nodes(G, pos, 
                          node_color=[node_colors.get(node, 'gray') for node in G.nodes()],
                          node_size=[node_sizes.get(node, 100) for node in G.nodes()])
    
    nx.draw_networkx_edges(G, pos, 
                          edge_color=[edge_colors.get(edge, 'black') for edge in G.edges()])
    
    nx.draw_networkx_labels(G, pos, font_size=8)
    
    plt.title(title)
    plt.axis('off')
    
    return plt.gcf()

def differential_network_analysis(G1: nx.Graph, G2: nx.Graph) -> Dict:
    """
    Perform differential network analysis between two networks.
    
    Parameters:
    -----------
    G1 : nx.Graph
        First network
    G2 : nx.Graph
        Second network
        
    Returns:
    --------
    Dict
        Dictionary of differential network metrics
    """
    results = {}
    
    # Common nodes and edges
    common_nodes = set(G1.nodes()).intersection(set(G2.nodes()))
    results['common_nodes'] = len(common_nodes)
    results['nodes_only_in_G1'] = len(set(G1.nodes()) - set(G2.nodes()))
    results['nodes_only_in_G2'] = len(set(G2.nodes()) - set(G1.nodes()))
    
    common_edges = set(G1.edges()).intersection(set(G2.edges()))
    results['common_edges'] = len(common_edges)
    results['edges_only_in_G1'] = len(set(G1.edges()) - set(G2.edges()))
    results['edges_only_in_G2'] = len(set(G2.edges()) - set(G1.edges()))
    
    # Compute differential node metrics for common nodes
    if common_nodes:
        # Degree difference
        degree_diff = {}
        for node in common_nodes:
            degree_diff[node] = G2.degree(node) - G1.degree(node)
        
        results['degree_diff'] = degree_diff
        
        # Betweenness centrality difference
        bc1 = nx.betweenness_centrality(G1)
        bc2 = nx.betweenness_centrality(G2)
        bc_diff = {}
        for node in common_nodes:
            bc_diff[node] = bc2.get(node, 0) - bc1.get(node, 0)
        
        results['betweenness_diff'] = bc_diff
        
        # Closeness centrality difference
        cc1 = nx.closeness_centrality(G1)
        cc2 = nx.closeness_centrality(G2)
        cc_diff = {}
        for node in common_nodes:
            cc_diff[node] = cc2.get(node, 0) - cc1.get(node, 0)
        
        results['closeness_diff'] = cc_diff
    
    return results

def enrichment_analysis(node_set: List[str], 
                      background_set: List[str], 
                      annotations: Dict[str, List[str]],
                      method: str = 'hypergeometric') -> pd.DataFrame:
    """
    Perform enrichment analysis on a set of nodes.
    
    Parameters:
    -----------
    node_set : List[str]
        Set of nodes to analyze
    background_set : List[str]
        Background set of nodes
    annotations : Dict[str, List[str]]
        Dictionary mapping annotation terms to lists of nodes
    method : str
        Statistical method ('hypergeometric' or 'fisher')
        
    Returns:
    --------
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
        # Calculate adjusted p-values (Benjamini-Hochberg)
        result_df = result_df.sort_values('P_value')
        result_df['FDR'] = result_df['P_value'] * len(result_df) / (result_df.index + 1)
        result_df = result_df.sort_values('FDR')
    
    return result_df

def integrate_experimental_data(G: nx.Graph, 
                               df: pd.DataFrame, 
                               id_column: str,
                               value_column: str,
                               p_value_column: Optional[str] = None,
                               p_value_threshold: float = 0.05) -> nx.Graph:
    """
    Integrate experimental data into a network.
    
    Parameters:
    -----------
    G : nx.Graph
        NetworkX graph
    df : pd.DataFrame
        Dataframe with experimental data
    id_column : str
        Column name with node identifiers
    value_column : str
        Column name with values (e.g., fold change)
    p_value_column : str, optional
        Column name with p-values
    p_value_threshold : float
        P-value threshold for significance
        
    Returns:
    --------
    nx.Graph
        Graph with integrated experimental data
    """
    # Create a copy of the graph to avoid modifying the original
    G_copy = G.copy()
    
    # Create a mapping of node IDs to data values
    data_dict = {}
    for _, row in df.iterrows():
        node_id = str(row[id_column])
        value = row[value_column]
        
        # Skip if p-value is significant (if provided)
        if p_value_column and pd.notna(row[p_value_column]):
            if row[p_value_column] > p_value_threshold:
                continue
        
        data_dict[node_id] = value
    
    # Add data as node attributes
    for node in G_copy.nodes():
        node_str = str(node)
        if node_str in data_dict:
            G_copy.nodes[node]['experimental_value'] = data_dict[node_str]
            # Add a flag for significant nodes
            G_copy.nodes[node]['significant'] = True
        else:
            G_copy.nodes[node]['experimental_value'] = None
            G_copy.nodes[node]['significant'] = False
    
    return G_copy

def find_active_modules(G: nx.Graph, 
                       score_attr: str = 'experimental_value',
                       algorithm: str = 'greedy',
                       n_modules: int = 5) -> List[List[str]]:
    """
    Find active modules (subnetworks) based on node scores.
    
    Parameters:
    -----------
    G : nx.Graph
        NetworkX graph with node scores
    score_attr : str
        Node attribute containing scores
    algorithm : str
        Algorithm for finding modules ('greedy' or 'simulated_annealing')
    n_modules : int
        Number of modules to return
        
    Returns:
    --------
    List[List[str]]
        List of modules (each module is a list of node IDs)
    """
    # Convert absolute scores to Z-scores
    scores = []
    for node in G.nodes():
        score = G.nodes[node].get(score_attr)
        if score is not None:
            scores.append(score)
    
    if not scores:
        logger.error("No scores found in the network")
        return []
    
    mean_score = np.mean(scores)
    std_score = np.std(scores)
    
    if std_score == 0:
        logger.error("All scores are identical, can't compute Z-scores")
        return []
    
    for node in G.nodes():
        score = G.nodes[node].get(score_attr)
        if score is not None:
            G.nodes[node]['z_score'] = (score - mean_score) / std_score
        else:
            G.nodes[node]['z_score'] = 0.0
    
    if algorithm == 'greedy':
        return _find_modules_greedy(G, n_modules)
    elif algorithm == 'simulated_annealing':
        return _find_modules_sa(G, n_modules)
    else:
        logger.error(f"Unknown algorithm: {algorithm}")
        return []

def _find_modules_greedy(G: nx.Graph, n_modules: int) -> List[List[str]]:
    """
    Find active modules using a greedy algorithm.
    
    Parameters:
    -----------
    G : nx.Graph
        NetworkX graph with z_score node attribute
    n_modules : int
        Number of modules to return
        
    Returns:
    --------
    List[List[str]]
        List of modules (each module is a list of node IDs)
    """
    modules = []
    remaining_graph = G.copy()
    
    for _ in range(n_modules):
        if remaining_graph.number_of_nodes() == 0:
            break
        
        # Start from the node with the highest score
        seed_node = max(remaining_graph.nodes(), 
                       key=lambda n: remaining_graph.nodes[n].get('z_score', 0))
        
        current_module = [seed_node]
        current_score = remaining_graph.nodes[seed_node].get('z_score', 0)
        
        # Grow the module greedily
        while True:
            # Get neighbors of the current module
            neighbors = set()
            for node in current_module:
                neighbors.update(set(remaining_graph.neighbors(node)))
            
            # Remove nodes already in the module
            neighbors -= set(current_module)
            
            if not neighbors:
                break
            
            # Find the neighbor that increases the score the most
            best_node = None
            best_score_increase = 0
            
            for neighbor in neighbors:
                new_score = current_score + remaining_graph.nodes[neighbor].get('z_score', 0)
                score_increase = new_score - current_score
                
                if score_increase > best_score_increase:
                    best_score_increase = score_increase
                    best_node = neighbor
            
            # Stop if no neighbor improves the score
            if best_score_increase <= 0:
                break
            
            # Add the best neighbor to the module
            current_module.append(best_node)
            current_score += best_score_increase
        
        # Add the module to the results
        modules.append(current_module)
        
        # Remove the module nodes from the remaining graph
        remaining_graph.remove_nodes_from(current_module)
    
    return modules

def _find_modules_sa(G: nx.Graph, n_modules: int, 
                    max_iter: int = 1000, 
                    temp_init: float = 1.0,
                    temp_final: float = 0.01) -> List[List[str]]:
    """
    Find active modules using simulated annealing.
    
    Parameters:
    -----------
    G : nx.Graph
        NetworkX graph with z_score node attribute
    n_modules : int
        Number of modules to return
    max_iter : int
        Maximum number of iterations
    temp_init : float
        Initial temperature
    temp_final : float
        Final temperature
        
    Returns:
    --------
    List[List[str]]
        List of modules (each module is a list of node IDs)
    """
    import random
    import math
    
    # Initialize with random modules
    all_nodes = list(G.nodes())
    modules = []
    
    # Start with n_modules random seed nodes
    if len(all_nodes) < n_modules:
        seeds = all_nodes
    else:
        seeds = random.sample(all_nodes, n_modules)
    
    for seed in seeds:
        modules.append([seed])
    
    # Calculate initial score
    def module_score(module):
        return sum(G.nodes[n].get('z_score', 0) for n in module) / len(module) if module else 0
    
    def total_score(modules):
        return sum(module_score(m) for m in modules)
    
    current_score = total_score(modules)
    best_modules = [list(m) for m in modules]
    best_score = current_score
    
    # Simulated annealing
    for i in range(max_iter):
        # Calculate current temperature
        temp = temp_init * (temp_final / temp_init) ** (i / max_iter)
        
        # Choose a random modification
        r = random.random()
        
        if r < 0.3 and modules:  # Remove a node
            # Choose a random module
            if not any(modules):  # Skip if all modules are empty
                continue
                
            module_idx = random.choice([i for i, m in enumerate(modules) if m])
            if not modules[module_idx]:
                continue
                
            # Remove a random node
            node_idx = random.randrange(len(modules[module_idx]))
            removed_node = modules[module_idx].pop(node_idx)
            
            # Calculate new score
            new_score = total_score(modules)
            
            # Accept or reject
            if new_score > current_score or random.random() < math.exp((new_score - current_score) / temp):
                current_score = new_score
                if new_score > best_score:
                    best_score = new_score
                    best_modules = [list(m) for m in modules]
            else:
                # Undo change
                modules[module_idx].append(removed_node)
        
        elif r < 0.7:  # Add a node
            # Choose a random module
            module_idx = random.randrange(len(modules))
            
            # Find candidate nodes (neighbors of the module)
            module_nodes = set(modules[module_idx])
            candidates = set()
            
            for node in module_nodes:
                candidates.update(G.neighbors(node))
            
            # Remove nodes already in any module
            all_module_nodes = set().union(*modules)
            candidates -= all_module_nodes
            
            if not candidates:
                continue
                
            # Add a random candidate
            new_node = random.choice(list(candidates))
            modules[module_idx].append(new_node)
            
            # Calculate new score
            new_score = total_score(modules)
            
            # Accept or reject
            if new_score > current_score or random.random() < math.exp((new_score - current_score) / temp):
                current_score = new_score
                if new_score > best_score:
                    best_score = new_score
                    best_modules = [list(m) for m in modules]
            else:
                # Undo change
                modules[module_idx].pop()
        
        else:  # Move a node between modules
            if len(modules) < 2:
                continue
                
            # Choose source and destination modules
            src_candidates = [i for i, m in enumerate(modules) if m]
            if not src_candidates:
                continue
                
            src_idx = random.choice(src_candidates)
            dst_idx = random.choice([i for i in range(len(modules)) if i != src_idx])
            
            # Move a random node
            if not modules[src_idx]:
                continue
                
            node_idx = random.randrange(len(modules[src_idx]))
            moved_node = modules[src_idx].pop(node_idx)
            modules[dst_idx].append(moved_node)
            
            # Calculate new score
            new_score = total_score(modules)
            
            # Accept or reject
            if new_score > current_score or random.random() < math.exp((new_score - current_score) / temp):
                current_score = new_score
                if new_score > best_score:
                    best_score = new_score
                    best_modules = [list(m) for m in modules]
            else:
                # Undo change
                modules[dst_idx].pop()
                modules[src_idx].append(moved_node)
    
    # Filter out empty modules and sort by score
    best_modules = [m for m in best_modules if m]
    best_modules.sort(key=module_score, reverse=True)
    
    return best_modules[:n_modules]
    