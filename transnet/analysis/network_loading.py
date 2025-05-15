"""
Network loading functionality for transnet networks.
"""

import os
import pandas as pd
import networkx as nx
from typing import Dict, Optional, Tuple, List, Any

def load_network_from_dataframe(
    edge_df: pd.DataFrame,
    node_df: Optional[pd.DataFrame] = None,
    source_col: str = "SourceNode",
    target_col: str = "TargetNode",
    weight_col: Optional[str] = "Weight",
    source_layer_col: Optional[str] = "SourceLayer",
    target_layer_col: Optional[str] = "TargetLayer",
    edge_type_col: Optional[str] = "Comment"
) -> nx.Graph:
    """
    Load a network from a DataFrame into a NetworkX graph.
    
    Parameters:
    -----------
    edge_df : pd.DataFrame
        DataFrame containing edge information
    node_df : pd.DataFrame, optional
        DataFrame containing node information (if available)
    source_col : str
        Column name for source nodes
    target_col : str
        Column name for target nodes
    weight_col : str, optional
        Column name for edge weights
    source_layer_col : str, optional
        Column name for source node layers
    target_layer_col : str, optional
        Column name for target node layers
    edge_type_col : str, optional
        Column name for edge types
        
    Returns:
    --------
    nx.Graph
        NetworkX graph representation of the network
    """
    G = nx.Graph()
    
    # Add edges with attributes
    for _, row in edge_df.iterrows():
        source = row[source_col]
        target = row[target_col]
        
        # Collect edge attributes
        attrs = {}
        if weight_col and weight_col in edge_df.columns:
            attrs['weight'] = row[weight_col]
        
        if edge_type_col and edge_type_col in edge_df.columns:
            attrs['type'] = row[edge_type_col]
            
        # Add edge with attributes
        G.add_edge(source, target, **attrs)
        
        # Add node layer information
        if source_layer_col and source_layer_col in edge_df.columns:
            G.nodes[source]['layer'] = row[source_layer_col]
            
        if target_layer_col and target_layer_col in edge_df.columns:
            G.nodes[target]['layer'] = row[target_layer_col]
    
    # Add additional node attributes from node_df if available
    if node_df is not None:
        node_id_col = node_df.columns[0]  # Assume first column is node ID
        
        for _, row in node_df.iterrows():
            node_id = row[node_id_col]
            if node_id in G:
                for col in node_df.columns:
                    if col != node_id_col:
                        G.nodes[node_id][col] = row[col]
    
    return G

def load_network_from_files(
    edge_file: str,
    node_file: Optional[str] = None,
    delimiter: str = '\t',
    **kwargs
) -> nx.Graph:
    """
    Load a network from files into a NetworkX graph.
    
    Parameters:
    -----------
    edge_file : str
        Path to edge list file
    node_file : str, optional
        Path to node attribute file
    delimiter : str
        Delimiter used in files
    **kwargs
        Additional arguments to pass to load_network_from_dataframe
        
    Returns:
    --------
    nx.Graph
        NetworkX graph representation of the network
    """
    edge_df = pd.read_csv(edge_file, delimiter=delimiter)
    
    node_df = None
    if node_file and os.path.exists(node_file):
        node_df = pd.read_csv(node_file, delimiter=delimiter)
    
    return load_network_from_dataframe(edge_df, node_df, **kwargs)

def load_network_from_transnet(transnet_obj: Any) -> nx.Graph:
    """
    Load a network from a Transnet object into a NetworkX graph.
    
    Parameters:
    -----------
    transnet_obj : Transnet
        Transnet object
        
    Returns:
    --------
    nx.Graph
        NetworkX graph representation of the network
    """
    # Generate the interaction DataFrame
    interactions = transnet_obj.generate_interaction_df()
    
    # Create an empty graph
    G = nx.Graph()
    
    # Add edges with attributes
    for _, row in interactions.iterrows():
        source = row["SourceNode"]
        target = row["TargetNode"]
        weight = row["Weight"]
        
        # Add edge with attributes
        G.add_edge(source, target, weight=weight, type=row.get("Comment", ""))
        
        # Add node layer information
        G.nodes[source]['layer'] = row["SourceLayer"]
        G.nodes[target]['layer'] = row["TargetLayer"]
    
    return G

def integrate_experimental_data_from_files(
    G: nx.Graph,
    data_files: Dict[str, Tuple[str, str, Optional[str]]],
    p_value_threshold: float = 0.05
) -> nx.Graph:
    """
    Integrate experimental data from multiple files into a network.
    
    Parameters:
    -----------
    G : nx.Graph
        NetworkX graph
    data_files : Dict[str, Tuple[str, str, Optional[str]]]
        Dictionary mapping data types to tuples of 
        (file_path, value_column, p_value_column)
    p_value_threshold : float
        P-value threshold for significance
        
    Returns:
    --------
    nx.Graph
        NetworkX graph with integrated data
    """
    from transnet.analysis.network_analysis import integrate_experimental_data
    
    # Create a copy of the graph
    G_copy = G.copy()
    
    # Integrate each data file
    for data_type, (file_path, value_col, p_value_col) in data_files.items():
        df = pd.read_csv(file_path)
        
        # Determine ID column by checking first column
        id_col = df.columns[0]
        
        # Integrate the data
        G_copy = integrate_experimental_data(
            G_copy, 
            df, 
            id_column=id_col,
            value_column=value_col,
            p_value_column=p_value_col,
            p_value_threshold=p_value_threshold
        )
        
        # Tag nodes with the data type
        for node in G_copy.nodes():
            if G_copy.nodes[node].get('experimental_value') is not None:
                G_copy.nodes[node][f'has_{data_type}_data'] = True
    
    return G_copy
    