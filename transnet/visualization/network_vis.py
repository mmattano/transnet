"""
Network visualization methods for TransNet multi-layer networks.

This module provides a range of visualization functions for TransNet networks, from
basic static plots to advanced interactive visualizations.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx
import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Union, Any
import logging
import io
import base64
from matplotlib.colors import LinearSegmentedColormap, Normalize

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Default color schemes for different network layers
DEFAULT_LAYER_COLORS = {
    'Transcriptome': '#1f77b4',   # Blue
    'Proteome': '#ff7f0e',        # Orange
    'Metabolome': '#2ca02c',      # Green
    'Reactions': '#d62728',       # Red
    'Pathways': '#9467bd',        # Purple
    'Unknown': '#7f7f7f'          # Gray
}

# Default node shapes for different network layers
DEFAULT_LAYER_SHAPES = {
    'Transcriptome': 'o',     # Circle
    'Proteome': 's',          # Square
    'Metabolome': 'h',        # Hexagon
    'Reactions': 'd',         # Diamond
    'Pathways': 'p',          # Pentagon
    'Unknown': 'o'            # Circle
}

def plot_network_basic(
    G: nx.Graph,
    node_colors: Optional[Dict[str, str]] = None,
    node_sizes: Optional[Dict[str, float]] = None,
    edge_colors: Optional[Dict[Tuple[str, str], str]] = None,
    title: str = "Network Visualization",
    layout: str = "spring",
    figsize: Tuple[int, int] = (12, 10),
    show_labels: bool = True,
    label_font_size: int = 8,
    save_path: Optional[str] = None,
    dpi: int = 300,
    node_alpha: float = 0.7,
    edge_alpha: float = 0.5,
    edge_width: float = 1.0,
    node_border_width: float = 0.5
) -> plt.Figure:
    """
    Create a basic static network visualization.
    
    Parameters:
    -----------
    G : nx.Graph
        NetworkX graph to visualize
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
    show_labels : bool
        Whether to show node labels
    label_font_size : int
        Font size for node labels
    save_path : str, optional
        Path to save the visualization
    dpi : int
        DPI for saved image
    node_alpha : float
        Transparency of nodes
    edge_alpha : float
        Transparency of edges
    edge_width : float
        Width of edges
    node_border_width : float
        Width of node borders
        
    Returns:
    --------
    plt.Figure
        Matplotlib figure
    """
    plt.figure(figsize=figsize)
    
    # Choose layout
    if layout == 'spring':
        pos = nx.spring_layout(G, k=0.3, iterations=50)
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
        pos = nx.spring_layout(G, k=0.3, iterations=50)
    
    # Set default node colors (by layer)
    if node_colors is None:
        node_colors = {}
        for node in G.nodes():
            layer = G.nodes[node].get('layer', 'Unknown')
            node_colors[node] = DEFAULT_LAYER_COLORS.get(layer, DEFAULT_LAYER_COLORS['Unknown'])
    
    # Set default node sizes (by degree)
    if node_sizes is None:
        degrees = dict(G.degree())
        max_degree = max(degrees.values()) if degrees else 1
        node_sizes = {node: 50 + 150 * (degrees.get(node, 0) / max_degree) for node in G.nodes()}
    
    # Set default edge colors
    if edge_colors is None:
        edge_colors = {edge: '#888888' for edge in G.edges()}
    
    # Draw the network
    nx.draw_networkx_nodes(
        G, pos, 
        node_color=[node_colors.get(node, DEFAULT_LAYER_COLORS['Unknown']) for node in G.nodes()],
        node_size=[node_sizes.get(node, 50) for node in G.nodes()],
        alpha=node_alpha,
        linewidths=node_border_width,
        edgecolors='black'
    )
    
    nx.draw_networkx_edges(
        G, pos, 
        edge_color=[edge_colors.get(edge, '#888888') for edge in G.edges()],
        width=edge_width,
        alpha=edge_alpha
    )
    
    if show_labels:
        # Truncate long node labels
        labels = {}
        for node in G.nodes():
            if len(str(node)) > 15:
                labels[node] = str(node)[:12] + '...'
            else:
                labels[node] = str(node)
        
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=label_font_size)
    
    plt.title(title)
    plt.axis('off')
    
    # Create legend for layers
    legend_elements = []
    layers_seen = set()
    
    for node in G.nodes():
        layer = G.nodes[node].get('layer', 'Unknown')
        if layer not in layers_seen:
            layers_seen.add(layer)
            color = node_colors.get(node, DEFAULT_LAYER_COLORS.get(layer, DEFAULT_LAYER_COLORS['Unknown']))
            legend_elements.append(
                mpatches.Patch(color=color, alpha=0.7, label=layer)
            )
    
    plt.legend(handles=legend_elements, loc='upper right')
    
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
    
    return plt.gcf()

def plot_multilayer_network(
    G: nx.Graph,
    layer_vertical_spacing: float = 2.0,
    layer_order: Optional[List[str]] = None,
    node_color_by: str = 'layer',
    node_size_by: str = 'degree',
    edge_color_by: str = 'type',
    title: str = "Multi-layer Network Visualization",
    figsize: Tuple[int, int] = (14, 10),
    show_labels: bool = True,
    label_font_size: int = 8,
    save_path: Optional[str] = None,
    dpi: int = 300,
    layer_colors: Optional[Dict[str, str]] = None,
    experimental_attr: Optional[str] = None,
    color_map: str = 'coolwarm',
    node_alpha: float = 0.7,
    edge_alpha: float = 0.4
) -> plt.Figure:
    """
    Create a visualization of a multi-layer network with layers positioned vertically.
    
    Parameters:
    -----------
    G : nx.Graph
        NetworkX graph to visualize
    layer_vertical_spacing : float
        Vertical spacing between layers
    layer_order : List[str], optional
        Order of layers from bottom to top
    node_color_by : str
        Attribute to color nodes by ('layer', 'experimental', or node attribute name)
    node_size_by : str
        Attribute to size nodes by ('degree', 'betweenness', or node attribute name)
    edge_color_by : str
        Attribute to color edges by ('type', 'weight', or edge attribute name)
    title : str
        Plot title
    figsize : Tuple[int, int]
        Figure size
    show_labels : bool
        Whether to show node labels
    label_font_size : int
        Font size for node labels
    save_path : str, optional
        Path to save the visualization
    dpi : int
        DPI for saved image
    layer_colors : Dict[str, str], optional
        Dictionary mapping layer names to colors
    experimental_attr : str, optional
        Node attribute containing experimental values (e.g., 'fc' for fold change)
    color_map : str
        Matplotlib colormap name for experimental data
    node_alpha : float
        Transparency of nodes
    edge_alpha : float
        Transparency of edges
        
    Returns:
    --------
    plt.Figure
        Matplotlib figure
    """
    plt.figure(figsize=figsize)
    
    # Get all layers in the network
    all_layers = set()
    for node in G.nodes():
        layer = G.nodes[node].get('layer', 'Unknown')
        all_layers.add(layer)
    
    # Determine layer order if not provided
    if layer_order is None:
        # Default order: Transcriptome -> Proteome -> Reactions -> Metabolome
        default_order = ['Transcriptome', 'Proteome', 'Reactions', 'Metabolome', 'Pathways']
        layer_order = [layer for layer in default_order if layer in all_layers]
        # Add any layers not in the default order
        for layer in all_layers:
            if layer not in layer_order:
                layer_order.append(layer)
    
    # Use default layer colors if not provided
    if layer_colors is None:
        layer_colors = DEFAULT_LAYER_COLORS
    
    # Create a position dictionary with each layer on a different y-level
    pos = {}
    layer_positions = {}
    
    for i, layer in enumerate(layer_order):
        # Get nodes in this layer
        layer_nodes = [node for node in G.nodes() if G.nodes[node].get('layer', 'Unknown') == layer]
        
        # Skip empty layers
        if not layer_nodes:
            continue
        
        # Create a position for nodes in this layer
        y_pos = i * layer_vertical_spacing
        layer_positions[layer] = y_pos
        
        # Position nodes in a horizontal line with some randomness
        for j, node in enumerate(layer_nodes):
            x_pos = j - len(layer_nodes) / 2  # Center horizontally
            x_pos += np.random.uniform(-0.2, 0.2)  # Add some randomness
            pos[node] = np.array([x_pos, y_pos])
    
    # Spread nodes horizontally within each layer using spring layout
    for layer, y_pos in layer_positions.items():
        layer_nodes = [node for node in G.nodes() if G.nodes[node].get('layer', 'Unknown') == layer]
        
        # Skip empty layers
        if not layer_nodes:
            continue
            
        # Create a subgraph for this layer
        layer_subgraph = G.subgraph(layer_nodes)
        
        # Get the positions for nodes in the layer
        try:
            layer_pos = nx.spring_layout(layer_subgraph, k=0.3, iterations=50, dim=1)
            
            # Update the positions (x only, keep y the same)
            for node, position in layer_pos.items():
                pos[node] = np.array([position[0] * 5, y_pos])  # Scale x for better spacing
        except:
            # If spring layout fails (e.g., for disconnected graphs), use a simple line layout
            for j, node in enumerate(layer_nodes):
                x_pos = (j - len(layer_nodes) / 2) * 2  # Scale for better spacing
                pos[node] = np.array([x_pos, y_pos])
    
    # Prepare node colors
    node_colors = []
    
    if node_color_by == 'layer':
        # Color by layer
        for node in G.nodes():
            layer = G.nodes[node].get('layer', 'Unknown')
            node_colors.append(layer_colors.get(layer, DEFAULT_LAYER_COLORS['Unknown']))
    
    elif node_color_by == 'experimental' and experimental_attr is not None:
        # Color by experimental value
        values = []
        for node in G.nodes():
            value = G.nodes[node].get(experimental_attr)
            if value is not None:
                values.append(float(value))
        
        if values:
            vmin, vmax = min(values), max(values)
            
            # Create colormap
            cmap = plt.get_cmap(color_map)
            norm = Normalize(vmin=vmin, vmax=vmax)
            
            for node in G.nodes():
                value = G.nodes[node].get(experimental_attr)
                if value is not None:
                    node_colors.append(cmap(norm(float(value))))
                else:
                    node_colors.append('#7f7f7f')  # Gray for nodes without data
            
            # Add a colorbar
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=plt.gca(), orientation='vertical', pad=0.05)
            cbar.set_label(experimental_attr)
        else:
            # Fallback to layer colors if no experimental values
            for node in G.nodes():
                layer = G.nodes[node].get('layer', 'Unknown')
                node_colors.append(layer_colors.get(layer, DEFAULT_LAYER_COLORS['Unknown']))
    
    else:
        # Color by some other attribute
        for node in G.nodes():
            attr_value = G.nodes[node].get(node_color_by)
            if attr_value is not None:
                # If the attribute is a string, use it as a color
                if isinstance(attr_value, str) and attr_value.startswith('#'):
                    node_colors.append(attr_value)
                else:
                    # Otherwise, color by layer
                    layer = G.nodes[node].get('layer', 'Unknown')
                    node_colors.append(layer_colors.get(layer, DEFAULT_LAYER_COLORS['Unknown']))
            else:
                # Fallback to layer color
                layer = G.nodes[node].get('layer', 'Unknown')
                node_colors.append(layer_colors.get(layer, DEFAULT_LAYER_COLORS['Unknown']))
    
    # Prepare node sizes
    node_sizes = []
    
    if node_size_by == 'degree':
        # Size by degree
        degrees = dict(G.degree())
        max_degree = max(degrees.values()) if degrees else 1
        for node in G.nodes():
            node_sizes.append(50 + 150 * (degrees.get(node, 0) / max_degree))
    else:
        # Size by some other attribute
        attr_values = []
        for node in G.nodes():
            attr_value = G.nodes[node].get(node_size_by)
            if attr_value is not None and isinstance(attr_value, (int, float)):
                attr_values.append(float(attr_value))
        
        if attr_values:
            max_attr = max(attr_values)
            for node in G.nodes():
                attr_value = G.nodes[node].get(node_size_by)
                if attr_value is not None and isinstance(attr_value, (int, float)):
                    node_sizes.append(50 + 150 * (float(attr_value) / max_attr))
                else:
                    node_sizes.append(50)  # Default size
        else:
            # Default sizes
            node_sizes = [50] * len(G.nodes())
    
    # Prepare edge colors
    edge_colors = []
    
    if edge_color_by == 'type':
        # Color by edge type
        edge_types = {}
        for u, v, data in G.edges(data=True):
            edge_type = data.get('type', 'default')
            edge_types[edge_type] = edge_types.get(edge_type, 0) + 1
        
        # Generate colors for each edge type
        type_colors = {}
        for i, edge_type in enumerate(edge_types.keys()):
            type_colors[edge_type] = plt.cm.tab10(i % 10)
        
        for u, v, data in G.edges(data=True):
            edge_type = data.get('type', 'default')
            edge_colors.append(type_colors.get(edge_type, '#888888'))
    
    elif edge_color_by == 'weight':
        # Color by edge weight
        weights = []
        for u, v, data in G.edges(data=True):
            weight = data.get('weight', 1.0)
            weights.append(float(weight))
        
        if weights:
            min_weight, max_weight = min(weights), max(weights)
            
            # Create colormap
            cmap = plt.cm.viridis
            norm = Normalize(vmin=min_weight, vmax=max_weight)
            
            for u, v, data in G.edges(data=True):
                weight = data.get('weight', 1.0)
                edge_colors.append(cmap(norm(float(weight))))
        else:
            # Default color
            edge_colors = ['#888888'] * len(G.edges())
    
    else:
        # Color by some other attribute or default
        edge_colors = ['#888888'] * len(G.edges())
    
    # Draw the network
    nx.draw_networkx_nodes(
        G, pos, 
        node_color=node_colors,
        node_size=node_sizes,
        alpha=node_alpha,
        linewidths=0.5,
        edgecolors='black'
    )
    
    nx.draw_networkx_edges(
        G, pos, 
        edge_color=edge_colors,
        width=1.0,
        alpha=edge_alpha,
        connectionstyle='arc3,rad=0.1'  # Curved edges
    )
    
    if show_labels:
        # Truncate long node labels and only show for bigger nodes
        labels = {}
        for node in G.nodes():
            if len(str(node)) > 15:
                labels[node] = str(node)[:12] + '...'
            else:
                labels[node] = str(node)
        
        node_idx = {node: i for i, node in enumerate(G.nodes())}
        show_label_indices = [i for i, size in enumerate(node_sizes) if size > 100]
        labels_to_show = {node: labels[node] for node in G.nodes() if node_idx[node] in show_label_indices}
        
        nx.draw_networkx_labels(G, pos, labels=labels_to_show, font_size=label_font_size)
    
    # Add layer labels
    for layer, y_pos in layer_positions.items():
        plt.text(
            -8, y_pos, 
            layer, 
            fontsize=12, 
            ha='right', 
            va='center',
            bbox=dict(facecolor=layer_colors.get(layer, DEFAULT_LAYER_COLORS['Unknown']), alpha=0.3, boxstyle='round,pad=0.5')
        )
    
    plt.title(title)
    plt.axis('off')
    
    # Create legend for layers
    legend_elements = []
    for layer in layer_positions.keys():
        color = layer_colors.get(layer, DEFAULT_LAYER_COLORS['Unknown'])
        legend_elements.append(
            mpatches.Patch(color=color, alpha=0.7, label=layer)
        )
    
    plt.legend(handles=legend_elements, loc='upper right')
    
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
    
    return plt.gcf()

def plot_subnetwork_heatmap(
    G: nx.Graph,
    node_list: List[str],
    attribute: str,
    figsize: Tuple[int, int] = (12, 10),
    title: str = "Node Attribute Heatmap",
    color_map: str = 'viridis',
    save_path: Optional[str] = None,
    dpi: int = 300,
    label_fontsize: int = 10,
    cluster: bool = True
) -> plt.Figure:
    """
    Create a heatmap visualization of node attributes for a subset of nodes.
    
    Parameters:
    -----------
    G : nx.Graph
        NetworkX graph
    node_list : List[str]
        List of nodes to include in the heatmap
    attribute : str
        Node attribute to visualize
    figsize : Tuple[int, int]
        Figure size
    title : str
        Plot title
    color_map : str
        Matplotlib colormap name
    save_path : str, optional
        Path to save the visualization
    dpi : int
        DPI for saved image
    label_fontsize : int
        Font size for node labels
    cluster : bool
        Whether to cluster nodes by similarity
        
    Returns:
    --------
    plt.Figure
        Matplotlib figure
    """
    # Filter nodes that exist in the graph and have the attribute
    valid_nodes = [node for node in node_list if node in G.nodes() and attribute in G.nodes[node]]
    
    if not valid_nodes:
        logger.warning(f"No nodes with attribute '{attribute}' found in the given node list")
        return plt.figure()
    
    # Create a matrix of attribute values
    attr_values = [G.nodes[node][attribute] for node in valid_nodes]
    
    if not all(isinstance(val, (int, float)) for val in attr_values):
        logger.warning(f"Attribute '{attribute}' contains non-numeric values")
        return plt.figure()
    
    # Create a DataFrame for the heatmap
    df = pd.DataFrame({
        'Node': valid_nodes,
        'Layer': [G.nodes[node].get('layer', 'Unknown') for node in valid_nodes],
        'Value': attr_values
    })
    
    # Sort by layer and value
    df = df.sort_values(['Layer', 'Value'])
    
    # Plot the heatmap
    plt.figure(figsize=figsize)
    
    # Create heatmap data
    heatmap_data = df['Value'].values.reshape(-1, 1)
    
    if cluster and len(valid_nodes) > 2:
        # Cluster nodes by attribute similarity
        from scipy.cluster.hierarchy import linkage, dendrogram
        try:
            # Compute linkage matrix
            Z = linkage(heatmap_data, 'ward')
            
            # Compute dendrogram
            dendro = dendrogram(Z, no_plot=True)
            
            # Reorder based on dendrogram
            df = df.iloc[dendro['leaves']]
            heatmap_data = df['Value'].values.reshape(-1, 1)
        except Exception as e:
            logger.warning(f"Clustering failed: {e}")
    
    # Plot heatmap
    plt.imshow(
        heatmap_data,
        aspect='auto',
        cmap=color_map,
        interpolation='nearest'
    )
    
    # Add colorbar
    plt.colorbar(label=attribute)
    
    # Add node labels
    plt.yticks(
        range(len(df)), 
        df['Node'].values,
        fontsize=label_fontsize
    )
    
    # Add layer markers
    layer_colors = {}
    for i, (_, row) in enumerate(df.iterrows()):
        layer = row['Layer']
        if layer not in layer_colors:
            layer_colors[layer] = DEFAULT_LAYER_COLORS.get(layer, DEFAULT_LAYER_COLORS['Unknown'])
        
        plt.gca().add_patch(
            plt.Rectangle(
                (-1.5, i - 0.5),
                1,
                1,
                color=layer_colors[layer],
                alpha=0.7
            )
        )
    
    # Add legend for layers
    legend_elements = []
    for layer, color in layer_colors.items():
        legend_elements.append(
            mpatches.Patch(color=color, alpha=0.7, label=layer)
        )
    
    plt.legend(
        handles=legend_elements,
        loc='upper left',
        bbox_to_anchor=(1.05, 1),
        title="Layers"
    )
    
    plt.title(title)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
    
    return plt.gcf()

def plot_network_metrics(
    G: nx.Graph,
    metrics: List[str] = ['degree', 'betweenness', 'closeness', 'eigenvector'],
    top_n: int = 10,
    figsize: Tuple[int, int] = (14, 10),
    title: str = "Network Metrics",
    save_path: Optional[str] = None,
    dpi: int = 300,
    color_by_layer: bool = True
) -> plt.Figure:
    """
    Create a visualization of various network metrics.
    
    Parameters:
    -----------
    G : nx.Graph
        NetworkX graph
    metrics : List[str]
        List of metrics to visualize ('degree', 'betweenness', 'closeness', 'eigenvector')
    top_n : int
        Number of top nodes to show for each metric
    figsize : Tuple[int, int]
        Figure size
    title : str
        Plot title
    save_path : str, optional
        Path to save the visualization
    dpi : int
        DPI for saved image
    color_by_layer : bool
        Whether to color bars by node layer
        
    Returns:
    --------
    plt.Figure
        Matplotlib figure
    """
    # Calculate metrics
    metric_functions = {
        'degree': nx.degree_centrality,
        'betweenness': nx.betweenness_centrality,
        'closeness': nx.closeness_centrality,
        'eigenvector': nx.eigenvector_centrality_numpy
    }
    
    metric_values = {}
    for metric in metrics:
        if metric in metric_functions:
            try:
                metric_values[metric] = metric_functions[metric](G)
            except Exception as e:
                logger.warning(f"Error calculating {metric} centrality: {e}")
                metric_values[metric] = {}
        else:
            logger.warning(f"Unknown metric: {metric}")
    
    if not metric_values:
        logger.warning("No metrics could be calculated")
        return plt.figure()
    
    # Create subplot layout
    num_metrics = len(metric_values)
    fig, axes = plt.subplots(num_metrics, 1, figsize=figsize)
    
    if num_metrics == 1:
        axes = [axes]
    
    # Plot each metric
    for i, (metric, values) in enumerate(metric_values.items()):
        ax = axes[i]
        
        # Get top nodes
        top_nodes = sorted(values.items(), key=lambda x: x[1], reverse=True)[:top_n]
        
        if not top_nodes:
            ax.text(0.5, 0.5, f"No data for {metric}", ha='center', va='center')
            continue
        
        nodes, values = zip(*top_nodes)
        
        # Create bars with colors by layer if requested
        bar_colors = []
        if color_by_layer:
            for node in nodes:
                layer = G.nodes[node].get('layer', 'Unknown')
                bar_colors.append(DEFAULT_LAYER_COLORS.get(layer, DEFAULT_LAYER_COLORS['Unknown']))
        else:
            bar_colors = DEFAULT_LAYER_COLORS['Unknown']
        
        # Plot bars
        bars = ax.barh(range(len(nodes)), values, color=bar_colors, alpha=0.7)
        
        # Add node labels
        node_labels = []
        for node in nodes:
            if len(str(node)) > 20:
                node_labels.append(str(node)[:17] + '...')
            else:
                node_labels.append(str(node))
        
        ax.set_yticks(range(len(nodes)))
        ax.set_yticklabels(node_labels)
        
        # Add layer markers if coloring by layer
        if color_by_layer:
            # Create legend for layers
            legend_elements = []
            layers_seen = set()
            
            for node in nodes:
                layer = G.nodes[node].get('layer', 'Unknown')
                if layer not in layers_seen:
                    layers_seen.add(layer)
                    color = DEFAULT_LAYER_COLORS.get(layer, DEFAULT_LAYER_COLORS['Unknown'])
                    legend_elements.append(
                        mpatches.Patch(color=color, alpha=0.7, label=layer)
                    )
            
            ax.legend(handles=legend_elements, loc='lower right')
        
        ax.set_title(f"{metric.capitalize()} Centrality - Top {len(nodes)} Nodes")
        ax.set_xlabel("Centrality Value")
        ax.grid(axis='x', linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    fig.suptitle(title, fontsize=16, y=1.02)
    
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
    
    return fig

def plot_community_network(
    G: nx.Graph,
    communities: Dict[str, int],
    figsize: Tuple[int, int] = (12, 10),
    title: str = "Network Communities",
    layout: str = "spring",
    save_path: Optional[str] = None,
    dpi: int = 300,
    show_labels: bool = True,
    label_font_size: int = 8,
    node_alpha: float = 0.7,
    edge_alpha: float = 0.3,
    show_legend: bool = True,
    layer_markers: bool = True
) -> plt.Figure:
    """
    Create a visualization of a network with nodes colored by community.
    
    Parameters:
    -----------
    G : nx.Graph
        NetworkX graph
    communities : Dict[str, int]
        Dictionary mapping nodes to community IDs
    figsize : Tuple[int, int]
        Figure size
    title : str
        Plot title
    layout : str
        Layout algorithm ('spring', 'circular', 'kamada_kawai', 'spectral', 'random', 'shell')
    save_path : str, optional
        Path to save the visualization
    dpi : int
        DPI for saved image
    show_labels : bool
        Whether to show node labels
    label_font_size : int
        Font size for node labels
    node_alpha : float
        Transparency of nodes
    edge_alpha : float
        Transparency of edges
    show_legend : bool
        Whether to show the community legend
    layer_markers : bool
        Whether to add layer markers to nodes
        
    Returns:
    --------
    plt.Figure
        Matplotlib figure
    """
    plt.figure(figsize=figsize)
    
    # Choose layout
    if layout == 'spring':
        pos = nx.spring_layout(G, k=0.3, iterations=50)
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
        pos = nx.spring_layout(G, k=0.3, iterations=50)
    
    # Get all community IDs
    community_ids = set(communities.values())
    
    # Create a colormap for communities
    cmap = plt.cm.tab20
    
    # Prepare node colors and sizes
    node_colors = []
    node_sizes = []
    node_shapes = []
    
    # Degree-based sizes
    degrees = dict(G.degree())
    max_degree = max(degrees.values()) if degrees else 1
    
    for node in G.nodes():
        # Color by community
        community_id = communities.get(node, -1)
        if community_id >= 0:
            color_idx = community_id % 20  # Use tab20 colormap (20 colors)
            node_colors.append(cmap(color_idx))
        else:
            node_colors.append('#7f7f7f')  # Gray for nodes without community
        
        # Size by degree
        node_sizes.append(50 + 150 * (degrees.get(node, 0) / max_degree))
        
        # Shape by layer if requested
        if layer_markers:
            layer = G.nodes[node].get('layer', 'Unknown')
            node_shapes.append(DEFAULT_LAYER_SHAPES.get(layer, DEFAULT_LAYER_SHAPES['Unknown']))
        else:
            node_shapes.append('o')  # Default shape
    
    # Draw nodes by shape group if using layer markers
    if layer_markers:
        # Group nodes by shape
        shape_groups = {}
        for i, node in enumerate(G.nodes()):
            shape = node_shapes[i]
            if shape not in shape_groups:
                shape_groups[shape] = []
            shape_groups[shape].append(i)
        
        # Draw each shape group separately
        for shape, indices in shape_groups.items():
            node_list = [list(G.nodes())[i] for i in indices]
            node_color_list = [node_colors[i] for i in indices]
            node_size_list = [node_sizes[i] for i in indices]
            
            nx.draw_networkx_nodes(
                G, pos, 
                nodelist=node_list,
                node_color=node_color_list,
                node_size=node_size_list,
                alpha=node_alpha,
                linewidths=0.5,
                edgecolors='black',
                node_shape=shape
            )
    else:
        # Draw all nodes together
        nx.draw_networkx_nodes(
            G, pos, 
            node_color=node_colors,
            node_size=node_sizes,
            alpha=node_alpha,
            linewidths=0.5,
            edgecolors='black'
        )
    
    # Draw edges
    nx.draw_networkx_edges(
        G, pos, 
        width=1.0,
        alpha=edge_alpha
    )
    
    if show_labels:
        # Truncate long node labels and only show for bigger nodes
        labels = {}
        for node in G.nodes():
            if len(str(node)) > 15:
                labels[node] = str(node)[:12] + '...'
            else:
                labels[node] = str(node)
        
        # Only show labels for larger nodes to avoid clutter
        node_idx = {node: i for i, node in enumerate(G.nodes())}
        show_label_indices = [i for i, size in enumerate(node_sizes) if size > 100]
        labels_to_show = {node: labels[node] for node in G.nodes() if node_idx[node] in show_label_indices}
        
        nx.draw_networkx_labels(G, pos, labels=labels_to_show, font_size=label_font_size)
    
    plt.title(title)
    plt.axis('off')
    
    # Create legends
    if show_legend:
        # Count nodes in each community
        community_counts = {}
        for _, comm_id in communities.items():
            community_counts[comm_id] = community_counts.get(comm_id, 0) + 1
        
        # Get top communities by size
        top_communities = sorted(community_counts.items(), key=lambda x: x[1], reverse=True)[:10]
        
        # Create legend for communities
        legend_elements = []
        for comm_id, count in top_communities:
            color_idx = comm_id % 20
            legend_elements.append(
                mpatches.Patch(
                    color=cmap(color_idx), 
                    alpha=0.7, 
                    label=f"Community {comm_id} ({count} nodes)"
                )
            )
        
        plt.legend(
            handles=legend_elements,
            loc='upper left',
            bbox_to_anchor=(1.05, 1),
            title="Top Communities"
        )
    
    # Add layer legend if using layer markers
    if layer_markers:
        # Get all layers in the network
        all_layers = set()
        for node in G.nodes():
            layer = G.nodes[node].get('layer', 'Unknown')
            all_layers.add(layer)
        
        # Create legend for layers
        layer_elements = []
        for layer in all_layers:
            layer_elements.append(
                plt.Line2D(
                    [0], [0],
                    marker=DEFAULT_LAYER_SHAPES.get(layer, DEFAULT_LAYER_SHAPES['Unknown']),
                    color='w',
                    markerfacecolor='black',
                    markersize=10,
                    label=layer
                )
            )
        
        if layer_elements:
            second_legend = plt.legend(
                handles=layer_elements,
                loc='lower left',
                bbox_to_anchor=(1.05, 0),
                title="Layers"
            )
            
            # Add the second legend manually if showing both legends
            if show_legend:
                plt.gca().add_artist(second_legend)
    
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
    
    return plt.gcf()

def create_interactive_html_network(
    G: nx.Graph,
    title: str = "Interactive Network Visualization",
    width: int = 1000,
    height: int = 800,
    save_path: Optional[str] = None,
    color_by: str = 'layer',
    size_by: str = 'degree',
    experimental_attr: Optional[str] = None,
    node_hover_info: List[str] = None,
    show_legend: bool = True
) -> str:
    """
    Create an interactive HTML network visualization using D3.js.
    
    Parameters:
    -----------
    G : nx.Graph
        NetworkX graph
    title : str
        Visualization title
    width : int
        Width of the visualization in pixels
    height : int
        Height of the visualization in pixels
    save_path : str, optional
        Path to save the HTML file
    color_by : str
        Node attribute to color by ('layer', 'experimental', or node attribute name)
    size_by : str
        Node attribute to size by ('degree', 'betweenness', or node attribute name)
    experimental_attr : str, optional
        Node attribute containing experimental values (e.g., 'fc' for fold change)
    node_hover_info : List[str], optional
        List of node attributes to show on hover
    show_legend : bool
        Whether to show the legend
        
    Returns:
    --------
    str
        HTML content for the visualization
    """
    # Convert graph to JSON format
    try:
        import json
    except ImportError:
        logger.error("json module not available")
        return "<p>Error: json module not available</p>"
    
    # Prepare node data
    nodes_data = []
    for node in G.nodes():
        node_data = {
            "id": str(node),
            "label": str(node)
        }
        
        # Add all node attributes
        for attr, value in G.nodes[node].items():
            if isinstance(value, (int, float, str, bool)):
                node_data[attr] = value
            elif value is None:
                node_data[attr] = None
            else:
                # Convert complex objects to string
                node_data[attr] = str(value)
        
        # Set node color
        if color_by == 'layer':
            layer = G.nodes[node].get('layer', 'Unknown')
            node_data["color"] = DEFAULT_LAYER_COLORS.get(layer, DEFAULT_LAYER_COLORS['Unknown'])
        elif color_by == 'experimental' and experimental_attr is not None:
            value = G.nodes[node].get(experimental_attr)
            if value is not None:
                # Will be colored in the HTML/JS
                node_data["value"] = float(value)
            else:
                node_data["color"] = '#7f7f7f'  # Gray for nodes without data
        else:
            # Try to use the specified attribute
            attr_value = G.nodes[node].get(color_by)
            if attr_value is not None:
                if isinstance(attr_value, str) and attr_value.startswith('#'):
                    node_data["color"] = attr_value
                else:
                    # Default to layer color
                    layer = G.nodes[node].get('layer', 'Unknown')
                    node_data["color"] = DEFAULT_LAYER_COLORS.get(layer, DEFAULT_LAYER_COLORS['Unknown'])
            else:
                # Default to layer color
                layer = G.nodes[node].get('layer', 'Unknown')
                node_data["color"] = DEFAULT_LAYER_COLORS.get(layer, DEFAULT_LAYER_COLORS['Unknown'])
        
        # Set node size
        if size_by == 'degree':
            node_data["size"] = G.degree(node) * 2 + 5
        else:
            # Try to use the specified attribute
            attr_value = G.nodes[node].get(size_by)
            if attr_value is not None and isinstance(attr_value, (int, float)):
                node_data["size"] = float(attr_value) * 5 + 5
            else:
                node_data["size"] = 10  # Default size
        
        nodes_data.append(node_data)
    
    # Prepare edge data
    edges_data = []
    for u, v, data in G.edges(data=True):
        edge_data = {
            "source": str(u),
            "target": str(v),
            "value": data.get('weight', 1.0)
        }
        
        # Add edge type if available
        if 'type' in data:
            edge_data["type"] = data['type']
        
        edges_data.append(edge_data)
    
    # Create JSON data
    graph_data = {
        "nodes": nodes_data,
        "links": edges_data
    }
    
    # Get all layers for legend
    layers = set()
    for node in G.nodes():
        layer = G.nodes[node].get('layer', 'Unknown')
        layers.add(layer)
    
    # Determine layer colors for legend
    layer_colors = {}
    for layer in layers:
        layer_colors[layer] = DEFAULT_LAYER_COLORS.get(layer, DEFAULT_LAYER_COLORS['Unknown'])
    
    # Create HTML content with proper escaping of JavaScript code
    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <title>{title}</title>
        <script src="https://d3js.org/d3.v7.min.js"></script>
        <style>
            body {{
                font-family: Arial, sans-serif;
                margin: 0;
                padding: 0;
                background-color: #f5f5f5;
            }}
            
            .container {{
                max-width: {width}px;
                margin: 0 auto;
                padding: 20px;
            }}
            
            h1 {{
                text-align: center;
                margin-bottom: 20px;
            }}
            
            #network {{
                width: 100%;
                height: {height}px;
                border: 1px solid #ddd;
                background-color: white;
                border-radius: 4px;
                box-shadow: 0 1px 3px rgba(0,0,0,0.1);
            }}
            
            .tooltip {{
                position: absolute;
                background-color: rgba(255, 255, 255, 0.9);
                padding: 10px;
                border-radius: 4px;
                box-shadow: 0 1px 3px rgba(0,0,0,0.2);
                pointer-events: none;
                max-width: 300px;
                z-index: 999;
            }}
            
            .tooltip p {{
                margin: 0;
                margin-bottom: 5px;
                font-size: 12px;
            }}
            
            .tooltip p.title {{
                font-weight: bold;
                font-size: 14px;
                margin-bottom: 8px;
                border-bottom: 1px solid #ddd;
                padding-bottom: 5px;
            }}
            
            .legend {{
                position: absolute;
                top: 20px;
                right: 20px;
                background-color: rgba(255, 255, 255, 0.9);
                padding: 10px;
                border-radius: 4px;
                box-shadow: 0 1px 3px rgba(0,0,0,0.2);
                max-width: 200px;
                font-size: 12px;
            }}
            
            .legend-title {{
                font-weight: bold;
                margin-bottom: 5px;
                border-bottom: 1px solid #ddd;
                padding-bottom: 3px;
            }}
            
            .legend-item {{
                display: flex;
                align-items: center;
                margin-bottom: 3px;
            }}
            
            .legend-color {{
                width: 12px;
                height: 12px;
                margin-right: 5px;
                border-radius: 3px;
            }}
            
            .controls {{
                margin-top: 10px;
                display: flex;
                justify-content: center;
            }}
            
            .controls button {{
                margin: 0 5px;
                padding: 5px 10px;
                background-color: #4CAF50;
                color: white;
                border: none;
                border-radius: 4px;
                cursor: pointer;
            }}
            
            .controls button:hover {{
                background-color: #45a049;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>{title}</h1>
            <div id="network"></div>
            <div class="controls">
                <button id="zoom-in">Zoom In</button>
                <button id="zoom-out">Zoom Out</button>
                <button id="reset">Reset</button>
            </div>
        </div>
        
        <script>
            // Network data
            const graphData = {json.dumps(graph_data)};
            
            // Layer colors
            const layerColors = {json.dumps(layer_colors)};
            
            // Create the network visualization
            const width = {width};
            const height = {height};
            
            // Create SVG
            const svg = d3.select("#network")
                .append("svg")
                .attr("width", "100%")
                .attr("height", "100%")
                .attr("viewBox", [0, 0, width, height]);
            
            // Create group for zoom/pan
            const g = svg.append("g");
            
            // Create tooltip
            const tooltip = d3.select("body")
                .append("div")
                .attr("class", "tooltip")
                .style("opacity", 0);
            
            // Create legend if enabled
            {f'const legend = svg.append("g").attr("class", "legend");' if show_legend else '// Legend disabled'}
            
            // Create color scale for experimental values if needed
            {f'const colorScale = d3.scaleSequential(d3.interpolateRdBu).domain([-1, 1]);' if color_by == 'experimental' and experimental_attr is not None else '// No experimental color scale needed'}
            
            // Create force simulation
            const simulation = d3.forceSimulation(graphData.nodes)
                .force("link", d3.forceLink(graphData.links).id(d => d.id).distance(100))
                .force("charge", d3.forceManyBody().strength(-200))
                .force("center", d3.forceCenter(width / 2, height / 2));
            
            // Create links
            const link = g.append("g")
                .selectAll("line")
                .data(graphData.links)
                .join("line")
                .attr("stroke", "#999")
                .attr("stroke-opacity", 0.6)
                .attr("stroke-width", d => Math.sqrt(d.value));
            
            // Create nodes
            const node = g.append("g")
                .selectAll("circle")
                .data(graphData.nodes)
                .join("circle")
                .attr("r", d => d.size)
                .attr("fill", d => {{
                    {f'if (d.hasOwnProperty("value")) {{ return colorScale(d.value); }}' if color_by == 'experimental' and experimental_attr is not None else ''}
                    return d.color;
                }})
                .attr("stroke", "#fff")
                .attr("stroke-width", 1.5)
                .call(drag(simulation));
            
            // Add node labels
            const label = g.append("g")
                .selectAll("text")
                .data(graphData.nodes)
                .join("text")
                .text(d => d.label)
                .attr("font-size", 10)
                .attr("dx", d => d.size + 5)
                .attr("dy", 4)
                .style("pointer-events", "none")
                .style("visibility", d => d.size > 15 ? "visible" : "hidden");
            
            // Set up node hover behavior
            node.on("mouseover", function(event, d) {{
                // Highlight node
                d3.select(this).attr("stroke", "#000").attr("stroke-width", 2);
                
                // Show tooltip
                tooltip.transition()
                    .duration(200)
                    .style("opacity", 0.9);
                
                let tooltipContent = `<p class="title">${{d.label}}</p>`;
                tooltipContent += `<p>Layer: ${{d.layer || 'Unknown'}}</p>`;
                
                {f'tooltipContent += `<p>{experimental_attr}: ${{d.value || "N/A"}}</p>`;' if experimental_attr else '// No experimental attribute to display'}
                
                // Add additional hover info if specified
                {f'const hoverAttrs = {json.dumps(node_hover_info) if node_hover_info else []};' if node_hover_info else '// No additional hover info specified'}
                {f'hoverAttrs.forEach(attr => {{' if node_hover_info else ''}
                    {f'if (d.hasOwnProperty(attr)) {{ tooltipContent += `<p>${{attr}}: ${{d[attr]}}</p>`; }}' if node_hover_info else ''}
                {f'}});' if node_hover_info else ''}
                
                tooltip.html(tooltipContent)
                    .style("left", (event.pageX + 10) + "px")
                    .style("top", (event.pageY - 10) + "px");
            }})
            .on("mouseout", function() {{
                // Remove highlight
                d3.select(this).attr("stroke", "#fff").attr("stroke-width", 1.5);
                
                // Hide tooltip
                tooltip.transition()
                    .duration(500)
                    .style("opacity", 0);
            }});
            
            // Create the legend if enabled
            if ({str(show_legend).lower()}) {{
                legend.append("text")
                    .attr("class", "legend-title")
                    .attr("x", width - 180)
                    .attr("y", 30)
                    .text("Network Layers");
                
                const layers = Object.keys(layerColors);
                
                layers.forEach((layer, i) => {{
                    const g = legend.append("g")
                        .attr("transform", `translate(${{width - 180}}, ${{45 + i * 20}})`);
                    
                    g.append("rect")
                        .attr("class", "legend-color")
                        .attr("width", 12)
                        .attr("height", 12)
                        .attr("fill", layerColors[layer]);
                    
                    g.append("text")
                        .attr("x", 20)
                        .attr("y", 10)
                        .attr("font-size", 12)
                        .text(layer);
                }});
                
                {f'// Add experimental color scale legend' if color_by == 'experimental' and experimental_attr is not None else '// No experimental legend needed'}
                {f'if ({str(color_by == "experimental" and experimental_attr is not None).lower()}) {{' if color_by == 'experimental' and experimental_attr is not None else ''}
                    {f'const legendY = 45 + layers.length * 20 + 20;' if color_by == 'experimental' and experimental_attr is not None else ''}
                    {f'legend.append("text")' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("class", "legend-title")' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("x", width - 180)' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("y", legendY)' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.text("{experimental_attr}");' if color_by == 'experimental' and experimental_attr is not None else ''}
                    
                    {f'// Create gradient' if color_by == 'experimental' and experimental_attr is not None else ''}
                    {f'const gradient = legend.append("linearGradient")' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("id", "gradient")' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("x1", "0%")' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("y1", "0%")' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("x2", "100%")' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("y2", "0%");' if color_by == 'experimental' and experimental_attr is not None else ''}
                    
                    {f'gradient.append("stop")' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("offset", "0%")' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("stop-color", colorScale(-1));' if color_by == 'experimental' and experimental_attr is not None else ''}
                    
                    {f'gradient.append("stop")' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("offset", "50%")' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("stop-color", colorScale(0));' if color_by == 'experimental' and experimental_attr is not None else ''}
                    
                    {f'gradient.append("stop")' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("offset", "100%")' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("stop-color", colorScale(1));' if color_by == 'experimental' and experimental_attr is not None else ''}
                    
                    {f'// Add gradient rectangle' if color_by == 'experimental' and experimental_attr is not None else ''}
                    {f'legend.append("rect")' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("x", width - 180)' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("y", legendY + 10)' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("width", 100)' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("height", 10)' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.style("fill", "url(#gradient)");' if color_by == 'experimental' and experimental_attr is not None else ''}
                    
                    {f'// Add min and max labels' if color_by == 'experimental' and experimental_attr is not None else ''}
                    {f'legend.append("text")' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("x", width - 180)' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("y", legendY + 35)' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("font-size", 10)' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.text("-1");' if color_by == 'experimental' and experimental_attr is not None else ''}
                    
                    {f'legend.append("text")' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("x", width - 130)' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("y", legendY + 35)' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("font-size", 10)' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.text("0");' if color_by == 'experimental' and experimental_attr is not None else ''}
                    
                    {f'legend.append("text")' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("x", width - 80)' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("y", legendY + 35)' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.attr("font-size", 10)' if color_by == 'experimental' and experimental_attr is not None else ''}
                        {f'.text("1");' if color_by == 'experimental' and experimental_attr is not None else ''}
                {f'}}' if color_by == 'experimental' and experimental_attr is not None else ''}
            }}
            
            // Set up zoom
            const zoom = d3.zoom()
                .scaleExtent([0.1, 4])
                .on("zoom", zoomed);
            
            svg.call(zoom);
            
            function zoomed(event) {{
                g.attr("transform", event.transform);
            }}
            
            // Drag function for nodes
            function drag(simulation) {{
                function dragstarted(event, d) {{
                    if (!event.active) simulation.alphaTarget(0.3).restart();
                    d.fx = d.x;
                    d.fy = d.y;
                }}
                
                function dragged(event, d) {{
                    d.fx = event.x;
                    d.fy = event.y;
                }}
                
                function dragended(event, d) {{
                    if (!event.active) simulation.alphaTarget(0);
                    d.fx = null;
                    d.fy = null;
                }}
                
                return d3.drag()
                    .on("start", dragstarted)
                    .on("drag", dragged)
                    .on("end", dragended);
            }}
            
            // Update simulation on each tick
            simulation.on("tick", () => {{
                link
                    .attr("x1", d => d.source.x)
                    .attr("y1", d => d.source.y)
                    .attr("x2", d => d.target.x)
                    .attr("y2", d => d.target.y);
                
                node
                    .attr("cx", d => d.x)
                    .attr("cy", d => d.y);
                
                label
                    .attr("x", d => d.x)
                    .attr("y", d => d.y);
            }});
            
            // Control buttons
            document.getElementById("zoom-in").addEventListener("click", () => {{
                svg.transition().duration(500).call(zoom.scaleBy, 1.5);
            }});
            
            document.getElementById("zoom-out").addEventListener("click", () => {{
                svg.transition().duration(500).call(zoom.scaleBy, 0.75);
            }});
            
            document.getElementById("reset").addEventListener("click", () => {{
                svg.transition().duration(500).call(
                    zoom.transform,
                    d3.zoomIdentity.translate(width / 2, height / 2).scale(1)
                );
            }});
        </script>
    </body>
    </html>
    """
    
    # Save to file if path provided
    if save_path:
        with open(save_path, 'w', encoding='utf-8') as f:
            f.write(html)
    
    return html

def export_network_visualization(G: nx.Graph, method: str = 'static', **kwargs) -> Any:
    """
    Export a network visualization using various methods.
    
    Parameters:
    -----------
    G : nx.Graph
        NetworkX graph to visualize
    method : str
        Visualization method ('static', 'multilayer', 'interactive', 'base64')
    **kwargs : dict
        Additional parameters for the specific visualization method
        
    Returns:
    --------
    Any
        Visualization output (Figure, HTML string, or base64 string)
    """
    if method == 'static':
        return plot_network_basic(G, **kwargs)
    
    elif method == 'multilayer':
        return plot_multilayer_network(G, **kwargs)
    
    elif method == 'interactive':
        return create_interactive_html_network(G, **kwargs)
    
    elif method == 'base64':
        # Create figure
        fig = plot_network_basic(G, **kwargs)
        
        # Convert to base64
        buffer = io.BytesIO()
        fig.savefig(buffer, format='png', dpi=kwargs.get('dpi', 300), bbox_inches='tight')
        buffer.seek(0)
        
        # Encode as base64
        image_base64 = base64.b64encode(buffer.getvalue()).decode('utf-8')
        plt.close(fig)
        
        return f"data:image/png;base64,{image_base64}"
    
    else:
        logger.error(f"Unknown visualization method: {method}")
        return None