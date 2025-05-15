"""
Example of using TransNet's visualization module to create network visualizations.
"""

import os
import sys
import logging
import matplotlib.pyplot as plt
import networkx as nx
from pathlib import Path

# Add the parent directory to the path
sys.path.append(str(Path(__file__).parent.parent))

from transnet.biology.transnet import Transnet
from transnet.biology.layers import Pathways, Proteome, Metabolome
from transnet.visualization.network_vis import (
    plot_network_basic,
    plot_multilayer_network,
    plot_community_network,
    plot_network_metrics,
    plot_subnetwork_heatmap,
    create_interactive_html_network
)
from transnet.analysis.network_analysis import (
    compute_network_statistics, 
    identify_hubs,
    compute_centrality_measures,
    detect_communities
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def build_sample_network():
    """Build a small sample network for visualization."""
    logger.info("Building sample network")
    
    # Create a pathways layer with just a few pathways
    pathways = Pathways()
    pathways.kegg_organism = "mmu"  # Mouse
    pathways.populate()
    # Keep only 2 pathways for demonstration
    pathways.pathways = pathways.pathways[:2]
    pathways.fill_pathways()
    
    # Create a simplified proteome layer
    proteome = Proteome()
    proteome.ncbi_organism = "10090"  # Mouse
    proteome.kegg_organism = "mmu"
    proteome.populate(uniprot=True)
    # Keep only a small subset of proteins
    proteome.proteins = proteome.proteins[:20]
    
    # Add some experimental data to proteins
    import numpy as np
    for i, protein in enumerate(proteome.proteins):
        # Add some synthetic experimental values
        protein.fc = np.sin(i * 0.5) * 2  # Values between -2 and 2
        protein.adj_p_value = 0.01 + 0.04 * np.random.random()  # Values between 0.01 and 0.05
    
    # Create a small metabolome layer
    metabolome = Metabolome()
    metabolome.populate()
    metabolome.metabolites = metabolome.metabolites[:15]
    
    # Add synthetic experimental data to metabolites
    for i, metabolite in enumerate(metabolome.metabolites):
        metabolite.fc = np.cos(i * 0.4) * 1.5  # Values between -1.5 and 1.5
        metabolite.adj_p_value = 0.005 + 0.04 * np.random.random()  # Values between 0.005 and 0.045
    
    # Create the integrated network
    transnet = Transnet(
        name="visualization_demo_network",
        pathways=pathways,
        proteome=proteome,
        metabolome=metabolome
    )
    
    logger.info("Generating network graph")
    G = transnet.generate_graph()
    
    # Add some additional edge types for visualization examples
    for u, v, attrs in list(G.edges(data=True)):
        if G.nodes[u].get('layer') == 'Proteome' and G.nodes[v].get('layer') == 'Proteome':
            attrs['type'] = 'protein-protein'
        elif G.nodes[u].get('layer') == 'Proteome' and G.nodes[v].get('layer') == 'Metabolome':
            attrs['type'] = 'protein-metabolite'
        elif G.nodes[u].get('layer') == 'Pathways' and G.nodes[v].get('layer') == 'Proteome':
            attrs['type'] = 'pathway-protein'
        else:
            attrs['type'] = 'other'
    
    return G, transnet

def create_visualizations(G, output_dir):
    """Create various visualizations and save them to the output directory."""
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Basic network visualization
    logger.info("Creating basic network visualization")
    fig1 = plot_network_basic(
        G,
        title="Basic Network Visualization",
        figsize=(10, 8),
        save_path=os.path.join(output_dir, "basic_network.png")
    )
    plt.close(fig1)
    
    # 2. Multi-layer network visualization
    logger.info("Creating multi-layer network visualization")
    fig2 = plot_multilayer_network(
        G,
        title="Multi-layer Network Visualization",
        figsize=(12, 10),
        layer_vertical_spacing=4.0,
        save_path=os.path.join(output_dir, "multilayer_network.png")
    )
    plt.close(fig2)
    
    # 3. Network with experimental data
    logger.info("Creating network visualization with experimental data")
    fig3 = plot_multilayer_network(
        G,
        title="Network with Fold Change Data",
        figsize=(12, 10),
        node_color_by='experimental',
        experimental_attr='fc',
        color_map='coolwarm',
        save_path=os.path.join(output_dir, "network_with_experimental_data.png")
    )
    plt.close(fig3)
    
    # 4. Community detection and visualization
    logger.info("Detecting communities and creating visualization")
    try:
        communities = detect_communities(G, method='louvain')
        fig4 = plot_community_network(
            G,
            communities=communities,
            title="Network Communities",
            figsize=(12, 10),
            save_path=os.path.join(output_dir, "community_network.png")
        )
        plt.close(fig4)
    except Exception as e:
        logger.error(f"Error in community detection: {e}")
    
    # 5. Network metrics visualization
    logger.info("Creating network metrics visualization")
    fig5 = plot_network_metrics(
        G,
        metrics=['degree', 'betweenness'],
        top_n=10,
        figsize=(12, 8),
        title="Network Centrality Metrics",
        save_path=os.path.join(output_dir, "network_metrics.png"),
        color_by_layer=True
    )
    plt.close(fig5)
    
    # 6. Heatmap of experimental data
    logger.info("Creating experimental data heatmap")
    # Get nodes with experimental data
    nodes_with_data = [
        node for node in G.nodes() 
        if G.nodes[node].get('fc') is not None
    ]
    
    if nodes_with_data:
        fig6 = plot_subnetwork_heatmap(
            G,
            node_list=nodes_with_data[:15],  # Limit to 15 nodes for readability
            attribute='fc',
            figsize=(10, 8),
            title="Fold Change Heatmap",
            color_map='coolwarm',
            save_path=os.path.join(output_dir, "experimental_heatmap.png"),
            cluster=True
        )
        plt.close(fig6)
    
    # 7. Interactive HTML visualization
    logger.info("Creating interactive HTML visualization")
    html = create_interactive_html_network(
        G,
        title="Interactive Network Visualization",
        width=1000,
        height=800,
        color_by='layer',
        size_by='degree',
        experimental_attr='fc',
        node_hover_info=['layer', 'fc', 'adj_p_value'],
        save_path=os.path.join(output_dir, "interactive_network.html")
    )
    
    logger.info(f"All visualizations saved to {output_dir}")

def main():
    """Main function to demonstrate network visualizations."""
    # Build sample network
    G, transnet = build_sample_network()
    
    # Create output directory
    output_dir = "visualization_examples"
    
    # Create and save visualizations
    create_visualizations(G, output_dir)
    
    logger.info(f"Network visualization examples saved to {output_dir}")
    logger.info("Open interactive_network.html in a web browser to view the interactive visualization")

if __name__ == "__main__":
    main()
    