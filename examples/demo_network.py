"""
Demo script to showcase how to build and analyze a network using transnet.
This script builds a small network for the specified organism and demonstrates
basic analysis functionality.
"""

import os
import argparse
import logging
import matplotlib.pyplot as plt
from datetime import datetime

from transnet.biology.transnet import Transnet
from transnet.biology.layers import Pathways, Reactions, Proteome, Metabolome, Transcriptome
from transnet.analysis.network_analysis import (
    compute_network_statistics, 
    identify_hubs,
    compute_centrality_measures,
    detect_communities,
    plot_network
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Organism configuration - simplified for demo
ORGANISM_CONFIG = {
    "human": {
        "kegg_org": "hsa",
        "ncbi_org": "9606",
        "ensembl_org": "homo_sapiens"
    },
    "mouse": {
        "kegg_org": "mmu",
        "ncbi_org": "10090",
        "ensembl_org": "mus_musculus"
    },
    "yeast": {
        "kegg_org": "sce",
        "ncbi_org": "4932",
        "ensembl_org": "saccharomyces_cerevisiae"
    },
    "ecoli": {
        "kegg_org": "eco",
        "ncbi_org": "511145",
        "ensembl_org": "escherichia_coli_str_k_12_substr_mg1655"
    }
}

def build_small_demo_network(organism):
    """
    Build a small demo network for the specified organism.
    
    Parameters:
    -----------
    organism : str
        Organism identifier (human, mouse, yeast, ecoli)
        
    Returns:
    --------
    Transnet
        A small network for the specified organism
    """
    if organism not in ORGANISM_CONFIG:
        raise ValueError(f"Unknown organism: {organism}. Available options: {list(ORGANISM_CONFIG.keys())}")
    
    config = ORGANISM_CONFIG[organism]
    logger.info(f"Building demo network for {organism}")
    
    # Create network layers - limited scope for demo
    
    # 1. Pathways layer - limit to 5 pathways for demo
    logger.info("Building pathways layer")
    pathways = Pathways()
    pathways.kegg_organism = config["kegg_org"]
    pathways.populate()
    pathways.pathways = pathways.pathways[:5]  # Just use first 5 pathways
    pathways.fill_pathways()
    logger.info(f"Created pathways layer with {len(pathways.pathways)} pathways")
    
    # 2. Proteome layer - limit to proteins in selected pathways
    logger.info("Building proteome layer")
    proteome = Proteome()
    proteome.ncbi_organism = config["ncbi_org"]
    proteome.kegg_organism = config["kegg_org"]
    proteome.populate(uniprot=True)
    
    # Filter proteins to only those in the selected pathways
    genes_in_pathways = set()
    for pathway in pathways.pathways:
        genes_in_pathways.update(pathway.genes)
    
    # Keep only 50 proteins max for the demo
    proteome.proteins = proteome.proteins[:50]
    logger.info(f"Populated proteome with {len(proteome.proteins)} proteins")
    
    # Get protein-protein interactions
    proteome.get_interaction_partners()
    logger.info("Added protein-protein interactions")
    
    # 3. Metabolome layer - just use a few metabolites
    logger.info("Building metabolome layer")
    metabolome = Metabolome()
    metabolome.populate()
    metabolome.metabolites = metabolome.metabolites[:30]  # Just use 30 metabolites
    logger.info(f"Created metabolome layer with {len(metabolome.metabolites)} metabolites")
    
    # Create the integrated network
    logger.info("Creating integrated network")
    transnet = Transnet(
        name=f"{organism}_demo_network",
        pathways=pathways,
        proteome=proteome,
        metabolome=metabolome
    )
    
    return transnet

def analyze_network(transnet):
    """
    Perform basic analysis on the network.
    
    Parameters:
    -----------
    transnet : Transnet
        The network to analyze
    """
    logger.info("Generating network graph")
    G = transnet.generate_graph()
    
    # Basic statistics
    logger.info("Computing network statistics")
    stats = compute_network_statistics(G)
    print("\nNetwork Statistics:")
    for key, value in stats.items():
        print(f"  {key}: {value}")
    
    # Hub identification
    logger.info("Identifying hub nodes")
    hubs = identify_hubs(G, top_n=5)
    print("\nTop 5 Hub Nodes:")
    for node, degree in hubs:
        print(f"  {node}: {degree} connections")
    
    # Centrality measures
    logger.info("Computing centrality measures")
    centrality = compute_centrality_measures(G, top_n=5)
    print("\nTop Nodes by Betweenness Centrality:")
    for node, value in centrality['betweenness_centrality']:
        print(f"  {node}: {value:.4f}")
    
    # Community detection
    logger.info("Detecting communities")
    try:
        communities = detect_communities(G, method='louvain')
        community_sizes = {}
        for node, community_id in communities.items():
            if community_id not in community_sizes:
                community_sizes[community_id] = 0
            community_sizes[community_id] += 1
        
        print("\nCommunity Distribution:")
        for community_id, size in sorted(community_sizes.items(), key=lambda x: x[1], reverse=True)[:5]:
            print(f"  Community {community_id}: {size} nodes")
    except Exception as e:
        logger.error(f"Error in community detection: {e}")
        print("\nCommunity detection failed. This often requires additional packages.")
    
    # Network visualization
    logger.info("Creating network visualization")
    try:
        logger.info("Plotting network")
        plt.figure(figsize=(10, 8))
        plot_network(
            G, 
            title=f"{transnet.name} Visualization",
            layout='spring'
        )
        
        # Save the plot
        output_dir = "output"
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, f"{transnet.name}_visualization.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Network visualization saved to {output_file}")
        
        # Also save the network
        network_file = os.path.join(output_dir, f"{transnet.name}")
        os.makedirs(network_file, exist_ok=True)
        transnet.save_network(network_file)
        logger.info(f"Network saved to {network_file}")
        
    except Exception as e:
        logger.error(f"Error in network visualization: {e}")
        print("\nNetwork visualization failed.")

def main():
    parser = argparse.ArgumentParser(description="Build and analyze a demo network")
    parser.add_argument("--organism", default="mouse", choices=list(ORGANISM_CONFIG.keys()),
                        help="Organism to build network for")
    
    args = parser.parse_args()
    
    logger.info(f"Starting demo for {args.organism}")
    
    # Build the network
    transnet = build_small_demo_network(args.organism)
    
    # Analyze the network
    analyze_network(transnet)
    
    logger.info("Demo completed successfully")

if __name__ == "__main__":
    main()
    