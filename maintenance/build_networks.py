"""
Build pre-made networks for common model organisms.
This script is intended to be run periodically via GitHub Actions to maintain
up-to-date network databases for common model organisms.
"""

import os
import argparse
import logging
from datetime import datetime
from transnet.biology.transnet import Transnet
from transnet.biology.layers import Pathways, Reactions, Proteome, Metabolome, Transcriptome
from transnet.api.brenda import BRENDA_api

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("network_build.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def build_organism_network(organism, output_dir):
    """
    Build a network for a specific organism.
    
    Parameters:
    -----------
    organism : str
        Organism identifier (human, mouse, yeast, ecoli)
    output_dir : str
        Directory to save the network
    """
    # Organism config
    organism_config = {
        "human": {
            "kegg_org": "hsa",
            "ncbi_org": "9606",
            "ensembl_org": "homo_sapiens",
            "genome_chip": "hg38"
        },
        "mouse": {
            "kegg_org": "mmu",
            "ncbi_org": "10090",
            "ensembl_org": "mus_musculus",
            "genome_chip": "mm10"
        },
        "yeast": {
            "kegg_org": "sce",
            "ncbi_org": "4932",
            "ensembl_org": "saccharomyces_cerevisiae",
            "genome_chip": "sacCer3"
        },
        "ecoli": {
            "kegg_org": "eco",
            "ncbi_org": "511145",
            "ensembl_org": "escherichia_coli",
            "genome_chip": "eschColi_K12"
        }
    }
    
    if organism not in organism_config:
        logger.error(f"Unknown organism: {organism}")
        return False
    
    config = organism_config[organism]
    logger.info(f"Building network for {organism}")
    
    try:
        # Initialize BRENDA API
        brenda_api = None  # Replace with your BRENDA API credentials
        
        # Create network layers
        pathways = Pathways()
        pathways.kegg_organism = config["kegg_org"]
        pathways.populate()
        pathways.fill_pathways()
        
        reactions = Reactions()
        reactions.populate(from_api=True)
        
        proteome = Proteome()
        proteome.ncbi_organism = config["ncbi_org"]
        proteome.kegg_organism = config["kegg_org"]
        proteome.populate(uniprot=True)
        proteome.get_interaction_partners()
        proteome.get_metabolites()
        if brenda_api:
            proteome.get_activators_and_inhibitors(brenda_api)
        
        metabolome = Metabolome()
        # Populate with all KEGG compounds
        metabolome.populate()
        
        transcriptome = Transcriptome()
        transcriptome.kegg_organism = config["kegg_org"]
        transcriptome.organism_full = config["ensembl_org"]
        transcriptome.populate(ensembl=True)
        transcriptome.fill_gene_info()
        
        # Create the integrated network
        transnet = Transnet(
            name=f"{organism}_network",
            pathways=pathways,
            reactions=reactions,
            proteome=proteome,
            metabolome=metabolome,
            transcriptome=transcriptome
        )
        
        # Add transcription factors to proteins
        proteome.get_transcription_factor_targets(
            genome_ChIP=config["genome_chip"],
            distance_ChIP=5
        )
        
        # Save the network
        timestamp = datetime.now().strftime("%Y%m%d")
        org_dir = os.path.join(output_dir, organism, timestamp)
        os.makedirs(org_dir, exist_ok=True)
        
        transnet.save_network(org_dir)
        
        # Create a latest symlink
        latest_dir = os.path.join(output_dir, organism, "latest")
        if os.path.exists(latest_dir) and os.path.islink(latest_dir):
            os.unlink(latest_dir)
        
        os.symlink(timestamp, latest_dir, target_is_directory=True)
        
        logger.info(f"Successfully built network for {organism}")
        return True
    
    except Exception as e:
        logger.error(f"Error building network for {organism}: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Build integrated networks for model organisms")
    parser.add_argument("--organisms", nargs="+", default=["human", "mouse", "yeast", "ecoli"],
                        help="Organisms to build networks for")
    parser.add_argument("--output-dir", default="data",
                        help="Directory to save the networks")
    
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    for organism in args.organisms:
        build_organism_network(organism, args.output_dir)

if __name__ == "__main__":
    main()
