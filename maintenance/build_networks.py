"""Fix for the import error in build_networks.py.

This script adds the parent directory to the Python path before importing transnet modules.
"""

import os
import sys
import argparse
import logging
import traceback
from datetime import datetime
import pandas as pd

# Add the parent directory to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Now we can import from transnet
from transnet.biology.transnet import Transnet
from transnet.biology.layers import Pathways, Reactions, Proteome, Metabolome, Transcriptome

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

# Organism configuration
ORGANISM_CONFIG = {
    "human": {
        "kegg_org": "hsa",
        "ncbi_org": "9606",
        "ensembl_org": "homo_sapiens",
        "ensembl_release": 109,
        "genome_chip": "hg38"
    },
    "mouse": {
        "kegg_org": "mmu",
        "ncbi_org": "10090",
        "ensembl_org": "mus_musculus",
        "ensembl_release": 109,
        "genome_chip": "mm10"
    },
    "yeast": {
        "kegg_org": "sce",
        "ncbi_org": "4932",
        "ensembl_org": "saccharomyces_cerevisiae",
        "ensembl_release": 109,
        "genome_chip": "sacCer3"
    },
    "ecoli": {
        "kegg_org": "eco",
        "ncbi_org": "511145",
        "ensembl_org": "escherichia_coli_str_k_12_substr_mg1655",
        "ensembl_release": 109,
        "genome_chip": "eschColi_K12"
    }
}

def build_organism_network(organism, output_dir, debug=False):
    """
    Build a network for a specific organism.
    
    Parameters:
    -----------
    organism : str
        Organism identifier (human, mouse, yeast, ecoli)
    output_dir : str
        Directory to save the network
    debug : bool
        If True, run in debug mode with smaller network components
    
    Returns:
    --------
    bool
        Whether the network build was successful
    """
    if organism not in ORGANISM_CONFIG:
        logger.error(f"Unknown organism: {organism}")
        return False
    
    config = ORGANISM_CONFIG[organism]
    logger.info(f"Building network for {organism}")
    
    try:
        # Create network layers
        # 1. Pathways layer
        logger.info("Building pathways layer")
        pathways = Pathways()
        pathways.kegg_organism = config["kegg_org"]
        pathways.populate()
        pathways.fill_pathways()
        logger.info(f"Created pathways layer with {len(pathways.pathways)} pathways")
        
        # 2. Reactions layer
        logger.info("Building reactions layer")
        reactions = Reactions()
        reactions.populate(from_api=True)
        logger.info(f"Created reactions layer with {len(reactions.reactions)} reactions")
        
        # 3. Proteome layer
        logger.info("Building proteome layer")
        proteome = Proteome()
        proteome.ncbi_organism = config["ncbi_org"]
        proteome.kegg_organism = config["kegg_org"]
        proteome.populate(uniprot=True)
        
        if debug:
            # In debug mode, limit to 100 proteins for faster testing
            proteome.proteins = proteome.proteins[:100] if len(proteome.proteins) > 100 else proteome.proteins
            
        logger.info(f"Populated proteome with {len(proteome.proteins)} proteins")
        
        logger.info("Getting protein-protein interactions")
        proteome.get_interaction_partners()
        
        logger.info("Getting protein-metabolite associations")
        proteome.get_metabolites()
        
        logger.info(f"Completed proteome layer for {organism}")
        
        # 4. Metabolome layer
        logger.info("Building metabolome layer")
        metabolome = Metabolome()
        # Populate with all KEGG compounds
        metabolome.populate()
        logger.info(f"Created metabolome layer with {len(metabolome.metabolites)} metabolites")
        
        # 5. Transcriptome layer
        logger.info("Building transcriptome layer")
        transcriptome = Transcriptome()
        transcriptome.kegg_organism = config["kegg_org"]
        transcriptome.organism_full = config["ensembl_org"]
        
        if debug:
            # In debug mode, use simplified transcriptome population
            logger.info("Running in debug mode - using simplified transcriptome population")
            # Just populate with a minimal set using KEGG instead of Ensembl
            kegg_genes = pd.DataFrame()
            transcriptome.genes = []
        else:
            try:
                transcriptome.populate(
                    ensembl=True,
                    ensembl_release=config["ensembl_release"]
                )
                transcriptome.fill_gene_info()
                logger.info(f"Created transcriptome layer with {len(transcriptome.genes)} genes")
            except Exception as e:
                logger.error(f"Error populating transcriptome: {e}")
                logger.error(traceback.format_exc())
                logger.warning("Will continue with empty transcriptome layer")
                transcriptome.genes = []
        
        # Create the integrated network
        logger.info("Creating integrated network")
        transnet = Transnet(
            name=f"{organism}_network",
            pathways=pathways,
            reactions=reactions,
            proteome=proteome,
            metabolome=metabolome,
            transcriptome=transcriptome
        )
        
        # Add transcription factors to proteins
        if not debug:
            logger.info("Getting transcription factor targets")
            try:
                proteome.get_transcription_factor_targets(
                    genome_ChIP=config["genome_chip"],
                    distance_ChIP=5
                )
            except Exception as e:
                logger.error(f"Error getting transcription factor targets: {e}")
                logger.error(traceback.format_exc())
                logger.warning("Will continue without transcription factor targets")
        
        # Save the network
        timestamp = datetime.now().strftime("%Y%m%d")
        org_dir = os.path.join(output_dir, organism, timestamp)
        os.makedirs(org_dir, exist_ok=True)
        
        logger.info(f"Saving network to {org_dir}")
        transnet.save_network(org_dir)
        
        # Create a latest symlink
        latest_dir = os.path.join(output_dir, organism, "latest")
        if os.path.exists(latest_dir):
            if os.path.islink(latest_dir):
                logger.info(f"Removing existing symlink at {latest_dir}")
                os.unlink(latest_dir)
            else:
                logger.warning(f"{latest_dir} exists but is not a symlink. Removing...")
                import shutil
                shutil.rmtree(latest_dir)
        
        # Create relative symlink
        logger.info(f"Creating 'latest' symlink to {timestamp}")
        os.symlink(timestamp, latest_dir, target_is_directory=True)
        
        # Create a simple summary file
        summary = {
            "organism": organism,
            "build_date": timestamp,
            "pathways_count": len(pathways.pathways),
            "reactions_count": len(reactions.reactions),
            "proteins_count": len(proteome.proteins),
            "metabolites_count": len(metabolome.metabolites),
            "genes_count": len(transcriptome.genes) if hasattr(transcriptome, 'genes') else 0
        }
        
        with open(os.path.join(org_dir, "summary.txt"), "w") as f:
            for key, value in summary.items():
                f.write(f"{key}: {value}\n")
        
        logger.info(f"Successfully built network for {organism}")
        return True
    
    except Exception as e:
        logger.error(f"Error building network for {organism}: {e}")
        logger.error(traceback.format_exc())
        return False

def main():
    parser = argparse.ArgumentParser(description="Build integrated networks for model organisms")
    parser.add_argument("--organisms", nargs="+", default=["human", "mouse", "yeast", "ecoli"],
                        help="Organisms to build networks for")
    parser.add_argument("--output-dir", default="data",
                        help="Directory to save the networks")
    parser.add_argument("--debug", action="store_true",
                        help="Run in debug mode with smaller network components")
    
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    results = {}
    for organism in args.organisms:
        logger.info(f"======= Starting build for {organism} =======")
        success = build_organism_network(organism, args.output_dir, args.debug)
        results[organism] = "Success" if success else "Failed"
        logger.info(f"======= Completed build for {organism}: {results[organism]} =======")
    
    # Print summary
    logger.info("===== Build Summary =====")
    for organism, result in results.items():
        logger.info(f"{organism}: {result}")
    logger.info("========================")

if __name__ == "__main__":
    main()
    