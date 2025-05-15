"""
Further modified version of build_networks.py with better error handling and dependency handling.
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

# Check for required packages and install if missing
try:
    import pyensembl
except ImportError:
    logging.warning("pyensembl package not found. Attempting to install it.")
    import subprocess
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pyensembl"])
        import pyensembl
    except Exception as e:
        logging.error(f"Failed to install pyensembl: {e}")
        logging.error("Will continue with limited functionality")

# Now we can import from transnet
try:
    from transnet.biology.transnet import Transnet
    from transnet.biology.layers import Pathways, Reactions, Proteome, Metabolome, Transcriptome
except ImportError as e:
    logging.error(f"Import error: {e}")
    logging.error("This may be due to missing dependencies. Make sure all dependencies are installed.")
    sys.exit(1)

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
        try:
            pathways = Pathways()
            pathways.kegg_organism = config["kegg_org"]
            pathways.populate()
            pathways.fill_pathways()
            logger.info(f"Created pathways layer with {len(pathways.pathways)} pathways")
        except Exception as e:
            logger.error(f"Error creating pathways layer: {e}")
            logger.error(traceback.format_exc())
            pathways = Pathways()
            pathways.pathways = []
            logger.warning("Continuing with empty pathways layer")
        
        # 2. Reactions layer
        logger.info("Building reactions layer")
        try:
            reactions = Reactions()
            reactions.populate(from_api=True)
            logger.info(f"Created reactions layer with {len(reactions.reactions)} reactions")
        except Exception as e:
            logger.error(f"Error creating reactions layer: {e}")
            logger.error(traceback.format_exc())
            reactions = Reactions()
            reactions.reactions = []
            logger.warning("Continuing with empty reactions layer")
        
        # 3. Proteome layer
        logger.info("Building proteome layer")
        try:
            proteome = Proteome()
            proteome.ncbi_organism = config["ncbi_org"]
            proteome.kegg_organism = config["kegg_org"]
            proteome.populate(uniprot=True)
            
            if debug:
                # In debug mode, limit to 100 proteins for faster testing
                proteome.proteins = proteome.proteins[:100] if len(proteome.proteins) > 100 else proteome.proteins
                
            logger.info(f"Populated proteome with {len(proteome.proteins)} proteins")
            
            # Skip interaction retrieval in debug mode to make it faster
            if not debug:
                logger.info("Getting protein-protein interactions")
                try:
                    proteome.get_interaction_partners()
                except Exception as e:
                    logger.error(f"Error getting protein-protein interactions: {e}")
                    logger.error(traceback.format_exc())
                    logger.warning("Continuing without protein-protein interactions")
            else:
                logger.info("Skipping protein-protein interactions in debug mode")
            
            logger.info("Getting protein-metabolite associations")
            try:
                proteome.get_metabolites()
            except Exception as e:
                logger.error(f"Error getting protein-metabolite associations: {e}")
                logger.error(traceback.format_exc())
                logger.warning("Continuing without protein-metabolite associations")
            
            logger.info(f"Completed proteome layer for {organism}")
        except Exception as e:
            logger.error(f"Error creating proteome layer: {e}")
            logger.error(traceback.format_exc())
            proteome = Proteome()
            proteome.proteins = []
            logger.warning("Continuing with empty proteome layer")
        
        # 4. Metabolome layer
        logger.info("Building metabolome layer")
        try:
            metabolome = Metabolome()
            # Populate with all KEGG compounds
            metabolome.populate()
            logger.info(f"Created metabolome layer with {len(metabolome.metabolites)} metabolites")
        except Exception as e:
            logger.error(f"Error creating metabolome layer: {e}")
            logger.error(traceback.format_exc())
            metabolome = Metabolome()
            metabolome.metabolites = []
            logger.warning("Continuing with empty metabolome layer")
        
        # 5. Transcriptome layer
        logger.info("Building transcriptome layer")
        transcriptome = Transcriptome()
        transcriptome.kegg_organism = config["kegg_org"]
        transcriptome.organism_full = config["ensembl_org"]
        
        if debug:
            # In debug mode, use simplified transcriptome population
            logger.info("Running in debug mode - using simplified transcriptome population")
            # Just populate with a minimal set using KEGG
            try:
                import numpy as np
                from transnet.api.kegg import kegg_list_genes
                
                # Get a small set of genes from KEGG
                kegg_genes = kegg_list_genes(config["kegg_org"])
                if len(kegg_genes) > 0:
                    # Sample a small number of genes
                    indices = np.random.choice(len(kegg_genes), min(50, len(kegg_genes)), replace=False)
                    sampled_genes = kegg_genes.iloc[indices]
                    
                    # Create Gene objects
                    transcriptome.genes = []
                    for _, row in sampled_genes.iterrows():
                        from transnet.biology.elements import Gene
                        gene = Gene(
                            kegg_id=row['gene_id'], 
                            name=row['name(s)'], 
                            description=row['description'],
                            kegg_organism=config["kegg_org"]
                        )
                        transcriptome.genes.append(gene)
                    
                    logger.info(f"Created transcriptome layer with {len(transcriptome.genes)} genes from KEGG")
                else:
                    logger.warning("No genes found in KEGG for this organism")
                    transcriptome.genes = []
            except Exception as e:
                logger.error(f"Error creating simplified transcriptome layer: {e}")
                logger.error(traceback.format_exc())
                transcriptome.genes = []
                logger.warning("Continuing with empty transcriptome layer")
        else:
            try:
                # Skip pyensembl in GitHub Actions to avoid installation issues
                if 'GITHUB_ACTIONS' in os.environ:
                    logger.info("Running in GitHub Actions - skipping Ensembl transcriptome population")
                    transcriptome.genes = []
                else:
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
        if not debug and len(proteome.proteins) > 0:
            logger.info("Getting transcription factor targets")
            try:
                proteome.get_transcription_factor_targets(
                    genome_ChIP=config["genome_chip"],
                    distance_ChIP=5
                )
            except Exception as e:
                logger.error(f"Error getting transcription factor targets: {e}")
                logger.error(traceback.format_exc())
                logger.warning("Continuing without transcription factor targets")
        else:
            logger.info("Skipping transcription factor targets in debug mode")
        
        # Save the network
        timestamp = datetime.now().strftime("%Y%m%d")
        org_dir = os.path.join(output_dir, organism, timestamp)
        os.makedirs(org_dir, exist_ok=True)
        
        logger.info(f"Saving network to {org_dir}")
        try:
            transnet.save_network(org_dir)
        except Exception as e:
            logger.error(f"Error saving network: {e}")
            logger.error(traceback.format_exc())
            
            # Try to save minimal information
            logger.info("Attempting to save minimal network information")
            try:
                # Save minimal information
                with open(os.path.join(org_dir, "network_info.txt"), "w") as f:
                    f.write(f"Network for {organism}\n")
                    f.write(f"Build date: {timestamp}\n")
                    f.write(f"Pathways: {len(pathways.pathways)}\n")
                    f.write(f"Reactions: {len(reactions.reactions)}\n")
                    f.write(f"Proteins: {len(proteome.proteins)}\n")
                    f.write(f"Metabolites: {len(metabolome.metabolites)}\n")
                    f.write(f"Genes: {len(transcriptome.genes)}\n")
            except Exception as e2:
                logger.error(f"Error saving minimal information: {e2}")
        
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
        try:
            logger.info(f"Creating 'latest' symlink to {timestamp}")
            os.symlink(timestamp, latest_dir, target_is_directory=True)
        except Exception as e:
            logger.error(f"Error creating symlink: {e}")
            logger.warning("Continuing without symlink")
        
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
    
    logger.info(f"Starting build with: organisms={args.organisms}, output_dir={args.output_dir}, debug={args.debug}")
    logger.info(f"Python version: {sys.version}")
    logger.info(f"Current directory: {os.getcwd()}")
    
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
    