"""Biological layers (omics), layers in the network
"""

from typing import List, Dict, Any, Optional, Union
import pandas as pd
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

__all__ = [
    "Reactions",
    "Pathways",
    "Transcriptome",
    "Proteome",
    "Metabolome",
]

# Import elements after defining __all__ to avoid circular imports
from transnet.biology.elements import Reaction, Pathway, Metabolite, Gene, Protein
from transnet.api.kegg import (
    kegg_create_reaction_table, kegg_list_pathways, kegg_link_pathway, 
    kegg_link_ec, kegg_list_genes, kegg_conv_ncbi_idtable, kegg_ec_to_cpds, 
    kegg_list_compounds, kegg_to_chebi, chebi_to_kegg
)
from transnet.api.uniprot import uniprot_list_proteins, uniprot_add_entrez_id
from transnet.api.string import string_map_identifiers, string_get_interactions, translate_string_dict
from transnet.api.ensembl import ensembl_download_release, ensembl_get_transcripts, ensembl_expander
from transnet.api.chip_atlas import get_ChIP_data, get_ChIP_exps
from transnet.api.chem_info import chemical_info_converter

class OmicsLayer:
    """Base class for all omics layers."""
    
    def __init__(self, name: str = None):
        self.name = name
        self.elements = []
        
    def add_experimental_data(self, 
                             input_data: pd.DataFrame = None, 
                             id_column: str = None, 
                             id_type: str = None):
        """Add experimental data to this layer."""
        pass
        
    def populate(self):
        """Populate this layer with data from databases."""
        pass
        
    def __repr__(self):
        return f"<{self.__class__.__name__}: {len(self.elements)} elements>"


class Reactions(OmicsLayer):
    """
    Reactions layer representing biochemical reactions.
    """
    
    def __init__(self):
        super().__init__(name="Reactions")
        self.reactions = []
    
    def __repr__(self):
        return f"<Reactions: {len(self.reactions)} reactions in total>"
    
    def populate(
            self,
            from_api: bool = True,
            df = None,
            ):
        """
        Populate reactions from KEGG API or from a dataframe.
        
        Parameters:
        -----------
        from_api : bool
            Whether to populate reactions from KEGG API
        df : pd.DataFrame
            Dataframe containing reaction information
        """
        if from_api:
            logger.info("Populating reactions from KEGG API")
            reaction_table = kegg_create_reaction_table()
            for _, row in reaction_table.iterrows():
                new_reaction = Reaction(
                    id = row['reaction'],
                    name = row['name'],
                    equation = row['equation'],
                    definition = row['definition'],
                    enzyme = row['enzyme'],
                    substrates = row['substrates'],
                    products = row['products'],
                    stoichiometry_substrates = row['stoichiometry_substrates'],
                    stoichiometry_products = row['stoichiometry_products'],
                )
                self.reactions.append(new_reaction)
            logger.info(f"Populated {len(self.reactions)} reactions from KEGG API")
        elif df is not None:
            logger.info("Populating reactions from dataframe")
            reaction_table = df
            for _, row in reaction_table.iterrows():
                new_reaction = Reaction(
                    id = row['reaction'],
                    name = row['name'],
                    equation = row['equation'],
                    definition = row['definition'],
                    enzyme = row['enzyme'],
                    substrates = row['substrates'],
                    products = row['products'],
                    stoichiometry_substrates = row['stoichiometry_substrates'],
                    stoichiometry_products = row['stoichiometry_products'],
                )
                self.reactions.append(new_reaction)
            logger.info(f"Populated {len(self.reactions)} reactions from dataframe")
        else:
            logger.warning("No data source provided for reaction population")


class Pathways(OmicsLayer):
    """
    Pathways layer representing biological pathways.
    """

    def __init__(self):
        super().__init__(name="Pathways")
        self.kegg_organism = None
        self.pathways = []
        self.organism_full = None

    def __repr__(self):
        return f"<Pathways of: {self.kegg_organism}, {len(self.pathways)} pathways in total>"

    def populate(self):
        """
        Populate pathways from KEGG API.
        """
        if not self.kegg_organism:
            logger.error("KEGG organism code must be set before populating pathways")
            return
            
        logger.info(f"Populating pathways for {self.kegg_organism} from KEGG API")
        for _, row in kegg_list_pathways(self.kegg_organism).iterrows():
            new_pathway = Pathway(
                id = row[0],
                name = row[1],
                kegg_organism = self.kegg_organism,
            )
            self.pathways.append(new_pathway)
        logger.info(f"Populated {len(self.pathways)} pathways from KEGG API")

    def fill_pathways(self):
        """
        Fill in additional pathway information from KEGG.
        """
        if not self.kegg_organism:
            logger.error("KEGG organism code must be set before filling pathways")
            return
            
        logger.info(f"Filling pathway information for {self.kegg_organism}")
        kegg_pathway_links = kegg_link_pathway(self.kegg_organism)
        kegg_pathway_ecs = kegg_link_ec(self.kegg_organism)
        kegg_pathway_list = kegg_list_pathways(self.kegg_organism)

        for pathway in self.pathways:
            # Checks which info is missing and fills it
            # based on what is available
            # Start with KEGG ID
            if pathway.id is not None:
                # Check if name is there
                if pathway.name is None:
                    try:
                        pathway.name = kegg_pathway_list[
                            "description"
                        ][
                            kegg_pathway_list["pathways_id"]
                            == pathway.id
                        ].to_string(
                            index=False
                        )
                    except Exception as e:
                        logger.warning(f"Error getting pathway name for {pathway.id}: {e}")
                
                # Check if genes are there
                if len(pathway.genes) == 0:
                    try:
                        pathway.genes = kegg_pathway_links[
                            kegg_pathway_links['pathway'] == pathway.id][
                            "kegg_gene_id"
                        ].to_list()
                    except Exception as e:
                        logger.warning(f"Error getting pathway genes for {pathway.id}: {e}")
                
                # Check if EC numbers are there
                if len(pathway.ecs) == 0:
                    protein_ecs = []
                    for gene in kegg_pathway_links[
                        kegg_pathway_links['pathway'] == pathway.id][
                        "kegg_gene_id"
                    ].to_list():
                        if gene in kegg_pathway_ecs["kegg_gene_id"].to_list():
                            try:
                                ecs_tmp = kegg_pathway_ecs[
                                    kegg_pathway_ecs["kegg_gene_id"] == gene
                                ]["ec_number"]
                                if len(ecs_tmp) == 1:
                                    protein_ecs.append(ecs_tmp.to_string(index=False))
                                else:
                                    for ec in ecs_tmp.to_list():
                                        protein_ecs.append(ec)
                            except Exception as e:
                                logger.warning(f"Error getting EC numbers for {gene}: {e}")
                    pathway.ecs = list(set(protein_ecs))
        
        logger.info(f"Successfully filled pathway information")


class Metabolome(OmicsLayer):
    """
    Metabolome layer representing metabolites.
    """

    def __init__(self):
        super().__init__(name="Metabolome")
        self.metabolites = []
        self.input_id_type = None
        self.input_ids = []
        self.input_data = None

    def __repr__(self):
        return f"<Metabolome, {len(self.metabolites)} metabolites in total>"
    
    def add_experimental_data(
        self,
        input_data: pd.DataFrame = None,
        metabolite_column_name: str = None,
        id_type: str = None,
    ):
        """
        Add experimental metabolomics data.
        
        Parameters:
        -----------
        input_data : pd.DataFrame
            Dataframe containing experimental data
        metabolite_column_name : str
            Column name containing metabolite identifiers
        id_type : str
            Type of identifier (chebi, kegg, pubchem, inchi, inchikey, smile)
        """
        self.input_data = input_data
        self.input_ids = input_data[metabolite_column_name].tolist()
        
        valid_id_types = ["chebi", "kegg", "pubchem", "inchi", "inchikey", "smile"]
        if id_type in valid_id_types:
            self.input_id_type = id_type
        else:
            logger.warning(f"Unknown ID type: {id_type}. Valid types are: {valid_id_types}")

    def populate(
        self,
        metabolite_column_name: str = None,
        input_data_value_column_name: str = None,
        input_data_p_value_column_name: str = None,
    ):
        """
        Populate metabolome with data.
        
        Parameters:
        -----------
        metabolite_column_name : str
            Column name containing metabolite identifiers
        input_data_value_column_name : str
            Column name containing fold change values
        input_data_p_value_column_name : str
            Column name containing p-values
        """
        logger.info("Populating metabolome layer")
        self.metabolites = []
        
        # Get KEGG compounds
        try:
            kegg_compounds_df = kegg_list_compounds()
            kegg_compounds_df_sample = kegg_compounds_df.head(500)  # Use just a sample for testing
            
            # Convert KEGG to ChEBI
            try:
                kegg_converted = kegg_to_chebi(kegg_compounds_df_sample["kegg_compounds"].to_list())
                # Get chemical info
                converted_ids = chemical_info_converter(kegg_converted["chebi_compounds"].dropna().to_list())
                
                # Create metabolites
                for _, row in converted_ids.iterrows():
                    fc = None
                    adj_p_value = None
                    
                    try:
                        # Create metabolite object
                        new_metabolite = Metabolite(
                            pubchem_id = row.get('pubchem_ids'),
                            kegg_name = None,  # Would need to join with KEGG data
                            kegg_compound_id = None,  # Would need to join with KEGG data
                            inchi = row.get('inchi'),
                            inchikey = row.get('inchikeys'),
                            chebi_id = row.get('chebi_ids'),
                            smile = row.get('smiles'),
                            fc = fc,
                            adj_p_value = adj_p_value,
                        )
                        self.metabolites.append(new_metabolite)
                    except Exception as e:
                        logger.error(f"Error creating metabolite: {e}")
                        continue
            except Exception as e:
                logger.error(f"Error converting KEGG to ChEBI: {e}")
                # Fallback: just use KEGG compounds
                for _, row in kegg_compounds_df_sample.iterrows():
                    try:
                        new_metabolite = Metabolite(
                            kegg_compound_id=row.get("kegg_compounds"),
                            kegg_name=row.get("name")
                        )
                        self.metabolites.append(new_metabolite)
                    except Exception as e:
                        logger.error(f"Error creating metabolite from KEGG: {e}")
                        continue
        except Exception as e:
            logger.error(f"Error getting KEGG compounds: {e}")
            logger.warning("Using empty metabolome")
            
        logger.info(f"Populated metabolome with {len(self.metabolites)} metabolites")


class Transcriptome(OmicsLayer):
    """
    Transcriptome layer representing genes and transcripts.
    """

    def __init__(self):
        super().__init__(name="Transcriptome")
        self.kegg_organism = None
        self.organism_full = None
        self.input_id_type = None
        self.input_ids = []
        self.input_data = None
        self.genes = []

    def __repr__(self):
        return f"<Transcriptome of: {self.kegg_organism}, {len(self.genes)} genes in total>"

    def add_experimental_data(
        self,
        input_data: pd.DataFrame = None,
        transcript_column_name: str = None,
        id_type: str = None,
    ):
        """
        Add experimental transcriptomics data.
        
        Parameters:
        -----------
        input_data : pd.DataFrame
            Dataframe containing experimental data
        transcript_column_name : str
            Column name containing transcript identifiers
        id_type : str
            Type of identifier (ensembl, etc.)
        """
        self.input_data = input_data
        self.input_ids = input_data[transcript_column_name].tolist()
        if id_type == "ensembl":
            self.input_id_type = "ensembl"
        else:
            logger.warning(f"Unknown ID type: {id_type}. Currently only 'ensembl' is supported.")

    def populate(
        self,
        kegg_organism: str = None,
        organism_full: str = None,
        ensembl: bool = False,
        ensembl_release: int = 109,
        kegg_ftp: bool = False,
        kegg_api: bool = False,
        ftp_base_path: str = None,
        biotype: str = None,
        transcript_column_name: str = None,
        input_data_value_column_name: str = None,
        input_data_p_value_column_name: str = None,
    ):
        """
        Populate transcriptome with data.
        
        Parameters:
        -----------
        kegg_organism : str
            KEGG organism code
        organism_full : str
            Full organism name for Ensembl
        ensembl : bool
            Whether to use Ensembl as data source
        ensembl_release : int
            Ensembl release version
        biotype : str
            Filter genes by biotype
        transcript_column_name : str
            Column name containing transcript identifiers
        input_data_value_column_name : str
            Column name containing fold change values
        input_data_p_value_column_name : str
            Column name containing p-values
        """
        self.kegg_organism = kegg_organism or self.kegg_organism
        self.organism_full = organism_full or self.organism_full
        self.genes = []
        
        if not ensembl and not kegg_api and not kegg_ftp:
            # Default to KEGG API if no source specified
            kegg_api = True
            
        if kegg_api:
            logger.info(f"Populating genes from KEGG API for {self.kegg_organism}")
            try:
                kegg_genes = kegg_list_genes(self.kegg_organism)
                # Use only a sample for testing
                kegg_genes_sample = kegg_genes.head(100)
                
                for _, row in kegg_genes_sample.iterrows():
                    try:
                        new_gene = Gene(
                            kegg_id=row.get('gene_id'),
                            name=row.get('name(s)'),
                            description=row.get('description'),
                            kegg_organism=self.kegg_organism
                        )
                        self.genes.append(new_gene)
                    except Exception as e:
                        logger.error(f"Error creating gene from KEGG: {e}")
                        continue
                logger.info(f"Populated transcriptome with {len(self.genes)} genes from KEGG API")
            except Exception as e:
                logger.error(f"Error getting genes from KEGG API: {e}")
                logger.warning("Using empty transcriptome")
                
        elif ensembl:
            logger.info(f"Populating genes from Ensembl for {self.organism_full}")
            try:
                ensembl_download_release(
                    release_number=ensembl_release,
                    organism_full=self.organism_full,
                )
                transcripts_df = ensembl_get_transcripts(
                    release_number=ensembl_release,
                    organism_full=self.organism_full,
                )
                
                if biotype:
                    logger.info(f"Filtering transcripts by biotype: {biotype}")
                    if not isinstance(biotype, list):
                        biotype = [biotype]
                    transcripts_df = transcripts_df[
                        transcripts_df["biotype"].isin(biotype)
                    ]
                
                # Use only a sample for testing
                transcripts_df_sample = transcripts_df.head(100)
                
                ensembl_df = ensembl_expander(transcripts_df_sample)
                
                for _, row in ensembl_df.iterrows():
                    try:
                        fc = None
                        adj_p_value = None
                        
                        new_gene = Gene(
                            kegg_id=row.get("entrez_id"),
                            ncbi_id=row.get("entrez_id"),
                            ensembl_id=row.get("ensembl_gene_id"),
                            transcript_id=row.get("ensembl_transcript_id"),
                            uniprot_id=row.get("uniprot_id"),
                            type=row.get("biotype"),
                            name=row.get("gene_description"),
                            description=row.get("gene_name"),
                            kegg_organism=self.kegg_organism,
                            fc=fc,
                            adj_p_value=adj_p_value,
                        )
                        self.genes.append(new_gene)
                    except Exception as e:
                        logger.error(f"Error creating gene from Ensembl: {e}")
                        continue
                logger.info(f"Populated transcriptome with {len(self.genes)} genes from Ensembl")
            except Exception as e:
                logger.error(f"Error getting genes from Ensembl: {e}")
                logger.warning("Using empty transcriptome")
        
        elif kegg_ftp:
            logger.warning("KEGG FTP not implemented yet")
            

    def fill_gene_info(self):
        """
        Fill in additional gene information from KEGG.
        """
        if not self.kegg_organism:
            logger.error("KEGG organism code must be set before filling gene info")
            return
            
        logger.info(f"Filling gene information for {self.kegg_organism}")

        try:
            kegg_gene_list = kegg_list_genes(self.kegg_organism)
            kegg_ncbi_idtable = kegg_conv_ncbi_idtable(self.kegg_organism)
            kegg_ec_link = kegg_link_ec(self.kegg_organism)

            for gene in self.genes:
                # Start with KEGG ID
                if gene.kegg_id is not None:
                    # Check if name is there
                    if gene.name is None:
                        try:
                            kegg_info = kegg_gene_list[
                                kegg_gene_list["gene_id"] == gene.kegg_id
                            ]
                            if not kegg_info.empty:
                                gene.name = kegg_info["name(s)"].to_string(index=False)
                        except Exception as e:
                            logger.warning(f"Error getting gene name for {gene.kegg_id}: {e}")
                    
                    # Check if description is there
                    if gene.description is None:
                        try:
                            kegg_info = kegg_gene_list[
                                kegg_gene_list["gene_id"] == gene.kegg_id
                            ]
                            if not kegg_info.empty:
                                gene.description = kegg_info["description"].to_string(
                                    index=False
                                )
                        except Exception as e:
                            logger.warning(f"Error getting gene description for {gene.kegg_id}: {e}")
                    
                    # Check if ncbi_id is there
                    if gene.ncbi_id is None:
                        try:
                            ncbi_info = kegg_ncbi_idtable[
                                kegg_ncbi_idtable["kegg_id"] == gene.kegg_id
                            ]
                            if not ncbi_info.empty:
                                gene.ncbi_id = ncbi_info["ncbi_id"].to_string(index=False)
                        except Exception as e:
                            logger.warning(f"Error getting NCBI ID for {gene.kegg_id}: {e}")
                    
                    # Check if related EC numbers are there
                    if len(gene.related_ecs) == 0:
                        try:
                            protein_ecs = []
                            if gene.kegg_id in kegg_ec_link["kegg_gene_id"].to_list():
                                ecs_tmp = kegg_ec_link[
                                    kegg_ec_link["kegg_gene_id"] == gene.kegg_id
                                ]["ec_number"]
                                if len(ecs_tmp) == 1:
                                    protein_ecs.append(ecs_tmp.to_string(index=False))
                                else:
                                    for ec in ecs_tmp.to_list():
                                        protein_ecs.append(ec)
                            gene.related_ecs = list(set(protein_ecs))
                        except Exception as e:
                            logger.warning(f"Error getting EC numbers for {gene.kegg_id}: {e}")
            
            logger.info(f"Successfully filled gene information")
        except Exception as e:
            logger.error(f"Error filling gene information: {e}")


class Proteome(OmicsLayer):
    """
    Proteome layer representing proteins.
    """

    def __init__(self):
        super().__init__(name="Proteome")
        self.proteins = []
        self.organism_full = None
        self.input_id_type = None
        self.input_ids = []
        self.input_data = None
        self.ncbi_organism = None
        self.kegg_organism = None

    def __repr__(self):
        return f"<Proteome of: {self.ncbi_organism or self.kegg_organism}, {len(self.proteins)} proteins in total>"

    def add_experimental_data(
        self,
        input_data: pd.DataFrame = None,
        protein_column_name: str = None,
        id_type: str = None,
    ):
        """
        Add experimental proteomics data.
        
        Parameters:
        -----------
        input_data : pd.DataFrame
            Dataframe containing experimental data
        protein_column_name : str
            Column name containing protein identifiers
        id_type : str
            Type of identifier (uniprot, etc.)
        """
        self.input_data = input_data
        self.input_ids = input_data[protein_column_name].tolist()
        if id_type == "uniprot":
            self.input_id_type = "uniprot"
        else:
            logger.warning(f"Unknown ID type: {id_type}. Currently only 'uniprot' is supported.")

    def populate(
        self,
        kegg_organism: str = None,
        ncbi_organism: str = None,
        organism_full: str = None,
        ensembl: bool = False,
        uniprot: bool = True,
        ensembl_release: int = 109,
        protein_column_name: str = None,
        input_data_value_column_name: str = None,
        input_data_p_value_column_name: str = None,
    ):
        """
        Populate proteome with data.
        
        Parameters:
        -----------
        kegg_organism : str
            KEGG organism code
        ncbi_organism : str
            NCBI organism code
        organism_full : str
            Full organism name for Ensembl
        ensembl : bool
            Whether to use Ensembl as data source
        uniprot : bool
            Whether to use UniProt as data source
        ensembl_release : int
            Ensembl release version
        protein_column_name : str
            Column name containing protein identifiers
        input_data_value_column_name : str
            Column name containing fold change values
        input_data_p_value_column_name : str
            Column name containing p-values
        """
        self.kegg_organism = kegg_organism or self.kegg_organism
        self.ncbi_organism = ncbi_organism or self.ncbi_organism
        self.organism_full = organism_full or self.organism_full
        
        if not ensembl and not uniprot:
            logger.warning("Neither Ensembl nor UniProt selected as data source. Using UniProt as default.")
            uniprot = True
            
        if uniprot:
            logger.info(f"Populating proteins from UniProt for {self.ncbi_organism}")
            try:
                uniprot_proteins = uniprot_list_proteins(self.ncbi_organism)
                
                # Use only a sample for testing
                uniprot_proteins_sample = uniprot_proteins.head(50)
                
                # Add Entrez IDs
                uniprot_proteins_sample = uniprot_add_entrez_id(uniprot_proteins_sample)
                
                for _, row in uniprot_proteins_sample.iterrows():
                    try:
                        # Filter by input IDs if specified
                        if self.input_id_type == "uniprot" and row["Entry"] not in self.input_ids:
                            continue
                        
                        fc = None
                        adj_p_value = None
                        
                        # Parse EC numbers
                        ec_numbers = []
                        try:
                            if not pd.isna(row["EC number"]):
                                ec_strings = str(row["EC number"])
                                if "; " in ec_strings:
                                    ec_numbers = ec_strings.split("; ")
                                elif "  " in ec_strings:
                                    ec_numbers = ec_strings.split("  ")
                                else:
                                    ec_numbers = [ec_strings]
                        except Exception as e:
                            logger.error(f"Error parsing EC numbers for {row['Entry']}: {e}")
                        
                        new_protein = Protein(
                            uniprot_id=row.get("Entry"),
                            uniprot_name=row.get("Entry Name"),
                            review_status=row.get("Reviewed"),
                            name=row.get("Protein names"),
                            gene=row.get("Gene Names"),
                            organism_full=row.get("Organism"),
                            length=row.get("Length"),
                            ec_number=ec_numbers,
                            ensembl_id=str(row.get("Ensembl")).split(";") if not pd.isna(row.get("Ensembl")) else [],
                            entrez_id=row.get("Entrez"),
                            ncbi_organism=self.ncbi_organism,
                            fc=fc,
                            adj_p_value=adj_p_value,
                        )
                        self.proteins.append(new_protein)
                    except Exception as e:
                        logger.error(f"Error creating protein from UniProt: {e}")
                        continue
                
                logger.info(f"Populated proteome with {len(self.proteins)} proteins from UniProt")
            except Exception as e:
                logger.error(f"Error getting proteins from UniProt: {e}")
                logger.warning("Using empty proteome")
        
        elif ensembl:
            logger.info("Populating proteins from Ensembl - limited functionality")
            # Simplified implementation for now
            logger.warning("Ensembl protein population not fully implemented")
            self.proteins = []

    def get_interaction_partners(self):
        """
        Get protein-protein interaction partners from STRING database.
        """
        if not self.ncbi_organism or not self.proteins:
            logger.error("NCBI organism code and proteins must be set before getting interactions")
            return
            
        logger.info(f"Getting protein-protein interactions from STRING for {self.ncbi_organism}")
        
        try:
            # Extract protein IDs
            proteins = [protein.uniprot_id for protein in self.proteins]
            proteins_flattened = []
            
            # Flatten list of protein IDs if some are lists
            for protein_entries in proteins:
                if isinstance(protein_entries, list):
                    for gene_entry in protein_entries:
                        proteins_flattened.append(gene_entry)
                else:
                    proteins_flattened.append(protein_entries)
            
            # Limit to 20 proteins for testing
            proteins_limited = proteins_flattened[:20]
            
            # Map to STRING identifiers
            string_mapping_dict = string_map_identifiers(
                protein_list=proteins_limited,
                species=self.ncbi_organism,
            )
            
            # Get interactions
            string_interactions_dict = string_get_interactions(
                protein_list=list(string_mapping_dict.values()),
                species=self.ncbi_organism,
                cutoff_score=700,
            )
            
            # Translate STRING identifiers back to UniProt
            translated_interactions_dict = translate_string_dict(
                mapping_dict=string_mapping_dict,
                interactions_dict=string_interactions_dict,
            )
            
            # Assign interaction partners to proteins
            for protein in self.proteins:
                if isinstance(protein.uniprot_id, list):
                    for protein_id in protein.uniprot_id:
                        try:
                            protein.interaction_partners = translated_interactions_dict[protein_id]
                            break
                        except KeyError:
                            continue
                else:
                    try:
                        protein.interaction_partners = translated_interactions_dict[protein.uniprot_id]
                    except KeyError:
                        continue
            
            logger.info(f"Successfully mapped protein-protein interactions")
        except Exception as e:
            logger.error(f"Error getting protein-protein interactions: {e}")
            logger.warning("Continuing without protein-protein interactions")

    def get_transcription_factor_targets(
        self,
        genome_ChIP="mm10",
        distance_ChIP=5,
        cell_type_class_ChIP=None,
        cell_type_ChIP=None,
    ):
        """
        Get transcription factor targets from ChIP-Atlas.
        
        Parameters:
        -----------
        genome_ChIP : str
            Genome assembly in ChIP-Atlas
        distance_ChIP : int
            Distance from transcription start site in kb
        cell_type_class_ChIP : str
            Cell type class in ChIP-Atlas
        cell_type_ChIP : str
            Specific cell type in ChIP-Atlas
        """
        if not self.proteins:
            logger.error("Proteins must be set before getting transcription factor targets")
            return
            
        logger.info(f"Getting transcription factor targets from ChIP-Atlas for {genome_ChIP}")
        
        try:
            # Get binding data from ChIP-Atlas
            scores_unsort, gene_to_experiment = get_ChIP_data(
                genome=genome_ChIP, distance=distance_ChIP
            )
            
            # Get experiment info
            exp_info = get_ChIP_exps(
                genome=genome_ChIP,
                cell_type_class=cell_type_class_ChIP,
                cell_type=cell_type_ChIP,
            )
            
            # Assign targets to transcription factors
            for protein in self.proteins:
                if isinstance(protein.gene, list):
                    for gene in protein.gene:
                        if gene in gene_to_experiment.keys():
                            experiments = gene_to_experiment[gene]
                            for experiment in experiments:
                                if experiment in exp_info[0].to_list():
                                    # Get targets with non-zero score
                                    try:
                                        exp_col = pd.to_numeric(scores_unsort.loc[:, experiment])
                                        mask = [val != 0 for val in exp_col.to_list()]
                                        exp_col = exp_col[mask]
                                        protein.transcription_factor_targets = exp_col.index.to_list()
                                    except Exception as e:
                                        logger.error(f"Error processing ChIP-Atlas data: {e}")
            
            logger.info(f"Successfully mapped transcription factor targets")
        except Exception as e:
            logger.error(f"Error getting transcription factor targets: {e}")
            logger.warning("Continuing without transcription factor targets")

    def get_metabolites(self):
        """
        Get metabolites associated with enzymes from KEGG.
        """
        if not self.proteins:
            logger.error("Proteins must be set before getting metabolites")
            return
            
        logger.info("Getting enzyme-associated metabolites")
        
        try:
            # Get enzyme to compound mapping from KEGG
            kegg_ec_to_cpds_data = kegg_ec_to_cpds()
            kegg_compounds_list = kegg_list_compounds()
            
            # Connect metabolites to proteins based on EC numbers
            for protein in self.proteins:
                if protein.ec_number and len(protein.ec_number) > 0:
                    metabolites = []
                    for ec in protein.ec_number:
                        try:
                            # Get compounds associated with the EC number
                            if ec in kegg_ec_to_cpds_data["ec_number"].values:
                                compounds = kegg_ec_to_cpds_data[
                                    kegg_ec_to_cpds_data["ec_number"] == ec
                                ]["kegg_compounds"].tolist()
                                
                                # Add compound IDs directly
                                for compound_id in compounds:
                                    metabolites.append(compound_id)
                        except Exception as e:
                            logger.warning(f"Error getting KEGG metabolites for {ec}: {e}")
                    
                    protein.metabolites = list(set(metabolites))
            
            logger.info(f"Successfully mapped enzyme-associated metabolites")
        except Exception as e:
            logger.error(f"Error getting enzyme-associated metabolites: {e}")
            logger.warning("Continuing without enzyme-associated metabolites")
            