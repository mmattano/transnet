"""Biological layers (omics), layers in the network
"""

from typing import List, Dict, Any, Optional, Union
import pandas as pd
from ..api.kegg import *
from ..api.uniprot import *
from ..api.brenda import *
from ..api.string import *
from ..api.ensembl import *
from ..api.chip_atlas import *
from ..api.chem_info import *
from .elements import *
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
                    pathway.name = kegg_pathway_list[
                        "description"
                    ][
                        kegg_pathway_list["pathways_id"]
                        == pathway.id
                    ].to_string(
                        index=False
                    )
                # Check if genes are there
                if len(pathway.genes) == 0:
                    pathway.genes = kegg_pathway_links[
                        kegg_pathway_links['pathway'] == pathway.id][
                        "kegg_gene_id"
                    ].to_list()
                # Check if EC numbers are there
                if len(pathway.ecs) == 0:
                    protein_ecs = []
                    for gene in kegg_pathway_links[
                        kegg_pathway_links['pathway'] == pathway.id][
                        "kegg_gene_id"
                    ].to_list():
                        if gene in kegg_pathway_ecs["kegg_gene_id"].to_list():
                            ecs_tmp = kegg_pathway_ecs[
                                kegg_pathway_ecs["kegg_gene_id"] == gene
                            ]["ec_number"]
                            if len(ecs_tmp) == 1:
                                protein_ecs.append(ecs_tmp.to_string(index=False))
                            else:
                                for ec in ecs_tmp.to_list():
                                    protein_ecs.append(ec)
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
        # Currently only works if you provide input data
        if self.input_id_type == "kegg":
            kegg_compounds_df = kegg_list_compounds()
            kegg_compounds_df = kegg_compounds_df[
                kegg_compounds_df["compound_id"].isin(self.input_ids)
            ]
            kegg_converted = kegg_to_chebi(kegg_compounds_df.kegg_compounds)
            converted_ids = chemical_info_converter(
                kegg_converted.chebi_compounds.to_list()
                )
            compound_id_table = pd.merge(
                pd.merge(
                kegg_compounds_df, kegg_converted, on="kegg_compounds"
                ),
                converted_ids, 
                left_on="chebi_compounds", right_on="chebi_ids",
                how="right")
            input_id_column = 'kegg_compounds'
        elif self.input_id_type:
            converted_ids = chemical_info_converter(self.input_ids)
            chebi_kegg_conversion = chebi_to_kegg(converted_ids.chebi_ids.to_list())
            kegg_all_compounds = kegg_list_compounds()
            kegg_names = []
            for kegg_compound in chebi_kegg_conversion.kegg_compounds.to_list():
                if kegg_compound is None:
                    kegg_name = None
                else:
                    try:
                        kegg_name = kegg_all_compounds[
                            kegg_all_compounds.kegg_compounds == kegg_compound
                        ]['name'].to_list()[0]
                    except IndexError:
                        kegg_name = None
                kegg_names.append(kegg_name)
            
            compound_id_table = converted_ids
            compound_id_table['kegg_compounds'] = chebi_kegg_conversion.kegg_compounds
            compound_id_table['kegg_names'] = kegg_names
            if self.input_id_type == "chebi":
                input_id_column = 'chebi_ids'            
            elif self.input_id_type == "pubchem":
                input_id_column = 'pubchem_ids'
            elif self.input_id_type == "inchi":
                input_id_column = 'inchi'
            elif self.input_id_type == "inchikey":
                input_id_column = 'inchikeys'
            elif self.input_id_type == "smile":
                input_id_column = 'smiles'
        else:
            logger.warning("No ID type specified, cannot populate metabolome")
            return
        
        # Create metabolite objects
        for _, row in compound_id_table.iterrows():
            try:
                fc = None
                adj_p_value = None
                
                if self.input_data is not None and metabolite_column_name and input_id_column in row:
                    matches = self.input_data[
                        self.input_data[metabolite_column_name] == row[input_id_column]
                    ]
                    if len(matches) > 0 and input_data_value_column_name:
                        fc = matches[input_data_value_column_name].iloc[0]
                    if len(matches) > 0 and input_data_p_value_column_name:
                        adj_p_value = matches[input_data_p_value_column_name].iloc[0]
                
                new_metabolite = Metabolite(
                    pubchem_id = row.get('pubchem_ids'),
                    kegg_name = row.get('kegg_names'),
                    kegg_compound_id = row.get('kegg_compounds'),
                    inchi = row.get('inchi'),
                    inchikey = row.get('inchikeys'),
                    chebi_id = row.get('chebi_ids'),
                    smile = row.get('smiles'),
                    fc = fc,
                    adj_p_value = adj_p_value,
                )
                self.metabolites.append(new_metabolite)
            except Exception as e:
                logger.error(f"Error creating metabolite from row: {row}")
                logger.error(f"Error: {e}")
        
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
        self.kegg_organism = kegg_organism
        self.organism_full = organism_full
        self.genes = []
        
        if ensembl:
            logger.info(f"Populating genes from Ensembl (release {ensembl_release})")
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
            if self.input_id_type == "ensembl":
                logger.info("Filtering genes by input data")
                transcripts_df = transcripts_df[
                    transcripts_df["ensembl_gene_id"].isin(self.input_ids)
                ]
                logger.info(f"Found {len(transcripts_df)} genes in input data")
            
            ensembl_df = ensembl_expander(transcripts_df)
            
            for _, row in ensembl_df.iterrows():
                fc = None
                adj_p_value = None
                
                if self.input_id_type == "ensembl" and self.input_data is not None:
                    matches = self.input_data[
                        self.input_data[transcript_column_name] == row["ensembl_gene_id"]
                    ]
                    if len(matches) > 0 and input_data_value_column_name:
                        fc = matches[input_data_value_column_name].iloc[0]
                    if len(matches) > 0 and input_data_p_value_column_name:
                        adj_p_value = matches[input_data_p_value_column_name].iloc[0]
                
                new_gene = Gene(
                    kegg_id = row.get("entrez_id"),
                    ncbi_id = row.get("entrez_id"),
                    ensembl_id = row.get("ensembl_gene_id"),
                    transcript_id = row.get("ensembl_transcript_id"),
                    uniprot_id = row.get("uniprot_id"),
                    type = row.get("biotype"),
                    name = row.get("gene_description"),
                    description = row.get("gene_name"),
                    kegg_organism = self.kegg_organism,
                    fc = fc,
                    adj_p_value = adj_p_value,
                )
                self.genes.append(new_gene)
        
        logger.info(f"Populated transcriptome with {len(self.genes)} genes")

    def fill_gene_info(self):
        """
        Fill in additional gene information from KEGG.
        """
        if not self.kegg_organism:
            logger.error("KEGG organism code must be set before filling gene info")
            return
            
        logger.info(f"Filling gene information for {self.kegg_organism}")

        kegg_gene_list = kegg_list_genes(self.kegg_organism)
        kegg_ncbi_idtable = kegg_conv_ncbi_idtable(self.kegg_organism)
        kegg_ec_link = kegg_link_ec(self.kegg_organism)

        for gene in self.genes:
            # Start with KEGG ID
            if gene.kegg_id is not None:
                # Check if name is there
                if gene.name is None:
                    kegg_info = kegg_gene_list[
                        kegg_gene_list["gene_id"] == gene.kegg_id
                    ]
                    if not kegg_info.empty:
                        gene.name = kegg_info["name(s)"].to_string(index=False)
                
                # Check if description is there
                if gene.description is None:
                    kegg_info = kegg_gene_list[
                        kegg_gene_list["gene_id"] == gene.kegg_id
                    ]
                    if not kegg_info.empty:
                        gene.description = kegg_info["description"].to_string(
                            index=False
                        )
                
                # Check if ncbi_id is there
                if gene.ncbi_id is None:
                    ncbi_info = kegg_ncbi_idtable[
                        kegg_ncbi_idtable["kegg_id"] == gene.kegg_id
                    ]
                    if not ncbi_info.empty:
                        gene.ncbi_id = ncbi_info["ncbi_id"].to_string(index=False)
                
                # Check if related EC numbers are there
                if len(gene.related_ecs) == 0:
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
        
        logger.info(f"Successfully filled gene information")


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
        uniprot: bool = False,
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
        self.kegg_organism = kegg_organism
        self.ncbi_organism = ncbi_organism
        self.organism_full = organism_full
        
        if not ensembl and not uniprot:
            logger.warning("Neither Ensembl nor UniProt selected as data source")
            
        if ensembl:
            logger.info("Populating proteins from Ensembl (currently not implemented)")
            # Implementation for Ensembl population would go here
            
        elif uniprot:
            if self.input_id_type == "uniprot":
                logger.info("Populating proteins from UniProt filtered by input data")
            else:
                logger.info("Populating all proteins from UniProt")
                
            uniprot_proteins = uniprot_list_proteins(self.ncbi_organism)
            uniprot_proteins = uniprot_add_entrez_id(uniprot_proteins)
            
            for _, row in uniprot_proteins.iterrows():
                if self.input_id_type == "uniprot" and row["Entry"] not in self.input_ids:
                    continue
                
                fc = None
                adj_p_value = None
                
                if self.input_id_type == "uniprot" and self.input_data is not None:
                    matches = self.input_data[
                        self.input_data[protein_column_name] == row["Entry"]
                    ]
                    if len(matches) > 0 and input_data_value_column_name:
                        fc = matches[input_data_value_column_name].iloc[0]
                    if len(matches) > 0 and input_data_p_value_column_name:
                        adj_p_value = matches[input_data_p_value_column_name].iloc[0]
                
                # Parse EC numbers
                ec_numbers = []
                try:
                    if not pd.isna(row["EC number"]):
                        ec_numbers = str(row["EC number"]).split("; ")
                except Exception as e:
                    logger.error(f"Error parsing EC numbers for {row['Entry']}: {e}")
                
                new_protein = Protein(
                    uniprot_id = row.get("Entry"),
                    uniprot_name = row.get("Entry Name"),
                    review_status = row.get("Reviewed"),
                    name = row.get("Protein names"),
                    gene = row.get("Gene Names"),
                    organism_full = row.get("Organism"),
                    length = row.get("Length"),
                    ec_number = ec_numbers,
                    ensembl_id = row.get("Ensembl"),
                    entrez_id = row.get("Entrez"),
                    ncbi_organism = self.ncbi_organism,
                    fc = fc,
                    adj_p_value = adj_p_value,
                )
                self.proteins.append(new_protein)
        
        logger.info(f"Populated proteome with {len(self.proteins)} proteins")

    def get_interaction_partners(self):
        """
        Get protein-protein interaction partners from STRING database.
        """
        if not self.ncbi_organism:
            logger.error("NCBI organism code must be set before getting interactions")
            return
            
        logger.info(f"Getting protein-protein interactions from STRING for {self.ncbi_organism}")
        
        proteins = [protein.uniprot_id for protein in self.proteins]
        # Flatten list of protein IDs if some are lists
        proteins_flattened = []
        for protein_entries in proteins:
            if isinstance(protein_entries, list):
                for gene_entry in protein_entries:
                    proteins_flattened.append(gene_entry)
            else:
                proteins_flattened.append(protein_entries)
        
        # Limit it to 1000 proteins at a time
        start_slice = 0
        string_mapping_dict = {}
        while start_slice < len(proteins_flattened):
            end_slice = min(start_slice + 1000, len(proteins_flattened))
            try:
                temp = string_map_identifiers(
                    protein_list=proteins_flattened[start_slice:end_slice],
                    species=self.ncbi_organism,
                )
                string_mapping_dict.update(temp)
            except Exception as e:
                logger.error(f"Error mapping STRING identifiers: {e}")
            start_slice = end_slice
        
        # Get interactions in batches
        start_slice = 0
        string_interactions_dict = {}
        while start_slice < len(list(string_mapping_dict.values())):
            end_slice = min(start_slice + 1000, len(list(string_mapping_dict.values())))
            try:
                temp = string_get_interactions(
                    protein_list=list(string_mapping_dict.values())[
                        start_slice:end_slice
                    ],
                    species=self.ncbi_organism,
                    cutoff_score=700,
                )
                string_interactions_dict.update(temp)
            except Exception as e:
                logger.error(f"Error getting STRING interactions: {e}")
            start_slice = end_slice
        
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
                        protein.interaction_partners = translated_interactions_dict[
                            protein_id
                        ]
                        break
                    except KeyError:
                        continue
            else:
                try:
                    protein.interaction_partners = translated_interactions_dict[
                        protein.uniprot_id
                    ]
                except KeyError:
                    continue
        
        logger.info(f"Successfully mapped protein-protein interactions")

    def get_transcription_factor_targets(
        self,
        genome_ChIP="mm10",
        distance_ChIP=5,
        cell_type_class_ChIP="Adipocyte",
        cell_type_ChIP="Brown preadipocytes",
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
        logger.info(f"Getting transcription factor targets from ChIP-Atlas for {genome_ChIP}")
        
        # Get binding data from ChIP-Atlas
        scores_unsort, gene_to_experiment = get_ChIP_data(
            genome=genome_ChIP, distance=distance_ChIP
        )
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
                                exp_col = pd.to_numeric(
                                    scores_unsort.loc[:, experiment]
                                )
                                mask = [val != 0 for val in exp_col.to_list()]
                                exp_col = exp_col[mask]
                                protein.transcription_factor_targets = (
                                    exp_col.index.to_list()
                                )
        
        logger.info(f"Successfully mapped transcription factor targets")

    def get_activators_and_inhibitors(self, brenda_api=None):
        """
        Get activators and inhibitors of enzymes from BRENDA.
        
        Parameters:
        -----------
        brenda_api : BRENDA_api
            BRENDA API object for querying
        """
        if not brenda_api:
            logger.error("BRENDA API object must be provided")
            return
            
        logger.info("Getting activators and inhibitors from BRENDA")
        
        # For each protein with EC numbers, get activators and inhibitors
        for protein in self.proteins:
            if protein.ec_number and len(protein.ec_number) > 0:
                activators = []
                inhibitors = []
                
                for ec in protein.ec_number:
                    try:
                        # Get activators
                        act = brenda_api.get_activators(ec_number=ec)
                        if act:
                            activators.extend([(a, None) for a in act])
                        
                        # Get inhibitors
                        inh = brenda_api.get_inhibitors(ec_number=ec)
                        if inh:
                            inhibitors.extend([(i, None) for i in inh])
                    except Exception as e:
                        logger.warning(f"Error getting BRENDA data for {ec}: {e}")
                
                protein.activators = activators
                protein.inhibitors = inhibitors
        
        logger.info(f"Successfully mapped enzyme activators and inhibitors")

    def get_metabolites(self, brenda_api=None):
        """
        Get metabolites associated with enzymes from KEGG and BRENDA.
        
        Parameters:
        -----------
        brenda_api : BRENDA_api
            BRENDA API object for querying
        """
        logger.info("Getting enzyme-associated metabolites")
        
        # Get enzyme to compound mapping from KEGG
        database = {
            "kegg_ec_to_cpds": kegg_ec_to_cpds(),
            "kegg_compounds_list": kegg_list_compounds(),
        }
        
        # If BRENDA API is available, use it for more detailed substrate/product info
        if brenda_api:
            logger.info("Using BRENDA for detailed substrate/product information")
            for protein in self.proteins:
                if protein.ec_number and len(protein.ec_number) > 0:
                    substrates = []
                    products = []
                    
                    for ec in protein.ec_number:
                        try:
                            # Get substrates
                            subs = brenda_api.get_substrates(ec_number=ec)
                            if subs:
                                substrates.extend(subs)
                            
                            # Get products
                            prods = brenda_api.get_products(ec_number=ec)
                            if prods:
                                products.extend(prods)
                        except Exception as e:
                            logger.warning(f"Error getting BRENDA data for {ec}: {e}")
                    
                    protein.substrates = substrates
                    protein.products = products
        else:
            # Use KEGG for basic metabolite associations
            for protein in self.proteins:
                if protein.ec_number and len(protein.ec_number) > 0:
                    kegg_ec_to_cpds_list = database["kegg_ec_to_cpds"]
                    kegg_compounds_list = database["kegg_compounds_list"]
                    
                    metabolites = []
                    for ec in protein.ec_number:
                        try:
                            # Get compounds associated with the EC number
                            compounds = kegg_ec_to_cpds_list[
                                kegg_ec_to_cpds_list["ec_number"] == ec
                            ]["kegg_compounds"].tolist()
                            
                            # Get names for the compounds
                            for compound in compounds:
                                name = kegg_compounds_list[
                                    kegg_compounds_list["kegg_compounds"] == compound
                                ]["name"].to_string(index=False)
                                if name:
                                    metabolites.append(name)
                        except Exception as e:
                            logger.warning(f"Error getting KEGG data for {ec}: {e}")
                    
                    protein.metabolites = list(set(metabolites))
        
        logger.info(f"Successfully mapped enzyme-associated metabolites")
        