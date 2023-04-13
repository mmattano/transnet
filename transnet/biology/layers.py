"""Biological layers (omics), layers in the network
"""

__all__ = [
    "Reactions",
    "Pathways",
    "Transcriptome",
    "Proteome",
    "Metabolome",
]

from transnet.api.kegg import *
from transnet.api.uniprot import *
from transnet.api.brenda import *
from transnet.api.string import *
from transnet.api.ensembl import *
from transnet.api.chip_atlas import *
from transnet.api.chem_info import *
from transnet.biology.elements import *
import pandas as pd


class Reactions:
    """
    Reactions
    """
    
    def __init__(
        self,
    ):
        self.reactions = []
    
    def __repr__(self):
        return f"<Reactions: {len(self.reactions)} reactions in total>"
    
    def populate(
            self,
            from_api: bool = True,
            df = None,
            ):
        if from_api:
            print("-- Populate reactions from API")
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
            print("-- Populate reactions: done")
        if not from_api and df is not None:
            print("-- Populate reactions from file")
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
            print("-- Populate reactions: done")


class Pathways:
    """
    Pathways
    """

    def __init__(
        self,
    ):
        self.kegg_organism = None
        self.pathways = []
        self.organism_full = None

    def __repr__(self):
        return f"<Pathways of: {self.kegg_organism}, {len(self.pathways)} pathways in total>"

    def populate(self,):
        print("-- Populate pathways from API")
        for _, row in kegg_list_pathways(self.kegg_organism).iterrows():
            new_pathway = Pathway(
                id = row[0],
                name = row[1],
                kegg_organism = self.kegg_organism,
            )
            self.pathways.append(new_pathway)
        print("-- Populate pathways: done")

    def fill_pathways(self):
        """
        Checks which info is missing and fills it
        based on what is available
        """
        print("-- Fill pathways")
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
                            #for ec in ecs_tmp.to_string(index=False).split("; "):
                            #    protein_ecs.append(ec)
                            if len(ecs_tmp) == 1:
                                protein_ecs.append(ecs_tmp.to_string(index=False))
                            else:
                                for ec in ecs_tmp.to_list():
                                    protein_ecs.append(ec)
                    pathway.ecs = list(set(protein_ecs))
        print("-- Fill pathways: done")


class Metabolome():
    """
    Metabolome
    """

    def __init__(self):
        self.metabolites = list
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
        self.input_data = input_data
        self.input_ids = input_data[metabolite_column_name].tolist()
        if id_type == "chebi":
            self.input_id_type = "chebi"
        elif id_type == "kegg":
            self.input_id_type = "kegg"
        elif id_type == "pubchem":
            self.input_id_type = "pubchem"
        elif id_type == "inchi":
            self.input_id_type = "inchi"
        elif id_type == "inchikey":
            self.input_id_type = "inchikey"
        elif id_type == "smile":
            self.input_id_type = "smile"

    def populate(
        self,
        metabolite_column_name: str = None,
        input_data_value_column_name: str = None,
        input_data_p_value_column_name: str = None,
    ):
        print("-- Populate metabolome, currently only works if you provide input data")
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
            #compound_id_table = pd.merge(
            #    converted_ids, chebi_kegg_conversion,
            #    left_on="chebi_ids", right_on="chebi_compounds",
            #    how="inner")
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
        # create metabolites
        for _, row in compound_id_table.iterrows():
            try:
                new_metabolite = Metabolite(
                    pubchem_id = row['pubchem_ids'],
                    kegg_name = row['kegg_names'],
                    kegg_compound_id = row['kegg_compounds'],
                    inchi = row['inchi'],
                    inchikey = row['inchikeys'],
                    chebi_id = row['chebi_ids'],
                    smile = row['smiles'],
                    fc = self.input_data[
                        self.input_data[metabolite_column_name] == row[input_id_column]
                    ][input_data_value_column_name].to_list()[0],
                    adj_p_value = self.input_data[
                        self.input_data[metabolite_column_name] == row[input_id_column]
                    ][input_data_p_value_column_name].to_list()[0],
                )
            except:
                print(row['pubchem_ids'], row['kegg_names'], row['kegg_compounds'], row['inchi'], row['inchikeys'], row['chebi_ids'], row['smiles'])
            self.metabolites.append(new_metabolite)
        print("-- Populate metabolome: done")

class Transcriptome:
    """
    Transcriptome
    """

    # def __new__(self, *args, **kwargs):
    #    return super().__new__(self)

    def __init__(self):
        self.kegg_organism = None
        self.organism_full = None
        self.input_id_type = None
        self.input_ids = []
        self.input_data = None

    def __repr__(self):
        return f"<Transcriptome of: {self.kegg_organism}, {len(self.genes)} genes in total>"

    def add_experimental_data(
        self,
        input_data: pd.DataFrame = None,
        transcript_column_name: str = None,
        id_type: str = None,
    ):
        self.input_data = input_data
        input_transcripts = input_data[transcript_column_name].tolist()
        if id_type == "ensembl":
            self.input_id_type = "ensembl"
            self.input_ids = input_transcripts

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
        self.kegg_organism = kegg_organism
        self.organism_full = organism_full
        self.genes = []
        # Currently commented out because we should probably just provide Ensembl and get KEGG info from there
        # if kegg_ftp:
        #    print('-- Populate genes from KEGG FTP')
        #    for _, row in kegg_ftp_list_genes(self.kegg_organism, ftp_base_path).iterrows():
        #        new_gene = Gene()
        #        new_gene.kegg_id = row[0]
        #        new_gene.ncbi_id = row[0]
        #        new_gene.type = row[1]
        #        new_gene.name = row[2]
        #        new_gene.description = row[3]
        #        new_gene.kegg_organism = self.kegg_organism
        #        self.genes.append(new_gene)
        # elif kegg_api:
        #    print('-- Populate genes from KEGG API')
        #    for _, row in kegg_list_genes(self.kegg_organism).iterrows():
        #        new_gene = Gene()
        #        new_gene.kegg_id = row[0]
        #        new_gene.ncbi_id = row[0]
        #        new_gene.type = row[1]
        #        new_gene.name = row[2]
        #        new_gene.description = row[3]
        #        new_gene.kegg_organism = self.kegg_organism
        #        self.genes.append(new_gene)
        if ensembl:
            ensembl_download_release(
                release_number=ensembl_release,
                organism_full=self.organism_full,
            )
            transcripts_df = ensembl_get_transcripts(
                release_number=ensembl_release,
                organism_full=self.organism_full,
            )
            if biotype:
                print(f"-- Filter transcripts by biotype: {biotype}")
                if not isinstance(biotype, list):
                    biotype = [biotype]
                # mask = [True if transcript_type in biotype else False for transcript_type in transcripts_df['biotype'].to_list()]
                # transcripts_df = transcripts_df[mask]
                transcripts_df = transcripts_df[
                    transcripts_df["biotype"].isin(biotype)
                ]
            if self.input_id_type == "ensembl":
                print("-- Populate genes from Ensembl reduced by input data")
                # mask = [True if transcript in self.input_ids else False for transcript in transcripts_df['ensembl_gene_id'].to_list()]
                # transcripts_df = transcripts_df[mask]
                transcripts_df = transcripts_df[
                    transcripts_df["ensembl_gene_id"].isin(self.input_ids)
                ]
                print(f"--- {len(transcripts_df)} genes found in input data")
            else:
                print("-- Populate genes from Ensembl (full)")
            ensembl_df = ensembl_expander(transcripts_df)
            for _, row in ensembl_df.iterrows():
                print(row)
                new_gene = Gene(
                    kegg_id = row["entrez_id"],
                    ncbi_id = row["entrez_id"],
                    ensembl_id = row["ensembl_gene_id"],
                    transcript_id = row["ensembl_transcript_id"],
                    uniprot_id = row["uniprot_id"],
                    type = row["biotype"],
                    name = row["gene_description"],
                    description = row["gene_name"],
                    kegg_organism = self.kegg_organism,
                    fc = self.input_data[
                        self.input_data[transcript_column_name] == row["ensembl_gene_id"]
                        ][input_data_value_column_name].to_list()[0],
                    adj_p_value = self.input_data[
                        self.input_data[transcript_column_name] == row["ensembl_gene_id"]
                        ][input_data_p_value_column_name].to_list()[0],
                )
                self.genes.append(new_gene)
        print("-- Populate genes: done")

    def fill_gene_info(self):
        """
        Checks which info is missing and fills it
        based on what is available
        """
        print("-- Fill genes")

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
                    gene.name = kegg_info["name(s)"].to_string(index=False)
                # Check if description is there
                if gene.description is None:
                    kegg_info = kegg_gene_list[
                        kegg_gene_list["gene_id"] == gene.kegg_id
                    ]
                    gene.description = kegg_info["description"].to_string(
                        index=False
                    )
                # Check if ncbi_id is there
                if gene.ncbi_id is None:
                    gene.ncbi_id = kegg_ncbi_idtable[
                        kegg_ncbi_idtable["kegg_id"] == gene.kegg_id
                    ]["ncbi_id"].to_string(index=False)
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
        print("-- Fill genes: done")


class Proteome:
    """
    Proteome
    """

    def __init__(self,):
        self.proteins = []
        self.organism_full = None
        self.input_id_type = None
        self.input_ids = []
        self.input_data = None

    def __repr__(self):
        return f"<Proteome of: {self.ncbi_organism}, {len(self.proteins)} proteins in total>"

    def add_experimental_data(
        self,
        input_data: pd.DataFrame = None,
        protein_column_name: str = None,
        id_type: str = None,
    ):

        self.input_data = input_data
        input_poteins = input_data[protein_column_name].tolist()
        if id_type == "uniprot":
            self.input_id_type = "uniprot"
            self.input_ids = input_poteins

    def populate(
        self,
        kegg_organism: str = None,
        ncbi_organism: str = None,
        organism_full: str = None,
        ensembl: bool = False,  # print a message if both are false
        uniprot: bool = False,
        ensembl_release: int = 109,
        protein_column_name: str = None,
        input_data_value_column_name: str = None,
        input_data_p_value_column_name: str = None,
    ):

        self.kegg_organism = kegg_organism
        self.ncbi_organism = ncbi_organism
        self.organism_full = organism_full
        if ensembl:
            print("-- Populate proteins from Ensembl, currently commented out")
        #    ensembl_download_release(release_number=ensembl_release, organism_full=self.organism_full)
        #    transcripts_df = ensembl_get_transcripts(release_number=ensembl_release, organism_full=self.organism_full)
        #    ensembl_df = ensembl_expander(transcripts_df)
        #    # and then set up a similar loop as for uniprot
        #    for _, row in ensembl_df.iterrows():
        #        new_protein = Protein()
        #        new_protein.ensembl_id = row['ensembl_gene_id']
        #        new_protein.uniprot_id = row['uniprot_id']
        #        new_protein.gene = row['ncbi_gene_id']
        #        new_protein.ncbi_organism = self.ncbi_organism
        #        self.proteins.append(new_protein)
        elif uniprot:
            if self.input_id_type == "uniprot":
                print(
                    "-- Populate proteins from UniProt reduced by input data"
                )
            else:
                print("-- Populate proteins from UniProt (full)")
            uniprot_proteins = uniprot_list_proteins(self.ncbi_organism)
            uniprot_proteins = uniprot_add_entrez_id(uniprot_proteins)
            for _, row in uniprot_proteins.iterrows():
                if self.input_id_type == "uniprot":
                    if row["Entry"] not in self.input_ids:
                        continue
                try:
                    print(row["EC number"])
                    print(str(row["EC number"]).split("; "))
                except:
                    print("No EC number")
                new_protein = Protein(
                    uniprot_id = row["Entry"],
                    uniprot_name = row["Entry Name"],
                    review_status = row["Reviewed"],
                    name = row["Protein names"],
                    gene = row["Gene Names"],
                    organism_full = row["Organism"],
                    length = row["Length"],
                    ec_number = str(row["EC number"]).split("; "),
                    ensembl_id = row["Ensembl"],
                    entrez_id = row["Entrez"],
                    ncbi_organism = self.ncbi_organism,
                    fc = self.input_data[
                        self.input_data[protein_column_name] == row["Entry"]
                        ][input_data_value_column_name].to_list()[0],
                    adj_p_value = self.input_data[
                        self.input_data[protein_column_name] == row["Entry"]
                        ][input_data_p_value_column_name].to_list()[0],
                )
                self.proteins.append(new_protein)
        print("-- Populate proteins: done")

    def get_interaction_partners(self):
        proteins = [protein.uniprot_id for protein in self.proteins]
        proteins_flattened = []
        for protein_entries in proteins:
            if isinstance(protein_entries, list):
                for gene_entry in protein_entries:
                    proteins_flattened.append(gene_entry)
            else:
                proteins_flattened.append(protein_entries)
        # limit it to 1000 genes at a time
        start_slice = 0
        string_mapping_dict = {}
        while start_slice < len(proteins_flattened):
            end_slice = start_slice + 1000
            temp = string_map_identifiers(
                protein_list=proteins_flattened[start_slice:end_slice],
                species=self.ncbi_organism,
            )
            string_mapping_dict.update(temp)
            start_slice = end_slice
        # limit it to 1000 genes at a time
        start_slice = 0
        string_interactions_dict = {}
        while start_slice < len(list(string_mapping_dict.values())):
            end_slice = start_slice + 1000
            temp = string_get_interactions(
                protein_list=list(string_mapping_dict.values())[
                    start_slice:end_slice
                ],
                species=self.ncbi_organism,
                cutoff_score=700,
            )
            string_interactions_dict.update(temp)
            start_slice = end_slice
        translated_interactions_dict = translate_string_dict(
            mapping_dict=string_mapping_dict,
            interactions_dict=string_interactions_dict,
        )
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

    def get_transcription_factor_targets(
        self,
        genome_ChIP="mm10",
        distance_ChIP=5,
        cell_type_class_ChIP="Adipocyte",
        cell_type_ChIP="Brown preadipocytes",
    ):
        # access info
        scores_unsort, gene_to_experiment = get_ChIP_data(
            genome=genome_ChIP, distance=distance_ChIP
        )
        exp_info = get_ChIP_exps(
            genome=genome_ChIP,
            cell_type_class=cell_type_class_ChIP,
            cell_type=cell_type_ChIP,
        )
        for protein in self.proteins:
            if isinstance(protein.gene, list):
                for gene in protein.gene:
                    if gene in gene_to_experiment.keys():
                        experiments = gene_to_experiment[gene]
                        for experiment in experiments:
                            if experiment in exp_info[0].to_list():
                                exp_col = pd.to_numeric(
                                    scores_unsort.loc[:, experiment]
                                )
                                mask = [val != 0 for val in exp_col.to_list()]
                                exp_col = exp_col[mask]
                                protein.transcription_factor_targets = (
                                    exp_col.index.to_list()
                                )

    def get_metabolites(self, brenda_api: BRENDA_api = None):
        """
        Checks which info is missing and fills it
        based on what is available
        """
        print("-- Fill proteins")

        database = {
            "kegg_ec_to_cpds": kegg_ec_to_cpds(),
            "kegg_compounds_list": kegg_list_compounds(),
        }

        if brenda_api:
            print("--- Fill proteins: fill metabolites and brenda info")
            for protein in self.proteins:
                if protein.ec_number is not None:
                    substrates = []
                    products = []
                    for ec in protein.ec_number.split("; "):
                        try:
                            substrates.extend(
                                brenda_api.get_substrates(ec_number=ec,)
                            )
                            products.extend(
                                brenda_api.get_products(ec_number=ec,)
                            )
                        except Exception as e:
                            print(e)
                            print(
                                f"Could not find substrates/products for {ec}"
                            )
                    protein.substrates = substrates
                    protein.products = products
                else:
                    pass
        else:
            for protein in self.proteins:
                if database:
                    kegg_ec_to_cpds_list = database["kegg_ec_to_cpds"]
                    kegg_compounds_list = database["kegg_compounds_list"]
                    if protein.ec_number is not None:
                        compounds = []
                        for ec in protein.ec_number:
                            try:
                                compounds.extend(
                                    kegg_ec_to_cpds_list[
                                        kegg_ec_to_cpds_list["ec_number"] == ec
                                    ]["kegg_compounds"].tolist()
                                )
                            except Exception as e:
                                print(e)
                                print(
                                    f"Could not find substrates/products for {ec}"
                                )
                        metabolites = [
                            kegg_compounds_list[
                                kegg_compounds_list["kegg_compounds"]
                                == compound
                            ]["name"].to_string(index=False)
                            for compound in compounds
                        ]
                        protein.metabolites = list(set(metabolites))

        print("-- Fill proteins: done")
