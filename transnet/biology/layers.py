"""Biological layers (omics), layers in the network
"""

__all__ = [
    "Pathways",
    "Transcriptome",
    "Proteome",
    # "Fluxome",
    # "Metabolome",
    # "Genome",
]

from transnet.api.kegg import *
from transnet.api.uniprot import *
from transnet.api.brenda import *
from transnet.api.string import *
from transnet.biology.elements import *
import pandas as pd

class Pathways():
    """
    Pathways
    """

    def __init__(
        self,
        kegg_organism: str = None,
        names: list = [],
        ids: list = [],
        pathways: list = [],
        ):
        self.kegg_organism = kegg_organism
        self.names = names
        self.ids = ids
        self.pathways = pathways
    
    def __repr__(self):
        return f"<Pathways of: {self.kegg_organism}, {len(self.ids)} pathways in total>"
    
    def populate(
        self,
        ):
        print('-- Populate pathways')
        for _, row in kegg_list_pathways(self.kegg_organism).iterrows():
            self.ids.append(row[0])
            self.names.append(row[1])
            self.pathways.append(
                Pathway(
                    id = row[0],
                    name = row[1],
                    kegg_organism = self.kegg_organism,
                    )
                )
        print('-- Populate pathways: done')
    
    def fill_pathways(self):
        """
        Checks which info is missing and fills it
        based on what is available
        """
        print('-- Fill pathways')
        for pathway in self.pathways:
            pathway.fill()
        print('-- Fill pathways: done')


class Transcriptome():
    """
    Transcriptome
    """

    def __init__(
        self,
        kegg_organism: str = None,
        ids: list = [],
        names: list = [],
        descriptions: list = [],
        genes: list = [],
        ):
        self.kegg_organism = kegg_organism
        self.ids = ids
        self.names = names
        self.descriptions = descriptions
        self.genes = genes

    def __repr__(self):
        return f"<Transcriptome of: {self.kegg_organism}, {len(self.ids)} genes in total>"

    def populate(
        self,
        ):
        print('-- Populate genes')
        for _, row in kegg_list_genes(self.kegg_organism).iterrows():
            self.ids.append(row[0])
            self.names.append(row[1])
            self.descriptions.append(row[2])
            self.genes.append(
                Gene(
                    kegg_id = row[0],
                    name = row[1],
                    description = row[2],
                    kegg_organism = self.kegg_organism,
                    )
                )
        print('-- Populate genes: done')
    
    def fill_genes(self):
        """
        Checks which info is missing and fills it
        based on what is available
        """
        print('-- Fill genes')
        database = {
            "kegg_gene_list": kegg_list_genes(self.kegg_organism),
            "kegg_ncbi_idtable": kegg_conv_ncbi_idtable(self.kegg_organism),
            "kegg_ec_link": kegg_link_ec(self.kegg_organism),
        }
        for gene in self.genes:
            gene.fill(database)
        print('-- Fill genes: done')


class Proteome():
    """
    Proteome
    """

    def __init__(
        self,
        ncbi_organism: str = None,
        uniprot_ids: list = [],
        proteins: list = [],
        brenda_api: BRENDA_api = None,
        ):
        self.ncbi_organism = ncbi_organism
        self.uniprot_ids = uniprot_ids
        self.proteins = proteins
        self.brenda_api = brenda_api

    def __repr__(self):
        return f"<Proteome of: {self.ncbi_organism}, {len(self.uniprot_ids)} proteins in total>"

    def populate(self):
        print('-- Populate proteins')
        for _, row in uniprot_list_proteins(self.ncbi_organism).iterrows():
            self.uniprot_ids.append(row[0])
            if not pd.isna(row[7]):
                _ec_number = row[7]
            else:
                _ec_number = None
            self.proteins.append(
                Protein(
                    uniprot_id = row[0],
                    uniprot_name = row[1],
                    review_status = row[2],
                    name = row[3],
                    gene = row[4],
                    organism_full = row[5],
                    length = row[6],
                    ec_number = _ec_number,
                    ncbi_organism = self.ncbi_organism,
                    )
                )
        print('-- Populate proteins: done')
    
    def get_interaction_partners(self):
        genes = [protein.gene for protein in self.proteins]
        
        genes_flattened = []
        for gene_entries in genes:
            if isinstance(gene_entries, list):
                for gene_entry in gene_entries:
                    genes_flattened.append(gene_entry)
        
        # limit it to 1000 genes at a time
        start_slice = 0
        string_mapping_dict = {}
        while start_slice < len(genes_flattened):
            end_slice = start_slice + 1000
            temp = string_map_identifiers(
                protein_list=genes_flattened[start_slice:end_slice],
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
                protein_list=list(string_mapping_dict.values())[start_slice:end_slice],
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
            if isinstance(protein.gene, list):
                for gene in protein.gene:
                    try:
                        protein.interaction_partners = translated_interactions_dict[gene]
                        #print(translated_interactions_dict[gene])
                        break
                    except KeyError:
                        continue

    def fill_proteins(self):
        """
        Checks which info is missing and fills it
        based on what is available
        """
        print('-- Fill proteins')
        print('--- Fill proteins: get interaction partners')
        self.get_interaction_partners()
        
        print('--- Fill proteins: fill brenda')
        for protein in self.proteins:
            protein.fill(brenda_api=self.brenda_api)
        print('-- Fill proteins: done')


#class Fluxome():
#    """
#    Fluxome
#    """
#
#    def __init__(self, reactions: list = None, fluxes: dict = None):
#        self.reactions = reactions
#        self.fluxes = fluxes
#
#    def __repr__(self):
#        return f"<Fluxomics>"
#
#
#class Metabolome():
#    """
#    Metabolome
#    """
#
#    def __init__(self):
#        self.metabolites = list
#        self.abundances = float
#
#    def __repr__(self):
#        return f"<Metabolome>"
#
#
#class Genome():
#    """
#    Genome
#    """
#
#    def __init__(self):
#        self.genes = list
#
#    def __repr__(self):
#        return f"<Genome>"