"""Biological elements, nodes in the network
"""

__all__ = [
    # "Metabolite",
    # "Reaction",
    "Protein",
    "Gene",
    "Pathway",
]

from transnet.api.kegg import *
from transnet.api.uniprot import *
from transnet.api.brenda import *
from transnet.api.string import *
from typing import Union

#class Metabolite():
#    """
#    Metabolite
#    """
#
#    def __init__(
#        self, id: str = None, name: str = None, kegg_id: str = None,
#        pathways: list = None, proteins: list = None,
#        reactions: list = None, abundance: float = None,
#        abundance_std: float = None, diff_exp: bool = None,
#        ):
#        self.id = id
#        self.name = name
#        self.kegg_id = kegg_id
#        # self.formula = str
#        # self.inchi = str
#        # self.inchikey = str
#        self.pathways = pathways
#        # self.genes = list
#        self.proteins = proteins
#        self.reactions = reactions
#        # self.interaction_metabolites = list
#        self.abundance = abundance
#        self.abundance_std = abundance_std
#        self.diff_exp = diff_exp
#
#    def __repr__(self):
#        return f"<Metabolite: {self.name}, ID: {self.id}>"
#
#
#class Reaction():
#    """
#    Reaction
#    """
#
#    def __init__(
#        self, # ID, substrates, products
#        ):
#        #self.
#
#    def __repr__(self):
#        return f"<Reaction: {self.name}, ID: {self.id}>"
#    
#    def print_function(self):
#        func = ""
#        first_negative = True
#        for i, value in enumerate(self.stoichiometry):
#            if value < 0 & first_negative:
#                func += " -> "
#                first_negative = False
#            elif value > 0 & value != len(self.stoichiometry) - 1:
#                func += " + "
#            func += str(value) + "*" + str(self.metabolites[i])


class Protein():
    """
    Protein
    """

    def __init__(
        self,
        uniprot_id: str = None,
        uniprot_name: str = None,
        review_status: str = None,
        name: str = None,
        gene: str = None,
        organism_full: str = None,
        length: str = None,
        ncbi_organism: str = None,
        ec_number : str = None,
        pathways: list = None,
        reactions: list = None,
        substrates: list = [],
        products: list = [],
        interaction_partners: list = [],
        ):
        self.uniprot_id = uniprot_id
        self.uniprot_name = uniprot_name
        self.review_status = review_status
        self.name = name
        self.gene = gene
        self.organism_full = organism_full
        self.length = length
        self.ncbi_organism = ncbi_organism
        self.ec_number = ec_number
        self.pathways = pathways
        self.reactions = reactions
        self.substrates = substrates
        self.products = products
        self.interaction_partners = interaction_partners

    def __repr__(self):
        return f"<Protein: {self.name}, ID: {self.uniprot_id}>"

    def fill(self, brenda_api : BRENDA_api = None):
        """
        Checks which info is missing and fills it
        based on what is available
        """

        if brenda_api:
            if self.ec_number is not None:
                substrates = []
                products = []
                for ec in self.ec_number.split('; '):
                    try:
                        substrates.extend(brenda_api.get_substrates(
                            ec_number=ec,
                            ))
                        products.extend(brenda_api.get_products(
                            ec_number=ec,
                            ))
                    except Exception as e:
                        print(e)
                        print(f'Could not find substrates/products for {ec}')
                self.substrates = substrates
                self.products = products
            else:
                pass
        else:
            pass


# change `fill` so that is complements it if there's alredy something in the list but more in the DB
class Gene():
    """
    Gene
    """

    def __init__(
        self,
        name: Union[str, list] = None,
        ncbi_id: Union[str, list] = None,
        kegg_id: str = None,
        related_ecs: list = [],
        description: str = None,
        kegg_organism: str = None,
    ):
        self.name = name
        self.ncbi_id = ncbi_id
        self.kegg_id = kegg_id
        self.related_ecs = related_ecs
        self.description = description
        self.kegg_organism = kegg_organism

    def __repr__(self):
        if isinstance(self.name, list):
            return f"<Gene: {self.name[0]}, ID: {self.kegg_id}>"
        else:
            return f"<Gene: {self.name}, ID: {self.kegg_id}>"
    
    def fill(self, database: dict = None):
        """
        Checks which info is missing and fills it
        based on what is available
        """

        if database:
            kegg_gene_list = database["kegg_gene_list"]
            kegg_ncbi_idtable = database["kegg_ncbi_idtable"]
            kegg_ec_link = database["kegg_ec_link"]

        # Start with KEGG ID
        if self.kegg_id is not None:
            # Check if name is there
            if self.name is None:
                kegg_info_table = kegg_gene_list#kegg_list_genes(self.kegg_organism)
                kegg_info = kegg_info_table[kegg_info_table['gene_id'] == self.kegg_id]
                self.name = kegg_info['name(s)'].to_string(index=False)
            # Check if description is there
            if self.description is None:
                kegg_info_table = kegg_gene_list#kegg_list_genes(self.kegg_organism)
                kegg_info = kegg_info_table[kegg_info_table['gene_id'] == self.kegg_id]
                self.description = kegg_info['description'].to_string(index=False)
            # Check if ncbi_id is there
            if self.ncbi_id is None:
                kegg_ncbi_conversion_table = kegg_ncbi_idtable#kegg_conv_ncbi_idtable(self.kegg_organism)
                self.ncbi_id = kegg_ncbi_conversion_table[
                    kegg_ncbi_conversion_table['kegg_id'] == self.kegg_id
                    ]['ncbi_id'].to_string(index=False)
            # Check if related EC numbers are there
            if len(self.related_ecs) == 0:
                conversion_table = kegg_ec_link#kegg_link_ec(self.kegg_organism)
                protein_ecs = []
                if self.kegg_id in conversion_table['kegg_gene_id'].to_list():
                    ecs_tmp = conversion_table[
                        conversion_table['kegg_gene_id'] == self.kegg_id
                        ]['ec_number']
                    if len(ecs_tmp) == 1:
                        protein_ecs.append(
                            ecs_tmp.to_string(index=False)
                            )
                    else:
                        for ec in ecs_tmp.to_list():
                            protein_ecs.append(ec)
                self.related_ecs = list(set(protein_ecs))


# change `fill` so that is complements it if there's alredy something in the list but more in the DBv
class Pathway():
    """
    Pathway
    """

    def __init__(
        self,
        name: str = None,
        genes: list = [],
        id: str = None,
        ecs: list = [],
        kegg_organism: str = None,
    ):
        self.name = name
        self.id = id
        self.genes = genes
        self.ecs = ecs
        self.kegg_organism = kegg_organism
    
    def __repr__(self):
        return f"<Pathway: {self.name}, ID: {self.id}>"
    
    def fill(self):
        """
        Checks which info is missing and fills it
        based on what is available
        """

        # Start with KEGG ID
        if self.id is not None:
            # Check if name is there
            if self.name is None:
                self.name = kegg_list_pathways(
                    self.kegg_organism
                    )['description'][
                        kegg_list_pathways(
                            self.kegg_organism
                            )['pathways_id'] == self.id
                            ].to_string(index=False)
            # Check if genes are there
            if len(self.genes) == 0:
                self.genes = kegg_link_pathway(
                    self.kegg_organism, self.id
                    )['kegg_gene_id'].to_list()
            # Check if EC numbers are there
            if len(self.ecs) == 0:
                conversion_table = kegg_link_ec(self.kegg_organism)
                protein_ecs = []
                for gene in kegg_link_pathway(
                    self.kegg_organism, self.id
                    )['kegg_gene_id'].to_list():
                    if gene in conversion_table['kegg_gene_id'].to_list():
                        ecs_tmp = conversion_table[
                            conversion_table['kegg_gene_id'] == gene
                            ]['ec_number']
                        if len(ecs_tmp) == 1:
                            protein_ecs.append(
                                ecs_tmp.to_string(index=False)
                                )
                        else:
                            for ec in ecs_tmp.to_list():
                                protein_ecs.append(ec)
                self.ecs = list(set(protein_ecs))
