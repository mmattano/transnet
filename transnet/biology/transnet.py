"""Connecting transomics network
"""

__all__ = [
    "Transnet",
]

from transnet.api.kegg import *
from transnet.api.uniprot import *
from transnet.api.brenda import *
from transnet.api.string import *
from transnet.biology.elements import *
from transnet.biology.layers import *
import pandas as pd


class Transnet:
    """
    Transnet
    """

    def __init__(
            self,
            #name: str = None,
            #pathways: Pathways = None,
            #transcriptome: Transcriptome = None,
            #proteome: Proteome = None,
            #metabolome: Metabolome = None,
            #reactions: Reactions = None,
            ):
        self.name = None
        self.pathways = None
        self.transcriptome = None
        self.proteome = None
        self.metabolome = None
        self.reactions = None

    def __repr__(self):
        return f"<Transnet {self.name}>"

    def protein_metabolite_interaction(self, weight):
        interactions = []
        #######################################
        # only to reduce the network size
        new_proteome = []
        #######################################
        present_metabolites = [
            metabolite.kegg_name for metabolite in self.metabolome.metabolites
        ]
        for protein in self.proteome.proteins:
            for metabolite in protein.metabolites:
                if metabolite in present_metabolites:
                    new_interaction = [
                        protein.uniprot_id,
                        metabolite,
                        weight,
                        "Proteome",
                        "Metabolome",
                    ]
                    interactions.append(new_interaction)
        #######################################
                    new_proteome.append(protein)
        self.proteome.proteins = new_proteome
        #######################################
        return interactions

    def protein_protein_interaction(self, weight):
        interactions = []
        present_proteins = []
        for protein in self.proteome.proteins:
            if isinstance(protein.uniprot_id, list):
                present_proteins += [
                    protein_id for protein_id in protein.uniprot_id
                ]
            else:
                present_proteins.append(protein.uniprot_id)
        for protein in self.proteome.proteins:
            for interacting_protein in protein.interaction_partners:
                if interacting_protein in present_proteins:
                    new_interaction = [
                        protein.uniprot_id,
                        interacting_protein,
                        weight,
                        "Proteome",
                        "Proteome",
                    ]
                    interactions.append(new_interaction)
        return interactions

    def gene_protein_interaction(self, weight):
        interactions = []
        present_genes = [gene.ncbi_id for gene in self.transcriptome.genes]
        for i, protein in enumerate(self.proteome.proteins):
            if isinstance(protein.entrez_id, list):
                for gene in protein.entrez_id:
                    if gene in present_genes:
                        new_interaction = [
                            gene,
                            protein.uniprot_id,
                            weight,
                            "Transcriptome",
                            "Proteome",
                        ]
                        interactions.append(new_interaction)
        return interactions

    def transcription_factor_interaction(self, weight):
        interactions = []
        present_genes = [gene.name for gene in self.transcriptome.genes]
        for protein in self.proteome.proteins:
            for (
                transcription_factor_substrate
            ) in protein.transcription_factor_substrates:
                if transcription_factor_substrate in present_genes:
                    new_interaction = [
                        protein.uniprot_id,
                        transcription_factor_substrate,
                        weight,
                        "Proteome",
                        "Transcriptome",
                    ]
                    interactions.append(new_interaction)
        return interactions

    def reactions_interactions(self, weight):
        interactions = []
        present_metabolites = [
            metabolite.kegg_compound_id for metabolite in self.metabolome.metabolites
        ]
        ec_uniprot = {}
        for enzyme in self.proteome.proteins:
            for ec in enzyme.ec_number:
                ec_uniprot[ec] = enzyme.uniprot_id
        kegg_met_name = {}
        for metabolite in self.metabolome.metabolites:
            kegg_met_name[metabolite.kegg_compound_id] = metabolite.kegg_name
        for reaction in self.reactions.reactions:
            participating_metabolites = [
                metabolite for metabolite in reaction.substrates + reaction.products
            ]
            if reaction.enzyme in ec_uniprot:
                for metabolite in participating_metabolites:
                    if metabolite in present_metabolites:
                        new_interaction = [
                            ec_uniprot[reaction.enzyme],
                            reaction.id,
                            weight,
                            "Proteome",
                            "Reactions",
                        ]
                        interactions.append(new_interaction)
                        new_interaction = [
                            reaction.id,
                            kegg_met_name[metabolite],
                            weight,
                            "Reactions",
                            "Metabolome",
                        ]
                        interactions.append(new_interaction)
        return interactions

    def move_transcription_factor(self, interactions):
        present_genes = [gene.name for gene in self.transcriptome.genes]
        for protein in self.proteome.proteins:
            for (
                transcription_factor_substrate
            ) in protein.transcription_factor_substrates:
                if transcription_factor_substrate in present_genes:
                    for i, interaction in enumerate(interactions):
                        if interaction[0] == protein.uniprot_id:
                            temp_interaction = interaction
                            temp_interaction[3] = "Transcription factors"
                            interactions[i] = temp_interaction
        return interactions

    # def ptm_regulation(self):
    #    pass

    def generate_interaction_df(self):
        column_names = [
            "SourceNode",
            "TargetNode",
            "Weight",
            "SourceLayer",
            "TargetLayer",
        ]
        weight = 1
        interactions = []
        interactions += self.protein_metabolite_interaction(weight)
        interactions += self.protein_protein_interaction(weight)
        interactions += self.gene_protein_interaction(weight)
        interactions += self.transcription_factor_interaction(weight)
        interactions += self.reactions_interactions(weight)
        interactions = self.move_transcription_factor(interactions)

        interaction_df = pd.DataFrame(interactions, columns=column_names)
        interaction_df = interaction_df.drop_duplicates().reset_index(drop=True)

        return interaction_df
    
    def generate_adjacency_matrix(self):
        # This potentially will have to be changed to separate elements that are in more than one layer
        # This could be done by adding a prefix to the element name, e.g. 'm_' for metabolite, 'p_' for protein, etc.

        interaction_df = self.generate_interaction_df()
        # Create a list of all elements in the network
        interacting_elements = list(set(interaction_df['SourceNode'].tolist() + interaction_df['TargetNode'].tolist()))
        adjacency_matrix = pd.DataFrame(index=interacting_elements, columns=interacting_elements)
        adjacency_matrix = adjacency_matrix.fillna(0)

        # Fill the adjacency matrix, currently undirected
        for index, row in interaction_df.iterrows():
            adjacency_matrix.loc[row['SourceNode'], row['TargetNode']] = 1
            adjacency_matrix.loc[row['TargetNode'], row['SourceNode']] = 1

        return adjacency_matrix
