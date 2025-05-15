"""Connecting transomics network
"""

from typing import List, Dict, Any, Optional, Union
import pandas as pd
import networkx as nx
import logging
from ..api.kegg import *
from ..api.uniprot import *
from ..api.brenda import *
from ..api.string import *
from .elements import *
from .layers import *

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

__all__ = [
    "Transnet",
]


class Transnet:
    """
    Transnet class for connecting transomics networks.
    """

    def __init__(
            self,
            name: str = "Transnet",
            pathways: Pathways = None,
            transcriptome: Transcriptome = None,
            proteome: Proteome = None,
            metabolome: Metabolome = None,
            reactions: Reactions = None,
            ):
        self.name = name
        self.pathways = pathways
        self.transcriptome = transcriptome
        self.proteome = proteome
        self.metabolome = metabolome
        self.reactions = reactions
        self.graph = None

    def __repr__(self):
        return f"<Transnet {self.name}>"

    def protein_metabolite_interaction(self, weight: float = 1.0):
        """
        Create protein-metabolite interactions.
        
        Parameters:
        -----------
        weight : float
            Weight of the interaction edges
            
        Returns:
        --------
        List of interactions [source, target, weight, source_layer, target_layer]
        """
        interactions = []
        
        # Check if required layers exist
        if not self.proteome or not self.metabolome:
            logger.warning("Proteome or Metabolome layer missing")
            return interactions
            
        # Get metabolites in the network
        present_metabolites = [
            metabolite.kegg_compound_id for metabolite in self.metabolome.metabolites
        ]
        
        # Create interactions
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
        
        logger.info(f"Created {len(interactions)} protein-metabolite interactions")
        return interactions

    def protein_protein_interaction(self, weight: float = 1.0):
        """
        Create protein-protein interactions.
        
        Parameters:
        -----------
        weight : float
            Weight of the interaction edges
            
        Returns:
        --------
        List of interactions [source, target, weight, source_layer, target_layer]
        """
        interactions = []
        
        # Check if required layer exists
        if not self.proteome:
            logger.warning("Proteome layer missing")
            return interactions
            
        # Get proteins in the network
        present_proteins = []
        for protein in self.proteome.proteins:
            if isinstance(protein.uniprot_id, list):
                present_proteins += [
                    protein_id for protein_id in protein.uniprot_id
                ]
            else:
                present_proteins.append(protein.uniprot_id)
        
        # Create interactions
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
        
        logger.info(f"Created {len(interactions)} protein-protein interactions")
        return interactions

    def gene_protein_interaction(self, weight: float = 1.0):
        """
        Create gene-protein interactions.
        
        Parameters:
        -----------
        weight : float
            Weight of the interaction edges
            
        Returns:
        --------
        List of interactions [source, target, weight, source_layer, target_layer]
        """
        interactions = []
        
        # Check if required layers exist
        if not self.transcriptome or not self.proteome:
            logger.warning("Transcriptome or Proteome layer missing")
            return interactions
            
        # Get genes in the network
        present_genes = [gene.ncbi_id for gene in self.transcriptome.genes]
        
        # Create interactions
        for protein in self.proteome.proteins:
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
        
        logger.info(f"Created {len(interactions)} gene-protein interactions")
        return interactions

    def transcription_factor_interaction(self, weight: float = 1.0):
        """
        Create transcription factor-gene interactions.
        
        Parameters:
        -----------
        weight : float
            Weight of the interaction edges
            
        Returns:
        --------
        List of interactions [source, target, weight, source_layer, target_layer]
        """
        interactions = []
        
        # Check if required layers exist
        if not self.transcriptome or not self.proteome:
            logger.warning("Transcriptome or Proteome layer missing")
            return interactions
            
        # Get genes in the network
        present_genes = [gene.name for gene in self.transcriptome.genes]
        
        # Create interactions
        for protein in self.proteome.proteins:
            for target in protein.transcription_factor_targets:
                if target in present_genes:
                    new_interaction = [
                        protein.uniprot_id,
                        target,
                        weight,
                        "Proteome",
                        "Transcriptome",
                    ]
                    interactions.append(new_interaction)
        
        logger.info(f"Created {len(interactions)} transcription factor interactions")
        return interactions

    def enzyme_reaction_interaction(self, weight: float = 1.0):
        """
        Create enzyme-reaction-metabolite interactions.
        
        Parameters:
        -----------
        weight : float
            Weight of the interaction edges
            
        Returns:
        --------
        List of interactions [source, target, weight, source_layer, target_layer]
        """
        interactions = []
        
        # Check if required layers exist
        if not self.proteome or not self.metabolome or not self.reactions:
            logger.warning("Proteome, Metabolome, or Reactions layer missing")
            return interactions
            
        # Get metabolites in the network
        present_metabolites = [
            metabolite.kegg_compound_id for metabolite in self.metabolome.metabolites
        ]
        
        # Create EC number to UniProt mapping
        ec_uniprot = {}
        for enzyme in self.proteome.proteins:
            for ec in enzyme.ec_number:
                ec_uniprot[ec] = enzyme.uniprot_id
        
        # Create metabolite ID to name mapping
        kegg_met_name = {}
        for metabolite in self.metabolome.metabolites:
            kegg_met_name[metabolite.kegg_compound_id] = metabolite.kegg_name
        
        # Create interactions
        for reaction in self.reactions.reactions:
            participating_metabolites = [
                metabolite for metabolite in reaction.substrates + reaction.products
            ]
            
            if reaction.enzyme in ec_uniprot:
                for metabolite in participating_metabolites:
                    if metabolite in present_metabolites:
                        # Enzyme to reaction interaction
                        new_interaction = [
                            ec_uniprot[reaction.enzyme],
                            reaction.id,
                            weight,
                            "Proteome",
                            "Reactions",
                        ]
                        interactions.append(new_interaction)
                        
                        # Reaction to metabolite interaction
                        new_interaction = [
                            reaction.id,
                            kegg_met_name[metabolite],
                            weight,
                            "Reactions",
                            "Metabolome",
                        ]
                        interactions.append(new_interaction)
        
        logger.info(f"Created {len(interactions)} enzyme-reaction-metabolite interactions")
        return interactions

    def activator_inhibitor_interactions(self, weight: float = 1.0):
        """
        Create enzyme activator/inhibitor interactions.
        
        Parameters:
        -----------
        weight : float
            Weight of the interaction edges
            
        Returns:
        --------
        List of interactions [source, target, weight, source_layer, target_layer, interaction_type]
        """
        interactions = []
        
        # Check if required layers exist
        if not self.proteome or not self.metabolome or not self.reactions:
            logger.warning("Proteome, Metabolome, or Reactions layer missing")
            return interactions
            
        # Get metabolites in the network
        present_metabolites = [
            metabolite.kegg_compound_id for metabolite in self.metabolome.metabolites
        ]
        
        # Get reactions
        present_reactions = []
        ec_reaction = {}
        for reaction in self.reactions.reactions:
            if '   ' in str(reaction.enzyme):
                present_reactions.append(reaction.enzyme.split('   ')[0])
                ec_reaction[reaction.enzyme.split('   ')[0]] = reaction.id
            else:
                present_reactions.append(reaction.enzyme)
                ec_reaction[reaction.enzyme] = reaction.id
        
        # Create EC number to UniProt mapping
        ec_uniprot = {}
        for enzyme in self.proteome.proteins:
            for ec in enzyme.ec_number:
                ec_uniprot[ec] = enzyme.uniprot_id
        
        # Create interactions for activators
        for enzyme in self.proteome.proteins:
            if len(enzyme.activators) > 0:
                for ec in enzyme.ec_number:
                    if ec in present_reactions:
                        for activator in enzyme.activators:
                            if isinstance(activator, tuple):
                                activator_id = activator[1]  # Assuming tuple (name, id)
                            else:
                                activator_id = activator
                                
                            if activator_id in present_metabolites:
                                # Protein to reaction interaction (activation)
                                new_interaction = [
                                    enzyme.uniprot_id,
                                    ec_reaction[ec],
                                    weight,
                                    "Proteome",
                                    "Reactions",
                                    "Activator"
                                ]
                                interactions.append(new_interaction)
                                
                                # Reaction to metabolite interaction (activation)
                                new_interaction = [
                                    ec_reaction[ec],
                                    activator_id,
                                    weight,
                                    "Reactions",
                                    "Metabolome",
                                    "Activator"
                                ]
                                interactions.append(new_interaction)
        
        # Create interactions for inhibitors
        for enzyme in self.proteome.proteins:
            if len(enzyme.inhibitors) > 0:
                for ec in enzyme.ec_number:
                    if ec in present_reactions:
                        for inhibitor in enzyme.inhibitors:
                            if isinstance(inhibitor, tuple):
                                inhibitor_id = inhibitor[1]  # Assuming tuple (name, id)
                            else:
                                inhibitor_id = inhibitor
                                
                            if inhibitor_id in present_metabolites:
                                # Protein to reaction interaction (inhibition)
                                new_interaction = [
                                    enzyme.uniprot_id,
                                    ec_reaction[ec],
                                    weight,
                                    "Proteome",
                                    "Reactions",
                                    "Inhibitor"
                                ]
                                interactions.append(new_interaction)
                                
                                # Reaction to metabolite interaction (inhibition)
                                new_interaction = [
                                    ec_reaction[ec],
                                    inhibitor_id,
                                    weight,
                                    "Reactions",
                                    "Metabolome",
                                    "Inhibitor"
                                ]
                                interactions.append(new_interaction)
        
        logger.info(f"Created {len(interactions)} activator/inhibitor interactions")
        return interactions

    def generate_interaction_df(self):
        """
        Generate a dataframe of all interactions in the network.
        
        Returns:
        --------
        pd.DataFrame
            Dataframe with columns [SourceNode, TargetNode, Weight, SourceLayer, TargetLayer]
        """
        column_names = [
            "SourceNode",
            "TargetNode",
            "Weight",
            "SourceLayer",
            "TargetLayer",
        ]
        
        weight = 1.0
        interactions = []
        
        # Generate all types of interactions
        interactions.extend(self.protein_metabolite_interaction(weight))
        interactions.extend(self.protein_protein_interaction(weight))
        interactions.extend(self.gene_protein_interaction(weight))
        interactions.extend(self.transcription_factor_interaction(weight))
        interactions.extend(self.enzyme_reaction_interaction(weight))
        
        # Create dataframe and remove duplicates
        interaction_df = pd.DataFrame(interactions, columns=column_names)
        interaction_df = interaction_df.drop_duplicates().reset_index(drop=True)
        
        logger.info(f"Generated interaction dataframe with {len(interaction_df)} interactions")
        return interaction_df
    
    def generate_graph(self):
        """
        Generate a NetworkX graph from the network.
        
        Returns:
        --------
        nx.Graph
            NetworkX graph of the network
        """
        # Get interaction dataframe
        df = self.generate_interaction_df()
        
        # Create graph
        G = nx.Graph()
        
        # Add nodes with attributes
        node_layers = {}
        for _, row in df.iterrows():
            node_layers[row["SourceNode"]] = row["SourceLayer"]
            node_layers[row["TargetNode"]] = row["TargetLayer"]
        
        for node, layer in node_layers.items():
            G.add_node(node, layer=layer)
        
        # Add edges with weights
        for _, row in df.iterrows():
            G.add_edge(row["SourceNode"], row["TargetNode"], weight=row["Weight"])
        
        self.graph = G
        logger.info(f"Generated graph with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
        return G
    
    def generate_adjacency_matrix(self):
        """
        Generate an adjacency matrix from the network.
        
        Returns:
        --------
        pd.DataFrame
            Adjacency matrix of the network
        """
        # This potentially will have to be changed to separate elements that are in more than one layer
        # This could be done by adding a prefix to the element name, e.g. 'm_' for metabolite, 'p_' for protein, etc.

        interaction_df = self.generate_interaction_df()
        # Create a list of all elements in the network
        interacting_elements = list(set(interaction_df['SourceNode'].tolist() + interaction_df['TargetNode'].tolist()))
        adjacency_matrix = pd.DataFrame(index=interacting_elements, columns=interacting_elements)
        adjacency_matrix = adjacency_matrix.fillna(0)

        # Fill the adjacency matrix, currently undirected
        for index, row in interaction_df.iterrows():
            adjacency_matrix.loc[row['SourceNode'], row['TargetNode']] = row['Weight']
            adjacency_matrix.loc[row['TargetNode'], row['SourceNode']] = row['Weight']

        logger.info(f"Generated adjacency matrix of size {adjacency_matrix.shape}")
        return adjacency_matrix
    
    def save_network(self, output_dir: str):
        """
        Save the network to CSV files.
        
        Parameters:
        -----------
        output_dir : str
            Directory to save the network files
        """
        import os
        
        # Create directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Generate interaction dataframe
        df = self.generate_interaction_df()
        
        # Save interaction dataframe
        df.to_csv(os.path.join(output_dir, "interactions.csv"), index=False)
        
        # Save node information
        node_data = []
        for layer_name, layer in [
            ("Transcriptome", self.transcriptome), 
            ("Proteome", self.proteome), 
            ("Metabolome", self.metabolome), 
            ("Reactions", self.reactions)
        ]:
            if layer is None:
                continue
                
            if layer_name == "Transcriptome" and layer.genes:
                for gene in layer.genes:
                    node_data.append({
                        "ID": gene.ncbi_id or gene.ensembl_id,
                        "Name": gene.name,
                        "Type": "Gene",
                        "Layer": layer_name,
                        "FC": gene.fc,
                        "P_value": gene.adj_p_value
                    })
            elif layer_name == "Proteome" and layer.proteins:
                for protein in layer.proteins:
                    node_data.append({
                        "ID": protein.uniprot_id,
                        "Name": protein.name,
                        "Type": "Protein",
                        "Layer": layer_name,
                        "FC": protein.fc,
                        "P_value": protein.adj_p_value
                    })
            elif layer_name == "Metabolome" and layer.metabolites:
                for metabolite in layer.metabolites:
                    node_data.append({
                        "ID": metabolite.kegg_compound_id,
                        "Name": metabolite.kegg_name,
                        "Type": "Metabolite",
                        "Layer": layer_name,
                        "FC": metabolite.fc,
                        "P_value": metabolite.adj_p_value
                    })
            elif layer_name == "Reactions" and layer.reactions:
                for reaction in layer.reactions:
                    node_data.append({
                        "ID": reaction.id,
                        "Name": reaction.name,
                        "Type": "Reaction",
                        "Layer": layer_name,
                        "FC": None,
                        "P_value": None
                    })
        
        # Save node data
        node_df = pd.DataFrame(node_data)
        node_df.to_csv(os.path.join(output_dir, "nodes.csv"), index=False)
        
        # Save adjacency matrix
        adj_matrix = self.generate_adjacency_matrix()
        adj_matrix.to_csv(os.path.join(output_dir, "adjacency_matrix.csv"))
        
        logger.info(f"Saved network to {output_dir}")
    
    @classmethod
    def load_network(cls, input_dir: str):
        """
        Load a network from CSV files.
        
        Parameters:
        -----------
        input_dir : str
            Directory containing the network files
            
        Returns:
        --------
        Transnet
            Loaded network
        """
        import os
        
        # Create new Transnet object
        transnet = cls()
        
        # Load interaction dataframe
        df = pd.read_csv(os.path.join(input_dir, "interactions.csv"))
        
        # Extract layers from interaction dataframe
        layers = set(df["SourceLayer"].tolist() + df["TargetLayer"].tolist())
        
        # Initialize layers
        if "Transcriptome" in layers:
            transnet.transcriptome = Transcriptome()
            transnet.transcriptome.genes = []
        
        if "Proteome" in layers:
            transnet.proteome = Proteome()
            transnet.proteome.proteins = []
        
        if "Metabolome" in layers:
            transnet.metabolome = Metabolome()
            transnet.metabolome.metabolites = []
        
        if "Reactions" in layers:
            transnet.reactions = Reactions()
            transnet.reactions.reactions = []
        
        # Load node information
        node_df = pd.read_csv(os.path.join(input_dir, "nodes.csv"))
        
        # Create nodes based on their types
        for _, row in node_df.iterrows():
            if row["Type"] == "Gene":
                gene = Gene(
                    ncbi_id=row["ID"],
                    name=row["Name"],
                    fc=row["FC"],
                    adj_p_value=row["P_value"]
                )
                transnet.transcriptome.genes.append(gene)
            elif row["Type"] == "Protein":
                protein = Protein(
                    uniprot_id=row["ID"],
                    name=row["Name"],
                    fc=row["FC"],
                    adj_p_value=row["P_value"]
                )
                transnet.proteome.proteins.append(protein)
            elif row["Type"] == "Metabolite":
                metabolite = Metabolite(
                    kegg_compound_id=row["ID"],
                    kegg_name=row["Name"],
                    fc=row["FC"],
                    adj_p_value=row["P_value"]
                )
                transnet.metabolome.metabolites.append(metabolite)
            elif row["Type"] == "Reaction":
                reaction = Reaction(
                    id=row["ID"],
                    name=row["Name"]
                )
                transnet.reactions.reactions.append(reaction)
        
        logger.info(f"Loaded network from {input_dir}")
        return transnet
        