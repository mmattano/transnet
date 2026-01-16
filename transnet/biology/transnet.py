"""Connecting transomics network
"""

from typing import List, Dict, Any, Optional, Union, Tuple
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
        
        self.cross_layer_edges = pd.DataFrame()
        self.edge_types = {}
        self._edge_index_built = False

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
    
    def build_cross_layer_edge_index(self) -> pd.DataFrame:
        """
        Build queryable edge structure across all layers.
        
        Creates a comprehensive DataFrame of all cross-layer interactions with metadata
        about source/target layers, edge types, and biological evidence.
        
        Returns:
        --------
        pd.DataFrame
            DataFrame with columns: source, target, source_layer, target_layer,
            edge_type, weight, evidence
        """
        logger.info("Building cross-layer edge index")
        edges = []
        
        # 1. Protein-Metabolite edges from enzymatic reactions
        if self.proteome and self.metabolome:
            logger.info("Adding protein-metabolite edges")
            for protein in self.proteome.proteins:
                if hasattr(protein, 'metabolites') and protein.metabolites:
                    for metabolite in protein.metabolites:
                        edges.append({
                            'source': str(protein.uniprot_id),
                            'target': str(metabolite),
                            'source_layer': 'Proteome',
                            'target_layer': 'Metabolome',
                            'edge_type': 'enzymatic',
                            'weight': 1.0,
                            'evidence': f"EC:{','.join(protein.ec_number) if protein.ec_number else 'unknown'}"
                        })
        
        # 2. Gene-Protein edges (translation)
        if self.transcriptome and self.proteome:
            logger.info("Adding gene-protein edges")
            gene_protein_interactions = self.gene_protein_interaction()
            for interaction in gene_protein_interactions:
                edges.append({
                    'source': str(interaction[0]),  # gene
                    'target': str(interaction[1]),  # protein
                    'source_layer': 'Transcriptome',
                    'target_layer': 'Proteome',
                    'edge_type': 'translation',
                    'weight': float(interaction[2]),
                    'evidence': 'gene_to_protein'
                })
        
        # 3. TF-Gene edges (transcriptional regulation)
        if self.transcriptome and self.proteome:
            logger.info("Adding transcription factor-gene edges")
            tf_interactions = self.transcription_factor_interaction()
            for interaction in tf_interactions:
                edges.append({
                    'source': str(interaction[0]),  # TF protein
                    'target': str(interaction[1]),  # target gene
                    'source_layer': 'Proteome',
                    'target_layer': 'Transcriptome',
                    'edge_type': 'transcriptional_regulation',
                    'weight': float(interaction[2]),
                    'evidence': 'ChIP-Atlas'
                })
        
        # 4. Protein-Protein edges (PPI from STRING)
        if self.proteome:
            logger.info("Adding protein-protein interaction edges")
            ppi_interactions = self.protein_protein_interaction()
            for interaction in ppi_interactions:
                edges.append({
                    'source': str(interaction[0]),
                    'target': str(interaction[1]),
                    'source_layer': 'Proteome',
                    'target_layer': 'Proteome',
                    'edge_type': 'protein_interaction',
                    'weight': float(interaction[2]),
                    'evidence': 'STRING'
                })
        
        # 5. Enzyme-Reaction-Metabolite edges
        if self.proteome and self.metabolome and self.reactions:
            logger.info("Adding enzyme-reaction-metabolite edges")
            enzyme_rxn_interactions = self.enzyme_reaction_interaction()
            for interaction in enzyme_rxn_interactions:
                # These create two edges: enzyme->reaction and reaction->metabolite
                if interaction[3] == 'Proteome' and interaction[4] == 'Reactions':
                    edges.append({
                        'source': str(interaction[0]),
                        'target': str(interaction[1]),
                        'source_layer': 'Proteome',
                        'target_layer': 'Reactions',
                        'edge_type': 'catalysis',
                        'weight': float(interaction[2]),
                        'evidence': 'KEGG_reactions'
                    })
                elif interaction[3] == 'Reactions' and interaction[4] == 'Metabolome':
                    edges.append({
                        'source': str(interaction[0]),
                        'target': str(interaction[1]),
                        'source_layer': 'Reactions',
                        'target_layer': 'Metabolome',
                        'edge_type': 'reaction_participant',
                        'weight': float(interaction[2]),
                        'evidence': 'KEGG_reactions'
                    })
        
        # Create DataFrame
        self.cross_layer_edges = pd.DataFrame(edges)
        
        # Remove duplicates (keep first occurrence)
        if not self.cross_layer_edges.empty:
            self.cross_layer_edges = self.cross_layer_edges.drop_duplicates(
                subset=['source', 'target', 'edge_type'],
                keep='first'
            ).reset_index(drop=True)
        
        # Create edge type indices for fast lookup
        if not self.cross_layer_edges.empty:
            self.edge_types = {
                edge_type: self.cross_layer_edges[
                    self.cross_layer_edges['edge_type'] == edge_type
                ].copy()
                for edge_type in self.cross_layer_edges['edge_type'].unique()
            }
        else:
            self.edge_types = {}
        
        self._edge_index_built = True
        
        logger.info(
            f"Built cross-layer edge index with {len(self.cross_layer_edges)} edges "
            f"across {len(self.edge_types)} edge types"
        )
        
        return self.cross_layer_edges
    
    def get_neighbors(
        self, 
        node_id: str, 
        layer: Optional[str] = None, 
        edge_type: Optional[str] = None,
        direction: str = 'both'
    ) -> List[str]:
        """
        Get all neighbors of a node, optionally filtered by layer/type.
        
        Parameters:
        -----------
        node_id : str
            Node identifier
        layer : str, optional
            Filter neighbors by layer (e.g., 'Proteome', 'Metabolome')
        edge_type : str, optional
            Filter by edge type (e.g., 'enzymatic', 'translation')
        direction : str
            Direction of edges: 'outgoing', 'incoming', or 'both' (default)
        
        Returns:
        --------
        List[str]
            List of neighbor node IDs
        """
        if not self._edge_index_built:
            logger.warning("Edge index not built. Building now...")
            self.build_cross_layer_edge_index()
        
        if self.cross_layer_edges.empty:
            logger.warning("No edges in cross-layer edge index")
            return []
        
        neighbors = []
        
        # Get outgoing edges
        if direction in ['outgoing', 'both']:
            edges_from = self.cross_layer_edges[
                self.cross_layer_edges['source'] == str(node_id)
            ]
            
            if layer:
                edges_from = edges_from[edges_from['target_layer'] == layer]
            
            if edge_type:
                edges_from = edges_from[edges_from['edge_type'] == edge_type]
            
            neighbors.extend(edges_from['target'].tolist())
        
        # Get incoming edges
        if direction in ['incoming', 'both']:
            edges_to = self.cross_layer_edges[
                self.cross_layer_edges['target'] == str(node_id)
            ]
            
            if layer:
                edges_to = edges_to[edges_to['source_layer'] == layer]
            
            if edge_type:
                edges_to = edges_to[edges_to['edge_type'] == edge_type]
            
            neighbors.extend(edges_to['source'].tolist())
        
        return list(set(neighbors))
    
    def find_paths(
        self, 
        source: str, 
        target: str, 
        max_length: int = 3,
        allowed_layers: Optional[List[str]] = None,
        allowed_edge_types: Optional[List[str]] = None
    ) -> List[List[str]]:
        """
        Find all paths between two nodes up to max_length.
        
        This is the core method for cross-layer discovery. Given a changed
        gene and a changed metabolite, find mechanistic paths connecting them.
        
        Parameters:
        -----------
        source : str
            Source node ID
        target : str
            Target node ID
        max_length : int
            Maximum path length (number of edges)
        allowed_layers : List[str], optional
            Only use nodes from these layers
        allowed_edge_types : List[str], optional
            Only use edges of these types
        
        Returns:
        --------
        List[List[str]]
            List of paths, where each path is a list of node IDs
        
        Examples:
        ---------
        # Find paths from a gene to a metabolite
        >>> paths = transnet.find_paths('ENSG00000123456', 'C00002', max_length=4)
        >>> for path in paths:
        ...     print(' -> '.join(path))
        """
        if not self._edge_index_built:
            logger.warning("Edge index not built. Building now...")
            self.build_cross_layer_edge_index()
        
        # Build filtered graph if needed
        edges_to_use = self.cross_layer_edges.copy()
        
        if allowed_edge_types:
            edges_to_use = edges_to_use[
                edges_to_use['edge_type'].isin(allowed_edge_types)
            ]
        
        if allowed_layers:
            edges_to_use = edges_to_use[
                edges_to_use['source_layer'].isin(allowed_layers) &
                edges_to_use['target_layer'].isin(allowed_layers)
            ]
        
        # Build NetworkX graph
        G = nx.DiGraph()
        
        for _, row in edges_to_use.iterrows():
            G.add_edge(
                row['source'], 
                row['target'],
                weight=row['weight'],
                edge_type=row['edge_type']
            )
        
        # Find paths
        try:
            paths = list(nx.all_simple_paths(
                G, 
                str(source), 
                str(target), 
                cutoff=max_length
            ))
            logger.info(f"Found {len(paths)} paths from {source} to {target}")
            return paths
        except (nx.NodeNotFound, nx.NetworkXNoPath) as e:
            logger.warning(f"Could not find paths: {e}")
            return []
    
    def get_path_annotations(
        self, 
        path: List[str]
    ) -> List[Dict[str, Any]]:
        """
        Get annotations for edges in a path.
        
        Parameters:
        -----------
        path : List[str]
            List of node IDs representing a path
        
        Returns:
        --------
        List[Dict[str, Any]]
            List of edge annotations, one for each edge in the path
        """
        if not self._edge_index_built:
            self.build_cross_layer_edge_index()
        
        annotations = []
        
        for i in range(len(path) - 1):
            source = str(path[i])
            target = str(path[i + 1])
            
            # Find edge(s) between these nodes
            edge_matches = self.cross_layer_edges[
                (self.cross_layer_edges['source'] == source) &
                (self.cross_layer_edges['target'] == target)
            ]
            
            if not edge_matches.empty:
                # Take first match if multiple
                edge_data = edge_matches.iloc[0].to_dict()
                annotations.append(edge_data)
            else:
                # Edge not found (shouldn't happen if path is valid)
                annotations.append({
                    'source': source,
                    'target': target,
                    'edge_type': 'unknown',
                    'weight': 0.0,
                    'evidence': 'not_found'
                })
        
        return annotations
    
    def query_cross_layer_relationships(
        self,
        changed_nodes: Dict[str, List[str]],
        p_value_threshold: float = 0.05,
        path_max_length: int = 3
    ) -> pd.DataFrame:
        """
        Query cross-layer relationships between changed nodes.
        
        This is the main entry point for cross-layer discovery analysis.
        Given sets of changed nodes from different omics layers, find
        mechanistic paths connecting them.
        
        Parameters:
        -----------
        changed_nodes : Dict[str, List[str]]
            Dictionary mapping layer names to lists of changed node IDs
            Example: {'Transcriptome': ['gene1', 'gene2'], 
                     'Metabolome': ['met1', 'met2']}
        p_value_threshold : float
            P-value threshold for considering a node as changed
        path_max_length : int
            Maximum path length to search
        
        Returns:
        --------
        pd.DataFrame
            DataFrame with columns: source, target, path_length, path, 
            source_layer, target_layer, edge_types
        """
        if not self._edge_index_built:
            self.build_cross_layer_edge_index()
        
        logger.info(f"Querying cross-layer relationships among {sum(len(v) for v in changed_nodes.values())} changed nodes")
        
        results = []
        
        # For each pair of layers
        layer_pairs = [
            (l1, l2) for l1 in changed_nodes.keys() 
            for l2 in changed_nodes.keys() 
            if l1 != l2
        ]
        
        for source_layer, target_layer in layer_pairs:
            source_nodes = changed_nodes[source_layer]
            target_nodes = changed_nodes[target_layer]
            
            logger.info(f"Finding paths from {source_layer} to {target_layer}")
            
            # Find paths between all pairs
            for source in source_nodes:
                for target in target_nodes:
                    paths = self.find_paths(
                        source, 
                        target, 
                        max_length=path_max_length
                    )
                    
                    for path in paths:
                        # Get edge types in path
                        annotations = self.get_path_annotations(path)
                        edge_types = [a['edge_type'] for a in annotations]
                        
                        results.append({
                            'source': source,
                            'target': target,
                            'source_layer': source_layer,
                            'target_layer': target_layer,
                            'path_length': len(path) - 1,
                            'path': ' -> '.join(path),
                            'edge_types': ', '.join(edge_types),
                            'evidence': ', '.join([a['evidence'] for a in annotations])
                        })
        
        results_df = pd.DataFrame(results)
        
        if not results_df.empty:
            # Sort by path length
            results_df = results_df.sort_values('path_length').reset_index(drop=True)
            logger.info(f"Found {len(results_df)} cross-layer paths")
        else:
            logger.warning("No cross-layer paths found")
        
        return results_df
        