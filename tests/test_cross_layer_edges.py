"""
Tests for cross-layer edge functionality.
"""

import pytest
import pandas as pd
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from transnet.biology.transnet import Transnet
from transnet.biology.layers import Proteome, Metabolome, Transcriptome
from transnet.biology.elements import Protein, Metabolite, Gene


@pytest.fixture
def simple_network():
    """Create a simple test network."""
    # Create a simple proteome
    proteome = Proteome()
    proteome.ncbi_organism = "9606"
    proteome.kegg_organism = "hsa"
    
    # Add a few proteins with metabolites
    protein1 = Protein(
        uniprot_id="P12345",
        name="Test Protein 1",
        ec_number=["1.1.1.1"]
    )
    protein1.metabolites = ["C00002", "C00008"]  # ATP, ADP
    
    protein2 = Protein(
        uniprot_id="P67890",
        name="Test Protein 2",
        ec_number=["2.7.1.1"]
    )
    protein2.metabolites = ["C00008", "C00020"]  # ADP, AMP
    
    protein1.interaction_partners = ["P67890"]
    
    proteome.proteins = [protein1, protein2]
    
    # Create a simple metabolome
    metabolome = Metabolome()
    met1 = Metabolite(kegg_compound_id="C00002", kegg_name="ATP")
    met2 = Metabolite(kegg_compound_id="C00008", kegg_name="ADP")
    met3 = Metabolite(kegg_compound_id="C00020", kegg_name="AMP")
    metabolome.metabolites = [met1, met2, met3]
    
    # Create a simple transcriptome
    transcriptome = Transcriptome()
    gene1 = Gene(
        ensembl_id="ENSG00001",
        name="Gene1",
        kegg_id="hsa:12345"
    )
    gene2 = Gene(
        ensembl_id="ENSG00002",
        name="Gene2",
        kegg_id="hsa:67890"
    )
    transcriptome.genes = [gene1, gene2]
    
    # Create network
    transnet = Transnet(
        name="test_network",
        proteome=proteome,
        metabolome=metabolome,
        transcriptome=transcriptome
    )
    
    return transnet


def test_build_edge_index(simple_network):
    """Test building the cross-layer edge index."""
    edges = simple_network.build_cross_layer_edge_index()
    
    assert isinstance(edges, pd.DataFrame)
    assert not edges.empty
    
    # Check required columns
    required_cols = ['source', 'target', 'source_layer', 'target_layer', 
                     'edge_type', 'weight', 'evidence']
    for col in required_cols:
        assert col in edges.columns
    
    # Check that we have protein-metabolite edges
    enzymatic_edges = edges[edges['edge_type'] == 'enzymatic']
    assert len(enzymatic_edges) > 0
    
    # Check that edges have correct layers
    assert all(enzymatic_edges['source_layer'] == 'Proteome')
    assert all(enzymatic_edges['target_layer'] == 'Metabolome')


def test_edge_types_index(simple_network):
    """Test that edge type indices are created correctly."""
    simple_network.build_cross_layer_edge_index()
    
    assert isinstance(simple_network.edge_types, dict)
    assert len(simple_network.edge_types) > 0
    
    # Check that enzymatic edges are indexed
    if 'enzymatic' in simple_network.edge_types:
        enzymatic_df = simple_network.edge_types['enzymatic']
        assert isinstance(enzymatic_df, pd.DataFrame)
        assert all(enzymatic_df['edge_type'] == 'enzymatic')


def test_get_neighbors_basic(simple_network):
    """Test basic neighbor retrieval."""
    simple_network.build_cross_layer_edge_index()
    
    # Get neighbors of a protein
    neighbors = simple_network.get_neighbors("P12345")
    
    assert isinstance(neighbors, list)
    assert len(neighbors) > 0
    
    # Should include metabolites and other proteins
    assert "C00002" in neighbors or "C00008" in neighbors


def test_get_neighbors_filtered_by_layer(simple_network):
    """Test neighbor retrieval filtered by layer."""
    simple_network.build_cross_layer_edge_index()
    
    # Get only metabolite neighbors of a protein
    neighbors = simple_network.get_neighbors("P12345", layer="Metabolome")
    
    assert isinstance(neighbors, list)
    
    # All neighbors should be metabolites (C-prefixed)
    for n in neighbors:
        assert n.startswith("C")


def test_get_neighbors_filtered_by_edge_type(simple_network):
    """Test neighbor retrieval filtered by edge type."""
    simple_network.build_cross_layer_edge_index()
    
    # Get neighbors via enzymatic edges
    neighbors = simple_network.get_neighbors("P12345", edge_type="enzymatic")
    
    assert isinstance(neighbors, list)


def test_get_neighbors_direction(simple_network):
    """Test directional neighbor retrieval."""
    simple_network.build_cross_layer_edge_index()
    
    # Outgoing edges
    outgoing = simple_network.get_neighbors("P12345", direction="outgoing")
    
    # Incoming edges
    incoming = simple_network.get_neighbors("P12345", direction="incoming")
    
    # Both directions
    both = simple_network.get_neighbors("P12345", direction="both")
    
    assert isinstance(outgoing, list)
    assert isinstance(incoming, list)
    assert isinstance(both, list)
    
    # 'both' should contain all neighbors from outgoing and incoming
    assert set(both) == set(outgoing + incoming)


def test_find_paths(simple_network):
    """Test finding paths between nodes."""
    simple_network.build_cross_layer_edge_index()
    
    # Try to find paths from protein to metabolite
    paths = simple_network.find_paths("P12345", "C00002", max_length=2)
    
    assert isinstance(paths, list)
    
    # If paths exist, they should be lists of node IDs
    if len(paths) > 0:
        assert isinstance(paths[0], list)
        assert len(paths[0]) >= 2  # At least source and target


def test_find_paths_not_exist(simple_network):
    """Test finding paths that don't exist."""
    simple_network.build_cross_layer_edge_index()
    
    # Try to find paths between unconnected nodes
    paths = simple_network.find_paths("P12345", "NONEXISTENT", max_length=3)
    
    assert isinstance(paths, list)
    assert len(paths) == 0


def test_get_path_annotations(simple_network):
    """Test getting annotations for a path."""
    simple_network.build_cross_layer_edge_index()
    
    # Create a simple path
    path = ["P12345", "C00002"]
    
    annotations = simple_network.get_path_annotations(path)
    
    assert isinstance(annotations, list)
    assert len(annotations) == len(path) - 1  # One annotation per edge
    
    # Check annotation structure
    if len(annotations) > 0:
        ann = annotations[0]
        assert 'source' in ann
        assert 'target' in ann
        assert 'edge_type' in ann
        assert 'weight' in ann
        assert 'evidence' in ann


def test_query_cross_layer_relationships(simple_network):
    """Test querying cross-layer relationships."""
    simple_network.build_cross_layer_edge_index()
    
    # Create a set of changed nodes
    changed_nodes = {
        'Proteome': ['P12345'],
        'Metabolome': ['C00002', 'C00008']
    }
    
    results = simple_network.query_cross_layer_relationships(
        changed_nodes,
        path_max_length=3
    )
    
    assert isinstance(results, pd.DataFrame)
    
    # If results exist, check structure
    if not results.empty:
        required_cols = ['source', 'target', 'source_layer', 'target_layer',
                        'path_length', 'path', 'edge_types', 'evidence']
        for col in required_cols:
            assert col in results.columns


def test_edge_index_persistence(simple_network):
    """Test that edge index is built only once."""
    # First build
    simple_network.build_cross_layer_edge_index()
    edges_count_1 = len(simple_network.cross_layer_edges)
    
    # Second call should use cached version
    simple_network.build_cross_layer_edge_index()
    edges_count_2 = len(simple_network.cross_layer_edges)
    
    assert edges_count_1 == edges_count_2
    assert simple_network._edge_index_built is True


def test_empty_network():
    """Test edge index with empty network."""
    transnet = Transnet(name="empty")
    
    edges = transnet.build_cross_layer_edge_index()
    
    assert isinstance(edges, pd.DataFrame)
    assert edges.empty
    assert len(transnet.edge_types) == 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])