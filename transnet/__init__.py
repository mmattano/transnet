"""TransNet: trans-omics network reconstruction and analysis.

TransNet builds a trans-omic network, in which molecules from different omic
layers are joined by typed, directed and signed regulatory edges::

    signal -> TF -> gene -> enzyme protein -> REACTION <- metabolite

and analyses it. Each edge records its relationship (``edge_type``), whether
it increases or decreases its target (``sign``), and its source.

::

    from transnet import load_example_network, load_example_omics
    from transnet import map_omics_to_network, reaction_regulation_table

    graph = load_example_network()
    map_omics_to_network(graph, load_example_omics(), id_column="id",
                         log2fc_column="log2FC", qvalue_column="padj")
    table = reaction_regulation_table(graph)
    table[table["controversial"]]

The analyses are in :mod:`transnet.analysis.transomics`, factor models read
through the network in :mod:`transnet.analysis.factors`, figures in
:mod:`transnet.visualization` and export formats in :mod:`transnet.io`.
"""

try:                                    # one version, from the installed metadata
    from importlib.metadata import version as _version

    __version__ = _version("transnet")
except Exception:                       # pragma: no cover - source checkout
    __version__ = "0.2.0"

# --- the network -----------------------------------------------------------
from transnet.biology.transnet import Transnet, to_simple_graph
from transnet.biology.layers import (
    Pathways, Reactions, Transcriptome, Proteome, Metabolome, Signaling,
)
from transnet.biology.elements import (
    Reaction, Pathway, Metabolite, Gene, Protein, SignalingProtein,
)
from transnet.biology.schema import (
    EDGE_TYPES, LAYERS, LAYER_HIERARCHY, DISPLAY_ORDER,
    available_layers, available_edge_types, top_layer_present,
)

# --- trans-omics analysis (the headline API) -------------------------------
from transnet.analysis.transomics import (
    # data mapping
    map_omics_to_network, map_modification_sites, MappingReport,
    regulated_nodes, is_responsive,
    build_name_id_map, build_alias_id_map,
    # network reconstruction reconstruction / responsive subnetwork
    responsive_subnetwork, compare_transomic_networks,
    # reaction regulation axes and per-pathway balance
    reaction_regulation_table, regulation_axis_summary,
    transcription_factor_activity, expression_concordance,
    metabolite_regulatory_roles, regulatory_role_enrichment,
    # signed regulatory paths signed paths
    trace_regulatory_paths, path_consistency_summary,
    regulatory_motifs, structural_vulnerability, convergence_significance,
    # cross-layer connectivity, trans-omic hubs layer-aware topology
    cross_layer_connectivity, layer_coverage, transomic_hubs,
    # response timing temporal structure
    assign_temporal_parameters, temporal_network_structure,
    split_by_response_class,
    # regulatory-flow propagation
    hierarchical_propagation, downstream_influence,
)

# --- example data ----------------------------------------------------------
from transnet.data import (
    load_example_network,
    load_example_omics,
    load_example_pathways,
    load_example_phosphoproteomics,
    load_example_timecourse,
)

# --- database clients ------------------------------------------------------
from transnet.api.brenda import BrendaClient, brenda_enrich_proteins
from transnet.api.reactome import (
    reactome_get_pathways,
    reactome_over_representation,
)
from transnet.api.hmdb import (
    hmdb_get_metabolite,
    hmdb_search_metabolites,
    hmdb_enrich_metabolites,
)

__all__ = [
    # network
    'Transnet', 'to_simple_graph',
    'Pathways', 'Reactions', 'Transcriptome', 'Proteome', 'Metabolome',
    'Signaling',
    'Reaction', 'Pathway', 'Metabolite', 'Gene', 'Protein', 'SignalingProtein',
    'EDGE_TYPES', 'LAYERS', 'LAYER_HIERARCHY', 'DISPLAY_ORDER',
    'available_layers', 'available_edge_types', 'top_layer_present',
    # trans-omics analysis
    'map_omics_to_network', 'map_modification_sites', 'MappingReport', 'regulated_nodes', 'is_responsive', 'build_name_id_map', 'build_alias_id_map',
    'responsive_subnetwork', 'compare_transomic_networks',
    'reaction_regulation_table', 'regulation_axis_summary',
    'transcription_factor_activity', 'expression_concordance',
    'metabolite_regulatory_roles', 'regulatory_role_enrichment',
    'trace_regulatory_paths', 'path_consistency_summary',
    'regulatory_motifs', 'structural_vulnerability', 'convergence_significance',
    'cross_layer_connectivity', 'layer_coverage', 'transomic_hubs',
    'assign_temporal_parameters', 'temporal_network_structure',
    'split_by_response_class',
    'hierarchical_propagation', 'downstream_influence',
    # example data
    'load_example_network', 'load_example_omics',
    'load_example_pathways', 'load_example_phosphoproteomics', 'load_example_timecourse',
    # database clients
    'BrendaClient', 'brenda_enrich_proteins',
    'reactome_get_pathways', 'reactome_over_representation',
    'hmdb_get_metabolite', 'hmdb_search_metabolites', 'hmdb_enrich_metabolites',
]
