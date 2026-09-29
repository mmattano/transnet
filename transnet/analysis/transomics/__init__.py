"""Trans-omics network analysis -- TransNet's headline API.

A trans-omic network is a typed, directed, signed regulatory hierarchy:

    signal -> TF -> gene -> enzyme protein -> REACTION <- metabolite

with the metabolic reaction as the convergence point.  Every analysis here is
meaningful only because of that structure -- remove the edges and none of it can
be computed.  For generic multi-omics factor analysis, which does not use the
network, see :mod:`transnet.analysis.factor`.

The documented analysis catalogue
---------------------------------

=====  ============================================  ===========================================
Code   Analysis                                      Function
=====  ============================================  ===========================================
network reconstruction     Trans-omic network reconstruction             :meth:`~transnet.Transnet.generate_graph`,
                                                     :func:`responsive_subnetwork`
reaction regulation axes     Reaction regulation-axis attribution          :func:`reaction_regulation_table`
per-pathway balance     Per-pathway regulation balance                :func:`regulation_axis_summary`
signed regulatory paths     Signed regulatory-path tracing                :func:`trace_regulatory_paths`
cross-layer connectivity     Cross-layer connectivity                      :func:`cross_layer_connectivity`
trans-omic hubs     Trans-omic hub identification                 :func:`transomic_hubs`
response timing     Temporal / dose structure on the network      :func:`assign_temporal_parameters`,
                                                     :func:`temporal_network_structure`
=====  ============================================  ===========================================

Layers are optional
-------------------

No analysis here requires a particular layer.  Each one discovers what the
network actually contains and degrades by reporting weaker evidence rather than
by failing: :func:`reaction_regulation_table` records which chain of layers
supported its gene-axis call, and :func:`trace_regulatory_paths` infers its own
starting layer.  A transcriptome + metabolome network gets the full trans-omics
story for the layers it has.
"""

from transnet.analysis.transomics.mapping import (
    MappingReport,
    build_alias_id_map,
    build_name_id_map,
    is_responsive,
    map_modification_sites,
    map_omics_to_network,
    regulated_nodes,
)
from transnet.analysis.transomics.subnetwork import (
    compare_transomic_networks,
    responsive_subnetwork,
)
from transnet.analysis.transomics.regulation import (
    metabolite_regulatory_roles,
    reaction_regulation_table,
    regulation_axis_summary,
    regulatory_role_enrichment,
)
from transnet.analysis.transomics.paths import (
    path_consistency_summary,
    trace_regulatory_paths,
)
from transnet.analysis.transomics.structure import (
    convergence_significance,
    regulatory_motifs,
    structural_vulnerability,
)
from transnet.analysis.transomics.topology import (
    cross_layer_connectivity,
    layer_coverage,
    transomic_hubs,
)
from transnet.analysis.transomics.temporal import (
    assign_temporal_parameters,
    split_by_response_class,
    temporal_network_structure,
)
from transnet.analysis.transomics.propagation import (
    downstream_influence,
    hierarchical_propagation,
)

# Re-exported so callers can introspect a network without reaching into the
# biology package.
from transnet.analysis.transomics.gene_axis import (
    expression_concordance,
    transcription_factor_activity,
)

from transnet.biology.schema import (
    DISPLAY_ORDER,
    EDGE_TYPES,
    LAYER_HIERARCHY,
    available_edge_types,
    available_layers,
    top_layer_present,
)

__all__ = [
    # network reconstruction -- reconstruction and responsive subnetwork
    "responsive_subnetwork",
    "compare_transomic_networks",
    # data mapping
    "map_omics_to_network",
    "map_modification_sites",
    "MappingReport",
    "regulated_nodes",
    "is_responsive",
    "build_name_id_map",
    "build_alias_id_map",
    # reaction regulation axes and per-pathway balance
    "reaction_regulation_table",
    "regulation_axis_summary",
    "metabolite_regulatory_roles",
    "regulatory_role_enrichment",
    # signed regulatory paths -- signed paths
    # transcription-factor activity, transcript-protein concordance -- the gene-expression axis
    "transcription_factor_activity",
    "expression_concordance",
    "trace_regulatory_paths",
    "path_consistency_summary",
    # cross-layer connectivity, trans-omic hubs -- layer-aware topology
    "cross_layer_connectivity",
    "layer_coverage",
    "transomic_hubs",
    "regulatory_motifs",
    "structural_vulnerability",
    "convergence_significance",
    # response timing -- temporal structure
    "assign_temporal_parameters",
    "temporal_network_structure",
    "split_by_response_class",
    # regulatory-flow propagation
    "hierarchical_propagation",
    "downstream_influence",
    # schema introspection
    "EDGE_TYPES",
    "LAYER_HIERARCHY",
    "DISPLAY_ORDER",
    "available_layers",
    "available_edge_types",
    "top_layer_present",
]
