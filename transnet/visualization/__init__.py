"""TransNet visualization.

The trans-omics figures -- :func:`plot_transomic_network`,
:func:`plot_regulation_axes` and :func:`plot_layer_connectivity` -- draw the
network as a layered, typed, signed object.  Beside them are the
ordinary graph figures -- the census, communities, a values heatmap -- and an
interactive plotly view that writes a standalone HTML file.
"""

from .transomics_vis import (
    plot_transomic_network,
    transomic_backbone,
    plot_regulation_axes,
    plot_layer_connectivity,
    LAYER_COLORS,
    EDGE_STYLES,
)
from .findings import (
    plot_axis_composition,
    plot_controversial_reactions,
    plot_metabolite_regulators,
    plot_regulatory_paths,
    plot_transomic_hubs,
    plot_tf_activity,
    plot_expression_concordance,
    plot_downstream_influence,
    plot_layer_changes,
    plot_similarity_heatmap,
    plot_temporal_structure,
    plot_condition_comparison,
    plot_convergence_null,
    plot_modification_axis,
    plot_regulatory_motifs,
    plot_structural_vulnerability,
)
from .factor_vis import (
    plot_factor_network,
    plot_factor_overview,
    plot_factor_scores,
)
from .network_vis import (
    plot_community_network,
    plot_network_metrics,
    plot_values_heatmap,
)
from .interactive import plot_transomic_network_interactive

__all__ = [
    # trans-omics figures
    'plot_transomic_network',
    'transomic_backbone',
    'plot_regulation_axes',
    'plot_layer_connectivity',
    'LAYER_COLORS',
    'EDGE_STYLES',
    # what each analysis finds
    'plot_axis_composition',
    'plot_controversial_reactions',
    'plot_metabolite_regulators',
    'plot_regulatory_paths',
    'plot_transomic_hubs',
    'plot_tf_activity',
    'plot_expression_concordance',
    'plot_downstream_influence',
    'plot_layer_changes',
    'plot_similarity_heatmap',
    'plot_temporal_structure',
    'plot_condition_comparison',
    'plot_convergence_null',
    'plot_structural_vulnerability',
    'plot_regulatory_motifs',
    'plot_modification_axis',
    'plot_factor_network',
    'plot_factor_overview',
    'plot_factor_scores',
    # general graph figures
    'plot_community_network',
    'plot_network_metrics',
    'plot_values_heatmap',
    'plot_transomic_network_interactive',
]
