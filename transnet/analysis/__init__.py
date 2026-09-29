"""TransNet analysis.

The headline API is :mod:`transnet.analysis.transomics` -- the analyses that
are meaningful only because the network's edges are typed, directed and
signed. Its names are re-exported here.

Beside it, three modules that read the same network in ordinary ways:

:mod:`transnet.analysis.network_analysis`
    Graph statistics, centrality, communities and active modules. Useful as
    the contrast to the trans-omic readings: degree centrality finds the
    best-connected molecule *within* the graph, while
    :func:`~transnet.analysis.transomics.transomic_hubs` finds the one that
    joins the most layers, and they rarely agree.
:mod:`transnet.analysis.network_propagation`
    Undirected random-walk-with-restart and a permutation test for
    connectivity. For regulatory flow that respects direction and sign, use
    :func:`~transnet.analysis.transomics.hierarchical_propagation`.
:mod:`transnet.analysis.data_integration`
    Preprocessing: differential expression, time-course statistics,
    trajectory clustering and identifier mapping.

:mod:`transnet.analysis.factors` fits NMF factors and reads them through the
network; import it explicitly. Other decompositions come from scikit-learn.
"""

from transnet.analysis import (
    data_integration,
    network_analysis,
    network_propagation,
    summaries,
    transomics,
)

from transnet.analysis.transomics import *          # noqa: F401,F403
from transnet.analysis.transomics import __all__ as _transomics_all

from transnet.analysis.data_integration import (
    cluster_temporal_trajectories,
    compute_differential_expression,
    compute_timecourse_statistics,
    id_mapping,
)
from transnet.analysis.network_analysis import (
    compute_centrality_measures,
    compute_network_statistics,
    detect_communities,
    enrichment_analysis,
    find_active_modules,
    identify_hubs,
)
from transnet.analysis.network_propagation import (
    network_enrichment_permutation,
    random_walk_with_restart,
)
from transnet.analysis.summaries import (
    hub_rankings,
    path_verdicts,
    versus_chance,
)

__all__ = [
    "data_integration",
    "network_analysis",
    "network_propagation",
    "summaries",
    "transomics",
    *_transomics_all,
    "cluster_temporal_trajectories",
    "compute_differential_expression",
    "compute_timecourse_statistics",
    "id_mapping",
    "compute_centrality_measures",
    "compute_network_statistics",
    "detect_communities",
    "enrichment_analysis",
    "find_active_modules",
    "identify_hubs",
    "network_enrichment_permutation",
    "random_walk_with_restart",
    "hub_rankings",
    "path_verdicts",
    "versus_chance",
]
