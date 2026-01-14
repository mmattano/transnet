"""Network analysis and data integration tools"""

from transnet.analysis import (
    network_analysis, 
    network_loading, 
    data_integration,
    multi_omics_integration
)

# Import key classes for convenience
from transnet.analysis.multi_omics_integration import (
    MultiOmicsIntegrator,
    integrate_with_transnet,
    compute_factor_network_enrichment,
    perform_factor_enrichment
)

__all__ = [
    'network_analysis', 
    'network_loading', 
    'data_integration',
    'multi_omics_integration',
    'MultiOmicsIntegrator',
    'integrate_with_transnet',
    'compute_factor_network_enrichment',
    'perform_factor_enrichment'
]