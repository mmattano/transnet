"""API modules for external database access"""

from transnet.api import (
    kegg, uniprot, string, ensembl, chip_atlas, chem_info, brenda, reactome, hmdb
)

# KEGG: the signalling relations that build the Signaling layer. The rest of the
# KEGG client is reached through the module, since it is builder machinery.
from transnet.api.kegg import kegg_signaling_relations

# BRENDA
from transnet.api.brenda import (
    BrendaClient,
    brenda_get_km_values,
    brenda_get_kcat_values,
    brenda_get_ki_values,
    brenda_get_inhibitors,
    brenda_get_activators,
    brenda_get_substrates,
    brenda_get_products,
    brenda_get_cofactors,
    brenda_enrich_proteins,
)

# Reactome
from transnet.api.reactome import (
    reactome_get_pathways,
    reactome_get_pathway_reactions,
    reactome_get_pathway_entities,
    reactome_map_ids_to_pathways,
    reactome_over_representation,
    reactome_get_pathway_hierarchy,
)

# HMDB
from transnet.api.hmdb import (
    hmdb_get_metabolite,
    hmdb_search_metabolites,
    hmdb_get_pathways,
    hmdb_get_diseases,
    hmdb_map_ids,
    hmdb_enrich_metabolites,
)

__all__ = [
    # sub-modules
    'kegg', 'uniprot', 'string', 'ensembl', 'chip_atlas', 'chem_info',
    'kegg_signaling_relations',
    'brenda', 'reactome', 'hmdb',
    # BRENDA
    'BrendaClient',
    'brenda_get_km_values', 'brenda_get_kcat_values', 'brenda_get_ki_values',
    'brenda_get_inhibitors', 'brenda_get_activators', 'brenda_get_substrates',
    'brenda_get_products', 'brenda_get_cofactors', 'brenda_enrich_proteins',
    # Reactome
    'reactome_get_pathways', 'reactome_get_pathway_reactions',
    'reactome_get_pathway_entities', 'reactome_map_ids_to_pathways',
    'reactome_over_representation', 'reactome_get_pathway_hierarchy',
    # HMDB
    'hmdb_get_metabolite', 'hmdb_search_metabolites',
    'hmdb_get_pathways', 'hmdb_get_diseases',
    'hmdb_map_ids', 'hmdb_enrich_metabolites',
]
