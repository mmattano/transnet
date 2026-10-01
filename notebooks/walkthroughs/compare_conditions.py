# %% [markdown]
# # Two conditions, compared as networks
#
# Comparing two networks usually means listing the nodes and edges that
# differ. Because every edge here has a type, the comparison can also say
# **which kinds of regulation** differ. A condition that loses its
# transcriptional regulation but keeps its allosteric regulation is a
# different biological situation from the reverse, and only a comparison that
# knows the edge types can tell them apart.
#
# The example compares the insulin-sensitive and insulin-resistant data on the
# same network. `compare_transomic_networks` takes the responsive subnetwork of
# each and reports shared and unique edges per edge type, and nodes per layer.

# %%
import pandas as pd

from transnet import (
    compare_transomic_networks,
    load_example_network,
    load_example_omics,
    map_omics_to_network,
    metabolite_regulatory_roles,
    reaction_regulation_table,
    regulatory_role_enrichment,
    responsive_subnetwork,
)

ARROW = {1: "activated", -1: "inhibited", 0: "-"}


def mapped(condition):
    graph = load_example_network()
    map_omics_to_network(graph, load_example_omics(condition), id_column="id",
                         log2fc_column="log2FC", qvalue_column="padj")
    return graph


sensitive, resistant = mapped("insulin_sensitive"), mapped("insulin_resistant")

comparison = compare_transomic_networks(
    responsive_subnetwork(sensitive), responsive_subnetwork(resistant),
    "sensitive", "resistant",
)
comparison["edges_by_type"]

# %% [markdown]
# The Jaccard index is the number of shared edges divided by the number of
# edges in either network: 1 means identical, 0 means nothing in common.

# %%
print(f"edge Jaccard: {comparison['summary']['edge_jaccard']:.2f}")
comparison["nodes_by_layer"]

# %% [markdown]
# ## Molecules that reversed direction
#
# These molecules changed significantly in both conditions, but in opposite
# directions. Comparing two lists of changed genes would count them as
# shared, and so as agreement.

# %%
shifts = comparison["regulation_shifts"]
shifts.assign(sensitive=shifts["sensitive"].map(ARROW),
              resistant=shifts["resistant"].map(ARROW)) if not shifts.empty else "none"

# %% [markdown]
# ## Which axis drives each reaction, in each condition
#
# The regulation-axis table from the *which axis regulates each reaction*
# walkthrough, computed for both conditions. Only reactions whose calls differ
# are shown.

# %%
axes = pd.concat(
    [reaction_regulation_table(g).set_index("reaction")[["name", "gene_axis", "metabolite_axis"]]
       .rename(columns={"gene_axis": f"{label} gene", "metabolite_axis": f"{label} metabolite"})
     for label, g in [("sensitive", sensitive), ("resistant", resistant)]],
    axis=1,
)
axes = axes.loc[:, ~axes.columns.duplicated()]
changed = axes[(axes["sensitive gene"] != axes["resistant gene"])
               | (axes["sensitive metabolite"] != axes["resistant metabolite"])]
changed

# %% [markdown]
# ## Do the changed metabolites regulate anything?
#
# A metabolite that changed may simply be carried along by the change in flux,
# or it may itself act on an enzyme as an allosteric regulator. The network
# holds the known regulators from BRENDA as signed allosteric edges.
# `metabolite_regulatory_roles` lists, for each changed metabolite, the
# reactions it activates or inhibits. `regulatory_role_enrichment` tests
# whether changed metabolites are regulators more often than measured
# metabolites in general (Fisher's exact test).

# %%
for label, graph in [("sensitive", sensitive), ("resistant", resistant)]:
    result = regulatory_role_enrichment(metabolite_regulatory_roles(graph))
    counts = result["counts"]
    row = result["enrichment"].set_index("role").loc["any"]
    print(f"{label:<11}{counts['n_differential_regulators']} of {counts['n_differential']} "
          f"changed metabolites regulate an enzyme "
          f"({counts['fraction_differential_regulators']:.0%} against a background of "
          f"{counts['fraction_background_regulators']:.0%}; q = {row['q_value']:.2g})")

# %%
roles = regulatory_role_enrichment(metabolite_regulatory_roles(sensitive))
roles["regulators"][["name", "log2fc", "reactions_activated", "reactions_inhibited"]]

# %% [markdown]
# `plot_metabolite_regulators` puts the share against its background beside the
# regulators themselves, so an unremarkable enrichment and a useful list of
# named metabolites can both be read off one figure.

# %%
from transnet.visualization import plot_metabolite_regulators

plot_metabolite_regulators(roles)

# %% [markdown]
# `plot_condition_comparison` summarises the comparison. The left panel shows,
# per edge type, the edges kept in both conditions and those gained or lost; the
# right panel lists the molecules that reversed direction.

# %%
import matplotlib.pyplot as plt

from transnet.visualization import plot_condition_comparison

plot_condition_comparison(comparison, "sensitive", "resistant")
plt.show()
