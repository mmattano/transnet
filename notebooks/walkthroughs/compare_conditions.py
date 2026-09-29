# %% [markdown]
# # Two conditions, compared as networks
#
# A differential-network comparison reports which nodes and edges differ. A
# *trans-omic* comparison asks the sharper question: **which kinds of
# regulation** changed. A condition that loses its transcriptional arm but
# keeps its allosteric one is a different story from the reverse, and only an
# edge-type-aware comparison can tell them apart.

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

# %%
print(f"edge Jaccard: {comparison['summary']['edge_jaccard']:.2f}")
comparison["nodes_by_layer"]

# %% [markdown]
# ## Molecules that reversed direction
#
# Not "changed in both": changed *the other way*, which a shared gene list
# reports as agreement.

# %%
shifts = comparison["regulation_shifts"]
shifts.assign(sensitive=shifts["sensitive"].map(ARROW),
              resistant=shifts["resistant"].map(ARROW)) if not shifts.empty else "none"

# %% [markdown]
# ## Which axis drives each reaction, in each condition

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
# A metabolite whose concentration moved is either a regulator acting back on an
# enzyme, which is a mechanistic hypothesis, or a passenger carried along by
# flux. BRENDA annotations separate the two. They are already in the network as
# signed allosteric edges.

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

# %%
import matplotlib.pyplot as plt

from transnet.visualization import plot_condition_comparison

plot_condition_comparison(comparison, "sensitive", "resistant")
plt.show()
