# %% [markdown]
# # Which axis regulates each reaction
#
# This is the type of question trans-omics is designed to answer. A metabolic reaction is
# driven from two directions at once:
#
# * the **enzyme axis**: how much enzyme there is, set by transcription and
#   translation;
# * the **metabolite axis**: how hard that enzyme works, set by its
#   substrates and its allosteric regulators.
#
# A list of changed molecules cannot separate them. The network can, per
# reaction, and can say when the two **disagree**.

# %%
import pandas as pd

from transnet import (
    load_example_network,
    load_example_omics,
    map_omics_to_network,
    reaction_regulation_table,
    regulation_axis_summary,
)

graph = load_example_network()
map_omics_to_network(
    graph, load_example_omics("insulin_resistant"),
    id_column="id", log2fc_column="log2FC", qvalue_column="padj",
)

table = reaction_regulation_table(graph)
regulated = table[(table["gene_axis"] != 0) | (table["metabolite_axis"] != 0)]
print(f"{len(regulated)} of {len(table)} reactions are regulated")
table.head(6)[["reaction", "name", "gene_axis", "metabolite_axis",
               "gene_axis_evidence", "controversial"]]

# %% [markdown]
# ## What the evidence actually was
#
# `gene_axis_evidence` records how far up the hierarchy the call reached. A
# reaction supported by `tf_gene_protein` is a stronger claim than one
# supported by `protein` alone, and a missing layer weakens a claim visibly
# instead of silently.

# %%
table["gene_axis_evidence"].value_counts(dropna=False)

# %% [markdown]
# ## Controversial reactions
#
# The two axes point opposite ways: more enzyme with less substrate, or less
# enzyme with more activator. Reading one layer alone gives the incomplete answers for
# these reactions.

# %%
controversial = table[table["controversial"]]
controversial[["reaction", "name", "gene_axis", "metabolite_axis", "allosteric_regulators"]]

# %% [markdown]
# `plot_controversial_reactions` draws each one as a tug of war: every bar is one
# measured molecule's push on the reaction, coloured by the axis it belongs to.
# The enzymes and the metabolites pulling against each other are named, which is
# what turns a flag into a hypothesis.

# %%
from transnet.visualization import plot_controversial_reactions

plot_controversial_reactions(graph, table)

# %% [markdown]
# ## Per-pathway balance
#
# For each pathway, how much of it is driven by which axis.

# %%
summary = regulation_axis_summary(table)
summary

# %%
import matplotlib.pyplot as plt

from transnet.visualization import plot_regulation_axes

plot_regulation_axes(summary)
plt.show()
