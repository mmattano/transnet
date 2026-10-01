# %% [markdown]
# # Which axis regulates each reaction
#
# This is the central question trans-omics is designed to answer. The rate of
# a metabolic reaction can change for two different reasons:
#
# * the **enzyme axis**: the amount of enzyme changed, because its gene or
#   protein went up or down;
# * the **metabolite axis**: the same amount of enzyme works faster or slower,
#   because its substrates, products or allosteric regulators changed.
#
# A list of changed molecules cannot tell these apart. The network can,
# because it knows which enzyme catalyses each reaction and which metabolites
# act on it. For every reaction, `reaction_regulation_table` reports the
# direction each axis pushes, and flags the reactions where the two axes push in
# **opposite** directions. These are called *controversial* reactions.
#
# This walkthrough goes through the table, the per-pathway summary, the figures,
# and two additions: a third axis from phosphorylation sites, and whether the
# protein changes follow their transcripts.

# %%
import matplotlib.pyplot as plt
import pandas as pd

from transnet import (
    expression_concordance,
    load_example_network,
    load_example_omics,
    load_example_pathways,
    map_modification_sites,
    map_omics_to_network,
    reaction_regulation_table,
    regulation_axis_summary,
)
from transnet.visualization import (
    plot_controversial_reactions,
    plot_expression_concordance,
    plot_regulation_axes,
)

graph = load_example_network()
map_omics_to_network(
    graph, load_example_omics("insulin_sensitive"),
    id_column="id", log2fc_column="log2FC", qvalue_column="padj",
)

table = reaction_regulation_table(graph)
regulated = table[(table["gene_axis"] != 0) | (table["metabolite_axis"] != 0)]
print(f"{len(regulated)} of {len(table)} reactions are regulated on at least one axis")
table.head(6)[["reaction", "name", "gene_axis", "metabolite_axis",
               "gene_axis_evidence", "controversial"]]

# %% [markdown]
# `gene_axis` and `metabolite_axis` are +1 (the axis speeds the reaction up),
# -1 (slows it down) or 0 (no change). For the metabolite axis, more substrate
# or activator counts as +1, and more product or inhibitor as -1.
#
# ## What the enzyme-axis call is based on
#
# `gene_axis_evidence` records which measurements the enzyme-axis call used:
#
# | Value | Meaning |
# |---|---|
# | `protein` | the enzyme protein changed |
# | `gene_protein` | the protein changed, and its transcript changed the same way |
# | `gene` | only the transcript was measured or changed |
# | `tf_gene_protein` | as `gene_protein`, and a transcription factor of the gene also changed |
# | empty | no enzyme change |
#
# A study without a proteome still gets an enzyme axis, from transcripts, and
# this column makes the weaker evidence visible.

# %%
table["gene_axis_evidence"].value_counts(dropna=False).to_frame("reactions")

# %% [markdown]
# ## Controversial reactions
#
# In a controversial reaction the two axes point opposite ways, for example
# more enzyme but also more of an inhibitor. Looking at one layer alone would
# give the wrong answer for these reactions: the transcriptome says the
# reaction speeds up, the metabolome says it slows down, and only together do
# they show a conflict.

# %%
controversial = table[table["controversial"]]
controversial[["reaction", "name", "gene_axis", "metabolite_axis", "allosteric_regulators"]]

# %% [markdown]
# `plot_controversial_reactions` draws each controversial reaction as a tug of
# war. Every bar is one measured molecule's push on the reaction: bars to the
# right speed it up, bars to the left slow it down, and the colour shows the
# axis. Naming the enzymes and metabolites that pull against each other turns a
# flag into a testable hypothesis. In hexokinase, for example, more Hk2 enzyme
# meets more glucose-6-phosphate, the enzyme's own product and a known
# inhibitor.

# %%
plot_controversial_reactions(graph, table)
plt.show()

# %% [markdown]
# ## Per-pathway balance
#
# `regulation_axis_summary` counts, for each pathway, how many regulated
# reactions each axis activates or inhibits and how many are controversial. It
# needs to know which pathway each reaction belongs to. The example data comes
# with that map. For an organism-wide network, use
# `transnet.api.kegg_reaction_pathways`, which reads it from KEGG. A reaction in
# several pathways is counted in each.

# %%
pathways = load_example_pathways()
summary = regulation_axis_summary(table, pathway_map=pathways)
summary[["pathway", "n_reactions", "gene_activated", "gene_inhibited",
         "metabolite_activated", "metabolite_inhibited", "n_controversial"]]

# %% [markdown]
# `plot_regulation_axes` draws the same table. The left panel is the enzyme
# axis and the middle panel the metabolite axis, with activation to the right
# (red) and inhibition to the left (blue). The right panel is the fraction of
# each pathway's regulated reactions that are controversial.

# %%
plot_regulation_axes(summary)
plt.show()

# %% [markdown]
# ## A third axis: phosphorylation of the enzyme
#
# An enzyme's activity can also change through phosphorylation, without any
# change in its amount. Phosphoproteomics measures individual sites, so
# `map_modification_sites` adds each site as its own node in the Signaling
# layer, with an edge to the protein it sits on. `reaction_regulation_table`
# then reports a `phospho_axis`: whether a changed site sits on the reaction's
# enzyme.
#
# The table below is made up for illustration: one pyruvate kinase site goes
# up, another does not change, and a site on phosphoenolpyruvate carboxykinase
# goes down. Whether phosphorylation at a given site activates or inhibits the
# enzyme is rarely known, so `phospho_axis_effect` stays 0 unless the edge was
# given a sign.

# %%
sites = pd.DataFrame({
    "protein": ["P53657", "P53657", "Q9Z2V4"],     # Pklr, Pklr, Pck1
    "site": ["S12", "S99", "S458"],
    "log2FC": [1.4, 0.1, -1.1],
    "padj": [0.001, 0.8, 0.01],
})
print(map_modification_sites(graph, sites))

with_sites = reaction_regulation_table(graph)
with_sites[with_sites["n_phosphosites_changed"] > 0][
    ["reaction", "name", "gene_axis", "phospho_axis", "phospho_axis_effect",
     "phospho_axis_via", "metabolite_axis"]]

# %% [markdown]
# ## Do the proteins follow their transcripts?
#
# The enzyme axis combines two steps, transcription and translation.
# `expression_concordance` checks the second step: for every gene measured as
# both transcript and protein, it asks whether the two changed together. The
# categories are `concordant` (both changed the same way), `protein_only`
# (a sign of regulation after transcription), `transcript_only`, `discordant`
# (opposite directions) and `unchanged`.

# %%
concordance = expression_concordance(graph)
print(concordance["counts"])
print(f"correlation of transcript and protein fold changes: "
      f"{concordance['correlation']:.2f}")
concordance["table"][["name", "gene_log2fc", "protein_log2fc", "category"]].head(8)

# %% [markdown]
# `plot_expression_concordance` plots each gene's protein change against its
# transcript change, coloured by category. Points on the diagonal changed
# equally in both layers. In this small example every changed protein follows
# its transcript. In the studies, the off-diagonal points are the interesting
# ones.

# %%
plot_expression_concordance(concordance)
plt.show()
