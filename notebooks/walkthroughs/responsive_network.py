# %% [markdown]
# # Mapping data onto the network
#
# Each omics table identifies its molecules in its own way: Entrez or Ensembl
# gene ids, UniProt accessions, KEGG compound ids or plain names. Mapping the
# tables onto the network is where data are most easily lost without anyone
# noticing, so `map_omics_to_network` returns a report of how much of each
# table it could place.

# %%
import pandas as pd

from transnet import (
    layer_coverage,
    load_example_network,
    load_example_omics,
    map_omics_to_network,
    regulated_nodes,
    responsive_subnetwork,
)

graph = load_example_network()
tables = load_example_omics("insulin_sensitive")
{layer: frame.shape for layer, frame in tables.items()}

# %%
report = map_omics_to_network(
    graph, tables, id_column="id", log2fc_column="log2FC", qvalue_column="padj",
    qvalue_threshold=0.05,
)
report.per_layer

# %% [markdown]
# `match_fraction` is the share of a table's rows that found a node. A value
# below about 0.5 usually means the table uses a different kind of identifier
# from the network, for example Ensembl ids against a network keyed by Entrez
# ids. `unmatched` lists the rows that could not be placed, and `match_rate` is
# the share over all tables.

# %%
print(f"overall match rate {report.match_rate:.0%}")
for layer, missing in report.unmatched.items():
    if missing:
        print(f"{layer}: {len(missing)} unmatched, e.g. {', '.join(missing[:4])}")

# %% [markdown]
# `layer_coverage` asks the same question from the network's side: of the nodes
# in each layer, how many were measured, and how many of those changed. The
# metabolome is usually the limiting layer.

# %%
layer_coverage(graph)

# %% [markdown]
# `plot_layer_changes` draws what an analysis of each layer on its own would
# report: the number of molecules up and down per layer, and nothing about how
# they are connected.

# %%
from transnet.visualization import plot_layer_changes

coverage = layer_coverage(graph)
plot_layer_changes(
    coverage.assign(group="insulin-sensitive")[["group", "layer", "n_up", "n_down"]],
    title="Changed molecules per layer",
)

# %% [markdown]
# ## The responsive subnetwork
#
# `responsive_subnetwork` keeps the molecules that changed and the edges
# between them. This is the part of the network that responded, and it is what
# most trans-omic analyses and figures work on.

# %%
responsive = responsive_subnetwork(graph)
print(f"{responsive.number_of_nodes()} molecules, {responsive.number_of_edges()} relationships")
pd.Series({d["edge_type"]: 1 for _, _, d in responsive.edges(data=True)}).index.tolist()

# %% [markdown]
# Reactions are kept even though a reaction is never measured itself. Without
# them an enzyme and the metabolites of its reaction would not be connected, and
# nothing could be said across layers. `regulated_nodes` lists the changed
# molecules.

# %%
regulated = regulated_nodes(graph)
pd.DataFrame(
    [{"node": n, "name": graph.nodes[n].get("name"), "layer": graph.nodes[n].get("layer"),
      "log2FC": graph.nodes[n].get("log2fc")} for n in regulated[:10]]
)

# %% [markdown]
# ## Responsive, with or without a direction
#
# Each mapped node has a `regulated` attribute: +1 (up), -1 (down), or 0. A
# molecule can be significant without a usable fold change, for example from a
# test across several time points; it then has `regulated` 0 but still counts
# as responsive. `is_responsive` checks for either case, and it is the test
# `responsive_subnetwork` uses, so such a molecule still joins its reaction.

# %%
from transnet import is_responsive

states = pd.Series(
    [graph.nodes[n].get("regulated") for n in graph if is_responsive(graph.nodes[n])]
).value_counts().sort_index()
states.rename("nodes").to_frame().rename_axis("regulated")

# %% [markdown]
# ## What it looks like

# %%
import matplotlib.pyplot as plt

from transnet.visualization import plot_transomic_network

figure = plot_transomic_network(responsive, title="Responsive trans-omic network")
plt.show()
