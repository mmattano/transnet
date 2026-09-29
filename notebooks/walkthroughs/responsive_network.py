# %% [markdown]
# # Mapping data onto the network
#
# Measurements arrive keyed by whatever the instrument produced. Mapping them
# onto the network is where studies lose data without noticing, so
# `map_omics_to_network` returns a coverage report rather than mapping quietly.

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
# A `match_fraction` below about 0.5 usually means an identifier-type mismatch,
# such as Ensembl ids against a network keyed by Entrez. The report names the
# features it could not place, and `match_rate` summarises the whole mapping.

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
# `plot_layer_changes` draws what a per-layer analysis would report: counts up
# and down per layer, and nothing about how they connect.

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
# The molecules that changed, plus the relationships that join them. This is
# the object the trans-omics literature works on: the part of the hierarchy
# that moved, with its wiring.

# %%
responsive = responsive_subnetwork(graph)
print(f"{responsive.number_of_nodes()} molecules, {responsive.number_of_edges()} relationships")
pd.Series({d["edge_type"]: 1 for _, _, d in responsive.edges(data=True)}).index.tolist()

# %% [markdown]
# Reactions are kept as connectors even though a reaction is never itself
# measured. Without them the enzyme and the metabolite it acts on fall into
# separate components and no cross-layer statement is possible.

# %%
regulated = regulated_nodes(graph)
pd.DataFrame(
    [{"node": n, "name": graph.nodes[n].get("name"), "layer": graph.nodes[n].get("layer"),
      "log2FC": graph.nodes[n].get("log2fc")} for n in regulated[:10]]
)

# %% [markdown]
# ## Responsive, with or without a direction
#
# `regulated` carries a direction: +1, -1, or 0 for a molecule that changed
# without a usable fold change. `is_responsive` asks the weaker question, which
# is the one the subnetwork uses, so a metabolite measured as changed but
# without a sign still joins its reaction.

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
