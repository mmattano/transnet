# %% [markdown]
# # Hubs across layers, and timing
#
# Two readings that need the whole hierarchy: which molecules *join* the
# layers, and whether the wiring explains the order in which things respond.

# %%
import pandas as pd

from transnet import (
    assign_temporal_parameters,
    load_example_network,
    load_example_omics,
    load_example_timecourse,
    map_omics_to_network,
    split_by_response_class,
    temporal_network_structure,
    transomic_hubs,
)
from transnet.analysis import compute_centrality_measures

graph = load_example_network()
map_omics_to_network(
    graph, load_example_omics(), id_column="id", log2fc_column="log2FC", qvalue_column="padj",
)

# %% [markdown]
# ## A hub within a layer is not a trans-omic hub
#
# Degree centrality finds whatever is best annotated: a spliceosome protein,
# a ribosomal subunit. Cross-layer degree finds the molecules through which
# one layer's response reaches another.

# %%
hubs = transomic_hubs(graph, top_percent=10)
hubs.head(10)[["node", "name", "layer", "degree", "cross_layer_degree",
               "n_layers_touched", "versatility"]]

# %%
centrality = compute_centrality_measures(graph, top_n=5)
pd.DataFrame({
    "by cross-layer degree": [graph.nodes[n].get("name", n) for n in hubs["node"].head(5)],
    "by degree centrality": [graph.nodes[n].get("name", n)
                             for n, _ in centrality["degree_centrality"]],
    "by betweenness": [graph.nodes[n].get("name", n)
                       for n, _ in centrality["betweenness_centrality"]],
})

# %% [markdown]
# ## Response times on the network
#
# A time course gives each molecule a half-response time,
# which is the time it takes for the molecule to reach half of its maximum response.
# As is presented in Morita et al. (2025). The question is
# whether the best-connected molecules respond *fastest*.

# %%
TIMEPOINTS = [0, 5, 15, 30, 60]          # minutes after the glucose bolus

parameters = assign_temporal_parameters(
    graph, {"Metabolome": load_example_timecourse()},
    time_columns=TIMEPOINTS, id_column="id",
)
ranked = parameters.dropna(subset=["t_half"]).sort_values("t_half")
ranked = ranked.assign(name=[graph.nodes[n].get("name", n) for n in ranked["node"]])
pd.concat([ranked.head(3), ranked.tail(3)])[["node", "name", "layer", "t_half", "degree"]]

# %%
structure = temporal_network_structure(graph)
timing = structure["degree_vs_thalf"]
print(timing["interpretation"])
structure["per_layer_thalf"]

# %% [markdown]
# The interpretation states a direction only when the correlation is
# significant. On this small example it is not, and the notebook says so.
# a reported trend that a test does not support is the failure mode this
# wording exists to prevent.

# %%
split = split_by_response_class(graph)
pd.DataFrame([{"class": name, "molecules": sub.number_of_nodes()}
              for name, sub in split["subnetworks"].items()])

# %%
import matplotlib.pyplot as plt

from transnet.visualization import plot_temporal_structure, plot_transomic_hubs

plot_transomic_hubs(hubs)
plot_temporal_structure(graph, structure, time_unit="min")
plt.show()
