# %% [markdown]
# # Hubs across layers, and timing
#
# Two questions that need the whole network: which molecules *connect* the
# layers, and whether the network's structure explains the order in which
# molecules respond over time.

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
# Ranking molecules by their number of connections finds whatever is best
# annotated, such as spliceosome or ribosomal proteins with hundreds of
# interaction partners. `transomic_hubs` instead counts connections to *other
# layers* (`cross_layer_degree`) and how many layers a molecule touches. These
# are the molecules through which one layer's response reaches another.
# `versatility` combines both into one score.

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
# A time course gives each molecule a half-response time (`t_half`): the time
# it takes to reach half of its largest change. `assign_temporal_parameters`
# computes it for every measured molecule and stores it on the node. The
# question is then whether the best-connected molecules respond fastest.

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
# `temporal_network_structure` correlates each molecule's number of
# connections with its half-response time, and compares the layers. The
# `interpretation` text states a direction only when the correlation is
# statistically significant. On this small example it is not, so no trend is
# claimed.
#
# `split_by_response_class` divides the molecules into fast and slow
# responders (and, when dose data are present, sensitive and insensitive ones),
# split at the median by default. It returns the subnetwork of each class, so
# the groups can be analysed separately.

# %%
split = split_by_response_class(graph)
pd.DataFrame([{"class": name, "molecules": sub.number_of_nodes()}
              for name, sub in split["subnetworks"].items()])

# %% [markdown]
# `plot_transomic_hubs` ranks the hubs by cross-layer connections, coloured by
# layer. `plot_temporal_structure` plots half-response time against the number
# of connections for each layer, with the result of the correlation test.

# %%
import matplotlib.pyplot as plt

from transnet.visualization import plot_temporal_structure, plot_transomic_hubs

plot_transomic_hubs(hubs)
plt.show()

# %%
plot_temporal_structure(graph, structure, time_unit="min")
plt.show()
