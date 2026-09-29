# %% [markdown]
# # The trans-omic network
#
# A trans-omic network is a typed, directed, signed regulatory hierarchy:
#
# ```
# signal -> TF -> gene -> enzyme protein -> REACTION <- metabolite
# ```
#
# Every edge records which relationship it represents, which way the regulation
# runs, and what evidence put it there. This walkthrough loads the bundled
# example network and reads those properties off it.

# %%
import pandas as pd

from transnet import (
    EDGE_TYPES,
    LAYERS,
    available_edge_types,
    available_layers,
    cross_layer_connectivity,
    load_example_network,
    to_simple_graph,
)
from transnet.analysis import compute_network_statistics

graph = load_example_network()
graph

# %% [markdown]
# ## The schema
#
# `LAYERS` and `EDGE_TYPES` describe what the package can represent, before any
# data is involved. Each edge type declares the layers it runs between, its
# default sign and the database it comes from, so a network can be queried for
# its own structure rather than documented separately.

# %%
pd.DataFrame(
    [
        {"edge_type": name, "from": kind.source_layer, "to": kind.target_layer,
         "default sign": kind.sign, "source": kind.source_db}
        for name, kind in EDGE_TYPES.items()
    ]
).set_index("edge_type")

# %% [markdown]
# `LAYERS` is the full set; `available_layers` reports the ones a given network
# has. Analyses work on whichever layers a study measured, so a study without a
# proteome still gets an enzyme axis, inferred from transcripts.

# %%
print("schema:", ", ".join(LAYERS))
print("this network:", " -> ".join(available_layers(graph)))

# %% [markdown]
# Two orderings exist and they are not the same. `LAYER_HIERARCHY` is the
# direction regulation flows, which is what path tracing and propagation follow.
# `DISPLAY_ORDER` is the order figures stack the layers.

# %%
from transnet import DISPLAY_ORDER, LAYER_HIERARCHY

print("regulation flows:", " -> ".join(LAYER_HIERARCHY))
print("figures stack:   ", " -> ".join(DISPLAY_ORDER))

pd.Series(available_edge_types(graph), name="edges").sort_values(ascending=False).to_frame()

# %% [markdown]
# ## Signs
#
# Each relationship carries a sign. `allosteric_inhibition` is `-1`, so a rise
# in that metabolite predicts a fall in the reaction it regulates. Signed path
# tracing depends on this.

# %%
pd.DataFrame(
    [
        {"source": u, "target": v, "edge_type": d["edge_type"], "sign": d.get("sign")}
        for u, v, d in graph.edges(data=True)
        if d["edge_type"] in {"allosteric_activation", "allosteric_inhibition", "catalysis"}
    ]
).head(8)

# %% [markdown]
# ## Direction and parallel edges
#
# The network is a `MultiDiGraph`. A metabolite can be both a substrate of a
# reaction and its allosteric regulator. Those are two different claims, stored
# as two different edges.

# %%
for u, v in [("C00668", "R00299"), ("R00299", "C00668")]:
    records = graph.get_edge_data(u, v) or {}
    print(f"{u} -> {v}: " + (", ".join(r["edge_type"] for r in records.values()) or "nothing"))

# %% [markdown]
# Some graph algorithms are defined only for simple undirected graphs.
# `to_simple_graph` collapses the multigraph for those, keeping the node
# attributes. The typed edges are what the trans-omic analyses need, so they
# never use it.

# %%
simple = to_simple_graph(graph)
print(f"{graph.number_of_edges()} typed edges collapse to "
      f"{simple.number_of_edges()} undirected connections")

# %% [markdown]
# ## Cross-layer connectivity
#
# A network whose edges sit inside single layers is a stack of separate
# single-omics networks. This fraction measures how much of it crosses.

# %%
connectivity = cross_layer_connectivity(graph)
print(f"{connectivity['cross_layer_fraction']:.0%} of edges cross between layers")
connectivity["matrix"]

# %% [markdown]
# The same matrix as a figure. Off-diagonal cells are the trans-omic edges; the
# diagonal is within-layer wiring, which for a real organism is dominated by the
# protein interactome.

# %%
from transnet.visualization import plot_layer_connectivity

plot_layer_connectivity(connectivity, title="Example network: edges between layers")

# %% [markdown]
# ## The ordinary census
#
# Size, density and components, answered here so the trans-omic readings that
# follow have a baseline. `plot_network_metrics` draws the degree distribution
# and component sizes beside it.

# %%
statistics = compute_network_statistics(graph)
{k: statistics[k] for k in ("nodes", "edges", "density", "connected_components",
                            "largest_component_ratio", "avg_degree", "max_degree")}

# %%
from transnet.visualization import plot_network_metrics

plot_network_metrics(graph)

# %% [markdown]
# ## Writing it out and reading it back
#
# A network is saved as two CSVs, `nodes.csv` and `interactions.csv`. The
# studies read a prebuilt network rather than rebuilding it, so the round trip
# has to preserve everything an analysis reads: layers, node types, gene
# symbols, edge types, signs and stoichiometry.

# %%
import tempfile
from pathlib import Path

from transnet.io import read_network, write_network

with tempfile.TemporaryDirectory() as directory:
    write_network(graph, directory)
    print(", ".join(sorted(p.name for p in Path(directory).iterdir())))
    again = read_network(str(Path(directory) / "interactions.csv"),
                         nodes_file=str(Path(directory) / "nodes.csv"))

print(f"{graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges "
      f"-> {again.number_of_nodes()} nodes, {again.number_of_edges()} edges")
print("edge types preserved:",
      available_edge_types(graph) == available_edge_types(again))

# %% [markdown]
# ## The object API
#
# The graph is the working representation, and `Transnet` is the builder that
# produces it: layers hold typed elements, and `generate_graph` turns them into
# the `MultiDiGraph` every analysis takes. `maintenance/build_networks.py` uses
# this side. A study that loads a prebuilt network does not need it.

# %%
from transnet import Metabolite, Protein, Reaction, Transnet
from transnet.biology.layers import Metabolome, Proteome, Reactions

network = Transnet(name="one reaction")

enzyme = Protein()
enzyme.uniprot_id = "P52789"
enzyme.name = "Hexokinase-2"
enzyme.gene = ["Hk2"]
enzyme.ec_number = ["2.7.1.1"]
proteome = Proteome()
proteome.proteins = [enzyme]
network.proteome = proteome

glucose = Metabolite()
glucose.kegg_compound_id = "C00031"
glucose.kegg_name = "D-Glucose"
glucose.pubchem_id = "5793"
g6p = Metabolite()
g6p.kegg_compound_id = "C00092"
g6p.kegg_name = "D-Glucose 6-phosphate"
g6p.pubchem_id = "5958"
metabolome = Metabolome()
metabolome.metabolites = [glucose, g6p]
network.metabolome = metabolome

reaction = Reaction(
    id="R00299", name="ATP:D-glucose 6-phosphotransferase",
    enzyme=["2.7.1.1"], substrates=["C00031"], products=["C00092"],
    reversible=False,
)
reactions = Reactions()
reactions.reactions = [reaction]
network.reactions = reactions

built = network.generate_graph()
print(f"{built.number_of_nodes()} nodes, {built.number_of_edges()} edges, "
      f"layers {', '.join(available_layers(built))}")
pd.DataFrame([{"source": u, "target": v, "edge_type": d["edge_type"]}
              for u, v, d in built.edges(data=True)])

# %% [markdown]
# The Signaling layer is the one layer the bundled network leaves out by default,
# because most studies do not measure a phosphoproteome. Asking for it adds the
# insulin cascade as `SignalingProtein` nodes, and every analysis then starts one
# layer higher.

# %%
from transnet import SignalingProtein

with_signalling = load_example_network(signaling=True)
print(" -> ".join(available_layers(with_signalling)))
pd.DataFrame([
    {"node": n, "name": d.get("name"), "node_type": d.get("node_type")}
    for n, d in with_signalling.nodes(data=True) if d.get("layer") == "Signaling"
])

# %%
kinase = SignalingProtein(id="P31750", name="Akt1", uniprot_id="P31750", sign=1)
print(kinase)

# %% [markdown]
# ## Identifier maps
#
# Metabolomics reports PubChem or plain names while the network holds KEGG
# compound identifiers, so a table keyed the other way matches nothing.
# `build_alias_id_map` reads the translation out of the layer itself, with no
# database call, and the result goes straight to `map_omics_to_network`.

# %%
from transnet import build_alias_id_map

print(build_alias_id_map(network, "Metabolome", "pubchem_id"))
