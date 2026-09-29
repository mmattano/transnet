# %% [markdown]
# # Signed regulatory paths
#
# A stimulus reaches a metabolite through a chain of
# regulatory steps. The signal transduction is testable:
# in a given network, multiply the signs along the
# chain, multiply by the direction the source moved, and
# then compare it with the
# metabolite's measurement.

# %%
import matplotlib.pyplot as plt
import pandas as pd

from transnet import (
    load_example_network,
    load_example_omics,
    map_omics_to_network,
    path_consistency_summary,
    responsive_subnetwork,
    top_layer_present,
    trace_regulatory_paths,
)
from transnet.analysis import versus_chance

graph = load_example_network()
map_omics_to_network(
    graph, load_example_omics(), id_column="id", log2fc_column="log2FC", qvalue_column="padj",
)

paths = trace_regulatory_paths(graph, target_layer="Metabolome", max_length=6)
print(f"traced from {top_layer_present(graph)}, the highest layer present")
paths.head(6)[["path", "edge_types", "sign", "source_regulated", "predicted", "observed",
               "consistent"]]

# %% [markdown]
# ## Two checks on the paths
#
# * A path through a molecule that **was measured and did not change** is
#   dropped: the data contradict it. Unmeasured intermediates are kept, since
#   unmeasured is unknown, not unchanged.
# * The rate is quoted **per target molecule**, not per path. One metabolite
#   reached by thirty overlapping paths is one piece of evidence, not thirty.

# %%
verdicts = path_consistency_summary(paths)
tested = verdicts[(verdicts["observed"] != 0) & (verdicts["predicted"] != 0)]
agree = int(tested["agrees"].sum())
print(f"{agree} of {len(tested)} changed metabolites predicted correctly: "
      f"{versus_chance(agree, len(tested))}")
verdicts[["target", "observed", "predicted", "agrees", "n_paths", "shortest_consistent_path"]]

# %% [markdown]
# ## Sign propagation, step by step
#
# ATP inhibits pyruvate kinase, so a rise in ATP predicts a fall in pyruvate.
# ATP is a currency metabolite, which tracing routes around by default; naming
# it as the source keeps it. This is a structural demonstration, so
# intermediates that did not change are kept too.

# %%
inhibited = trace_regulatory_paths(
    graph, sources=["C00002"], target_layer="Metabolome", max_length=3,
    regulated_only=False, allow_unchanged_intermediates=True,
)
for _, row in inhibited[inhibited["edge_types"].str.contains("allosteric_inhibition")].head(3).iterrows():
    names = " -> ".join(str(graph.nodes[n].get("name", n)) for n in row["path"].split(" -> "))
    print(f"{names}\n    {row['edge_types']}  predicts {row['sign']:+d}")

# %% [markdown]
# ## Layer-optional, demonstrated
#
# No analysis requires a particular layer. With a Signaling layer present,
# paths start there; without one they start at the highest layer that exists,
# and the result says which.

# %%
for label, network in [("no signaling", load_example_network()),
                       ("with signaling", load_example_network(signaling=True))]:
    map_omics_to_network(network, load_example_omics(), id_column="id",
                         log2fc_column="log2FC", qvalue_column="padj")
    traced = trace_regulatory_paths(network, target_layer="Metabolome", max_length=6)
    print(f"{label:<16} starts at {top_layer_present(network):<13} {len(traced):>4} paths")

# %% [markdown]
# ## Signal flow, from receptor to metabolite
#
# The Signaling layer is the top of the hierarchy, and with it the network spans
# every step from a receptor to a metabolite concentration:
#
# ```
# Insr -> Irs1 -> Pik3r1 -> Akt1 -| Foxo1 -> gene -> enzyme -> REACTION -> metabolite
# ```
#
# The bundled example carries that cascade as `phosphorylation` and `kinase_tf`
# edges from KEGG's insulin pathway, and a phosphoproteomics table to map onto
# it. Signalling molecules are measured by phosphorylation rather than by
# abundance, which is why they arrive in their own table.

# %%
from transnet import load_example_phosphoproteomics

signalling = load_example_network(signaling=True)
phospho = load_example_phosphoproteomics()
tables = dict(load_example_omics())
tables["Signaling"] = phospho

report = map_omics_to_network(
    signalling, tables, id_column="id", log2fc_column="log2FC",
    qvalue_column="padj", qvalue_threshold=0.05,
)
report.per_layer[["layer", "n_supplied", "n_matched", "n_up", "n_down"]]

# %% [markdown]
# The cascade as the network holds it. `kinase_tf` edges carry a sign because
# KEGG annotates the effect: Akt1 activating Srebf1 is +1, Akt1 inhibiting Foxo1
# is -1. The unsigned `phosphorylation` edges up the chain are the usual case,
# since a site's effect on activity is rarely recorded.

# %%
pd.DataFrame([
    {"from": signalling.nodes[u].get("symbol", u),
     "to": signalling.nodes[v].get("symbol", v),
     "edge_type": d["edge_type"], "sign": d.get("sign"),
     "from log2FC": signalling.nodes[u].get("log2fc")}
    for u, v, d in signalling.edges(data=True)
    if d["edge_type"] in ("phosphorylation", "kinase_tf")
])

# %% [markdown]
# Now trace from the receptor. The chain of edge types is the hierarchy itself,
# end to end: phosphorylation up the cascade, `kinase_tf` into a transcription
# factor, transcription, translation, catalysis, and finally the reaction that
# changes a metabolite.

# %%
from_receptor = trace_regulatory_paths(signalling, target_layer="Metabolome",
                                       max_length=7)
print(f"{len(from_receptor)} paths from {top_layer_present(signalling)}")
from_receptor["edge_types"].value_counts().head(5).to_frame("paths")

# %% [markdown]
# Fourteen of the twenty-six pass through a `kinase_tf` step, and all of them
# reach the same transcription factor. None goes through Foxo1, even though Akt1
# inhibits it and that edge is in the network.
#
# The reason is the drop rule from the top of this notebook. Foxo1 was measured
# and did not change (log2FC -0.11, not significant), so a path claiming Akt1
# reached a metabolite through Foxo1 is contradicted by the data and is dropped.
# The one signed branch of this cascade that the data refuse is removed without
# anyone having to notice it by hand.

# %%
through_tf = from_receptor[from_receptor["edge_types"].str.contains("kinase_tf")]

# The factor is the target of the kinase_tf step, which is the node one place
# after the kinase in the path.
factors = set()
for path, kinds in zip(through_tf["path"], through_tf["edge_types"]):
    nodes, steps = path.split(" -> "), kinds.split(" -> ")
    for index, step in enumerate(steps):
        if step == "kinase_tf":
            factors.add(signalling.nodes[nodes[index + 1]].get("name",
                                                               nodes[index + 1]))
print(f"{len(through_tf)} paths through a kinase_tf step; "
      f"factor reached: {', '.join(sorted(factors))}")

foxo1 = signalling.nodes["Q9R1E0"]
print(f"Foxo1: log2FC {foxo1['log2fc']}, regulated {foxo1['regulated']} "
      f"-> paths through it are dropped")
through_tf.head(4)[["path", "sign", "predicted", "observed", "consistent"]]

# %% [markdown]
# And the same question by propagation: push the measured phosphorylation
# forward from the receptor and read what arrives at each layer. This is the
# reading a phosphoproteomics study wants, and it needs the signed Signaling
# edges to mean anything.

# %%
from transnet import hierarchical_propagation

receptor_seeds = {
    node: float(data["log2fc"])
    for node, data in signalling.nodes(data=True)
    if data.get("layer") == "Signaling" and data.get("regulated")
    and data.get("log2fc") is not None
}
arrived = hierarchical_propagation(signalling, receptor_seeds)
arrived.groupby("layer")["score"].agg(["count", "mean"]).round(3)

# %%
arrived[arrived["layer"] == "Metabolome"].head(8)[
    ["name", "score", "predicted_direction"]]

# %% [markdown]
# The same cascade as a figure. With a Signaling layer present the drawing spans
# five rows, and the receptor sits at the top of the hierarchy rather than the
# transcriptome.

# %%
from transnet.visualization import plot_regulatory_paths, plot_transomic_network

plot_transomic_network(
    responsive_subnetwork(signalling),
    title="Insulin signal, from receptor to metabolite",
)
plt.show()

# %% [markdown]
# And the paths themselves, drawn from the receptor down.

# %%
plot_regulatory_paths(from_receptor, signalling, top_n=6)
plt.show()

# %%
from transnet.visualization import plot_regulatory_paths

plot_regulatory_paths(paths, graph, top_n=6)
plt.show()

# %% [markdown]
# ## Propagation, the other way to ask
#
# Path tracing enumerates routes to a target. Propagation pushes the changed
# molecules outwards along signed edges and reads off where the signal arrives,
# so a target reached by many weak routes and one reached by a single strong
# route are scored differently.
#
# `hierarchical_propagation` returns every node it reaches.

# %%
from transnet import downstream_influence, hierarchical_propagation

seeds = {
    node: float(data["log2fc"])
    for node, data in graph.nodes(data=True)
    if data.get("layer") in ("Transcriptome", "Proteome")
    and data.get("regulated") and data.get("log2fc") is not None
}
scores = hierarchical_propagation(graph, seeds)
scores.head(8)[["node", "layer", "name", "score", "predicted_direction"]]

# %% [markdown]
# `downstream_influence` restricts that to one layer and adds the measured
# direction, so the prediction can be scored. `agrees` is the comparison, and
# `versus_chance` states it against a coin flip.

# %%
influence = downstream_influence(graph, seeds, target_layer="Metabolome")
measured = influence[influence["observed"].fillna(0) != 0]
agree = int(measured["agrees"].astype("boolean").fillna(False).sum())
print(f"{agree} of {len(measured)} changed metabolites predicted correctly: "
      f"{versus_chance(agree, len(measured))}")
measured[["name", "score", "predicted_direction", "observed", "agrees"]]

# %% [markdown]
# The figure puts the predicted score against the measured direction for each
# metabolite, so the ones the hierarchy gets wrong are visible individually.

# %%
from transnet.visualization import plot_downstream_influence

plot_downstream_influence(influence)
plt.show()
