# %% [markdown]
# # Signed regulatory paths
#
# A change in one molecule reaches other molecules through a chain of
# regulatory steps: a transcription factor changes a gene, the gene changes the
# amount of an enzyme, the enzyme changes how fast a reaction runs, and the
# reaction changes the level of a metabolite. Each step either passes the
# change on in the same direction or reverses it.
#
# Because every edge in a TransNet network records which of the two it does,
# the network can *predict* the direction a metabolite should move, and that
# prediction can be checked against the measurement. This walkthrough shows how
# the prediction is made, how it is scored, and a second method,
# propagation, that asks the same question differently.
#
# | Part | Question | Functions |
# |---|---|---|
# | 1 | How is a direction predicted along a path? | `trace_regulatory_paths` |
# | 2 | Do the predictions match the data? | `path_consistency_summary`, `plot_regulatory_paths` |
# | 3 | Can a signal be followed from a receptor to a metabolite? | Signaling layer, `map_omics_to_network` |
# | 4 | What does propagation add? | `hierarchical_propagation`, `downstream_influence`, `plot_downstream_influence` |
# | 5 | What happens when a layer is missing? | `top_layer_present` |

# %%
import matplotlib.pyplot as plt
import pandas as pd

from transnet import (
    downstream_influence,
    hierarchical_propagation,
    load_example_network,
    load_example_omics,
    load_example_phosphoproteomics,
    map_omics_to_network,
    path_consistency_summary,
    responsive_subnetwork,
    top_layer_present,
    trace_regulatory_paths,
)
from transnet.analysis import versus_chance
from transnet.visualization import (
    plot_downstream_influence,
    plot_regulatory_paths,
    plot_transomic_network,
)

graph = load_example_network()
map_omics_to_network(
    graph, load_example_omics(), id_column="id", log2fc_column="log2FC", qvalue_column="padj",
)


def readable(path):
    """A path of node ids, written with the molecules' names."""
    return " -> ".join(str(graph.nodes[n].get("name", n)) for n in path.split(" -> "))


# %% [markdown]
# ## 1. How signs combine along a path
#
# Each edge carries a sign: +1 if more of the upstream molecule means more of
# the downstream one, and -1 if it means less. The sign of a whole path is the
# product of its edge signs.
#
# For example, ATP is an allosteric inhibitor of pyruvate kinase (sign -1), and
# pyruvate kinase produces pyruvate (sign +1). The path ATP -> pyruvate kinase ->
# pyruvate therefore has sign -1 x +1 = -1: more ATP should mean less pyruvate.
#
# The cell below finds that path. Two default settings are switched off so that
# it appears:
#
# * Path tracing normally skips very common cofactors such as ATP, ADP and
#   water, called *currency metabolites*. They take part in so many reactions
#   that almost any two reactions are connected through them, which produces
#   meaningless paths. Naming ATP as the starting point (`sources`) overrides
#   the rule for that one molecule.
# * Path tracing normally drops paths that pass through a molecule that was
#   measured and did not change (see part 2).
#   `allow_unchanged_intermediates=True` keeps them. We use it here because
#   this cell only illustrates how signs combine, not whether the data support
#   the path.

# %%
from_atp = trace_regulatory_paths(
    graph, sources=["C00002"], target_layer="Metabolome", max_length=3,
    regulated_only=False, allow_unchanged_intermediates=True,
)
for _, row in from_atp[from_atp["edge_types"].str.contains("allosteric_inhibition")].iterrows():
    print(f"{readable(row['path'])}\n    edges: {row['edge_types']}\n"
          f"    path sign: {row['sign']:+d}")

# %% [markdown]
# ## 2. Tracing the data
#
# Called with defaults, `trace_regulatory_paths` starts from every changed
# molecule in the highest layer the network has and follows directed edges down
# to the metabolome. Each row is one path:
#
# | Column | Meaning |
# |---|---|
# | `sign` | product of the edge signs along the path |
# | `source_regulated` | measured direction of the starting molecule (+1 up, -1 down) |
# | `predicted` | `sign` x `source_regulated`: the direction the path predicts for the target |
# | `observed` | measured direction of the target (0 = no significant change) |
# | `consistent` | whether `predicted` equals `observed` |
# | `unsigned_steps` | edges with unknown sign; a path with any is only a tentative prediction |
#
# A decreased enzyme on a +1 path therefore predicts a *decrease*.

# %%
paths = trace_regulatory_paths(graph, target_layer="Metabolome", max_length=6)
print(f"{len(paths)} paths, starting in the {top_layer_present(graph)} layer")
paths.head(6)[["path", "edge_types", "sign", "source_regulated", "predicted",
               "observed", "consistent", "unsigned_steps"]]

# %% [markdown]
# ### Two rules that keep the score meaningful
#
# * **Paths through unchanged molecules are dropped.** If an intermediate
#   molecule was measured and did not change, it cannot have passed a change
#   on, so the data contradict the path. Molecules that were *not measured* are
#   kept, because for them there is no evidence either way.
# * **One verdict per target molecule, not per path.** A metabolite reached
#   by thirty paths that share most of their steps is still one observation.
#   Counting each path would inflate the evidence thirtyfold.
#
# `path_consistency_summary` applies the second rule. It gives one row per
# metabolite: the direction most of its paths predict, whether that matches the
# measurement (`agrees`), and the shortest path that makes the correct
# prediction. `versus_chance` compares the number of correct predictions with
# what guessing up or down at random would achieve.

# %%
verdicts = path_consistency_summary(paths)
tested = verdicts[(verdicts["observed"] != 0) & (verdicts["predicted"] != 0)]
agree = int(tested["agrees"].sum())
print(f"{agree} of {len(tested)} changed metabolites predicted correctly: "
      f"{versus_chance(agree, len(tested))}")
verdicts.assign(target=[graph.nodes[n].get("name", n) for n in verdicts["target"]])[
    ["target", "observed", "predicted", "agrees", "n_paths", "shortest_consistent_path"]]

# %% [markdown]
# The one metabolite predicted wrongly is glucose (C00031). The network predicts
# that more glucose-6-phosphatase makes more glucose, but glucose went down.
#
# `plot_regulatory_paths` draws the best-supported paths, one row per path.
# Circles are molecules and squares are reactions. A molecule's fill shows its
# measured change (red up, blue down, grey unchanged, hollow not measured) and
# its outline shows its layer. Each arrow is labelled with its edge type, and a
# red or blue arrow is a signed regulatory step. The text on the right compares
# the predicted and observed direction of the last molecule.

# %%
plot_regulatory_paths(paths, graph, top_n=6)
plt.show()

# %% [markdown]
# ## 3. From a receptor to a metabolite
#
# The Signaling layer sits at the top of the hierarchy. With it, a single path
# can run from a hormone receptor all the way to a metabolite:
#
# ```
# Insr -> Irs1 -> Pik3r1 -> Akt1 -| Foxo1 -> gene -> enzyme -> REACTION -> metabolite
# ```
#
# The bundled example carries this insulin cascade from KEGG, as
# `phosphorylation` edges (kinase to substrate) and `kinase_tf` edges (kinase to
# transcription factor). Signalling proteins are measured by phosphoproteomics,
# which reports how much of a protein is phosphorylated rather than how much of
# it there is, so their data arrive as a separate table and are mapped onto the
# `Signaling` layer.

# %%
signalling = load_example_network(signaling=True)
tables = dict(load_example_omics())
tables["Signaling"] = load_example_phosphoproteomics()

report = map_omics_to_network(
    signalling, tables, id_column="id", log2fc_column="log2FC",
    qvalue_column="padj", qvalue_threshold=0.05,
)
report.per_layer[["layer", "n_supplied", "n_matched", "n_up", "n_down"]]

# %% [markdown]
# The signalling edges as the network stores them. KEGG annotates the effect of
# a kinase on a transcription factor, so `kinase_tf` edges carry a sign: Akt1
# activates Srebf1 (+1) and inhibits Foxo1 (-1). KEGG rarely records whether a
# phosphorylation activates or inhibits its target, so most `phosphorylation`
# edges have sign 0 (unknown).

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
# Tracing now starts at the receptor. The most common chains of edge types show
# the whole hierarchy in order: phosphorylation down the cascade, a `kinase_tf`
# step into a transcription factor, transcription, translation, catalysis, and
# the reaction that makes the metabolite.

# %%
from_receptor = trace_regulatory_paths(signalling, target_layer="Metabolome", max_length=7)
print(f"{len(from_receptor)} paths, starting in the {top_layer_present(signalling)} layer")
from_receptor["edge_types"].value_counts().head(5).to_frame("paths")

# %% [markdown]
# Which transcription factor do the paths go through? Akt1 has signed edges to
# two factors, Srebf1 and Foxo1, but every path that survives goes through
# Srebf1.
#
# This is the drop rule from part 2 at work. Foxo1 was measured and did not
# change (log2FC -0.11, not significant). A path claiming that Akt1 changed a
# metabolite *through* Foxo1 is therefore contradicted by the data, and is
# removed automatically.

# %%
through_tf = from_receptor[from_receptor["edge_types"].str.contains("kinase_tf")]

factors = set()
for path, kinds in zip(through_tf["path"], through_tf["edge_types"]):
    nodes, steps = path.split(" -> "), kinds.split(" -> ")
    for index, step in enumerate(steps):
        if step == "kinase_tf":                # the factor is the node after this step
            factors.add(signalling.nodes[nodes[index + 1]].get("name", nodes[index + 1]))
print(f"{len(through_tf)} paths through a kinase_tf step; "
      f"factors reached: {', '.join(sorted(factors))}")

foxo1 = signalling.nodes["Q9R1E0"]
print(f"Foxo1: log2FC {foxo1['log2fc']}, regulated {foxo1['regulated']}, "
      f"so paths through it are dropped")

# %% [markdown]
# The same network as a figure, with the Signaling layer as the top row, and
# then the paths from the receptor.

# %%
plot_transomic_network(responsive_subnetwork(signalling),
                       title="Insulin signal, from receptor to metabolite")
plt.show()

# %%
plot_regulatory_paths(from_receptor, signalling, top_n=6)
plt.show()

# %% [markdown]
# ## 4. Propagation: the same question, asked differently
#
# Path tracing lists individual routes. Propagation instead starts from the
# measured changes and pushes them forward along the signed edges, adding up
# everything that arrives at each molecule. A metabolite reached by many weak
# routes and one reached by a single strong route then get different scores,
# which path tracing cannot express.
#
# Use path tracing when you want to name the mechanism behind one change. Use
# propagation when you want a single predicted score for every molecule in a
# layer.
#
# `hierarchical_propagation` returns every molecule the signal reaches. Here
# the seeds are the changed transcripts and proteins, weighted by their fold
# changes.

# %%
seeds = {
    node: float(data["log2fc"])
    for node, data in graph.nodes(data=True)
    if data.get("layer") in ("Transcriptome", "Proteome")
    and data.get("regulated") and data.get("log2fc") is not None
}
scores = hierarchical_propagation(graph, seeds)
scores.groupby("layer")["score"].agg(["count", "mean"]).round(3)

# %% [markdown]
# `downstream_influence` keeps one layer and adds the measured direction, so
# each prediction can be scored (`agrees`). As in part 2, `versus_chance`
# compares the number of correct predictions with random guessing.

# %%
influence = downstream_influence(graph, seeds, target_layer="Metabolome")
measured = influence[influence["observed"].fillna(0) != 0]
agree = int(measured["agrees"].astype("boolean").fillna(False).sum())
print(f"{agree} of {len(measured)} changed metabolites predicted correctly: "
      f"{versus_chance(agree, len(measured))}")
measured[["name", "score", "predicted_direction", "observed", "agrees"]]

# %% [markdown]
# `plot_downstream_influence` plots each metabolite's predicted score against
# its measured direction, so the metabolites the network gets wrong stand out
# individually.

# %%
plot_downstream_influence(influence)
plt.show()

# %% [markdown]
# Propagation works from the signalling layer too. Seeding with the measured
# phosphorylation changes shows how far the insulin signal is predicted to
# reach in each layer.

# %%
receptor_seeds = {
    node: float(data["log2fc"])
    for node, data in signalling.nodes(data=True)
    if data.get("layer") == "Signaling" and data.get("regulated")
    and data.get("log2fc") is not None
}
arrived = hierarchical_propagation(signalling, receptor_seeds)
arrived.groupby("layer")["score"].agg(["count", "mean"]).round(3)

# %% [markdown]
# ## 5. When a layer is missing
#
# No analysis requires a particular layer. Path tracing starts at the highest
# layer the network contains and reports which one that was. With the Signaling
# layer, paths start at the receptor; without it, they start at the proteome.

# %%
with_phospho = dict(load_example_omics())
with_phospho["Signaling"] = load_example_phosphoproteomics()

for label, network, data in [
        ("without signalling", load_example_network(), load_example_omics()),
        ("with signalling", load_example_network(signaling=True), with_phospho)]:
    map_omics_to_network(network, data, id_column="id",
                         log2fc_column="log2FC", qvalue_column="padj")
    traced = trace_regulatory_paths(network, target_layer="Metabolome", max_length=6)
    print(f"{label:<19} starts at {top_layer_present(network):<10} {len(traced):>4} paths")
