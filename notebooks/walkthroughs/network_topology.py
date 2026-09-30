# %% [markdown]
# # The network as a graph
#
# The previous walkthroughs read the network as a regulatory hierarchy, where
# edges have a direction and a sign. This one uses standard graph methods
# instead: summary statistics, centrality, communities, active modules,
# diffusion, recurring wiring patterns (motifs), and weak points. Each section
# says what the method answers and how to read its output on a trans-omic
# network.

# %%
import matplotlib.pyplot as plt
import pandas as pd

from transnet import (
    convergence_significance,
    load_example_network,
    load_example_omics,
    load_example_pathways,
    map_omics_to_network,
    regulated_nodes,
    regulatory_motifs,
    structural_vulnerability,
    transomic_hubs,
)
from transnet.analysis import (
    compute_centrality_measures,
    compute_network_statistics,
    detect_communities,
    enrichment_analysis,
    find_active_modules,
    identify_hubs,
    network_enrichment_permutation,
    random_walk_with_restart,
)
from transnet.visualization import (
    plot_community_network,
    plot_convergence_null,
    plot_network_metrics,
    plot_regulatory_motifs,
    plot_structural_vulnerability,
    transomic_backbone,
)

graph = load_example_network()
map_omics_to_network(
    graph, load_example_omics(), id_column="id", log2fc_column="log2FC", qvalue_column="padj",
)


def name(node, limit=40):
    return str(graph.nodes[node].get("name", node))[:limit]


statistics = compute_network_statistics(graph)
{k: statistics[k] for k in ("nodes", "edges", "connected_components",
                            "largest_component_ratio", "avg_degree", "avg_clustering")}

# %% [markdown]
# Average path length and diameter take time proportional to the square of the
# number of nodes. They are computed only when the largest component is small,
# and reported as `None` otherwise, because on an organism-wide network they
# would take hours.
#
# `plot_network_metrics` draws three panels: the degree distribution (how many
# molecules have 1, 2, 3 ... connections), the sizes of the connected
# components, and the number of molecules per layer.

# %%
plot_network_metrics(graph)
plt.show()

# %% [markdown]
# ## Centrality is not the same as connecting layers
#
# Degree centrality ranks molecules by their number of connections, and
# betweenness by how many shortest paths run through them. `transomic_hubs`
# ranks by something else: how many connections a molecule has to *other
# layers*. On an organism-wide network the first two favour well-annotated
# proteins such as ribosomal subunits, while the third finds the molecules
# through which one layer's response reaches the next.

# %%
centrality = compute_centrality_measures(graph, top_n=5)
hubs = transomic_hubs(graph, top_percent=10)

pd.DataFrame({
    "degree centrality": [name(n) for n, _ in centrality["degree_centrality"]],
    "betweenness": [name(n) for n, _ in centrality["betweenness_centrality"]],
    "most connections": [name(n) for n, _ in identify_hubs(graph, top_n=5)],
    "most cross-layer connections": [name(n) for n in hubs["node"].head(5)],
})

# %% [markdown]
# ## Do communities span layers?
#
# `detect_communities` splits the network into groups of nodes that are more
# densely connected to each other than to the rest (Louvain method by default).
# A community that contains genes, proteins, reactions *and* metabolites is a
# functional unit across layers. A community confined to one layer usually
# reflects how that layer was annotated, for example a cluster of protein
# interactions.

# %%
communities = detect_communities(graph, method="louvain")
membership = pd.Series(communities, name="community")
layers = pd.Series({n: graph.nodes[n].get("layer") for n in communities}, name="layer")
composition = pd.crosstab(membership, layers)
print(f"{len(composition)} communities, "
      f"{int((composition > 0).sum(axis=1).gt(1).sum())} spanning more than one layer")
composition.head(8)

# %% [markdown]
# `plot_community_network` shows the layer make-up of the largest communities
# (left) and draws one of them (right). In the drawing, the fill colour is the
# layer and squares are reactions. A red or blue ring marks a molecule that went
# up or down. Solid, dashed and dotted lines are mass flow, catalysis or
# regulation, and allosteric or interaction edges. The best-connected and the
# changed molecules are labelled.

# %%
plot_community_network(graph, communities)
plt.show()

# %% [markdown]
# ## Active modules
#
# `responsive_subnetwork` keeps every molecule that passes the significance
# cut-off, whether or not it is connected to the others. An active module is
# built differently. Each molecule gets a score for how strongly it changed.
# The module starts at the strongest change and grows one neighbour at a time,
# staying connected, for as long as each addition raises the module's total
# score. Reactions have no measurement of their own, so a reaction is added
# together with the molecule behind it. This lets a module run from an enzyme,
# through its reaction, to a metabolite. The result is a small, connected set of
# molecules that changed together, which is easier to interpret than a
# scattered list.

# %%
modules = find_active_modules(graph, score_attr="log2fc", n_modules=3)
for i, module in enumerate(modules, start=1):
    layer_counts = pd.Series([graph.nodes[n].get("layer") for n in module]).value_counts()
    print(f"module {i}: {len(module)} molecules "
          f"({', '.join(f'{c} {l}' for l, c in layer_counts.items())})")
    print("   ", ", ".join(sorted({name(n, 22) for n in module})[:12]), "...")

# %% [markdown]
# To see what a module is about, test which pathways are over-represented
# among its reactions. `enrichment_analysis` runs a hypergeometric test of each
# pathway's reactions in the module against all reactions in the network, and
# reports the p-value and the false discovery rate (FDR). In this small
# example the module takes in most of the network, so no pathway stands out.
# On an organism-wide network a module is a small fraction of the whole, and
# this test names what it is about.

# %%
reactions = [n for n, d in graph.nodes(data=True) if d.get("layer") == "Reactions"]
in_module = [n for n in modules[0] if n in reactions]
pathway_members = {}
for reaction, names in load_example_pathways().items():
    for pathway in names:
        pathway_members.setdefault(pathway, []).append(reaction)
enrichment_analysis(in_module, reactions, pathway_members)

# %% [markdown]
# ## Motifs: recurring wiring patterns
#
# Some regulatory patterns can be recognised from the signed wiring alone,
# before looking at any data. `regulatory_motifs` searches for three:
#
# * **product inhibition**: a reaction's product inhibits the same reaction,
#   the simplest negative feedback;
# * **allosteric feedback**: a metabolite made by one reaction regulates
#   another reaction that feeds it;
# * **feed-forward**: a transcription factor controls a gene whose enzyme
#   catalyses a reaction that the factor's other targets also affect.
#
# By default only molecules that changed are considered.

# %%
motifs = regulatory_motifs(graph)
print(motifs["counts"])
plot_regulatory_motifs(motifs)
plt.show()
motifs["motifs"][["motif", "reaction_name", "metabolite_name", "enzyme", "sign_product"]]

# %% [markdown]
# Product inhibition is one explanation for a *controversial* reaction in the
# regulation-axis analysis: the enzyme goes up, makes more product, and the
# product slows the reaction down again. Hexokinase, inhibited by its product
# glucose-6-phosphate, is both a motif here and a controversial reaction in the
# reaction-regulation walkthrough.

# %% [markdown]
# ## Weak points: molecules that hold the response together
#
# A *cut molecule* is a molecule whose removal splits the responsive network
# into separate pieces. It is the only connection between the pieces, so if it
# were missing, or its enzyme inhibited, the response on one side could no
# longer reach the other. `structural_vulnerability` lists these molecules:
# `fragments` is the number of pieces left after removing it, and
# `largest_loss` is the share of the response that is cut off.

# %%
vulnerability = structural_vulnerability(graph)
vulnerability

# %%
plot_structural_vulnerability(vulnerability)
plt.show()

# %% [markdown]
# ## Is the convergence between layers more than chance?
#
# Many analyses here rest on reactions where a changed enzyme and a changed
# metabolite meet. Some such meetings are expected by chance alone: if many
# molecules change, and some reactions have many connections, changed molecules
# will sometimes land on the same reaction.
#
# To check, `convergence_significance` keeps the network and the number of
# changed molecules in each layer fixed, reassigns at random which molecules
# count as changed, and counts the meetings again. Repeating this 500 times
# gives the number expected by chance. If the real count is well above that
# range, the layers converge more than chance would produce.

# %%
convergence = convergence_significance(graph, n_randomisations=500)
print(f"{convergence['observed']} reactions where both axes changed; "
      f"{convergence['null_mean']:.1f} expected by chance "
      f"(z = {convergence['z']:+.1f}, p = {convergence['p_value']:.3g})")

# %% [markdown]
# The figure shows the counts from the 500 random reassignments as a histogram,
# with the real count as a vertical line. In this small example the real count
# lies inside the random range, so the convergence could be chance. With 89
# molecules there is little room to detect it; the studies test it on
# organism-wide networks.

# %%
plot_convergence_null(convergence)
plt.show()

# %% [markdown]
# ## Diffusion: what lies close to the response
#
# Random walk with restart spreads a starting score over the network. A walker
# repeatedly steps to a random neighbour, and at each step jumps back to one of
# the starting nodes with probability 1 - `alpha`. The time spent at each node
# is its score. Nodes close to many changed molecules score higher than nodes
# close to one. The factor analyses in the studies use this to find the network
# neighbourhood of each factor. Here the changed molecules are the starting
# nodes, and the table lists the unchanged molecules closest to them.

# %%
seeds = {node: 1.0 for node in regulated_nodes(graph)}
settled = random_walk_with_restart(seeds, graph, alpha=0.85)
nearest = (pd.Series(settled).drop(labels=list(seeds), errors="ignore")
           .sort_values(ascending=False).head(8))
pd.DataFrame({
    "name": [name(n, 44) for n in nearest.index],
    "layer": [graph.nodes[n].get("layer") for n in nearest.index],
    "score": nearest.values,
})

# %% [markdown]
# A diffusion score is only meaningful against a baseline, because a large set
# of well-connected nodes is close to everything.
# `network_enrichment_permutation` asks whether the changed molecules are more
# connected to each other than chance: it counts the edges among them, then
# counts the edges among random sets of the same size.

# %%
enrichment = network_enrichment_permutation(list(seeds), graph, n_permutations=500)
print(f"{enrichment['observed_edges']} edges among the {enrichment['n_features_in_graph']} "
      f"changed molecules; random sets of that size have "
      f"{enrichment['mean_null_edges']:.1f} +/- {enrichment['std_null_edges']:.1f} "
      f"(p = {enrichment['p_value']:.3g})")

# %% [markdown]
# ## Choosing what to draw
#
# An organism-wide network is far too large to draw. Keeping the
# highest-degree nodes does not help, because on a real network those are
# interactome hubs such as ribosomal proteins, and no reactions survive.
# `transomic_backbone` selects the reactions where changed enzymes and changed
# metabolites meet, and keeps the molecules around them. Every network figure in
# the studies is drawn this way.

# %%
backbone = transomic_backbone(graph, max_reactions=8)
print(f"{backbone.number_of_nodes()} nodes selected from {graph.number_of_nodes()}")
pd.Series([d.get("layer") for _, d in backbone.nodes(data=True)]).value_counts().to_frame("nodes")

# %% [markdown]
# BRENDA lists every compound shown to affect an enzyme in a test tube,
# including drugs and laboratory reagents that no cell contains.
# `metabolic_effectors_only=True` (the default) draws only regulators that also
# occur in the organism's metabolism.

# %%
for only_metabolic in (True, False):
    selected = transomic_backbone(graph, max_reactions=8,
                                  metabolic_effectors_only=only_metabolic)
    print(f"metabolic_effectors_only={only_metabolic!s:<5} "
          f"{selected.number_of_nodes():>3} nodes, {selected.number_of_edges():>3} edges")

# %% [markdown]
# To save the network or send it to another program, see the walkthrough
# *Saving, exporting and sharing a network*.
