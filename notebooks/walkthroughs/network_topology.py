# %% [markdown]
# # The network as a network
#
# Everything so far read the network as a *regulatory hierarchy*. This
# notebook reads it as an ordinary graph: components, centrality, communities
# and modules.

# %%
import matplotlib.pyplot as plt
import pandas as pd

from transnet import (
    load_example_network,
    load_example_omics,
    map_omics_to_network,
    regulated_nodes,
    responsive_subnetwork,
    transomic_hubs,
)
from transnet.analysis import (
    compute_centrality_measures,
    compute_network_statistics,
    detect_communities,
    enrichment_analysis,
    find_active_modules,
    identify_hubs,
)
from transnet.visualization import (
    plot_community_network,
    plot_network_metrics,
    plot_transomic_network_interactive,
)

graph = load_example_network()
map_omics_to_network(
    graph, load_example_omics(), id_column="id", log2fc_column="log2FC", qvalue_column="padj",
)

statistics = compute_network_statistics(graph)
{k: statistics[k] for k in ("nodes", "edges", "connected_components",
                            "largest_component_ratio", "avg_degree", "avg_clustering")}

# %% [markdown]
# Path length and diameter are quadratic in the node count, so they are
# computed only when the largest component is small enough and reported as
# `None` otherwise, since an interactome would take hours.

# %%
plot_network_metrics(graph)
plt.show()

# %% [markdown]
# ## Centrality is not trans-omic hubness
#
# The two rankings answer different questions. Degree centrality asks which
# molecule is best connected; `transomic_hubs` asks which molecule joins the
# *layers* the most.

# %%
centrality = compute_centrality_measures(graph, top_n=5)
hubs = transomic_hubs(graph, top_percent=10)

pd.DataFrame({
    "degree centrality": [graph.nodes[n].get("name", n) for n, _ in centrality["degree_centrality"]],
    "betweenness": [graph.nodes[n].get("name", n) for n, _ in centrality["betweenness_centrality"]],
    "plain degree hubs": [graph.nodes[n].get("name", n) for n, _ in identify_hubs(graph, top_n=5)],
    "cross-layer hubs": [graph.nodes[n].get("name", n) for n in hubs["node"].head(5)],
})

# %% [markdown]
# ## Do communities span layers?
#
# A community inside one layer says something about how well that layer is
# annotated. A community holding genes, proteins *and* metabolites is a
# candidate module.

# %%
communities = detect_communities(graph, method="louvain")
membership = pd.Series(communities, name="community")
layers = pd.Series({n: graph.nodes[n].get("layer") for n in communities}, name="layer")
composition = pd.crosstab(membership, layers)
print(f"{len(composition)} communities, "
      f"{int((composition > 0).sum(axis=1).gt(1).sum())} spanning more than one layer")
composition.head(8)

# %%
plot_community_network(graph, communities)
plt.show()

# %% [markdown]
# ## Active modules
#
# `responsive_subnetwork` takes everything past a significance threshold.
# Active modules instead grow the highest-scoring *connected* subnetwork, so a
# molecule just under the threshold is kept when its neighbours carry the
# module (Ideker et al. 2002).

# %%
modules = find_active_modules(graph, score_attr="log2fc", n_modules=3)
for i, module in enumerate(modules, start=1):
    names = [str(graph.nodes[n].get("name", n))[:22] for n in module]
    print(f"module {i}: {len(module):>3} molecules  {', '.join(names[:6])}")

# %%
background = [n for n, d in graph.nodes(data=True) if d.get("layer") == "Reactions"]
largest = [n for n in max(modules, key=len) if n in background]
pathways = {
    "glycolysis": ["R00299", "R00756", "R01786", "R00200"],
    "TCA cycle": ["R00351", "R01325", "R00709"],
}
enrichment_analysis(largest, background, pathways)

# %% [markdown]
# ## Motifs
#
# A pathway diagram tells you that hexokinase is inhibited by its own product.
# A typed, signed network lets you *find* this pattern, and others like it, without
# being told where to look.

# %%
from transnet import (
    convergence_significance,
    regulatory_motifs,
    structural_vulnerability,
)

motifs = regulatory_motifs(graph)
print(motifs["counts"])

from transnet.visualization import plot_regulatory_motifs

plot_regulatory_motifs(motifs)
plt.show()
motifs["motifs"][["motif", "reaction_name", "metabolite_name", "enzyme", "sign_product"]]

# %% [markdown]
# These are the reactions that can slow down while their enzyme rises, which is the
# mechanism behind a "controversial" call in the regulation axes, found here
# from the wiring alone.

# %% [markdown]
# ## What holds the response together?
#
# A cut molecule is one whose removal splits the response into disconnected
# pieces: the single route by which one layer's response reaches another.
# Morita et al. call this the structural half of a network's robustness.

# %%
vulnerability = structural_vulnerability(graph)
vulnerability

# %%
from transnet.visualization import plot_structural_vulnerability

plot_structural_vulnerability(vulnerability)
plt.show()

# %% [markdown]
# ## Do the layers converge more than chance?
#
# "A changed enzyme and a changed metabolite meet at this reaction" is the claim
# the catalogue rests on. Some convergence happens for free, given how many
# molecules changed and how the degrees are distributed. The null keeps the
# network and the number of changes per layer, and shuffles which molecules
# changed.

# %%
convergence = convergence_significance(graph, n_randomisations=500)
print(f"{convergence['observed']} reactions where both axes changed; "
      f"null expects {convergence['null_mean']:.1f} "
      f"(z = {convergence['z']:+.1f}, p = {convergence['p_value']:.3g})")

# %%
from transnet.visualization import plot_convergence_null

plot_convergence_null(convergence)
plt.show()

# %% [markdown]
# ## Diffusion
#
# Random walk with restart spreads a seed score over the network until it
# settles, so a node close to many seeds scores higher than one near a single
# seed. It is the assumption behind the factor readings in the studies, and 
# it also answers a simpler question: which molecules sit near the response.

# %%
from transnet.analysis import network_enrichment_permutation, random_walk_with_restart

seeds = {node: 1.0 for node in regulated_nodes(graph)}
settled = random_walk_with_restart(seeds, graph, alpha=0.85)
nearest = (pd.Series(settled).drop(labels=list(seeds), errors="ignore")
           .sort_values(ascending=False).head(8))
pd.DataFrame({
    "name": [graph.nodes[n].get("name", n)[:44] for n in nearest.index],
    "layer": [graph.nodes[n].get("layer") for n in nearest.index],
    "score": nearest.values,
})

# %% [markdown]
# A diffusion score means nothing on its own, because a large set of
# well-connected nodes is close to everything. `network_enrichment_permutation`
# supplies the baseline: it counts the edges among a feature set, then draws
# random sets of the same size and counts again.

# %%
enrichment = network_enrichment_permutation(list(seeds), graph, n_permutations=500)
print(f"{enrichment['observed_edges']} edges among the {enrichment['n_features_in_graph']} "
      f"responsive molecules; random sets of that size give "
      f"{enrichment['mean_null_edges']:.1f} +/- {enrichment['std_null_edges']:.1f} "
      f"(p = {enrichment['p_value']:.3g})")

# %% [markdown]
# ## Choosing what to draw
#
# A trans-omic network trimmed by degree keeps the interactome hubs, ribosomal
# and splicing proteins, and drops every reaction. `transomic_backbone` selects
# by where regulation converges instead, and it is the step behind every study
# figure in the documentation.

# %%
from transnet.visualization import transomic_backbone

backbone = transomic_backbone(graph, max_reactions=8)
print(f"{backbone.number_of_nodes()} nodes selected from {graph.number_of_nodes()}")
pd.Series([d.get("layer") for _, d in backbone.nodes(data=True)]).value_counts().to_frame("nodes")

# %% [markdown]
# `metabolic_effectors_only` decides whether laboratory reagents count as
# regulators. BRENDA records what has been shown to affect an enzyme in vitro,
# including compounds no cell contains.

# %%
for only_metabolic in (True, False):
    selected = transomic_backbone(graph, max_reactions=8,
                                  metabolic_effectors_only=only_metabolic)
    print(f"metabolic_effectors_only={only_metabolic!s:<5} "
          f"{selected.number_of_nodes():>3} nodes, {selected.number_of_edges():>3} edges")

# %% [markdown]
# ## Out of Python
#
# Four ways to hand the network to someone else. The plotly file is
# self-contained: it opens in any browser with no Python at all, which is what
# to send a collaborator.

# %%
from pathlib import Path

from transnet.io import to_arena3d, to_cytoscape_json, to_transomics2cytoscape

out = Path("exports")
out.mkdir(exist_ok=True)

responsive = responsive_subnetwork(graph)
figure = plot_transomic_network_interactive(responsive, layout="layered",
                                            title="Responsive trans-omic network")
figure.write_html(out / "network.html")

to_cytoscape_json(responsive, str(out / "network_cytoscape.json"))
to_arena3d(responsive, str(out / "network_arena3d"))
to_transomics2cytoscape(responsive, zip_path=out / "network_transomics2cytoscape.zip")

print("wrote:", ", ".join(sorted(p.name for p in out.iterdir())))
print("network.html is self-contained: open it in a browser, no Python needed.")

# %% [markdown]
# The same figure, inline. Hover a node for what it is and what it did; drag
# to zoom; click a legend entry to hide a layer or a relationship class.

# %%
import plotly.io as pio

pio.renderers.default = "notebook_connected"
figure
