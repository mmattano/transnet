# %% [markdown]
# # Brown adipocytes: thermogenic lipolysis
#
# **Question.** When brown fat cells are stimulated to produce heat, which
# reactions change, through which regulatory mechanism, and in what order?
#
# **Data.** Immortalised mouse brown adipocytes stimulated with
# norepinephrine, which triggers heat production without shivering. Six
# replicates per time point; transcriptome and proteome at 0, 4 and 24 h;
# metabolome at seven time points between 0 and 24 h. The data and the
# published analysis are
# [Anagho-Mattanovich *et al.*, *iScience* 28(9):113382, 2025](https://doi.org/10.1016/j.isci.2025.113382),
# which compared five ways of integrating them. This notebook runs the whole
# TransNet catalogue on the data.
#
# Needs the mouse network: `python maintenance/build_networks.py --organisms
# mouse --brenda`.

# %%
import matplotlib.pyplot as plt
import pandas as pd

from transnet import (
    build_name_id_map,
    assign_temporal_parameters,
    compare_transomic_networks,
    cross_layer_connectivity,
    downstream_influence,
    temporal_network_structure,
    expression_concordance,
    layer_coverage,
    map_omics_to_network,
    metabolite_regulatory_roles,
    reaction_regulation_table,
    regulation_axis_summary,
    regulatory_role_enrichment,
    responsive_subnetwork,
    trace_regulatory_paths,
    transcription_factor_activity,
    transomic_hubs,
)
from transnet.analysis.factors import (
    matched_transcript_protein,
    pairing_robustness,
    sample_pairing_check,
    factor_cross_layer_agreement,
    factor_design_association,
    factor_network_propagation,
    factor_network_coherence,
    factor_variance_explained,
    fit_factors,
)
from transnet.analysis import (
    cluster_temporal_trajectories,
    compute_timecourse_statistics,
    id_mapping,
    path_verdicts,
    versus_chance,
)
from transnet.datasets import (
    BROWN_ADIPOCYTE_HOURS,
    DATA_DIR,
    brown_adipocyte_contrast,
    load_brown_adipocyte_omics,
)
from transnet.io import read_network

OUT = DATA_DIR / "brown_adipocyte_results"
OUT.mkdir(parents=True, exist_ok=True)

omics = load_brown_adipocyte_omics()
pd.DataFrame({layer: {"features": frame.shape[0], "samples": frame.shape[1]}
              for layer, frame in omics.items()}).T

# %% [markdown]
# ## The network
#
# The organism-wide mouse network: KEGG reactions, UniProt proteins, STRING
# interactions, ChIP-Atlas binding and BRENDA allosteric effectors.

# %%
source = DATA_DIR / "mouse" / "latest"
graph = read_network(str(source / "interactions.csv"), nodes_file=str(source / "nodes.csv"))
print(f"{graph.number_of_nodes():,} molecules, {graph.number_of_edges():,} relationships")

# %% [markdown]
# ## How much of this network crosses layers
#
# If most edges run within a layer, the network is in effect a set of separate
# single-omics networks. The share of edges that cross layers measures how much
# of it is genuinely trans-omic.

# %%
from transnet.visualization import plot_layer_connectivity

connectivity = cross_layer_connectivity(graph)
print(f"{connectivity['cross_layer_fraction']:.0%} of edges cross between layers")
layer_figure = plot_layer_connectivity(connectivity,
                                       title="Mouse network: edges between layers")
plt.show()

# %% [markdown]
# ## Identifiers
#
# The transcripts are identified by Ensembl ids and the network's genes by
# Entrez ids; the metabolites by name and the network's by KEGG compound ids.
# Both translations are built once and cached.

# %%
idmap_file = OUT / "ensembl_to_entrez.json"
if idmap_file.exists():
    gene_map = pd.read_json(idmap_file, typ="series").astype(str).to_dict()
else:
    mapped = id_mapping(
        pd.DataFrame({"feature": omics["Transcriptome"].index}),
        id_col="feature", from_type="ensembl.gene", to_type="entrezgene",
        organism="mouse",
    ).dropna(subset=["entrezgene_id"])
    gene_map = dict(zip(mapped["feature"].astype(str),
                        mapped["entrezgene_id"].astype(str)))
    pd.Series(gene_map).to_json(idmap_file)

metabolite_map = build_name_id_map(graph, list(omics["Metabolome"].index))
id_map = {"Transcriptome": gene_map, "Metabolome": metabolite_map}
print(f"{len(gene_map):,} transcripts and {len(metabolite_map)} metabolites resolved")

# %% [markdown]
# ## The contrast
#
# The main contrast is 24 h of stimulation against unstimulated cells. Every
# layer is on a log2 scale, so a fold change is a difference of means. Its
# standard error is kept too, which later allows a protein's change to be
# compared with its transcript's without a significance cut-off.

# %%
CONTROL, TREATMENT = "0h", "24h"
tables = {layer: brown_adipocyte_contrast(frame, CONTROL, TREATMENT)
          for layer, frame in omics.items()}

report = map_omics_to_network(
    graph, tables, id_column="feature", log2fc_column="log2FC",
    qvalue_column="padj", se_column="se", id_map=id_map, qvalue_threshold=0.05,
)
report.per_layer

# %% [markdown]
# `brown_adipocyte_contrast` is a thin wrapper around
# `compute_differential_expression`, which tests every feature of a
# features x samples table between two sets of sample columns and corrects for
# multiple testing (Benjamini-Hochberg). Called directly on the metabolome, it
# also shows how much the choice of test matters: a Welch t-test (used
# throughout this study) against the rank-based Wilcoxon test.

# %%
from transnet.analysis import compute_differential_expression
from transnet.datasets import timepoint_columns

metabolome_table = omics["Metabolome"].reset_index()
metabolome_table = metabolome_table.rename(columns={metabolome_table.columns[0]: "feature"})
pd.DataFrame({
    method: {"changed at q <= 0.05": int((compute_differential_expression(
        metabolome_table,
        control_cols=timepoint_columns(omics["Metabolome"], CONTROL),
        treatment_cols=timepoint_columns(omics["Metabolome"], TREATMENT),
        id_col="feature", method=method, data_is_log=True,
    )["adj_p_value"] <= 0.05).sum())}
    for method in ("t-test", "wilcoxon")
})

# %%
coverage = layer_coverage(graph)
coverage

# %% [markdown]
# The metabolome covers a small fraction of the network's compounds, and that
# limit runs through everything below: a reaction is judged on the metabolites
# that were measured, not on the ones that exist.

# %% [markdown]
# ## The responsive network
#
# The molecules that changed and the edges between them.

# %%
responsive = responsive_subnetwork(graph)
print(f"{responsive.number_of_nodes():,} responsive molecules, "
      f"{responsive.number_of_edges():,} relationships between them")

# %% [markdown]
# `detect_communities` splits the responsive network into densely connected
# groups. Most of the large ones here are a transcription factor with the
# hundreds of genes it binds, so they contain transcripts and one protein. The
# figure draws the community with the most reactions, the metabolic part of the
# response, with the most connected and the changed molecules named.

# %%
from transnet.analysis import detect_communities
from transnet.visualization import plot_community_network

communities = detect_communities(responsive, method="louvain")
membership = pd.DataFrame({"community": pd.Series(communities)})
membership["layer"] = [responsive.nodes[n].get("layer") for n in membership.index]
reactions_per_community = (membership["layer"] == "Reactions").groupby(
    membership["community"]).sum()
metabolic = int(reactions_per_community.idxmax())
print(f"{membership['community'].nunique()} communities; community {metabolic} holds "
      f"{int(reactions_per_community.max())} reactions")
community_figure = plot_community_network(
    responsive, communities, community=metabolic, label_top=20,
    title="Brown adipocytes: communities in the responsive network")
plt.show()

# %% [markdown]
# ## Which axis regulates each reaction

# %%
regulation = reaction_regulation_table(graph)
regulated = regulation[(regulation["gene_axis"] != 0) | (regulation["metabolite_axis"] != 0)]
enzyme_only = regulated[(regulated["gene_axis"] != 0) & (regulated["metabolite_axis"] == 0)]
metabolite_only = regulated[(regulated["gene_axis"] == 0) & (regulated["metabolite_axis"] != 0)]
both = regulated[(regulated["gene_axis"] != 0) & (regulated["metabolite_axis"] != 0)]
controversial = regulated[regulated["controversial"]]

print(f"{len(regulated):,} reactions regulated: {len(enzyme_only):,} through enzyme amount "
      f"only, {len(metabolite_only):,} through metabolites only, {len(both):,} through both")
print(f"{len(controversial):,} are controversial -- the two axes point opposite ways")
regulation.to_csv(OUT / "reaction_regulation.csv", index=False)
regulated["gene_axis_evidence"].value_counts(dropna=False).to_frame("reactions")

# %% [markdown]
# `gene_axis_evidence` shows what each enzyme-axis call rests on: the protein
# alone, the transcript alone, or both moving together (`gene_protein`).
#
# ### The controversial reactions
#
# In these reactions the enzyme axis and the metabolite axis point in opposite
# directions. Each bar in the figure is one measured molecule's push on the
# reaction: to the right it speeds the reaction up, to the left it slows it
# down.

# %%
from transnet.visualization import plot_axis_composition, plot_controversial_reactions

controversial.head(10)[["reaction", "name", "gene_axis", "metabolite_axis",
                        "allosteric_regulators"]]

# %%
controversial_figure = plot_controversial_reactions(graph, regulation, max_enzymes=6,
                                                    title="Brown adipocytes, 24 h: "
                                                          "controversial reactions")
plt.show()

# %%
composition = pd.DataFrame([{
    "contrast": f"{TREATMENT} vs {CONTROL}",
    "enzyme_axis_only": len(enzyme_only),
    "metabolite_axis_only": len(metabolite_only),
    "both_axes": len(both),
    "controversial": len(controversial),
}])
composition_figure = plot_axis_composition(composition)
plt.show()

# %% [markdown]
# ## Per-pathway balance
#
# The same attribution per pathway: which pathways change through enzyme
# amount, which through their metabolites, and where the two disagree.
# `kegg_reaction_pathways` assigns each reaction to its mouse KEGG pathways.

# %%
from transnet.api import kegg_reaction_pathways
from transnet.visualization import plot_regulation_axes

pathways = kegg_reaction_pathways(regulation["reaction"], organism="mmu")
pathway_balance = regulation_axis_summary(regulation, pathway_map=pathways)
pathway_balance.head(12)[["pathway", "n_reactions", "gene_activated", "gene_inhibited",
                          "metabolite_activated", "metabolite_inhibited", "n_controversial"]]

# %%
axes_figure = plot_regulation_axes(
    pathway_balance[pathway_balance["pathway"] != "unassigned"].head(15),
    title="Brown adipocytes, 24 h: regulation by pathway")
plt.show()

# %% [markdown]
# ## Is each protein change transcriptional?
#
# For every gene measured as both transcript and protein,
# `expression_concordance` asks whether the protein followed its transcript.
# The table below shows a few genes of lipid metabolism and thermogenesis.

# %%
concordance = expression_concordance(graph)
counts = concordance["counts"]
print(f"{counts['concordant']} concordant, {counts['protein_only']} protein-only, "
      f"{counts['transcript_only']} transcript-only, {counts['discordant']} against "
      f"their transcript")
print(f"{counts.get('protein_beyond_transcript')} proteins moved significantly further "
      f"than their transcript (of {counts.get('tested_difference')} testable)")
concordance["table"].to_csv(OUT / "concordance.csv", index=False)

named = concordance["table"].assign(
    symbol=lambda d: [graph.nodes[p].get("symbol") or graph.nodes[p].get("name")
                      for p in d["protein"]])
PAPER_GENES = ["Sqle", "Fdft1", "Pnpla2", "Plin2", "Pdk4", "Ucp1"]
columns = [c for c in ("symbol", "gene_log2fc", "protein_log2fc", "category",
                       "difference", "difference_q", "protein_beyond_transcript")
           if c in named.columns]
named[named["symbol"].isin(PAPER_GENES)][columns]

# %%
from transnet.visualization import plot_expression_concordance

concordance_figure = plot_expression_concordance(concordance)
plt.show()

# %% [markdown]
# ## Which factors drive the responsive genes?

# %%
factors = transcription_factor_activity(graph, min_targets=5)
implicated = factors[factors["q_value"] <= 0.05]
print(f"{len(implicated)} of {len(factors)} transcription factors implicated")
implicated.head(10)[["name", "n_responsive_targets", "n_up", "n_down", "q_value",
                     "inferred_activity", "factor_regulated"]]

# %%
from transnet.visualization import plot_tf_activity

tf_figure = plot_tf_activity(factors, title="Brown adipocytes, 24 h: transcription factors")
plt.show()

# %% [markdown]
# ## Do the changed metabolites regulate anything?

# %%
roles = regulatory_role_enrichment(metabolite_regulatory_roles(graph))
row = roles["enrichment"].set_index("role").loc["any"]
counts = roles["counts"]
print(f"{counts['n_differential_regulators']} of {counts['n_differential']} changed "
      f"metabolites act on an enzyme ({counts['fraction_differential_regulators']:.0%} "
      f"against {counts['fraction_background_regulators']:.0%} of all measured; "
      f"q = {row['q_value']:.2g})")
roles["regulators"].head(10)[["name", "log2fc", "reactions_activated", "reactions_inhibited"]]

# %%
from transnet.visualization import plot_metabolite_regulators

regulator_figure = plot_metabolite_regulators(roles)
plt.show()

# %% [markdown]
# ## Signed paths to the metabolome

# %%
paths = trace_regulatory_paths(graph, target_layer="Metabolome", max_length=4, max_paths=20000)
verdicts, agree, tested = path_verdicts(paths)
print(f"{agree} of {tested} changed metabolites predicted in the right direction: "
      f"{versus_chance(agree, tested)}")
verdicts.head(10)[["target", "observed", "predicted", "agrees", "n_paths"]]

# %% [markdown]
# The best-supported paths, each coloured by the measured change of its
# molecules. When the overall rate is close to chance, the figure shows whether
# the wrong predictions are spread across many metabolites or concentrated on a
# few.

# %%
from transnet.visualization import plot_regulatory_paths

paths_figure = plot_regulatory_paths(paths, graph, top_n=8)
plt.show()

# %% [markdown]
# ## Hubs across layers

# %%
hubs = transomic_hubs(responsive, top_percent=2)
hubs.head(10)[["name", "layer", "cross_layer_degree", "n_layers_touched"]]

# %%
from transnet.visualization import plot_transomic_hubs

hub_figure = plot_transomic_hubs(hubs, title="Brown adipocytes: molecules connecting layers")
plt.show()

# %% [markdown]
# ## Does the upper hierarchy predict the metabolites?
#
# Propagating the changed genes and proteins forward along signed edges gives
# each metabolite a predicted direction. Comparing it with the measurement is
# the most demanding test in the catalogue, because every step of the network
# has to be right. On real data it often fails, and the result is reported
# either way.

# %%
seeds = {n: float(d["log2fc"]) for n, d in graph.nodes(data=True)
         if d.get("layer") in ("Transcriptome", "Proteome")
         and d.get("regulated") and d.get("log2fc") is not None}
influence = downstream_influence(graph, seeds, target_layer="Metabolome")
measured = influence[influence["observed"].fillna(0) != 0] if "observed" in influence else influence
if not measured.empty and "agrees" in measured:
    agree = int(measured["agrees"].astype("boolean").fillna(False).sum())
    print(f"{agree} of {len(measured)} changed metabolites predicted correctly: "
          f"{versus_chance(agree, len(measured))}")
measured.head(10)

# %%
from transnet.visualization import plot_downstream_influence

influence_figure = plot_downstream_influence(influence)
plt.show()

# %% [markdown]
# ## Timing on the network
#
# The seven-point metabolome gives each metabolite a half-response time: the
# time it takes to reach half of its largest change. Morita *et al.* asked
# whether the best-connected molecules respond fastest.

# %%
metabolome_labels = [str(c).split("_")[0] for c in omics["Metabolome"].columns]
metabolome_means = pd.DataFrame({
    timepoint: omics["Metabolome"].loc[:, [c for c, label in
                                           zip(omics["Metabolome"].columns, metabolome_labels)
                                           if label == timepoint]].mean(axis=1)
    for timepoint in dict.fromkeys(metabolome_labels)
})
# the network is keyed by KEGG compound id and the table by name; the columns
# are named for the timepoint, so rename them to the hours they stand for
hours = [BROWN_ADIPOCYTE_HOURS[label] for label in metabolome_means.columns]
timecourse = (metabolome_means.rename(index=metabolite_map)
              .set_axis(hours, axis=1).reset_index())
timecourse = timecourse.rename(columns={timecourse.columns[0]: "feature"})
timecourse = timecourse[timecourse["feature"].isin(graph.nodes)]
parameters = assign_temporal_parameters(
    graph, {"Metabolome": timecourse}, id_column="feature", time_columns=hours,
)
print(f"{parameters['t_half'].notna().sum()} metabolites have a half-response time")
structure = temporal_network_structure(graph)
print(structure["degree_vs_thalf"]["interpretation"])
structure["per_layer_thalf"]

# %%
from transnet.visualization import plot_temporal_structure

timing_figure = plot_temporal_structure(graph, structure, time_unit="h")
plt.show()

# %% [markdown]
# ## 4 h against 24 h
#
# The early (4 h) and late (24 h) responses of the same cells, compared as
# networks. The question is not only which molecules differ, but which *kinds*
# of regulation each response uses.

# %%
early = graph.copy()
early_tables = {layer: brown_adipocyte_contrast(frame, CONTROL, "4h")
                for layer, frame in omics.items()}
map_omics_to_network(early, early_tables, id_column="feature", log2fc_column="log2FC",
                     qvalue_column="padj", se_column="se", id_map=id_map,
                     qvalue_threshold=0.05)
comparison = compare_transomic_networks(
    responsive_subnetwork(early), responsive_subnetwork(graph), "4h", "24h")
print(f"edge Jaccard between the two responses: "
      f"{comparison['summary']['edge_jaccard']:.2f}")
comparison["edges_by_type"]

# %%
from transnet.visualization import plot_condition_comparison

comparison_figure = plot_condition_comparison(comparison, "4h", "24h", graph=graph)
plt.show()

# %% [markdown]
# ## Factors, read through the network
#
# The published analysis fitted MOFA factors to these data. A factor model
# finds molecules that vary together across samples. It cannot tell whether a
# factor's strongest molecules are also connected biochemically; the network
# can.

# %%
shared_samples = sorted(
    set(omics["Transcriptome"].columns)
    & set(omics["Proteome"].columns)
    & set(omics["Metabolome"].columns)
)
matrices = {layer: frame[shared_samples].T for layer, frame in omics.items()}

# A few metabolites are missing in some samples. network_guided_imputation fills
# each gap from the same sample's values for the metabolite's network
# neighbours, which is a better guess than the average over all samples.
from transnet.analysis.factors import network_guided_imputation

missing = int(matrices["Metabolome"].isna().sum().sum())
matrices["Metabolome"] = network_guided_imputation(matrices["Metabolome"], graph,
                                                   id_map=metabolite_map)
print(f"{missing} missing metabolite values filled from network neighbours")
design = pd.DataFrame(
    {"timepoint": [str(s).split("_")[0] for s in shared_samples]}, index=shared_samples)
print(f"{len(shared_samples)} samples measured in all three layers: "
      + ", ".join(sorted(set(design['timepoint']))))

# %% [markdown]
# **First, are the layers measured on the same cultures?** A joint factor
# model, like MOFA and DIABLO, treats `0h_01` as the same sample in every
# layer. The files share column *names*, but whether they share *samples* can
# be tested: within a time point, a culture with a high transcript level of a
# gene should, if it is the same culture, also tend to have a high level of the
# protein.

# %%
transcripts, proteins = matched_transcript_protein(
    omics["Transcriptome"], omics["Proteome"], graph, id_maps=id_map)
groups = design["timepoint"]
pairing = sample_pairing_check(transcripts, proteins, groups, n_permutations=500)
print(f"{transcripts.shape[1]:,} genes measured as transcript and protein")
print(f"agreement as listed {pairing['observed']:.3f}, under shuffled pairings "
      f"{pairing['null_mean']:.3f} (p = {pairing['p_value']:.2f}) -> "
      f"{'paired' if pairing['paired'] else 'pairing not confirmed'}")

# %% [markdown]
# The pairing is confirmed: the layers come from the same cultures, so the
# joint model is sound. A factor carried by differences between cultures within
# a time point is then as real as one carried by the time points themselves.
# The robustness check below still reports which kind each factor is.

# %%
factorisation = fit_factors(matrices, n_components=5)
association = factor_design_association(factorisation.factors_, design, ["timepoint"])
variance = factor_variance_explained(factorisation)
coherence = factor_network_coherence(graph, factorisation.loadings_, id_maps=id_map,
                                     top_n=100, n_permutations=500)

(association.pivot(index="factor", columns="term", values="partial_eta_squared")
 .join(variance[variance["layer"] == "all layers"].set_index("factor")["variance_accounted"])
 .join(coherence["table"].set_index("factor")[["fold_enrichment", "q_value"]])
 .round(3))

# %% [markdown]
# **Every factor across the design**, and **every factor on the network**: the
# top features the network actually joins, drawn with the reactions between
# them.

# %%
from transnet.visualization import plot_factor_overview, plot_factor_scores

FACTOR_OUT = OUT / "factors"
FACTOR_OUT.mkdir(parents=True, exist_ok=True)

plot_factor_scores(factorisation.factors_, design, x="timepoint",
                   x_order=["0h", "4h", "24h"], association=association,
                   title="Brown adipocytes: factor scores across the design"
                   ).savefig(FACTOR_OUT / "factor_scores.png", dpi=150, bbox_inches="tight")
plt.show()

plot_factor_overview(association, variance, coherence["table"],
                     title="Brown adipocytes: what each factor follows, explains and connects"
                     ).savefig(FACTOR_OUT / "factor_overview.png", dpi=150, bbox_inches="tight")
plt.show()

# %%
from transnet.visualization import plot_factor_network

for factor in factorisation.factors_.columns:
    linked = coherence["nodes"].get(factor, [])
    if not linked:
        print(f"{factor}: none of its top features are joined on the network")
        continue
    figure = plot_factor_network(graph, factorisation.loadings_, factor, id_maps=id_map,
                                 nodes=linked, max_nodes=60)
    figure.savefig(FACTOR_OUT / f"network_{factor}.png", dpi=150, bbox_inches="tight")
    plt.show()

# %% [markdown]
# Three further readings, each asking something of a factor that the factor
# model cannot ask of itself.
#
# **Do the layers agree?** A joint factor is fitted on all layers at once, but
# one layer alone can carry it. `factor_layer_scores` projects the samples onto
# each factor using one layer at a time, and `factor_cross_layer_agreement`
# correlates those projections between layers.

# %%
from transnet.analysis.factors import factor_layer_scores

layer_scores = factor_layer_scores(factorisation)
print({layer: scores.shape for layer, scores in layer_scores.items()})
agreement = factor_cross_layer_agreement(factorisation)
agreement.pivot_table(index="factor", columns=["layer_a", "layer_b"],
                      values="correlation").round(2)

# %% [markdown]
# **Do a factor's layers land in the same place on the network?** A joint
# factor includes transcripts, proteins and metabolites, but that makes it
# trans-omic only if those molecules are related biochemically. The strongest
# molecules of each layer are spread over the network separately (random walk
# with restart), and the overlap of the resulting profiles is compared with
# that of random molecules *from the same layers*. Drawing the random molecules
# from all layers together would give mostly transcripts, because the
# transcriptome is a hundred times larger than the metabolome.

# %%
propagation = factor_network_propagation(graph, factorisation.loadings_,
                                         id_maps=id_map, top_n=50, n_permutations=100)
propagation["table"].round(3)

# %% [markdown]
# **Which factors to trust.** All the readings above in one table. A factor
# worth interpreting follows the design, is carried by the layers *together*
# (whatever the pairing of replicates), and has its molecules converge on the
# network. Each column rules out a different way a factor can mislead.

# %%
robustness = pairing_robustness(factorisation, groups, n_shuffles=200)

synthesis = (
    association.pivot(index="factor", columns="term", values="partial_eta_squared")
    .rename(columns={"timepoint": "timepoint eta2"})
    .join(robustness.set_index("factor")[["agreement", "retained", "verdict"]])
    .join(coherence["table"].set_index("factor")[["fold_enrichment"]]
          .rename(columns={"fold_enrichment": "direct links vs chance"}))
    .join(propagation["table"].set_index("factor")[["overlap", "null_overlap", "q_value"]]
          .rename(columns={"q_value": "overlap q"}))
)
# a factor is read further when the layers carry it together, its joint
# structure does not rest on an unconfirmed pairing, and its layers converge on
# the network
joint = synthesis["verdict"] != "single-layer"
trusted_pairing = (synthesis["verdict"] == "between-group") | pairing["paired"]
synthesis["interpret?"] = (joint & trusted_pairing & (synthesis["overlap q"] <= 0.05)
                           ).map({True: "yes", False: "no"})
synthesis.round(3)

# %% [markdown]
# For the factor that passes every reading, the molecules that receive the
# most signal, *relative to random starting molecules from the same layers*,
# show what it is about. Comparing with random starts stops the network's
# best-connected molecules from ranking first for every factor.

# %%
trusted = synthesis[synthesis["interpret?"] == "yes"]
leading = (trusted if not trusted.empty else synthesis).sort_values(
    "overlap", ascending=False).index[0]
print(f"{leading}: where its transcripts, proteins and metabolites meet")
propagation["top_nodes"][leading].head(12)[["name", "layer", "enrichment"]]

# %% [markdown]
# **What is the factor about?** Its strongest molecules, drawn where they sit
# on the network rather than as a bar chart of loadings.

# %%
plot_factor_network(graph, factorisation.loadings_, leading, id_maps=id_map, top_n=12)
plt.show()

# %% [markdown]
# ## Timing: three metabolic states
#
# Do the metabolites follow a few distinct time courses? Clustering the
# metabolite time courses over all seven time points groups metabolites that
# change together.

# %%
metabolome = omics["Metabolome"]
labels = [str(c).split("_")[0] for c in metabolome.columns]
# one x-value per group, in the order the groups appear
hours = [BROWN_ADIPOCYTE_HOURS[label] for label in dict.fromkeys(labels)]

omnibus = compute_timecourse_statistics(metabolome, groups=labels, method="anova",
                                        trend_test=True, continuous_x=hours)
print(f"{(omnibus['fdr'] <= 0.05).sum()} of {len(omnibus)} metabolites change across "
      f"the time course (one-way ANOVA over all seven points, BH-corrected); "
      f"{(omnibus['trend_pval'] <= 0.05).sum()} follow a monotone trend")
omnibus.sort_values("fdr").head(8)

# %%
clusters, profiles = cluster_temporal_trajectories(metabolome, timepoints=labels,
                                                   n_clusters=3)
clusters.value_counts().rename("metabolites").to_frame()

# %%
ordered = [t for t in BROWN_ADIPOCYTE_HOURS if t in profiles.columns]
profiles[ordered].T.plot(marker="o", figsize=(7, 3.4))
plt.xlabel("timepoint")
plt.ylabel("mean log2 intensity (scaled)")
plt.title("Metabolite trajectories, clustered")
plt.show()

# %% [markdown]
# ## The wiring itself: motifs, bottlenecks, and a null
#
# The analyses above read the data through the network. These three look at
# the structure of the responsive network itself, which a list of molecules
# cannot provide at all.

# %%
from transnet import (
    convergence_significance,
    regulatory_motifs,
    structural_vulnerability,
)

motifs = regulatory_motifs(graph, responsive_only=True)
print(motifs["counts"])
motifs["motifs"].head(12)[["motif", "reaction_name", "metabolite_name", "enzyme",
                           "sign_product"]]

# %% [markdown]
# Each product-inhibition row is a reaction slowed down by its own product,
# with both the enzyme and the product measured here. This is one way a
# reaction can slow down while its enzyme increases.

# %%
vulnerability = structural_vulnerability(responsive)
vulnerability.head(10)

# %%
from transnet.visualization import (
    plot_convergence_null,
    plot_regulatory_motifs,
    plot_structural_vulnerability,
)

vulnerability_figure = plot_structural_vulnerability(vulnerability)
plt.show()

# %% [markdown]
# These are the molecules the response runs through: removing one would
# disconnect part of the response from the rest.

# %%
convergence = convergence_significance(graph, n_randomisations=200)
print(f"{convergence['observed']:,} reactions have both a changed enzyme and a changed "
      f"metabolite; shuffling which molecules changed gives "
      f"{convergence['null_mean']:.0f} (z = {convergence['z']:+.1f}, "
      f"p = {convergence['p_value']:.3g})")

# %%
convergence_figure = plot_convergence_null(convergence)
plt.show()

# %%
motif_figure = plot_regulatory_motifs(motifs)
plt.show()

# %% [markdown]
# ## Glycolysis enzymes
#
# The transcript and protein changes of the glycolysis enzymes measured in both
# layers.

# %%
GLYCOLYSIS = ["Pfkl", "Pfkp", "Pfkm", "Pklr", "Pkm", "Hk1", "Hk2", "Gpi1", "Eno1"]
named[named["symbol"].isin(GLYCOLYSIS)][
    [c for c in ("symbol", "gene_log2fc", "protein_log2fc", "category") if c in named.columns]
]

# %% [markdown]
# ## Figures and exports
#
# Every figure is saved, and the responsive network is exported in every format
# TransNet writes (see the *saving, exporting and sharing* walkthrough).

# %%
from transnet.io import to_arena3d, to_cytoscape_json, to_transomics2cytoscape, write_network
from transnet.visualization import (
    plot_transomic_network,
    plot_transomic_network_interactive,
    transomic_backbone,
)

for figure, name in [
    (plot_transomic_network(responsive, title=f"Brown adipocytes, {TREATMENT} vs {CONTROL}"),
     "network"),
    (layer_figure, "layer_connectivity"),
    (community_figure, "communities"),
    (composition_figure, "axis_composition"),
    (axes_figure, "regulation_axes"),
    (controversial_figure, "controversial_reactions"),
    (concordance_figure, "concordance"),
    (tf_figure, "tf_activity"),
    (regulator_figure, "metabolite_regulators"),
    (hub_figure, "transomic_hubs"),
    (influence_figure, "downstream_influence"),
    (timing_figure, "temporal_structure"),
    (comparison_figure, "early_vs_late"),
    (convergence_figure, "convergence_null"),
    (vulnerability_figure, "structural_vulnerability"),
    (motif_figure, "regulatory_motifs"),
    (paths_figure, "regulatory_paths"),
]:
    figure.savefig(OUT / f"{name}.png", dpi=150, bbox_inches="tight")
plt.show()

exports = OUT / "exports"
exports.mkdir(exist_ok=True)
write_network(responsive, str(exports / "csv"))
to_cytoscape_json(responsive, exports / "network_cytoscape.json")
to_arena3d(responsive, exports / "network_arena3d.json")
to_transomics2cytoscape(responsive, zip_path=exports / "network_transomics2cytoscape.zip")
# a browser draws a few hundred nodes well; export the backbone the figures use
html_network = transomic_backbone(graph, max_reactions=40)
plot_transomic_network_interactive(html_network, layout="layered",
                                   title="Brown adipocytes, 24 h").write_html(
    exports / "network.html")
print(f"figures and tables in {OUT}")
print("exports:", ", ".join(sorted(p.name for p in exports.iterdir())))
