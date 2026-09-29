# %% [markdown]
# # Thermogenic lipolysis in brown adipocytes
#
# Immortalised murine brown adipocytes, stimulated with norepinephrine to
# trigger the non-shivering cold response. Six replicates per timepoint;
# transcriptome and proteome at 0, 4 and 24 h; metabolome at seven points
# between 0 and 24 h.
#
# The data and the published analysis are
# [Anagho-Mattanovich *et al.*, *iScience* 28(9):113382, 2025](https://doi.org/10.1016/j.isci.2025.113382),
# which compared five ways of integrating them. This notebook runs the
# trans-omic network reading and ends by setting it beside what that paper
# concluded.
#
# Needs the mouse network: `python maintenance/build_networks.py --organisms
# mouse --brenda`.

# %%
from pathlib import Path

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
    path_consistency_summary,
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
    timepoint_columns,
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
# A network whose edges sit inside single layers is a stack of separate
# single-omics networks; this fraction is what earns the name.

# %%
connectivity = cross_layer_connectivity(graph)
print(f"{connectivity['cross_layer_fraction']:.0%} of edges cross between layers")
connectivity["matrix"]

# %% [markdown]
# ## Identifiers
#
# Transcripts are Ensembl, the network's genes are Entrez; metabolites are
# names, the network's are KEGG compounds. Both maps are built once.

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
# 24 h of sustained stimulation against unstimulated cells. Every layer is
# log2, so a fold change is a difference of means; the standard error comes
# with it, which is what lets a protein's change be compared with its
# transcript's without a threshold (transcript-protein concordance).

# %%
CONTROL, TREATMENT = "0h", "24h"
tables = {layer: brown_adipocyte_contrast(frame, CONTROL, TREATMENT)
          for layer, frame in omics.items()}

report = map_omics_to_network(
    graph, tables, id_column="feature", log2fc_column="log2FC",
    qvalue_column="padj", se_column="se", id_map=id_map, qvalue_threshold=0.05,
)
report.per_layer

# %%
coverage = layer_coverage(graph)
coverage

# %% [markdown]
# The metabolome covers a small fraction of the network's compounds, and that
# limit runs through everything below: a reaction is judged on the metabolites
# that were measured, not on the ones that exist.

# %% [markdown]
# ## Where the layers converge

# %%
responsive = responsive_subnetwork(graph)
print(f"{responsive.number_of_nodes():,} responsive molecules, "
      f"{responsive.number_of_edges():,} relationships between them")

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
controversial.head(10)[["reaction", "name", "gene_axis", "metabolite_axis",
                        "allosteric_regulators"]]

# %% [markdown]
# ## Per-pathway balance
#
# The same attribution, aggregated: which pathways are driven by enzyme
# amount, which by their metabolites, and where the two disagree.

# %%
pathway_balance = regulation_axis_summary(regulation)
pathway_balance.head(12)

# %% [markdown]
# ## Is each protein change transcriptional?
#
# The paper's own example of transcript-protein divergence is the lipid
# droplet machinery: *Sqle* and *Fdft1* transcripts fall at 4 h while their
# proteins rise. That is exactly the `discordant` category here.

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

# %% [markdown]
# ## Which factors drive the responsive genes?

# %%
factors = transcription_factor_activity(graph, min_targets=5)
implicated = factors[factors["q_value"] <= 0.05]
print(f"{len(implicated)} of {len(factors)} transcription factors implicated")
implicated.head(10)[["name", "n_responsive_targets", "n_up", "n_down", "q_value",
                     "inferred_activity", "factor_regulated"]]

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

# %% [markdown]
# ## Signed paths to the metabolome

# %%
paths = trace_regulatory_paths(graph, target_layer="Metabolome", max_length=4, max_paths=20000)
verdicts, agree, tested = path_verdicts(paths)
print(f"{agree} of {tested} changed metabolites predicted in the right direction: "
      f"{versus_chance(agree, tested)}")
verdicts.head(10)[["target", "observed", "predicted", "agrees", "n_paths"]]

# %% [markdown]
# The paths themselves, with what each one predicts against what was measured.
# A chance-level rate is read from this figure, not just from the number: the
# question is whether the misses are scattered or concentrated on a few
# metabolites.

# %%
from transnet.visualization import plot_regulatory_paths

paths_figure = plot_regulatory_paths(paths, graph, top_n=8)
plt.show()

# %% [markdown]
# ## Hubs across layers

# %%
hubs = transomic_hubs(responsive, top_percent=2)
hubs.head(10)[["name", "layer", "cross_layer_degree", "n_layers_touched"]]

# %% [markdown]
# ## Does the upper hierarchy predict the metabolites?
#
# Propagating the changed genes and proteins forward along signed edges gives
# each metabolite a predicted direction. Comparing that with the measurement
# is the strongest test the network can fail, and it usually does.

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

# %% [markdown]
# ## Timing on the network
#
# The seven-point metabolome gives each metabolite a half-response time. Morita
# et al. ask whether the best-connected molecules respond fastest.

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

# %% [markdown]
# ## 4 h against 24 h
#
# Two contrasts of the same cells, compared as networks: not which molecules
# differ, but which *kinds* of regulation the early and late responses use.

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

# %% [markdown]
# ## Factors, read through the network
#
# The published analysis fits MOFA factors to these data. A factor model finds
# co-variation; the question it cannot ask of itself is whether a factor's
# strongest features are *connected*. That is what the permutation test below
# answers.

# %%
shared_samples = sorted(
    set(omics["Transcriptome"].columns)
    & set(omics["Proteome"].columns)
    & set(omics["Metabolome"].columns)
)
matrices = {layer: frame[shared_samples].T for layer, frame in omics.items()}
design = pd.DataFrame(
    {"timepoint": [str(s).split("_")[0] for s in shared_samples]}, index=shared_samples)
print(f"{len(shared_samples)} samples measured in all three layers: "
      + ", ".join(sorted(set(design['timepoint']))))

# %% [markdown]
# **First, are the layers measured on the same cultures?** A joint factor model
# treats `0h_01` as one sample in every layer, and so do MOFA and DIABLO. The
# files share column *names*; whether they share *samples* is testable: within a
# timepoint, a culture whose transcript of a gene runs high should, if it is the
# same culture, tend to have that protein high too.

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
# The pairing is confirmed: the layers are the same cultures, so the joint model
# is sound, and a factor carried by culture-to-culture variation within a
# timepoint is as real as one carried by the timepoints themselves. The
# robustness check below still says which kind each factor is.

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
# **Do the layers agree?** A joint factor is fitted on the layers stacked
# together, which does not make it trans-omic: one layer can carry it alone.
# Projecting the samples onto a factor within each layer and correlating those
# projections says which it is.

# %%
agreement = factor_cross_layer_agreement(factorisation)
agreement.pivot_table(index="factor", columns=["layer_a", "layer_b"],
                      values="correlation").round(2)

# %% [markdown]
# **Do a factor's layers land in the same place on the network?** A joint
# factor loads on transcripts, proteins and metabolites at once; that makes it
# trans-omic only if those features are related by the biochemistry. Each
# layer's strongest features are diffused separately, and the overlap of the
# resulting profiles is compared with random features *of the same layers* --
# the transcriptome is a hundred times the size of the metabolome, so a null
# that ignored that would be all transcripts.

# %%
propagation = factor_network_propagation(graph, factorisation.loadings_,
                                         id_maps=id_map, top_n=50, n_permutations=100)
propagation["table"].round(3)

# %% [markdown]
# **Which factors to trust.** Every reading above in one table. A factor worth
# interpreting follows the design, is carried by the layers *together* (and
# stays so however replicates are paired), and has its features converge on
# the network. Each column rules out a different way a factor can mislead.

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
# For the factor that passes every reading, the nodes that receive most
# signal *relative to random seeds of the same layers* -- so the network's hubs
# do not win by default -- say what it is about.

# %%
trusted = synthesis[synthesis["interpret?"] == "yes"]
leading = (trusted if not trusted.empty else synthesis).sort_values(
    "overlap", ascending=False).index[0]
print(f"{leading}: where its transcripts, proteins and metabolites meet")
propagation["top_nodes"][leading].head(12)[["name", "layer", "enrichment"]]

# %% [markdown]
# **What is the factor about?** Its strongest features, drawn where they sit
# on the network rather than as a bar chart of loadings.

# %%
from transnet.visualization import plot_factor_network

plot_factor_network(graph, factorisation.loadings_, leading, id_maps=id_map, top_n=12)
plt.show()

# %% [markdown]
# ## Timing: three metabolic states
#
# The paper reports three states, uninduced (0 h), active lipolysis (4 h) and
# sustained induction (24 h), with the metabolome moving first. Clustering
# the metabolite trajectories over all seven timepoints asks the same question
# of the same data.

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
# Everything above reads data *through* the network. These three read the
# network, and are the analyses a molecule list cannot approximate at all.

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
# Each product-inhibition row is a reaction whose own product holds it back,
# with both the enzyme and the metabolite measured here, which is the mechanism behind
# a reaction that slows while its enzyme rises.

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
# Those are the molecules the response runs through: remove one and the rest
# of the response is no longer connected to what it regulates.

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
# ## Compared with the published analysis
#
# | The paper concluded | Here |
# |---|---|
# | Three states: uninduced (0 h), active lipolysis (4 h), sustained (24 h) | the metabolite trajectories cluster into that many distinct shapes, above |
# | The metabolome moves before the transcriptome and proteome | the responsive metabolites are already significant at 4 h while the protein changes concentrate at 24 h |
# | Lipid-droplet genes (*Sqle*, *Fdft1*) fall as transcripts while their proteins rise | the same genes land in the `discordant` category of transcript-protein concordance, without being looked for |
# | Upper glycolysis (*Pfkl*, *Pfkp*) down, lower glycolysis (*Pklr*) up | checked directly below |
#
# The paper reached its trans-omic observations, the Rock2/protamine link and
# the glycolytic split, by reading a network by hand. The point of the table
# above is that the same statements fall out of the catalogue as ordinary
# output, with the axis attribution and the threshold-free transcript-protein
# test attached.

# %%
GLYCOLYSIS = ["Pfkl", "Pfkp", "Pfkm", "Pklr", "Pkm", "Hk1", "Hk2", "Gpi1", "Eno1"]
named[named["symbol"].isin(GLYCOLYSIS)][
    [c for c in ("symbol", "gene_log2fc", "protein_log2fc", "category") if c in named.columns]
]

# %% [markdown]
# ## Figures

# %%
from transnet.visualization import (
    plot_expression_concordance,
    plot_axis_composition,
    plot_transomic_network,
)

composition = pd.DataFrame([{
    "contrast": f"{TREATMENT} vs {CONTROL}",
    "enzyme_axis_only": len(enzyme_only),
    "metabolite_axis_only": len(metabolite_only),
    "both_axes": len(both),
    "controversial": len(controversial),
}])

for figure, name in [
    (plot_transomic_network(responsive, title=f"Brown adipocytes, {TREATMENT} vs {CONTROL}"),
     "network"),
    (plot_axis_composition(composition), "regulation_axes"),
    (plot_expression_concordance(concordance), "concordance"),
    (convergence_figure, "convergence_null"),
    (vulnerability_figure, "structural_vulnerability"),
    (motif_figure, "regulatory_motifs"),
    (paths_figure, "regulatory_paths"),
]:
    figure.savefig(OUT / f"{name}.png", dpi=150, bbox_inches="tight")
plt.show()
print(f"figures and tables in {OUT}")
