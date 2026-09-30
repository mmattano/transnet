# %% [markdown]
# # Obese liver after a glucose load: a genome-wide time course
#
# **Question.** How does the liver's response to glucose differ between lean
# and obese mice, genome-wide and over time, and which regulatory mechanisms
# carry each response?
#
# **Data.** Wild-type and leptin-deficient obese (ob/ob) mouse liver, sampled
# at 0, 20, 60, 120 and 240 minutes after an oral glucose load. The
# transcriptome covers 14,292 genes and the metabolome 162 compounds, each with
# the authors' own fold changes, q-values and half-response times. The
# measurements were published by Kokaji *et al.*, *Science Signaling*
# 13(660):eaaz1236, 2020. They are downloaded when the notebook runs and are
# never redistributed.
#
# The panel study (`obese_liver.py`) could check three claims of this paper only
# partly, on 19 genes at one time point. Here they are tested on the data they
# were made from. Two analyses appear only in this study: the
# transcription-factor inference is compared with the authors' own inference on
# the same data, and measured enzyme affinities show whether a metabolite's
# change can matter to its enzyme at all.

# %%
import json
import warnings

import matplotlib.pyplot as plt
import pandas as pd

from transnet import (
    compare_transomic_networks,
    convergence_significance,
    cross_layer_connectivity,
    downstream_influence,
    layer_coverage,
    map_omics_to_network,
    metabolite_regulatory_roles,
    reaction_regulation_table,
    regulation_axis_summary,
    regulatory_motifs,
    regulatory_role_enrichment,
    responsive_subnetwork,
    structural_vulnerability,
    trace_regulatory_paths,
    transcription_factor_activity,
)
from transnet.analysis import hub_rankings, path_verdicts, versus_chance
from transnet.datasets import (
    DATA_DIR,
    KOKAJI_TIMEPOINTS,
    fetch_kokaji_tables,
    kokaji_contrast,
    load_kokaji_table,
)
from transnet.io import read_network

warnings.filterwarnings("ignore", category=RuntimeWarning)

OUT = DATA_DIR / "published_results" / "liver_timecourse"
OUT.mkdir(parents=True, exist_ok=True)
TIMEPOINT = 240
LOG2FC_THRESHOLD = 0.585        # 1.5-fold, as the Kuroda laboratory defines it
QVALUE_THRESHOLD = 0.1

fetch_kokaji_tables()
source = DATA_DIR / "mouse" / "latest"
network = read_network(str(source / "interactions.csv"),
                       nodes_file=str(source / "nodes.csv"))
print(f"{network.number_of_nodes():,} molecules, "
      f"{network.number_of_edges():,} relationships")

# %%
from transnet.visualization import plot_layer_connectivity

connectivity = cross_layer_connectivity(network)
print(f"{connectivity['cross_layer_fraction']:.0%} of edges cross between layers")
layer_figure = plot_layer_connectivity(connectivity,
                                       title="Mouse network: edges between layers")
plt.show()

# %% [markdown]
# ## The contrasts
#
# No statistics are recomputed. The supplement gives a fold change, p-value
# and q-value per time point and genotype, so each contrast is read from it and
# the ratio converted to log2. The analysis differs from the paper only in what
# is done with the network, not in the statistics. The main contrast is 240
# minutes against 0 minutes.
#
# The metabolites come with KEGG compound ids, which the network also uses, so
# no name matching is needed.

# %%
contrasts = {}
for genotype in ("WT", "ob/ob"):
    contrasts[genotype] = {
        "Metabolome": kokaji_contrast("S1", genotype, TIMEPOINT),
        "Transcriptome": kokaji_contrast("S4", genotype, TIMEPOINT),
    }

pd.DataFrame([
    {"genotype": genotype, "layer": layer, "features": len(table),
     "q <= 0.1": int((table["padj"] <= QVALUE_THRESHOLD).sum()),
     "and 1.5-fold": int(((table["padj"] <= QVALUE_THRESHOLD)
                          & (table["log2FC"].abs() >= LOG2FC_THRESHOLD)).sum())}
    for genotype, tables in contrasts.items() for layer, table in tables.items()
]).set_index(["genotype", "layer"])

# %% [markdown]
# The first result is visible before any network is used: at 240 minutes more
# metabolites change in wild-type than in obese liver, and more genes change in
# obese than in wild-type liver. This is the paper's main finding, and it is
# already in the raw numbers.
#
# ## The response over time
#
# The same comparison at every time point. Each column is one genotype at one
# time point; each row is a metabolite that changed at some time point. Rows are
# sorted by their largest change.

# %%
from transnet.visualization import plot_similarity_heatmap, plot_values_heatmap

over_time, changed_sets = {}, {}
for genotype in ("WT", "ob/ob"):
    for minutes in KOKAJI_TIMEPOINTS:
        column = f"{genotype} {minutes} min"
        for layer, sheet in (("Metabolome", "S1"), ("Transcriptome", "S4")):
            contrast = kokaji_contrast(sheet, genotype, minutes)
            hit = contrast[(contrast["padj"] <= QVALUE_THRESHOLD)
                           & (contrast["log2FC"].abs() >= LOG2FC_THRESHOLD)]
            changed_sets.setdefault(layer, {})[column] = set(hit["feature"])
            if layer == "Metabolome":
                over_time[column] = contrast.set_index("feature")["log2FC"]

values = pd.DataFrame(over_time)
values = values.loc[values.index.isin(set().union(*changed_sets["Metabolome"].values()))]
timecourse_figure = plot_values_heatmap(values, graph=network, max_rows=35,
                                        title="Changed metabolites, by genotype and time")
plt.show()

# %% [markdown]
# How similar are the responses? The Jaccard index of two sets of changed genes
# is the number they share divided by the number in either: 1 means the same
# genes changed, 0 means none in common.

# %%
columns = list(changed_sets["Transcriptome"])
jaccard = pd.DataFrame(
    [[len(changed_sets["Transcriptome"][a] & changed_sets["Transcriptome"][b])
      / max(len(changed_sets["Transcriptome"][a] | changed_sets["Transcriptome"][b]), 1)
      for b in columns] for a in columns],
    index=columns, columns=columns)
similarity_figure = plot_similarity_heatmap(
    jaccard, title="Overlap of the changed genes between time points and genotypes")
plt.show()

# %% [markdown]
# ## Identifiers
#
# Genes are identified by Ensembl ids and the network uses Entrez ids.
# Translating them takes thousands of lookups, so the result is cached.

# %%
def cached_map(name, features, from_type, to_type):
    path = OUT / f"idmap_{name}.json"
    if path.exists():
        return json.loads(path.read_text())
    from transnet.analysis import id_mapping
    mapped = id_mapping(pd.DataFrame({"feature": sorted(features)}), id_col="feature",
                        from_type=from_type, to_type=to_type, organism="mouse")
    column = f"{to_type}_id"
    result = {k: str(v) for k, v in zip(mapped["feature"], mapped[column])
              if pd.notna(v)}
    path.write_text(json.dumps(result))
    return result


genes = set(contrasts["WT"]["Transcriptome"]["feature"])
metabolite_ids = set(contrasts["WT"]["Metabolome"]["feature"])
id_map = {
    "Transcriptome": cached_map("ensembl_to_entrez", genes, "ensembl.gene", "entrezgene"),
    "Metabolome": {i: i for i in metabolite_ids},      # already KEGG compounds
}
{layer: len(m) for layer, m in id_map.items()}

# %% [markdown]
# ## A network per genotype

# %%
graphs, reports = {}, []
for genotype, tables in contrasts.items():
    graph = network.copy()
    report = map_omics_to_network(
        graph, tables, id_column="feature", log2fc_column="log2FC",
        qvalue_column="padj", id_map=id_map,
        qvalue_threshold=QVALUE_THRESHOLD, log2fc_threshold=LOG2FC_THRESHOLD,
    )
    graphs[genotype] = graph
    reports.append(report.per_layer.assign(genotype=genotype))

pd.concat(reports)[["genotype", "layer", "n_supplied", "n_matched",
                    "match_fraction", "n_up", "n_down"]].reset_index(drop=True)

# %%
layer_coverage(graphs["WT"])

# %% [markdown]
# ## Which axis regulates each reaction
#
# There is no proteome here, so the enzyme axis rests on transcripts alone.
# `gene_axis_evidence` records this (`gene` rather than `gene_protein`), and it
# is why these counts cannot be compared directly with a study that measured
# proteins.

# %%
regulation, rows = {}, []
for genotype, graph in graphs.items():
    table = reaction_regulation_table(graph)
    regulation[genotype] = table
    regulated = table[(table["gene_axis"] != 0) | (table["metabolite_axis"] != 0)]
    rows.append({
        "genotype": genotype,
        "regulated": len(regulated),
        "enzyme only": int(((regulated["gene_axis"] != 0)
                            & (regulated["metabolite_axis"] == 0)).sum()),
        "metabolite only": int(((regulated["gene_axis"] == 0)
                                & (regulated["metabolite_axis"] != 0)).sum()),
        "both": int(((regulated["gene_axis"] != 0)
                     & (regulated["metabolite_axis"] != 0)).sum()),
        "controversial": int(regulated["controversial"].sum()),
    })
    table.to_csv(OUT / f"reaction_regulation_{genotype.replace('/', '')}.csv",
                 index=False)

axes = pd.DataFrame(rows).set_index("genotype")
axes

# %% [markdown]
# This is the claim the paper is built on: healthy liver responds to glucose
# through its metabolites, obese liver through gene expression. The 19-gene
# panel could only test it partly. Here the whole liver is available and the two
# axes can be counted in each genotype.

# %%
from transnet.visualization import plot_axis_composition

composition = pd.DataFrame([
    {"contrast": genotype,
     "enzyme_axis_only": row["enzyme only"],
     "metabolite_axis_only": row["metabolite only"],
     "both_axes": row["both"],
     "controversial": row["controversial"]}
    for genotype, row in axes.iterrows()
])
plot_axis_composition(composition)
plt.show()

# %% [markdown]
# ### The controversial reactions
#
# Where the two axes point in opposite directions. Each bar is one measured
# molecule's push on the reaction: to the right it speeds the reaction up, to
# the left it slows it down.

# %%
from transnet.visualization import plot_controversial_reactions

controversial_figures = {}
for genotype, graph in graphs.items():
    table = regulation[genotype]
    print(f"{genotype}: {int(table['controversial'].sum())} controversial reactions")
    controversial_figures[genotype] = plot_controversial_reactions(
        graph, table, max_enzymes=5,
        title=f"{genotype} liver: controversial reactions")
    plt.show()

# %% [markdown]
# ## Per-pathway balance, and the metabolites doing the regulating
#
# `kegg_reaction_pathways` assigns each reaction to its mouse KEGG pathways, and
# `regulation_axis_summary` counts, per pathway, the reactions each axis
# activates or inhibits. The figures show the pathways with the most regulated
# reactions in each genotype.

# %%
from transnet.api import kegg_reaction_pathways
from transnet.visualization import plot_regulation_axes

reactions = [n for n, d in network.nodes(data=True) if d.get("layer") == "Reactions"]
pathways = kegg_reaction_pathways(reactions, organism="mmu")
balance = {genotype: regulation_axis_summary(table, pathway_map=pathways)
           for genotype, table in regulation.items()}
pd.concat([b.assign(genotype=g) for g, b in balance.items()]).set_index("genotype").head(12)

# %%
axis_figures = {}
for genotype, summary in balance.items():
    top = summary[summary["pathway"] != "unassigned"].head(12)
    axis_figures[genotype] = plot_regulation_axes(
        top, title=f"{genotype} liver, 240 min after glucose: regulation by pathway")
    plt.show()

# %%
role_rows = []
for genotype, graph in graphs.items():
    result = regulatory_role_enrichment(metabolite_regulatory_roles(graph))
    counts = result["counts"]
    row = result["enrichment"].set_index("role").loc["any"]
    role_rows.append({
        "genotype": genotype,
        "changed metabolites": counts["n_differential"],
        "of those regulators": counts["n_differential_regulators"],
        "background share": round(counts["fraction_background_regulators"], 3),
        "q": round(float(row["q_value"]), 3),
    })
pd.DataFrame(role_rows).set_index("genotype")

# %%
from transnet.visualization import plot_metabolite_regulators

regulator_figure = plot_metabolite_regulators(
    regulatory_role_enrichment(metabolite_regulatory_roles(graphs["WT"])))
plt.show()

# %% [markdown]
# ## Saturation, and the currency metabolites
#
# A metabolite can be an annotated regulator of an enzyme and still be irrelevant
# to it if its concentration sits far from the affinity. Table S13 gives measured
# Km and Ki values with a saturation index per genotype, for ATP and NADP+.
#
# Those two are exactly what TransNet treats as *currency* metabolites and
# excludes from the metabolite axis, on the grounds that everything uses them and
# a change in ATP is not a regulatory statement about any particular reaction.
# The saturation index tests that decision rather than assuming it: if these
# enzymes are already near-saturated for ATP, then ATP moving changes nothing,
# and excluding it is right.

# %%
kinetics = load_kokaji_table("S13", sheet="Km, Ki (ATP, NADP+)", header_rows=1)
kinetics = kinetics.assign(
    wt=pd.to_numeric(kinetics["Saturation index of WT mice at 0min"], errors="coerce"),
    obese=pd.to_numeric(kinetics["Saturation index of ob/ob mice at 0min"],
                        errors="coerce"),
).dropna(subset=["wt", "obese"])

print(f"{len(kinetics)} enzyme-metabolite pairs with a measured constant")
kinetics.groupby(["Km or Ki", "Metabolite"])[["wt", "obese"]].mean().round(3)

# %% [markdown]
# A saturation index near 1 means the enzyme is saturated, so a further rise in
# that metabolite changes nothing. Comparing the genotypes says whether obesity
# moves an enzyme into or out of the range where its regulator still has an
# effect.

# %%
shifted = kinetics.assign(shift=kinetics["obese"] - kinetics["wt"])
shifted.nlargest(8, "shift")[["Km or Ki", "Metabolite", "EC number",
                              "Geometric mean", "wt", "obese", "shift"]]

# %% [markdown]
# Two aggregate means would hide the shape of this. Every pair, with the
# diagonal: a point above the line is an enzyme sitting closer to saturation in
# obese liver, where a further change in that metabolite has less effect.

# %%
from transnet.visualization.palette import DOWN, INK, MUTED, SURFACE, UP, style_axes

saturation_figure, saturation_axes = plt.subplots(figsize=(5.2, 5))
saturation_axes.plot([0, 1], [0, 1], color=MUTED, linewidth=1, zorder=1)
for constant, colour in (("Km", UP), ("Ki", DOWN)):
    subset = shifted[shifted["Km or Ki"] == constant]
    saturation_axes.scatter(subset["wt"], subset["obese"], s=34, alpha=0.75,
                            color=colour, edgecolor=SURFACE, linewidth=0.6,
                            label=f"{constant} (n = {len(subset)})", zorder=2)
above = int((shifted["shift"] > 0).sum())
saturation_axes.set_xlabel("saturation index, wild-type liver")
saturation_axes.set_ylabel("saturation index, obese liver")
saturation_axes.set_title(
    f"ATP and NADP+ saturation\n{above} of {len(shifted)} pairs sit closer to "
    f"saturation in obese liver", loc="left", fontsize=11, fontweight="bold",
    color=INK)
saturation_axes.legend(frameon=False, fontsize=9, loc="lower right")
style_axes(saturation_axes)
saturation_figure.tight_layout()
plt.show()

# %% [markdown]
# The same pairs, matched to the reactions the network regulated, so the
# saturation index attaches to a regulation call rather than to an EC number on
# its own.

# %%
# The EC number sits on the catalysis edge, not on the reaction node, because one
# reaction can be catalysed by enzymes of more than one EC class.
reaction_ec = {}
for enzyme, reaction, data in graphs["WT"].edges(data=True):
    if data.get("edge_type") != "catalysis":
        continue
    for ec in str(data.get("ec") or "").split(","):
        ec = ec.strip()
        if ec:
            reaction_ec.setdefault(reaction, set()).add(ec)

ec_saturation = (shifted.groupby(shifted["EC number"].astype(str).str.strip())
                 [["wt", "obese", "shift"]].mean().round(3))

rows = []
for reaction, ecs in reaction_ec.items():
    matched = [ec for ec in ecs if ec in ec_saturation.index]
    if not matched:
        continue
    row = regulation["WT"].loc[regulation["WT"]["reaction"] == reaction]
    if row.empty:
        continue
    row = row.iloc[0]
    if row["gene_axis"] == 0 and row["metabolite_axis"] == 0:
        continue
    values = ec_saturation.loc[matched].mean()
    rows.append({"reaction": reaction, "name": str(row["name"])[:44],
                 "ec": ", ".join(sorted(matched)),
                 "gene_axis": row["gene_axis"],
                 "metabolite_axis": row["metabolite_axis"],
                 "wt": round(values["wt"], 3), "obese": round(values["obese"], 3),
                 "shift": round(values["shift"], 3)})

regulated_with_kinetics = pd.DataFrame(rows)
print(f"{len(regulated_with_kinetics)} regulated reactions have a measured "
      f"saturation index")
regulated_with_kinetics.sort_values("shift", key=abs, ascending=False).head(10)

# %% [markdown]
# ## Transcription factors, against a published inference
#
# Elsewhere in this documentation there is nothing to check the
# transcription-factor ranking against. Here the authors inferred factors from
# the same transcriptome by motif enrichment over gene clusters (Table S7), so
# the two inferences can be compared.
#
# The methods differ. Motif enrichment asks which DNA sequence motifs are
# over-represented in the promoters of a gene cluster.
# `transcription_factor_activity` asks which factors' ChIP-Atlas target genes
# are over-represented among the genes that changed. When two different
# methods agree on the same data, that agreement is stronger evidence than
# either result alone.

# %%
activity = {}
for genotype, graph in graphs.items():
    table = transcription_factor_activity(graph, min_targets=5)
    activity[genotype] = table
    implicated = table[table["q_value"] <= 0.05] if not table.empty else table
    print(f"{genotype:<6} {len(implicated)} of {len(table)} factors implicated "
          f"at q <= 0.05")

activity["WT"].head(12)[["name", "n_targets", "n_responsive_targets", "n_up",
                         "n_down", "q_value", "inferred_activity",
                         "factor_regulated"]]

# %%
from transnet.visualization import plot_tf_activity

tf_figure = plot_tf_activity(activity["WT"],
                             title="Transcription factors, wild-type liver")
plt.show()

# %%
motifs_published = load_kokaji_table("S7", sheet="Motif")
enriched_columns = [c for c in motifs_published.columns if "isEnrichedMotif" in c]
enriched_mask = motifs_published[enriched_columns].apply(
    pd.to_numeric, errors="coerce").fillna(0).gt(0).any(axis=1)
their_factors = set(motifs_published.loc[enriched_mask, "Name"].astype(str).str.upper())
print(f"the authors call {len(their_factors)} motifs enriched in at least one "
      f"gene cluster, out of {len(motifs_published)} tested")

ours = activity["WT"]
ours_implicated = set(ours.loc[ours["q_value"] <= 0.05, "name"]
                      .astype(str).str.upper())
ours_tested = set(ours["name"].astype(str).str.upper())
shared = their_factors & ours_tested

pd.DataFrame([
    {"in both analyses": len(shared),
     "they call enriched, we test": len(their_factors & ours_tested),
     "we implicate and they do too": len(ours_implicated & their_factors),
     "we implicate, they do not": len(ours_implicated - their_factors),
     "they implicate, we do not": len(their_factors & ours_tested
                                      - ours_implicated)}
])

# %%
sorted(ours_implicated & their_factors)

# %% [markdown]
# ## Timing, against the authors' own half-response times
#
# The other studies compute half-response times with TransNet. This supplement
# provides the authors' own, per genotype and for both layers, so the question
# of whether well-connected molecules respond first can be asked with the
# paper's numbers.

# %%
timing = []
for genotype in ("WT", "ob/ob"):
    for layer, table in (("Metabolome", "S1"), ("Transcriptome", "S4")):
        contrast = kokaji_contrast(table, genotype, TIMEPOINT)
        if "t_half" not in contrast:
            continue
        timed = contrast.dropna(subset=["t_half"])
        timing.append({"genotype": genotype, "layer": layer,
                       "with a t-half": len(timed),
                       "median t-half (min)": round(timed["t_half"].median(), 1)})
pd.DataFrame(timing).set_index(["genotype", "layer"])

# %% [markdown]
# Morita *et al.* report that in healthy liver the best-connected molecules
# respond first, and that this relationship is lost in obesity. The same test,
# on the authors' half-response times and this network's connections:

# %%
from scipy import stats

degree_timing = []
for genotype, graph in graphs.items():
    for layer, table in (("Metabolome", "S1"), ("Transcriptome", "S4")):
        contrast = kokaji_contrast(table, genotype, TIMEPOINT)
        if "t_half" not in contrast:
            continue
        mapped = contrast.assign(
            node=contrast["feature"].map(lambda f: id_map[layer].get(f, f)))
        mapped = mapped[mapped["node"].isin(graph)].dropna(subset=["t_half"])
        if len(mapped) < 10:
            continue
        degrees = mapped["node"].map(dict(graph.degree()))
        rho, p = stats.spearmanr(degrees, mapped["t_half"])
        degree_timing.append({"genotype": genotype, "layer": layer,
                              "n": len(mapped), "rho": round(rho, 3),
                              "p": float(f"{p:.3g}")})
pd.DataFrame(degree_timing).set_index(["genotype", "layer"])

# %% [markdown]
# `temporal_network_structure` runs the same test once the half-response times
# are stored on the network's nodes, and `plot_temporal_structure` draws it.

# %%
from transnet import temporal_network_structure
from transnet.visualization import plot_temporal_structure

timed_graph = graphs["WT"].copy()
for layer, sheet in (("Metabolome", "S1"), ("Transcriptome", "S4")):
    contrast = kokaji_contrast(sheet, "WT", TIMEPOINT)
    if "t_half" not in contrast:
        continue
    for feature, t_half in zip(contrast["feature"], contrast["t_half"]):
        node = id_map[layer].get(feature, feature)
        if node in timed_graph and pd.notna(t_half):
            timed_graph.nodes[node]["t_half"] = float(t_half)
timing_structure = temporal_network_structure(timed_graph)
print(timing_structure["degree_vs_thalf"]["interpretation"])
timing_figure = plot_temporal_structure(timed_graph, timing_structure, time_unit="min",
                                        title="Wild-type liver: connections against timing")
plt.show()

# %% [markdown]
# ## Hubs

# %%
hub_tables = {}
for genotype, graph in graphs.items():
    hubs, without_binding, binding_share = hub_rankings(graph, top_percent=2)
    hub_tables[genotype] = hubs
    print(f"{genotype:<6} " + ", ".join(hubs.head(3)["name"].str[:34]))

hub_tables["WT"].head(10)[["name", "layer", "cross_layer_degree",
                           "n_layers_touched"]]

# %%
from transnet.visualization import plot_transomic_hubs

hub_figure = plot_transomic_hubs(hub_tables["WT"], title="Wild-type liver: molecules "
                                                          "connecting layers")
plt.show()

# %% [markdown]
# ## Signed paths and propagation

# %%
traced, traced_paths, influences, path_rows = {}, {}, {}, []
for genotype, graph in graphs.items():
    paths = trace_regulatory_paths(graph, target_layer="Metabolome", max_length=4,
                                   max_paths=20000)
    summary, agree_paths, tested_paths = path_verdicts(paths)
    traced[genotype] = (summary, agree_paths, tested_paths)
    traced_paths[genotype] = paths
    seeds = {n: float(d["log2fc"]) for n, d in graph.nodes(data=True)
             if d.get("layer") == "Transcriptome" and d.get("regulated")
             and d.get("log2fc") is not None}
    influence = downstream_influence(graph, seeds, target_layer="Metabolome")
    influences[genotype] = influence
    measured = influence[influence["observed"].fillna(0) != 0]
    agree = int(measured["agrees"].astype("boolean").fillna(False).sum())
    path_rows.append({
        "genotype": genotype,
        "paths": len(paths),
        "metabolites the paths score": tested_paths,
        "paths predict correctly": agree_paths,
        "propagation reaches": len(measured),
        "propagation correct": agree,
    })
pd.DataFrame(path_rows).set_index("genotype")

# %%
for genotype, (summary, agree_paths, tested_paths) in traced.items():
    if not tested_paths:
        print(f"{genotype}: no path reaches a changed metabolite")
        continue
    print(f"{genotype}: {agree_paths} of {tested_paths} metabolites predicted "
          f"correctly: {versus_chance(agree_paths, tested_paths)}")

# %% [markdown]
# Thousands of paths reach the metabolome, but only 162 metabolites were
# measured, so few paths end at a metabolite whose change can be checked. The
# path figure therefore mostly shows predictions without a measurement to test
# them against. Propagation gives every measured metabolite a predicted
# direction, and `plot_downstream_influence` compares it with the measurement.

# %%
from transnet.visualization import plot_downstream_influence, plot_regulatory_paths

paths_figure = plot_regulatory_paths(traced_paths["WT"], graphs["WT"], top_n=8)
plt.show()

# %%
influence_figure = plot_downstream_influence(
    influences["WT"], title="Wild-type liver: metabolites predicted from the transcripts")
plt.show()

# %% [markdown]
# ## The wiring

# %%
wiring = []
for genotype, graph in graphs.items():
    motifs = regulatory_motifs(graph, responsive_only=True)
    convergence = convergence_significance(graph, n_randomisations=200)
    vulnerability = structural_vulnerability(responsive_subnetwork(graph), top_n=5)
    wiring.append({
        "genotype": genotype,
        **{k: v for k, v in motifs["counts"].items()},
        "convergent reactions": convergence["observed"],
        "expected": round(convergence["null_mean"], 1),
        "z": round(convergence["z"], 1),
        "p": float(f"{convergence['p_value']:.3g}"),
        "holds it together": ", ".join(vulnerability["name"].head(2))
        if not vulnerability.empty else "",
    })
pd.DataFrame(wiring).set_index("genotype")

# %% [markdown]
# The convergence test counts reactions where a changed *protein* and a
# changed metabolite meet. This study has no proteome, so the count is 0 by
# construction, not because the layers fail to converge. The figure is shown
# for completeness.

# %%
from transnet.visualization import (
    plot_convergence_null,
    plot_regulatory_motifs,
    plot_structural_vulnerability,
)

convergence_figure = plot_convergence_null(
    convergence_significance(graphs["WT"], n_randomisations=200),
    title="Wild-type liver: convergence against the shuffled null")
plt.show()

# %%
motif_figure = plot_regulatory_motifs(
    regulatory_motifs(graphs["WT"], responsive_only=True),
    title="Wild-type liver: regulatory motifs")
plt.show()

# %% [markdown]
# The molecules whose removal would split the wild-type response into separate
# pieces.

# %%
vulnerability_figure = plot_structural_vulnerability(
    structural_vulnerability(responsive_subnetwork(graphs["WT"]), top_n=10),
    title="Wild-type liver: molecules that hold the response together")
plt.show()

# %% [markdown]
# ## The two genotypes as networks

# %%
comparison = compare_transomic_networks(
    responsive_subnetwork(graphs["WT"]), responsive_subnetwork(graphs["ob/ob"]),
    "WT", "ob/ob",
)
print(f"edge Jaccard between the two responses: "
      f"{comparison['summary']['edge_jaccard']:.2f}")
comparison["edges_by_type"]

# %%
comparison["nodes_by_layer"]

# %%
from transnet.visualization import plot_condition_comparison

comparison_figure = plot_condition_comparison(comparison, "WT", "ob/ob")
plt.show()

# %% [markdown]
# ## The published claims, checked
#
# The three claims from this paper that the 19-gene panel could only partly
# test, now tested on the data they were made from.

# %%
wt, obese = axes.loc["WT"], axes.loc["ob/ob"]
metabolite_share = {g: (row["metabolite only"] + row["both"]) / max(row["regulated"], 1)
                    for g, row in axes.iterrows()}
enzyme_share = {g: (row["enzyme only"] + row["both"]) / max(row["regulated"], 1)
                for g, row in axes.iterrows()}

claims = pd.DataFrame([
    {"claim": "Healthy hepatic glucose responses rely on regulation by metabolites",
     "evidence": f"{metabolite_share['WT']:.0%} of regulated reactions in WT "
                 f"carry a changed metabolite"},
    {"claim": "In ob/ob liver, regulation by metabolites is lost",
     "evidence": f"{metabolite_share['ob/ob']:.0%} in ob/ob against "
                 f"{metabolite_share['WT']:.0%} in WT; "
                 f"{obese['metabolite only'] + obese['both']} reactions against "
                 f"{wt['metabolite only'] + wt['both']}"},
    {"claim": "ob/ob glucose responses depend instead on slow gene expression",
     "evidence": f"{enzyme_share['ob/ob']:.0%} of ob/ob regulated reactions "
                 f"carry a changed transcript against {enzyme_share['WT']:.0%} "
                 f"in WT, over {len(contrasts['ob/ob']['Transcriptome']):,} genes"},
])
claims.to_csv(OUT / "published_claims.csv", index=False)
claims

# %% [markdown]
# ## Figures and exports
#
# Every figure is saved, and the wild-type responsive network is exported in
# every format TransNet writes (see the *saving, exporting and sharing*
# walkthrough).

# %%
from transnet.io import to_arena3d, to_cytoscape_json, to_transomics2cytoscape, write_network
from transnet.visualization import (
    plot_transomic_network,
    plot_transomic_network_interactive,
    transomic_backbone,
)

for figure, name in [
    (layer_figure, "layer_connectivity"),
    (timecourse_figure, "metabolites_over_time"),
    (similarity_figure, "gene_overlap"),
    (controversial_figures["WT"], "controversial_WT"),
    (controversial_figures["ob/ob"], "controversial_obob"),
    (axis_figures["WT"], "regulation_axes_WT"),
    (axis_figures["ob/ob"], "regulation_axes_obob"),
    (timing_figure, "temporal_structure"),
    (hub_figure, "transomic_hubs"),
    (influence_figure, "downstream_influence"),
    (vulnerability_figure, "structural_vulnerability"),
    (plot_transomic_network(responsive_subnetwork(graphs["WT"]),
                            title="Wild-type liver, 240 min after glucose"),
     "network_WT"),
    (plot_transomic_network(responsive_subnetwork(graphs["ob/ob"]),
                            title="ob/ob liver, 240 min after glucose"),
     "network_obob"),
    (plot_axis_composition(composition), "axes_by_genotype"),
    (tf_figure, "tf_activity"),
    (regulator_figure, "metabolite_regulators"),
    (comparison_figure, "genotypes_compared"),
    (convergence_figure, "convergence_null"),
    (motif_figure, "regulatory_motifs"),
    (saturation_figure, "saturation"),
    (paths_figure, "regulatory_paths"),
]:
    figure.savefig(OUT / f"{name}.png", dpi=150, bbox_inches="tight")
plt.show()

responsive = responsive_subnetwork(graphs["WT"])
exports = OUT / "exports"
exports.mkdir(exist_ok=True)
write_network(responsive, str(exports / "csv"))
to_cytoscape_json(responsive, exports / "network_cytoscape.json")
to_arena3d(responsive, exports / "network_arena3d.json")
to_transomics2cytoscape(responsive, zip_path=exports / "network_transomics2cytoscape.zip")
# a browser draws a few hundred nodes well; export the backbone the figures use
html_network = transomic_backbone(graphs["WT"], max_reactions=40)
plot_transomic_network_interactive(html_network, layout="layered",
                                   title="Wild-type liver, 240 min").write_html(
    exports / "network.html")
print(f"figures and tables in {OUT}")
print("exports:", ", ".join(sorted(p.name for p in exports.iterdir())))
