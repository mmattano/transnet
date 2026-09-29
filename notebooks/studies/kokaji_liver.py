# %% [markdown]
# # Kokaji liver: a genome-wide glucose response, in two genotypes
#
# Wild-type and leptin-deficient obese (ob/ob) mouse liver, sampled at 0, 20,
# 60, 120 and 240 minutes after an oral glucose bolus. The transcriptome covers
# 14,292 genes and the metabolome 162 compounds, both with the authors' own fold
# changes, q-values and half-response times, and eleven signalling proteins were
# measured by western blot.
#
# `obese_liver.py` scores three claims from this study on a 19-gene panel at one
# timepoint and has to mark them "not reproduced here". This notebook uses the
# study's own data, so the same claims can be tested with the right instrument.
#
# Two things here exist nowhere else in this documentation: the transcription
# factor inference is scored against a published inference on the same data, and
# the metabolite axis is made quantitative with measured enzyme affinities.
#
# The measurements belong to Kokaji *et al.*, *Science Signaling*
# 13(660):eaaz1236, 2020. They are downloaded at run time into a directory that
# is not part of this repository and never redistributed.

# %%
import json
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from transnet import (
    build_name_id_map,
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
    transomic_hubs,
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

OUT = DATA_DIR / "published_results" / "kokaji_liver"
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
connectivity = cross_layer_connectivity(network)
print(f"{connectivity['cross_layer_fraction']:.0%} of edges cross between layers")
connectivity["matrix"]

# %% [markdown]
# ## The contrasts
#
# Nothing is recomputed. The supplement gives a fold change, p-value and q-value
# per timepoint per genotype, so the contrast is read off it and converted from a
# ratio to log2. What differs between this analysis and the paper is the network
# reading, not the statistics.
#
# Metabolites arrive with KEGG compound ids, which is unusual and useful: no name
# matching is needed for the layer that name matching usually limits.

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
# The first result is visible before any network: at 240 minutes the wild-type
# metabolome moves more than the obese one (36 metabolites against 11), while the
# obese transcriptome moves more than the wild-type one (552 genes against 425).
# That is the paper's headline, and it comes out of the raw contrast.

# %% [markdown]
# ## Identifiers
#
# Genes are Ensembl and the network is keyed by Entrez. The map is thousands of
# lookups, so it is cached.

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
# `gene_axis_evidence` records that, and it is why these counts are not
# comparable with a study that measured proteins.

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
# This is the claim the paper is built on, and the 19-gene panel in
# `obese_liver.py` could not test it: healthy liver answers glucose through
# metabolites, obese liver through gene expression. Here the whole liver is
# available and the two axes can be counted in each genotype.

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
# ## Per-pathway balance, and the metabolites doing the regulating

# %%
pd.concat([regulation_axis_summary(t).assign(genotype=g)
           for g, t in regulation.items()]).set_index("genotype")

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
# Every other transcription-factor result in this documentation is unvalidated:
# the ranking looks plausible and nothing says whether it is right. Here the
# authors inferred factors from the same transcriptome by motif enrichment over
# gene clusters (Table S7), so the two inferences can be compared.
#
# They are not the same method. Motif enrichment asks which sequences are
# over-represented in a cluster's promoters; `transcription_factor_activity` asks
# which ChIP-Atlas-bound target sets are over-represented among responsive genes.
# Agreement between two different methods on one dataset is worth more than
# either alone.

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
# Every other timing result here uses half-response times this package computed.
# The supplement distributes them, per genotype, for both layers, so the
# degree-versus-timing question can be asked with the paper's own numbers.

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
# respond first, and that the relationship is lost in obesity. The same test, on
# the authors' t-half values and this network's degrees:

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
# ## Hubs

# %%
hub_tables = {}
for genotype, graph in graphs.items():
    hubs, without_binding, binding_share = hub_rankings(graph, top_percent=2)
    hub_tables[genotype] = hubs
    print(f"{genotype:<6} " + ", ".join(hubs.head(3)["name"].str[:34]))

hub_tables["WT"].head(10)[["name", "layer", "cross_layer_degree",
                           "n_layers_touched"]]

# %% [markdown]
# ## Signed paths and propagation

# %%
traced, traced_paths, path_rows = {}, {}, []
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
# 5,788 paths reach the metabolome in wild-type liver and none of them scores a
# changed metabolite, so the figure shows what the paths predict with no
# measurement to check it against. That is the honest picture of a layer whose
# coverage runs out.

# %%
from transnet.visualization import plot_regulatory_paths

paths_figure = plot_regulatory_paths(traced_paths["WT"], graphs["WT"], top_n=8)
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
# The convergence null, drawn. Wild-type liver has no proteome here, so the
# enzyme side of a convergence is never measured and the count is 0 by
# construction rather than by measurement.

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
# ## The published claims, scored
#
# The three claims `obese_liver.py` had to leave open, tested on the data they
# were made from.

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
# ## Figures

# %%
from transnet.visualization import plot_transomic_network

for figure, name in [
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
print(f"figures and tables in {OUT}")
