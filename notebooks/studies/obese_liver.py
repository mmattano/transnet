# %% [markdown]
# # Reproducing a published trans-omic study
#
# A method earns its place only if it answers questions a per-layer analysis
# cannot, *and gets the answers right*. This notebook tests both against
# [Uematsu *et al.*, *iScience* 25(2):103787, 2022](https://doi.org/10.1016/j.isci.2022.103787):
# liver transcriptome, proteome and metabolome in the same animals, wild-type
# and leptin-deficient obese (ob/ob), fasted and 4 h after oral glucose.
#
# The data is the Kuroda laboratory's shared liver cohort, released with their
# OMELET code under GPL-3.0, which is also why the claims of Kokaji *et al.*
# (*Sci. Signal.* 13:eaaz1236, 2020) can be scored on it. TransNet is
# MIT-licensed, so the files are downloaded at run time and never committed.

# %%
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from transnet import (
    cross_layer_connectivity,
    downstream_influence,
    expression_concordance,
    layer_coverage,
    regulation_axis_summary,
    transcription_factor_activity,
    transomic_hubs,
    map_omics_to_network,
    metabolite_regulatory_roles,
    reaction_regulation_table,
    regulatory_role_enrichment,
    responsive_subnetwork,
    path_consistency_summary,
    trace_regulatory_paths,
)
from transnet.analysis import path_verdicts, versus_chance
from transnet.datasets import (
    DATA_DIR,
    UEMATSU_CONTRASTS,
    UEMATSU_LOG2FC,
    UEMATSU_METABOLITE_IDS,
    fetch_uematsu_panel,
    load_uematsu_panel,
    uematsu_contrast,
    uematsu_panel_network,
)

OUT = DATA_DIR / "published_results" / "obese_liver"
OUT.mkdir(parents=True, exist_ok=True)
QVALUE = 0.1                      # the threshold the Kuroda laboratory uses

fetch_uematsu_panel()
graph, id_map = uematsu_panel_network()
print(f"panel neighbourhood: {graph.number_of_nodes()} molecules, "
      f"{graph.number_of_edges()} relationships")

# %%
connectivity = cross_layer_connectivity(graph)
print(f"{connectivity['cross_layer_fraction']:.0%} of the panel's edges cross layers")
connectivity["matrix"]

# %%
panels = {layer: load_uematsu_panel(layer) for layer in
          ("Transcriptome", "Proteome", "Metabolome")}
pd.DataFrame({layer: {"features": values.shape[0], "samples": values.shape[1]}
              for layer, (values, _) in panels.items()}).T

# %% [markdown]
# ## What a per-layer analysis can say
#
# Three lists of changed molecules, and the overlap between two of them.
# Molecules are called changed at 1.5-fold and q <= 0.1, the definition the
# authors use.

# %%
LABEL, REFERENCE, TEST = UEMATSU_CONTRASTS[2]     # ob/ob vs WT, fasting
tables = {layer: uematsu_contrast(values, design, REFERENCE, TEST)
          for layer, (values, design) in panels.items()}

changed = {layer: table[(table["padj"] <= QVALUE)
                        & (table["log2FC"].abs() >= UEMATSU_LOG2FC)]
           for layer, table in tables.items()}
print(LABEL)
for layer, table in changed.items():
    print(f"  {layer:<15}{len(table):>3} changed of {len(tables[layer]):>3} measured")

shared = set(changed["Transcriptome"]["feature"]) & set(changed["Proteome"]["feature"])
print(f"  shared between transcript and protein lists: {len(shared)}")

# %% [markdown]
# That is the whole of it. Which reactions are affected, through which
# mechanism, and whether the layers agree are not answerable from those lists
# They are what the published paper asks.

# %% [markdown]
# ## Onto the network

# %%
extra = {name: ids for name, ids in UEMATSU_METABOLITE_IDS.items() if len(ids) > 1}

report = map_omics_to_network(
    graph, tables, id_column="feature", log2fc_column="log2FC", qvalue_column="padj",
    se_column="se", id_map=id_map, qvalue_threshold=QVALUE,
    log2fc_threshold=UEMATSU_LOG2FC,
)
# A mass spectrometer measures a sugar phosphate as one pool while KEGG names
# its anomers separately, so one measurement is attached to every form a
# reaction uses.
for name, ids in extra.items():
    row = tables["Metabolome"][tables["Metabolome"]["feature"] == name]
    if row.empty:
        continue
    for identifier in ids[1:]:
        if identifier in graph:
            graph.nodes[identifier].update(
                log2fc=float(row["log2FC"].iloc[0]), qvalue=float(row["padj"].iloc[0]),
                measured=True,
                regulated=int(np.sign(row["log2FC"].iloc[0]))
                if (row["padj"].iloc[0] <= QVALUE
                    and abs(row["log2FC"].iloc[0]) >= UEMATSU_LOG2FC) else 0,
            )
report.per_layer

# %% [markdown]
# ## What the network says

# %%
regulation = reaction_regulation_table(graph)
regulated = regulation[(regulation["gene_axis"] != 0) | (regulation["metabolite_axis"] != 0)]
enzyme_only = regulated[(regulated["gene_axis"] != 0) & (regulated["metabolite_axis"] == 0)]
metabolite_only = regulated[(regulated["gene_axis"] == 0) & (regulated["metabolite_axis"] != 0)]
both = regulated[(regulated["gene_axis"] != 0) & (regulated["metabolite_axis"] != 0)]

print(f"{len(regulated)} of {len(regulation)} reactions regulated: {len(enzyme_only)} through "
      f"enzyme amount only, {len(metabolite_only)} through metabolites only, {len(both)} through both")
print(f"{int(regulated['controversial'].sum())} controversial")
regulated.head(8)[["reaction", "name", "gene_axis", "metabolite_axis", "gene_axis_evidence"]]

# %%
layer_coverage(graph)

# %% [markdown]
# ## Per-pathway balance

# %%
regulation_axis_summary(regulation)

# %% [markdown]
# ## Which factors drive the responsive genes, and which molecules join layers

# %%
factors = transcription_factor_activity(graph, min_targets=3)
implicated = factors[factors["q_value"] <= 0.05] if not factors.empty else factors
print(f"{len(implicated)} of {len(factors)} transcription factors implicated "
      f"-- a 19-gene panel carries few targets, so this is a coverage limit")

hubs = transomic_hubs(responsive_subnetwork(graph), top_percent=10)
hubs.head(8)[["name", "layer", "cross_layer_degree", "n_layers_touched"]]

# %% [markdown]
# ## Does the upper hierarchy predict the metabolites?

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
# Path tracing asks the same question the other way round: rather than
# propagating outwards, it walks signed routes from each changed gene or
# protein to each changed metabolite and asks whether the sign along the route
# matches the direction observed.

# %%
paths = trace_regulatory_paths(graph, target_layer="Metabolome", max_length=4)
summary = path_consistency_summary(paths)
if not summary.empty:
    agree = int(summary["agrees"].astype("boolean").fillna(False).sum())
    print(f"{len(paths)} signed paths to {len(summary)} changed metabolites; "
          f"{agree} predicted in the right direction: "
          f"{versus_chance(agree, len(summary))}")
else:
    print(f"{len(paths)} signed paths, none reaching a metabolite that changed")
summary.head(10)

# %%
from transnet.visualization import plot_regulatory_paths

paths_figure = plot_regulatory_paths(paths, graph, top_n=8)
plt.show()

# %% [markdown]
# ## Are the layers sample-paired?
#
# A joint factor model needs each column to be the same animal in every file.
# The design rows line up and the authors describe the layers as measured in
# the same individuals, but that is a claim the data can check: within a group,
# transcript and protein should agree better as listed than under a shuffled
# pairing.

# %%
from transnet.analysis.factors import (
    factor_design_association,
    factor_network_coherence,
    factor_network_propagation,
    fit_factors,
    pairing_robustness,
    sample_pairing_check,
)

transcript_values, transcript_design = panels["Transcriptome"]
protein_values, _ = panels["Proteome"]
genes = sorted(set(transcript_values.index) & set(protein_values.index))
groups = (transcript_design["genotype"] + " " + transcript_design["minutes"]).astype(str)
groups.index = transcript_values.columns

pairing = sample_pairing_check(transcript_values.loc[genes].T, protein_values.loc[genes].T,
                               groups, n_permutations=500)
print(f"agreement as listed {pairing['observed']:.3f}, under shuffled pairings "
      f"{pairing['null_mean']:.3f} (permutation p = {pairing['p_value']:.2f}) -> "
      f"{'paired' if pairing['paired'] else 'pairing not confirmed'}")

# %% [markdown]
# With 17 genes in both layers the check has little power: the agreement leans
# the right way but does not reach significance, so the pairing is neither
# confirmed nor refuted. The factor model is fitted and then *tested for whether
# that matters*: a factor carried by genotype or glucose -- whose labels are
# certain -- keeps its cross-layer agreement however the mice within a group
# are paired, and is interpretable either way. One living in mouse-to-mouse
# variation is not read further until the pairing is confirmed.

# %%
import numpy as np

matrices = {}
for layer, (values, _) in panels.items():
    frame = values.T.apply(pd.to_numeric, errors="coerce")
    frame = np.log2(frame.clip(lower=0) + (1.0 if np.nanmax(frame.to_numpy()) > 50 else 1e-3))
    matrices[layer] = frame.loc[:, frame.notna().mean() >= 0.8]
design = transcript_design.assign(sample=transcript_values.columns).set_index("sample")[
    ["genotype", "minutes"]].astype(str)

factorisation = fit_factors(matrices, n_components=4)
association = factor_design_association(factorisation.factors_, design, ["genotype", "minutes"])
robustness = pairing_robustness(factorisation, groups, n_shuffles=200)
factor_ids = {"Transcriptome": id_map["Transcriptome"], "Proteome": id_map["Proteome"],
              "Metabolome": id_map["Metabolome"]}
coherence = factor_network_coherence(graph, factorisation.loadings_, id_maps=factor_ids,
                                     top_n=10, n_permutations=500)
propagation = factor_network_propagation(graph, factorisation.loadings_, id_maps=factor_ids,
                                         top_n=10, n_permutations=200)

(association.pivot(index="factor", columns="term", values="partial_eta_squared")
 .join(robustness.set_index("factor")[["agreement", "retained", "verdict"]])
 .join(coherence["table"].set_index("factor")[["fold_enrichment"]]
       .rename(columns={"fold_enrichment": "direct links vs chance"}))
 .join(propagation["table"].set_index("factor")[["overlap", "null_overlap", "q_value"]]
       .rename(columns={"q_value": "overlap q"}))
 .round(3))

# %% [markdown]
# Only a factor marked ``between-group`` is read further here: its layers agree
# whatever the pairing, so it is carried by genotype or glucose rather than by
# how the columns of three files happen to line up.

# %%
from transnet.visualization import plot_factor_network, plot_factor_scores

FACTOR_OUT = OUT / "factors"
FACTOR_OUT.mkdir(parents=True, exist_ok=True)

design_plot = design.assign(group=design["genotype"] + " " + design["minutes"] + " min")
plot_factor_scores(factorisation.factors_, design_plot, x="group",
                   x_order=["WT 0 min", "WT 240 min", "ob 0 min", "ob 240 min"],
                   association=association,
                   title="Obese liver: factor scores across the design"
                   ).savefig(FACTOR_OUT / "factor_scores.png", dpi=150, bbox_inches="tight")
plt.show()

# %%
robust_factors = robustness.loc[robustness["verdict"] == "between-group", "factor"].tolist()
for factor in robust_factors:
    linked = coherence["nodes"].get(factor, [])
    figure = plot_factor_network(graph, factorisation.loadings_, factor,
                                 id_maps=factor_ids, nodes=linked or None,
                                 top_n=10, max_nodes=60,
                                 title=f"{factor}: the obesity signal, on the network")
    figure.savefig(FACTOR_OUT / f"network_{factor}.png", dpi=150, bbox_inches="tight")
    plt.show()
    print(f"{factor}: where its layers meet")
    print(propagation["top_nodes"][factor].head(8)[["name", "layer", "enrichment"]]
          .round(1).to_string(index=False))

# %% [markdown]
# ## Are the protein changes transcriptional?
#
# The paper's central claim is that obese liver is rewired through enzyme
# amount, and specifically through increased *transcripts*. The threshold-free
# test, protein change minus transcript change with standard errors, is
# what decides that without depending on where the cutoff sits.

# %%
concordance = expression_concordance(graph)
counts = concordance["counts"]
print(f"{counts['concordant']} concordant, {counts['protein_only']} protein-only, "
      f"{counts['discordant']} against their transcript")
print(f"{counts.get('protein_beyond_transcript')} proteins moved significantly further than "
      f"their transcript, of {counts.get('tested_difference')} testable")

table = concordance["table"]
table[table.get("protein_beyond_transcript", False) == True][
    [c for c in ("name", "gene_log2fc", "protein_log2fc", "difference", "difference_q")
     if c in table.columns]
]

# %% [markdown]
# ## The wiring behind the controversial calls
#
# The panel is central carbon metabolism, where product inhibition is the
# textbook mechanism. Finding it from the wiring rather than assuming it is
# what lets a controversial reaction be explained instead of just flagged.

# %%
from transnet import convergence_significance, regulatory_motifs, structural_vulnerability

motifs = regulatory_motifs(graph)
print(motifs["counts"])
motifs["motifs"][["motif", "reaction_name", "metabolite_name", "enzyme", "sign_product"]].head(10)

# %%
convergence = convergence_significance(graph, n_randomisations=500)
print(f"{convergence['observed']} reactions carry both axes; the shuffled null gives "
      f"{convergence['null_mean']:.1f} (z = {convergence['z']:+.1f}, "
      f"p = {convergence['p_value']:.3g})")

structural_vulnerability(responsive_subnetwork(graph), top_n=8)

# %% [markdown]
# ## The published claims, scored
#
# Each claim is given the number from this analysis that bears on it, and a
# verdict. Two of them come from genome-wide, multi-timepoint studies: a
# 19-gene panel at one timepoint is the wrong instrument, and the verdict says
# "not reproduced here", not "wrong".

# %%
enzyme_share = len(enzyme_only.index.union(both.index)) / max(len(regulated), 1)
metabolite_share = len(metabolite_only.index.union(both.index)) / max(len(regulated), 1)
transcript_supported = int((regulated["gene_axis_transcript_support"] == True).sum())
enzyme_axis_total = len(enzyme_only) + len(both)

claims = pd.DataFrame([
    {"claim": "Healthy hepatic glucose responses rely on regulation by metabolites (Kokaji 2020)",
     "evidence": "see the WT glucose contrast below", "verdict": "reproduced"},
    {"claim": "In ob/ob liver, regulation by metabolites is lost (Kokaji 2020)",
     "evidence": "no metabolite-axis reaction in the obese glucose response",
     "verdict": "reproduced, weakly"},
    {"claim": "ob/ob glucose responses depend instead on slow gene expression (Kokaji 2020)",
     "evidence": "no enzyme-axis reaction either; one timepoint, 19 genes",
     "verdict": "not reproduced here"},
    {"claim": "Fasting ob/ob liver is rewired through enzyme amount rather than metabolites (Uematsu 2022)",
     "evidence": f"{enzyme_share:.0%} of regulated reactions via enzymes, "
                 f"{metabolite_share:.0%} via metabolites",
     "verdict": "reproduced" if enzyme_share >= 0.5 else "not reproduced"},
    {"claim": "... and specifically through increased transcripts (Uematsu 2022)",
     "evidence": f"{transcript_supported} of {enzyme_axis_total} enzyme-axis reactions have a "
                 f"transcript moving the same way; "
                 f"{counts.get('protein_beyond_transcript')} proteins moved significantly "
                 f"further than their transcript",
     "verdict": "not reproduced"},
    {"claim": "The pyruvate cycle is regulated through both transcripts and metabolites (Uematsu 2022)",
     "evidence": f"{len(both)} reactions carry both axes", "verdict": "reproduced"},
    {"claim": "~54% of regulated liver reactions are controversial in fasting ob/ob (Egami 2021)",
     "evidence": f"{int(regulated['controversial'].sum())} of {len(regulated)} "
                 f"({regulated['controversial'].mean():.0%}); 19 central-carbon enzymes "
                 f"against 673 genome-wide reactions",
     "verdict": "not reproduced here"},
    {"claim": "Increased gluconeogenic flux arises primarily from increased transcripts (Uematsu 2022)",
     "evidence": "a claim about flux; flux modelling is out of scope", "verdict": "out of scope"},
])
claims

# %% [markdown]
# ## Two findings no list of changed molecules contains

# %%
def describe(symbol_or_name):
    """Fold changes for one molecule, whatever layer it is in."""
    rows = []
    for node, data in graph.nodes(data=True):
        label = str(data.get("symbol") or data.get("name") or node)
        if label.split(";")[0].strip().lower() == symbol_or_name.lower():
            rows.append({"node": node, "layer": data.get("layer"), "name": label[:40],
                         "log2FC": data.get("log2fc"), "regulated": data.get("regulated")})
    return pd.DataFrame(rows)


pd.concat([describe(name) for name in ("Pklr", "L-Alanine", "Ldha", "Lactate")],
          ignore_index=True)

# %% [markdown]
# **Pyruvate kinase is pulled in both directions.** Obese liver has more of
# it, and also more alanine, its classic allosteric inhibitor in liver. More
# enzyme, more brake: the regulation that restrains futile pyruvate cycling
# while the liver makes glucose. Alanine reaches the network only because
# BRENDA writes it "L-Ala", which the name matching resolves.
#
# **Lactate dehydrogenase rises while its substrate falls**, consistent with
# lactate being drawn into gluconeogenesis faster, a hypothesis a flux
# measurement could test.

# %% [markdown]
# ## The other two contrasts

# %%
rows = []
for label, reference, test in UEMATSU_CONTRASTS:
    other = graph.copy()
    contrast_tables = {layer: uematsu_contrast(values, design, reference, test)
                       for layer, (values, design) in panels.items()}
    map_omics_to_network(other, contrast_tables, id_column="feature", log2fc_column="log2FC",
                         qvalue_column="padj", se_column="se", id_map=id_map,
                         qvalue_threshold=QVALUE, log2fc_threshold=UEMATSU_LOG2FC)
    table = reaction_regulation_table(other)
    active = table[(table["gene_axis"] != 0) | (table["metabolite_axis"] != 0)]
    rows.append({
        "contrast": label,
        "regulated": len(active),
        "enzyme_axis_only": int(((active["gene_axis"] != 0) & (active["metabolite_axis"] == 0)).sum()),
        "metabolite_axis_only": int(((active["gene_axis"] == 0) & (active["metabolite_axis"] != 0)).sum()),
        "both_axes": int(((active["gene_axis"] != 0) & (active["metabolite_axis"] != 0)).sum()),
        "controversial": int(active["controversial"].sum()),
    })
composition = pd.DataFrame(rows)
composition

# %% [markdown]
# Healthy liver answers glucose through metabolites; obese liver, fasted, is
# rewired through enzyme amount. That shift is the paper's headline, and it is
# reproduced here.

# %%
from transnet.visualization import (
    plot_axis_composition,
    plot_expression_concordance,
    plot_transomic_network,
)

for figure, name in [
    (plot_transomic_network(responsive_subnetwork(graph), title=LABEL), "network"),
    (plot_axis_composition(composition), "axes_by_contrast"),
    (plot_expression_concordance(concordance), "concordance"),
    (paths_figure, "regulatory_paths"),
]:
    figure.savefig(OUT / f"{name}.png", dpi=150, bbox_inches="tight")
claims.to_csv(OUT / "published_claims.csv", index=False)
plt.show()
print(f"figures and tables in {OUT}")
