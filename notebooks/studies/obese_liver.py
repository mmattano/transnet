# %% [markdown]
# # Obese liver after a glucose load: a 19-gene metabolic panel
#
# **Question.** Which reactions change in fasting obese liver, and through
# which regulatory mechanism? What can a trans-omic analysis say about a small
# targeted panel that separate lists of changed molecules cannot?
#
# **Data.** [Uematsu *et al.*, *iScience* 25(2):103787, 2022](https://doi.org/10.1016/j.isci.2022.103787)
# measured the liver transcriptome, proteome and metabolome in the same mice:
# wild-type and leptin-deficient obese (ob/ob) animals, fasted and 4 hours
# after an oral glucose load. The panel covers central carbon metabolism: 19
# enzyme genes, their proteins, and 32 metabolites.
#
# The data files are downloaded when the notebook runs.

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from transnet import (
    convergence_significance,
    cross_layer_connectivity,
    downstream_influence,
    expression_concordance,
    layer_coverage,
    map_omics_to_network,
    metabolite_regulatory_roles,
    path_consistency_summary,
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
from transnet.analysis import versus_chance
from transnet.api import kegg_reaction_pathways
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
from transnet.io import to_arena3d, to_cytoscape_json, to_transomics2cytoscape, write_network
from transnet.visualization import (
    plot_axis_composition,
    plot_controversial_reactions,
    plot_convergence_null,
    plot_downstream_influence,
    plot_expression_concordance,
    plot_layer_connectivity,
    plot_metabolite_regulators,
    plot_regulation_axes,
    plot_regulatory_motifs,
    plot_regulatory_paths,
    plot_structural_vulnerability,
    plot_transomic_hubs,
    plot_transomic_network,
    plot_transomic_network_interactive,
)

OUT = DATA_DIR / "published_results" / "obese_liver"
OUT.mkdir(parents=True, exist_ok=True)
QVALUE = 0.1                      # the threshold the Kuroda laboratory uses


def save(figure, name):
    figure.savefig(OUT / f"{name}.png", dpi=150, bbox_inches="tight")
    plt.show()


fetch_uematsu_panel()
graph, id_map = uematsu_panel_network()
print(f"panel neighbourhood: {graph.number_of_nodes()} molecules, "
      f"{graph.number_of_edges()} relationships")

# %% [markdown]
# The network is the part of the organism-wide mouse network around the panel:
# its genes, their proteins, the reactions those enzymes catalyse, and every
# metabolite on or regulating those reactions. With no protein interactions
# among 19 enzymes, every edge connects two different layers.

# %%
connectivity = cross_layer_connectivity(graph)
print(f"{connectivity['cross_layer_fraction']:.0%} of the panel's edges cross layers")
save(plot_layer_connectivity(connectivity, title="Obese liver panel: edges between layers"),
     "layer_connectivity")

# %%
panels = {layer: load_uematsu_panel(layer) for layer in
          ("Transcriptome", "Proteome", "Metabolome")}
pd.DataFrame({layer: {"features": values.shape[0], "samples": values.shape[1]}
              for layer, (values, _) in panels.items()}).T

# %% [markdown]
# ## What a per-layer analysis can say
#
# Three lists of changed molecules, and the overlap between two of them. A
# molecule counts as changed at 1.5-fold and q <= 0.1, the definition the
# authors use. The main contrast is obese against lean liver, fasted.

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
# That is all the separate lists can say. They cannot say which reactions are
# affected, through which mechanism, or whether the layers agree. The rest of
# this notebook answers those questions on the network.
#
# ## Mapping the data onto the network

# %%
extra = {name: ids for name, ids in UEMATSU_METABOLITE_IDS.items() if len(ids) > 1}

report = map_omics_to_network(
    graph, tables, id_column="feature", log2fc_column="log2FC", qvalue_column="padj",
    se_column="se", id_map=id_map, qvalue_threshold=QVALUE,
    log2fc_threshold=UEMATSU_LOG2FC,
)
# Mass spectrometry measures a sugar phosphate as one pool, while KEGG lists its
# anomers separately, so the measurement is attached to every form a reaction
# uses.
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
# `layer_coverage` shows the limit of this panel: nearly every transcript and
# protein is measured, but only a quarter of the metabolites that the panel's
# reactions touch.

# %%
layer_coverage(graph)

# %% [markdown]
# ## Which axis regulates each reaction
#
# For every reaction, `reaction_regulation_table` asks whether it changed
# through the amount of its enzyme (the enzyme or gene axis), through its
# substrates, products and allosteric regulators (the metabolite axis), or both,
# and whether the two axes agree.

# %%
regulation = reaction_regulation_table(graph)
regulated = regulation[(regulation["gene_axis"] != 0) | (regulation["metabolite_axis"] != 0)]
enzyme_only = regulated[(regulated["gene_axis"] != 0) & (regulated["metabolite_axis"] == 0)]
metabolite_only = regulated[(regulated["gene_axis"] == 0) & (regulated["metabolite_axis"] != 0)]
both = regulated[(regulated["gene_axis"] != 0) & (regulated["metabolite_axis"] != 0)]
controversial = regulated[regulated["controversial"]]

print(f"{len(regulated)} of {len(regulation)} reactions regulated: {len(enzyme_only)} through "
      f"enzyme amount only, {len(metabolite_only)} through metabolites only, "
      f"{len(both)} through both")
print(f"{len(controversial)} controversial (the two axes point opposite ways)")
regulation.to_csv(OUT / "reaction_regulation.csv", index=False)
regulated.head(8)[["reaction", "name", "gene_axis", "metabolite_axis", "gene_axis_evidence"]]

# %% [markdown]
# `gene_axis_evidence` shows what each enzyme-axis call rests on. With both
# transcript and protein measured for almost every enzyme, most calls are
# supported by both layers (`gene_protein`).

# %%
regulated["gene_axis_evidence"].value_counts(dropna=False).to_frame("reactions")

# %% [markdown]
# ### The controversial reactions
#
# Each bar is one measured molecule's push on the reaction: to the right it
# speeds the reaction up, to the left it slows it down. Yellow bars are the
# enzyme axis, pink bars the metabolite axis.

# %%
controversial[["reaction", "name", "gene_axis", "metabolite_axis", "allosteric_regulators"]]

# %%
save(plot_controversial_reactions(graph, regulation), "controversial_reactions")

# %% [markdown]
# ### Per-pathway balance
#
# `kegg_reaction_pathways` assigns each reaction to its KEGG pathways (mouse
# pathways only, overview maps excluded). `regulation_axis_summary` then counts,
# per pathway, how many reactions each axis activates or inhibits, and
# `plot_regulation_axes` draws it. Reactions that KEGG places in no mouse
# pathway are counted as `unassigned` and left out of the figure.

# %%
pathways = kegg_reaction_pathways(regulation["reaction"], organism="mmu")
pathway_balance = regulation_axis_summary(regulation, pathway_map=pathways)
pathway_balance.head(10)[["pathway", "n_reactions", "gene_activated", "gene_inhibited",
                          "metabolite_activated", "metabolite_inhibited", "n_controversial"]]

# %%
assigned = pathway_balance[(pathway_balance["pathway"] != "unassigned")
                           & (pathway_balance["n_reactions"] >= 2)]
save(plot_regulation_axes(assigned,
                          title="Obese against lean liver, fasted: regulation by pathway"),
     "regulation_axes")

# %% [markdown]
# ### Do the changed metabolites regulate enzymes?
#
# A changed metabolite may simply follow the change in flux, or it may itself
# act on an enzyme. `metabolite_regulatory_roles` lists, for each metabolite,
# the reactions it activates or inhibits according to BRENDA, and
# `regulatory_role_enrichment` tests whether changed metabolites are
# regulators more often than measured metabolites in general.

# %%
roles = regulatory_role_enrichment(metabolite_regulatory_roles(graph))
counts = roles["counts"]
print(f"{counts['n_differential_regulators']} of {counts['n_differential']} changed metabolites "
      f"regulate an enzyme ({counts['fraction_differential_regulators']:.0%}, against "
      f"{counts['fraction_background_regulators']:.0%} of all measured metabolites)")
roles["regulators"][["name", "log2fc", "reactions_activated", "reactions_inhibited"]]

# %%
save(plot_metabolite_regulators(roles), "metabolite_regulators")

# %% [markdown]
# ### Are the protein changes transcriptional?
#
# Most enzyme-axis reactions change through the amount of enzyme. Is that
# change transcriptional? For every gene measured in both layers,
# `expression_concordance` asks whether the protein
# followed its transcript. Because standard errors were mapped, it also tests
# directly whether the protein changed *more* than its transcript, which does
# not depend on where a significance cut-off falls.

# %%
concordance = expression_concordance(graph)
concordance_counts = concordance["counts"]
print(f"{concordance_counts['concordant']} concordant, "
      f"{concordance_counts['protein_only']} protein-only, "
      f"{concordance_counts['discordant']} against their transcript")
print(f"{concordance_counts.get('protein_beyond_transcript')} proteins moved significantly "
      f"further than their transcript, of {concordance_counts.get('tested_difference')} testable")

table = concordance["table"]
table[table.get("protein_beyond_transcript", False) == True][  # noqa: E712
    [c for c in ("name", "gene_log2fc", "protein_log2fc", "difference", "difference_q")
     if c in table.columns]
]

# %%
save(plot_expression_concordance(concordance), "concordance")

# %% [markdown]
# ## Transcription factors, and the molecules that connect layers
#
# The panel network has no transcription-factor edges that reach its 19
# genes, so the factor inference cannot run here; this is a limit of the panel,
# not a finding. The liver time-course study runs it genome-wide. The
# cross-layer hubs are informative: they are the molecules through which the
# response in one layer reaches another.

# %%
factors = transcription_factor_activity(graph, min_targets=3)
implicated = factors[factors["q_value"] <= 0.05] if not factors.empty else factors
print(f"{len(factors)} transcription factors testable, {len(implicated)} implicated")

hubs = transomic_hubs(responsive_subnetwork(graph), top_percent=10)
hubs.head(8)[["name", "layer", "cross_layer_degree", "n_layers_touched"]]

# %%
save(plot_transomic_hubs(hubs, title="Obese liver: molecules connecting layers"),
     "transomic_hubs")

# %% [markdown]
# ## Do the enzyme changes predict the metabolite changes?
#
# If the network is right, the changes in transcripts and proteins should
# predict which way the metabolites moved. Two methods ask this (see the
# *signed regulatory paths* walkthrough). Propagation pushes the enzyme changes
# forward along signed edges and gives each metabolite a predicted direction.

# %%
seeds = {n: float(d["log2fc"]) for n, d in graph.nodes(data=True)
         if d.get("layer") in ("Transcriptome", "Proteome")
         and d.get("regulated") and d.get("log2fc") is not None}
influence = downstream_influence(graph, seeds, target_layer="Metabolome")
measured = influence[influence["observed"].fillna(0) != 0]
if not measured.empty:
    agree = int(measured["agrees"].astype("boolean").fillna(False).sum())
    print(f"{agree} of {len(measured)} changed metabolites predicted correctly: "
          f"{versus_chance(agree, len(measured))}")
measured.head(10)[["name", "score", "predicted_direction", "observed", "agrees"]]

# %%
save(plot_downstream_influence(influence), "downstream_influence")

# %% [markdown]
# Path tracing lists the individual routes from each changed gene or protein to
# each changed metabolite, and checks whether the sign along each route matches
# the measured direction. `path_consistency_summary` gives one verdict per
# metabolite.

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
save(plot_regulatory_paths(paths, graph, top_n=8), "regulatory_paths")

# %% [markdown]
# ## Recurring wiring patterns, weak points, and a null model
#
# `regulatory_motifs` searches the signed wiring among the changed molecules
# for product inhibition, allosteric feedback and feed-forward patterns, which
# can explain a controversial reaction rather than just flag it. On this panel
# it finds none: the one allosteric regulator that changed, alanine, acts on
# pyruvate kinase, which does not produce it.

# %%
motifs = regulatory_motifs(graph)
print(motifs["counts"])
save(plot_regulatory_motifs(motifs), "regulatory_motifs")
motifs["motifs"][["motif", "reaction_name", "metabolite_name", "enzyme", "sign_product"]].head(10)

# %% [markdown]
# Is the number of reactions where both axes changed more than chance would
# produce? `convergence_significance` reassigns at random which molecules
# changed, keeping the number per layer fixed, and counts again.

# %%
convergence = convergence_significance(graph, n_randomisations=500)
print(f"{convergence['observed']} reactions with both axes changed; "
      f"{convergence['null_mean']:.1f} expected by chance "
      f"(z = {convergence['z']:+.1f}, p = {convergence['p_value']:.3g})")
save(plot_convergence_null(convergence), "convergence_null")

# %% [markdown]
# `structural_vulnerability` lists the molecules whose removal would split the
# responsive network into separate pieces.

# %%
vulnerability = structural_vulnerability(responsive_subnetwork(graph), top_n=8)
save(plot_structural_vulnerability(vulnerability), "structural_vulnerability")
vulnerability

# %% [markdown]
# ## Are the layers from the same animals?
#
# A joint factor model assumes that column *i* is the same mouse in every file.
# The authors describe the layers as measured in the same individuals, and the
# data can check this: within a group, a mouse whose transcript of a gene is
# high should also tend to have more of that protein, and this agreement
# should be stronger in the listed pairing than in shuffled ones.

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
# With only 17 genes in both layers the check has little power: the agreement
# points the right way but is not significant, so the pairing is neither
# confirmed nor refuted.
#
# The factor model is therefore fitted and then tested for whether the pairing
# matters. `pairing_robustness` shuffles the mice within each group and checks
# whether each factor keeps its agreement across layers. A factor marked
# `between-group` is carried by genotype or glucose, whose labels are certain,
# and can be interpreted whatever the pairing. A factor that depends on
# mouse-to-mouse variation is not interpreted until the pairing is confirmed.

# %%
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
# In the table, the genotype and minutes columns give the share of each factor's
# variation explained by that part of the design. `direct links vs chance` says
# how much more often the factor's top molecules are directly connected on the
# network than random molecules, and `overlap` how much the network
# neighbourhoods of its top molecules in different layers overlap.

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
# ## Two findings that no list of changed molecules contains

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
# **Pyruvate kinase is pushed in both directions.** Obese liver has more of the
# enzyme, and also more alanine, its classic allosteric inhibitor in liver. The
# extra enzyme is held back by the extra inhibitor. This is the regulation that
# limits wasteful cycling between pyruvate and phosphoenolpyruvate while the
# liver makes glucose. Alanine is found as a regulator only because the name
# matching resolves BRENDA's spelling, "L-Ala".
#
# **Lactate dehydrogenase rises while its substrate falls.** This is consistent
# with lactate being used faster for glucose production, a hypothesis that a
# flux measurement could test.

# %% [markdown]
# ## The other two contrasts
#
# The same analysis for the glucose response of lean and of obese liver.

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
# Healthy liver responds to glucose through its metabolites; fasted obese liver
# is changed through the amount of its enzymes.

# %%
save(plot_axis_composition(composition), "axes_by_contrast")

# %% [markdown]
# ## The responsive network

# %%
responsive = responsive_subnetwork(graph)
save(plot_transomic_network(responsive, title=LABEL), "network")

# %% [markdown]
# ## Exports
#
# The responsive network in every format TransNet writes, for use in other
# tools (see the *saving, exporting and sharing* walkthrough).

# %%
exports = OUT / "exports"
exports.mkdir(exist_ok=True)
write_network(responsive, str(exports / "csv"))
to_cytoscape_json(responsive, exports / "network_cytoscape.json")
to_arena3d(responsive, exports / "network_arena3d.json")
to_transomics2cytoscape(responsive, zip_path=exports / "network_transomics2cytoscape.zip")
plot_transomic_network_interactive(responsive, layout="layered", title=LABEL).write_html(
    exports / "network.html")
print(f"figures and tables in {OUT}")
print("exports:", ", ".join(sorted(p.name for p in exports.iterdir())))
