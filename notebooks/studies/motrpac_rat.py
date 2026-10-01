# %% [markdown]
# # MoTrPAC: endurance training in six rat tissues
#
# **Question.** Endurance training changes many tissues at once. Do they
# respond through the same regulatory mechanisms, and how much of the response
# lies in enzyme modification rather than enzyme amount?
#
# **Data.** The MoTrPAC consortium trained rats on a treadmill for 1, 2, 4 and 8
# weeks and measured transcriptome, proteome and metabolome in many tissues,
# with phosphorylation sites in every tissue used here and acetylation and
# ubiquitination sites in heart and liver
# ([*Nature* 629:174-183, 2024](https://doi.org/10.1038/s41586-023-06877-w)).
# Six tissues are mapped onto one rat network, so their responses can be
# compared directly: **one network, six different responses on it**.
#
# MoTrPAC distributes a time-course test (an ANOVA) whose statistic says that a
# molecule changed but not in which direction. The direction is therefore
# computed here from the normalised data: each time point against its sedentary
# controls, as the mean of the differences within each sex, with a Welch test on
# values centred per sex.
#
# Needs the rat network: `python maintenance/build_networks.py --organisms rat
# --brenda`.

# %%
import json
import matplotlib.pyplot as plt
import pandas as pd

from transnet import (
    build_name_id_map,
    compare_transomic_networks,
    expression_concordance,
    map_modification_sites,
    map_omics_to_network,
    reaction_regulation_table,
    responsive_subnetwork,
    trace_regulatory_paths,
    transcription_factor_activity,
)
from transnet import (                      # noqa: E402  (grouped for readability)
    cross_layer_connectivity,
    downstream_influence,
    layer_coverage,
    metabolite_regulatory_roles,
    regulation_axis_summary,
    regulatory_role_enrichment,
    temporal_network_structure,
    assign_temporal_parameters,
)
from transnet.analysis import hub_rankings, path_verdicts, versus_chance
from transnet.datasets import (
    DATA_DIR,
    MOTRPAC_PTM_ASSAYS,
    load_motrpac_contrast,
    load_motrpac_ptm,
    motrpac_tissues,
)
from transnet.io import read_network

OUT = DATA_DIR / "motrpac_results" / "transomics"
OUT.mkdir(parents=True, exist_ok=True)
TIMEPOINT = "8w"
TISSUE_FOCUS = "SKM_GN"

source = DATA_DIR / "rat" / "latest"
network = read_network(str(source / "interactions.csv"), nodes_file=str(source / "nodes.csv"))
tissues = motrpac_tissues()
print(f"{network.number_of_nodes():,} molecules, {network.number_of_edges():,} relationships")
print("tissues:", ", ".join(tissues))

# %%
from transnet.visualization import plot_layer_connectivity

connectivity = cross_layer_connectivity(network)
print(f"{connectivity['cross_layer_fraction']:.0%} of edges cross between layers")
layer_figure = plot_layer_connectivity(connectivity, title="Rat network: edges between layers")
plt.show()

# %% [markdown]
# ## Identifiers
#
# The transcripts are identified by Ensembl ids and the proteins by RefSeq
# accessions, while the network uses Entrez ids and UniProt accessions. The
# translations take thousands of lookups, so they are cached.

# %%
def cached_map(name, features, from_type, to_type):
    path = OUT / f"idmap_{name}.json"
    if path.exists():
        return json.loads(path.read_text())
    from transnet.analysis import id_mapping
    mapped = id_mapping(pd.DataFrame({"feature": sorted(features)}), id_col="feature",
                        from_type=from_type, to_type=to_type, organism="rat")
    column = f"{to_type}_id"
    result = {k: str(v) for k, v in zip(mapped["feature"], mapped[column]) if pd.notna(v)}
    path.write_text(json.dumps(result))
    return result


contrasts = {tissue: load_motrpac_contrast(tissue, TIMEPOINT) for tissue in tissues}

transcripts = {f for t in contrasts.values() for f in t.get("Transcriptome", pd.DataFrame()).get("feature", [])}
proteins = {f for t in contrasts.values() for f in t.get("Proteome", pd.DataFrame()).get("feature", [])}
metabolites = {f for t in contrasts.values() for f in t.get("Metabolome", pd.DataFrame()).get("feature", [])}

id_map = {
    "Transcriptome": cached_map("ensembl_to_entrez", transcripts, "ensembl.gene", "entrezgene"),
    "Proteome": cached_map("refseq_to_uniprot", proteins, "refseq", "uniprot"),
    "Metabolome": build_name_id_map(network, sorted(metabolites)),
}
{layer: len(m) for layer, m in id_map.items()}

# %% [markdown]
# ## Each tissue gets its own copy
#
# Every tissue is mapped onto its own copy of the same network: the network is
# shared, the measured response is not.

# %%
graphs, mapped = {}, []
for tissue, tables in contrasts.items():
    graph = network.copy()
    report = map_omics_to_network(
        graph, tables, id_column="feature", log2fc_column="log2FC", qvalue_column="padj",
        se_column="se", id_map=id_map, qvalue_threshold=0.05,
    )
    graphs[tissue] = graph
    row = report.per_layer.assign(tissue=tissue)
    mapped.append(row[["tissue", "layer", "n_supplied", "n_matched", "n_up", "n_down"]])

pd.concat(mapped).reset_index(drop=True)

# %% [markdown]
# `layer_coverage` shows, for the focus tissue (gastrocnemius muscle), how
# much of each network layer was measured and how much of it changed.

# %%
layer_coverage(graphs[TISSUE_FOCUS])

# %% [markdown]
# ## Modification sites
#
# Besides protein amounts, MoTrPAC measures three protein modifications, and
# each says something different about an enzyme. Phosphorylation often switches
# an enzyme's activity; it was measured in every tissue here. Acetylation and
# ubiquitination (which marks proteins for degradation) were measured in heart
# and liver only.
#
# These assays measure individual sites, not proteins. `map_modification_sites`
# therefore adds each site as its own node in the Signaling layer, with an edge
# to its protein; several sites on one protein stay separate because they can
# move in opposite directions. The edge has no sign: whether phosphorylation at
# a given site raises or lowers an enzyme's activity is rarely known. So the
# direction the site moved and its effect on the reaction are reported
# separately.

# %%
ptm_tables, ptm_reports = {}, []
for assay in ("PHOSPHO", "ACETYL", "UBIQ"):
    for tissue in tissues:
        if tissue not in MOTRPAC_PTM_ASSAYS[assay]:
            continue
        table = load_motrpac_ptm(tissue, assay, TIMEPOINT)
        if table.empty:
            continue
        ptm_tables[(tissue, assay)] = table
        changed = table[table["padj"] <= 0.05]
        ptm_reports.append({
            "tissue": tissue, "assay": assay, "sites": len(table),
            "proteins": table["protein"].nunique(),
            "changed sites": len(changed),
            "changed proteins": changed["protein"].nunique(),
        })

ptm_summary = pd.DataFrame(ptm_reports).set_index(["assay", "tissue"])
ptm_summary

# %% [markdown]
# Phosphorylation changes in every tissue. Ubiquitination hardly changes: 43
# sites in heart and 9 in liver at q <= 0.05, against 490 and 489
# phosphorylation sites. Whatever protein degradation contributes to the
# training response, it does not show in these ubiquitination sites at eight
# weeks.
#
# The sites are attached to the proteins through the same RefSeq-to-UniProt
# translation as the proteome, extended to the proteins that only the
# modification assays measured.

# %%
ptm_proteins = {p for table in ptm_tables.values() for p in table["protein"]}
site_id_map = dict(id_map["Proteome"])
missing = sorted(ptm_proteins - set(site_id_map))
if missing:
    site_id_map.update(cached_map("ptm_refseq_to_uniprot", missing, "refseq", "uniprot"))

attached = []
for (tissue, assay), table in ptm_tables.items():
    if assay != "PHOSPHO":
        continue        # only phosphorylation feeds the phospho axis; see below
    report = map_modification_sites(graphs[tissue], table, id_map=site_id_map)
    attached.append({"tissue": tissue, **report})

pd.DataFrame(attached).set_index("tissue")

# %% [markdown]
# The three modifications side by side, per tissue. Bars to the right count
# sites that went up, bars to the left sites that went down.

# %%
from transnet.visualization import plot_layer_changes

ptm_changes = pd.DataFrame([
    {"group": tissue, "layer": assay,
     "n_up": int((table.loc[table["padj"] <= 0.05, "log2FC"] > 0).sum()),
     "n_down": int((table.loc[table["padj"] <= 0.05, "log2FC"] < 0).sum())}
    for (tissue, assay), table in ptm_tables.items()
])
modification_figure = plot_layer_changes(
    ptm_changes, title="Changed modification sites per tissue")
plt.show()

# %% [markdown]
# ## Which axis regulates each reaction

# %%
regulation, rows = {}, []
for tissue, graph in graphs.items():
    table = reaction_regulation_table(graph)
    regulation[tissue] = table
    regulated = table[(table["gene_axis"] != 0) | (table["metabolite_axis"] != 0)]
    rows.append({
        "tissue": tissue,
        "regulated": len(regulated),
        "enzyme only": int(((regulated["gene_axis"] != 0) & (regulated["metabolite_axis"] == 0)).sum()),
        "metabolite only": int(((regulated["gene_axis"] == 0) & (regulated["metabolite_axis"] != 0)).sum()),
        "both": int(((regulated["gene_axis"] != 0) & (regulated["metabolite_axis"] != 0)).sum()),
        "controversial": int(regulated["controversial"].sum()),
        "transcript supported": int((regulated["gene_axis_transcript_support"] == True).sum()),
    })
    (OUT / "regulation").mkdir(parents=True, exist_ok=True)
    table.to_csv(OUT / "regulation" / f"{tissue}_reaction_regulation.csv", index=False)

axes = pd.DataFrame(rows).set_index("tissue")
axes

# %% [markdown]
# The tissues differ in the *kind* of regulation, not only in its amount.
# Tissues whose metabolome changed a lot regulate many reactions through
# metabolites; in tissues whose metabolome barely changed, the enzyme axis
# dominates by default.
#
# `plot_axis_composition` draws the table, one bar per tissue.

# %%
from transnet.visualization import plot_axis_composition

composition = axes.reset_index().rename(columns={
    "tissue": "contrast", "enzyme only": "enzyme_axis_only",
    "metabolite only": "metabolite_axis_only", "both": "both_axes",
})
composition_figure = plot_axis_composition(composition, label="contrast")
plt.show()

# %% [markdown]
# In gastrocnemius muscle, `gene_axis_evidence` shows what each enzyme-axis
# call rests on, and `plot_controversial_reactions` shows the reactions where
# the two axes pull in opposite directions.

# %%
focus_regulated = regulation[TISSUE_FOCUS][(regulation[TISSUE_FOCUS]["gene_axis"] != 0)
                                           | (regulation[TISSUE_FOCUS]["metabolite_axis"] != 0)]
focus_regulated["gene_axis_evidence"].value_counts(dropna=False).to_frame("reactions")

# %%
from transnet.visualization import plot_controversial_reactions

controversial_figure = plot_controversial_reactions(
    graphs[TISSUE_FOCUS], regulation[TISSUE_FOCUS], max_enzymes=6,
    title=f"{TISSUE_FOCUS}: controversial reactions after {TIMEPOINT} of training")
plt.show()

# %% [markdown]
# ## Enzyme amount against enzyme modification
#
# With the sites mapped, the regulation table has a phosphorylation axis beside
# the gene axis. The important reactions are those whose enzyme amount did not
# change while its phosphorylation did: an analysis of protein amounts alone
# would call them unregulated.

# %%
phospho_rows = []
for tissue, table in regulation.items():
    if "phospho_axis" not in table:
        continue
    with_sites = table[table["n_phosphosites_changed"] > 0]
    steady = with_sites[with_sites["gene_axis"] == 0]
    phospho_rows.append({
        "tissue": tissue,
        "reactions with a changed site": len(with_sites),
        "amount steady, site moved": len(steady),
        "one direction": int((with_sites["phospho_axis"] != 0).sum()),
        "sites disagree": int(((with_sites["phospho_axis"] == 0)
                               & (with_sites["n_phosphosites_changed"] > 1)).sum()),
    })

phospho_summary = pd.DataFrame(phospho_rows)
phospho_summary.set_index("tissue")

# %% [markdown]
# The left part of each bar is what a proteome alone cannot show: reactions
# whose enzyme amount stayed the same while its phosphorylation changed.

# %%
from transnet.visualization import plot_modification_axis

phospho_axis_figure = plot_modification_axis(
    phospho_summary.assign(
        steady=phospho_summary["amount steady, site moved"],
        also_moved=(phospho_summary["reactions with a changed site"]
                    - phospho_summary["amount steady, site moved"]),
    ).rename(columns={"tissue": "group"})[["group", "steady", "also_moved"]],
    title="Enzyme amount against enzyme phosphorylation")
plt.show()

# %% [markdown]
# `phospho_axis_effect` stays 0 throughout, because the site edges have no
# sign. The analysis reports that these enzymes are regulated by
# phosphorylation and which way the sites moved, but not whether that speeds
# the reaction up or slows it down. That would need site-level annotation that
# the network does not have, and guessing it would put an invented direction
# into every path through the reaction.
#
# The reactions in gastrocnemius muscle:

# %%
focus_table = regulation[TISSUE_FOCUS]
focus_phospho = focus_table[(focus_table["n_phosphosites_changed"] > 0)
                            & (focus_table["gene_axis"] == 0)]
focus_phospho.head(12)[["reaction", "name", "phospho_axis",
                        "n_phosphosites_changed", "metabolite_axis"]]

# %% [markdown]
# ## Per-pathway balance, and the metabolites doing the regulating
#
# `kegg_reaction_pathways` assigns each reaction to its rat KEGG pathways, and
# `regulation_axis_summary` counts, per pathway, the reactions each axis
# activates or inhibits in gastrocnemius muscle.

# %%
from transnet.api import kegg_reaction_pathways
from transnet.visualization import plot_regulation_axes

pathways = kegg_reaction_pathways(regulation[TISSUE_FOCUS]["reaction"], organism="rno")
balance = regulation_axis_summary(regulation[TISSUE_FOCUS], pathway_map=pathways)
balance.head(10)[["pathway", "n_reactions", "gene_activated", "gene_inhibited",
                  "metabolite_activated", "metabolite_inhibited", "n_controversial"]]

# %%
axes_figure = plot_regulation_axes(
    balance[balance["pathway"] != "unassigned"].head(15),
    title=f"{TISSUE_FOCUS} after {TIMEPOINT}: regulation by pathway")
plt.show()

# %% [markdown]
# Per tissue: how many changed metabolites are known allosteric regulators, and
# whether that is more than among all measured metabolites.

# %%
rows = []
for tissue, graph in graphs.items():
    result = regulatory_role_enrichment(metabolite_regulatory_roles(graph))
    counts = result["counts"]
    rows.append({
        "tissue": tissue,
        "changed metabolites": counts["n_differential"],
        "of those, regulators": counts["n_differential_regulators"],
        "background rate": round(counts["fraction_background_regulators"], 2),
        "q": round(float(result["enrichment"].set_index("role").loc["any", "q_value"]), 3),
    })
pd.DataFrame(rows).set_index("tissue")

# %%
from transnet.visualization import plot_metabolite_regulators

regulator_figure = plot_metabolite_regulators(
    regulatory_role_enrichment(metabolite_regulatory_roles(graphs[TISSUE_FOCUS])),
    title=f"{TISSUE_FOCUS}: changed metabolites that regulate enzymes")
plt.show()

# %% [markdown]
# ## Transcription factors
#
# ChIP-Atlas has little binding data for the rat genome, so few factors have
# enough measured targets to be tested at all. The result is reported for the
# focus tissue with that limit in mind.

# %%
from transnet.visualization import plot_tf_activity

tf_table = transcription_factor_activity(graphs[TISSUE_FOCUS], min_targets=5)
print(f"{len(tf_table)} factors testable, "
      f"{int((tf_table['q_value'] <= 0.05).sum()) if len(tf_table) else 0} implicated "
      f"at q <= 0.05")
tf_figure = plot_tf_activity(tf_table, title=f"{TISSUE_FOCUS}: transcription factors")
plt.show()

# %% [markdown]
# ## Signed paths, and whether the upper layers predict the metabolites
#
# For each tissue, path tracing and propagation predict the direction of each
# changed metabolite from the changed transcripts and proteins (see the
# *signed regulatory paths* walkthrough). Paths start at the proteome: the
# mapped phosphosites form a Signaling layer above it, but their edges have no
# sign, so they cannot predict a direction. The table counts how many
# metabolites each method scores and how many it predicts correctly.

# %%
rows, traced_paths, influences = [], {}, {}
for tissue, graph in graphs.items():
    # The phosphosites are a Signaling layer, which paths would otherwise start
    # from; the question here is what the enzyme layers predict.
    paths = trace_regulatory_paths(graph, source_layer="Proteome",
                                   target_layer="Metabolome", max_length=4,
                                   max_paths=20000)
    traced_paths[tissue] = paths
    verdicts, agree, tested = path_verdicts(paths)

    seeds = {n: float(d["log2fc"]) for n, d in graph.nodes(data=True)
             if d.get("layer") in ("Transcriptome", "Proteome")
             and d.get("regulated") and d.get("log2fc") is not None}
    influence = downstream_influence(graph, seeds, target_layer="Metabolome")
    influences[tissue] = influence
    measured = influence[influence["observed"].fillna(0) != 0] if "observed" in influence else influence
    predicted = int(measured["agrees"].astype("boolean").fillna(False).sum()) \
        if "agrees" in measured and not measured.empty else 0

    rows.append({"tissue": tissue, "paths": len(paths),
                 "metabolites tested": tested, "predicted by paths": agree,
                 "metabolites reached": len(measured), "predicted by propagation": predicted,
                 "paths vs chance": versus_chance(agree, tested) if tested else "-"})
pd.DataFrame(rows).set_index("tissue")

# %% [markdown]
# The paths in trained muscle, each coloured by the measured change of its
# molecules. When the overall rate is close to chance, the figure shows the
# wrong predictions one by one instead of hiding them in a percentage.

# %%
from transnet.visualization import plot_regulatory_paths

paths_figure = plot_regulatory_paths(traced_paths[TISSUE_FOCUS], graphs[TISSUE_FOCUS],
                                     top_n=8)
plt.show()

# %%
from transnet.visualization import plot_downstream_influence

influence_figure = plot_downstream_influence(
    influences[TISSUE_FOCUS], title=f"{TISSUE_FOCUS}: metabolites predicted from the "
                                    f"enzyme layers")
plt.show()

# %% [markdown]
# The predictions are scored once per metabolite. Scored per path, the rate
# would look overwhelming, because one well-connected metabolite is reached by
# dozens of overlapping paths and each would count as a separate success.

# %% [markdown]
# ## Is each protein change transcriptional?

# %%
concordance = {}
rows = []
for tissue, graph in graphs.items():
    result = expression_concordance(graph)
    concordance[tissue] = result
    counts = result["counts"]
    rows.append({
        "tissue": tissue,
        "concordant": counts.get("concordant", 0),
        "protein only": counts.get("protein_only", 0),
        "against transcript": counts.get("discordant", 0),
        "beyond transcript": counts.get("protein_beyond_transcript"),
        "rho": result.get("correlation"),
    })
pd.DataFrame(rows).set_index("tissue")

# %% [markdown]
# The size of the "protein only" class depends partly on statistical power:
# on how many transcripts passed the significance threshold at all. The
# column that supports a claim is **beyond transcript**: a direct test of
# whether the protein changed more than its transcript, using the standard
# errors, with no threshold. It assumes both platforms report fold changes on
# comparable scales. The isobaric labelling used for the proteome compresses
# ratios, which makes a positive result conservative.

# %% [markdown]
# ## Does ubiquitination explain the proteins that changed alone?
#
# "Protein changed, transcript did not" is the largest class in every tissue.
# The usual explanations are changes in translation rate or in protein
# degradation. MoTrPAC measures ubiquitination, the signal for degradation, in
# heart and liver, so the question can be tested: of the proteins that changed
# without their transcript, how many carry a changed ubiquitination site?

# %%
ubiquitin_rows = []
for tissue in MOTRPAC_PTM_ASSAYS["UBIQ"]:
    table = ptm_tables.get((tissue, "UBIQ"))
    genes = concordance.get(tissue, {}).get("table")
    if table is None or genes is None or genes.empty:
        continue

    ubiquitinated = {site_id_map.get(p, p)
                     for p in table.loc[table["padj"] <= 0.05, "protein"]}
    measured = {site_id_map.get(p, p) for p in table["protein"]}
    discordant = genes[genes["category"].isin(["protein_only", "discordant"])]
    concordant = genes[genes["category"] == "concordant"]

    def share(frame):
        proteins = set(frame["protein"].dropna()) & measured
        if not proteins:
            return None, 0
        return len(proteins & ubiquitinated) / len(proteins), len(proteins)

    discordant_share, n_discordant = share(discordant)
    concordant_share, n_concordant = share(concordant)
    ubiquitin_rows.append({
        "tissue": tissue,
        "changed ubiquitin sites": int((table["padj"] <= 0.05).sum()),
        "protein-only measured": n_discordant,
        "protein-only ubiquitinated": (None if discordant_share is None
                                       else round(discordant_share, 3)),
        "concordant measured": n_concordant,
        "concordant ubiquitinated": (None if concordant_share is None
                                     else round(concordant_share, 3)),
    })

pd.DataFrame(ubiquitin_rows).set_index("tissue") if ubiquitin_rows else "no ubiquitin data"

# %% [markdown]
# With 43 changed ubiquitination sites in heart and 9 in liver, this assay
# cannot explain a class of several hundred proteins. The conclusion is that
# ubiquitination sites do not account for the protein-only class at eight
# weeks, not that degradation plays no part.

# %% [markdown]
# ## Acetylation in heart and liver
#
# The activity of many mitochondrial enzymes is regulated by acetylation,
# which MoTrPAC measured in heart and liver.

# %%
acetyl_rows = []
for tissue in MOTRPAC_PTM_ASSAYS["ACETYL"]:
    table = ptm_tables.get((tissue, "ACETYL"))
    if table is None:
        continue
    changed = table[table["padj"] <= 0.05]
    up = int((changed["log2FC"] > 0).sum())
    acetyl_rows.append({
        "tissue": tissue,
        "changed sites": len(changed),
        "up": up, "down": len(changed) - up,
        "proteins": changed["protein"].nunique(),
    })

pd.DataFrame(acetyl_rows).set_index("tissue") if acetyl_rows else "no acetyl data"

# %% [markdown]
# ## Hubs, within and across tissues

# %%
hub_tables = {}
for tissue, graph in graphs.items():
    hubs, without_binding, binding_share = hub_rankings(graph, top_percent=2)
    hub_tables[tissue] = hubs
    top = hubs.head(3)["name"].tolist()
    print(f"{tissue:<8}{', '.join(str(n)[:34] for n in top)}")

shared = pd.concat(
    [table.head(20).assign(tissue=tissue) for tissue, table in hub_tables.items()]
)
recurring = (shared.groupby("name")["tissue"].nunique().sort_values(ascending=False)
             .rename("tissues").to_frame().query("tissues > 1"))
recurring.head(10)

# %%
from transnet.visualization import plot_transomic_hubs

hub_figure = plot_transomic_hubs(hub_tables[TISSUE_FOCUS],
                                 title=f"{TISSUE_FOCUS}: molecules connecting layers")
plt.show()

# %% [markdown]
# ## Tissues compared as networks
#
# `compare_transomic_networks` compares the responsive networks of two tissues.
# The Jaccard index of their edges is the number of edges both contain divided
# by the number in either: 1 means the same edges responded, 0 means none in
# common.

# %%
jaccard = pd.DataFrame(index=tissues, columns=tissues, dtype=float)
for a in tissues:
    for b in tissues:
        if a == b:
            jaccard.loc[a, b] = 1.0
        else:
            comparison = compare_transomic_networks(
                responsive_subnetwork(graphs[a]), responsive_subnetwork(graphs[b]), a, b)
            jaccard.loc[a, b] = comparison["summary"]["edge_jaccard"]
jaccard.round(2)

# %% [markdown]
# The comparison of the two most similar tissues in detail: which edge types
# they share, and which molecules changed in opposite directions.

# %%
from transnet.visualization import plot_condition_comparison

pairs = [(a, b) for i, a in enumerate(tissues) for b in tissues[i + 1:]]
first, second = max(pairs, key=lambda pair: jaccard.loc[pair[0], pair[1]])
closest = compare_transomic_networks(responsive_subnetwork(graphs[first]),
                                     responsive_subnetwork(graphs[second]), first, second)
comparison_figure = plot_condition_comparison(closest, first, second, graph=network)
plt.show()

# %% [markdown]
# ## Timing across the training weeks
#
# The four contrasts of the same tissue (1, 2, 4 and 8 weeks against the same
# controls) give each molecule a time course, and so a half-response time.

# %%
weeks = {"1w": 1.0, "2w": 2.0, "4w": 4.0, "8w": 8.0}
trajectories = {}
for label in weeks:
    for layer, table in load_motrpac_contrast(TISSUE_FOCUS, label).items():
        keyed = table.assign(node=table["feature"].map(id_map.get(layer, {})))
        # several features can map to one node (two probes, one gene); average
        # them rather than letting the last one win
        keyed = keyed.dropna(subset=["node"]).groupby("node")["log2FC"].mean()
        trajectories.setdefault(layer, {})[weeks[label]] = keyed

timing = graphs[TISSUE_FOCUS].copy()
tables = {}
for layer, series in trajectories.items():
    frame = pd.DataFrame(series)
    frame.insert(0, "feature", frame.index)
    tables[layer] = frame.reset_index(drop=True)

parameters = assign_temporal_parameters(timing, tables, id_column="feature",
                                        time_columns=list(weeks.values()))
print(f"{parameters['t_half'].notna().sum():,} molecules have a half-response time")
structure = temporal_network_structure(timing)
print(structure["degree_vs_thalf"]["interpretation"])
structure["per_layer_thalf"]

# %%
from transnet.visualization import plot_temporal_structure

timing_figure = plot_temporal_structure(timing, structure, time_unit="weeks",
                                        title=f"{TISSUE_FOCUS}: connections against timing")
plt.show()

# %% [markdown]
# ## The wiring, per tissue
#
# Because every tissue is mapped onto the same network, the wiring patterns
# each response uses and the molecules each response depends on can be compared
# directly.

# %%
from transnet import (
    convergence_significance,
    regulatory_motifs,
    structural_vulnerability,
)

rows = []
for tissue, graph in graphs.items():
    motifs = regulatory_motifs(graph, responsive_only=True)["counts"]
    convergence = convergence_significance(graph, n_randomisations=100)
    cut = structural_vulnerability(responsive_subnetwork(graph), top_n=3)
    rows.append({
        "tissue": tissue,
        "product inhibition": motifs.get("product_inhibition", 0),
        "product activation": motifs.get("product_activation", 0),
        "feed-forward": motifs.get("feed_forward", 0),
        "convergent reactions": convergence["observed"],
        "null": round(convergence["null_mean"], 1),
        "z": round(convergence["z"], 1),
        "holds it together": ", ".join(cut["name"].head(2)) if not cut.empty else "-",
    })
pd.DataFrame(rows).set_index("tissue")

# %% [markdown]
# Read the convergence column against the chance expectation beside it. A
# tissue with many changed molecules has many convergent reactions by chance
# alone; the z-score says how far the real count exceeds that.
#
# The three figures below show gastrocnemius muscle.

# %%
from transnet.visualization import (
    plot_convergence_null,
    plot_regulatory_motifs,
    plot_structural_vulnerability,
)

motif_figure = plot_regulatory_motifs(
    regulatory_motifs(graphs[TISSUE_FOCUS], responsive_only=True),
    title=f"{TISSUE_FOCUS}: regulatory motifs")
plt.show()

# %%
convergence_figure = plot_convergence_null(
    convergence_significance(graphs[TISSUE_FOCUS], n_randomisations=200),
    title=f"{TISSUE_FOCUS}: convergence against chance")
plt.show()

# %%
vulnerability_figure = plot_structural_vulnerability(
    structural_vulnerability(responsive_subnetwork(graphs[TISSUE_FOCUS]), top_n=10),
    title=f"{TISSUE_FOCUS}: molecules that hold the response together")
plt.show()

# %% [markdown]
# ## Mitochondrial enzymes across tissues
#
# For each tissue, how many TCA-cycle and oxidative-phosphorylation enzymes were
# measured, and how many went up or down.

# %%
TCA_ENZYMES = ["Cs", "Aco2", "Idh2", "Idh3a", "Ogdh", "Sdha", "Sdhb", "Fh", "Mdh1", "Mdh2",
               "Got1", "Got2", "Pdha1", "Dlat"]
OXPHOS = ["Ndufa1", "Ndufs1", "Sdhc", "Uqcrc1", "Uqcrc2", "Cox4i1", "Cox5a", "Atp5f1a",
          "Atp5f1b"]

rows = []
for tissue, graph in graphs.items():
    symbols = {str(d.get("symbol") or d.get("name")): n for n, d in graph.nodes(data=True)
               if d.get("layer") == "Proteome"}
    for label, panel in [("TCA cycle", TCA_ENZYMES), ("oxidative phosphorylation", OXPHOS)]:
        measured = [symbols[s] for s in panel if s in symbols
                    and graph.nodes[symbols[s]].get("log2fc") is not None]
        up = sum(1 for n in measured if (graph.nodes[n].get("regulated") or 0) > 0)
        down = sum(1 for n in measured if (graph.nodes[n].get("regulated") or 0) < 0)
        rows.append({"tissue": tissue, "panel": label, "measured": len(measured),
                     "up": up, "down": down})
mitochondrial = pd.DataFrame(rows).pivot(index="tissue", columns="panel",
                                         values=["measured", "up", "down"])
mitochondrial

# %% [markdown]
# ## Factors, read through the network
#
# Sex and training are hard to separate in this cohort: every factor mixes the
# two.

# %%
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
from transnet.datasets import load_motrpac_matrices

matrices, design = load_motrpac_matrices(TISSUE_FOCUS)
print(f"{len(design)} animals with every layer; "
      + ", ".join(f"{l} {m.shape[1]:,}" for l, m in matrices.items()))
design.groupby(["timepoint", "sex"]).size().rename("animals").to_frame()

# %% [markdown]
# **Are the layers from the same animals?** MoTrPAC labels every sample by
# animal, so they should be, and the data can confirm it: within a group, an
# animal with a high transcript level of a gene should also tend to have a high
# level of that protein. This also validates the check itself, because a
# dataset paired by design has to pass it.

# %%
transcripts, proteins = matched_transcript_protein(
    matrices["Transcriptome"].T, matrices["Proteome"].T, graphs[TISSUE_FOCUS], id_maps=id_map)
groups = design["timepoint"].astype(str) + " " + design["sex"].astype(str)
pairing = sample_pairing_check(transcripts, proteins, groups, n_permutations=500)
print(f"{transcripts.shape[1]:,} genes in both layers; agreement as listed "
      f"{pairing['observed']:.3f}, shuffled {pairing['null_mean']:.3f} "
      f"(p = {pairing['p_value']:.3f}) -> "
      f"{'paired' if pairing['paired'] else 'pairing not confirmed'}")

# %%
factorisation = fit_factors(matrices, n_components=6)
association = factor_design_association(factorisation.factors_, design, ["timepoint", "sex"])
variance = factor_variance_explained(factorisation)
coherence = factor_network_coherence(
    graphs[TISSUE_FOCUS], factorisation.loadings_, id_maps=id_map, top_n=100, n_permutations=500)

summary = (association.pivot(index="factor", columns="term", values="partial_eta_squared")
           .join(variance[variance["layer"] == "all layers"].set_index("factor")["variance_accounted"])
           .join(coherence["table"].set_index("factor")[["fold_enrichment", "q_value"]]))
summary.round(3)

# %% [markdown]
# A factor that follows the design but whose top molecules are *not*
# connected on the network shows molecules that vary together without a known
# mechanism linking them. A factor model alone cannot make this distinction.

# %% [markdown]
# **Every factor across the design.** Each panel is one factor's scores by
# training week, split by sex (filled: female, hollow: male), headed by the
# design term it follows most strongly. Read with the table above: a factor
# whose sexes separate at every timepoint is a sex factor, one whose sexes
# diverge only late in training is the sex-specific response MoTrPAC reports.

# %%
from transnet.visualization import (
    plot_factor_network,
    plot_factor_overview,
    plot_factor_scores,
)

FACTOR_OUT = OUT / "factors" / TISSUE_FOCUS
FACTOR_OUT.mkdir(parents=True, exist_ok=True)

scores_figure = plot_factor_scores(
    factorisation.factors_, design, x="timepoint",
    x_order=["control", "1w", "2w", "4w", "8w"], hue="sex", association=association,
    title=f"{TISSUE_FOCUS}: factor scores across the design")
scores_figure.savefig(FACTOR_OUT / "factor_scores.png", dpi=150, bbox_inches="tight")
plt.show()

# %%
overview_figure = plot_factor_overview(association, variance, coherence["table"],
                                       title=f"{TISSUE_FOCUS}: what each factor follows, "
                                             f"explains and connects")
overview_figure.savefig(FACTOR_OUT / "factor_overview.png", dpi=150, bbox_inches="tight")
plt.show()

# %% [markdown]
# **Every factor on the network.** For each factor, the top molecules that are
# connected on the network, drawn in layers with the reactions between them.
# This turns a list of loadings into biochemistry: which enzymes, which of
# their transcripts, and which metabolites of the reactions they catalyse.

# %%
for factor in factorisation.factors_.columns:
    linked = coherence["nodes"].get(factor, [])
    if not linked:
        print(f"{factor}: none of its top features are joined on the network")
        continue
    figure = plot_factor_network(
        graphs[TISSUE_FOCUS], factorisation.loadings_, factor, id_maps=id_map,
        nodes=linked, max_nodes=60,
        title=f"{TISSUE_FOCUS} {factor}: its top features, connected on the network")
    figure.savefig(FACTOR_OUT / f"network_{factor}.png", dpi=150, bbox_inches="tight")
    plt.show()

# %%
factorisation.factors_.to_csv(FACTOR_OUT / "scores.csv")
association.to_csv(FACTOR_OUT / "design_association.csv", index=False)
variance.to_csv(FACTOR_OUT / "variance_accounted.csv", index=False)
coherence["table"].to_csv(FACTOR_OUT / "network_coherence.csv", index=False)
for layer, frame in factorisation.loadings_.items():
    frame.to_csv(FACTOR_OUT / f"loadings_{layer.lower()}.csv")
print(f"factor tables and figures in {FACTOR_OUT}")

# %% [markdown]
# **Do the layers agree about each factor?** A factor fitted on all layers at
# once can still be carried by one layer alone. Projecting the animals onto a
# factor using each layer separately, and correlating those projections, says
# which.

# %%
agreement = factor_cross_layer_agreement(factorisation)
agreement.pivot_table(index="factor", columns=["layer_a", "layer_b"],
                      values="correlation").round(2)

# %% [markdown]
# **Do a factor's layers land in the same place on the network?** The
# strongest molecules of each layer are spread over the network separately. If
# the resulting profiles overlap more than those of random molecules from the
# same layers, the factor's transcripts, proteins and metabolites are related
# biochemically, not only statistically.

# %%
propagation = factor_network_propagation(graphs[TISSUE_FOCUS], factorisation.loadings_,
                                         id_maps=id_map, top_n=50, n_permutations=100)
propagation["table"].round(3)

# %% [markdown]
# **Which factors to trust.** All readings in one table: the part of the
# design a factor follows; whether the layers agree about it, and whether that
# agreement comes from the design (*between-group*) or from differences between
# animals (*within-group*, which is real here because the pairing is
# confirmed); whether its molecules are directly linked; and whether its layers
# converge on the network.

# %%
robustness = pairing_robustness(factorisation, groups, n_shuffles=200)
synthesis = (
    association.pivot(index="factor", columns="term", values="partial_eta_squared")
    .join(robustness.set_index("factor")[["agreement", "retained", "verdict"]])
    .join(coherence["table"].set_index("factor")[["fold_enrichment"]]
          .rename(columns={"fold_enrichment": "direct links vs chance"}))
    .join(propagation["table"].set_index("factor")[["overlap", "null_overlap", "q_value"]]
          .rename(columns={"q_value": "overlap q"}))
)
synthesis.round(3)


# %%
joint = synthesis["verdict"] != "single-layer"
trusted_pairing = (synthesis["verdict"] == "between-group") | pairing["paired"]
robust = synthesis[joint & trusted_pairing & (synthesis["overlap q"] <= 0.05)]
leading = (robust if not robust.empty else synthesis).sort_values(
    "overlap", ascending=False).index[0]
print(f"{leading}: where its transcripts, proteins and metabolites meet")
propagation["top_nodes"][leading].head(12)[["name", "layer", "enrichment"]]

# %%
plot_factor_network(graphs[TISSUE_FOCUS], factorisation.loadings_, leading,
                    id_maps=id_map, top_n=12)
plt.show()

# %% [markdown]
# ## Figures and exports
#
# Every figure is saved, and the responsive network of gastrocnemius muscle is
# exported in every format TransNet writes (see the *saving, exporting and
# sharing* walkthrough).

# %%
from transnet.io import to_arena3d, to_cytoscape_json, to_transomics2cytoscape, write_network
from transnet.visualization import (
    plot_expression_concordance,
    plot_similarity_heatmap,
    plot_transomic_network,
    plot_transomic_network_interactive,
    plot_values_heatmap,
    transomic_backbone,
)

FOCUS = TISSUE_FOCUS

changes = pd.DataFrame({
    tissue: {n: d.get("log2fc") for n, d in graph.nodes(data=True)
             if d.get("log2fc") is not None}
    for tissue, graph in graphs.items()
})

figures = {
    f"network_{FOCUS}": plot_transomic_network(
        responsive_subnetwork(graphs[FOCUS]),
        title=f"{FOCUS}: most-regulated reactions after {TIMEPOINT} of training"),
    "axes_by_tissue": composition_figure,
    "layer_connectivity": layer_figure,
    f"regulation_axes_{FOCUS}": axes_figure,
    f"controversial_{FOCUS}": controversial_figure,
    f"metabolite_regulators_{FOCUS}": regulator_figure,
    f"tf_activity_{FOCUS}": tf_figure,
    f"downstream_influence_{FOCUS}": influence_figure,
    f"transomic_hubs_{FOCUS}": hub_figure,
    "closest_tissues": comparison_figure,
    f"temporal_structure_{FOCUS}": timing_figure,
    f"regulatory_motifs_{FOCUS}": motif_figure,
    f"convergence_null_{FOCUS}": convergence_figure,
    f"structural_vulnerability_{FOCUS}": vulnerability_figure,
    f"concordance_{FOCUS}": plot_expression_concordance(concordance[FOCUS]),
    "modification_sites": modification_figure,
    "phospho_axis": phospho_axis_figure,
    f"paths_{FOCUS}": paths_figure,
    "edge_jaccard": plot_similarity_heatmap(
        jaccard.astype(float), title="Tissues compared by shared regulatory edges"),
    "factor_overview": plot_factor_overview(association, variance, coherence["table"]),
    "cross_tissue_changes": plot_values_heatmap(
        changes.dropna(thresh=3), graph=network, max_rows=30,
        title="The same molecules across tissues"),
}
for name, figure in figures.items():
    figure.savefig(OUT / f"{name}.png", dpi=150, bbox_inches="tight")
plt.show()

responsive = responsive_subnetwork(graphs[FOCUS])
exports = OUT / "exports"
exports.mkdir(exist_ok=True)
write_network(responsive, str(exports / "csv"))
to_cytoscape_json(responsive, exports / "network_cytoscape.json")
to_arena3d(responsive, exports / "network_arena3d.json")
to_transomics2cytoscape(responsive, zip_path=exports / "network_transomics2cytoscape.zip")
# a browser draws a few hundred nodes well; export the backbone the figures use
html_network = transomic_backbone(graphs[FOCUS], max_reactions=40)
plot_transomic_network_interactive(html_network, layout="layered",
                                   title=f"{FOCUS}, {TIMEPOINT}").write_html(
    exports / "network.html")
print(f"figures in {OUT}")
print("exports:", ", ".join(sorted(p.name for p in exports.iterdir())))
