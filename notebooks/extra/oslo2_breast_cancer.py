# %% [markdown]
# # Trans-omic differences between clinical outcome groups
#
# Worked on the Oslo2 breast cancer cohort of
# [Sharma *et al.*, *Oncogenesis* 13:22, 2024](https://doi.org/10.1038/s41389-024-00521-6):
# 335 untreated patients with transcriptome, proteome and metabolome from the
# same tumours, integrated there with MOFA+ into three prognostic clusters
# (MOC1 basal-enriched, MOC2 luminal B with the poorest long-term survival,
# MOC3 luminal A with the best).
#
# That analysis answers *which molecules separate the groups*. This one asks
# the next question: **through which mechanism** — is a reaction pushed by the
# amount of enzyme present, or by the metabolites acting on it, and do the two
# agree? Nothing here is specific to the MOC clusters: set `GROUPING` to any
# clinical variable and the same path runs.
#
# ---
#
# ## What you need
#
# A directory (set `DATA` below) holding four CSVs:
#
# | File | Rows | Columns | Notes |
# |---|---|---|---|
# | `transcriptome.csv` | genes | sample ids | log2 expression; first column gene symbol or Agilent probe id |
# | `proteome.csv` | antibodies | sample ids | RPPA signal; first column the antibody name |
# | `metabolome.csv` | metabolites | sample ids | HR-MAS intensities; first column the compound name |
# | `clinical.csv` | samples | variables | must include the grouping variable; `survival_months` and `event` enable the survival section |
# | `antibody_map.csv` | antibodies | `antibody,uniprot` | which protein each antibody reports; phospho antibodies may be left blank |
#
# Sample ids must agree across files. You also need the human network built
# once: `python maintenance/build_networks.py --organisms human --brenda`.

# %%
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from transnet.datasets import DATA_DIR      # the repository's data directory

DATA = Path(os.environ.get("OSLO2_DATA", DATA_DIR / "oslo2"))
OUT = Path(os.environ.get("OSLO2_OUT", DATA_DIR / "oslo2_results"))

GROUPING = "moc_cluster"        # any column of clinical.csv
REFERENCE, TEST = "MOC3", "MOC2"        # best against worst long-term survival
QVALUE = 0.05

REQUIRED = ["transcriptome.csv", "proteome.csv", "metabolome.csv", "clinical.csv",
            "antibody_map.csv"]
missing = [name for name in REQUIRED if not (DATA / name).exists()]
if missing:
    raise SystemExit(
        f"Missing from {DATA.resolve()}:\n  " + "\n  ".join(missing)
        + "\n\nSee the table above for the expected format, or set OSLO2_DATA to "
          "the directory holding these files."
    )

OUT.mkdir(parents=True, exist_ok=True)

# %%
from transnet import (
    build_name_id_map,
    compare_transomic_networks,
    expression_concordance,
    map_omics_to_network,
    metabolite_regulatory_roles,
    reaction_regulation_table,
    regulatory_role_enrichment,
    responsive_subnetwork,
    trace_regulatory_paths,
    transcription_factor_activity,
    transomic_hubs,
)
from transnet.analysis import compute_differential_expression, id_mapping, path_verdicts, versus_chance
from transnet.io import read_network

transcriptome = pd.read_csv(DATA / "transcriptome.csv", index_col=0)
proteome = pd.read_csv(DATA / "proteome.csv", index_col=0)
metabolome = pd.read_csv(DATA / "metabolome.csv", index_col=0)
clinical = pd.read_csv(DATA / "clinical.csv", index_col=0)
antibodies = pd.read_csv(DATA / "antibody_map.csv")

omics = {"Transcriptome": transcriptome, "Proteome": proteome, "Metabolome": metabolome}
complete = set(clinical.index.astype(str))
for frame in omics.values():
    complete &= set(frame.columns.astype(str))

pd.DataFrame({
    layer: {"features": frame.shape[0], "samples": frame.shape[1]}
    for layer, frame in omics.items()
} | {"clinical": {"features": clinical.shape[1], "samples": clinical.shape[0]}}).T \
    .assign(complete_cases=len(complete))

# %% [markdown]
# ## The network
#
# The organism-wide human network: KEGG reactions, UniProt proteins, STRING
# interactions, ChIP-Atlas binding and BRENDA allosteric effectors.

# %%
source = DATA_DIR / "human" / "latest"
if not (source / "interactions.csv").exists():
    raise SystemExit(
        f"No human network at {source.resolve()}. Build it once with:\n"
        f"  python maintenance/build_networks.py --organisms human --brenda"
    )

network = read_network(str(source / "interactions.csv"), nodes_file=str(source / "nodes.csv"))
print(f"{network.number_of_nodes():,} molecules, {network.number_of_edges():,} relationships")

# %% [markdown]
# ## Identifiers
#
# Three different problems. Transcripts are gene symbols or Agilent probes
# against a network keyed by Entrez. Antibodies name a protein *and sometimes
# a phosphosite*: the total-protein antibodies become the Proteome layer,
# while the phospho ones are kept aside rather than silently merged, because
# a phosphosite changing is not the same claim as the protein changing.
# Metabolites are compound names, which resolve against the network's synonyms
# offline.

# %%
symbol_to_node = {
    str(data.get("name")): node for node, data in network.nodes(data=True)
    if data.get("layer") == "Transcriptome" and data.get("name")
}
features = [str(f) for f in transcriptome.index]
gene_map = {f: symbol_to_node[f] for f in features if f in symbol_to_node}

if len(gene_map) < 0.3 * len(features):        # probe ids rather than symbols
    mapped = id_mapping(pd.DataFrame({"feature": features}), id_col="feature",
                        from_type="reporter", to_type="entrezgene", organism="human")
    gene_map = {k: str(v) for k, v in zip(mapped["feature"], mapped["entrezgene_id"])
                if pd.notna(v)}

phospho = antibodies["antibody"].str.contains(r"_p[STY]|phospho", case=False, regex=True,
                                              na=False)
total_protein = antibodies[~phospho].dropna(subset=["uniprot"])
protein_map = dict(zip(total_protein["antibody"].astype(str),
                       total_protein["uniprot"].astype(str)))
metabolite_map = build_name_id_map(network, [str(m) for m in metabolome.index])

print(f"{len(gene_map):,} of {len(features):,} transcripts, "
      f"{len(protein_map)} of {len(antibodies)} antibodies "
      f"({int(phospho.sum())} phospho held aside), "
      f"{len(metabolite_map)} of {len(metabolome)} metabolites resolved")
id_map = {"Transcriptome": gene_map, "Proteome": protein_map, "Metabolome": metabolite_map}

# %% [markdown]
# The metabolite count is the limiting number for everything that follows: a
# reaction is judged on the metabolites that were measured, and an 18-compound
# panel constrains how much of the metabolite axis is visible.

# %% [markdown]
# ## One contrast between outcome groups

# %%
def contrast(frame, reference_samples, test_samples):
    """Welch's t-test per feature, BH within the layer, with standard errors."""
    tidy = frame.reset_index()
    tidy = tidy.rename(columns={tidy.columns[0]: "feature"})
    result = compute_differential_expression(
        tidy, control_cols=list(reference_samples), treatment_cols=list(test_samples),
        id_col="feature", method="t-test", data_is_log=True,
    ).rename(columns={"log2_fold_change": "log2FC", "adj_p_value": "padj",
                      "p_value": "pvalue"})
    a = frame[list(reference_samples)].apply(pd.to_numeric, errors="coerce")
    b = frame[list(test_samples)].apply(pd.to_numeric, errors="coerce")
    se = np.sqrt(b.var(axis=1, ddof=1) / b.notna().sum(axis=1)
                 + a.var(axis=1, ddof=1) / a.notna().sum(axis=1))
    result["se"] = result["feature"].map(se)
    return result


def group_samples(frame, group):
    members = clinical.index[clinical[GROUPING].astype(str) == group].astype(str)
    return [s for s in frame.columns.astype(str) if s in set(members)]


tables = {
    layer: contrast(frame, group_samples(frame, REFERENCE), group_samples(frame, TEST))
    for layer, frame in omics.items()
}
{layer: int((table["padj"] <= QVALUE).sum()) for layer, table in tables.items()}

# %%
graph = network.copy()
report = map_omics_to_network(
    graph, tables, id_column="feature", log2fc_column="log2FC", qvalue_column="padj",
    se_column="se", id_map=id_map, qvalue_threshold=QVALUE,
)
report.per_layer

# %% [markdown]
# ## Which axis separates the groups
#
# The pathways the published analysis found separating the clusters are
# D-glutamine/D-glutamate and alanine/aspartate/glutamate metabolism — and the
# metabolite panel contains glutamate, glutamine, alanine and glycine. Those
# are exactly the reactions where the enzyme axis (GLS, GLUD1, GOT1/GOT2,
# ASNS) meets the metabolite axis, so they are worth reading reaction by
# reaction rather than as an enrichment score.

# %%
regulation = reaction_regulation_table(graph)
regulated = regulation[(regulation["gene_axis"] != 0) | (regulation["metabolite_axis"] != 0)]
print(f"{len(regulated):,} of {len(regulation):,} reactions regulated; "
      f"{int(regulated['controversial'].sum())} controversial")
regulation.to_csv(OUT / f"reaction_regulation_{TEST}_vs_{REFERENCE}.csv", index=False)

GLUTAMINE_ENZYMES = ["GLS", "GLS2", "GLUD1", "GOT1", "GOT2", "ASNS", "GPT", "GPT2",
                     "ASS1", "CAD", "PPAT"]
# gene_axis_via lists the enzymes (and genes) that carried the call
focus = regulated[regulated["gene_axis_via"].astype(str).str.contains(
    "|".join(GLUTAMINE_ENZYMES), case=False, na=False)]
if focus.empty:
    print("no reaction of the glutamine/glutamate panel is regulated in this contrast")
focus[["reaction", "name", "gene_axis", "gene_axis_via", "metabolite_axis",
       "controversial", "allosteric_regulators"]]

# %% [markdown]
# ## Is each protein change transcriptional?
#
# RPPA measures ~150 proteins against a genome-wide transcriptome, so this is
# a small but directly paired comparison: for each antibody's protein, did its
# transcript move the same way?

# %%
concordance = expression_concordance(graph)
counts = concordance["counts"]
print(f"{counts['concordant']} concordant, {counts['protein_only']} protein-only, "
      f"{counts['discordant']} against their transcript; "
      f"{counts.get('protein_beyond_transcript')} moved significantly further than "
      f"their transcript")
concordance["table"].to_csv(OUT / "concordance.csv", index=False)
concordance["table"].head(15)

# %% [markdown]
# ## Transcription factors behind the responsive genes

# %%
factors = transcription_factor_activity(graph, min_targets=10)
implicated = factors[factors["q_value"] <= 0.05]
print(f"{len(implicated)} of {len(factors)} factors implicated")
implicated.head(15)[["name", "n_responsive_targets", "n_up", "n_down", "q_value",
                     "inferred_activity", "factor_regulated"]]

# %% [markdown]
# ## Do the changed metabolites regulate anything?

# %%
roles = regulatory_role_enrichment(metabolite_regulatory_roles(graph))
roles["regulators"][["name", "log2fc", "reactions_activated", "reactions_inhibited"]]

# %% [markdown]
# ## Signed paths and cross-layer hubs

# %%
paths = trace_regulatory_paths(graph, target_layer="Metabolome", max_length=4, max_paths=20000)
verdicts, agree, tested = path_verdicts(paths)
print(f"{agree} of {tested} changed metabolites predicted in the right direction: "
      f"{versus_chance(agree, tested)}")

responsive = responsive_subnetwork(graph)
hubs = transomic_hubs(responsive, top_percent=2)
hubs.head(10)[["name", "layer", "cross_layer_degree", "n_layers_touched"]]

# %% [markdown]
# ## The same question of several clinical groupings
#
# The point of doing this per outcome rather than once: whether the kind of
# regulation that separates poor from good prognosis is the same kind that
# separates ER+ from ER−, or something else entirely.

# %%
GROUPINGS = [
    ("moc_cluster", "MOC3", "MOC2"),
    ("er_status", "positive", "negative"),
    ("grade", "1", "3"),
]

summaries, networks = [], {}
for column, reference, test in GROUPINGS:
    if column not in clinical.columns:
        continue
    available = set(clinical[column].astype(str))
    if not {reference, test} <= available:
        print(f"skipping {column}: needs {reference} and {test}, has {sorted(available)[:6]}")
        continue

    GROUPING = column
    layer_tables = {
        layer: contrast(frame, group_samples(frame, reference), group_samples(frame, test))
        for layer, frame in omics.items()
    }
    other = network.copy()
    map_omics_to_network(other, layer_tables, id_column="feature", log2fc_column="log2FC",
                         qvalue_column="padj", se_column="se", id_map=id_map,
                         qvalue_threshold=QVALUE)
    table = reaction_regulation_table(other)
    active = table[(table["gene_axis"] != 0) | (table["metabolite_axis"] != 0)]
    networks[f"{test} vs {reference}"] = other
    summaries.append({
        "contrast": f"{column}: {test} vs {reference}",
        "regulated": len(active),
        "enzyme_axis_only": int(((active["gene_axis"] != 0) & (active["metabolite_axis"] == 0)).sum()),
        "metabolite_axis_only": int(((active["gene_axis"] == 0) & (active["metabolite_axis"] != 0)).sum()),
        "both_axes": int(((active["gene_axis"] != 0) & (active["metabolite_axis"] != 0)).sum()),
        "controversial": int(active["controversial"].sum()),
    })

composition = pd.DataFrame(summaries)
composition.to_csv(OUT / "regulation_by_grouping.csv", index=False)
composition

# %%
labels = list(networks)
shared = pd.DataFrame(index=labels, columns=labels, dtype=float)
for a in labels:
    for b in labels:
        shared.loc[a, b] = 1.0 if a == b else compare_transomic_networks(
            responsive_subnetwork(networks[a]), responsive_subnetwork(networks[b]), a, b
        )["summary"]["edge_jaccard"]
shared.round(2)

# %% [markdown]
# ## Beside the published MOFA+ analysis
#
# Their factors separate the prognostic groups; the question they do not ask
# is whether a factor's strongest features are **connected** on the network.
# A factor that follows the design but scatters across the network is
# co-variation without a mechanism.

# %%
from transnet.analysis.factors import (
    factor_design_association,
    factor_network_coherence,
    factor_variance_explained,
    fit_factors,
)

samples = sorted(complete)
matrices = {layer: frame[samples].T.apply(pd.to_numeric, errors="coerce")
            for layer, frame in omics.items()}
design = clinical.loc[[s for s in samples]].astype(str)

factorisation = fit_factors(matrices, n_components=6)
association = factor_design_association(factorisation.factors_, design, [GROUPING])
variance = factor_variance_explained(factorisation)
coherence = factor_network_coherence(graph, factorisation.loadings_, id_maps=id_map,
                                     top_n=100, n_permutations=500)

(association.pivot(index="factor", columns="term", values="partial_eta_squared")
 .join(variance[variance["layer"] == "all layers"].set_index("factor")["variance_accounted"])
 .join(coherence["table"].set_index("factor")[["fold_enrichment", "q_value"]])
 .round(3))

# %% [markdown]
# For reference, the same matrices through PCA — from scikit-learn, since
# there is no reason for this package to wrap it.

# %%
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

joint = pd.concat([m.fillna(m.mean()) for m in matrices.values()], axis=1)
components = PCA(n_components=5).fit(StandardScaler().fit_transform(joint))
pd.Series(components.explained_variance_ratio_,
          index=[f"PC{i+1}" for i in range(5)], name="variance explained").to_frame()

# %% [markdown]
# ## Does the network reading carry prognostic information?
#
# A Cox model on a network-derived summary — how much of each patient's
# response runs through the metabolite axis — beside the clusters themselves.
# Needs `survival_months` and `event` (1 = death) in `clinical.csv`.

# %%
if {"survival_months", "event"} <= set(clinical.columns):
    from statsmodels.duration.hazard_regression import PHReg

    changed_metabolites = {n for n, d in graph.nodes(data=True)
                           if d.get("layer") == "Metabolome" and d.get("regulated")}
    measured = [name for name, node in metabolite_map.items() if node in changed_metabolites]

    if not measured:
        print("No metabolite changed between these groups, so there is no burden to "
              "test. Widen the contrast or raise the q-value threshold.")
    else:
        burden = matrices["Metabolome"].reindex(columns=measured).abs().mean(axis=1)
        survival = clinical.reindex(burden.index)[["survival_months", "event"]].astype(float)
        frame = survival.join(burden.rename("metabolite_burden")).dropna()
        print(f"{len(frame)} patients, {int(frame['event'].sum())} events, "
              f"{len(measured)} metabolites in the burden score")
        model = PHReg(frame["survival_months"], frame[["metabolite_burden"]],
                      status=frame["event"]).fit()
        print(model.summary())
else:
    print("No survival columns in clinical.csv; skipping the Cox model.")

# %% [markdown]
# ## Hand-off
#
# The interactive file opens in any browser with no Python installed.

# %%
from transnet.visualization import plot_transomic_network_interactive

if responsive.number_of_nodes():
    figure = plot_transomic_network_interactive(
        responsive, layout="layered", max_nodes=2000,
        title=f"Trans-omic response: {TEST} against {REFERENCE}",
    )
    figure.write_html(OUT / "network.html")
else:
    print("Nothing responded at this threshold, so there is no network to draw.")
print(f"results in {OUT.resolve()}")
sorted(p.name for p in OUT.iterdir())
