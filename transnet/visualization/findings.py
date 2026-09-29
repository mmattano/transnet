"""Figures for what each trans-omic analysis finds -- the part only the network shows.

Every analysis in the catalogue returns a table; these draw the result a
per-layer analysis could not have produced, so a finding reads at a glance:

========================================  ==========================================
:func:`plot_axis_composition`             which layer regulates the reactions (reaction regulation axes)
:func:`plot_controversial_reactions`      enzyme and metabolite axes pulling apart (reaction regulation axes)
:func:`plot_metabolite_regulators`        changed metabolites that act on enzymes (metabolite regulatory roles)
:func:`plot_regulatory_paths`             signed paths that explain a change (signed regulatory paths)
:func:`plot_transomic_hubs`               molecules connecting layers (trans-omic hubs)
:func:`plot_tf_activity`                  transcription factors behind the genes (transcription-factor activity)
:func:`plot_expression_concordance`       transcriptional vs post-transcriptional (transcript-protein concordance)
:func:`plot_downstream_influence`         predicted against observed metabolites
========================================  ==========================================

Colour follows :mod:`transnet.visualization.palette`: red/blue only ever mean
direction; layers and axes use their own validated hues.
"""

from typing import Dict, List, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from transnet.biology.schema import CURRENCY_METABOLITES
from transnet.visualization.palette import (
    AXIS_COLORS, BASELINE, CONCORDANCE_COLORS, DOWN, INK, LAYER_COLORS, MUTED,
    SECONDARY, SIGN_COLORS, SURFACE, UNCHANGED, UNMEASURED, UP, style_axes,
)

__all__ = [
    "reaction_contributions",
    "plot_axis_composition",
    "plot_controversial_reactions",
    "plot_metabolite_regulators",
    "plot_regulatory_paths",
    "plot_transomic_hubs",
    "plot_tf_activity",
    "plot_expression_concordance",
    "plot_downstream_influence",
    "plot_layer_changes",
    "plot_similarity_heatmap",
    "plot_temporal_structure",
    "plot_condition_comparison",
    "plot_convergence_null",
    "plot_structural_vulnerability",
    "plot_regulatory_motifs",
    "plot_modification_axis",
]

_EFFECT = {"allosteric_activation": 1, "allosteric_inhibition": -1,
           "substrate": 1, "product": -1}


def _name(graph, node, limit: int = 22) -> str:
    data = graph.nodes.get(node, {}) if graph is not None else {}
    text = data.get("symbol") or str(data.get("name") or node).split(";")[0]
    if " (" in text:
        text = text.split(" (")[0]
    return text[:limit]


def _title(axes, text, subtitle=None):
    """Bold title, with an optional subtitle on its own line beneath it."""
    if subtitle:
        axes.set_title(f"{text}\n", loc="left", fontsize=11.5, fontweight="bold",
                       color=INK, pad=6)
        axes.text(0, 1.012, subtitle, transform=axes.transAxes, fontsize=8.5,
                  color=SECONDARY, va="bottom", ha="left")
    else:
        axes.set_title(text, loc="left", fontsize=11.5, fontweight="bold", color=INK, pad=8)


def _empty(message: str, figsize=(6, 2)):
    figure, axes = plt.subplots(figsize=figsize)
    axes.axis("off")
    axes.text(0.5, 0.5, message, ha="center", va="center", color=SECONDARY, fontsize=10)
    return figure


# ------------------------------------------------------------------------------
# reaction regulation axes -- regulation axes
# ------------------------------------------------------------------------------


def plot_axis_composition(summary: pd.DataFrame, label: str = "contrast",
                          title: str = "Which layer regulates each reaction"):
    """Stacked bars: regulated reactions by the axis that carries them (reaction regulation axes).

    ``summary`` needs ``enzyme_axis_only``, ``both_axes`` and
    ``metabolite_axis_only`` per row, and optionally ``controversial``.
    """
    parts = [("enzyme_axis_only", "enzyme amount only", AXIS_COLORS["enzyme"]),
             ("both_axes", "both axes", AXIS_COLORS["both"]),
             ("metabolite_axis_only", "metabolites only", AXIS_COLORS["metabolite"])]
    labels = summary[label].astype(str).tolist()
    figure, axes = plt.subplots(figsize=(8.5, 0.62 * len(labels) + 1.9))
    style_axes(axes, grid_axis="x")
    left = np.zeros(len(labels))
    for column, name, color in parts:
        values = summary[column].to_numpy(dtype=float)
        axes.barh(labels, values, left=left, color=color, label=name,
                  edgecolor=SURFACE, linewidth=2, height=0.62)
        for i, (start, value) in enumerate(zip(left, values)):
            if value and value >= 0.06 * max(1.0, summary[[c for c, _, _ in parts]].sum(axis=1).max()):
                axes.text(start + value / 2, i, f"{int(value)}", ha="center", va="center",
                          color=INK, fontsize=8.5)
        left += values
    if "controversial" in summary.columns:
        for i, (total, n) in enumerate(zip(left, summary["controversial"])):
            if n:
                axes.text(total + 0.01 * max(left.max(), 1), i, f"{int(n)} controversial",
                          va="center", fontsize=8.5, color=SECONDARY)
    axes.invert_yaxis()
    axes.set_xlabel("regulated reactions")
    axes.set_xlim(0, max(left.max(), 1) * 1.28)
    axes.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3, frameon=False,
                fontsize=9, labelcolor=SECONDARY)
    _title(axes, title, "enzyme = transcript/protein amount; metabolites = substrates, "
                        "products, allosteric effectors")
    figure.tight_layout()
    return figure


def plot_controversial_reactions(graph, table: pd.DataFrame, max_enzymes: int = 6,
                                 title: str = "Controversial reactions: the axes pull apart"):
    """Tug-of-war per controversial enzyme (reaction regulation axes).

    Each bar is one measured molecule's push on the reaction -- right drives it
    up, left drives it down -- coloured by the axis it belongs to. A
    controversial reaction is the one where the two colours point opposite
    ways, which neither a gene list nor a metabolite list can show.
    """
    rows = table[table.get("controversial", False) == True]  # noqa: E712
    if rows.empty:
        return _empty("No controversial reactions: enzyme and metabolite axes agree")

    groups = []
    for via, group in rows.groupby(rows["gene_axis_via"].fillna("?")):
        first = group.iloc[0]
        contributions = reaction_contributions(graph, first["reaction"])
        contributions = contributions[contributions["regulated"] != 0].copy()
        contributions = contributions.drop_duplicates("molecule")
        if contributions.empty:
            continue
        groups.append((via, group, contributions))
    groups = sorted(groups, key=lambda g: -g[2]["contribution"].abs().sum())[:max_enzymes]

    height = sum(len(c) + 1.4 for _, _, c in groups)
    figure, axes = plt.subplots(figsize=(8.5, 0.38 * height + 1.6))
    style_axes(axes, grid_axis="x")
    # start below the top edge so the first heading clears the subtitle
    y, ticks, ticklabels = 0.9, [], []
    for via, group, contributions in groups:
        enzymes = [e for e in str(via).split(";") if e]
        heading = ", ".join(_name(graph, e) for e in enzymes) or "?"
        reactions = len(group)
        axes.text(0, y - 0.15, f"{heading}   ({reactions} reaction{'s' if reactions > 1 else ''}, "
                               f"e.g. {group.iloc[0]['reaction']})",
                  transform=axes.get_yaxis_transform(), fontsize=9.5, fontweight="bold",
                  color=INK, va="bottom", ha="left")
        y += 0.55
        for _, row in contributions.sort_values(["axis", "contribution"]).iterrows():
            axes.barh(y, row["contribution"], color=AXIS_COLORS[row["axis"]], height=0.62,
                      edgecolor=SURFACE, linewidth=1)
            role = {"enzyme": "protein", "transcript": "transcript"}.get(row["role"], row["role"])
            ticks.append(y)
            ticklabels.append(f"{row['name']} ({role}, {row['log2fc']:+.2f})")
            y += 1
        y += 0.85
    axes.set_yticks(ticks)
    axes.set_yticklabels(ticklabels, fontsize=8.5)
    axes.invert_yaxis()
    axes.axvline(0, color=BASELINE, linewidth=1)
    limit = max(1.0, max(abs(c["contribution"]).max() for _, _, c in groups) * 1.15)
    axes.set_xlim(-limit, limit)
    axes.set_xlabel("push on the reaction  (log2 FC, signed by effect)   <- down   |   up ->")
    handles = [plt.Rectangle((0, 0), 1, 1, color=AXIS_COLORS["enzyme"], label="enzyme axis"),
               plt.Rectangle((0, 0), 1, 1, color=AXIS_COLORS["metabolite"], label="metabolite axis")]
    axes.legend(handles=handles, loc="lower right", frameon=False, fontsize=9, labelcolor=SECONDARY)
    _title(axes, title, "inhibitors and products push against the reaction; enzymes, "
                        "substrates and activators push with it")
    figure.tight_layout()
    return figure


# ------------------------------------------------------------------------------
# metabolite regulatory roles -- metabolite regulators
# ------------------------------------------------------------------------------

def plot_metabolite_regulators(outcome: Dict, top_n: int = 15,
                               title: str = "Changed metabolites that regulate enzymes"):
    """Share of changed metabolites that are regulators, and who they are (metabolite regulatory roles).

    ``outcome`` is :func:`~transnet.regulatory_role_enrichment`'s result.
    """
    counts = outcome.get("counts") or {}
    regulators = outcome.get("regulators", pd.DataFrame())
    if not counts:
        return _empty("No measured metabolites to classify")

    figure, (left, right) = plt.subplots(
        1, 2, figsize=(10, max(3.2, 0.34 * min(top_n, max(len(regulators), 1)) + 1.8)),
        gridspec_kw={"width_ratios": [1, 2.2]})
    style_axes(left, grid_axis="y")
    diff_n = counts.get("n_differential", 0)
    diff_r = counts.get("n_differential_regulators", 0)
    back_n = counts.get("n_background", 0)
    back_r = counts.get("n_background_regulators", 0)
    shares = [diff_r / diff_n if diff_n else 0, back_r / back_n if back_n else 0]
    left.bar(["changed", "all measured"], shares, color=[LAYER_COLORS["Metabolome"], UNCHANGED],
             width=0.55, edgecolor=SURFACE, linewidth=2)
    for i, (share, (r, n)) in enumerate(zip(shares, [(diff_r, diff_n), (back_r, back_n)])):
        left.text(i, share + 0.02, f"{r}/{n}", ha="center", fontsize=9, color=INK)
    left.set_ylim(0, 1.05)
    left.set_ylabel("share that regulate an enzyme")
    enrichment = outcome.get("enrichment", pd.DataFrame())
    if not enrichment.empty and "role" in enrichment.columns:
        any_row = enrichment[enrichment["role"] == "any"]
        if not any_row.empty and pd.notna(any_row.iloc[0]["p_value"]):
            left.text(0.5, 1.0, f"enrichment p = {any_row.iloc[0]['p_value']:.2g}",
                      transform=left.transAxes, ha="center", va="bottom", fontsize=8.5,
                      color=SECONDARY)
    _title(left, "Regulators among changes")

    style_axes(right, grid_axis="x")
    if regulators.empty:
        right.axis("off")
        right.text(0.5, 0.5, "none of the changed metabolites\nhas an annotated regulatory target",
                   ha="center", va="center", color=SECONDARY, fontsize=9.5)
    else:
        shown = regulators.assign(
            magnitude=regulators["log2fc"].abs()).nlargest(top_n, "magnitude").iloc[::-1]
        names = [str(n or m).split(";")[0][:26] for n, m in zip(shown["name"], shown["metabolite"])]
        colors = [UP if v > 0 else DOWN for v in shown["log2fc"]]
        right.barh(names, shown["log2fc"], color=colors, height=0.6, edgecolor=SURFACE, linewidth=1)
        right.axvline(0, color=BASELINE, linewidth=1)
        span = max(1.0, shown["log2fc"].abs().max())
        right.set_xlim(-span * 1.9, span * 1.9)
        for i, (_, row) in enumerate(shown.iterrows()):
            text = (f"activates {int(row['n_reactions_activated'])} · "
                    f"inhibits {int(row['n_reactions_inhibited'])}")
            x = row["log2fc"]
            right.text(x + (0.05 * span if x >= 0 else -0.05 * span), i, text,
                       ha="left" if x >= 0 else "right", va="center", fontsize=8, color=SECONDARY)
        right.set_xlabel("log2 fold change")
        _title(right, "Which ones, and how many reactions they act on")
    figure.suptitle(title, x=0.01, ha="left", fontsize=12.5, fontweight="bold", color=INK)
    figure.tight_layout()
    return figure


# ------------------------------------------------------------------------------
# signed regulatory paths -- signed paths
# ------------------------------------------------------------------------------

def plot_regulatory_paths(paths: pd.DataFrame, graph, top_n: int = 6,
                          title: str = "Signed regulatory paths that explain a change"):
    """The shortest fully signed paths consistent with measurement (signed regulatory paths).

    One row per path, drawn left to right from regulator to metabolite: node
    fill is the measured change, ring the layer, edge colour the sign. The
    right-hand column compares the sign the path predicts with what was seen.
    """
    if paths is None or paths.empty:
        return _empty("No regulatory paths reach a measured metabolite")
    chosen = paths[(paths["consistent"] == True) & (paths["unsigned_steps"] == 0)]  # noqa: E712
    if chosen.empty:
        chosen = paths[paths["consistent"] == True]  # noqa: E712
    # Prefer what only a trans-omic network can show: paths through regulation
    # (allosteric, transcriptional, phosphorylation) and across more layers,
    # before the shortest enzyme -> product chains.
    chosen = chosen.assign(
        _regulatory=chosen["edge_types"].str.contains(
            "allosteric|transcriptional|phosphorylation|kinase").astype(int),
        _layers=chosen["layers"].map(lambda s: len(set(x.strip() for x in str(s).split("->")))),
    ).sort_values(["_regulatory", "_layers", "length"], ascending=[False, False, True])
    chosen = chosen.drop_duplicates("target").head(top_n)
    if chosen.empty:
        return _empty("No path is consistent with the measured direction")

    longest = int(chosen["length"].max()) + 1
    figure, axes = plt.subplots(figsize=(2.1 * longest + 3.2, 0.95 * len(chosen) + 1.4))
    axes.axis("off")
    for row_index, (_, row) in enumerate(chosen.iterrows()):
        nodes = [n.strip() for n in str(row["path"]).split("->")]
        kinds = [k.strip() for k in str(row["edge_types"]).split("->")]
        y = -row_index
        for i, node in enumerate(nodes):
            data = graph.nodes.get(node, {})
            state = int(data.get("regulated") or 0)
            measured = data.get("log2fc") is not None or data.get("measured")
            fill = SIGN_COLORS[state] if state else (UNCHANGED if measured else UNMEASURED)
            marker = "s" if data.get("layer") == "Reactions" else "o"
            axes.scatter(i * 2.1, y, s=260, c=fill, marker=marker, zorder=3,
                         edgecolors=LAYER_COLORS.get(data.get("layer"), MUTED), linewidths=2)
            axes.text(i * 2.1, y - 0.33, _name(graph, node, 16), ha="center", va="top",
                      fontsize=8, color=INK)
            if i < len(nodes) - 1:
                kind = kinds[i] if i < len(kinds) else ""
                sign = -1 if "inhibition" in kind else (1 if "activation" in kind else 0)
                color = SIGN_COLORS[sign] if sign else MUTED
                axes.annotate("", xy=(i * 2.1 + 1.82, y), xytext=(i * 2.1 + 0.28, y),
                              arrowprops=dict(arrowstyle="-|>", color=color, linewidth=1.4))
                axes.text(i * 2.1 + 1.05, y + 0.12, kind.replace("allosteric_", "").replace("_", " "),
                          ha="center", va="bottom", fontsize=7, color=SECONDARY)
        predicted = int(row["predicted"]) if "predicted" in row else int(row["sign"])
        observed = int(row["observed"])
        arrow = {1: "up", -1: "down", 0: "?"}
        axes.text(longest * 2.1 + 0.2, y,
                  f"predicts {arrow[predicted]}, observed {arrow[observed]}",
                  va="center", fontsize=8.5, color=INK)
    axes.set_xlim(-0.8, longest * 2.1 + 3.6)
    axes.set_ylim(-len(chosen) + 0.3, 0.9)
    total = len(paths)
    consistent = int((paths["consistent"] == True).sum())  # noqa: E712
    axes.set_title(f"{title}\n{consistent} of {total} traced paths agree with the measured "
                   f"direction", loc="left", fontsize=11.5, fontweight="bold", color=INK)
    figure.tight_layout()
    return figure


# ------------------------------------------------------------------------------
# trans-omic hubs -- hubs
# ------------------------------------------------------------------------------

def plot_transomic_hubs(hubs: pd.DataFrame, top_n: int = 20,
                        title: str = "Molecules that connect layers"):
    """Top molecules by cross-layer degree, coloured by their layer (trans-omic hubs).

    Ranking by *cross-layer* rather than total degree is what separates a
    trans-omic hub from an interactome hub: a ribosomal protein with hundreds
    of protein partners connects nothing across layers.
    """
    if hubs is None or hubs.empty:
        return _empty("No hubs")
    shown = hubs.sort_values("cross_layer_degree", ascending=False).head(top_n).iloc[::-1]
    figure, axes = plt.subplots(figsize=(8.5, 0.3 * len(shown) + 1.6))
    style_axes(axes, grid_axis="x")
    labels = [str(n or node).split(";")[0][:26] for n, node in zip(shown["name"], shown["node"])]
    colors = [LAYER_COLORS.get(l, MUTED) for l in shown["layer"]]
    axes.barh(labels, shown["cross_layer_degree"], color=colors, height=0.65,
              edgecolor=SURFACE, linewidth=1)
    span = shown["cross_layer_degree"].max()
    for i, (_, row) in enumerate(shown.iterrows()):
        state = int(row.get("regulated") or 0)
        mark = {1: " up", -1: " down"}.get(state, "")
        axes.text(row["cross_layer_degree"] + 0.01 * span, i,
                  f"{int(row['n_layers_touched'])} layer(s){mark}", va="center",
                  fontsize=8, color=UP if state > 0 else DOWN if state < 0 else SECONDARY)
    axes.set_xlim(0, span * 1.25)
    axes.set_xlabel("edges to other layers")
    present = [l for l in ("Signaling", "Transcriptome", "Proteome", "Reactions", "Metabolome")
               if l in set(shown["layer"])]
    handles = [plt.Rectangle((0, 0), 1, 1, color=LAYER_COLORS[l], label=l) for l in present]
    axes.legend(handles=handles, loc="lower right", frameon=False, fontsize=8.5, labelcolor=SECONDARY)
    _title(axes, title, "ranked by edges that cross into another layer, not by total degree")
    figure.tight_layout()
    return figure


# ------------------------------------------------------------------------------
# transcription-factor activity -- transcription factors
# ------------------------------------------------------------------------------

def plot_tf_activity(table: pd.DataFrame, top_n: int = 15, q_threshold: float = 0.05,
                     title: str = "Transcription factors behind the responsive genes"):
    """Responsive targets per factor, split by direction (transcription-factor activity).

    Bars to the right count targets that rose, to the left targets that fell.
    Factors whose targets are enriched among responsive genes are named in ink,
    the rest in grey; the factor's own measured change is noted beside it --
    a factor active through its targets but unchanged itself is regulated
    post-translationally.
    """
    if table is None or table.empty:
        return _empty("No transcription factor could be tested")
    shown = table.sort_values(["q_value", "p_value"]).head(top_n).iloc[::-1]
    figure, axes = plt.subplots(figsize=(9, 0.34 * len(shown) + 1.8))
    style_axes(axes, grid_axis="x")
    y = np.arange(len(shown))
    axes.barh(y, shown["n_up"], color=UP, height=0.62, edgecolor=SURFACE, linewidth=1,
              label="targets up")
    axes.barh(y, -shown["n_down"], color=DOWN, height=0.62, edgecolor=SURFACE, linewidth=1,
              label="targets down")
    axes.axvline(0, color=BASELINE, linewidth=1)
    labels = []
    for _, row in shown.iterrows():
        own = row.get("factor_regulated")
        own_text = ("" if own is None or pd.isna(own) else
                    {1: ", itself up", -1: ", itself down", 0: ", itself unchanged"}[int(own)])
        labels.append(f"{row['name']}  (q={row['q_value']:.1g}{own_text})")
    axes.set_yticks(y)
    axes.set_yticklabels(labels, fontsize=8.5)
    for tick, (_, row) in zip(axes.get_yticklabels(), shown.iterrows()):
        tick.set_color(INK if row["q_value"] <= q_threshold else MUTED)
    span = max(shown["n_up"].max(), shown["n_down"].max(), 1)
    axes.set_xlim(-span * 1.15, span * 1.15)
    axes.set_xlabel("responsive target genes   <- down   |   up ->")
    axes.legend(loc="lower right", frameon=False, fontsize=8.5, labelcolor=SECONDARY)
    _title(axes, title, f"named in ink: targets enriched among responsive genes "
                        f"(q <= {q_threshold})")
    figure.tight_layout()
    return figure


# ------------------------------------------------------------------------------
# transcript-protein concordance -- transcript vs protein
# ------------------------------------------------------------------------------

def plot_expression_concordance(result: Dict, label_top: int = 10,
                                title: str = "Is each protein change transcriptional?"):
    """Transcript against protein fold change, by concordance class (transcript-protein concordance)."""
    table = result.get("table", pd.DataFrame()) if result else pd.DataFrame()
    if table.empty:
        return _empty("No gene-protein pair with both layers measured")
    data = table.dropna(subset=["gene_log2fc", "protein_log2fc"])
    figure, axes = plt.subplots(figsize=(7.2, 6.4))
    style_axes(axes, grid_axis="both")
    names = {"concordant": "concordant", "protein_only": "protein only (post-transcriptional)",
             "transcript_only": "transcript only", "discordant": "discordant",
             "unchanged": "unchanged"}
    counts = result.get("counts", {})
    for category in ("unchanged", "transcript_only", "protein_only", "concordant", "discordant"):
        part = data[data["category"] == category]
        if part.empty:
            continue
        axes.scatter(part["gene_log2fc"], part["protein_log2fc"], s=26 if category != "unchanged" else 12,
                     c=CONCORDANCE_COLORS[category], edgecolors=SURFACE, linewidths=0.6,
                     label=f"{names[category]} ({counts.get(category, len(part))})",
                     zorder=2 if category == "unchanged" else 3)
    span = float(np.nanmax(np.abs(data[["gene_log2fc", "protein_log2fc"]].to_numpy()))) * 1.1 or 1
    axes.plot([-span, span], [-span, span], color=BASELINE, linewidth=0.8, zorder=1)
    axes.axhline(0, color=BASELINE, linewidth=0.8)
    axes.axvline(0, color=BASELINE, linewidth=0.8)
    axes.set_xlim(-span, span)
    axes.set_ylim(-span, span)
    interesting = data[data["category"].isin(["discordant", "protein_only"])]
    interesting = interesting.assign(
        gap=(interesting["protein_log2fc"] - interesting["gene_log2fc"]).abs()).nlargest(label_top, "gap")
    for _, row in interesting.iterrows():
        axes.annotate(str(row["name"])[:14], (row["gene_log2fc"], row["protein_log2fc"]),
                      xytext=(4, 3), textcoords="offset points", fontsize=7.5, color=INK)
    axes.set_xlabel("transcript log2 fold change")
    axes.set_ylabel("protein log2 fold change")
    axes.legend(loc="upper left", frameon=False, fontsize=8, labelcolor=SECONDARY)
    rho = result.get("correlation")
    _title(axes, title, f"Spearman rho = {rho:.2f} across {len(data)} gene-protein pairs"
                        if rho is not None and not pd.isna(rho) else None)
    figure.tight_layout()
    return figure


# ------------------------------------------------------------------------------
# influence
# ------------------------------------------------------------------------------

def plot_downstream_influence(table: pd.DataFrame, top_n: int = 25,
                              title: str = "Predicted from the upper layers, checked against measurement"):
    """Propagated influence on metabolites against their measured direction."""
    if table is None or table.empty:
        return _empty("No influence scores")
    measured = table[table["observed"].fillna(0) != 0]
    if measured.empty:
        return _empty("No measured metabolite to check a prediction against")
    shown = measured.assign(magnitude=measured["score"].abs()).nlargest(top_n, "magnitude")
    shown = shown.sort_values("score")
    figure, axes = plt.subplots(figsize=(8, 0.3 * len(shown) + 1.8))
    style_axes(axes, grid_axis="x")
    y = np.arange(len(shown))
    axes.hlines(y, 0, shown["score"], color=BASELINE, linewidth=1)
    for i, (_, row) in enumerate(shown.iterrows()):
        agrees = bool(row.get("agrees"))
        color = UP if row["observed"] > 0 else DOWN
        axes.scatter(row["score"], i, s=70, c=color if agrees else SURFACE,
                     edgecolors=color, linewidths=1.6, zorder=3)
    axes.set_yticks(y)
    axes.set_yticklabels([str(n)[:28] for n in shown["name"]], fontsize=8.5)
    axes.axvline(0, color=BASELINE, linewidth=1)
    axes.set_xlabel("predicted influence   <- down   |   up ->")
    agree = int(measured["agrees"].astype("boolean").fillna(False).astype(bool).sum())
    handles = [plt.Line2D([], [], marker="o", linestyle="", markersize=8, color=UP, label="observed up"),
               plt.Line2D([], [], marker="o", linestyle="", markersize=8, color=DOWN, label="observed down"),
               plt.Line2D([], [], marker="o", linestyle="", markersize=8, markerfacecolor=SURFACE,
                          markeredgecolor=SECONDARY, label="prediction disagrees")]
    axes.legend(handles=handles, loc="lower right", frameon=False, fontsize=8, labelcolor=SECONDARY)
    _title(axes, title, f"{agree} of {len(measured)} measured metabolites move the predicted way")
    figure.tight_layout()
    return figure


# ------------------------------------------------------------------------------
# baselines and comparisons
# ------------------------------------------------------------------------------

def plot_layer_changes(table: pd.DataFrame, group: str = "group",
                       title: str = "Changed molecules per layer"):
    """Up and down counts per layer and group -- what a per-layer analysis reports.

    ``table`` has ``group``, ``layer``, ``n_up`` and ``n_down`` columns. Bars
    to the right count increases, to the left decreases; one row per layer
    within each group, layers in display order.
    """
    from transnet.biology.schema import order_layers

    if table is None or table.empty:
        return _empty("Nothing changed")
    groups = list(dict.fromkeys(table[group]))
    layers = order_layers(table["layer"].unique())
    rows = [(g, l) for g in groups for l in layers
            if not table[(table[group] == g) & (table["layer"] == l)].empty]
    figure, axes = plt.subplots(figsize=(8.5, 0.3 * len(rows) + 0.35 * len(groups) + 1.6))
    style_axes(axes, grid_axis="x")
    y, ticks, labels, heads = 0.0, [], [], []
    for g in groups:
        heads.append((y, g))
        y += 0.7
        for l in layers:
            part = table[(table[group] == g) & (table["layer"] == l)]
            if part.empty:
                continue
            up, down = int(part["n_up"].sum()), int(part["n_down"].sum())
            axes.barh(y, up, color=UP, height=0.7, edgecolor=SURFACE, linewidth=1)
            axes.barh(y, -down, color=DOWN, height=0.7, edgecolor=SURFACE, linewidth=1)
            for value, sign in ((up, 1), (down, -1)):
                if value:
                    axes.text(sign * value, y, f" {value:,} " if sign > 0 else f" {value:,} ",
                              ha="left" if sign > 0 else "right", va="center", fontsize=8,
                              color=SECONDARY)
            ticks.append(y)
            labels.append(l)
            y += 0.85
        y += 0.4
    axes.set_yticks(ticks)
    axes.set_yticklabels(labels, fontsize=8.5)
    for position, name in heads:
        axes.text(0.0, position + 0.2, name, transform=axes.get_yaxis_transform(),
                  fontsize=9.5, fontweight="bold", color=INK, ha="left", va="center")
    axes.invert_yaxis()
    span = max(1, table["n_up"].max(), table["n_down"].max())
    axes.set_xlim(-span * 1.3, span * 1.3)
    axes.axvline(0, color=BASELINE, linewidth=1)
    axes.set_xlabel("molecules   <- decreased   |   increased ->")
    _title(axes, title)
    figure.tight_layout()
    return figure


def plot_similarity_heatmap(matrix: pd.DataFrame, title: str, value_label: str = "Jaccard"):
    """A symmetric similarity matrix with values printed in each cell."""
    if matrix is None or matrix.empty:
        return _empty("Nothing to compare")
    figure, axes = plt.subplots(figsize=(0.75 * len(matrix) + 2.6, 0.7 * len(matrix) + 1.8))
    values = matrix.to_numpy(dtype=float)
    image = axes.imshow(values, cmap="Purples", vmin=0, vmax=max(float(np.nanmax(values)), 1e-9))
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            axes.text(j, i, f"{values[i, j]:.2f}", ha="center", va="center", fontsize=8.5,
                      color=SURFACE if values[i, j] > 0.6 * np.nanmax(values) else INK)
    axes.set_xticks(range(len(matrix)))
    axes.set_xticklabels(matrix.columns, rotation=40, ha="right", fontsize=8.5, color=SECONDARY)
    axes.set_yticks(range(len(matrix)))
    axes.set_yticklabels(matrix.index, fontsize=8.5, color=SECONDARY)
    for side in axes.spines.values():
        side.set_visible(False)
    colorbar = figure.colorbar(image, ax=axes, fraction=0.046, pad=0.04)
    colorbar.set_label(value_label, color=SECONDARY)
    colorbar.outline.set_visible(False)
    axes.set_title(title, loc="left", fontsize=11.5, fontweight="bold", color=INK)
    figure.tight_layout()
    return figure


# ------------------------------------------------------------------------------
# response timing -- timing on the network
# ------------------------------------------------------------------------------

def plot_temporal_structure(graph, structure: Dict, time_unit: str = "",
                            title: str = "Does the wiring explain the timing?"):
    """Degree against response time, per layer, with the Spearman test (response timing).

    Morita et al. report that in healthy liver the best-connected molecules
    respond first. A direction is only stated when the correlation is
    significant, as :func:`~transnet.temporal_network_structure` decides.
    """
    from transnet.biology.schema import order_layers

    nodes = [(n, d) for n, d in graph.nodes(data=True) if d.get("t_half") is not None]
    if len(nodes) < 3:
        return _empty("Too few molecules carry a response time")
    figure, (left, right) = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [1.5, 1]})
    style_axes(left, grid_axis="both")
    layers = order_layers({d.get("layer") for _, d in nodes})
    for layer in layers:
        part = [(graph.degree(n), float(d["t_half"])) for n, d in nodes if d.get("layer") == layer]
        left.scatter([p[0] for p in part], [p[1] for p in part], s=22, alpha=0.8,
                     color=LAYER_COLORS.get(layer, MUTED), edgecolors=SURFACE, linewidths=0.5,
                     label=f"{layer} ({len(part)})")
    left.set_xscale("symlog")
    left.set_xlabel("degree on the network")
    left.set_ylabel(f"response half-time{f' ({time_unit})' if time_unit else ''}")
    left.legend(frameon=False, fontsize=8, labelcolor=SECONDARY)
    relation = structure.get("degree_vs_thalf", {})
    rho = relation.get("spearman_rho")
    subtitle = (f"Spearman rho = {rho:+.2f}, p = {relation['p_value']:.2g}, n = {relation['n']} -- "
                f"{relation.get('interpretation', '')}" if rho is not None else None)
    _title(left, title, subtitle)

    style_axes(right, grid_axis="x")
    per_layer = structure.get("per_layer_thalf", pd.DataFrame())
    if isinstance(per_layer, pd.DataFrame) and not per_layer.empty:
        ordered = per_layer.set_index("layer").reindex(order_layers(per_layer["layer"])).dropna(how="all")
        right.barh(ordered.index, ordered["median_t_half"], height=0.6, edgecolor=SURFACE, linewidth=1,
                   color=[LAYER_COLORS.get(l, MUTED) for l in ordered.index])
        for i, (layer, row) in enumerate(ordered.iterrows()):
            right.text(row["median_t_half"], i, f"  {row['median_t_half']:.1f} (n={int(row['n'])})",
                       va="center", fontsize=8.5, color=SECONDARY)
        right.invert_yaxis()
        right.set_xlim(0, ordered["median_t_half"].max() * 1.45)
    right.set_xlabel(f"median half-time{f' ({time_unit})' if time_unit else ''}")
    _title(right, "Which layer responds first")
    figure.tight_layout()
    return figure


# ------------------------------------------------------------------------------
# comparing conditions
# ------------------------------------------------------------------------------

def plot_condition_comparison(result: Dict, name1: str, name2: str, graph=None, top_n: int = 20,
                              title: Optional[str] = None):
    """Relationships kept, gained and lost, and molecules that changed direction.

    ``result`` is :func:`~transnet.compare_transomic_networks`'s output. The
    left panel is what an edge-type-aware comparison adds over a node overlap:
    *which kinds* of regulation a condition keeps. The right lists molecules
    whose direction differs between the two, which a structural comparison
    alone would call identical.
    """
    edges = result.get("edges_by_type", pd.DataFrame())
    shifts = result.get("regulation_shifts", pd.DataFrame())
    figure, (left, right) = plt.subplots(1, 2, figsize=(12, max(3.8, 0.28 * max(len(edges), min(top_n, len(shifts))) + 2)))
    style_axes(left, grid_axis="x")
    if not edges.empty:
        data = edges.iloc[::-1]
        y = np.arange(len(data))
        shared = data["n_shared"].to_numpy()
        gained = data[f"gained_in_{name2}"].to_numpy()
        lost = data[f"lost_in_{name2}"].to_numpy()
        left.barh(y, shared, color=UNCHANGED, height=0.62, label="in both", edgecolor=SURFACE, linewidth=1)
        left.barh(y, gained, left=shared, color=LAYER_COLORS["Reactions"], height=0.62,
                  label=f"only in {name2}", edgecolor=SURFACE, linewidth=1)
        left.barh(y, -lost, color=LAYER_COLORS["Metabolome"], height=0.62, label=f"only in {name1}",
                  edgecolor=SURFACE, linewidth=1)
        left.set_yticks(y)
        left.set_yticklabels([e.replace("_", " ") for e in data["edge_type"]], fontsize=8.5)
        left.axvline(0, color=BASELINE, linewidth=1)
        left.set_xlabel("responsive edges")
        left.legend(frameon=False, fontsize=8, labelcolor=SECONDARY, loc="lower right")
    summary = result.get("summary", {})
    _title(left, "Relationships kept and changed",
           f"edge Jaccard {summary.get('edge_jaccard', float('nan')):.2f}, "
           f"node Jaccard {summary.get('node_jaccard', float('nan')):.2f}")

    style_axes(right, grid_axis="x")
    if shifts.empty:
        right.axis("off")
        right.text(0.5, 0.5, "no molecule changed direction", ha="center", va="center", color=SECONDARY)
    else:
        data = shifts.copy()
        data["reversed"] = (data[name1] * data[name2]) < 0
        data["magnitude"] = (data[f"log2fc_{name2}"].fillna(0) - data[f"log2fc_{name1}"].fillna(0)).abs()
        data = data.sort_values(["reversed", "magnitude"], ascending=False).head(top_n).iloc[::-1]
        y = np.arange(len(data))
        for i, (_, row) in enumerate(data.iterrows()):
            a, b = row.get(f"log2fc_{name1}"), row.get(f"log2fc_{name2}")
            a = 0.0 if a is None or pd.isna(a) else float(a)
            b = 0.0 if b is None or pd.isna(b) else float(b)
            right.plot([a, b], [i, i], color=BASELINE, linewidth=1.2, zorder=1)
            right.scatter(a, i, s=36, color=SURFACE, edgecolors=SECONDARY, linewidths=1.2, zorder=2)
            right.scatter(b, i, s=36, color=UP if b > 0 else DOWN, zorder=3)
        labels = []
        for _, row in data.iterrows():
            label = _name(graph, row["node"]) if graph is not None else str(row.get("name") or row["node"])[:22]
            labels.append(f"{label} ({row['layer']})")
        right.set_yticks(y)
        right.set_yticklabels(labels, fontsize=8)
        right.axvline(0, color=BASELINE, linewidth=1)
        right.set_xlabel(f"log2 fold change: hollow = {name1}, filled = {name2}")
        n_reversed = int(((shifts[name1] * shifts[name2]) < 0).sum())
        _title(right, "Molecules whose regulation changed",
               f"{len(shifts)} changed call, {n_reversed} reversed direction")
    figure.suptitle(title or f"{name1} vs {name2}", x=0.01, ha="left", fontsize=12.5,
                    fontweight="bold", color=INK)
    figure.tight_layout()
    return figure


# ------------------------------------------------------------------------------
# the wiring: motifs, vulnerability, and the convergence null
# ------------------------------------------------------------------------------


def plot_convergence_null(result: Dict, title: str = "Convergence against the shuffled null"):
    """The null distribution with the observed count marked (convergence null).

    "A changed enzyme and a changed metabolite meet at this reaction" is the
    claim the rest of the catalogue rests on, and some convergence happens for
    free given how many molecules changed and how the degrees are distributed.
    This draws what the shuffled null produced against what was observed, which
    is the only honest way to read a convergence count.

    ``result`` is :func:`~transnet.convergence_significance`'s output.
    """
    null = np.asarray(result.get("null", []), dtype=float)
    observed = result.get("observed", 0)
    if null.size == 0:
        return _empty("No randomisations to compare against")

    figure, axes = plt.subplots(figsize=(6.4, 3.2))
    highest = int(max(null.max(), observed))
    lowest = int(min(null.min(), observed))
    bins = range(lowest, highest + 2) if highest > lowest else 10
    axes.hist(null, bins=bins, color=UNCHANGED, edgecolor=SURFACE, linewidth=0.6)
    axes.axvline(observed, color=UP, linewidth=2)
    axes.annotate(f"observed {observed:,}", xy=(observed, axes.get_ylim()[1]),
                  xytext=(-6, -10), textcoords="offset points",
                  ha="right", va="top", fontsize=9, color=UP, fontweight="bold")

    mean, z = result.get("null_mean"), result.get("z")
    subtitle = None
    if mean is not None:
        subtitle = f"null mean {mean:,.1f}"
        if z is not None:
            subtitle += f", z = {z:+.1f}"
        if result.get("p_value") is not None:
            subtitle += f", p = {result['p_value']:.3g}"
    _title(axes, title, subtitle)
    axes.set_xlabel("reactions with a changed enzyme and a changed metabolite")
    axes.set_ylabel("randomisations")
    style_axes(axes)
    figure.tight_layout()
    return figure


def plot_structural_vulnerability(table: pd.DataFrame, top_n: int = 10,
                                  title: str = "Molecules that hold the response together"):
    """How far the response falls apart when one molecule is removed.

    ``table`` is :func:`~transnet.structural_vulnerability`'s output. Bars are
    the share of the response disconnected by removing that molecule
    (``largest_loss``); the count of fragments it leaves is annotated, since a
    molecule can disconnect a large share as one piece or shatter it into many.
    """
    if table is None or table.empty:
        return _empty("No molecule splits this response: no articulation point")

    shown = table.head(top_n).iloc[::-1]
    labels = [str(row.get("name") or row.get("node"))[:28] for _, row in shown.iterrows()]
    colours = [LAYER_COLORS.get(row.get("layer"), MUTED) for _, row in shown.iterrows()]

    figure, axes = plt.subplots(figsize=(7, 0.42 * len(shown) + 1.6))
    positions = np.arange(len(shown))
    axes.barh(positions, shown["largest_loss"].astype(float), color=colours,
              edgecolor=SURFACE, linewidth=0.6)
    for position, (_, row) in zip(positions, shown.iterrows()):
        axes.text(float(row["largest_loss"]), position,
                  f"  {int(row['fragments'])} fragments", va="center",
                  fontsize=8.5, color=SECONDARY)

    axes.set_yticks(positions)
    axes.set_yticklabels(labels, fontsize=9)
    axes.set_xlabel("share of the response disconnected")
    _title(axes, title, "bar colour is the molecule's layer")
    style_axes(axes)
    figure.tight_layout()
    return figure


def plot_regulatory_motifs(result: Dict, title: str = "Regulatory motifs found in the wiring"):
    """Counts per motif type (regulatory motifs).

    ``result`` is :func:`~transnet.regulatory_motifs`'s output. The motifs are
    found from the network rather than assumed, so an empty result is a
    statement about the annotation and is drawn as one.
    """
    counts = (result or {}).get("counts") or {}
    counts = {name: int(value) for name, value in counts.items() if value}
    if not counts:
        return _empty("No regulatory motif survives the responsiveness filter")

    order = sorted(counts, key=counts.get)
    figure, axes = plt.subplots(figsize=(6.4, 0.5 * len(order) + 1.5))
    positions = np.arange(len(order))
    axes.barh(positions, [counts[name] for name in order],
              color=AXIS_COLORS["both"], edgecolor=SURFACE, linewidth=0.6)
    for position, name in zip(positions, order):
        axes.text(counts[name], position, f"  {counts[name]:,}", va="center",
                  fontsize=8.5, color=SECONDARY)

    axes.set_yticks(positions)
    axes.set_yticklabels([name.replace("_", " ") for name in order], fontsize=9)
    axes.set_xlabel("instances")
    _title(axes, title)
    style_axes(axes)
    figure.tight_layout()
    return figure


def plot_modification_axis(summary: pd.DataFrame, group: str = "group",
                           title: str = "Enzyme amount against enzyme modification"):
    """Reactions regulated by a modification, split by whether the amount moved too.

    The column that carries the claim is the one an abundance-only reading
    cannot produce: reactions whose enzyme held steady while its modification
    state changed. Those are called unregulated by a proteome alone.

    ``summary`` needs, per row, ``group``, ``steady`` (amount steady, site moved)
    and ``also_moved`` (both changed).
    """
    if summary is None or summary.empty:
        return _empty("No reaction carries a changed modification site")

    frame = summary.set_index(group)
    frame = frame.loc[frame[["steady", "also_moved"]].sum(axis=1).sort_values().index]
    positions = np.arange(len(frame))

    figure, axes = plt.subplots(figsize=(7, 0.46 * len(frame) + 1.8))
    axes.barh(positions, frame["steady"], color=LAYER_COLORS["Signaling"],
              edgecolor=SURFACE, linewidth=0.6, label="amount steady, site moved")
    axes.barh(positions, frame["also_moved"], left=frame["steady"],
              color=LAYER_COLORS["Proteome"], edgecolor=SURFACE, linewidth=0.6,
              label="amount moved as well")

    axes.set_yticks(positions)
    axes.set_yticklabels(frame.index, fontsize=9)
    axes.set_xlabel("reactions with a changed modification site")
    _title(axes, title,
           "the left block is what a proteome alone would call unregulated")
    axes.legend(loc="lower right", frameon=False, fontsize=8.5)
    style_axes(axes)
    figure.tight_layout()
    return figure
