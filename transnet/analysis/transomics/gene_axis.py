"""The top of the hierarchy: transcription factors, genes and proteins.

* :func:`transcription_factor_activity` infers which transcription factors
  drove the genes that changed, and in which direction.
* :func:`expression_concordance` compares each gene with the protein it
  encodes, to separate transcriptional from post-transcriptional changes.
"""

import logging
from typing import Optional

import numpy as np
import pandas as pd

from transnet.analysis.transomics.mapping import is_responsive

logger = logging.getLogger(__name__)

__all__ = ["transcription_factor_activity", "expression_concordance"]


from transnet.analysis.transomics.mapping import is_measured as _measured


def _bh(p_values: pd.Series) -> pd.Series:
    """Benjamini-Hochberg q-values, NaN-safe, monotone."""
    q = pd.Series(np.nan, index=p_values.index)
    valid = p_values.dropna()
    if valid.empty:
        return q
    order = valid.sort_values()
    n = len(order)
    raw = order.to_numpy() * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(raw[::-1])[::-1].clip(max=1.0)
    q.loc[order.index] = adjusted
    return q


def transcription_factor_activity(
    graph,
    min_targets: int = 5,
    min_confidence: Optional[float] = None,
) -> pd.DataFrame:
    """Infer which transcription factors drove the genes that changed.

    For each factor, its measured target genes are tested for enrichment among
    the changed genes (one-sided Fisher's exact test, Benjamini-Hochberg across
    factors). If the changed targets moved mostly one way, that is the factor's
    inferred activity (binomial test). The factor's own measured change is
    reported beside it: a factor inferred active whose level did not change is
    usually controlled by modification, ligand binding or location.

    Parameters
    ----------
    graph : networkx.MultiDiGraph
        A mapped trans-omic network with ``transcriptional_regulation`` edges.
    min_targets : int
        Skip factors with fewer measured targets than this.
    min_confidence : float, optional
        Ignore target edges with a ChIP-Atlas binding score below this.

    Returns
    -------
    pandas.DataFrame
        One row per tested factor, most significant first: ``factor``, ``name``,
        ``n_targets``, ``n_responsive_targets``, ``n_up``, ``n_down``,
        ``odds_ratio``, ``p_value``, ``q_value``, ``direction_p_value``,
        ``inferred_activity`` (+1, -1 or 0) and ``factor_regulated``.
    """
    columns = [
        "factor", "name", "n_targets", "n_responsive_targets", "n_up", "n_down",
        "odds_ratio", "p_value", "q_value", "direction_p_value",
        "inferred_activity", "factor_regulated",
    ]

    genes = {
        n: d for n, d in graph.nodes(data=True)
        if d.get("layer") == "Transcriptome" and _measured(d)
    }
    if not genes:
        logger.warning(
            "No measured Transcriptome nodes; transcription factor activity "
            "needs gene-level measurements"
        )
        return pd.DataFrame(columns=columns)

    responsive = {n for n, d in genes.items() if is_responsive(d)}
    n_measured, n_responsive = len(genes), len(responsive)

    targets = {}
    for u, v, data in graph.edges(data=True):
        if data.get("edge_type") != "transcriptional_regulation" or v not in genes:
            continue
        if min_confidence is not None:
            confidence = data.get("confidence")
            try:
                score = float(confidence)
            except (TypeError, ValueError):
                continue        # unscored edge: cannot meet a threshold
            # A network read back from CSV carries missing scores as NaN, and
            # every comparison against NaN is False, so the plain `<` test let
            # unscored edges through any threshold.
            if not score == score or score < min_confidence:
                continue
        targets.setdefault(u, set()).add(v)

    if not targets:
        logger.warning(
            "No transcriptional_regulation edges reach a measured gene; "
            "transcription factor activity cannot be inferred"
        )
        return pd.DataFrame(columns=columns)

    from scipy.stats import binomtest, fisher_exact

    rows = []
    for factor, factor_targets in targets.items():
        if len(factor_targets) < min_targets:
            continue
        hits = factor_targets & responsive
        a = len(hits)
        b = len(factor_targets) - a
        c = n_responsive - a
        d = (n_measured - len(factor_targets)) - c
        odds, p_value = fisher_exact([[a, b], [c, max(d, 0)]], alternative="greater")

        up = sum(1 for g in hits if (genes[g].get("regulated") or 0) > 0)
        down = sum(1 for g in hits if (genes[g].get("regulated") or 0) < 0)
        signed = up + down
        direction_p = (binomtest(up, signed, 0.5).pvalue if signed else float("nan"))

        attributes = graph.nodes[factor]
        rows.append({
            "factor": factor,
            "name": attributes.get("symbol") or attributes.get("name") or factor,
            "n_targets": len(factor_targets),
            "n_responsive_targets": a,
            "n_up": up,
            "n_down": down,
            "odds_ratio": float(odds),
            "p_value": float(p_value),
            "direction_p_value": float(direction_p),
            "factor_regulated": (
                int(attributes.get("regulated") or 0) if _measured(attributes) else None
            ),
        })

    if not rows:
        logger.warning(
            f"No transcription factor has {min_targets} or more measured targets"
        )
        return pd.DataFrame(columns=columns)

    table = pd.DataFrame(rows)
    table["q_value"] = _bh(table["p_value"])
    direction_q = _bh(table["direction_p_value"])
    table["inferred_activity"] = np.where(
        (table["q_value"] <= 0.05) & (direction_q <= 0.05),
        np.sign(table["n_up"] - table["n_down"]), 0,
    ).astype(int)

    table = table.sort_values(["q_value", "p_value"]).reset_index(drop=True)
    logger.info(
        f"Transcription factor activity: {len(table)} factors tested against "
        f"{n_responsive} of {n_measured} responsive genes; "
        f"{int((table['q_value'] <= 0.05).sum())} implicated, "
        f"{int((table['inferred_activity'] != 0).sum())} with a direction"
    )
    return table[columns]


def expression_concordance(graph, alpha: float = 0.05) -> dict:
    """Compare each measured gene with the protein it encodes.

    Each gene-protein pair measured in both layers is classified as
    ``concordant`` (both changed the same way), ``protein_only``,
    ``transcript_only``, ``discordant`` (opposite directions) or ``unchanged``.

    These categories depend on each layer's statistical power. When standard
    errors are mapped, each pair is also tested directly for whether the protein
    changed more than its transcript, ``z = (protein - transcript) /
    sqrt(se_p^2 + se_t^2)`` with Benjamini-Hochberg correction; pairs where it did
    are marked ``protein_beyond_transcript``. The test assumes both platforms
    report fold changes on comparable scales.

    Parameters
    ----------
    graph : networkx.MultiDiGraph
        A mapped trans-omic network with ``translation`` edges.
    alpha : float
        False discovery rate for the direct test.

    Returns
    -------
    dict
        ``table``: one row per pair, with both fold changes, both calls, the
        category and, when tested, the difference and its q-value.
        ``counts``: pairs per category, plus ``tested_difference`` and
        ``protein_beyond_transcript``. ``correlation``: Spearman correlation of
        the gene and protein fold changes.
    """
    rows = []
    for gene, protein, data in graph.edges(data=True):
        if data.get("edge_type") != "translation":
            continue
        g, p = graph.nodes[gene], graph.nodes[protein]
        if not (_measured(g) and _measured(p)):
            continue
        gr, pr = int(g.get("regulated") or 0), int(p.get("regulated") or 0)
        if gr and pr:
            category = "concordant" if gr == pr else "discordant"
        elif pr:
            category = "protein_only"
        elif gr:
            category = "transcript_only"
        else:
            category = "unchanged"
        gene_fc, protein_fc = g.get("log2fc"), p.get("log2fc")
        gene_se, protein_se = g.get("se"), p.get("se")
        difference, difference_p = np.nan, np.nan
        if None not in (gene_fc, protein_fc, gene_se, protein_se):
            spread = float(np.hypot(protein_se, gene_se))
            if spread > 0 and np.isfinite(spread):
                from scipy.stats import norm
                difference = float(protein_fc) - float(gene_fc)
                difference_p = float(2 * norm.sf(abs(difference) / spread))
        rows.append({
            "gene": gene,
            "protein": protein,
            "name": p.get("symbol") or g.get("name") or protein,
            "gene_log2fc": gene_fc,
            "protein_log2fc": protein_fc,
            "gene_regulated": gr,
            "protein_regulated": pr,
            "category": category,
            "difference": difference,
            "difference_p": difference_p,
        })

    categories = ["concordant", "protein_only", "transcript_only", "discordant", "unchanged"]
    table = pd.DataFrame(rows, columns=[
        "gene", "protein", "name", "gene_log2fc", "protein_log2fc",
        "gene_regulated", "protein_regulated", "category", "difference", "difference_p",
    ])
    counts = {c: int((table["category"] == c).sum()) for c in categories}

    table["difference_q"] = _bh(table["difference_p"]) if not table.empty else pd.Series(dtype=float)
    tested = table["difference_p"].notna() if not table.empty else pd.Series(dtype=bool)
    table["protein_beyond_transcript"] = (
        tested
        & (table["difference_q"] <= alpha)
        & (table["protein_regulated"] != 0)
        & (np.sign(table["difference"]) == table["protein_regulated"])
    ) if not table.empty else pd.Series(dtype=bool)
    counts["tested_difference"] = int(tested.sum()) if not table.empty else 0
    counts["protein_beyond_transcript"] = int(table["protein_beyond_transcript"].sum()) if not table.empty else 0

    correlation = float("nan")
    paired = table.dropna(subset=["gene_log2fc", "protein_log2fc"])
    if len(paired) >= 3:
        from scipy.stats import spearmanr
        correlation = float(spearmanr(
            paired["gene_log2fc"].astype(float), paired["protein_log2fc"].astype(float)
        ).correlation)

    if table.empty:
        logger.warning(
            "No translation edge joins a measured gene to a measured protein; "
            "expression concordance needs both layers"
        )
    order = {c: i for i, c in enumerate(categories)}
    table = table.sort_values(
        "category", key=lambda s: s.map(order)
    ).reset_index(drop=True)
    return {"table": table, "counts": counts, "correlation": correlation}
