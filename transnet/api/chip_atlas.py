"""Functions to access ChIP-Atlas database.

ChIP-Atlas (https://chip-atlas.dbcls.jp) aggregates publicly available
ChIP-seq data and provides pre-computed target-gene scores for thousands of
transcription factors across multiple genomes.

Two main use-cases
------------------
1. ``get_chip_tf_targets`` — download binding-score tables for a set of TFs
   and return a clean (TF, gene, score) table suitable for Transnet.
2. ``get_ChIP_data`` / ``get_ChIP_exps`` — legacy bulk-download helpers kept
   for backward compatibility.
"""

__all__ = [
    "get_chip_tf_targets",
    "list_chip_tfs",
    "get_ChIP_data",
    "list_ChIP_cell_type_classes",
    "list_ChIP_cell_types",
    "get_ChIP_exps",
]

import logging
from io import StringIO
from typing import List, Optional

import pandas as pd
import requests


logger = logging.getLogger(__name__)

# ChIP-Atlas serves large per-factor tables and intermittently times out.
# A dropped factor is invisible in the finished network -- it simply has no
# targets -- so transient failures are retried rather than logged and skipped.
_RETRY_STATUS = frozenset({429, 500, 502, 503, 504})


def _get_with_retry(
    url: str,
    timeout: int = 60,
    retries: int = 3,
    backoff: float = 2.0,
    label: str = "",
):
    """GET *url*, retrying timeouts, connection errors and 5xx responses.

    Returns the final :class:`requests.Response`, or raises the last exception
    if every attempt failed. A 404 is returned immediately: it means the
    resource genuinely does not exist, and retrying cannot change that.
    """
    import time as _time

    last_error = None
    for attempt in range(1, retries + 1):
        try:
            response = requests.get(url, timeout=timeout)
            if response.status_code == 404 or response.status_code not in _RETRY_STATUS:
                return response
            last_error = f"HTTP {response.status_code}"
        except (requests.Timeout, requests.ConnectionError) as exc:
            last_error = f"{type(exc).__name__}: {exc}"
        except requests.RequestException:
            raise

        if attempt < retries:
            delay = backoff ** (attempt - 1)
            logger.warning(
                f"ChIP-Atlas request failed for {label or url} "
                f"({last_error}); retry {attempt}/{retries - 1} in {delay:.0f}s"
            )
            _time.sleep(delay)

    raise requests.RequestException(
        f"ChIP-Atlas request for {label or url} failed after {retries} "
        f"attempts; last error: {last_error}"
    )

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_ANALYSIS_LIST_URL = (
    "http://dbarchive.biosciencedbc.jp/kyushu-u/metadata/analysisList.tab"
)
_EXPERIMENT_LIST_URL = (
    "http://dbarchive.biosciencedbc.jp/kyushu-u/metadata/experimentList.tab"
)
_TARGET_URL = (
    "https://chip-atlas.dbcls.jp/data"
    "/{genome}/target/{protein}.{distance}.tsv"
)

# Column names for experimentList.tab
_EXP_COLS = [
    "exp_id",        # SRX/ERX accession
    "genome",        # e.g. mm10, hg38
    "exp_class",     # antigen class: TFs and others, Histone, etc.
    "antigen",       # the factor or mark that was immunoprecipitated
    "cell_class",    # broad cell-type class
    "cell_type",     # specific cell type
    "cell_description",  # free-text description of the cell type
]


# ---------------------------------------------------------------------------
# Analysis list helper (cached)
# ---------------------------------------------------------------------------


def _fetch_analysis_list() -> pd.DataFrame:
    """Download and parse the ChIP-Atlas analysis list (all TFs × genomes)."""
    r = _get_with_retry(_ANALYSIS_LIST_URL, timeout=60, label="analysis list")
    r.raise_for_status()
    rows = []
    for line in r.text.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) >= 4:
            rows.append({"protein": parts[0], "cell_type": parts[1],
                         "strand": parts[2], "genome": parts[3]})
    return pd.DataFrame(rows)


def chip_cache_dir():
    """Directory holding the downloaded ChIP-Atlas experiment list.

    Override with the ``TRANSNET_CHIP_CACHE`` environment variable.
    """
    import os
    from pathlib import Path

    configured = os.environ.get("TRANSNET_CHIP_CACHE")
    return Path(configured) if configured else Path.home() / ".cache" / "transnet" / "chip_atlas"


def _fetch_experiment_list() -> pd.DataFrame:
    """The ChIP-Atlas experiment list, downloaded once (about 200 MB) and cached."""
    path = chip_cache_dir() / "experimentList.tab"
    if path.exists():
        text = path.read_text()
    else:
        logger.info("Downloading the ChIP-Atlas experiment list (about 200 MB, once)")
        r = _get_with_retry(_EXPERIMENT_LIST_URL, timeout=600, label="experiment list")
        r.raise_for_status()
        text = r.text
        path.parent.mkdir(parents=True, exist_ok=True)
        partial = path.with_suffix(".part")
        partial.write_text(text)
        partial.replace(path)            # never leave a half-written cache behind
    rows = []
    for line in text.strip().split("\n"):
        parts = line.split("\t")
        row = {}
        for i, col in enumerate(_EXP_COLS):
            row[col] = parts[i] if i < len(parts) else None
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Primary public function
# ---------------------------------------------------------------------------

def get_chip_tf_targets(
    proteins: List[str],
    genome: str = "mm10",
    distance: int = 5,
    score_col: str = "average",
    min_score: float = 0.0,
) -> pd.DataFrame:
    """Download ChIP-Atlas target-gene scores for a set of TF proteins.

    For each TF name in *proteins* that has data in ChIP-Atlas for the given
    *genome*, the pre-computed ``{TF}|Average`` score column is downloaded
    and stacked into a tidy DataFrame.

    Parameters
    ----------
    proteins : list of str
        TF names **exactly as they appear in ChIP-Atlas** (case-sensitive).
        Use ``list_chip_tfs()`` to see what is available for a genome.
    genome : str
        UCSC genome assembly code, e.g. ``'mm10'``, ``'hg38'``.
    distance : int
        Distance threshold (kb) for assigning peaks to genes.
        Must be one of 1, 5, or 10.
    score_col : str
        Which score column to return — ``'average'`` for the
        ``{TF}|Average`` column, ``'all'`` to keep all per-experiment
        columns.
    min_score : float
        Minimum average score to include in the result.  Setting this
        > 0 reduces the table size significantly.

    Returns
    -------
    pd.DataFrame
        Tidy table with columns ``['tf', 'gene', 'score']``
        (or ``['tf', 'gene', 'exp_id', 'score']`` when
        ``score_col='all'``).  Returns an empty DataFrame if no data
        was found.
    """
    if distance not in (1, 5, 10):
        raise ValueError(f"distance must be 1, 5, or 10; got {distance}")

    # Validate which TFs are available for this genome
    try:
        analysis = _fetch_analysis_list()
        available = set(
            analysis[analysis["genome"] == genome]["protein"].tolist()
        )
    except Exception as exc:
        logger.warning(
            f"Could not fetch ChIP-Atlas analysis list: {exc}. "
            "Proceeding without validation."
        )
        available = set(proteins)

    records = []
    missing = []

    for tf in proteins:
        if tf not in available:
            missing.append(tf)
            continue

        url = _TARGET_URL.format(genome=genome, protein=tf, distance=distance)
        try:
            r = _get_with_retry(url, timeout=60, label=tf)
            if r.status_code == 404:
                missing.append(tf)
                continue
            r.raise_for_status()
            text = r.text.strip()
            if not text:
                missing.append(tf)
                continue

            df = pd.read_csv(StringIO(text), sep="\t", index_col=0)
            # The first column after the index is '{TF}|Average'
            avg_col = [c for c in df.columns if "|Average" in c]
            if not avg_col:
                missing.append(tf)
                continue

            if score_col == "average":
                for gene, row in df[avg_col[0]].items():
                    try:
                        score = float(row)
                    except (TypeError, ValueError):
                        continue
                    if score >= min_score:
                        records.append(
                            {"tf": tf, "gene": str(gene), "score": score}
                        )
            else:
                # All per-experiment columns
                exp_cols = [
                    c for c in df.columns
                    if "|" in c and "|Average" not in c
                ]
                for gene, gene_row in df.iterrows():
                    for col in exp_cols:
                        try:
                            score = float(gene_row[col])
                        except (TypeError, ValueError):
                            continue
                        if score >= min_score:
                            exp_id = col.split("|")[0]
                            records.append({
                                "tf": tf, "gene": str(gene),
                                "exp_id": exp_id, "score": score,
                            })

            logger.debug(f"ChIP-Atlas: loaded {len(df)} target genes for {tf}")

        except requests.RequestException as exc:
            logger.warning(f"ChIP-Atlas request failed for {tf}: {exc}")
            missing.append(tf)

    if missing:
        logger.info(
            f"ChIP-Atlas: no data for {len(missing)} TF(s): "
            f"{missing[:10]}{'...' if len(missing) > 10 else ''}"
        )

    if not records:
        logger.warning("ChIP-Atlas: no target-gene data retrieved.")
        return pd.DataFrame(columns=["tf", "gene", "score"])

    result = pd.DataFrame(records)
    logger.info(
        f"ChIP-Atlas: {len(result)} TF→gene associations "
        f"for {result['tf'].nunique()} TF(s)"
    )
    return result


def list_chip_tfs(genome: str = "mm10") -> List[str]:
    """Return the list of TF names available in ChIP-Atlas for *genome*.

    Parameters
    ----------
    genome : str
        UCSC genome assembly code.

    Returns
    -------
    list of str
        Sorted list of TF names.
    """
    analysis = _fetch_analysis_list()
    return sorted(analysis[analysis["genome"] == genome]["protein"].tolist())


# ---------------------------------------------------------------------------
# Legacy helpers (kept for backward compatibility)
# ---------------------------------------------------------------------------

def get_ChIP_data(
    genome: str = "hg38",
    distance: int = 5,
    proteins: Optional[List[str]] = None,
):
    """Download ChIP-Atlas target-gene data.

    .. deprecated::
        Use :func:`get_chip_tf_targets` instead.  This wrapper now accepts
        an optional *proteins* filter so it no longer downloads all ~700 TFs
        by default (which was impractical).

    Parameters
    ----------
    genome : str
        Genome assembly code.
    distance : int
        Distance threshold in kb (1, 5, or 10).
    proteins : list of str, optional
        Restrict to these TF names.  If *None*, all available TFs for
        *genome* are used (may be slow).

    Returns
    -------
    scores_unsort : pd.DataFrame
        Wide matrix with genes as rows and experiment IDs as columns.
    gene_to_experiment : dict
        Mapping ``{TF: [experiment_id, ...]}``.
    failed : list of str
        Factors whose download failed after retries. These are *absent* from
        ``scores_unsort``, which is otherwise indistinguishable from a factor
        that genuinely has no targets, so callers should surface them.
    """
    try:
        analysis = _fetch_analysis_list()
        available = analysis[analysis["genome"] == genome]["protein"].tolist()
    except Exception as exc:
        logger.error(f"Could not fetch ChIP-Atlas analysis list: {exc}")
        return pd.DataFrame(), {}, []

    if proteins is not None:
        targets = [p for p in proteins if p in available]
        if not targets:
            logger.warning(
                "None of the requested proteins found in ChIP-Atlas."
            )
            return pd.DataFrame(), {}, []
    else:
        targets = available
        logger.warning(
            f"Downloading ChIP-Atlas data for all {len(targets)} TFs "
            f"in {genome}. Pass proteins=[...] to restrict."
        )

    all_dfs = []
    gene_to_experiment: dict = {}
    failed: List[str] = []

    for protein in targets:
        url = _TARGET_URL.format(
            genome=genome, protein=protein, distance=distance
        )
        try:
            r = _get_with_retry(url, timeout=60, label=protein)
            if r.status_code != 200:
                logger.debug(f"No data for {protein} ({r.status_code})")
                continue
            text = r.text.strip()
            if not text:
                continue
            df = pd.read_csv(StringIO(text), sep="\t", index_col=0)
            if df.empty:
                continue
            all_dfs.append(df)
            # Build gene_to_experiment mapping
            avg_col = [c for c in df.columns if "|Average" in c]
            exp_cols = [c.split("|")[0] for c in df.columns
                        if "|" in c and "|Average" not in c]
            if avg_col:
                tf_name = avg_col[0].split("|")[0]
                gene_to_experiment[tf_name] = exp_cols
        except Exception as exc:
            # Exhausted retries. Record it: a factor lost here is
            # indistinguishable from a factor with no targets once the
            # network is built.
            failed.append(protein)
            logger.error(
                f"ChIP-Atlas: giving up on {protein} after retries: {exc}"
            )

    if failed:
        logger.error(
            f"ChIP-Atlas: {len(failed)} of {len(targets)} factor(s) could not "
            f"be downloaded and are MISSING from the result: "
            f"{failed[:10]}{'...' if len(failed) > 10 else ''}"
        )

    if not all_dfs:
        return pd.DataFrame(), gene_to_experiment, failed

    # Concatenate all TF data along columns (outer join on gene index)
    merged = pd.concat(all_dfs, axis=1, join="outer")
    # Keep only the per-experiment score columns (not Average columns)
    exp_cols_all = [
        c for c in merged.columns if "|" in c and "|Average" not in c
    ]
    # Deduplicate columns
    seen = set()
    keep = []
    for c in exp_cols_all:
        key = c.split("|")[0]
        if key not in seen:
            seen.add(key)
            keep.append(c)
    scores_unsort = merged[keep].copy()
    scores_unsort.columns = [c.split("|")[0] for c in keep]
    scores_unsort.dropna(axis=0, how="all", inplace=True)

    return scores_unsort, gene_to_experiment, failed


def _get_ChIP_experiment_info(genome: str = "hg38") -> pd.DataFrame:
    """Download and parse the ChIP-Atlas experiment list for *genome*."""
    df = _fetch_experiment_list()
    return df[df["genome"] == genome].reset_index(drop=True)


def list_ChIP_cell_type_classes(genome: str = "hg38") -> list:
    """Return unique cell-type classes for *genome* from ChIP-Atlas."""
    return sorted(
        set(_get_ChIP_experiment_info(genome)["cell_class"].dropna())
    )


def list_ChIP_cell_types(
    genome: str = "hg38", cell_type_class: str = None
) -> list:
    """Return unique cell types for *genome* (optionally filtered by class)."""
    exp_info = _get_ChIP_experiment_info(genome)
    if cell_type_class is None:
        return sorted(set(exp_info["cell_type"].dropna()))
    return sorted(set(
        exp_info[
            exp_info["cell_class"] == cell_type_class
        ]["cell_type"].dropna()
    ))


def get_ChIP_exps(
    genome: str = "hg38",
    cell_type_class: str = None,
    cell_type: str = None,
) -> pd.DataFrame:
    """Return ChIP-Atlas experiments for *genome*, optionally filtered."""
    exp_info = _get_ChIP_experiment_info(genome)
    if cell_type_class is not None:
        exp_info = exp_info[exp_info["cell_class"] == cell_type_class]
    if cell_type is not None:
        exp_info = exp_info[exp_info["cell_type"] == cell_type]
    return exp_info.reset_index(drop=True)
