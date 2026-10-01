"""Datasets the studies in ``notebooks/`` run on.

Loading is not analysis, but it is where most of the work hides: MoTrPAC
distributes an unsigned ANOVA and its own phenotype table, the Uematsu panel
arrives as three headerless CSVs from a GPL-licensed repository, and the brown
adipocyte matrices are semicolon-separated with decimal commas. Keeping that
here means a notebook opens with a question rather than with parsing.

Each loader returns tidy frames with the identifiers the network uses:
Ensembl or Entrez for genes, UniProt for proteins, KEGG ``C#####`` (or
PubChem, mapped later) for metabolites.
"""

from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import logging
import warnings

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

__all__ = [
    "DATA_DIR",
    "read_master_reactions",
    "load_brown_adipocyte_omics",
    "brown_adipocyte_contrast",
    "load_motrpac_contrast",
    "load_motrpac_matrices",
    "load_motrpac_ptm",
    "motrpac_objects",
    "fetch_motrpac_raw",
    "motrpac_tissues",
    "MOTRPAC_PTM_ASSAYS",
    "MOTRPAC_TISSUES",
    "MOTRPAC_TIMEPOINTS",
    "fetch_uematsu_panel",
    "load_uematsu_panel",
    "uematsu_panel_network",
    "uematsu_contrast",
    "UEMATSU_METABOLITE_IDS",
    "UEMATSU_CONTRASTS",
    "fetch_kokaji_tables",
    "load_kokaji_table",
    "kokaji_contrast",
    "KOKAJI_TABLES",
    "KOKAJI_TIMEPOINTS",
    "KOKAJI_GENOTYPES",
]

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def read_master_reactions(path: Optional[Path] = None) -> pd.DataFrame:
    """The checked-in KEGG reaction table, with its list columns intact.

    ``data/master_reactions.csv`` exists so that a network build does not make
    ~12,000 REST calls for data that changes rarely. Its list-like columns need
    converters rather than pandas' type inference, which would give a
    semicolon-joined column a numeric dtype whenever every row holds one bare
    number, and would turn the literal ``None`` that marks variable
    stoichiometry (KEGG's ``n`` for a polymer) into NaN.
    """
    path = Path(path or DATA_DIR / "master_reactions.csv")
    if not path.exists():
        raise FileNotFoundError(f"No reaction table at {path}")

    def split(value):
        return [v for v in str(value).split(";") if v] if value else []

    def stoichiometry(value):
        parsed = []
        for item in split(value):
            try:
                parsed.append(None if item in ("None", "nan") else float(item))
            except ValueError:
                parsed.append(None)          # variable stoichiometry, e.g. 'n'
        return parsed

    frame = pd.read_csv(path, converters={
        "enzyme": split, "substrates": split, "products": split,
        "stoichiometry_substrates": stoichiometry,
        "stoichiometry_products": stoichiometry,
    })
    if "reversible" not in frame.columns:     # older tables predate it
        frame["reversible"] = True
    logger.info(f"Loaded {len(frame):,} reactions from {path.name}")
    return frame


# ---------------------------------------------------------------------------
# Brown adipocytes — Anagho-Mattanovich et al., iScience 28(9):113382, 2025
# ---------------------------------------------------------------------------

#: Immortalised murine brown adipocytes, norepinephrine-stimulated, n = 6 per
#: timepoint. Transcriptome and proteome at 0/4/24 h; the metabolome adds
#: 0.5, 1, 2 and 8 h.
BROWN_ADIPOCYTE_TIMEPOINTS = {
    "Transcriptome": ["0h", "4h", "24h"],
    "Proteome": ["0h", "4h", "24h"],
    "Metabolome": ["0h", "05h", "1h", "2h", "4h", "8h", "24h"],
}

#: Hours, for the time-course statistics.
BROWN_ADIPOCYTE_HOURS = {"0h": 0.0, "05h": 0.5, "1h": 1.0, "2h": 2.0,
                         "4h": 4.0, "8h": 8.0, "24h": 24.0}


def load_brown_adipocyte_omics(source: Optional[Path] = None) -> Dict[str, pd.DataFrame]:
    """The three matrices, features x samples, on a log2 scale.

    Transcripts are CPM: genes below one count in half the samples are
    dropped before the log, since a log-CPM of a near-absent gene is noise
    with a fold change attached. Proteins arrive imputed and already log2.
    Metabolites are intensities: features missing in more than half the
    samples are dropped and the rest filled with the feature median, which is
    what the published analysis did.

    Columns are named ``<timepoint>_<replicate>``.
    """
    source = Path(source or DATA_DIR / "sample_data" / "experimental_data")
    if not (source / "gene_data_cpm.csv").exists():
        raise FileNotFoundError(
            f"No brown adipocyte data at {source}. Expected gene_data_cpm.csv, "
            f"proteomics_imputed_log2.csv and metabolomics_062023_intensities.csv."
        )

    counts = pd.read_csv(source / "gene_data_cpm.csv", index_col=0)
    expressed = (counts > 1).sum(axis=1) >= counts.shape[1] // 2
    transcriptome = np.log2(counts.loc[expressed] + 1)

    proteome = pd.read_csv(source / "proteomics_imputed_log2.csv",
                           sep=";", index_col=0, decimal=",")

    raw = pd.read_csv(source / "metabolomics_062023_intensities.csv",
                      sep=";", index_col=0, decimal=",")
    timepoints = BROWN_ADIPOCYTE_TIMEPOINTS["Metabolome"]
    columns = [c for c in raw.columns if any(str(c).startswith(f"{t}_") for t in timepoints)]
    metabolome = raw[columns].copy()
    # Keyed by compound name: the network is keyed by KEGG compound id, and
    # names resolve against its synonyms offline (build_name_id_map), while
    # PubChem ids would need a lookup service.
    metabolome.index = raw["Name"].astype(str)
    metabolome = metabolome.loc[metabolome.isna().mean(axis=1) <= 0.5]
    metabolome = metabolome.apply(lambda row: row.fillna(row.median()), axis=1)
    metabolome = np.log2(metabolome + 1)

    logger.info(f"Brown adipocytes: {transcriptome.shape[0]:,} transcripts, "
                f"{proteome.shape[0]:,} proteins, {metabolome.shape[0]:,} metabolites")
    return {"Transcriptome": transcriptome, "Proteome": proteome,
            "Metabolome": metabolome}


def timepoint_columns(frame: pd.DataFrame, timepoint: str) -> List[str]:
    """Sample columns belonging to one timepoint."""
    return [c for c in frame.columns
            if str(c).startswith(f"{timepoint}_") or str(c).split("_")[0] == timepoint]


def brown_adipocyte_contrast(frame: pd.DataFrame, control: str = "0h",
                             treatment: str = "24h") -> pd.DataFrame:
    """One timepoint against another: ``feature, log2FC, pvalue, padj, se``.

    The matrices are already log2, so a fold change is a difference of means.
    The standard error comes back too, because comparing a protein's change
    with its transcript's needs both, without a significance threshold.
    """
    from transnet.analysis.data_integration import compute_differential_expression

    control_columns = timepoint_columns(frame, control)
    treatment_columns = timepoint_columns(frame, treatment)
    if not control_columns or not treatment_columns:
        raise ValueError(f"no samples for {treatment} vs {control}")

    tidy = frame.reset_index()
    tidy = tidy.rename(columns={tidy.columns[0]: "feature"})
    result = compute_differential_expression(
        tidy, control_cols=control_columns, treatment_cols=treatment_columns,
        id_col="feature", method="t-test", data_is_log=True,
    ).rename(columns={"log2_fold_change": "log2FC", "adj_p_value": "padj",
                      "p_value": "pvalue"})

    numeric = frame.apply(pd.to_numeric, errors="coerce")
    a, b = numeric[control_columns], numeric[treatment_columns]
    se = np.sqrt(b.var(axis=1, ddof=1) / b.notna().sum(axis=1)
                 + a.var(axis=1, ddof=1) / a.notna().sum(axis=1))
    result["se"] = result["feature"].map(se)
    return result


# ---------------------------------------------------------------------------
# MoTrPAC — endurance training, Nature 629:174-183, 2024
# ---------------------------------------------------------------------------

MOTRPAC_RAW = DATA_DIR / "raw" / "motrpac" / "rda_cache"

#: Tissue as the result files name it -> (rda file label, label inside data).
MOTRPAC_TISSUES = {
    "CORTEX": ("CORTEX", "CORTEX"), "HEART": ("HEART", "HEART"),
    "KIDNEY": ("KIDNEY", "KIDNEY"), "LIVER": ("LIVER", "LIVER"),
    "LUNG": ("LUNG", "LUNG"), "SKM_GN": ("SKMGN", "SKM-GN"),
    "WAT_SC": ("WATSC", "WAT-SC"),
}
MOTRPAC_TIMEPOINTS = ("1w", "2w", "4w", "8w")
MOTRPAC_LAYERS = {"Transcriptome": "tran", "Proteome": "prot", "Metabolome": "meta"}

#: Post-translational modification assays -> the tissues that have them. All
#: three share one feature format, ``<RefSeq protein>_<site>`` (for example
#: ``NP_001030329.2_S14s``), and the same animal columns as the proteome, so one
#: loader covers them. Each answers a different question about a protein whose
#: amount did or did not change: PHOSPHO its activity state, UBIQ whether it is
#: being degraded, ACETYL the mitochondrial enzyme regulation exercise acts on.
MOTRPAC_PTM_ASSAYS = {
    "PHOSPHO": ("CORTEX", "HEART", "KIDNEY", "LIVER", "LUNG", "SKM_GN", "WAT_SC"),
    "ACETYL": ("HEART", "LIVER"),
    "UBIQ": ("HEART", "LIVER"),
}

#: Where the MotrpacRatTraining6moData objects are published. Nothing is
#: redistributed here: the .rda files are fetched into a gitignored cache.
MOTRPAC_SOURCE = (
    "https://github.com/MoTrPAC/MotrpacRatTraining6moData/raw/main/data/{object}.rda"
)

_PHENO: Optional[pd.DataFrame] = None


def motrpac_tissues(results_dir: Optional[Path] = None) -> List[str]:
    """Tissues with data present, from the distributed result tables."""
    results_dir = Path(results_dir or DATA_DIR / "motrpac_results")
    suffix = "_tran_timecourse_anova.csv"
    found = {p.name[: -len(suffix)] for p in results_dir.glob(f"*{suffix}")}
    return sorted(found & set(MOTRPAC_TISSUES))


def _motrpac_pheno() -> pd.DataFrame:
    """Sample id -> training group and sex.

    Keyed by *both* identifiers MoTrPAC uses: per-tissue datasets name columns
    by vial label, the flat metabolomics file by animal (``pid``). Keying on
    vial label alone matched no metabolomics column at all.
    """
    global _PHENO
    if _PHENO is None:
        import pyreadr

        frame = next(iter(pyreadr.read_r(str(MOTRPAC_RAW / "PHENO.rda")).values()))
        by_vial = frame.assign(sample=frame["viallabel"].astype(str))
        by_animal = frame.drop_duplicates("pid").assign(sample=lambda d: d["pid"].astype(str))
        _PHENO = (pd.concat([by_vial, by_animal]).drop_duplicates("sample")
                    .set_index("sample")[["group", "sex"]])
    return _PHENO


def _motrpac_pheno_raw() -> pd.DataFrame:
    """MoTrPAC's phenotype table as distributed, one row per vial."""
    import pyreadr

    return next(iter(pyreadr.read_r(str(MOTRPAC_RAW / "PHENO.rda")).values()))


def motrpac_objects(tissues: Sequence[str] = tuple(MOTRPAC_TISSUES),
                    assays: Sequence[str] = ("TRNSCRPT", "PROT", "PHOSPHO"),
                    ) -> List[str]:
    """The MotrpacRatTraining6moData object names a study needs.

    ``PHENO`` and the flat metabolomics file are always included, since every
    analysis needs the design and the metabolome is not per-tissue.
    """
    names = ["PHENO", "METAB_NORM_DATA_FLAT"]
    for tissue in tissues:
        file_label, _ = MOTRPAC_TISSUES.get(tissue, (tissue, tissue))
        for assay in assays:
            if assay in MOTRPAC_PTM_ASSAYS and tissue not in MOTRPAC_PTM_ASSAYS[assay]:
                continue
            names.append(f"{assay}_{file_label}_NORM_DATA")
    return names


def fetch_motrpac_raw(objects: Sequence[str], target: Optional[Path] = None,
                      refresh: bool = False) -> Path:
    """Download MoTrPAC data objects into the local cache.

    The consortium publishes the processed data as R objects in the
    ``MotrpacRatTraining6moData`` package, which is also a plain file listing on
    GitHub, so each object is one HTTPS fetch. They are cached under
    ``data/raw/motrpac/`` and never committed.

    Parameters
    ----------
    objects : sequence of str
        Object names without the ``.rda`` suffix, as :func:`motrpac_objects`
        produces them.
    target : pathlib.Path, optional
        Cache directory; defaults to the one the loaders read.
    refresh : bool
        Re-download objects already present.

    Returns
    -------
    pathlib.Path
        The cache directory.
    """
    import urllib.error
    import urllib.request

    target = Path(target or MOTRPAC_RAW)
    target.mkdir(parents=True, exist_ok=True)
    for name in objects:
        path = target / f"{name}.rda"
        if path.exists() and not refresh:
            continue
        url = MOTRPAC_SOURCE.format(object=name)
        logger.info(f"fetching {name}.rda")
        try:
            urllib.request.urlretrieve(url, path)
        except urllib.error.HTTPError as exc:
            path.unlink(missing_ok=True)
            logger.error(f"{name}: {exc}. Check the object name against "
                         f"https://motrpac.github.io/MotrpacRatTraining6moData/")
    return target


def _motrpac_matrix(tissue: str, layer: str) -> Optional[pd.DataFrame]:
    """Log-scale normalised data for one tissue and layer, or None."""
    import pyreadr

    file_label, data_label = MOTRPAC_TISSUES.get(tissue, (tissue, tissue))
    prefix = MOTRPAC_LAYERS.get(layer, layer)
    if layer in MOTRPAC_PTM_ASSAYS:                 # PHOSPHO, ACETYL, UBIQ
        path = MOTRPAC_RAW / f"{layer}_{file_label}_NORM_DATA.rda"
    elif prefix == "tran":
        path = MOTRPAC_RAW / f"TRNSCRPT_{file_label}_NORM_DATA.rda"
    elif prefix == "prot":
        path = MOTRPAC_RAW / f"PROT_{file_label}_NORM_DATA.rda"
    else:
        path = MOTRPAC_RAW / "METAB_NORM_DATA_FLAT.rda"
    if not path.exists():
        return None
    frame = next(iter(pyreadr.read_r(str(path)).values()))
    if "tissue" in frame.columns:
        frame = frame[frame["tissue"] == data_label]
    return frame if len(frame) else None


def load_motrpac_contrast(tissue: str, timepoint: str = "8w",
                          layers: Sequence[str] = ("Transcriptome", "Proteome", "Metabolome"),
                          ) -> Dict[str, pd.DataFrame]:
    """Signed contrast per layer: ``timepoint`` against sedentary controls.

    MoTrPAC distributes a time-course ANOVA, whose F statistic has no sign, so
    every molecule reached the network as "changed, direction unknown" and
    nothing needing a direction could run. This computes the direction from
    the normalised data instead.

    Sex is a large effect in these animals, so it is adjusted for rather than
    pooled over: ``log2FC`` is the mean of the within-sex differences, and the
    Welch test runs on sex-centred values. The data are already log scale, so
    a difference of means is a log2 fold change.

    Returns ``{layer: DataFrame}`` with ``feature, log2FC, pvalue, se, padj``.
    """
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    if timepoint not in MOTRPAC_TIMEPOINTS:
        raise ValueError(f"timepoint must be one of {MOTRPAC_TIMEPOINTS}")

    pheno = _motrpac_pheno()
    tables: Dict[str, pd.DataFrame] = {}
    for layer in layers:
        frame = _motrpac_matrix(tissue, layer)
        if frame is None:
            continue

        samples = [c for c in frame.columns if str(c).isdigit() and str(c) in pheno.index]
        info = pheno.loc[samples]
        keep = info["group"].isin(["control", timepoint]).to_numpy()
        samples = [s for s, k in zip(samples, keep) if k]
        info = info[keep]
        treated = (info["group"] == timepoint).to_numpy()
        if treated.sum() < 2 or (~treated).sum() < 2:
            continue

        values = frame[samples].apply(pd.to_numeric, errors="coerce").to_numpy(float)
        centred, differences = values.copy(), []
        # A feature can be missing from one sex entirely, which makes numpy warn
        # about empty slices for that row. The result is NaN and is dropped
        # below, so the warning is noise in every notebook that loads this data.
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        for sex in info["sex"].unique():
            in_sex = (info["sex"] == sex).to_numpy()
            if (in_sex & treated).sum() and (in_sex & ~treated).sum():
                differences.append(np.nanmean(values[:, in_sex & treated], axis=1)
                                   - np.nanmean(values[:, in_sex & ~treated], axis=1))
            centred[:, in_sex] -= np.nanmean(values[:, in_sex], axis=1, keepdims=True)

        with np.errstate(invalid="ignore"):
            log2fc = np.nanmean(np.vstack(differences), axis=0)
            treated_values, control_values = centred[:, treated], centred[:, ~treated]
            se = np.sqrt(
                np.nanvar(treated_values, axis=1, ddof=1) / np.sum(~np.isnan(treated_values), axis=1)
                + np.nanvar(control_values, axis=1, ddof=1) / np.sum(~np.isnan(control_values), axis=1))
            pvalue = stats.ttest_ind(treated_values, control_values, axis=1,
                                     equal_var=False, nan_policy="omit").pvalue

        table = pd.DataFrame({
            "feature": frame["feature_ID"].astype(str).to_numpy(),
            "log2FC": log2fc, "pvalue": np.asarray(pvalue, dtype=float), "se": se,
        }).dropna(subset=["log2FC", "pvalue"])
        # A metabolite measured on two platforms appears twice; keep the
        # measurement with the strongest evidence rather than an arbitrary one.
        table = (table.sort_values("pvalue").drop_duplicates("feature", keep="first")
                      .reset_index(drop=True))
        table["padj"] = multipletests(table["pvalue"], method="fdr_bh")[1]
        tables[layer] = table

    logger.info(f"MoTrPAC {tissue} {timepoint} vs control: "
                + ", ".join(f"{l} {len(t):,}" for l, t in tables.items()))
    return tables


def load_motrpac_ptm(tissue: str, assay: str = "PHOSPHO", timepoint: str = "8w",
                     ) -> pd.DataFrame:
    """Signed contrast for one post-translational modification assay.

    The same statistic as :func:`load_motrpac_contrast`: ``timepoint`` against
    sedentary controls, as the mean of the within-sex differences, Welch test on
    sex-centred values, Benjamini-Hochberg across sites.

    The unit is a modification site, not a protein. ``feature`` is MoTrPAC's
    ``<RefSeq protein>_<site>`` identifier, so mapping onto a network needs the
    protein part resolved to the identifier the network uses, and several sites
    on one protein stay separate rows.

    Parameters
    ----------
    tissue : str
        As :data:`MOTRPAC_TISSUES` names it.
    assay : {"PHOSPHO", "ACETYL", "UBIQ"}
        Which modification. Not every tissue has every assay; see
        :data:`MOTRPAC_PTM_ASSAYS`.
    timepoint : str
        One of :data:`MOTRPAC_TIMEPOINTS`.

    Returns
    -------
    pandas.DataFrame
        ``feature, protein, site, log2FC, pvalue, se, padj``. Empty if the
        object for that tissue and assay is not in the cache; fetch it with
        :func:`fetch_motrpac_raw`.
    """
    if assay not in MOTRPAC_PTM_ASSAYS:
        raise ValueError(f"assay must be one of {sorted(MOTRPAC_PTM_ASSAYS)}")
    if tissue not in MOTRPAC_PTM_ASSAYS[assay]:
        raise ValueError(
            f"{assay} was not measured in {tissue}; it covers "
            f"{', '.join(MOTRPAC_PTM_ASSAYS[assay])}"
        )

    tables = load_motrpac_contrast(tissue, timepoint, layers=(assay,))
    table = tables.get(assay)
    if table is None or table.empty:
        logger.warning(
            f"no {assay} data for {tissue}; fetch "
            f"{assay}_{MOTRPAC_TISSUES[tissue][0]}_NORM_DATA with fetch_motrpac_raw()"
        )
        return pd.DataFrame(columns=["feature", "protein", "site", "log2FC",
                                     "pvalue", "se", "padj"])

    # "NP_001030329.2_S14s" -> protein "NP_001030329.2", site "S14s". A site
    # string can name several positions ("S14sS15s"), which stays as one row
    # because that is what was measured.
    split = table["feature"].astype(str).str.rsplit("_", n=1, expand=True)
    table = table.assign(protein=split[0], site=split[1])
    changed = int((table["padj"] <= 0.05).sum())
    logger.info(f"MoTrPAC {tissue} {assay} {timepoint} vs control: "
                f"{len(table):,} sites, {changed:,} at q <= 0.05, "
                f"{table['protein'].nunique():,} proteins")
    return table[["feature", "protein", "site", "log2FC", "pvalue", "se", "padj"]]


def load_motrpac_matrices(tissue: str, min_complete: float = 0.8
                          ) -> Tuple[Dict[str, pd.DataFrame], pd.DataFrame]:
    """Per-animal matrices and design, for a joint factor model.

    Returns ``({layer: animals x features}, design)`` restricted to animals
    measured in every layer. Columns are keyed by animal (``pid``), not by
    vial: the per-tissue files name columns by vial label while the flat
    metabolomics file names them by animal, so factorising without this
    translation finds no shared samples at all.

    A metabolite measured on several platforms appears more than once; the
    most variable measurement is kept.
    """
    pheno_raw = _motrpac_pheno_raw()
    vial_to_animal = dict(zip(pheno_raw["viallabel"].astype(str),
                              pheno_raw["pid"].astype(str)))

    matrices: Dict[str, pd.DataFrame] = {}
    for layer in MOTRPAC_LAYERS:
        frame = _motrpac_matrix(tissue, layer)
        if frame is None:
            continue
        samples = [c for c in frame.columns if str(c).isdigit()]
        values = frame[samples].apply(pd.to_numeric, errors="coerce")
        values.index = frame["feature_ID"].astype(str).to_numpy()
        values.columns = [vial_to_animal.get(str(c), str(c)) for c in samples]
        values = values.T.groupby(level=0).mean()                 # animals x features
        order = values.var().sort_values(ascending=False).index
        values = values.loc[:, order]
        values = values.loc[:, ~values.columns.duplicated()]
        matrices[layer] = values.loc[:, values.notna().mean() >= min_complete]

    if len(matrices) < 2:
        raise ValueError(f"{tissue}: fewer than two layers available")

    animals = sorted(set.intersection(*(set(m.index) for m in matrices.values())))
    matrices = {layer: m.loc[animals] for layer, m in matrices.items()}
    design = (pheno_raw.drop_duplicates("pid").assign(pid=lambda d: d["pid"].astype(str))
              .set_index("pid").loc[animals, ["group", "sex"]]
              .rename(columns={"group": "timepoint"}))
    logger.info(f"MoTrPAC {tissue}: {len(animals)} animals with every layer")
    return matrices, design


# ---------------------------------------------------------------------------
# Uematsu et al., iScience 25(2):103787, 2022 — fetched, never redistributed
# ---------------------------------------------------------------------------

UEMATSU_DIR = DATA_DIR / "published" / "uematsu2022"
UEMATSU_SOURCE = "https://raw.githubusercontent.com/usa0ri/OMELET/{branch}/data/{name}.csv"
UEMATSU_FILES = {"Transcriptome": "transcriptome", "Proteome": "proteome",
                 "Metabolome": "metabolome"}

#: Differential calls as the Kuroda laboratory defines them: 1.5-fold, q <= 0.1.
UEMATSU_LOG2FC = 0.585

#: (label, reference, test), each group as (genotype, minutes after glucose).
UEMATSU_CONTRASTS: List[Tuple[str, Tuple[str, str], Tuple[str, str]]] = [
    ("WT glucose response", ("WT", "0"), ("WT", "240")),
    ("ob/ob glucose response", ("ob", "0"), ("ob", "240")),
    ("ob/ob vs WT, fasting", ("WT", "0"), ("ob", "0")),
]

#: The panel's abbreviations, curated to KEGG compounds. A mass spectrometer
#: measures a sugar phosphate as one pool while KEGG names anomers separately,
#: so one measurement attaches to every form a reaction uses. G3P here is
#: glycerol 3-phosphate -- the panel's Gpd1 reaction.
UEMATSU_METABOLITE_IDS: Dict[str, List[str]] = {
    "G6P": ["C00092", "C00668", "C01172"], "F6P": ["C00085", "C05345"],
    "F1,6P": ["C00354", "C05378"], "DHAP": ["C00111"], "3PG": ["C00197"],
    "2PG": ["C00631"], "PEP": ["C00074"], "G3P": ["C00093"], "G1P": ["C00103"],
    "Lactate": ["C00186"], "Ala": ["C00041"], "Glu": ["C00025"],
    "Asp": ["C00049"], "Citrate": ["C00158"], "Succinate": ["C00042"],
    "Fumarate": ["C00122"], "Malate": ["C00149"], "Acetyl-CoA": ["C00024"],
    "CoA": ["C00010"], "ATP": ["C00002"], "ADP": ["C00008"], "AMP": ["C00020"],
    "NADH": ["C00004"], "NAD+": ["C00003"], "NADPH": ["C00005"],
    "NADP+": ["C00006"], "GTP": ["C00044"], "GDP": ["C00035"],
    "FAD": ["C00016"], "Leu": ["C00123"], "Phe": ["C00079"],
    "Glycogen": ["C00182"],
}


def fetch_uematsu_panel(target: Optional[Path] = None, refresh: bool = False) -> Path:
    """Download the published data into a gitignored directory.

    The authors release it with their OMELET code under GPL-3.0; TransNet is
    MIT-licensed, so the files are fetched at run time and never committed.
    """
    import urllib.request

    target = Path(target or UEMATSU_DIR)
    target.mkdir(parents=True, exist_ok=True)
    for name in UEMATSU_FILES.values():
        path = target / f"{name}.csv"
        if path.exists() and not refresh:
            continue
        for branch in ("main", "master"):
            try:
                url = UEMATSU_SOURCE.format(branch=branch, name=name)
                urllib.request.urlretrieve(url, path)
                logger.info(f"Fetched {name}.csv from OMELET ({branch})")
                break
            except Exception:                      # pragma: no cover - network
                continue
        else:
            raise RuntimeError(f"could not download {name}.csv from the OMELET repository")
    return target


def uematsu_panel_network(mouse_network: Optional[Path] = None,
                          source: Optional[Path] = None):
    """The panel's neighbourhood of the mouse network.

    Returns ``(graph, id_map)``: the published panel is 19 enzyme genes and
    32 metabolites of central carbon metabolism. This walks out from those genes -- gene to protein by
    translation, protein to reaction by catalysis, then every metabolite on or
    regulating those reactions -- so the analysis runs on the same network as
    every other study rather than on a hand-drawn pathway.
    """
    from transnet.io import read_network

    mouse_network = Path(mouse_network or DATA_DIR / "mouse" / "latest")
    if not (mouse_network / "interactions.csv").exists():
        raise FileNotFoundError(
            f"No mouse network at {mouse_network}. Build one with: "
            f"python maintenance/build_networks.py --organisms mouse --brenda"
        )

    edges = pd.read_csv(mouse_network / "interactions.csv", dtype=str, low_memory=False)
    nodes = pd.read_csv(mouse_network / "nodes.csv", dtype=str, low_memory=False)

    symbols = set(load_uematsu_panel("Transcriptome", source)[0].index) | \
        set(load_uematsu_panel("Proteome", source)[0].index)
    genes = nodes[(nodes["Layer"] == "Transcriptome") & nodes["Name"].isin(symbols)]
    translation = edges[(edges["edge_type"] == "translation")
                        & edges["source"].isin(set(genes["ID"]))]
    catalysis = edges[(edges["edge_type"] == "catalysis")
                      & edges["source"].isin(set(translation["target"]))]
    reactions = set(catalysis["target"])
    chemistry = edges[
        edges["edge_type"].isin(["substrate", "product",
                                 "allosteric_activation", "allosteric_inhibition"])
        & (edges["source"].isin(reactions) | edges["target"].isin(reactions))
    ]
    kept = pd.concat([translation, catalysis, chemistry])
    ids = set(kept["source"]) | set(kept["target"])

    scratch = DATA_DIR / "published_results" / "obese_liver" / "network"
    scratch.mkdir(parents=True, exist_ok=True)
    kept.to_csv(scratch / "edges.csv", index=False)
    nodes[nodes["ID"].isin(ids)].to_csv(scratch / "nodes.csv", index=False)

    missing = sorted(symbols - set(genes["Name"]))
    if missing:
        logger.info(f"Panel genes absent from the mouse network: {missing}")

    graph = read_network(str(scratch / "edges.csv"), nodes_file=str(scratch / "nodes.csv"))

    # The panel is keyed by gene symbol; the network by Entrez and UniProt.
    # Walk translation once to key both layers by symbol, and label the
    # proteins with it -- UniProt's stored name is a full description.
    gene_ids, protein_ids = {}, {}
    for node, data in graph.nodes(data=True):
        if data.get("layer") == "Transcriptome" and data.get("name") in symbols:
            gene_ids[data["name"]] = node
            for _, protein, edge in graph.out_edges(node, data=True):
                if edge.get("edge_type") == "translation":
                    protein_ids.setdefault(data["name"], protein)
    for symbol, protein in protein_ids.items():
        graph.nodes[protein]["symbol"] = symbol

    id_map = {"Transcriptome": gene_ids, "Proteome": protein_ids,
              "Metabolome": {name: ids[0] for name, ids in UEMATSU_METABOLITE_IDS.items()}}
    return graph, id_map


def load_uematsu_panel(layer: str, source: Optional[Path] = None
                       ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """One layer as ``(values, design)``: features x samples, and the groups."""
    source = Path(source or UEMATSU_DIR)
    raw = pd.read_csv(source / f"{UEMATSU_FILES[layer]}.csv", header=None)
    design = raw.iloc[:3, 1:].T
    design.columns = ["tissue", "genotype", "minutes"]
    design = design.reset_index(drop=True)
    values = raw.iloc[3:].set_index(0)
    values = values[values.index != "Index"].apply(pd.to_numeric, errors="coerce")
    values.columns = range(values.shape[1])
    return values, design


def uematsu_contrast(values: pd.DataFrame, design: pd.DataFrame,
                     reference: Tuple[str, str], test: Tuple[str, str]) -> pd.DataFrame:
    """Welch's t-test on log2 values, Benjamini-Hochberg within the layer.

    Transcript counts get a pseudocount of 1; relative protein and metabolite
    amounts a small one, so zeros stay finite without distorting ratios.
    """
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    pseudocount = 1.0 if np.nanmax(values.to_numpy()) > 50 else 1e-3
    logged = np.log2(values.clip(lower=0) + pseudocount)

    def columns(group):
        return design.index[(design["genotype"] == group[0])
                            & (design["minutes"] == group[1])]

    a, b = logged[columns(reference)], logged[columns(test)]
    with np.errstate(invalid="ignore"):
        pvalue = stats.ttest_ind(b, a, axis=1, equal_var=False, nan_policy="omit").pvalue
    table = pd.DataFrame({
        "feature": values.index,
        "log2FC": (b.mean(axis=1) - a.mean(axis=1)).to_numpy(),
        "pvalue": np.asarray(pvalue, dtype=float),
        "se": np.sqrt(b.var(axis=1, ddof=1) / b.notna().sum(axis=1)
                      + a.var(axis=1, ddof=1) / a.notna().sum(axis=1)).to_numpy(),
    }).dropna(subset=["log2FC", "pvalue"])
    table["padj"] = multipletests(table["pvalue"], method="fdr_bh")[1]
    return table.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Kokaji et al., Science Signaling 13(660):eaaz1236, 2020 — fetched at run
# time, never redistributed
# ---------------------------------------------------------------------------

KOKAJI_DIR = DATA_DIR / "published" / "kokaji2020"

#: The preprint of this study (bioRxiv 653758) carries the same supplementary
#: tables as the journal version, openly downloadable, so a reader without a
#: subscription can still run the analysis.
KOKAJI_SOURCE = (
    "https://www.biorxiv.org/content/biorxiv/early/2019/05/31/653758/"
    "DC{number}/embed/media-{number}.xlsx?download=true"
)

#: Supplementary table -> (bioRxiv media number, what it holds). The journal and
#: the preprint number their media files differently from their tables, which is
#: why this map exists rather than a filename pattern.
KOKAJI_TABLES = {
    "S1": (7, "liver metabolome: 162 metabolites with KEGG ids, per-timepoint "
              "fold changes, q-values and t-half, in both genotypes"),
    "S4": (10, "liver transcriptome: 14,292 genes, the same columns"),
    "S7": (13, "transcription factors inferred by the authors, with the motif "
               "enrichment behind them"),
    "S10": (2, "western-blot phosphorylation of 11 signalling proteins"),
    "S13": (5, "metabolite amounts against enzyme Km and Ki"),
}

#: Minutes after the oral glucose bolus.
KOKAJI_TIMEPOINTS = (20, 60, 120, 240)

#: Genotype as the tables spell it.
KOKAJI_GENOTYPES = ("WT", "ob/ob")


def fetch_kokaji_tables(tables: Sequence[str] = tuple(KOKAJI_TABLES),
                        target: Optional[Path] = None,
                        refresh: bool = False) -> Path:
    """Download the supplementary tables into a gitignored directory.

    The measurements belong to Kokaji *et al.*; TransNet fetches them at run
    time and never commits them, as for the Uematsu panel.

    Files that cannot be fetched are reported rather than raising, since a
    publisher can change a supplement URL at any time. :func:`load_kokaji_table`
    then says exactly which file is missing and where to get it by hand.
    """
    import urllib.error
    import urllib.request

    target = Path(target or KOKAJI_DIR)
    target.mkdir(parents=True, exist_ok=True)
    for name in tables:
        if name not in KOKAJI_TABLES:
            raise ValueError(f"unknown table {name!r}; expected one of "
                             f"{sorted(KOKAJI_TABLES)}")
        number, _ = KOKAJI_TABLES[name]
        path = target / f"Table_{name}.xlsx"
        if path.exists() and not refresh:
            continue
        logger.info(f"fetching Kokaji Table {name}")
        request = urllib.request.Request(
            KOKAJI_SOURCE.format(number=number),
            headers={"User-Agent": "Mozilla/5.0 (TransNet dataset fetch)"},
        )
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                path.write_bytes(response.read())
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as exc:
            path.unlink(missing_ok=True)
            logger.error(f"Table {name}: {exc}")
    return target


def _kokaji_path(table: str, directory: Optional[Path] = None) -> Path:
    """Where a table should be, with a precise message when it is not."""
    directory = Path(directory or KOKAJI_DIR)
    path = directory / f"Table_{table}.xlsx"
    if path.exists():
        return path
    number, what = KOKAJI_TABLES[table]
    raise FileNotFoundError(
        f"Kokaji Table {table} not found at {path}.\n"
        f"It holds: {what}.\n"
        f"Either run transnet.datasets.fetch_kokaji_tables(), or download it by "
        f"hand and save it as {path.name} in {directory}:\n"
        f"  - the journal version: Science Signaling 13(660):eaaz1236, 2020, "
        f"supplementary Table {table}\n"
        f"  - the preprint, openly available: "
        f"https://www.biorxiv.org/content/10.1101/653758v1.supplementary-material"
        f" (media-{number}.xlsx)"
    )


def load_kokaji_table(table: str, sheet: Optional[str] = None,
                      directory: Optional[Path] = None,
                      header_rows: int = 2) -> pd.DataFrame:
    """One supplementary table, with its two header rows joined.

    Every analysis sheet in this supplement has the same shape: two header rows,
    where the upper names a comparison ("Difference between 0 min and 20 min in
    WT mice") and the lower a quantity ("Fold change", "P value", "Q value").
    They are joined into single column names, so a column reads
    ``"Difference between 0 min and 20 min in WT mice | Fold change"``.

    Parameters
    ----------
    table : str
        ``"S1"``, ``"S4"``, ``"S7"``, ``"S10"`` or ``"S13"``.
    sheet : str, optional
        Which sheet. Defaults to the analysis sheet for the glucose arm, which
        is the one every contrast comes from.
    directory : pathlib.Path, optional
        Where the files are; defaults to ``data/published/kokaji2020/``.
    header_rows : int
        How many header rows the sheet has. Two for the analysis sheets; the
        kinetics sheet of Table S13 has one, and reading it as two would turn
        its first enzyme into a column name.
    """
    path = _kokaji_path(table, directory)
    book = pd.ExcelFile(path)
    if sheet is None:
        candidates = [s for s in book.sheet_names
                      if "Analysis" in s and "water" not in s.lower()]
        if not candidates:
            candidates = [s for s in book.sheet_names if s != "Readme"]
        sheet = candidates[0]

    raw = book.parse(sheet, header=None)
    if header_rows == 1:
        frame = raw.iloc[1:].copy()
        frame.columns = [str(name) for name in raw.iloc[0]]
        return frame.reset_index(drop=True)

    upper, lower = raw.iloc[0].ffill(), raw.iloc[1]
    names = [
        str(low) if pd.isna(up) or str(up) == "nan" else f"{up} | {low}"
        for up, low in zip(upper, lower)
    ]
    frame = raw.iloc[2:].copy()
    frame.columns = names
    return frame.reset_index(drop=True)


def kokaji_contrast(table: str = "S1", genotype: str = "WT",
                    timepoint: int = 240,
                    directory: Optional[Path] = None) -> pd.DataFrame:
    """One timepoint against 0 min, for one genotype, as the authors computed it.

    This study distributes its own fold changes, p-values and q-values per
    timepoint, so nothing is recomputed here: the contrast is read off the
    supplement. That makes the comparison against the paper's conclusions a
    comparison of *network readings* rather than of statistics.

    The fold changes are ratios, not logs, so they are converted to log2 for the
    network, which is what every analysis expects.

    Parameters
    ----------
    table : {"S1", "S4", "S10"}
        S1 the metabolome, S4 the transcriptome, S10 the western blots.
    genotype : {"WT", "ob/ob"}
    timepoint : int
        One of :data:`KOKAJI_TIMEPOINTS`.
    directory : pathlib.Path, optional

    Returns
    -------
    pandas.DataFrame
        ``feature, name, log2FC, pvalue, padj``, plus ``t_half`` and
        ``responded`` from the authors' own response call. ``feature`` is the
        KEGG compound id for S1, the Ensembl gene id for S4 and the protein name
        for S10.
    """
    if timepoint not in KOKAJI_TIMEPOINTS:
        raise ValueError(f"timepoint must be one of {KOKAJI_TIMEPOINTS}")
    if genotype not in KOKAJI_GENOTYPES:
        raise ValueError(f"genotype must be one of {KOKAJI_GENOTYPES}")

    frame = load_kokaji_table(table, directory=directory)
    id_column = next((c for c in frame.columns
                      if c in ("KEGG ID", "Ensembl ID", "Name")), frame.columns[0])
    name_column = "Name" if "Name" in frame.columns else id_column

    mice = "WT mice" if genotype == "WT" else "ob/ob mice"
    prefix = f"Difference between 0 min and {timepoint} min in {mice}"
    wanted = {"log2FC": f"{prefix} | Fold change", "pvalue": f"{prefix} | P value",
              "padj": f"{prefix} | Q value"}
    missing = [c for c in wanted.values() if c not in frame.columns]
    if missing:
        raise KeyError(
            f"Table {table} has no column {missing[0]!r}. Its comparisons are: "
            + "; ".join(sorted({c.split(" | ")[0] for c in frame.columns if " | " in c}))
        )

    out = pd.DataFrame({
        "feature": frame[id_column].astype(str).str.strip(),
        "name": frame[name_column].astype(str).str.strip(),
    })
    for new, old in wanted.items():
        out[new] = pd.to_numeric(frame[old], errors="coerce")
    # The supplement gives ratios; every analysis here expects log2.
    out["log2FC"] = np.log2(out["log2FC"].where(out["log2FC"] > 0))

    response = f"Response of {mice}"
    if f"{response} | T1/2" in frame.columns:
        out["t_half"] = pd.to_numeric(frame[f"{response} | T1/2"], errors="coerce")
    if f"{response} | Response" in frame.columns:
        out["responded"] = frame[f"{response} | Response"].astype(str).str.upper().eq("TRUE")

    out = out[out["feature"].ne("nan") & out["feature"].ne("")]
    changed = int((out["padj"] <= 0.1).sum())
    logger.info(f"Kokaji Table {table} {genotype} {timepoint} min vs 0: "
                f"{len(out):,} features, {changed:,} at q <= 0.1")
    return out.reset_index(drop=True)
