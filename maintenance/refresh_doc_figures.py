#!/usr/bin/env python3
"""Copy the figures the documentation shows out of the analysis outputs.

The study pages show real figures from real runs, but ``data/*_results/`` is
run output: regenerable, large, and not in the repository. So the handful of
figures the docs actually reference live under ``docs/source/figures/``, and
this script refreshes them after a rerun.

::

    make studies
    python maintenance/refresh_doc_figures.py

Anything missing is reported rather than silently skipped, since a missing
figure means the page it belongs to is describing a run that never happened.
"""

import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEST = ROOT / "docs" / "source" / "figures"

#: published path -> where the analysis writes it.
FIGURES = {
    # notebooks/studies/brown_adipocytes.py
    "brown_adipocytes/network.png": "data/brown_adipocyte_results/network.png",
    "brown_adipocytes/regulation_axes.png": "data/brown_adipocyte_results/regulation_axes.png",
    "brown_adipocytes/concordance.png": "data/brown_adipocyte_results/concordance.png",
    "brown_adipocytes/regulatory_paths.png":
        "data/brown_adipocyte_results/regulatory_paths.png",
    "brown_adipocytes/regulatory_motifs.png":
        "data/brown_adipocyte_results/regulatory_motifs.png",
    "brown_adipocytes/convergence_null.png":
        "data/brown_adipocyte_results/convergence_null.png",
    "brown_adipocytes/structural_vulnerability.png":
        "data/brown_adipocyte_results/structural_vulnerability.png",
    "brown_adipocytes/factor_scores.png": "data/brown_adipocyte_results/factors/factor_scores.png",
    "brown_adipocytes/factor_overview.png":
        "data/brown_adipocyte_results/factors/factor_overview.png",
    "brown_adipocytes/network_Factor3.png":
        "data/brown_adipocyte_results/factors/network_Factor3.png",
    # notebooks/studies/motrpac_rat.py
    "motrpac/network_SKM_GN.png":
        "data/motrpac_results/transomics/network_SKM_GN.png",
    "motrpac/axes_by_tissue.png": "data/motrpac_results/transomics/axes_by_tissue.png",
    "motrpac/modification_sites.png":
        "data/motrpac_results/transomics/modification_sites.png",
    "motrpac/regulatory_paths.png":
        "data/motrpac_results/transomics/paths_SKM_GN.png",
    "motrpac/phospho_axis.png": "data/motrpac_results/transomics/phospho_axis.png",
    "kokaji/saturation.png": "data/published_results/kokaji_liver/saturation.png",
    "motrpac/concordance_SKM_GN.png":
        "data/motrpac_results/transomics/concordance_SKM_GN.png",
    "motrpac/edge_jaccard.png": "data/motrpac_results/transomics/edge_jaccard.png",
    "motrpac/factor_overview.png":
        "data/motrpac_results/transomics/factors/SKM_GN/factor_overview.png",
    "motrpac/factor_scores.png":
        "data/motrpac_results/transomics/factors/SKM_GN/factor_scores.png",
    "motrpac/network_Factor5.png":
        "data/motrpac_results/transomics/factors/SKM_GN/network_Factor5.png",
    "motrpac/network_Factor6.png":
        "data/motrpac_results/transomics/factors/SKM_GN/network_Factor6.png",
    "motrpac/cross_tissue_changes.png":
        "data/motrpac_results/transomics/cross_tissue_changes.png",
    # notebooks/studies/obese_liver.py -- our own figures of our own computed
    # results; the Uematsu data they are derived from is never redistributed.
    "published/concordance.png": "data/published_results/obese_liver/concordance.png",
    "published/regulatory_paths.png":
        "data/published_results/obese_liver/regulatory_paths.png",
    "published/network.png": "data/published_results/obese_liver/network.png",
    "published/axes_by_contrast.png": "data/published_results/obese_liver/axes_by_contrast.png",
    "published/factor_scores.png":
        "data/published_results/obese_liver/factors/factor_scores.png",
    "published/network_Factor1.png":
        "data/published_results/obese_liver/factors/network_Factor1.png",
    "published/network_Factor4.png":
        "data/published_results/obese_liver/factors/network_Factor4.png",
    # notebooks/studies/kokaji_liver.py -- our figures of our own results; the
    # Kokaji supplement they derive from is never redistributed.
    "kokaji/network_WT.png": "data/published_results/kokaji_liver/network_WT.png",
    "kokaji/network_obob.png": "data/published_results/kokaji_liver/network_obob.png",
    "kokaji/axes_by_genotype.png":
        "data/published_results/kokaji_liver/axes_by_genotype.png",
    "kokaji/regulatory_paths.png":
        "data/published_results/kokaji_liver/regulatory_paths.png",
    "kokaji/convergence_null.png":
        "data/published_results/kokaji_liver/convergence_null.png",
    "kokaji/regulatory_motifs.png":
        "data/published_results/kokaji_liver/regulatory_motifs.png",
    "kokaji/metabolite_regulators.png":
        "data/published_results/kokaji_liver/metabolite_regulators.png",
    "kokaji/tf_activity.png": "data/published_results/kokaji_liver/tf_activity.png",
    "kokaji/genotypes_compared.png":
        "data/published_results/kokaji_liver/genotypes_compared.png",
}


def main() -> int:
    copied, missing = 0, []
    for published, source in sorted(FIGURES.items()):
        origin = ROOT / source
        if not origin.exists():
            missing.append((published, source))
            continue
        target = DEST / published
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(origin, target)
        copied += 1
        print(f"  {published}")

    print(f"\n{copied} figure(s) refreshed into {DEST.relative_to(ROOT)}")
    if missing:
        print(f"\n{len(missing)} missing -- rerun the analysis that writes them:")
        for published, source in missing:
            print(f"  {published:<38} <- {source}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
