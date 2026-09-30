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

#: Per study: the folder under docs/source/figures, the folder the notebook
#: writes to, and the figures the study page shows (published name, file
#: relative to the output folder).
STUDIES = {
    # notebooks/studies/brown_adipocytes.py
    "brown_adipocytes": ("data/brown_adipocyte_results", [
        "network.png", "layer_connectivity.png", "communities.png",
        "axis_composition.png", "regulation_axes.png", "controversial_reactions.png",
        "concordance.png", "tf_activity.png", "metabolite_regulators.png",
        "transomic_hubs.png", "downstream_influence.png", "temporal_structure.png",
        "early_vs_late.png", "regulatory_paths.png", "regulatory_motifs.png",
        "convergence_null.png", "structural_vulnerability.png",
        ("factor_scores.png", "factors/factor_scores.png"),
        ("factor_overview.png", "factors/factor_overview.png"),
        ("network_Factor3.png", "factors/network_Factor3.png"),
    ]),
    # notebooks/studies/motrpac_rat.py
    "motrpac": ("data/motrpac_results/transomics", [
        "network_SKM_GN.png", "layer_connectivity.png", "axes_by_tissue.png",
        "modification_sites.png", "phospho_axis.png",
        ("regulatory_paths.png", "paths_SKM_GN.png"),
        "concordance_SKM_GN.png", "edge_jaccard.png", "cross_tissue_changes.png",
        "regulation_axes_SKM_GN.png", "controversial_SKM_GN.png",
        "metabolite_regulators_SKM_GN.png", "tf_activity_SKM_GN.png",
        "downstream_influence_SKM_GN.png", "transomic_hubs_SKM_GN.png",
        "closest_tissues.png", "temporal_structure_SKM_GN.png",
        "regulatory_motifs_SKM_GN.png", "convergence_null_SKM_GN.png",
        "structural_vulnerability_SKM_GN.png",
        ("factor_overview.png", "factors/SKM_GN/factor_overview.png"),
        ("factor_scores.png", "factors/SKM_GN/factor_scores.png"),
        ("network_Factor5.png", "factors/SKM_GN/network_Factor5.png"),
        ("network_Factor6.png", "factors/SKM_GN/network_Factor6.png"),
    ]),
    # notebooks/studies/obese_liver.py -- figures of results computed here; the
    # Uematsu data they derive from is never redistributed.
    "obese_liver_panel": ("data/published_results/obese_liver", [
        "network.png", "layer_connectivity.png", "regulation_axes.png",
        "controversial_reactions.png", "metabolite_regulators.png", "concordance.png",
        "transomic_hubs.png", "downstream_influence.png", "regulatory_paths.png",
        "axes_by_contrast.png",
        ("factor_scores.png", "factors/factor_scores.png"),
        ("network_Factor1.png", "factors/network_Factor1.png"),
        ("network_Factor4.png", "factors/network_Factor4.png"),
    ]),
    # notebooks/studies/liver_timecourse.py -- figures of results computed here;
    # the Kokaji supplement they derive from is never redistributed.
    "liver_timecourse": ("data/published_results/liver_timecourse", [
        "network_WT.png", "network_obob.png", "layer_connectivity.png",
        "metabolites_over_time.png", "gene_overlap.png", "axes_by_genotype.png",
        "controversial_WT.png", "regulation_axes_WT.png", "regulation_axes_obob.png",
        "metabolite_regulators.png", "saturation.png", "tf_activity.png",
        "temporal_structure.png", "transomic_hubs.png", "regulatory_paths.png",
        "downstream_influence.png", "convergence_null.png", "regulatory_motifs.png",
        "structural_vulnerability.png", "genotypes_compared.png",
    ]),
}

#: published path -> where the analysis writes it.
FIGURES = {
    f"{folder}/{entry if isinstance(entry, str) else entry[0]}":
        f"{source}/{entry if isinstance(entry, str) else entry[1]}"
    for folder, (source, entries) in STUDIES.items()
    for entry in entries
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
