"""The walkthroughs must run, offline, from a clean checkout.

They are the documentation's executable half: if one breaks, the docs build
fails and the README's first instruction is wrong. Marked slow, since running
all of them takes about a minute.
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
WALKTHROUGHS = sorted((ROOT / "notebooks" / "walkthroughs").glob("*.py"))
STUDIES = sorted((ROOT / "notebooks" / "studies").glob("*.py"))

#: The teaching order, which the filenames no longer carry. The Makefile and
#: the documentation toctree state the same order, and this test is what keeps
#: a new walkthrough from being added to one and forgotten in the others.
WALKTHROUGH_ORDER = [
    "build_network", "responsive_network", "reaction_regulation",
    "regulatory_paths", "temporal_and_hubs", "compare_conditions",
    "network_topology", "export_network", "transcription_factors",
    "external_annotation",
]

#: Walkthroughs that call live databases, so they cannot run offline.
NEEDS_NETWORK = {"external_annotation"}


def test_every_walkthrough_has_a_place_in_the_order():
    assert sorted(p.stem for p in WALKTHROUGHS) == sorted(WALKTHROUGH_ORDER)


def test_the_documentation_lists_them_in_that_order():
    index = (ROOT / "docs" / "source" / "index.rst").read_text()
    listed = [line.strip().split("/")[-1] for line in index.splitlines()
              if line.strip().startswith("notebooks/walkthroughs/")]
    assert listed == WALKTHROUGH_ORDER


def test_the_studies_are_all_present():
    assert sorted(p.stem for p in STUDIES) == [
        "brown_adipocytes", "liver_timecourse", "motrpac_rat", "obese_liver",
    ]


@pytest.mark.slow
@pytest.mark.parametrize("notebook",
                         [p for p in WALKTHROUGHS if p.stem not in NEEDS_NETWORK],
                         ids=lambda p: p.stem)
def test_walkthrough_runs_offline(notebook, tmp_path):
    environment = {**os.environ, "MPLBACKEND": "Agg"}
    result = subprocess.run([sys.executable, str(notebook)], cwd=tmp_path,
                            capture_output=True, text=True, env=environment, timeout=600)
    assert result.returncode == 0, result.stderr[-2000:]


@pytest.mark.parametrize("notebook", WALKTHROUGHS + STUDIES, ids=lambda p: p.stem)
def test_notebook_parses_as_a_notebook(notebook):
    """jupytext must be able to read it, since that is how the docs render it."""
    jupytext = pytest.importorskip("jupytext")
    parsed = jupytext.reads(notebook.read_text(), fmt="py:light")
    assert parsed.cells, f"{notebook.name} has no cells"
    assert any(cell.cell_type == "markdown" for cell in parsed.cells), \
        f"{notebook.name} has no prose -- results need stating, not just printing"

