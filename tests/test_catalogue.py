"""The documented analysis catalogue must match the actual API.

The catalogue table is the spine of the README and the docs. If a
function is renamed or removed, these tests fail rather than letting the
documentation quietly describe an API that no longer exists.
"""

import os
import re

import pytest

import transnet

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: Every analysis in the catalogue: the heading it carries in the docs and the
#: README, and the functions it promises. Keyed by name -- the A1..A10 codes
#: are gone, because a reader should not have to hold a lookup table.
CATALOGUE = {
    "Trans-omic network reconstruction": ["responsive_subnetwork"],
    "Reaction regulation-axis attribution": ["reaction_regulation_table"],
    "Per-pathway regulation balance": ["regulation_axis_summary"],
    "Metabolite regulatory roles": ["metabolite_regulatory_roles",
                                    "regulatory_role_enrichment"],
    "Signed regulatory-path tracing": ["trace_regulatory_paths",
                                       "path_consistency_summary"],
    "Cross-layer connectivity": ["cross_layer_connectivity", "layer_coverage"],
    "Trans-omic hub identification": ["transomic_hubs"],
    "Temporal and dose structure on the network": [
        "assign_temporal_parameters", "temporal_network_structure",
        "split_by_response_class"],
    "Responsive transcription-factor inference": ["transcription_factor_activity"],
    "Transcript–protein concordance": ["expression_concordance"],
    "Regulatory motifs": ["regulatory_motifs"],
    "Structural vulnerability": ["structural_vulnerability"],
    "Convergence against a null model": ["convergence_significance"],
}

SUPPORTING = [
    "map_omics_to_network", "regulated_nodes", "compare_transomic_networks",
    "hierarchical_propagation", "downstream_influence",
    "load_example_network", "load_example_omics", "load_example_timecourse",
    "to_simple_graph", "available_layers", "available_edge_types",
    "top_layer_present",
]

ALL_NAMES = [name for names in CATALOGUE.values() for name in names] + SUPPORTING


@pytest.mark.parametrize("name", ALL_NAMES)
def test_catalogue_name_is_importable_from_the_package_root(name):
    assert hasattr(transnet, name), (
        f"{name} is documented but not exported from `transnet`"
    )
    assert name in transnet.__all__, f"{name} is missing from transnet.__all__"


@pytest.mark.parametrize("name", ALL_NAMES)
def test_catalogue_function_is_documented(name):
    function = getattr(transnet, name)
    assert function.__doc__, f"{name} has no docstring"
    assert len(function.__doc__.strip()) > 60, f"{name} has a stub docstring"


def _section(heading):
    """The text of one catalogue section, from its heading to the next one."""
    docs = _read("docs/source/transomics_analyses.rst")
    match = re.search(rf"^{re.escape(heading)}\n-+\n(.*?)(?=^\S[^\n]*\n-{{3,}}\n|\Z)",
                      docs, re.MULTILINE | re.DOTALL)
    return match.group(1) if match else ""


@pytest.mark.parametrize("analysis", sorted(CATALOGUE))
def test_each_analysis_names_its_source(analysis):
    """Provenance lives in the catalogue, one :Reference: field per analysis."""
    assert ":Reference:" in _section(analysis), (
        f"'{analysis}' in docs/source/transomics_analyses.rst names no reference"
    )


def _read(path):
    with open(os.path.join(ROOT, path)) as handle:
        return handle.read()


@pytest.mark.parametrize("analysis", sorted(CATALOGUE))
def test_readme_lists_every_analysis(analysis):
    readme = _read("README.md")
    assert f"| {analysis}" in readme, f"{analysis} is missing from the README table"


@pytest.mark.parametrize("analysis", sorted(CATALOGUE))
def test_docs_have_a_section_per_analysis(analysis):
    docs = _read("docs/source/transomics_analyses.rst")
    assert re.search(rf"^{re.escape(analysis)}\s*$", docs, re.MULTILINE), (
        f"'{analysis}' has no section in docs/source/transomics_analyses.rst"
    )


def test_no_analysis_is_referred_to_by_a_code():
    """The A1..A10 numbering is gone: names only, everywhere a reader looks."""
    for path in ("README.md", "docs/source/transomics_analyses.rst",
                 "docs/source/brown_adipocytes.rst", "docs/source/motrpac_study.rst",
                 "docs/source/obese_liver_panel.rst"):
        found = re.findall(r"\bA\d+b?\b", _read(path))
        assert not found, f"{path} still refers to {sorted(set(found))}"


@pytest.mark.parametrize("name", [n for names in CATALOGUE.values() for n in names])
def test_docs_reference_each_function(name):
    docs = _read("docs/source/transomics_analyses.rst")
    assert name in docs, f"{name} is not mentioned in the analysis catalogue"


def test_every_documented_edge_type_exists():
    """The edge-type tables in the README and docs must match the schema."""
    from transnet.biology.schema import EDGE_TYPES

    for path in ("README.md", "docs/source/network_model.rst"):
        text = _read(path)
        for edge_type in EDGE_TYPES:
            assert f"`{edge_type}`" in text or f"``{edge_type}``" in text, (
                f"edge type {edge_type!r} is missing from {path}"
            )


def test_citations_file_covers_the_catalogue():
    """Every analysis names its source; the bibliography must hold them all."""
    citations = _read("CITATIONS.md")
    for names in CATALOGUE.values():
        for name in names:
            assert name in citations, f"{name} is missing from CITATIONS.md"


