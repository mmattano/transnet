"""KEGG's REST endpoints, exercised against the live service.

Marked ``network`` and deselected by default. These exist because KEGG changes
endpoint spellings without notice -- `/link/ec/<organism>` began returning 400
and silently broke every organism network build. A unit test with a mocked
response would not have caught that; only a live call does.

Run with:  pytest -m network tests/test_kegg_endpoints.py
"""

import pytest

pytestmark = pytest.mark.network


def test_gene_to_ec_links_resolve():
    """`/link/ec/<organism>` now 400s; `/link/enzyme/<organism>` is the spelling."""
    from transnet.api.kegg import kegg_link_ec

    table = kegg_link_ec("sce")
    assert not table.empty
    assert list(table.columns) == ["kegg_gene_id", "ec_number"]
    # EC numbers, not gene ids, must be in the second column.
    assert table["ec_number"].str.match(r"^\d+\.").any()


def test_gene_to_pathway_links_resolve():
    from transnet.api.kegg import kegg_link_pathway

    table = kegg_link_pathway("sce")
    assert not table.empty
    assert list(table.columns) == ["pathway", "kegg_gene_id"]


def test_pathway_listing_resolves():
    from transnet.api.kegg import kegg_list_pathways

    table = kegg_list_pathways("sce")
    assert not table.empty
    assert list(table.columns) == ["pathways_id", "description"]


def test_compound_listing_resolves():
    from transnet.api.kegg import kegg_list_compounds

    table = kegg_list_compounds()
    assert len(table) > 10000, "KEGG should list tens of thousands of compounds"


def test_ncbi_conversion_resolves():
    from transnet.api.kegg import kegg_conv_ncbi_idtable

    table = kegg_conv_ncbi_idtable("sce")
    assert not table.empty
    assert list(table.columns) == ["ncbi_id", "kegg_id"]


def test_ec_to_compound_links_resolve():
    from transnet.api.kegg import kegg_ec_to_cpds

    table = kegg_ec_to_cpds()
    assert not table.empty
