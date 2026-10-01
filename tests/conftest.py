"""Shared fixtures: small trans-omic networks in several layer configurations.

The point of the parametrised ``configured_network`` fixture is that every
trans-omics analysis must run on all of them.  A study with no phosphoproteomics
and no proteomics is a normal study, not a broken input.
"""

import pandas as pd
import pytest

from transnet import Transnet
from transnet.biology.elements import (
    Gene, Metabolite, Protein, Reaction, SignalingProtein,
)
from transnet.biology.layers import (
    Metabolome, Proteome, Reactions, Signaling, Transcriptome,
)


def _reactions():
    """Hexokinase and phosphofructokinase, both feedback-inhibited."""
    layer = Reactions()
    layer.reactions = [
        Reaction(
            id="R00299", name="hexokinase", enzyme=["2.7.1.1"],
            substrates=["C00031", "C00002"], products=["C00668", "C00008"],
            stoichiometry_substrates=[1.0, 1.0],
            stoichiometry_products=[1.0, 1.0],
            reversible=False,
        ),
        Reaction(
            id="R00756", name="6-phosphofructokinase", enzyme=["2.7.1.11"],
            substrates=["C00085", "C00002"], products=["C00354", "C00008"],
            stoichiometry_substrates=[1.0, 1.0],
            stoichiometry_products=[1.0, 1.0],
            reversible=False,
        ),
    ]
    return layer


def _metabolome():
    layer = Metabolome()
    layer.metabolites = [
        Metabolite(kegg_compound_id=cid, kegg_name=name)
        for cid, name in [
            ("C00031", "D-Glucose"), ("C00002", "ATP"), ("C00668", "G6P"),
            ("C00008", "ADP"), ("C00085", "F6P"), ("C00354", "F1,6BP"),
        ]
    ]
    return layer


def _transcriptome():
    layer = Transcriptome()
    hk1 = Gene(ncbi_id="3098", name="HK1")
    pfkm = Gene(ncbi_id="5213", name="PFKM")
    srebf1 = Gene(ncbi_id="6720", name="SREBF1")
    # EC annotations let a transcriptome reach the reaction layer even when no
    # Proteome is available.
    hk1.related_ecs = ["2.7.1.1"]
    pfkm.related_ecs = ["2.7.1.11"]
    layer.genes = [hk1, pfkm, srebf1]
    return layer


def _proteome(with_tf=True, with_allostery=True):
    layer = Proteome()
    hk1 = Protein(uniprot_id="P19367", name="HK1",
                  ec_number=["2.7.1.1"], entrez_id=["3098"])
    pfkm = Protein(uniprot_id="P08237", name="PFKM",
                   ec_number=["2.7.1.11"], entrez_id=["5213"])
    srebf1 = Protein(uniprot_id="P36956", name="SREBF1", entrez_id=["6720"])

    if with_allostery:
        hk1.inhibitors = ["C00668"]        # G6P inhibits hexokinase
        pfkm.activators = ["C00354"]       # F1,6BP activates PFK
        pfkm.inhibitors = ["C00002"]       # ATP inhibits PFK

    if with_tf:
        srebf1.transcription_factor_targets = ["HK1", "PFKM"]

    hk1.interaction_partners = ["P08237"]
    layer.proteins = [hk1, pfkm, srebf1]
    return layer


def _signaling():
    layer = Signaling()
    layer.populate_from_table(
        pd.DataFrame({
            "kinase": ["P31749", "P31749"],
            "substrate": ["P36956", "P08237"],
            "sign": [1, 1],
            "target_type": ["tf", "protein"],
        }),
        sign_column="sign",
        target_type_column="target_type",
        evidence="test",
    )
    layer.nodes[0].name = "AKT1"
    return layer


def build_network(layers="full"):
    """Build one of the standard test configurations.

    Parameters
    ----------
    layers : {"full", "no_signaling", "minimal", "no_reactions"}
        ``full`` has every layer including Signaling; ``no_signaling`` is the
        common case; ``minimal`` is transcriptome + metabolome + reactions with
        no Proteome, so no catalysis or allosteric edges exist at all;
        ``no_reactions`` drops the reaction layer entirely.
    """
    kwargs = {}
    if layers in ("full", "no_signaling", "no_reactions"):
        kwargs["proteome"] = _proteome()
        kwargs["transcriptome"] = _transcriptome()
    if layers == "minimal":
        kwargs["transcriptome"] = _transcriptome()
    if layers != "no_reactions":
        kwargs["reactions"] = _reactions()
    kwargs["metabolome"] = _metabolome()
    if layers == "full":
        kwargs["signaling"] = _signaling()

    return Transnet(name=f"test-{layers}", **kwargs)


@pytest.fixture
def network_builder():
    """The :func:`build_network` factory, for tests that need a specific shape."""
    return build_network


@pytest.fixture(params=["full", "no_signaling", "minimal"])
def configured_network(request):
    """A trans-omic graph in each supported layer configuration."""
    transnet = build_network(request.param)
    graph = transnet.generate_graph()
    graph.graph["configuration"] = request.param
    return graph


@pytest.fixture
def full_network():
    """The complete six-layer network."""
    return build_network("full").generate_graph()


@pytest.fixture
def standard_network():
    """The common case: every layer except Signaling."""
    return build_network("no_signaling").generate_graph()


@pytest.fixture
def omics_tables():
    """Differential results for each layer, matching the fixture identifiers."""
    return {
        "Signaling": pd.DataFrame({
            "id": ["P31749"], "log2FC": [1.8], "padj": [0.001],
        }),
        "Transcriptome": pd.DataFrame({
            "id": ["3098", "5213", "6720"],
            "log2FC": [1.5, -1.2, 0.9],
            "padj": [0.01, 0.02, 0.03],
        }),
        "Proteome": pd.DataFrame({
            "id": ["P19367", "P08237"],
            "log2FC": [1.1, -0.9],
            "padj": [0.01, 0.04],
        }),
        "Metabolome": pd.DataFrame({
            "id": ["C00668", "C00031", "C00354", "C00002"],
            "log2FC": [2.0, -0.8, 1.4, 0.2],
            "padj": [0.001, 0.02, 0.01, 0.6],   # ATP not significant
        }),
    }
