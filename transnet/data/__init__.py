"""The bundled example trans-omic network.

A hand-curated slice of mouse hepatic glucose metabolism -- glycolysis,
gluconeogenesis, the pentose-phosphate entry point, and the insulin signalling
branch above them -- with real KEGG, UniProt and EC identifiers.  Everything
loads from disk in milliseconds and needs no network access or credentials, so
the examples and tests are reproducible anywhere.

The contrast in the bundled omics data is fed-versus-fasted-like: insulin
signalling raises glycolytic enzyme expression while product and allosteric
feedback pushes several of the same reactions the other way, which is what makes
the regulation-axis analysis interesting rather than trivially concordant.

::

    from transnet import load_example_network, load_example_omics

    graph = load_example_network()                 # no Signaling layer
    graph = load_example_network(signaling=True)   # with the kinase layer
    tables = load_example_omics()

Regenerate the files with ``python maintenance/build_example_data.py``.
"""

import logging
import os
from typing import Dict

import pandas as pd

logger = logging.getLogger(__name__)

__all__ = [
    "EXAMPLE_DIR",
    "load_example_network",
    "load_example_omics",
    "load_example_timecourse",
]

#: Directory holding the bundled example files.
EXAMPLE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "example")


def _path(name: str) -> str:
    path = os.path.join(EXAMPLE_DIR, name)
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Example file {name!r} not found at {path}. Regenerate the bundled "
            f"data with: python maintenance/build_example_data.py"
        )
    return path


def load_example_network(signaling: bool = False):
    """Load the example trans-omic network.

    Parameters
    ----------
    signaling : bool
        Add the optional Signaling layer (insulin receptor -> IRS1 -> PI3K ->
        AKT1, and AKT1's activating and inhibitory edges onto transcription
        factors).  ``False`` by default, because most studies do not measure
        the phosphoproteome -- and every analysis works either way.

    Returns
    -------
    networkx.MultiDiGraph
        Typed, directed, signed.  Roughly 60 nodes and 200 edges spanning
        Transcriptome, Proteome, Reactions and Metabolome (plus Signaling when
        requested).
    """
    from transnet.io.network_io import read_network

    graph = read_network(EXAMPLE_DIR)

    if signaling:
        edges = pd.read_csv(_path("signaling_edges.csv"))
        for _, row in edges.iterrows():
            kinase = str(row["kinase"])
            target = str(row["substrate"])
            if kinase not in graph:
                graph.add_node(
                    kinase, layer="Signaling", node_type="SignalingProtein",
                    name=row.get("kinase_name", kinase),
                )
            else:
                # A kinase measured in the proteome is promoted to the
                # signalling layer, where its regulatory role lives.
                graph.nodes[kinase]["layer"] = "Signaling"
                graph.nodes[kinase]["node_type"] = "SignalingProtein"

            edge_type = "kinase_tf" if row["target_type"] == "tf" else "phosphorylation"
            graph.add_edge(
                kinase, target, key=edge_type,
                edge_type=edge_type, role=edge_type,
                sign=int(row["sign"]), directed=True, weight=1.0,
                stoichiometry=None, ec="", source_db="KEGG",
                evidence="KEGG hsa04910 (insulin signaling)", confidence=None,
            )

    return graph


def load_example_omics(condition: str = "insulin_sensitive") -> Dict[str, pd.DataFrame]:
    """Load the example differential-expression tables.

    Parameters
    ----------
    condition : {"insulin_sensitive", "insulin_resistant"}
        Which contrast to load.  The two conditions are what
        :func:`~transnet.compare_transomic_networks` compares in example 06:
        in the resistant condition the signalling-driven arm is blunted while
        the allosteric arm still responds.

    Returns
    -------
    dict of str to pandas.DataFrame
        ``{layer: table}`` ready to hand to
        :func:`~transnet.map_omics_to_network`.  Each table has ``id``,
        ``log2FC`` and ``padj`` columns.
    """
    valid = {"insulin_sensitive", "insulin_resistant"}
    if condition not in valid:
        raise ValueError(f"condition must be one of {sorted(valid)}, got {condition!r}")

    if condition == "insulin_sensitive":
        return {
            "Transcriptome": pd.read_csv(_path("transcriptomics.csv"), dtype={"id": str}),
            "Proteome": pd.read_csv(_path("proteomics.csv")),
            "Metabolome": pd.read_csv(_path("metabolomics.csv")),
        }

    return {
        "Transcriptome": pd.read_csv(
            _path("transcriptomics_insulin_resistant.csv"), dtype={"id": str}),
        "Proteome": pd.read_csv(_path("proteomics_insulin_resistant.csv")),
        "Metabolome": pd.read_csv(_path("metabolomics_insulin_resistant.csv")),
    }


def load_example_phosphoproteomics() -> pd.DataFrame:
    """Phosphoproteomics for the optional Signaling layer."""
    return pd.read_csv(_path("phosphoproteomics.csv"))


def load_example_timecourse() -> pd.DataFrame:
    """Metabolite time course, in minutes after a glucose bolus.

    Returns
    -------
    pandas.DataFrame
        Columns ``id``, ``"0"``, ``"5"``, ``"15"``, ``"30"``, ``"60"``, ready
        for :func:`~transnet.assign_temporal_parameters`.
    """
    return pd.read_csv(_path("metabolome_timecourse.csv"))


def load_example_pathways() -> Dict[str, list]:
    """Pathway membership of the example reactions.

    Returns
    -------
    dict
        ``{reaction_id: [pathway, ...]}``, the ``pathway_map`` for
        :func:`~transnet.regulation_axis_summary`.
    """
    table = pd.read_csv(_path("reaction_pathways.csv"))
    return table.groupby("reaction")["pathway"].apply(list).to_dict()


__all__.append("load_example_phosphoproteomics")
__all__.append("load_example_pathways")
