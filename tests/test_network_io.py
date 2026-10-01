"""Writing a network to CSV and reading it back, which had no test at all.

The study analyses read ``data/<organism>/<date>/`` rather than rebuilding, so
whatever these files drop is lost to every one of them. The gene symbol was
dropped this way: proteins came back labelled with their full UniProt
description, and every figure built from a loaded network showed
"Aryl hydrocarbon receptor nuclear translocator (ARNT protein) (Dioxin
receptor..." where it had shown "Arnt".
"""

import networkx as nx
import pytest

from transnet.io import read_network, write_network


@pytest.fixture
def graph():
    graph = nx.MultiDiGraph()
    graph.add_node("24379", layer="Transcriptome", node_type="Gene",
                   name="Gad1", symbol="Gad1", log2fc=1.2, qvalue=0.01)
    graph.add_node("P97875", layer="Proteome", node_type="Protein",
                   name="Aryl hydrocarbon receptor nuclear translocator "
                        "(ARNT protein) (Dioxin receptor)",
                   symbol="Arnt")
    graph.add_node("R00299", layer="Reactions", node_type="Reaction",
                   name="ATP:D-glucose 6-phosphotransferase", reversible=False)
    graph.add_node("C00031", layer="Metabolome", node_type="Metabolite",
                   name="D-Glucose; Grape sugar")
    graph.add_edge("24379", "P97875", edge_type="translation", sign=1)
    graph.add_edge("P97875", "R00299", edge_type="catalysis", sign=1, ec="2.7.1.1")
    graph.add_edge("C00031", "R00299", edge_type="substrate", sign=1, stoichiometry=1)
    return graph


class TestRoundTrip:

    def test_the_gene_symbol_survives(self, graph, tmp_path):
        write_network(graph, str(tmp_path))
        back = read_network(str(tmp_path / "interactions.csv"),
                            nodes_file=str(tmp_path / "nodes.csv"))
        assert back.nodes["P97875"]["symbol"] == "Arnt"
        assert back.nodes["P97875"]["name"].startswith("Aryl hydrocarbon")

    def test_layers_types_and_names_survive(self, graph, tmp_path):
        write_network(graph, str(tmp_path))
        back = read_network(str(tmp_path / "interactions.csv"),
                            nodes_file=str(tmp_path / "nodes.csv"))
        assert back.number_of_nodes() == graph.number_of_nodes()
        for node, data in graph.nodes(data=True):
            assert back.nodes[node]["layer"] == data["layer"]
            assert back.nodes[node]["node_type"] == data["node_type"]
            assert back.nodes[node]["name"] == data["name"]

    def test_edge_types_and_signs_survive(self, graph, tmp_path):
        write_network(graph, str(tmp_path))
        back = read_network(str(tmp_path / "interactions.csv"),
                            nodes_file=str(tmp_path / "nodes.csv"))
        kinds = sorted(d.get("edge_type") for _, _, d in back.edges(data=True))
        assert kinds == ["catalysis", "substrate", "translation"]
        assert back.number_of_edges() == graph.number_of_edges()

    def test_a_table_without_the_symbol_column_still_loads(self, graph, tmp_path):
        """Networks exported before the column existed must keep working."""
        import pandas as pd

        write_network(graph, str(tmp_path))
        nodes = pd.read_csv(tmp_path / "nodes.csv").drop(columns=["Symbol"])
        nodes.to_csv(tmp_path / "nodes.csv", index=False)

        back = read_network(str(tmp_path / "interactions.csv"),
                            nodes_file=str(tmp_path / "nodes.csv"))
        assert "symbol" not in back.nodes["P97875"]
        assert back.nodes["P97875"]["name"].startswith("Aryl hydrocarbon")


def test_the_builder_writes_the_same_columns_the_reader_expects(tmp_path):
    """``Transnet.save_network`` and ``write_network`` are separate writers;
    a network saved by the builder must read back with its symbols too."""
    import pandas as pd

    from transnet import Transnet
    from transnet.biology.layers import Proteome
    from transnet.biology.elements import Protein

    network = Transnet(name="tiny")
    protein = Protein()
    protein.uniprot_id = "P97875"
    protein.name = "Aryl hydrocarbon receptor nuclear translocator (ARNT protein)"
    protein.gene = ["Arnt"]
    proteome = Proteome()
    proteome.proteins = [protein]
    network.proteome = proteome

    network.save_network(str(tmp_path))
    nodes = pd.read_csv(tmp_path / "nodes.csv")
    assert nodes.loc[nodes["ID"] == "P97875", "Symbol"].iloc[0] == "Arnt"
