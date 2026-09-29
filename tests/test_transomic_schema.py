"""The typed network: edge vocabulary, direction, sign, and identifier integrity."""

import networkx as nx
import pytest

from transnet.biology.schema import (
    EDGE_TYPES,
    INTERACTION_COLUMNS,
    available_edge_types,
    available_layers,
    top_layer_present,
)
from transnet.biology.transnet import to_simple_graph


class TestInteractionSchema:
    def test_every_edge_declares_a_known_type(self, standard_network):
        for _, _, data in standard_network.edges(data=True):
            assert data["edge_type"] in EDGE_TYPES

    def test_interaction_frame_has_canonical_columns(self, network_builder):
        df = network_builder("no_signaling").generate_interaction_df()
        assert list(df.columns) == INTERACTION_COLUMNS

    def test_graph_is_directed_multigraph_by_default(self, standard_network):
        assert isinstance(standard_network, nx.MultiDiGraph)

    def test_signs_match_the_edge_type_definitions(self, standard_network):
        for _, _, data in standard_network.edges(data=True):
            expected = EDGE_TYPES[data["edge_type"]].sign
            assert data["sign"] == expected


class TestReactionEdges:
    """Regression tests for the bug that disconnected the reaction layer."""

    def test_reaction_metabolite_edges_use_kegg_compound_ids(self, standard_network):
        # Reaction edges once carried compound *names* as node ids, so every
        # metabolite existed twice and the reaction layer was disconnected.
        for u, v, data in standard_network.edges(data=True):
            if data["edge_type"] == "substrate":
                assert u.startswith("C"), f"substrate node {u!r} is not a KEGG id"
            if data["edge_type"] == "product":
                assert v.startswith("C"), f"product node {v!r} is not a KEGG id"

    def test_every_reaction_edge_endpoint_is_a_real_node(self, standard_network):
        nodes = set(standard_network.nodes)
        for u, v, data in standard_network.edges(data=True):
            if data["edge_type"] in ("substrate", "product", "catalysis"):
                assert u in nodes and v in nodes

    def test_substrate_and_product_are_distinct_relationships(self, standard_network):
        types = available_edge_types(standard_network)
        assert types.get("substrate", 0) > 0
        assert types.get("product", 0) > 0

    def test_substrate_points_into_the_reaction_and_product_out(self, standard_network):
        for u, v, data in standard_network.edges(data=True):
            if data["edge_type"] == "substrate":
                assert standard_network.nodes[v]["layer"] == "Reactions"
            if data["edge_type"] == "product":
                assert standard_network.nodes[u]["layer"] == "Reactions"

    def test_stoichiometry_is_carried_onto_the_edge(self, standard_network):
        coefficients = [
            data.get("stoichiometry")
            for _, _, data in standard_network.edges(data=True)
            if data["edge_type"] in ("substrate", "product")
        ]
        assert coefficients and all(c == 1.0 for c in coefficients)

    def test_ec_number_survives_onto_catalysis_edges(self, standard_network):
        ecs = [
            data.get("ec")
            for _, _, data in standard_network.edges(data=True)
            if data["edge_type"] == "catalysis"
        ]
        assert ecs and all(ec for ec in ecs)

    def test_reversibility_reaches_the_node(self, standard_network):
        assert standard_network.nodes["R00299"]["reversible"] is False


class TestAllostericEdges:
    def test_allosteric_edges_run_metabolite_to_reaction(self, standard_network):
        found = False
        for u, v, data in standard_network.edges(data=True):
            if data["edge_type"].startswith("allosteric"):
                found = True
                assert standard_network.nodes[u]["layer"] == "Metabolome"
                assert standard_network.nodes[v]["layer"] == "Reactions"
        assert found, "fixture should produce allosteric edges"

    def test_inhibition_is_negative_and_activation_positive(self, standard_network):
        for _, _, data in standard_network.edges(data=True):
            if data["edge_type"] == "allosteric_inhibition":
                assert data["sign"] == -1
            if data["edge_type"] == "allosteric_activation":
                assert data["sign"] == 1

    def test_one_metabolite_can_hold_two_roles_on_one_reaction(self, standard_network):
        # G6P is both the product of hexokinase and its allosteric inhibitor.
        # A simple graph cannot represent this; the multigraph must.
        roles = {
            data["edge_type"]
            for _, _, data in standard_network.edges(data=True)
            if "C00668" in (_, _) or True
        }
        parallel = standard_network.get_edge_data("C00668", "R00299")
        assert parallel is not None
        assert "allosteric_inhibition" in parallel
        product = standard_network.get_edge_data("R00299", "C00668")
        assert product is not None and "product" in product


class TestSimpleGraphProjection:
    def test_projection_is_undirected_and_simple(self, standard_network):
        simple = to_simple_graph(standard_network)
        assert not simple.is_directed()
        assert not simple.is_multigraph()

    def test_projection_keeps_every_node_and_its_attributes(self, standard_network):
        simple = to_simple_graph(standard_network)
        assert set(simple.nodes) == set(standard_network.nodes)
        for node in simple.nodes:
            assert simple.nodes[node]["layer"] == standard_network.nodes[node]["layer"]

    def test_projection_records_the_collapsed_relationships(self, standard_network):
        simple = to_simple_graph(standard_network)
        collapsed = simple["C00668"]["R00299"]["edge_types"]
        assert set(collapsed) == {"product", "allosteric_inhibition"}

    def test_projection_works_with_standard_networkx_algorithms(self, standard_network):
        simple = to_simple_graph(standard_network)
        centrality = nx.degree_centrality(simple)
        assert len(centrality) == simple.number_of_nodes()


class TestLayerIntrospection:
    def test_available_layers_returns_hierarchy_order(self, full_network):
        layers = available_layers(full_network)
        assert layers[0] == "Signaling"
        assert layers.index("Proteome") < layers.index("Metabolome")

    def test_top_layer_tracks_what_is_present(self, network_builder):
        full = network_builder("full").generate_graph()
        standard = network_builder("no_signaling").generate_graph()
        minimal = network_builder("minimal").generate_graph()

        assert top_layer_present(full) == "Signaling"
        assert top_layer_present(standard) == "Proteome"
        # No Proteome, so the hierarchy starts one layer lower.
        assert top_layer_present(minimal) == "Transcriptome"

    def test_minimal_network_has_no_proteome_derived_edges(self, network_builder):
        minimal = network_builder("minimal").generate_graph()
        types = available_edge_types(minimal)
        assert "catalysis" not in types
        assert "allosteric_inhibition" not in types
        # but the metabolic skeleton is still there ...
        assert types.get("substrate", 0) > 0
        assert types.get("product", 0) > 0
        # ... and the transcriptome still reaches it, via EC annotations
        assert types.get("gene_catalysis", 0) > 0
        assert "Transcriptome" in available_layers(minimal)


class TestSignalingLayerIsOptional:
    def test_signaling_edges_appear_only_when_the_layer_is_present(self, network_builder):
        full = network_builder("full").generate_graph()
        standard = network_builder("no_signaling").generate_graph()

        assert "kinase_tf" in available_edge_types(full)
        assert "kinase_tf" not in available_edge_types(standard)

    def test_network_without_signaling_is_still_well_formed(self, standard_network):
        assert standard_network.number_of_edges() > 0
        assert "Signaling" not in available_layers(standard_network)


class TestNoReactionsLayer:
    def test_enzymatic_shortcut_replaces_the_reaction_chain(self, network_builder):
        graph = network_builder("no_reactions").generate_graph()
        types = available_edge_types(graph)
        assert "catalysis" not in types
        assert "substrate" not in types
        # Nothing to shortcut to in this fixture (no protein.metabolites set),
        # but the network must still build and stay usable.
        assert "Reactions" not in available_layers(graph)
        assert graph.number_of_nodes() > 0


class TestBrendaNameResolution:
    """BRENDA reports effectors by name; the metabolome is keyed by KEGG id.

    Matching the two directly produced no allosteric edges at all, which is
    silent rather than loud -- the network simply came out with no metabolite
    regulation axis.
    """

    @staticmethod
    def _network(effector_names):
        from transnet import Transnet
        from transnet.biology.elements import Metabolite, Protein, Reaction
        from transnet.biology.layers import Metabolome, Proteome, Reactions

        reactions = Reactions()
        reactions.reactions = [Reaction(
            id="R00756", name="PFK", enzyme=["2.7.1.11"],
            substrates=["C00085"], products=["C00354"],
            stoichiometry_substrates=[1.0], stoichiometry_products=[1.0],
        )]
        metabolome = Metabolome()
        metabolome.metabolites = [
            Metabolite(kegg_compound_id="C00002",
                       kegg_name="ATP; Adenosine 5'-triphosphate"),
            Metabolite(kegg_compound_id="C00158",
                       kegg_name="Citrate; Citric acid"),
            Metabolite(kegg_compound_id="C00085", kegg_name="D-Fructose 6-phosphate"),
            Metabolite(kegg_compound_id="C00354",
                       kegg_name="D-Fructose 1,6-bisphosphate"),
        ]
        proteome = Proteome()
        enzyme = Protein(uniprot_id="P12382", name="Pfkl", ec_number=["2.7.1.11"])
        enzyme.inhibitors = list(effector_names)
        proteome.proteins = [enzyme]

        return Transnet(
            reactions=reactions, metabolome=metabolome, proteome=proteome
        ).generate_graph()

    def _allosteric(self, graph):
        return {
            (u, v) for u, v, d in graph.edges(data=True)
            if d["edge_type"] == "allosteric_inhibition"
        }

    def test_a_common_name_resolves_to_its_kegg_id(self):
        graph = self._network(["ATP"])
        assert ("C00002", "R00756") in self._allosteric(graph)

    def test_a_synonym_resolves_too(self):
        # "citrate" is the second synonym in the KEGG name string
        graph = self._network(["citrate"])
        assert ("C00158", "R00756") in self._allosteric(graph)

    def test_matching_is_case_insensitive(self):
        assert self._allosteric(self._network(["atp"])) == \
               self._allosteric(self._network(["ATP"]))

    def test_a_kegg_id_still_works_directly(self):
        graph = self._network(["C00002"])
        assert ("C00002", "R00756") in self._allosteric(graph)

    def test_an_unknown_effector_is_skipped_not_fatal(self):
        graph = self._network(["ATP", "not-a-real-compound"])
        edges = self._allosteric(graph)
        assert ("C00002", "R00756") in edges
        assert len(edges) == 1


class TestBrendaEffectorQuality:
    """BRENDA lists every compound that changes an enzyme's activity in vitro.

    On mouse pyruvate kinase that meant erythropoietin, interleukin-2,
    kaempferol, luteolin and curcumin beside fructose-1,6-bisphosphate --
    while L-alanine, its classic liver inhibitor, was lost because BRENDA
    writes it "L-Ala".
    """

    @staticmethod
    def _network(effectors, n_filler=250, endogenous_only=False):
        from transnet import Transnet
        from transnet.biology.elements import Metabolite, Protein, Reaction
        from transnet.biology.layers import Metabolome, Proteome, Reactions

        reactions = Reactions()
        # pyruvate kinase, plus enough organism-catalysed reactions that the
        # endogenous filter has a representative set to judge against
        reactions.reactions = [Reaction(
            id="R00200", name="pyruvate kinase", enzyme=["2.7.1.40"],
            substrates=["C00074", "C00008"], products=["C00022", "C00002"],
            stoichiometry_substrates=[1.0, 1.0], stoichiometry_products=[1.0, 1.0],
        ), Reaction(
            id="R00258", name="alanine transaminase", enzyme=["2.6.1.2"],
            substrates=["C00041", "C00026"], products=["C00022", "C00025"],
            stoichiometry_substrates=[1.0, 1.0], stoichiometry_products=[1.0, 1.0],
        )] + [Reaction(
            id=f"R9{i:04d}", name="filler", enzyme=["2.7.1.40"],
            substrates=["C00354"], products=["C00085"],
            stoichiometry_substrates=[1.0], stoichiometry_products=[1.0],
        ) for i in range(n_filler)]

        metabolome = Metabolome()
        metabolome.metabolites = [
            Metabolite(kegg_compound_id=c, kegg_name=n) for c, n in [
                ("C00074", "Phosphoenolpyruvate; PEP"), ("C00008", "ADP"),
                ("C00022", "Pyruvate"), ("C00002", "ATP"),
                ("C00041", "L-Alanine; L-2-Aminopropionic acid"),
                ("C00026", "2-Oxoglutarate"), ("C00025", "L-Glutamate"),
                ("C00354", "D-Fructose 1,6-bisphosphate"),
                ("C00085", "D-Fructose 6-phosphate"),
                ("C10107", "Kaempferol"),        # plant flavonoid
                ("C00284", "EDTA"),              # assay chelator
            ]
        ]
        proteome = Proteome()
        enzyme = Protein(uniprot_id="P53657", name="Pklr", ec_number=["2.7.1.40"])
        enzyme.activators = [e for e, kind in effectors if kind == "+"]
        enzyme.inhibitors = [e for e, kind in effectors if kind == "-"]
        proteome.proteins = [enzyme, Protein(uniprot_id="P25409", name="Gpt",
                                             ec_number=["2.6.1.2"])]
        network = Transnet(reactions=reactions, metabolome=metabolome, proteome=proteome)
        return network.allosteric_interactions(endogenous_only=endogenous_only)

    @staticmethod
    def _pairs(edges, reaction="R00200"):
        return {(e["source"], e["edge_type"]) for e in edges if e["target"] == reaction}

    def test_three_letter_amino_acid_codes_resolve(self):
        pairs = self._pairs(self._network([("L-Ala", "-")]))
        assert ("C00041", "allosteric_inhibition") in pairs

    def test_diphosphate_resolves_to_bisphosphate(self):
        pairs = self._pairs(self._network([("D-fructose 1,6-diphosphate", "+")]))
        assert ("C00354", "allosteric_activation") in pairs

    def test_by_default_every_resolved_effector_is_kept(self):
        """No silent curation: BRENDA's data enter the network as reported."""
        pairs = self._pairs(self._network([
            ("D-fructose 1,6-bisphosphate", "+"), ("kaempferol", "-"), ("EDTA", "+"),
        ]))
        assert {("C00354", "allosteric_activation"), ("C10107", "allosteric_inhibition"),
                ("C00284", "allosteric_activation")} <= pairs

    def test_xenobiotics_are_dropped_on_request(self):
        pairs = self._pairs(self._network([
            ("D-fructose 1,6-bisphosphate", "+"), ("kaempferol", "-"), ("EDTA", "+"),
        ], endogenous_only=True))
        assert ("C00354", "allosteric_activation") in pairs
        assert not any(source in ("C10107", "C00284") for source, _ in pairs)

    def test_a_small_network_keeps_everything_even_when_filtering(self):
        """Too few reactions to judge: ATP must not be dropped as foreign."""
        pairs = self._pairs(self._network([("kaempferol", "-")], n_filler=0,
                                          endogenous_only=True))
        assert ("C10107", "allosteric_inhibition") in pairs


class TestTranslationFallsBackToGeneSymbol:
    """~5% of proteins get no Entrez id; Pklr and Pcx were among them."""

    @staticmethod
    def _network(entrez, genes):
        from transnet import Transnet
        from transnet.biology.elements import Gene, Protein
        from transnet.biology.layers import Proteome, Transcriptome

        transcriptome = Transcriptome()
        transcriptome.genes = [Gene(ncbi_id=i, name=n) for i, n in genes]
        proteome = Proteome()
        proteome.proteins = [Protein(uniprot_id="P53657", gene=["Pklr"], entrez_id=entrez)]
        return Transnet(transcriptome=transcriptome, proteome=proteome)

    def test_entrez_is_used_when_present(self):
        edges = self._network(["18770"], [("18770", "Pklr")]).gene_protein_interaction()
        assert [(e["source"], e["evidence"]) for e in edges] == [("18770", "NCBI/UniProt")]

    def test_symbol_links_a_protein_without_entrez(self):
        edges = self._network([], [("18770", "Pklr")]).gene_protein_interaction()
        assert [(e["source"], e["target"], e["evidence"]) for e in edges] == [
            ("18770", "P53657", "gene symbol")
        ]

    def test_an_ambiguous_symbol_is_not_guessed(self):
        edges = self._network([], [("18770", "Pklr"), ("99999", "Pklr")]).gene_protein_interaction()
        assert edges == []
