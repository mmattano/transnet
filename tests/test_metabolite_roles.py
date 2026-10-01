"""Which differential metabolites carry a BRENDA regulatory function?

The question this answers is: of the metabolites that changed between two
conditions, how many act back on an enzyme rather than merely flowing through
the pathway.
"""

import pandas as pd
import pytest

from transnet.analysis.transomics import (
    map_omics_to_network,
    metabolite_regulatory_roles,
    regulatory_role_enrichment,
)


@pytest.fixture
def mapped(standard_network, omics_tables):
    map_omics_to_network(
        standard_network,
        {k: v for k, v in omics_tables.items() if k != "Signaling"},
        id_column="id", log2fc_column="log2FC", qvalue_column="padj",
    )
    return standard_network


class TestRoles:
    def test_an_allosteric_inhibitor_is_labelled_as_one(self, mapped):
        roles = metabolite_regulatory_roles(mapped).set_index("metabolite")
        # G6P inhibits hexokinase in the fixture
        assert roles.loc["C00668", "role"] == "inhibitor"
        assert bool(roles.loc["C00668", "is_allosteric_regulator"])
        assert roles.loc["C00668", "n_reactions_inhibited"] == 1

    def test_an_allosteric_activator_is_labelled_as_one(self, mapped):
        roles = metabolite_regulatory_roles(mapped).set_index("metabolite")
        # F1,6BP activates PFK
        assert roles.loc["C00354", "role"] == "activator"
        assert roles.loc["C00354", "n_reactions_activated"] == 1

    def test_a_pure_participant_is_not_called_a_regulator(self, mapped):
        roles = metabolite_regulatory_roles(mapped).set_index("metabolite")
        # Glucose is only a substrate
        assert roles.loc["C00031", "role"] == "substrate/product only"
        assert not bool(roles.loc["C00031", "is_allosteric_regulator"])
        assert bool(roles.loc["C00031", "is_substrate"])

    def test_dual_role_metabolites_are_marked_both(self, mapped):
        roles = metabolite_regulatory_roles(mapped).set_index("metabolite")
        # ATP inhibits PFK and is a substrate of two reactions
        assert roles.loc["C00002", "role"] == "inhibitor"
        assert bool(roles.loc["C00002", "is_substrate"])

    def test_the_regulated_reactions_are_named(self, mapped):
        roles = metabolite_regulatory_roles(mapped).set_index("metabolite")
        assert "R00299" in roles.loc["C00668", "reactions_inhibited"]

    def test_measured_values_are_carried_through(self, mapped):
        roles = metabolite_regulatory_roles(mapped).set_index("metabolite")
        assert roles.loc["C00668", "regulated"] == 1
        assert roles.loc["C00668", "log2fc"] == 2.0

    def test_works_without_any_data_mapped(self, standard_network):
        roles = metabolite_regulatory_roles(standard_network)
        assert not roles.empty
        assert (roles["regulated"] == 0).all()
        assert roles["is_allosteric_regulator"].any()

    def test_network_without_allosteric_edges_reports_no_regulators(
        self, network_builder
    ):
        graph = network_builder("minimal").generate_graph()
        roles = metabolite_regulatory_roles(graph)
        assert not roles.empty
        assert not roles["is_allosteric_regulator"].any()

    def test_no_metabolome_layer_returns_an_empty_typed_frame(self, network_builder):
        graph = network_builder("no_reactions").generate_graph()
        # still has a Metabolome, but no reactions to regulate
        roles = metabolite_regulatory_roles(graph)
        assert "is_allosteric_regulator" in roles.columns


class TestEnrichment:
    def test_counts_answer_the_question_directly(self, mapped):
        result = regulatory_role_enrichment(metabolite_regulatory_roles(mapped))
        counts = result["counts"]

        assert counts["n_differential"] > 0
        assert counts["n_differential_regulators"] <= counts["n_differential"]
        assert 0.0 <= counts["fraction_differential_regulators"] <= 1.0

    def test_regulators_shortlist_only_contains_differential_regulators(self, mapped):
        result = regulatory_role_enrichment(metabolite_regulatory_roles(mapped))
        shortlist = result["regulators"]
        assert (shortlist["regulated"] != 0).all()
        assert shortlist["is_allosteric_regulator"].all()

    def test_enrichment_has_a_row_per_role(self, mapped):
        result = regulatory_role_enrichment(metabolite_regulatory_roles(mapped))
        assert set(result["enrichment"]["role"]) == {"any", "activator", "inhibitor"}

    def test_two_by_two_counts_are_internally_consistent(self, mapped):
        roles = metabolite_regulatory_roles(mapped)
        result = regulatory_role_enrichment(roles)
        counts = result["counts"]

        row = result["enrichment"].set_index("role").loc["any"]
        assert (
            row["n_differential_with_role"] + row["n_differential_without_role"]
            == counts["n_differential"]
        )
        assert row["n_differential_with_role"] == counts["n_differential_regulators"]

    def test_background_defaults_to_measured_metabolites(self, mapped):
        roles = metabolite_regulatory_roles(mapped)
        result = regulatory_role_enrichment(roles)
        assert result["counts"]["n_background"] == int(roles["log2fc"].notna().sum())

    def test_explicit_background_is_honoured(self, mapped):
        roles = metabolite_regulatory_roles(mapped)
        subset = ["C00668", "C00031", "C00354"]
        result = regulatory_role_enrichment(roles, background=subset)
        assert result["counts"]["n_background"] == 3

    def test_enrichment_is_detected_when_regulators_are_the_changed_ones(
        self, standard_network
    ):
        # Only the two allosteric regulators move.
        map_omics_to_network(
            standard_network,
            {"Metabolome": pd.DataFrame({
                "id": ["C00668", "C00354", "C00031", "C00085", "C00008"],
                "log2FC": [2.0, 1.8, 0.01, 0.01, 0.01],
                "padj": [0.001, 0.001, 0.9, 0.9, 0.9],
            })},
            id_column="id", log2fc_column="log2FC", qvalue_column="padj",
        )
        result = regulatory_role_enrichment(
            metabolite_regulatory_roles(standard_network)
        )
        row = result["enrichment"].set_index("role").loc["any"]
        assert result["counts"]["fraction_differential_regulators"] == 1.0
        assert row["p_value"] < 0.5

    def test_empty_roles_table_is_handled(self):
        result = regulatory_role_enrichment(pd.DataFrame(columns=[
            "metabolite", "regulated", "is_allosteric_regulator", "role",
            "log2fc", "qvalue",
        ]))
        assert result["counts"] == {} or result["counts"]["n_background"] == 0

    def test_q_values_are_bounded(self, mapped):
        result = regulatory_role_enrichment(metabolite_regulatory_roles(mapped))
        q = result["enrichment"]["q_value"].dropna()
        assert ((q >= 0) & (q <= 1)).all()


class TestIdentifierTranslation:
    """Real omics tables are keyed by whatever the instrument produced."""

    def test_alias_map_translates_ensembl_to_node_ids(self, network_builder):
        from transnet import build_alias_id_map

        network = network_builder("no_signaling")
        network.transcriptome.genes[0].ensembl_id = "ENSMUSG00000000001"

        mapping = build_alias_id_map(network, "Transcriptome", "ensembl_id")
        assert mapping["ENSMUSG00000000001"] == network.transcriptome.genes[0].ncbi_id

    def test_alias_map_expands_list_valued_attributes(self, network_builder):
        from transnet import build_alias_id_map

        network = network_builder("no_signaling")
        network.proteome.proteins[0].ensembl_id = ["ENSMUSP1", "ENSMUSP2"]

        mapping = build_alias_id_map(network, "Proteome", "ensembl_id")
        assert mapping["ENSMUSP1"] == "P19367"
        assert mapping["ENSMUSP2"] == "P19367"

    def test_alias_map_tolerates_float_formatted_ids(self, network_builder):
        from transnet import build_alias_id_map

        network = network_builder("no_signaling")
        network.metabolome.metabolites[0].pubchem_id = "5793.0"

        mapping = build_alias_id_map(network, "Metabolome", "pubchem_id")
        assert mapping["5793"] == "C00031"
        assert mapping["5793.0"] == "C00031"

    def test_unknown_layer_is_rejected(self, network_builder):
        from transnet import build_alias_id_map

        with pytest.raises(ValueError, match="Unknown layer"):
            build_alias_id_map(network_builder("no_signaling"), "Nope", "x")

    def test_absent_layer_gives_an_empty_map(self, network_builder):
        from transnet import build_alias_id_map

        network = network_builder("minimal")     # no Proteome
        assert build_alias_id_map(network, "Proteome", "ensembl_id") == {}

    def test_alias_map_feeds_straight_into_mapping(self, standard_network,
                                                   network_builder):
        from transnet import build_alias_id_map, map_omics_to_network

        network = network_builder("no_signaling")
        network.metabolome.metabolites[0].pubchem_id = "5793"

        report = map_omics_to_network(
            standard_network,
            {"Metabolome": pd.DataFrame({
                "id": ["5793"], "log2FC": [1.5], "padj": [0.01]})},
            id_column="id", log2fc_column="log2FC", qvalue_column="padj",
            id_map={"Metabolome": build_alias_id_map(
                network, "Metabolome", "pubchem_id")},
        )
        assert report.per_layer.iloc[0]["n_matched"] == 1
        assert standard_network.nodes["C00031"]["regulated"] == 1


class TestChebiPrefixHandling:
    """KEGG writes `chebi:10`, ChEBI clients write `CHEBI:10`."""

    @pytest.mark.network
    def test_prefixed_and_bare_ids_resolve_alike(self):
        from transnet.api.kegg import chebi_to_kegg

        result = chebi_to_kegg(["CHEBI:4167", "chebi:4167", "4167"]).set_index(
            "chebi_compounds")["kegg_compounds"]
        assert set(result) == {"C00031"}, (
            "all three spellings of a ChEBI id must resolve to the same compound"
        )
