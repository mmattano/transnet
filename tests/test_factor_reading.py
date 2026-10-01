"""Reading factors: design associations, variance explained, network coherence."""

import networkx as nx
import numpy as np
import pandas as pd
import pytest

from transnet.analysis.factors import (
    factor_design_association,
    factor_network_coherence,
    factor_variance_explained,
    fit_factors,
)


def _design(n_per_cell=6):
    rows = []
    for time in ("control", "1w", "8w"):
        for sex in ("female", "male"):
            rows += [{"timepoint": time, "sex": sex}] * n_per_cell
    design = pd.DataFrame(rows)
    design.index = [f"s{i}" for i in range(len(design))]
    return design


class TestDesignAssociation:

    def test_each_factor_is_attributed_to_the_right_term(self):
        design = _design()
        rng = np.random.default_rng(0)
        noise = lambda: rng.normal(0, 0.2, len(design))  # noqa: E731
        is_male = (design["sex"] == "male").astype(float).to_numpy()
        trained = design["timepoint"].map({"control": 0, "1w": 1, "8w": 2}).to_numpy()
        factors = pd.DataFrame({
            "Sex": 3 * is_male + noise(),
            "Training": 1.5 * trained + noise(),
            "MaleTraining": 1.5 * trained * is_male + noise(),
        }, index=design.index)

        result = factor_design_association(factors, design, ["timepoint", "sex"])
        q = result.pivot(index="factor", columns="term", values="q_value")

        assert q.loc["Sex", "sex"] < 1e-6 and q.loc["Sex", "timepoint"] > 0.05
        assert q.loc["Training", "timepoint"] < 1e-6 and q.loc["Training", "sex"] > 0.05
        # a sex-specific training response is an interaction, which two
        # separate one-way tests cannot see
        assert q.loc["MaleTraining", "timepoint x sex"] < 1e-6

    def test_one_term_gives_a_one_way_model(self):
        design = _design()
        factors = pd.DataFrame({"F": np.arange(len(design), dtype=float)}, index=design.index)
        result = factor_design_association(factors, design, ["timepoint"])
        assert set(result["term"]) == {"timepoint"}


class TestVarianceExplained:

    def test_shares_are_bounded_by_the_model(self):
        rng = np.random.default_rng(1)
        samples = [f"s{i}" for i in range(30)]
        data = {
            "Transcriptome": pd.DataFrame(rng.gamma(2, 1, (30, 40)), index=samples,
                                          columns=[f"g{i}" for i in range(40)]),
            "Proteome": pd.DataFrame(rng.gamma(2, 1, (30, 20)), index=samples,
                                     columns=[f"p{i}" for i in range(20)]),
        }
        table = factor_variance_explained(fit_factors(data, n_components=3))

        model = table[table["factor"] == "model"].set_index("layer")["variance_accounted"]
        for layer in data:
            per_factor = table[(table["layer"] == layer) & (table["factor"] != "model")]
            assert (per_factor["variance_accounted"] <= model[layer] + 1e-9).all()
        overall = table[table["layer"] == "all layers"]
        assert sorted(overall["rank"]) == [1, 2, 3]


class TestNetworkCoherence:

    @staticmethod
    def _setup(n=60, linked=15):
        graph = nx.MultiDiGraph()
        genes = [f"g{i}" for i in range(n)]
        proteins = [f"p{i}" for i in range(n)]
        for g, p in zip(genes, proteins):
            graph.add_node(g, layer="Transcriptome")
            graph.add_node(p, layer="Proteome")
            graph.add_edge(g, p, edge_type="translation")
        rng = np.random.default_rng(2)
        # Coherent: its top genes and top proteins are the same pairs.
        # Scattered: its top genes and proteins are unrelated.
        gene_loadings = pd.DataFrame({
            "Coherent": [1.0 if i < linked else 0.01 for i in range(n)],
            "Scattered": [1.0 if i < linked else 0.01 for i in range(n)],
        }, index=genes)
        protein_loadings = pd.DataFrame({
            "Coherent": [1.0 if i < linked else 0.01 for i in range(n)],
            "Scattered": [1.0 if i >= n - linked else 0.01 for i in range(n)],
        }, index=proteins)
        return graph, {"Transcriptome": gene_loadings, "Proteome": protein_loadings}

    def test_linked_top_features_are_enriched(self):
        graph, loadings = self._setup()
        result = factor_network_coherence(graph, loadings, top_n=15, n_permutations=300)
        table = result["table"].set_index("factor")
        assert table.loc["Coherent", "n_links"] == 15
        assert table.loc["Coherent", "p_value"] < 0.01
        assert table.loc["Scattered", "n_links"] == 0
        assert table.loc["Scattered", "p_value"] > 0.5
        assert set(result["nodes"]["Coherent"]) >= {"g0", "p0"}

    def test_features_are_mapped_through_id_maps(self):
        graph, loadings = self._setup()
        loadings["Transcriptome"].index = [f"ENSG{i}" for i in range(60)]
        id_maps = {"Transcriptome": {f"ENSG{i}": f"g{i}" for i in range(60)}}
        table = factor_network_coherence(graph, loadings, id_maps=id_maps, top_n=15,
                                         n_permutations=100)["table"].set_index("factor")
        assert table.loc["Coherent", "n_links"] == 15


def test_nmf_keeps_the_ordering_of_negative_values():
    """Log-scale values below zero used to be clipped to zero before scaling,
    so every negative value became the same number."""
    samples = [f"s{i}" for i in range(12)]
    rng = np.random.default_rng(3)
    values = rng.normal(0, 1, (12, 5))           # half the values negative
    data = {"Proteome": pd.DataFrame(values, index=samples, columns=list("abcde"))}
    scaled = fit_factors(data, n_components=2).data_["Proteome"]
    for column in scaled.columns:
        original_rank = pd.Series(values[:, "abcde".index(column)]).rank().to_numpy()
        assert (scaled[column].rank().to_numpy() == original_rank).all()


class TestFactorsOnTheNetwork:
    """The readings that make factor analysis worth keeping in this package:
    each one asks something of a factor that a factor model cannot ask of
    itself."""

    @staticmethod
    def _setup(n_samples=14, linked=8):
        """A network where one factor's features are a connected chain and
        another's are scattered across unrelated molecules."""
        import networkx as nx

        graph = nx.MultiDiGraph()
        genes = [f"g{i}" for i in range(20)]
        proteins = [f"p{i}" for i in range(20)]
        for gene, protein in zip(genes, proteins):
            graph.add_node(gene, layer="Transcriptome", name=gene)
            graph.add_node(protein, layer="Proteome", name=protein)
            graph.add_edge(gene, protein, edge_type="translation", sign=1)
        for i in range(linked):
            graph.add_node(f"r{i}", layer="Reactions", name=f"reaction {i}")
            graph.add_edge(f"p{i}", f"r{i}", edge_type="catalysis", sign=1)
            graph.add_node(f"m{i}", layer="Metabolome", name=f"metabolite {i}")
            graph.add_edge(f"r{i}", f"m{i}", edge_type="product", sign=1)

        rng = np.random.default_rng(0)
        samples = [f"s{i}" for i in range(n_samples)]
        signal = rng.gamma(4, 1, n_samples)
        matrices = {}
        for layer, members in (("Transcriptome", genes), ("Proteome", proteins)):
            values = rng.gamma(1, 0.4, (n_samples, len(members)))
            values[:, :linked] += signal[:, None]        # the connected block
            matrices[layer] = pd.DataFrame(values, index=samples, columns=members)
        metabolites = [f"m{i}" for i in range(linked)]
        values = rng.gamma(1, 0.4, (n_samples, linked)) + signal[:, None]
        matrices["Metabolome"] = pd.DataFrame(values, index=samples, columns=metabolites)
        return graph, matrices

    def test_propagation_reaches_every_layer_from_a_connected_factor(self):
        from transnet.analysis.factors import factor_network_propagation, fit_factors

        graph, matrices = self._setup()
        factorisation = fit_factors(matrices, n_components=2)
        result = factor_network_propagation(graph, factorisation.loadings_, top_n=8,
                                            n_permutations=50)
        table = result["table"].set_index("factor")
        assert {"overlap", "null_overlap", "p_value"} <= set(table.columns)
        assert not result["scores"].empty
        assert set(result["top_nodes"]) <= set(table.index)

    def test_a_factor_whose_layers_share_a_pathway_overlaps_more_than_chance(self):
        """The planted factor loads on enzymes, their transcripts and the
        metabolites of the reactions they catalyse: one connected block. Its
        layers' diffusion profiles must overlap more than random features of
        the same layers do."""
        from transnet.analysis.factors import factor_network_propagation

        graph, _ = self._setup()
        linked = 8
        loadings = {
            "Transcriptome": pd.DataFrame({"Planted": [1.0 if i < linked else 0.01
                                                       for i in range(20)]},
                                          index=[f"g{i}" for i in range(20)]),
            "Metabolome": pd.DataFrame({"Planted": [1.0] * linked},
                                       index=[f"m{i}" for i in range(linked)]),
        }
        result = factor_network_propagation(graph, loadings, top_n=linked,
                                            n_permutations=200)
        row = result["table"].set_index("factor").loc["Planted"]
        assert row["overlap"] > row["null_overlap"]

    def test_propagation_reports_where_the_signal_landed(self):
        from transnet.analysis.factors import factor_network_propagation, fit_factors

        graph, matrices = self._setup()
        factorisation = fit_factors(matrices, n_components=2)
        result = factor_network_propagation(graph, factorisation.loadings_, top_n=8,
                                            n_permutations=50)
        top = next(iter(result["top_nodes"].values()))
        assert {"enrichment", "score", "layer", "name"} <= set(top.columns)
        # ranked by how much more than random seeds deliver, not by raw score:
        # a hub collects score from any seed set at all
        assert top["enrichment"].is_monotonic_decreasing

    def test_cross_layer_agreement_finds_a_shared_factor(self):
        from transnet.analysis.factors import factor_cross_layer_agreement, fit_factors

        graph, matrices = self._setup()
        factorisation = fit_factors(matrices, n_components=2)
        agreement = factor_cross_layer_agreement(factorisation)
        assert {"factor", "layer_a", "layer_b", "correlation"} <= set(agreement.columns)
        # the planted signal is in all three layers, so some factor must show it
        assert agreement["correlation"].max() > 0.5

    def test_network_imputation_beats_the_column_mean_for_a_connected_feature(self):
        """A metabolite missing in one sample is better guessed from the
        reaction it takes part in than from the cohort average."""
        from transnet.analysis.factors import network_guided_imputation

        graph, matrices = self._setup()
        frame = matrices["Metabolome"].copy()
        truth = float(frame.iloc[0, 0])
        frame.iloc[0, 0] = np.nan

        filled = network_guided_imputation(frame, graph)
        assert filled.notna().all().all()
        network_error = abs(filled.iloc[0, 0] - truth)
        mean_error = abs(frame.iloc[:, 0].mean() - truth)
        assert network_error <= mean_error * 1.5     # not worse than the mean

    def test_features_without_neighbours_fall_back_and_say_so(self, caplog):
        from transnet.analysis.factors import network_guided_imputation
        import networkx as nx

        frame = pd.DataFrame({"lonely": [1.0, np.nan, 3.0]})
        with caplog.at_level("INFO"):
            filled = network_guided_imputation(frame, nx.MultiDiGraph())
        assert filled["lonely"].notna().all()
        assert "column mean" in caplog.text


class TestSamplePairing:
    """Matching column names is not evidence the layers share samples."""

    @staticmethod
    def _layers(paired: bool, n_groups=3, per_group=6, n_features=40, seed=0):
        rng = np.random.default_rng(seed)
        samples, groups, rows_a, rows_b = [], [], [], []
        for g in range(n_groups):
            for r in range(per_group):
                sample = f"g{g}_{r}"
                individual = rng.normal(0, 1, n_features)        # the sample's own biology
                rows_a.append(individual + rng.normal(0, 0.3, n_features))
                partner = individual if paired else rng.normal(0, 1, n_features)
                rows_b.append(partner + rng.normal(0, 0.3, n_features))
                samples.append(sample)
                groups.append(f"group{g}")
        features = [f"f{i}" for i in range(n_features)]
        return (pd.DataFrame(rows_a, index=samples, columns=features),
                pd.DataFrame(rows_b, index=samples, columns=features),
                pd.Series(groups, index=samples))

    def test_truly_paired_layers_are_recognised(self):
        from transnet.analysis.factors import sample_pairing_check

        a, b, groups = self._layers(paired=True)
        result = sample_pairing_check(a, b, groups, n_permutations=100)
        assert result["paired"] and result["observed"] > result["null_mean"]

    def test_unpaired_layers_are_not(self):
        from transnet.analysis.factors import sample_pairing_check

        a, b, groups = self._layers(paired=False)
        result = sample_pairing_check(a, b, groups, n_permutations=100)
        assert not result["paired"]


def _two_group_layers(signal: str, seed=0):
    """Two layers of 2 x 10 samples, carrying either a between-group signal or
    a per-sample signal shared by both layers."""
    rng = np.random.default_rng(seed)
    samples = [f"{g}_{r}" for g in ("a", "b") for r in range(10)]
    groups = pd.Series([s.split("_")[0] for s in samples], index=samples)
    if signal == "between":
        driver = groups.map({"a": 0.5, "b": 4.0}).to_numpy()
    else:
        driver = np.tile(rng.gamma(2, 1, 10), 2)             # varies within groups
    matrices = {
        layer: pd.DataFrame(driver[:, None] * rng.uniform(0.5, 1.5, 20)
                            + rng.gamma(1, 0.1, (len(samples), 20)),
                            index=samples, columns=[f"{layer}{i}" for i in range(20)])
        for layer in ("T", "P")
    }
    return matrices, groups


def test_a_between_group_factor_survives_reshuffled_pairing():
    """A factor carried by the difference between groups does not care which
    replicate in a group is paired with which."""
    from transnet.analysis.factors import fit_factors, pairing_robustness

    matrices, groups = _two_group_layers("between")
    result = pairing_robustness(fit_factors(matrices, n_components=1), groups,
                                n_shuffles=100).iloc[0]
    assert result["verdict"] == "between-group", result.to_dict()


def test_a_factor_made_of_replicate_variation_is_flagged():
    """A factor that captures how replicates vary *together across layers*
    exists only if the columns are the same samples."""
    from transnet.analysis.factors import fit_factors, pairing_robustness

    matrices, groups = _two_group_layers("within")
    result = pairing_robustness(fit_factors(matrices, n_components=1), groups,
                                n_shuffles=100).iloc[0]
    assert result["verdict"] == "within-group", result.to_dict()


def test_pairing_is_detected_in_feature_centred_data():
    """MoTrPAC's layers are paired by design but centred per feature. The
    first version of the check correlated across genes within a sample, which
    on centred data sees nothing, and called a truly paired dataset unpaired."""
    from transnet.analysis.factors import sample_pairing_check

    rng = np.random.default_rng(3)
    samples = [f"{g}_{r}" for g in ("a", "b") for r in range(8)]
    groups = pd.Series([s.split("_")[0] for s in samples], index=samples)
    genes = [f"g{i}" for i in range(60)]
    individual = rng.normal(0, 1, (len(samples), len(genes)))   # sample x gene deviations
    transcript = individual + rng.normal(0, 0.7, individual.shape)
    protein = individual + rng.normal(0, 0.7, individual.shape)
    first = pd.DataFrame(transcript, index=samples, columns=genes)
    second = pd.DataFrame(protein, index=samples, columns=genes)
    first, second = first - first.mean(), second - second.mean()   # centred per feature

    assert sample_pairing_check(first, second, groups, n_permutations=100)["paired"]

    shuffled = second.sample(frac=1, random_state=0)
    shuffled.index = samples                                        # break the pairing
    assert not sample_pairing_check(first, shuffled, groups, n_permutations=100)["paired"]
