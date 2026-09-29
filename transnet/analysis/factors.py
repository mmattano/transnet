"""Multi-omics factors, read through the trans-omic network.

Factor analysis by itself is not trans-omics: it asks which molecules co-vary,
not how a signal propagates through a biochemical network. It earns a place
here only because of the three readings at the bottom of this module, which
put a factor *on the network* and ask whether its strongest features are
connected, which part of the design it follows, and how much it reconstructs.

Only NMF is implemented, because only NMF is used: its parts-based, additive
factors are what the coherence test is interpretable on. For PCA, factor
analysis or CCA, use scikit-learn directly --
``sklearn.decomposition.PCA(n_components=5).fit_transform(matrix)`` is one
line, always current, and needs no wrapper here.

References
----------
Argelaguet R, et al. MOFA+: a statistical framework for comprehensive
integration of multi-modal single-cell data. *Genome Biology* 21:111, 2020.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import logging

import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import MinMaxScaler

logger = logging.getLogger(__name__)

__all__ = [
    "Factorisation",
    "fit_factors",
    "factor_design_association",
    "factor_variance_explained",
    "factor_network_coherence",
]


@dataclass
class Factorisation:
    """The result of :func:`fit_factors`.

    ``factors_`` is samples x factors (the scores); ``loadings_`` maps each
    layer to a features x factors frame; ``data_`` holds the scaled matrices
    the fit actually saw, which is what :func:`factor_variance_explained`
    reconstructs against.
    """

    factors_: pd.DataFrame
    loadings_: Dict[str, pd.DataFrame] = field(default_factory=dict)
    data_: Dict[str, pd.DataFrame] = field(default_factory=dict)

    @property
    def n_components(self) -> int:
        return self.factors_.shape[1]

    def top_features(self, layer: str, factor: str, n: int = 20) -> pd.Series:
        """The ``n`` features loading most strongly on one factor."""
        return self.loadings_[layer][factor].sort_values(ascending=False).head(n)


def fit_factors(
    matrices: Dict[str, pd.DataFrame],
    n_components: int = 6,
    random_state: int = 0,
    impute: bool = True,
) -> Factorisation:
    """Fit NMF factors jointly across layers (early fusion).

    Parameters
    ----------
    matrices : dict
        ``{layer: DataFrame}``, samples x features, sharing an index. Layers
        are concatenated feature-wise, so every layer contributes to every
        factor -- which is the point: a factor that loads on transcripts and
        metabolites together is the one worth taking to the network.
    n_components : int
        Number of factors.
    random_state : int
        Seed, so a rerun gives the same factors.
    impute : bool
        Replace missing values with the feature mean before fitting. NMF has
        no missing-value handling of its own.

    Returns
    -------
    Factorisation

    Notes
    -----
    Each feature is scaled to [0, 1] before fitting, not clipped at zero. Log
    ratios and log-CPM are negative for a large share of values (47% of the
    MoTrPAC proteome); clipping first collapsed all of them onto zero and
    discarded that half of the signal. The clip after scaling absorbs
    floating-point error only.
    """
    if not matrices:
        raise ValueError("no matrices to factorise")

    # Keep the first layer's sample order rather than sorting: "s10" sorts
    # before "s2", which would silently reorder every row against the design.
    shared = set.intersection(*(set(m.index) for m in matrices.values()))
    samples = [s for s in next(iter(matrices.values())).index if s in shared]
    if not samples:
        raise ValueError("the layers share no samples; align them on sample id first")

    scaled: Dict[str, pd.DataFrame] = {}
    for layer in sorted(matrices):
        frame = matrices[layer].loc[samples]
        if impute and frame.isnull().any().any():
            frame = pd.DataFrame(
                SimpleImputer(strategy="mean").fit_transform(frame),
                index=frame.index, columns=frame.columns,
            )
        values = np.clip(MinMaxScaler().fit_transform(frame), 0.0, 1.0)
        scaled[layer] = pd.DataFrame(values, index=frame.index, columns=frame.columns)

    joint = pd.concat([scaled[layer] for layer in sorted(scaled)], axis=1)
    logger.info(f"Fitting {n_components} NMF factors on {joint.shape[0]} samples "
                f"x {joint.shape[1]:,} features across {len(scaled)} layers")

    model = NMF(n_components=n_components, random_state=random_state, max_iter=500)
    scores = model.fit_transform(joint.to_numpy())
    names = [f"Factor{i + 1}" for i in range(n_components)]

    loadings, start = {}, 0
    for layer in sorted(scaled):
        width = scaled[layer].shape[1]
        loadings[layer] = pd.DataFrame(
            model.components_[:, start:start + width].T,
            index=scaled[layer].columns, columns=names,
        )
        start += width

    return Factorisation(
        factors_=pd.DataFrame(scores, index=joint.index, columns=names),
        loadings_=loadings,
        data_=scaled,
    )


def factor_design_association(
    factors: pd.DataFrame,
    design: pd.DataFrame,
    terms: List[str],
    interaction: bool = True,
) -> pd.DataFrame:
    """Test every factor against the study design in one linear model.

    Testing timepoint and sex separately, as one-way ANOVAs, asks the wrong
    question of a factorial design: a factor that carries a *sex-specific*
    training response looks like a sex factor with a weak timepoint effect.
    Here each factor's scores are modelled as ``score ~ term1 * term2``
    (main effects plus their interaction when two terms are given), with
    Type II sums of squares so each main effect is adjusted for the other.

    Parameters
    ----------
    factors : pandas.DataFrame
        Samples x factors, e.g. ``MultiOmicsIntegrator.factors_``.
    design : pandas.DataFrame
        Samples x covariates, indexed like ``factors``.
    terms : list of str
        One or two categorical columns of ``design``.
    interaction : bool
        Include the interaction when two terms are given.

    Returns
    -------
    pandas.DataFrame
        One row per factor and model term: ``factor``, ``term``, ``df``,
        ``F``, ``p_value``, ``partial_eta_squared`` (share of the variance not
        explained by other terms that this term explains) and ``q_value``
        (Benjamini-Hochberg across every factor and term).
    """
    import statsmodels.formula.api as smf
    from statsmodels.stats.anova import anova_lm
    from statsmodels.stats.multitest import multipletests

    common = factors.index.intersection(design.index)
    if len(common) < 4:
        raise ValueError("fewer than four samples shared by factors and design")
    data = design.loc[common, terms].astype(str).copy()
    safe = {term: f"t{i}" for i, term in enumerate(terms)}
    data = data.rename(columns=safe)

    formula_terms = [f"C({safe[t]})" for t in terms]
    rhs = " * ".join(formula_terms) if (interaction and len(terms) == 2) else " + ".join(formula_terms)
    labels = {f"C({safe[t]})": t for t in terms}
    if interaction and len(terms) == 2:
        labels[f"C({safe[terms[0]]}):C({safe[terms[1]]})"] = f"{terms[0]} x {terms[1]}"

    rows = []
    for factor in factors.columns:
        frame = data.assign(score=factors.loc[common, factor].astype(float).to_numpy())
        table = anova_lm(smf.ols(f"score ~ {rhs}", data=frame).fit(), typ=2)
        residual = table.loc["Residual", "sum_sq"]
        for key, label in labels.items():
            if key not in table.index:
                continue
            ss = table.loc[key, "sum_sq"]
            rows.append({
                "factor": factor, "term": label, "df": int(table.loc[key, "df"]),
                "F": float(table.loc[key, "F"]), "p_value": float(table.loc[key, "PR(>F)"]),
                "partial_eta_squared": float(ss / (ss + residual)) if (ss + residual) > 0 else np.nan,
            })
    result = pd.DataFrame(rows)
    valid = result["p_value"].notna()
    result["q_value"] = np.nan
    if valid.any():
        result.loc[valid, "q_value"] = multipletests(result.loc[valid, "p_value"], method="fdr_bh")[1]
    return result


def factor_variance_explained(factorisation: "Factorisation") -> pd.DataFrame:
    """How much of each layer every factor reconstructs.

    NMF factors, unlike principal components, come in no order of importance:
    "Factor 1" is a label, not a rank. This measures each one, so factors can
    be compared on what they explain rather than on their number. For a factor
    k the variance accounted for in layer L is

        VAF(k, L) = 1 - ||X_L - w_k h_{k,L}||^2 / ||X_L||^2

    on the scaled data the model was fitted to, alongside the whole model's.

    Returns
    -------
    pandas.DataFrame
        ``factor``, ``layer`` (and ``"all layers"``), ``variance_accounted``,
        plus a ``rank`` of factors by their all-layer share.
    """
    if factorisation.factors_ is None or not factorisation.loadings_:
        raise ValueError("fit the integrator before measuring variance")
    weights = factorisation.factors_.to_numpy()
    rows, totals = [], {}
    for layer, loadings in factorisation.loadings_.items():
        data = factorisation.data_[layer].loc[factorisation.factors_.index, loadings.index].to_numpy()
        norm = float((data ** 2).sum())
        full = weights @ loadings.to_numpy().T
        rows.append({"factor": "model", "layer": layer,
                     "variance_accounted": 1 - float(((data - full) ** 2).sum()) / norm})
        for k, factor in enumerate(loadings.columns):
            part = np.outer(weights[:, k], loadings[factor].to_numpy())
            residual = float(((data - part) ** 2).sum())
            rows.append({"factor": factor, "layer": layer,
                         "variance_accounted": 1 - residual / norm})
            totals.setdefault(factor, [0.0, 0.0])
            totals[factor][0] += residual
            totals[factor][1] += norm
    for factor, (residual, norm) in totals.items():
        rows.append({"factor": factor, "layer": "all layers",
                     "variance_accounted": 1 - residual / norm})
    table = pd.DataFrame(rows)
    overall = table[(table["layer"] == "all layers")].sort_values(
        "variance_accounted", ascending=False)
    rank = {f: i + 1 for i, f in enumerate(overall["factor"])}
    table["rank"] = table["factor"].map(rank)
    return table


def factor_network_coherence(
    graph,
    loadings: Dict[str, pd.DataFrame],
    id_maps: Optional[Dict[str, Dict[str, str]]] = None,
    top_n: int = 100,
    n_permutations: int = 1000,
    random_state: int = 0,
) -> Dict[str, object]:
    """Do a factor's top features connect *across layers* on the network?

    A factor built from co-variation alone can mix unrelated genes, proteins
    and metabolites. If its strongest features also sit on connected parts of
    the trans-omic network -- a gene and the protein it encodes, a factor and
    its target gene, an enzyme and a metabolite of the reaction it catalyses
    -- the factor describes one regulatory programme rather than a
    coincidence. This counts those links among each factor's top features and
    compares the count with random feature sets of the same size drawn from
    the same mappable features.

    Links counted: ``translation`` (gene-protein), ``transcriptional_regulation``
    (factor-gene) and protein-metabolite pairs joined through a reaction the
    protein catalyses (substrate, product or allosteric regulator).

    Parameters
    ----------
    graph : networkx.MultiDiGraph
    loadings : dict
        ``{layer: features x factors}`` -- ``MultiOmicsIntegrator.loadings_``.
    id_maps : dict, optional
        ``{layer: {feature id: node id}}`` for features not keyed like nodes.
    top_n : int
        Top features per layer, by loading.
    n_permutations : int
    random_state : int

    Returns
    -------
    dict
        ``table`` : one row per factor -- ``n_links`` observed, ``null_mean``,
        ``fold_enrichment``, ``p_value`` (one-sided permutation), ``q_value``.
        ``nodes`` : ``{factor: [network nodes of its linked top features]}``,
        for drawing the factor's subnetwork.
    """
    from statsmodels.stats.multitest import multipletests

    id_maps = id_maps or {}
    layer_of = {n: d.get("layer") for n, d in graph.nodes(data=True)}

    def to_nodes(layer, features):
        mapping = id_maps.get(layer, {})
        out = []
        for feature in features:
            node = mapping.get(str(feature), str(feature))
            if node in layer_of and layer_of[node] == layer:
                out.append(node)
        return out

    translation, regulates, reach = set(), set(), {}
    catalyses = {}
    reaction_metabolites = {}
    for u, v, data in graph.edges(data=True):
        kind = data.get("edge_type")
        if kind == "translation":
            translation.add((u, v))
        elif kind == "transcriptional_regulation":
            regulates.add((u, v))
        elif kind == "catalysis":
            catalyses.setdefault(u, set()).add(v)
        elif kind in ("substrate", "allosteric_activation", "allosteric_inhibition"):
            reaction_metabolites.setdefault(v, set()).add(u)
        elif kind == "product":
            reaction_metabolites.setdefault(u, set()).add(v)
    for protein, reactions in catalyses.items():
        reach[protein] = set().union(*(reaction_metabolites.get(r, set()) for r in reactions))

    def count(genes, proteins, metabolites):
        genes, proteins, metabolites = set(genes), set(proteins), set(metabolites)
        links = sum(1 for g in genes for p in proteins if (g, p) in translation)
        links += sum(1 for p in proteins for g in genes if (p, g) in regulates)
        links += sum(len(reach.get(p, set()) & metabolites) for p in proteins)
        return links

    pools = {layer: to_nodes(layer, frame.index) for layer, frame in loadings.items()}
    rng = np.random.default_rng(random_state)
    factors = next(iter(loadings.values())).columns
    rows, members = [], {}
    for factor in factors:
        chosen = {}
        for layer, frame in loadings.items():
            ranked = frame[factor].abs().sort_values(ascending=False).index
            chosen[layer] = to_nodes(layer, ranked)[:top_n]
        observed = count(chosen.get("Transcriptome", []), chosen.get("Proteome", []),
                         chosen.get("Metabolome", []))
        null = np.empty(n_permutations)
        for i in range(n_permutations):
            sample = {
                layer: rng.choice(pool, size=min(len(chosen[layer]), len(pool)), replace=False)
                if len(pool) else []
                for layer, pool in pools.items()
            }
            null[i] = count(sample.get("Transcriptome", []), sample.get("Proteome", []),
                            sample.get("Metabolome", []))
        rows.append({
            "factor": factor, "n_links": observed,
            "null_mean": float(null.mean()),
            "fold_enrichment": observed / null.mean() if null.mean() > 0 else np.nan,
            "p_value": (1 + int((null >= observed).sum())) / (1 + n_permutations),
        })
        genes, proteins, metabolites = (set(chosen.get(l, [])) for l in
                                        ("Transcriptome", "Proteome", "Metabolome"))
        linked = set()
        for g in genes:
            for p in proteins:
                if (g, p) in translation or (p, g) in regulates:
                    linked.update((g, p))
        for p in proteins:
            hit = reach.get(p, set()) & metabolites
            if hit:
                linked.add(p)
                linked.update(hit)
                linked.update(r for r in catalyses.get(p, set())
                              if reaction_metabolites.get(r, set()) & hit)
        members[factor] = sorted(linked)
    table = pd.DataFrame(rows)
    table["q_value"] = multipletests(table["p_value"], method="fdr_bh")[1]
    return {"table": table, "nodes": members}


# ---------------------------------------------------------------------------
# Factors, read through the network
# ---------------------------------------------------------------------------

def _mapped_top_features(loadings: Dict[str, pd.DataFrame], factor: str, graph,
                         id_maps: Optional[Dict[str, Dict[str, str]]], top_n: int):
    """The factor's strongest features per layer, as node ids in ``graph``.

    Returns ``{layer: {node: weight}}`` -- keyed by layer because a null that
    ignores which layers a factor loads on is not a null for that factor: the
    transcriptome is ten times the size of the metabolome, so uniform draws
    are effectively all transcripts.
    """
    seeds: Dict[str, Dict[str, float]] = {}
    for layer, frame in loadings.items():
        if factor not in frame.columns:
            continue
        ranked = frame[factor].sort_values(ascending=False).head(top_n)
        mapping = (id_maps or {}).get(layer, {})
        for feature, weight in ranked.items():
            node = mapping.get(str(feature), str(feature))
            if node in graph and weight > 0:
                layer_seeds = seeds.setdefault(layer, {})
                layer_seeds[node] = max(layer_seeds.get(node, 0.0), float(weight))
    return seeds


def _propagation_graph(graph, exclude_edge_types, exclude_nodes):
    """The graph diffusion should actually run on.

    Two classes of edge make propagation meaningless on a real network, and
    they are the same two that make path tracing meaningless: STRING
    associations, which connect almost any protein to almost any other, and
    currency metabolites, which sit on half the reactions. Left in, every
    factor's signal ends up on water and ADP and no seed set can be told from
    another.
    """
    import networkx as nx

    keep_edges = [
        (u, v) for u, v, data in graph.edges(data=True)
        if (data.get("edge_type") or "unknown") not in set(exclude_edge_types)
    ]
    trimmed = nx.Graph()
    trimmed.add_nodes_from(graph.nodes(data=True))
    trimmed.add_edges_from(keep_edges)
    trimmed.remove_nodes_from(set(exclude_nodes) & set(trimmed.nodes))
    trimmed.remove_edges_from(nx.selfloop_edges(trimmed))
    return trimmed


def _propagator(graph, alpha: float):
    """Build the RWR operator once, to be reused for many seed vectors.

    Running :func:`random_walk_with_restart` per seed set rebuilds the sparse
    matrix every time, which on an interactome is minutes per factor and made
    the permutation test unusable. The matrix depends only on the network, so
    it is built once and applied to every seed vector -- factors and null
    draws alike -- in one power iteration.
    """
    import networkx as nx
    from scipy import sparse

    simple = graph
    nodes = list(simple.nodes())
    index = {node: i for i, node in enumerate(nodes)}

    adjacency = nx.to_scipy_sparse_array(simple, nodelist=nodes, format="csr", weight=None)
    degree = np.asarray(adjacency.sum(axis=0)).ravel()
    degree[degree == 0] = 1.0
    transition = adjacency @ sparse.diags(1.0 / degree)

    def propagate(seed_matrix: np.ndarray, max_iter: int = 100, tol: float = 1e-8):
        """seed_matrix: nodes x k, column-normalised restart vectors."""
        scores = seed_matrix.copy()
        for _ in range(max_iter):
            updated = alpha * (transition @ scores) + (1 - alpha) * seed_matrix
            if np.abs(updated - scores).max() < tol:
                scores = updated
                break
            scores = updated
        return scores

    return nodes, index, propagate


def factor_network_propagation(
    graph,
    loadings: Dict[str, pd.DataFrame],
    id_maps: Optional[Dict[str, Dict[str, str]]] = None,
    top_n: int = 50,
    alpha: float = 0.85,
    n_permutations: int = 100,
    random_state: int = 0,
    exclude_edge_types: Sequence[str] = ("protein_interaction",),
    exclude_nodes: Optional[Sequence[str]] = None,
) -> Dict[str, object]:
    """Do a factor's features from different layers land in the same place?

    A joint factor loads on transcripts, proteins and metabolites at once, but
    that makes it trans-omic only if those features are related by the
    biochemistry -- if the transcripts it picks encode the enzymes acting on
    the metabolites it picks. The network says whether they are.

    Each layer's strongest features are seeded separately and diffused with a
    random walk with restart. If the layers' profiles overlap, the factor's
    features from different layers sit in the same region of the network: a
    candidate mechanism. If they do not, the factor is co-variation across
    layers with nothing joining it.

    The overlap is the mean cosine similarity between the layers' diffusion
    profiles. It is compared with the same statistic for random features of
    the same layers and the same number per layer -- a null that matters,
    because the transcriptome is ten times the size of the metabolome and
    uniform draws would be almost all transcripts.

    Parameters
    ----------
    graph : networkx.Graph
        The trans-omic network.
    loadings : dict
        ``{layer: features x factors}``, as :class:`Factorisation` holds.
    id_maps : dict, optional
        ``{layer: {feature id: node id}}`` when the tables are keyed by
        something other than the network's identifiers.
    top_n : int
        Seeds per layer per factor.
    alpha : float
        Restart parameter of the random walk.
    n_permutations : int
        Layer-matched random seed sets for the null. They are propagated
        together with the observed seeds, one factor at a time.
    random_state : int
    exclude_edge_types, exclude_nodes : sequence, optional
        Relationships and molecules diffusion routes around. By default STRING
        associations and currency metabolites, for the same reason path
        tracing excludes them: left in, every factor's signal ends up on water
        and ADP.

    Returns
    -------
    dict
        ``table`` : one row per factor -- seeds per layer, the observed
        overlap, its null mean, p and q.
        ``scores`` : nodes x factors, the combined diffusion profile.
        ``top_nodes`` : ``{factor: DataFrame}`` of the nodes receiving most
        signal *relative to layer-matched random seeds*, which is where to look
        for what the factor is about. Ranking by raw score instead returns the
        network's hubs for every factor alike.

    References
    ----------
    Cowen L, et al. Network propagation: a universal amplifier of genetic
    associations. *Nature Reviews Genetics* 18:551-562, 2017.
    """
    from itertools import combinations

    from transnet.biology.schema import CURRENCY_METABOLITES

    if exclude_nodes is None:
        exclude_nodes = CURRENCY_METABOLITES
    diffusion_graph = _propagation_graph(graph, exclude_edge_types, exclude_nodes)
    nodes, index, propagate = _propagator(diffusion_graph, alpha)
    layer_of = {node: str(graph.nodes[node].get("layer") or "Unknown") for node in nodes}
    by_layer: Dict[str, List[int]] = {}
    for node in nodes:
        by_layer.setdefault(layer_of[node], []).append(index[node])

    rng = np.random.default_rng(random_state)
    factors = list(next(iter(loadings.values())).columns)

    def normalise(matrix):
        totals = matrix.sum(axis=0)
        totals[totals == 0] = 1.0
        return matrix / totals

    def overlap(profiles):
        """Mean pairwise cosine similarity between layer profiles."""
        values = []
        for a, b in combinations(range(profiles.shape[1]), 2):
            x, y = profiles[:, a], profiles[:, b]
            norm = np.linalg.norm(x) * np.linalg.norm(y)
            values.append(float(x @ y / norm) if norm else 0.0)
        return float(np.mean(values)) if values else np.nan

    rows, score_columns, top_nodes = [], {}, {}
    for factor in factors:
        seeds = _mapped_top_features(loadings, factor, diffusion_graph, id_maps, top_n)
        layers = sorted(seeds)
        seed_counts = {layer: len(seeds[layer]) for layer in layers}
        row = {"factor": factor, **{f"seeds_{layer}": seed_counts[layer] for layer in layers}}

        if len(layers) < 2:
            rows.append({**row, "overlap": np.nan, "null_overlap": np.nan, "p_value": np.nan})
            continue

        # columns: [observed per layer] then [null draw k, per layer] for each k
        width = len(layers)
        matrix = np.zeros((len(nodes), width * (1 + n_permutations)))
        for j, layer in enumerate(layers):
            for node, weight in seeds[layer].items():
                matrix[index[node], j] = weight
        for k in range(n_permutations):
            for j, layer in enumerate(layers):
                pool = by_layer.get(layer, [])
                size = min(seed_counts[layer], len(pool))
                if size:
                    chosen = rng.choice(len(pool), size=size, replace=False)
                    matrix[[pool[c] for c in chosen], width * (1 + k) + j] = 1.0

        diffused = propagate(normalise(matrix))
        observed_profiles = diffused[:, :width]
        observed = overlap(observed_profiles)
        null = np.array([overlap(diffused[:, width * (1 + k): width * (2 + k)])
                         for k in range(n_permutations)])

        combined = observed_profiles.sum(axis=1)
        null_combined = np.mean(
            [diffused[:, width * (1 + k): width * (2 + k)].sum(axis=1)
             for k in range(n_permutations)], axis=0) if n_permutations else combined
        floor = np.median(null_combined[null_combined > 0]) if n_permutations else 1e-12
        enrichment = pd.Series(combined / np.maximum(null_combined, floor), index=nodes)
        series = pd.Series(combined, index=nodes)
        score_columns[factor] = series

        ranked = (enrichment[series > floor].sort_values(ascending=False).head(25)
                  .rename("enrichment").to_frame())
        ranked["score"] = series.reindex(ranked.index)
        ranked["layer"] = [layer_of.get(node, "Unknown") for node in ranked.index]
        ranked["name"] = [str(graph.nodes[node].get("symbol") or graph.nodes[node].get("name")
                              or node).split(";")[0][:48] for node in ranked.index]
        top_nodes[factor] = ranked

        rows.append({
            **row, "overlap": observed,
            "null_overlap": float(null.mean()) if len(null) else np.nan,
            "p_value": float((np.sum(null >= observed) + 1) / (len(null) + 1))
            if len(null) else np.nan,
        })

    table = pd.DataFrame(rows)
    if not table.empty and table["p_value"].notna().any():
        from statsmodels.stats.multitest import multipletests

        testable = table["p_value"].notna()
        table.loc[testable, "q_value"] = multipletests(
            table.loc[testable, "p_value"], method="fdr_bh")[1]
    scores_frame = pd.DataFrame(score_columns) if score_columns else pd.DataFrame()
    return {"table": table, "scores": scores_frame, "top_nodes": top_nodes}


def factor_layer_scores(factorisation: "Factorisation") -> Dict[str, pd.DataFrame]:
    """Each layer's own projection of the samples onto the joint factors."""
    return {
        layer: pd.DataFrame(
            factorisation.data_[layer].to_numpy() @ loadings.to_numpy(),
            index=factorisation.data_[layer].index, columns=loadings.columns,
        )
        for layer, loadings in factorisation.loadings_.items()
    }


def factor_cross_layer_agreement(factorisation: "Factorisation",
                                 method: str = "spearman") -> pd.DataFrame:
    """Do the layers agree about each factor?

    A joint factor is fitted on the layers stacked together, which does not
    make it trans-omic: one layer can carry it alone. Projecting the samples
    onto a factor *within each layer* and correlating those projections says
    which it is -- and a factor whose transcriptome and metabolome projections
    agree is the kind worth tracing on the network.

    Returns
    -------
    pandas.DataFrame
        One row per factor and layer pair: ``factor``, ``layer_a``,
        ``layer_b``, ``correlation``, ``p_value``.
    """
    from itertools import combinations

    from scipy import stats

    scores = factor_layer_scores(factorisation)
    rows = []
    for factor in factorisation.factors_.columns:
        for layer_a, layer_b in combinations(sorted(scores), 2):
            a = scores[layer_a][factor]
            b = scores[layer_b][factor].reindex(a.index)
            usable = a.notna() & b.notna()
            if usable.sum() < 3:
                continue
            if method == "pearson":
                correlation, p_value = stats.pearsonr(a[usable], b[usable])
            else:
                correlation, p_value = stats.spearmanr(a[usable], b[usable])
            rows.append({"factor": factor, "layer_a": layer_a, "layer_b": layer_b,
                         "correlation": float(correlation), "p_value": float(p_value)})
    return pd.DataFrame(rows)


def network_guided_imputation(frame: pd.DataFrame, graph,
                              id_map: Optional[Dict[str, str]] = None,
                              max_distance: int = 2) -> pd.DataFrame:
    """Fill missing values from a feature's network neighbours.

    The column mean says "this feature, on average"; the network says "this
    feature's neighbours, in this sample". For a metabolite missing in one
    animal, the substrates and products of the reactions it takes part in are
    a better guess than the cohort mean, and they are already in the graph.

    Features absent from the network, or with no measured neighbour, fall
    back to the column mean, and the fallback is logged rather than silent.
    """
    import networkx as nx

    filled = frame.copy()
    mapping = id_map or {}
    undirected = graph.to_undirected() if graph.is_directed() else graph
    fallbacks = 0

    for feature in frame.columns[frame.isna().any()]:
        node = mapping.get(str(feature), str(feature))
        neighbours = []
        if node in undirected:
            reachable = nx.single_source_shortest_path_length(
                undirected, node, cutoff=max_distance)
            inverse = {v: k for k, v in mapping.items()}
            for neighbour in reachable:
                if neighbour == node:
                    continue
                column = inverse.get(neighbour, neighbour)
                if column in frame.columns and frame[column].notna().any():
                    neighbours.append(column)

        missing = frame[feature].isna()
        if neighbours:
            # neighbours are on their own scales; standardise before averaging
            block = frame[neighbours].apply(
                lambda column: (column - column.mean()) / (column.std(ddof=0) or 1.0))
            predicted = block.mean(axis=1)
            scale = frame[feature].std(ddof=0) or 1.0
            filled.loc[missing, feature] = (
                frame[feature].mean() + predicted[missing] * scale)
        else:
            fallbacks += 1
            filled.loc[missing, feature] = frame[feature].mean()

    if fallbacks:
        logger.info(f"Network imputation: {fallbacks} feature(s) had no measured "
                    f"neighbour within {max_distance} steps and used the column mean")
    return filled


__all__ += ["factor_network_propagation", "factor_layer_scores",
            "factor_cross_layer_agreement", "network_guided_imputation"]


def sample_pairing_check(first: pd.DataFrame, second: pd.DataFrame,
                         groups: pd.Series, n_permutations: int = 200,
                         min_samples: int = 4, random_state: int = 0) -> Dict[str, object]:
    """Are the layers really measured on the same samples?

    A joint factor model treats column ``i`` of every layer as one sample.
    When the files merely share column *names*, that is an assumption, and the
    data can test it where it bites: within a group of replicates, a sample
    whose transcript of a gene is unusually high should, if it is the same
    sample, tend to have that protein unusually high too.

    So the statistic is computed **per gene, across the samples of a group**:
    the Spearman correlation between a gene's transcript and protein over the
    replicates, averaged over genes and groups. Shuffling which replicate is
    paired with which destroys it only if the pairing was real. Correlating
    across genes *within* a sample instead -- the obvious first attempt --
    measures gene abundance, which is the same whatever the pairing, and says
    nothing at all when values are centred per feature, as MoTrPAC's are.

    Parameters
    ----------
    first, second : pandas.DataFrame
        samples x features, with the features already matched between them
        (the same gene in both) -- see :func:`matched_transcript_protein`.
    groups : pandas.Series
        Group label per sample; shuffling happens only within a group, so a
        treatment effect cannot pass for pairing.
    n_permutations : int
    min_samples : int
        Replicates a group needs to contribute.

    Returns
    -------
    dict
        ``observed`` mean per-gene within-group correlation as listed,
        ``null_mean`` under shuffled pairings, ``p_value``, ``n_genes``, and
        ``paired``: whether the check supports treating the columns as the
        same samples.
    """
    rng = np.random.default_rng(random_state)
    features = [f for f in first.columns if f in set(second.columns)]

    def ranked(values):
        """Column-wise ranks, NaN-safe, centred and scaled for fast Spearman."""
        frame = pd.DataFrame(values).rank(axis=0)
        frame = frame - frame.mean(axis=0)
        norm = np.sqrt((frame ** 2).sum(axis=0)).replace(0, np.nan)
        return (frame / norm).to_numpy()

    blocks = []
    for _, members in groups.groupby(groups):
        samples = [s for s in members.index if s in first.index and s in second.index]
        if len(samples) < min_samples:
            continue
        a = first.loc[samples, features].to_numpy(float)
        b = second.loc[samples, features].to_numpy(float)
        usable = ~(np.isnan(a).any(axis=0) | np.isnan(b).any(axis=0))
        usable &= (np.nanstd(a, axis=0) > 0) & (np.nanstd(b, axis=0) > 0)
        if usable.sum() == 0:
            continue
        blocks.append((a[:, usable], b[:, usable]))

    if not blocks:
        return {"observed": np.nan, "null_mean": np.nan, "p_value": np.nan,
                "n_genes": 0, "paired": False}

    def statistic(orders):
        values = []
        for (a, b), order in zip(blocks, orders):
            ra, rb = ranked(a), ranked(b[order])
            values.append(np.nanmean(np.nansum(ra * rb, axis=0)))
        return float(np.mean(values))

    observed = statistic([np.arange(len(a)) for a, _ in blocks])
    null = np.array([statistic([rng.permutation(len(a)) for a, _ in blocks])
                     for _ in range(n_permutations)])
    p_value = float((np.sum(null >= observed) + 1) / (len(null) + 1))
    return {"observed": observed, "null_mean": float(null.mean()), "p_value": p_value,
            "n_genes": int(max(a.shape[1] for a, _ in blocks)),
            "paired": bool(p_value <= 0.05)}


__all__ += ["sample_pairing_check"]


def pairing_robustness(factorisation: "Factorisation", groups: pd.Series,
                       n_shuffles: int = 200, random_state: int = 0) -> pd.DataFrame:
    """Is a factor carried by the design, or by variation between samples?

    A factor on which the layers agree can owe that agreement to two different
    things, and they need different evidence. When :func:`sample_pairing_check`
    cannot confirm that the layers share samples, only the first is safe; when
    it can, both are real. A joint factor is usable regardless of the pairing
    *if the layers agree about it whatever the pairing*. That is the case for a factor carried by
    differences between groups -- timepoints, genotypes -- because the group
    labels are known: every replicate of a timepoint scores alike, so it does
    not matter which replicate is paired with which. A factor that lives in
    replicate-to-replicate variation is joint only if the pairing is real.

    The factor's loadings are held fixed; each layer's samples are projected
    onto them, the rows of every layer but the first are shuffled within each
    group, and the agreement between the layers' projections is recomputed.

    Returns
    -------
    pandas.DataFrame
        One row per factor: ``agreement`` (mean cross-layer Spearman of the
        per-layer projections as listed), ``shuffled_agreement`` (mean under
        within-group shuffles), ``retained`` (their ratio) and ``verdict``:

        ``between-group``
            the layers agree just as well under any within-group pairing, so
            the factor is carried by differences between groups -- timepoints,
            genotypes -- whose labels are certain. Interpretable whatever the
            pairing.
        ``within-group``
            the agreement collapses when replicates are re-paired, so the
            factor lives in sample-to-sample variation shared across layers.
            Real and often the interesting part -- individual animals that
            respond more strongly in every layer -- but only if the samples
            really are paired; :func:`sample_pairing_check` decides that.
        ``single-layer``
            the layers do not agree about it even as listed: one layer's view
            of the design wearing a joint label.
    """
    from itertools import combinations

    from scipy import stats

    rng = np.random.default_rng(random_state)
    projections = factor_layer_scores(factorisation)
    layers = sorted(projections)
    samples = list(projections[layers[0]].index)
    groups = groups.reindex(samples)

    def agreement(frames, factor):
        values = []
        for a, b in combinations(layers, 2):
            x, y = frames[a][factor].to_numpy(), frames[b][factor].to_numpy()
            if np.std(x) > 0 and np.std(y) > 0:
                values.append(stats.spearmanr(x, y).statistic)
        return float(np.mean(values)) if values else np.nan

    positions = {group: [samples.index(s) for s in members.index]
                 for group, members in groups.groupby(groups)}
    rows = []
    for factor in factorisation.factors_.columns:
        observed = agreement(projections, factor)
        shuffled = []
        for _ in range(n_shuffles):
            frames = {layers[0]: projections[layers[0]]}
            for layer in layers[1:]:
                order = np.arange(len(samples))
                for members in positions.values():
                    order[members] = rng.permutation(members)
                frame = projections[layer].iloc[order].copy()
                frame.index = samples
                frames[layer] = frame
            shuffled.append(agreement(frames, factor))
        shuffled_mean = float(np.nanmean(shuffled))
        retained = shuffled_mean / observed if observed and observed > 0 else np.nan
        if not observed or observed < 0.3:
            verdict = "single-layer"
        elif retained >= 0.8:
            verdict = "between-group"
        else:
            verdict = "within-group"
        rows.append({"factor": factor, "agreement": observed,
                     "shuffled_agreement": shuffled_mean, "retained": retained,
                     "verdict": verdict})
    return pd.DataFrame(rows)


__all__ += ["pairing_robustness"]


def matched_transcript_protein(transcriptome: pd.DataFrame, proteome: pd.DataFrame,
                               graph, id_maps: Optional[Dict[str, Dict[str, str]]] = None
                               ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Transcript and protein matrices for the same genes, keyed by gene.

    The pairing check needs each gene in both layers. Transcripts and proteins
    are keyed differently (Ensembl or symbols against UniProt), and the network
    already knows which transcript encodes which protein: its ``translation``
    edges. Both inputs are features x samples; both outputs are samples x
    genes with identical columns.
    """
    maps = id_maps or {}
    transcript_map = maps.get("Transcriptome", {})
    protein_map = maps.get("Proteome", {})
    protein_node_to_feature = {protein_map.get(str(f), str(f)): f for f in proteome.index}

    translation = {u: v for u, v, d in graph.edges(data=True)
                   if d.get("edge_type") == "translation"}
    pairs = {}
    for feature in transcriptome.index:
        gene = transcript_map.get(str(feature), str(feature))
        protein = translation.get(gene)
        if protein in protein_node_to_feature:
            pairs[feature] = protein_node_to_feature[protein]

    samples = [s for s in transcriptome.columns if s in set(proteome.columns)]
    transcripts = transcriptome.loc[list(pairs), samples].T
    proteins = proteome.loc[list(pairs.values()), samples].T
    proteins.columns = transcripts.columns
    return transcripts, proteins


__all__ += ["matched_transcript_protein"]
