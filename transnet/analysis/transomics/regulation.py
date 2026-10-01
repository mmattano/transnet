"""Which axis regulates each reaction: enzyme amount, metabolites, or both.

A reaction can change through the *gene-expression axis* (the amount of its
enzyme, set by transcription and translation) or through the *metabolite axis*
(its substrates, products and allosteric regulators). Reactions where the two
point in opposite directions are flagged as ``controversial``.

Both axes use whatever layers the network has. ``gene_axis_evidence`` records
which layers supported each enzyme-axis call, so a network with fewer layers
gives visibly weaker calls rather than wrong ones.
"""

from typing import Dict, Optional, Sequence

import logging

import numpy as np
import pandas as pd

from transnet.biology.schema import (
    CURRENCY_METABOLITES, available_edge_types, available_layers,
)

logger = logging.getLogger(__name__)

__all__ = [
    "reaction_regulation_table",
    "regulation_axis_summary",
    "metabolite_regulatory_roles",
    "regulatory_role_enrichment",
]


#: Edge types that regulate a reaction through metabolite concentration.
#: Substrate increase pushes a reaction forward; product increase pushes back.
_METABOLITE_EFFECT = {
    "allosteric_activation": 1,
    "allosteric_inhibition": -1,
    "substrate": 1,
    "product": -1,
}


def _sign_label(value: int) -> str:
    return {1: "activated", -1: "inhibited", 0: "unchanged"}[int(value)]


def _incoming(graph, node):
    """Predecessor edges of ``node``, working for directed and undirected graphs."""
    if graph.is_directed():
        if graph.is_multigraph():
            return [(u, d) for u, _, d in graph.in_edges(node, data=True)]
        return [(u, d) for u, _, d in graph.in_edges(node, data=True)]
    if graph.is_multigraph():
        return [(u, d) for u, _, d in graph.edges(node, data=True)]
    return [(u, d) for u, _, d in graph.edges(node, data=True)]


def _outgoing(graph, node):
    if graph.is_directed():
        return [(v, d) for _, v, d in graph.out_edges(node, data=True)]
    if graph.is_multigraph():
        return [(v, d) for _, v, d in graph.edges(node, data=True)]
    return [(v, d) for _, v, d in graph.edges(node, data=True)]


def _enzyme_phospho_axis(graph, enzyme) -> Dict[str, object]:
    """Phosphorylation feeding one enzyme, kept apart from the gene axis.

    A ``phosphorylation`` edge from a changed Signaling node says the enzyme's
    modification state moved. Whether that raises or lowers catalytic activity is
    site-specific and usually unrecorded: the schema gives ``phosphorylation`` a
    sign of 0, and KEGG only signs a relation when it annotates it as activation
    or inhibition. So two things are reported separately:

    ``direction``
        Which way the site itself moved, from the Signaling node's ``regulated``.
    ``effect``
        What that does to the reaction, which is ``direction * sign`` where the
        edge carries a sign, and 0 (unknown) where it does not.

    Folding ``effect`` into the gene axis would assert a direction the data does
    not have, so :func:`reaction_regulation_table` keeps it in its own columns.
    """
    directions, effects, via = [], set(), []
    for source, data in _incoming(graph, enzyme):
        if data.get("edge_type") != "phosphorylation":
            continue
        state = graph.nodes[source].get("regulated", 0) or 0
        if not state:
            continue
        via.append(source)
        directions.append(int(state))
        try:
            sign = int(data.get("sign") or 0)
        except (TypeError, ValueError):
            sign = 0
        effects.add(int(state) * sign)

    if not via:
        return {"direction": 0, "effect": 0, "via": [], "n_sites": 0}

    signed = {e for e in effects if e}
    distinct = set(directions)
    return {
        # Sites on one enzyme can move opposite ways, and then no single
        # direction is honest. ``n_sites`` is what separates that from no
        # site having changed at all.
        "direction": distinct.pop() if len(distinct) == 1 else 0,
        "effect": signed.pop() if len(signed) == 1 else 0,
        "via": sorted(set(via)),
        "n_sites": len(via),
    }


def _enzyme_gene_axis(graph, enzyme) -> Dict[str, object]:
    """Direction of the gene-expression axis feeding one enzyme.

    Prefers the enzyme's own measured change; falls back to the transcript that
    encodes it, and then to the transcription factors regulating that
    transcript.  Returns the direction plus the evidence level it came from.
    """
    # gene -> protein (translation)
    genes = [
        u for u, data in _incoming(graph, enzyme)
        if data.get("edge_type") == "translation"
    ]
    gene_states = [
        (g, graph.nodes[g].get("regulated", 0) or 0) for g in genes
    ]
    measured_genes = [
        g for g in genes
        if graph.nodes[g].get("measured")
        or graph.nodes[g].get("log2fc") is not None
        or graph.nodes[g].get("qvalue") is not None
    ]

    protein_state = graph.nodes[enzyme].get("regulated", 0) or 0
    if protein_state:
        # ``evidence`` names the strongest single source, so a changed protein
        # reads "protein" whether or not its transcript moved too. Whether it
        # did is a separate question -- transcriptional or post-transcriptional
        # regulation -- answered here.
        support = (
            any(s == protein_state for _, s in gene_states)
            if measured_genes else None
        )
        return {"direction": int(protein_state), "evidence": "protein",
                "via": [enzyme], "transcript_support": support}
    regulated_genes = [(g, s) for g, s in gene_states if s]
    if regulated_genes:
        directions = {s for _, s in regulated_genes}
        direction = regulated_genes[0][1] if len(directions) == 1 else 0
        return {
            "direction": int(direction),
            "evidence": "gene_protein",
            "via": [g for g, _ in regulated_genes],
            "transcript_support": True,
        }

    # TF -> gene -> protein
    tf_hits = []
    for gene in genes:
        for tf, data in _incoming(graph, gene):
            if data.get("edge_type") != "transcriptional_regulation":
                continue
            state = graph.nodes[tf].get("regulated", 0) or 0
            if state:
                tf_hits.append((tf, state))
    if tf_hits:
        directions = {s for _, s in tf_hits}
        # ChIP binding does not tell us whether a factor activates or represses,
        # so a TF-only chain gives presence of regulation, not its direction.
        direction = tf_hits[0][1] if len(directions) == 1 else 0
        return {
            "direction": int(direction),
            "evidence": "tf_gene_protein",
            "via": [t for t, _ in tf_hits],
        }

    return {"direction": 0, "evidence": None, "via": []}


def reaction_regulation_table(
    graph,
    reactions: Optional[Sequence[str]] = None,
    include_substrate_product: bool = True,
    exclude_currency_metabolites: bool = True,
    currency_metabolites: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """Attribute each reaction's regulation to the gene and metabolite axes.

    Parameters
    ----------
    graph : networkx.Graph
        A trans-omic network with omics data mapped onto it.
    reactions : sequence of str, optional
        Reactions to include. Default: the whole Reactions layer.
    include_substrate_product : bool
        Count substrate and product changes in the metabolite axis. False
        restricts the axis to allosteric regulators.
    exclude_currency_metabolites : bool
        Ignore cofactors such as ATP, NADH and water as substrates and products
        (default), since they take part in hundreds of reactions. Their allosteric
        roles are always kept.
    currency_metabolites : sequence of str, optional
        Replaces the default set, ``transnet.biology.schema.CURRENCY_METABOLITES``.

    Returns
    -------
    pandas.DataFrame
        One row per reaction:

        ``reaction``, ``name``, ``ec``, ``reversible``
            The reaction.
        ``gene_axis``
            Direction of enzyme-level regulation: +1, -1 or 0.
        ``gene_axis_evidence``
            Layers that supported it: ``"protein"``, ``"gene_protein"``,
            ``"gene"``, ``"tf_gene_protein"`` or None.
        ``gene_axis_via``
            The molecules behind the call.
        ``gene_axis_transcript_support``
            True if a transcript of the enzyme changed the same way, False if the
            enzyme changed and its measured transcript did not, None if no
            transcript was measured.
        ``phospho_axis``, ``phospho_axis_effect``, ``n_phosphosites_changed``, ``phospho_axis_via``
            Phosphorylation of the enzyme, when sites are mapped:
            the direction the sites moved, its effect on the reaction (0 unless
            the edge has a sign), and how many sites changed.
        ``metabolite_axis``, ``n_activators``, ``n_inhibitors``
            Direction of metabolite regulation and the counts behind it.
        ``allosteric_regulators``, ``substrates_changed``, ``products_changed``
            The metabolites behind the metabolite axis.
        ``consensus``
            Overall direction when the axes agree or only one is informative.
        ``controversial``
            True when the two axes point in opposite directions.
    """
    columns = [
        "reaction", "name", "ec", "reversible",
        "gene_axis", "gene_axis_evidence", "gene_axis_via",
        "gene_axis_transcript_support",
        "phospho_axis", "phospho_axis_effect", "n_phosphosites_changed",
        "phospho_axis_via",
        "metabolite_axis", "n_activators", "n_inhibitors",
        "allosteric_regulators", "substrates_changed", "products_changed",
        "consensus", "controversial",
    ]

    if "Reactions" not in available_layers(graph):
        logger.warning(
            "No Reactions layer in this network; reaction regulation cannot be "
            "attributed. Present layers: %s", available_layers(graph)
        )
        return pd.DataFrame(columns=columns)

    edge_types = available_edge_types(graph)
    if not any(t in edge_types for t in _METABOLITE_EFFECT):
        logger.warning(
            "Network has a Reactions layer but no substrate/product/allosteric "
            "edges; the metabolite axis will be empty for every reaction."
        )

    if reactions is None:
        reactions = [
            n for n, d in graph.nodes(data=True) if d.get("layer") == "Reactions"
        ]

    effect_types = dict(_METABOLITE_EFFECT)
    if not include_substrate_product:
        effect_types.pop("substrate", None)
        effect_types.pop("product", None)

    if currency_metabolites is not None:
        currency = set(currency_metabolites)
    elif exclude_currency_metabolites:
        currency = set(CURRENCY_METABOLITES)
    else:
        currency = set()

    #: Mass-action roles only; an allosteric edge is never suppressed.
    mass_action = {"substrate", "product"}

    rows = []
    for reaction in reactions:
        attrs = graph.nodes[reaction]

        # --- gene-expression axis -------------------------------------------
        catalysts = _incoming(graph, reaction)
        enzymes = [
            u for u, data in catalysts if data.get("edge_type") == "catalysis"
        ]
        # Networks built without a Proteome link genes straight to reactions;
        # the transcript is then the only enzyme-level evidence available.
        genes_as_catalysts = [
            u for u, data in catalysts if data.get("edge_type") == "gene_catalysis"
        ]

        axis_results = [_enzyme_gene_axis(graph, e) for e in enzymes]
        phospho_results = [_enzyme_phospho_axis(graph, e) for e in enzymes]
        for gene in genes_as_catalysts:
            state = graph.nodes[gene].get("regulated", 0) or 0
            axis_results.append({
                "direction": int(state),
                "evidence": "gene" if state else None,
                "via": [gene] if state else [],
                "transcript_support": True if state else None,
            })
        informative = [r for r in axis_results if r["direction"] != 0]

        if informative:
            directions = {r["direction"] for r in informative}
            gene_axis = informative[0]["direction"] if len(directions) == 1 else 0
            # Report the strongest evidence level that contributed.
            order = ["protein", "gene_protein", "gene", "tf_gene_protein"]
            evidence = sorted(
                {r["evidence"] for r in informative if r["evidence"]},
                key=order.index,
            )[0]
            via = sorted({m for r in informative for m in r["via"]})
            agreeing = [r for r in informative if r["direction"] == gene_axis]
            supports = [r.get("transcript_support") for r in agreeing]
            if any(s is True for s in supports):
                transcript_support = True
            elif any(s is False for s in supports):
                transcript_support = False
            else:
                transcript_support = None
        else:
            gene_axis = 0
            evidence = None
            via = []
            transcript_support = None
            # Even with no change, record how far the chain reaches, so a
            # missing layer is distinguishable from a measured non-response.
            reachable = [r["evidence"] for r in axis_results if r["evidence"]]
            if enzymes and not reachable:
                evidence = None

        # --- phosphorylation, reported beside the gene axis rather than in it -
        phospho_directions = {r["direction"] for r in phospho_results if r["direction"]}
        phospho_effects = {r["effect"] for r in phospho_results if r["effect"]}
        phospho_via = sorted({m for r in phospho_results for m in r["via"]})
        phospho_axis = phospho_directions.pop() if len(phospho_directions) == 1 else 0
        phospho_effect = phospho_effects.pop() if len(phospho_effects) == 1 else 0
        n_phosphosites = sum(r["n_sites"] for r in phospho_results)

        # --- metabolite axis -------------------------------------------------
        activators, inhibitors = [], []
        allosteric, subs_changed, prods_changed = [], [], []

        for neighbour, data in _incoming(graph, reaction):
            etype = data.get("edge_type")
            if etype not in effect_types:
                continue
            if etype in mass_action and str(neighbour) in currency:
                continue
            state = graph.nodes[neighbour].get("regulated", 0) or 0
            if not state:
                continue
            effect = effect_types[etype] * state
            if etype in ("allosteric_activation", "allosteric_inhibition"):
                allosteric.append(f"{neighbour}({_sign_label(effect)})")
            else:
                subs_changed.append(neighbour)
            (activators if effect > 0 else inhibitors).append(neighbour)

        # Products leave the reaction, so they are outgoing edges.
        for neighbour, data in _outgoing(graph, reaction):
            if data.get("edge_type") != "product" or "product" not in effect_types:
                continue
            if str(neighbour) in currency:
                continue
            state = graph.nodes[neighbour].get("regulated", 0) or 0
            if not state:
                continue
            effect = effect_types["product"] * state
            prods_changed.append(neighbour)
            (activators if effect > 0 else inhibitors).append(neighbour)

        if activators and not inhibitors:
            metabolite_axis = 1
        elif inhibitors and not activators:
            metabolite_axis = -1
        elif activators and inhibitors:
            metabolite_axis = int(np.sign(len(activators) - len(inhibitors)))
        else:
            metabolite_axis = 0

        # --- consensus -------------------------------------------------------
        controversial = bool(
            gene_axis and metabolite_axis and gene_axis != metabolite_axis
        )
        if controversial:
            consensus = 0
        else:
            consensus = gene_axis or metabolite_axis

        rows.append({
            "reaction": reaction,
            "name": attrs.get("name"),
            "ec": attrs.get("ec", ""),
            "reversible": attrs.get("reversible", True),
            "gene_axis": gene_axis,
            "gene_axis_evidence": evidence,
            "gene_axis_via": ";".join(via),
            "gene_axis_transcript_support": transcript_support,
            "phospho_axis": phospho_axis,
            "phospho_axis_effect": phospho_effect,
            "n_phosphosites_changed": n_phosphosites,
            "phospho_axis_via": ";".join(phospho_via),
            "metabolite_axis": metabolite_axis,
            "n_activators": len(activators),
            "n_inhibitors": len(inhibitors),
            "allosteric_regulators": ";".join(allosteric),
            "substrates_changed": ";".join(subs_changed),
            "products_changed": ";".join(prods_changed),
            "consensus": consensus,
            "controversial": controversial,
        })

    table = pd.DataFrame(rows, columns=columns)

    if not table.empty:
        responsive = table[(table["gene_axis"] != 0) | (table["metabolite_axis"] != 0)]
        logger.info(
            f"{len(responsive)}/{len(table)} reactions differentially regulated; "
            f"{int(table['controversial'].sum())} controversial "
            f"({table['controversial'].mean():.0%} of all reactions)"
        )

    return table


def regulation_axis_summary(
    table: pd.DataFrame,
    pathway_map: Optional[Dict[str, str]] = None,
    responsive_only: bool = True,
) -> pd.DataFrame:
    """Summarise a regulation table per pathway.

    Parameters
    ----------
    table : pandas.DataFrame
        Output of :func:`reaction_regulation_table`.
    pathway_map : dict, optional
        ``{reaction: pathway}`` or ``{reaction: [pathway, ...]}``, for example
        from :func:`transnet.api.kegg_reaction_pathways`. A reaction in several
        pathways is counted in each. Without a map all reactions form one
        ``"all"`` row.
    responsive_only : bool
        Count only reactions regulated on at least one axis (default).

    Returns
    -------
    pandas.DataFrame
        One row per pathway: reaction count, reactions activated and inhibited
        on each axis, and the number and fraction of controversial reactions.
    """
    columns = [
        "pathway", "n_reactions",
        "gene_activated", "gene_inhibited",
        "metabolite_activated", "metabolite_inhibited",
        "n_controversial", "fraction_controversial",
        "fraction_gene_driven", "fraction_metabolite_driven",
    ]
    if table.empty:
        return pd.DataFrame(columns=columns)

    data = table.copy()
    if responsive_only:
        data = data[(data["gene_axis"] != 0) | (data["metabolite_axis"] != 0)]
    if data.empty:
        return pd.DataFrame(columns=columns)

    if pathway_map:
        def _pathways(reaction):
            value = pathway_map.get(reaction)
            if value is None or (isinstance(value, float) and np.isnan(value)):
                return ["unassigned"]
            return [value] if isinstance(value, str) else list(value) or ["unassigned"]
        data["pathway"] = data["reaction"].map(_pathways)
        data = data.explode("pathway")
    else:
        data["pathway"] = "all"

    rows = []
    for pathway, group in data.groupby("pathway"):
        n = len(group)
        gene_driven = int(((group["gene_axis"] != 0)).sum())
        met_driven = int(((group["metabolite_axis"] != 0)).sum())
        rows.append({
            "pathway": pathway,
            "n_reactions": n,
            "gene_activated": int((group["gene_axis"] == 1).sum()),
            "gene_inhibited": int((group["gene_axis"] == -1).sum()),
            "metabolite_activated": int((group["metabolite_axis"] == 1).sum()),
            "metabolite_inhibited": int((group["metabolite_axis"] == -1).sum()),
            "n_controversial": int(group["controversial"].sum()),
            "fraction_controversial": float(group["controversial"].mean()),
            "fraction_gene_driven": gene_driven / n,
            "fraction_metabolite_driven": met_driven / n,
        })

    return (
        pd.DataFrame(rows, columns=columns)
        .sort_values("n_reactions", ascending=False)
        .reset_index(drop=True)
    )


def metabolite_regulatory_roles(
    graph,
    metabolites: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """Classify each metabolite as an allosteric regulator or not.

    The allosteric edges come from BRENDA. A changed metabolite that regulates an
    enzyme is a candidate mechanism; one that is only a substrate or product often
    just follows the flux.

    Parameters
    ----------
    graph : networkx.Graph
        A trans-omic network, ideally with omics data mapped onto it.
    metabolites : sequence of str, optional
        Metabolites to include. Default: the whole Metabolome layer.

    Returns
    -------
    pandas.DataFrame
        One row per metabolite: ``metabolite``, ``name``; the measured
        ``regulated``, ``log2fc``, ``qvalue``; ``role`` (``"activator"``,
        ``"inhibitor"``, ``"both"``, ``"substrate/product only"`` or ``"none"``);
        ``is_allosteric_regulator``; the number and ids of the reactions it
        activates and inhibits; ``n_enzymes_regulated``; and ``is_substrate``,
        ``is_product``.
    """
    columns = [
        "metabolite", "name", "regulated", "log2fc", "qvalue",
        "role", "is_allosteric_regulator",
        "n_reactions_activated", "n_reactions_inhibited",
        "n_reactions_regulated", "n_enzymes_regulated",
        "reactions_activated", "reactions_inhibited",
        "is_substrate", "is_product",
    ]

    if "Metabolome" not in available_layers(graph):
        logger.warning(
            "No Metabolome layer in this network; there are no metabolite "
            "roles to report. Present layers: %s", available_layers(graph)
        )
        return pd.DataFrame(columns=columns)

    edge_types = available_edge_types(graph)
    if not any(t.startswith("allosteric") for t in edge_types):
        logger.warning(
            "This network has no allosteric edges, so every metabolite will be "
            "reported as a non-regulator. Populate BRENDA effectors with "
            "Proteome.get_brenda_kinetics() before building the network."
        )

    if metabolites is None:
        metabolites = [
            n for n, d in graph.nodes(data=True) if d.get("layer") == "Metabolome"
        ]

    rows = []
    for metabolite in metabolites:
        if metabolite not in graph:
            continue
        attrs = graph.nodes[metabolite]

        activated, inhibited, enzymes = [], [], set()
        is_substrate = is_product = False

        for target, data in _outgoing(graph, metabolite):
            edge_type = data.get("edge_type")
            if edge_type == "allosteric_activation":
                activated.append(target)
            elif edge_type == "allosteric_inhibition":
                inhibited.append(target)
            elif edge_type == "substrate":
                is_substrate = True
                continue
            else:
                continue
            for ec in str(data.get("ec", "") or "").split(","):
                if ec:
                    enzymes.add(ec)

        # Products point reaction -> metabolite, so they arrive as in-edges.
        for _, data in _incoming(graph, metabolite):
            if data.get("edge_type") == "product":
                is_product = True
                break

        if activated and inhibited:
            role = "both"
        elif activated:
            role = "activator"
        elif inhibited:
            role = "inhibitor"
        elif is_substrate or is_product:
            role = "substrate/product only"
        else:
            role = "none"

        rows.append({
            "metabolite": metabolite,
            "name": attrs.get("name"),
            "regulated": attrs.get("regulated", 0) or 0,
            "log2fc": attrs.get("log2fc"),
            "qvalue": attrs.get("qvalue"),
            "role": role,
            "is_allosteric_regulator": bool(activated or inhibited),
            "n_reactions_activated": len(set(activated)),
            "n_reactions_inhibited": len(set(inhibited)),
            "n_reactions_regulated": len(set(activated) | set(inhibited)),
            "n_enzymes_regulated": len(enzymes),
            "reactions_activated": ";".join(sorted(set(activated))),
            "reactions_inhibited": ";".join(sorted(set(inhibited))),
            "is_substrate": is_substrate,
            "is_product": is_product,
        })

    table = pd.DataFrame(rows, columns=columns)

    if not table.empty:
        regulators = int(table["is_allosteric_regulator"].sum())
        logger.info(
            f"{regulators}/{len(table)} metabolites act as allosteric regulators"
        )
    return table


def regulatory_role_enrichment(
    roles: pd.DataFrame,
    background: Optional[Sequence[str]] = None,
) -> Dict[str, object]:
    """Test whether changed metabolites are regulators more often than expected.

    Fisher's exact test of changed (yes/no) against regulator (yes/no), overall
    and per role, with Benjamini-Hochberg correction.

    Parameters
    ----------
    roles : pandas.DataFrame
        Output of :func:`metabolite_regulatory_roles`.
    background : sequence of str, optional
        The measured metabolites to compare against. Default: every metabolite
        in ``roles`` with a measurement.

    Returns
    -------
    dict
        ``counts``: numbers of background, changed, and changed-and-regulator
        metabolites, with the fractions. ``enrichment``: one row per role
        (``any``, ``activator``, ``inhibitor``) with odds ratio, p- and q-value.
        ``regulators``: the changed metabolites that are regulators.
    """
    empty_enrichment = pd.DataFrame(columns=[
        "role", "n_differential_with_role", "n_differential_without_role",
        "n_background_with_role", "n_background_without_role",
        "odds_ratio", "p_value", "q_value",
    ])

    if roles.empty:
        return {
            "counts": {}, "enrichment": empty_enrichment,
            "regulators": pd.DataFrame(columns=roles.columns),
        }

    measured = roles[roles["log2fc"].notna() | roles["qvalue"].notna()]
    if background is not None:
        pool = roles[roles["metabolite"].isin(set(background))]
    elif not measured.empty:
        pool = measured
    else:
        logger.warning(
            "No metabolite carries a measurement; using every metabolite in the "
            "network as the background, which will overstate the background size."
        )
        pool = roles

    differential = pool[pool["regulated"] != 0]
    non_differential = pool[pool["regulated"] == 0]

    role_masks = {
        "any": lambda frame: frame["is_allosteric_regulator"],
        "activator": lambda frame: frame["role"].isin(["activator", "both"]),
        "inhibitor": lambda frame: frame["role"].isin(["inhibitor", "both"]),
    }

    from scipy.stats import fisher_exact

    rows = []
    for role, mask in role_masks.items():
        a = int(mask(differential).sum())              # differential, has role
        b = int(len(differential) - a)                 # differential, no role
        c = int(mask(non_differential).sum())          # other, has role
        d = int(len(non_differential) - c)             # other, no role

        if (a + b) == 0 or (c + d) == 0:
            odds, p_value = float("nan"), float("nan")
        else:
            odds, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")

        rows.append({
            "role": role,
            "n_differential_with_role": a,
            "n_differential_without_role": b,
            "n_background_with_role": c,
            "n_background_without_role": d,
            "odds_ratio": float(odds),
            "p_value": float(p_value),
        })

    enrichment = pd.DataFrame(rows)
    valid = enrichment["p_value"].notna()
    enrichment["q_value"] = np.nan
    if valid.any():
        ranked = enrichment.loc[valid].sort_values("p_value")
        n = len(ranked)
        raw = ranked["p_value"].to_numpy() * n / np.arange(1, n + 1)
        ranked = ranked.assign(
            q_value=np.minimum.accumulate(raw[::-1])[::-1].clip(max=1.0)
        )
        enrichment.loc[ranked.index, "q_value"] = ranked["q_value"].to_numpy()

    n_differential = len(differential)
    n_differential_regulators = int(differential["is_allosteric_regulator"].sum())
    n_pool_regulators = int(pool["is_allosteric_regulator"].sum())

    counts = {
        "n_background": len(pool),
        "n_differential": n_differential,
        "n_differential_regulators": n_differential_regulators,
        "fraction_differential_regulators": (
            n_differential_regulators / n_differential if n_differential else 0.0
        ),
        "n_background_regulators": n_pool_regulators,
        "fraction_background_regulators": (
            n_pool_regulators / len(pool) if len(pool) else 0.0
        ),
    }

    logger.info(
        f"{n_differential_regulators}/{n_differential} differential metabolites "
        f"act as allosteric regulators "
        f"({counts['fraction_differential_regulators']:.0%}; background "
        f"{counts['fraction_background_regulators']:.0%})"
    )

    regulators = differential[differential["is_allosteric_regulator"]].sort_values(
        "n_reactions_regulated", ascending=False
    ).reset_index(drop=True)

    return {
        "counts": counts,
        "enrichment": enrichment[[
            "role", "n_differential_with_role", "n_differential_without_role",
            "n_background_with_role", "n_background_without_role",
            "odds_ratio", "p_value", "q_value",
        ]],
        "regulators": regulators,
    }
