"""The structure of the responsive network: motifs, weak points, a null model.

* :func:`regulatory_motifs` finds signed wiring patterns such as product
  inhibition.
* :func:`structural_vulnerability` finds the molecules whose removal would
  split the response into separate pieces.
* :func:`convergence_significance` tests whether the layers converge on
  shared reactions more than chance would produce.
"""

from typing import Dict, List, Optional, Sequence

import logging

import networkx as nx
import numpy as np
import pandas as pd

from transnet.analysis.transomics.mapping import is_responsive
from transnet.biology.schema import CURRENCY_METABOLITES, order_layers

logger = logging.getLogger(__name__)

__all__ = [
    "regulatory_motifs",
    "structural_vulnerability",
    "convergence_significance",
]

#: Relationships that carry a regulatory claim rather than mass flow.
_REGULATORY = frozenset({
    "allosteric_activation", "allosteric_inhibition", "transcriptional_regulation",
    "phosphorylation", "kinase_tf",
})


def _edge_records(graph, u, v) -> List[dict]:
    data = graph.get_edge_data(u, v)
    if data is None:
        return []
    return list(data.values()) if graph.is_multigraph() else [data]


def _first(records, predicate):
    for record in records:
        if predicate(record):
            return record
    return None


def regulatory_motifs(graph, responsive_only: bool = False,
                      exclude_currency: bool = True) -> Dict[str, object]:
    """Find signed wiring patterns among the molecules.

    ``product_inhibition`` / ``product_activation``
        A reaction regulated by its own product.
    ``allosteric_feedback``
        A metabolite made by one reaction regulates another reaction that
        feeds it.
    ``feed_forward``
        A transcription factor controls two enzymes of the same reaction.

    Parameters
    ----------
    graph : networkx.MultiDiGraph
    responsive_only : bool
        Consider only molecules that changed (and factors that changed).
    exclude_currency : bool
        Skip currency metabolites such as ATP and water.

    Returns
    -------
    dict
        ``motifs``: one row per instance, with the molecules, signs and sign
        product (negative means the loop opposes its input). ``counts``:
        instances per motif type.
    """
    currency = set(CURRENCY_METABOLITES) if exclude_currency else set()
    responsive = None
    if responsive_only:
        responsive = {n for n, d in graph.nodes(data=True) if is_responsive(d)}

    def keep(node) -> bool:
        return node not in currency and (responsive is None or node in responsive)

    def name(node) -> str:
        data = graph.nodes[node]
        return str(data.get("symbol") or data.get("name") or node).split(";")[0][:40]

    rows = []
    reactions = [n for n, d in graph.nodes(data=True) if d.get("layer") == "Reactions"]

    for reaction in reactions:
        products = [v for _, v, d in graph.out_edges(reaction, data=True)
                    if d.get("edge_type") == "product" and keep(v)]
        regulators = {u: d for u, _, d in graph.in_edges(reaction, data=True)
                      if str(d.get("edge_type", "")).startswith("allosteric") and keep(u)}
        enzymes = [u for u, _, d in graph.in_edges(reaction, data=True)
                   if d.get("edge_type") == "catalysis" and keep(u)]

        # product inhibition: this reaction makes what regulates it
        for metabolite, record in regulators.items():
            if metabolite in products:
                sign = int(record.get("sign", 0) or 0)
                rows.append({
                    "motif": "product_inhibition" if sign < 0 else "product_activation",
                    "reaction": reaction, "reaction_name": name(reaction),
                    "metabolite": metabolite, "metabolite_name": name(metabolite),
                    "enzyme": name(enzymes[0]) if enzymes else None,
                    "sign_product": sign,
                })

        # allosteric feedback at range: a product of another reaction regulates
        # this one, and this one feeds that reaction's substrates
        for metabolite, record in regulators.items():
            if metabolite in products:
                continue
            sources = [u for u, _, d in graph.in_edges(metabolite, data=True)
                       if d.get("edge_type") == "product" and u != reaction]
            for upstream in sources:
                upstream_substrates = {u for u, _, d in graph.in_edges(upstream, data=True)
                                       if d.get("edge_type") == "substrate" and keep(u)}
                if upstream_substrates & set(products):
                    rows.append({
                        "motif": "allosteric_feedback",
                        "reaction": reaction, "reaction_name": name(reaction),
                        "metabolite": metabolite, "metabolite_name": name(metabolite),
                        "enzyme": name(upstream), "sign_product": int(record.get("sign", 0) or 0),
                    })
                    break

    # feed-forward: TF -> gene -> protein -> reaction, with the TF also
    # regulating a second gene whose protein hits the same reaction
    for factor, targets in _transcription_targets(graph, keep).items():
        reached: Dict[str, List[str]] = {}
        for gene in targets:
            for _, protein, data in graph.out_edges(gene, data=True):
                if data.get("edge_type") != "translation" or not keep(protein):
                    continue
                for _, reaction, catalysis in graph.out_edges(protein, data=True):
                    if catalysis.get("edge_type") == "catalysis":
                        reached.setdefault(reaction, []).append(gene)
        for reaction, genes in reached.items():
            if len(set(genes)) >= 2:
                rows.append({
                    "motif": "feed_forward",
                    "reaction": reaction, "reaction_name": name(reaction),
                    "metabolite": None, "metabolite_name": None,
                    "enzyme": name(factor),
                    "sign_product": 0,       # ChIP-Atlas binding carries no sign
                })

    motifs = pd.DataFrame(rows, columns=["motif", "reaction", "reaction_name",
                                         "metabolite", "metabolite_name", "enzyme",
                                         "sign_product"])
    counts = motifs["motif"].value_counts().to_dict() if not motifs.empty else {}
    logger.info(f"Motifs found: {counts}")
    return {"motifs": motifs, "counts": counts}


def _transcription_targets(graph, keep) -> Dict[str, List[str]]:
    """Factors with at least two targets that pass ``keep``.

    The factor must pass it too. ChIP-Atlas gives a well-studied factor
    thousands of targets, so counting feed-forward motifs through factors that
    did not themselves respond produces tens of thousands of instances and no
    information -- on the brown adipocyte network, 10,272 of them.
    """
    targets: Dict[str, List[str]] = {}
    for factor, gene, data in graph.edges(data=True):
        if (data.get("edge_type") == "transcriptional_regulation"
                and keep(gene) and keep(factor)):
            targets.setdefault(factor, []).append(gene)
    return {factor: genes for factor, genes in targets.items() if len(genes) >= 2}


def structural_vulnerability(graph, responsive_only: bool = True,
                             top_n: int = 25) -> pd.DataFrame:
    """Find the molecules whose removal splits the responsive network.

    These cut molecules are the single routes by which the response in one part
    of the network reaches another.

    Parameters
    ----------
    graph : networkx.Graph
    responsive_only : bool
        Analyse the responsive subnetwork (default).
    top_n : int, optional
        Return only the most damaging molecules.

    Returns
    -------
    pandas.DataFrame
        One row per cut molecule, most damaging first: ``node``, ``name``,
        ``layer``, ``fragments`` (pieces left after removal), ``largest_loss``
        (share of the response cut off) and ``cross_layer_degree``.
    """
    scope = graph
    if responsive_only:
        nodes = [n for n, d in graph.nodes(data=True) if is_responsive(d)]
        scope = graph.subgraph(nodes)
    simple = nx.Graph()
    simple.add_nodes_from(scope.nodes(data=True))
    simple.add_edges_from((u, v) for u, v in scope.edges())

    if simple.number_of_nodes() == 0:
        return pd.DataFrame(columns=["node", "name", "layer", "fragments",
                                     "largest_loss", "cross_layer_degree"])

    component = max(nx.connected_components(simple), key=len)
    largest = simple.subgraph(component).copy()

    rows = []
    for node in nx.articulation_points(largest):
        remaining = largest.copy()
        remaining.remove_node(node)
        pieces = sorted((len(c) for c in nx.connected_components(remaining)), reverse=True)
        layers = {str(graph.nodes[neighbour].get("layer")) for neighbour in largest[node]}
        data = graph.nodes[node]
        rows.append({
            "node": node,
            "name": str(data.get("symbol") or data.get("name") or node).split(";")[0][:40],
            "layer": data.get("layer"),
            "fragments": len(pieces),
            "largest_loss": sum(pieces[1:]) / max(largest.number_of_nodes() - 1, 1),
            "cross_layer_degree": len(layers - {str(data.get("layer"))}),
        })

    table = pd.DataFrame(rows)
    if table.empty:
        logger.info("No cut molecules: the response has no single point of failure")
        return table
    table = table.sort_values(["largest_loss", "fragments"], ascending=False)
    return table.head(top_n).reset_index(drop=True)


def convergence_significance(graph, n_randomisations: int = 100,
                             random_state: int = 0) -> Dict[str, object]:
    """Test whether changed enzymes and metabolites meet more than chance.

    Counts the reactions with both a changed enzyme and a changed metabolite,
    then repeats the count after reassigning at random which molecules changed,
    keeping the network and the number of changed molecules per layer fixed.

    Parameters
    ----------
    graph : networkx.MultiDiGraph
    n_randomisations : int
        Number of random reassignments.
    random_state : int, optional
        Seed for reproducibility.

    Returns
    -------
    dict
        ``observed``, ``null_mean``, ``z``, ``p_value``, and ``null`` (the
        random counts, for plotting).
    """
    rng = np.random.default_rng(random_state)

    reactions = [n for n, d in graph.nodes(data=True) if d.get("layer") == "Reactions"]
    enzymes_of, metabolites_of = {}, {}
    for reaction in reactions:
        enzymes_of[reaction] = [u for u, _, d in graph.in_edges(reaction, data=True)
                                if d.get("edge_type") == "catalysis"]
        metabolites_of[reaction] = [
            u for u, _, d in graph.in_edges(reaction, data=True)
            if d.get("edge_type") in ("substrate", "allosteric_activation",
                                      "allosteric_inhibition")
        ] + [v for _, v, d in graph.out_edges(reaction, data=True)
             if d.get("edge_type") == "product"]

    by_layer: Dict[str, List[str]] = {}
    for node, data in graph.nodes(data=True):
        layer = str(data.get("layer"))
        if layer != "Reactions":
            by_layer.setdefault(layer, []).append(node)

    def convergent(responsive: set) -> int:
        return sum(
            1 for reaction in reactions
            if any(e in responsive for e in enzymes_of[reaction])
            and any(m in responsive for m in metabolites_of[reaction])
        )

    responsive = {n for n, d in graph.nodes(data=True) if is_responsive(d)}

    # The enzyme side of a convergence is a catalysis edge, which comes from the
    # Proteome. A study without a proteome therefore scores zero convergence
    # whatever its data says, and that zero must not read as a result.
    if not any(e in responsive for enzymes in enzymes_of.values() for e in enzymes):
        logger.warning(
            "No responsive catalyst reaches any reaction: convergence is 0 by "
            "construction, not by measurement. This network has no measured "
            "Proteome layer, or none of its enzymes changed."
        )

    observed = convergent(responsive)

    null = []
    for _ in range(n_randomisations):
        shuffled = set()
        for layer, members in by_layer.items():
            n_changed = sum(1 for node in members if node in responsive)
            if n_changed:
                shuffled |= {members[i] for i in
                             rng.choice(len(members), size=n_changed, replace=False)}
        null.append(convergent(shuffled))

    null = np.array(null, dtype=float)
    spread = null.std(ddof=1) or 1.0
    result = {
        "observed": observed,
        "null_mean": float(null.mean()),
        "z": float((observed - null.mean()) / spread),
        "p_value": float((np.sum(null >= observed) + 1) / (len(null) + 1)),
        "null": null,
    }
    logger.info(f"Convergence: {observed} reactions, null {null.mean():.1f} "
                f"(z = {result['z']:.1f}, p = {result['p_value']:.3g})")
    return result
