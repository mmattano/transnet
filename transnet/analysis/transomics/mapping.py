"""Map measured omics data onto trans-omic network nodes.

One entry point, one set of node attributes, and an explicit report of what did
and did not match.  Identifier mismatch between an omics table and a network is
the most common failure in trans-omic analysis and the easiest to miss, so
:func:`map_omics_to_network` always tells you how many features landed.
"""

from typing import Any, Dict, List, Mapping, Optional, Sequence

import logging

import numpy as np
import pandas as pd

from transnet.biology.schema import available_layers

logger = logging.getLogger(__name__)

__all__ = [
    "MappingReport",
    "map_omics_to_network",
    "map_modification_sites",
    "regulated_nodes",
    "is_responsive",
    "is_measured",
    "build_name_id_map",
    "build_alias_id_map",
]


def is_responsive(node_data: Mapping[str, Any]) -> bool:
    """Did this molecule change, regardless of whether a direction is known?

    Analyses that need only "something happened here" -- subnetwork extraction,
    hub ranking, coverage -- use this, so an omnibus test such as a time-course
    ANOVA counts. Analyses that need a direction read ``regulated`` instead.
    """
    if "responsive" in node_data:
        return bool(node_data["responsive"])
    return bool(node_data.get("regulated", 0))


def is_measured(node_data: Mapping[str, Any]) -> bool:
    """Was this molecule measured at all, whatever columns its table had?

    "Measured and unchanged" is a different statement from "not measured", and
    several analyses turn on the difference -- a path through an unchanged
    molecule is contradicted, a path through an unmeasured one is untested.
    """
    return bool(
        node_data.get("measured")
        or node_data.get("responsive")
        or any(node_data.get(key) is not None
               for key in ("value", "log2fc", "qvalue"))
    )


class MappingReport:
    """What happened when omics tables were mapped onto a network.

    Attributes
    ----------
    per_layer : pandas.DataFrame
        One row per layer: features supplied, matched, unmatched, and how many
        of the matched were called regulated up or down.
    unmatched : dict
        ``{layer: [identifier, ...]}`` for features that found no node, capped
        at ``max_unmatched_examples`` per layer.
    """

    def __init__(self, per_layer: pd.DataFrame, unmatched: Dict[str, List[str]]):
        self.per_layer = per_layer
        self.unmatched = unmatched

    @property
    def match_rate(self) -> float:
        """Overall fraction of supplied features that matched a node."""
        supplied = self.per_layer["n_supplied"].sum()
        if not supplied:
            return 0.0
        return float(self.per_layer["n_matched"].sum() / supplied)

    def __repr__(self):
        return (
            f"<MappingReport {self.per_layer['n_matched'].sum()}/"
            f"{self.per_layer['n_supplied'].sum()} features mapped "
            f"({self.match_rate:.1%})>"
        )

    def __str__(self):
        lines = ["Omics -> network mapping", "-" * 60]
        lines.append(self.per_layer.to_string(index=False))

        poor = self.per_layer[self.per_layer["match_fraction"] < 0.5]
        for _, row in poor.iterrows():
            layer = row["layer"]
            examples = ", ".join(self.unmatched.get(layer, [])[:5])
            covered = row.get("layer_covered", 0.0)

            # A low match rate has two very different causes, and conflating
            # them sends people to debug the wrong thing.
            if covered >= 0.5:
                lines.append(
                    f"\n  - {layer}: {row['match_fraction']:.0%} of your features "
                    f"matched, but {covered:.0%} of the network's {layer} nodes "
                    f"were covered."
                )
                lines.append(
                    "    The identifiers are fine; the network is simply smaller "
                    "than the assay."
                )
            else:
                lines.append(
                    f"\n  ! {layer}: only {row['match_fraction']:.0%} of features "
                    f"matched, covering {covered:.0%} of the network's {layer} "
                    f"nodes. This usually means an identifier-type mismatch."
                )
                lines.append(
                    "    See build_alias_id_map() / build_name_id_map() to "
                    "translate."
                )
            if examples:
                lines.append(f"    unmatched examples: {examples}")
        return "\n".join(lines)


def _coerce_float(value) -> Optional[float]:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return None if np.isnan(result) else result


def map_omics_to_network(
    graph,
    tables: Mapping[str, pd.DataFrame],
    *,
    id_column: str = None,
    value_column: str = None,
    log2fc_column: str = None,
    qvalue_column: str = None,
    qvalue_threshold: float = 0.05,
    log2fc_threshold: float = 0.0,
    id_map: Optional[Mapping[str, Mapping[str, str]]] = None,
    max_unmatched_examples: int = 20,
    se_column: str = None,
) -> MappingReport:
    """Write measured values onto the nodes of ``graph``, in place.

    Every node that matches a supplied feature gains ``value``, ``log2fc``,
    ``qvalue`` and ``regulated``; ``regulated`` is ``+1`` / ``-1`` / ``0`` and
    is what the rest of the trans-omics analyses read.

    Parameters
    ----------
    graph : networkx.Graph
        A trans-omic network.  Nodes must carry a ``layer`` attribute.
    tables : mapping of str to pandas.DataFrame
        ``{layer_name: table}``, e.g.
        ``{"Transcriptome": deg_table, "Metabolome": dem_table}``.  Layers not
        present in the graph are reported and skipped -- supplying a
        phosphoproteomics table to a network with no Signaling layer is a
        warning, not an error.
    id_column, value_column, log2fc_column, qvalue_column : str, optional
        Column names, applied to every table.  ``id_column`` defaults to the
        first column.  If ``log2fc_column`` is absent but ``value_column`` is
        given, the value is used for direction as well.
    qvalue_threshold : float
        A feature is called regulated only if its q-value is at or below this.
        Tables with no q-value column are called on fold change alone.
    log2fc_threshold : float
        Minimum absolute log2 fold change for a feature to count as regulated.
    id_map : mapping, optional
        ``{layer: {table_identifier: node_identifier}}`` for layers whose
        table uses a different identifier space than the network.
    max_unmatched_examples : int
        How many unmatched identifiers to keep per layer for the report.

    Returns
    -------
    MappingReport
        Print it -- it names any layer where most features failed to match.

    Examples
    --------
    >>> report = map_omics_to_network(
    ...     G,
    ...     {"Transcriptome": degs, "Metabolome": dems},
    ...     id_column="id", log2fc_column="log2FC", qvalue_column="padj",
    ... )
    >>> print(report)
    """
    nodes_by_layer: Dict[str, Dict[str, Any]] = {}
    for node, data in graph.nodes(data=True):
        layer = data.get("layer")
        if layer:
            nodes_by_layer.setdefault(str(layer), {})[str(node)] = node

    graph_layers = set(nodes_by_layer)
    rows = []
    unmatched: Dict[str, List[str]] = {}

    for layer, table in tables.items():
        layer = str(layer)
        if table is None or len(table) == 0:
            continue

        if layer not in graph_layers:
            logger.warning(
                f"No '{layer}' layer in this network; skipping its {len(table)} "
                f"features. Present layers: {sorted(graph_layers)}"
            )
            rows.append({
                "layer": layer, "n_supplied": len(table), "n_matched": 0,
                "match_fraction": 0.0, "n_layer_nodes": 0, "layer_covered": 0.0,
                "n_up": 0, "n_down": 0, "n_undirected": 0,
                "note": "layer absent from network",
            })
            unmatched[layer] = [str(v) for v in table.iloc[:, 0].head(max_unmatched_examples)]
            continue

        layer_nodes = nodes_by_layer[layer]
        translation = dict(id_map.get(layer, {})) if id_map else {}

        ids = id_column if id_column and id_column in table.columns else table.columns[0]

        matched = 0
        n_up = n_down = n_undirected = 0
        matched_nodes = set()
        misses: List[str] = []

        for _, row in table.iterrows():
            raw_id = str(row[ids])
            node_id = translation.get(raw_id, raw_id)
            node = layer_nodes.get(node_id)
            if node is None:
                if len(misses) < max_unmatched_examples:
                    misses.append(raw_id)
                continue

            value = _coerce_float(row[value_column]) if value_column and value_column in table.columns else None
            log2fc = _coerce_float(row[log2fc_column]) if log2fc_column and log2fc_column in table.columns else None
            qvalue = _coerce_float(row[qvalue_column]) if qvalue_column and qvalue_column in table.columns else None

            if log2fc is None:
                log2fc = value

            significant = True if qvalue is None else qvalue <= qvalue_threshold
            if significant and log2fc is not None:
                significant = abs(log2fc) > log2fc_threshold

            direction = 0
            if significant and log2fc is not None:
                direction = 1 if log2fc > 0 else -1

            attrs = graph.nodes[node]
            attrs["value"] = value
            attrs["log2fc"] = log2fc
            attrs["qvalue"] = qvalue
            if se_column and se_column in table.columns:
                # standard error of log2fc; lets analyses compare changes
                # across layers without a significance threshold
                attrs["se"] = _coerce_float(row[se_column])
            # Explicit, because none of the three above is reliably set: a
            # table with only log2FC and q-value leaves `value` empty, and
            # figures inferring "measured" from it drew every measured but
            # unchanged molecule as unmeasured.
            attrs["measured"] = True
            # `responsive` is "this molecule changed"; `regulated` is "...and in
            # this direction". An ANOVA or other omnibus test gives the first
            # without the second, which is a normal result, not a missing one.
            attrs["responsive"] = bool(significant)
            attrs["regulated"] = direction

            matched += 1
            matched_nodes.add(str(node))
            if direction > 0:
                n_up += 1
            elif direction < 0:
                n_down += 1
            elif significant:
                n_undirected += 1

        if misses:
            unmatched[layer] = misses

        rows.append({
            "layer": layer,
            "n_supplied": len(table),
            "n_matched": matched,
            "match_fraction": matched / len(table) if len(table) else 0.0,
            "n_layer_nodes": len(layer_nodes),
            # Distinct nodes, not rows: two aliases can resolve to one node,
            # which would otherwise push coverage above 100%.
            "layer_covered": (
                len(matched_nodes) / len(layer_nodes) if layer_nodes else 0.0
            ),
            "n_up": n_up,
            "n_down": n_down,
            "n_undirected": n_undirected,
            "note": "",
        })

    per_layer = pd.DataFrame(rows, columns=[
        "layer", "n_supplied", "n_matched", "match_fraction",
        "n_layer_nodes", "layer_covered",
        "n_up", "n_down", "n_undirected", "note",
    ])
    report = MappingReport(per_layer, unmatched)

    logger.info(str(report))
    return report


def regulated_nodes(
    graph,
    layers: Optional[Sequence[str]] = None,
    direction: Optional[int] = None,
    include_undirected: bool = False,
) -> List[str]:
    """Nodes called differentially regulated by :func:`map_omics_to_network`.

    Parameters
    ----------
    graph : networkx.Graph
    layers : sequence of str, optional
        Restrict to these layers.  Defaults to every layer present.
    direction : int, optional
        ``1`` for up, ``-1`` for down, ``None`` (default) for either.
    include_undirected : bool
        Also return molecules called significant by a test that gives no
        direction (an ANOVA, say). Ignored when ``direction`` is given.

    Returns
    -------
    list of str
    """
    wanted = set(layers) if layers is not None else set(available_layers(graph))
    result = []
    for node, data in graph.nodes(data=True):
        if data.get("layer") not in wanted:
            continue
        state = data.get("regulated", 0) or 0
        if direction is not None:
            if state != direction:
                continue
        elif state == 0 and not (include_undirected and is_responsive(data)):
            continue
        result.append(node)
    return result


def build_name_id_map(
    graph,
    names: Sequence[str],
    layer: str = "Metabolome",
) -> Dict[str, str]:
    """Resolve a table's compound *names* to the identifiers the network uses.

    Real metabolomics tables are keyed by name ("Isocitrate", "alpha-Aminoadipic
    acid"), while the network is keyed by KEGG compound id. This builds the
    ``id_map`` that :func:`map_omics_to_network` needs, matching against every
    synonym in each node's ``name`` (KEGG stores them ``;``-separated),
    case-insensitively and ignoring hyphen/space differences.

    Parameters
    ----------
    graph : networkx.Graph
    names : sequence of str
        The identifiers appearing in your table.
    layer : str
        Layer to resolve against.

    Returns
    -------
    dict
        ``{name: node_id}`` for the names that resolved. Names absent from the
        result did not match any node -- pass the dict straight to
        ``map_omics_to_network(..., id_map={layer: mapping})`` and the mapping
        report will count the rest as unmatched.

    Examples
    --------
    >>> mapping = build_name_id_map(graph, metabolomics["feature"])
    >>> map_omics_to_network(graph, {"Metabolome": metabolomics},
    ...                      id_map={"Metabolome": mapping})
    """
    def _normalise(value: str) -> str:
        text = str(value).strip().lower()
        for char in "-_ ":
            text = text.replace(char, "")
        return text

    index: Dict[str, str] = {}
    for node, data in graph.nodes(data=True):
        if data.get("layer") != layer:
            continue
        index.setdefault(_normalise(node), str(node))
        for synonym in str(data.get("name") or "").split(";"):
            key = _normalise(synonym)
            # First synonym wins, so a common name is not shadowed by an
            # obscure one belonging to a different compound.
            if key and key not in index:
                index[key] = str(node)

    mapping = {
        str(name): index[_normalise(name)]
        for name in names if _normalise(name) in index
    }

    logger.info(
        f"Resolved {len(mapping)}/{len(list(names))} '{layer}' names to node ids"
    )
    return mapping


#: How each layer's elements are stored, and which attribute becomes the node
#: id.  Mirrors ``Transnet._node_metadata``.
_LAYER_ELEMENTS = {
    "Transcriptome": ("genes", ("ncbi_id", "ensembl_id")),
    "Proteome": ("proteins", ("uniprot_id",)),
    "Metabolome": ("metabolites", ("kegg_compound_id",)),
    "Reactions": ("reactions", ("id",)),
    "Signaling": ("nodes", ("id",)),
}


def build_alias_id_map(
    network,
    layer: str,
    alias_attribute: str,
) -> Dict[str, str]:
    """Map a layer's *other* identifiers to the ones its nodes use.

    A network node is keyed by one identifier -- Entrez for genes, KEGG for
    metabolites -- but the element behind it usually knows several. An omics
    table keyed by Ensembl or PubChem therefore matches nothing, even though
    the translation is sitting in the layer already. This extracts it, with no
    database call.

    Parameters
    ----------
    network : transnet.Transnet
        A network whose layers are populated.
    layer : str
        Layer to read, e.g. ``"Metabolome"``.
    alias_attribute : str
        The element attribute holding your table's identifier, e.g.
        ``"pubchem_id"`` or ``"ensembl_id"``. List-valued attributes are
        expanded, so every alias maps to the node.

    Returns
    -------
    dict
        ``{alias: node_id}``, ready to pass as
        ``map_omics_to_network(..., id_map={layer: mapping})``.

    Examples
    --------
    >>> id_map = {
    ...     "Transcriptome": build_alias_id_map(net, "Transcriptome", "ensembl_id"),
    ...     "Metabolome": build_alias_id_map(net, "Metabolome", "pubchem_id"),
    ... }
    >>> map_omics_to_network(graph, tables, id_map=id_map)
    """
    if layer not in _LAYER_ELEMENTS:
        raise ValueError(
            f"Unknown layer {layer!r}; expected one of {sorted(_LAYER_ELEMENTS)}"
        )

    collection_name, node_attributes = _LAYER_ELEMENTS[layer]
    layer_object = getattr(network, layer.lower(), None)
    if layer_object is None:
        logger.warning(f"Network has no {layer} layer; alias map is empty")
        return {}

    elements = getattr(layer_object, collection_name, None) or []

    mapping: Dict[str, str] = {}
    for element in elements:
        node_id = None
        for attribute in node_attributes:
            value = getattr(element, attribute, None)
            if value is not None and str(value) not in ("", "nan", "None"):
                node_id = str(value)
                break
        if node_id is None:
            continue

        alias = getattr(element, alias_attribute, None)
        aliases = alias if isinstance(alias, (list, tuple, set)) else [alias]
        for item in aliases:
            if item is None:
                continue
            text = str(item).strip()
            if text and text not in ("nan", "None"):
                mapping.setdefault(text, node_id)
                # Numeric identifiers (PubChem CIDs, Entrez) appear with and
                # without a trailing ".0" depending on how they were read.
                if text.endswith(".0"):
                    mapping.setdefault(text[:-2], node_id)

    logger.info(
        f"Built {layer} alias map from '{alias_attribute}': "
        f"{len(mapping):,} identifiers -> node ids"
    )
    return mapping


def map_modification_sites(
    graph,
    table: pd.DataFrame,
    protein_column: str = "protein",
    site_column: str = "site",
    log2fc_column: str = "log2FC",
    qvalue_column: str = "padj",
    qvalue_threshold: float = 0.05,
    log2fc_threshold: float = 0.0,
    id_map: Optional[Mapping[str, str]] = None,
    edge_type: str = "phosphorylation",
    sign: int = 0,
    evidence: str = "measured modification site",
) -> Dict[str, int]:
    """Attach measured modification sites to the proteins they sit on.

    Phosphoproteomics measures sites, not proteins, and a site is a regulator of
    its protein's activity rather than a measure of how much of it there is. So
    each site becomes its own ``Signaling`` node, ``<protein>_<site>``, with an
    edge into the protein. Several sites on one protein stay separate nodes,
    because they can move in opposite directions.

    The edge is unsigned by default. Whether more phosphorylation at a given site
    raises or lowers catalytic activity is site-specific and usually unrecorded,
    and :func:`~transnet.reaction_regulation_table` reports the site's direction
    and its effect on the reaction in separate columns for that reason. Pass
    ``sign`` only for a set of sites whose effect is actually known.

    Parameters
    ----------
    graph : networkx.MultiDiGraph
        A network to add to, modified in place.
    table : pandas.DataFrame
        One row per site, as :func:`transnet.datasets.load_motrpac_ptm` returns.
    protein_column, site_column : str
        Columns naming the protein and the site.
    log2fc_column, qvalue_column : str
        Columns holding the measured change and its adjusted p-value.
    qvalue_threshold, log2fc_threshold : float
        A site counts as changed at or below the q-value and at or above the
        absolute fold change.
    id_map : mapping, optional
        ``{protein in the table: node in the graph}``, for the usual case where
        the assay reports RefSeq accessions and the network is keyed by UniProt.
        Proteins absent from the map are looked up directly.
    edge_type : {"phosphorylation", ...}
        The relationship to record. Acetylation and ubiquitination have no edge
        type of their own in the schema, so they are attached as
        ``phosphorylation`` only if you mean them to feed the phospho axis;
        otherwise keep them out of the graph and analyse the table directly.
    sign : int
        Edge sign: 0 (unknown) unless the effect on activity is known.
    evidence : str
        Provenance recorded on each edge.

    Returns
    -------
    dict
        ``n_sites`` supplied, ``n_mapped`` attached to a protein in the graph,
        ``n_changed`` of those that passed the thresholds, and ``n_proteins``
        carrying at least one attached site.
    """
    proteins = {
        node for node, data in graph.nodes(data=True)
        if data.get("layer") == "Proteome"
    }
    if not proteins:
        logger.warning("No Proteome layer: modification sites have nothing to "
                       "attach to")
        return {"n_sites": len(table), "n_mapped": 0, "n_changed": 0,
                "n_proteins": 0}

    mapped, changed, touched = 0, 0, set()
    for row in table.itertuples(index=False):
        source = str(getattr(row, protein_column))
        target = (id_map or {}).get(source, source)
        if target not in proteins:
            continue

        site = str(getattr(row, site_column))
        node = f"{target}_{site}"
        log2fc = _coerce_float(getattr(row, log2fc_column, None))
        qvalue = _coerce_float(getattr(row, qvalue_column, None))

        is_changed = (
            qvalue is not None and qvalue <= qvalue_threshold
            and log2fc is not None and abs(log2fc) >= log2fc_threshold
        )
        graph.add_node(
            node, layer="Signaling", node_type="SignalingProtein",
            name=f"{graph.nodes[target].get('symbol', target)} {site}",
            symbol=f"{graph.nodes[target].get('symbol', target)}_{site}",
            measured=True, log2fc=log2fc, qvalue=qvalue,
            regulated=int(np.sign(log2fc)) if is_changed and log2fc else 0,
        )
        graph.add_edge(
            node, target, key=edge_type, edge_type=edge_type, role=edge_type,
            sign=sign, directed=True, weight=1.0, stoichiometry=None, ec="",
            source_db="measured", evidence=evidence, confidence=None,
        )
        mapped += 1
        changed += bool(is_changed)
        touched.add(target)

    logger.info(
        f"{mapped:,} of {len(table):,} sites attached to {len(touched):,} "
        f"proteins; {changed:,} changed at q <= {qvalue_threshold}"
    )
    return {"n_sites": len(table), "n_mapped": mapped, "n_changed": changed,
            "n_proteins": len(touched)}
