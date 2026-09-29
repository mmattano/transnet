"""Schema for the trans-omic network: layers, edge types, and attribute keys.

This module is the single source of truth for the vocabulary used across the
whole package.  Builders emit these ``edge_type`` values, the analysis modules
dispatch on them, and the visualisation and export code style them.  Nothing
else in TransNet should hard-code an edge type or a layer name.

The taxonomy follows the trans-omic framework of Yugi & Kuroda, in which a
biochemical network is reconstructed by connecting omic layers through a small
number of well-defined regulatory relationships.

References
----------
Yugi K, Kubota H, Hatano A, Kuroda S. Trans-Omics: How To Reconstruct
Biochemical Networks Across Multiple 'Omic' Layers. *Trends in Biotechnology*
34(4):276-290, 2016.
"""

from typing import Dict, Iterable, List, Optional, Set

__all__ = [
    "LAYERS",
    "LAYER_HIERARCHY",
    "EDGE_TYPES",
    "EdgeType",
    "INTERACTION_COLUMNS",
    "NODE_DATA_ATTRS",
    "GENE_AXIS_EDGE_TYPES",
    "METABOLITE_AXIS_EDGE_TYPES",
    "CURRENCY_METABOLITES",
    "metabolic_pool",
    "available_layers",
    "available_edge_types",
    "edge_type_info",
    "is_cross_layer",
    "top_layer_present",
]


# --------------------------------------------------------------------------
# Layers
# --------------------------------------------------------------------------

#: All layer names TransNet understands.
LAYERS = (
    "Signaling",
    "Transcriptome",
    "Proteome",
    "Reactions",
    "Metabolome",
    "Pathways",
)

#: Regulatory hierarchy, from the top of the trans-omic network downwards.
#: Used to infer a sensible default source layer for path tracing when the
#: caller does not name one.  Layers absent from a given network are simply
#: skipped -- see :func:`available_layers`.
LAYER_HIERARCHY = (
    "Signaling",
    "Proteome",
    "Transcriptome",
    "Reactions",
    "Metabolome",
)
"""Regulatory order, used to decide where a regulatory path *starts*.

Proteome sits above Transcriptome here because transcription factors are
proteins that regulate genes. This is logic, not presentation -- see
:data:`DISPLAY_ORDER` for how layers are shown.
"""

DISPLAY_ORDER = (
    "Signaling",
    "Transcriptome",
    "Proteome",
    "Reactions",
    "Metabolome",
)
"""Order in which layers are presented: figures, tables and reports.

The trans-omic convention reads down the flow of information -- gene, then
protein, then the reaction the protein catalyses, then the metabolite -- so
this is the order readers expect. It differs from :data:`LAYER_HIERARCHY` only
in putting Transcriptome above Proteome; a protein appears twice in the
regulatory chain (as transcription factor and as enzyme), and a single linear
order cannot show both.
"""


class EdgeType:
    """Definition of one inter- or intra-layer regulatory relationship."""

    __slots__ = (
        "name",
        "source_layer",
        "target_layer",
        "sign",
        "directed",
        "source_db",
        "description",
    )

    def __init__(
        self,
        name: str,
        source_layer: str,
        target_layer: str,
        sign: int = 0,
        directed: bool = True,
        source_db: str = "",
        description: str = "",
    ):
        self.name = name
        self.source_layer = source_layer
        self.target_layer = target_layer
        #: +1 activating, -1 inhibiting, 0 sign unknown / not applicable.
        self.sign = sign
        self.directed = directed
        self.source_db = source_db
        self.description = description

    def __repr__(self):
        return (
            f"<EdgeType {self.name}: {self.source_layer}->{self.target_layer} "
            f"sign={self.sign:+d}>"
        )


#: The fixed edge-type vocabulary.  ``sign`` is the *default* for the type;
#: individual edges may override it (e.g. a phosphorylation known to be
#: inhibitory).
EDGE_TYPES: Dict[str, EdgeType] = {
    et.name: et
    for et in [
        EdgeType(
            "phosphorylation",
            "Signaling",
            "Proteome",
            sign=0,
            source_db="KEGG / user phosphoproteomics",
            description="Kinase phosphorylates a substrate protein.",
        ),
        EdgeType(
            "kinase_tf",
            "Signaling",
            "Proteome",
            sign=0,
            source_db="KEGG signaling pathways",
            description="Signaling protein regulates a transcription factor.",
        ),
        EdgeType(
            "transcriptional_regulation",
            "Proteome",
            "Transcriptome",
            sign=0,
            source_db="ChIP-Atlas",
            description="Transcription factor binds and regulates a target gene.",
        ),
        EdgeType(
            "translation",
            "Transcriptome",
            "Proteome",
            sign=1,
            source_db="NCBI/UniProt identifier join",
            description="Gene is translated into its protein product.",
        ),
        EdgeType(
            "protein_interaction",
            "Proteome",
            "Proteome",
            sign=0,
            directed=False,
            source_db="STRING",
            description="Physical or functional protein-protein association.",
        ),
        EdgeType(
            "catalysis",
            "Proteome",
            "Reactions",
            sign=1,
            source_db="KEGG (EC number)",
            description="Enzyme catalyses a metabolic reaction.",
        ),
        EdgeType(
            "gene_catalysis",
            "Transcriptome",
            "Reactions",
            sign=1,
            source_db="KEGG (EC number)",
            description=(
                "Gene encodes an enzyme catalysing this reaction. Used when no "
                "Proteome layer is available, so that a transcriptome + "
                "metabolome study still reaches the reaction layer."
            ),
        ),
        EdgeType(
            "substrate",
            "Metabolome",
            "Reactions",
            sign=1,
            source_db="KEGG reaction equation",
            description="Metabolite is consumed by a reaction.",
        ),
        EdgeType(
            "product",
            "Reactions",
            "Metabolome",
            sign=1,
            source_db="KEGG reaction equation",
            description="Metabolite is produced by a reaction.",
        ),
        EdgeType(
            "allosteric_activation",
            "Metabolome",
            "Reactions",
            sign=1,
            source_db="BRENDA",
            description="Metabolite allosterically activates the reaction's enzyme.",
        ),
        EdgeType(
            "allosteric_inhibition",
            "Metabolome",
            "Reactions",
            sign=-1,
            source_db="BRENDA",
            description="Metabolite allosterically inhibits the reaction's enzyme.",
        ),
        EdgeType(
            "enzymatic",
            "Proteome",
            "Metabolome",
            sign=0,
            source_db="KEGG (EC number)",
            description=(
                "Enzyme acts on a metabolite; a reaction-free shortcut retained "
                "for networks built without a Reactions layer."
            ),
        ),
    ]
}


#: Column order of the canonical interaction DataFrame written to disk.
INTERACTION_COLUMNS = [
    "source",
    "target",
    "source_layer",
    "target_layer",
    "edge_type",
    "role",
    "sign",
    "directed",
    "weight",
    "stoichiometry",
    "ec",
    "source_db",
    "evidence",
    "confidence",
]


#: Node attribute keys that carry measured or derived experimental data.
NODE_DATA_ATTRS = (
    "value",
    "log2fc",
    "qvalue",
    "regulated",
    "t_half",
    "ec50",
)


#: Cofactors and ubiquitous small molecules that participate in hundreds of
#: reactions.  As *substrates and products* they carry no specificity: a change
#: in NADH marks every dehydrogenase in the network as metabolite-regulated,
#: which is an artefact of connectivity rather than a finding.  Excluding them
#: from mass-action contributions is standard practice in metabolic network
#: analysis.
#:
#: Their *allosteric* roles are a different matter and are never excluded -- AMP
#: activating phosphofructokinase is exactly the kind of specific regulation
#: this package exists to find.
CURRENCY_METABOLITES = frozenset({
    "C00001",  # H2O
    "C00002",  # ATP
    "C00003",  # NAD+
    "C00004",  # NADH
    "C00005",  # NADPH
    "C00006",  # NADP+
    "C00007",  # O2
    "C00008",  # ADP
    "C00009",  # Orthophosphate
    "C00010",  # CoA
    "C00011",  # CO2
    "C00013",  # Diphosphate
    "C00014",  # NH3
    "C00020",  # AMP
    "C00027",  # H2O2
    "C00080",  # H+
    "C00016",  # FAD
    "C01352",  # FADH2
    "C00035",  # GDP
    "C00044",  # GTP
    "C00015",  # UDP
    "C00075",  # UTP
    "C00063",  # CTP
    "C00112",  # CDP
    "C00019",  # S-adenosyl-L-methionine
    "C00021",  # S-adenosyl-L-homocysteine
})


#: Edge types that carry regulation of a reaction through gene expression /
#: enzyme abundance -- the "gene regulation axis".
GENE_AXIS_EDGE_TYPES = frozenset({
    "catalysis", "gene_catalysis", "translation",
    "transcriptional_regulation", "kinase_tf", "phosphorylation",
})

#: Edge types that carry regulation of a reaction by metabolites -- the
#: "allosteric / metabolic regulation axis".
METABOLITE_AXIS_EDGE_TYPES = frozenset(
    {"allosteric_activation", "allosteric_inhibition", "substrate", "product"}
)


# --------------------------------------------------------------------------
# Introspection helpers
#
# Every analysis in TransNet discovers what a network actually contains rather
# than assuming a fixed set of layers.  A network without a Signaling layer, or
# without a Proteome, is a first-class trans-omic network -- analyses degrade by
# reporting weaker evidence, never by failing.
# --------------------------------------------------------------------------


def order_layers(layers) -> List[str]:
    """Sort layer names into :data:`DISPLAY_ORDER`, unknown layers last.

    For anything that presents layers -- figures, tables, menus -- so none of
    them falls back to alphabetical order, which puts Metabolome first.
    """
    layers = list(dict.fromkeys(layers))
    known = [layer for layer in DISPLAY_ORDER if layer in layers]
    return known + sorted(layer for layer in layers if layer not in DISPLAY_ORDER)


def available_layers(graph) -> List[str]:
    """Return the layers actually present in ``graph``, in display order.

    Ordered by :data:`DISPLAY_ORDER` (Transcriptome, Proteome, Reactions,
    Metabolome). Layers not named there are appended alphabetically so that
    user-defined layers are never silently dropped. Code that needs the
    *regulatory* order -- where a path starts -- uses :data:`LAYER_HIERARCHY`
    or :func:`top_layer_present`, not this.

    Parameters
    ----------
    graph : networkx.Graph
        Any NetworkX graph whose nodes carry a ``layer`` attribute.

    Returns
    -------
    list of str
    """
    present: Set[str] = set()
    for _, data in graph.nodes(data=True):
        layer = data.get("layer")
        if layer:
            present.add(str(layer))

    return order_layers(present)


def available_edge_types(graph) -> Dict[str, int]:
    """Return a count of each ``edge_type`` present in ``graph``.

    Edges written before the typed schema, or loaded from a legacy file, carry
    no ``edge_type``; those are counted under ``"unknown"``.

    Parameters
    ----------
    graph : networkx.Graph

    Returns
    -------
    dict
        Mapping of edge type name to the number of edges of that type.
    """
    counts: Dict[str, int] = {}
    for _, _, data in graph.edges(data=True):
        etype = data.get("edge_type") or "unknown"
        counts[etype] = counts.get(etype, 0) + 1
    return dict(sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])))


def edge_type_info(name: str) -> Optional[EdgeType]:
    """Look up an :class:`EdgeType` definition by name, or ``None``."""
    return EDGE_TYPES.get(name)


def is_cross_layer(edge_type: str) -> bool:
    """True if ``edge_type`` connects two different layers.

    Unknown edge types are treated as cross-layer, since the common case for an
    untyped edge in this package is an inter-layer link from a legacy file.
    """
    info = EDGE_TYPES.get(edge_type)
    if info is None:
        return True
    return info.source_layer != info.target_layer


def top_layer_present(graph, candidates: Optional[Iterable[str]] = None) -> Optional[str]:
    """Return the highest layer of :data:`LAYER_HIERARCHY` present in ``graph``.

    This is how path-tracing picks a starting layer when the caller does not
    name one: a network with phosphoproteomics starts at ``Signaling``, one
    without starts at ``Proteome``, and a transcriptome-only network starts at
    ``Transcriptome``.

    Returns
    -------
    str or None
        ``None`` if the graph has no recognised layer at all.
    """
    present = set(available_layers(graph))
    order = list(candidates) if candidates is not None else list(LAYER_HIERARCHY)
    for layer in order:
        if layer in present:
            return layer
    return None


def metabolic_pool(graph) -> set:
    """Metabolites that take part in a reaction of this organism's network.

    BRENDA reports every effector an enzyme was ever tested with, including
    laboratory reagents and inhibitors from other species -- p-chloromercuri-
    benzoate, trichloroethene conjugates. They are real measurements, so the
    allosteric edges keep them; but a molecule the organism's own metabolic
    network neither makes nor consumes cannot be part of its regulation, so
    figures select regulators from this pool rather than from every effector.

    The pool is read off the graph, not from a curated list: a metabolite
    qualifies when some reaction has it as a substrate or a product.
    """
    pool = set()
    for u, v, data in graph.edges(data=True):
        edge_type = data.get("edge_type")
        if edge_type == "substrate":
            pool.add(u)
        elif edge_type == "product":
            pool.add(v)
    return pool
