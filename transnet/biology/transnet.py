"""Connecting trans-omic networks.

The network TransNet builds is a *typed, directed, signed* regulatory
hierarchy, following the trans-omic framework of Yugi & Kuroda:

    signal -> TF -> gene -> enzyme protein -> REACTION <- metabolite

with the metabolic reaction as the convergence point.  Every edge carries the
relationship it represents (``edge_type``), the direction of its regulatory
effect (``sign``), and the evidence it came from.  See
:mod:`transnet.biology.schema` for the vocabulary.

Layers are optional.  A network built from transcriptomics and metabolomics
alone is a first-class trans-omic network -- the builders for layers that are
absent simply contribute no edges.

References
----------
Yugi K, Kubota H, Hatano A, Kuroda S. Trans-Omics: How To Reconstruct
Biochemical Networks Across Multiple 'Omic' Layers. *Trends in Biotechnology*
34(4):276-290, 2016.
"""

from typing import List, Dict, Any, Optional
import numpy as np
import pandas as pd
import networkx as nx
import logging
from .elements import Reaction, Metabolite, Gene, Protein
from .layers import (
    Reactions, Pathways, Transcriptome, Proteome, Metabolome, Signaling,
)
from .schema import EDGE_TYPES, INTERACTION_COLUMNS

logger = logging.getLogger(__name__)

__all__ = [
    "Transnet",
    "to_simple_graph",
]


def _as_str_list(val) -> List[str]:
    """Normalise ``None`` / scalar / list / ndarray into a list of strings."""
    if val is None:
        return []
    if isinstance(val, np.ndarray):
        return [str(v) for v in val.flat if v is not None and str(v) != "nan"]
    if isinstance(val, (list, tuple, set)):
        return [str(v) for v in val if v is not None and str(v) != "nan"]
    if isinstance(val, float) and np.isnan(val):
        return []
    return [str(val)]


#: Three-letter amino acid codes as BRENDA writes them ("L-Ala").
_AMINO_ACIDS = {
    "ala": "alanine", "arg": "arginine", "asn": "asparagine", "asp": "aspartate",
    "cys": "cysteine", "gln": "glutamine", "glu": "glutamate", "gly": "glycine",
    "his": "histidine", "ile": "isoleucine", "leu": "leucine", "lys": "lysine",
    "met": "methionine", "phe": "phenylalanine", "pro": "proline",
    "ser": "serine", "thr": "threonine", "trp": "tryptophan", "tyr": "tyrosine",
    "val": "valine",
}


def _name_variants(name: str) -> List[str]:
    """Spellings of a compound name that KEGG might use instead.

    BRENDA and KEGG disagree in predictable ways: amino acids by three-letter
    code ("L-Ala" / "L-Alanine"), "diphosphate" for "bisphosphate", and acids
    against their anions ("glutamic acid" / "glutamate"). Order matters: the
    name as given is always tried first.
    """
    variants = [name]
    lowered = name.strip().lower()
    stem = lowered[2:] if lowered[:2] in ("l-", "d-") else lowered
    prefix = lowered[:2] if lowered[:2] in ("l-", "d-") else ""
    if stem in _AMINO_ACIDS:
        full = _AMINO_ACIDS[stem]
        variants += [f"{prefix or 'l-'}{full}", full]
        if full.endswith("ate"):
            variants.append(f"{prefix or 'l-'}{full[:-3]}ic acid")
    if "diphosphate" in lowered:
        variants.append(lowered.replace("diphosphate", "bisphosphate"))
    if lowered.endswith("ic acid"):
        variants.append(lowered[: -len("ic acid")] + "ate")
    elif lowered.endswith("ate") and not lowered.endswith("phosphate"):
        variants.append(lowered[:-3] + "ic acid")
    return list(dict.fromkeys(variants))


def _edge(
    source,
    target,
    edge_type: str,
    *,
    role: str = "",
    sign: Optional[int] = None,
    weight: float = 1.0,
    stoichiometry=None,
    ec: str = "",
    evidence: str = "",
    confidence=None,
    source_layer: Optional[str] = None,
    target_layer: Optional[str] = None,
) -> Dict[str, Any]:
    """Build one canonical edge record.

    Layers, sign and directedness default to the definition of ``edge_type`` in
    :data:`transnet.biology.schema.EDGE_TYPES`, so callers only state what is
    specific to the individual edge.
    """
    info = EDGE_TYPES.get(edge_type)
    if info is None:
        raise ValueError(
            f"Unknown edge_type {edge_type!r}. "
            f"Valid types: {sorted(EDGE_TYPES)}"
        )
    return {
        "source": str(source),
        "target": str(target),
        "source_layer": source_layer or info.source_layer,
        "target_layer": target_layer or info.target_layer,
        "edge_type": edge_type,
        "role": role or edge_type,
        "sign": info.sign if sign is None else int(sign),
        "directed": info.directed,
        "weight": float(weight),
        "stoichiometry": stoichiometry,
        "ec": ec,
        "source_db": info.source_db,
        "evidence": evidence or info.source_db,
        "confidence": confidence,
    }


def to_simple_graph(graph: nx.Graph) -> nx.Graph:
    """Project a typed :class:`networkx.MultiDiGraph` to a simple undirected graph.

    Centrality, community detection and random-walk algorithms generally expect
    a simple graph.  This collapses parallel edges (keeping the largest weight)
    and drops direction, while preserving all node attributes.

    Parameters
    ----------
    graph : networkx.Graph
        Typically the ``MultiDiGraph`` returned by
        :meth:`Transnet.generate_graph`.

    Returns
    -------
    networkx.Graph
        Undirected simple graph with ``weight`` and a ``edge_types`` list
        recording which relationships were collapsed into each edge.
    """
    simple = nx.Graph()
    simple.add_nodes_from(graph.nodes(data=True))

    for u, v, data in graph.edges(data=True):
        weight = float(data.get("weight", 1.0) or 1.0)
        etype = data.get("edge_type", "unknown")
        if simple.has_edge(u, v):
            existing = simple[u][v]
            existing["weight"] = max(existing["weight"], weight)
            if etype not in existing["edge_types"]:
                existing["edge_types"].append(etype)
        else:
            simple.add_edge(u, v, weight=weight, edge_types=[etype])

    return simple


class Transnet:
    """
    Transnet class for connecting trans-omic networks.

    Parameters
    ----------
    name : str
        Label for this network.
    pathways, transcriptome, proteome, metabolome, reactions, signaling
        Layer objects.  All are optional; builders skip layers that are
        ``None``.  ``signaling`` supplies the phosphoproteome / kinase layer at
        the top of the trans-omic hierarchy and is only needed for studies that
        measured it.
    """

    def __init__(
            self,
            name: str = "Transnet",
            pathways: Pathways = None,
            transcriptome: Transcriptome = None,
            proteome: Proteome = None,
            metabolome: Metabolome = None,
            reactions: Reactions = None,
            signaling: Signaling = None,
            ):
        self.name = name
        self.pathways = pathways
        self.transcriptome = transcriptome
        self.proteome = proteome
        self.metabolome = metabolome
        self.reactions = reactions
        self.signaling = signaling
        self.graph = None

        self.cross_layer_edges = pd.DataFrame()
        self.edge_types = {}
        self._edge_index_built = False

    def __repr__(self):
        return f"<Transnet {self.name}>"

    #: Attributes added after the first release, with the value a network that
    #: predates them should take.
    _DEFAULTS = {
        "name": "Transnet",
        "pathways": None,
        "transcriptome": None,
        "proteome": None,
        "metabolome": None,
        "reactions": None,
        "signaling": None,
        "graph": None,
        "edge_types": None,
        "_edge_index_built": False,
    }

    def __setstate__(self, state):
        """Restore a pickled network, filling in attributes it predates.

        Unpickling bypasses ``__init__``, so a network pickled before a layer
        existed comes back without that attribute and every builder that checks
        for it raises. Backfilling here means an old cache keeps working: its
        layers are intact, so calling ``generate_graph()`` regenerates the graph
        under the current schema.
        """
        self.__dict__.update(state)
        for attribute, default in self._DEFAULTS.items():
            if not hasattr(self, attribute):
                setattr(self, attribute, default)
        if getattr(self, "edge_types", None) is None:
            self.edge_types = {}
        if not hasattr(self, "cross_layer_edges"):
            self.cross_layer_edges = pd.DataFrame()

    def layer_map(self) -> Dict[str, Any]:
        """Return ``{layer_name: layer_object}`` for the layers that are present."""
        candidates = {
            "Signaling": self.signaling,
            "Transcriptome": self.transcriptome,
            "Proteome": self.proteome,
            "Metabolome": self.metabolome,
            "Reactions": self.reactions,
            "Pathways": self.pathways,
        }
        return {name: layer for name, layer in candidates.items() if layer is not None}

    def protein_metabolite_interaction(self, weight: float = 1.0):
        """Enzyme -> metabolite edges (``enzymatic``).

        A reaction-free shortcut linking an enzyme to the compounds its EC
        number acts on.  Useful when no Reactions layer is available; when one
        is, prefer the explicit ``catalysis`` / ``substrate`` / ``product``
        edges built by :meth:`enzyme_reaction_interaction` and
        :meth:`reaction_metabolite_interaction`.

        Returns
        -------
        list of dict
            Canonical edge records.
        """
        edges = []

        if not self.proteome or not self.metabolome:
            logger.info("Skipping enzymatic edges: Proteome or Metabolome layer absent")
            return edges

        present_metabolites = {
            str(m.kegg_compound_id) for m in self.metabolome.metabolites
        }

        for protein in self.proteome.proteins:
            ec = ",".join(_as_str_list(protein.ec_number))
            for metabolite in _as_str_list(protein.metabolites):
                if metabolite in present_metabolites:
                    edges.append(_edge(
                        protein.uniprot_id, metabolite, "enzymatic",
                        weight=weight, ec=ec,
                        evidence=f"EC:{ec}" if ec else "KEGG",
                    ))

        logger.info(f"Created {len(edges)} enzymatic protein-metabolite edges")
        return edges

    def protein_protein_interaction(self, weight: float = 1.0):
        """Protein <-> protein edges (``protein_interaction``, undirected).

        One of the five trans-omic connection technologies.  Sourced from
        STRING; the association score is carried on the edge as ``confidence``
        where the layer recorded it.

        Returns
        -------
        list of dict
        """
        edges = []

        if not self.proteome:
            logger.info("Skipping protein-protein edges: Proteome layer absent")
            return edges

        present_proteins = set()
        for protein in self.proteome.proteins:
            present_proteins.update(_as_str_list(protein.uniprot_id))

        scores = getattr(self.proteome, "interaction_scores", None) or {}

        for protein in self.proteome.proteins:
            source = str(protein.uniprot_id)
            for partner in _as_str_list(protein.interaction_partners):
                if partner in present_proteins:
                    edges.append(_edge(
                        source, partner, "protein_interaction",
                        weight=weight,
                        confidence=scores.get((source, partner)),
                        evidence="STRING",
                    ))

        logger.info(f"Created {len(edges)} protein-protein edges")
        return edges

    def gene_protein_interaction(self, weight: float = 1.0):
        """Gene -> protein edges (``translation``, sign +1).

        Matched on NCBI/Entrez identifiers, normalised so that ``17309.0`` and
        ``17309`` join correctly.

        Returns
        -------
        list of dict
        """
        edges = []

        if not self.transcriptome or not self.proteome:
            logger.info("Skipping translation edges: Transcriptome or Proteome layer absent")
            return edges

        def _norm_ncbi(val):
            """Normalise an NCBI/Entrez ID to a plain integer string, e.g. 17309.0 -> '17309'."""
            try:
                return str(int(float(val)))
            except (TypeError, ValueError):
                return None

        # Build a set of normalised NCBI IDs present in the transcriptome,
        # mapping normalised ID -> gene identifier used as graph node
        ncbi_to_gene_node = {}
        for gene in self.transcriptome.genes:
            norm = _norm_ncbi(gene.ncbi_id)
            if norm:
                ncbi_to_gene_node[norm] = gene.ncbi_id

        # Gene symbol -> gene node, only where the symbol is unambiguous.
        symbol_to_nodes: Dict[str, set] = {}
        for gene in self.transcriptome.genes:
            node = gene.ncbi_id or gene.ensembl_id or gene.name
            for symbol in _as_str_list(gene.name):
                if symbol and node:
                    symbol_to_nodes.setdefault(symbol, set()).add(str(node))

        by_symbol = 0
        for protein in self.proteome.proteins:
            entrez_ids = protein.entrez_id if isinstance(protein.entrez_id, (list, np.ndarray)) else []
            linked = False
            for raw_id in entrez_ids:
                norm = _norm_ncbi(raw_id)
                if norm and norm in ncbi_to_gene_node:
                    edges.append(_edge(
                        ncbi_to_gene_node[norm], protein.uniprot_id, "translation",
                        weight=weight, evidence="NCBI/UniProt",
                    ))
                    linked = True
            if linked:
                continue
            # About 5% of proteins get no Entrez id from mygene -- in mouse,
            # pyruvate kinase L (Pklr) and pyruvate carboxylase (Pcx) among
            # them -- and were silently left without a gene. Their UniProt
            # gene symbol still names the gene; use it when it points to
            # exactly one gene node.
            for symbol in _as_str_list(protein.gene):
                nodes = symbol_to_nodes.get(symbol, set())
                if len(nodes) == 1:
                    edges.append(_edge(
                        next(iter(nodes)), protein.uniprot_id, "translation",
                        weight=weight, evidence="gene symbol",
                    ))
                    by_symbol += 1
                    break

        logger.info(
            f"Created {len(edges)} translation edges"
            + (f" ({by_symbol} by gene symbol, for proteins with no Entrez id)"
               if by_symbol else "")
        )
        return edges

    def transcription_factor_interaction(self, weight: float = 1.0):
        """Transcription factor -> gene edges (``transcriptional_regulation``).

        Sourced from ChIP-Atlas.  Binding evidence does not tell us whether the
        factor activates or represses its target, so ``sign`` is 0 (unknown);
        downstream sign propagation treats these edges as sign-preserving but
        flags the path as carrying an unsigned step.

        Returns
        -------
        list of dict
        """
        edges = []

        if not self.transcriptome or not self.proteome:
            logger.info(
                "Skipping transcriptional regulation edges: "
                "Transcriptome or Proteome layer absent"
            )
            return edges

        # ChIP-Atlas reports targets by gene symbol, so match on gene name and
        # emit the identifier actually used as the graph node.
        name_to_node = {}
        for gene in self.transcriptome.genes:
            if gene.name:
                name_to_node[str(gene.name)] = gene.ncbi_id or gene.ensembl_id or gene.name

        for protein in self.proteome.proteins:
            # Mean ChIP-Atlas binding score per target, when the proteome was
            # populated by get_transcription_factor_targets. Older pickles and
            # hand-built proteomes have no scores; the edge then carries none.
            scores = getattr(protein, "transcription_factor_target_scores", None) or {}
            for target in _as_str_list(protein.transcription_factor_targets):
                if target in name_to_node:
                    edges.append(_edge(
                        protein.uniprot_id, name_to_node[target],
                        "transcriptional_regulation",
                        weight=weight, evidence="ChIP-Atlas",
                        confidence=scores.get(target),
                    ))

        logger.info(f"Created {len(edges)} transcriptional regulation edges")
        return edges

    def _ec_to_uniprot(self) -> Dict[str, List[str]]:
        """Map each EC number in the proteome to the UniProt accessions carrying it."""
        ec_uniprot: Dict[str, List[str]] = {}
        if not self.proteome:
            return ec_uniprot
        for enzyme in self.proteome.proteins:
            for ec in _as_str_list(enzyme.ec_number):
                if ec:
                    ec_uniprot.setdefault(ec, []).append(str(enzyme.uniprot_id))
        return ec_uniprot

    def _catalysed_reactions(self) -> Dict[str, List[str]]:
        """Map reaction id -> EC numbers with a known catalyst in this organism.

        Catalysts come from the Proteome where there is one, and otherwise from
        EC numbers annotated on the transcriptome, so a study without
        proteomics still gets an organism-specific reaction set rather than all
        of KEGG.
        """
        known_ecs = set(self._ec_to_uniprot())
        if self.transcriptome:
            for gene in self.transcriptome.genes:
                known_ecs.update(_as_str_list(getattr(gene, "related_ecs", None)))

        catalysed: Dict[str, List[str]] = {}
        if not self.reactions:
            return catalysed
        for reaction in self.reactions.reactions:
            hits = [ec for ec in _as_str_list(reaction.enzyme) if ec in known_ecs]
            if hits:
                catalysed[str(reaction.id)] = hits
        return catalysed

    def enzyme_reaction_interaction(self, weight: float = 1.0):
        """Enzyme -> reaction edges (``catalysis``, sign +1).

        Joined on EC number, which is carried onto the edge in the ``ec``
        field rather than being discarded after the join.  This is the
        "metabolic regulation" connection technology of Yugi & Kuroda (2016).

        Returns
        -------
        list of dict
        """
        edges = []

        if not self.proteome or not self.reactions:
            logger.info("Skipping catalysis edges: Proteome or Reactions layer absent")
            return edges

        ec_uniprot = self._ec_to_uniprot()

        for reaction in self.reactions.reactions:
            for ec in _as_str_list(reaction.enzyme):
                for uniprot_id in ec_uniprot.get(ec, []):
                    edges.append(_edge(
                        uniprot_id, reaction.id, "catalysis",
                        weight=weight, ec=ec, evidence=f"KEGG EC:{ec}",
                    ))

        logger.info(f"Created {len(edges)} catalysis edges")
        return edges

    def gene_reaction_interaction(self, weight: float = 1.0):
        """Gene -> reaction edges (``gene_catalysis``, sign +1).

        Joined on the EC numbers annotated on each gene.  Only built when there
        is no Proteome layer: with proteomics the enzyme itself is the better
        evidence, and the ``translation`` + ``catalysis`` chain says so
        explicitly.  Without it, this is what keeps a transcriptome +
        metabolome study connected to its reaction layer instead of leaving the
        transcriptome as an island.

        Returns
        -------
        list of dict
        """
        edges = []

        if not self.transcriptome or not self.reactions:
            return edges

        ec_genes: Dict[str, List[str]] = {}
        for gene in self.transcriptome.genes:
            node = gene.ncbi_id or gene.ensembl_id or gene.name
            for ec in _as_str_list(getattr(gene, "related_ecs", None)):
                if ec:
                    ec_genes.setdefault(ec, []).append(str(node))

        if not ec_genes:
            logger.info(
                "No gene carries an EC annotation; skipping gene->reaction edges. "
                "Populate the transcriptome with KEGG EC links to connect the "
                "transcriptome to the reaction layer without proteomics."
            )
            return edges

        for reaction in self.reactions.reactions:
            for ec in _as_str_list(reaction.enzyme):
                for gene_node in ec_genes.get(ec, []):
                    edges.append(_edge(
                        gene_node, reaction.id, "gene_catalysis",
                        weight=weight, ec=ec, evidence=f"KEGG EC:{ec}",
                    ))

        logger.info(f"Created {len(edges)} gene-reaction edges")
        return edges

    def reaction_metabolite_interaction(
        self,
        weight: float = 1.0,
        require_enzyme: bool = True,
    ):
        """Metabolite -> reaction (``substrate``) and reaction -> metabolite (``product``) edges.

        Substrate and product roles are kept distinct -- they are what makes a
        reaction node directional, and what lets the metabolite regulation axis
        be read off the network.  Stoichiometric coefficients from the KEGG
        equation are carried on the edge.

        Metabolites are addressed by their KEGG compound identifier
        (``C#####``), the same identifier used everywhere else in the network,
        so the reaction layer connects to the metabolome rather than forming a
        parallel set of name-keyed nodes.

        Parameters
        ----------
        weight : float
            Edge weight.
        require_enzyme : bool
            When True (default) only reactions with at least one catalysing
            enzyme present in the Proteome contribute edges, which keeps the
            network organism-specific.  Set False to include every reaction in
            the Reactions layer -- the right choice when no Proteome is
            available.

        Returns
        -------
        list of dict
        """
        edges = []

        if not self.reactions or not self.metabolome:
            logger.info(
                "Skipping substrate/product edges: Reactions or Metabolome layer absent"
            )
            return edges

        present_metabolites = {
            str(m.kegg_compound_id) for m in self.metabolome.metabolites
        }

        catalysed = self._catalysed_reactions() if require_enzyme else None
        if require_enzyme and not catalysed:
            logger.info(
                "No reaction has a catalysing enzyme in the Proteome; "
                "building substrate/product edges for all reactions instead"
            )
            catalysed = None

        def _pairs(ids, stoichiometry):
            """Zip compound ids with their coefficients, tolerating length mismatch."""
            coefficients = _as_str_list(stoichiometry)
            for index, compound in enumerate(_as_str_list(ids)):
                coefficient = None
                if index < len(coefficients):
                    try:
                        coefficient = float(coefficients[index])
                    except (TypeError, ValueError):
                        coefficient = None  # variable stoichiometry, e.g. 'n'
                yield compound, coefficient

        for reaction in self.reactions.reactions:
            reaction_id = str(reaction.id)
            if catalysed is not None and reaction_id not in catalysed:
                continue

            ec = ",".join(_as_str_list(reaction.enzyme))

            for compound, coefficient in _pairs(
                reaction.substrates, reaction.stoichiometry_substrates
            ):
                if compound in present_metabolites:
                    edges.append(_edge(
                        compound, reaction_id, "substrate",
                        weight=weight, stoichiometry=coefficient, ec=ec,
                        evidence="KEGG reaction equation",
                    ))

            for compound, coefficient in _pairs(
                reaction.products, reaction.stoichiometry_products
            ):
                if compound in present_metabolites:
                    edges.append(_edge(
                        reaction_id, compound, "product",
                        weight=weight, stoichiometry=coefficient, ec=ec,
                        evidence="KEGG reaction equation",
                    ))

            # A reversible reaction can run either way, so each participant is
            # both consumed and produced.  Reversibility is recorded on the
            # reaction node rather than duplicating edges here; analyses that
            # care read it from the node attribute.

        logger.info(f"Created {len(edges)} substrate/product edges")
        return edges

    def _metabolite_name_index(self):
        """Return a callable resolving a compound id *or name* to a KEGG id.

        Databases disagree about how to name a compound: KEGG lists synonyms
        separated by ``;``, BRENDA reports a single common name. This indexes
        every synonym the metabolome knows, case- and whitespace-insensitively,
        so an annotation keyed by name still lands on the right node.
        """
        by_id = {}
        by_name = {}

        def _normalise(value):
            return str(value).strip().lower()

        for metabolite in self.metabolome.metabolites:
            compound = str(metabolite.kegg_compound_id)
            by_id[compound] = compound
            for synonym in _as_str_list(metabolite.kegg_name):
                for part in synonym.split(";"):
                    part = _normalise(part)
                    # First synonym wins, so a common name is not overwritten
                    # by an obscure one from another compound.
                    if part and part not in by_name:
                        by_name[part] = compound

        def resolve(value):
            if value is None:
                return None
            text = str(value).strip()
            if text in by_id:
                return by_id[text]
            for candidate in _name_variants(text):
                found = by_name.get(_normalise(candidate))
                if found:
                    return found
            return None

        return resolve

    #: Below this many organism-catalysed reactions the endogenous-effector
    #: filter is skipped: too small a reaction set to decide what is foreign.
    MIN_REACTIONS_FOR_ENDOGENOUS_FILTER = 200

    def allosteric_interactions(self, weight: float = 1.0, endogenous_only: bool = False):
        """Metabolite -> reaction allosteric edges (``allosteric_activation`` / ``allosteric_inhibition``).

        Allosteric regulation is the connection technology that makes a
        trans-omic network more than a pathway map: a metabolite modulates a
        reaction without being consumed by it, and its sign is known.  Sourced
        from BRENDA, which reports activators and inhibitors per EC number.

        The edge runs metabolite -> reaction, matching the direction of the
        regulatory effect.

        Parameters
        ----------
        weight : float
        endogenous_only : bool
            Off by default: BRENDA's effectors enter the network as BRENDA
            reports them. When True, keep only effectors whose KEGG compound is
            a substrate or product of a reaction this organism catalyses (EC
            number present in its proteome or transcriptome).

            BRENDA records every compound shown to change an enzyme's activity
            in vitro, so it lists drugs, flavonoids and assay reagents -- on
            mouse pyruvate kinase, luteolin and curcumin next to
            fructose-1,6-bisphosphate. The rule removes most of those, but it
            is a proxy and it misclassifies: on the mouse network it also
            removes physiological ions (Zn2+, Ca2+, Mg2+, Mn2+, K+), heparin,
            and deoxycholate -- an endogenous bile acid that metabolomics
            measures -- and it keeps any drug a mouse enzyme metabolises.

            In practice it changes little: the metabolite axis and
            :func:`~transnet.metabolite_regulatory_roles` only use metabolites
            that were *measured*, and assay reagents are not. Use it to clean
            figures and edge counts, not to change a result.

        Returns
        -------
        list of dict

        References
        ----------
        Kokaji T, et al. Transomics analysis reveals allosteric and gene
        regulation axes for altered hepatic glucose-responsive metabolism in
        obesity. *Science Signaling* 13:eaaz1236, 2020.
        """
        edges = []

        if not self.proteome or not self.metabolome or not self.reactions:
            logger.info(
                "Skipping allosteric edges: Proteome, Metabolome or Reactions layer absent"
            )
            return edges

        # BRENDA reports effectors by compound *name* ("ATP", "citrate"), while
        # the metabolome is keyed by KEGG compound id. Resolving names here is
        # what lets BRENDA annotations reach the graph at all; matching the two
        # directly would silently produce no allosteric edges.
        resolver = self._metabolite_name_index()

        # EC number -> reactions it catalyses
        ec_reactions: Dict[str, List[str]] = {}
        for reaction in self.reactions.reactions:
            for ec in _as_str_list(reaction.enzyme):
                if ec:
                    ec_reactions.setdefault(ec, []).append(str(reaction.id))

        # "Endogenous" means a substrate or product of a reaction this organism
        # can catalyse -- not of any KEGG reaction, because KEGG includes plant
        # flavonoid biosynthesis, so kaempferol and luteolin would pass. The
        # filter needs a representative reaction set to judge against; on a
        # handful of reactions, ATP would look foreign simply because no
        # reaction present uses it.
        endogenous = set()
        catalysed = self._catalysed_reactions() if endogenous_only else {}
        if endogenous_only and len(catalysed) >= self.MIN_REACTIONS_FOR_ENDOGENOUS_FILTER:
            for reaction in self.reactions.reactions:
                if str(reaction.id) in catalysed:
                    endogenous.update(_as_str_list(reaction.substrates))
                    endogenous.update(_as_str_list(reaction.products))
        elif endogenous_only:
            logger.info(
                f"Only {len(catalysed)} organism-catalysed reactions; too few to "
                f"judge which BRENDA effectors are endogenous, so all are kept"
            )

        unresolved = set()
        xenobiotic = set()
        for enzyme in self.proteome.proteins:
            regulators = [
                ("allosteric_activation", _as_str_list(enzyme.activators)),
                ("allosteric_inhibition", _as_str_list(enzyme.inhibitors)),
            ]
            for edge_type, metabolites in regulators:
                if not metabolites:
                    continue
                for ec in _as_str_list(enzyme.ec_number):
                    for reaction_id in ec_reactions.get(ec, []):
                        for metabolite in metabolites:
                            compound = resolver(metabolite)
                            if compound is None:
                                unresolved.add(metabolite)
                                continue
                            if endogenous and compound not in endogenous:
                                xenobiotic.add(metabolite)
                                continue
                            edges.append(_edge(
                                compound, reaction_id, edge_type,
                                weight=weight, ec=ec,
                                evidence=f"BRENDA:{metabolite}",
                            ))

        if unresolved:
            logger.info(
                f"{len(unresolved)} BRENDA effector name(s) matched no metabolite "
                f"in this network and were skipped, e.g. "
                f"{sorted(unresolved)[:5]}"
            )
        if xenobiotic:
            logger.info(
                f"{len(xenobiotic)} BRENDA effector(s) are not metabolites of this "
                f"network (drugs, chelators, signalling proteins) and were "
                f"skipped, e.g. {sorted(xenobiotic)[:6]}; pass "
                f"endogenous_only=False to keep them"
            )
        logger.info(f"Created {len(edges)} allosteric regulation edges")
        return edges

    #: Retained for backwards compatibility; prefer
    #: :meth:`allosteric_interactions`, whose name says what the edges are.
    activator_inhibitor_interactions = allosteric_interactions

    def signaling_interactions(self, weight: float = 1.0):
        """Signaling -> protein edges (``phosphorylation``, ``kinase_tf``).

        The top of the trans-omic hierarchy.  Returns an empty list when no
        Signaling layer is present, which is the common case -- a network
        without phosphoproteomics is complete in its own right, it simply
        starts one layer lower.

        Returns
        -------
        list of dict
        """
        edges = []

        if not self.signaling:
            return edges

        present_proteins = set()
        if self.proteome:
            for protein in self.proteome.proteins:
                present_proteins.update(_as_str_list(protein.uniprot_id))

        for node in self.signaling.nodes:
            source = str(node.id)
            for target in _as_str_list(node.substrates):
                if not present_proteins or target in present_proteins:
                    edges.append(_edge(
                        source, target, "phosphorylation",
                        weight=weight, sign=node.sign,
                        evidence=node.evidence or "KEGG",
                    ))
            for target in _as_str_list(node.tf_targets):
                if not present_proteins or target in present_proteins:
                    edges.append(_edge(
                        source, target, "kinase_tf",
                        weight=weight, sign=node.sign,
                        evidence=node.evidence or "KEGG",
                    ))

        logger.info(f"Created {len(edges)} signaling edges")
        return edges

    def generate_interaction_df(self, require_enzyme: bool = True) -> pd.DataFrame:
        """Assemble every edge in the network into one canonical DataFrame.

        This is the single source of truth for the network's edges: the graph,
        the adjacency matrix, the saved CSV and the cross-layer edge index are
        all derived from it, so an edge type cannot exist in one view and be
        missing from another.

        Parameters
        ----------
        require_enzyme : bool
            Passed through to :meth:`reaction_metabolite_interaction`.

        Returns
        -------
        pandas.DataFrame
            Columns as in
            :data:`transnet.biology.schema.INTERACTION_COLUMNS`.
        """
        weight = 1.0
        records: List[Dict[str, Any]] = []

        records.extend(self.signaling_interactions(weight))
        records.extend(self.transcription_factor_interaction(weight))
        records.extend(self.gene_protein_interaction(weight))
        records.extend(self.protein_protein_interaction(weight))
        records.extend(self.enzyme_reaction_interaction(weight))
        records.extend(self.reaction_metabolite_interaction(weight, require_enzyme=require_enzyme))
        records.extend(self.allosteric_interactions(weight))

        # Without proteomics the transcriptome would be an island, so link
        # genes straight to the reactions their EC numbers catalyse.
        if not self.proteome:
            records.extend(self.gene_reaction_interaction(weight))

        # The reaction-free enzyme->metabolite shortcut would duplicate the
        # catalysis/substrate/product chain, so only include it when there is no
        # Reactions layer to carry that structure.
        if not self.reactions:
            records.extend(self.protein_metabolite_interaction(weight))

        if not records:
            logger.warning("No interactions generated; are any layers populated?")
            return pd.DataFrame(columns=INTERACTION_COLUMNS)

        interaction_df = pd.DataFrame(records, columns=INTERACTION_COLUMNS)
        interaction_df = interaction_df.drop_duplicates(
            subset=["source", "target", "edge_type"]
        ).reset_index(drop=True)

        logger.info(
            f"Generated interaction dataframe with {len(interaction_df)} edges "
            f"across {interaction_df['edge_type'].nunique()} edge types"
        )
        return interaction_df

    def _node_metadata(self) -> Dict[str, Dict[str, Any]]:
        """Collect node attributes (layer, type, name, reversibility) from the layers."""
        metadata: Dict[str, Dict[str, Any]] = {}

        def _put(node_id, **attrs):
            if node_id is None or str(node_id) == "nan":
                return
            metadata.setdefault(str(node_id), {}).update(attrs)

        if self.transcriptome:
            for gene in self.transcriptome.genes:
                _put(gene.ncbi_id or gene.ensembl_id,
                     layer="Transcriptome", node_type="Gene", name=gene.name)
        if self.proteome:
            for protein in self.proteome.proteins:
                # UniProt names are full descriptions ("Large ribosomal
                # subunit protein uL2"); the gene symbol is what a reader
                # recognises in a figure, so keep it alongside.
                symbols = _as_str_list(protein.gene)
                _put(protein.uniprot_id,
                     layer="Proteome", node_type="Protein", name=protein.name,
                     symbol=symbols[0] if symbols else "",
                     ec=",".join(_as_str_list(protein.ec_number)))
        if self.metabolome:
            for metabolite in self.metabolome.metabolites:
                _put(metabolite.kegg_compound_id,
                     layer="Metabolome", node_type="Metabolite",
                     name=metabolite.kegg_name)
        if self.reactions:
            for reaction in self.reactions.reactions:
                _put(reaction.id,
                     layer="Reactions", node_type="Reaction", name=reaction.name,
                     reversible=getattr(reaction, "reversible", True),
                     ec=",".join(_as_str_list(reaction.enzyme)))
        if self.signaling:
            for node in self.signaling.nodes:
                _put(node.id, layer="Signaling", node_type="SignalingProtein",
                     name=node.name)

        return metadata

    def generate_graph(self, directed: bool = True, require_enzyme: bool = True):
        """Build the trans-omic graph.

        Parameters
        ----------
        directed : bool
            When True (default) returns a :class:`networkx.MultiDiGraph` that
            preserves edge direction and allows two molecules to be connected by
            more than one relationship (e.g. a metabolite that is both a
            substrate of a reaction and an allosteric inhibitor of it).  Set
            False for an undirected simple graph -- equivalent to calling
            :func:`to_simple_graph` on the directed result.
        require_enzyme : bool
            Passed through to :meth:`generate_interaction_df`.

        Returns
        -------
        networkx.MultiDiGraph or networkx.Graph
            Nodes carry ``layer``, ``node_type`` and ``name``; edges carry the
            full schema (``edge_type``, ``role``, ``sign``, ``weight``,
            ``stoichiometry``, ``ec``, ``evidence``, ``confidence``).
        """
        df = self.generate_interaction_df(require_enzyme=require_enzyme)
        metadata = self._node_metadata()

        G = nx.MultiDiGraph()

        # Layer assignment comes from the layer objects where possible, so a
        # node's layer does not depend on which edge row happened to mention it
        # last.  Edge-derived layers are the fallback for nodes we have no
        # element object for.
        for _, row in df.iterrows():
            for node, layer in (
                (row["source"], row["source_layer"]),
                (row["target"], row["target_layer"]),
            ):
                attrs = metadata.get(str(node))
                if attrs is None:
                    attrs = {"layer": layer, "node_type": layer, "name": str(node)}
                    metadata[str(node)] = attrs
                G.add_node(str(node), **attrs)

        for _, row in df.iterrows():
            G.add_edge(
                str(row["source"]),
                str(row["target"]),
                key=row["edge_type"],
                **{col: row[col] for col in INTERACTION_COLUMNS
                   if col not in ("source", "target", "source_layer", "target_layer")},
            )

        if not directed:
            G = to_simple_graph(G)

        self.graph = G
        logger.info(
            f"Generated {'directed multigraph' if directed else 'simple graph'} "
            f"with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges"
        )
        return G

    def generate_adjacency_matrix(self, symmetric: bool = True) -> pd.DataFrame:
        """Adjacency matrix of the network.

        Parameters
        ----------
        symmetric : bool
            When True (default) the matrix is symmetrised, matching the
            undirected projection used by most centrality algorithms.  Set
            False to keep edge direction, in which case ``A[i, j]`` is the
            weight of the edge from ``i`` to ``j``.
        """
        interaction_df = self.generate_interaction_df()
        if interaction_df.empty:
            return pd.DataFrame()

        elements = sorted(
            set(interaction_df["source"]).union(interaction_df["target"])
        )
        adjacency = pd.DataFrame(0.0, index=elements, columns=elements)

        for _, row in interaction_df.iterrows():
            adjacency.loc[row["source"], row["target"]] = row["weight"]
            if symmetric:
                adjacency.loc[row["target"], row["source"]] = row["weight"]

        logger.info(f"Generated adjacency matrix of size {adjacency.shape}")
        return adjacency

    def save_network(self, output_dir: str, save_adjacency: bool = False):
        """Write the network to ``interactions.csv`` and ``nodes.csv``.

        Parameters
        ----------
        output_dir : str
            Directory to write into; created if missing.
        save_adjacency : bool
            Also write a dense ``adjacency_matrix.csv``.  Off by default: the
            matrix is quadratic in node count and unusable at genome scale.
        """
        import os

        os.makedirs(output_dir, exist_ok=True)

        df = self.generate_interaction_df()
        df.to_csv(os.path.join(output_dir, "interactions.csv"), index=False)

        node_data = []
        if self.transcriptome and self.transcriptome.genes:
            for gene in self.transcriptome.genes:
                node_data.append({
                    "ID": gene.ncbi_id or gene.ensembl_id,
                    "Name": gene.name, "Symbol": gene.name, "Type": "Gene",
                    "Layer": "Transcriptome",
                    "FC": gene.fc, "P_value": gene.adj_p_value, "Reversible": None,
                })
        if self.proteome and self.proteome.proteins:
            for protein in self.proteome.proteins:
                # UniProt's name is a full description; the gene symbol is what
                # a reader recognises, and without this column a network read
                # back from CSV labels every protein with the description.
                symbols = _as_str_list(protein.gene)
                node_data.append({
                    "ID": protein.uniprot_id,
                    "Name": protein.name,
                    "Symbol": symbols[0] if symbols else None,
                    "Type": "Protein", "Layer": "Proteome",
                    "FC": protein.fc, "P_value": protein.adj_p_value, "Reversible": None,
                })
        if self.metabolome and self.metabolome.metabolites:
            for metabolite in self.metabolome.metabolites:
                node_data.append({
                    "ID": metabolite.kegg_compound_id,
                    "Name": metabolite.kegg_name, "Type": "Metabolite",
                    "Layer": "Metabolome",
                    "FC": metabolite.fc, "P_value": metabolite.adj_p_value,
                    "Reversible": None,
                })
        if self.reactions and self.reactions.reactions:
            for reaction in self.reactions.reactions:
                node_data.append({
                    "ID": reaction.id, "Name": reaction.name, "Type": "Reaction",
                    "Layer": "Reactions", "FC": None, "P_value": None,
                    "Reversible": getattr(reaction, "reversible", True),
                })
        if self.signaling and self.signaling.nodes:
            for node in self.signaling.nodes:
                node_data.append({
                    "ID": node.id, "Name": node.name, "Type": "SignalingProtein",
                    "Layer": "Signaling", "FC": node.fc, "P_value": node.adj_p_value,
                    "Reversible": None,
                })

        pd.DataFrame(node_data).to_csv(
            os.path.join(output_dir, "nodes.csv"), index=False
        )

        if save_adjacency:
            self.generate_adjacency_matrix().to_csv(
                os.path.join(output_dir, "adjacency_matrix.csv")
            )

        logger.info(f"Saved network to {output_dir}")

    @classmethod
    def load_network(cls, input_dir: str):
        """Load a network previously written by :meth:`save_network`.

        Parameters
        ----------
        input_dir : str
            Directory containing ``interactions.csv`` and ``nodes.csv``.

        Returns
        -------
        Transnet
        """
        import os

        transnet = cls()

        df = pd.read_csv(os.path.join(input_dir, "interactions.csv"))
        node_df = pd.read_csv(os.path.join(input_dir, "nodes.csv"))

        layers = set()
        for col in ("source_layer", "target_layer"):
            if col in df.columns:
                layers.update(df[col].dropna().tolist())
        if "Layer" in node_df.columns:
            layers.update(node_df["Layer"].dropna().tolist())

        if "Transcriptome" in layers:
            transnet.transcriptome = Transcriptome()
            transnet.transcriptome.genes = []
        if "Proteome" in layers:
            transnet.proteome = Proteome()
            transnet.proteome.proteins = []
        if "Metabolome" in layers:
            transnet.metabolome = Metabolome()
            transnet.metabolome.metabolites = []
        if "Reactions" in layers:
            transnet.reactions = Reactions()
            transnet.reactions.reactions = []
        if "Signaling" in layers:
            transnet.signaling = Signaling()
            transnet.signaling.nodes = []

        def _value(row, column):
            return row[column] if column in row.index else None

        for _, row in node_df.iterrows():
            node_type = row.get("Type")
            if node_type == "Gene":
                transnet.transcriptome.genes.append(Gene(
                    ncbi_id=row["ID"], name=row["Name"],
                    fc=_value(row, "FC"), adj_p_value=_value(row, "P_value"),
                ))
            elif node_type == "Protein":
                transnet.proteome.proteins.append(Protein(
                    uniprot_id=row["ID"], name=row["Name"],
                    fc=_value(row, "FC"), adj_p_value=_value(row, "P_value"),
                ))
            elif node_type == "Metabolite":
                transnet.metabolome.metabolites.append(Metabolite(
                    kegg_compound_id=row["ID"], kegg_name=row["Name"],
                    fc=_value(row, "FC"), adj_p_value=_value(row, "P_value"),
                ))
            elif node_type == "Reaction":
                reversible = _value(row, "Reversible")
                transnet.reactions.reactions.append(Reaction(
                    id=row["ID"], name=row["Name"],
                    reversible=True if pd.isna(reversible) else bool(reversible),
                ))
            elif node_type == "SignalingProtein":
                transnet.signaling.nodes.append(SignalingProtein(
                    id=row["ID"], name=row["Name"],
                    fc=_value(row, "FC"), adj_p_value=_value(row, "P_value"),
                ))

        transnet._loaded_interactions = df
        logger.info(f"Loaded network from {input_dir} ({len(df)} edges)")
        return transnet

    def build_cross_layer_edge_index(self, require_enzyme: bool = True) -> pd.DataFrame:
        """Build the queryable edge index used by the graph-query helpers.

        Derived from :meth:`generate_interaction_df`, so it can never disagree
        with the graph about which edges exist -- the previous implementation
        rebuilt the edges independently and silently omitted allosteric
        regulation.

        Returns
        -------
        pandas.DataFrame
            The canonical edge table; ``source``, ``target``, ``source_layer``,
            ``target_layer``, ``edge_type``, ``weight`` and ``evidence`` are
            always present.
        """
        logger.info("Building cross-layer edge index")

        self.cross_layer_edges = self.generate_interaction_df(
            require_enzyme=require_enzyme
        )

        if not self.cross_layer_edges.empty:
            self.edge_types = {
                edge_type: group.copy()
                for edge_type, group in self.cross_layer_edges.groupby("edge_type")
            }
        else:
            self.edge_types = {}

        self._edge_index_built = True

        logger.info(
            f"Built cross-layer edge index with {len(self.cross_layer_edges)} edges "
            f"across {len(self.edge_types)} edge types"
        )
        return self.cross_layer_edges

    def get_neighbors(
        self, 
        node_id: str, 
        layer: Optional[str] = None, 
        edge_type: Optional[str] = None,
        direction: str = 'both'
    ) -> List[str]:
        """
        Get all neighbors of a node, optionally filtered by layer/type.
        
        Parameters:
        -----------
        node_id : str
            Node identifier
        layer : str, optional
            Filter neighbors by layer (e.g., 'Proteome', 'Metabolome')
        edge_type : str, optional
            Filter by edge type (e.g., 'enzymatic', 'translation')
        direction : str
            Direction of edges: 'outgoing', 'incoming', or 'both' (default)
        
        Returns:
        --------
        List[str]
            List of neighbor node IDs
        """
        if not self._edge_index_built:
            logger.warning("Edge index not built. Building now...")
            self.build_cross_layer_edge_index()
        
        if self.cross_layer_edges.empty:
            logger.warning("No edges in cross-layer edge index")
            return []
        
        neighbors = []
        
        # Get outgoing edges
        if direction in ['outgoing', 'both']:
            edges_from = self.cross_layer_edges[
                self.cross_layer_edges['source'] == str(node_id)
            ]
            
            if layer:
                edges_from = edges_from[edges_from['target_layer'] == layer]
            
            if edge_type:
                edges_from = edges_from[edges_from['edge_type'] == edge_type]
            
            neighbors.extend(edges_from['target'].tolist())
        
        # Get incoming edges
        if direction in ['incoming', 'both']:
            edges_to = self.cross_layer_edges[
                self.cross_layer_edges['target'] == str(node_id)
            ]
            
            if layer:
                edges_to = edges_to[edges_to['source_layer'] == layer]
            
            if edge_type:
                edges_to = edges_to[edges_to['edge_type'] == edge_type]
            
            neighbors.extend(edges_to['source'].tolist())
        
        return list(set(neighbors))
    
    def find_paths(
        self, 
        source: str, 
        target: str, 
        max_length: int = 3,
        allowed_layers: Optional[List[str]] = None,
        allowed_edge_types: Optional[List[str]] = None
    ) -> List[List[str]]:
        """
        Find all paths between two nodes up to max_length.
        
        This is the core method for cross-layer discovery. Given a changed
        gene and a changed metabolite, find mechanistic paths connecting them.
        
        Parameters:
        -----------
        source : str
            Source node ID
        target : str
            Target node ID
        max_length : int
            Maximum path length (number of edges)
        allowed_layers : List[str], optional
            Only use nodes from these layers
        allowed_edge_types : List[str], optional
            Only use edges of these types
        
        Returns:
        --------
        List[List[str]]
            List of paths, where each path is a list of node IDs
        
        Examples:
        ---------
        # Find paths from a gene to a metabolite
        >>> paths = transnet.find_paths('ENSG00000123456', 'C00002', max_length=4)
        >>> for path in paths:
        ...     print(' -> '.join(path))
        """
        if not self._edge_index_built:
            logger.warning("Edge index not built. Building now...")
            self.build_cross_layer_edge_index()
        
        # Build filtered graph if needed
        edges_to_use = self.cross_layer_edges.copy()
        
        if allowed_edge_types:
            edges_to_use = edges_to_use[
                edges_to_use['edge_type'].isin(allowed_edge_types)
            ]
        
        if allowed_layers:
            edges_to_use = edges_to_use[
                edges_to_use['source_layer'].isin(allowed_layers) &
                edges_to_use['target_layer'].isin(allowed_layers)
            ]
        
        # Build NetworkX graph
        G = nx.DiGraph()
        
        for _, row in edges_to_use.iterrows():
            G.add_edge(
                row['source'], 
                row['target'],
                weight=row['weight'],
                edge_type=row['edge_type']
            )
        
        # Find paths
        try:
            paths = list(nx.all_simple_paths(
                G, 
                str(source), 
                str(target), 
                cutoff=max_length
            ))
            logger.info(f"Found {len(paths)} paths from {source} to {target}")
            return paths
        except (nx.NodeNotFound, nx.NetworkXNoPath) as e:
            logger.warning(f"Could not find paths: {e}")
            return []
    
    def get_path_annotations(
        self, 
        path: List[str]
    ) -> List[Dict[str, Any]]:
        """
        Get annotations for edges in a path.
        
        Parameters:
        -----------
        path : List[str]
            List of node IDs representing a path
        
        Returns:
        --------
        List[Dict[str, Any]]
            List of edge annotations, one for each edge in the path
        """
        if not self._edge_index_built:
            self.build_cross_layer_edge_index()
        
        annotations = []
        
        for i in range(len(path) - 1):
            source = str(path[i])
            target = str(path[i + 1])
            
            # Find edge(s) between these nodes
            edge_matches = self.cross_layer_edges[
                (self.cross_layer_edges['source'] == source) &
                (self.cross_layer_edges['target'] == target)
            ]
            
            if not edge_matches.empty:
                # Take first match if multiple
                edge_data = edge_matches.iloc[0].to_dict()
                annotations.append(edge_data)
            else:
                # Edge not found (shouldn't happen if path is valid)
                annotations.append({
                    'source': source,
                    'target': target,
                    'edge_type': 'unknown',
                    'weight': 0.0,
                    'evidence': 'not_found'
                })
        
        return annotations
    
    def query_cross_layer_relationships(
        self,
        changed_nodes: Dict[str, List[str]],
        p_value_threshold: float = 0.05,
        path_max_length: int = 3
    ) -> pd.DataFrame:
        """
        Query cross-layer relationships between changed nodes.
        
        This is the main entry point for cross-layer discovery analysis.
        Given sets of changed nodes from different omics layers, find
        mechanistic paths connecting them.
        
        Parameters:
        -----------
        changed_nodes : Dict[str, List[str]]
            Dictionary mapping layer names to lists of changed node IDs,
            e.g. ``{'Transcriptome': ['gene1'], 'Metabolome': ['met1']}``
        p_value_threshold : float
            P-value threshold for considering a node as changed
        path_max_length : int
            Maximum path length to search
        
        Returns:
        --------
        pd.DataFrame
            DataFrame with columns: source, target, path_length, path, 
            source_layer, target_layer, edge_types
        """
        if not self._edge_index_built:
            self.build_cross_layer_edge_index()
        
        logger.info(f"Querying cross-layer relationships among {sum(len(v) for v in changed_nodes.values())} changed nodes")
        
        results = []
        
        # For each pair of layers
        layer_pairs = [
            (l1, l2) for l1 in changed_nodes.keys() 
            for l2 in changed_nodes.keys() 
            if l1 != l2
        ]
        
        for source_layer, target_layer in layer_pairs:
            source_nodes = changed_nodes[source_layer]
            target_nodes = changed_nodes[target_layer]
            
            logger.info(f"Finding paths from {source_layer} to {target_layer}")
            
            # Find paths between all pairs
            for source in source_nodes:
                for target in target_nodes:
                    paths = self.find_paths(
                        source, 
                        target, 
                        max_length=path_max_length
                    )
                    
                    for path in paths:
                        # Get edge types in path
                        annotations = self.get_path_annotations(path)
                        edge_types = [a['edge_type'] for a in annotations]
                        
                        results.append({
                            'source': source,
                            'target': target,
                            'source_layer': source_layer,
                            'target_layer': target_layer,
                            'path_length': len(path) - 1,
                            'path': ' -> '.join(path),
                            'edge_types': ', '.join(edge_types),
                            'evidence': ', '.join([a['evidence'] for a in annotations])
                        })
        
        results_df = pd.DataFrame(results)

        if not results_df.empty:
            # Sort by path length
            results_df = results_df.sort_values('path_length').reset_index(drop=True)
            logger.info(f"Found {len(results_df)} cross-layer paths")
        else:
            logger.warning("No cross-layer paths found")

        return results_df

    def enrich_with_brenda(
        self,
        organism: Optional[str] = None,
        email: Optional[str] = None,
        password: Optional[str] = None,
        fields: Optional[List[str]] = None,
    ) -> None:
        """Enrich the Proteome layer with kinetic data from BRENDA.

        Delegates to :func:`~transnet.api.brenda.brenda_enrich_proteins`.
        After calling this method, each :class:`~transnet.biology.elements.Protein`
        that has EC numbers will have its ``activators``, ``inhibitors``,
        ``substrates``, ``products``, and (optionally) ``cofactors`` lists
        populated with compound names from BRENDA.

        These data are then used by :meth:`activator_inhibitor_interactions`
        when building the full multi-layer graph.

        Parameters
        ----------
        organism : str, optional
            Restrict BRENDA queries to this organism (e.g. ``"Homo sapiens"``).
        email, password : str, optional
            BRENDA credentials.  Fall back to ``BRENDA_EMAIL`` /
            ``BRENDA_PASSWORD`` environment variables or an interactive prompt.
        fields : list of str, optional
            Subset of ``['inhibitors', 'activators', 'substrates', 'products',
            'cofactors']``.  All fields are fetched by default.
        """
        if self.proteome is None or not self.proteome.proteins:
            logger.warning(
                "enrich_with_brenda: Proteome layer is empty. "
                "Call proteome.populate() before enriching with BRENDA."
            )
            return

        logger.info(
            f"Enriching proteome ({len(self.proteome.proteins)} proteins) with BRENDA"
        )
        from ..api.brenda import brenda_enrich_proteins
        brenda_enrich_proteins(
            proteins=self.proteome.proteins,
            organism=organism,
            email=email,
            password=password,
            fields=fields,
        )
