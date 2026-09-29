"""Biological elements (genes, proteins, etc.) that form the basis of a biological network
"""

from typing import List, Dict, Any, Optional, Union

__all__ = [
    "Reaction",
    "Pathway",
    "Metabolite",
    "Gene",
    "Protein",
    "SignalingProtein",
]

class Reaction:
    """
    Reaction class representing a biochemical reaction.
    """
    
    def __init__(
            self,
            id: str = None,
            name: str = None,
            equation: str = None,
            definition: str = None,
            enzyme: str = None,
            substrates: List[str] = None,
            products: List[str] = None,
            stoichiometry_substrates: List[float] = None,
            stoichiometry_products: List[float] = None,
            reversible: bool = True,
            ):
        self.id = id
        self.name = name
        self.equation = equation
        self.definition = definition
        self.enzyme = enzyme
        self.substrates = substrates if substrates is not None else []
        self.products = products if products is not None else []
        self.stoichiometry_substrates = stoichiometry_substrates if stoichiometry_substrates is not None else []
        self.stoichiometry_products = stoichiometry_products if stoichiometry_products is not None else []
        #: Whether the reaction can run in both directions.  KEGG writes most
        #: equations with a reversible arrow, so this defaults to True; an
        #: irreversible arrow in the source equation sets it to False.
        self.reversible = reversible

    def __repr__(self):
        return f"<Reaction {self.id}: {self.name}>"


class Pathway:
    """
    Pathway class representing a biological pathway.
    """
    
    def __init__(
            self,
            id: str = None,
            name: str = None,
            kegg_organism: str = None,
            genes: List[str] = None,
            ecs: List[str] = None,
            ):
        self.id = id
        self.name = name
        self.kegg_organism = kegg_organism
        self.genes = genes if genes is not None else []
        self.ecs = ecs if ecs is not None else []
        
    def __repr__(self):
        return f"<Pathway {self.id}: {self.name}>"


class Metabolite:
    """
    Metabolite class representing a small molecule.
    """
    
    def __init__(
            self,
            pubchem_id: str = None,
            kegg_name: str = None,
            kegg_compound_id: str = None,
            inchi: str = None,
            inchikey: str = None,
            chebi_id: str = None,
            smile: str = None,
            fc: float = None,
            adj_p_value: float = None,
            data: Any = None,
            ):
        self.pubchem_id = pubchem_id
        self.kegg_name = kegg_name
        self.kegg_compound_id = kegg_compound_id
        self.inchi = inchi
        self.inchikey = inchikey
        self.chebi_id = chebi_id
        self.smile = smile
        self.fc = fc
        self.adj_p_value = adj_p_value
        self.data = data
        
    def __repr__(self):
        return f"<Metabolite {self.kegg_compound_id}: {self.kegg_name}>"


class Gene:
    """
    Gene class representing a gene or transcript.
    """
    
    def __init__(
            self,
            kegg_id: str = None,
            ncbi_id: str = None,
            ensembl_id: str = None,
            transcript_id: str = None,
            uniprot_id: str = None,
            type: str = None,
            name: str = None,
            description: str = None,
            kegg_organism: str = None,
            fc: float = None,
            adj_p_value: float = None,
            data: Any = None,
            ):
        self.kegg_id = kegg_id
        self.ncbi_id = ncbi_id
        self.ensembl_id = ensembl_id
        self.transcript_id = transcript_id
        self.uniprot_id = uniprot_id
        self.type = type
        self.name = name
        self.description = description
        self.kegg_organism = kegg_organism
        self.fc = fc
        self.adj_p_value = adj_p_value
        self.data = data
        self.related_ecs = []
        
    def __repr__(self):
        return f"<Gene {self.kegg_id or self.ensembl_id}: {self.name}>"


class Protein:
    """
    Protein class representing a protein.
    """
    
    def __init__(
            self,
            uniprot_id: str = None,
            uniprot_name: str = None,
            review_status: str = None,
            name: str = None,
            gene: str = None,
            organism_full: str = None,
            length: int = None,
            ec_number: List[str] = None,
            ensembl_id: List[str] = None,
            entrez_id: List[str] = None,
            ncbi_organism: str = None,
            fc: float = None,
            adj_p_value: float = None,
            data: Any = None,
            ):
        self.uniprot_id = uniprot_id
        self.uniprot_name = uniprot_name
        self.review_status = review_status
        self.name = name
        self.gene = gene
        self.organism_full = organism_full
        self.length = length
        self.ec_number = ec_number if ec_number is not None else []
        self.ensembl_id = ensembl_id if ensembl_id is not None else []
        self.entrez_id = entrez_id if entrez_id is not None else []
        self.ncbi_organism = ncbi_organism
        self.fc = fc
        self.adj_p_value = adj_p_value
        self.data = data
        
        # Interactions and related elements
        self.interaction_partners = []
        self.transcription_factor_targets = []
        # {target gene symbol: mean ChIP-Atlas binding score}; populated
        # alongside transcription_factor_targets and carried onto the
        # transcriptional_regulation edges as `confidence`.
        self.transcription_factor_target_scores = {}
        self.activators = []
        self.inhibitors = []
        self.substrates = []
        self.products = []
        self.metabolites = []
        
    def __repr__(self):
        return f"<Protein {self.uniprot_id}: {self.name}>"


class SignalingProtein:
    """A node in the signaling layer: a kinase, phosphatase or phosphosite.

    The signaling layer sits at the top of the trans-omic hierarchy and is
    optional -- many studies do not measure the phosphoproteome.  When it is
    present it supplies the ``phosphorylation`` and ``kinase_tf`` edges that
    let regulatory paths be traced from a stimulus down to metabolites.

    Parameters
    ----------
    id : str
        Node identifier, normally a UniProt accession, optionally suffixed with
        a phosphosite (e.g. ``"P31749_S473"``).
    name : str
        Display name.
    uniprot_id : str
        Accession of the parent protein, used to link a phosphosite back to the
        Proteome layer.
    site : str
        Modified residue, e.g. ``"S473"``.
    sign : int
        Direction of this node's regulatory effect on its targets: ``+1``
        activating, ``-1`` inhibiting, ``0`` unknown.
    """

    def __init__(
            self,
            id: str = None,
            name: str = None,
            uniprot_id: str = None,
            site: str = None,
            sign: int = 0,
            pathway: str = None,
            evidence: str = None,
            fc: float = None,
            adj_p_value: float = None,
            data: Any = None,
            ):
        self.id = id
        self.name = name
        self.uniprot_id = uniprot_id
        self.site = site
        self.sign = sign
        self.pathway = pathway
        self.evidence = evidence
        self.fc = fc
        self.adj_p_value = adj_p_value
        self.data = data

        #: Proteins this node phosphorylates.
        self.substrates = []
        #: Transcription factors this node regulates.
        self.tf_targets = []

    def __repr__(self):
        return f"<SignalingProtein {self.id}: {self.name}>"
