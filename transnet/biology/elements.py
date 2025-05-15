"""Biological elements (genes, proteins, etc.) that form the basis of a biological network
"""

from typing import List, Dict, Any, Optional, Union

__all__ = [
    "Reaction",
    "Pathway",
    "Metabolite",
    "Gene",
    "Protein",
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
        self.activators = []
        self.inhibitors = []
        self.substrates = []
        self.products = []
        self.metabolites = []
        
    def __repr__(self):
        return f"<Protein {self.uniprot_id}: {self.name}>"
        