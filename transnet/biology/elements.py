"""Biological elements, nodes in the network
"""

__all__ = [
    "Metabolite",
    "Reaction",
    "Protein",
    "Gene",
    "Pathway",
]


class Metabolite:
    """
    Metabolite
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

    def __repr__(self):
        return f"<Metabolite: {self.kegg_name}>"


class Reaction:
    """
    Reaction
    """

    def __init__(
            self,
            id: str = None,
            name: str = None,
            equation: str = None,
            definition: str = None,
            enzyme: str = None,
            substrates: list = [],
            products: list = [],
            stoichiometry_substrates: list = [],
            stoichiometry_products: list = [],
            ):
        self.id = id
        self.name = name
        self.equation = equation
        self.definition = definition
        self.enzyme = enzyme
        self.substrates = substrates
        self.products = products
        self.stoichiometry_substrates = stoichiometry_substrates
        self.stoichiometry_products = stoichiometry_products

    def __repr__(self):
        return f"<Reaction: {self.name}, ID: {self.id}>"


class Protein:
    """
    Protein
    """

    def __init__(
            self,
            ensembl_id: str = None,
            uniprot_id: str = None,
            uniprot_name: str = None,
            review_status: str = None,
            entrez_id: str = None,
            name: str = None,
            gene: str = None,
            organism_full: str = None,
            length: str = None,
            ncbi_organism: str = None,
            ec_number: list = [],
            pathways: str = None,
            reactions: str = None,
            substrates: list = [],
            products: list = [],
            interaction_partners: list = [],
            metabolites: list = [],
            transcription_factor_substrates: list = [],
            fc: float = None,
            adj_p_value: float = None,
            ):
        self.ensembl_id = ensembl_id
        self.uniprot_id = uniprot_id
        self.uniprot_name = uniprot_name
        self.review_status = review_status
        self.entrez_id = entrez_id
        self.name = name
        self.gene = gene
        self.organism_full = organism_full
        self.length = length
        self.ncbi_organism = ncbi_organism
        self.ec_number = ec_number
        self.pathways = pathways
        self.reactions = reactions
        self.substrates = substrates
        self.products = products
        self.interaction_partners = interaction_partners
        self.metabolites = metabolites
        self.transcription_factor_substrates = transcription_factor_substrates
        self.fc = fc
        self.adj_p_value = adj_p_value

    def __repr__(self):
        return f"<Protein: {self.name}, ID: {self.uniprot_id}>"


class Gene:
    """
    Gene
    """

    def __init__(
            self,
            name: str = None,
            ncbi_id: str = None,
            kegg_id: str = None,
            ensembl_id: str = None,
            transcript_id: str = None,
            uniprot_id: str = None,
            type: str = None,
            related_ecs: list = [],
            description: str = None,
            kegg_organism: str = None,
            fc: float = None,
            adj_p_value: float = None,
            ):
        self.name = name
        self.ncbi_id = ncbi_id
        self.kegg_id = kegg_id
        self.ensembl_id = ensembl_id
        self.transcript_id = transcript_id
        self.uniprot_id = uniprot_id
        self.type = type
        self.related_ecs = related_ecs
        self.description = description
        self.kegg_organism = kegg_organism
        self.fc = fc
        self.adj_p_value = adj_p_value

    def __repr__(self):
        if isinstance(self.name, list):
            return f"<Gene: {self.name[0]}, ID: {self.kegg_id}>"
        else:
            return f"<Gene: {self.name}, ID: {self.kegg_id}>"


class Pathway:
    """
    Pathway
    """

    def __init__(
            self,
            name: str = None,
            id: str = None,
            genes: list = [],
            ecs: list = [],
            kegg_organism: str = None,
            ):
        self.name = name
        self.id = id
        self.genes = genes
        self.ecs = ecs
        self.kegg_organism = kegg_organism

    def __repr__(self):
        return f"<Pathway: {self.name}, ID: {self.id}>"
