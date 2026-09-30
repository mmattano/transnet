"""Biological layers (omics), layers in the network
"""

from typing import List, Dict, Any, Optional, Union
import pandas as pd
import logging

logger = logging.getLogger(__name__)

__all__ = [
    "Reactions",
    "Pathways",
    "Transcriptome",
    "Proteome",
    "Metabolome",
    "Signaling",
]

# Import elements after defining __all__ to avoid circular imports
from transnet.biology.elements import (
    Reaction, Pathway, Metabolite, Gene, Protein, SignalingProtein,
)
from transnet.api.kegg import (
    kegg_create_reaction_table, kegg_list_pathways, kegg_link_pathway, 
    kegg_link_ec, kegg_list_genes, kegg_conv_ncbi_idtable, kegg_ec_to_cpds, 
    kegg_list_compounds, kegg_to_chebi, chebi_to_kegg
)
from transnet.api.uniprot import uniprot_list_proteins, uniprot_add_entrez_id
from transnet.api.string import string_map_identifiers, string_get_interactions, translate_string_dict
from transnet.api.ensembl import ensembl_download_release, ensembl_get_transcripts, ensembl_expander
from transnet.api.chip_atlas import get_ChIP_data, get_ChIP_exps
from transnet.api.chem_info import chemical_info_converter

class OmicsLayer:
    """Base class for all omics layers."""
    
    def __init__(self, name: str = None):
        self.name = name
        self.elements = []
        
    def add_experimental_data(self,
                              input_data: pd.DataFrame = None,
                              id_column: str = None,
                              id_type: str = None):
        """Add experimental data to this layer."""
        raise NotImplementedError(
            f"{self.__class__.__name__} does not take experimental data; "
            "map measurements onto the built graph with map_omics_to_network")

    def populate(self):
        """Populate this layer with data from databases."""
        raise NotImplementedError(f"{self.__class__.__name__} has no populate method")
        
    def __repr__(self):
        return f"<{self.__class__.__name__}: {len(self.elements)} elements>"


class Reactions(OmicsLayer):
    """
    Reactions layer representing biochemical reactions.
    """
    
    def __init__(self):
        super().__init__(name="Reactions")
        self.reactions = []
    
    def __repr__(self):
        return f"<Reactions: {len(self.reactions)} reactions in total>"
    
    def populate(
            self,
            from_api: bool = True,
            df = None,
            ):
        """
        Populate reactions from KEGG API or from a dataframe.
        
        Parameters
        ----------
        from_api : bool
            Whether to populate reactions from KEGG API
        df : pd.DataFrame
            Dataframe containing reaction information
        """
        if from_api:
            logger.info("Populating reactions from KEGG API")
            reaction_table = kegg_create_reaction_table()
            for _, row in reaction_table.iterrows():
                new_reaction = Reaction(
                    id = row['reaction'],
                    name = row['name'],
                    equation = row['equation'],
                    definition = row['definition'],
                    enzyme = row['enzyme'],
                    substrates = row['substrates'],
                    products = row['products'],
                    stoichiometry_substrates = row['stoichiometry_substrates'],
                    stoichiometry_products = row['stoichiometry_products'],
                )
                self.reactions.append(new_reaction)
            logger.info(f"Populated {len(self.reactions)} reactions from KEGG API")
        elif df is not None:
            logger.info("Populating reactions from dataframe")
            reaction_table = df
            for _, row in reaction_table.iterrows():
                new_reaction = Reaction(
                    id = row['reaction'],
                    name = row['name'],
                    equation = row['equation'],
                    definition = row['definition'],
                    enzyme = row['enzyme'],
                    substrates = row['substrates'],
                    products = row['products'],
                    stoichiometry_substrates = row['stoichiometry_substrates'],
                    stoichiometry_products = row['stoichiometry_products'],
                    # Tables written before reversibility was recorded default
                    # to True, matching KEGG's usual reversible arrow.
                    reversible = bool(row['reversible'])
                        if 'reversible' in row.index and pd.notna(row['reversible'])
                        else True,
                )
                self.reactions.append(new_reaction)
            logger.info(f"Populated {len(self.reactions)} reactions from dataframe")
        else:
            logger.warning("No data source provided for reaction population")


class Pathways(OmicsLayer):
    """
    Pathways layer representing biological pathways.
    """

    def __init__(self):
        super().__init__(name="Pathways")
        self.kegg_organism = None
        self.pathways = []
        self.organism_full = None

    def __repr__(self):
        return f"<Pathways of: {self.kegg_organism}, {len(self.pathways)} pathways in total>"

    def populate(self):
        """
        Populate pathways from KEGG API.
        """
        if not self.kegg_organism:
            logger.error("KEGG organism code must be set before populating pathways")
            return
            
        logger.info(f"Populating pathways for {self.kegg_organism} from KEGG API")
        for _, row in kegg_list_pathways(self.kegg_organism).iterrows():
            new_pathway = Pathway(
                id = row[0],
                name = row[1],
                kegg_organism = self.kegg_organism,
            )
            self.pathways.append(new_pathway)
        logger.info(f"Populated {len(self.pathways)} pathways from KEGG API")

    def fill_pathways(self):
        """
        Fill in additional pathway information from KEGG.
        """
        if not self.kegg_organism:
            logger.error("KEGG organism code must be set before filling pathways")
            return
            
        logger.info(f"Filling pathway information for {self.kegg_organism}")
        kegg_pathway_links = kegg_link_pathway(self.kegg_organism)
        kegg_pathway_ecs = kegg_link_ec(self.kegg_organism)
        kegg_pathway_list = kegg_list_pathways(self.kegg_organism)

        for pathway in self.pathways:
            # Checks which info is missing and fills it
            # based on what is available
            # Start with KEGG ID
            if pathway.id is not None:
                # Check if name is there
                if pathway.name is None:
                    try:
                        pathway.name = kegg_pathway_list[
                            "description"
                        ][
                            kegg_pathway_list["pathways_id"]
                            == pathway.id
                        ].to_string(
                            index=False
                        )
                    except Exception as e:
                        logger.warning(f"Error getting pathway name for {pathway.id}: {e}")
                
                # Check if genes are there
                if len(pathway.genes) == 0:
                    try:
                        pathway.genes = kegg_pathway_links[
                            kegg_pathway_links['pathway'] == pathway.id][
                            "kegg_gene_id"
                        ].to_list()
                    except Exception as e:
                        logger.warning(f"Error getting pathway genes for {pathway.id}: {e}")
                
                # Check if EC numbers are there
                if len(pathway.ecs) == 0:
                    protein_ecs = []
                    for gene in kegg_pathway_links[
                        kegg_pathway_links['pathway'] == pathway.id][
                        "kegg_gene_id"
                    ].to_list():
                        if gene in kegg_pathway_ecs["kegg_gene_id"].to_list():
                            try:
                                ecs_tmp = kegg_pathway_ecs[
                                    kegg_pathway_ecs["kegg_gene_id"] == gene
                                ]["ec_number"]
                                if len(ecs_tmp) == 1:
                                    protein_ecs.append(ecs_tmp.to_string(index=False))
                                else:
                                    for ec in ecs_tmp.to_list():
                                        protein_ecs.append(ec)
                            except Exception as e:
                                logger.warning(f"Error getting EC numbers for {gene}: {e}")
                    pathway.ecs = list(set(protein_ecs))
        
        logger.info(f"Successfully filled pathway information")

    def populate_from_reactome(self, species="9606"):
        """Populate pathways from Reactome instead of (or in addition to) KEGG.

        Queries the Reactome ContentService for all top-level pathways for
        the given species and appends them to ``self.pathways``.

        Parameters
        ----------
        species : str
            NCBI taxonomy ID (e.g. ``"9606"``) or Reactome species name.
        """
        from transnet.api.reactome import reactome_get_pathways
        logger.info(f"Populating pathways from Reactome for species {species}")
        try:
            df = reactome_get_pathways(species)
            if df.empty:
                logger.warning("Reactome returned no pathways")
                return
            for _, row in df.iterrows():
                new_pathway = Pathway(
                    id=row.get("stId"),
                    name=row.get("displayName") or row.get("name"),
                )
                self.pathways.append(new_pathway)
            logger.info(
                f"Added {len(df)} Reactome pathways "
                f"(total: {len(self.pathways)})"
            )
        except Exception as exc:
            logger.error(f"Reactome pathway population failed: {exc}")

    def over_representation_analysis(self, gene_ids, species="9606", p_value=0.05):
        """Run Reactome pathway over-representation analysis.

        Parameters
        ----------
        gene_ids : list
            Gene symbols, UniProt IDs, or Ensembl IDs.
        species : str
            NCBI taxonomy ID or Reactome species name.
        p_value : float
            P-value cut-off.

        Returns
        -------
        pd.DataFrame
            Enriched pathways with FDR-corrected p-values.
        """
        from transnet.api.reactome import reactome_over_representation
        logger.info(
            f"Running Reactome ORA for {len(gene_ids)} genes "
            f"(p<={p_value})"
        )
        try:
            return reactome_over_representation(
                gene_ids, species=species, p_value=p_value
            )
        except Exception as exc:
            logger.error(f"Reactome ORA failed: {exc}")
            return pd.DataFrame()


class Metabolome(OmicsLayer):
    """
    Metabolome layer representing metabolites.
    """

    def __init__(self):
        super().__init__(name="Metabolome")
        self.metabolites = []
        self.input_id_type = None
        self.input_ids = []
        self.input_data = None

    def __repr__(self):
        return f"<Metabolome, {len(self.metabolites)} metabolites in total>"
    
    def add_experimental_data(
        self,
        input_data: pd.DataFrame = None,
        metabolite_column_name: str = None,
        id_type: str = None,
    ):
        """
        Add experimental metabolomics data.
        
        Parameters
        ----------
        input_data : pd.DataFrame
            Dataframe containing experimental data
        metabolite_column_name : str
            Column name containing metabolite identifiers
        id_type : str
            Type of identifier (chebi, kegg, pubchem, inchi, inchikey, smile)
        """
        self.input_data = input_data
        self.input_ids = input_data[metabolite_column_name].tolist()
        
        valid_id_types = ["chebi", "kegg", "pubchem", "inchi", "inchikey", "smile"]
        if id_type in valid_id_types:
            self.input_id_type = id_type
        else:
            logger.warning(f"Unknown ID type: {id_type}. Valid types are: {valid_id_types}")

    def populate(
        self,
        metabolite_column_name: str = None,
        input_data_value_column_name: str = None,
        input_data_p_value_column_name: str = None,
    ):
        """
        Populate metabolome with all KEGG compounds.

        If experimental PubChem IDs have been registered via
        ``add_experimental_data(..., id_type='pubchem')``, those IDs are
        mapped to KEGG compound IDs through mychem.info (ChEBI) and the KEGG
        REST API, and the matching metabolite objects are flagged with
        ``data={'in_experiment': True}``.

        Parameters
        ----------
        metabolite_column_name : str
            Unused – kept for API compatibility.
        input_data_value_column_name : str
            Unused – kept for API compatibility.
        input_data_p_value_column_name : str
            Unused – kept for API compatibility.
        """
        logger.info("Populating metabolome layer from KEGG compound list")
        self.metabolites = []

        try:
            kegg_compounds_df = kegg_list_compounds()
        except Exception as e:
            logger.error(f"Error fetching KEGG compound list: {e}")
            logger.warning("Using empty metabolome")
            return

        # Map experimental PubChem IDs -> KEGG via mychem.info + chebi_to_kegg.
        # The PubChem -> KEGG correspondence is kept, not just a boolean flag:
        # an omics table keyed by PubChem needs it to reach the KEGG-keyed
        # nodes at all (see transnet.build_alias_id_map).
        exp_kegg_ids: set = set()
        kegg_to_pubchem: dict = {}
        if self.input_ids and self.input_id_type == "pubchem":
            input_ids = [str(p) for p in self.input_ids]
            logger.info(f"Mapping {len(input_ids)} experimental PubChem IDs to KEGG")
            try:
                chem_df = chemical_info_converter(input_ids)
                # Try every ChEBI id known for a compound, not just the first.
                # A compound's neutral form and its zwitterion are separate
                # ChEBI entries and KEGG cross-references only some of them,
                # so taking one at random loses ~20% of a typical panel.
                if "chebi_ids_all" in chem_df.columns:
                    candidates = [
                        list(c) if isinstance(c, (list, tuple)) else ([c] if c else [])
                        for c in chem_df["chebi_ids_all"]
                    ]
                else:  # older chemical_info_converter
                    candidates = [
                        [c] if c else [] for c in chem_df["chebi_ids"]
                    ]
                valid_chebi = sorted({c for group in candidates for c in group if c})
                if valid_chebi:
                    kegg_df = chebi_to_kegg(valid_chebi)
                    chebi_kegg_map = {
                        chebi: kegg
                        for chebi, kegg in zip(
                            kegg_df["chebi_compounds"], kegg_df["kegg_compounds"]
                        )
                        if kegg
                    }
                    for pubchem, group in zip(input_ids, candidates):
                        kid = next(
                            (chebi_kegg_map[c] for c in group if c in chebi_kegg_map),
                            None,
                        )
                        if kid:
                            exp_kegg_ids.add(kid)
                            kegg_to_pubchem.setdefault(kid, pubchem)
                logger.info(
                    f"Mapped {len(exp_kegg_ids)}/{len(input_ids)} experimental "
                    f"metabolites to KEGG IDs"
                )
                if not exp_kegg_ids:
                    logger.warning(
                        "No experimental metabolite mapped to a KEGG compound. "
                        "Check that the supplied IDs really are PubChem CIDs."
                    )
            except Exception as e:
                logger.error(f"PubChem->KEGG mapping failed: {e}")

        # Build Metabolite objects from full KEGG compound list
        for _, row in kegg_compounds_df.iterrows():
            kid = str(row.get("kegg_compounds", "")).strip()
            try:
                self.metabolites.append(Metabolite(
                    kegg_compound_id=kid,
                    kegg_name=str(row.get("name", "")),
                    pubchem_id=kegg_to_pubchem.get(kid),
                    data={"in_experiment": kid in exp_kegg_ids} if exp_kegg_ids else None,
                ))
            except Exception as e:
                logger.error(f"Error creating metabolite {kid}: {e}")

        logger.info(f"Populated metabolome with {len(self.metabolites)} metabolites"
                    + (f" ({len(exp_kegg_ids)} from experimental data)" if exp_kegg_ids else ""))

    def enrich_with_hmdb(self, fields=None):
        """Enrich metabolites with data from HMDB (Human Metabolome Database).

        For each metabolite in the layer, looks up the HMDB record using
        available identifiers (HMDB ID > ChEBI > KEGG > name) and populates
        missing fields such as InChI, InChIKey, SMILES, formula, molecular
        weight, pathway associations, and disease links.

        Parameters
        ----------
        fields : list of str, optional
            Subset of ``['inchi', 'inchikey', 'smiles', 'formula',
            'molecular_weight', 'chebi_id', 'pubchem_id',
            'pathways', 'diseases']``.
            All fields are enriched by default.
        """
        if not self.metabolites:
            logger.error("Metabolites must be populated before HMDB enrichment")
            return
        from transnet.api.hmdb import hmdb_enrich_metabolites
        logger.info(
            f"Enriching {len(self.metabolites)} metabolites with HMDB"
        )
        try:
            hmdb_enrich_metabolites(self.metabolites, fields=fields)
        except Exception as exc:
            logger.error(f"HMDB enrichment failed: {exc}")
            logger.warning("Continuing without HMDB data")


class Transcriptome(OmicsLayer):
    """
    Transcriptome layer representing genes and transcripts.
    """

    def __init__(self):
        super().__init__(name="Transcriptome")
        self.kegg_organism = None
        self.organism_full = None
        self.input_id_type = None
        self.input_ids = []
        self.input_data = None
        self.genes = []

    def __repr__(self):
        return f"<Transcriptome of: {self.kegg_organism}, {len(self.genes)} genes in total>"

    def add_experimental_data(
        self,
        input_data: pd.DataFrame = None,
        transcript_column_name: str = None,
        id_type: str = None,
    ):
        """
        Add experimental transcriptomics data.
        
        Parameters
        ----------
        input_data : pd.DataFrame
            Dataframe containing experimental data
        transcript_column_name : str
            Column name containing transcript identifiers
        id_type : str
            Type of identifier (ensembl, etc.)
        """
        self.input_data = input_data
        self.input_ids = input_data[transcript_column_name].tolist()
        if id_type == "ensembl":
            self.input_id_type = "ensembl"
        else:
            logger.warning(f"Unknown ID type: {id_type}. Currently only 'ensembl' is supported.")

    def populate(
        self,
        kegg_organism: str = None,
        organism_full: str = None,
        ensembl: bool = False,
        ensembl_release: int = 109,
        kegg_api: bool = False,
        biotype: str = None,
        transcript_column_name: str = None,
        input_data_value_column_name: str = None,
        input_data_p_value_column_name: str = None,
    ):
        """
        Populate transcriptome with data.
        
        Parameters
        ----------
        kegg_organism : str
            KEGG organism code
        organism_full : str
            Full organism name for Ensembl
        ensembl : bool
            Whether to use Ensembl as data source
        ensembl_release : int
            Ensembl release version
        biotype : str
            Filter genes by biotype
        transcript_column_name : str
            Column name containing transcript identifiers
        input_data_value_column_name : str
            Column name containing fold change values
        input_data_p_value_column_name : str
            Column name containing p-values
        """
        self.kegg_organism = kegg_organism or self.kegg_organism
        self.organism_full = organism_full or self.organism_full
        self.genes = []
        
        if not ensembl and not kegg_api:
            # Default to KEGG API if no source specified
            kegg_api = True
            
        if kegg_api:
            # If experimental ENSEMBL IDs are available, build from those with
            # mygene ENSEMBL→NCBI conversion.  Otherwise fall back to the full
            # KEGG gene list for the organism.
            if self.input_ids and self.input_id_type == "ensembl":
                logger.info(
                    f"Populating transcriptome from {len(self.input_ids)} experimental ENSEMBL IDs"
                )
                try:
                    import mygene as _mygene
                    mg = _mygene.MyGeneInfo()
                    mg_out = mg.querymany(
                        self.input_ids,
                        scopes="ensembl.gene",
                        fields=["entrezgene", "symbol"],
                        species={"mmu": "mouse", "hsa": "human", "rno": "rat",
                                 "sce": "yeast", "eco": "fruitfly"}.get(
                                    (self.kegg_organism or ""), "human"),
                        returnall=True,
                    )
                    ensembl_to_ncbi   = {}
                    ensembl_to_symbol = {}
                    for entry in mg_out.get("out", []):
                        eid = entry.get("query")
                        if "entrezgene" in entry:
                            try:
                                ensembl_to_ncbi[eid] = str(int(float(entry["entrezgene"])))
                            except (TypeError, ValueError):
                                pass
                        if "symbol" in entry:
                            ensembl_to_symbol[eid] = entry["symbol"]
                    for eid in self.input_ids:
                        try:
                            self.genes.append(Gene(
                                ensembl_id=eid,
                                ncbi_id=ensembl_to_ncbi.get(eid),
                                name=ensembl_to_symbol.get(eid, eid),
                                kegg_organism=self.kegg_organism,
                            ))
                        except Exception as e:
                            logger.error(f"Error creating gene {eid}: {e}")
                    n_ncbi = sum(1 for g in self.genes if g.ncbi_id)
                    logger.info(
                        f"Populated transcriptome with {len(self.genes)} genes "
                        f"({n_ncbi} with NCBI IDs)"
                    )
                except Exception as e:
                    logger.error(f"ENSEMBL→NCBI conversion failed: {e}")
                    logger.warning("Using empty transcriptome")
            else:
                logger.info(f"Populating genes from KEGG API for {self.kegg_organism}")
                try:
                    kegg_genes = kegg_list_genes(self.kegg_organism)
                    # Resolve NCBI/Entrez ids up front. Without them a gene
                    # cannot be joined to its protein, so the Transcriptome
                    # layer ends up with no `translation` edges at all -- the
                    # layer is present, connected downward by
                    # transcriptional_regulation, and silently missing its
                    # link up to the Proteome. One conversion call, then a
                    # dict lookup per gene.
                    kegg_to_ncbi = {}
                    try:
                        idtable = kegg_conv_ncbi_idtable(self.kegg_organism)
                        kegg_to_ncbi = dict(
                            zip(
                                idtable["kegg_id"].astype(str).str.strip(),
                                idtable["ncbi_id"].astype(str).str.strip(),
                            )
                        )
                        logger.info(
                            f"KEGG->NCBI: {len(kegg_to_ncbi):,} gene id mappings"
                        )
                    except Exception as e:
                        logger.warning(
                            f"KEGG->NCBI gene id conversion failed: {e}. "
                            f"Genes will have no NCBI id, so no translation "
                            f"edges will be built."
                        )

                    for _, row in kegg_genes.iterrows():
                        try:
                            self.genes.append(Gene(
                                kegg_id=row.get("gene_id"),
                                ncbi_id=kegg_to_ncbi.get(
                                    str(row.get("gene_id")).strip()
                                ),
                                name=row.get("name(s)"),
                                description=row.get("description"),
                                kegg_organism=self.kegg_organism,
                            ))
                        except Exception as e:
                            logger.error(f"Error creating gene from KEGG: {e}")
                    logger.info(f"Populated transcriptome with {len(self.genes)} genes from KEGG API")
                except Exception as e:
                    logger.error(f"Error getting genes from KEGG API: {e}")
                    logger.warning("Using empty transcriptome")
                
        elif ensembl:
            logger.info(f"Populating genes from Ensembl for {self.organism_full}")
            try:
                ensembl_download_release(
                    release_number=ensembl_release,
                    organism_full=self.organism_full,
                )
                transcripts_df = ensembl_get_transcripts(
                    release_number=ensembl_release,
                    organism_full=self.organism_full,
                )
                
                if biotype:
                    logger.info(f"Filtering transcripts by biotype: {biotype}")
                    if not isinstance(biotype, list):
                        biotype = [biotype]
                    transcripts_df = transcripts_df[
                        transcripts_df["biotype"].isin(biotype)
                    ]
                
                ensembl_df = ensembl_expander(transcripts_df)
                
                for _, row in ensembl_df.iterrows():
                    try:
                        fc = None
                        adj_p_value = None
                        
                        new_gene = Gene(
                            kegg_id=row.get("entrez_id"),
                            ncbi_id=row.get("entrez_id"),
                            ensembl_id=row.get("ensembl_gene_id"),
                            transcript_id=row.get("ensembl_transcript_id"),
                            uniprot_id=row.get("uniprot_id"),
                            type=row.get("biotype"),
                            name=row.get("gene_description"),
                            description=row.get("gene_name"),
                            kegg_organism=self.kegg_organism,
                            fc=fc,
                            adj_p_value=adj_p_value,
                        )
                        self.genes.append(new_gene)
                    except Exception as e:
                        logger.error(f"Error creating gene from Ensembl: {e}")
                        continue
                logger.info(f"Populated transcriptome with {len(self.genes)} genes from Ensembl")
            except Exception as e:
                logger.error(f"Error getting genes from Ensembl: {e}")
                logger.warning("Using empty transcriptome")

            

    def fill_gene_info(self):
        """
        Fill in additional gene information from KEGG.
        """
        if not self.kegg_organism:
            logger.error("KEGG organism code must be set before filling gene info")
            return
            
        logger.info(f"Filling gene information for {self.kegg_organism}")

        try:
            kegg_gene_list = kegg_list_genes(self.kegg_organism)
            kegg_ncbi_idtable = kegg_conv_ncbi_idtable(self.kegg_organism)
            kegg_ec_link = kegg_link_ec(self.kegg_organism)

            for gene in self.genes:
                # Start with KEGG ID
                if gene.kegg_id is not None:
                    # Check if name is there
                    if gene.name is None:
                        try:
                            kegg_info = kegg_gene_list[
                                kegg_gene_list["gene_id"] == gene.kegg_id
                            ]
                            if not kegg_info.empty:
                                gene.name = kegg_info["name(s)"].to_string(index=False)
                        except Exception as e:
                            logger.warning(f"Error getting gene name for {gene.kegg_id}: {e}")
                    
                    # Check if description is there
                    if gene.description is None:
                        try:
                            kegg_info = kegg_gene_list[
                                kegg_gene_list["gene_id"] == gene.kegg_id
                            ]
                            if not kegg_info.empty:
                                gene.description = kegg_info["description"].to_string(
                                    index=False
                                )
                        except Exception as e:
                            logger.warning(f"Error getting gene description for {gene.kegg_id}: {e}")
                    
                    # Check if ncbi_id is there
                    if gene.ncbi_id is None:
                        try:
                            ncbi_info = kegg_ncbi_idtable[
                                kegg_ncbi_idtable["kegg_id"] == gene.kegg_id
                            ]
                            if not ncbi_info.empty:
                                # .to_string() pads the value; the id is later
                                # compared as a string, so take the cell.
                                gene.ncbi_id = str(
                                    ncbi_info["ncbi_id"].iloc[0]
                                ).strip()
                        except Exception as e:
                            logger.warning(f"Error getting NCBI ID for {gene.kegg_id}: {e}")
                    
                    # Check if related EC numbers are there
                    if len(gene.related_ecs) == 0:
                        try:
                            protein_ecs = []
                            if gene.kegg_id in kegg_ec_link["kegg_gene_id"].to_list():
                                ecs_tmp = kegg_ec_link[
                                    kegg_ec_link["kegg_gene_id"] == gene.kegg_id
                                ]["ec_number"]
                                if len(ecs_tmp) == 1:
                                    protein_ecs.append(ecs_tmp.to_string(index=False))
                                else:
                                    for ec in ecs_tmp.to_list():
                                        protein_ecs.append(ec)
                            gene.related_ecs = list(set(protein_ecs))
                        except Exception as e:
                            logger.warning(f"Error getting EC numbers for {gene.kegg_id}: {e}")
            
            logger.info(f"Successfully filled gene information")
        except Exception as e:
            logger.error(f"Error filling gene information: {e}")


class Proteome(OmicsLayer):
    """
    Proteome layer representing proteins.
    """

    def __init__(self):
        super().__init__(name="Proteome")
        self.proteins = []
        self.organism_full = None
        self.input_id_type = None
        self.input_ids = []
        self.input_data = None
        self.ncbi_organism = None
        self.kegg_organism = None
        # Factors ChIP-Atlas could not be reached for, after retries. Empty
        # until get_transcription_factor_targets runs.
        self.chip_atlas_failed_factors = []

    def __repr__(self):
        return f"<Proteome of: {self.ncbi_organism or self.kegg_organism}, {len(self.proteins)} proteins in total>"

    def add_experimental_data(
        self,
        input_data: pd.DataFrame = None,
        protein_column_name: str = None,
        id_type: str = None,
    ):
        """
        Add experimental proteomics data.
        
        Parameters
        ----------
        input_data : pd.DataFrame
            Dataframe containing experimental data
        protein_column_name : str
            Column name containing protein identifiers
        id_type : str
            Type of identifier (uniprot, etc.)
        """
        self.input_data = input_data
        self.input_ids = input_data[protein_column_name].tolist()
        if id_type == "uniprot":
            self.input_id_type = "uniprot"
        else:
            logger.warning(f"Unknown ID type: {id_type}. Currently only 'uniprot' is supported.")

    def populate(
        self,
        kegg_organism: str = None,
        ncbi_organism: str = None,
        organism_full: str = None,
        ensembl: bool = False,
        uniprot: bool = True,
        ensembl_release: int = 109,
        protein_column_name: str = None,
        input_data_value_column_name: str = None,
        input_data_p_value_column_name: str = None,
        reviewed_only: bool = False,
        uniprot_organism: str = None,
    ):
        """
        Populate proteome with data.
        
        Parameters
        ----------
        kegg_organism : str
            KEGG organism code
        ncbi_organism : str
            NCBI organism code
        organism_full : str
            Full organism name for Ensembl
        ensembl : bool
            Whether to use Ensembl as data source
        uniprot : bool
            Whether to use UniProt as data source
        reviewed_only : bool
            Restrict to reviewed (Swiss-Prot) entries. Recommended for
            reference networks: the unreviewed TrEMBL bulk is mostly predicted
            isoforms, and it is what makes STRING and ChIP-Atlas lookups
            unreliable at this scale.
        uniprot_organism : str, optional
            Taxon to query UniProt with, when it differs from
            ``ncbi_organism``. UniProt files the curated entries of some
            model organisms under a reference *strain* rather than the
            species: yeast under S288c (559292, not 4932) and E. coli under
            K-12 (83333, not 511145). Queried at species level, reviewed-only
            returns 43 yeast proteins and 0 for E. coli. ``ncbi_organism``
            is left untouched because STRING uses the species-level ids.
        ensembl_release : int
            Ensembl release version
        protein_column_name : str
            Column name containing protein identifiers
        input_data_value_column_name : str
            Column name containing fold change values
        input_data_p_value_column_name : str
            Column name containing p-values
        """
        self.kegg_organism = kegg_organism or self.kegg_organism
        self.ncbi_organism = ncbi_organism or self.ncbi_organism
        self.organism_full = organism_full or self.organism_full
        
        if not ensembl and not uniprot:
            logger.warning("Neither Ensembl nor UniProt selected as data source. Using UniProt as default.")
            uniprot = True
            
        if uniprot:
            taxon = uniprot_organism or self.ncbi_organism
            logger.info(
                f"Populating proteins from UniProt for {taxon}"
                + (f" (strain taxon for {self.ncbi_organism})"
                   if uniprot_organism and uniprot_organism != self.ncbi_organism
                   else "")
            )
            try:
                uniprot_proteins = uniprot_list_proteins(
                    taxon, reviewed_only=reviewed_only
                )
                
                # Filter to experimental IDs if set, otherwise use full table
                if self.input_id_type == "uniprot" and self.input_ids:
                    input_id_set = set(self.input_ids)
                    uniprot_proteins = uniprot_proteins[
                        uniprot_proteins["Entry"].isin(input_id_set)
                    ]
                    logger.info(f"Filtered to {len(uniprot_proteins)} experimental proteins")

                # Add Entrez IDs
                uniprot_proteins = uniprot_add_entrez_id(uniprot_proteins)

                for _, row in uniprot_proteins.iterrows():
                    try:
                        # Parse EC numbers
                        ec_numbers = []
                        ec_raw = row.get("EC number")
                        try:
                            if ec_raw and not pd.isna(ec_raw):
                                ec_str = str(ec_raw).strip()
                                if "; " in ec_str:
                                    ec_numbers = [e.strip() for e in ec_str.split("; ") if e.strip()]
                                elif "  " in ec_str:
                                    ec_numbers = [e.strip() for e in ec_str.split("  ") if e.strip()]
                                else:
                                    ec_numbers = [ec_str]
                        except Exception as e:
                            logger.error(f"Error parsing EC numbers for {row.get('Entry')}: {e}")

                        # Normalise Entrez IDs to plain integer strings
                        raw_entrez = row.get("Entrez", [])
                        entrez_ids = []
                        if isinstance(raw_entrez, list):
                            for e in raw_entrez:
                                try:
                                    entrez_ids.append(str(int(float(e))))
                                except (TypeError, ValueError):
                                    pass

                        new_protein = Protein(
                            uniprot_id=row.get("Entry"),
                            uniprot_name=row.get("Entry Name"),
                            review_status=row.get("Reviewed"),
                            name=row.get("Protein names"),
                            gene=row.get("Gene Names"),
                            organism_full=row.get("Organism"),
                            length=row.get("Length"),
                            ec_number=ec_numbers,
                            ensembl_id=(
                                str(row.get("Ensembl")).split(";")
                                if row.get("Ensembl") and not pd.isna(row.get("Ensembl"))
                                else []
                            ),
                            entrez_id=entrez_ids,
                            ncbi_organism=self.ncbi_organism,
                        )
                        self.proteins.append(new_protein)
                    except Exception as e:
                        logger.error(f"Error creating protein from UniProt: {e}")
                        continue

                logger.info(f"Populated proteome with {len(self.proteins)} proteins from UniProt")
            except Exception as e:
                logger.error(f"Error getting proteins from UniProt: {e}")
                logger.warning("Using empty proteome")
        
        elif ensembl:
            logger.info("Populating proteins from Ensembl - limited functionality")
            # Simplified implementation for now
            logger.warning("Ensembl protein population not fully implemented")
            self.proteins = []

    def get_interaction_partners(self):
        """
        Get protein-protein interaction partners from STRING database.
        """
        if not self.ncbi_organism or not self.proteins:
            logger.error("NCBI organism code and proteins must be set before getting interactions")
            return
            
        logger.info(f"Getting protein-protein interactions from STRING for {self.ncbi_organism}")
        
        try:
            # Extract protein IDs
            proteins = [protein.uniprot_id for protein in self.proteins]
            proteins_flattened = []
            
            # Flatten list of protein IDs if some are lists
            for protein_entries in proteins:
                if isinstance(protein_entries, list):
                    for gene_entry in protein_entries:
                        proteins_flattened.append(gene_entry)
                else:
                    proteins_flattened.append(protein_entries)
            
            # STRING API accepts at most 2000 identifiers per request — batch.
            _BATCH = 2000

            def _map(identifiers, label):
                mapped = {}
                for _i in range(0, len(identifiers), _BATCH):
                    _batch = identifiers[_i: _i + _BATCH]
                    _m = string_map_identifiers(
                        protein_list=_batch,
                        species=self.ncbi_organism,
                    )
                    mapped.update(_m)
                    logger.info(
                        f"STRING mapping batch {_i // _BATCH + 1} ({label}): "
                        f"{len(_batch)} in → {len(_m)} mapped"
                    )
                return mapped

            string_mapping_dict = _map(proteins_flattened, "UniProt")

            # STRING resolves UniProt accessions for some organisms but not
            # others -- yeast and E. coli return nothing for them while
            # resolving gene symbols fine. Fall back rather than silently
            # producing a network with no protein interactions.
            symbol_to_uniprot = {}
            if not string_mapping_dict:
                logger.info(
                    f"STRING resolved no UniProt accessions for species="
                    f"{self.ncbi_organism}; retrying with gene symbols"
                )
                def _symbols(value):
                    """Gene symbols from a field that may be None, NaN, str or list."""
                    if value is None or isinstance(value, float):
                        return []          # NaN arrives as a float
                    if isinstance(value, str):
                        return [value]
                    try:
                        return [v for v in value if v is not None]
                    except TypeError:
                        return []

                symbols = []
                for protein in self.proteins:
                    uniprot = protein.uniprot_id
                    uniprot = uniprot[0] if isinstance(uniprot, list) and uniprot else uniprot
                    for symbol in _symbols(protein.gene):
                        symbol = str(symbol).strip()
                        if symbol and symbol not in symbol_to_uniprot:
                            symbol_to_uniprot[symbol] = uniprot
                            symbols.append(symbol)
                if symbols:
                    string_mapping_dict = _map(symbols, "gene symbol")

            if not string_mapping_dict:
                logger.warning(
                    "STRING identifier mapping returned no results for "
                    f"{len(proteins_flattened)} proteins (species={self.ncbi_organism}), "
                    "by UniProt accession or gene symbol. Check that "
                    "ncbi_organism is the NCBI taxon ID STRING expects."
                )
                return

            # Fetch interactions in batches too
            from collections import defaultdict as _dd
            string_interactions_dict = _dd(list)
            ensp_ids = list(string_mapping_dict.values())
            for _i in range(0, len(ensp_ids), _BATCH):
                _batch = ensp_ids[_i: _i + _BATCH]
                _ints = string_get_interactions(
                    protein_list=_batch,
                    species=self.ncbi_organism,
                    cutoff_score=700,
                )
                for _k, _v in _ints.items():
                    string_interactions_dict[_k].extend(_v)
                logger.info(
                    f"STRING interactions batch {_i // _BATCH + 1}: "
                    f"{len(_batch)} queried → {len(_ints)} with partners"
                )

            # Translate STRING identifiers back to the queried identifiers
            translated_interactions_dict = translate_string_dict(
                mapping_dict=string_mapping_dict,
                interactions_dict=string_interactions_dict,
            )

            # If the query went through gene symbols, map both sides back to
            # UniProt so the graph stays keyed by one identifier.
            if symbol_to_uniprot:
                translated_interactions_dict = {
                    symbol_to_uniprot.get(symbol, symbol): [
                        symbol_to_uniprot.get(partner, partner) for partner in partners
                    ]
                    for symbol, partners in translated_interactions_dict.items()
                }
            
            # Assign interaction partners to proteins
            for protein in self.proteins:
                if isinstance(protein.uniprot_id, list):
                    for protein_id in protein.uniprot_id:
                        try:
                            protein.interaction_partners = translated_interactions_dict[protein_id]
                            break
                        except KeyError:
                            continue
                else:
                    try:
                        protein.interaction_partners = translated_interactions_dict[protein.uniprot_id]
                    except KeyError:
                        continue
            
            n_with_partners = sum(
                1 for protein in self.proteins if protein.interaction_partners
            )
            logger.info(
                f"Successfully mapped protein-protein interactions "
                f"({n_with_partners:,} proteins with partners)"
            )
            return True
        except Exception as e:
            logger.error(f"Error getting protein-protein interactions: {e}")
            logger.warning("Continuing without protein-protein interactions")
            # Returning False rather than None lets a caller distinguish "the
            # lookup failed and should be retried" from "it ran and found
            # nothing" -- the difference between a transient timeout and a real
            # absence of data.
            return False

    def get_transcription_factor_targets(
        self,
        genome_ChIP="mm10",
        distance_ChIP=5,
        cell_type_class_ChIP=None,
        cell_type_ChIP=None,
        transcription_factors=None,
        restrict_to_proteome=True,
        score_threshold: float = 100.0,
        max_targets_per_tf: int = None,
    ):
        """
        Get transcription factor targets from ChIP-Atlas.

        Parameters
        ----------
        genome_ChIP : str
            Genome assembly in ChIP-Atlas
        distance_ChIP : int
            Distance from transcription start site in kb
        cell_type_class_ChIP : str
            Cell type class in ChIP-Atlas
        cell_type_ChIP : str
            Specific cell type in ChIP-Atlas
        transcription_factors : list of str, optional
            Fetch only these factors, by gene symbol.
        score_threshold : float
            Minimum mean ChIP-Atlas binding score, averaged over the selected
            experiments, for a gene to count as a target. The score is a MACS2
            peak significance derived from -log10(Q), so larger means stronger
            and more significant binding. Empirically the median listed gene
            scores 10-60 and the default of 100 keeps roughly the top 15-20%
            per factor. Pass 0 for the old behaviour of accepting every gene
            ChIP-Atlas lists (~15,000 per factor, which swamps the rest of
            the network).
        max_targets_per_tf : int, optional
            Additionally cap each factor at its N highest-scoring targets.
            Off by default; useful for architectural factors such as CTCF
            that bind most promoters genuinely and strongly.
        restrict_to_proteome : bool
            When True (default) and no explicit list is given, fetch only the
            factors whose gene symbol appears in this proteome. ChIP-Atlas
            holds ~700 factors per genome and downloading all of them takes
            hours, so fetching only measured proteins is the difference between
            a usable build and an overnight one. Set False for the full set.
        """
        if not self.proteins:
            logger.error("Proteins must be set before getting transcription factor targets")
            return

        if transcription_factors is None and restrict_to_proteome:
            def _symbols(protein):
                """Gene symbols, flattened. `gene` is a list; str() on it would
                yield "['Actb']", which matches nothing."""
                value = protein.gene
                if value is None or isinstance(value, float):
                    return []
                if isinstance(value, str):
                    return [value]
                try:
                    return [str(v).strip() for v in value if v is not None]
                except TypeError:
                    return []

            transcription_factors = sorted({
                symbol for protein in self.proteins
                for symbol in _symbols(protein) if symbol
            })
            logger.info(
                f"Restricting ChIP-Atlas to the {len(transcription_factors)} "
                f"gene symbols present in this proteome; pass "
                f"restrict_to_proteome=False for all factors"
            )

        logger.info(f"Getting transcription factor targets from ChIP-Atlas for {genome_ChIP}")
        total_targets = 0
        n_factors = 0

        try:
            # Get binding data from ChIP-Atlas.
            # get_ChIP_data returns (scores, gene_to_experiment, failed); the
            # two-value form is still accepted so older callers and test
            # doubles keep working.
            result = get_ChIP_data(
                genome=genome_ChIP, distance=distance_ChIP,
                proteins=transcription_factors,
            )
            if len(result) == 3:
                scores_unsort, gene_to_experiment, failed_factors = result
            else:
                scores_unsort, gene_to_experiment = result
                failed_factors = []
            # Keep the casualties on the layer: once the network is built, a
            # factor lost to a timeout looks exactly like one with no targets.
            self.chip_atlas_failed_factors = list(failed_factors)
            
            import numpy as _np
            # The experiment list only matters when the caller asked for a
            # cell-type restriction; without one it selects every experiment
            # for the genome, which is what we already have. Fetching it
            # anyway means downloading and parsing the full cross-genome
            # experimentList.tab on every build for no effect.
            scores_columns = set(scores_unsort.columns)
            if cell_type_class_ChIP is None and cell_type_ChIP is None:
                exp_info_set = scores_columns
            else:
                exp_info = get_ChIP_exps(
                    genome=genome_ChIP,
                    cell_type_class=cell_type_class_ChIP,
                    cell_type=cell_type_ChIP,
                )
                # get_ChIP_exps returns named columns; indexing it with the
                # integer 0 raised KeyError and lost every TF target silently.
                experiment_column = (
                    "exp_id" if "exp_id" in exp_info.columns
                    else exp_info.columns[0]
                )
                exp_info_set = set(exp_info[experiment_column].to_list())
                if not exp_info_set:
                    logger.warning(
                        f"No ChIP-Atlas experiments match cell_type_class="
                        f"{cell_type_class_ChIP!r}, cell_type={cell_type_ChIP!r} "
                        f"for {genome_ChIP}; no transcription factor targets "
                        f"will be assigned"
                    )

            # Assign targets to transcription factors
            for protein in self.proteins:
                # protein.gene may be str, list, or numpy array — normalise to list
                raw_gene = protein.gene
                if raw_gene is None:
                    gene_names = []
                elif isinstance(raw_gene, _np.ndarray):
                    gene_names = [str(g) for g in raw_gene.flat if g is not None]
                elif isinstance(raw_gene, list):
                    gene_names = [str(g) for g in raw_gene if g is not None]
                else:
                    gene_names = [str(raw_gene)]

                # Score each candidate target by its MEAN binding score across
                # the selected experiments, and keep only the ones that clear
                # score_threshold.
                #
                # The previous rule took the union of every gene with a
                # non-zero score in any experiment. ChIP-Atlas only lists
                # genes that show some binding, so that rule kept essentially
                # the whole file -- ~15,000 targets per factor, which made
                # transcriptional regulation 91% of the mouse network. The
                # mean is ChIP-Atlas's own consensus score (it is exactly the
                # "{TF}|Average" column when no cell-type filter is applied),
                # so a single strong peak in 1 of 193 experiments is
                # correctly downweighted instead of promoting a target.
                experiments = [
                    experiment
                    for gene in gene_names
                    for experiment in gene_to_experiment.get(gene, [])
                    if experiment in exp_info_set
                    and experiment in scores_columns
                ]
                # dict.fromkeys: a gene symbol can appear twice on one protein
                experiments = list(dict.fromkeys(experiments))
                if not experiments:
                    continue

                try:
                    block = scores_unsort.loc[:, experiments].apply(
                        pd.to_numeric, errors='coerce'
                    )
                    mean_score = block.mean(axis=1).dropna()
                    hits = mean_score[mean_score >= score_threshold]
                    if max_targets_per_tf is not None:
                        hits = hits.nlargest(max_targets_per_tf)
                except Exception as e:
                    logger.error(f"Error processing ChIP-Atlas data: {e}")
                    continue

                if len(hits):
                    protein.transcription_factor_targets = (
                        hits.sort_values(ascending=False).index.tolist()
                    )
                    protein.transcription_factor_target_scores = (
                        hits.round(2).to_dict()
                    )
                    total_targets += len(hits)
                    n_factors += 1
            
            if failed_factors:
                logger.warning(
                    f"{len(failed_factors)} transcription factor(s) are "
                    f"missing from this network because ChIP-Atlas could not "
                    f"be reached for them: "
                    f"{failed_factors[:10]}"
                    f"{'...' if len(failed_factors) > 10 else ''}. "
                    f"Rerun to retry them (see "
                    f"self.chip_atlas_failed_factors)."
                )

            logger.info(
                f"Successfully mapped transcription factor targets: "
                f"{n_factors} factors, {total_targets:,} target links "
                f"(mean score >= {score_threshold}"
                + (f", top {max_targets_per_tf} per factor"
                   if max_targets_per_tf else "") + ")"
            )
            return True
        except Exception as e:
            logger.error(f"Error getting transcription factor targets: {e}")
            logger.warning("Continuing without transcription factor targets")
            return False

    def get_metabolites(self):
        """
        Get metabolites associated with enzymes from KEGG.
        """
        if not self.proteins:
            logger.error("Proteins must be set before getting metabolites")
            return
            
        logger.info("Getting enzyme-associated metabolites")
        
        try:
            # Get enzyme to compound mapping from KEGG
            kegg_ec_to_cpds_data = kegg_ec_to_cpds()
            kegg_compounds_list = kegg_list_compounds()
            
            # Connect metabolites to proteins based on EC numbers
            for protein in self.proteins:
                if protein.ec_number and len(protein.ec_number) > 0:
                    metabolites = []
                    for ec in protein.ec_number:
                        try:
                            # Get compounds associated with the EC number
                            if ec in kegg_ec_to_cpds_data["ec_number"].values:
                                compounds = kegg_ec_to_cpds_data[
                                    kegg_ec_to_cpds_data["ec_number"] == ec
                                ]["kegg_compounds"].tolist()
                                
                                # Add compound IDs directly
                                for compound_id in compounds:
                                    metabolites.append(compound_id)
                        except Exception as e:
                            logger.warning(f"Error getting KEGG metabolites for {ec}: {e}")
                    
                    protein.metabolites = list(set(metabolites))
            
            logger.info(f"Successfully mapped enzyme-associated metabolites")
        except Exception as e:
            logger.error(f"Error getting enzyme-associated metabolites: {e}")
            logger.warning("Continuing without enzyme-associated metabolites")

    def get_brenda_kinetics(
        self,
        organism: str = None,
        email: str = None,
        password: str = None,
        fields: list = None,
    ):
        """Enrich proteins with kinetic data from BRENDA.

        Queries BRENDA for inhibitors, activators, substrates, products, and
        cofactors for each unique EC number in the proteome and populates the
        corresponding attributes on each :class:`~transnet.biology.elements.Protein`.

        Parameters
        ----------
        organism : str, optional
            Restrict BRENDA queries to this organism (e.g. ``"Homo sapiens"``).
            Leave as ``None`` to retrieve data from all organisms.
        email, password : str, optional
            BRENDA credentials.  Fall back to ``BRENDA_EMAIL`` / ``BRENDA_PASSWORD``
            environment variables, then to an interactive prompt.
        fields : list of str, optional
            Subset of ``['inhibitors', 'activators', 'substrates', 'products',
            'cofactors']``.  Fetches all by default.
        """
        if not self.proteins:
            logger.error("Proteins must be populated before querying BRENDA")
            return

        try:
            from transnet.api.brenda import brenda_enrich_proteins
        except ImportError as exc:
            logger.error(f"BRENDA import failed: {exc}")
            return

        logger.info(
            f"Fetching BRENDA kinetics for {len(self.proteins)} proteins "
            f"(organism filter: {organism or 'all'})"
        )
        try:
            brenda_enrich_proteins(
                proteins=self.proteins,
                organism=organism,
                email=email,
                password=password,
                fields=fields,
            )
            logger.info("BRENDA kinetics enrichment complete")
        except Exception as exc:
            logger.error(f"BRENDA kinetics enrichment failed: {exc}")
            logger.warning("Continuing without BRENDA kinetics data")


class Signaling(OmicsLayer):
    """Signaling layer: kinases, phosphatases and phosphosites.

    The top of the trans-omic hierarchy.  This layer is **optional** -- most
    studies do not measure the phosphoproteome, and a network without it is a
    complete trans-omic network that simply starts at the transcription factor
    or transcript level.  When present, it supplies the ``phosphorylation`` and
    ``kinase_tf`` edges that let a regulatory path be traced from a stimulus
    down through gene expression to a metabolic reaction.

    Because kinase-substrate resources vary widely in coverage and licensing,
    the layer is designed to be filled either from KEGG signaling pathways or
    from a user-supplied table (PhosphoSitePlus, an in-house phosphoproteomics
    result, or any two-column kinase/substrate list).
    """

    def __init__(self):
        super().__init__(name="Signaling")
        self.nodes = []
        self.input_id_type = None
        self.input_ids = []
        self.input_data = None

    def __repr__(self):
        return f"<Signaling, {len(self.nodes)} signaling nodes in total>"

    def add_experimental_data(
        self,
        input_data: pd.DataFrame = None,
        id_column: str = None,
        id_type: str = "uniprot",
        site_column: str = None,
    ):
        """Register phosphoproteomics measurements for this layer.

        Parameters
        ----------
        input_data : pandas.DataFrame
            Measured phosphosites, one row per site.
        id_column : str
            Column holding the protein identifier.
        id_type : str
            Identifier type; ``"uniprot"`` (default) or ``"symbol"``.
        site_column : str, optional
            Column holding the modified residue (e.g. ``"S473"``).  When given,
            node ids become ``"<protein>_<site>"`` so that distinct sites on the
            same protein stay distinct nodes.
        """
        self.input_data = input_data
        self.input_ids = input_data[id_column].tolist() if id_column else []
        self.input_id_type = id_type
        self._site_column = site_column
        self._id_column = id_column

    def populate_from_table(
        self,
        edges: pd.DataFrame,
        source_column: str = "kinase",
        target_column: str = "substrate",
        sign_column: str = None,
        target_type_column: str = None,
        evidence: str = "user",
    ):
        """Build the layer from a kinase-substrate table.

        Parameters
        ----------
        edges : pandas.DataFrame
            One row per kinase-target relationship.
        source_column, target_column : str
            Columns naming the kinase and its target.
        sign_column : str, optional
            Column holding ``+1`` / ``-1`` for activating / inhibiting effects.
            Absent means the sign is unknown (0).
        target_type_column : str, optional
            Column whose value ``"tf"`` marks a target as a transcription
            factor, producing a ``kinase_tf`` edge instead of a
            ``phosphorylation`` edge.
        evidence : str
            Provenance string recorded on every edge built from this table.

        Returns
        -------
        Signaling
            ``self``, so the call can be chained.
        """
        by_id = {node.id: node for node in self.nodes}

        for _, row in edges.iterrows():
            source = str(row[source_column])
            target = str(row[target_column])
            if not source or not target or source == "nan" or target == "nan":
                continue

            node = by_id.get(source)
            if node is None:
                sign = 0
                if sign_column and sign_column in row.index:
                    try:
                        sign = int(row[sign_column])
                    except (TypeError, ValueError):
                        sign = 0
                node = SignalingProtein(
                    id=source, name=source, uniprot_id=source,
                    sign=sign, evidence=evidence,
                )
                by_id[source] = node
                self.nodes.append(node)

            is_tf = False
            if target_type_column and target_type_column in row.index:
                is_tf = str(row[target_type_column]).lower() == "tf"

            bucket = node.tf_targets if is_tf else node.substrates
            if target not in bucket:
                bucket.append(target)

        logger.info(
            f"Populated signaling layer with {len(self.nodes)} nodes "
            f"from {len(edges)} kinase-target rows"
        )
        return self

    def populate(self, kegg_organism: str = None, pathway_ids: List[str] = None):
        """Populate kinase-substrate relationships from KEGG signaling pathways.

        Parameters
        ----------
        kegg_organism : str
            KEGG organism code, e.g. ``"hsa"``.
        pathway_ids : list of str, optional
            Restrict to these KEGG pathway ids.  Defaults to the signaling
            pathways KEGG classifies under "Signal transduction".

        Notes
        -----
        KEGG's coverage of phosphorylation relationships is partial.  For a
        study with its own phosphoproteomics, :meth:`populate_from_table` gives
        a more complete layer.
        """
        if not kegg_organism:
            logger.error("kegg_organism is required to populate the signaling layer")
            return self

        try:
            from transnet.api.kegg import kegg_signaling_relations
        except ImportError as exc:
            logger.error(f"KEGG signaling import failed: {exc}")
            return self

        try:
            relations = kegg_signaling_relations(
                kegg_organism, pathway_ids=pathway_ids
            )
        except Exception as exc:
            logger.error(f"Failed to fetch KEGG signaling relations: {exc}")
            logger.warning("Signaling layer left empty; the network is still valid without it")
            return self

        if relations is None or relations.empty:
            logger.warning("KEGG returned no signaling relations for this organism")
            return self

        return self.populate_from_table(
            relations,
            source_column="source",
            target_column="target",
            sign_column="sign",
            target_type_column="target_type",
            evidence="KEGG",
        )
