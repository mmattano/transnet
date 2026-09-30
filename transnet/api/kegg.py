"""Functions to access KEGG database."""

__all__ = [
    "kegg_conv_ncbi_idtable",
    "kegg_link_pathway",
    "kegg_link_ec",
    "kegg_pathwaymap_entries",
    "kegg_list_pathways",
    "kegg_list_genes",
    "kegg_ftp_list_genes",
    "kegg_list_organisms",
    "kegg_get_organism_info",
    "kegg_ec_to_cpds",
    "kegg_list_compounds",
    "kegg_create_reaction_table",
    "kegg_to_chebi",
    "chebi_to_kegg",
    "kegg_signaling_relations",
    "kegg_reaction_pathways",
]

import pandas as pd
from urllib.error import HTTPError
import time
import requests
import logging

logger = logging.getLogger(__name__)



def kegg_conv_ncbi_idtable(org_abb: str = "hsa",) -> pd.DataFrame:
    """
    Conversion table for NCBI IDs to KEGG IDs.

    Parameters
    ----------
    org_abb : str
        KEGG organism abbreviation.

    Returns
    -------
    organism_id_conversion : pd.DataFrame
        Dataframe with NCBI IDs and KEGG IDs.

    """

    url = f"http://rest.kegg.jp/conv/{org_abb}/ncbi-geneid"

    organism_id_conversion = pd.read_table(url, sep="\t", header=None)
    cols = organism_id_conversion.columns
    for col in cols:
        organism_id_conversion[col] = (
            organism_id_conversion[col].str.split(":", expand=True).iloc[:, 1]
        )
    organism_id_conversion.columns = ["ncbi_id", "kegg_id"]
    return organism_id_conversion


def kegg_link_pathway(
    org_abb: str = "hsa", pathway_id: str = None,
) -> pd.DataFrame:
    """
    Link table for genes in pathways.

    Parameters
    ----------
    org_abb : str
        KEGG organism abbreviation.
    pathway_id : str
        KEGG pathway ID.

    Returns
    -------
    organism_pathway_genes : pd.DataFrame
        Dataframe with genes and the pathways they are present in.

    """
    if pathway_id:
        url = f"http://rest.kegg.jp/link/{org_abb}/{pathway_id}"
    else:
        url = f"http://rest.kegg.jp/link/{org_abb}/pathway"

    organism_pathway_genes = pd.read_table(url, sep="\t", header=None)
    cols = organism_pathway_genes.columns
    for col in cols:
        organism_pathway_genes[col] = (
            organism_pathway_genes[col].str.split(":", expand=True).iloc[:, 1]
        )
    organism_pathway_genes.columns = ["pathway", "kegg_gene_id"]
    return organism_pathway_genes


def kegg_link_ec(org_abb: str = "hsa",) -> pd.DataFrame:
    """
    Link table for genes and their corresponding enzyme comission numbers.

    Parameters
    ----------
    org_abb : str
        KEGG organism abbreviation.

    Returns
    -------
    organism_ec : pd.DataFrame
        Dataframe with genes and their corresponding enzyme comission numbers.

    """
    # KEGG rejects "ec" as the target database when the source is an organism
    # (`/link/ec/sce` returns 400), though it still accepts it for pathway maps.
    # "enzyme" is the spelling that works, and returns the same two columns in
    # the same order: organism gene, then EC number.
    url = f"https://rest.kegg.jp/link/enzyme/{org_abb}"

    organism_ec = pd.read_table(url, sep="\t", header=None)
    cols = organism_ec.columns
    for col in cols:
        organism_ec[col] = (
            organism_ec[col].str.split(":", expand=True).iloc[:, 1]
        )
    organism_ec.columns = ["kegg_gene_id", "ec_number"]
    return organism_ec


def kegg_pathwaymap_entries(
    pathway_map: str = "map00010", entry_type: str = "",
) -> pd.DataFrame:
    """
    Shows info on elements in pathway map.

    Parameters
    ----------
    pathway_map : str
        KEGG map ID.
    entry_type : str
        Type of information requested: 'rn' = Reaction, 'ec' = EC number,
        'cpd' = compound.

    Returns
    -------
    pathwaymap_entries : pd.DataFrame
        Dataframe with genes and the pathways they are present in.

    Raises
    ------
    ValueError
        If entry_type is not 'rn', 'ec', or 'cpd'.
    """

    if entry_type in ["rn", "ec", "cpd"]:
        url = f"http://rest.kegg.jp/link/{entry_type}/{pathway_map}"

        pathwaymap_entries = pd.read_table(url, sep="\t", header=None)
        cols = pathwaymap_entries.columns
        for col in cols:
            pathwaymap_entries[col] = (
                pathwaymap_entries[col].str.split(":", expand=True).iloc[:, 1]
            )
        pathwaymap_entries.columns = ["pathway", f"{entry_type}_entries"]
        return pathwaymap_entries
    else:
        raise ValueError(f"'{entry_type}' is not a valid entry type.")


def kegg_list_pathways(org_abb: str = "hsa",) -> pd.DataFrame:
    """
    Reference table for KEGG pathways.

    Parameters
    ----------
    org_abb : str
        KEGG organism abbreviation.

    Returns
    -------
    organism_pathways : pd.DataFrame
        Dataframe with KEGG pathways and descriptions.

    """

    url = f"http://rest.kegg.jp/list/pathway/{org_abb}"

    organism_pathways = pd.read_table(url, sep="\t", header=None)
    cols = organism_pathways.columns
    # organism_pathways[cols[0]] = (
    #    #organism_pathways[cols[0]].str.split(":", expand=True).iloc[:, 1]
    #    organism_pathways[cols[0]].str.split(org_abb, expand=True).iloc[:, 1]
    # )
    organism_pathways[cols[1]] = (
        organism_pathways[cols[1]].str.split(" - ", expand=True).iloc[:, 0]
    )
    organism_pathways.columns = ["pathways_id", "description"]
    return organism_pathways


def kegg_list_genes(org_abb: str = "hsa",) -> pd.DataFrame:
    """
    Reference table for KEGG genes.

    Parameters
    ----------
    org_abb : str
        KEGG organism abbreviation.

    Returns
    -------
    organism_genes : pd.DataFrame
        Dataframe with KEGG genes, alternative names and descriptions.
    """

    url = f"http://rest.kegg.jp/list/{org_abb}"

    organism_genes = pd.read_table(url, sep="\t", header=None)
    cols = organism_genes.columns
    organism_genes.iloc[:, cols[0]] = (
        organism_genes[cols[0]].str.split(":", expand=True).iloc[:, 1]
    )
    organism_genes_info = organism_genes[cols[3]].str.split("; ", expand=True)
    organism_genes.iloc[:, cols[2]] = organism_genes_info.iloc[:, 0]
    organism_genes.iloc[:, cols[3]] = organism_genes_info.iloc[:, 1]
    organism_genes.columns = ["gene_id", "type", "name(s)", "description"]
    return organism_genes


def kegg_ftp_list_genes(
    org_abb: str = "hsa",
    ftp_base_path: str = "",
    path_to_genome_file: str = None,
) -> pd.DataFrame:
    """
    Reference table for KEGG genes.

    Parameters
    ----------
    org_abb : str
        KEGG organism abbreviation.
    ftp_base_path : str
        Base path to the folder containing the KEGG FTP files.
    path_to_genome_file : str, optional
        The path to the KEGG genome file. If None, the name will be pulled by the API, but the file will have to be in the base folder.

    Returns
    -------
    organism_genes_ftp : pd.DataFrame
        Dataframe with KEGG genes, alternative names and descriptions.
    """

    if path_to_genome_file is None:
        url = f"https://rest.kegg.jp/get/genome:{org_abb}"
        genome_info = pd.read_table(url, sep="\t", header=None)
        genome_info_sep = [" ".join(info.split()) for info in genome_info[0]]
        path_to_genome_file = [
            info.split(" ")[1] for info in genome_info_sep if "ENTRY " in info
        ][0]

    path = ftp_base_path + "/" + path_to_genome_file + ".kff"

    organism_genes = pd.read_table(path, sep="\t", header=None)
    cols = organism_genes.columns
    organism_genes_ftp = organism_genes.loc[
        :, [cols[5], cols[1], cols[7], cols[8]]
    ]
    organism_genes_ftp.columns = ["gene_id", "type", "name(s)", "description"]
    return organism_genes_ftp


def kegg_list_organisms() -> pd.DataFrame:
    """
    Reference table for KEGG organisms.

    Returns
    -------
    organisms : pd.DataFrame
        Dataframe with the organisms in the KEGG database
        and associated information.
    """

    url = "http://rest.kegg.jp/list/organism"

    organisms = pd.read_table(url, sep="\t", header=None)
    organisms.columns = ["kegg_id", "kegg_name", "name", "taxonomy"]
    return organisms


def kegg_get_organism_info(org_abb: str = "hsa",) -> pd.DataFrame:
    """
    Get full organism name and NCBI ID.

    Returns
    -------
    organisms : pd.DataFrame
        Dataframe with the organisms in the KEGG database
        and associated information.

    """
    url = f"https://rest.kegg.jp/get/genome:{org_abb}"

    genome_info = pd.read_table(
        url, sep="\t", header=None, skiprows=1, skipfooter=1
    )
    genome_info_sep = [" ".join(info.split()) for info in genome_info[0]]
    name = [
        info.split("NAME ")[1] for info in genome_info_sep if "NAME " in info
    ][0]
    taxonomy = [
        info.split("TAXONOMY TAX:")[1]
        for info in genome_info_sep
        if "TAXONOMY " in info
    ][0]
    return [name, taxonomy]


def kegg_ec_to_cpds() -> pd.DataFrame:
    """
    Produces lookup table of EC numbers and interacting metabolites.

    Returns
    -------
    ec_compounds : pd.DataFrame
        Dataframe with EC numbers of proteins and the metabolites that interact with them.

    """
    url = "https://rest.kegg.jp/link/compound/ec"

    ec_compounds = pd.read_table(url, sep="\t", header=None)
    cols = ec_compounds.columns
    ec_compounds[cols[0]] = (
        ec_compounds[cols[0]].str.split(":", expand=True).iloc[:, 1]
    )
    ec_compounds[cols[1]] = (
        ec_compounds[cols[1]].str.split(":", expand=True).iloc[:, 1]
    )
    ec_compounds.columns = ["ec_number", "kegg_compounds"]
    return ec_compounds


def kegg_list_compounds() -> pd.DataFrame:
    """
    Reference table for KEGG compounds.

    Returns
    -------
    compounds : pd.DataFrame
        Dataframe with KEGG compounds and their (alternative) names.

    """

    url = "https://rest.kegg.jp/list/compound"

    compounds = pd.read_table(url, sep="\t", header=None)
    cols = compounds.columns
    compounds.columns = ["kegg_compounds", "name"]
    return compounds


def _kegg_get_equation(reaction):
    """
    Get the equation for a KEGG reaction.

    Parameters
    ----------
    reaction : str
        KEGG reaction ID.
    
    Returns
    -------
    equation : str
        KEGG reaction equation.
    defintion : str
        KEGG reaction definition.
    enzyme : str
        EC number of the enzyme executing the reaction.
    repeat : bool
        True if the function has to repeat the request to the server.
    """

    reaction_url = "https://rest.kegg.jp/get/" + reaction
    repeat = False
    try:
        r = requests.get(reaction_url)

        info_string = r.content.decode('utf-8')

        info_lines = info_string.split('\n')

        info_dict = {}

        for line in info_lines:
            if line.startswith('NAME'):
                info_dict['NAME'] = line.split('NAME        ')[1]
            elif line.startswith('DEFINITION'):
                info_dict['DEFINITION'] = line.split('DEFINITION  ')[1]
            elif line.startswith('EQUATION'):
                info_dict['EQUATION'] = line.split('EQUATION    ')[1]
            elif line.startswith('ENZYME'):
                info_dict['ENZYME'] = line.split('ENZYME      ')[1]

        try:
            equation = info_dict['EQUATION']
            defintion = info_dict['DEFINITION']
        except KeyError:
            equation = None
            defintion = None
            repeat = True
        try:
            enzyme = info_dict['ENZYME']
        except KeyError:
            enzyme = None
    except requests.exceptions.RequestException:
        equation = None
        defintion = None
        enzyme = None
        repeat = True
    return equation, defintion, enzyme, repeat


def kegg_create_reaction_table(print_max_repeats_needed=False):
    """
    Create a table with all KEGG reactions and their associated information.

    Returns
    -------
    reactions : pd.DataFrame
        Dataframe with all KEGG reactions and their associated information.
    """

    url = "https://rest.kegg.jp/list/reaction"
    r = requests.get(url)
    kegg_reaction_info = [x.split("\t") for x in r.content.decode('utf-8').split("\n") if x != '']
    reactions = pd.DataFrame(kegg_reaction_info, columns=['reaction', 'name'])
    equations = []
    definitions = []
    enzymes = []
    substrates = []
    stoichiometry_substrates = []
    products = []
    stoichiometry_products = []
    reversibilities = []
    max_repeats_needed = 0

    # KEGG writes equations with one of three arrows.  Only "<=>" asserts
    # reversibility; the directed arrows mark an irreversible reaction, which
    # matters for tracing regulatory paths through the metabolic layer.
    arrows = [(" <=> ", True), (" => ", False), (" <= ", False)]

    for i, reaction in enumerate(reactions["reaction"].to_list()):
        equation, defintion, enzyme, repeat = _kegg_get_equation(reaction)
        if repeat:
            counter = 0
            while repeat:
                time.sleep(0.01)
                equation, defintion, enzyme, repeat = _kegg_get_equation(
                    reaction
                )
                counter += 1
                if counter > max_repeats_needed:
                    max_repeats_needed = counter
                time.sleep(0.1)
        if enzyme is not None:
            enzyme = enzyme.split('        ')

        reversible = True
        eq_parts = [equation]
        for arrow, arrow_reversible in arrows:
            if arrow in equation:
                eq_parts = equation.split(arrow)
                reversible = arrow_reversible
                # "A <= B" reads right-to-left; normalise so that eq_parts[0]
                # is always the substrate side.
                if arrow == " <= ":
                    eq_parts = eq_parts[::-1]
                break
        if len(eq_parts) < 2:
            logger.warning(f"Could not parse equation for {reaction}: {equation!r}")
            eq_parts = [eq_parts[0], ""]
        reversibilities.append(reversible)

        temp_substrates = eq_parts[0].split(" + ")
        stoichiometry_substrate = []
        for i, substrate in enumerate(temp_substrates):
            if len(substrate.split(" ")) > 1:
                stoic = substrate.split(" ")[0]
                # Handle non-numeric stoichiometry (e.g., 'n' for variable)
                try:
                    stoichiometry_substrate.append(float(stoic))
                except ValueError:
                    stoichiometry_substrate.append(None)  # Variable stoichiometry
                temp_substrates[i] = substrate.split(" ")[1]
            else:
                stoichiometry_substrate.append(1.0)
        temp_products = eq_parts[1].split(" + ")
        stoichiometry_product = []
        for i, product in enumerate(temp_products):
            if len(product.split(" ")) > 1:
                stoic = product.split(" ")[0]
                # Handle non-numeric stoichiometry (e.g., 'n' for variable)
                try:
                    stoichiometry_product.append(float(stoic))
                except ValueError:
                    stoichiometry_product.append(None)  # Variable stoichiometry
                temp_products[i] = product.split(" ")[1]
            else:
                stoichiometry_product.append(1.0)
        substrates.append(temp_substrates)
        products.append(temp_products)
        stoichiometry_substrates.append(stoichiometry_substrate)
        stoichiometry_products.append(stoichiometry_product)
        equations.append(equation)
        definitions.append(defintion)
        enzymes.append(enzyme)

    reactions["equation"] = equations
    reactions["definition"] = definitions
    reactions["enzyme"] = enzymes
    reactions["substrates"] = substrates
    reactions["products"] = products
    reactions["stoichiometry_substrates"] = stoichiometry_substrates
    reactions["stoichiometry_products"] = stoichiometry_products
    reactions["reversible"] = reversibilities

    if print_max_repeats_needed:
        print("Maximum number of repeats needed:", max_repeats_needed)

    return reactions


def kegg_to_chebi(kegg_compound_ids):
    """Generates conversion table from KEGG to ChEBI.
    
    Parameters
    ----------
    kegg_compound_ids : list
        A list of KEGG compound IDs.

    Returns
    -------
    kegg_chebi_conversion = pandas.DataFrame
        A dataframe with KEGG compound IDs and ChEBI IDs.
    """
    chebi_ids = []

    from bioservices import KEGG as _KEGG
    kegg_bio = _KEGG(verbose=False)
    map_kegg_chebi = kegg_bio.conv("chebi", "compound")

    for compound in kegg_compound_ids:
        cpd_id = f"cpd:{compound}"
        if cpd_id in map_kegg_chebi:
            chebi_id = map_kegg_chebi[cpd_id].upper()
            chebi_ids.append(chebi_id)
        else:
            chebi_ids.append(None)
            continue

    kegg_chebi_conversion = pd.DataFrame(
        {"kegg_compounds": kegg_compound_ids, "chebi_compounds": chebi_ids,},
        columns=["kegg_compounds", "chebi_compounds"],
    )

    return kegg_chebi_conversion

def chebi_to_kegg(chebi_compound_ids):
    """Convert ChEBI IDs to KEGG compound IDs using the KEGG REST API.

    Parameters
    ----------
    chebi_compound_ids : list
        A list of ChEBI compound IDs (e.g. ``['CHEBI:15422', 'CHEBI:16761']``).

    Returns
    -------
    chebi_kegg_conversion : pandas.DataFrame
        DataFrame with columns ``chebi_compounds`` and ``kegg_compounds``.
    """
    def _bare(value):
        """The numeric part of a ChEBI id.

        KEGG's conversion table writes ``chebi:10`` while ChEBI itself, and
        every client that talks to it, writes ``CHEBI:10``. Comparing the two
        verbatim never matches, so both sides are reduced to the number.
        """
        text = str(value).strip().lower()
        if text.startswith("chebi:"):
            text = text[len("chebi:"):]
        return text.strip()

    # Download the full ChEBI->KEGG mapping from KEGG REST in one request
    r = requests.get("https://rest.kegg.jp/conv/compound/chebi", timeout=120)
    r.raise_for_status()
    chebi_kegg_map = {}
    for line in r.text.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) == 2:
            chebi_kegg_map[_bare(parts[0])] = parts[1].replace("cpd:", "").strip()

    kegg_ids = [
        chebi_kegg_map.get(_bare(c)) if c is not None else None
        for c in chebi_compound_ids
    ]
    matched = sum(1 for k in kegg_ids if k)
    logger.info(
        f"chebi_to_kegg: {matched}/{len(kegg_ids)} ChEBI ids resolved to KEGG "
        f"compounds (table holds {len(chebi_kegg_map):,} mappings)"
    )
    return pd.DataFrame(
        {"chebi_compounds": chebi_compound_ids, "kegg_compounds": kegg_ids},
        columns=["chebi_compounds", "kegg_compounds"],
    )


#: KEGG relation subtypes that correspond to a phosphorylation-style
#: regulatory effect, with the sign of that effect.
_KEGG_RELATION_SIGNS = {
    "activation": 1,
    "expression": 1,
    "phosphorylation": 0,
    "inhibition": -1,
    "repression": -1,
    "dephosphorylation": 0,
}


def kegg_signaling_relations(organism, pathway_ids=None, max_pathways=None):
    """Extract kinase-substrate relations from KEGG signaling pathway KGML.

    Parses the KGML of an organism's signal-transduction pathways and returns
    the protein-protein relations KEGG annotates as activation, inhibition,
    phosphorylation or dephosphorylation.  These become the ``phosphorylation``
    and ``kinase_tf`` edges of the optional Signaling layer.

    Parameters
    ----------
    organism : str
        KEGG organism code, e.g. ``"hsa"``.
    pathway_ids : list of str, optional
        Restrict to these pathway ids (with or without the organism prefix).
        Defaults to every pathway KEGG lists for the organism whose name
        mentions signaling.
    max_pathways : int, optional
        Stop after this many pathways.  Useful for demos; KGML is fetched one
        pathway at a time.

    Returns
    -------
    pandas.DataFrame
        Columns ``source``, ``target``, ``sign``, ``subtype``, ``pathway``.
        Empty if nothing could be retrieved.

    Notes
    -----
    KEGG's phosphorylation annotation is partial.  A study with its own
    phosphoproteomics will get a fuller layer from
    :meth:`transnet.biology.layers.Signaling.populate_from_table`.
    """
    import xml.etree.ElementTree as ET

    if pathway_ids is None:
        try:
            pathways = kegg_list_pathways(organism)
        except Exception as exc:
            logger.error(f"Could not list KEGG pathways for {organism}: {exc}")
            return pd.DataFrame(
                columns=["source", "target", "sign", "subtype", "pathway"]
            )
        # kegg_list_pathways returns columns ["pathways_id", "description"].
        mask = pathways["description"].str.contains(
            "signal|signaling|MAPK|PI3K|insulin|AMPK|mTOR", case=False, na=False
        )
        pathway_ids = pathways.loc[mask, "pathways_id"].tolist()

    if max_pathways:
        pathway_ids = pathway_ids[:max_pathways]

    rows = []
    for pathway_id in pathway_ids:
        pid = str(pathway_id).replace("path:", "")
        if not pid.startswith(organism):
            # Accept a bare map number ("04010") as well as a prefixed id.
            pid = f"{organism}{pid.lstrip('map')}"
        try:
            response = requests.get(f"https://rest.kegg.jp/get/{pid}/kgml", timeout=30)
            response.raise_for_status()
            root = ET.fromstring(response.content)
        except Exception as exc:
            logger.debug(f"Skipping KGML for {pid}: {exc}")
            continue

        # entry id -> gene identifiers it represents
        entries = {}
        for entry in root.findall("entry"):
            if entry.get("type") not in ("gene", "ortholog"):
                continue
            names = [n.replace(f"{organism}:", "") for n in (entry.get("name") or "").split()]
            entries[entry.get("id")] = names

        for relation in root.findall("relation"):
            if relation.get("type") not in ("PPrel", "GErel"):
                continue
            sources = entries.get(relation.get("entry1"), [])
            targets = entries.get(relation.get("entry2"), [])
            if not sources or not targets:
                continue
            for subtype in relation.findall("subtype"):
                name = subtype.get("name")
                if name not in _KEGG_RELATION_SIGNS:
                    continue
                sign = _KEGG_RELATION_SIGNS[name]
                is_tf = relation.get("type") == "GErel"
                for source in sources:
                    for target in targets:
                        rows.append({
                            "source": source,
                            "target": target,
                            "sign": sign,
                            "subtype": name,
                            "target_type": "tf" if is_tf else "protein",
                            "pathway": pid,
                        })
        time.sleep(0.1)

    relations = pd.DataFrame(
        rows,
        columns=["source", "target", "sign", "subtype", "target_type", "pathway"],
    )
    if not relations.empty:
        relations = relations.drop_duplicates(
            subset=["source", "target", "subtype"]
        ).reset_index(drop=True)

    logger.info(
        f"Extracted {len(relations)} signaling relations from "
        f"{len(pathway_ids)} KEGG pathways"
    )
    return relations


#: KEGG's global and overview maps (map011xx, map012xx) contain most of
#: metabolism, so a reaction's membership in them says nothing about which
#: pathway it belongs to.
_OVERVIEW_MAP_PREFIXES = ("map011", "map012")


def kegg_reaction_pathways(reactions=None, organism: str = None,
                           exclude_overview: bool = True, cache: bool = True):
    """KEGG pathways each reaction belongs to.

    The result is the ``pathway_map`` that
    :func:`~transnet.regulation_axis_summary` needs for a per-pathway summary.
    A reaction usually sits in several pathways and is listed under each.

    Parameters
    ----------
    reactions : iterable of str, optional
        Reaction ids (``R00200``) to keep. Default: all KEGG reactions.
    organism : str, optional
        KEGG organism code (``"mmu"``). Keeps only pathways that exist in that
        organism, which removes, for example, antibiotic biosynthesis from a
        mouse analysis.
    exclude_overview : bool
        Leave out the global and overview maps such as "Metabolic pathways".
    cache : bool
        Store the KEGG answer in ``~/.cache/transnet/kegg`` and reuse it.

    Returns
    -------
    dict
        ``{reaction_id: [pathway name, ...]}``.
    """
    import io
    import os
    from pathlib import Path

    folder = Path(os.environ.get("TRANSNET_KEGG_CACHE",
                                 Path.home() / ".cache" / "transnet" / "kegg"))
    links_file, names_file = folder / "link_pathway_reaction.tsv", folder / "list_pathway.tsv"

    def fetch(url, path):
        if cache and path.exists():
            return path.read_text()
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        if cache:
            folder.mkdir(parents=True, exist_ok=True)
            path.write_text(response.text)
        return response.text

    links = pd.read_table(io.StringIO(fetch("https://rest.kegg.jp/link/pathway/reaction",
                                            links_file)),
                          header=None, names=["reaction", "pathway"])
    names = pd.read_table(io.StringIO(fetch("https://rest.kegg.jp/list/pathway", names_file)),
                          header=None, names=["pathway", "name"])
    links["reaction"] = links["reaction"].str.replace("rn:", "", regex=False)
    links["pathway"] = links["pathway"].str.replace("path:", "", regex=False)
    links = links[links["pathway"].str.startswith("map")]
    if exclude_overview:
        links = links[~links["pathway"].str.startswith(_OVERVIEW_MAP_PREFIXES)]
    if reactions is not None:
        links = links[links["reaction"].isin(set(reactions))]
    if organism:
        present = pd.read_table(
            io.StringIO(fetch(f"https://rest.kegg.jp/list/pathway/{organism}",
                              folder / f"list_pathway_{organism}.tsv")),
            header=None, names=["pathway", "name"])
        kept = {"map" + str(p)[len(organism):] for p in present["pathway"]}
        links = links[links["pathway"].isin(kept)]
    links = links.copy()
    links["name"] = links["pathway"].map(dict(zip(names["pathway"], names["name"])))
    links["name"] = links["name"].fillna(links["pathway"])
    return links.groupby("reaction")["name"].apply(sorted).to_dict()
