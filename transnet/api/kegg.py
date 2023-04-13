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
]

import pandas as pd
from urllib.error import HTTPError
from bioservices import KEGG, ChEBI
import time
import requests


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
    url = f"http://rest.kegg.jp/link/ec/{org_abb}"

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
    max_repeats_needed = 0

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
        eq_parts = equation.split(" <=> ")
        temp_substrates = eq_parts[0].split(" + ")
        stoichiometry_substrate = []
        for i, substrate in enumerate(temp_substrates):
            if len(substrate.split(" ")) > 1:
                stoic = substrate.split(" ")[0]
                stoichiometry_substrate.append(stoic)
                temp_substrates[i] = substrate.split(" ")[1]
            else:
                stoichiometry_substrate.append(1)
        temp_products = eq_parts[1].split(" + ")
        stoichiometry_product = []
        for i, product in enumerate(temp_products):
            if len(product.split(" ")) > 1:
                stoic = product.split(" ")[0]
                stoichiometry_product.append(stoic)
                temp_products[i] = product.split(" ")[1]
            else:
                stoichiometry_product.append(1)
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

    kegg_bio = KEGG(verbose=False)
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
    """Inverse of kegg_to_chebi.
    
    Parameters
    ----------
    chebi_compound_ids : list
        A list of ChEBI compound IDs.
        
    Returns
    -------
    chebi_kegg_conversion = pandas.DataFrame
        A dataframe with ChEBI compound IDs and KEGG IDs.
    """
    kegg_ids = []

    chebi_bio = ChEBI()
    for compound in chebi_compound_ids:
        if compound is None:
            kegg_ids.append(None)
            continue
        chebi_entry = chebi_bio.getCompleteEntity(compound)
        try:
            found_kegg = False
            for db_links in chebi_entry.DatabaseLinks:
                if db_links.type == "KEGG COMPOUND accession":
                    kegg_ids.append(db_links.data)
                    found_kegg = True
                    break
            if not found_kegg:
                kegg_ids.append(None)
        except AttributeError:
            kegg_ids.append(None)
            continue

    chebi_kegg_conversion = pd.DataFrame(
        {"chebi_compounds": chebi_compound_ids, "kegg_compounds": kegg_ids,},
        columns=["chebi_compounds", "kegg_compounds"],
    )

    return chebi_kegg_conversion
