"""Functions to access KEGG database."""

__all__ = [
    "kegg_conv_ncbi_idtable",
    "kegg_link_pathway",
    "kegg_link_ec",
    "kegg_pathwaymap_entries",
    "kegg_list_pathways",
    "kegg_list_genes",
    "kegg_list_organisms",
]

import pandas as pd


def kegg_conv_ncbi_idtable(
    org_abb: str = "hsa",
) -> pd.DataFrame:
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
    org_abb: str = "hsa",
    pathway_id: str = None,
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


def kegg_link_ec(
    org_abb: str = "hsa",
) -> pd.DataFrame:
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
    pathway_map: str = "map00010",
    entry_type: str = "",
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


def kegg_list_pathways(
    org_abb: str = "hsa",
) -> pd.DataFrame:
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
    organism_pathways[cols[0]] = (
        organism_pathways[cols[0]].str.split(":", expand=True).iloc[:, 1]
    )
    organism_pathways[cols[1]] = (
        organism_pathways[cols[1]].str.split(" - ", expand=True).iloc[:, 0]
    )
    organism_pathways.columns = ["pathways_id", "description"]
    return organism_pathways


def kegg_list_genes(
    org_abb: str = "hsa",
) -> pd.DataFrame:
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
    organism_genes[cols[0]] = (
        organism_genes[cols[0]].str.split(":", expand=True).iloc[:, 1]
    )
    organism_genes_info = organism_genes[cols[1]].str.split("; ", expand=True)
    organism_genes[cols[1]] = organism_genes_info.iloc[:, 0]
    organism_genes[2] = organism_genes_info.iloc[:, 1]
    organism_genes.columns = ["gene_id", "name(s)", "description"]
    return organism_genes


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
