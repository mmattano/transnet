"""Functions to access UniProt database."""

__all__ = [
    "uniprot_list_proteins",
]

import pandas as pd


def uniprot_list_proteins(
    org_ncbi_num: str = 9606,
) -> pd.DataFrame:
    """
    Reference table for KEGG genes.

    Parameters
    ----------
    org_ncbi_num : str
        NCBI organism number.

    Returns
    -------
    organism_proteins : pd.DataFrame
        Dataframe with UniProt proteins and additional information.

    """
    url = (
        f"https://www.uniprot.org/uniprot/?query=reviewed:"
        f"yes+AND+organism:{org_ncbi_num}&format=tab&"
        f"columns=id,entry%20name,reviewed,protein%20names,"
        f"genes,organism,length,ec"
    )

    organism_proteins = pd.read_table(url, sep="\t", header=0)
    cols = organism_proteins.columns
    organism_proteins[cols[4]] = (
        organism_proteins[cols[4]].str.split(" ")
    )
    return organism_proteins
