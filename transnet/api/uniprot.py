"""Functions to access UniProt database."""

__all__ = ["uniprot_list_proteins", "uniprot_add_entrez_id"]

import pandas as pd
import mygene


def uniprot_list_proteins(org_ncbi_num: str = 9606,) -> pd.DataFrame:
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
        f"https://rest.uniprot.org/uniprotkb/stream?fields="
        f"accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2C"
        f"organism_name%2Clength%2Cxref_string%2Cxref_brenda%2C"
        f"xref_reactome%2Cxref_ensembl%2Cxref_kegg%2Cxref_geneid%2C"
        f"ec&format=tsv&query=(organism_id:{org_ncbi_num})"
    )

    organism_proteins = pd.read_table(url, sep="\t", header=0)
    cols = organism_proteins.columns
    organism_proteins[cols[4]] = organism_proteins[cols[4]].str.split(" ")
    return organism_proteins


def uniprot_add_entrez_id(uniprot_df):
    """
    Add Entrez ID to uniprot dataframe.

    Parameters
    ----------
    uniprot_df : pd.DataFrame
        Dataframe with UniProt proteins and additional information.
    
    Returns
    -------
    uniprot_df : pd.DataFrame
        Dataframe with UniProt proteins and additional information.
    """
    uniprot_df = uniprot_df.copy(deep=True)
    uniprot_entries = uniprot_df["Entry"].to_list()
    mg = mygene.MyGeneInfo()
    # Query uniprot entries names using mygene
    pinfo = mg.querymany(
        uniprot_entries,
        scopes="uniprot",
        returnall=True,
        fields=["query", "name", "symbol", "entrezgene", "uniprot"],
    )
    entrez_id_list = [[] for i in range(len(uniprot_entries))]
    for entry in pinfo["out"]:
        try:
            entrez_id_list[uniprot_entries.index(entry["query"])].append(
                entry["entrezgene"]
            )
        except KeyError:
            pass
    uniprot_df["Entrez"] = entrez_id_list
    return uniprot_df
