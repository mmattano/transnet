"""Functions to access Ensembl database."""

__all__ = [
    "ensembl_download_release",
    "ensembl_get_transcripts",
    "ensembl_expander",
]

import subprocess
from pyensembl import EnsemblRelease
import pandas as pd
import mygene


def ensembl_download_release(release_number, organism_full):
    """Download Ensembl release.

    Parameters
    ----------
    release_number : int
        Ensembl release number.
    organism_full : str
        Organism name.
    """

    cmd_str = f"pyensembl install --release {release_number} --species {organism_full}"
    subprocess.run(cmd_str, shell=True)


def ensembl_get_transcripts(
    release_number, organism_full, drop_duplicates=True
):
    """Get transcripts from Ensembl.
    
    Parameters
    ----------
    organism_full : str
        Organism name.
    release_number : int
        Ensembl release number.
    
    Returns
    -------
    transcript_df : pd.DataFrame
        List of transcripts.
    """

    data = EnsemblRelease(release_number, species=f"{organism_full}")
    transcript_df = pd.DataFrame(
        [
            [transcript.gene_id, transcript.transcript_id, transcript.biotype]
            for transcript in data.transcripts()
        ],
        columns=["ensembl_gene_id", "ensembl_transcript_id", "biotype"],
    )

    return transcript_df


def ensembl_expander(transcript_df):
    """Expand transcript dataframe with gene names and descriptions.

    Parameters
    ----------
    transcript_df : pd.DataFrame
        Data frame of transcripts.
    
    Returns
    -------
    transcript_df : pd.DataFrame
        Data frame of transcripts with gene names and descriptions.
    """
    transcript_df = transcript_df.copy(deep=True)

    # Get gene names
    gene_ids_list = transcript_df["ensembl_transcript_id"].to_list()
    mg = mygene.MyGeneInfo()
    # Query gene names using mygene
    ginfo = mg.querymany(
        gene_ids_list,
        scopes="ensembltranscript",
        returnall=True,
        fields=["query", "name", "symbol", "entrezgene", "uniprot"],
    )
    # Set up lists to write into
    # I'm doing it like this to potentially account for missing genes
    # Like this I also don't need to run it twice to pick out info from the
    # 'missing' entry
    gene_name_list = [None] * len(gene_ids_list)
    gene_description_list = [None] * len(gene_ids_list)
    entrez_id_list = [None] * len(gene_ids_list)
    uniprot_id_list = [None] * len(gene_ids_list)

    for out_data in ginfo["out"]:
        to_index = gene_ids_list.index(out_data["query"])
        try:
            gene_name_list[to_index] = out_data["name"]
        except KeyError:
            gene_name_list[to_index] = None
        try:
            gene_description_list[to_index] = out_data["symbol"]
        except KeyError:
            gene_description_list[to_index] = None
        try:
            entrez_id_list[to_index] = out_data["entrezgene"]
        except KeyError:
            entrez_id_list[to_index] = None
        try:
            uniprot_id_list[to_index] = out_data["uniprot"]["Swiss-Prot"]
        except KeyError:
            uniprot_id_list[to_index] = None

    for duplicate in ginfo["dup"]:
        indices = [i for i, x in enumerate(gene_ids_list) if x == duplicate[0]]
        dup_info = [
            entry for entry in ginfo["out"] if entry["query"] == duplicate[0]
        ][0]
        for dup_index in indices:
            try:
                gene_name_list[dup_index] = dup_info["name"]
            except KeyError:
                gene_name_list[dup_index] = None
            try:
                gene_description_list[dup_index] = dup_info["symbol"]
            except KeyError:
                gene_description_list[dup_index] = None
            try:
                entrez_id_list[dup_index] = dup_info["entrezgene"]
            except KeyError:
                entrez_id_list[dup_index] = None
            try:
                uniprot_id_list[dup_index] = dup_info["uniprot"]["Swiss-Prot"]
            except KeyError:
                uniprot_id_list[dup_index] = None
    if len(ginfo["dup"]) > 0:
        print("Duplicates were reasigned to the first entry in the list.")

    transcript_df.loc[:, "gene_name"] = gene_name_list
    transcript_df.loc[:, "gene_description"] = gene_description_list
    transcript_df.loc[:, "entrez_id"] = entrez_id_list
    transcript_df.loc[:, "uniprot_id"] = uniprot_id_list

    return transcript_df
