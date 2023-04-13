"""Functions to access ChIP-Atlas database."""

__all__ = [
    "get_ChIP_data",
    "list_ChIP_cell_type_classes",
    "list_ChIP_cell_types",
    "get_ChIP_exps",
]

import pandas as pd
import numpy as np
import tempfile
import requests
from functools import reduce


def get_ChIP_data(genome: str = "hg38", distance: int = 5):
    """Get ChIP-Atlas data for a given genome and distance.
    
    Parameters
    ----------
    genome : str, optional
        Genome to use, by default human, 'hg38'. Check the available genomes in the ChIP-Atlas website.
    distance : int, optional
        Distance in kb, by default 5. Can be 1, 5 or 10.
    
    Returns
    -------
    scores_unsort : pd.DataFrame
        Dataframe with the scores for each gene and each experiment.
    gene_to_experiment : dict
        Dictionary with the gene as key and the experiments that were performed for the corresponding genes as values.
    """

    scores_unsort = pd.DataFrame()
    gene_to_experiment = {}

    # Creating a temporary directory to store the downloaded files
    # with tempfile.TemporaryDirectory() as tmpdirname:
    # Accessing ChIP-Atlas data to get all the gene interactions for the given genome
    url_genome = (
        "http://dbarchive.biosciencedbc.jp/kyushu-u/metadata/analysisList.tab"
    )
    genome_data = requests.get(url_genome)
    analysisList = pd.DataFrame(
        [
            x.split("\t")
            for x in genome_data.content.decode("utf-8").split("\n")
        ]
    )
    reduced_analysisList = analysisList[
        analysisList[analysisList.columns[-1]] == genome
    ]
    proteins = reduced_analysisList[analysisList.columns[0]].to_list()
    print(len(proteins), "transcription factors found for", genome)
    # Now that we have the list of proteins, we can download the data for each gene
    all_dfs = []
    for i, protein in enumerate(proteins):
        # modified version of the established version
        # Downloading the data for each gene
        url_exp_data = f"https://chip-atlas.dbcls.jp/data/{genome}/target/{protein}.{distance}.tsv"
        exp_data = requests.get(url_exp_data)
        gene_exp_data = pd.DataFrame(
            [
                x.split("\t")
                for x in exp_data.content.decode("utf-8").split("\n")
            ]
        )
        gene_exp_data.columns = gene_exp_data.iloc[0]
        gene_exp_data = gene_exp_data.drop(gene_exp_data.index[0])
        if gene_exp_data.columns[0] == "Target_genes":
            # I am writing it to a list to merge them later (this is a bit faster than merging them one by one)
            gene_exp_data.set_index("Target_genes", inplace=True)
            all_dfs.append(gene_exp_data)
        else:
            # ChIP-Atlas doesn't have data for some genes but still returns an empty file
            print(
                f"ChIP-Atlas does not have data for {protein}"
            )  # This is just to keep track of the genes that don't have data
    for i, df in enumerate(all_dfs):
        if i == 0:
            df_all_ChIP_data = df
        else:
            df_all_ChIP_data = df_all_ChIP_data.join(
                df, how="left", lsuffix=f"{i-1}"
            )

    # Now we have to clean the data
    # First, we have to remove some of the columns that are not experiments
    experiments_original = list(df_all_ChIP_data.columns)
    experiments_ids = []
    exp_indices = []
    current_gene = ""
    gene_to_experiment = {}
    for i, x in enumerate(experiments_original):
        if "|Average" in x:
            current_gene = x.split("|Average")[0]
        if "|" in x and "|Average" not in x:
            exp = x.split("|")[0]
            experiments_ids.append(exp)
            if current_gene in gene_to_experiment:
                gene_to_experiment[current_gene].append(exp)
            else:
                gene_to_experiment[current_gene] = [exp]
            # We need to keep the indices of the experiments to filter the scores later
            exp_indices.append(i)
    # Now we have to remove the duplicates
    experiments_unique_ids, indices_unique = np.unique(
        experiments_ids, return_index=True
    )
    experiments_unique_ids = list(experiments_unique_ids)
    indices_unique = list(indices_unique)
    # Now we're tracing back the indices of the experiments to keep
    indices_to_keep = [i for i in exp_indices if i in indices_unique]
    experiments_ids_to_keep = [
        x.split("|")[0]
        for i, x in enumerate(experiments_original)
        if i in indices_to_keep
    ]

    # Now we can filter the scores
    scores_unsort_pre = df_all_ChIP_data.iloc[:, 2:]
    scores_unsort = scores_unsort_pre.iloc[:, indices_to_keep]
    scores_unsort.columns = experiments_ids_to_keep
    scores_unsort.dropna(axis=0, how="all", inplace=True)
    # And now we have a new dataframe with the scores and the corresponding experiments
    # The rows still correspond to the genes

    return scores_unsort, gene_to_experiment


def _get_ChIP_experiment_info(genome: str = "hg38"):
    """Get the list of experiments for a given genome.
    
    Parameters
    ----------
    genome : str, optional
        Genome to use, by default human, 'hg38'. Check the available genomes in the ChIP-Atlas website.
    
    Returns
    -------
    genome_specific_exps : pd.DataFrame
        Dataframe with the experiments for the given genome and associted information.
    """

    # Getting the list of experiments registered in ChIP-Atlas
    url_exp_list = "http://dbarchive.biosciencedbc.jp/kyushu-u/metadata/experimentList.tab"
    exp_list = requests.get(url_exp_list)
    exp_list_df = pd.DataFrame(
        [x.split("\t") for x in exp_list.content.decode("utf-8").split("\n")]
    )
    # Getting the experiments for the given genome
    genome_specific_exps = exp_list_df[
        exp_list_df[exp_list_df.columns[1]] == genome
    ]

    return genome_specific_exps


def list_ChIP_cell_type_classes(genome: str = "hg38"):
    """Get the cell type classes for a given genome.

    Parameters
    ----------
    genome : str, optional
        Genome to use, by default human, 'hg38'. Check the available genomes in the ChIP-Atlas website.
    
    Returns
    -------
    list
        List with the cell type classes for the given genome.
    """

    return list(set(_get_ChIP_experiment_info(genome)[4]))


def list_ChIP_cell_types(genome: str = "hg38", cell_type_class: str = None):
    """Get the cell types for a given genome.

    Parameters
    ----------
    genome : str, optional
        Genome to use, by default human, 'hg38'. Check the available genomes in the ChIP-Atlas website.
    cell_type_class : str, optional
        Cell type class to use, by default None. Check list_ChIP_cell_type_classes() for the available cell type classes.
    
    Returns
    -------
    list
        List with the cell types for the given genome.
    """
    exp_info = _get_ChIP_experiment_info(genome)
    if cell_type_class is None:
        return list(set(exp_info[5]))
    else:
        return list(set(exp_info[exp_info[4] == cell_type_class][5]))


def get_ChIP_exps(
    genome: str = "hg38", cell_type_class: str = None, cell_type: str = None
):
    """Get the experiments for a given genome, cell type class and cell type.

    Parameters
    ----------
    genome : str, optional
        Genome to use, by default human, 'hg38'. Check the available genomes in the ChIP-Atlas website.
    cell_type_class : str, optional
        Cell type class to use, by default None. Check list_ChIP_cell_type_classes() for the available cell type classes.
    cell_type : str, optional
        Cell type to use, by default None. Check list_ChIP_cell_types() for the available cell types.
    
    Returns
    -------
    pd.DataFrame
        List with the experiments for the given genome, optionally filtered by cell type class and cell type.
    """
    exp_info = _get_ChIP_experiment_info(genome)
    if cell_type_class is None:
        return exp_info
    elif cell_type is None:
        return exp_info[exp_info[4] == cell_type_class]
    else:
        return exp_info[
            (exp_info[4] == cell_type_class) & (exp_info[5] == cell_type)
        ]
