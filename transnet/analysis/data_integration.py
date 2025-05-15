"""
Functionality for integrating experimental data into networks.
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple, Union, Any
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def preprocess_transcriptomics(
    df: pd.DataFrame,
    ensembl_id_col: str = "Ensembl",
    value_cols: List[str] = None,
    log2_transform: bool = False,
    min_exp_threshold: float = 0.0
) -> pd.DataFrame:
    """
    Preprocess transcriptomics data.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe
    ensembl_id_col : str
        Column with Ensembl IDs
    value_cols : List[str], optional
        Columns with expression values
    log2_transform : bool
        Whether to log2 transform the data
    min_exp_threshold : float
        Minimum expression threshold
        
    Returns:
    --------
    pd.DataFrame
        Preprocessed dataframe
    """
    # Create a copy to avoid modifying the original
    result = df.copy()
    
    # Ensure Ensembl IDs are strings
    result[ensembl_id_col] = result[ensembl_id_col].astype(str)
    
    # If value columns not specified, use all numeric columns except the ID column
    if value_cols is None:
        value_cols = [col for col in result.columns 
                      if col != ensembl_id_col and pd.api.types.is_numeric_dtype(result[col])]
    
    # Apply minimum expression threshold
    if min_exp_threshold > 0:
        for col in value_cols:
            result = result[result[col] >= min_exp_threshold]
    
    # Apply log2 transformation if requested
    if log2_transform:
        for col in value_cols:
            # Avoid log2(0)
            result[col] = np.log2(result[col] + 1)
    
    return result

def compute_differential_expression(
    df: pd.DataFrame, 
    control_cols: List[str],
    treatment_cols: List[str],
    id_col: str,
    method: str = "t-test"
) -> pd.DataFrame:
    """
    Compute differential expression between conditions.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe
    control_cols : List[str]
        Columns for control condition
    treatment_cols : List[str]
        Columns for treatment condition
    id_col : str
        Column with IDs
    method : str
        Statistical method ('t-test', 'wilcoxon', or 'fold-change')
        
    Returns:
    --------
    pd.DataFrame
        Dataframe with differential expression results
    """
    from scipy import stats
    
    # Create result dataframe
    result = pd.DataFrame({id_col: df[id_col]})
    
    # Compute fold change
    control_mean = df[control_cols].mean(axis=1)
    treatment_mean = df[treatment_cols].mean(axis=1)
    
    result['fold_change'] = treatment_mean / control_mean
    result['log2_fold_change'] = np.log2(result['fold_change'])
    
    # Compute p-values based on the method
    p_values = []
    
    if method == "t-test":
        for idx, row in df.iterrows():
            control_values = row[control_cols].values
            treatment_values = row[treatment_cols].values
            
            # Skip if there's not enough data
            if len(control_values) < 2 or len(treatment_values) < 2:
                p_values.append(np.nan)
                continue
            
            t_stat, p_val = stats.ttest_ind(treatment_values, control_values)
            p_values.append(p_val)
            
    elif method == "wilcoxon":
        for idx, row in df.iterrows():
            control_values = row[control_cols].values
            treatment_values = row[treatment_cols].values
            
            # Skip if there's not enough data
            if len(control_values) < 2 or len(treatment_values) < 2:
                p_values.append(np.nan)
                continue
            
            try:
                _, p_val = stats.ranksums(treatment_values, control_values)
                p_values.append(p_val)
            except ValueError:
                p_values.append(np.nan)
                
    elif method == "fold-change":
        # No p-values for simple fold change
        p_values = [np.nan] * len(df)
        
    else:
        logger.error(f"Unknown method: {method}")
        p_values = [np.nan] * len(df)
    
    result['p_value'] = p_values
    
    # Compute adjusted p-values (Benjamini-Hochberg)
    valid_p = result['p_value'].dropna()
    
    if len(valid_p) > 0:
        rank = valid_p.rank()
        fdr = valid_p * len(valid_p) / rank
        fdr[fdr > 1] = 1
        
        # Add back to result
        result.loc[valid_p.index, 'adj_p_value'] = fdr
    else:
        result['adj_p_value'] = np.nan
    
    return result

def id_mapping(
    df: pd.DataFrame,
    id_col: str,
    from_type: str,
    to_type: str,
    organism: str = "human"
) -> pd.DataFrame:
    """
    Map IDs from one type to another.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe
    id_col : str
        Column with IDs to map
    from_type : str
        Source ID type
    to_type : str
        Target ID type
    organism : str
        Organism name
        
    Returns:
    --------
    pd.DataFrame
        Dataframe with mapped IDs
    """
    try:
        import mygene
        mg = mygene.MyGeneInfo()
        
        # Extract IDs to map
        ids = df[id_col].tolist()
        
        # Perform mapping
        result = mg.querymany(ids, scopes=from_type, fields=to_type, species=organism)
        
        # Create mapping dictionary
        id_map = {}
        for item in result:
            query_id = item.get('query', '')
            mapped_id = item.get(to_type, None)
            
            if mapped_id:
                id_map[query_id] = mapped_id
        
        # Add mapped IDs to dataframe
        new_col = f"{to_type}_id"
        df[new_col] = df[id_col].map(id_map)
        
        return df
        
    except ImportError:
        logger.error("mygene package is required for ID mapping")
        df[f"{to_type}_id"] = None
        return df

def integrate_experimental_data_into_transnet(
    transnet_obj: Any,
    transcriptomics_df: Optional[pd.DataFrame] = None,
    proteomics_df: Optional[pd.DataFrame] = None,
    metabolomics_df: Optional[pd.DataFrame] = None
) -> None:
    """
    Integrate experimental data into a Transnet object.
    
    Parameters:
    -----------
    transnet_obj : Transnet
        Transnet object
    transcriptomics_df : pd.DataFrame, optional
        Transcriptomics data
    proteomics_df : pd.DataFrame, optional
        Proteomics data
    metabolomics_df : pd.DataFrame, optional
        Metabolomics data
    """
    # Integrate transcriptomics data
    if transcriptomics_df is not None:
        transnet_obj.transcriptome.add_experimental_data(
            input_data=transcriptomics_df,
            transcript_column_name=transcriptomics_df.columns[0],
            id_type="ensembl"
        )
    
    # Integrate proteomics data
    if proteomics_df is not None:
        transnet_obj.proteome.add_experimental_data(
            input_data=proteomics_df,
            protein_column_name=proteomics_df.columns[0],
            id_type="uniprot"
        )
    
    # Integrate metabolomics data
    if metabolomics_df is not None:
        transnet_obj.metabolome.add_experimental_data(
            input_data=metabolomics_df,
            metabolite_column_name=metabolomics_df.columns[0],
            id_type="pubchem"
        )
        