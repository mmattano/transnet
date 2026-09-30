"""
Functionality for integrating experimental data into networks.
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple, Union, Any
import logging

logger = logging.getLogger(__name__)


def _finite(series) -> np.ndarray:
    """Finite float values from one row of measurements."""
    values = series.to_numpy(dtype=float)
    return values[np.isfinite(values)]


def compute_differential_expression(
    df: pd.DataFrame, 
    control_cols: List[str],
    treatment_cols: List[str],
    id_col: str,
    method: str = "t-test",
    data_is_log: bool = False,
) -> pd.DataFrame:
    """
    Compute differential expression between conditions.
    
    Parameters
    ----------
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
        
    Returns
    -------
    pd.DataFrame
        Dataframe with differential expression results
    """
    from scipy import stats
    
    # Create result dataframe
    result = pd.DataFrame({id_col: df[id_col]})
    
    # Coerce once: a frame carrying an identifier column is mixed dtype, and
    # `iterrows` then yields object Series that scipy cannot operate on.
    values = df[list(control_cols) + list(treatment_cols)].apply(
        pd.to_numeric, errors="coerce"
    )
    control_values_all = values[list(control_cols)]
    treatment_values_all = values[list(treatment_cols)]

    control_mean = control_values_all.mean(axis=1)
    treatment_mean = treatment_values_all.mean(axis=1)

    if data_is_log:
        result['log2_fold_change'] = treatment_mean - control_mean
        result['fold_change'] = np.power(2.0, result['log2_fold_change'])
    else:
        with np.errstate(divide='ignore', invalid='ignore'):
            fold_change = treatment_mean / control_mean
            result['fold_change'] = fold_change
            result['log2_fold_change'] = np.log2(fold_change.where(fold_change > 0))
    
    # Compute p-values based on the method
    p_values = []
    
    if method == "t-test":
        for position in range(len(df)):
            control_values = _finite(control_values_all.iloc[position])
            treatment_values = _finite(treatment_values_all.iloc[position])

            # Skip if there's not enough data
            if len(control_values) < 2 or len(treatment_values) < 2:
                p_values.append(np.nan)
                continue

            _, p_val = stats.ttest_ind(treatment_values, control_values)
            p_values.append(p_val)
            
    elif method == "wilcoxon":
        for position in range(len(df)):
            control_values = _finite(control_values_all.iloc[position])
            treatment_values = _finite(treatment_values_all.iloc[position])

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
    
    Parameters
    ----------
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
        
    Returns
    -------
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


def compute_timecourse_statistics(
    df: pd.DataFrame,
    groups: List[str],
    method: str = 'anova',
    trend_test: bool = False,
    continuous_x: Optional[List[float]] = None,
) -> pd.DataFrame:
    """
    Compute per-feature time-course statistics across ≥2 groups.

    Parameters
    ----------
    df : pd.DataFrame
        Features × samples (rows = features, columns = samples).
    groups : List[str]
        Group label for each sample; length must equal ``len(df.columns)``.
    method : str
        ``'anova'`` (one-way F-test) or ``'kruskal'`` (Kruskal-Wallis H-test).
    trend_test : bool
        If True, also run a linear trend test via ``scipy.stats.linregress``
        using ``continuous_x`` as the x-values.  Adds ``slope``, ``r_value``,
        and ``trend_pval`` columns to the output.
    continuous_x : List[float], optional
        Numeric value assigned to each *unique* group (order matters).
        Required when ``trend_test=True``.  Example: ``[0, 1, 2, 4, 8]`` for
        five timepoints.  Must align with the sorted unique values of ``groups``.

    Returns
    -------
    pd.DataFrame
        Columns: ``feature``, ``stat``, ``pval``, ``fdr``, and optionally
        ``slope``, ``r_value``, ``trend_pval``.
    """
    from scipy import stats
    try:
        from statsmodels.stats.multitest import multipletests
        _have_sm = True
    except ImportError:
        _have_sm = False

    groups_arr = np.asarray(groups)
    if len(groups_arr) != len(df.columns):
        raise ValueError(
            f"Length of groups ({len(groups_arr)}) must equal number of "
            f"samples ({len(df.columns)})"
        )

    unique_groups = list(dict.fromkeys(groups_arr))  # preserve insertion order

    if trend_test:
        if continuous_x is None:
            raise ValueError("continuous_x is required when trend_test=True")
        if len(continuous_x) != len(unique_groups):
            raise ValueError(
                f"continuous_x length ({len(continuous_x)}) must equal number "
                f"of unique groups ({len(unique_groups)})"
            )
        group_x = dict(zip(unique_groups, continuous_x))

    # Use positional numpy indexing to avoid issues when df.index contains
    # duplicate labels (df.loc[label] returns a 2-D DataFrame for duplicates,
    # making boolean masks and .values 2-D).
    _data = df.to_numpy().astype(float)   # (n_features, n_samples)
    _group_masks = {g: groups_arr == g for g in unique_groups}

    if trend_test:
        x_vals_full = np.array([group_x[g] for g in groups_arr])

    rows = []
    for feat_i, feature in enumerate(df.index):
        _row_data = _data[feat_i]   # 1-D, length n_samples
        vals_by_group = []
        for g in unique_groups:
            v = _row_data[_group_masks[g]]
            v = v[~np.isnan(v)]
            if len(v) > 0:
                vals_by_group.append(v)

        if len(vals_by_group) < 2:
            stat, pval = np.nan, np.nan
        else:
            try:
                if method == 'anova':
                    stat, pval = stats.f_oneway(*vals_by_group)
                elif method == 'kruskal':
                    stat, pval = stats.kruskal(*vals_by_group)
                else:
                    raise ValueError(f"Unknown method: {method}")
            except Exception:
                stat, pval = np.nan, np.nan

        row: Dict[str, Any] = {'feature': feature, 'stat': stat, 'pval': pval}

        if trend_test:
            y_vals = _row_data.copy()
            mask = ~np.isnan(y_vals)
            if mask.sum() > 2:
                res = stats.linregress(x_vals_full[mask], y_vals[mask])
                row['slope'] = res.slope
                row['r_value'] = res.rvalue
                row['trend_pval'] = res.pvalue
            else:
                row['slope'] = np.nan
                row['r_value'] = np.nan
                row['trend_pval'] = np.nan

        rows.append(row)

    result = pd.DataFrame(rows)

    # BH-FDR correction on the omnibus p-values
    valid_mask = result['pval'].notna()
    if valid_mask.sum() > 0:
        pvals_valid = result.loc[valid_mask, 'pval'].values
        if _have_sm:
            _, fdr_vals, _, _ = multipletests(pvals_valid, method='fdr_bh')
        else:
            # Manual BH fallback
            n = len(pvals_valid)
            order = np.argsort(pvals_valid)
            ranks = np.empty(n)
            ranks[order] = np.arange(1, n + 1)
            fdr_vals = np.minimum(pvals_valid * n / ranks, 1.0)
            fdr_sorted = fdr_vals[order]
            for i in range(n - 2, -1, -1):
                fdr_sorted[i] = min(fdr_sorted[i], fdr_sorted[i + 1])
            fdr_vals[order] = fdr_sorted
        result.loc[valid_mask, 'fdr'] = fdr_vals
    else:
        result['fdr'] = np.nan

    if trend_test and valid_mask.sum() > 0:
        trend_valid = result['trend_pval'].notna()
        if trend_valid.sum() > 0:
            tp_valid = result.loc[trend_valid, 'trend_pval'].values
            if _have_sm:
                _, trend_fdr, _, _ = multipletests(tp_valid, method='fdr_bh')
            else:
                n = len(tp_valid)
                order = np.argsort(tp_valid)
                ranks = np.empty(n)
                ranks[order] = np.arange(1, n + 1)
                trend_fdr = np.minimum(tp_valid * n / ranks, 1.0)
            result.loc[trend_valid, 'trend_fdr'] = trend_fdr

    return result.sort_values('pval').reset_index(drop=True)


def cluster_temporal_trajectories(
    df: pd.DataFrame,
    timepoints: List[str],
    n_clusters: int = 5,
    method: str = 'kmeans',
    random_state: int = 42,
) -> Tuple[pd.Series, pd.DataFrame]:
    """
    Cluster features by their mean expression profile across ordered timepoints.

    Parameters
    ----------
    df : pd.DataFrame
        Features × samples (rows = features, columns = samples).
    timepoints : List[str]
        Timepoint label for each sample; length must equal ``len(df.columns)``.
    n_clusters : int
        Number of trajectory clusters.
    method : str
        Currently only ``'kmeans'`` is supported.
    random_state : int
        Random seed.

    Returns
    -------
    assignments : pd.Series
        Cluster id (0-indexed) indexed by feature name.
    centroids : pd.DataFrame
        Mean z-scored profiles for each cluster (clusters × unique timepoints).
    """
    from sklearn.cluster import KMeans
    from sklearn.preprocessing import StandardScaler

    tp_arr = np.asarray(timepoints)
    if len(tp_arr) != len(df.columns):
        raise ValueError(
            f"Length of timepoints ({len(tp_arr)}) must equal number of "
            f"samples ({len(df.columns)})"
        )

    unique_tp = list(dict.fromkeys(tp_arr))

    # Mean per feature per timepoint, then z-score across timepoints
    mean_profiles = np.vstack([
        df.loc[:, tp_arr == tp].mean(axis=1).values
        for tp in unique_tp
    ]).T  # (n_features, n_timepoints)

    scaler = StandardScaler()
    profiles_scaled = scaler.fit_transform(mean_profiles)

    if method == 'kmeans':
        km = KMeans(n_clusters=n_clusters, random_state=random_state, n_init=10)
        labels = km.fit_predict(profiles_scaled)
        centroid_data = scaler.inverse_transform(km.cluster_centers_)
    else:
        raise ValueError(f"Unknown clustering method: {method}")

    assignments = pd.Series(labels, index=df.index, name='cluster')

    centroids = pd.DataFrame(
        centroid_data,
        columns=unique_tp,
        index=[f'Cluster{i}' for i in range(n_clusters)],
    )

    return assignments, centroids
