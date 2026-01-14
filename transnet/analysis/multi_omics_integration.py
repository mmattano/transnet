"""
Statistical multi-omics integration methods for TransNet.

Implements dimensionality reduction and integration approaches using scikit-learn:
- Non-negative Matrix Factorization (NMF)
- Factor Analysis (FA)
- Principal Component Analysis (PCA)
- Missing data imputation strategies
- Cross-omics factor correlation analysis
- Network-aware enrichment testing
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Union
from sklearn.decomposition import NMF, FactorAnalysis, PCA
from sklearn.impute import SimpleImputer, KNNImputer
from sklearn.preprocessing import StandardScaler
from scipy import stats
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class MultiOmicsIntegrator:
    """
    Integrates multi-omics data using dimensionality reduction.
    
    Supports three integration strategies:
    1. Early fusion: Concatenate features then reduce
    2. Late fusion: Reduce each omics separately, then combine factors
    3. Concatenation: Stack features without reduction
    """
    
    def __init__(
        self,
        n_components: int = 10,
        method: str = 'nmf',
        integration_strategy: str = 'early',
        random_state: int = 42
    ):
        """
        Parameters:
        -----------
        n_components : int
            Number of latent factors/components
        method : str
            'nmf', 'pca', or 'fa' (factor analysis)
        integration_strategy : str
            'early', 'late', or 'concatenation'
        random_state : int
            Random seed for reproducibility
        """
        self.n_components = n_components
        self.method = method.lower()
        self.integration_strategy = integration_strategy
        self.random_state = random_state
        
        self.models_ = {}
        self.scalers_ = {}
        self.feature_names_ = {}
        self.sample_names_ = None
        self.factors_ = None
        self.loadings_ = {}
        self.variance_explained_ = {}
        
    def fit_transform(
        self,
        omics_data: Dict[str, pd.DataFrame],
        imputation_strategy: str = 'knn'
    ) -> pd.DataFrame:
        """
        Fit the integration model and transform data.
        
        Parameters:
        -----------
        omics_data : Dict[str, pd.DataFrame]
            Dictionary mapping omics layer names to DataFrames
            Rows = samples, Columns = features
        imputation_strategy : str
            'knn', 'mean', 'median', or 'none'
            
        Returns:
        --------
        pd.DataFrame
            Integrated factors (samples x components)
        """
        
        if not omics_data:
            raise ValueError("omics_data cannot be empty")
        
        # Validate sample alignment
        sample_names = None
        for layer_name, df in omics_data.items():
            if sample_names is None:
                sample_names = df.index
            elif not sample_names.equals(df.index):
                raise ValueError(f"Sample names not aligned in {layer_name}")
        
        self.sample_names_ = sample_names
        
        # Preprocess each omics layer
        processed_data = {}
        for layer_name, df in omics_data.items():
            logger.info(f"Processing {layer_name}: {df.shape}")
            
            self.feature_names_[layer_name] = df.columns.tolist()
            
            # Handle missing values
            if imputation_strategy != 'none':
                df = self._impute_missing(df, imputation_strategy)
            
            # Ensure non-negative for NMF
            if self.method == 'nmf':
                df = df.clip(lower=0)
            
            # Standardize
            scaler = StandardScaler()
            data_scaled = scaler.fit_transform(df)
            self.scalers_[layer_name] = scaler
            
            processed_data[layer_name] = pd.DataFrame(
                data_scaled,
                index=df.index,
                columns=df.columns
            )
        
        # Perform integration
        if self.integration_strategy == 'early':
            self.factors_ = self._early_fusion(processed_data)
        elif self.integration_strategy == 'late':
            self.factors_ = self._late_fusion(processed_data)
        elif self.integration_strategy == 'concatenation':
            self.factors_ = self._concatenation(processed_data)
        else:
            raise ValueError(f"Unknown strategy: {self.integration_strategy}")
        
        return self.factors_
    
    def transform(
        self,
        omics_data: Dict[str, pd.DataFrame]
    ) -> pd.DataFrame:
        """
        Transform new data using fitted models.
        
        Parameters:
        -----------
        omics_data : Dict[str, pd.DataFrame]
            New omics data to transform
            
        Returns:
        --------
        pd.DataFrame
            Transformed factors
        """
        
        if not self.models_:
            raise ValueError("Model not fitted yet. Call fit_transform() first.")
        
        # Preprocess using fitted scalers
        processed_data = {}
        for layer_name, df in omics_data.items():
            if layer_name not in self.scalers_:
                raise ValueError(f"Layer {layer_name} not in fitted model")
            
            data_scaled = self.scalers_[layer_name].transform(df)
            processed_data[layer_name] = pd.DataFrame(
                data_scaled,
                index=df.index,
                columns=df.columns
            )
        
        # Transform based on strategy
        if self.integration_strategy == 'early':
            concat_df = pd.concat(
                [processed_data[k] for k in sorted(processed_data.keys())],
                axis=1
            )
            factors = self.models_['integrated'].transform(concat_df.values)
            
            return pd.DataFrame(
                factors,
                index=concat_df.index,
                columns=[f'Factor{i+1}' for i in range(factors.shape[1])]
            )
        
        elif self.integration_strategy == 'late':
            all_factors = []
            
            for layer_name, df in processed_data.items():
                factors = self.models_[layer_name].transform(df.values)
                factor_df = pd.DataFrame(
                    factors,
                    index=df.index,
                    columns=[f'{layer_name}_Factor{i+1}' for i in range(factors.shape[1])]
                )
                all_factors.append(factor_df)
            
            return pd.concat(all_factors, axis=1)
        
        else:  # concatenation
            concat_df = pd.concat(
                [processed_data[k] for k in sorted(processed_data.keys())],
                axis=1
            )
            return concat_df
    
    def _impute_missing(
        self,
        df: pd.DataFrame,
        strategy: str
    ) -> pd.DataFrame:
        """Impute missing values."""
        
        if not df.isnull().any().any():
            return df
        
        logger.info(f"Imputing with strategy: {strategy}")
        
        if strategy == 'knn':
            imputer = KNNImputer(n_neighbors=5)
        elif strategy in ['mean', 'median']:
            imputer = SimpleImputer(strategy=strategy)
        else:
            raise ValueError(f"Unknown imputation strategy: {strategy}")
        
        data_imputed = imputer.fit_transform(df)
        
        return pd.DataFrame(
            data_imputed,
            index=df.index,
            columns=df.columns
        )
    
    def _get_model(self, n_features: int):
        """Initialize dimensionality reduction model."""
        
        n_comp = min(self.n_components, n_features)
        
        if self.method == 'nmf':
            return NMF(
                n_components=n_comp,
                init='nndsvda',
                random_state=self.random_state,
                max_iter=500
            )
        elif self.method == 'fa':
            return FactorAnalysis(
                n_components=n_comp,
                random_state=self.random_state,
                max_iter=500
            )
        elif self.method == 'pca':
            return PCA(
                n_components=n_comp,
                random_state=self.random_state
            )
        else:
            raise ValueError(f"Unknown method: {self.method}")
    
    def _early_fusion(
        self,
        processed_data: Dict[str, pd.DataFrame]
    ) -> pd.DataFrame:
        """Early fusion: concatenate features then reduce."""
        
        concat_df = pd.concat(
            [processed_data[k] for k in sorted(processed_data.keys())],
            axis=1
        )
        
        logger.info(f"Early fusion shape: {concat_df.shape}")
        
        model = self._get_model(concat_df.shape[1])
        factors = model.fit_transform(concat_df.values)
        
        self.models_['integrated'] = model
        
        # Track variance explained
        if hasattr(model, 'explained_variance_ratio_'):
            self.variance_explained_['integrated'] = model.explained_variance_ratio_
        
        # Store loadings for each omics layer
        start_idx = 0
        for layer_name in sorted(processed_data.keys()):
            n_features = len(self.feature_names_[layer_name])
            end_idx = start_idx + n_features
            
            layer_loadings = model.components_[:, start_idx:end_idx]
            self.loadings_[layer_name] = pd.DataFrame(
                layer_loadings.T,
                index=self.feature_names_[layer_name],
                columns=[f'Factor{i+1}' for i in range(layer_loadings.shape[0])]
            )
            
            start_idx = end_idx
        
        return pd.DataFrame(
            factors,
            index=self.sample_names_,
            columns=[f'Factor{i+1}' for i in range(factors.shape[1])]
        )
    
    def _late_fusion(
        self,
        processed_data: Dict[str, pd.DataFrame]
    ) -> pd.DataFrame:
        """Late fusion: reduce each omics separately, then concatenate factors."""
        
        all_factors = []
        
        for layer_name, df in processed_data.items():
            logger.info(f"Late fusion for {layer_name}")
            
            model = self._get_model(df.shape[1])
            factors = model.fit_transform(df.values)
            
            self.models_[layer_name] = model
            
            # Track variance explained
            if hasattr(model, 'explained_variance_ratio_'):
                self.variance_explained_[layer_name] = model.explained_variance_ratio_
            
            # Store loadings
            self.loadings_[layer_name] = pd.DataFrame(
                model.components_.T,
                index=self.feature_names_[layer_name],
                columns=[f'{layer_name}_Factor{i+1}' for i in range(factors.shape[1])]
            )
            
            factor_df = pd.DataFrame(
                factors,
                index=df.index,
                columns=[f'{layer_name}_Factor{i+1}' for i in range(factors.shape[1])]
            )
            
            all_factors.append(factor_df)
        
        combined = pd.concat(all_factors, axis=1)
        logger.info(f"Late fusion combined shape: {combined.shape}")
        
        return combined
    
    def _concatenation(
        self,
        processed_data: Dict[str, pd.DataFrame]
    ) -> pd.DataFrame:
        """Concatenation: stack features without dimensionality reduction."""
        
        concat_df = pd.concat(
            [processed_data[k] for k in sorted(processed_data.keys())],
            axis=1
        )
        
        logger.info(f"Concatenation shape: {concat_df.shape}")
        
        return concat_df
    
    def compute_cross_omics_correlations(
        self,
        omics1: str,
        omics2: str,
        method: str = 'pearson'
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Compute correlations between factors across omics layers.
        
        Only works with late fusion strategy where each omics has separate factors.
        
        Parameters:
        -----------
        omics1, omics2 : str
            Names of omics layers to compare
        method : str
            'pearson' or 'spearman'
            
        Returns:
        --------
        Tuple[pd.DataFrame, pd.DataFrame]
            (correlation_matrix, pvalue_matrix)
        """
        
        if self.integration_strategy != 'late':
            raise ValueError("Cross-omics correlations only available with late fusion")
        
        if self.factors_ is None:
            raise ValueError("Model not fitted yet")
        
        # Get factors for each omics
        omics1_cols = [c for c in self.factors_.columns if c.startswith(f'{omics1}_')]
        omics2_cols = [c for c in self.factors_.columns if c.startswith(f'{omics2}_')]
        
        if not omics1_cols or not omics2_cols:
            raise ValueError(f"Could not find factors for {omics1} or {omics2}")
        
        factors1 = self.factors_[omics1_cols]
        factors2 = self.factors_[omics2_cols]
        
        # Compute correlations
        n_factors1 = len(omics1_cols)
        n_factors2 = len(omics2_cols)
        
        corr_matrix = np.zeros((n_factors1, n_factors2))
        pval_matrix = np.zeros((n_factors1, n_factors2))
        
        for i in range(n_factors1):
            for j in range(n_factors2):
                if method == 'pearson':
                    r, p = stats.pearsonr(factors1.iloc[:, i], factors2.iloc[:, j])
                elif method == 'spearman':
                    r, p = stats.spearmanr(factors1.iloc[:, i], factors2.iloc[:, j])
                else:
                    raise ValueError(f"Unknown method: {method}")
                
                corr_matrix[i, j] = r
                pval_matrix[i, j] = p
        
        corr_df = pd.DataFrame(
            corr_matrix,
            index=omics1_cols,
            columns=omics2_cols
        )
        
        pval_df = pd.DataFrame(
            pval_matrix,
            index=omics1_cols,
            columns=omics2_cols
        )
        
        return corr_df, pval_df
    
    def get_top_features(
        self,
        layer_name: str,
        factor_idx: int,
        n_top: int = 20,
        return_weights: bool = True
    ) -> Union[List[str], pd.DataFrame]:
        """
        Get top contributing features for a factor.
        
        Parameters:
        -----------
        layer_name : str
            Omics layer name
        factor_idx : int
            Factor index (0-indexed)
        n_top : int
            Number of top features to return
        return_weights : bool
            If True, return DataFrame with weights
            
        Returns:
        --------
        List[str] or pd.DataFrame
            Top features (and weights if requested)
        """
        
        if layer_name not in self.loadings_:
            raise ValueError(f"Layer {layer_name} not found in loadings")
        
        loadings = self.loadings_[layer_name]
        
        if factor_idx >= len(loadings.columns):
            raise ValueError(f"Factor index {factor_idx} out of range")
        
        factor_col = loadings.columns[factor_idx]
        
        abs_loadings = loadings[factor_col].abs()
        top_features = abs_loadings.nlargest(n_top)
        
        if return_weights:
            result = pd.DataFrame({
                'Feature': top_features.index,
                'Loading': loadings.loc[top_features.index, factor_col],
                'AbsLoading': top_features.values
            })
            return result.reset_index(drop=True)
        else:
            return top_features.index.tolist()
    
    def explained_variance_ratio(self) -> Dict[str, np.ndarray]:
        """
        Get explained variance ratio for each model.
        Only works for PCA models.
        """
        
        if self.method != 'pca':
            logger.warning("Explained variance only available for PCA")
            return {}
        
        return self.variance_explained_


def integrate_with_transnet(
    transnet_obj,
    omics_data: Dict[str, pd.DataFrame],
    integrator: Optional[MultiOmicsIntegrator] = None,
    **integrator_kwargs
) -> Tuple[pd.DataFrame, MultiOmicsIntegrator]:
    """
    Integrate multi-omics data and map factors to TransNet graph.
    
    Parameters:
    -----------
    transnet_obj : Transnet
        TransNet object with network structure
    omics_data : Dict[str, pd.DataFrame]
        Dictionary of omics DataFrames (samples x features)
    integrator : MultiOmicsIntegrator, optional
        Pre-configured integrator
    **integrator_kwargs
        Arguments for MultiOmicsIntegrator if not provided
        
    Returns:
    --------
    Tuple[pd.DataFrame, MultiOmicsIntegrator]
        (factors, fitted_integrator)
    """
    
    if integrator is None:
        integrator = MultiOmicsIntegrator(**integrator_kwargs)
    
    # Fit and transform
    factors = integrator.fit_transform(omics_data)
    
    logger.info(f"Generated {factors.shape[1]} factors for {factors.shape[0]} samples")
    
    # Map loadings to network nodes
    G = transnet_obj.generate_graph()
    
    for layer_name, loadings_df in integrator.loadings_.items():
        logger.info(f"Mapping {layer_name} loadings to network")
        
        for feature in loadings_df.index:
            if feature in G.nodes():
                for col in loadings_df.columns:
                    G.nodes[feature][f'loading_{col}'] = loadings_df.loc[feature, col]
    
    transnet_obj.graph = G
    
    return factors, integrator


def compute_factor_network_enrichment(
    G,
    factor_scores: pd.Series,
    pathway_attribute: str = 'pathway',
    top_percentile: float = 0.1
) -> pd.DataFrame:
    """
    Test for enrichment of high-scoring nodes in network pathways/modules.
    
    Uses hypergeometric test to identify pathways enriched in top factor-loading nodes.
    
    Parameters:
    -----------
    G : nx.Graph
        NetworkX graph with pathway annotations
    factor_scores : pd.Series
        Factor loadings for nodes (index=node_id, values=loadings)
    pathway_attribute : str
        Node attribute containing pathway membership
    top_percentile : float
        Percentile threshold for "high-scoring" nodes
        
    Returns:
    --------
    pd.DataFrame
        Enrichment results for each pathway
    """
    from scipy.stats import hypergeom
    
    # Define high-scoring nodes
    threshold = np.percentile(factor_scores.dropna().abs(), (1 - top_percentile) * 100)
    high_nodes = set(factor_scores[factor_scores.abs() > threshold].index)
    
    # Get pathway membership
    pathways = {}
    for node in G.nodes():
        pw = G.nodes[node].get(pathway_attribute)
        if pw:
            if isinstance(pw, list):
                for p in pw:
                    if p not in pathways:
                        pathways[p] = set()
                    pathways[p].add(node)
            else:
                if pw not in pathways:
                    pathways[pw] = set()
                pathways[pw].add(node)
    
    # Hypergeometric test
    N = len(factor_scores)
    K = len(high_nodes)
    
    results = []
    for pathway, pathway_nodes in pathways.items():
        n = len(pathway_nodes)
        k = len(high_nodes & pathway_nodes)
        
        if k == 0:
            continue
        
        # P(X >= k)
        pval = hypergeom.sf(k - 1, N, K, n)
        
        results.append({
            'pathway': pathway,
            'pathway_size': n,
            'overlap': k,
            'expected': (K * n) / N,
            'enrichment': k / ((K * n) / N) if (K * n) > 0 else 0,
            'pvalue': pval
        })
    
    if not results:
        return pd.DataFrame()
    
    df = pd.DataFrame(results)
    
    # FDR correction
    from statsmodels.stats.multitest import multipletests
    if len(df) > 0:
        _, df['fdr'], _, _ = multipletests(df['pvalue'], method='fdr_bh')
    
    return df.sort_values('pvalue')


def perform_factor_enrichment(
    integrator: MultiOmicsIntegrator,
    layer_name: str,
    factor_idx: int,
    annotations: Dict[str, List[str]],
    n_top: int = 100
) -> pd.DataFrame:
    """
    Perform pathway/gene set enrichment on factor loadings.
    
    Parameters:
    -----------
    integrator : MultiOmicsIntegrator
        Fitted integrator
    layer_name : str
        Omics layer name
    factor_idx : int
        Factor index (0-indexed)
    annotations : Dict[str, List[str]]
        Mapping of pathway/gene set names to feature lists
    n_top : int
        Number of top features to consider
        
    Returns:
    --------
    pd.DataFrame
        Enrichment results
    """
    
    from transnet.analysis.network_analysis import enrichment_analysis
    
    # Get top features
    top_features_df = integrator.get_top_features(
        layer_name, factor_idx, n_top=n_top, return_weights=True
    )
    top_features = top_features_df['Feature'].tolist()
    
    # Get all features as background
    all_features = integrator.feature_names_[layer_name]
    
    # Perform enrichment
    enrichment_results = enrichment_analysis(
        node_set=top_features,
        background_set=all_features,
        annotations=annotations,
        method='hypergeometric'
    )
    
    return enrichment_results