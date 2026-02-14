"""
Statistical multi-omics integration methods for TransNet.

Implements dimensionality reduction and integration approaches using scikit-learn:
- Non-negative Matrix Factorization (NMF)
- Factor Analysis (FA)
- Principal Component Analysis (PCA)
- Canonical Correlation Analysis (CCA)
- Missing data imputation strategies (including network-based)
- Cross-omics factor correlation analysis
- Network-aware enrichment testing
- Network-guided correlation and propagation (TransNet's unique capability)
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Union
from sklearn.decomposition import NMF, FactorAnalysis, PCA
from sklearn.cross_decomposition import CCA
from sklearn.impute import SimpleImputer, KNNImputer
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from scipy import stats
from scipy.stats import pearsonr, spearmanr
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
    
    Also provides network-aware methods that leverage biological network topology.
    """
    
    def __init__(
        self,
        n_components: int = 10,
        method: str = 'nmf',
        integration_strategy: str = 'early',
        random_state: int = 42,
        network=None
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
        network : Transnet, optional
            Pre-built network for network-guided methods
        """
        self.n_components = n_components
        self.method = method.lower()
        self.integration_strategy = integration_strategy
        self.random_state = random_state
        self.network = network
        
        self.models_ = {}
        self.scalers_ = {}
        self.feature_names_ = {}
        self.sample_names_ = None
        self.factors_ = None
        self.loadings_ = {}
        self.variance_explained_ = {}
        
        # Build network edge index if network provided
        if network is not None:
            logger.info(f"Initialized with network: {network.name}")
            if not hasattr(network, '_edge_index_built') or not network._edge_index_built:
                logger.info("Building cross-layer edge index")
                network.build_cross_layer_edge_index()
        
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
            'knn', 'mean', 'median', 'network', or 'none'
            
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
                df = self._impute_missing(df, imputation_strategy, layer_name)
            
            # Scale data based on method requirements
            if self.method == 'nmf':
                # NMF requires non-negative values
                # Ensure non-negative first, then scale to [0,1]
                df = df.clip(lower=0)
                scaler = MinMaxScaler()
                data_scaled = scaler.fit_transform(df)
            else:
                # PCA and FA work with standardized data
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
        strategy: str,
        layer_name: str = None
    ) -> pd.DataFrame:
        """Impute missing values."""
        
        if not df.isnull().any().any():
            return df
        
        logger.info(f"Imputing with strategy: {strategy}")
        
        if strategy == 'network' and self.network is not None:
            return self._network_imputation_single_layer(df, layer_name)
        elif strategy == 'knn':
            imputer = KNNImputer(n_neighbors=5)
        elif strategy == 'iterative':
            imputer = IterativeImputer(random_state=self.random_state)
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
    
    def _network_imputation_single_layer(
        self,
        df: pd.DataFrame,
        layer_name: str
    ) -> pd.DataFrame:
        """
        Network-based imputation for a single layer.
        Uses connected features to impute missing values.
        """
        if self.network is None:
            raise ValueError("Network required for network-based imputation")
        
        logger.info(f"Performing network-based imputation for {layer_name}")
        
        import networkx as nx
        G = self.network.generate_graph()
        
        df_imputed = df.copy()
        
        for feature in df.columns:
            missing_idx = df[feature].isna()
            
            if not missing_idx.any():
                continue
            
            feature_str = str(feature)
            if feature_str not in G:
                # Fall back to mean imputation
                df_imputed.loc[missing_idx, feature] = df[feature].mean()
                continue
            
            # Get network neighbors
            try:
                neighbors = []
                for node in nx.single_source_shortest_path_length(G, feature_str, cutoff=2):
                    if node != feature_str:
                        neighbors.append(node)
                
                if not neighbors:
                    df_imputed.loc[missing_idx, feature] = df[feature].mean()
                    continue
                
                # Use neighbor data
                neighbor_data = []
                neighbor_weights = []
                
                for neighbor in neighbors[:10]:  # Limit to 10 neighbors
                    if neighbor in df.columns:
                        dist = nx.shortest_path_length(G, feature_str, neighbor)
                        weight = 1.0 / (dist ** 2)
                        
                        neighbor_data.append(df[neighbor])
                        neighbor_weights.append(weight)
                
                if not neighbor_data:
                    df_imputed.loc[missing_idx, feature] = df[feature].mean()
                    continue
                
                # Weighted average
                neighbor_matrix = np.column_stack(neighbor_data)
                weights_array = np.array(neighbor_weights)
                weights_array = weights_array / weights_array.sum()
                
                imputed_values = np.average(neighbor_matrix, axis=1, weights=weights_array)
                df_imputed.loc[missing_idx, feature] = imputed_values[missing_idx]
                
            except Exception as e:
                logger.warning(f"Network imputation failed for {feature}: {e}")
                df_imputed.loc[missing_idx, feature] = df[feature].mean()
        
        return df_imputed
    
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
    
    def _early_fusion(self, processed_data: Dict[str, pd.DataFrame]) -> pd.DataFrame:
        """Concatenate all features, then reduce."""
        
        logger.info("Performing early fusion")
        
        # Sort keys for consistency
        sorted_keys = sorted(processed_data.keys())
        
        # Concatenate all features
        concat_df = pd.concat(
            [processed_data[k] for k in sorted_keys],
            axis=1
        )
        
        logger.info(f"Combined shape: {concat_df.shape}")
        
        # Fit dimensionality reduction
        model = self._get_model(concat_df.shape[1])
        factors = model.fit_transform(concat_df.values)
        
        # Store model and loadings
        self.models_['integrated'] = model
        
        # Split loadings back to layers
        col_idx = 0
        for layer_name in sorted_keys:
            n_features = len(processed_data[layer_name].columns)
            
            if hasattr(model, 'components_'):
                layer_loadings = model.components_[:, col_idx:col_idx + n_features].T
                self.loadings_[layer_name] = pd.DataFrame(
                    layer_loadings,
                    index=processed_data[layer_name].columns,
                    columns=[f'Factor{i+1}' for i in range(self.n_components)]
                )
            
            col_idx += n_features
        
        # Store variance explained for PCA
        if self.method == 'pca':
            self.variance_explained_['integrated'] = model.explained_variance_ratio_
        
        return pd.DataFrame(
            factors,
            index=concat_df.index,
            columns=[f'Factor{i+1}' for i in range(factors.shape[1])]
        )
    
    def _late_fusion(self, processed_data: Dict[str, pd.DataFrame]) -> pd.DataFrame:
        """Reduce each omics separately, then concatenate factors."""
        
        logger.info("Performing late fusion")
        
        all_factors = []
        
        for layer_name, df in processed_data.items():
            logger.info(f"Reducing {layer_name}")
            
            model = self._get_model(df.shape[1])
            factors = model.fit_transform(df.values)
            
            self.models_[layer_name] = model
            
            # Store loadings
            if hasattr(model, 'components_'):
                self.loadings_[layer_name] = pd.DataFrame(
                    model.components_.T,
                    index=df.columns,
                    columns=[f'{layer_name}_Factor{i+1}' for i in range(factors.shape[1])]
                )
            
            # Store variance explained for PCA
            if self.method == 'pca':
                self.variance_explained_[layer_name] = model.explained_variance_ratio_
            
            # Create factor DataFrame
            factor_df = pd.DataFrame(
                factors,
                index=df.index,
                columns=[f'{layer_name}_Factor{i+1}' for i in range(factors.shape[1])]
            )
            all_factors.append(factor_df)
        
        # Concatenate all factors
        return pd.concat(all_factors, axis=1)
    
    def _concatenation(self, processed_data: Dict[str, pd.DataFrame]) -> pd.DataFrame:
        """Simply concatenate all features without reduction."""
        
        logger.info("Performing concatenation (no reduction)")
        
        sorted_keys = sorted(processed_data.keys())
        concat_df = pd.concat(
            [processed_data[k] for k in sorted_keys],
            axis=1
        )
        
        return concat_df
    
    def compute_cross_omics_correlation(
        self,
        omics1: str,
        omics2: str,
        method: str = 'pearson',
        min_samples: int = 3
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Compute pairwise correlations between factors of two omics layers.
        
        Parameters:
        -----------
        omics1, omics2 : str
            Names of omics layers
        method : str
            'pearson' or 'spearman'
        min_samples : int
            Minimum samples required
            
        Returns:
        --------
        Tuple[pd.DataFrame, pd.DataFrame]
            (correlation_matrix, pvalue_matrix)
        """
        
        if self.factors_ is None:
            raise ValueError("No factors computed. Run fit_transform() first.")
        
        # Get factor columns for each omics
        if self.integration_strategy == 'late':
            omics1_cols = [c for c in self.factors_.columns if c.startswith(f'{omics1}_')]
            omics2_cols = [c for c in self.factors_.columns if c.startswith(f'{omics2}_')]
        else:
            # For early fusion, use loadings
            if omics1 not in self.loadings_ or omics2 not in self.loadings_:
                raise ValueError(f"Loadings not available for {omics1} or {omics2}")
            
            omics1_cols = self.loadings_[omics1].columns.tolist()
            omics2_cols = self.loadings_[omics2].columns.tolist()
        
        # Compute correlations
        corr_matrix = np.zeros((len(omics1_cols), len(omics2_cols)))
        pval_matrix = np.zeros((len(omics1_cols), len(omics2_cols)))
        
        for i, col1 in enumerate(omics1_cols):
            for j, col2 in enumerate(omics2_cols):
                vals1 = self.factors_[col1].dropna()
                vals2 = self.factors_[col2].dropna()
                
                common_idx = vals1.index.intersection(vals2.index)
                
                if len(common_idx) < min_samples:
                    corr_matrix[i, j] = np.nan
                    pval_matrix[i, j] = np.nan
                    continue
                
                v1 = vals1.loc[common_idx]
                v2 = vals2.loc[common_idx]
                
                if method == 'pearson':
                    r, p = pearsonr(v1, v2)
                elif method == 'spearman':
                    r, p = spearmanr(v1, v2)
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
    
    # =============================================================================
    # Network-Aware Methods (TransNet's Key Differentiator)
    # =============================================================================
    
    def network_aware_correlation(
        self,
        data1: pd.DataFrame,
        data2: pd.DataFrame,
        layer1: str,
        layer2: str,
        alpha: float = 0.5,
        network_weight: str = 'uniform',
        correlation_method: str = 'pearson',
        min_samples: int = 3
    ) -> pd.DataFrame:
        """
        Correlate features between layers with network prior.
        
        TransNet's key innovation: combines statistical correlation with biological connection.
        
        Parameters:
        -----------
        data1, data2 : pd.DataFrame
            Data matrices (samples x features)
        layer1, layer2 : str
            Layer names ('Transcriptome', 'Proteome', 'Metabolome')
        alpha : float
            Balance between correlation (0) and network prior (1)
        network_weight : str
            How to weight network edges ('uniform', 'degree', 'centrality')
        correlation_method : str
            'pearson' or 'spearman'
        min_samples : int
            Minimum non-missing samples required
        
        Returns:
        --------
        pd.DataFrame
            Columns: feature1, feature2, correlation, p_value, network_prior,
                    combined_score, network_path
        """
        if self.network is None:
            raise ValueError("Network required for network-aware methods")
        
        logger.info(f"Computing network-aware correlations between {layer1} and {layer2}")
        
        # Align samples
        common_samples = set(data1.index).intersection(set(data2.index))
        if len(common_samples) < min_samples:
            raise ValueError(f"Need at least {min_samples} common samples")
        
        data1_aligned = data1.loc[list(common_samples)]
        data2_aligned = data2.loc[list(common_samples)]
        
        # Compute correlation matrix
        logger.info("Computing pairwise correlations")
        corr_results = self._compute_correlation_matrix(
            data1_aligned, 
            data2_aligned,
            method=correlation_method,
            min_samples=min_samples
        )
        
        # Build network prior matrix
        logger.info("Building network prior matrix")
        network_prior = self._build_network_prior_matrix(
            data1.columns.tolist(),
            data2.columns.tolist(),
            layer1,
            layer2,
            network_weight
        )
        
        # Combine correlation and network prior
        logger.info("Combining statistical and network evidence")
        results = []
        
        for _, row in corr_results.iterrows():
            feat1 = row['feature1']
            feat2 = row['feature2']
            corr = row['correlation']
            pval = row['p_value']
            
            # Get network prior
            prior = network_prior.loc[feat1, feat2] if (feat1 in network_prior.index and feat2 in network_prior.columns) else 0.0
            
            # Combined score: balance statistics and biology
            combined = (1 - alpha) * abs(corr) + alpha * prior
            
            # Find network path if connected
            path = None
            path_length = None
            if prior > 0:
                try:
                    paths = self.network.find_paths(feat1, feat2, max_length=3)
                    if paths:
                        path = ' -> '.join(paths[0])
                        path_length = len(paths[0]) - 1
                except Exception as e:
                    logger.debug(f"Could not find path between {feat1} and {feat2}: {e}")
            
            results.append({
                'feature1': feat1,
                'feature2': feat2,
                'layer1': layer1,
                'layer2': layer2,
                'correlation': corr,
                'p_value': pval,
                'network_prior': prior,
                'combined_score': combined,
                'network_path': path,
                'path_length': path_length
            })
        
        if not results:
            logger.warning("No correlation results found")
            return pd.DataFrame(columns=[
                'feature1', 'feature2', 'layer1', 'layer2', 
                'correlation', 'p_value', 'network_prior', 
                'combined_score', 'network_path', 'path_length'
            ])

        results_df = pd.DataFrame(results)
        if not results_df.empty:
            results_df = results_df.sort_values('combined_score', ascending=False)
        
        logger.info(f"Found {len(results_df)} feature pairs with combined evidence")
        logger.info(f"  - {(results_df['network_prior'] > 0).sum()} pairs with network connection")
        logger.info(f"  - {(results_df['p_value'] < 0.05).sum()} pairs with significant correlation")
        
        return results_df
    
    def _compute_correlation_matrix(
        self,
        data1: pd.DataFrame,
        data2: pd.DataFrame,
        method: str = 'pearson',
        min_samples: int = 3
    ) -> pd.DataFrame:
        """Compute pairwise correlations between features."""
        results = []
        
        for feat1 in data1.columns:
            for feat2 in data2.columns:
                vals1 = data1[feat1].dropna()
                vals2 = data2[feat2].dropna()
                
                common_idx = vals1.index.intersection(vals2.index)
                
                if len(common_idx) < min_samples:
                    continue
                
                vals1_aligned = vals1.loc[common_idx]
                vals2_aligned = vals2.loc[common_idx]
                
                try:
                    if method == 'pearson':
                        corr, pval = pearsonr(vals1_aligned, vals2_aligned)
                    elif method == 'spearman':
                        corr, pval = spearmanr(vals1_aligned, vals2_aligned)
                    else:
                        raise ValueError(f"Unknown correlation method: {method}")
                    
                    results.append({
                        'feature1': feat1,
                        'feature2': feat2,
                        'correlation': corr,
                        'p_value': pval,
                        'n_samples': len(common_idx)
                    })
                except Exception as e:
                    logger.debug(f"Could not compute correlation for {feat1}-{feat2}: {e}")
        
        return pd.DataFrame(results)
    
    def _build_network_prior_matrix(
        self,
        features1: List[str],
        features2: List[str],
        layer1: str,
        layer2: str,
        weight_type: str
    ) -> pd.DataFrame:
        """Build prior knowledge matrix from network."""
        import networkx as nx
        
        if not hasattr(self.network, '_edge_index_built') or not self.network._edge_index_built:
            self.network.build_cross_layer_edge_index()
        
        G = self.network.generate_graph()
        prior = pd.DataFrame(0.0, index=features1, columns=features2)
        
        # Pre-compute centralities if needed
        centrality = None
        if weight_type == 'centrality':
            logger.info("Computing betweenness centrality")
            centrality = nx.betweenness_centrality(G)
        
        for f1 in features1:
            f1_str = str(f1)
            
            if f1_str not in G:
                continue
            
            for f2 in features2:
                f2_str = str(f2)
                
                if f2_str not in G:
                    continue
                
                # Direct edge
                if G.has_edge(f1_str, f2_str) or G.has_edge(f2_str, f1_str):
                    if weight_type == 'uniform':
                        prior.loc[f1, f2] = 1.0
                    elif weight_type == 'degree':
                        prior.loc[f1, f2] = 1.0 / np.sqrt(G.degree(f1_str) * G.degree(f2_str))
                    elif weight_type == 'centrality':
                        prior.loc[f1, f2] = (centrality[f1_str] + centrality[f2_str]) / 2
                
                # Indirect connection
                elif nx.has_path(G, f1_str, f2_str):
                    try:
                        path_length = nx.shortest_path_length(G, f1_str, f2_str)
                        if path_length <= 3:
                            prior.loc[f1, f2] = 1.0 / (path_length ** 2)
                    except nx.NetworkXNoPath:
                        pass
        
        return prior
    
    def network_propagation(
        self,
        changed_features: Dict[str, float],
        layer: str,
        diffusion_alpha: float = 0.5,
        n_steps: int = 3,
        score_threshold: float = 0.01
    ) -> pd.DataFrame:
        """
        Propagate signal from changed features through network.
        
        Use case: If metabolites X,Y,Z changed significantly, which proteins
        are likely affected based on network topology?
        
        Parameters:
        -----------
        changed_features : dict
            {feature_id: fold_change} for significantly changed features
        layer : str
            Source layer
        diffusion_alpha : float
            Restart probability
        n_steps : int
            Number of diffusion steps
        score_threshold : float
            Minimum propagation score to report
        
        Returns:
        --------
        pd.DataFrame
            Columns: feature_id, layer, propagation_score, distance_to_source
        """
        if self.network is None:
            raise ValueError("Network required for network propagation")
        
        logger.info(f"Propagating signal from {len(changed_features)} changed features in {layer}")
        
        import networkx as nx
        G = self.network.generate_graph()
        
        changed_features = {str(k): v for k, v in changed_features.items()}
        
        # Initialize signal
        signal = {str(node): 0.0 for node in G.nodes()}
        for feature, fc in changed_features.items():
            if feature in signal:
                signal[feature] = abs(fc)
        
        # Normalize
        total_signal = sum(signal.values())
        if total_signal > 0:
            signal = {k: v / total_signal for k, v in signal.items()}
        else:
            logger.warning("No signal to propagate")
            return pd.DataFrame()
        
        # Diffusion
        logger.info(f"Running diffusion for {n_steps} steps")
        for step in range(n_steps):
            new_signal = {}
            
            for node in G.nodes():
                node_str = str(node)
                
                restart = changed_features.get(node_str, 0.0)
                if restart > 0:
                    restart = abs(restart) / total_signal if total_signal > 0 else 0.0
                
                neighbors = list(G.neighbors(node_str))
                neighbor_avg = np.mean([signal[str(n)] for n in neighbors]) if neighbors else 0.0
                
                new_signal[node_str] = (
                    diffusion_alpha * restart +
                    (1 - diffusion_alpha) * neighbor_avg
                )
            
            signal = new_signal
        
        # Calculate distances
        logger.info("Computing distances to source features")
        results = []
        
        for node in G.nodes():
            node_str = str(node)
            
            if signal[node_str] > score_threshold:
                min_dist = float('inf')
                for source in changed_features.keys():
                    if source in G:
                        try:
                            dist = nx.shortest_path_length(G, source, node_str)
                            min_dist = min(min_dist, dist)
                        except nx.NetworkXNoPath:
                            pass
                
                node_layer = G.nodes[node_str].get('layer', 'Unknown')
                
                results.append({
                    'feature_id': node_str,
                    'layer': node_layer,
                    'propagation_score': signal[node_str],
                    'distance_to_source': min_dist if min_dist != float('inf') else None,
                    'is_source': node_str in changed_features
                })
        
        results_df = pd.DataFrame(results).sort_values('propagation_score', ascending=False)
        
        logger.info(f"Found {len(results_df)} features with propagation score > {score_threshold}")
        
        return results_df
    
    def canonical_correlation_analysis(
        self,
        data1: pd.DataFrame,
        data2: pd.DataFrame,
        n_components: int = 5
    ) -> Dict[str, Union[pd.DataFrame, List[float]]]:
        """
        Canonical Correlation Analysis between two omics layers.
        
        Parameters:
        -----------
        data1, data2 : pd.DataFrame
            Data matrices (samples x features)
        n_components : int
            Number of canonical components
        
        Returns:
        --------
        Dict
            Dictionary with transformed data and correlations
        """
        logger.info(f"Running CCA with {n_components} components")
        
        # Align samples
        common_samples = set(data1.index).intersection(set(data2.index))
        
        if len(common_samples) == 0:
            raise ValueError("No common samples between data matrices")
        
        data1_aligned = data1.loc[list(common_samples)].fillna(0)
        data2_aligned = data2.loc[list(common_samples)].fillna(0)
        
        # Standardize
        scaler1 = StandardScaler()
        scaler2 = StandardScaler()
        
        data1_scaled = scaler1.fit_transform(data1_aligned)
        data2_scaled = scaler2.fit_transform(data2_aligned)
        
        # Run CCA
        model = CCA(n_components=n_components)
        X_c, Y_c = model.fit_transform(data1_scaled, data2_scaled)
        
        # Calculate canonical correlations
        correlations = [pearsonr(X_c[:, i], Y_c[:, i])[0] for i in range(n_components)]
        
        results = {
            'X_canonical': pd.DataFrame(
                X_c,
                index=data1_aligned.index,
                columns=[f'CC{i+1}_X' for i in range(n_components)]
            ),
            'Y_canonical': pd.DataFrame(
                Y_c,
                index=data2_aligned.index,
                columns=[f'CC{i+1}_Y' for i in range(n_components)]
            ),
            'correlations': correlations,
            'model': model
        }
        
        logger.info(f"CCA completed. Canonical correlations: {correlations}")
        
        return results


# =============================================================================
# Utility Functions
# =============================================================================

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
        integrator = MultiOmicsIntegrator(network=transnet_obj, **integrator_kwargs)
    
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
    
    Parameters:
    -----------
    G : nx.Graph
        NetworkX graph with pathway annotations
    factor_scores : pd.Series
        Factor loadings for nodes
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
    try:
        from statsmodels.stats.multitest import multipletests
        if len(df) > 0:
            _, df['fdr'], _, _ = multipletests(df['pvalue'], method='fdr_bh')
    except ImportError:
        logger.warning("statsmodels not available for FDR correction")
    
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
        Mapping of pathway names to feature lists
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