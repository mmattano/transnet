"""
Tests for multi-omics integration module.
"""

import pytest
import numpy as np
import pandas as pd
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from transnet.analysis.multi_omics_integration import MultiOmicsIntegrator


@pytest.fixture
def synthetic_omics_data():
    """Generate synthetic multi-omics data for testing."""
    np.random.seed(42)
    n_samples = 30
    n_proteins = 50
    n_metabolites = 40
    
    # Generate with positive values for NMF compatibility
    proteomics = pd.DataFrame(
        np.random.rand(n_samples, n_proteins) * 10,
        columns=[f'P{i}' for i in range(n_proteins)],
        index=[f'S{i}' for i in range(n_samples)]
    )
    
    metabolomics = pd.DataFrame(
        np.random.rand(n_samples, n_metabolites) * 10,
        columns=[f'M{i}' for i in range(n_metabolites)],
        index=[f'S{i}' for i in range(n_samples)]
    )
    
    return {
        'proteomics': proteomics,
        'metabolomics': metabolomics
    }


@pytest.fixture
def synthetic_omics_data_with_missing():
    """Generate synthetic data with missing values."""
    np.random.seed(42)
    n_samples = 30
    n_proteins = 50
    n_metabolites = 40
    
    proteomics = pd.DataFrame(
        np.random.rand(n_samples, n_proteins) * 10,
        columns=[f'P{i}' for i in range(n_proteins)],
        index=[f'S{i}' for i in range(n_samples)]
    )
    
    metabolomics = pd.DataFrame(
        np.random.rand(n_samples, n_metabolites) * 10,
        columns=[f'M{i}' for i in range(n_metabolites)],
        index=[f'S{i}' for i in range(n_samples)]
    )
    
    # Introduce missing values
    proteomics.iloc[0:5, 0:10] = np.nan
    metabolomics.iloc[10:15, 5:15] = np.nan
    
    return {
        'proteomics': proteomics,
        'metabolomics': metabolomics
    }


# Test Early Fusion
def test_early_fusion_pca(synthetic_omics_data):
    """Test early fusion with PCA."""
    integrator = MultiOmicsIntegrator(
        n_components=5,
        method='pca',
        integration_strategy='early',
        random_state=42
    )
    
    factors = integrator.fit_transform(synthetic_omics_data, imputation_strategy='none')
    
    assert factors.shape == (30, 5)
    assert integrator.factors_ is not None
    assert 'proteomics' in integrator.feature_names_
    assert 'metabolomics' in integrator.feature_names_


def test_early_fusion_fa(synthetic_omics_data):
    """Test early fusion with Factor Analysis."""
    integrator = MultiOmicsIntegrator(
        n_components=5,
        method='fa',
        integration_strategy='early',
        random_state=42
    )
    
    factors = integrator.fit_transform(synthetic_omics_data, imputation_strategy='none')
    
    assert factors.shape == (30, 5)
    assert integrator.factors_ is not None


def test_early_fusion_nmf(synthetic_omics_data):
    """Test early fusion with NMF."""
    integrator = MultiOmicsIntegrator(
        n_components=5,
        method='nmf',
        integration_strategy='early',
        random_state=42
    )
    
    factors = integrator.fit_transform(synthetic_omics_data, imputation_strategy='none')
    
    assert factors.shape == (30, 5)
    assert integrator.factors_ is not None
    # NMF should produce non-negative factors
    assert (factors >= 0).all().all()


# Test Late Fusion
def test_late_fusion_pca(synthetic_omics_data):
    """Test late fusion with PCA."""
    integrator = MultiOmicsIntegrator(
        n_components=5,
        method='pca',
        integration_strategy='late',
        random_state=42
    )
    
    factors = integrator.fit_transform(synthetic_omics_data, imputation_strategy='none')
    
    # Late fusion creates separate factors per layer
    assert factors.shape[0] == 30
    assert factors.shape[1] == 10  # 5 factors x 2 layers


def test_late_fusion_nmf(synthetic_omics_data):
    """Test late fusion with NMF."""
    integrator = MultiOmicsIntegrator(
        n_components=5,
        method='nmf',
        integration_strategy='late',
        random_state=42
    )
    
    factors = integrator.fit_transform(synthetic_omics_data, imputation_strategy='none')
    
    assert factors.shape[0] == 30
    assert factors.shape[1] == 10
    # All factors should be non-negative
    assert (factors >= 0).all().all()


# Test Concatenation
def test_concatenation(synthetic_omics_data):
    """Test concatenation strategy."""
    integrator = MultiOmicsIntegrator(
        n_components=5,
        method='pca',
        integration_strategy='concatenation',
        random_state=42
    )
    
    factors = integrator.fit_transform(synthetic_omics_data, imputation_strategy='none')
    
    # Concatenation should combine all features
    assert factors.shape == (30, 90)  # 50 + 40 features


# Test Imputation Strategies
def test_knn_imputation(synthetic_omics_data_with_missing):
    """Test KNN imputation."""
    integrator = MultiOmicsIntegrator(
        n_components=5,
        method='pca',
        integration_strategy='early',
        random_state=42
    )
    
    factors = integrator.fit_transform(
        synthetic_omics_data_with_missing,
        imputation_strategy='knn'
    )
    
    assert factors.shape == (30, 5)
    assert not factors.isna().any().any()


def test_mean_imputation(synthetic_omics_data_with_missing):
    """Test mean imputation."""
    integrator = MultiOmicsIntegrator(
        n_components=5,
        method='pca',
        integration_strategy='early',
        random_state=42
    )
    
    factors = integrator.fit_transform(
        synthetic_omics_data_with_missing,
        imputation_strategy='mean'
    )
    
    assert factors.shape == (30, 5)
    assert not factors.isna().any().any()


# Test Error Handling
def test_misaligned_samples():
    """Test error handling for misaligned samples."""
    proteomics = pd.DataFrame(
        np.random.rand(30, 50),
        columns=[f'P{i}' for i in range(50)],
        index=[f'S{i}' for i in range(30)]
    )
    
    metabolomics = pd.DataFrame(
        np.random.rand(25, 40),  # Different number of samples
        columns=[f'M{i}' for i in range(40)],
        index=[f'S{i}' for i in range(25)]
    )
    
    omics_data = {
        'proteomics': proteomics,
        'metabolomics': metabolomics
    }
    
    integrator = MultiOmicsIntegrator(
        n_components=5,
        method='pca',
        integration_strategy='early'
    )
    
    with pytest.raises(ValueError, match="Sample names not aligned"):
        integrator.fit_transform(omics_data)


def test_empty_data():
    """Test error handling for empty data."""
    integrator = MultiOmicsIntegrator(
        n_components=5,
        method='pca',
        integration_strategy='early'
    )
    
    with pytest.raises(ValueError, match="cannot be empty"):
        integrator.fit_transform({})


def test_invalid_method():
    """Test error handling for invalid method."""
    integrator = MultiOmicsIntegrator(
        n_components=5,
        method='invalid_method',
        integration_strategy='early'
    )
    
    omics_data = {
        'proteomics': pd.DataFrame(np.random.rand(30, 50)),
        'metabolomics': pd.DataFrame(np.random.rand(30, 40))
    }
    
    with pytest.raises(ValueError, match="Unknown method"):
        integrator.fit_transform(omics_data)


def test_invalid_strategy():
    """Test error handling for invalid strategy."""
    integrator = MultiOmicsIntegrator(
        n_components=5,
        method='pca',
        integration_strategy='invalid_strategy'
    )
    
    omics_data = {
        'proteomics': pd.DataFrame(np.random.rand(30, 50)),
        'metabolomics': pd.DataFrame(np.random.rand(30, 40))
    }
    
    with pytest.raises(ValueError, match="Unknown strategy"):
        integrator.fit_transform(omics_data)


# Test Transform on New Data
def test_transform_new_data(synthetic_omics_data):
    """Test transforming new data with fitted model."""
    integrator = MultiOmicsIntegrator(
        n_components=5,
        method='pca',
        integration_strategy='early',
        random_state=42
    )
    
    # Fit on training data
    integrator.fit_transform(synthetic_omics_data, imputation_strategy='none')
    
    # Generate new test data with same features
    test_proteomics = pd.DataFrame(
        np.random.rand(10, 50) * 10,
        columns=[f'P{i}' for i in range(50)],
        index=[f'T{i}' for i in range(10)]
    )
    
    test_metabolomics = pd.DataFrame(
        np.random.rand(10, 40) * 10,
        columns=[f'M{i}' for i in range(40)],
        index=[f'T{i}' for i in range(10)]
    )
    
    test_data = {
        'proteomics': test_proteomics,
        'metabolomics': test_metabolomics
    }
    
    # Transform new data
    test_factors = integrator.transform(test_data)
    
    assert test_factors.shape == (10, 5)


# Test Loadings
def test_get_loadings(synthetic_omics_data):
    """Test retrieving factor loadings."""
    integrator = MultiOmicsIntegrator(
        n_components=5,
        method='pca',
        integration_strategy='early',
        random_state=42
    )
    
    integrator.fit_transform(synthetic_omics_data, imputation_strategy='none')
    
    assert 'proteomics' in integrator.loadings_
    assert 'metabolomics' in integrator.loadings_
    
    # Check loadings shape
    prot_loadings = integrator.loadings_['proteomics']
    assert prot_loadings.shape[0] == 50  # Number of features
    assert prot_loadings.shape[1] == 5   # Number of factors


def test_get_top_features(synthetic_omics_data):
    """Test getting top contributing features."""
    integrator = MultiOmicsIntegrator(
        n_components=5,
        method='pca',
        integration_strategy='early',
        random_state=42
    )
    
    integrator.fit_transform(synthetic_omics_data, imputation_strategy='none')
    
    top_features = integrator.get_top_features(
        layer_name='proteomics',
        factor_idx=0,
        n_top=10,
        return_weights=True
    )
    
    assert isinstance(top_features, pd.DataFrame)
    assert len(top_features) == 10
    # Check for actual column names (capitalized)
    assert 'Feature' in top_features.columns or 'feature' in top_features.columns
    assert 'Loading' in top_features.columns or 'weight' in top_features.columns


# Test Cross-Omics Correlation
def test_cross_omics_correlation(synthetic_omics_data):
    """Test cross-omics factor correlation analysis."""
    integrator = MultiOmicsIntegrator(
        n_components=5,
        method='pca',
        integration_strategy='late',
        random_state=42
    )
    
    integrator.fit_transform(synthetic_omics_data, imputation_strategy='none')
    
    corr_matrix, pval_matrix = integrator.compute_cross_omics_correlation(
        'proteomics',
        'metabolomics',
        method='pearson'
    )
    
    assert isinstance(corr_matrix, pd.DataFrame)
    assert isinstance(pval_matrix, pd.DataFrame)
    assert corr_matrix.shape == pval_matrix.shape
    assert corr_matrix.shape == (5, 5)  # 5 factors per layer


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
