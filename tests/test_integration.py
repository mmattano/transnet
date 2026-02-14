#!/usr/bin/env python3
"""
TransNet Multi-Omics Integration Demo
Uses uploaded transcriptomics, proteomics, metabolomics data
"""

import pandas as pd
import numpy as np
from transnet.analysis.multi_omics_integration import MultiOmicsIntegrator
import matplotlib.pyplot as plt
import seaborn as sns

# Load data
transcriptomics = pd.read_csv('/Users/krv114/Desktop/Data/Multi/transcriptomics_log2.csv', index_col=0)
proteomics = pd.read_csv('/Users/krv114/Desktop/Data/Multi/proteomics_log2.csv', index_col=0)
metabolomics = pd.read_csv('/Users/krv114/Desktop/Data/Multi/metabolomics_filtered_log2.csv', index_col=0)

# Transpose: samples as rows, features as columns
transcriptomics = transcriptomics.T
proteomics = proteomics.T
metabolomics = metabolomics.T

print(f"Transcriptomics: {transcriptomics.shape} (samples × genes)")
print(f"Proteomics: {proteomics.shape} (samples × proteins)")
print(f"Metabolomics: {metabolomics.shape} (samples × metabolites)")

# Create omics dictionary
omics_data = {
    'transcriptomics': transcriptomics,
    'proteomics': proteomics,
    'metabolomics': metabolomics
}

# Initialize integrator
integrator = MultiOmicsIntegrator(
    n_components=5,
    method='nmf',  # or 'pca', 'fa'
    integration_strategy='early',  # or 'late'
    random_state=42
)

# Fit and transform
factors = integrator.fit_transform(omics_data, imputation_strategy='knn')
print(f"\nIntegrated factors: {factors.shape}")
print(factors.head())

# Access loadings for each layer
for layer in ['transcriptomics', 'proteomics', 'metabolomics']:
    print(f"\n{layer.capitalize()} loadings: {integrator.loadings_[layer].shape}")
    
    # Get top features for Factor 1
    top_features = integrator.get_top_features(layer, factor_idx=0, n_top=10)
    print(f"Top 10 features in Factor 1:\n{top_features}")

# Visualize factors
fig, ax = plt.subplots(figsize=(10, 6))
scatter = ax.scatter(factors.iloc[:, 0], factors.iloc[:, 1], 
                     c=range(len(factors)), cmap='viridis', s=100)
ax.set_xlabel('Factor 1')
ax.set_ylabel('Factor 2')
ax.set_title('Multi-Omics Integration: Sample Projection')

# Add sample labels
for idx, label in enumerate(factors.index):
    ax.annotate(label, (factors.iloc[idx, 0], factors.iloc[idx, 1]), 
                fontsize=8, alpha=0.7)

plt.colorbar(scatter, label='Sample Order')
plt.tight_layout()
plt.savefig('results/factor_projection.png', dpi=300)
print("\nSaved factor projection to results/factor_projection.png")

# Cross-layer correlation analysis
transcriptomics_factors = integrator.factors_by_layer_['transcriptomics']
proteomics_factors = integrator.factors_by_layer_['proteomics']

# Compute factor correlations
factor_corr = pd.DataFrame(
    np.corrcoef(transcriptomics_factors.T, proteomics_factors.T)[:5, 5:],
    index=[f'RNA_F{i+1}' for i in range(5)],
    columns=[f'Protein_F{i+1}' for i in range(5)]
)

print("\nCross-layer factor correlations:")
print(factor_corr)

# Heatmap
plt.figure(figsize=(8, 6))
sns.heatmap(factor_corr, annot=True, fmt='.2f', cmap='coolwarm', center=0)
plt.title('Transcriptome-Proteome Factor Correlations')
plt.tight_layout()
plt.savefig('results/cross_layer_correlation.png', dpi=300)
print("Saved correlation heatmap to results/cross_layer_correlation.png")

# Compare integration strategies
strategies = ['early', 'late']
methods = ['nmf', 'pca']

results = {}
for strategy in strategies:
    for method in methods:
        key = f"{strategy}_{method}"
        integrator_test = MultiOmicsIntegrator(
            n_components=5,
            method=method,
            integration_strategy=strategy,
            random_state=42
        )
        factors_test = integrator_test.fit_transform(omics_data, imputation_strategy='knn')
        results[key] = factors_test

# Plot comparison
fig, axes = plt.subplots(2, 2, figsize=(14, 12))
for idx, (key, factors_test) in enumerate(results.items()):
    ax = axes[idx // 2, idx % 2]
    ax.scatter(factors_test.iloc[:, 0], factors_test.iloc[:, 1], alpha=0.6, s=100)
    ax.set_xlabel('Factor 1')
    ax.set_ylabel('Factor 2')
    ax.set_title(key.replace('_', ' + ').upper())
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('results/strategy_comparison.png', dpi=300)
print("Saved strategy comparison to results/strategy_comparison.png")

print("\n=== Integration Complete ===")