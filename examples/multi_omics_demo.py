"""
Demonstration of multi-omics integration with TransNet.

This example shows how to:
1. Load multi-omics data (simulated TCGA-style)
2. Integrate using different strategies
3. Map factors to biological networks
4. Perform network-aware enrichment analysis
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from transnet.analysis.multi_omics_integration import (
    MultiOmicsIntegrator,
    integrate_with_transnet,
    compute_factor_network_enrichment
)
from transnet.biology.transnet import Transnet
from transnet.biology.layers import Proteome, Metabolome

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def generate_synthetic_omics_data(
    n_samples: int = 50,
    n_proteins: int = 200,
    n_metabolites: int = 150,
    n_factors: int = 5,
    noise_level: float = 0.3
):
    """
    Generate synthetic multi-omics data with shared latent structure.
    
    Simulates TCGA-style data with:
    - Proteomics: log2(intensity)
    - Metabolomics: log2(concentration)
    - Shared latent factors driving both layers
    """
    
    logger.info("Generating synthetic multi-omics data")
    
    np.random.seed(42)
    
    # Generate shared latent factors
    factors = np.random.randn(n_samples, n_factors)
    
    # Generate protein loadings
    protein_loadings = np.random.randn(n_proteins, n_factors)
    protein_loadings = protein_loadings / np.linalg.norm(protein_loadings, axis=0)
    
    # Generate metabolite loadings
    metabolite_loadings = np.random.randn(n_metabolites, n_factors)
    metabolite_loadings = metabolite_loadings / np.linalg.norm(metabolite_loadings, axis=0)
    
    # Generate data
    proteomics = factors @ protein_loadings.T
    metabolomics = factors @ metabolite_loadings.T
    
    # Add noise
    proteomics += np.random.randn(*proteomics.shape) * noise_level
    metabolomics += np.random.randn(*metabolomics.shape) * noise_level
    
    # Add missing values (10%)
    proteomics_mask = np.random.rand(*proteomics.shape) > 0.1
    proteomics = np.where(proteomics_mask, proteomics, np.nan)
    
    metabolomics_mask = np.random.rand(*metabolomics.shape) > 0.1
    metabolomics = np.where(metabolomics_mask, metabolomics, np.nan)
    
    # Convert to DataFrames
    sample_names = [f"Sample_{i+1}" for i in range(n_samples)]
    protein_names = [f"Protein_{i+1}" for i in range(n_proteins)]
    metabolite_names = [f"Metabolite_{i+1}" for i in range(n_metabolites)]
    
    proteomics_df = pd.DataFrame(
        proteomics,
        index=sample_names,
        columns=protein_names
    )
    
    metabolomics_df = pd.DataFrame(
        metabolomics,
        index=sample_names,
        columns=metabolite_names
    )
    
    return {
        'proteomics': proteomics_df,
        'metabolomics': metabolomics_df,
        'true_factors': factors
    }


def compare_integration_strategies(omics_data, output_dir):
    """Compare different integration strategies."""
    
    logger.info("Comparing integration strategies")
    
    strategies = ['early', 'late']
    methods = ['nmf', 'pca', 'fa']
    
    fig, axes = plt.subplots(len(strategies), len(methods), figsize=(15, 10))
    
    for i, strategy in enumerate(strategies):
        for j, method in enumerate(methods):
            logger.info(f"Testing {strategy} fusion with {method}")
            
            integrator = MultiOmicsIntegrator(
                n_components=5,
                method=method,
                integration_strategy=strategy,
                random_state=42
            )
            
            try:
                factors = integrator.fit_transform(omics_data)
                
                # Plot first two factors
                ax = axes[i, j]
                ax.scatter(factors.iloc[:, 0], factors.iloc[:, 1], alpha=0.6)
                ax.set_xlabel('Factor 1')
                ax.set_ylabel('Factor 2')
                ax.set_title(f'{strategy.capitalize()} + {method.upper()}')
                ax.grid(True, alpha=0.3)
                
            except Exception as e:
                logger.error(f"Error with {strategy}/{method}: {e}")
                axes[i, j].text(0.5, 0.5, f'Failed\n{str(e)[:30]}', 
                               ha='center', va='center', transform=axes[i, j].transAxes)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'integration_comparison.png'), dpi=300)
    logger.info(f"Saved comparison plot to {output_dir}")
    plt.close()


def analyze_factor_loadings(integrator, omics_data, output_dir):
    """Analyze and visualize factor loadings."""
    
    logger.info("Analyzing factor loadings")
    
    # Get top features for each factor
    for layer_name in omics_data.keys():
        logger.info(f"Top features for {layer_name}")
        
        # Create heatmap of top features across factors
        n_factors = len(integrator.loadings_[layer_name].columns)
        top_features_per_factor = []
        
        for factor_idx in range(min(3, n_factors)):  # First 3 factors
            top_df = integrator.get_top_features(
                layer_name, factor_idx, n_top=20, return_weights=True
            )
            print(f"\nFactor {factor_idx + 1} top features:")
            print(top_df.head(10))
            
            top_features_per_factor.extend(top_df['Feature'].tolist())
        
        # Plot heatmap
        if top_features_per_factor:
            top_features_unique = list(dict.fromkeys(top_features_per_factor))[:30]
            
            heatmap_data = integrator.loadings_[layer_name].loc[
                top_features_unique,
                integrator.loadings_[layer_name].columns[:min(5, n_factors)]
            ]
            
            plt.figure(figsize=(10, 12))
            sns.heatmap(
                heatmap_data,
                cmap='RdBu_r',
                center=0,
                cbar_kws={'label': 'Loading'},
                yticklabels=True
            )
            plt.title(f'{layer_name.capitalize()} Factor Loadings')
            plt.xlabel('Factor')
            plt.ylabel('Feature')
            plt.tight_layout()
            plt.savefig(
                os.path.join(output_dir, f'{layer_name}_loadings_heatmap.png'),
                dpi=300,
                bbox_inches='tight'
            )
            plt.close()


def demonstrate_cross_omics_correlation(integrator, output_dir):
    """Demonstrate cross-omics factor correlation analysis."""
    
    if integrator.integration_strategy != 'late':
        logger.info("Skipping cross-omics correlation (not late fusion)")
        return
    
    logger.info("Computing cross-omics correlations")
    
    try:
        corr_matrix, pval_matrix = integrator.compute_cross_omics_correlations(
            'proteomics', 'metabolomics', method='pearson'
        )
        
        print("\nCross-omics factor correlations:")
        print(corr_matrix)
        
        # Plot correlation heatmap
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        # Correlation values
        sns.heatmap(
            corr_matrix,
            annot=True,
            fmt='.2f',
            cmap='RdBu_r',
            center=0,
            vmin=-1,
            vmax=1,
            ax=axes[0],
            cbar_kws={'label': 'Correlation'}
        )
        axes[0].set_title('Proteomics-Metabolomics Factor Correlations')
        axes[0].set_xlabel('Metabolomics Factors')
        axes[0].set_ylabel('Proteomics Factors')
        
        # Significance (log p-values)
        log_pval = -np.log10(pval_matrix + 1e-300)
        sns.heatmap(
            log_pval,
            annot=True,
            fmt='.1f',
            cmap='YlOrRd',
            ax=axes[1],
            cbar_kws={'label': '-log10(p-value)'}
        )
        axes[1].set_title('Statistical Significance')
        axes[1].set_xlabel('Metabolomics Factors')
        axes[1].set_ylabel('Proteomics Factors')
        
        plt.tight_layout()
        plt.savefig(
            os.path.join(output_dir, 'cross_omics_correlations.png'),
            dpi=300,
            bbox_inches='tight'
        )
        plt.close()
        
    except Exception as e:
        logger.error(f"Error computing cross-omics correlations: {e}")


def demonstrate_network_integration(omics_data, output_dir):
    """Demonstrate integration with biological networks."""
    
    logger.info("Demonstrating network integration")
    
    # Create a simplified network
    proteome = Proteome()
    proteome.ncbi_organism = "10090"  # Mouse
    proteome.kegg_organism = "mmu"
    
    # For demo, just create a few fake proteins matching our data
    from transnet.biology.elements import Protein
    for protein_name in omics_data['proteomics'].columns[:20]:
        protein = Protein(
            uniprot_id=protein_name,
            name=protein_name,
            ncbi_organism="10090"
        )
        proteome.proteins.append(protein)
    
    metabolome = Metabolome()
    from transnet.biology.elements import Metabolite
    for met_name in omics_data['metabolomics'].columns[:20]:
        metabolite = Metabolite(
            kegg_compound_id=met_name,
            kegg_name=met_name
        )
        metabolome.metabolites.append(metabolite)
    
    # Create transnet object
    transnet = Transnet(
        name="demo_network",
        proteome=proteome,
        metabolome=metabolome
    )
    
    # Integrate omics data with network
    logger.info("Integrating multi-omics data with network")
    
    factors, integrator = integrate_with_transnet(
        transnet,
        omics_data,
        n_components=5,
        method='nmf',
        integration_strategy='early'
    )
    
    logger.info(f"Generated factors shape: {factors.shape}")
    logger.info(f"Network has {transnet.graph.number_of_nodes() if transnet.graph else 0} nodes")
    
    # Save factors
    factors.to_csv(os.path.join(output_dir, 'integrated_factors.csv'))
    logger.info("Saved integrated factors")
    
    return factors, integrator, transnet


def main():
    """Run the complete multi-omics integration demo."""
    
    # Create output directory
    output_dir = "multi_omics_output"
    os.makedirs(output_dir, exist_ok=True)
    
    logger.info("=" * 60)
    logger.info("TransNet Multi-Omics Integration Demo")
    logger.info("=" * 60)
    
    # Generate synthetic data
    data = generate_synthetic_omics_data(
        n_samples=50,
        n_proteins=200,
        n_metabolites=150,
        n_factors=5
    )
    
    omics_data = {
        'proteomics': data['proteomics'],
        'metabolomics': data['metabolomics']
    }
    
    logger.info(f"Generated data shapes:")
    for name, df in omics_data.items():
        logger.info(f"  {name}: {df.shape}")
    
    # Compare integration strategies
    compare_integration_strategies(omics_data, output_dir)
    
    # Detailed analysis with one strategy
    logger.info("\nDetailed analysis with early fusion + NMF")
    integrator = MultiOmicsIntegrator(
        n_components=5,
        method='nmf',
        integration_strategy='early',
        random_state=42
    )
    
    factors = integrator.fit_transform(omics_data)
    logger.info(f"Generated {factors.shape[1]} factors for {factors.shape[0]} samples")
    
    # Analyze loadings
    analyze_factor_loadings(integrator, omics_data, output_dir)
    
    # Test late fusion for cross-omics correlation
    logger.info("\nTesting late fusion for cross-omics analysis")
    integrator_late = MultiOmicsIntegrator(
        n_components=5,
        method='pca',
        integration_strategy='late',
        random_state=42
    )
    
    factors_late = integrator_late.fit_transform(omics_data)
    demonstrate_cross_omics_correlation(integrator_late, output_dir)
    
    # Network integration
    factors_net, integrator_net, transnet = demonstrate_network_integration(
        omics_data, output_dir
    )
    
    logger.info("\n" + "=" * 60)
    logger.info(f"Demo complete! Results saved to {output_dir}/")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()