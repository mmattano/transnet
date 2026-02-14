#!/usr/bin/env python3
"""
TransNet Web Interface - Multi-Omics Integration
Streamlit app for interactive analysis
"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from io import BytesIO
import sys
sys.path.insert(0, '/path/to/transnet')  # Adjust path

from transnet.analysis.multi_omics_integration import MultiOmicsIntegrator

st.set_page_config(page_title="TransNet Multi-Omics", layout="wide")

st.title("🧬 TransNet: Multi-Omics Integration")
st.markdown("Statistical integration of transcriptomics, proteomics, and metabolomics data")

# Sidebar configuration
st.sidebar.header("Configuration")

method = st.sidebar.selectbox(
    "Integration Method",
    options=['nmf', 'pca', 'fa'],
    help="NMF: Non-negative Matrix Factorization\nPCA: Principal Component Analysis\nFA: Factor Analysis"
)

strategy = st.sidebar.selectbox(
    "Integration Strategy",
    options=['early', 'late'],
    help="Early: concatenate features then reduce\nLate: reduce each layer separately"
)

n_components = st.sidebar.slider("Number of Components", 2, 20, 5)

imputation = st.sidebar.selectbox(
    "Missing Data Imputation",
    options=['knn', 'iterative', 'mean', 'median'],
    help="Strategy for handling missing values"
)

# File upload
st.header("1. Upload Data")

col1, col2, col3 = st.columns(3)

with col1:
    st.subheader("Transcriptomics")
    rna_file = st.file_uploader(
        "Upload log2 RNA-seq data (CSV)",
        type=['csv'],
        key='rna',
        help="Features as rows, samples as columns"
    )
    
with col2:
    st.subheader("Proteomics")
    protein_file = st.file_uploader(
        "Upload log2 protein data (CSV)",
        type=['csv'],
        key='protein',
        help="Features as rows, samples as columns"
    )
    
with col3:
    st.subheader("Metabolomics")
    metabolite_file = st.file_uploader(
        "Upload log2 metabolite data (CSV)",
        type=['csv'],
        key='metabolite',
        help="Features as rows, samples as columns"
    )

# Load example data button
if st.sidebar.button("Load Example Data"):
    st.session_state['example_data'] = True

# Process uploaded or example data
omics_data = {}
layer_names = []

if 'example_data' in st.session_state or (rna_file and protein_file and metabolite_file):
    
    if 'example_data' in st.session_state:
        # Load example data from mounted files
        try:
            rna_df = pd.read_csv('/mnt/user-data/uploads/transcriptomics_log2.csv', index_col=0).T
            protein_df = pd.read_csv('/mnt/user-data/uploads/proteomics_log2.csv', index_col=0).T
            metabolite_df = pd.read_csv('/mnt/user-data/uploads/metabolomics_filtered_log2.csv', index_col=0).T
            st.info("✓ Loaded example data (mouse time-course)")
        except:
            st.error("Example data files not found")
            st.stop()
    else:
        # Load uploaded files
        rna_df = pd.read_csv(rna_file, index_col=0).T
        protein_df = pd.read_csv(protein_file, index_col=0).T
        metabolite_df = pd.read_csv(metabolite_file, index_col=0).T
    
    # Build omics dictionary
    omics_data = {
        'transcriptomics': rna_df,
        'proteomics': protein_df,
        'metabolomics': metabolite_df
    }
    layer_names = ['transcriptomics', 'proteomics', 'metabolomics']
    
    # Display data summary
    st.header("2. Data Summary")
    
    summary_data = []
    for layer, df in omics_data.items():
        summary_data.append({
            'Layer': layer.capitalize(),
            'Features': df.shape[1],
            'Samples': df.shape[0],
            'Missing (%)': f"{(df.isna().sum().sum() / df.size * 100):.1f}"
        })
    
    st.dataframe(pd.DataFrame(summary_data), use_container_width=True)
    
    # Run integration
    if st.button("🚀 Run Integration", type="primary", use_container_width=True):
        
        with st.spinner("Integrating omics layers..."):
            
            try:
                # Initialize integrator
                integrator = MultiOmicsIntegrator(
                    n_components=n_components,
                    method=method,
                    integration_strategy=strategy,
                    imputation_strategy=imputation,
                    random_state=42
                )
                
                # Fit transform
                factors = integrator.fit_transform(omics_data)
                
                st.success(f"✓ Integration complete: {factors.shape[0]} samples × {factors.shape[1]} factors")
                
                # Store in session state
                st.session_state['factors'] = factors
                st.session_state['integrator'] = integrator
                st.session_state['omics_data'] = omics_data
                
            except Exception as e:
                st.error(f"Integration failed: {str(e)}")
                st.stop()
    
    # Display results if integration completed
    if 'factors' in st.session_state:
        
        factors = st.session_state['factors']
        integrator = st.session_state['integrator']
        
        st.header("3. Results")
        
        # Tabs for different views
        tab1, tab2, tab3, tab4 = st.tabs([
            "📊 Factor Projection", 
            "🔥 Loadings", 
            "🔗 Cross-Layer Correlation",
            "💾 Export"
        ])
        
        with tab1:
            st.subheader("Sample Projection onto Integrated Factors")
            
            # Select factors to plot
            col1, col2 = st.columns(2)
            with col1:
                x_factor = st.selectbox("X-axis Factor", range(1, n_components + 1), index=0)
            with col2:
                y_factor = st.selectbox("Y-axis Factor", range(1, n_components + 1), index=1)
            
            # Plot
            fig, ax = plt.subplots(figsize=(10, 6))
            scatter = ax.scatter(
                factors.iloc[:, x_factor - 1], 
                factors.iloc[:, y_factor - 1],
                c=range(len(factors)), 
                cmap='viridis', 
                s=150,
                alpha=0.7,
                edgecolors='black',
                linewidth=0.5
            )
            
            # Add labels
            for idx, label in enumerate(factors.index):
                ax.annotate(
                    label, 
                    (factors.iloc[idx, x_factor - 1], factors.iloc[idx, y_factor - 1]),
                    fontsize=9,
                    alpha=0.8
                )
            
            ax.set_xlabel(f'Factor {x_factor}', fontsize=12)
            ax.set_ylabel(f'Factor {y_factor}', fontsize=12)
            ax.set_title(f'{strategy.capitalize()} Fusion + {method.upper()}', fontsize=14)
            ax.grid(True, alpha=0.3)
            plt.colorbar(scatter, ax=ax, label='Sample Index')
            
            st.pyplot(fig)
            
            # Download button
            buf = BytesIO()
            fig.savefig(buf, format='png', dpi=300, bbox_inches='tight')
            buf.seek(0)
            st.download_button(
                "Download Plot",
                data=buf,
                file_name="factor_projection.png",
                mime="image/png"
            )
        
        with tab2:
            st.subheader("Top Features by Factor Loading")
            
            # Select layer and factor
            col1, col2 = st.columns(2)
            with col1:
                selected_layer = st.selectbox(
                    "Omics Layer",
                    options=layer_names,
                    format_func=lambda x: x.capitalize()
                )
            with col2:
                selected_factor = st.selectbox(
                    "Factor",
                    range(n_components),
                    format_func=lambda x: f"Factor {x + 1}"
                )
            
            n_top = st.slider("Number of top features", 10, 100, 20)
            
            # Get top features
            top_features = integrator.get_top_features(
                selected_layer, 
                selected_factor, 
                n_top=n_top,
                return_weights=True
            )
            
            # Display table
            st.dataframe(top_features, use_container_width=True)
            
            # Heatmap of loadings
            fig, ax = plt.subplots(figsize=(8, max(6, n_top * 0.3)))
            
            heatmap_data = integrator.loadings_[selected_layer].loc[
                top_features['Feature'].tolist(),
                [f'Factor{i+1}' for i in range(min(5, n_components))]
            ]
            
            sns.heatmap(
                heatmap_data,
                cmap='RdBu_r',
                center=0,
                annot=False,
                cbar_kws={'label': 'Loading'},
                ax=ax
            )
            ax.set_title(f'{selected_layer.capitalize()} Loadings', fontsize=14)
            ax.set_ylabel('Feature', fontsize=12)
            ax.set_xlabel('Factor', fontsize=12)
            
            st.pyplot(fig)
        
        with tab3:
            st.subheader("Cross-Layer Factor Correlations")
            
            # Select layers to compare
            col1, col2 = st.columns(2)
            with col1:
                layer1 = st.selectbox("Layer 1", layer_names, index=0, key='layer1')
            with col2:
                layer2 = st.selectbox("Layer 2", layer_names, index=1, key='layer2')
            
            if layer1 == layer2:
                st.warning("Select different layers to compare")
            else:
                # Compute correlations
                factors1 = integrator.factors_by_layer_[layer1]
                factors2 = integrator.factors_by_layer_[layer2]
                
                corr_matrix = pd.DataFrame(
                    np.corrcoef(factors1.T, factors2.T)[:n_components, n_components:],
                    index=[f'{layer1[:5].upper()}_F{i+1}' for i in range(n_components)],
                    columns=[f'{layer2[:5].upper()}_F{i+1}' for i in range(n_components)]
                )
                
                # Heatmap
                fig, ax = plt.subplots(figsize=(8, 6))
                sns.heatmap(
                    corr_matrix,
                    annot=True,
                    fmt='.2f',
                    cmap='coolwarm',
                    center=0,
                    vmin=-1,
                    vmax=1,
                    ax=ax,
                    cbar_kws={'label': 'Correlation'}
                )
                ax.set_title(f'{layer1.capitalize()} vs {layer2.capitalize()} Factor Correlations')
                
                st.pyplot(fig)
                
                # Display correlation table
                st.dataframe(corr_matrix.style.background_gradient(cmap='coolwarm', vmin=-1, vmax=1))
        
        with tab4:
            st.subheader("Export Results")
            
            # Export factors
            st.markdown("**Integrated Factors**")
            csv_factors = factors.to_csv()
            st.download_button(
                "Download Factor Scores (CSV)",
                data=csv_factors,
                file_name="integrated_factors.csv",
                mime="text/csv"
            )
            
            # Export loadings
            st.markdown("**Feature Loadings**")
            for layer in layer_names:
                csv_loadings = integrator.loadings_[layer].to_csv()
                st.download_button(
                    f"Download {layer.capitalize()} Loadings (CSV)",
                    data=csv_loadings,
                    file_name=f"{layer}_loadings.csv",
                    mime="text/csv",
                    key=f"download_{layer}"
                )

else:
    st.info("👆 Upload data files or click 'Load Example Data' to begin")
    
    # Display instructions
    st.markdown("""
    ### Getting Started
    
    1. **Upload your data**: Provide CSV files with log2-transformed values
       - Features should be in rows, samples in columns
       - First column should contain feature IDs
    
    2. **Configure parameters**: Use the sidebar to select:
       - Integration method (NMF, PCA, or Factor Analysis)
       - Strategy (early or late fusion)
       - Number of components to extract
    
    3. **Run analysis**: Click "Run Integration" to perform statistical integration
    
    4. **Explore results**: View factor projections, feature loadings, and cross-layer correlations
    
    ### Data Format Example
    
    ```
    ,Sample1,Sample2,Sample3
    Gene1,5.2,6.1,5.8
    Gene2,4.3,4.5,4.7
    ```
    """)

# Footer
st.markdown("---")
st.markdown("**TransNet** | Multi-Omics Network Integration | [GitHub](https://github.com/mmattano/transnet)")