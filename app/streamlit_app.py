import streamlit as st
import pandas as pd
import sys
from pathlib import Path
import networkx as nx

# Add transnet to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from transnet.biology.transnet import Transnet
from transnet.analysis.multi_omics_integration import MultiOmicsIntegrator
from transnet.visualization.network_vis import plot_network_basic, plot_multilayer_network

st.set_page_config(
    page_title="TransNet - Trans-Omics Network Integration",
    layout="wide"
)

st.title("🧬 TransNet: Trans-Omics Network Integration")

# Sidebar
st.sidebar.header("Configuration")
organism = st.sidebar.selectbox("Organism", ["mouse", "human"])

# Initialize network as None
if 'network' not in st.session_state:
    st.session_state['network'] = None
    st.session_state['network_loaded'] = False

# Network loading section in sidebar
with st.sidebar:
    st.markdown("---")
    st.subheader("Network")
    
    if st.button("Load Pre-built Network", type="primary"):
        with st.spinner(f"Loading {organism} network..."):
            try:
                network = Transnet.load_network(f"data/{organism}/latest")
                st.session_state['network'] = network
                st.session_state['network_loaded'] = True
                st.success(f"✅ Loaded {organism} network")
            except Exception as e:
                st.error(f"Failed to load network: {e}")
                st.session_state['network_loaded'] = False
    
    if st.session_state['network_loaded']:
        st.success(f"Network: {st.session_state['network'].name}")
        # Show basic stats
        try:
            if st.session_state['network'].proteome:
                st.write(f"Proteins: {len(st.session_state['network'].proteome.proteins)}")
            if st.session_state['network'].metabolome:
                st.write(f"Metabolites: {len(st.session_state['network'].metabolome.metabolites)}")
        except:
            pass
    else:
        st.info("No network loaded. Analysis will use statistics only.")

# Tabs
tab1, tab2, tab3, tab4 = st.tabs(["📊 Upload Data", "🔬 Cross-Layer Discovery", "📈 Visualization", "📋 Results"])

# Tab 1: Upload
with tab1:
    st.header("Upload Multi-Omics Data")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.subheader("Transcriptomics")
        rna_file = st.file_uploader("RNA-seq (CSV)", type=['csv'], key='rna')
        if st.button("Load Sample RNA"):
            try:
                st.session_state['rna_data'] = pd.read_csv('data/sample_data/mouse/transcriptomics.csv')
                st.success("Loaded sample transcriptomics")
            except Exception as e:
                st.error(f"Failed to load: {e}")
        
        if 'rna_data' in st.session_state:
            st.write(f"{len(st.session_state['rna_data'])} genes")
            st.dataframe(st.session_state['rna_data'].head())
    
    with col2:
        st.subheader("Proteomics")
        prot_file = st.file_uploader("Proteomics (CSV)", type=['csv'], key='prot')
        if st.button("Load Sample Protein"):
            try:
                st.session_state['prot_data'] = pd.read_csv('data/sample_data/mouse/proteomics.csv')
                st.success("Loaded sample proteomics")
            except Exception as e:
                st.error(f"Failed to load: {e}")
        
        if 'prot_data' in st.session_state:
            st.write(f"{len(st.session_state['prot_data'])} proteins")
            st.dataframe(st.session_state['prot_data'].head())
    
    with col3:
        st.subheader("Metabolomics")
        met_file = st.file_uploader("Metabolomics (CSV)", type=['csv'], key='met')
        if st.button("Load Sample Metabolite"):
            try:
                st.session_state['met_data'] = pd.read_csv('data/sample_data/mouse/metabolomics.csv')
                st.success("Loaded sample metabolomics")
            except Exception as e:
                st.error(f"Failed to load: {e}")
        
        if 'met_data' in st.session_state:
            st.write(f"{len(st.session_state['met_data'])} metabolites")
            st.dataframe(st.session_state['met_data'].head())

# Tab 2: Analysis
with tab2:
    st.header("🔬 Cross-Layer Discovery")
    
    st.markdown("""
    Find mechanistically connected changes across omics layers using **network-aware correlation**.
    
    TransNet combines statistical correlation with biological network topology.
    """)
    
    if 'prot_data' not in st.session_state or 'met_data' not in st.session_state:
        st.warning("Upload proteomics and metabolomics data first")
    else:
        col1, col2 = st.columns(2)
        
        with col1:
            alpha = st.slider("Network weight (0=stats, 1=network)", 0.0, 1.0, 0.5)
            if not st.session_state['network_loaded']:
                st.info("No network loaded. Setting to 0 (pure statistics).")
                alpha = 0.0
        with col2:
            min_corr = st.slider("Min correlation", 0.0, 1.0, 0.3)
        
        if st.button("Run Analysis", type="primary"):
            with st.spinner("Running network-aware correlation..."):
                # Prepare data matrices (samples x features)
                prot = st.session_state['prot_data']
                met = st.session_state['met_data']

                # Extract sample columns (exclude ID, name, log2FC, p_value)
                sample_cols = [c for c in prot.columns if c.startswith('Sample_')]

                if not sample_cols:
                    st.error("Data must have Sample_1, Sample_2, ... columns")
                else:
                    prot_df = prot[sample_cols].T  # Transpose: samples as rows
                    prot_df.columns = prot['protein_id'].values
                    
                    met_df = met[sample_cols].T
                    met_df.columns = met['metabolite_id'].values
                    
                    # Use network if loaded, else None
                    network = st.session_state['network'] if st.session_state['network_loaded'] else None
                    integrator = MultiOmicsIntegrator(network=network)
                    
                    try:
                        results = integrator.network_aware_correlation(
                            data1=met_df,
                            data2=prot_df,
                            layer1='Metabolome',
                            layer2='Proteome',
                            alpha=alpha,
                            min_samples=3
                        )
                    
                        results = results[abs(results['correlation']) >= min_corr]
                        st.session_state['results'] = results
                        
                        st.success(f"Found {len(results)} connections")
                    except Exception as e:
                        st.error(f"Analysis failed: {e}")
                        st.exception(e)

# Tab 3: Visualization
with tab3:
    st.header("📈 Network Visualization")
    
    if not st.session_state['network_loaded']:
        st.warning("Load a pre-built network first (see sidebar)")
    else:
        network = st.session_state['network']
        
        if st.button("Generate Network Graph"):
            with st.spinner("Building network graph..."):
                try:
                    G = network.generate_graph()
                    st.session_state['graph'] = G
                    st.success(f"Generated graph with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
                except Exception as e:
                    st.error(f"Failed to generate graph: {e}")
                    st.exception(e)
        
        if 'graph' in st.session_state:
            G = st.session_state['graph']
            
            if G.number_of_nodes() == 0:
                st.warning("Network graph is empty")
            else:
                st.subheader("Network Overview")
                
                col1, col2, col3 = st.columns(3)
                col1.metric("Nodes", G.number_of_nodes())
                col2.metric("Edges", G.number_of_edges())
                col3.metric("Density", f"{nx.density(G):.4f}")
                
                # Visualization options
                viz_type = st.selectbox(
                    "Visualization Type",
                    ["Basic Network", "Multi-layer Network"]
                )
                
                show_labels = st.checkbox("Show labels", value=False)
                
                if st.button("Generate Visualization"):
                    with st.spinner("Creating visualization..."):
                        try:
                            if viz_type == "Basic Network":
                                st.subheader("Basic Network View")
                                fig = plot_network_basic(
                                    G, 
                                    figsize=(12, 10),
                                    show_labels=show_labels,
                                    title=f"{organism.capitalize()} Network"
                                )
                                st.pyplot(fig)
                            
                            elif viz_type == "Multi-layer Network":
                                st.subheader("Multi-layer Network View")
                                
                                color_by = st.selectbox(
                                    "Color nodes by",
                                    ["layer", "degree"]
                                )
                                
                                fig = plot_multilayer_network(
                                    G,
                                    node_color_by=color_by,
                                    show_labels=show_labels,
                                    title=f"{organism.capitalize()} Multi-layer Network",
                                    figsize=(14, 12)
                                )
                                st.pyplot(fig)
                        except Exception as e:
                            st.error(f"Visualization failed: {e}")
                            st.exception(e)

# Tab 4: Results
with tab4:
    st.header("📋 Results")
    
    if 'results' in st.session_state:
        results = st.session_state['results']
        
        col1, col2, col3 = st.columns(3)
        col1.metric("Total Connections", len(results))
        col2.metric("Network Supported", sum(results['network_prior'] > 0))
        col3.metric("High Confidence", sum(results['combined_score'] > 0.7))
        
        st.subheader("Cross-Layer Connections")
        st.dataframe(
            results[[
                'feature1', 'feature2', 'correlation', 
                'network_prior', 'combined_score', 'network_path'
            ]],
            height=400
        )
        
        # Add interpretation
        st.subheader("Interpretation")
        
        if not st.session_state['network_loaded']:
            st.info("ℹ️ Analysis performed without network (pure statistics)")
        elif sum(results['network_prior'] > 0) > 0:
            st.success(f"✅ Found {sum(results['network_prior'] > 0)} connections with network support")
            st.write("These connections are supported by known biological relationships in the network.")
        else:
            st.info("ℹ️ No network-supported connections found. This may indicate:")
            st.write("- Novel relationships not yet in databases")
            st.write("- Sample IDs don't match network IDs")
            st.write("- Statistical correlations without direct biological connection")
        
        # Export
        st.subheader("Export Results")
        csv = results.to_csv(index=False)
        st.download_button(
            "📥 Download Results (CSV)", 
            csv, 
            "transnet_results.csv", 
            "text/csv",
            help="Download all cross-layer connections as CSV"
        )
    else:
        st.info("Run cross-layer discovery analysis first")