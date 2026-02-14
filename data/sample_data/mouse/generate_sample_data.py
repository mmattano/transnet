import pandas as pd
import numpy as np

np.random.seed(42)

# Generate 10 samples instead of just values
n_samples = 10

# Transcriptomics: 50 genes x 10 samples
genes = [f"ENSMUSG{i:011d}" for i in range(1, 51)]
gene_names = [f"Gene{i}" for i in range(1, 51)]

# Create sample matrix
rna_data = np.random.randn(50, n_samples) * 0.5
rna_data[:10] = np.random.randn(10, n_samples) * 2 + 1.5

# Create DataFrame with samples as columns
rna_df = pd.DataFrame(rna_data, columns=[f'Sample_{i+1}' for i in range(n_samples)])
rna_df.insert(0, 'gene_id', genes)
rna_df.insert(1, 'gene_name', gene_names)

# Add fold change and p-value as aggregate stats
rna_df['log2FC'] = rna_df[[f'Sample_{i+1}' for i in range(n_samples)]].mean(axis=1)
rna_df['p_value'] = np.random.uniform(0.001, 0.9, 50)
rna_df.loc[:9, 'p_value'] = np.random.uniform(0.001, 0.04, 10)

# Proteomics: 20 proteins x 10 samples
proteins = [f"P{i:05d}" for i in range(1, 21)]
protein_names = [f"Protein{i}" for i in range(1, 21)]

prot_data = np.random.randn(20, n_samples) * 0.7
# Correlate first 5 proteins with first 5 genes
prot_data[:5] = rna_data[:5] * 0.6 + np.random.randn(5, n_samples) * 0.3

prot_df = pd.DataFrame(prot_data, columns=[f'Sample_{i+1}' for i in range(n_samples)])
prot_df.insert(0, 'protein_id', proteins)
prot_df.insert(1, 'protein_name', protein_names)
prot_df['log2FC'] = prot_df[[f'Sample_{i+1}' for i in range(n_samples)]].mean(axis=1)
prot_df['p_value'] = np.random.uniform(0.05, 0.9, 20)
prot_df.loc[:4, 'p_value'] = np.random.uniform(0.001, 0.04, 5)

# Metabolomics: 10 metabolites x 10 samples
metabolites = [f"C{i:05d}" for i in range(1, 11)]
met_names = [f"Metabolite{i}" for i in range(1, 11)]

met_data = np.random.randn(10, n_samples) * 0.5
# Correlate first 3 metabolites with first 3 proteins
met_data[:3] = prot_data[:3] * 0.7 + np.random.randn(3, n_samples) * 0.2

met_df = pd.DataFrame(met_data, columns=[f'Sample_{i+1}' for i in range(n_samples)])
met_df.insert(0, 'metabolite_id', metabolites)
met_df.insert(1, 'metabolite_name', met_names)
met_df['log2FC'] = met_df[[f'Sample_{i+1}' for i in range(n_samples)]].mean(axis=1)
met_df['p_value'] = np.random.uniform(0.05, 0.9, 10)
met_df.loc[:2, 'p_value'] = np.random.uniform(0.001, 0.04, 3)

# Save
rna_df.to_csv('transcriptomics.csv', index=False)
prot_df.to_csv('proteomics.csv', index=False)
met_df.to_csv('metabolomics.csv', index=False)

print(f"Generated:\n  {len(rna_df)} genes x {n_samples} samples")
print(f"  {len(prot_df)} proteins x {n_samples} samples")
print(f"  {len(met_df)} metabolites x {n_samples} samples")