# TransNet: Trans-Omics Network Integration and Analysis

[![Python package](https://github.com/mmattano/transnet/actions/workflows/python-package.yml/badge.svg)](https://github.com/mmattano/transnet/actions/workflows/python-package.yml)
[![Update Networks](https://github.com/mmattano/transnet/actions/workflows/update.yml/badge.svg)](https://github.com/mmattano/transnet/actions/workflows/update.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

TransNet is a Python package designed for the generation, integration, and analysis of multi-layer biological networks from trans-omics data. It enables researchers to:

1. Automatically build comprehensive biological networks by pulling data from various databases
2. Integrate experimental data from different omics layers (transcriptomics, proteomics, metabolomics)
3. Analyze and visualize integrated networks
4. Access pre-built networks for common model organisms

## Features

- **Database Integration**: Pull data from multiple biological databases (KEGG, UniProt, Ensembl, STRING, ChIP-Atlas)
- **Multi-omics Layers**: Build networks with transcriptome, proteome, metabolome, pathways, and reaction layers
- **Pre-built Networks**: Access regularly updated networks for human, mouse, yeast, and E. coli
- **Experimental Data Integration**: Map your experimental data onto networks
- **Network Analysis**: Perform various analyses including:
  - Centrality measures
  - Community detection
  - Differential network analysis
  - Enrichment analysis
  - Active module identification
- **Data Export**: Save networks and results in CSV format for easy sharing and further analysis

## Installation

```bash
# Install from PyPI
pip install transnet

# Or install from source
git clone https://github.com/mmattano/transnet.git
cd transnet
pip install -e .