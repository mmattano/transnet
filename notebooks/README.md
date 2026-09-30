# Notebooks

Every runnable thing in this repository. They are
[jupytext](https://jupytext.readthedocs.io/) scripts: plain `.py` files that are
also notebooks, so they run as scripts, open in Jupyter, and render in the
documentation without a third copy of the content to keep in step.

The filenames carry no numbers. The order the walkthroughs are meant to be read
in is stated in `docs/source/index.rst` and in the Makefile, and
`tests/test_notebooks.py` fails if those two disagree.

```bash
python notebooks/walkthroughs/reaction_regulation.py            # as a script
jupytext --to notebook notebooks/walkthroughs/reaction_regulation.py   # as a notebook
make notebooks                                        # every offline walkthrough
```

## Walkthroughs (offline, seconds)

Run against `transnet/data/example/`: a curated slice of mouse hepatic glucose
metabolism with real KEGG, UniProt and EC identifiers. No credentials, no
network access.

| | Covers |
|---|---|
| `build_network` | the typed network: layers, edge types, signs, the builder API, and the organism registry |
| `responsive_network` | mapping data onto the network, the coverage report, and the responsive subnetwork |
| `reaction_regulation` | which axis regulates each reaction, where the axes disagree, per-pathway balance, the phospho axis, transcript–protein concordance |
| `regulatory_paths` | predicting a metabolite's direction along signed paths, from receptor to metabolite, and by propagation |
| `temporal_and_hubs` | molecules that connect layers, and whether the network explains response timing |
| `compare_conditions` | two conditions compared as typed networks; which changed metabolites regulate enzymes |
| `network_topology` | statistics, centrality, communities, active modules, motifs, weak points, a null model, diffusion |
| `export_network` | every export format: CSV, adjacency matrix, Cytoscape, Arena3D, transomics2cytoscape, HTML, static figures |
| `transcription_factors` | inferring which factors drove the responsive genes, and where the annotation runs out |
| `external_annotation` | Reactome, HMDB, BRENDA, ChIP-Atlas and KEGG called directly (**needs network access**) |

The walkthroughs work on differential tables, which is what the bundled example
holds. The analyses that need **sample-level** matrices (factors read through
the network, network-guided imputation, cross-layer agreement) run on the
studies below, where the data has samples; the rendered study notebooks in the
documentation show every one of them with its output.

## Studies (need a built network, minutes)

```bash
python maintenance/build_networks.py --organisms mouse --brenda   # once
make studies
```

| | Data |
|---|---|
| `brown_adipocytes` | norepinephrine-stimulated brown adipocytes (Anagho-Mattanovich *et al.*, iScience 2025) |
| `motrpac_rat` | MoTrPAC endurance training, six tissues on one rat network |
| `obese_liver` | lean and obese mouse liver after glucose, 19-enzyme panel (Uematsu *et al.* 2022), fetched at run time |
| `liver_timecourse` | lean and obese mouse liver after glucose, genome-wide over four hours (Kokaji *et al.* 2020), fetched at run time |

Each writes CSVs, figures and network exports under `data/<study>_results/`. The narrative
versions, with the figures and the comparison against each paper's own
conclusions, are in `docs/source/`.

## Conventions

* Loading data is not analysis: the study loaders live in `transnet.datasets`,
  so a notebook opens with a question rather than with parsing.
* Prose goes in markdown cells, not `print()` calls.
* Results are stated with their baseline (`versus_chance` for directional
  predictions, permutation tests for network coherence), and failures are
  reported as failures.
