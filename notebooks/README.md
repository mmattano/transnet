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

Run against `data/example/`: a curated slice of mouse hepatic glucose
metabolism with real KEGG, UniProt and EC identifiers. No credentials, no
network access.

| | Covers |
|---|---|
| `build_network` | the typed network: layers, edge types, signs, and how much of it crosses layers |
| `responsive_network` | mapping data on, the coverage report, and the responsive subnetwork |
| `reaction_regulation` | the flagship analysis: which axis regulates each reaction, and where they disagree |
| `regulatory_paths` | signed path tracing, and what happens when a layer is missing |
| `temporal_and_hubs` | hubs that join layers, and whether the wiring explains the timing |
| `compare_conditions` | two conditions as typed networks; which metabolites act back on enzymes |
| `network_topology` | statistics, centrality, communities, active modules, diffusion, and four ways to export |
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
| `obese_liver` | Uematsu *et al.* 2022, fetched at run time, published claims scored |
| `kokaji_liver` | Kokaji *et al.* 2020: genome-wide liver time course, WT and ob/ob, fetched at run time |

Each writes CSVs and figures under `data/<study>_results/`. The narrative
versions, with the figures and the comparison against each paper's own
conclusions, are in `docs/source/`.

## extra/

Analyses written for a specific collaboration rather than for the
documentation. `extra/oslo2_breast_cancer.py` compares clinical outcome groups
in a breast cancer cohort; its README states what data it expects.

## Conventions

* Loading data is not analysis: the study loaders live in `transnet.datasets`,
  so a notebook opens with a question rather than with parsing.
* Prose goes in markdown cells, not `print()` calls.
* Results are stated with their baseline (`versus_chance` for directional
  predictions, permutation tests for network coherence), and failures are
  reported as failures.
