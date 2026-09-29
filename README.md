# TransNet: Trans-Omics Network Reconstruction and Analysis

[![Python package](https://github.com/mmattano/transnet/actions/workflows/python-package.yml/badge.svg)](https://github.com/mmattano/transnet/actions/workflows/python-package.yml)
[![Update Networks](https://github.com/mmattano/transnet/actions/workflows/update.yml/badge.svg)](https://github.com/mmattano/transnet/actions/workflows/update.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

TransNet builds a **trans-omic network**, a typed, directed, signed regulatory
hierarchy connecting omic layers, and analyses it as a network.

```
signal  ->  TF  ->  gene  ->  enzyme protein  ->  REACTION  <-  metabolite
                                                                (substrate,
                                                                 product,
                                                                 allosteric
                                                                 regulator)
```

Every edge records **which** relationship it represents, **which direction** its
regulatory effect runs in, and **what evidence** it came from. That is what lets
TransNet answer the question trans-omics exists to ask:

> Is this reaction regulated through gene expression, or through its allosteric
> effectors, and do the two agree?

**Trans-omics is not multi-omics with more layers.** Multi-omics integration asks
which molecules co-vary across data types. Trans-omics asks how a signal
*propagates* through a biochemical network: which regulatory relationship carried
it, in which direction, and whether the observed changes are consistent with that
mechanism. Every analysis below is meaningful only because of the network: remove
the edges and none of it can be computed.

## Install

```bash
pip install -e .          # from a clone
pip install transnet      # from PyPI
```

## Quickstart

```python
from transnet import (
    load_example_network, load_example_omics,
    map_omics_to_network, reaction_regulation_table,
)

graph = load_example_network()          # bundled, offline, no credentials
tables = load_example_omics()

report = map_omics_to_network(
    graph, tables,
    id_column="id", log2fc_column="log2FC", qvalue_column="padj",
)
print(report)     # how many features matched a node, per layer

table = reaction_regulation_table(graph)
table[table["controversial"]]           # reactions whose two axes disagree
```

```
reaction  name                          gene_axis   metabolite_axis  controversial
R00299    hexokinase                    activated   inhibited        True
R00756    6-phosphofructokinase         activated   inhibited        True
R00835    glucose-6-phosphate dehydrog. activated   inhibited        True
R00341    phosphoenolpyruvate carboxyk. inhibited   activated        True
R00303    glucose-6-phosphatase         inhibited   activated        True
```

More enzyme meeting more inhibitor; less enzyme meeting more substrate. Neither
axis alone would have shown the conflict.

## The analysis catalogue

These analyses define the field. Each maps to one function and one paper; full
details in [`docs/source/transomics_analyses.rst`](docs/source/transomics_analyses.rst),
bibliography in [`CITATIONS.md`](CITATIONS.md).

| Analysis | Function | Source |
|---|---|---|
| Trans-omic network reconstruction | `Transnet.generate_graph`, `responsive_subnetwork` | Yugi *et al.* 2016 |
| Reaction regulation-axis attribution | `reaction_regulation_table` | Kokaji *et al.* 2020; Egami *et al.* 2021 |
| Per-pathway regulation balance | `regulation_axis_summary` | Egami *et al.* 2021 |
| Metabolite regulatory roles (which changed metabolites regulate an enzyme) | `metabolite_regulatory_roles`, `regulatory_role_enrichment` | Kokaji *et al.* 2020 |
| Signed regulatory-path tracing | `trace_regulatory_paths` | Kawata *et al.* 2018 |
| Cross-layer connectivity | `cross_layer_connectivity` | Sugimoto *et al.* 2024 |
| Trans-omic hub identification | `transomic_hubs` | Morita *et al.* 2025 |
| Temporal and dose structure on the network | `assign_temporal_parameters`, `temporal_network_structure` | Morita 2025; Kawata 2018 |
| Responsive transcription-factor inference | `transcription_factor_activity` | Kokaji *et al.* 2022; Maehara *et al.* 2025 |
| Transcript–protein concordance (post-transcriptional regulation) | `expression_concordance` | Liu, Beyer & Aebersold 2016 |
| Regulatory motifs (product inhibition, feedback, feed-forward) | `regulatory_motifs` | Milo *et al.* 2002 |
| Structural vulnerability (what holds the response together) | `structural_vulnerability` | Morita *et al.* 2025 |
| Convergence against a null model | `convergence_significance` | Maslov & Sneppen 2002 |

### Does it reproduce published trans-omics?

`notebooks/studies/obese_liver.py` runs the catalogue on the liver data of
Uematsu *et al.* (*iScience* 2022; wild-type and ob/ob mice, fasting and after
glucose) and checks the published claims one by one. It reproduces 4 of 6 it can
test -- healthy liver responds to glucose through metabolites, obese liver is
rewired through enzyme amount -- reports the two it does not, and sets every
answer beside what a per-layer analysis of the same data can say. See
[docs/source/published_study.rst](docs/source/published_study.rst).

Two further studies are written up the same way, each set beside the
conclusions its own authors published: brown adipocytes under norepinephrine
([docs/source/brown_adipocytes.rst](docs/source/brown_adipocytes.rst), the data
of Anagho-Mattanovich *et al.* 2025) and MoTrPAC endurance training across six
rat tissues ([docs/source/motrpac_study.rst](docs/source/motrpac_study.rst),
against the consortium's *Nature* 2024 paper). Both report what the catalogue
fails to show on that data as plainly as what it shows.

## Layers are optional

**No analysis requires a particular layer.** A study with only transcriptomics
and metabolomics is a normal study, not a degraded input. Each function
discovers what the network actually contains and degrades by reporting *weaker
evidence*, never by failing:

- `reaction_regulation_table` records which chain of layers supported each
  gene-axis call in `gene_axis_evidence`: `"protein"`, `"gene_protein"`,
  `"gene"`, `"tf_gene_protein"`, or `None`;
- `trace_regulatory_paths` infers its own starting layer and reports the
  hierarchy it used;
- an omics table supplied for a layer the network lacks is *reported* in the
  mapping report, not raised.

`notebooks/walkthroughs/regulatory_paths.py` demonstrates this by running the same analysis
twice, with and without the optional Signaling layer.

## Notebooks

Everything runnable lives in [`notebooks/`](notebooks/), as jupytext scripts:
plain `.py` files that are also notebooks. Open them in Jupyter
(`jupytext --to notebook notebooks/walkthroughs/build_network.py`), run them as scripts,
or read them rendered in the documentation.

The walkthroughs run offline against `data/example/`, a curated slice of mouse
hepatic glucose metabolism with real KEGG, UniProt and EC identifiers. They
finish in seconds and need no credentials:

```bash
make notebooks          # or: python notebooks/walkthroughs/build_network.py
```

| Notebook | Covers |
|---|---|
| `build_network` | the typed network: layers, edge types, signs, and how much of it crosses layers |
| `responsive_network` | mapping data on, the coverage report, and the responsive subnetwork |
| `reaction_regulation` | the flagship analysis: which axis regulates each reaction, and where they disagree |
| `regulatory_paths` | signed path tracing, and what happens when a layer is missing |
| `temporal_and_hubs` | hubs that join layers, and whether the wiring explains the timing |
| `compare_conditions` | two conditions as typed networks; which metabolites act back on enzymes |
| `network_topology` | statistics, centrality, communities, active modules, diffusion, and four ways to export |
| `transcription_factors` | which factors drove the responsive genes, and where the annotation runs out |
| `external_annotation` | Reactome, HMDB, BRENDA, ChIP-Atlas and KEGG called directly (needs network access) |

The studies need a built network and take minutes:

```bash
make studies            # brown adipocytes, MoTrPAC, Uematsu, Kokaji
```

| Notebook | Data |
|---|---|
| `brown_adipocytes` | norepinephrine-stimulated brown adipocytes, three layers, seven timepoints |
| `motrpac_rat` | MoTrPAC endurance training, six tissues on one rat network, with phospho-, acetyl- and ubiquitin-proteomics |
| `obese_liver` | Uematsu *et al.* 2022, fetched at run time, claims scored |
| `kokaji_liver` | Kokaji *et al.* 2020, a genome-wide liver time course in two genotypes |

`notebooks/extra/` holds analyses written for collaborators rather than for
this documentation.

### BRENDA credentials

The allosteric edges come from BRENDA. Register free at
[brenda-enzymes.org](https://www.brenda-enzymes.org/), then either export the
credentials, pass them inline, or let the scripts prompt:

```bash
export BRENDA_EMAIL='you@example.com' BRENDA_PASSWORD='...'
python maintenance/build_networks.py --organisms mouse --brenda
```

`Proteome.get_brenda_kinetics()` reads the same two variables and prompts if
they are unset. Without BRENDA the metabolite regulation axis is limited to
substrates and products, and the analyses say so.

## The network model

A `networkx.MultiDiGraph`. Both properties are load-bearing:

- **Directed**, because a substrate flows into a reaction and a product flows
  out; a transcription factor regulates its target and not the reverse.
- **Multigraph**, because glucose-6-phosphate is both the *product* of
  hexokinase and its *allosteric inhibitor*: two edges with opposite signs that
  a simple graph collapses into one meaningless edge.

```python
>>> graph.get_edge_data("C00668", "R00299")     # G6P -> hexokinase
{'allosteric_inhibition': {..., 'sign': -1, ...}}
>>> graph.get_edge_data("R00299", "C00668")     # hexokinase -> G6P
{'product': {..., 'sign': 1, ...}}
```

| Edge type | Direction | Sign | Source |
|---|---|---|---|
| `phosphorylation` | Signaling → Proteome | ± | KEGG / user phosphoproteomics |
| `kinase_tf` | Signaling → Proteome (TF) | ± | KEGG signaling pathways |
| `transcriptional_regulation` | Proteome (TF) → Transcriptome | 0 | ChIP-Atlas |
| `translation` | Transcriptome → Proteome | +1 | NCBI/UniProt identifier join |
| `protein_interaction` | Proteome ↔ Proteome | 0 | STRING |
| `catalysis` | Proteome → Reactions | +1 | KEGG (EC number) |
| `gene_catalysis` | Transcriptome → Reactions | +1 | KEGG (EC number) |
| `substrate` | Metabolome → Reactions | +1 | KEGG reaction equation |
| `product` | Reactions → Metabolome | +1 | KEGG reaction equation |
| `allosteric_activation` | Metabolome → Reactions | +1 | BRENDA |
| `allosteric_inhibition` | Metabolome → Reactions | −1 | BRENDA |
| `enzymatic` | Proteome → Metabolome | 0 | KEGG (EC number) |

Sign `0` means **unknown**, not "no effect": ChIP-Atlas says a factor binds a
promoter, not whether it activates or represses. Analyses count those steps
separately, so a tentative prediction is never presented as a confident one.

Where an algorithm needs a simple undirected graph, `to_simple_graph()` projects
one and records what it collapsed rather than discarding it silently.

Full details: [`docs/source/network_model.rst`](docs/source/network_model.rst).

## Building from databases

```python
from transnet import Transnet, Proteome, Metabolome, Reactions

proteome = Proteome()
proteome.populate(kegg_organism="mmu", ncbi_organism="10090")
proteome.get_interaction_partners()          # STRING
proteome.get_transcription_factor_targets()  # ChIP-Atlas
proteome.get_brenda_kinetics()               # BRENDA: allosteric effectors

network = Transnet(proteome=proteome, metabolome=metabolome, reactions=reactions)
graph = network.generate_graph()
network.save_network("data/mouse/latest")
```

ChIP-Atlas lists every gene with any detectable binding, so
`get_transcription_factor_targets()` keeps only targets whose *mean* binding
score across the factor's experiments reaches `score_threshold` (default 100).
Without it a factor brings in ~15,000 targets and transcriptional regulation
becomes 91% of the network. Pass `score_threshold=0` for the permissive
behaviour; see [the network model docs](docs/source/network_model.rst) for the
score scale.

Networks for human, mouse, rat, yeast and *E. coli* are built by
`maintenance/build_networks.py`, each written to `data/<organism>/<date>/` as
`nodes.csv` plus `interactions.csv`. Rat is the network behind the MoTrPAC
study, mouse the one behind the Uematsu study. They are large (the human one is
~330 MB) and regenerable, so they are not kept in the repository. Build them
with `make networks`.

## Visualization

```python
from transnet.visualization import (
    plot_transomic_network,     # 2.5D stacked layers (transomics2cytoscape style)
    plot_regulation_axes,       # per-pathway axis balance
    plot_layer_connectivity,    # layer x layer edge heatmap
)
```

Networks also export to Cytoscape, Arena3D and transomics2cytoscape via
`transnet.io`.

## Factors, read through the network

`transnet.analysis.factors` fits NMF factors and then reads them *through the
network*: which part of the design each factor follows, how much it
reconstructs, and whether
its strongest features are actually connected. Other decompositions are not
reimplemented here; `sklearn.decomposition.PCA` is one line and always
current.

## Tests

```bash
pytest                  # offline, ~10 s
pytest -m network       # additionally hits live KEGG/UniProt
```

## Citing

TransNet implements published methods; each carries its reference in its
docstring. The framework as a whole is due to Yugi K, Kubota H, Hatano A,
Kuroda S. *Trans-Omics: How To Reconstruct Biochemical Networks Across Multiple
'Omic' Layers.* Trends in Biotechnology 34(4):276-290, 2016. Full bibliography
in [`CITATIONS.md`](CITATIONS.md).

## License

MIT; see [LICENSE](LICENSE).
