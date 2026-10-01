# TransNet: Trans-Omics Network Reconstruction and Analysis

[![Tests](https://github.com/mmattano/transnet/actions/workflows/tests.yml/badge.svg)](https://github.com/mmattano/transnet/actions/workflows/tests.yml)
[![Docs](https://github.com/mmattano/transnet/actions/workflows/docs.yml/badge.svg)](https://mmattano.github.io/transnet/)
[![PyPI](https://img.shields.io/pypi/v/transnet.svg)](https://pypi.org/project/transnet/)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue.svg)](pyproject.toml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Documentation: [mmattano.github.io/transnet](https://mmattano.github.io/transnet/)**

TransNet builds a **trans-omic network**, in which molecules from different
omic layers are joined by typed, directed and signed regulatory edges, and
analyses it.

```
signal  ->  TF  ->  gene  ->  enzyme protein  ->  REACTION  <-  metabolite
                                                                (substrate,
                                                                 product,
                                                                 allosteric
                                                                 regulator)
```

Every edge records **which** relationship it represents, whether it increases
or decreases its target, and **which database** it comes from. This lets
TransNet answer the central question of trans-omics:

> Is this reaction regulated through the amount of its enzyme, or through the
> metabolites acting on it, and do the two agree?

Multi-omics integration usually asks which molecules change together across
data types. Trans-omics asks how a change travels through the biochemical
network: which regulatory relationship carried it, in which direction, and
whether the measured changes are consistent with that mechanism. None of the
analyses below can be computed without the network's edges.

## Install

```bash
pip install transnet      # from PyPI
pip install -e .          # from a clone
```

## Quickstart

```python
from transnet import (
    load_example_network, load_example_omics,
    map_omics_to_network, reaction_regulation_table,
)

graph = load_example_network()          # bundled, works offline
tables = load_example_omics()

report = map_omics_to_network(
    graph, tables,
    id_column="id", log2fc_column="log2FC", qvalue_column="padj",
)
print(report)     # how many rows of each table matched a node

table = reaction_regulation_table(graph)
table[table["controversial"]]           # reactions whose two axes disagree
```

```
reaction                              name  gene_axis  metabolite_axis  controversial
  R00299                        hexokinase          1               -1           True
  R00756             6-phosphofructokinase          1               -1           True
  R00835 glucose-6-phosphate dehydrogenase          1               -1           True
  R00341 phosphoenolpyruvate carboxykinase         -1                1           True
  R00303             glucose-6-phosphatase         -1                1           True
```

In the first three reactions the enzyme increased (`gene_axis` +1) while the
metabolites acting on it pushed the reaction down (`metabolite_axis` -1); in
the last two it is the other way round. Neither layer alone would show the
conflict.

## The analysis catalogue

Each analysis maps to one function and one publication. Details are in the
[analysis catalogue](https://mmattano.github.io/transnet/transomics_analyses.html),
and the bibliography in [`CITATIONS.md`](CITATIONS.md).

| Analysis | Function | Source |
|---|---|---|
| Trans-omic network reconstruction | `Transnet.generate_graph`, `responsive_subnetwork` | Yugi *et al.* 2016 |
| Reaction regulation-axis attribution | `reaction_regulation_table` | Kokaji *et al.* 2020; Egami *et al.* 2021 |
| Per-pathway regulation balance | `regulation_axis_summary`, `kegg_reaction_pathways` | Egami *et al.* 2021 |
| Metabolite regulatory roles (which changed metabolites regulate an enzyme) | `metabolite_regulatory_roles`, `regulatory_role_enrichment` | Kokaji *et al.* 2020 |
| Signed regulatory-path tracing | `trace_regulatory_paths`, `path_consistency_summary` | Kawata *et al.* 2018 |
| Signal flow and propagation | `hierarchical_propagation`, `downstream_influence` | Yugi *et al.* 2016 |
| Cross-layer connectivity | `cross_layer_connectivity` | Sugimoto *et al.* 2024 |
| Trans-omic hub identification | `transomic_hubs` | Morita *et al.* 2025 |
| Temporal and dose structure on the network | `assign_temporal_parameters`, `temporal_network_structure` | Morita 2025; Kawata 2018 |
| Responsive transcription-factor inference | `transcription_factor_activity` | Kokaji *et al.* 2022; Maehara *et al.* 2025 |
| Transcript–protein concordance (post-transcriptional regulation) | `expression_concordance` | Liu, Beyer & Aebersold 2016 |
| Regulatory motifs (product inhibition, feedback, feed-forward) | `regulatory_motifs` | Milo *et al.* 2002 |
| Structural vulnerability (what holds the response together) | `structural_vulnerability` | Morita *et al.* 2025 |
| Convergence against a null model | `convergence_significance` | Maslov & Sneppen 2002 |
| Comparing conditions | `compare_transomic_networks` | Egami *et al.* 2021 |

Every analysis has a figure in `transnet.visualization`, and factor models can
be read through the network with `transnet.analysis.factors`.

## Studies on published data

Four studies run the whole catalogue on published data. Each has a
[summary page](https://mmattano.github.io/transnet/) and a notebook with every
table and figure, and each reports what the analysis fails to show as plainly
as what it shows.

| Study | Data | Question |
|---|---|---|
| [Brown adipocytes](docs/source/brown_adipocytes.rst) | mouse brown fat cells stimulated with norepinephrine; three layers, up to seven time points (Anagho-Mattanovich *et al.* 2025) | which reactions carry heat production, and through which mechanism |
| [MoTrPAC](docs/source/motrpac_study.rst) | endurance training in six rat tissues, with phospho-, acetyl- and ubiquitin-proteomics (MoTrPAC 2024) | do tissues respond through the same mechanisms; how much lies in enzyme modification |
| [Obese liver: metabolic panel](docs/source/obese_liver_panel.rst) | lean and obese mouse liver, fasted and after glucose; 19 enzymes (Uematsu *et al.* 2022) | which reactions change in obese liver, and through which mechanism |
| [Obese liver: time course](docs/source/liver_timecourse.rst) | the same comparison genome-wide, over four hours (Kokaji *et al.* 2020) | how the two genotypes differ as networks over time |

## Layers are optional

**No analysis requires a particular layer.** A study with only
transcriptomics and metabolomics is a normal study. Each function uses the
layers the network has and reports weaker evidence, rather than failing, when
a layer is missing:

- `reaction_regulation_table` records in `gene_axis_evidence` which layers
  supported each enzyme-axis call: `"protein"`, `"gene_protein"`, `"gene"`,
  `"tf_gene_protein"`, or `None`;
- `trace_regulatory_paths` starts from the highest layer present and reports
  which one that was;
- a table supplied for a layer the network lacks is reported in the mapping
  report, not raised as an error.

## Notebooks

Everything runnable is in [`notebooks/`](notebooks/), as jupytext scripts:
plain `.py` files that are also notebooks. Open them in Jupyter
(`jupytext --to notebook notebooks/walkthroughs/build_network.py`), run them as
scripts, or read them rendered in the documentation.

The walkthroughs run offline on `transnet/data/example/`, a small curated part of mouse
liver glucose metabolism with real KEGG, UniProt and EC identifiers. They finish
in seconds and need no credentials:

```bash
make notebooks
```

| Walkthrough | Covers |
|---|---|
| `build_network` | the network: layers, edge types, signs, the builder API, the organism registry |
| `responsive_network` | mapping data onto the network, the coverage report, the responsive subnetwork |
| `reaction_regulation` | which axis regulates each reaction, controversial reactions, per-pathway balance, the phospho axis, transcript–protein concordance |
| `regulatory_paths` | predicting a metabolite's direction along signed paths, from a receptor to a metabolite, and by propagation |
| `temporal_and_hubs` | molecules that connect layers, and whether the network explains response timing |
| `compare_conditions` | two conditions compared as networks; which changed metabolites regulate enzymes |
| `network_topology` | statistics, centrality, communities, active modules, motifs, weak points, a null model, diffusion |
| `export_network` | every export format: CSV, adjacency matrix, Cytoscape, Arena3D, transomics2cytoscape, HTML, static figures |
| `transcription_factors` | which factors drove the changed genes, and where the annotation runs out |
| `external_annotation` | Reactome, HMDB, BRENDA, ChIP-Atlas and KEGG called directly (needs network access) |

The studies need a built network and take minutes to hours:

```bash
make studies
```

## Organism networks

Reference networks for human, mouse, rat, yeast and *E. coli* are built from
KEGG, UniProt, Ensembl, STRING, ChIP-Atlas and BRENDA by
`maintenance/build_networks.py`, and written to `data/<organism>/<date>/` as
`nodes.csv` and `interactions.csv`. They are several hundred megabytes and can
be rebuilt at any time, so they are not kept in the repository.

**The networks are updated by hand**, not on a schedule, because a rebuild
changes the numbers in every study and should be checked:

```bash
export BRENDA_EMAIL='you@example.com' BRENDA_PASSWORD='...'   # free registration
make networks ORGANISMS=mouse        # 30-90 min; resumes if interrupted
make studies                         # rerun the studies, refresh the doc figures
```

Read the build report printed at the end: a missing layer or edge type is
reported, with a non-zero exit status. Without BRENDA credentials the network
has no allosteric edges, and the metabolite axis is limited to substrates and
products. Adding an organism takes one call to
`transnet.organisms.register_organism`; see
[organism networks](https://mmattano.github.io/transnet/organism_networks.html).

## The network model

A `networkx.MultiDiGraph`. Both properties matter:

- **Directed**, because a substrate flows into a reaction and a product flows
  out, and a transcription factor regulates its target, not the reverse.
- **Multigraph**, because glucose 6-phosphate is both the *product* of
  hexokinase and its *allosteric inhibitor*: two edges with opposite signs,
  which a simple graph would merge into one.

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
promoter, not whether it activates or represses it. Analyses count such steps
separately, so a tentative prediction is never presented as a confident one.

Full details: [the network model](https://mmattano.github.io/transnet/network_model.html).

## Export

A network, or any part of it, can be written for other tools
(see the [export walkthrough](notebooks/walkthroughs/export_network.py)):

| Format | Function |
|---|---|
| CSV (`nodes.csv`, `interactions.csv`) | `transnet.io.write_network` / `read_network` |
| Adjacency matrix | `Transnet.generate_adjacency_matrix` |
| Cytoscape JSON | `transnet.io.to_cytoscape_json` |
| Arena3D Web | `transnet.io.to_arena3d` |
| transomics2cytoscape (R) | `transnet.io.to_transomics2cytoscape` |
| Self-contained interactive HTML | `plot_transomic_network_interactive(...).write_html(...)` |

## Tests

```bash
make test               # offline unit tests, about 15 s
pytest -m slow          # also runs every offline walkthrough
pytest -m network       # also calls the live databases
```

## Citing

A TransNet manuscript is in preparation. Each analysis names the publication
it comes from in the [analysis catalogue](https://mmattano.github.io/transnet/transomics_analyses.html).
The framework as a whole builds on Yugi K, Kubota H, Hatano A, Kuroda S.
*Trans-Omics: How To Reconstruct Biochemical Networks Across Multiple 'Omic'
Layers.* Trends in Biotechnology 34(4):276-290, 2016. The full bibliography is
in [`CITATIONS.md`](CITATIONS.md).

## License

MIT; see [LICENSE](LICENSE).
