# Citations

The bibliography lives in [`docs/source/references.rst`](docs/source/references.rst)
and is rendered as the **References** page of the documentation. Every paper is
written out once there, so a correction is made in one place, and
`python maintenance/verify_references.py` re-queries Crossref for each DOI and
fails if an author, journal, volume or year has drifted.

The analysis catalogue (`docs/source/transomics_analyses.rst`) names the source
of each analysis in a `:Reference:` field, and `tests/test_catalogue.py` fails if
an analysis has none. The docstrings describe what the code does and do not
repeat the references.

The framework TransNet implements, a biochemical network reconstructed across
omic layers through a defined set of connection technologies, is due to the
Kuroda group: **Yugi K, Kubota H, Hatano A, Kuroda S.** Trans-Omics: How To
Reconstruct Biochemical Networks Across Multiple 'Omic' Layers. *Trends in
Biotechnology* 34(4):276-290, 2016.
[doi:10.1016/j.tibtech.2015.12.013](https://doi.org/10.1016/j.tibtech.2015.12.013)

## Which reference belongs to which analysis

| Analysis | Implemented by | Reference |
|----------|----------------|-----------|
| Trans-omic network reconstruction | `responsive_subnetwork` | Yugi 2016 |
| Reaction regulation-axis attribution | `reaction_regulation_table` | Kokaji 2020; Egami 2021 |
| Per-pathway regulation balance | `regulation_axis_summary` | Egami 2021 |
| Metabolite regulatory roles | `metabolite_regulatory_roles`, `regulatory_role_enrichment` | Kokaji 2020 |
| Signed regulatory-path tracing | `trace_regulatory_paths`, `path_consistency_summary` | Kawata 2018 |
| Signal flow through the hierarchy | `kegg_signaling_relations`, `map_modification_sites`, `hierarchical_propagation`, `downstream_influence` | Yugi 2016; Kawata 2018; Ohno 2020 |
| Cross-layer connectivity | `cross_layer_connectivity`, `layer_coverage` | Sugimoto 2024 (iTraNet) |
| Trans-omic hub identification | `transomic_hubs` | Morita 2025; De Domenico 2015 |
| Temporal and dose structure | `assign_temporal_parameters`, `temporal_network_structure`, `split_by_response_class` | Morita 2025; Kawata 2018 |
| Responsive transcription-factor inference | `transcription_factor_activity` | Kokaji 2022; Maehara 2025 |
| Regulatory motifs | `regulatory_motifs` | Milo 2002 |
| Structural vulnerability | `structural_vulnerability` | Morita 2025 |
| Convergence against a null model | `convergence_significance` | Maslov & Sneppen 2002 |
| Active modules | `find_active_modules` | Ideker 2002 |
| Transcript–protein concordance | `expression_concordance` | Liu 2016; Egami 2021 |
| Factors read through the network | `fit_factors` and the `factor_*` readings | Argelaguet 2020; Cowen 2017 |

## Databases

The network is assembled from these, and they should be cited alongside it:
[KEGG](https://www.kegg.jp/), [UniProt](https://www.uniprot.org/),
[BRENDA](https://www.brenda-enzymes.org/), [STRING](https://string-db.org/),
[ChIP-Atlas](https://chip-atlas.org/), [Ensembl](https://www.ensembl.org/),
[Reactome](https://reactome.org/) and [HMDB](https://hmdb.ca/).
