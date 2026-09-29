# Trans-omic analysis of the Oslo2 breast cancer cohort

This folder contains one analysis notebook. It takes transcriptome, proteome
and metabolome from the same tumours and asks, per metabolic reaction,
**through which mechanism two clinical groups differ**: whether a reaction is
pushed by the amount of enzyme present, or by the metabolites acting on it,
and whether the two agree.

It is written against the data described in Sharma *et al.*, *Comprehensive
multi-omics analysis of breast cancer reveals distinct long-term prognostic
subtypes*, **Oncogenesis 13:22 (2024)**, and is complementary to the MOFA+
analysis published there: that analysis identifies the molecules separating
the prognostic clusters, this one asks which reactions they act on.

## What you need

Python 3.10+, and TransNet:

```bash
pip install -e /path/to/transnet          # or: pip install transnet
python maintenance/build_networks.py --organisms human --brenda
```

The last step builds the human network (KEGG reactions, UniProt proteins,
STRING interactions, ChIP-Atlas binding, BRENDA allosteric effectors). It
takes 30–90 minutes and needs free [BRENDA](https://www.brenda-enzymes.org/)
credentials in `BRENDA_EMAIL` and `BRENDA_PASSWORD`. It is done once.

## Your data

Put five CSVs in one directory and point `OSLO2_DATA` at it:

| File | Rows | Columns |
|---|---|---|
| `transcriptome.csv` | genes (symbol or Agilent probe id in the first column) | sample ids |
| `proteome.csv` | antibodies (name in the first column) | sample ids |
| `metabolome.csv` | metabolites (compound name in the first column) | sample ids |
| `clinical.csv` | sample ids in the first column | clinical variables |
| `antibody_map.csv` | one row per antibody | `antibody,uniprot` |

Notes:

* **Sample ids must agree** across all five files. Samples missing from a
  layer are simply not used for that layer; the notebook reports how many
  have all three.
* **Expression and RPPA values should be on a log scale** already (as
  distributed). Metabolite intensities may be raw.
* `antibody_map.csv` says which protein each antibody reports. Leave
  `uniprot` blank for phosphosite antibodies: the notebook keeps them aside
  rather than treating a phosphosite change as a change in protein amount.
* `clinical.csv` needs the grouping variable (`moc_cluster` by default). Add
  `survival_months` and `event` (1 = death) to enable the Cox section.

Nothing is uploaded anywhere and no data leaves your machine.

## Running it

```bash
export OSLO2_DATA=/path/to/your/csvs
export OSLO2_OUT=/path/to/results          # optional
jupytext --to notebook oslo2_breast_cancer.py   # optional: as a .ipynb
python oslo2_breast_cancer.py                   # or run the notebook
```

If a file is missing it stops on the first cell and names what it needs.

## What it produces

* `reaction_regulation_<test>_vs_<reference>.csv`: every reaction, the axis
  that carries it, the enzymes involved, the allosteric regulators, and
  whether the two axes disagree.
* `concordance.csv`: for each antibody's protein, whether its transcript
  moved the same way, plus a threshold-free test of whether the protein moved
  *further* than its transcript.
* `regulation_by_grouping.csv`: the same analysis across several clinical
  groupings, so that "what distinguishes poor from good prognosis" can be
  compared with "what distinguishes ER+ from ER−".
* `network.html`: the responsive network, interactive and self-contained. Open
  it in any browser; no Python needed.

## Changing the question

The first cell holds everything you are likely to change:

```python
GROUPING = "moc_cluster"      # any column of clinical.csv
REFERENCE, TEST = "MOC3", "MOC2"
QVALUE = 0.05
```

and the `GROUPINGS` list further down controls the multi-contrast comparison.

## Caveats worth knowing before reading the output

* **The metabolite panel is the limiting factor.** With 18 named compounds,
  the metabolite axis is visible only where those compounds act. A reaction
  reported as "enzyme axis only" may simply have no measured metabolite.
* **RPPA covers ~150 proteins**, so the transcript-protein comparison is
  precise but narrow.
* **"Protein only" depends on transcriptome power**, which is why the
  threshold-free difference test is reported beside it.
* This is a cross-sectional cohort, so everything here is an association
  between groups of patients, not a time course: the direction of an edge is
  the direction of the annotated regulation, not evidence of causation in
  these tumours.

Questions about the method: see `docs/source/transomics_analyses.rst` in the
TransNet repository, which documents each analysis and the paper it comes
from.
