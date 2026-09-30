"""Generate the bundled example trans-omic network in ``transnet/data/example/``.

A hand-curated slice of mouse hepatic glucose metabolism -- glycolysis,
gluconeogenesis, the pentose-phosphate entry point and the insulin signalling
branch above them -- with real KEGG, UniProt and EC identifiers.  Small enough
that every example script runs offline in seconds, and rich enough to show every
edge type, including reactions where the gene-expression and allosteric axes
disagree.

The network is produced by the ordinary :class:`~transnet.Transnet` builders, so
the bundled files cannot drift from the live schema.

Run::

    python maintenance/build_example_data.py
"""

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from transnet import Transnet                                    # noqa: E402
from transnet.biology.elements import (                          # noqa: E402
    Gene, Metabolite, Protein, Reaction,
)
from transnet.biology.layers import (                            # noqa: E402
    Metabolome, Proteome, Reactions, Signaling, Transcriptome,
)

OUTPUT_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "transnet", "data", "example"
)

# --------------------------------------------------------------------------
# Metabolites (KEGG compound ids)
# --------------------------------------------------------------------------
METABOLITES = {
    "C00031": "D-Glucose",
    "C00668": "alpha-D-Glucose 6-phosphate",
    "C00085": "D-Fructose 6-phosphate",
    "C00354": "D-Fructose 1,6-bisphosphate",
    "C00665": "beta-D-Fructose 2,6-bisphosphate",
    "C00111": "Glycerone phosphate",
    "C00118": "D-Glyceraldehyde 3-phosphate",
    "C00236": "3-Phospho-D-glyceroyl phosphate",
    "C00197": "3-Phospho-D-glycerate",
    "C00631": "2-Phospho-D-glycerate",
    "C00074": "Phosphoenolpyruvate",
    "C00022": "Pyruvate",
    "C00186": "(S)-Lactate",
    "C00024": "Acetyl-CoA",
    "C00158": "Citrate",
    "C00036": "Oxaloacetate",
    "C00041": "L-Alanine",
    "C00002": "ATP",
    "C00008": "ADP",
    "C00020": "AMP",
    "C00003": "NAD+",
    "C00004": "NADH",
    "C00006": "NADP+",
    "C00005": "NADPH",
    "C00009": "Orthophosphate",
    "C00001": "H2O",
    "C01236": "D-Glucono-1,5-lactone 6-phosphate",
}

# --------------------------------------------------------------------------
# Reactions: (id, name, EC, substrates, products, reversible)
# --------------------------------------------------------------------------
REACTIONS = [
    ("R00299", "hexokinase", "2.7.1.1",
     ["C00031", "C00002"], ["C00668", "C00008"], False),
    ("R00771", "glucose-6-phosphate isomerase", "5.3.1.9",
     ["C00668"], ["C00085"], True),
    ("R00756", "6-phosphofructokinase", "2.7.1.11",
     ["C00085", "C00002"], ["C00354", "C00008"], False),
    ("R00762", "fructose-1,6-bisphosphatase", "3.1.3.11",
     ["C00354", "C00001"], ["C00085", "C00009"], False),
    ("R01070", "fructose-bisphosphate aldolase", "4.1.2.13",
     ["C00354"], ["C00111", "C00118"], True),
    ("R01015", "triose-phosphate isomerase", "5.3.1.1",
     ["C00111"], ["C00118"], True),
    ("R01061", "glyceraldehyde-3-phosphate dehydrogenase", "1.2.1.12",
     ["C00118", "C00009", "C00003"], ["C00236", "C00004"], True),
    ("R01512", "phosphoglycerate kinase", "2.7.2.3",
     ["C00236", "C00008"], ["C00197", "C00002"], True),
    ("R01518", "phosphoglycerate mutase", "5.4.2.11",
     ["C00197"], ["C00631"], True),
    ("R00658", "enolase", "4.2.1.11",
     ["C00631"], ["C00074", "C00001"], True),
    ("R00200", "pyruvate kinase", "2.7.1.40",
     ["C00074", "C00008"], ["C00022", "C00002"], False),
    ("R00703", "L-lactate dehydrogenase", "1.1.1.27",
     ["C00022", "C00004"], ["C00186", "C00003"], True),
    ("R00209", "pyruvate dehydrogenase", "1.2.4.1",
     ["C00022", "C00003"], ["C00024", "C00004"], False),
    ("R00835", "glucose-6-phosphate dehydrogenase", "1.1.1.49",
     ["C00668", "C00006"], ["C01236", "C00005"], False),
    ("R00341", "phosphoenolpyruvate carboxykinase", "4.1.1.32",
     ["C00036", "C00002"], ["C00074", "C00008"], False),
    ("R00303", "glucose-6-phosphatase", "3.1.3.9",
     ["C00668", "C00001"], ["C00031", "C00009"], False),
]

# --------------------------------------------------------------------------
# Enzymes: (UniProt, symbol, Entrez, EC, activators, inhibitors)
# Allosteric effectors follow BRENDA's annotation for the mouse enzymes.
# --------------------------------------------------------------------------
ENZYMES = [
    ("P17710", "Hk1", "15275", "2.7.1.1", [], ["C00668"]),
    ("O08528", "Hk2", "15277", "2.7.1.1", [], ["C00668"]),
    ("P06745", "Gpi1", "14751", "5.3.1.9", [], []),
    ("P12382", "Pfkl", "18641", "2.7.1.11",
     ["C00665", "C00020"], ["C00002", "C00158"]),
    ("Q9QXD6", "Fbp1", "14121", "3.1.3.11", [], ["C00020", "C00665"]),
    ("P05064", "Aldoa", "11674", "4.1.2.13", [], []),
    ("P17751", "Tpi1", "21991", "5.3.1.1", [], []),
    ("P16858", "Gapdh", "14433", "1.2.1.12", [], []),
    ("P09411", "Pgk1", "18655", "2.7.2.3", [], []),
    ("Q9DBJ1", "Pgam1", "18648", "5.4.2.11", [], []),
    ("P17182", "Eno1", "13806", "4.2.1.11", [], []),
    ("P53657", "Pklr", "18770", "2.7.1.40", ["C00354"], ["C00002", "C00041"]),
    ("P06151", "Ldha", "16828", "1.1.1.27", [], []),
    ("P35486", "Pdha1", "18597", "1.2.4.1", [], ["C00024", "C00004"]),
    ("Q00612", "G6pdx", "14381", "1.1.1.49", [], ["C00005"]),
    ("Q9Z2V4", "Pck1", "18534", "4.1.1.32", [], []),
    ("P35576", "G6pc", "14377", "3.1.3.9", [], []),
]

# Transcription factors: (UniProt, symbol, Entrez, target gene symbols)
TRANSCRIPTION_FACTORS = [
    ("Q9WTN3", "Srebf1", "20787", ["Pklr", "Ldha", "Hk2"]),
    ("Q9R1E0", "Foxo1", "56458", ["G6pc", "Pck1", "Fbp1"]),
    ("P49698", "Hnf4a", "15378", ["Pck1", "G6pc", "Aldoa"]),
    ("Q99MZ3", "Mlxipl", "58805", ["Pklr", "Fbp1", "Ldha"]),
]

# Protein-protein interactions (STRING-style), (a, b, confidence)
INTERACTIONS = [
    ("P17710", "P06745", 0.90),
    ("P12382", "P05064", 0.85),
    ("P05064", "P17751", 0.93),
    ("P17751", "P16858", 0.88),
    ("P16858", "P09411", 0.91),
    ("P53657", "P06151", 0.82),
    ("Q9Z2V4", "P35576", 0.87),
]

# Signaling: insulin -> AKT -> transcription factors.
# sign -1 marks an inhibitory phosphorylation: AKT phosphorylating FOXO1
# excludes it from the nucleus, so an active AKT *represses* FOXO1 targets.
SIGNALING = [
    ("P15208", "Insr", "P35569", 1, "protein"),
    ("P35569", "Irs1", "P26450", 1, "protein"),
    ("P26450", "Pik3r1", "P31750", 1, "protein"),
    ("P31750", "Akt1", "Q9R1E0", -1, "tf"),
    ("P31750", "Akt1", "Q9WTN3", 1, "tf"),
    ("P31750", "Akt1", "P12382", 1, "protein"),
]
SIGNALING_NAMES = {
    "P15208": "Insr", "P35569": "Irs1",
    "P26450": "Pik3r1", "P31750": "Akt1",
}


def build():
    metabolome = Metabolome()
    metabolome.metabolites = [
        Metabolite(kegg_compound_id=cid, kegg_name=name)
        for cid, name in METABOLITES.items()
    ]

    reactions = Reactions()
    reactions.reactions = [
        Reaction(
            id=rid, name=name, enzyme=[ec],
            equation=f"{' + '.join(subs)} <=> {' + '.join(prods)}",
            substrates=list(subs), products=list(prods),
            stoichiometry_substrates=[1.0] * len(subs),
            stoichiometry_products=[1.0] * len(prods),
            reversible=reversible,
        )
        for rid, name, ec, subs, prods, reversible in REACTIONS
    ]

    proteome = Proteome()
    transcriptome = Transcriptome()
    proteins, genes = [], []

    for uniprot, symbol, entrez, ec, activators, inhibitors in ENZYMES:
        protein = Protein(
            uniprot_id=uniprot, name=symbol, gene=symbol,
            ec_number=[ec], entrez_id=[entrez], ncbi_organism="10090",
        )
        protein.activators = list(activators)
        protein.inhibitors = list(inhibitors)
        proteins.append(protein)

        gene = Gene(ncbi_id=entrez, name=symbol)
        gene.related_ecs = [ec]
        genes.append(gene)

    for uniprot, symbol, entrez, targets in TRANSCRIPTION_FACTORS:
        protein = Protein(
            uniprot_id=uniprot, name=symbol, gene=symbol,
            entrez_id=[entrez], ncbi_organism="10090",
        )
        protein.transcription_factor_targets = list(targets)
        proteins.append(protein)
        genes.append(Gene(ncbi_id=entrez, name=symbol))

    for uniprot, symbol in SIGNALING_NAMES.items():
        proteins.append(Protein(
            uniprot_id=uniprot, name=symbol, gene=symbol,
            ncbi_organism="10090",
        ))

    by_id = {p.uniprot_id: p for p in proteins}
    for a, b, confidence in INTERACTIONS:
        by_id[a].interaction_partners.append(b)
        by_id[b].interaction_partners.append(a)

    proteome.proteins = proteins
    proteome.interaction_scores = {
        (a, b): c for a, b, c in INTERACTIONS
    } | {(b, a): c for a, b, c in INTERACTIONS}
    transcriptome.genes = genes

    signaling = Signaling()
    signaling.populate_from_table(
        pd.DataFrame(
            [{"kinase": src, "substrate": tgt, "sign": sign, "target_type": kind}
             for src, _, tgt, sign, kind in SIGNALING],
        ),
        sign_column="sign", target_type_column="target_type",
        evidence="KEGG hsa04910 (insulin signaling)",
    )
    for node in signaling.nodes:
        node.name = SIGNALING_NAMES.get(node.id, node.id)

    return Transnet(
        name="example-mouse-hepatic-glucose-metabolism",
        reactions=reactions, metabolome=metabolome,
        proteome=proteome, transcriptome=transcriptome,
    ), signaling


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    network, signaling = build()

    # The default example network has no Signaling layer, because most studies
    # do not measure the phosphoproteome. Example 04 adds it from a side file.
    network.save_network(OUTPUT_DIR)

    signaling_rows = [
        {"kinase": src, "kinase_name": SIGNALING_NAMES[src],
         "substrate": tgt, "sign": sign, "target_type": kind}
        for src, _, tgt, sign, kind in SIGNALING
    ]
    pd.DataFrame(signaling_rows).to_csv(
        os.path.join(OUTPUT_DIR, "signaling_edges.csv"), index=False
    )

    _write_omics()
    _write_timecourse()
    _write_pathways()

    graph = network.generate_graph()
    print(f"Wrote example network to {OUTPUT_DIR}")
    print(f"  {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")


def _write_omics():
    """Differential results for a fed-vs-fasted style contrast.

    Chosen so that the two regulation axes genuinely disagree on some
    reactions: insulin signalling raises glycolytic enzyme expression while
    product/allosteric feedback pushes the other way.
    """
    transcriptome = pd.DataFrame([
        ("15275", "Hk1", 0.35, 0.21), ("15277", "Hk2", 1.62, 0.001),
        ("14751", "Gpi1", 0.12, 0.64), ("18641", "Pfkl", 1.05, 0.004),
        ("14121", "Fbp1", -1.48, 0.002), ("11674", "Aldoa", 0.44, 0.09),
        ("21991", "Tpi1", 0.08, 0.81), ("14433", "Gapdh", 0.21, 0.33),
        ("18655", "Pgk1", 0.31, 0.18), ("18648", "Pgam1", 0.05, 0.90),
        ("13806", "Eno1", 0.28, 0.24), ("18770", "Pklr", 1.85, 0.0005),
        ("16828", "Ldha", 1.24, 0.003), ("18597", "Pdha1", 0.52, 0.06),
        ("14381", "G6pdx", 0.95, 0.01), ("18534", "Pck1", -2.10, 0.0001),
        ("14377", "G6pc", -1.76, 0.0008), ("20787", "Srebf1", 1.32, 0.002),
        ("56458", "Foxo1", -0.18, 0.55), ("15378", "Hnf4a", -0.42, 0.12),
        ("58805", "Mlxipl", 1.10, 0.006),
    ], columns=["id", "symbol", "log2FC", "padj"])

    proteome = pd.DataFrame([
        ("P17710", "Hk1", 0.22, 0.44), ("O08528", "Hk2", 1.18, 0.008),
        ("P06745", "Gpi1", 0.09, 0.77), ("P12382", "Pfkl", 0.87, 0.02),
        ("Q9QXD6", "Fbp1", -1.05, 0.01), ("P05064", "Aldoa", 0.33, 0.19),
        ("P17751", "Tpi1", 0.04, 0.93), ("P16858", "Gapdh", 0.15, 0.51),
        ("P09411", "Pgk1", 0.26, 0.29), ("Q9DBJ1", "Pgam1", 0.02, 0.96),
        ("P17182", "Eno1", 0.19, 0.42), ("P53657", "Pklr", 1.41, 0.001),
        ("P06151", "Ldha", 0.98, 0.01), ("P35486", "Pdha1", 0.38, 0.15),
        ("Q00612", "G6pdx", 0.71, 0.03), ("Q9Z2V4", "Pck1", -1.62, 0.0006),
        ("P35576", "G6pc", -1.33, 0.003), ("Q9WTN3", "Srebf1", 1.09, 0.007),
        ("Q9R1E0", "Foxo1", -0.11, 0.68), ("P49698", "Hnf4a", -0.29, 0.31),
        ("Q99MZ3", "Mlxipl", 0.84, 0.02),
    ], columns=["id", "symbol", "log2FC", "padj"])

    metabolome = pd.DataFrame([
        ("C00031", "D-Glucose", -1.15, 0.004),
        ("C00668", "G6P", 1.42, 0.001),
        ("C00085", "F6P", 0.88, 0.02),
        ("C00354", "F1,6BP", 1.67, 0.0007),
        ("C00665", "F2,6BP", 1.95, 0.0003),
        ("C00111", "DHAP", 0.51, 0.08),
        ("C00118", "G3P", 0.46, 0.11),
        ("C00236", "1,3BPG", 0.22, 0.47),
        ("C00197", "3PG", 0.34, 0.20),
        ("C00631", "2PG", 0.29, 0.28),
        ("C00074", "PEP", -0.62, 0.04),
        ("C00022", "Pyruvate", 1.08, 0.006),
        ("C00186", "Lactate", 1.55, 0.001),
        ("C00024", "Acetyl-CoA", 0.79, 0.03),
        ("C00158", "Citrate", 0.94, 0.01),
        ("C00036", "Oxaloacetate", -0.55, 0.06),
        ("C00041", "L-Alanine", -0.71, 0.03),
        ("C00002", "ATP", 0.14, 0.62),
        ("C00008", "ADP", -0.09, 0.75),
        ("C00020", "AMP", -0.83, 0.02),
        ("C00003", "NAD+", -0.18, 0.49),
        ("C00004", "NADH", 0.66, 0.04),
        ("C00006", "NADP+", -0.24, 0.38),
        ("C00005", "NADPH", 0.58, 0.05),
        ("C00009", "Pi", -0.12, 0.71),
        ("C01236", "6-phosphoglucono-lactone", 0.61, 0.04),
    ], columns=["id", "name", "log2FC", "padj"])

    # A second condition, for the two-condition comparison example: the same
    # animals made insulin resistant, so the signalling-driven arm is blunted
    # while the allosteric arm still responds.
    #
    # Significance has to follow the effect size. Copying `padj` across while
    # scaling `log2FC` -- which this used to do -- leaves both conditions with
    # the *same set* of significant features, so every set-based statistic
    # (the metabolite regulatory-role enrichment in example 06, for one) comes
    # out bit-identical in both conditions by construction. Instead the
    # p-value is pushed back through the normal quantile it implies, the
    # underlying z is scaled by the same factor as the effect, and a new
    # p-value is read off. A blunted response then loses significance and an
    # amplified one gains it, as it would in a real experiment.
    def _rescale_padj(frame, scaled):
        """Recompute `padj` for effect sizes that were scaled by row."""
        from scipy.stats import norm

        original = frame["log2FC"].astype(float).to_numpy()
        new_effect = scaled["log2FC"].astype(float).to_numpy()
        padj = frame["padj"].astype(float).to_numpy()

        # two-sided p -> |z|
        z = norm.isf(np.clip(padj, 1e-12, 1 - 1e-12) / 2.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            factor = np.where(
                np.abs(original) > 1e-9,
                np.abs(new_effect) / np.abs(original),
                1.0,
            )
        factor = np.nan_to_num(factor, nan=1.0, posinf=1.0)
        return np.clip(2.0 * norm.sf(z * factor), 1e-6, 0.999).round(3)

    # The metabolome is blunted only mildly: the point of the comparison is
    # that the *allosteric* arm still responds while the signalling-driven
    # transcriptional arm collapses. Blunting it as hard as the transcriptome
    # leaves too few differential metabolites for the regulatory-role
    # enrichment in example 06 to say anything.
    resistant = metabolome.copy()
    resistant["log2FC"] = resistant["log2FC"] * 0.75
    resistant.loc[resistant["id"] == "C00668", "log2FC"] = 1.85
    resistant.loc[resistant["id"] == "C00031", "log2FC"] = 0.92
    resistant["padj"] = _rescale_padj(metabolome, resistant)

    # Insulin-driven transcription is blunted, but the gluconeogenic genes
    # insulin normally suppresses are now de-repressed and go the other way.
    def _blunt(frame, key_ids, key_values, key_padj):
        blunted = frame.copy()
        blunted["log2FC"] = blunted["log2FC"] * 0.25
        blunted.loc[blunted["id"].isin(key_ids), "log2FC"] = key_values
        blunted["padj"] = _rescale_padj(frame, blunted)
        # The de-repressed genes are measured as strongly as before, so their
        # significance is stated outright rather than derived.
        blunted.loc[blunted["id"].isin(key_ids), "padj"] = key_padj
        return blunted

    resistant_transcriptome = _blunt(
        transcriptome, ["18534", "14377"], [0.85, 0.62], [0.01, 0.02])
    resistant_proteome = _blunt(
        proteome, ["Q9Z2V4", "P35576"], [0.71, 0.55], [0.02, 0.03])

    for name, frame in [
        ("transcriptomics.csv", transcriptome),
        ("proteomics.csv", proteome),
        ("metabolomics.csv", metabolome),
        ("metabolomics_insulin_resistant.csv", resistant),
        ("transcriptomics_insulin_resistant.csv", resistant_transcriptome),
        ("proteomics_insulin_resistant.csv", resistant_proteome),
    ]:
        frame.to_csv(os.path.join(OUTPUT_DIR, name), index=False)

    phospho = pd.DataFrame([
        ("P15208", "Insr", 2.10, 0.0002),
        ("P35569", "Irs1", 1.74, 0.001),
        ("P26450", "Pik3r1", 1.32, 0.004),
        ("P31750", "Akt1", 2.35, 0.0001),
    ], columns=["id", "symbol", "log2FC", "padj"])
    phospho.to_csv(os.path.join(OUTPUT_DIR, "phosphoproteomics.csv"), index=False)


def _write_timecourse():
    """Metabolite time courses, in minutes after the glucose bolus.

    Hub metabolites (ATP, G6P) are given fast responses and peripheral ones
    slow responses, so the degree-versus-t-half relationship is visible.
    """
    rows = [
        ("C00668", 0.00, 1.05, 1.38, 1.42, 1.40),
        ("C00002", 0.00, 0.09, 0.13, 0.14, 0.14),
        ("C00085", 0.00, 0.48, 0.79, 0.88, 0.86),
        ("C00354", 0.00, 0.71, 1.35, 1.67, 1.64),
        ("C00665", 0.00, 1.42, 1.85, 1.95, 1.90),
        ("C00031", 0.00, -0.82, -1.08, -1.15, -1.12),
        ("C00111", 0.00, 0.11, 0.29, 0.44, 0.51),
        ("C00118", 0.00, 0.08, 0.24, 0.39, 0.46),
        ("C00197", 0.00, 0.05, 0.14, 0.26, 0.34),
        ("C00631", 0.00, 0.04, 0.11, 0.22, 0.29),
        ("C00074", 0.00, -0.08, -0.22, -0.46, -0.62),
        ("C00022", 0.00, 0.19, 0.55, 0.89, 1.08),
        ("C00186", 0.00, 0.21, 0.68, 1.24, 1.55),
        ("C00024", 0.00, 0.12, 0.38, 0.64, 0.79),
        ("C00158", 0.00, 0.14, 0.44, 0.77, 0.94),
        ("C00020", 0.00, -0.58, -0.77, -0.83, -0.81),
        ("C00004", 0.00, 0.28, 0.51, 0.62, 0.66),
        ("C00005", 0.00, 0.16, 0.38, 0.52, 0.58),
    ]
    pd.DataFrame(
        rows, columns=["id", "0", "5", "15", "30", "60"]
    ).to_csv(os.path.join(OUTPUT_DIR, "metabolome_timecourse.csv"), index=False)


#: Which pathway each example reaction belongs to. The reversible steps shared
#: by glycolysis and gluconeogenesis are listed under both.
PATHWAYS = {
    "Glycolysis": ["R00299", "R00771", "R00756", "R01070", "R01015", "R01061",
                   "R01512", "R01518", "R00658", "R00200"],
    "Gluconeogenesis": ["R00341", "R00762", "R00303", "R00771", "R01070", "R01015",
                        "R01061", "R01512", "R01518", "R00658"],
    "Pyruvate fate": ["R00703", "R00209"],
    "Pentose phosphate pathway": ["R00835"],
}


def _write_pathways():
    rows = [{"reaction": reaction, "pathway": pathway}
            for pathway, reactions in PATHWAYS.items() for reaction in reactions]
    pd.DataFrame(rows).to_csv(os.path.join(OUTPUT_DIR, "reaction_pathways.csv"), index=False)


if __name__ == "__main__":
    main()
