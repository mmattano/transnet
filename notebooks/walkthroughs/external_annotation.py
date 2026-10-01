# %% [markdown]
# # External annotation
#
# The network is assembled from databases, and the same clients that build it can
# be called on a result. This walkthrough shows each one on the responsive
# subnetwork of the bundled example: Reactome for pathway context, HMDB for the
# metabolites, BRENDA for the kinetics behind the allosteric edges, ChIP-Atlas
# for the transcription-factor edges, and KEGG for the signalling layer.
#
# **This notebook needs network access**, which is why `make notebooks` leaves it
# out and the documentation renders it without executing. BRENDA also needs
# credentials in `BRENDA_EMAIL` and `BRENDA_PASSWORD`; the cells that use it are
# skipped when those are absent.

# %%
import os

import pandas as pd

from transnet import (
    load_example_network,
    load_example_omics,
    map_omics_to_network,
    regulated_nodes,
    responsive_subnetwork,
)

graph = load_example_network()
map_omics_to_network(graph, load_example_omics(), id_column="id",
                     log2fc_column="log2FC", qvalue_column="padj")
responsive = responsive_subnetwork(graph)

genes = [n for n in regulated_nodes(graph)
         if graph.nodes[n].get("layer") == "Transcriptome"]
metabolites = [n for n in regulated_nodes(graph)
               if graph.nodes[n].get("layer") == "Metabolome"]
print(f"{len(genes)} changed genes, {len(metabolites)} changed metabolites")

# %% [markdown]
# ## Reactome: what pathways a result sits in
#
# Pathway over-representation is the standard first question about a gene list,
# and it is the comparison trans-omics is measured against: it names pathways,
# not reactions, and it cannot say through which axis a pathway was regulated.
# Running it here makes that difference concrete rather than asserted.
#
# The example is mouse data; Reactome's species id for mouse is 10090.

# %%
from transnet.api import reactome_map_ids_to_pathways, reactome_over_representation

pathways = reactome_over_representation(genes, species="10090")
pathways.head(10)

# %% [markdown]
# `reactome_map_ids_to_pathways` answers the reverse question: which pathways
# contain a given molecule. It is useful for putting a name to a single
# molecule that an analysis singled out, such as a cross-layer hub. It takes
# UniProt accessions, so it is given the Proteome nodes rather than the
# Entrez-keyed Transcriptome ones.

# %%
proteins = [n for n in regulated_nodes(graph)
            if graph.nodes[n].get("layer") == "Proteome"]
reactome_map_ids_to_pathways(proteins[:5], species="10090").head(10)

# %% [markdown]
# `reactome_get_pathways` lists an organism's top-level pathways and
# `reactome_get_pathway_reactions` the reactions inside one, which is how a
# pathway name becomes a set of reaction identifiers to match against the
# network's own Reactions layer.

# %%
from transnet.api import reactome_get_pathway_reactions, reactome_get_pathways

catalogue = reactome_get_pathways(species="10090")
print(f"{len(catalogue)} top-level mouse pathways")
catalogue.head(5)

# %%
if not catalogue.empty:
    inside = reactome_get_pathway_reactions(catalogue.iloc[0]["stId"])
    print(f"{len(inside)} reactions in {catalogue.iloc[0]['displayName']}")
    inside.head(5)

# %% [markdown]
# `reactome_get_pathway_entities` lists the molecules in a pathway (proteins,
# small molecules and complexes), and `reactome_get_pathway_hierarchy` returns
# the tree that nests pathways inside larger ones.

# %%
from transnet.api import reactome_get_pathway_entities, reactome_get_pathway_hierarchy

if not catalogue.empty:
    entities = reactome_get_pathway_entities(catalogue.iloc[0]["stId"])
    print(f"{len(entities)} molecules in {catalogue.iloc[0]['displayName']}")
    entities.head(5)

# %%
tree = reactome_get_pathway_hierarchy(species="10090")
print(f"{len(tree)} top-level branches in the mouse pathway hierarchy")

# %% [markdown]
# ## HMDB: the metabolite layer
#
# Metabolomics reports names, and names are ambiguous. `hmdb_search_metabolites`
# resolves one, and `hmdb_get_diseases` gives the clinical associations HMDB
# records for it.
#
# hmdb.ca currently answers programmatic XML queries with HTTP 403, so these
# calls return empty frames here. They log the refusal and do not raise, which is
# what keeps a network build from dying on one blocked database. If HMDB
# annotation matters to your work, download their release and read it locally.

# %%
from transnet.api import hmdb_get_diseases, hmdb_search_metabolites

found = hmdb_search_metabolites("pyruvate", max_results=5)
found[["hmdb_id", "name", "formula"]] if not found.empty else found

# %%
if not found.empty:
    diseases = hmdb_get_diseases(found.iloc[0]["hmdb_id"])
    print(f"{len(diseases)} disease associations for {found.iloc[0]['name']}")
    diseases.head(8)

# %%
from transnet.api.hmdb import hmdb_get_metabolite, hmdb_get_pathways

if not found.empty:
    record = hmdb_get_metabolite(found.iloc[0]["hmdb_id"])
    print(f"{len(record)} fields in the HMDB record")
    hmdb_get_pathways(found.iloc[0]["hmdb_id"]).head(5)

# %%
from transnet.api import hmdb_map_ids

hmdb_map_ids(["C00022", "C00031"], from_type="kegg_id")

# %% [markdown]
# `hmdb_enrich_metabolites` works on `Metabolite` objects rather than on the
# graph, filling in formula, InChI, ChEBI and PubChem identifiers in place. This
# is how the builder gives the metabolome layer the cross-references that later
# let an omics table keyed by PubChem be mapped on.

# %%
from transnet.api import hmdb_enrich_metabolites
from transnet.biology.elements import Metabolite

pyruvate = Metabolite()
pyruvate.kegg_compound_id = "C00022"
pyruvate.kegg_name = "Pyruvate"

hmdb_enrich_metabolites([pyruvate], fields=["chebi_id", "pubchem_id", "inchikey"])
{"chebi": pyruvate.chebi_id, "pubchem": pyruvate.pubchem_id,
 "inchikey": pyruvate.inchikey}

# %% [markdown]
# ## BRENDA: where the allosteric edges come from
#
# `allosteric_activation` and `allosteric_inhibition` edges are BRENDA records.
# Reading them back for one enzyme shows what the network is asserting, and the
# kinetic constants add something the graph does not hold: whether a measured
# concentration is anywhere near the one that matters.
#
# A Ki beside a measured metabolite change turns "this metabolite inhibits this
# enzyme" into "this metabolite inhibits this enzyme at concentrations like the
# ones measured", which is the difference between an annotation and a mechanism.

# %%
have_credentials = bool(os.environ.get("BRENDA_EMAIL")
                        and os.environ.get("BRENDA_PASSWORD"))
print("BRENDA credentials present" if have_credentials
      else "no BRENDA credentials: the next cells are skipped")

# %%
if have_credentials:
    from transnet.api import (
        brenda_get_activators,
        brenda_get_inhibitors,
        brenda_get_ki_values,
        brenda_get_km_values,
    )

    EC = "2.7.1.40"        # pyruvate kinase, the example's controversial reaction

    inhibitors = brenda_get_inhibitors(EC, organism="Mus musculus")
    activators = brenda_get_activators(EC, organism="Mus musculus")
    print(f"EC {EC}: {len(inhibitors)} inhibitors, {len(activators)} activators "
          f"recorded for mouse")
    inhibitors.head(8)

# %%
if have_credentials:
    from transnet.api import brenda_get_kcat_values

    km = brenda_get_km_values(EC, organism="Mus musculus")
    ki = brenda_get_ki_values(EC, organism="Mus musculus")
    kcat = brenda_get_kcat_values(EC, organism="Mus musculus")
    print(f"{len(km)} Km, {len(ki)} Ki and {len(kcat)} kcat values")
    pd.concat([km.head(4).assign(constant="Km"), ki.head(4).assign(constant="Ki"),
               kcat.head(4).assign(constant="kcat")])

# %% [markdown]
# BRENDA also records what a reaction consumes and produces and what it needs as
# a cofactor. The network takes its substrates and products from KEGG reaction
# equations, so these are a cross-check rather than a source, and disagreements
# between the two are worth looking at before trusting a reaction.

# %%
if have_credentials:
    from transnet.api import (
        brenda_get_cofactors,
        brenda_get_products,
        brenda_get_substrates,
    )

    for label, fetch in (("substrates", brenda_get_substrates),
                         ("products", brenda_get_products),
                         ("cofactors", brenda_get_cofactors)):
        frame = fetch(EC, organism="Mus musculus")
        print(f"{label:<11} {len(frame)}")

# %% [markdown]
# `brenda_enrich_proteins` works on `Protein` objects, as the network builder
# does: it fetches the kinetic data for each protein's EC numbers and stores it
# on the object. `BrendaClient` is the authenticated connection all of these
# functions share; use it directly for a BRENDA query the helpers do not cover.

# %%
if have_credentials:
    from transnet import Protein
    from transnet.api import brenda_enrich_proteins

    kinase = Protein(uniprot_id="P53657", name="Pklr", ec_number=[EC])
    brenda_enrich_proteins([kinase], organism="Mus musculus")
    print({field: len(getattr(kinase, field, None) or [])
           for field in ("inhibitors", "activators", "substrates", "products", "cofactors")})

# %% [markdown]
# ## ChIP-Atlas: the transcription-factor edges
#
# `list_chip_tfs` is the set of factors with data for a genome, and
# `get_chip_tf_targets` returns the target genes for the ones you ask about,
# with the binding score that becomes the edge's `confidence`. `min_score`
# controls how much of the transcriptome each factor appears to own, which is the
# limit `transcription_factors.py` discusses.

# %%
from transnet.api.chip_atlas import get_chip_tf_targets, list_chip_tfs

available = list_chip_tfs(genome="mm10")
print(f"{len(available)} factors with mm10 data; Foxo1 present: "
      f"{'Foxo1' in available}")

# %%
targets = get_chip_tf_targets(["Foxo1"], genome="mm10", distance=5, min_score=100)
print(f"{len(targets)} Foxo1 target genes at a binding score of 100 or more")
targets.head(8)

# %% [markdown]
# `get_ChIP_data` downloads the raw target table for a factor, one column per
# ChIP-seq experiment, which `get_chip_tf_targets` summarises into one score per
# gene. It returns the table, the experiments per factor, and the factors whose
# download failed.

# %%
from transnet.api.chip_atlas import get_ChIP_data

raw, experiments, failed = get_ChIP_data(genome="mm10", distance=5, proteins=["Foxo1"])
print(f"{raw.shape[0]} genes, {raw.shape[1]} columns; "
      f"{len(experiments.get('Foxo1', []))} Foxo1 experiments")
raw.iloc[:5, :6]

# %% [markdown]
# ChIP-Atlas also records the cell type of every experiment.
# `list_ChIP_cell_type_classes` and `list_ChIP_cell_types` list them, and
# `get_ChIP_exps` lists the experiments, so the evidence for a factor's targets
# can be restricted to a relevant tissue. All three read ChIP-Atlas's full
# experiment list. It is a 200 MB download, stored in `~/.cache/transnet` and
# reused afterwards, so the cell below runs only when the environment variable
# `TRANSNET_LARGE_DOWNLOADS` is set.

# %%
if os.environ.get("TRANSNET_LARGE_DOWNLOADS"):
    from transnet.api.chip_atlas import (
        get_ChIP_exps,
        list_ChIP_cell_type_classes,
        list_ChIP_cell_types,
    )

    classes = list_ChIP_cell_type_classes(genome="mm10")
    print(f"{len(classes)} cell-type classes, e.g. {', '.join(sorted(classes)[:5])}")
    print(f"{len(list_ChIP_cell_types(genome='mm10', cell_type_class='Liver'))} "
          f"liver cell types")
    get_ChIP_exps(genome="mm10", cell_type_class="Liver").head(5)
else:
    print("set TRANSNET_LARGE_DOWNLOADS=1 to list ChIP-Atlas cell types and experiments")

# %% [markdown]
# ## KEGG: the signalling layer
#
# `kegg_signaling_relations` parses the KGML of an organism's signal-transduction
# pathways into kinase-substrate relations. These become the `phosphorylation`
# and `kinase_tf` edges of the Signaling layer, which
# `maintenance/build_networks.py --signaling` adds to a network and the
# phosphoproteomics studies map sites onto.
#
# KEGG's own subtype annotation supplies the sign: `activation` is +1,
# `inhibition` -1, and a plain `phosphorylation` is 0, since whether a site
# activates or inhibits is not recorded.

# %%
from transnet.api import kegg_signaling_relations

relations = kegg_signaling_relations("mmu", pathway_ids=["mmu04910"])
print(f"{len(relations)} relations from the insulin signalling pathway")
relations.head(10)

# %%
relations["subtype"].value_counts().to_frame("relations")
