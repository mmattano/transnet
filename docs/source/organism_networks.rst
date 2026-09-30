.. _organism-networks:

Building a network for your organism
====================================

Every study in this documentation maps its measurements onto a **reference
network**: an organism-wide network assembled from public databases, built once
and reused. The study supplies the data; the network supplies the connections.
This page describes how the reference networks are built and updated, and how to
build one for an organism TransNet does not yet include.

Five are configured and built: human, mouse, rat, yeast and E. coli.

What a reference network is assembled from
------------------------------------------

Each database contributes one part of the hierarchy, and the builder asks each
one separately.

.. list-table::
   :header-rows: 1
   :widths: 22 30 48

   * - Database
     - Contributes
     - Becomes
   * - KEGG
     - Reactions, compounds, EC numbers, pathways
     - The Reactions and Metabolome layers, and ``substrate`` / ``product`` edges
   * - UniProt
     - Reviewed proteins with EC annotation
     - The Proteome layer and ``catalysis`` edges
   * - Ensembl
     - Transcript models and identifiers
     - The Transcriptome layer and ``translation`` edges
   * - STRING
     - Protein-protein interactions
     - ``protein_interaction`` edges
   * - ChIP-Atlas
     - Transcription factor binding
     - ``transcriptional_regulation`` edges
   * - BRENDA
     - Allosteric activators and inhibitors
     - ``allosteric_activation`` / ``allosteric_inhibition`` edges
   * - KEGG KGML
     - Signal-transduction relations
     - The optional Signaling layer

BRENDA and the Signaling layer are optional and off by default, because both are
slow and neither is needed by every analysis. What they add is described in
:ref:`brenda-effectors` and :ref:`signal-flow`.

The organism registry
---------------------

The databases name organisms differently: KEGG uses a three-letter code, UniProt
a taxon id that is sometimes a strain rather than a species, Ensembl a lowercase
species name, ChIP-Atlas a genome assembly. An organism is therefore defined by
one entry in :data:`transnet.organisms.ORGANISMS` saying what to ask each
database for.

.. list-table::
   :header-rows: 1
   :widths: 26 40 34

   * - Key
     - What it is
     - Where to look it up
   * - ``kegg_org``
     - Three-letter KEGG organism code
     - :func:`transnet.api.kegg.kegg_list_organisms`
   * - ``organism_full``
     - Species name as BRENDA writes it
     - brenda-enzymes.org
   * - ``ncbi_org``
     - NCBI taxonomy id, used for STRING
     - NCBI Taxonomy
   * - ``uniprot_org``
     - UniProt taxon, where reviewed entries sit under a strain
     - uniprot.org (optional)
   * - ``ensembl_org``, ``ensembl_release``
     - Ensembl species name and release
     - ensembl.org
   * - ``genome_chip``
     - ChIP-Atlas genome assembly
     - chip-atlas.org
   * - ``transcriptome_source``
     - ``"kegg"`` where Ensembl has no data for the organism
     - optional
   * - ``expected_absent``
     - Edge types no database can supply for this organism
     - optional

Rat, end to end
---------------

Rat is a good example because it is ordinary in every respect but one.

Its registry entry:

.. code-block:: python

    "rat": {
        "kegg_org": "rno",
        "organism_full": "Rattus norvegicus",
        "ncbi_org": "10116",
        "ensembl_org": "rattus_norvegicus",
        "ensembl_release": 109,
        # rn6 is the only rat assembly ChIP-Atlas carries (52 factors).
        "genome_chip": "rn6",
    },

Building it:

.. code-block:: bash

    python maintenance/build_networks.py --organisms rat --brenda

That writes ``data/rat/<date>/`` with ``nodes.csv``, ``interactions.csv`` and a
summary, and points ``data/rat/latest`` at it. The result holds **74,363
molecules and 258,015 typed relationships, 55% of which cross between layers**,
over 12,384 reactions, 8,233 proteins, 19,656 metabolites and 54,991 genes. This
is the network :doc:`motrpac_study` maps six tissues onto.

The one quirk is the comment in the entry. ChIP-Atlas carries only the ``rn6``
assembly for rat, with data for only 52 factors. This is why the rat network has
few transcription-factor edges, and why the MoTrPAC study can test only 32
factors and finds none in four of its six tissues. Limits like this come from
the databases, not from the analysis, so it is worth checking which of your
organism's results they affect.

Adding another organism
-----------------------

Fill in one entry by the same pattern and run the same command. From your own
code, without editing anything in the repository:

.. code-block:: python

    from transnet.organisms import register_organism

    register_organism(
        "zebrafish",
        kegg_org="dre",
        organism_full="Danio rerio",
        ncbi_org="7955",
        ensembl_org="danio_rerio",
        ensembl_release=109,
        genome_chip="danRer11",
    )

Then ``--organisms zebrafish`` builds it. Every analysis in
:ref:`the catalogue <transomics-analyses>` works on the result without
modification, and any layer or edge type the databases could not supply is
simply absent rather than fatal.

Two keys exist for organisms that do not fit the common case.
``transcriptome_source="kegg"`` builds the transcriptome from KEGG genes where
Ensembl has no data, which is how E. coli is handled. ``expected_absent`` lists
edge types no database supplies for the organism, so their absence is reported
as expected rather than as a gap: E. coli declares
``transcriptional_regulation``, because ChIP-Atlas carries no bacterial genome.

.. _build-report:

Reading the build report
------------------------

A build that finishes is not necessarily a build that worked. An unreachable
database leaves a layer empty, and the network is then quietly missing a whole
class of edge. The builder therefore ends by listing what it made and naming
anything missing:

.. code-block:: text

    --- rat build report ---
      74,363 nodes, 258,015 edges
      layers: ['Transcriptome', 'Proteome', 'Reactions', 'Metabolome']
        translation                        54,991
        catalysis                          32,899
        ...

Anything absent that was expected is reported as a gap, and the exit status is
non-zero, so a partial build cannot be mistaken for a complete one. One check
worth knowing about: a reviewed proteome of fewer than a thousand proteins is
flagged, because it almost always means the UniProt taxon is wrong rather than
that the organism is small. Yeast once built with 43 proteins for exactly that
reason, which is why its entry carries a ``uniprot_org`` override.

.. _updating-networks:

Updating the networks
---------------------

The organism networks are rebuilt by hand, not on a schedule. The databases
change a few times a year, a build takes 30-90 minutes per organism, and a
rebuilt network changes the numbers in every study, so an update should be a
deliberate step that someone checks.

To update one organism:

1. Set the BRENDA credentials (free registration at brenda-enzymes.org):

   .. code-block:: bash

       export BRENDA_EMAIL='you@example.com' BRENDA_PASSWORD='...'

2. Build it. The build saves checkpoints, so an interrupted build resumes where
   it stopped:

   .. code-block:: bash

       make networks ORGANISMS=mouse
       # the same as: python maintenance/build_networks.py --organisms mouse --brenda

3. Read the build report printed at the end, and ``data/mouse/latest/summary.txt``.
   A non-zero exit status means a layer or edge type is missing (see
   :ref:`build-report` below).
4. Rerun the studies that use the organism and refresh the documentation
   figures:

   .. code-block:: bash

       make studies

5. Check the numbers quoted on the study pages against the new notebook
   output before committing.

The networks are several hundred megabytes and can be rebuilt at any time, so
they are not committed to the repository. ``data/<organism>/latest`` points to
the most recent dated build, and older builds can be deleted.

Assembling a network step by step
---------------------------------

``maintenance/build_networks.py`` is a script around the layer objects, which
can also be used directly, for example to build a network that covers only part
of an organism. Each step below calls one database and needs network access.

.. code-block:: python

    from transnet import Transnet
    from transnet.biology.layers import (
        Metabolome, Pathways, Proteome, Reactions, Signaling, Transcriptome,
    )

    # Reactions and compounds from KEGG
    reactions = Reactions()
    reactions.populate(from_api=True)          # every KEGG reaction equation

    metabolome = Metabolome()
    metabolome.populate()                      # every KEGG compound
    metabolome.enrich_with_hmdb()              # optional: HMDB cross-references

    # Proteins from UniProt, then their interactions and regulators
    proteome = Proteome()
    proteome.populate(kegg_organism="mmu", ncbi_organism="10090")
    proteome.get_interaction_partners()        # STRING protein interactions
    proteome.get_metabolites()                 # KEGG enzyme-compound links
    proteome.get_transcription_factor_targets(genome_ChIP="mm10")   # ChIP-Atlas
    proteome.get_brenda_kinetics(organism="Mus musculus")           # BRENDA effectors

    # Genes from Ensembl, or from KEGG with kegg_api=True
    transcriptome = Transcriptome()
    transcriptome.populate(kegg_organism="mmu", organism_full="mus_musculus",
                           ensembl=True)
    transcriptome.fill_gene_info()             # KEGG gene names, when built from KEGG

    # Optional layers
    pathways = Pathways()
    pathways.kegg_organism = "mmu"
    pathways.populate()                        # KEGG pathways
    pathways.fill_pathways()                   # their genes and compounds

    signaling = Signaling()
    signaling.populate(kegg_organism="mmu")    # kinase relations from KEGG KGML

    network = Transnet(name="mouse", reactions=reactions, metabolome=metabolome,
                       proteome=proteome, transcriptome=transcriptome,
                       pathways=pathways, signaling=signaling)
    graph = network.generate_graph()
    network.save_network("data/mouse/custom")

A Signaling layer can also be built from your own kinase-substrate table with
``Signaling.populate_from_table``, and Reactome pathways can replace KEGG's with
``Pathways.populate_from_reactome``. ``Pathways.over_representation_analysis``
runs Reactome's enrichment test on a gene list.

:doc:`notebooks/walkthroughs/build_network` shows the same object API on a
network small enough to read, without any database calls.

.. seealso::

   :doc:`network_model` for what the resulting graph contains, including
   :ref:`brenda-effectors` and :ref:`chip-atlas-threshold`, the two places where
   a database's conventions most affect the network.
