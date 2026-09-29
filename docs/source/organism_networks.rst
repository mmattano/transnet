.. _organism-networks:

Building a network for your organism
====================================

Every study in this documentation maps its measurements onto a **reference
network**: an organism-wide, typed, directed, signed graph assembled from public
databases, built once and reused. The study supplies the data; the network
supplies the wiring. This page describes how that network is made and how to
make one for an organism TransNet does not already carry.

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
assembly for rat, and only 52 factors on it. That single line is why the rat
network's transcription-factor coverage is thin, and why the MoTrPAC study can
test only 32 factors and finds none in four of its six tissues. A configuration
choice made once, for a reason external to TransNet, shows up as a limit on a
result four pages away. It is worth knowing which of your organism's numbers are
like this.

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

Assembling a network in your own code
-------------------------------------

The builder is a script around the layer objects, and those can be used
directly when you want a network that is not organism-wide:

.. code-block:: python

    from transnet import Transnet
    from transnet.biology.layers import Proteome, Metabolome

    proteome = Proteome()
    proteome.populate(kegg_organism="mmu", ncbi_organism="10090")
    proteome.get_interaction_partners()          # STRING
    proteome.get_transcription_factor_targets()  # ChIP-Atlas
    proteome.get_brenda_kinetics()               # BRENDA effectors

    metabolome = Metabolome()
    metabolome.populate(kegg_organism="mmu", ncbi_organism="10090")

    network = Transnet(proteome=proteome, metabolome=metabolome,
                       reactions=reactions)
    graph = network.generate_graph()

``notebooks/walkthroughs/build_network.py`` shows the same object API on a
network small enough to read.

.. seealso::

   :doc:`network_model` for what the resulting graph contains, including
   :ref:`brenda-effectors` and :ref:`chip-atlas-threshold`, the two places where
   a database's conventions most affect the network.
