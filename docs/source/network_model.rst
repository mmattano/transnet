.. _network-model:

The trans-omic network model
============================

TransNet's network is a :class:`networkx.MultiDiGraph`. Both properties are
load-bearing:

**Directed**, because direction *is* the biology. A substrate flows into a
reaction and a product flows out; a transcription factor regulates its target
and not the reverse. An undirected graph cannot express any of that, and every
analysis that follows regulatory flow becomes impossible.

**Multigraph**, because one pair of molecules can stand in more than one
relationship. Glucose-6-phosphate is the *product* of hexokinase and its
*allosteric inhibitor*. Those are two different edges with opposite signs; a
simple graph collapses them into one edge that means nothing.

.. code-block:: python

    >>> graph = load_example_network()
    >>> graph.get_edge_data("C00668", "R00299")     # G6P -> hexokinase
    {'allosteric_inhibition': {..., 'sign': -1, ...}}
    >>> graph.get_edge_data("R00299", "C00668")     # hexokinase -> G6P
    {'product': {..., 'sign': 1, ...}}

Most centrality and community-detection implementations need a simple
undirected graph, and :func:`transnet.to_simple_graph`
projects one, collapsing parallel edges and recording what was collapsed in an
``edge_types`` list.


Layers
------

.. list-table::
   :header-rows: 1
   :widths: 18 18 64

   * - Layer
     - Node identifier
     - Contents
   * - ``Signaling``
     - UniProt (optionally ``_site``)
     - Kinases, phosphatases, phosphosites.
   * - ``Proteome``
     - UniProt accession
     - Proteins, including enzymes and transcription factors
   * - ``Transcriptome``
     - NCBI/Entrez gene ID
     - Transcripts
   * - ``Reactions``
     - KEGG reaction ID (``R#####``)
     - Metabolic reactions; the convergence point of the hierarchy
   * - ``Metabolome``
     - KEGG compound ID (``C#####``)
     - Metabolites
   * - ``Pathways``
     - KEGG / Reactome ID
     - Pathway annotation

Layers are ordered by :data:`transnet.LAYER_HIERARCHY`, and
:func:`transnet.available_layers` reports which are actually present in a given
graph. Nothing in the package assumes a fixed set; see
:ref:`transomics-analyses`.


Edge types
----------

The edge-type vocabulary is fixed and defined once, in
:data:`transnet.EDGE_TYPES`. Each entry carries the layers it connects, its
default sign, and the database it comes from. The taxonomy implements the five
connection approaches of :ref:`Yugi et al. 2016 <ref-yugi2016>`.

.. list-table::
   :header-rows: 1
   :widths: 26 24 8 42

   * - ``edge_type``
     - Direction
     - Sign
     - Source
   * - ``phosphorylation``
     - Signaling → Proteome
     - ±
     - KEGG / user phosphoproteomics
   * - ``kinase_tf``
     - Signaling → Proteome (TF)
     - ±
     - KEGG signaling pathways
   * - ``transcriptional_regulation``
     - Proteome (TF) → Transcriptome
     - 0
     - ChIP-Atlas
   * - ``translation``
     - Transcriptome → Proteome
     - +1
     - NCBI/UniProt identifier join
   * - ``protein_interaction``
     - Proteome ↔ Proteome
     - 0
     - STRING (score in ``confidence``)
   * - ``catalysis``
     - Proteome → Reactions
     - +1
     - KEGG, joined on EC number
   * - ``gene_catalysis``
     - Transcriptome → Reactions
     - +1
     - KEGG EC; used when there is no Proteome layer
   * - ``substrate``
     - Metabolome → Reactions
     - +1
     - KEGG reaction equation
   * - ``product``
     - Reactions → Metabolome
     - +1
     - KEGG reaction equation
   * - ``allosteric_activation``
     - Metabolome → Reactions
     - +1
     - BRENDA
   * - ``allosteric_inhibition``
     - Metabolome → Reactions
     - −1
     - BRENDA
   * - ``enzymatic``
     - Proteome → Metabolome
     - 0
     - KEGG; shortcut used only without a Reactions layer

Sign ``0`` means **unknown**, not "no effect". For example, ChIP-Atlas evidence
says a factor binds a promoter, not whether it activates or represses it.
Analyses treat unsigned steps as sign-preserving but count them: a path with
``unsigned_steps > 0`` predicts a direction only tentatively, and
:func:`~transnet.trace_regulatory_paths` reports that column so the distinction
is never lost.


Node attributes
---------------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Key
     - Meaning
   * - ``layer``
     - Which layer the node belongs to
   * - ``node_type``
     - ``Gene``, ``Protein``, ``Metabolite``, ``Reaction``, ...
   * - ``name``
     - Display name (gene symbol, compound name, reaction name)
   * - ``reversible``
     - Reactions only; parsed from the KEGG equation arrow
   * - ``ec``
     - EC number(s), on enzymes and reactions
   * - ``value``, ``log2fc``, ``qvalue``
     - Written by :func:`~transnet.map_omics_to_network`
   * - ``regulated``
     - ``+1`` / ``-1`` / ``0``; what the analyses read
   * - ``t_half``, ``ec50``
     - Written by :func:`~transnet.assign_temporal_parameters`


Edge attributes
---------------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Key
     - Meaning
   * - ``edge_type``
     - One of the fixed vocabulary above
   * - ``sign``
     - ``+1`` activating, ``-1`` inhibiting, ``0`` unknown
   * - ``role``
     - Finer-grained role, defaulting to ``edge_type``
   * - ``directed``
     - Whether the relationship is inherently directional
   * - ``weight``
     - Edge weight
   * - ``stoichiometry``
     - Coefficient, on substrate and product edges
   * - ``ec``
     - EC number that produced the edge
   * - ``source_db``, ``evidence``
     - Provenance
   * - ``confidence``
     - Source score where one exists (e.g. STRING)


Reading and writing
-------------------

A saved network is ``interactions.csv`` plus ``nodes.csv``:

.. code-block:: python

    from transnet.io import read_network, write_network

    graph = read_network("data/example")
    write_network(subnetwork, "results/responsive")

:func:`~transnet.io.read_network` also reads networks written before the typed
schema, filling the missing fields with neutral defaults and warning that edge
signs and roles need a rebuild to recover.


Building from databases
-----------------------

A network is assembled from KEGG, UniProt, Ensembl, STRING, ChIP-Atlas and
optionally BRENDA, one organism at a time, and saved as the two CSVs above. What
each database contributes, how an organism is configured, and how to build one
for an organism TransNet does not already carry are covered in
:doc:`organism_networks`.

The two sections below are the places where a database's own conventions most
affect what ends up in the network, whichever organism it is built for.

.. _brenda-effectors:

BRENDA effectors: things to keep in mind
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BRENDA reports activators and inhibitors by *name*; the network needs KEGG
compound ids. Two rules translate one to the other, and neither alters the
data:

* **Names are matched case- and whitespace-insensitively against every KEGG
  synonym**, and a small set of spelling equivalences is tried when the name as
  given does not match: three-letter amino acid codes (``L-Ala`` = L-alanine),
  ``diphosphate`` = ``bisphosphate``, and acid = anion (``citric acid`` =
  ``citrate``). These are nomenclature, the same compound written two ways.
* **Every resolved effector is kept.** BRENDA lists any compound shown to change
  an enzyme's activity in vitro, including drugs, flavonoids and assay reagents.
  ``allosteric_interactions(endogenous_only=True)`` removes effectors that are
  not a substrate or product of a reaction the organism catalyses, but it is off
  by default because the rule misclassifies: on the mouse network it removes 542
  of 888 resolved compounds, most correctly (luteolin, EDTA, cyclosporin A,
  Triton X-100), but also physiological ions (Zn\ :sup:`2+`, Ca\ :sup:`2+`,
  Mg\ :sup:`2+`), heparin and deoxycholate, an endogenous bile acid. It barely
  changes results in any case: the metabolite axis and metabolite regulatory roles only use *measured*
  metabolites, and assay reagents are usually not measured.

Names BRENDA gives that match no KEGG compound -- mostly synthetic inhibitors --
are logged and skipped.

.. _chip-atlas-threshold:

ChIP-Atlas binding scores
~~~~~~~~~~~~~~~~~~~~~~~~~

ChIP-Atlas lists every gene showing *any* detectable binding for a factor, so
accepting them all gives roughly 15,000 targets per factor. In a mouse network
that is 4.3 million ``transcriptional_regulation`` edges -- 91% of the graph --
and every cross-layer analysis is then dominated by one edge type.

:meth:`~transnet.biology.layers.Proteome.get_transcription_factor_targets`
therefore scores each candidate target by its **mean binding score across the
factor's experiments** and keeps those at or above ``score_threshold``
(default 100). The mean is ChIP-Atlas's own consensus score -- it equals the
``{TF}|Average`` column when no cell-type restriction is applied -- so a
single strong peak in 1 of 193 experiments is downweighted rather than
promoted to a target. The score is a MACS2 peak significance derived from
:math:`-\log_{10}(Q)`, so larger means stronger and more significant binding.
Across mouse factors the median listed gene scores 10-60, and the default
threshold keeps roughly the top 15-20% per factor: for the whole mouse
proteome it yields 1.55 million target links across 716 factors, against
roughly 11 million with no threshold.

.. code-block:: python

    # permissive: every gene ChIP-Atlas lists
    proteome.get_transcription_factor_targets(score_threshold=0)

    # strict, and capped for architectural factors such as CTCF that
    # genuinely bind most promoters
    proteome.get_transcription_factor_targets(
        score_threshold=250, max_targets_per_tf=1000,
    )

The retained score is written onto each ``transcriptional_regulation`` edge as
``confidence``, so downstream filtering does not require a rebuild.
``maintenance/build_networks.py`` exposes the same knob as
``--chip-score-threshold``.

``generate_graph(require_enzyme=True)`` (the default) restricts the reaction
layer to reactions with a catalysing enzyme present in the organism, which keeps
the network organism-specific rather than importing all of KEGG.
