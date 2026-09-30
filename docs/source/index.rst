TransNet
========

**Trans-omics network reconstruction and analysis.**

TransNet builds a *trans-omic network*: a network that connects molecules
from different omic layers through the regulatory relationships between
them, and analyses it.

.. code-block:: text

    signal  ->  TF  ->  gene  ->  enzyme protein  ->  REACTION  <-  metabolite

Every edge records which relationship it represents (for example catalysis or
allosteric inhibition), whether it increases or decreases its target, and which
database it comes from. This lets TransNet answer questions that a separate
analysis of each layer cannot, such as: *is this reaction regulated through
the amount of its enzyme, or through the metabolites acting on it, and do the
two agree?*

Multi-omics integration usually asks which molecules change together across
data types. Trans-omics asks how a change travels through the biochemical
network from one layer to the next. :ref:`transomics-analyses` describes the
analyses that follow from this.


Quickstart
----------

.. code-block:: bash

   pip install transnet        # or, from a clone: pip install -e .

.. code-block:: python

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
   print(report)          # how many measurements matched a node, per layer

   table = reaction_regulation_table(graph)
   table[table["controversial"]]     # reactions whose two axes disagree

The walkthroughs continue from here, starting with
:doc:`notebooks/walkthroughs/build_network`.


Where to find what
------------------

* **The network model** and **organism networks** describe what a TransNet
  network contains and how the organism-wide networks are built and updated.
* **The analysis catalogue** explains each analysis, what it needs, and where
  it comes from.
* **The walkthroughs** introduce the functions one topic at a time. They run
  offline, in seconds, on a small bundled example of mouse liver glucose
  metabolism.
* **The studies** apply the whole catalogue to published data. Each has a
  summary page and a full notebook with every table and figure.

.. toctree::
   :maxdepth: 2

   network_model
   organism_networks
   transomics_analyses

.. toctree::
   :maxdepth: 1
   :caption: Walkthroughs

   notebooks/walkthroughs/build_network
   notebooks/walkthroughs/responsive_network
   notebooks/walkthroughs/reaction_regulation
   notebooks/walkthroughs/regulatory_paths
   notebooks/walkthroughs/temporal_and_hubs
   notebooks/walkthroughs/compare_conditions
   notebooks/walkthroughs/network_topology
   notebooks/walkthroughs/export_network
   notebooks/walkthroughs/transcription_factors
   notebooks/walkthroughs/external_annotation

The four studies, and what each one shows:

* :doc:`brown_adipocytes`: mouse brown fat cells stimulated with
  norepinephrine (Anagho-Mattanovich *et al.* 2025), three layers and up to seven time
  points; the full catalogue on one cell type.
* :doc:`motrpac_study`: endurance training in six rat tissues (MoTrPAC 2024),
  with phospho-, acetyl- and ubiquitin-proteomics; tissues compared on one
  network.
* :doc:`obese_liver_panel`: lean and obese mouse liver after a glucose load,
  measured on a 19-enzyme panel (Uematsu *et al.* 2022); published claims
  checked one by one.
* :doc:`liver_timecourse`: the same comparison genome-wide, over four hours
  (Kokaji *et al.* 2020); two genotypes compared as networks, and inferences
  checked against the authors' own.

.. toctree::
   :maxdepth: 2
   :caption: Studies

   brown_adipocytes
   motrpac_study
   obese_liver_panel
   liver_timecourse

.. toctree::
   :maxdepth: 1
   :caption: Study notebooks, with every output

   studies/brown_adipocytes
   studies/motrpac_rat
   studies/obese_liver
   studies/liver_timecourse

.. toctree::
   :maxdepth: 2

   api
   references
   author


Citing the methods
------------------

A TransNet manuscript is in preparation. Each analysis in
:ref:`transomics-analyses` names the publication it comes from. The framework
as a whole builds on :ref:`Yugi et al. 2016 <ref-yugi2016>`, and
:doc:`references` is the full bibliography.


Indices
-------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
