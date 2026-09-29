TransNet
========

**Trans-omics network reconstruction and analysis.**

TransNet builds a *trans-omic network*, a typed, directed, signed regulatory
hierarchy connecting omic layers, and analyses it as a network:

.. code-block:: text

    signal  ->  TF  ->  gene  ->  enzyme protein  ->  REACTION  <-  metabolite

Every edge records the relationship it represents, the direction of its
regulatory effect, and the evidence behind it. That is what lets TransNet answer
the questions such as: *is this reaction regulated through gene
expression or through its allosteric effectors, and do the two agree?*

Multi-omics integration asks which molecules co-vary; trans-omics asks
how a signal propagates through a biochemical network. Here we aim to
answer the questions that arise from the interplay between different
omic layers. See :ref:`transomics-analyses` for the analyses that
follow from that difference.


Quickstart
----------

.. code-block:: bash

   pip install -e .

.. code-block:: python

   from transnet import (
       load_example_network, load_example_omics,
       map_omics_to_network, reaction_regulation_table,
   )

   graph = load_example_network()
   tables = load_example_omics()

   report = map_omics_to_network(
       graph, tables,
       id_column="id", log2fc_column="log2FC", qvalue_column="padj",
   )
   print(report)          # how many features matched, per layer

   table = reaction_regulation_table(graph)
   table[table["controversial"]]     # reactions whose two axes disagree

Then you can work your way through the examples:

.. code-block:: bash

   python notebooks/walkthroughs/build_network.py
   python notebooks/walkthroughs/reaction_regulation.py

Contents
--------

Where to find what
------------------

* **The analysis catalogue** explains each analysis, where it comes from, and
  what it needs.
* **The walkthroughs** run offline (without accessing any external databases)
  on bundled differential tables.
* **The studies** run the whole catalogue on real data. Each has a summary page
  and a full notebook.

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
   notebooks/walkthroughs/transcription_factors
   notebooks/walkthroughs/external_annotation

.. toctree::
   :maxdepth: 2
   :caption: Studies

   brown_adipocytes
   motrpac_study
   published_study
   kokaji_study

.. toctree::
   :maxdepth: 1
   :caption: Study notebooks, with every output

   studies/brown_adipocytes
   studies/motrpac_rat
   studies/obese_liver
   studies/kokaji_liver

.. toctree::
   :maxdepth: 2

   api
   references
   author


Citing the methods
------------------

A TransNet manuscript is in preparation. Each analysis in
:ref:`transomics-analyses` names its source. The framework as a
whole is largely based on concepts described in
:ref:`Yugi et al. 2016 <ref-yugi2016>`; :doc:`references` is the full
bibliography.


Indices
-------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
