.. _transomics-analyses:

The trans-omics analysis catalogue
==================================

Multi-omics integration usually asks which molecules change together across
data types. Trans-omics asks how a change *travels* through the biochemical
network: which regulatory relationship carried it, in which direction, and
whether the measured changes are consistent with that mechanism.

The approach was developed largely by the Kuroda laboratory. All its analyses
work on one kind of network: layers of molecules joined by typed, directed and
signed regulatory edges, which converge on the metabolic reactions.

.. code-block:: text

    signal  ->  TF  ->  gene  ->  enzyme protein  ->  REACTION  <-  metabolite
                                                                    (substrate,
                                                                     product,
                                                                     allosteric
                                                                     regulator)

The analyses
------------

.. list-table::
   :header-rows: 1
   :widths: 30 44 26

   * - Analysis
     - What it answers
     - Function
   * - :ref:`Network reconstruction <network-reconstruction>`
     - Which molecules responded, and what joins them
     - ``responsive_subnetwork``
   * - :ref:`Reaction regulation axes <regulation-axes>`
     - Whether a reaction is driven by enzyme amount or by metabolite
       concentration, and whether the two agree
     - ``reaction_regulation_table``
   * - :ref:`Per-pathway balance <pathway-balance>`
     - Which axis carries each pathway
     - ``regulation_axis_summary``
   * - :ref:`Metabolite regulatory roles <metabolite-roles>`
     - Whether the changed metabolites regulate anything
     - ``metabolite_regulatory_roles``
   * - :ref:`Signal flow <signal-flow>`
     - Which metabolites a stimulus reaches, and with what sign
     - ``hierarchical_propagation``
   * - :ref:`Signed regulatory paths <signed-paths>`
     - Whether the sign along a chain matches the measurement
     - ``trace_regulatory_paths``
   * - :ref:`Cross-layer connectivity <cross-layer-connectivity>`
     - How much of the network crosses layers, and how well each is covered
     - ``cross_layer_connectivity``
   * - :ref:`Trans-omic hubs <transomic-hubs>`
     - Which molecules are hubs across layers rather than within one
     - ``transomic_hubs``
   * - :ref:`Response timing <response-timing>`
     - Whether well-connected molecules respond first
     - ``temporal_network_structure``
   * - :ref:`Transcription-factor inference <tf-activity>`
     - Which factors drove the responsive genes, and in which direction
     - ``transcription_factor_activity``
   * - :ref:`Transcript-protein concordance <concordance>`
     - Whether a protein change is transcriptional
     - ``expression_concordance``
   * - :ref:`Factors through the network <factors-on-the-network>`
     - Whether a factor's features are related by biochemistry or only by
       covariance
     - ``fit_factors``
   * - :ref:`Regulatory motifs <regulatory-motifs>`
     - Which signed wiring patterns the network contains
     - ``regulatory_motifs``
   * - :ref:`Structural vulnerability <structural-vulnerability>`
     - Which molecules hold the response together
     - ``structural_vulnerability``
   * - :ref:`Convergence null <convergence-null>`
     - Whether the layers converge more than chance allows
     - ``convergence_significance``
   * - :ref:`Comparing conditions <comparing-conditions>`
     - Which kinds of regulation were gained and lost
     - ``compare_transomic_networks``


Layers are optional
-------------------

**No analysis requires a particular layer to be present.** Each
function discovers what the network actually contains:

* :func:`~transnet.reaction_regulation_table` records which chain of layers
  supported each gene-axis call in ``gene_axis_evidence``;
* :func:`~transnet.trace_regulatory_paths` infers its own starting layer and
  reports the hierarchy it used;
* a table supplied for a layer the network lacks is reported in the mapping
  report, not raised.

:doc:`notebooks/walkthroughs/regulatory_paths` shows this by running the same
analysis with and without the optional Signaling layer.


Identifiers, and why the mapping report exists
----------------------------------------------

Each layer of a network uses one kind of identifier: Entrez ids for genes,
UniProt accessions for proteins, KEGG compound ids for metabolites. An omics
table uses whatever the instrument and pipeline produced, such as Ensembl gene
ids, RefSeq protein accessions, PubChem ids or metabolite names. A mismatch
raises no error; the analysis simply runs on almost no data.
:func:`~transnet.map_omics_to_network` therefore returns a report of how many
rows of each table found a node, and which did not.

Two helpers resolve the common cases without a database call:

:func:`~transnet.build_alias_id_map`
    Reads the translation off the layer objects, which usually know several
    identifiers per element even though the node carries one.

    .. code-block:: python

        id_map = {
            "Transcriptome": build_alias_id_map(net, "Transcriptome", "ensembl_id"),
            "Metabolome": build_alias_id_map(net, "Metabolome", "pubchem_id"),
        }
        map_omics_to_network(graph, tables, id_map=id_map)

:func:`~transnet.build_name_id_map`
    Matches compound *names* against every KEGG synonym on the network's nodes,
    case- and punctuation-insensitively -- for metabolomics tables that carry
    names rather than identifiers.

    .. code-block:: python

        mapping = build_name_id_map(graph, metabolomics["feature"])

A layer reported at 0 % matched is almost always caused by an identifier
mismatch, not by the biology.


.. _network-reconstruction:

Trans-omic network reconstruction
---------------------------------

:meth:`transnet.Transnet.generate_graph` assembles
the typed graph from the layer objects; :func:`~transnet.responsive_subnetwork`
extracts the part of it that responded: the molecules that changed, in every
layer, plus the regulatory edges joining them. Reaction nodes are retained as
connectors.

.. code-block:: python

    from transnet import load_example_network, map_omics_to_network
    from transnet import responsive_subnetwork

    graph = load_example_network()
    map_omics_to_network(graph, tables, log2fc_column="log2FC",
                         qvalue_column="padj")
    responsive = responsive_subnetwork(graph)

:Function: :meth:`transnet.Transnet.generate_graph`,
   :func:`~transnet.responsive_subnetwork`
:Needs: A network with omics mapped onto it. Any set of layers.
:Example: ``notebooks/walkthroughs/build_network.py``,
   ``notebooks/walkthroughs/responsive_network.py``
:Reference: :ref:`Yugi et al. 2016 <ref-yugi2016>`


.. _regulation-axes:

Reaction regulation-axis attribution
------------------------------------

The rate of a reaction can change for two reasons. On the *gene-expression
axis* (also called the enzyme axis) the amount of enzyme changes. On the
*metabolite axis* the same amount of enzyme works faster or slower, because its
substrates, products or allosteric regulators changed. The two axes are
measured in different layers, and they can point in opposite directions.
Separating them is what most clearly distinguishes trans-omics from pathway
enrichment, which cannot say how a pathway was regulated.
:ref:`Kokaji et al. <ref-kokaji2020>` found that roughly half of all regulated
reactions in obese mouse liver received *opposing* input from the two axes.
TransNet calls these **controversial** reactions.

:func:`~transnet.reaction_regulation_table` returns one row per reaction:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Column
     - Meaning
   * - ``gene_axis``
     - Direction (+1/-1/0) of enzyme-level regulation
   * - ``gene_axis_evidence``
     - Which chain supported it: ``"protein"``, ``"gene_protein"``,
       ``"gene"``, ``"tf_gene_protein"``, or ``None``
   * - ``gene_axis_via``
     - The molecules that carried the evidence
   * - ``metabolite_axis``
     - Direction of metabolite-mediated regulation
   * - ``allosteric_regulators``
     - Effectors driving it, each with its effect
   * - ``substrates_changed`` / ``products_changed``
     - Mass-action contributions
   * - ``consensus``
     - Overall direction when the axes agree
   * - ``controversial``
     - True when the two axes point opposite ways

Each axis uses whatever the network has. The metabolite axis needs only a
Metabolome and a Reactions layer. The gene axis takes the longest chain
available: TF → gene → protein where all three exist, gene → protein without TF
evidence, protein alone without a transcriptome, or the transcript itself where
there is no proteome at all. A row whose ``gene_axis_evidence`` is ``None`` was
assessed on the metabolite axis alone, so that claim rests on one axis.

Currency metabolites
~~~~~~~~~~~~~~~~~~~~

NADH, ATP, H\ :sub:`2`\ O and their kin take part in hundreds of reactions.
Counted as substrates and products, a change in one marks every dehydrogenase in
the network as metabolite-regulated, which is a statement about connectivity
rather than biology. They are therefore excluded from mass-action
contributions by default
(:data:`transnet.biology.schema.CURRENCY_METABOLITES`; override with
``currency_metabolites=`` or disable with
``exclude_currency_metabolites=False``).

Their *allosteric* roles are never excluded. AMP activating phosphofructokinase
is exactly the specific regulation this analysis exists to find, and the
distinction is only possible because the edges are typed: the same metabolite
can be suppressed in one role and kept in another.

On a real mouse dataset this moved the controversial fraction from 7 % to 3 %,
and replaced a list of NADH-driven oxidoreductases with specific enzymes --
glutamate dehydrogenase, 2-oxoglutarate dehydrogenase, carnitine
acetyltransferase.

A third axis, where phosphoproteomics exists
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An enzyme's activity is set by how much of it there is and by its modification
state. Where a study measures phosphosites,
:func:`~transnet.map_modification_sites` attaches each site to its protein as a
Signaling node, and the table gains ``phospho_axis``,
``phospho_axis_effect``, ``n_phosphosites_changed`` and ``phospho_axis_via``.

These are separate columns rather than part of ``gene_axis``. Whether more
phosphorylation at a given site raises or lowers catalytic activity is
site-specific and usually unrecorded: the schema gives ``phosphorylation`` a sign
of 0, and KEGG signs a relation only where it annotates activation or inhibition.
``phospho_axis`` therefore reports which way the site moved and
``phospho_axis_effect`` what that does to the reaction, which stays 0 unless the
edge carries a sign. Merging the two would assert a direction the data does not
have, and that direction would then propagate into every path traced through the
reaction.

Sites on one enzyme can move in opposite directions, which leaves
``phospho_axis`` at 0. ``n_phosphosites_changed`` is what separates that case
from no site having changed. The reaction class worth naming is the one where
``gene_axis`` is 0 and ``phospho_axis`` is not: the enzyme's amount held steady
while its modification state moved, and an abundance-only reading calls it
unregulated. :doc:`motrpac_study` finds 114 such reactions in heart and 69 in
trained muscle.

.. code-block:: python

    table = reaction_regulation_table(graph)
    table[table["controversial"]][
        ["reaction", "name", "gene_axis", "metabolite_axis",
         "allosteric_regulators"]
    ]

:Function: :func:`~transnet.reaction_regulation_table`
:Needs: A Reactions layer. The gene axis needs a Proteome or a Transcriptome,
   the metabolite axis a Metabolome, and the phospho axis a Signaling layer.
   Whichever are present is what gets used.
:Example: ``notebooks/walkthroughs/reaction_regulation.py``
:Reference: :ref:`Kokaji et al. 2020 <ref-kokaji2020>`,
   :ref:`Egami et al. 2021 <ref-egami2021>`


.. _pathway-balance:

Per-pathway regulation balance
------------------------------

Two pathways can change equally much, but through different mechanisms: one
through the amount of its enzymes, the other through its metabolites.
:func:`~transnet.regulation_axis_summary` summarises the regulation axes per
pathway: how many of a pathway's regulated reactions each axis activates or
inhibits, and what fraction are controversial.
:func:`~transnet.visualization.plot_regulation_axes` draws this as paired bar
charts, in the style of :ref:`Kokaji et al. <ref-kokaji2020>`

The summary needs to know which pathway each reaction belongs to.
:func:`transnet.api.kegg_reaction_pathways` reads this from KEGG, restricted
to the pathways present in one organism and without KEGG's overview maps. A
reaction in several pathways is counted in each. The bundled example has its
own map, :func:`~transnet.load_example_pathways`.

.. code-block:: python

    from transnet import regulation_axis_summary
    from transnet.api import kegg_reaction_pathways
    from transnet.visualization import plot_regulation_axes

    pathways = kegg_reaction_pathways(table["reaction"], organism="mmu")
    summary = regulation_axis_summary(table, pathway_map=pathways)
    plot_regulation_axes(summary)

:Function: :func:`~transnet.regulation_axis_summary`,
   :func:`transnet.api.kegg_reaction_pathways`
:Needs: A regulation table from :ref:`the axis attribution <regulation-axes>`
   and a reaction-to-pathway map. Without a map, all reactions are summarised
   in a single row.
:Example: ``notebooks/walkthroughs/reaction_regulation.py``, and every study
:Reference: :ref:`Egami et al. 2021 <ref-egami2021>`, Fig. 5


.. _metabolite-roles:

Metabolite regulatory roles
---------------------------

A metabolite can sit in the network two ways. It can be consumed and produced,
mass flowing through the pathway, or it can be an allosteric effector: changing
how fast an enzyme works without being changed itself. BRENDA annotates the
second kind. A differential metabolite that is a regulator is a mechanistic
hypothesis; a differential metabolite that is only a substrate is often a
consequence of flux rather than a cause of anything.

The question this answers is how many of the changed metabolites actually
regulate something, and whether that share is higher than the measured
background. :func:`~transnet.metabolite_regulatory_roles` classifies every
metabolite from the allosteric edges on the graph:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Column
     - Meaning
   * - ``role``
     - ``activator``, ``inhibitor``, ``both``, ``substrate/product only``,
       or ``none``
   * - ``is_allosteric_regulator``
     - True for the first three
   * - ``n_reactions_activated`` / ``n_reactions_inhibited``
     - How much it regulates
   * - ``reactions_activated`` / ``reactions_inhibited``
     - Which reactions, by id
   * - ``n_enzymes_regulated``
     - Distinct EC numbers affected
   * - ``is_substrate`` / ``is_product``
     - Whether it also participates in mass flow
   * - ``regulated``, ``log2fc``, ``qvalue``
     - The measured change, where data was mapped

:func:`~transnet.regulatory_role_enrichment` then answers the question with a
count and a test: ``counts["n_differential_regulators"]`` out of
``counts["n_differential"]``, plus Fisher's exact against the measured
background, run for regulators overall and for activators and inhibitors
separately. Its ``regulators`` frame is the shortlist of metabolites that both
changed and have somewhere to act.

.. code-block:: python

    from transnet import metabolite_regulatory_roles, regulatory_role_enrichment

    roles = metabolite_regulatory_roles(graph)
    result = regulatory_role_enrichment(roles)

    counts = result["counts"]
    print(f"{counts['n_differential_regulators']} of "
          f"{counts['n_differential']} differential metabolites regulate "
          f"an enzyme ({counts['fraction_differential_regulators']:.0%}; "
          f"background {counts['fraction_background_regulators']:.0%})")

    result["regulators"][
        ["name", "log2fc", "role", "reactions_inhibited"]
    ]

Without allosteric edges the function warns and reports every metabolite as a
non-regulator, rather than failing, so a network built without BRENDA gives an
answer that looks like a negative result and is not one.

.. note::

   BRENDA reports effectors by compound *name*, while the metabolome is keyed by
   KEGG compound id. TransNet resolves names (and KEGG synonyms) to ids when the
   edges are built, and logs any effector name it could not place. Matching the
   two directly would silently produce no allosteric edges at all.

:func:`~transnet.visualization.plot_metabolite_regulators` shows the share
of changed metabolites that are regulators beside the background share, and
lists the regulators.

:Function: :func:`~transnet.metabolite_regulatory_roles`,
   :func:`~transnet.regulatory_role_enrichment`
:Needs: A Metabolome and a Reactions layer, plus BRENDA allosteric edges on the
   Proteome (:meth:`~transnet.biology.layers.Proteome.get_brenda_kinetics`, or
   ``build_networks.py --brenda``).
:Example: ``notebooks/walkthroughs/compare_conditions.py``, and every study
:Reference: :ref:`Kokaji et al. 2020 <ref-kokaji2020>`


.. _signed-paths:

Signed regulatory-path tracing
------------------------------

Every edge in a TransNet network has a sign: +1 if more of the upstream
molecule means more of the downstream one, -1 if it means less, and 0 if the
effect is not known. Multiplying the signs along a chain of edges gives the
sign of the whole chain. For example, ATP inhibits pyruvate kinase (-1) and
pyruvate kinase produces pyruvate (+1), so the path ATP → pyruvate kinase →
pyruvate has sign -1: more ATP should mean less pyruvate.

This turns the network into something that makes predictions.
:func:`~transnet.trace_regulatory_paths` follows directed paths from each
changed molecule to a target layer, usually the metabolome, and predicts the
target's direction as the path sign times the measured direction of the
starting molecule. A *decreased* enzyme on a +1 path therefore predicts a
decrease. Each prediction is then compared with the measurement. Because the
prediction can be wrong, this is a test of the network, not an illustration.

Each row of the result is one path:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Column
     - Meaning
   * - ``sign``
     - product of the edge signs along the path
   * - ``source_regulated``
     - measured direction of the starting molecule
   * - ``predicted``
     - ``sign`` × ``source_regulated``: the direction predicted for the target
   * - ``observed``
     - measured direction of the target (0 = no significant change)
   * - ``consistent``
     - whether the prediction matches the measurement
   * - ``unsigned_steps``
     - edges with unknown sign, such as ChIP-Atlas binding edges; a path with
       any of these is only a tentative prediction

Two rules keep the result meaningful:

* **Paths through unchanged molecules are dropped.** If a molecule in the
  middle of a path was measured and did not change, it cannot have passed a
  change on, so the data contradict the path. Molecules that were not measured
  are kept, because there is no evidence either way.
  ``allow_unchanged_intermediates=True`` switches this off, which is useful
  only to show how signs combine.
* **Report one verdict per target molecule, not per path.**
  :func:`~transnet.path_consistency_summary` gives one row per molecule: the
  direction most of its paths predict (``predicted``), whether that matches
  (``agrees``), and the shortest path that predicts correctly. Paths are not
  independent observations. In MoTrPAC skeletal muscle one metabolite was
  reached by dozens of paths that shared most of their steps; counting each
  path separately turned a result indistinguishable from chance into
  p = 4e-17.

Two kinds of edge are skipped by default, because they connect almost
everything to everything:

* ``protein_interaction`` edges from STRING have neither direction nor sign,
  so they are not regulatory steps. Include them with
  ``exclude_edge_types=()``.
* **Currency metabolites** such as water, ATP and NAD take part in hundreds of
  reactions. A path through one of them joins two reactions that have nothing
  else to do with each other. Include them with ``exclude_nodes=()``, or name
  one in ``sources`` or ``targets`` to use it as an end point.

On the organism-wide mouse network, these two rules reduced 8,454 paths,
almost all of them passing through water, to 165. The shortest of those are the
glycolytic chain: hexokinase → fructose 6-phosphate → phosphofructokinase →
fructose 1,6-bisphosphate.

Paths start at the highest layer the network contains: Signaling if a
phosphoproteome was mapped, otherwise Proteome, otherwise Transcriptome.
``source_layer`` overrides this.

.. code-block:: python

    paths = trace_regulatory_paths(graph, target_layer="Metabolome")
    paths[paths["consistent"] & (paths["unsigned_steps"] == 0)]

    verdicts = path_consistency_summary(paths)      # one row per molecule
    verdicts["agrees"].sum(), (verdicts["predicted"] != 0).sum()

:func:`~transnet.visualization.plot_regulatory_paths` draws the best-supported
paths, one row each, with every molecule coloured by its measured change.

:Function: :func:`~transnet.trace_regulatory_paths`,
   :func:`~transnet.path_consistency_summary`
:Needs: Two measured layers with signed edges between them, and a target layer.
:Example: :doc:`notebooks/walkthroughs/regulatory_paths`, and every study
:Reference: :ref:`Kawata et al. 2018 <ref-kawata2018>`,
   :ref:`Yugi et al. 2016 <ref-yugi2016>`


.. _signal-flow:

Signal flow from a receptor, and propagation
--------------------------------------------

The Signaling layer sits above the Proteome. It holds kinases and other
signalling proteins, joined to the proteome by two edge types:
``phosphorylation`` (a kinase acting on another protein) and ``kinase_tf`` (a
signalling protein acting on a transcription factor). With this layer, a single
path can run from a hormone receptor to a metabolite:

.. code-block:: text

    Insr -> Irs1 -> Pik3r1 -> Akt1 -| Foxo1 -> gene -> enzyme -> REACTION -> metabolite

The layer can come from two sources:

* :func:`transnet.api.kegg_signaling_relations` reads KEGG's signalling
  pathway maps, and ``maintenance/build_networks.py --signaling`` adds the
  relations to a network. KEGG's annotation gives the sign: +1 for
  activation, -1 for inhibition, and 0 for a phosphorylation whose effect is
  not recorded.
* :func:`~transnet.map_modification_sites` adds a study's *measured*
  phosphorylation sites. Each site becomes its own Signaling node with an edge
  to its protein, because several sites on one protein can move in opposite
  directions. These edges also feed the phospho axis of the
  :ref:`regulation-axis analysis <regulation-axes>`.

On the bundled example, 26 paths run from the insulin receptor to a
metabolite, and 14 of them pass through a ``kinase_tf`` step. All 14 reach
Srebf1 and none reaches Foxo1, although Akt1 inhibits Foxo1 and the edge is in
the network. The reason is the first rule above: Foxo1 was measured and did not
change, so every path through it is contradicted by the data and dropped.

Path tracing lists individual routes. **Propagation** answers the same
question with a single score per molecule.
:func:`~transnet.hierarchical_propagation` starts from the measured changes and
pushes them forward along directed, signed edges, adding up what arrives at
each molecule; a score that passes an inhibiting edge arrives negative. A
molecule reached by many weak routes and one reached by a single strong route
then get different scores, which path tracing cannot express.
:func:`~transnet.downstream_influence` keeps one layer, adds the measured
direction and marks whether each prediction agrees;
:func:`~transnet.visualization.plot_downstream_influence` plots predicted
score against measured direction.

Use path tracing to name the mechanism behind one change, and propagation to
get a predicted direction for every molecule in a layer.
:func:`transnet.analysis.network_propagation.random_walk_with_restart` is a
different tool: it spreads scores without regard to direction or sign, to find
what lies *near* a set of molecules.

:func:`~transnet.visualization.plot_transomic_network` draws the Signaling
layer as the top row when the network has one.

:Function: :func:`transnet.api.kegg_signaling_relations`,
   :func:`~transnet.map_modification_sites`,
   :func:`~transnet.hierarchical_propagation`,
   :func:`~transnet.downstream_influence`
:Needs: For signal flow, a Signaling layer, from KEGG
   (``build_networks.py --signaling``) or from measured phosphorylation sites.
   Reaching the metabolome also needs a Proteome and a Reactions layer.
   Propagation works from any layer.
:Example: :doc:`notebooks/walkthroughs/regulatory_paths`;
   :doc:`motrpac_study` maps phosphorylation sites from six tissues.
:Reference: :ref:`Yugi et al. 2016 <ref-yugi2016>`,
   :ref:`Kawata et al. 2018 <ref-kawata2018>`


.. _cross-layer-connectivity:

Cross-layer connectivity
------------------------

If 95 % of a network's edges run within a layer, it is in effect a set of
separate single-omics networks, and trans-omic conclusions drawn from it are
weak. This analysis measures that before any conclusions are drawn.

:func:`~transnet.cross_layer_connectivity` returns
the layer × layer edge-count matrix, a per-relationship breakdown, per-layer
coverage (also available on its own as :func:`~transnet.layer_coverage`), and
the headline ``cross_layer_fraction``.

:func:`~transnet.layer_coverage` is the one to check straight after mapping
data: it reports, per layer, how many nodes carry a measurement and how many
were called regulated. A layer with near-zero ``measured_fraction`` usually
means an identifier mismatch rather than a biological result.

.. code-block:: python

    connectivity = cross_layer_connectivity(graph)
    connectivity["cross_layer_fraction"]
    connectivity["matrix"]

:Function: :func:`~transnet.cross_layer_connectivity`,
   :func:`~transnet.layer_coverage`
:Needs: Any network. Coverage additionally needs omics mapped onto it.
:Example: ``notebooks/walkthroughs/build_network.py``
:Reference: :ref:`Sugimoto et al. 2024 <ref-sugimoto2024>`


.. _transomic-hubs:

Trans-omic hub identification
-----------------------------

A hub in a single-omics network is whatever has most neighbours in that one
layer. The interesting molecule in a trans-omic network is the one that reaches
*several* layers, which a within-layer degree ranking cannot see.

:func:`~transnet.transomic_hubs` ranks molecules by
degree among the responsive nodes, where :ref:`Morita et al. <ref-morita2025>`
define hubs as the top 2 %, and adds two multi-layer quantities:
``n_layers_touched``, and a
``versatility`` score that weights degree by the number of layers a node
reaches. A molecule with 20 neighbours spread over three layers outranks one
with 20 neighbours inside a single layer.

In healthy mouse liver this surfaces ATP and AMP, at degrees above 200,
coordinating reactions throughout metabolism. A within-layer degree ranking
would never identify them as the coordinators they are.

.. code-block:: python

    hubs = transomic_hubs(graph, top_percent=2)
    hubs[hubs["is_hub"]]

:Function: :func:`~transnet.transomic_hubs`
:Needs: A mapped network with at least two layers.
:Example: ``notebooks/walkthroughs/temporal_and_hubs.py``
:Reference: :ref:`Morita et al. 2025 <ref-morita2025>`,
   :ref:`De Domenico et al. 2015 <ref-dedomenico2015>`


.. _response-timing:

Temporal and dose structure on the network
------------------------------------------

Timing is a property of the response, and connectivity a property of the
network. Putting them together asks whether the wiring explains the order in
which things moved.

:func:`~transnet.assign_temporal_parameters`
computes a half-response time (and EC50 from dose–response data) per molecule
and writes them onto nodes. :func:`~transnet.temporal_network_structure` then
tests the architecture:

* Spearman correlation of degree against t½. A negative value means hubs lead
  the response;
* the distribution of correlations between the time courses of *connected*
  molecules. Two peaks near ±0.75 mean that connected molecules move together
  (or exactly opposite); a single peak near zero means they do not;
* median response time per layer, showing the order in which layers move.

:func:`~transnet.split_by_response_class` splits the network into fast/slow and
sensitive/insensitive subnetworks.

:ref:`Morita et al. <ref-morita2025>` found the degree-versus-t½ relationship
intact in healthy liver
and destroyed in obese liver *while the network's structure was unchanged*:
structural robustness with temporal vulnerability. That finding is inaccessible
without putting the timing onto the graph.

.. code-block:: python

    assign_temporal_parameters(graph, {"Metabolome": timecourse},
                               time_columns=[0, 5, 15, 30, 60])
    structure = temporal_network_structure(graph)
    structure["degree_vs_thalf"]

:Function: :func:`~transnet.assign_temporal_parameters`,
   :func:`~transnet.temporal_network_structure`,
   :func:`~transnet.split_by_response_class`
:Needs: A time course of at least three points for one layer, and EC50 needs a
   dose-response series.
:Example: ``notebooks/walkthroughs/temporal_and_hubs.py``
:Reference: :ref:`Morita et al. 2025 <ref-morita2025>`,
   :ref:`Kawata et al. 2018 <ref-kawata2018>` (EC50 by t-half response classes)


.. _tf-activity:

Responsive transcription-factor inference
-----------------------------------------

ChIP-Atlas binding says a factor *can* regulate a
gene, not whether it activates or represses it, so ``transcriptional_regulation``
edges carry no sign. The responsive targets supply it. This is how the Kuroda
laboratory identifies differentially regulated transcription factors, and it is
the principled way to turn "regulation, sign unknown" edges into a direction.

:func:`~transnet.transcription_factor_activity` tests
each factor's measured targets for enrichment among responsive genes (one-sided
Fisher's exact test, Benjamini-Hochberg across factors), then asks whether its
responsive targets move one way more than chance (binomial test). Mostly-up
targets mean an activated activator *or a relieved repressor*; the factor's own
measured change is reported alongside, and a factor inferred active whose level
did not change is regulated post-translationally.

.. code-block:: python

    tf = transcription_factor_activity(graph, min_targets=10)
    tf[["name", "n_responsive_targets", "n_up", "n_down",
        "q_value", "inferred_activity", "factor_regulated"]]

:func:`~transnet.visualization.plot_tf_activity` ranks the factors by
significance and shows the direction of their targets.

:Function: :func:`~transnet.transcription_factor_activity`
:Needs: Transcriptome measurements and ``transcriptional_regulation`` edges,
   genome-wide. A targeted gene panel carries too few targets per factor to
   test, and ChIP-Atlas coverage varies sharply by organism.
:Example: ``notebooks/walkthroughs/transcription_factors.py``;
   :doc:`liver_timecourse` scores the inference against a published one.
:Reference: :ref:`Kokaji et al. 2022 <ref-kokaji2022>`,
   :ref:`Maehara et al. 2025 <ref-maehara2025>`


.. _concordance:

Transcript–protein concordance
------------------------------

A protein can change because its transcript did, or for a reason downstream of
transcription: translation rate, degradation, modification. Separating the two
needs both layers measured on the same samples.

:func:`~transnet.expression_concordance` classifies
every ``translation`` edge whose gene and protein were both measured as
*concordant*, *protein-only*, *transcript-only*, *discordant* or *unchanged*,
and reports the fold-change correlation. A protein-only change whose transcript
moved but missed the threshold is not evidence of post-transcriptional control
on its own.

The categories depend on each layer's statistical power: with a transcriptome
test that detects few changes, most changed proteins land in *protein-only*
whatever the biology. When standard errors are mapped
(``map_omics_to_network(..., se_column=)``), each pair is also tested directly
-- is the protein change larger than the transcript change? -- and
``protein_beyond_transcript`` marks the pairs where it is. That, not the
category, is the evidence for post-transcriptional regulation -- provided both
platforms report fold changes on comparable scales. Isobaric proteomics
compresses ratios, which makes a positive call conservative; other platform
pairings should be checked before the counts are compared across studies.

The same distinction is carried into :ref:`reaction regulation axes <regulation-axes>` as
``gene_axis_transcript_support``, so each enzyme-driven reaction says whether
its regulation was transcriptional.

.. code-block:: python

    result = expression_concordance(graph)
    result["counts"], result["correlation"]

mRNA levels explain only part of the variation in protein levels, so a large
protein-only class is normal biology. A transcriptome-only study cannot see
this class at all. :func:`~transnet.visualization.plot_expression_concordance`
plots protein against transcript change, coloured by category.

:Function: :func:`~transnet.expression_concordance`
:Needs: Transcriptome and Proteome measurements on the same samples. Standard
   errors (``map_omics_to_network(..., se_column=)``) enable the direct test.
:Example: ``notebooks/studies/brown_adipocytes.py``,
   ``notebooks/studies/motrpac_rat.py``
:Reference: :ref:`Liu et al. 2016 <ref-liu2016>`


.. _factors-on-the-network:

Factors, read through the network
---------------------------------

A factor model finds molecules that co-vary. Whether they are related by the
biochemistry or only by the numbers is a question the factorisation cannot ask
of itself.

A joint factor loads on transcripts, proteins
and metabolites at once, and nothing in the factorisation says whether those
are *the same biology* -- whether the transcripts it picks encode the enzymes
acting on the metabolites it picks. The network does.

Seven readings, each asking something of a factor that the factor model cannot
ask of itself:

:func:`~transnet.analysis.factors.factor_design_association`
    Every factor against the whole design in one linear model, with
    interactions -- a factor carrying a *sex-specific* training response
    otherwise looks like a sex factor with a weak timepoint effect.
:func:`~transnet.analysis.factors.factor_cross_layer_agreement`
    Project the samples onto a factor *within each layer* and correlate the
    projections. A factor fitted on stacked layers can be carried by one layer
    alone; this says whether it is.
:func:`~transnet.analysis.factors.factor_network_coherence`
    Are the factor's top features directly joined on the network, more than
    random features are?
:func:`~transnet.analysis.factors.factor_network_propagation`
    Diffuse each layer's top features separately and measure how much the
    profiles overlap, against random features **of the same layers** -- the
    transcriptome is a hundred times the size of the metabolome, so a null
    drawn uniformly would be all transcripts. The top nodes are ranked by how
    much more signal they receive than from random seeds, so the network's
    hubs do not win by default.
:func:`~transnet.analysis.factors.network_guided_imputation`
    Fill a missing value from the feature's network neighbours rather than the
    column mean: for a metabolite missing in one sample, the substrates and
    products of its reactions are a better guess than the cohort average.
:func:`~transnet.analysis.factors.sample_pairing_check`
    A joint model treats column *i* as one sample in every layer. Files that
    share column *names* need not share samples. Within a group, a sample
    whose transcript of a gene runs high should, if it is the same sample,
    tend to have that protein high too -- so the statistic is computed per
    gene across the replicates, and compared with shuffled pairings.
:func:`~transnet.analysis.factors.pairing_robustness`
    Whether a factor's cross-layer agreement comes from the design or from
    sample-to-sample variation. ``between-group`` factors are carried by
    group differences whose labels are certain and are interpretable whatever
    the pairing; ``within-group`` factors live in variation shared by the
    layers within a group, real when the pairing is confirmed; ``single-layer``
    ones are not joint at all.

The pairing check confirms the pairing of MoTrPAC, whose layers are keyed by
animal (p = 0.002), and of the brown adipocyte data (p = 0.002). On the Uematsu
panel, with only 17 genes measured in both layers, it points the right way
without reaching significance.

On the brown adipocyte data, four of five factors have a cross-layer overlap
two to three times the null (q = 0.025). The strongest lands on
branched-chain amino acid breakdown (*Bckdhb*, 3-hydroxyisobutyrate) beside
the complex III assembly factor *Uqcc4*. Branched-chain amino acids are a known
fuel for heat production in brown fat, and here they were found from the
network rather than from a pathway list.

:func:`~transnet.visualization.plot_factor_overview`,
:func:`~transnet.visualization.plot_factor_scores` and
:func:`~transnet.visualization.plot_factor_network` draw these readings.

:Function: :func:`~transnet.analysis.factors.fit_factors` and the
   ``factor_*`` readings above
:Needs: Sample-level matrices for at least two layers, measured on the same
   samples. Differential tables are not enough: these readings work on the
   sample axis, which is why they appear in the studies and not the
   walkthroughs.
:Example: ``notebooks/studies/brown_adipocytes.py``,
   ``notebooks/studies/motrpac_rat.py``
:Reference: :ref:`Argelaguet et al. 2020 <ref-argelaguet2020>`,
   :ref:`Cowen et al. 2017 <ref-cowen2017>`


.. _regulatory-motifs:

Regulatory motifs
-----------------

A pathway diagram tells you that hexokinase is inhibited by its own product.
In a typed, directed, signed network such patterns can be *found*, wherever
they occur, without knowing in advance where to look: any reaction whose
product has an inhibiting edge back to that reaction is product inhibition. On the bundled example the search returns
hexokinase inhibited by glucose-6-phosphate and pyruvate dehydrogenase
inhibited by acetyl-CoA, neither of them named anywhere in the code.

:func:`~transnet.regulatory_motifs` reports three patterns and their sign
products:

``product_inhibition`` / ``product_activation``
    A reaction regulated by the metabolite it makes. These are the reactions
    that can slow down while their enzyme rises -- the mechanism behind a
    *controversial* call in the regulation axes.
``allosteric_feedback``
    A metabolite made by one reaction regulating another reaction that feeds
    it: feedback at pathway range rather than at one step.
``feed_forward``
    A transcription factor whose targets include two enzymes of the same
    reaction -- control applied twice, on different timescales.

Currency metabolites are excluded by default (ATP is a product of half the
network, and counting it makes every reaction look self-regulating), and
``responsive_only=True`` turns a census of the network into a census of *this*
response. The transcription factor must itself have responded: without that
rule, ChIP-Atlas's thousands of targets per factor produced 10,272
feed-forward instances on one mouse dataset and no information.

:Function: :func:`~transnet.regulatory_motifs`
:Needs: A Reactions layer with allosteric edges. Feed-forward motifs also need
   ``transcriptional_regulation`` edges.
:Example: ``notebooks/walkthroughs/network_topology.py``
:Reference: :ref:`Milo et al. 2002 <ref-milo2002>`


.. _structural-vulnerability:

Structural vulnerability
------------------------

A *cut molecule* (in graph terms, a cut vertex) is a molecule whose removal
splits the network into separate pieces. In a trans-omic network these are the
places where the response in one layer reaches another through a single route:
if that molecule were missing, or its enzyme inhibited, the pieces would no
longer be connected.

:func:`~transnet.structural_vulnerability` ranks
the cut molecules of the responsive network by how much of it they strand:
``fragments`` counts the pieces left behind, ``largest_loss`` the share of the
response disconnected. On the bundled example the two most consequential are
the transcription factors ChREBP (*Mlxipl*) and SREBP-1 (*Srebf1*), each
holding together a quarter of the response.

:Function: :func:`~transnet.structural_vulnerability`
:Needs: A responsive subnetwork. A network whose response is already
   disconnected has no cut vertex to find, which is reported rather than
   treated as a null result.
:Example: ``notebooks/walkthroughs/network_topology.py``
:Reference: :ref:`Morita et al. 2025 <ref-morita2025>`


.. _convergence-null:

Convergence against a null model
--------------------------------

Many analyses in this catalogue rest on reactions where a changed enzyme and
a changed metabolite meet. Some such meetings are expected by chance alone: if
many molecules change, and some reactions have many connections, changed
molecules will sometimes land on the same reaction.

:func:`~transnet.convergence_significance` measures how many meetings chance
would produce. It keeps the network and the number of changed molecules in
each layer fixed, reassigns at random which molecules count as changed, and
counts again, many times. It returns the real count, the mean of the random
counts, a z-score and a p-value, and the random counts themselves, which
:func:`~transnet.visualization.plot_convergence_null` draws as a histogram.

On the brown adipocyte study, 170 reactions have both axes changed, against 19
expected by chance (z = +4.4, p = 0.015): the layers converge far more than
chance would produce. On the Uematsu panel the numbers are 14 against 6.3
(z = +1.9, p = 0.05), which is only just significant.

:Function: :func:`~transnet.convergence_significance`
:Needs: A Reactions layer with a measured enzyme layer and a Metabolome. The
   enzyme side reads ``catalysis`` edges, so a study with no proteome scores
   zero convergence by construction and is warned about it.
:Example: ``notebooks/walkthroughs/network_topology.py``
:Reference: :ref:`Maslov and Sneppen 2002 <ref-maslov2002>`


Figures: what each analysis finds
---------------------------------

Every analysis has a figure in :mod:`transnet.visualization`, drawing the
result a per-layer analysis could not have produced. The study notebooks save
one per analysis, beside the table it came from.

.. list-table::
   :header-rows: 1
   :widths: 46 54

   * - Figure
     - Shows
   * - :func:`~transnet.visualization.plot_transomic_network`
     - network reconstruction -- stacked layers, most-regulated reactions
   * - :func:`~transnet.visualization.plot_transomic_network_interactive`
     - network reconstruction -- the same, in a browser: hover, zoom, and a standalone HTML file
   * - :func:`~transnet.visualization.plot_layer_connectivity`
     - cross-layer connectivity -- cross-layer edge counts
   * - :func:`~transnet.visualization.plot_axis_composition`
     - reaction regulation axes -- which axis regulates the reactions
   * - :func:`~transnet.visualization.plot_controversial_reactions`
     - reaction regulation axes -- the tug-of-war, per enzyme
   * - :func:`~transnet.visualization.plot_regulation_axes`
     - per-pathway balance -- which axis activates or inhibits each pathway
   * - :func:`~transnet.visualization.plot_metabolite_regulators`
     - metabolite regulatory roles -- changed metabolites that act on an enzyme
   * - :func:`~transnet.visualization.plot_regulatory_paths`
     - signed regulatory paths -- signed paths, predicted against observed
   * - :func:`~transnet.visualization.plot_transomic_hubs`
     - trans-omic hubs -- cross-layer hubs
   * - :func:`~transnet.visualization.plot_temporal_structure`
     - response timing -- degree against response time
   * - :func:`~transnet.visualization.plot_tf_activity`
     - transcription-factor activity -- factors behind the responsive genes
   * - :func:`~transnet.visualization.plot_expression_concordance`
     - transcript-protein concordance -- transcript against protein
   * - :func:`~transnet.visualization.plot_downstream_influence`
     - predicted against measured metabolites
   * - :func:`~transnet.visualization.plot_condition_comparison`
     - relationships kept, lost and reversed between conditions
   * - :func:`~transnet.visualization.plot_values_heatmap`
     - the same molecules across tissues or clinical groups
   * - :func:`~transnet.visualization.plot_layer_changes`
     - changed molecules per layer, which is what a per-layer analysis reports
   * - :func:`~transnet.visualization.plot_regulatory_motifs`
     - regulatory motifs -- instances per motif type
   * - :func:`~transnet.visualization.plot_structural_vulnerability`
     - structural vulnerability -- how far the response falls apart per molecule
   * - :func:`~transnet.visualization.plot_convergence_null`
     - convergence null -- the shuffled null with the observed count marked
   * - :func:`~transnet.visualization.plot_factor_scores`
     - factors read through the network -- factor scores across the design
   * - :func:`~transnet.visualization.plot_factor_overview`
     - factors -- variance, design association and network coherence together
   * - :func:`~transnet.visualization.plot_factor_network`
     - factors -- a factor's top loadings placed on the network
   * - :func:`~transnet.visualization.plot_community_network`
     - communities, whether they span layers, and one community drawn with
       labels
   * - :func:`~transnet.visualization.plot_network_metrics`
     - degree distribution, component sizes and molecules per layer

To save any figure, or the network itself, in another format, see
:doc:`notebooks/walkthroughs/export_network`.

Colours come from one validated palette (:mod:`transnet.visualization.palette`):
red and blue only ever mean direction, layers and axes take hues that are neither.

Statistics reported against a baseline
--------------------------------------

Every result is reported against what chance would give:

* **Directional predictions** (signed paths, propagation): the share of
  correct predictions is compared with the 50 % expected from guessing up or
  down at random, using a binomial test (:func:`transnet.analysis.versus_chance`).
* **Enrichments** (metabolite regulatory roles, transcription-factor
  activity): the result states when the test finds no over-representation.
* **Hubs** are ranked among the molecules that responded, following
  :ref:`Morita et al. <ref-morita2025>` Over the whole network, transcription
  factors with thousands of ChIP-Atlas targets would rank highest whatever the
  data.

.. _comparing-conditions:

Comparing conditions
--------------------

:func:`~transnet.compare_transomic_networks` compares two networks
layer-by-layer and *relationship-by-relationship*. A generic differential
network comparison reports which nodes and edges differ; this reports which
*kinds of regulation* were gained and lost. A condition that loses its
transcriptional arm but keeps its allosteric one is a different biological story
from the reverse, and only an edge-type-aware comparison distinguishes them.

It also lists the molecules that changed in both conditions but in opposite
directions, which a comparison of two gene lists would count as agreement.
:func:`~transnet.visualization.plot_condition_comparison` draws both.

Tissues are compared this way, as separate networks, rather than connected.
Inter-organ networks (:ref:`Egami et al. 2021 <ref-egami2021>`) join tissues
through circulating metabolites, which TransNet does not do.

:Function: :func:`~transnet.compare_transomic_networks`
:Needs: Two networks built from the same reference, each with omics mapped on.
:Example: ``notebooks/walkthroughs/compare_conditions.py``;
   :doc:`motrpac_study` (six tissues) and :doc:`liver_timecourse` (two
   genotypes)
:Reference: :ref:`Egami et al. 2021 <ref-egami2021>`
