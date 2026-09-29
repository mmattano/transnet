.. _transomics-analyses:

The trans-omics analysis catalogue
==================================

Trans-omics is not multi-omics with more layers. Multi-omics integration asks
which molecules co-vary across data types; trans-omics asks how a signal
*propagates* through a biochemical network: which regulatory relationship
carried it, in which direction, and whether the observed changes are consistent
with that mechanism.

Trans-omics networks, largely defined by the Kuroda group, are built on one object: a typed,
directed, signed regulatory hierarchy converging on the metabolic reaction.

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

``notebooks/walkthroughs/regulatory_paths.py`` demonstrates this by running the same
analysis twice, with and without the optional Signaling layer.


Identifiers, and why the mapping report exists
----------------------------------------------

Every layer of a network is IDed by one identifier -- Entrez for genes, UniProt
for proteins, KEGG compound ids for metabolites -- while a real omics table is
IDed by whatever the instrument and pipeline produced: Ensembl gene ids, RefSeq
protein accessions, PubChem CIDs, RefMet names. Mismatch produces no error, only
an analysis quietly computed over nothing. Use :func:`~transnet.map_omics_to_network`
to get a coverage report and see which
names in a layer are mapped poorly.

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

A layer reported at 0 % is almost always an identifier-type mismatch, not due to
a biological result.


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

A reaction can be regulated through the *gene-expression axis*, meaning how much
enzyme there is, or through the *metabolite axis*, meaning how hard that enzyme
works. Both are regulation, they are measured in different layers, and they can
disagree. Separating them is what most clearly distinguishes trans-omics from
pathway enrichment. :ref:`Kokaji et al. <ref-kokaji2020>` found roughly half of all
differentially regulated reactions in obese mouse liver receiving *opposing*
input from the two axes. TransNet calls these **controversial** reactions.

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

Two pathways can be equally "changed" and yet be changed by entirely different
mechanisms, one transcriptionally and one allosterically. That distinction is
what this reports. :func:`~transnet.regulation_axis_summary` rolls the
regulation axes up per pathway: how many of a pathway's reactions each axis
activates or inhibits, and what fraction is controversial.
:func:`~transnet.visualization.plot_regulation_axes` draws it as the paired
bar figure of :ref:`Kokaji et al. <ref-kokaji2020>`

.. code-block:: python

    from transnet import regulation_axis_summary
    from transnet.visualization import plot_regulation_axes

    summary = regulation_axis_summary(table, pathway_map=reaction_to_pathway)
    plot_regulation_axes(summary)

:Function: :func:`~transnet.regulation_axis_summary`
:Needs: A regulation table from :ref:`the axis attribution <regulation-axes>`.
   A pathway map gives per-pathway rows; without one the whole network is
   summarised as a single row.
:Example: ``notebooks/walkthroughs/reaction_regulation.py``
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

Tissues are compared this way rather than connected. Inter-organ networks
(:ref:`Egami et al. 2021 <ref-egami2021>`) join tissues through circulating
metabolites, which TransNet does not do.

:Function: :func:`~transnet.metabolite_regulatory_roles`,
   :func:`~transnet.regulatory_role_enrichment`
:Needs: A Metabolome and a Reactions layer, plus BRENDA allosteric edges on the
   Proteome (:meth:`~transnet.biology.layers.Proteome.get_brenda_kinetics`, or
   ``build_networks.py --brenda``).
:Function: :func:`~transnet.compare_transomic_networks`
:Needs: Two networks built from the same reference, each with omics mapped on.
:Example: ``notebooks/walkthroughs/compare_conditions.py``
:Reference: :ref:`Egami et al. 2021 <ref-egami2021>`
:Reference: :ref:`Kokaji et al. 2020 <ref-kokaji2020>`


.. _signal-flow:

Signal flow through the hierarchy
---------------------------------

A stimulus arrives at a receptor. Which metabolite concentrations does it
reach, through which transcription factors and enzymes, and with what sign?

No other analysis here spans the whole hierarchy.
The Signaling layer sits above the Proteome and connects to it two ways:
``phosphorylation`` for a kinase acting on a substrate protein, and
``kinase_tf`` for a signalling protein acting on a transcription factor. With
those edges present, a traced path runs receptor to metabolite in one object:

.. code-block:: text

    Insr -> Irs1 -> Pik3r1 -> Akt1 -| Foxo1 -> gene -> enzyme -> REACTION -> metabolite

The layer comes from two sources, and they answer different questions.

* :func:`transnet.api.kegg_signaling_relations` parses the KGML of an organism's
  signal-transduction pathways into kinase-substrate and kinase-TF relations.
  ``maintenance/build_networks.py --signaling`` adds them to a network. KEGG's
  own subtype supplies the sign: ``activation`` is +1, ``inhibition`` -1, and a
  plain ``phosphorylation`` is 0.
* :func:`~transnet.map_modification_sites` attaches a study's *measured* sites.
  Phosphoproteomics measures sites, not proteins, and a site regulates its
  protein's activity rather than reporting how much of it there is, so each site
  becomes its own Signaling node with an edge into the protein. Several sites on
  one protein stay separate, because they can move in opposite directions. These
  edges feed the phospho axis of
  :ref:`reaction regulation axes <regulation-axes>`.

:func:`~transnet.hierarchical_propagation` pushes seed scores *forward along
directed, signed edges*, so a score arriving at a reaction through an inhibitor
arrives negative.
:func:`~transnet.downstream_influence` restricts the result to one layer and
compares the predicted direction against measurement. Path tracing
(:ref:`signed regulatory paths <signed-paths>`) answers the same question by
enumerating routes instead, and the two disagree in a useful way: a target
reached by many weak routes scores differently from one reached by a single
strong route.

Ordinary undirected random-walk-with-restart remains available in
:func:`transnet.analysis.network_propagation.random_walk_with_restart` for
diffusion-based similarity, which is a different question and has no direction.

The signs earn their keep on the bundled example: 26 paths run from the
insulin receptor to a metabolite and 14 pass through a ``kinase_tf`` step. All 14
reach Srebf1; none reaches Foxo1, although Akt1 inhibits Foxo1 and that edge is
in the network. Foxo1 was measured and did not change, so paths through it are
contradicted by the data and dropped. The signed branch the data refuse is
removed without anyone having to spot it.

:func:`~transnet.visualization.plot_transomic_network` draws the
Signaling layer as the top row when a network has one, and
:func:`~transnet.visualization.transomic_backbone` walks up from each selected
enzyme to the sites measured on it. Its ``max_signalling_per_enzyme`` governs how
many, and 0 leaves the layer out.

:Function: :func:`transnet.api.kegg_signaling_relations`,
   :func:`~transnet.map_modification_sites`,
   :func:`~transnet.hierarchical_propagation`,
   :func:`~transnet.downstream_influence`
:Needs: A Signaling layer, from KEGG KGML (``build_networks.py --signaling``)
   or from measured modification sites. Reaching the metabolome also needs a
   Proteome and a Reactions layer.
:Example: ``notebooks/walkthroughs/regulatory_paths.py``;
   ``notebooks/studies/motrpac_rat.py`` maps phosphosites from six tissues,
   and :doc:`kokaji_study` has a signalling readout but no proteome, which is
   what the layer needs to reach the enzyme axis.
:Reference: :ref:`Yugi et al. 2016 <ref-yugi2016>` (kinase-substrate as one of
   the five connection technologies), :ref:`Kawata et al. 2018 <ref-kawata2018>`


.. _signed-paths:

Signed regulatory-path tracing
------------------------------

Does the sign product along signal → TF → gene → enzyme → reaction → metabolite
match the observed change? A network that claims a mechanism can be asked to
predict, and the prediction can be wrong, which is what makes this a test rather
than an illustration.

:func:`~transnet.trace_regulatory_paths` enumerates
directed paths through the typed graph, multiplies the edge signs, and compares
the prediction with the target's measured direction. It reports
``unsigned_steps`` separately, because a path through a ChIP-Atlas binding edge
(where activation versus repression is unknown) predicts a direction only
tentatively. :func:`~transnet.path_consistency_summary` gives the shortest
consistent path per target, which is the most parsimonious explanation the network
offers.

A path's prediction is its sign multiplied by the source's measured direction:
a *decreased* enzyme on a +1 path predicts a decrease. (Earlier versions
compared the bare path sign with the target, which inverted the verdict for
every path starting from a decreased molecule.)

Two rules keep the consistency rate meaning what it appears to mean:

* A path through a molecule that *was measured and did not change* is dropped.
  The data it is scored against contradict it: an intermediate that did not
  move passed nothing on. Unmeasured intermediates are kept, since unmeasured
  is unknown rather than unchanged. Pass
  ``allow_unchanged_intermediates=True`` to trace structurally instead, as
  ``notebooks/walkthroughs/regulatory_paths.py`` does when demonstrating sign propagation.
* Quote the rate **per target molecule**, not per path.
  :func:`~transnet.path_consistency_summary` gives one row per molecule with
  the direction most of its paths predict (``predicted``) and whether that
  matches the measurement (``agrees``). Paths are not independent observations:
  in MoTrPAC skeletal muscle one hub metabolite was reached by dozens of paths
  sharing most of their steps, and counting each path separately turned a
  result that is not distinguishable from chance into "p = 4e-17".

Two classes of edge make path tracing meaningless on a real network, and both
are excluded by default:

* ``protein_interaction`` -- a STRING association is undirected and unsigned, so
  it is not a regulatory step, and tens of thousands of them connect almost any
  protein to almost any other. Override with ``exclude_edge_types=()``.
* **currency metabolites** -- a path hopping through water or ATP joins two
  reactions that have nothing to do with each other. Override with
  ``exclude_nodes=()``, or name one in ``sources`` / ``targets`` to keep it as
  an endpoint while still routing around it in the middle.

On a real mouse network these two took the result from 8,454 paths -- almost all
of them water-hops -- to 165, of which the shortest are the glycolytic chain:
hexokinase to F6P to phosphofructokinase to F1,6BP.

``source_layer=None`` infers the highest layer present:
Signaling when a phosphoproteome exists, otherwise Proteome, otherwise
Transcriptome. The result reports the hierarchy used.

.. code-block:: python

    paths = trace_regulatory_paths(graph, target_layer="Metabolome")
    paths[paths["consistent"] & (paths["unsigned_steps"] == 0)]

    verdicts = path_consistency_summary(paths)      # one row per molecule
    verdicts["agrees"].sum(), (verdicts["predicted"] != 0).sum()

:Function: :func:`~transnet.trace_regulatory_paths`,
   :func:`~transnet.path_consistency_summary`
:Needs: Two measured layers with signed edges between them, and a target layer
   to trace to. The starting layer is inferred from what is present.
:Example: ``notebooks/walkthroughs/regulatory_paths.py``
:Reference: :ref:`Kawata et al. 2018 <ref-kawata2018>`,
   :ref:`Yugi et al. 2016 <ref-yugi2016>`


.. _cross-layer-connectivity:

Cross-layer connectivity
------------------------

A network that is 95 % within-layer edges is a stack of separate single-omics
networks wearing one name. This is how to find that out before drawing
conclusions from it.

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
* the distribution of correlations between *connected* molecules, bimodal near
  ±0.75 says the network is coherently driven, unimodal near zero says it is
  not;
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

On the bundled mouse time course this finds Esrra (targets 499 down, 181 up,
the factor itself down), Rxra (targets down, factor level unchanged -- the
signature of a ligand-activated receptor) and the Polycomb components Suz12,
Eed and Rnf2 (targets up while Suz12 and Eed fall: loss of repression).

:Function: :func:`~transnet.transcription_factor_activity`
:Needs: Transcriptome measurements and ``transcriptional_regulation`` edges,
   genome-wide. A targeted gene panel carries too few targets per factor to
   test, and ChIP-Atlas coverage varies sharply by organism.
:Example: ``notebooks/walkthroughs/transcription_factors.py``;
   :doc:`kokaji_study` scores the inference against a published one.
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

mRNA explains only part of protein variation, so a large protein-only class is
expected biology -- and the class a transcriptome-only study attributes to
nothing. On the mouse time course gene and protein fold changes correlate at
rho = 0.12, and 212 pairs move in opposite directions.

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

The first version of the pairing check correlated transcript with protein
*across genes within a sample*. That measures gene abundance, which is the same
whatever the pairing, and on data centred per feature -- MoTrPAC's -- it sees
nothing: it called MoTrPAC, whose layers are keyed by animal, unpaired. The
per-gene version passes MoTrPAC (p = 0.002), which is its validation, and
confirms the brown adipocyte pairing too (p = 0.002). On the Uematsu panel,
with 17 genes in both layers, it leans the right way without reaching
significance.

Two earlier versions of the propagation reading were wrong in instructive
ways. Ranking nodes by raw diffused score returned water and ADP for every
factor, since a hub collects score from any seed set. And a "concentration"
statistic compared against uniformly drawn seeds measured how well connected
the seeds were rather than anything about the factor. Both are fixed by the
two rules above, and on the brown adipocyte data the result is informative:
four of five factors have cross-layer overlap two to three times the null
(q = 0.025), and the strongest lands on branched-chain amino acid catabolism
(*Bckdhb*, 3-hydroxyisobutyrate) beside the complex III assembly factor
*Uqcc4* -- a known thermogenic fuel in brown fat, found from the wiring rather
than from a pathway list.

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

A pathway diagram asserts that hexokinase is
inhibited by its own product. A typed, directed, signed graph lets you *find*
that, and everything shaped like it, without being told where to look: a
reaction whose product carries an inhibiting edge back to it is product
inhibition, wherever it occurs. On the bundled example the search returns
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

A cut vertex -- a molecule whose removal splits
its component -- is meaningless for a list and central for a network. In a
trans-omic network these are the places where one layer's response reaches
another through a single route.

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

"A changed enzyme and a changed metabolite meet
at this reaction" is the claim the whole catalogue rests on, and it needs a
baseline: a network with this many changed molecules produces some convergence
for free.

:func:`~transnet.convergence_significance` holds
the network and the number of changed molecules per layer fixed, shuffles
*which* molecules changed, and recounts. It returns the observed count, the
null mean, a z-score and a permutation p-value, plus the null distribution so
it can be drawn.

On the brown adipocyte study, 170 reactions carry both axes against a null of
19 (z = +4.4, p = 0.015): the convergence is a property of the response, not of
the network's degree distribution. On the Uematsu panel, 14 against 6.3
(z = +1.9, p = 0.05) -- real, but close enough to the boundary that the number
matters.

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
     - per-pathway balance -- per-pathway balance
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
     - communities, and whether they span layers
   * - :func:`~transnet.visualization.plot_network_metrics`
     - degree distribution and component sizes

Colours come from one validated palette (:mod:`transnet.visualization.palette`):
red and blue only ever mean direction, layers and axes take hues that are neither.

Statistics reported against a baseline
--------------------------------------

A share of correct directional predictions (signed regulatory paths, influence) is compared with the
50% a coin flip achieves, with a binomial test; an enrichment (metabolite regulatory roles, transcription-factor activity) is
reported as not over-represented when its test says so. Hubs (trans-omic hubs) are ranked by
degree *within the responsive network*, as :ref:`Morita et al. <ref-morita2025>` define them: over the
whole network, transcription factors with thousands of ChIP-Atlas targets
dominate regardless of the data.

.. _comparing-conditions:

Comparing conditions
--------------------

:func:`~transnet.compare_transomic_networks` compares two networks
layer-by-layer and *relationship-by-relationship*. A generic differential
network comparison reports which nodes and edges differ; this reports which
*kinds of regulation* were gained and lost. A condition that loses its
transcriptional arm but keeps its allosteric one is a different biological story
from the reverse, and only an edge-type-aware comparison distinguishes them.

:Example: ``notebooks/walkthroughs/compare_conditions.py``
