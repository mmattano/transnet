.. _brown-adipocytes:

Brown adipocytes: thermogenic lipolysis
=======================================

The non-shivering cold response, measured across three layers in the same
cells. This page reports what the catalogue produces on a study with fold
changes in every layer.

.. code-block:: bash

    python notebooks/studies/brown_adipocytes.py

Every number and figure below comes from that notebook; rerun it and
``python maintenance/refresh_doc_figures.py`` to bring this page up to date.

.. tip::

   This page is the summary. :doc:`The full notebook <studies/brown_adipocytes>`
   runs every analysis in the catalogue on this data, with all of its tables and
   figures: the factor analyses read through the network, the regulatory motifs,
   the structural vulnerability and the convergence null.


The data
--------

Immortalised murine brown adipocytes stimulated with norepinephrine, six
replicates per timepoint: transcriptome (13,182 genes) and proteome (3,321
proteins) at 0, 4 and 24 h, metabolome (171 features) at seven points between
0 and 24 h. These are the measurements of
:ref:`Anagho-Mattanovich et al. 2025 <ref-anaghomattanovich2025>`, whose own
analysis compared five ways of integrating them: post-hoc correlation, PLS
(DIABLO), factor analysis (MOFA), trans-omics, and knowledge graphs.

They are cultured cells, not animals. The contrast below is 24 h of sustained
stimulation against unstimulated cells at :math:`q \leq 0.05`.

The measurements map onto the organism-wide mouse network built by
``maintenance/build_networks.py --organisms mouse --brenda``: 106,859
molecules and 1.9 million typed relationships, **82% of which cross between
layers**. Transcripts reach it through Ensembl-to-Entrez mapping (12,012
resolved) and metabolites by name against the network's compound synonyms
(119 of 171).

What the network adds to the lists
----------------------------------

**6,568 molecules respond, joined by 21,383 relationships.** A per-layer
analysis produces the first number. The second counts the places where one
layer's response reaches another.

.. figure:: figures/brown_adipocytes/network.png
   :width: 100%

   The most-regulated reactions, drawn in the order regulation flows: gene,
   protein, reaction, metabolite. Node fill is the measured change; dotted
   edges are allosteric regulation.

Which axis regulates each reaction
----------------------------------

**2,393 reactions are regulated**: 1,902 through enzyme amount only, 229
through metabolites only, and 262 through both. In **154** the enzyme axis and
the metabolite axis point opposite ways, so reading either layer alone gives the
wrong answer.

.. figure:: figures/brown_adipocytes/regulation_axes.png
   :width: 70%

As proportions: **90% of the regulated reactions are gene-driven, 21%
metabolite-driven**, 6.4% controversial. This response runs through enzyme
amount. In the trained-rat tissues the two axes are near-balanced.

Of the changed metabolites, 21 of 38 act on an enzyme, against 52% of all
measured metabolites (q = 0.55). The named regulators are worth following
individually; the set as a whole shows no enrichment.

Coverage limits all of this, and the notebook prints it first: 13,002 of
57,540 transcripts on the network are measured and 6,882 respond, against 118 of
19,652 metabolites measured and 38 responding. Every metabolite-axis statement
rests on those 38.

Is each protein change transcriptional?
---------------------------------------

Of the genes measured in both layers, 165 are concordant, 170 change at the
protein level only, 746 at the transcript level only, and **126 move against
their transcript**. Tested directly, as protein change minus transcript change
with standard errors, **292 of 1,661 testable proteins moved significantly
further than their transcript**. That is evidence for post-transcriptional
regulation, provided the two platforms report fold changes on comparable
scales.

.. figure:: figures/brown_adipocytes/concordance.png
   :width: 100%

Which molecules hold the response together
------------------------------------------

Two questions about the same subnetwork, with different answers.

**Cross-layer hubs** are the molecules whose neighbours span most layers. Here
they are led by the RNA helicase Ddx21 and the transcriptional co-regulator
Hcfc1, each reaching two layers over thousands of edges, and by the
2-oxoglutarate dehydrogenase component Ogdh, which touches three: transcript,
protein and the TCA reactions it catalyses. Ogdh is the only metabolic enzyme in
the top ten, and the only one whose hub status comes from crossing layers rather
than from interactome density.

.. figure:: figures/brown_adipocytes/structural_vulnerability.png
   :width: 90%

**Structural vulnerability** asks which molecules, if removed, break the
response into disconnected pieces. Removing Ddx21 fragments it into 668 pieces
and disconnects 11% of it; Smc1a gives 131 pieces. Below those the drop is
steep: Hk2 and Ak1, the glycolytic and adenylate-kinase entries, split off two
or three molecules each. A handful of highly connected regulatory proteins hold
this response together, not its metabolic core.

Does the upper hierarchy predict the metabolites?
-------------------------------------------------

Propagating the changed genes and proteins forward along signed edges predicts
a direction for each metabolite they reach. **24 of 32 changed metabolites are
predicted correctly, 75%, binomial p = 0.007.** It is the one directional claim
in this study that beats chance, and path tracing over the same network does
not.

Timing on the network
---------------------

The seven-point metabolome gives 118 metabolites a half-response time. There is
**no association between connectivity and response time** (p = 0.227), where
:ref:`Morita et al. <ref-morita2025>` report one in liver. The layer ordering
holds, and that is the ordering the published analysis emphasises.

4 h against 24 h
----------------

The early and late responses share a fifth of their regulatory relationships
(edge Jaccard 0.20). A typed comparison says which kinds of regulation differ
between them, not only how many molecules moved.

Factors, read through the network
---------------------------------

The published analysis fits MOFA factors to these data. A joint factor model
treats ``0h_01`` as the same culture in every layer, and that is testable: within
a timepoint, a culture whose transcript of a gene runs high should tend to have
that protein high too. Per-gene agreement is 0.067 as listed against -0.001 under
shuffled pairings, p = 0.002, over 1,656 genes measured in both layers. The
layers are the same cultures, so the joint model is sound.

.. figure:: figures/brown_adipocytes/factor_scores.png
   :width: 100%

   Five factors across the 18 samples measured in all three layers, by
   timepoint.

Five NMF factors, read four ways:

.. list-table::
   :header-rows: 1
   :widths: 12 20 22 24 22

   * - Factor
     - Layers agree
     - Retained when re-paired
     - Carried by
     - Network overlap vs null
   * - Factor3
     - 0.58
     - 100%
     - between-group
     - 0.139 vs 0.052
   * - Factor5
     - 0.48
     - 92%
     - between-group
     - 0.112 vs 0.053
   * - Factor1
     - 0.20
     - --
     - single-layer
     - 0.060 vs 0.053
   * - Factor2
     - 0.23
     - --
     - single-layer
     - 0.147 vs 0.056
   * - Factor4
     - 0.03
     - --
     - single-layer
     - 0.119 vs 0.053

**Factors 3 and 5 pass every reading.** The layers agree about them, the
agreement comes from the timepoints, and their features from different layers
converge on the network more than random features of the same layers
(q = 0.025). Those two are the ones to interpret.

.. figure:: figures/brown_adipocytes/factor_overview.png
   :width: 100%

.. figure:: figures/brown_adipocytes/network_Factor3.png
   :width: 100%

   Factor3 placed on the network: where its transcripts, proteins and
   metabolites meet. The figure is the reason the factor is interpretable --
   the loadings are named network nodes, not feature indices.

Factor 2 is the case worth studying. Its features converge on the network more
strongly than any other, onto branched-chain amino acid catabolism (*Bckdhb*,
3-hydroxyisobutyrate) and the complex III assembly factor *Uqcc4*. Branched-chain
amino acids are a known thermogenic fuel in brown fat. But the layers do not
agree about which samples carry the factor, so it is one layer's view of coherent
biochemistry. Only the network reading separates those two cases.

Motifs, and whether the convergence is real
-------------------------------------------

Two readings of the wiring rather than of the data through it.

.. figure:: figures/brown_adipocytes/regulatory_motifs.png
   :width: 70%

**Motifs.** Searching the responsive network for regulatory shapes finds one
instance of product inhibition and five of allosteric feedback, all in methionine
and cysteine metabolism with L-methionine as the effector. The 228 feed-forward
motifs require both the transcription factor and its target enzyme to respond.
The TF inference is permissive on this data, with 225 factors implicated, so
those 228 are a candidate list.

.. figure:: figures/brown_adipocytes/convergence_null.png
   :width: 80%

**Is the convergence real?** 170 reactions carry both a changed enzyme and a
changed metabolite. Holding the network and the number of changed molecules per
layer fixed and shuffling which molecules changed gives 19 (z = +4.4, p = 0.015).
Every axis attribution above depends on that convergence being a property of the
response rather than of the network's shape, and this is the test of it.

Compared with the published analysis
------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 44 56

   * - The paper concluded
     - What the network reading shows
   * - Three metabolic states: uninduced (0 h), active lipolysis (4 h),
       sustained induction (24 h)
     - clustering the metabolite trajectories over all seven timepoints
       separates three shapes; 88 of 171 metabolites change across the course
       (one-way ANOVA, BH-corrected), 70 of them monotonically
   * - The metabolome moves before the transcriptome and proteome
     - the metabolite response is already significant at 4 h while the protein
       changes concentrate at 24 h
   * - Lipid-droplet genes (*Sqle*, *Fdft1*) fall as transcripts while their
       proteins rise
     - both land in the ``discordant`` category of :ref:`transcript-protein concordance <concordance>` without
       being looked for, as 126 genes do
   * - Upper glycolysis (*Pfkl*, *Pfkp*) down, lower glycolysis (*Pklr*) up
     - reported gene by gene in the notebook's glycolysis table
   * - Trans-omics "can highlight specific elements that might otherwise be
       overlooked", complementing the other four methods
     - 154 controversial reactions, each naming the enzymes and the metabolites
       pulling against each other

The paper reached its trans-omic observations, the Rock2/protamine link and the
glycolytic split, by reading a network by hand. The conclusions agree. Here they
arrive as ordinary output, with the axis attribution and the threshold-free
transcript-protein test attached.

What does not work on this data
-------------------------------

* **Signed paths** (:ref:`signed regulatory paths <signed-paths>`) do not beat
  chance: 19 of 28 changed metabolites predicted in the right direction, 68%,
  binomial p = 0.087. The figure below shows the paths one by one, which is how
  a chance-level rate should be judged.

  .. figure:: figures/brown_adipocytes/regulatory_paths.png
     :width: 100%

* **Changed metabolites are not enriched for regulators** (21 of 38 against a
  background of 52%, q = 0.55) -- the numbers are above, under the regulation
  axes. With BRENDA included, most measured metabolites are an annotated
  effector of something, so the comparison has little room to show enrichment.
* **Response time does not follow connectivity** (p = 0.227, n = 118).
* **Transcription-factor inference is too permissive here**: 225 of 703 factors
  reach q <= 0.05. With 6,175 responding transcripts and thousands of ChIP-Atlas
  targets per factor, the ranking is informative and the count is not. See
  :doc:`kokaji_study` for the same inference scored against a published one.

.. seealso::

   :doc:`motrpac_study` runs the same catalogue across six tissues.
   :doc:`published_study` and :doc:`kokaji_study` check it against published
   conclusions.
