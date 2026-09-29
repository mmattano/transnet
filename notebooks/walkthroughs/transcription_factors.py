# %% [markdown]
# # Transcription factors
#
# A changed transcript is an observation. The factor that changed it is an
# inference, and it is the one cross-layer question a transcriptome alone cannot
# answer. This walkthrough builds that inference from the network's
# `transcriptional_regulation` edges and shows what it can and cannot support.

# %%
import matplotlib.pyplot as plt
import pandas as pd

from transnet import (
    load_example_network,
    load_example_omics,
    map_omics_to_network,
    transcription_factor_activity,
)

graph = load_example_network()
map_omics_to_network(graph, load_example_omics(), id_column="id",
                     log2fc_column="log2FC", qvalue_column="padj")
graph

# %% [markdown]
# ## The edges
#
# A transcription factor is a Proteome node with outgoing
# `transcriptional_regulation` edges to Transcriptome nodes. The edges come from
# ChIP-Atlas, which records where a factor was observed bound, so `confidence`
# is a binding score and not a claim about activation or repression. The sign is
# therefore 0: a bound factor can do either.

# %%
pd.DataFrame(
    [
        {"factor": graph.nodes[u].get("symbol", u),
         "target": graph.nodes[v].get("symbol", v),
         "sign": d.get("sign"), "confidence": d.get("confidence")}
        for u, v, d in graph.edges(data=True)
        if d["edge_type"] == "transcriptional_regulation"
    ]
)

# %% [markdown]
# ## Two questions per factor
#
# `transcription_factor_activity` asks both:
#
# 1. **Is the factor implicated?** Its measured targets are tested for
#    enrichment among the responsive genes, one-sided Fisher, Benjamini-Hochberg
#    across factors.
# 2. **Which way?** If its responsive targets moved mostly one way, that is the
#    inferred activity, with a binomial test on the imbalance.
#
# `min_targets` sets how few measured targets is too few to test. Three is low
# for a real network and appropriate here, where the example has four factors
# with three targets each.

# %%
activity = transcription_factor_activity(graph, min_targets=2)
activity[["name", "n_targets", "n_responsive_targets", "n_up", "n_down",
          "odds_ratio", "q_value", "inferred_activity", "factor_regulated"]]

# %% [markdown]
# ## The factor that did not change
#
# `factor_regulated` is the factor's own measured direction. Foxo1 here has all
# three of its targets down while its own level is unchanged, so an
# expression-only reading of the proteome would call it inactive. Foxo1 is
# regulated by phosphorylation, which moves it out of the nucleus without
# changing how much of it there is. This is the usual case for signalling-
# controlled factors, and inference from targets is what finds it.

# %%
activity.set_index("name").loc[["Foxo1", "Srebf1"],
                               ["n_up", "n_down", "inferred_activity", "factor_regulated"]]

# %% [markdown]
# With four factors and three targets each, nothing survives correction
# (q = 0.12 at best). The direction is still readable, and the figure shows both:
# factors that pass are named in ink, the rest in grey.

# %%
from transnet.visualization import plot_tf_activity

plot_tf_activity(activity, q_threshold=0.05)
plt.show()

# %% [markdown]
# ## Stricter targets
#
# `min_confidence` raises the ChIP-Atlas binding score a target must have, which
# trades coverage for specificity. On an organism-wide network it is the main
# control over how much of the transcriptome a factor appears to own.
#
# The bundled example carries no binding scores, and an edge with no score
# cannot meet a threshold, so asking for one here removes every edge. That is
# the intended behaviour and worth seeing: a filtered result on an unscored
# network is empty rather than unfiltered.

# %%
for threshold in (None, 100, 500):
    stricter = transcription_factor_activity(graph, min_targets=2,
                                             min_confidence=threshold)
    targets = int(stricter["n_targets"].sum()) if len(stricter) else 0
    print(f"min_confidence={str(threshold):<5} {len(stricter)} factors testable, "
          f"{targets} measured targets between them")

# %% [markdown]
# ## Where this breaks
#
# On an organism-wide network the test becomes too permissive in one direction
# and too sparse in the other, and both failures appear in the studies:
#
# * **Mouse.** ChIP-Atlas lists thousands of targets for a well-studied factor.
#   When 6,175 transcripts respond, 225 of 703 factors reach q <= 0.05, which
#   says more about target-list size than about the biology. The ranking is
#   informative; the count is not.
# * **Rat.** Only 32 factors have enough `rn6` coverage to be testable at all,
#   and four of six MoTrPAC tissues implicate none.
#
# Neither is a property of the method. Both are properties of the annotation, so
# the honest reading states the coverage first. `notebooks/studies/kokaji_liver.py`
# scores this inference against a published one on the same data, which is the
# only way to know whether the ranking is right.
#
# The live ChIP-Atlas calls that build these edges, `get_chip_tf_targets` and
# `list_chip_tfs`, are in `external_annotation.py`, which needs the network.
