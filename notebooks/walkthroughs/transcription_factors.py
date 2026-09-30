# %% [markdown]
# # Transcription factors
#
# A changed transcript is measured. The transcription factor that changed it
# is not, and has to be inferred. The network's `transcriptional_regulation`
# edges say which factors bind near which genes, so a factor whose target genes
# changed more than expected is a candidate driver of the response. This
# walkthrough makes that inference and shows where it is reliable and where it
# is not.

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
# `min_targets` sets the smallest number of measured targets a factor needs to
# be tested. The example has four factors with three targets each, so it is set
# to 2 here; on an organism-wide network use a larger value.

# %%
activity = transcription_factor_activity(graph, min_targets=2)
activity[["name", "n_targets", "n_responsive_targets", "n_up", "n_down",
          "odds_ratio", "q_value", "inferred_activity", "factor_regulated"]]

# %% [markdown]
# ## The factor that did not change
#
# `factor_regulated` is the factor's own measured change. All three Foxo1
# targets went down while the amount of Foxo1 did not change, so looking at the
# proteome alone would suggest Foxo1 was not involved. But Foxo1 is controlled
# by phosphorylation, which moves it out of the nucleus without changing how
# much of it there is. This is typical of factors controlled by signalling, and
# inference from the targets is what detects it.

# %%
activity.set_index("name").loc[["Foxo1", "Srebf1"],
                               ["n_up", "n_down", "inferred_activity", "factor_regulated"]]

# %% [markdown]
# With four factors and three targets each, no factor is significant after
# correcting for multiple testing (the best q-value is 0.12). The inferred
# direction can still be read. In the figure, significant factors are labelled
# in black and the rest in grey.

# %%
from transnet.visualization import plot_tf_activity

plot_tf_activity(activity, q_threshold=0.05)
plt.show()

# %% [markdown]
# ## Stricter targets
#
# `min_confidence` sets the minimum ChIP-Atlas binding score for a target to
# count. A higher value keeps fewer, more reliable targets. On an
# organism-wide network it is the main control over how many genes each factor
# is assigned.
#
# The bundled example has no binding scores, and an edge without a score cannot
# pass a threshold, so any threshold removes every edge here. This is
# deliberate: filtering a network that has no scores gives an empty result, not
# an unfiltered one.

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
# On organism-wide networks the test fails in two opposite ways, and both
# appear in the studies:
#
# * **Too many hits (mouse).** ChIP-Atlas lists thousands of targets for a
#   well-studied factor. When 6,175 transcripts respond, 225 of 703 factors
#   reach q <= 0.05. That number reflects the size of the target lists more
#   than the biology. The ranking of factors is informative; the count is not.
# * **Too few testable factors (rat).** Only 32 factors have enough ChIP-Atlas
#   data for the rat genome (`rn6`) to be tested at all, and four of the six
#   MoTrPAC tissues implicate none.
#
# Both limits come from the annotation, not from the method, so report the
# coverage alongside any result. The liver time-course study
# (`notebooks/studies/liver_timecourse.py`) compares this inference with a
# published one on the same data, which is the only way to check whether the
# ranking is right.
#
# The live ChIP-Atlas calls that build these edges, `get_chip_tf_targets` and
# `list_chip_tfs`, are in `external_annotation.py`, which needs the network.
