"""Temporal and dose structure on the trans-omic network.

Time-course and dose-response data give every molecule two numbers -- how fast
it responds (t-half) and how sensitive it is (EC50).  Putting those on network
nodes turns them into statements about network architecture:

* do highly connected molecules respond *first*?  Morita et al. found a negative
  correlation between degree and t-half in healthy liver, meaning hubs lead the
  response, and found that correlation destroyed in obese liver even though the
  network's structure was unchanged -- structural robustness with temporal
  vulnerability;
* are neighbouring molecules co-regulated?  A bimodal distribution of neighbour
  correlations (peaks near +/-0.75) says the network is coherently driven; a
  unimodal distribution near zero says it is not;
* which subnetwork carries the fast, sensitive response?  Kawata et al. split
  the insulin trans-omic network by EC50 and t-half into selective subnetworks
  for induced versus basal insulin.

References
----------
Morita K, et al. Structural robustness and temporal vulnerability of the
starvation-responsive metabolic network in healthy and obese mouse liver.
*Science Signaling* 18, 2025.

Kawata K, et al. Trans-omic Analysis Reveals Selective Responses to Induced and
Basal Insulin across Signaling, Transcriptional, and Metabolic Networks.
*iScience* 7:212-229, 2018.
"""

from typing import Dict, Optional, Sequence

import logging

import numpy as np
import pandas as pd

from transnet.biology.schema import available_layers

logger = logging.getLogger(__name__)

__all__ = [
    "assign_temporal_parameters",
    "temporal_network_structure",
    "split_by_response_class",
]


def _resolve_columns(requested, available, layer: str, kind: str) -> list:
    """Match requested column labels against a table's actual columns.

    Callers naturally write ``time_columns=[0, 10, 30, 60]`` for a table whose
    headers are the strings ``"0"``, ``"10"``, ... .  Matching on ``str`` as
    well as identity means both spellings work, instead of silently indexing by
    position and producing nonsense.
    """
    lookup = {str(column): column for column in available}
    resolved, missing = [], []
    for column in requested:
        if column in available:
            resolved.append(column)
        elif str(column) in lookup:
            resolved.append(lookup[str(column)])
        else:
            missing.append(column)

    if missing:
        logger.warning(
            f"{layer}: {kind} columns {missing} are not in the table "
            f"(available: {list(available)}); they were skipped"
        )
    return resolved


def _numeric_axis(columns) -> np.ndarray:
    """Numeric axis for a set of column labels.

    Uses the labels themselves when they are numeric -- so unevenly spaced
    timepoints ("0", "10", "30", "60") are treated as unevenly spaced -- and
    falls back to position when they are not.
    """
    try:
        return np.asarray([float(str(c)) for c in columns], dtype=float)
    except ValueError:
        return np.arange(len(columns), dtype=float)


def _half_response_time(times: np.ndarray, values: np.ndarray) -> Optional[float]:
    """Time at which the response first reaches half its peak amplitude.

    Amplitude is measured from the first timepoint, so a molecule that falls is
    handled the same way as one that rises.  Linear interpolation between the
    bracketing timepoints.
    """
    if len(times) < 2:
        return None

    baseline = values[0]
    deviations = values - baseline
    peak_index = int(np.argmax(np.abs(deviations)))
    peak = deviations[peak_index]
    if peak == 0:
        return None

    half = peak / 2.0
    for i in range(1, peak_index + 1):
        previous, current = deviations[i - 1], deviations[i]
        if (previous < half <= current) or (previous > half >= current):
            if current == previous:
                return float(times[i])
            fraction = (half - previous) / (current - previous)
            return float(times[i - 1] + fraction * (times[i] - times[i - 1]))
    return float(times[peak_index])


def _ec50(doses: np.ndarray, values: np.ndarray) -> Optional[float]:
    """Dose producing half the maximal response, by linear interpolation."""
    if len(doses) < 2:
        return None
    baseline = values[0]
    deviations = values - baseline
    peak_index = int(np.argmax(np.abs(deviations)))
    peak = deviations[peak_index]
    if peak == 0:
        return None
    half = peak / 2.0
    for i in range(1, peak_index + 1):
        previous, current = deviations[i - 1], deviations[i]
        if (previous < half <= current) or (previous > half >= current):
            if current == previous:
                return float(doses[i])
            fraction = (half - previous) / (current - previous)
            return float(doses[i - 1] + fraction * (doses[i] - doses[i - 1]))
    return float(doses[peak_index])


def assign_temporal_parameters(
    graph,
    timecourse: Dict[str, pd.DataFrame],
    time_columns: Sequence[float] = None,
    id_column: str = None,
    dose_response: Optional[Dict[str, pd.DataFrame]] = None,
    dose_columns: Sequence[float] = None,
) -> pd.DataFrame:
    """Compute t-half (and optionally EC50) per molecule and write them onto nodes.

    Parameters
    ----------
    graph : networkx.Graph
    timecourse : dict of str to pandas.DataFrame
        ``{layer: table}`` where each table has one row per molecule and one
        column per timepoint.
    time_columns : sequence, optional
        The timepoint columns, in order.  Labels are matched against the
        table's headers by value or by string, so ``[0, 10, 30, 60]`` works for
        headers written as either integers or strings.  The labels also supply
        the time axis when they are numeric, so unevenly spaced timepoints are
        handled correctly.  Defaults to every column after ``id_column``.
    id_column : str, optional
        Identifier column; defaults to the first column.
    dose_response : dict of str to pandas.DataFrame, optional
        Same shape as ``timecourse`` but across doses, for EC50.
    dose_columns : sequence, optional
        The dose columns, in order; matched the same way as ``time_columns``.

    Returns
    -------
    pandas.DataFrame
        Columns ``node``, ``layer``, ``t_half``, ``ec50``, ``degree`` -- one row
        per molecule that matched a node.  Nodes also gain ``t_half`` and
        ``ec50`` attributes.
    """
    nodes_by_layer: Dict[str, set] = {}
    for node, data in graph.nodes(data=True):
        layer = data.get("layer")
        if layer:
            nodes_by_layer.setdefault(str(layer), set()).add(str(node))

    def _parameters(tables, columns, fn, key):
        results = {}
        if not tables:
            return results
        for layer, table in tables.items():
            layer = str(layer)
            if table is None or len(table) == 0:
                continue
            if layer not in nodes_by_layer:
                logger.warning(
                    f"No '{layer}' layer in this network; skipping its {key} data"
                )
                continue
            ids = id_column if id_column and id_column in table.columns else table.columns[0]
            requested = (
                list(columns) if columns is not None
                else [c for c in table.columns if c != ids]
            )
            value_columns = _resolve_columns(requested, table.columns, layer, key)
            if not value_columns:
                continue
            axis = _numeric_axis(value_columns)

            for _, row in table.iterrows():
                node = str(row[ids])
                if node not in nodes_by_layer[layer]:
                    continue
                try:
                    values = np.asarray(
                        [float(row[c]) for c in value_columns], dtype=float
                    )
                except (TypeError, ValueError, KeyError):
                    continue
                if np.isnan(values).any():
                    continue
                results[node] = fn(axis, values)
        return results

    t_halves = _parameters(timecourse, time_columns, _half_response_time, "time-course")
    ec50s = _parameters(dose_response, dose_columns, _ec50, "dose-response")

    rows = []
    for node in set(t_halves) | set(ec50s):
        t_half = t_halves.get(node)
        ec50 = ec50s.get(node)
        graph.nodes[node]["t_half"] = t_half
        graph.nodes[node]["ec50"] = ec50
        rows.append({
            "node": node,
            "layer": graph.nodes[node].get("layer"),
            "t_half": t_half,
            "ec50": ec50,
            "degree": graph.degree(node),
        })

    table = pd.DataFrame(rows, columns=["node", "layer", "t_half", "ec50", "degree"])
    logger.info(
        f"Assigned temporal parameters to {len(table)} nodes "
        f"({table['t_half'].notna().sum()} with t_half, "
        f"{table['ec50'].notna().sum()} with EC50)"
    )
    return table


def temporal_network_structure(
    graph,
    timecourse_matrix: Optional[pd.DataFrame] = None,
    max_neighbour_pairs: int = 50000,
    alpha: float = 0.05,
) -> Dict[str, object]:
    """Test whether the network's wiring explains its response timing.

    Parameters
    ----------
    graph : networkx.Graph
        A network with ``t_half`` on its nodes, from
        :func:`assign_temporal_parameters`.
    timecourse_matrix : pandas.DataFrame, optional
        Molecules x timepoints, indexed by node id.  When given, neighbour
        co-regulation is computed from it.
    max_neighbour_pairs : int
        Cap on how many neighbouring pairs to correlate.
    alpha : float
        Significance level above which no direction is claimed for the
        degree/t-half association.

    Returns
    -------
    dict
        ``degree_vs_thalf`` : dict
            Spearman rho, p-value and n between degree and t-half, plus
            ``significant`` and an ``interpretation``.  A negative rho means
            hubs respond first, but a direction is only stated when
            ``p < alpha``; otherwise the interpretation says so.
        ``neighbour_correlation`` : dict
            Mean and distribution summary of the correlation between connected
            molecules, with ``fraction_strong`` (``|r| > 0.9``) -- the coherence
            signature.  Absent when no matrix is supplied.
        ``per_layer_thalf`` : pandas.DataFrame
            Median t-half per layer, showing the order in which layers respond.
    """
    from scipy import stats

    nodes = [
        n for n, d in graph.nodes(data=True) if d.get("t_half") is not None
    ]
    result: Dict[str, object] = {}

    if len(nodes) >= 3:
        degrees = np.array([graph.degree(n) for n in nodes], dtype=float)
        t_halves = np.array([graph.nodes[n]["t_half"] for n in nodes], dtype=float)
        rho, pvalue = stats.spearmanr(degrees, t_halves)
        # State a direction only when the correlation supports one. Reading
        # "hubs respond earlier" off the sign of a rho that is not
        # distinguishable from zero is how a null result gets reported as a
        # finding -- and with the handful of nodes that usually carry a
        # t-half, that is the common case.
        significant = bool(pvalue < alpha)
        if not significant:
            interpretation = (
                f"no association between degree and response time "
                f"(p = {pvalue:.3g}, n = {len(nodes)})"
            )
        elif rho < 0:
            interpretation = "hubs respond earlier"
        else:
            interpretation = "hubs respond later"
        result["degree_vs_thalf"] = {
            "spearman_rho": float(rho),
            "p_value": float(pvalue),
            "n": len(nodes),
            "significant": significant,
            "alpha": float(alpha),
            "interpretation": interpretation,
        }
        logger.info(
            f"degree vs t_half: rho={rho:.3f}, p={pvalue:.3g}, n={len(nodes)}"
        )
    else:
        result["degree_vs_thalf"] = {
            "spearman_rho": None, "p_value": None, "n": len(nodes),
            "significant": False, "alpha": float(alpha),
            "interpretation": "too few nodes with t_half",
        }

    rows = []
    for layer in available_layers(graph):
        values = [
            graph.nodes[n]["t_half"] for n in nodes
            if graph.nodes[n].get("layer") == layer
        ]
        if values:
            rows.append({
                "layer": layer,
                "n": len(values),
                "median_t_half": float(np.median(values)),
                "mean_t_half": float(np.mean(values)),
            })
    result["per_layer_thalf"] = pd.DataFrame(
        rows, columns=["layer", "n", "median_t_half", "mean_t_half"]
    ).sort_values("median_t_half").reset_index(drop=True)

    if timecourse_matrix is not None and len(timecourse_matrix):
        profiles = timecourse_matrix.dropna()
        available = set(profiles.index.astype(str))
        correlations = []
        for u, v in graph.edges():
            if len(correlations) >= max_neighbour_pairs:
                break
            u, v = str(u), str(v)
            if u in available and v in available:
                a = profiles.loc[u].to_numpy(dtype=float)
                b = profiles.loc[v].to_numpy(dtype=float)
                if a.std() and b.std():
                    correlations.append(float(np.corrcoef(a, b)[0, 1]))

        if correlations:
            values = np.asarray(correlations)
            result["neighbour_correlation"] = {
                "n_pairs": len(values),
                "mean": float(values.mean()),
                "mean_absolute": float(np.abs(values).mean()),
                "fraction_strong": float((np.abs(values) > 0.9).mean()),
                "fraction_positive": float((values > 0).mean()),
                "values": values,
            }
            logger.info(
                f"neighbour co-regulation: mean |r| = "
                f"{np.abs(values).mean():.3f} over {len(values)} pairs"
            )

    return result


def split_by_response_class(
    graph,
    t_half_threshold: Optional[float] = None,
    ec50_threshold: Optional[float] = None,
) -> Dict[str, object]:
    """Split the network into fast/slow and sensitive/insensitive subnetworks.

    Thresholds default to the median of each parameter across the nodes that
    have it, which is the data-driven split used when no external cutoff
    applies.

    Parameters
    ----------
    graph : networkx.Graph
    t_half_threshold, ec50_threshold : float, optional

    Returns
    -------
    dict
        ``thresholds`` : dict
            The cutoffs actually used.
        ``classes`` : pandas.DataFrame
            Per node: its t-half, EC50 and assigned class.
        ``subnetworks`` : dict
            ``{class_name: subgraph}`` for each populated class.
    """
    nodes = [
        n for n, d in graph.nodes(data=True)
        if d.get("t_half") is not None or d.get("ec50") is not None
    ]
    if not nodes:
        logger.warning(
            "No node carries t_half or ec50; run assign_temporal_parameters() first"
        )
        return {"thresholds": {}, "classes": pd.DataFrame(), "subnetworks": {}}

    t_values = [graph.nodes[n]["t_half"] for n in nodes if graph.nodes[n].get("t_half") is not None]
    e_values = [graph.nodes[n]["ec50"] for n in nodes if graph.nodes[n].get("ec50") is not None]

    if t_half_threshold is None and t_values:
        t_half_threshold = float(np.median(t_values))
    if ec50_threshold is None and e_values:
        ec50_threshold = float(np.median(e_values))

    rows = []
    for node in nodes:
        t_half = graph.nodes[node].get("t_half")
        ec50 = graph.nodes[node].get("ec50")
        speed = None
        if t_half is not None and t_half_threshold is not None:
            speed = "fast" if t_half <= t_half_threshold else "slow"
        sensitivity = None
        if ec50 is not None and ec50_threshold is not None:
            sensitivity = "sensitive" if ec50 <= ec50_threshold else "insensitive"
        label = "_".join([p for p in (sensitivity, speed) if p]) or "unclassified"
        rows.append({
            "node": node,
            "layer": graph.nodes[node].get("layer"),
            "t_half": t_half,
            "ec50": ec50,
            "response_class": label,
        })

    classes = pd.DataFrame(rows)
    subnetworks = {
        label: graph.subgraph(group["node"].tolist()).copy()
        for label, group in classes.groupby("response_class")
    }

    logger.info(
        "Response classes: "
        + ", ".join(
            f"{label} ({sub.number_of_nodes()} nodes)"
            for label, sub in subnetworks.items()
        )
    )
    return {
        "thresholds": {
            "t_half": t_half_threshold,
            "ec50": ec50_threshold,
        },
        "classes": classes,
        "subnetworks": subnetworks,
    }
