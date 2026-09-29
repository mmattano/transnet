"""Reactome Pathway Database API client.

Reactome (https://reactome.org) is a curated, peer-reviewed database of
biological pathways covering reactions, complexes, and physical entities
across >20 species.

This module wraps the Reactome ContentService REST API and the Analysis
Service.

Key functions
-------------
reactome_get_pathways
    All top-level pathways for a species.
reactome_get_pathway_reactions
    Reactions (events) contained in a pathway.
reactome_get_pathway_entities
    Proteins and metabolites participating in a pathway.
reactome_map_ids_to_pathways
    Map a list of gene/protein/metabolite IDs to their Reactome pathways.
reactome_over_representation
    Pathway over-representation analysis for a gene list.
reactome_get_pathway_hierarchy
    Full pathway hierarchy tree for a species.

References
----------
https://reactome.org/ContentService/
https://reactome.org/AnalysisService/
"""

from __future__ import annotations

import logging
import time
from typing import Dict, List, Optional

import pandas as pd
import requests


logger = logging.getLogger(__name__)

__all__ = [
    "reactome_get_pathways",
    "reactome_get_pathway_reactions",
    "reactome_get_pathway_entities",
    "reactome_map_ids_to_pathways",
    "reactome_over_representation",
    "reactome_get_pathway_hierarchy",
]

_CONTENT_API = "https://reactome.org/ContentService/data"
_ANALYSIS_API = "https://reactome.org/AnalysisService"
_REQUEST_DELAY = 0.3  # seconds between calls

# Reactome uses full species names; map common NCBI IDs for convenience
_NCBI_TO_REACTOME_SPECIES: Dict[str, str] = {
    "9606": "Homo sapiens",
    "10090": "Mus musculus",
    "10116": "Rattus norvegicus",
    "4932": "Saccharomyces cerevisiae",
    "83333": "Escherichia coli",
    "6239": "Caenorhabditis elegans",
    "7227": "Drosophila melanogaster",
    "7955": "Danio rerio",
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _get(url: str, params: Optional[dict] = None, timeout: int = 30) -> dict | list:
    """GET a Reactome endpoint, raise on HTTP error, return JSON."""
    headers = {"Accept": "application/json"}
    response = requests.get(url, params=params, headers=headers, timeout=timeout)
    response.raise_for_status()
    time.sleep(_REQUEST_DELAY)
    return response.json()


def _post(url: str, data: str, params: Optional[dict] = None, timeout: int = 60) -> dict | list:
    """POST to a Reactome endpoint with plain-text body, return JSON."""
    headers = {"Accept": "application/json", "Content-Type": "text/plain"}
    response = requests.post(url, data=data.encode(), params=params, headers=headers, timeout=timeout)
    response.raise_for_status()
    time.sleep(_REQUEST_DELAY)
    return response.json()


def _resolve_species(species: str) -> str:
    """Return the Reactome species name for a NCBI taxon ID or species name."""
    return _NCBI_TO_REACTOME_SPECIES.get(str(species), species)


# ---------------------------------------------------------------------------
# Public API functions
# ---------------------------------------------------------------------------

def reactome_get_pathways(species: str = "9606") -> pd.DataFrame:
    """Retrieve all top-level pathways for a species from Reactome.

    Parameters
    ----------
    species : str
        NCBI taxonomy ID (e.g. ``"9606"``) or Reactome species name
        (e.g. ``"Homo sapiens"``).

    Returns
    -------
    pd.DataFrame
        Columns: stId, displayName, name, isInDisease, hasEHLD,
        releaseDate, schemaClass, isInferred, species.
    """
    species_name = _resolve_species(species)
    url = f"{_CONTENT_API}/pathways/top/{species_name}"
    logger.info(f"Fetching top-level Reactome pathways for {species_name}")
    try:
        data = _get(url)
        records = [
            {
                "stId": p.get("stId"),
                "displayName": p.get("displayName"),
                "name": p.get("name", [None])[0] if isinstance(p.get("name"), list) else p.get("name"),
                "isInDisease": p.get("isInDisease", False),
                "hasEHLD": p.get("hasEHLD", False),
                "releaseDate": p.get("releaseDate"),
                "schemaClass": p.get("schemaClass"),
                "isInferred": p.get("isInferred", False),
                "species": species_name,
            }
            for p in (data if isinstance(data, list) else [])
        ]
        df = pd.DataFrame(records)
        logger.info(f"Retrieved {len(df)} top-level pathways for {species_name}")
        return df
    except requests.HTTPError as exc:
        logger.error(f"Reactome pathways request failed: {exc}")
        return pd.DataFrame()


def reactome_get_pathway_reactions(pathway_id: str) -> pd.DataFrame:
    """Get all reactions/events contained in a Reactome pathway.

    Parameters
    ----------
    pathway_id : str
        Reactome stable identifier (e.g. ``"R-HSA-1643685"``).

    Returns
    -------
    pd.DataFrame
        Columns: stId, displayName, schemaClass, isInDisease,
        isInferred, releaseDate, order.
    """
    url = f"{_CONTENT_API}/pathway/{pathway_id}/containedEvents"
    logger.debug(f"Fetching reactions for pathway {pathway_id}")
    try:
        data = _get(url)
        records = []
        for i, event in enumerate(data if isinstance(data, list) else []):
            records.append({
                "stId": event.get("stId"),
                "displayName": event.get("displayName"),
                "schemaClass": event.get("schemaClass"),
                "isInDisease": event.get("isInDisease", False),
                "isInferred": event.get("isInferred", False),
                "releaseDate": event.get("releaseDate"),
                "order": i,
                "pathway_id": pathway_id,
            })
        return pd.DataFrame(records)
    except requests.HTTPError as exc:
        logger.warning(f"Reactome reactions request failed for {pathway_id}: {exc}")
        return pd.DataFrame()


def reactome_get_pathway_entities(pathway_id: str) -> pd.DataFrame:
    """Get the physical entities (proteins, metabolites, complexes) in a pathway.

    Parameters
    ----------
    pathway_id : str
        Reactome stable identifier.

    Returns
    -------
    pd.DataFrame
        Columns: stId, displayName, schemaClass, identifier,
        databaseName, pathway_id.
        ``identifier`` / ``databaseName`` are the cross-reference
        (e.g. UniProt / CHEBI) when available.
    """
    url = f"{_CONTENT_API}/pathway/{pathway_id}/participatingPhysicalEntities"
    logger.debug(f"Fetching entities for pathway {pathway_id}")
    try:
        data = _get(url)
        records = []
        for entity in (data if isinstance(data, list) else []):
            # Cross-reference to external DB (UniProt, ChEBI, etc.)
            xrefs = entity.get("crossReference", []) or []
            identifier = xrefs[0].get("identifier") if xrefs else None
            db_name = xrefs[0].get("databaseName") if xrefs else None

            records.append({
                "stId": entity.get("stId"),
                "displayName": entity.get("displayName"),
                "schemaClass": entity.get("schemaClass"),
                "identifier": identifier,
                "databaseName": db_name,
                "pathway_id": pathway_id,
            })
        return pd.DataFrame(records)
    except requests.HTTPError as exc:
        logger.warning(f"Reactome entities request failed for {pathway_id}: {exc}")
        return pd.DataFrame()


def reactome_map_ids_to_pathways(
    identifiers: list,
    species: str = "9606",
) -> pd.DataFrame:
    """Map a list of gene/protein/metabolite IDs to Reactome pathways.

    Uses the Reactome ``/data/mapping`` endpoint one ID at a time and
    aggregates the results.  For bulk ORA see :func:`reactome_over_representation`.

    Parameters
    ----------
    identifiers : list
        Gene symbols, UniProt IDs, CHEBI IDs, or other cross-reference IDs.
    species : str
        NCBI taxonomy ID or Reactome species name.

    Returns
    -------
    pd.DataFrame
        Columns: identifier, pathway_stId, pathway_name, pathway_top_level.
    """
    species_name = _resolve_species(species)
    logger.info(f"Mapping {len(identifiers)} identifiers to Reactome pathways")
    records = []

    for ident in identifiers:
        url = f"{_CONTENT_API}/mapping/UniProt/{ident}/pathways"
        try:
            data = _get(url)
            for pw in (data if isinstance(data, list) else []):
                records.append({
                    "identifier": ident,
                    "pathway_stId": pw.get("stId"),
                    "pathway_name": pw.get("displayName"),
                    "isInDisease": pw.get("isInDisease", False),
                    "isInferred": pw.get("isInferred", False),
                    "species": species_name,
                })
        except requests.HTTPError as exc:
            # 404 = no pathways for this identifier — normal, skip silently
            if exc.response is not None and exc.response.status_code != 404:
                logger.debug(f"No Reactome pathway mapping for {ident}: {exc}")

    df = pd.DataFrame(records)
    logger.info(
        f"Mapped {df['identifier'].nunique() if not df.empty else 0} identifiers "
        f"to {df['pathway_stId'].nunique() if not df.empty else 0} pathways"
    )
    return df


def reactome_over_representation(
    gene_ids: list,
    species: str = "9606",
    p_value: float = 0.05,
    include_disease: bool = False,
) -> pd.DataFrame:
    """Pathway over-representation analysis using the Reactome Analysis Service.

    Submits a list of identifiers to the Reactome ORA endpoint and returns
    significantly enriched pathways ranked by p-value.

    Parameters
    ----------
    gene_ids : list
        Gene symbols, UniProt IDs, ENSEMBL IDs, or other identifiers.
    species : str
        NCBI taxonomy ID or Reactome species name.
    p_value : float
        P-value cut-off for the returned results.
    include_disease : bool
        If False (default), exclude disease-specific pathways.

    Returns
    -------
    pd.DataFrame
        Columns: stId, name, entities_found, entities_total,
        entities_ratio, entities_pValue, entities_fdr,
        reactions_found, reactions_total, reactions_ratio,
        species_name.
    """
    species_name = _resolve_species(species)
    id_string = "\n".join(str(g) for g in gene_ids)

    url = f"{_ANALYSIS_API}/identifiers/"
    params = {
        "interactors": "false",
        "pageSize": 1000,
        "page": 1,
        "sortBy": "ENTITIES_PVALUE",
        "order": "ASC",
        "resource": "TOTAL",
        "pValue": p_value,
        "includeDisease": str(include_disease).lower(),
        "species": species_name,
    }

    logger.info(
        f"Running Reactome ORA for {len(gene_ids)} identifiers in {species_name}"
    )
    try:
        data = _post(url, id_string, params=params)
    except requests.HTTPError as exc:
        logger.error(f"Reactome ORA request failed: {exc}")
        return pd.DataFrame()

    pathways = data.get("pathways", []) if isinstance(data, dict) else []
    records = []
    for pw in pathways:
        ent = pw.get("entities", {})
        rxn = pw.get("reactions", {})
        sp_list = pw.get("species", {})
        records.append({
            "stId": pw.get("stId"),
            "name": pw.get("name"),
            "entities_found": ent.get("found"),
            "entities_total": ent.get("total"),
            "entities_ratio": ent.get("ratio"),
            "entities_pValue": ent.get("pValue"),
            "entities_fdr": ent.get("fdr"),
            "reactions_found": rxn.get("found"),
            "reactions_total": rxn.get("total"),
            "reactions_ratio": rxn.get("ratio"),
            "species_name": sp_list.get("name") if isinstance(sp_list, dict) else species_name,
        })

    df = pd.DataFrame(records)
    logger.info(
        f"ORA complete: {len(df)} pathways at p≤{p_value} "
        f"(token: {data.get('summary', {}).get('token', 'n/a') if isinstance(data, dict) else 'n/a'})"
    )
    return df


def reactome_get_pathway_hierarchy(species: str = "9606") -> list:
    """Return the full pathway hierarchy tree for a species.

    Parameters
    ----------
    species : str
        NCBI taxonomy ID or Reactome species name.

    Returns
    -------
    list of dict
        Nested list matching the Reactome hierarchy JSON structure.
        Each node has ``stId``, ``displayName``, ``children`` (list).
    """
    species_name = _resolve_species(species)
    url = f"{_CONTENT_API}/eventsHierarchy/{species_name}"
    logger.info(f"Fetching Reactome hierarchy for {species_name}")
    try:
        data = _get(url)
        return data if isinstance(data, list) else []
    except requests.HTTPError as exc:
        logger.error(f"Reactome hierarchy request failed: {exc}")
        return []
