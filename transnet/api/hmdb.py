"""Human Metabolome Database (HMDB) API client.

HMDB (https://hmdb.ca) is the most comprehensive freely available
database of human metabolites.  This module provides functions to look
up individual metabolites, search by name or formula, retrieve pathway
associations, and map between identifier types.

HMDB does not offer a formal REST API.  This module uses:

1. The HMDB XML record endpoint (``hmdb.ca/metabolites/{id}.xml``) for
   individual metabolite data.
2. The metabolite search endpoint
   (``hmdb.ca/metabolites.xml?search_query={q}&search_field={f}``)
   for name / formula / InChIKey searches.

All responses are parsed from XML using the Python standard library
``xml.etree.ElementTree`` — no extra dependencies required.

Key functions
-------------
hmdb_get_metabolite
    Full record for one HMDB ID.
hmdb_search_metabolites
    Search HMDB by name, formula, or InChIKey.
hmdb_get_pathways
    Pathway associations for a metabolite.
hmdb_get_diseases
    Disease associations for a metabolite.
hmdb_map_ids
    Batch ID mapping (name → HMDB, KEGG → HMDB, ChEBI → HMDB, etc.).
hmdb_enrich_metabolites
    Populate Metabolite objects with HMDB data in-place.
"""

from __future__ import annotations

import logging
import time
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional

import pandas as pd
import requests


logger = logging.getLogger(__name__)

__all__ = [
    "hmdb_get_metabolite",
    "hmdb_search_metabolites",
    "hmdb_get_pathways",
    "hmdb_get_diseases",
    "hmdb_map_ids",
    "hmdb_enrich_metabolites",
]

_BASE_URL = "https://hmdb.ca"
_REQUEST_DELAY = 0.5  # seconds between calls

# Namespaces present in HMDB XML
_NS = {"hmdb": "http://www.hmdb.ca"}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _get_xml(url: str, params: Optional[dict] = None, timeout: int = 30) -> ET.Element:
    """Fetch an HMDB XML endpoint and return the parsed root element."""
    headers = {"Accept": "application/xml, text/xml"}
    resp = requests.get(url, params=params, headers=headers, timeout=timeout)
    resp.raise_for_status()
    time.sleep(_REQUEST_DELAY)
    return ET.fromstring(resp.content)


def _text(element: ET.Element, tag: str, ns: str = "hmdb") -> Optional[str]:
    """Return stripped text of the first matching child element, or None."""
    child = element.find(f"{ns}:{tag}", _NS) if ns else element.find(tag)
    if child is not None and child.text:
        return child.text.strip() or None
    return None


def _text_list(element: ET.Element, path: str) -> List[str]:
    """Return a list of stripped text values for all matching elements."""
    results = []
    for child in element.findall(path, _NS):
        if child.text and child.text.strip():
            results.append(child.text.strip())
    return results


def _parse_metabolite_record(metabolite_el: ET.Element) -> dict:
    """Parse a single ``<metabolite>`` XML element into a flat dict."""
    t = lambda tag: _text(metabolite_el, tag)

    # Synonyms
    synonyms = _text_list(metabolite_el, "hmdb:synonyms/hmdb:synonym")

    # Pathways
    pathways = []
    for pw in metabolite_el.findall("hmdb:pathways/hmdb:pathway", _NS):
        name = _text(pw, "name")
        kegg_map = _text(pw, "kegg_map_id")
        smpdb = _text(pw, "smpdb_id")
        pathways.append({
            "name": name,
            "kegg_map_id": kegg_map,
            "smpdb_id": smpdb,
        })

    # Diseases
    diseases = []
    for dis in metabolite_el.findall("hmdb:diseases/hmdb:disease", _NS):
        dis_name = _text(dis, "name")
        omim = _text(dis, "omim_id")
        diseases.append({"name": dis_name, "omim_id": omim})

    # Normal concentrations (biofluid)
    concentrations = []
    for conc in metabolite_el.findall(
        "hmdb:normal_concentrations/hmdb:concentration", _NS
    ):
        concentrations.append({
            "biofluid": _text(conc, "biofluid"),
            "concentration_value": _text(conc, "concentration_value"),
            "concentration_units": _text(conc, "concentration_units"),
            "subject_age": _text(conc, "subject_age"),
            "subject_sex": _text(conc, "subject_sex"),
        })

    # External DB cross-references
    return {
        "hmdb_id": t("accession"),
        "name": t("name"),
        "description": t("description"),
        "chemical_formula": t("chemical_formula"),
        "monoisotopic_molecular_weight": t("monoisotopic_molecular_weight"),
        "average_molecular_weight": t("average_molecular_weight"),
        "smiles": t("smiles"),
        "inchi": t("inchi"),
        "inchikey": t("inchikey"),
        "cas_registry_number": t("cas_registry_number"),
        "pubchem_compound_id": t("pubchem_compound_id"),
        "chebi_id": t("chebi_id"),
        "kegg_id": t("kegg_id"),
        "synonyms": synonyms,
        "pathways": pathways,
        "diseases": diseases,
        "normal_concentrations": concentrations,
        "status": t("status"),
        "origin": t("origin"),
        "super_class": t("super_class"),
        "class": t("class"),
        "sub_class": t("sub_class"),
    }


# ---------------------------------------------------------------------------
# Public API functions
# ---------------------------------------------------------------------------

def hmdb_get_metabolite(hmdb_id: str) -> dict:
    """Retrieve the full HMDB record for a metabolite.

    Parameters
    ----------
    hmdb_id : str
        HMDB accession number, e.g. ``"HMDB0000001"`` or ``"HMDB0000122"``.
        Both the legacy 5-digit (``HMDB00001``) and current 7-digit
        (``HMDB0000001``) formats are accepted.

    Returns
    -------
    dict
        Keys: hmdb_id, name, chemical_formula, monoisotopic_molecular_weight,
        average_molecular_weight, smiles, inchi, inchikey, cas_registry_number,
        pubchem_compound_id, chebi_id, kegg_id, synonyms (list),
        pathways (list of dict), diseases (list of dict),
        normal_concentrations (list of dict), origin, super_class, class,
        sub_class, status.
        Returns an empty dict on error.
    """
    # Normalise to current 7-digit format
    accession = _normalise_hmdb_id(hmdb_id)
    url = f"{_BASE_URL}/metabolites/{accession}.xml"
    logger.debug(f"Fetching HMDB record for {accession}")
    try:
        root = _get_xml(url)
        # Root may be <metabolite> or <hmdb><metabolite>
        if root.tag.endswith("metabolite"):
            metabolite_el = root
        else:
            metabolite_el = root.find("hmdb:metabolite", _NS)
            if metabolite_el is None:
                metabolite_el = root.find("metabolite")
        if metabolite_el is None:
            logger.warning(f"Could not parse HMDB XML for {accession}")
            return {}
        return _parse_metabolite_record(metabolite_el)
    except requests.HTTPError as exc:
        logger.warning(f"HMDB fetch failed for {accession}: {exc}")
        return {}
    except ET.ParseError as exc:
        logger.warning(f"HMDB XML parse error for {accession}: {exc}")
        return {}


def hmdb_search_metabolites(
    query: str,
    search_field: str = "name",
    max_results: int = 20,
) -> pd.DataFrame:
    """Search HMDB metabolites by name, formula, InChIKey, or CAS number.

    Parameters
    ----------
    query : str
        The search term.
    search_field : str
        One of ``"name"``, ``"formula"``, ``"inchikey"``, ``"cas"``,
        ``"chebi_id"``, ``"kegg_id"``, ``"pubchem_id"``.
    max_results : int
        Maximum number of results to return.

    Returns
    -------
    pd.DataFrame
        Columns: hmdb_id, name, chemical_formula, average_molecular_weight,
        smiles, inchikey, chebi_id, kegg_id, pubchem_compound_id.
    """
    url = f"{_BASE_URL}/metabolites.xml"
    params = {"search_query": query, "search_field": search_field}
    logger.info(f"HMDB search: '{query}' in field '{search_field}'")
    try:
        root = _get_xml(url, params=params)
    except requests.HTTPError as exc:
        logger.error(f"HMDB search failed: {exc}")
        return pd.DataFrame()
    except ET.ParseError as exc:
        logger.error(f"HMDB search XML parse error: {exc}")
        return pd.DataFrame()

    records = []
    for met in root.findall("hmdb:metabolite", _NS):
        records.append({
            "hmdb_id": _text(met, "accession"),
            "name": _text(met, "name"),
            "chemical_formula": _text(met, "chemical_formula"),
            "average_molecular_weight": _text(met, "average_molecular_weight"),
            "smiles": _text(met, "smiles"),
            "inchikey": _text(met, "inchikey"),
            "chebi_id": _text(met, "chebi_id"),
            "kegg_id": _text(met, "kegg_id"),
            "pubchem_compound_id": _text(met, "pubchem_compound_id"),
        })
        if len(records) >= max_results:
            break

    df = pd.DataFrame(records)
    logger.info(f"HMDB search returned {len(df)} results for '{query}'")
    return df


def hmdb_get_pathways(hmdb_id: str) -> pd.DataFrame:
    """Return pathway associations for an HMDB metabolite.

    Parameters
    ----------
    hmdb_id : str
        HMDB accession number.

    Returns
    -------
    pd.DataFrame
        Columns: hmdb_id, pathway_name, kegg_map_id, smpdb_id.
    """
    record = hmdb_get_metabolite(hmdb_id)
    if not record:
        return pd.DataFrame()
    pathways = record.get("pathways", [])
    if not pathways:
        return pd.DataFrame(
            columns=["hmdb_id", "pathway_name", "kegg_map_id", "smpdb_id"]
        )
    df = pd.DataFrame(pathways)
    df.insert(0, "hmdb_id", record.get("hmdb_id", hmdb_id))
    df = df.rename(columns={"name": "pathway_name"})
    return df


def hmdb_get_diseases(hmdb_id: str) -> pd.DataFrame:
    """Return disease associations for an HMDB metabolite.

    Parameters
    ----------
    hmdb_id : str
        HMDB accession number.

    Returns
    -------
    pd.DataFrame
        Columns: hmdb_id, disease_name, omim_id.
    """
    record = hmdb_get_metabolite(hmdb_id)
    if not record:
        return pd.DataFrame()
    diseases = record.get("diseases", [])
    if not diseases:
        return pd.DataFrame(columns=["hmdb_id", "disease_name", "omim_id"])
    df = pd.DataFrame(diseases)
    df.insert(0, "hmdb_id", record.get("hmdb_id", hmdb_id))
    df = df.rename(columns={"name": "disease_name"})
    return df


def hmdb_map_ids(
    ids: list,
    from_type: str = "name",
) -> pd.DataFrame:
    """Map a list of metabolite identifiers to HMDB accessions.

    Performs one HMDB search per input identifier and returns the
    best match (first result by name similarity).

    Parameters
    ----------
    ids : list
        List of identifiers to map.
    from_type : str
        Source identifier type.  One of ``"name"``, ``"kegg_id"``,
        ``"chebi_id"``, ``"inchikey"``, ``"cas"``, ``"pubchem_id"``,
        ``"formula"``.

    Returns
    -------
    pd.DataFrame
        Columns: input_id, hmdb_id, name, chemical_formula,
        inchikey, kegg_id, chebi_id.
    """
    records = []
    logger.info(f"HMDB ID mapping: {len(ids)} {from_type} identifiers")
    for input_id in ids:
        try:
            results = hmdb_search_metabolites(str(input_id), search_field=from_type, max_results=1)
            if not results.empty:
                row = results.iloc[0].to_dict()
                row["input_id"] = input_id
                records.append(row)
            else:
                records.append({"input_id": input_id, "hmdb_id": None})
        except Exception as exc:
            logger.debug(f"HMDB mapping failed for {input_id}: {exc}")
            records.append({"input_id": input_id, "hmdb_id": None})

    df = pd.DataFrame(records)
    mapped = df["hmdb_id"].notna().sum()
    logger.info(f"HMDB ID mapping: {mapped}/{len(ids)} successfully mapped")
    return df


# ---------------------------------------------------------------------------
# Integration helper — enrich Metabolite objects
# ---------------------------------------------------------------------------

def hmdb_enrich_metabolites(
    metabolites: list,
    fields: Optional[List[str]] = None,
) -> None:
    """Enrich a list of Metabolite objects with HMDB data in-place.

    Attempts to look up each metabolite in HMDB using its existing
    identifiers (HMDB ID > ChEBI ID > KEGG ID > name, in that priority
    order) and populates any missing fields.

    Parameters
    ----------
    metabolites : list of Metabolite
        The metabolites to enrich.  Each should be a
        :class:`~transnet.biology.elements.Metabolite` instance with
        at least one of: ``hmdb_id``, ``chebi_id``, ``kegg_compound_id``,
        or ``kegg_name``.
    fields : list of str, optional
        Subset of ``['inchi', 'inchikey', 'smiles', 'formula',
        'molecular_weight', 'chebi_id', 'pubchem_id', 'pathways',
        'diseases']`` to populate.  All fields are enriched by default.

    Notes
    -----
    The HMDB accession is stored on the metabolite as ``metabolite.hmdb_id``
    if the attribute exists (it is added dynamically if not present).
    """
    if fields is None:
        fields = [
            "inchi", "inchikey", "smiles", "formula",
            "molecular_weight", "chebi_id", "pubchem_id",
            "pathways", "diseases",
        ]

    logger.info(f"Enriching {len(metabolites)} metabolites with HMDB data")
    enriched = 0

    for met in metabolites:
        record = _lookup_metabolite_record(met)
        if not record:
            continue

        enriched += 1

        # Store HMDB accession
        if not getattr(met, "hmdb_id", None) and record.get("hmdb_id"):
            met.hmdb_id = record["hmdb_id"]

        # Fill fields
        if "inchi" in fields and not getattr(met, "inchi", None):
            met.inchi = record.get("inchi")
        if "inchikey" in fields and not getattr(met, "inchikey", None):
            met.inchikey = record.get("inchikey")
        if "smiles" in fields and not getattr(met, "smile", None):
            met.smile = record.get("smiles")
        if "formula" in fields:
            if not hasattr(met, "formula") or not met.formula:
                met.formula = record.get("chemical_formula")
        if "molecular_weight" in fields:
            if not hasattr(met, "molecular_weight") or not met.molecular_weight:
                met.molecular_weight = record.get("average_molecular_weight")
        if "chebi_id" in fields and not getattr(met, "chebi_id", None):
            met.chebi_id = record.get("chebi_id")
        if "pubchem_id" in fields and not getattr(met, "pubchem_id", None):
            met.pubchem_id = record.get("pubchem_compound_id")

        # Pathways — stored as a list of dicts in met.hmdb_pathways
        if "pathways" in fields:
            if not hasattr(met, "hmdb_pathways"):
                met.hmdb_pathways = []
            met.hmdb_pathways = record.get("pathways", [])

        # Diseases — stored in met.hmdb_diseases
        if "diseases" in fields:
            if not hasattr(met, "hmdb_diseases"):
                met.hmdb_diseases = []
            met.hmdb_diseases = record.get("diseases", [])

    logger.info(f"HMDB enrichment complete: {enriched}/{len(metabolites)} metabolites matched")


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _normalise_hmdb_id(hmdb_id: str) -> str:
    """Convert legacy 5-digit HMDB IDs to current 7-digit format."""
    hmdb_id = hmdb_id.strip().upper()
    if hmdb_id.startswith("HMDB") and len(hmdb_id) == 9:
        # Legacy HMDB00001 → HMDB0000001
        return "HMDB" + hmdb_id[4:].zfill(7)
    return hmdb_id


def _lookup_metabolite_record(met) -> dict:
    """Try several identifiers to find the HMDB record for *met*."""
    # Priority 1: existing HMDB ID
    hmdb_id = getattr(met, "hmdb_id", None)
    if hmdb_id:
        record = hmdb_get_metabolite(hmdb_id)
        if record:
            return record

    # Priority 2: ChEBI ID
    chebi_id = getattr(met, "chebi_id", None)
    if chebi_id:
        try:
            results = hmdb_search_metabolites(str(chebi_id), search_field="chebi_id", max_results=1)
            if not results.empty and results.iloc[0].get("hmdb_id"):
                record = hmdb_get_metabolite(results.iloc[0]["hmdb_id"])
                if record:
                    return record
        except Exception:
            pass

    # Priority 3: KEGG compound ID
    kegg_id = getattr(met, "kegg_compound_id", None)
    if kegg_id:
        try:
            results = hmdb_search_metabolites(str(kegg_id), search_field="kegg_id", max_results=1)
            if not results.empty and results.iloc[0].get("hmdb_id"):
                record = hmdb_get_metabolite(results.iloc[0]["hmdb_id"])
                if record:
                    return record
        except Exception:
            pass

    # Priority 4: name search
    name = getattr(met, "kegg_name", None)
    if name:
        try:
            results = hmdb_search_metabolites(str(name), search_field="name", max_results=1)
            if not results.empty and results.iloc[0].get("hmdb_id"):
                record = hmdb_get_metabolite(results.iloc[0]["hmdb_id"])
                if record:
                    return record
        except Exception:
            pass

    return {}
