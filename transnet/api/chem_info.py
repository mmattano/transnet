"""Converts a list of chemical identifiers to a pandas dataframe.
"""

__all__ = [
    "chemical_info_converter",
]

import pandas as pd
import requests


def chemical_info_converter(input_ids):
    """
    Converts a list of chemical identifiers to a pandas dataframe with the following columns:
    chebi_ids: CHEBI identifier (the first, when several are known)
    chebi_ids_all: every CHEBI identifier known for the compound, as a list
    smiles: SMILES strings
    inchi: InChI strings
    inchikeys: InChIKeys
    pubchem_ids: PubChem identifiers
    Input can be any of these identifiers, but must be a list of strings.
    
    Parameters
    ----------
    input_ids: list
        List of chemical identifiers.
    
    Returns
    -------
    pandas.DataFrame
        Dataframe with the following columns:
        chebi_ids: CHEBI identifier (the first, when several are known)
        chebi_ids_all: every CHEBI identifier for the compound, as a list.
            Prefer this when translating onward to another database: the
            neutral and zwitterionic forms are separate CHEBI entries and
            only some are cross-referenced.
        smiles: SMILES strings
        inchi: InChI strings
        inchikeys: InChIKeys
        pubchem_ids: PubChem identifiers
    """

    MAX_SIMULTANEOUS_REQUESTS = 100

    chebi_ids = []
    chebi_ids_all = []
    smiles = []
    inchi = []
    inchikeys = []
    pubchem_ids = []

    for i in range(0, len(input_ids), MAX_SIMULTANEOUS_REQUESTS):
        start = i
        end = i + MAX_SIMULTANEOUS_REQUESTS
        if end > len(input_ids):
            end = len(input_ids)
        params = {
            "ids": ",".join(str(x) for x in input_ids[start:end]),
            "fields": "pubchem.cid, chebi.id, chebi.inchi, chebi.inchikey, chebi.smiles",
        }
        res = requests.post("http://mychem.info/v1/chem", params)
        con = res.json()
        for j in con:
            try:
                # mychem often returns several ChEBI entries for one compound
                # -- typically the neutral species and its zwitterion. Only
                # some of them are cross-referenced by KEGG, so keeping just
                # the first silently loses the mapping for the rest: on a
                # 171-metabolite mouse panel that was 106 compounds resolved
                # instead of 144. Every id is kept in `chebi_ids_all`, with
                # the first still in `chebi_ids` for callers that want one.
                entries = j["chebi"] if isinstance(j["chebi"], list) else [j["chebi"]]
                ids = [e["id"] for e in entries if isinstance(e, dict) and "id" in e]
                if not ids:
                    raise KeyError("chebi")
                chebi_ids.append(ids[0])
                chebi_ids_all.append(ids)
                first = entries[0] if isinstance(entries[0], dict) else {}
                smiles.append(first.get("smiles"))
                inchi.append(first.get("inchi"))
                inchikeys.append(first.get("inchikey"))
            except KeyError:
                chebi_ids.append(None)
                chebi_ids_all.append([])
                smiles.append(None)
                inchi.append(None)
                inchikeys.append(None)
            try:
                if isinstance(j["pubchem"], list):
                    pubchem_ids.append(str(j["pubchem"][0]["cid"]))
                else:
                    pubchem_ids.append(str(j["pubchem"]["cid"]))
            except KeyError:
                pubchem_ids.append(None)

    chem_info_df = pd.DataFrame(
        {
            "chebi_ids": chebi_ids,
            "chebi_ids_all": chebi_ids_all,
            "smiles": smiles,
            "inchi": inchi,
            "inchikeys": inchikeys,
            "pubchem_ids": pubchem_ids,
        }
    )

    return chem_info_df
