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
    chebi_ids: CHEBI identifiers
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
        chebi_ids: CHEBI identifiers
        smiles: SMILES strings
        inchi: InChI strings
        inchikeys: InChIKeys
        pubchem_ids: PubChem identifiers
    """

    MAX_SIMULTANEOUS_REQUESTS = 1000

    chebi_ids = []
    smiles = []
    inchi = []
    inchikeys = []
    pubchem_ids = []

    for i in range(0, len(input_ids), MAX_SIMULTANEOUS_REQUESTS):
        start = i
        end = i + 1000
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
                if isinstance(j["chebi"], list):
                    chebi_ids.append(j["chebi"][0]["id"])
                    try:
                        smiles.append(j["chebi"][0]["smiles"])
                    except KeyError:
                        smiles.append(None)
                    try:
                        inchi.append(j["chebi"][0]["inchi"])
                    except KeyError:
                        inchi.append(None)
                    try:
                        inchikeys.append(j["chebi"][0]["inchikey"])
                    except KeyError:
                        inchikeys.append(None)
                else:
                    chebi_ids.append(j["chebi"]["id"])
                    try:
                        smiles.append(j["chebi"]["smiles"])
                    except KeyError:
                        smiles.append(None)
                    try:
                        inchi.append(j["chebi"]["inchi"])
                    except KeyError:
                        inchi.append(None)
                    try:
                        inchikeys.append(j["chebi"]["inchikey"])
                    except KeyError:
                        inchikeys.append(None)
            except KeyError:
                chebi_ids.append(None)
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
            "smiles": smiles,
            "inchi": inchi,
            "inchikeys": inchikeys,
            "pubchem_ids": pubchem_ids,
        }
    )

    return chem_info_df
