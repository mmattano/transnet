"""STRING API"""

from collections import defaultdict

__all__ = [
    "string_map_identifiers",
    "string_get_interactions",
    "reverse_string_mapping",
    "translate_string_dict",
]

import requests
import time

_STRING_API = "https://string-db.org/api"


def string_map_identifiers(protein_list: list, species: str = "9606") -> dict:
    """Map UniProt / gene-name identifiers to STRING ENSP identifiers.

    Parameters
    ----------
    protein_list : list
        List of protein identifiers (UniProt IDs or gene symbols).
    species : str
        NCBI species taxon ID (e.g. ``"10090"`` for mouse).

    Returns
    -------
    dict
        ``{input_id: string_ensp_id}``
    """
    request_url = "/".join([_STRING_API, "tsv-no-header", "get_string_ids"])
    params = {
        "identifiers": "\r".join(protein_list),
        "species": species,
        "limit": 1,
        "echo_query": 1,
    }
    results = requests.post(request_url, data=params, timeout=120)
    time.sleep(1)

    string_map = {}
    for line in results.text.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        string_map[parts[0]] = parts[2]

    return string_map


def string_get_interactions(
    protein_list: list, species: str = "9606", cutoff_score: int = 700
) -> dict:
    """Fetch interaction partners for a list of STRING ENSP identifiers.

    Parameters
    ----------
    protein_list : list
        STRING ENSP identifiers (output of :func:`string_map_identifiers`).
    species : str
        NCBI species taxon ID.
    cutoff_score : int
        Minimum combined interaction score (0–1000).

    Returns
    -------
    dict
        ``{ensp_id: [partner_ensp_id, ...]}``
    """
    request_url = "/".join(
        [_STRING_API, "tsv-no-header", "interaction_partners"]
    )
    params = {
        "identifiers": "\r".join(protein_list),
        "species": species,
        "required_score": cutoff_score,
    }
    response = requests.post(request_url, data=params, timeout=120)
    time.sleep(1)

    if not response.ok:
        raise RuntimeError(
            f"STRING API error {response.status_code}: "
            f"{response.text[:200]}"
        )

    interactions_dict = defaultdict(list)
    for line in response.text.strip().split("\n"):
        parts = line.strip().split("\t")
        if len(parts) < 2:
            continue
        # Guard: STRING returns JSON error payloads on bad requests
        if parts[0].strip().lower() in ("error", "errormessage"):
            raise RuntimeError(
                f"STRING API returned error: {response.text[:200]}"
            )
        interactions_dict[parts[0]].append(parts[1])

    return interactions_dict


def reverse_string_mapping(mapping_dict: dict, string_identifier: str):
    key = list(mapping_dict.keys())[
        list(mapping_dict.values()).index(string_identifier)
    ]
    return key


def translate_string_dict(mapping_dict: dict, interactions_dict: dict) -> dict:
    """Translate STRING ENSP identifiers back to input identifiers.

    Parameters
    ----------
    mapping_dict : dict
        Output of :func:`string_map_identifiers`.
    interactions_dict : dict
        Output of :func:`string_get_interactions`.

    Returns
    -------
    dict
        ``{input_id: [partner_input_id, ...]}``
    """
    reverse = {v: k for k, v in mapping_dict.items()}
    translated_dict = {}
    for ensp_center, partners in interactions_dict.items():
        center_node = reverse.get(ensp_center)
        if center_node is None:
            continue
        translated_dict[center_node] = [
            reverse[p] for p in partners if p in reverse
        ]
    return translated_dict
