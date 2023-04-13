"""STRING API
"""
# Cheated in string_map_identifiers by just skipping the html part
# talk to string people about how to avoid that

__all__ = [
    "string_map_identifiers",
    "string_get_interactions",
    "reverse_string_mapping",
    "translate_string_dict",
]

import requests
import time


def string_map_identifiers(protein_list: list, species: str = "9606") -> dict:
    """
    Find the STRING identifiers for a given gene.

    Parameters
    ----------
    protein_list : list
        List of protein names.
    species : str
        Species identifier.
    
    Returns
    -------
    string_map : dict
        Dictionary with STRING identifiers.
    """

    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv-no-header"
    method = "get_string_ids"

    ## Set parameters

    params = {
        "identifiers": "\r".join(protein_list),  # your protein list
        "species": species,  # species NCBI identifier
        "limit": 1,  # only one (best) identifier per input protein
        "echo_query": 1,  # see your input identifiers in the output
        # "caller_identity" : "www.awesome_app.org" # your app name
    }

    ## Construct URL

    request_url = "/".join([string_api_url, output_format, method])

    ## Call STRING

    results = requests.post(request_url, data=params)

    ## Wait 1 second to not risk overloading the server

    time.sleep(1)

    ## Read and parse the results

    string_map = {}
    for line in results.text.strip().split("\n"):
        l = line.split("\t")
        try:
            input_identifier, string_identifier = l[0], l[2]
            # print("Input:", input_identifier, "STRING:", string_identifier, sep="\t")
            string_map[input_identifier] = string_identifier
        except IndexError:
            continue

    return string_map


def string_get_interactions(
    protein_list: list, species: str = "9606", cutoff_score: int = 700
) -> dict:
    """
    Find the STRING identifiers for a given gene.

    Parameters
    ----------
    protein_list : list
        List of protein names.
    species : str
        Species identifier.
    
    Returns
    -------
    interactions_dict : dict
        Dictionary with STRING identifiers of interaction partners.
    """

    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv-no-header"
    method = "interaction_partners"

    ## Construct the request

    request_url = "/".join([string_api_url, output_format, method])

    ## Set parameters

    params = {
        "identifiers": "%0d".join(protein_list),  # your protein
        "species": species,  # species NCBI identifier
        # "limit" : 5,
        # "caller_identity" : "www.awesome_app.org" # your app name
        "required_score": cutoff_score,
    }

    ## Call STRING

    response = requests.post(request_url, data=params)

    ## Wait 1 second to not risk overloading the server

    time.sleep(1)

    ## Read and parse the results

    from collections import defaultdict

    interactions_dict = defaultdict(list)

    for line in response.text.strip().split("\n"):

        l = line.strip().split("\t")
        query_ensp = l[0]
        # query_name = l[2]
        partner_ensp = l[1]
        # partner_name = l[3]
        # combined_score = l[5]

        interactions_dict[query_ensp].append(partner_ensp)

        ## print

        # print("\t".join([query_ensp, query_name, partner_ensp, partner_name, combined_score]))

    return interactions_dict


def reverse_string_mapping(mapping_dict: dict, string_identifier: str):
    key = list(mapping_dict.keys())[
        list(mapping_dict.values()).index(string_identifier)
    ]
    return key


def translate_string_dict(
    mapping_dict: dict, interactions_dict: dict,
):
    """
    Translate STRING identifiers to gene names.

    Parameters
    ----------
    mapping_dict : dict
        Dictionary with STRING identifier mapping.
    interactions_dict : dict
        Dictionary with STRING identifiers of interaction partners.

    Returns
    -------
    translated_dict : dict
        Dictionary with gene names.
    """

    translated_dict = {}
    for key, value in interactions_dict.items():
        try:
            center_node = reverse_string_mapping(
                mapping_dict=mapping_dict, string_identifier=key
            )

            interaction_nodes = []
            for node in value:
                try:
                    interaction_nodes.append(
                        reverse_string_mapping(
                            mapping_dict=mapping_dict, string_identifier=node
                        )
                    )
                except ValueError:
                    continue
            translated_dict[center_node] = interaction_nodes
        except ValueError:
            print("Did not work for", key)
            continue

    return translated_dict
