#!/usr/bin/env python3
"""
Utilities for loading AFluID metadata reference JSON trees.

Purpose
-------
This module converts the new JSON reference trees into the legacy dictionary
objects expected by final_report_utils.py.

Inputs
------
data/un_geodump/geo_tree.json
    Geographic hierarchy:
        region
            sub-region
                intermediate-region
                    country

data/ncbi_taxdump/host_taxonomy_tree.json
    Host taxonomy hierarchy:
        class
            order
                family
                    genus
                        species

Outputs
-------
geo_dict
    Legacy-compatible geographic dictionary:

        {
            "Portugal": ("Southern Europe", "Europe"),
            "Côte d’Ivoire": ("Western Africa", "Africa"),
        }

taxa_dict
    Legacy-compatible taxonomy dictionary:

        {
            "Homo sapiens": [
                "Homo sapiens",  # Species
                "Homo",          # Genus
                "",              # Subfamily
                "Hominidae",     # Family
                "",              # Division
                "Primates",      # Order
                "Mammalia",      # Class
            ]
        }

Notes
-----
This module does not apply synonym dictionaries yet.

Synonyms/exceptions should remain in structures.py and be applied later by
the parsing / cleaning layer.
"""

from pathlib import Path
import json
from typing import Any


# ============================================================
# Legacy level order expected by final_report_utils.py
# ============================================================

TAXA_LEGACY_LEVELS = [
    "Species",
    "Genus",
    "Subfamily",
    "Family",
    "Division",
    "Order",
    "Class",
]

TAXA_RANK_TO_INDEX = {
    "species": 0,
    "genus": 1,
    "subfamily": 2,
    "family": 3,
    "division": 4,
    "order": 5,
    "class": 6,
}

LEGACY_REGION_NAME_MAP = {
    "South-eastern Asia": "South-Eastern Asia",
}

LEGACY_AMERICAS_CONTINENT_BY_REGION = {
    "Northern America": "North America",
    "Central America": "North America",
    "Caribbean": "North America",
    "South America": "South America",
}


# ============================================================
# Basic JSON loading
# ============================================================

def load_json_payload(json_fp: str | Path) -> dict[str, Any]:
    """
    Load a JSON file and return its payload.

    Parameters
    ----------
    json_fp : str or pathlib.Path
        Path to the JSON file.

    Returns
    -------
    dict
        Parsed JSON payload.

    Raises
    ------
    FileNotFoundError
        If json_fp does not exist.

    ValueError
        If the JSON payload does not contain a dictionary.
    """
    json_fp = Path(json_fp)

    if not json_fp.exists():
        raise FileNotFoundError(f"Missing JSON file: {json_fp}")

    with json_fp.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)

    if not isinstance(payload, dict):
        raise ValueError(f"Invalid JSON payload in {json_fp}: expected object")

    return payload


def get_tree_from_payload(payload: dict[str, Any], source_name: str = "JSON") -> dict[str, Any]:
    """
    Extract the 'tree' object from a JSON payload.

    Parameters
    ----------
    payload : dict
        Parsed JSON payload.

    source_name : str
        Human-readable source name used in error messages.

    Returns
    -------
    dict
        Tree object.
    """
    tree = payload.get("tree")

    if not isinstance(tree, dict):
        raise ValueError(f"{source_name} does not contain a valid 'tree' object")

    return tree


# ============================================================
# Generic tree traversal
# ============================================================

def walk_named_tree(tree: dict[str, Any]):
    """
    Traverse a nominal tree and yield each node with its lineage.

    Expected input structure
    ------------------------
    {
        "Aves": {
            "rank": "class",
            "children": {
                ...
            }
        }
    }

    Yields
    ------
    tuple
        (node_name, node_data, lineage)

        where lineage is a list of dictionaries:

        [
            {"name": "Aves", "rank": "class"},
            {"name": "Anseriformes", "rank": "order"},
            ...
        ]
    """

    def _walk(node_name: str, node_data: dict[str, Any], lineage: list[dict[str, str]]):
        rank = str(node_data.get("rank", "")).strip()

        current = {
            "name": node_name,
            "rank": rank,
        }

        new_lineage = lineage + [current]

        yield node_name, node_data, new_lineage

        children = node_data.get("children", {})

        if not isinstance(children, dict):
            return

        for child_name, child_data in children.items():
            if isinstance(child_data, dict):
                yield from _walk(child_name, child_data, new_lineage)

    for root_name, root_data in tree.items():
        if isinstance(root_data, dict):
            yield from _walk(root_name, root_data, [])


# ============================================================
# Geographic JSON -> legacy geo_dict
# ============================================================

def lineage_to_geo_tuple(lineage: list[dict[str, str]]) -> tuple[str, str]:
    """
    Convert a geographic lineage into the legacy (region, continent) tuple.

    The UNSD-derived tree may look like:

        Africa > Sub-Saharan Africa > Western Africa > Côte d’Ivoire

    In the old geo_dict structure, this should become:

        "Côte d’Ivoire": ("Western Africa", "Africa")

    Rule
    ----
    - continent = top-level root node
    - region = nearest geographic parent before country
      usually intermediate-region if present, otherwise sub-region
    """
    if not lineage:
        return "Unknown", "Unknown"

    continent = lineage[0]["name"] or "Unknown"

    # The country is the final node, so its parent is the previous node.
    if len(lineage) >= 2:
        region = lineage[-2]["name"] or "Unknown"
    else:
        region = "Unknown"

    region = LEGACY_REGION_NAME_MAP.get(region, region)

    if continent == "Americas":
        continent = LEGACY_AMERICAS_CONTINENT_BY_REGION.get(region, continent)

    return region, continent


def build_geo_dict_from_tree(geo_tree: dict[str, Any]) -> dict[str, tuple[str, str]]:
    """
    Convert geo_tree.json tree object into legacy geo_dict.

    Parameters
    ----------
    geo_tree : dict
        Tree object from geo_tree.json.

    Returns
    -------
    dict
        Legacy-compatible dictionary:

            {
                country: (region, continent)
            }
    """
    geo_dict = {}

    for node_name, node_data, lineage in walk_named_tree(geo_tree):
        rank = str(node_data.get("rank", "")).strip().casefold()

        if rank != "country":
            continue

        region, continent = lineage_to_geo_tuple(lineage)
        geo_dict[node_name] = (region, continent)

    return geo_dict


def load_geo_dict_from_json(geo_json_fp: str | Path) -> dict[str, tuple[str, str]]:
    """
    Load geo_tree.json and convert it into legacy geo_dict.
    """
    payload = load_json_payload(geo_json_fp)
    tree = get_tree_from_payload(payload, source_name=str(geo_json_fp))
    return build_geo_dict_from_tree(tree)


# ============================================================
# Host taxonomy JSON -> legacy taxa_dict
# ============================================================

def lineage_to_taxa_list(lineage: list[dict[str, str]]) -> list[str]:
    """
    Convert a taxonomy lineage into the legacy taxa_dict list format.

    Legacy order expected by final_report_utils.py:

        0 Species
        1 Genus
        2 Subfamily
        3 Family
        4 Division
        5 Order
        6 Class

    Missing levels are stored as "X" because downstream parsing expects
    a non-empty placeholder.
    """
    taxa_list = ["X"] * len(TAXA_LEGACY_LEVELS)

    for item in lineage:
        rank = str(item.get("rank", "")).strip().casefold()
        name = str(item.get("name", "")).strip()

        if not rank or not name:
            continue

        idx = TAXA_RANK_TO_INDEX.get(rank)

        if idx is None:
            continue

        taxa_list[idx] = name

    # Preserve the legacy 7-slot contract even when the source tree omits
    # intermediate ranks such as subfamily/division.
    if taxa_list[TAXA_RANK_TO_INDEX["subfamily"]] == "X" and taxa_list[TAXA_RANK_TO_INDEX["family"]] != "X":
        taxa_list[TAXA_RANK_TO_INDEX["subfamily"]] = taxa_list[TAXA_RANK_TO_INDEX["family"]]

    if taxa_list[TAXA_RANK_TO_INDEX["division"]] == "X":
        order_name = taxa_list[TAXA_RANK_TO_INDEX["order"]]
        class_name = taxa_list[TAXA_RANK_TO_INDEX["class"]]

        if order_name != "X":
            taxa_list[TAXA_RANK_TO_INDEX["division"]] = order_name
        elif class_name != "X":
            taxa_list[TAXA_RANK_TO_INDEX["division"]] = class_name

    return taxa_list


def build_taxa_dict_from_tree(tax_tree: dict[str, Any]) -> dict[str, list[str]]:
    """
    Convert host_taxonomy_tree.json tree object into legacy taxa_dict.

    Parameters
    ----------
    tax_tree : dict
        Tree object from host_taxonomy_tree.json.

    Returns
    -------
    dict
        Legacy-compatible taxonomy dictionary.

    Notes
    -----
    Every named node is registered, not only species.

    This is important because valid host aliases may resolve to higher ranks,
    for example:

        "avian" -> "Aves"
        "duck"  -> "Anatidae"

    Therefore taxa_dict will contain entries such as:

        "Aves": ["", "", "", "", "", "", "Aves"]

        "Anatidae": ["", "", "", "Anatidae", "", "Anseriformes", "Aves"]

        "Anas platyrhynchos": [
            "Anas platyrhynchos",
            "Anas",
            "",
            "Anatidae",
            "",
            "Anseriformes",
            "Aves",
        ]
    """
    taxa_dict = {}

    for node_name, node_data, lineage in walk_named_tree(tax_tree):
        rank = str(node_data.get("rank", "")).strip().casefold()

        if rank not in TAXA_RANK_TO_INDEX:
            continue

        taxa_dict[node_name] = lineage_to_taxa_list(lineage)

    return taxa_dict


def load_taxa_dict_from_json(tax_json_fp: str | Path) -> dict[str, list[str]]:
    """
    Load host_taxonomy_tree.json and convert it into legacy taxa_dict.
    """
    payload = load_json_payload(tax_json_fp)
    tree = get_tree_from_payload(payload, source_name=str(tax_json_fp))
    return build_taxa_dict_from_tree(tree)


# ============================================================
# Combined loader for AFluID
# ============================================================

def load_reference_dicts(
    data_dir: str | Path = "data",
    geo_json_relpath: str | Path = "un_geodump/geo_tree.json",
    tax_json_relpath: str | Path = "ncbi_taxdump/host_taxonomy_tree.json",
) -> tuple[dict[str, tuple[str, str]], dict[str, list[str]]]:
    """
    Load both JSON reference trees and return legacy dictionaries.

    Parameters
    ----------
    data_dir : str or pathlib.Path, default "data"
        Base data directory.

    geo_json_relpath : str or pathlib.Path
        Relative path from data_dir to geo_tree.json.

    tax_json_relpath : str or pathlib.Path
        Relative path from data_dir to host_taxonomy_tree.json.

    Returns
    -------
    tuple
        (geo_dict, taxa_dict)
    """
    data_dir = Path(data_dir)

    geo_json_fp = data_dir / geo_json_relpath
    tax_json_fp = data_dir / tax_json_relpath

    geo_dict = load_geo_dict_from_json(geo_json_fp)
    taxa_dict = load_taxa_dict_from_json(tax_json_fp)

    return geo_dict, taxa_dict


# ============================================================
# Basic command-line check
# ============================================================

def main():
    """
    Basic command-line validation.

    Run from the AFluID repository root:

        python metadata_tree_utils.py
    """
    geo_dict, taxa_dict = load_reference_dicts("data")

    print("Loaded reference dictionaries")
    print(f"geo_dict entries:  {len(geo_dict)}")
    print(f"taxa_dict entries: {len(taxa_dict)}")

    print("\nGeo checks:")
    for country in ["Portugal", "United States of America", "Republic of Korea", "Taiwan", "Kosovo"]:
        print(f"  {country}: {geo_dict.get(country, 'NOT FOUND')}")

    print("\nTaxa checks:")
    for taxon in ["Aves", "Anatidae", "Homo sapiens", "Sus scrofa", "Gallus gallus"]:
        print(f"  {taxon}: {taxa_dict.get(taxon, 'NOT FOUND')}")


if __name__ == "__main__":
    main()
