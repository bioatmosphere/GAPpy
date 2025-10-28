"""
Genus Groups module for UVAFME vegetation model.

This module handles grouping of species by genus and manages
species/genus lists and their unique identifiers.

Translated from GenusGroups.f90
"""

from typing import List, Tuple, Optional
from .constants import MAX_NLEN
from .site import SiteData
from .species import SpeciesData
from .utilities import sort_strings


class Groups:
    """
    Data structure for storing genus groups and species information.

    Corresponds to the Groups type in the Fortran code.
    """

    def __init__(self):
        self.genusgroups: List[str] = []  # List of unique genus names
        self.spec_names: List[Tuple[str, str]] = []  # (genus_name, unique_id) pairs
        self.numgenera: int = 0
        self.numspecies: int = 0


def initialize_genus_groups(sites: List[SiteData]) -> Groups:
    """
    Initialize genus groups from site data.

    This function extracts all unique genera and species from all sites
    and creates organized lists for output purposes.

    Args:
        sites: List of site data structures

    Returns:
        Groups object containing organized genus/species information
    """
    group = Groups()

    # Collect all species data from all sites
    all_species_data = []

    for site in sites:
        for species in site.species:
            all_species_data.append(species)

    if not all_species_data:
        return group

    # Extract unique species IDs and genus names
    species_ids = [sp.unique_id for sp in all_species_data]
    genus_names = [sp.genus_name for sp in all_species_data]

    unique_species = get_unique_items(species_ids)
    unique_genera = get_unique_items(genus_names)

    group.numspecies = len(unique_species)
    group.numgenera = len(unique_genera)
    group.genusgroups = unique_genera

    # Create species names list with (genus_name, unique_id) pairs
    group.spec_names = []

    for unique_sp in unique_species:
        # Find the genus name for this species
        for sp_data in all_species_data:
            if sp_data.unique_id == unique_sp:
                group.spec_names.append((sp_data.genus_name, sp_data.unique_id))
                break

    return group


def get_unique_items(array: List[str]) -> List[str]:
    """
    Get unique items from a string array.

    This function sorts the array and removes duplicates,
    similar to the Fortran get_unique_items subroutine.

    Args:
        array: List of strings

    Returns:
        List of unique strings in sorted order
    """
    if not array:
        return []

    if len(array) == 1:
        return [array[0]]

    # Sort the array
    sorted_array = sort_strings(array)

    # Extract unique items
    unique_items = [sorted_array[0]]

    for i in range(1, len(sorted_array)):
        if sorted_array[i-1] != sorted_array[i]:
            unique_items.append(sorted_array[i].strip())

    return unique_items


def find_genus_index(groups: Groups, genus_name: str) -> Optional[int]:
    """
    Find the index of a genus in the genus groups list.

    Args:
        groups: Groups object
        genus_name: Name of genus to find

    Returns:
        Index of genus (0-based) or None if not found
    """
    try:
        return groups.genusgroups.index(genus_name)
    except ValueError:
        return None


def find_species_index(groups: Groups, unique_id: str) -> Optional[int]:
    """
    Find the index of a species in the species list.

    Args:
        groups: Groups object
        unique_id: Unique identifier of species to find

    Returns:
        Index of species (0-based) or None if not found
    """
    for i, (genus, species_id) in enumerate(groups.spec_names):
        if species_id == unique_id:
            return i
    return None


def get_species_by_genus(groups: Groups, genus_name: str) -> List[str]:
    """
    Get all species unique IDs that belong to a specific genus.

    Args:
        groups: Groups object
        genus_name: Name of genus

    Returns:
        List of unique IDs for species in the genus
    """
    species_in_genus = []

    for genus, species_id in groups.spec_names:
        if genus == genus_name:
            species_in_genus.append(species_id)

    return species_in_genus


def print_genus_summary(groups: Groups):
    """
    Print a summary of genus groups information.

    Args:
        groups: Groups object
    """
    print(f"Genus Groups Summary:")
    print(f"  Number of genera: {groups.numgenera}")
    print(f"  Number of species: {groups.numspecies}")
    print(f"  Genera: {', '.join(groups.genusgroups)}")

    print(f"  Species by genus:")
    for genus in groups.genusgroups:
        species = get_species_by_genus(groups, genus)
        print(f"    {genus}: {', '.join(species)}")


def validate_groups(groups: Groups) -> bool:
    """
    Validate the consistency of the Groups object.

    Args:
        groups: Groups object to validate

    Returns:
        True if valid, False otherwise
    """
    # Check that counts match actual list lengths
    if len(groups.genusgroups) != groups.numgenera:
        return False

    if len(groups.spec_names) != groups.numspecies:
        return False

    # Check that all species have valid genus names
    for genus, species_id in groups.spec_names:
        if genus not in groups.genusgroups:
            return False

    return True