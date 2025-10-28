"""
Sitelist module for UVAFME vegetation model.

This module handles initialization of site lists and manages
the attachment of species data to sites.

Translated from Sitelist.f90
"""

from typing import List, Optional, Tuple
import csv
import os
from .constants import INVALID, RNVALID
from .site import SiteData
from .species import SpeciesData
from .genus_groups import Groups, initialize_genus_groups
from .input_module import InputFileManager
from .io_utils import count_records, split_line, fatal_error, warning


# Global variable for number of sites
numsites = 0


def initialize_sitelist(sites: List[SiteData], species_present: Groups,
                       species_file: str, range_file: Optional[str] = None) -> None:
    """
    Initialize the complete site list with all associated data.

    This function follows the Fortran implementation more closely:
    1. Assumes sites are already loaded with basic site information
    2. Applies global parameter adjustments to sites
    3. Reads and attaches species data to sites
    4. Initializes genus groups

    Args:
        sites: List of SiteData objects (assumed to be pre-populated)
        species_present: Groups object to be populated with species info
        species_file: Path to species data file
        range_file: Optional path to species range file
    """
    global numsites

    # Set global numsites variable
    numsites = len(sites)
    print(f"Initializing sitelist with {numsites} sites")

    # Apply global parameter adjustments to sites
    for site in sites:
        if site.site_id != INVALID:
            # Initialize altitude to "no adjustment" value
            site.altitude = RNVALID

            # Apply site-specific adjustments based on global parameters
            apply_site_adjustments(site)

    # Initialize species lists for each site
    initialize_spec_list(sites, species_file, range_file)

    # Get the lists of genera and species present across all sites
    # This modifies species_present in place
    updated_groups = initialize_genus_groups(sites)

    # Copy the data from updated_groups to species_present
    if updated_groups:
        species_present.numspecies = updated_groups.numspecies
        species_present.spec_names = updated_groups.spec_names
        species_present.genusgroups = updated_groups.genusgroups

    print(f"Sitelist initialization complete. Species present: {species_present.numspecies}")


def apply_site_adjustments(site: SiteData) -> None:
    """
    Apply standard adjustments to site parameters.

    This mirrors the Fortran implementation that applies global parameter
    overrides and standard unit conversions to site data.

    Args:
        site: Site data to adjust
    """
    # Create parameters instance to access global overrides
    from .parameters import Parameters
    params = Parameters()

    # Apply global parameter overrides (if set to non-default values)
    if params.new_slope != RNVALID:
        site.slope = params.new_slope

    if params.fire_level != RNVALID:
        site.fire_prob = params.fire_level

    if params.wind_level != RNVALID:
        site.wind_prob = params.wind_level

    if params.SA_field_cap != RNVALID:
        site.soil.A_field_cap = params.SA_field_cap

    if params.A0_level_C != RNVALID:
        site.soil.A0_c0 = params.A0_level_C

    if params.A0_level_N != RNVALID:
        site.soil.A0_n0 = params.A0_level_N

    # Apply standard adjustments (these always happen)
    # Convert probabilities from per-1000-years to annual probabilities
    site.fire_prob = site.fire_prob / 1000.0
    site.wind_prob = site.wind_prob / 1000.0

    # Adjust soil parameters by root depth
    site.soil.A_w0 = site.soil.A_w0 * params.rootdepth
    site.soil.A_field_cap = site.soil.A_field_cap * params.rootdepth
    site.soil.A_perm_wp = site.soil.A_perm_wp * params.rootdepth

    # Initialize leaf area index
    site.leaf_area_ind = 1.0


def initialize_spec_list(sites: List[SiteData], species_file: str,
                        range_file: Optional[str] = None) -> None:
    """
    Initialize species lists for all sites.

    This function reads species data and optionally range data to determine
    which species are present at each site.

    Args:
        sites: List of site data structures
        species_file: Path to species data file
        range_file: Optional path to species range file
    """
    # Create input manager and read species data
    input_manager = InputFileManager()

    # Configure file paths
    input_manager.filenames['species'] = species_file
    if range_file:
        input_manager.filenames['range'] = range_file

    # Read all species data
    species_data = input_manager.read_species_data()

    # Read range data if available
    use_rangelist = False
    range_site_ids = []
    range_species_ids = []

    if range_file and os.path.exists(range_file):
        use_rangelist, range_site_ids, range_species_ids = read_rangelist(range_file)

    # Attach species to sites based on range data or use all species
    if use_rangelist:
        # Add range information to Site structures
        for site in sites:
            found_site = False
            for i, range_site_id in enumerate(range_site_ids):
                if site.site_id == range_site_id:
                    # Construct the list of species present for this site
                    found_site = True
                    attach_species_to_site(site, species_data, range_species_ids[i])
                    break

            # If site wasn't in the rangelist, add all species
            if not found_site:
                attach_species_to_site(site, species_data)
    else:
        # All species present in all sites
        for site in sites:
            attach_species_to_site(site, species_data)


def read_rangelist(range_file: str) -> Tuple[bool, List[int], List[List[str]]]:
    """
    Read species range list file.

    The range list file format:
    - First line: header with "site_id,lat,long,species1,species2,..."
    - Subsequent lines: site_id,lat,long,1,0,1,... (1=present, 0=absent)

    Args:
        range_file: Path to range list file

    Returns:
        Tuple of (use_rangelist, site_ids, species_lists)
        - use_rangelist: Whether range data was successfully read
        - site_ids: List of site IDs
        - species_lists: List of species ID lists for each site
    """
    if not os.path.exists(range_file):
        warning(f"Range file {range_file} not found, using all species for all sites")
        return False, [], []

    try:
        with open(range_file, 'r') as f:
            reader = csv.reader(f)

            # Read header line to get species names
            header = next(reader)
            if len(header) < 4:  # site_id, lat, long, at least one species
                warning("Invalid range file format")
                return False, [], []

            # Species names start from column 3 (after site_id, lat, long)
            species_names = header[3:]
            site_ids = []
            species_lists = []

            # Read data lines
            for row in reader:
                if len(row) < len(header):
                    continue  # Skip incomplete rows

                site_id = int(row[0])
                site_ids.append(site_id)

                # Extract species presence (1/0) and map to species names
                species_present = []
                for i, presence in enumerate(row[3:]):
                    if i < len(species_names) and int(presence) == 1:
                        species_present.append(species_names[i])

                species_lists.append(species_present)

        print(f"Range list read: {len(site_ids)} sites with species ranges")
        return True, site_ids, species_lists

    except Exception as e:
        warning(f"Error reading range file {range_file}: {e}")
        return False, [], []


def attach_species_to_site(site: SiteData, species_data: List[SpeciesData],
                          species_ids: Optional[List[str]] = None) -> None:
    """
    Attach species data to a specific site.

    Args:
        site: Site to attach species to
        species_data: List of all available species
        species_ids: Optional list of species IDs to include (if None, include all)
    """
    if species_ids is None:
        # Add all species to this site
        site.species = species_data.copy()
    else:
        # Add only specified species
        site.species = []
        for species in species_data:
            if species.unique_id in species_ids:
                site.species.append(species)

    print(f"Site {site.site_id}: attached {len(site.species)} species")


def validate_site_data(sites: List[SiteData]) -> bool:
    """
    Validate that site data is complete and consistent.

    Args:
        sites: List of sites to validate

    Returns:
        True if all sites are valid, False otherwise
    """
    for site in sites:
        # Check that site has valid ID
        if site.site_id == INVALID:
            continue

        # Check that site has climate data
        if site.site_wmo == RNVALID:
            print(f"Warning: No climate data for site {site.site_id}")
            return False

        # Check that site has species
        if len(site.species) == 0:
            print(f"Warning: No species present in site {site.site_id}")
            return False

    return True


def get_site_by_id(sites: List[SiteData], site_id: int) -> Optional[SiteData]:
    """
    Find a site by its ID.

    Args:
        sites: List of sites to search
        site_id: ID of site to find

    Returns:
        SiteData object if found, None otherwise
    """
    for site in sites:
        if site.site_id == site_id:
            return site
    return None


def count_total_species(sites: List[SiteData]) -> int:
    """
    Count total number of unique species across all sites.

    Args:
        sites: List of sites

    Returns:
        Number of unique species
    """
    all_species_ids = set()

    for site in sites:
        for species in site.species:
            all_species_ids.add(species.unique_id)

    return len(all_species_ids)


def print_site_summary(sites: List[SiteData]) -> None:
    """
    Print a summary of site information.

    Args:
        sites: List of sites
    """
    valid_sites = [s for s in sites if s.site_id != INVALID]

    print(f"Site Summary:")
    print(f"  Total sites: {len(sites)}")
    print(f"  Valid sites: {len(valid_sites)}")
    print(f"  Total unique species: {count_total_species(sites)}")

    for site in valid_sites:
        print(f"  Site {site.site_id} ({site.site_name}): {len(site.species)} species")