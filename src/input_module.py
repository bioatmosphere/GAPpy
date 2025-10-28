"""
Input module for UVAFME vegetation model.

This module handles reading and parsing of all input files including
parameters, site data, climate data, species data, and range data.

Translated from Input.f90
"""

import csv
import json
import os
from typing import List, Dict, Tuple, Optional, Any
import numpy as np

from .constants import *
from .parameters import Parameters
from .site import SiteData
from .species import SpeciesData
from .io_utils import count_records, split_line, quote_strip, fatal_error, warning


class InputFileManager:
    """
    Manages all input file operations for the UVAFME model.

    This class handles opening, reading, and parsing of various input files
    including parameters, site data, climate data, and species data.
    """

    def __init__(self):
        self.parameters = Parameters()
        self.file_handles = {}
        self.filenames = {}

    def initialize_input_files(self, filelist: str = ""):
        """
        Initialize input files and parameters.

        Args:
            filelist: Optional path to file containing list of input files
        """
        self.open_input_files(filelist)
        self.initialize_parameters()

    def open_input_files(self, filelist: str = ""):
        """
        Open all required input files.

        Args:
            filelist: Path to file listing input file locations
        """
        # In a full implementation, this would parse the filelist
        # For now, we'll use default file paths with fallback to ../input_data

        # Support multiple possible locations for input_data directory
        possible_paths = ['input_data', '../input_data', '../../input_data']
        input_base_path = None

        for path in possible_paths:
            if os.path.exists(path):
                input_base_path = path
                break

        if input_base_path is None:
            raise FileNotFoundError("Could not find input_data directory in any expected location")

        self.filenames = {
            'runtime': os.path.join(input_base_path, 'uvafme_config.json'),
            'sitelist': os.path.join(input_base_path, 'sites.csv'),
            'climate': os.path.join(input_base_path, 'climate.csv'),
            'species': os.path.join(input_base_path, 'species.csv'),
            'sites': os.path.join(input_base_path, 'sites.csv'),
            'range': os.path.join(input_base_path, 'range.csv'),  # Optional
            'altitude': os.path.join(input_base_path, 'altitude.csv'),  # Optional
            'climate_std': os.path.join(input_base_path, 'climate_std.csv')  # Optional
        }

    def initialize_parameters(self):
        """Read and initialize all model parameters."""

        # Set defaults first
        self.set_default_parameters()

        # Try to read from parameter file
        param_file = self.filenames.get('runtime')
        if param_file and os.path.exists(param_file):
            try:
                self.parameters.load_from_file(param_file)
                print(f"Loaded parameters from {param_file}")
            except Exception as e:
                print(f"Error reading parameter file: {e}")
                print("Using default values")
        else:
            print("Parameter file not found, using default values")

        # Validate and process climate change parameters
        self.process_climate_parameters()

    def set_default_parameters(self):
        """Set default parameter values."""
        # This mirrors the default values in the Fortran Input.f90
        self.parameters.debug = False
        self.parameters.spinup = False
        self.parameters.spinup_yrs = 500
        self.parameters.year_print_interval = 10
        self.parameters.begin_change_year = 300
        self.parameters.duration_of_change = 0
        self.parameters.with_clim_change = False
        self.parameters.plot_level_data = False
        self.parameters.tree_level_data = False
        self.parameters.adjust_for_elev = False
        self.parameters.same_climate = True
        self.parameters.fixed_seed = False
        self.parameters.use_gcm = False
        self.parameters.start_gcm = 0
        self.parameters.end_gcm = 100

        self.parameters.numyears = 1000
        self.parameters.numplots = 200
        self.parameters.maxtrees = 1000
        self.parameters.maxheight = 60

        self.parameters.incr_tmin_by = 0.0
        self.parameters.incr_tmax_by = 0.0
        self.parameters.incr_precip_by = 0.0
        self.parameters.decr_tmin_by = 0.0
        self.parameters.decr_tmax_by = 0.0
        self.parameters.decr_precip_by = 0.0

        self.parameters.plotsize = 500.0
        self.parameters.rootdepth = 0.8
        self.parameters.incr_or_decr = 'incr'

    def process_climate_parameters(self):
        """Process and validate climate change parameters."""
        if self.parameters.with_clim_change:
            if not self.parameters.use_gcm:
                self.parameters.linear_cc = True

                if self.parameters.duration_of_change == 0:
                    fatal_error("Inconsistent input for climate change: 0 years duration")

                # Check if any changes are specified
                if (self.parameters.incr_tmin_by == 0.0 and self.parameters.decr_tmax_by == 0.0 and
                    self.parameters.decr_tmin_by == 0.0 and self.parameters.decr_tmax_by == 0.0 and
                    self.parameters.incr_precip_by == 0.0 and self.parameters.decr_precip_by == 0.0):
                    fatal_error("Inconsistent climate change: no temp or precip changes")

                # Calculate per-year changes
                if self.parameters.incr_or_decr == 'incr':
                    self.parameters.tmin_change = self.parameters.incr_tmin_by / (self.parameters.duration_of_change + 1)
                    self.parameters.tmax_change = self.parameters.incr_tmax_by / (self.parameters.duration_of_change + 1)
                    self.parameters.precip_change = self.parameters.incr_precip_by / (self.parameters.duration_of_change + 1)
                elif self.parameters.incr_or_decr == 'decr':
                    # Handle negative values
                    if self.parameters.decr_tmin_by < 0.0:
                        print("Assuming decrease intended")
                        self.parameters.decr_tmin_by = -self.parameters.decr_tmin_by
                    if self.parameters.decr_tmax_by < 0.0:
                        print("Assuming decrease intended")
                        self.parameters.decr_tmax_by = -self.parameters.decr_tmax_by
                    if self.parameters.decr_precip_by < 0.0:
                        print("Assuming decrease intended")
                        self.parameters.decr_precip_by = -self.parameters.decr_precip_by

                    self.parameters.tmin_change = -self.parameters.decr_tmin_by / (self.parameters.duration_of_change + 1)
                    self.parameters.tmax_change = -self.parameters.decr_tmax_by / (self.parameters.duration_of_change + 1)
                    self.parameters.precip_change = -self.parameters.decr_precip_by / (self.parameters.duration_of_change + 1)
            else:
                # Using GCM data
                self.parameters.linear_cc = False
                self.parameters.tmin_change = 0.0
                self.parameters.tmax_change = 0.0
                self.parameters.precip_change = 0.0

                if self.parameters.duration_of_change == 0:
                    fatal_error("Inconsistent input for climate change: 0 years duration")
                else:
                    self.parameters.end_gcm = self.parameters.start_gcm + self.parameters.duration_of_change - 1

        # Set debug mode implications
        if self.parameters.debug:
            self.parameters.fixed_seed = True
            self.parameters.same_climate = True

    def read_sitelist(self) -> List[int]:
        """
        Read list of site IDs from sitelist file.

        Returns:
            List of site IDs
        """
        sitelist_file = self.filenames.get('sitelist')
        if not sitelist_file or not os.path.exists(sitelist_file):
            fatal_error('Unable to find site list file')

        site_ids = []

        try:
            with open(sitelist_file, 'r') as f:
                reader = csv.reader(f)
                for row in reader:
                    if row and row[0].strip().isdigit():
                        site_ids.append(int(row[0].strip()))
        except Exception as e:
            fatal_error(f"Error reading site list file: {e}")

        if not site_ids:
            fatal_error("No site data read")

        print(f'Site data initialized. Total read in: {len(site_ids)}')
        return site_ids

    def read_sites(self, site_ids: List[int]) -> List[SiteData]:
        """
        Read site data from sites file.

        Args:
            site_ids: List of site IDs to read

        Returns:
            List of SiteData objects
        """
        sites_file = self.filenames.get('sites')
        if not sites_file or not os.path.exists(sites_file):
            fatal_error('Unable to find site data file')

        sites = []

        try:
            with open(sites_file, 'r') as f:
                reader = csv.DictReader(f)

                for row in reader:
                    site_id = int(row['site'])  # 'site' maps to site_id

                    if site_id in site_ids:
                        site = SiteData()

                        # Parse site data from CSV row using UVAFME column names
                        site.initialize_site(
                            siteid=site_id,
                            sitename=row['name'],  # 'name' maps to site_name
                            siteregion=row['region'],  # 'region' maps to region
                            lat=float(row['latitude']),
                            long=float(row['longitude']),
                            wmo=float(row['wmo']),
                            elevation=float(row['elevation']),
                            slope=float(row['slope']),
                            Afc=float(row['soilA_field_cap']),  # 'soilA_field_cap' maps to A_field_cap
                            A_perm_wp=float(row['soilA_perm_wp']),  # 'soilA_perm_wp' maps to A_perm_wp
                            lai=float(row['lai']),
                            base_h=float(row['soil_base_h']),  # 'soil_base_h' maps to base_h
                            lai_w0=float(row['lai_w0']),
                            A0_w0=float(row['soilAO_w0']),  # 'soilAO_w0' maps to A0_w0
                            A_w0=float(row['soilA_w0']),  # 'soilA_w0' maps to A_w0
                            sbase_w0=float(row['sbase_w0']),
                            fire_prob=float(row['fire_prob']),
                            wind_prob=float(row['wind_prob']),
                            A0_c0=float(row['soilAO_c0']),  # 'soilAO_c0' maps to A0_c0
                            A0_n0=float(row['soilAO_n0']),  # 'soilAO_n0' maps to A0_n0
                            A_c0=float(row['soilA_c0']),  # 'soilA_c0' maps to A_c0
                            A_n0=float(row['soilA_n0']),  # 'soilA_n0' maps to A_n0
                            sbase_c0=float(row['sbase_c0']),
                            sbase_n0=float(row['sbase_n0']),
                            sigma=float(row['sigma']),
                            temp_lapse=[float(row[f'tmp_lapse_jan']), float(row[f'tmp_lapse_feb']), float(row[f'tmp_lapse_mar']),
                                      float(row[f'tmp_lapse_apr']), float(row[f'tmp_lapse_may']), float(row[f'tmp_lapse_jun']),
                                      float(row[f'tmp_lapse_jul']), float(row[f'tmp_lapse_aug']), float(row[f'tmp_lapse_sep']),
                                      float(row[f'tmp_lapse_oct']), float(row[f'tmp_lapse_nov']), float(row[f'tmp_lapse_dec'])],
                            prcp_lapse=[float(row[f'prcp_lapse_jan']), float(row[f'prcp_lapse_feb']), float(row[f'prcp_lapse_mar']),
                                       float(row[f'prcp_lapse_apr']), float(row[f'prcp_lapse_may']), float(row[f'prcp_lapse_jun']),
                                       float(row[f'prcp_lapse_jul']), float(row[f'prcp_lapse_aug']), float(row[f'prcp_lapse_sep']),
                                       float(row[f'prcp_lapse_oct']), float(row[f'prcp_lapse_nov']), float(row[f'prcp_lapse_dec'])]
                        )

                        sites.append(site)

        except Exception as e:
            fatal_error(f"Error reading site data file: {e}")

        return sites

    def read_climate(self, sites: List[SiteData]):
        """
        Read climate data and attach to sites.

        Args:
            sites: List of sites to attach climate data to
        """
        climate_file = self.filenames.get('climate')
        if not climate_file or not os.path.exists(climate_file):
            fatal_error('Unable to find climate file')

        try:
            with open(climate_file, 'r') as f:
                reader = csv.DictReader(f)

                for row in reader:
                    if not row['site']:  # Skip empty rows
                        continue
                    site_id = int(row['site'])  # 'site' maps to site_id

                    # Find matching site
                    site = next((s for s in sites if s.site_id == site_id), None)
                    if site:
                        # Read temperature and precipitation data using UVAFME column names
                        months = ['jan', 'feb', 'mar', 'apr', 'may', 'jun',
                                 'jul', 'aug', 'sep', 'oct', 'nov', 'dec']
                        tmin = [float(row[f'tmin_{month}']) for month in months]
                        tmax = [float(row[f'tmax_{month}']) for month in months]
                        prcp = [float(row[f'prcp_{month}']) for month in months]

                        self.attach_climate(site, tmin, tmax, prcp)

        except Exception as e:
            fatal_error(f"Error reading climate file: {e}")

        # Validate climate data
        for site in sites:
            if np.any(site.tmin == RNVALID) or site.site_wmo == RNVALID:
                print(f'No climate data for site number {site.site_id}')
                site.site_wmo = RNVALID

    def attach_climate(self, site: SiteData, tmin: List[float],
                      tmax: List[float], prcp: List[float]):
        """
        Attach climate data to a site.

        Args:
            site: Site to attach climate data to
            tmin: Monthly minimum temperatures
            tmax: Monthly maximum temperatures
            prcp: Monthly precipitation
        """
        site.tmin = np.array(tmin)
        site.tmax = np.array(tmax)
        site.precip = np.array(prcp)

    def read_climate_stds(self, sites: List[SiteData]):
        """
        Read climate standard deviations and attach to sites.

        Args:
            sites: List of sites to attach climate variability data to
        """
        climate_std_file = self.filenames.get('climate_std')
        if not climate_std_file or not os.path.exists(climate_std_file):
            print("Warning: Climate standard deviation file not found, using zero variability")
            # Set all std values to zero
            for site in sites:
                site.tmin_std = np.zeros(12)
                site.tmax_std = np.zeros(12)
                site.precip_std = np.zeros(12)
            return

        try:
            with open(climate_std_file, 'r') as f:
                reader = csv.DictReader(f)

                for row in reader:
                    if not row['site']:  # Skip empty rows
                        continue
                    site_id = int(row['site'])

                    # Find matching site
                    site = next((s for s in sites if s.site_id == site_id), None)
                    if site:
                        # Read temperature and precipitation standard deviations
                        months = ['jan', 'feb', 'mar', 'apr', 'may', 'jun',
                                 'jul', 'aug', 'sep', 'oct', 'nov', 'dec']
                        tmin_std = [float(row[f'tmn_std_{month}']) for month in months]
                        tmax_std = [float(row[f'tmx_std_{month}']) for month in months]
                        prcp_std = [float(row[f'prcp_std_{month}']) for month in months]

                        self.attach_climate_std(site, tmin_std, tmax_std, prcp_std)

        except Exception as e:
            print(f"Warning: Error reading climate std file: {e}")
            # Set all std values to zero on error
            for site in sites:
                site.tmin_std = np.zeros(12)
                site.tmax_std = np.zeros(12)
                site.precip_std = np.zeros(12)

    def attach_climate_std(self, site: SiteData, tmin_std: List[float],
                          tmax_std: List[float], prcp_std: List[float]):
        """
        Attach climate standard deviations to a site.

        Args:
            site: Site to attach climate variability data to
            tmin_std: Monthly minimum temperature standard deviations
            tmax_std: Monthly maximum temperature standard deviations
            prcp_std: Monthly precipitation standard deviations
        """
        site.tmin_std = np.array(tmin_std)
        site.tmax_std = np.array(tmax_std)
        site.precip_std = np.array(prcp_std) * 0.1  # Convert mm to cm

    def read_species_data(self) -> List[SpeciesData]:
        """
        Read species data from species file.

        Returns:
            List of SpeciesData objects
        """
        species_file = self.filenames.get('species')
        if not species_file or not os.path.exists(species_file):
            fatal_error('Unable to find species-list file')

        species_data = []

        try:
            with open(species_file, 'r') as f:
                reader = csv.DictReader(f)

                # Debug: print available columns
                sample_row = next(reader)
                f.seek(0)
                reader = csv.DictReader(f)
                print(f"Available CSV columns: {list(sample_row.keys())}")

                for row in reader:
                    species = SpeciesData()

                    # Helper function to safely get values with defaults
                    def safe_get(key, default_val=0, cast_func=str):
                        try:
                            if key in row and row[key].strip():
                                val = cast_func(row[key])
                                # Debug D_L reading
                                if key == 'D_L' and len(species_data) < 2:
                                    print(f"DEBUG: Reading D_L for species {len(species_data)+1}: raw='{row[key]}', parsed={val}")
                                return val
                            else:
                                print(f"Warning: Missing or empty column '{key}', using default {default_val}")
                                return default_val
                        except (ValueError, KeyError) as e:
                            print(f"Warning: Error reading column '{key}': {e}, using default {default_val}")
                            return default_val

                    # Map UVAFME CSV column names to expected names
                    # Initialize species from CSV data with mapped column access
                    species.initialize_species(
                        species_id=safe_get('Individual', 0, int),  # 'Individual' maps to species_id
                        genus_name=safe_get('Genus', 'Unknown'),  # 'Genus' maps to genus_name
                        taxonomic_name=safe_get('Scientific name', 'Unknown'),  # 'Scientific name' maps to taxonomic_name
                        unique_id=safe_get('Species_code', 'UNKNOWN'),  # 'Species_code' maps to unique_id
                        common_name=safe_get('Common name', 'Unknown'),  # 'Common name' maps to common_name
                        genus_id=safe_get('Group', 1, int),  # 'Group' maps to genus_id
                        shade_tol=safe_get('l', 3, int),  # 'l' (light tolerance) is shade_tol (1-5)
                        lownutr_tol=safe_get('n', 3, int),  # 'n' (nutrient tolerance) maps to lownutr_tol
                        stress_tol=safe_get('stress', 3, int),  # 'stress' maps to stress_tol
                        age_tol=safe_get('old', 3, int),  # 'old' (age tolerance) maps to age_tol
                        drought_tol=safe_get('d', 3, int),  # 'd' (drought tolerance) maps to drought_tol
                        flood_tol=safe_get('f', 3, int),  # 'f' (flood tolerance) maps to flood_tol
                        fire_tol=safe_get('fire', 3, int),  # 'fire' maps to fire_tol
                        max_age=safe_get('AGEmax', 300.0, float),  # 'AGEmax' maps to max_age
                        max_diam=safe_get('DBHmax', 100.0, float),  # 'DBHmax' maps to max_diam
                        max_ht=safe_get('Hmax', 30.0, float),  # 'Hmax' maps to max_ht
                        wood_bulk_dens=safe_get('bulk', 0.5, float),  # 'bulk' maps to wood_bulk_dens
                        rootdepth=self.parameters.rootdepth,
                        leafdiam_a=safe_get('D_L', 0.5, float),  # 'D_L' maps to leafdiam_a
                        leafarea_c=safe_get('L_C', 15.0, float),  # 'L_C' maps to leafarea_c
                        deg_day_min=safe_get('DEGDmin', 500.0, float),  # 'DEGDmin' maps to deg_day_min
                        deg_day_opt=safe_get('DEGDoptimum', 1200.0, float),  # 'DEGDoptimum' maps to deg_day_opt
                        deg_day_max=safe_get('DEGDmax', 2500.0, float),  # 'DEGDmax' maps to deg_day_max
                        seedling_lg=safe_get('s', 0.8, float),  # 's' is seedling growth parameter (arfa_0 shape param)
                        invader=safe_get('invader', 0.1, float),  # 'invader' maps to invader
                        seed_num=safe_get('seed', 100.0, float),  # 'seed' maps to seed_num
                        sprout_num=safe_get('sprout', 50.0, float),  # 'sprout' maps to sprout_num
                        seed_surv=safe_get('NDS', 0.3, float),  # 'NDS' (seedling survival) maps to seed_surv
                        arfa_0=safe_get('NDE', 0.8, float),  # 'NDE' maps to arfa_0
                        g=safe_get('g', 0.6, float),  # 'g' maps to g
                        conifer=(safe_get('evergreen', 0, int) == 1)  # 'evergreen' maps to conifer
                    )

                    # Debug: Print first two species leafdiam_a
                    if len(species_data) < 2:
                        print(f"DEBUG: Species {len(species_data)+1} ({species.genus_name}) initialized with leafdiam_a={species.leafdiam_a}")

                    species_data.append(species)

        except Exception as e:
            fatal_error(f"Error reading species data: {e}")

        print(f'Species data initialized. Total read in: {len(species_data)}')
        return species_data


# Global input manager instance
input_manager = InputFileManager()


# Convenience functions for compatibility
def initialize_input_files(filelist: str = ""):
    """Initialize input files and parameters."""
    input_manager.initialize_input_files(filelist)


def read_sitelist() -> List[int]:
    """Read list of site IDs."""
    return input_manager.read_sitelist()


def read_sites(site_ids: List[int]) -> List[SiteData]:
    """Read site data."""
    return input_manager.read_sites(site_ids)


def read_climate(sites: List[SiteData]):
    """Read and attach climate data to sites."""
    input_manager.read_climate(sites)


def read_species_data() -> List[SpeciesData]:
    """Read species data."""
    return input_manager.read_species_data()