"""
Input/Output utilities for UVAFME vegetation model.
Handles reading and writing of climate, species, and site data.
"""

import os
import csv
import json
import numpy as np
from typing import List, Dict, Tuple, Optional
from .constants import *
from .species import SpeciesData
from .site import SiteData
from .parameters import Parameters


def count_records(file_handle, nheaders: int = 0) -> int:
    """
    Count the number of records in a file, excluding headers.

    Args:
        file_handle: Open file handle or file path
        nheaders: Number of header lines to skip

    Returns:
        Number of data records
    """
    if isinstance(file_handle, str):
        with open(file_handle, 'r') as f:
            return count_records(f, nheaders)

    # Save current position
    current_pos = file_handle.tell()

    # Go to beginning
    file_handle.seek(0)

    # Count lines
    line_count = 0
    for line in file_handle:
        line_count += 1

    # Restore position
    file_handle.seek(current_pos)

    return max(0, line_count - nheaders)


def split_line(line: str) -> List[str]:
    """
    Split a line by commas and return fields.

    Args:
        line: Input line to split

    Returns:
        List of fields
    """
    # Strip whitespace and split by comma
    line = line.strip()
    if not line:
        return []

    # Simple CSV splitting - could be enhanced for quoted fields
    fields = [field.strip() for field in line.split(',')]
    return fields


def quote_strip(string: str) -> str:
    """
    Remove quotes (single or double) from a string.

    Args:
        string: Input string

    Returns:
        String with quotes removed
    """
    return string.strip().strip('"').strip("'")


def warning(message: str):
    """Print a warning message."""
    print(f"Warning: {message}")


def fatal_error(message: str):
    """Print an error message and exit."""
    print(f"Fatal Error: {message}")
    import sys
    sys.exit(1)


class UVAFMEReader:
    """Handles reading of UVAFME input files."""
    
    def __init__(self, base_path: str = "input_data"):
        self.base_path = base_path
        self.ensure_directories()
    
    def ensure_directories(self):
        """Ensure input directories exist."""
        os.makedirs(self.base_path, exist_ok=True)
    
    def read_species_file(self, filename: str = "species.csv") -> List[SpeciesData]:
        """Read species data from CSV file."""
        filepath = os.path.join(self.base_path, filename)
        species_list = []
        
        if not os.path.exists(filepath):
            print(f"Warning: Species file {filepath} not found, creating example")
            self.create_example_species_file(filepath)
        
        with open(filepath, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            
            # Peek at first row to determine format
            first_row = next(reader)
            csvfile.seek(0)
            reader = csv.DictReader(csvfile)
            
            if 'species_id' in first_row:
                # New format
                for row in reader:
                    species = self._parse_new_format(row)
                    species_list.append(species)
            else:
                # Original UVAFME format
                for row in reader:
                    species = self._parse_uvafme_format(row)
                    species_list.append(species)
        
        return species_list
    
    def _parse_new_format(self, row):
        """Parse species data from new format."""
        species = SpeciesData()
        species.initialize_species(
            species_id=int(row['species_id']),
            genus_name=row['genus_name'],
            taxonomic_name=row['taxonomic_name'],
            unique_id=row['unique_id'],
            common_name=row['common_name'],
            genus_id=int(row['genus_id']),
            shade_tol=int(row['shade_tol']),
            lownutr_tol=int(row['lownutr_tol']),
            stress_tol=int(row['stress_tol']),
            age_tol=int(row['age_tol']),
            drought_tol=int(row['drought_tol']),
            flood_tol=int(row['flood_tol']),
            fire_tol=int(row['fire_tol']),
            max_age=float(row['max_age']),
            max_diam=float(row['max_diam']),
            max_ht=float(row['max_ht']),
            wood_bulk_dens=float(row['wood_bulk_dens']),
            rootdepth=float(row['rootdepth']),
            leafdiam_a=float(row['leafdiam_a']),
            leafarea_c=float(row['leafarea_c']),
            deg_day_min=float(row['deg_day_min']),
            deg_day_opt=float(row['deg_day_opt']),
            deg_day_max=float(row['deg_day_max']),
            seedling_lg=float(row['seedling_lg']),
            invader=float(row['invader']),
            seed_num=float(row['seed_num']),
            sprout_num=float(row['sprout_num']),
            seed_surv=float(row['seed_surv']),
            arfa_0=float(row['arfa_0']),
            g=float(row['g']),
            conifer=row['conifer'].lower() == 'true'
        )
        return species
    
    def _parse_uvafme_format(self, row):
        """Parse species data from original UVAFME format."""
        species = SpeciesData()
        
        # Map UVAFME columns to our format
        genus_name = row['Genus'].strip("'")
        taxonomic_name = row['Scientific name'].strip("'")
        common_name = row['Common name'].strip("'")
        unique_id = row['Species_code']
        
        species.initialize_species(
            species_id=int(row['Individual']),
            genus_name=genus_name,
            taxonomic_name=taxonomic_name,
            unique_id=unique_id,
            common_name=common_name,
            genus_id=int(row['Group']),
            shade_tol=int(row['l']),  # light tolerance
            lownutr_tol=int(row['n']),  # nutrient tolerance
            stress_tol=int(row['stress']),
            age_tol=int(row['old']),
            drought_tol=int(row['d']),
            flood_tol=int(row['f']),
            fire_tol=int(row['fire']),
            max_age=float(row['AGEmax']),
            max_diam=float(row['DBHmax']),
            max_ht=float(row['Hmax']),
            wood_bulk_dens=float(row['bulk']),
            rootdepth=2.0,  # Default since not in original
            leafdiam_a=float(row['D_L']),
            leafarea_c=float(row['L_C']),
            deg_day_min=float(row['DEGDmin']),
            deg_day_opt=float(row['DEGDoptimum']),
            deg_day_max=float(row['DEGDmax']),
            seedling_lg=float(row['NDS']),
            invader=float(row['invader']),
            seed_num=float(row['seed']),
            sprout_num=float(row['sprout']),
            seed_surv=float(row['NDE']),
            arfa_0=float(row['s']),
            g=float(row['g']),
            conifer=int(row['evergreen']) == 1
        )
        return species
    
    def read_site_file(self, filename: str = "sites.csv") -> List[SiteData]:
        """Read site data from CSV file."""
        filepath = os.path.join(self.base_path, filename)
        sites_list = []
        
        if not os.path.exists(filepath):
            print(f"Warning: Site file {filepath} not found, creating example")
            self.create_example_site_file(filepath)
        
        with open(filepath, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            
            # Check format by looking at first row
            first_row = next(reader)
            csvfile.seek(0)
            reader = csv.DictReader(csvfile)
            
            if 'site_id' in first_row:
                # New format
                for row in reader:
                    sites_list.append(self._parse_new_site_format(row))
            else:
                # Original UVAFME format
                for row in reader:
                    if row['site']:  # Skip empty rows
                        sites_list.append(self._parse_uvafme_site_format(row))
        
        return sites_list
    
    def _parse_new_site_format(self, row):
        """Parse site data from new format."""
        site = SiteData()
        
        # Parse temperature and precipitation lapse rates
        temp_lapse = [float(x) for x in row['temp_lapse'].split(',')]
        prcp_lapse = [float(x) for x in row['prcp_lapse'].split(',')]
        
        site.initialize_site(
            siteid=int(row['site_id']),
            sitename=row['site_name'],
            siteregion=row['region'],
            lat=float(row['latitude']),
            long=float(row['longitude']),
            wmo=float(row['site_wmo']),
            elevation=float(row['elevation']),
            slope=float(row['slope']),
            Afc=float(row['A_field_cap']),
            A_perm_wp=float(row['A_perm_wp']),
            lai=float(row['leaf_area_ind']),
            base_h=float(row['base_h']),
            lai_w0=float(row['leaf_area_w0']),
            A0_w0=float(row['A0_w0']),
            A_w0=float(row['A_w0']),
            sbase_w0=float(row['sbase_w0']),
            fire_prob=float(row['fire_prob']),
            wind_prob=float(row['wind_prob']),
            A0_c0=float(row['A0_c0']),
            A0_n0=float(row['A0_n0']),
            A_c0=float(row['A_c0']),
            A_n0=float(row['A_n0']),
            sbase_c0=float(row['sbase_c0']),
            sbase_n0=float(row['sbase_n0']),
            sigma=float(row['sigma']),
            temp_lapse=temp_lapse,
            prcp_lapse=prcp_lapse
        )
        return site
    
    def _parse_uvafme_site_format(self, row):
        """Parse site data from original UVAFME format."""
        site = SiteData()
        
        # Parse temperature lapse rates (monthly)
        temp_lapse = [
            float(row['tmp_lapse_jan']), float(row['tmp_lapse_feb']), float(row['tmp_lapse_mar']),
            float(row['tmp_lapse_apr']), float(row['tmp_lapse_may']), float(row['tmp_lapse_jun']),
            float(row['tmp_lapse_jul']), float(row['tmp_lapse_aug']), float(row['tmp_lapse_sep']),
            float(row['tmp_lapse_oct']), float(row['tmp_lapse_nov']), float(row['tmp_lapse_dec'])
        ]
        
        # Parse precipitation lapse rates (monthly)
        prcp_lapse = [
            float(row['prcp_lapse_jan']), float(row['prcp_lapse_feb']), float(row['prcp_lapse_mar']),
            float(row['prcp_lapse_apr']), float(row['prcp_lapse_may']), float(row['prcp_lapse_jun']),
            float(row['prcp_lapse_jul']), float(row['prcp_lapse_aug']), float(row['prcp_lapse_sep']),
            float(row['prcp_lapse_oct']), float(row['prcp_lapse_nov']), float(row['prcp_lapse_dec'])
        ]
        
        site.initialize_site(
            siteid=int(row['site']),
            sitename=row['name'],
            siteregion=row['region'],
            lat=float(row['latitude']),
            long=float(row['longitude']),
            wmo=float(row['wmo']),
            elevation=float(row['elevation']),
            slope=float(row['slope']),
            Afc=float(row['soilA_field_cap']),
            A_perm_wp=float(row['soilA_perm_wp']),
            lai=float(row['lai']),
            base_h=float(row['soil_base_h']),
            lai_w0=float(row['lai_w0']),
            A0_w0=float(row['soilAO_w0']),
            A_w0=float(row['soilA_w0']),
            sbase_w0=float(row['sbase_w0']),
            fire_prob=float(row['fire_prob']),
            wind_prob=float(row['wind_prob']),
            A0_c0=float(row['soilAO_c0']),
            A0_n0=float(row['soilAO_n0']),
            A_c0=float(row['soilA_c0']),
            A_n0=float(row['soilA_n0']),
            sbase_c0=float(row['sbase_c0']),
            sbase_n0=float(row['sbase_n0']),
            sigma=float(row['sigma']),
            temp_lapse=temp_lapse,
            prcp_lapse=prcp_lapse
        )
        return site
    
    def read_climate_file(self, filename: str = "climate.csv") -> Dict[int, Dict[str, np.ndarray]]:
        """Read climate data from CSV file."""
        filepath = os.path.join(self.base_path, filename)
        climate_data = {}
        
        if not os.path.exists(filepath):
            print(f"Warning: Climate file {filepath} not found, creating example")
            self.create_example_climate_file(filepath)
        
        with open(filepath, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            
            # Check format by looking at first row
            first_row = next(reader)
            csvfile.seek(0)
            reader = csv.DictReader(csvfile)
            
            if 'site_id' in first_row:
                # New format
                for row in reader:
                    climate_data.update(self._parse_new_climate_format(row))
            else:
                # Original UVAFME format
                for row in reader:
                    if row['site']:  # Skip empty rows
                        climate_data.update(self._parse_uvafme_climate_format(row))
        
        return climate_data
    
    def _parse_new_climate_format(self, row):
        """Parse climate data from new format."""
        site_id = int(row['site_id'])
        
        # Parse monthly data
        tmin = [float(row[f'tmin_{i}']) for i in range(1, 13)]
        tmax = [float(row[f'tmax_{i}']) for i in range(1, 13)]
        precip = [float(row[f'precip_{i}']) for i in range(1, 13)]
        
        # Optional standard deviations
        tmin_std = [float(row.get(f'tmin_std_{i}', 0)) for i in range(1, 13)]
        tmax_std = [float(row.get(f'tmax_std_{i}', 0)) for i in range(1, 13)]
        precip_std = [float(row.get(f'precip_std_{i}', 0)) for i in range(1, 13)]
        
        return {site_id: {
            'tmin': np.array(tmin),
            'tmax': np.array(tmax),
            'precip': np.array(precip),
            'tmin_std': np.array(tmin_std),
            'tmax_std': np.array(tmax_std),
            'precip_std': np.array(precip_std)
        }}
    
    def _parse_uvafme_climate_format(self, row):
        """Parse climate data from original UVAFME format."""
        site_id = int(row['site'])
        
        # Parse monthly temperature data
        tmin = [
            float(row['tmin_jan']), float(row['tmin_feb']), float(row['tmin_mar']),
            float(row['tmin_apr']), float(row['tmin_may']), float(row['tmin_jun']),
            float(row['tmin_jul']), float(row['tmin_aug']), float(row['tmin_sep']),
            float(row['tmin_oct']), float(row['tmin_nov']), float(row['tmin_dec'])
        ]
        
        tmax = [
            float(row['tmax_jan']), float(row['tmax_feb']), float(row['tmax_mar']),
            float(row['tmax_apr']), float(row['tmax_may']), float(row['tmax_jun']),
            float(row['tmax_jul']), float(row['tmax_aug']), float(row['tmax_sep']),
            float(row['tmax_oct']), float(row['tmax_nov']), float(row['tmax_dec'])
        ]
        
        precip = [
            float(row['prcp_jan']), float(row['prcp_feb']), float(row['prcp_mar']),
            float(row['prcp_apr']), float(row['prcp_may']), float(row['prcp_jun']),
            float(row['prcp_jul']), float(row['prcp_aug']), float(row['prcp_sep']),
            float(row['prcp_oct']), float(row['prcp_nov']), float(row['prcp_dec'])
        ]
        
        # Standard deviations default to 0 for original format
        tmin_std = [0.0] * 12
        tmax_std = [0.0] * 12
        precip_std = [0.0] * 12
        
        return {site_id: {
            'tmin': np.array(tmin),
            'tmax': np.array(tmax),
            'precip': np.array(precip),
            'tmin_std': np.array(tmin_std),
            'tmax_std': np.array(tmax_std),
            'precip_std': np.array(precip_std)
        }}
    
    def read_filelist(self, filename: str = "filelist.txt") -> Dict[str, str]:
        """Read file list configuration."""
        filepath = os.path.join(self.base_path, filename)
        filelist = {}
        
        if not os.path.exists(filepath):
            print(f"Warning: Filelist {filepath} not found, creating example")
            self.create_example_filelist(filepath)
        
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    key, value = line.split('=', 1)
                    filelist[key.strip()] = value.strip()
        
        return filelist
    
    def create_example_species_file(self, filepath: str):
        """Create example species file."""
        example_species = [
            {
                'species_id': 1, 'genus_name': 'Quercus', 'taxonomic_name': 'Quercus alba',
                'unique_id': 'QUAL', 'common_name': 'White Oak', 'genus_id': 1,
                'shade_tol': 3, 'lownutr_tol': 2, 'stress_tol': 2, 'age_tol': 5,
                'drought_tol': 2, 'flood_tol': 2, 'fire_tol': 2, 'max_age': 300.0,
                'max_diam': 150.0, 'max_ht': 35.0, 'wood_bulk_dens': 0.68,
                'rootdepth': 2.0, 'leafdiam_a': 0.1, 'leafarea_c': 50.0,
                'deg_day_min': 1000.0, 'deg_day_opt': 2500.0, 'deg_day_max': 4000.0,
                'seedling_lg': 0.5, 'invader': 0.0, 'seed_num': 100.0,
                'sprout_num': 50.0, 'seed_surv': 0.1, 'arfa_0': 0.8, 'g': 0.3,
                'conifer': 'false'
            },
            {
                'species_id': 2, 'genus_name': 'Pinus', 'taxonomic_name': 'Pinus strobus',
                'unique_id': 'PIST', 'common_name': 'White Pine', 'genus_id': 2,
                'shade_tol': 2, 'lownutr_tol': 3, 'stress_tol': 3, 'age_tol': 4,
                'drought_tol': 3, 'flood_tol': 1, 'fire_tol': 1, 'max_age': 250.0,
                'max_diam': 120.0, 'max_ht': 40.0, 'wood_bulk_dens': 0.35,
                'rootdepth': 1.5, 'leafdiam_a': 0.15, 'leafarea_c': 40.0,
                'deg_day_min': 1200.0, 'deg_day_opt': 2800.0, 'deg_day_max': 4200.0,
                'seedling_lg': 0.3, 'invader': 0.0, 'seed_num': 200.0,
                'sprout_num': 0.0, 'seed_surv': 0.05, 'arfa_0': 0.9, 'g': 0.4,
                'conifer': 'true'
            }
        ]
        
        with open(filepath, 'w', newline='') as csvfile:
            fieldnames = list(example_species[0].keys())
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(example_species)
    
    def create_example_site_file(self, filepath: str):
        """Create example site file."""
        example_site = {
            'site_id': 1, 'site_name': 'Test Site', 'region': 'Virginia',
            'latitude': 37.0, 'longitude': -78.0, 'site_wmo': 123456,
            'elevation': 200.0, 'slope': 5.0, 'A_field_cap': 20.0,
            'A_perm_wp': 8.0, 'leaf_area_ind': 3.0, 'base_h': 70.0,
            'leaf_area_w0': 0.5, 'A0_w0': 2.0, 'A_w0': 15.0, 'sbase_w0': 40.0,
            'fire_prob': 0.01, 'wind_prob': 0.05, 'A0_c0': 10.0, 'A0_n0': 1.0,
            'A_c0': 50.0, 'A_n0': 5.0, 'sbase_c0': 100.0, 'sbase_n0': 10.0,
            'sigma': 0.1, 'temp_lapse': ','.join(['0.6'] * 12),
            'prcp_lapse': ','.join(['0.1'] * 12)
        }
        
        with open(filepath, 'w', newline='') as csvfile:
            fieldnames = list(example_site.keys())
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow(example_site)
    
    def create_example_climate_file(self, filepath: str):
        """Create example climate file."""
        # Example climate data for Virginia
        tmin = [-5, -2, 3, 8, 13, 18, 20, 19, 15, 9, 4, -1]
        tmax = [5, 8, 15, 22, 27, 32, 35, 33, 28, 22, 14, 7]
        precip = [80, 70, 90, 100, 120, 110, 90, 85, 75, 65, 70, 85]
        
        example_climate = {'site_id': 1}
        
        for i, (tmin_val, tmax_val, precip_val) in enumerate(zip(tmin, tmax, precip), 1):
            example_climate[f'tmin_{i}'] = tmin_val
            example_climate[f'tmax_{i}'] = tmax_val
            example_climate[f'precip_{i}'] = precip_val
        
        with open(filepath, 'w', newline='') as csvfile:
            fieldnames = list(example_climate.keys())
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow(example_climate)
    
    def create_example_filelist(self, filepath: str):
        """Create example filelist file."""
        example_filelist = """# UVAFME Input File List
# Format: key=value
species_file=species.csv
site_file=sites.csv
climate_file=climate.csv
config_file=uvafme_config.json
"""
        with open(filepath, 'w') as f:
            f.write(example_filelist)


class UVAFMEWriter:
    """Handles writing of UVAFME output files."""
    
    def __init__(self, base_path: str = None):
        # Support multiple possible locations for output_data directory
        if base_path is None:
            possible_output_paths = ['output_data', '../output_data', '../../output_data']
            self.base_path = None

            for path in possible_output_paths:
                if os.path.exists(path):
                    self.base_path = path
                    break

            # If none exist, use '../output_data' as the preferred default
            if self.base_path is None:
                self.base_path = '../output_data'
        else:
            self.base_path = base_path

        self.ensure_directories()
        self.output_files = {}
    
    def ensure_directories(self):
        """Ensure output directories exist."""
        os.makedirs(self.base_path, exist_ok=True)
    
    def initialize_output_files(self, species_present=None):
        """Initialize output files for writing."""
        # Species-level output
        species_file = os.path.join(self.base_path, "species_data.csv")
        self.output_files['species'] = open(species_file, 'w', newline='')
        species_writer = csv.writer(self.output_files['species'])
        species_writer.writerow(['year', 'site_id', 'species_id', 'genus_name', 
                                'biomass_c', 'biomass_n', 'basal_area', 'n_trees'])
        
        # Site-level output
        site_file = os.path.join(self.base_path, "site_data.csv")
        self.output_files['site'] = open(site_file, 'w', newline='')
        site_writer = csv.writer(self.output_files['site'])
        site_writer.writerow(['year', 'site_id', 'rain', 'pot_evap', 'act_evap',
                             'grow_days', 'deg_days', 'dry_days_upper', 'dry_days_base',
                             'flood_days'])
        
        # Soil-level output
        soil_file = os.path.join(self.base_path, "soil_data.csv")
        self.output_files['soil'] = open(soil_file, 'w', newline='')
        soil_writer = csv.writer(self.output_files['soil'])
        soil_writer.writerow(['year', 'site_id', 'A0_c0', 'A_c0', 'A0_n0', 'A_n0',
                             'BL_c0', 'BL_n0', 'total_C_rsp', 'avail_N'])
    
    def write_species_data(self, site, species_present, year):
        """Write species-level data."""
        if 'species' not in self.output_files:
            return
        
        writer = csv.writer(self.output_files['species'])
        for plot in site.plots:
            for species in plot.species:
                # Get species metrics from plot
                metrics = plot.sum_over_sg([species.unique_id], field='species')
                
                writer.writerow([
                    year, site.site_id, species.species_id, species.genus_name,
                    metrics['biomC'][0], metrics['biomN'][0], 
                    metrics['basal_area'][0], metrics['n_trees'][0]
                ])
    
    def write_site_data(self, site, year):
        """Write site-level data."""
        if 'site' not in self.output_files:
            return
        
        writer = csv.writer(self.output_files['site'])
        writer.writerow([
            year, site.site_id, site.rain, site.pot_evap_day, site.act_evap_day,
            site.grow_days, site.deg_days, site.dry_days_upper_layer,
            site.dry_days_base_layer, site.flood_days
        ])
    
    def write_soil_data(self, site, year):
        """Write soil-level data."""
        if 'soil' not in self.output_files:
            return
        
        writer = csv.writer(self.output_files['soil'])
        soil = site.soil
        writer.writerow([
            year, site.site_id, soil.A0_c0, soil.A_c0, soil.A0_n0, soil.A_n0,
            soil.BL_c0, soil.BL_n0, soil.total_C_rsp, soil.avail_N
        ])
    
    def write_tree_data(self, site, year):
        """Write tree-level data (if enabled)."""
        tree_file = os.path.join(self.base_path, f"trees_year_{year}.csv")
        
        with open(tree_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['year', 'site_id', 'plot_id', 'tree_id', 'species_id',
                            'genus_name', 'dbh', 'height', 'biomass_c', 'biomass_n',
                            'leaf_biomass', 'mort_marker'])
            
            for plot_id, plot in enumerate(site.plots):
                for tree_id, tree in enumerate(plot.trees):
                    writer.writerow([
                        year, site.site_id, plot_id, tree_id, tree.species_id,
                        tree.genus_name, tree.diam_bht, tree.forska_ht,
                        tree.biomC, tree.biomN, tree.leaf_bm, tree.mort_marker
                    ])
    
    def close_output_files(self):
        """Close all output files."""
        for file_handle in self.output_files.values():
            file_handle.close()
        self.output_files.clear()