"""
Output module for UVAFME vegetation model.

This module handles writing model output data to CSV files including
genus data, species data, site data, soil data, and tree data.

Translated from Output.f90
"""

import os
from typing import List, Dict, TextIO
import numpy as np

from .constants import *
from .site import SiteData
from .genus_groups import Groups
from .csv_file import CSVWriter, csv_write, register_csv_file
from .plot import PlotData
from .utilities import stddev


class OutputManager:
    """
    Manages all output file operations for the UVAFME model.

    This class handles opening, writing, and closing output files
    including genus data, species data, site data, and tree data.
    """

    def __init__(self, params=None):
        self.file_handles: Dict[str, TextIO] = {}
        self.csv_writers: Dict[str, CSVWriter] = {}
        self.params = params  # ModelParameters instance

        # Support multiple possible locations for output_data directory
        possible_output_paths = ['output_data', '../output_data', '../../output_data']
        self.output_dir = None

        for path in possible_output_paths:
            if os.path.exists(path):
                self.output_dir = path
                break

        # If none exist, use '../output_data' as the preferred default
        if self.output_dir is None:
            self.output_dir = '../output_data'

        # Plot adjustments (global variables from Fortran)
        self.plotscale = 0.0
        self.plotadj = 0.0
        self.plotrenorm = 0.0

    def initialize_output_files(self, species_present: Groups):
        """
        Initialize all output files and write headers.

        Args:
            species_present: Groups object containing species information
        """
        self.open_output_files()
        self.write_headers(species_present)

    def open_output_files(self):
        """Open all required output files."""
        # Ensure output directory exists
        os.makedirs(self.output_dir, exist_ok=True)

        # Define output files
        output_files = {
            'biom_by_g': 'genus_data.csv',
            'biom_by_s': 'species_data.csv',
            'clim_unit': 'site_data.csv',
            'c_and_n': 'soil_data.csv',
            'tld': 'tree_data.csv',
            'pl_biom_by_g': 'plot_genus_data.csv',
            'pl_biom_by_s': 'plot_species_data.csv'
        }

        # Open files and register with CSV writers
        for file_key, filename in output_files.items():
            filepath = os.path.join(self.output_dir, filename)
            self.file_handles[file_key] = open(filepath, 'w', newline='')
            self.csv_writers[file_key] = CSVWriter(self.file_handles[file_key])
            register_csv_file(file_key, self.file_handles[file_key])

    def write_headers(self, species_present: Groups):
        """
        Write header rows for all output files.

        Args:
            species_present: Groups object containing species information
        """
        # Genus data header (matching Fortran IO_Utils.f90 lines 224-245)
        genus_header = [
            'siteID', 'year', 'genus',
            '<0', '0-8', '8-28', '-48',
            '-68', '-88', '>88',
            'max_diam', 'max_hgt', 'leaf_area_ind',
            'basal_area', 'total_biomC', 'pl_biomC_std',
            'total_biomN', 'pl_biomN_std'
        ]
        csv_write(self.file_handles['biom_by_g'], genus_header)

        # Species data header (matching Fortran IO_Utils.f90 lines 248-266)
        species_header = [
            'siteID', 'year', 'genus', 'species',
            '<0', '0-8', '8-28', '-48',
            '-68', '-88', '>88',
            'max_diam', 'max_hgt', 'leaf_area_ind',
            'basal_area', 'total_biomC', 'pl_biomC_std',
            'total_biomN', 'pl_biomN_std'
        ]
        csv_write(self.file_handles['biom_by_s'], species_header)

        # Site data header
        site_header = ['site_id', 'year'] + self.get_site_csv_header()
        csv_write(self.file_handles['clim_unit'], site_header)

        # Soil data header
        soil_header = ['site_id', 'year'] + self.get_soil_csv_header()
        csv_write(self.file_handles['c_and_n'], soil_header)

        # Tree data header (if enabled)
        tree_header = ['site_id', 'year', 'plot', 'tree'] + self.get_tree_csv_header()
        csv_write(self.file_handles['tld'], tree_header)

    def get_site_csv_header(self) -> List[str]:
        """Get header for site data CSV."""
        return [
            'latitude', 'longitude', 'elevation', 'slope',
            'leaf_area_ind', 'grow_days', 'deg_days',
            'flood_days', 'dry_days_upper', 'dry_days_base',
            'pot_evap_day', 'act_evap_day', 'rain'
        ]

    def get_soil_csv_header(self) -> List[str]:
        """Get header for soil data CSV matching Fortran IO_Utils.f90 lines 198-213."""
        return [
            'a0c0', 'ac0', 'a0n0', 'an0', 'bc0', 'bn0',
            'soilresp', 'biomassC', 'C_into_A0', 'net_C_into_A0', 'net_prim_prodC',
            'biomassN', 'N_into_A0', 'net_N_into_A0', 'net_prim_prodN', 'avail_n'
        ]

    def get_tree_csv_header(self) -> List[str]:
        """Get header for tree data CSV."""
        return [
            'species_id', 'dbh', 'height', 'biomass_c',
            'biomass_n', 'leaf_biomass', 'age', 'crown_area'
        ]

    def write_genus_data(self, site: SiteData, species_pres: Groups, year: int):
        """
        Write genus-level biomass and demographic data.

        Args:
            site: Site data
            species_pres: Groups containing species information
            year: Current year
        """
        # Calculate plot scaling factors
        plotsize = self.params.plotsize if self.params else 500.0
        self.plotscale = HEC_TO_M2 / plotsize
        self.plotadj = self.plotscale / len(site.plots) if site.plots else self.plotscale
        self.plotrenorm = 1.0 / plotsize / len(site.plots) if site.plots else 1.0 / plotsize

        iwmo = site.site_id

        # Initialize arrays for plot-level data
        num_plots = len(site.plots)
        num_genera = species_pres.numgenera

        basal_area = np.zeros((num_plots, num_genera))
        biomC = np.zeros((num_plots, num_genera))
        biomN = np.zeros((num_plots, num_genera))
        leaf_bm = np.zeros((num_plots, num_genera))
        max_ht = np.zeros((num_plots, num_genera))
        max_diam = np.zeros((num_plots, num_genera))
        diam_categories = np.zeros((num_plots, num_genera, NHC), dtype=int)

        # Track which genera are present in site
        in_site = [False] * num_genera
        spec_index = [0] * num_genera

        # Check which genera are present and get their indices
        for is_gen in range(num_genera):
            for ip in range(num_plots):
                for ns, species in enumerate(site.plots[ip].species):
                    if species_pres.genusgroups[is_gen] == species.genus_name:
                        spec_index[is_gen] = ns
                        in_site[is_gen] = True
                        break
                if in_site[is_gen]:
                    break

        # Sum data over species/genera for each plot
        for ip in range(num_plots):
            plot_data = site.plots[ip].sum_over_sg(species_pres.genusgroups, 'genus')
            basal_area[ip, :] = plot_data['basal_area']
            leaf_bm[ip, :] = plot_data['leaf_bm']
            biomC[ip, :] = plot_data['biomC']
            biomN[ip, :] = plot_data['biomN']
            max_ht[ip, :] = plot_data['max_ht']
            max_diam[ip, :] = plot_data['max_diam']

            diam_cats = site.plots[ip].tree_dm_cats(species_pres.genusgroups, 'genus')
            diam_categories[ip, :, :] = diam_cats

        # Calculate totals and statistics
        total_biomc = np.zeros(num_genera)
        total_biomn = np.zeros(num_genera)
        total_basal = np.zeros(num_genera)
        total_laibm = np.zeros(num_genera)
        total_max_ht = np.zeros(num_genera)
        total_max_diam = np.zeros(num_genera)
        plot_mean_biomc = np.zeros(num_genera)
        plot_std_biomc = np.zeros(num_genera)
        plot_mean_biomn = np.zeros(num_genera)
        plot_std_biomn = np.zeros(num_genera)
        total_dm_cats = np.zeros((num_genera, NHC), dtype=int)

        for is_gen in range(num_genera):
            if in_site[is_gen]:
                ns = spec_index[is_gen]
                for ip in range(num_plots):
                    total_biomc[is_gen] += biomC[ip, is_gen]
                    total_biomn[is_gen] += biomN[ip, is_gen]
                    total_basal[is_gen] += basal_area[ip, is_gen]

                    # Calculate LAI biomass
                    if ns < len(site.plots[ip].species):
                        leafarea_c = site.plots[ip].species[ns].leafarea_c
                        total_laibm[is_gen] += leaf_bm[ip, is_gen] / (leafarea_c * 2.0)

                    total_max_ht[is_gen] = max(total_max_ht[is_gen], max_ht[ip, is_gen])
                    total_max_diam[is_gen] = max(total_max_diam[is_gen], max_diam[ip, is_gen])

                    total_dm_cats[is_gen, :] += diam_categories[ip, is_gen, :]

            # Calculate statistics
            plot_mean_biomc[is_gen], plot_std_biomc[is_gen] = stddev(biomC[:, is_gen].tolist(), RNVALID)
            plot_mean_biomn[is_gen], plot_std_biomn[is_gen] = stddev(biomN[:, is_gen].tolist(), RNVALID)

        # Apply scaling factors
        plot_mean_biomc *= self.plotscale
        plot_std_biomc *= self.plotscale
        plot_mean_biomn *= self.plotscale
        plot_std_biomn *= self.plotscale

        # Write data for each genus
        for is_gen in range(num_genera):
            row_data = [iwmo, year, species_pres.genusgroups[is_gen]]

            if in_site[is_gen]:
                # Apply scaling
                total_basal[is_gen] *= self.plotrenorm
                total_laibm[is_gen] *= self.plotrenorm
                total_biomc[is_gen] *= self.plotadj
                total_biomn[is_gen] *= self.plotrenorm * 10.0
                total_dm_cats[is_gen, :] = (total_dm_cats[is_gen, :] * self.plotadj).astype(int)

                row_data.extend([
                    total_dm_cats[is_gen, 0], total_dm_cats[is_gen, 1],
                    total_dm_cats[is_gen, 2], total_dm_cats[is_gen, 3],
                    total_dm_cats[is_gen, 4], total_dm_cats[is_gen, 5],
                    total_dm_cats[is_gen, 6],
                    total_max_diam[is_gen], total_max_ht[is_gen],
                    total_laibm[is_gen], total_basal[is_gen],
                    total_biomc[is_gen], plot_std_biomc[is_gen],
                    total_biomn[is_gen], plot_std_biomn[is_gen]
                ])
            else:
                # Fill with invalid values
                row_data.extend([RNVALID] * 14)

            csv_write(self.file_handles['biom_by_g'], row_data)

    def write_species_data(self, site: SiteData, species_pres: Groups, year: int):
        """
        Write species-level biomass and demographic data.

        Args:
            site: Site data
            species_pres: Groups containing species information
            year: Current year
        """
        # Calculate plot scaling factors (same as in write_genus_data)
        plotsize = self.params.plotsize if self.params else 500.0
        self.plotscale = HEC_TO_M2 / plotsize
        self.plotadj = self.plotscale / len(site.plots) if site.plots else self.plotscale
        self.plotrenorm = 1.0 / plotsize / len(site.plots) if site.plots else 1.0 / plotsize

        iwmo = site.site_id
        num_plots = len(site.plots)
        num_species = species_pres.numspecies

        # Arrays to hold per-plot data
        basal_area = np.zeros((num_plots, num_species))
        biomC = np.zeros((num_plots, num_species))
        biomN = np.zeros((num_plots, num_species))
        leaf_bm = np.zeros((num_plots, num_species))
        max_ht = np.zeros((num_plots, num_species))
        max_diam = np.zeros((num_plots, num_species))
        diam_categories = np.zeros((num_plots, num_species, NHC), dtype=int)

        # Track which species are present in site
        in_site = [False] * num_species
        spec_index = [0] * num_species

        # Check which species are present and get their indices
        for is_sp in range(num_species):
            genus_name, species_id = species_pres.spec_names[is_sp]
            for ip in range(num_plots):
                for ns, species in enumerate(site.plots[ip].species):
                    if species.unique_id == species_id:
                        spec_index[is_sp] = ns
                        in_site[is_sp] = True
                        break
                if in_site[is_sp]:
                    break

        # Sum data over species for each plot
        for ip in range(num_plots):
            plot_data = site.plots[ip].sum_over_sg(species_pres.spec_names, 'species')
            basal_area[ip, :] = plot_data['basal_area']
            leaf_bm[ip, :] = plot_data['leaf_bm']
            biomC[ip, :] = plot_data['biomC']
            biomN[ip, :] = plot_data['biomN']
            max_ht[ip, :] = plot_data['max_ht']
            max_diam[ip, :] = plot_data['max_diam']

            diam_cats = site.plots[ip].tree_dm_cats(species_pres.spec_names, 'species')
            diam_categories[ip, :, :] = diam_cats

        # Calculate totals and statistics
        total_biomc = np.zeros(num_species)
        total_biomn = np.zeros(num_species)
        total_basal = np.zeros(num_species)
        total_laibm = np.zeros(num_species)
        total_max_ht = np.zeros(num_species)
        total_max_diam = np.zeros(num_species)
        plot_mean_biomc = np.zeros(num_species)
        plot_std_biomc = np.zeros(num_species)
        plot_mean_biomn = np.zeros(num_species)
        plot_std_biomn = np.zeros(num_species)
        total_dm_cats = np.zeros((num_species, NHC), dtype=int)

        for is_sp in range(num_species):
            if in_site[is_sp]:
                ns = spec_index[is_sp]
                for ip in range(num_plots):
                    total_biomc[is_sp] += biomC[ip, is_sp]
                    total_biomn[is_sp] += biomN[ip, is_sp]
                    total_basal[is_sp] += basal_area[ip, is_sp]

                    # Calculate LAI biomass
                    if ns < len(site.plots[ip].species):
                        leafarea_c = site.plots[ip].species[ns].leafarea_c
                        total_laibm[is_sp] += leaf_bm[ip, is_sp] / (leafarea_c * 2.0)

                    total_max_ht[is_sp] = max(total_max_ht[is_sp], max_ht[ip, is_sp])
                    total_max_diam[is_sp] = max(total_max_diam[is_sp], max_diam[ip, is_sp])

                    total_dm_cats[is_sp, :] += diam_categories[ip, is_sp, :]

            # Calculate statistics
            plot_mean_biomc[is_sp], plot_std_biomc[is_sp] = stddev(biomC[:, is_sp].tolist(), RNVALID)
            plot_mean_biomn[is_sp], plot_std_biomn[is_sp] = stddev(biomN[:, is_sp].tolist(), RNVALID)

        # Apply scaling factors
        plot_mean_biomc *= self.plotscale
        plot_std_biomc *= self.plotscale
        plot_mean_biomn *= self.plotscale
        plot_std_biomn *= self.plotscale

        # Write data for each species
        for is_sp in range(num_species):
            genus_name, species_id = species_pres.spec_names[is_sp]
            row_data = [iwmo, year, genus_name, species_id]

            if in_site[is_sp]:
                # Apply scaling
                total_basal[is_sp] *= self.plotrenorm
                total_laibm[is_sp] *= self.plotrenorm
                total_biomc[is_sp] *= self.plotadj
                total_biomn[is_sp] *= self.plotrenorm * 10.0
                total_dm_cats[is_sp, :] = (total_dm_cats[is_sp, :] * self.plotadj).astype(int)

                row_data.extend([
                    total_dm_cats[is_sp, 0], total_dm_cats[is_sp, 1],
                    total_dm_cats[is_sp, 2], total_dm_cats[is_sp, 3],
                    total_dm_cats[is_sp, 4], total_dm_cats[is_sp, 5],
                    total_dm_cats[is_sp, 6],
                    total_max_diam[is_sp], total_max_ht[is_sp],
                    total_laibm[is_sp], total_basal[is_sp],
                    total_biomc[is_sp], plot_std_biomc[is_sp],
                    total_biomn[is_sp], plot_std_biomn[is_sp]
                ])
            else:
                # Fill with invalid values
                row_data.extend([RNVALID] * 14)

            csv_write(self.file_handles['biom_by_s'], row_data)

    def write_site_data(self, site: SiteData, year: int):
        """
        Write site-level climate and environmental data.

        Args:
            site: Site data
            year: Current year
        """
        row_data = [site.site_id, year]
        row_data.extend(self.extract_site_csv_data(site))
        csv_write(self.file_handles['clim_unit'], row_data)

    def extract_site_csv_data(self, site: SiteData) -> List[float]:
        """Extract site data for CSV output."""
        return [
            site.latitude, site.longitude, site.elevation, site.slope,
            site.leaf_area_ind, site.grow_days, site.deg_days,
            site.flood_days, site.dry_days_upper_layer, site.dry_days_base_layer,
            site.pot_evap_day, site.act_evap_day, site.rain
        ]

    def write_soil_data(self, site: SiteData, year: int):
        """
        Write soil carbon and nitrogen data.

        Args:
            site: Site data
            year: Current year
        """
        row_data = [site.site_id, year]
        row_data.extend(self.extract_soil_csv_data(site.soil))
        csv_write(self.file_handles['c_and_n'], row_data)

    def extract_soil_csv_data(self, soil) -> List[float]:
        """Extract soil data for CSV output matching Fortran Soil.f90 write_soil_csv."""
        return [
            soil.A0_c0, soil.A_c0, soil.A0_n0, soil.A_n0, soil.BL_c0, soil.BL_n0,
            soil.total_C_rsp, soil.biomC, soil.C_into_A0, soil.net_C_into_A0, soil.net_prim_prodC,
            soil.biomN, soil.N_into_A0, soil.net_N_into_A0, soil.net_prim_prodN, soil.avail_N
        ]

    def write_tree_data(self, site: SiteData, year: int):
        """
        Write individual tree data.

        Args:
            site: Site data
            year: Current year
        """
        iwmo = site.site_id

        for ip, plot in enumerate(site.plots):
            for it, tree in enumerate(plot.trees):
                row_data = [iwmo, year, ip + 1, it + 1]  # 1-based indexing
                row_data.extend(self.extract_tree_csv_data(tree))
                csv_write(self.file_handles['tld'], row_data)

    def extract_tree_csv_data(self, tree) -> List[float]:
        """Extract tree data for CSV output."""
        return [
            tree.species_id, tree.dbh, tree.height,
            tree.biomass_c, tree.biomass_n, tree.leaf_biomass,
            tree.age, tree.crown_area
        ]

    def close_output_files(self):
        """Close all output files."""
        for file_handle in self.file_handles.values():
            file_handle.close()
        self.file_handles.clear()
        self.csv_writers.clear()


# Global output manager instance
output_manager = OutputManager()


# Convenience functions for compatibility
def initialize_output_files(species_present: Groups):
    """Initialize output files and write headers."""
    output_manager.initialize_output_files(species_present)


def write_genus_data(site: SiteData, species_present: Groups, year: int):
    """Write genus data."""
    output_manager.write_genus_data(site, species_present, year)


def write_species_data(site: SiteData, species_present: Groups, year: int):
    """Write species data."""
    output_manager.write_species_data(site, species_present, year)


def write_site_data(site: SiteData, year: int):
    """Write site data."""
    output_manager.write_site_data(site, year)


def write_soil_data(site: SiteData, year: int):
    """Write soil data."""
    output_manager.write_soil_data(site, year)


def write_tree_data(site: SiteData, year: int):
    """Write tree data."""
    output_manager.write_tree_data(site, year)


def close_output_files():
    """Close output files."""
    output_manager.close_output_files()