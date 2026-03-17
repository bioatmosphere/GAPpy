"""
GAPpy (Gap Model - Python) main module.
Translated from UVAFME.f90
"""

import sys
import time
import os
import numpy as np
from .constants import *
from .parameters import Parameters
from .species import SpeciesData
from .site import SiteData
from .plot import PlotData
from .tree import TreeData
from .model import ForestModel
from .genus_groups import Groups, initialize_genus_groups
from .input_module import InputFileManager
from .output_module import OutputManager
from .sitelist import initialize_sitelist
from .climate import set_site_climate
from .dispersal import compute_distance_matrix, build_species_index_map, disperse_seeds


class GAPpyModel:
    """Main GAPpy model class."""
    
    def __init__(self):
        self.sites = []
        self.species_present = Groups()
        self.filelist = ""
        self.input_manager = InputFileManager()
        self.parameters = Parameters()
        self.output_manager = OutputManager(self.parameters)
        self.forest_model = ForestModel()
        
    def initialize_input_files(self, filelist=""):
        """Initialize input files and read data."""
        self.filelist = filelist
        
        # Use actual UVAFME input file names from input_data directory
        # Support multiple possible locations for input_data directory
        possible_paths = ['input_data', '../input_data', '../../input_data']
        input_base_path = None

        for path in possible_paths:
            if os.path.exists(path):
                input_base_path = path
                break

        if input_base_path is None:
            raise FileNotFoundError("Could not find input_data directory in any expected location")

        config_file = os.path.join(input_base_path, 'uvafme_config.json')
        species_file = os.path.join(input_base_path, 'UVAFME2012_specieslist.csv')
        site_file = os.path.join(input_base_path, 'UVAFME2012_site.csv')
        climate_file = os.path.join(input_base_path, 'UVAFME2012_climate.csv')
        climate_std_file = os.path.join(input_base_path, 'UVAFME2012_climate_stddev.csv')
        sitelist_file = os.path.join(input_base_path, 'UVAFME2012_sitelist.txt')
        
        # Load parameters
        try:
            self.parameters.load_from_file(config_file)
        except Exception as e:
            print(f"Using default parameters: {e}")
        
        # Setup input manager with file paths
        self.input_manager.filenames.update({
            'species': species_file,
            'sites': site_file,
            'climate': climate_file,
            'climate_std': climate_std_file,
            'sitelist': sitelist_file
        })

        # Read input data using input manager
        try:
            species_data = self.input_manager.read_species_data()
            print(f"Loaded {len(species_data)} species")

            # Create a Groups object from the species data
            self.species_present = Groups()
            self.species_present.numspecies = len(species_data)

            # Populate spec_names with (genus_name, unique_id) pairs
            self.species_present.spec_names = [(species.genus_name, species.unique_id) for species in species_data]

            # Extract unique genera
            unique_genera = list(set(species.genus_name for species in species_data))
            self.species_present.genusgroups = unique_genera
            self.species_present.numgenera = len(unique_genera)

            # Store the actual species data for later use
            self._species_data = species_data

            # Read site IDs and site data
            site_ids = [0]  # Use site ID 0 from the UVAFME data
            self.sites = self.input_manager.read_sites(site_ids)
            print(f"Loaded {len(self.sites)} sites")

            # Read and attach climate data
            self.input_manager.read_climate(self.sites)
            print("Climate data attached to sites")

            # Read and attach climate standard deviations
            self.input_manager.read_climate_stds(self.sites)
            print("Climate variability data attached to sites")

        except Exception as e:
            print(f"Warning: Error reading input files: {e}")
            # Initialize with empty data
            self.species_present = Groups()
            self.species_present.numspecies = 0
            self.species_present.spec_names = []
            self.species_present.genusgroups = []
            self.species_present.numgenera = 0
            self._species_data = []
            self.sites = []
    
    def initialize_sitelist(self):
        """Initialize site list and species data."""
        from .sitelist import apply_site_adjustments
        for site in self.sites:
            # Apply site adjustments (Fortran Sitelist.f90:33-61)
            # Must happen BEFORE attach_species so soil params are scaled
            apply_site_adjustments(site, self.parameters)

            # Attach species to site
            site.attach_species(self._species_data)

            # Initialize plots for each site
            site.numplots = self.parameters.numplots
            site.plots = []

            for plot_id in range(self.parameters.numplots):
                plot = PlotData()
                plot.initialize_plot(site.species, self.parameters.maxtrees, self.parameters.maxheight)
                site.plots.append(plot)

            # No initial trees - forest starts empty like Fortran
            # Trees establish through the renewal/recruitment process
    
    def initialize_output_files(self, species_present):
        """Initialize output files."""
        self.output_manager.initialize_output_files(species_present)
    
    def draw_banner(self, numsites, species_present):
        """Draw the startup banner."""
        print("=" * 80)
        print("                     GAPpy - Gap Model in Python")
        print("=" * 80)
        print()
        print("Running with parameters:")
        print(f"Number of sites: {numsites}")
        print(f"Number of years: {self.parameters.numyears}")
        print(f"Number of species: {species_present.numspecies if hasattr(species_present, 'numspecies') else len(species_present)}")
        print(f"Maximum number of trees: {self.parameters.maxtrees}")
        print(f"Maximum height of trees: {self.parameters.maxheight}")
        print(f"Plotsize: {self.parameters.plotsize}")
        print(f"Root depth: {self.parameters.rootdepth}")

        if self.parameters.with_clim_change:
            print("Running with climate change")
            print(f"Beginning at year: {self.parameters.begin_change_year}")
            print(f"Duration in years: {self.parameters.duration_of_change}")

        print(f"Printing interval in years: {self.parameters.year_print_interval}")
        print("=" * 80)
        print()
    
    def show_progress(self, site):
        """Show progress for current site."""
        num_site_species = len(site.species)
        print(f"Running for site {site.site_id} {site.site_name}")
        print(f"Number of species present {num_site_species}")
        
        if hasattr(site, 'altitude') and site.altitude != self.parameters.rnvalid:
            print(f"      Site altitude adjustment {site.altitude}")
        
        print(f"Site parameters: elevation {site.elevation:.2f}   slope     {site.slope:.2f}")
        print(f"                 fire/1000   {site.fire_prob:.2f}   wind/1000 {site.wind_prob:.2f}")
        print(f"                 SAFC*rtdpth {site.soil.A_field_cap:.2f}   A0_C      {site.soil.A0_c0:.2f}  A0_N {site.soil.A0_n0:.2f}")
        
    def run(self, filelist=""):
        """Main model execution loop."""
        
        # Handle command line arguments
        if len(sys.argv) > 1:
            filelist = sys.argv[1]
        
        # Prepare input files
        self.initialize_input_files(filelist)
        
        # Prepare site and species data
        self.initialize_sitelist()
        
        # Prepare output files
        self.initialize_output_files(self.species_present)
        
        # Write runtime vars to screen
        self.draw_banner(len(self.sites), self.species_present)
        
        # Start timing
        start_time = time.time()

        # Filter out invalid sites and prepare active site list
        active_sites = []
        for site in self.sites:
            if site.site_wmo == self.parameters.rnvalid:
                print(f"             No site or climate file for site {site.site_id}")
                print(f"             Skipping site {site.site_name}")
                print()
                continue
            if len(site.species) == 0:
                print(f"              No species present in site {site.site_id}")
                print(f"              Skipping site {site.site_name}")
                print()
                continue
            active_sites.append(site)

        # Initialize climate and RNG for each active site
        for site in active_sites:
            self.forest_model.set_site_rng_seed(self.parameters.fixed_seed)
            set_site_climate(self.parameters.same_climate, self.parameters.fixed_seed)
            self.show_progress(site)

        # Pre-compute inter-site dispersal structures if enabled
        use_dispersal = (self.parameters.seed_dispersal and len(active_sites) > 1)
        if use_dispersal:
            distance_matrix = compute_distance_matrix(active_sites)
            species_index_maps = build_species_index_map(active_sites)
            print(f"Seed dispersal enabled between {len(active_sites)} sites")

        # Main simulation loop: all sites advance in lockstep
        for year in range(self.parameters.numyears + 1):
            for site in active_sites:

                # Biogeochemical processes first
                self.forest_model.bio_geo_climate(site, year)

                # Write soil/CN/clim data after BioGeo but before tree dynamics
                self.output_manager.write_soil_data(site, year)
                self.output_manager.write_site_data(site, year)

                # Tree dynamics
                self.forest_model.canopy(site, year)
                self.forest_model.growth(site, year)
                self.forest_model.mortality(site)
                self.forest_model.renewal(site)
                self.forest_model.update_site_statistics(site)

                # Determine print interval
                if self.parameters.spinup:
                    if year < self.parameters.spinup_yrs:
                        print_interval = 10 * self.parameters.year_print_interval
                    else:
                        print_interval = self.parameters.year_print_interval
                else:
                    print_interval = self.parameters.year_print_interval

                # Write output data
                if (year % print_interval == 0) or (year == self.parameters.numyears):
                    self.output_manager.write_genus_data(site, self.species_present, year)
                    self.output_manager.write_species_data(site, self.species_present, year)
                    if self.parameters.tree_level_data:
                        self.output_manager.write_tree_data(site, year)

                # Progress indicator
                if year % 10 == 0:
                    if site.plots:
                        total_live = sum(len(plot.get_live_trees()) for plot in site.plots)
                        total_dead = sum(len(plot.get_dead_trees()) for plot in site.plots)
                        stats = site.plots[0].get_statistics() if site.plots else {}
                        print(f"  Site {site.site_id} Year {year}: {total_live} live, "
                              f"{total_dead} dead, "
                              f"BiomC={stats.get('total_biomass_c', 0):.6f}")
                    else:
                        print(f"  Site {site.site_id} Year {year}: No plots")

            # Inter-site seed dispersal after all sites complete this year
            if use_dispersal:
                disperse_seeds(active_sites, distance_matrix, species_index_maps)

            # Periodic timing update
            if year % 100 == 0 and year > 0:
                elapsed = time.time() - start_time
                print(f"Year {year} complete. Elapsed: {elapsed:.2f}s")

        # Report final timing
        total_time = time.time() - start_time
        print(f"Total simulation time: {total_time:.2f} seconds")
        print("=" * 80)

        # Close output files
        self.output_manager.close_output_files()


def main():
    """Main entry point for GAPpy model."""
    model = GAPpyModel()
    model.run()


if __name__ == "__main__":
    main()