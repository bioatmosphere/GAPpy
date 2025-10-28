"""
UVAFME (UVA Forest Model Enhanced) main module.
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


class UVAFMEModel:
    """Main UVAFME model class."""
    
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
        for site in self.sites:
            # Attach species to site
            site.attach_species(self._species_data)
            
            # Initialize plots for each site
            site.numplots = self.parameters.numplots
            site.plots = []
            
            for plot_id in range(self.parameters.numplots):
                plot = PlotData()
                plot.initialize_plot(site.species, self.parameters.maxtrees, self.parameters.maxheight)
                site.plots.append(plot)

            # NO initial trees - forest starts empty like Fortran
            # Trees will establish through renewal/recruitment process
            # self.forest_model.initialize_forest(site)  # REMOVED - not in Fortran
    
    def initialize_output_files(self, species_present):
        """Initialize output files."""
        self.output_manager.initialize_output_files(species_present)
    
    def draw_banner(self, numsites, species_present):
        """Draw the startup banner."""
        print("=" * 80)
        print("                       UVA Forest Model Enhanced")
        print("               Center For Regional Environmental Studies")
        print("                        University of Virginia")
        print("                   Department of Environmental Sciences")
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
        
    def bio_geo_climate(self, site, year):
        """Biogeochemical climate processing."""
        self.forest_model.bio_geo_climate(site, year)
    
    def write_soil_data(self, site, year):
        """Write soil data."""
        self.output_manager.write_soil_data(site, year)
    
    def write_site_data(self, site, year):
        """Write site data."""
        self.output_manager.write_site_data(site, year)
    
    def canopy(self, site):
        """Process canopy dynamics."""
        self.forest_model.canopy(site)
    
    def growth(self, site):
        """Process tree growth."""
        self.forest_model.growth(site)
    
    def mortality(self, site):
        """Process tree mortality."""
        self.forest_model.mortality(site)
    
    def renewal(self, site):
        """Process tree renewal."""
        self.forest_model.renewal(site)
    
    def write_genus_data(self, site, species_present, year):
        """Write genus data."""
        # Group species by genus
        # Get species data from the stored list
        species_data = self._species_data if hasattr(self, '_species_data') else []
        genera = list(set(species.genus_name for species in species_data))
        
        for plot in site.plots:
            genus_data = plot.sum_over_sg(genera, field='genus')
            
            # Write genus-level data (implementation depends on output format)
            # This is a placeholder for genus-specific output
            pass
    
    def write_species_data(self, site, species_present, year):
        """Write species data."""
        # Pass the Groups object which has the required numspecies and spec_names attributes
        self.output_manager.write_species_data(site, self.species_present, year)
    
    def write_tree_data(self, site, year):
        """Write tree data."""
        self.output_manager.write_tree_data(site, year)
    
    def close_output_files(self):
        """Close output files."""
        self.output_manager.close_output_files()
    
    def set_site_rng_seed(self, fixed_seed):
        """Set random number generator seed."""
        self.forest_model.set_site_rng_seed(fixed_seed)
    
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
        
        # Main site loop
        for sndx, site in enumerate(self.sites):
            
            # Check if site exists
            if site.site_wmo == self.parameters.rnvalid:
                print(f"             No site or climate file for site {site.site_id}")
                print(f"             Skipping site {site.site_name}")
                print()
                continue
            
            # Load climate and site specific vars, then adjust for altitude
            self.set_site_rng_seed(self.parameters.fixed_seed)
            set_site_climate(self.parameters.same_climate, self.parameters.fixed_seed)
            
            # Skip this site if no species are present
            if len(site.species) == 0:
                print(f"              No species present in site {site.site_id}")
                print(f"              Skipping site {site.site_name}")
                print()
                continue
            
            self.show_progress(site)
            
            # Run the model for this site
            for year in range(self.parameters.numyears + 1):
                
                # Run annual cycle
                self.forest_model.run_annual_cycle(site, year)
                
                # Write soil/CN/clim data for last year's trees but after BioGeo
                self.write_soil_data(site, year)
                self.write_site_data(site, year)
                
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
                    self.write_genus_data(site, self.species_present, year)
                    self.write_species_data(site, self.species_present, year)
                    if self.parameters.tree_level_data:
                        self.write_tree_data(site, year)
                
                # Progress indicator
                if year % 10 == 0:
                    if site.plots:
                        total_live = sum(len(plot.get_live_trees()) for plot in site.plots)
                        total_dead = sum(len(plot.get_dead_trees()) for plot in site.plots)
                        stats = site.plots[0].get_statistics() if site.plots else {}
                        print(f"  Year {year}: {total_live} live, {total_dead} dead, "
                              f"BiomC={stats.get('total_biomass_c', 0):.6f}")
                    else:
                        print(f"  Year {year}: No plots")
            
            # Report timing
            total_time = time.time()
            print(f"Cumulative time: {total_time - start_time:.2f} seconds")
            print("=" * 80)
        
        # Close output files
        self.close_output_files()


def main():
    """Main entry point for UVAFME model."""
    model = UVAFMEModel()
    model.run()


if __name__ == "__main__":
    main()