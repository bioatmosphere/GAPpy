"""
Parameters module for UVAFME vegetation model.
Translated from Parameters.f90
"""

import json
import os
from typing import Dict, Any
from .constants import *


class Parameters:
    """Global parameters for UVAFME model."""
    
    def __init__(self):
        # Basic parameters
        self.numyears = 100
        self.numplots = 1
        self.maxtrees = 10000
        self.maxheight = 50
        
        # Values for invalid/missing data
        self.rnvalid = -999.0
        self.invalid = -999
        
        # Variables determining whether to use explicit seeds for RNGs
        self.fixed_seed = False
        self.same_climate = False
        self.debug = False
        
        # Spinup
        self.spinup = False
        self.spinup_yrs = 50
        
        # Climate change variables
        self.incr_tmin_by = 0.0
        self.incr_tmax_by = 0.0
        self.incr_precip_by = 0.0
        self.decr_tmin_by = 0.0
        self.decr_tmax_by = 0.0
        self.decr_precip_by = 0.0
        self.tmin_change = 0.0
        self.tmax_change = 0.0
        self.precip_change = 0.0
        self.begin_change_year = 0
        self.start_gcm = 0
        self.end_gcm = 0
        self.duration_of_change = 0
        self.incr_or_decr = "incr"
        self.year_print_interval = 10
        self.linear_cc = False
        self.with_clim_change = False
        self.use_gcm = False
        self.adjust_for_elev = False
        self.plot_level_data = False
        self.tree_level_data = False
        
        # Plot parameters
        self.plotsize = 500.0  # square meters (m²) - matches UVAFME2012_runtime.txt
        self.rootdepth = 2.0  # meters
        
        # Variables for changes to all site values
        self.new_slope = 0.0
        self.fire_level = 0.0
        self.wind_level = 0.0
        self.SA_field_cap = 0.0
        self.A0_level_C = 0.0
        self.A0_level_N = 0.0
        
        # Counters
        self.clim_counter = 0
        self.rand_counter = 0
    
    def load_from_file(self, filename: str):
        """Load parameters from a JSON configuration file."""
        if os.path.exists(filename):
            with open(filename, 'r') as f:
                config = json.load(f)
                self.load_from_dict(config)
        elif os.path.exists(os.path.join("input_data/", filename)):
            alt_path = os.path.join("input_data/", filename)
            with open(alt_path, 'r') as f:
                config = json.load(f)
                self.load_from_dict(config)
        else:
            print(f"Warning: Parameter file {filename} not found, using defaults")
    
    def load_from_dict(self, config: Dict[str, Any]):
        """Load parameters from a dictionary."""
        for key, value in config.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                print(f"Warning: Unknown parameter {key}")
    
    def save_to_file(self, filename: str):
        """Save parameters to a JSON configuration file."""
        config = self.to_dict()
        with open(filename, 'w') as f:
            json.dump(config, f, indent=2)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert parameters to dictionary."""
        return {
            'numyears': self.numyears,
            'numplots': self.numplots,
            'maxtrees': self.maxtrees,
            'maxheight': self.maxheight,
            'rnvalid': self.rnvalid,
            'invalid': self.invalid,
            'fixed_seed': self.fixed_seed,
            'same_climate': self.same_climate,
            'debug': self.debug,
            'spinup': self.spinup,
            'spinup_yrs': self.spinup_yrs,
            'incr_tmin_by': self.incr_tmin_by,
            'incr_tmax_by': self.incr_tmax_by,
            'incr_precip_by': self.incr_precip_by,
            'decr_tmin_by': self.decr_tmin_by,
            'decr_tmax_by': self.decr_tmax_by,
            'decr_precip_by': self.decr_precip_by,
            'tmin_change': self.tmin_change,
            'tmax_change': self.tmax_change,
            'precip_change': self.precip_change,
            'begin_change_year': self.begin_change_year,
            'start_gcm': self.start_gcm,
            'end_gcm': self.end_gcm,
            'duration_of_change': self.duration_of_change,
            'incr_or_decr': self.incr_or_decr,
            'year_print_interval': self.year_print_interval,
            'linear_cc': self.linear_cc,
            'with_clim_change': self.with_clim_change,
            'use_gcm': self.use_gcm,
            'adjust_for_elev': self.adjust_for_elev,
            'plot_level_data': self.plot_level_data,
            'tree_level_data': self.tree_level_data,
            'plotsize': self.plotsize,
            'rootdepth': self.rootdepth,
            'new_slope': self.new_slope,
            'fire_level': self.fire_level,
            'wind_level': self.wind_level,
            'SA_field_cap': self.SA_field_cap,
            'A0_level_C': self.A0_level_C,
            'A0_level_N': self.A0_level_N
        }
    
    def create_default_config(self, filename: str = "uvafme_config.json"):
        """Create a default configuration file."""
        default_config = {
            "numyears": 100,
            "numplots": 1,
            "maxtrees": 10000,
            "maxheight": 50,
            "fixed_seed": False,
            "same_climate": False,
            "debug": False,
            "spinup": False,
            "spinup_yrs": 50,
            "year_print_interval": 10,
            "with_clim_change": False,
            "tree_level_data": False,
            "plot_level_data": False,
            "plotsize": 1.0,
            "rootdepth": 2.0,
            "fire_level": 0.01,
            "wind_level": 0.05
        }
        
        with open(filename, 'w') as f:
            json.dump(default_config, f, indent=2)
        
        print(f"Default configuration saved to {filename}")
    
    def validate(self):
        """Validate parameter values."""
        if self.numyears <= 0:
            raise ValueError("numyears must be positive")
        if self.numplots <= 0:
            raise ValueError("numplots must be positive")
        if self.maxtrees <= 0:
            raise ValueError("maxtrees must be positive")
        if self.maxheight <= 0:
            raise ValueError("maxheight must be positive")
        if self.plotsize <= 0:
            raise ValueError("plotsize must be positive")
        if self.rootdepth <= 0:
            raise ValueError("rootdepth must be positive")
        if self.year_print_interval <= 0:
            raise ValueError("year_print_interval must be positive")
        
        return True


# Global parameters instance
params = Parameters()