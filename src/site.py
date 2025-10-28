"""
Site module for UVAFME vegetation model.
Translated from Site.f90
"""

import numpy as np
from .constants import *
from .species import SpeciesData
from .soil import SoilData

# Placeholder for missing functions
def kron(x):
    """Kronecker delta function placeholder"""
    return 1.0 if x > 0 else 0.0


class SiteData:
    """Site data structure containing site attributes and species/plot data."""
    
    def __init__(self):
        self.soil = SoilData()
        self.plots = []
        self.species = []
        self.fc_flood = []
        self.fc_drought = []
        self.fc_degday = []
        
        # Site attributes
        self.region = ""
        self.site_name = ""
        self.site_id = 0
        self.numplots = 0
        self.site_wmo = 0.0
        self.latitude = 0.0
        self.longitude = 0.0
        self.elevation = 0.0
        self.altitude = 0.0
        self.slope = 0.0
        self.leaf_area_ind = 0.0
        self.leaf_area_w0 = 0.0
        self.sigma = 0.0
        
        # Climate arrays
        self.temp_lapse_r = np.zeros(NTEMPS)
        self.precip_lapse_r = np.zeros(NTEMPS)
        self.tmin = np.zeros(NTEMPS)
        self.tmax = np.zeros(NTEMPS)
        self.precip = np.zeros(NTEMPS)
        self.tmin_std = np.zeros(NTEMPS)
        self.tmax_std = np.zeros(NTEMPS)
        self.precip_std = np.zeros(NTEMPS)
        
        # Daily climate variables
        self.rain = 0.0
        self.pot_evap_day = 0.0
        self.act_evap_day = 0.0
        self.grow_days = 0.0
        self.deg_days = 0.0
        self.flood_days = 0.0
        self.dry_days_upper_layer = 0.0
        self.dry_days_base_layer = 0.0
        self.freeze = 0.0
        self.fire_prob = 0.0
        self.wind_prob = 0.0

        # Daily climate arrays (365 days)
        self.daily_tmin = np.zeros(365)
        self.daily_tmax = np.zeros(365)
        self.daily_precip = np.zeros(365)
    
    def initialize_site(self, siteid, sitename, siteregion, lat, long, wmo,
                       elevation, slope, Afc, A_perm_wp, lai, base_h,
                       lai_w0, A0_w0, A_w0, sbase_w0, fire_prob,
                       wind_prob, A0_c0, A0_n0, A_c0, A_n0, sbase_c0,
                       sbase_n0, sigma, temp_lapse, prcp_lapse):
        """Initialize site with all parameters from site file."""
        
        # Initialize properties from the site file
        self.site_id = siteid
        self.site_name = sitename
        self.region = siteregion
        self.latitude = lat
        self.longitude = long
        self.site_wmo = wmo
        self.elevation = elevation
        self.slope = slope
        self.sigma = sigma
        self.temp_lapse_r = np.array(temp_lapse)
        self.precip_lapse_r = np.array(prcp_lapse)
        self.leaf_area_w0 = lai_w0
        self.leaf_area_ind = lai
        # Convert fire/wind probabilities from per-1000-years to annual (like Fortran)
        self.fire_prob = fire_prob / 1000.0
        self.wind_prob = wind_prob / 1000.0
        
        # Initialize soil properties
        self.soil.A_field_cap = Afc
        self.soil.A_perm_wp = A_perm_wp
        self.soil.base_h = base_h
        self.soil.A0_w0 = A0_w0
        self.soil.A_w0 = A_w0
        self.soil.BL_w0 = sbase_w0
        self.soil.A0_c0 = A0_c0
        self.soil.A0_n0 = A0_n0
        self.soil.A_c0 = A_c0
        self.soil.A_n0 = A_n0
        self.soil.BL_c0 = sbase_c0
        self.soil.BL_n0 = sbase_n0
        
        self.adjust_for_altitude()
    
    def attach_climate(self, tmin, tmax, prcp):
        """Attach climate data to site."""
        self.tmin = np.array(tmin)
        self.tmax = np.array(tmax)
        self.precip = np.array(prcp) * MM_TO_CM
    
    def attach_climate_std(self, tmin_std, tmax_std, prcp_std):
        """Attach climate standard deviation data to site."""
        self.tmin_std = np.array(tmin_std)
        self.tmax_std = np.array(tmax_std)
        self.precip_std = np.array(prcp_std) * MM_TO_CM
    
    def attach_species(self, species_data, range_species_ids=None):
        """Attach species list for the current site."""
        
        if range_species_ids is not None:
            # Filter species based on range list
            num_all_species = len(species_data)
            num_range_species = len(range_species_ids)
            num_site_species = sum(1 for sid in range_species_ids if sid != 'NP')
            
            if num_site_species == 0:
                self.species = []
            else:
                self.species = []
                for species in species_data:
                    for range_id in range_species_ids:
                        if species.unique_id == range_id:
                            self.species.append(species)
                            break
        else:
            # No range list, all species in this site
            self.species = list(species_data)
        
        # Initialize plots based on species
        # This would need to be implemented based on Plot module
        # For now, we'll just set numplots from parameters
        # self.numplots = numplots  # This would come from parameters
    
    def adjust_for_altitude(self):
        """Adjust temperature and precipitation for altitude."""
        # Placeholder for rnvalid - would need to be defined in parameters
        rnvalid = -9999.0  # Common invalid value
        
        if self.altitude != rnvalid:
            tav = 0.0
            for z in range(12):
                self.tmax[z] = self.tmax[z] - (self.altitude - self.elevation) * \
                              self.temp_lapse_r[z] * 0.01
                self.tmin[z] = self.tmin[z] - (self.altitude - self.elevation) * \
                              self.temp_lapse_r[z] * 0.01
                tav = tav + 0.5 * (self.tmax[z] + self.tmin[z])
                self.precip[z] = max(self.precip[z] + 
                                   (self.altitude - self.elevation) * 
                                   self.precip_lapse_r[z] * 0.001, 0.0)
            
            tav = tav / 12.0
            self.freeze = kron(-tav - 2.0)
    
    def compute_clim_stds(self, tmin_stds, tmax_stds, prcp_stds):
        """Compute climate standard deviations (dummy for now)."""
        # Placeholder implementation
        pass