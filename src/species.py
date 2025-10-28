"""
Species module for UVAFME vegetation model.
Translated from Species.f90
"""

import math
import numpy as np
from .constants import *

class SpeciesData:
    """Species data structure containing attributes for all members of a species."""
    
    def __init__(self):
        # String attributes
        self.genus_name = ""
        self.taxonomic_name = ""
        self.unique_id = ""
        self.common_name = ""
        
        # Integer attributes
        self.genus_id = 0
        self.species_id = 0
        self.shade_tol = 0
        self.lownutr_tol = 0
        self.stress_tol = 0
        self.age_tol = 0
        self.drought_tol = 0
        self.flood_tol = 0
        self.fire_tol = 0
        
        # Real attributes
        self.max_age = 0.0
        self.max_diam = 0.0
        self.max_ht = 0.0
        self.wood_bulk_dens = 0.0
        self.rootdepth = 0.0
        self.leafdiam_a = 0.0
        self.leafarea_c = 0.0
        self.deg_day_min = 0.0
        self.deg_day_opt = 0.0
        self.deg_day_max = 0.0
        self.seed_surv = 0.0
        self.seedling_lg = 0.0
        self.invader = 0.0
        self.seed_num = 0.0
        self.sprout_num = 0.0
        self.arfa_0 = 0.0
        self.g = 0.0
        self.fc_fire = 0.0
        self.fc_wind = 0.0
        self.fc_degday = 0.0
        self.fc_drought = 0.0
        self.fc_flood = 0.0
        
        # Boolean attributes
        self.conifer = False
    
    def initialize_species(self, species_id, genus_name, taxonomic_name, unique_id, 
                          common_name, genus_id, shade_tol, lownutr_tol, stress_tol,
                          age_tol, drought_tol, flood_tol, fire_tol, max_age, max_diam,
                          max_ht, wood_bulk_dens, rootdepth, leafdiam_a, leafarea_c,
                          deg_day_min, deg_day_opt, deg_day_max, seedling_lg, invader,
                          seed_num, sprout_num, seed_surv, arfa_0, g, conifer):
        """Initialize species with all parameters."""
        
        ss = [1.1, 1.15, 1.2, 1.23, 1.25]
        adjust = [1.5, 1.55, 1.6, 1.65, 1.7]
        
        self.species_id = species_id
        self.genus_name = genus_name
        self.taxonomic_name = taxonomic_name
        self.common_name = common_name
        self.unique_id = unique_id
        self.genus_id = genus_id
        self.max_age = max_age
        self.max_diam = max_diam
        self.rootdepth = rootdepth
        self.wood_bulk_dens = wood_bulk_dens
        self.deg_day_min = deg_day_min
        self.deg_day_max = deg_day_max
        self.deg_day_opt = deg_day_opt
        self.shade_tol = shade_tol
        self.lownutr_tol = lownutr_tol
        self.drought_tol = drought_tol
        self.fire_tol = fire_tol
        self.flood_tol = flood_tol
        self.stress_tol = stress_tol
        self.age_tol = age_tol
        self.conifer = conifer
        self.invader = invader
        self.seed_num = seed_num
        self.sprout_num = sprout_num
        self.seed_surv = seed_surv
        self.seedling_lg = seedling_lg
        self.arfa_0 = arfa_0
        
        # Adjustments
        self.leafarea_c = leafarea_c / HEC_TO_M2
        self.max_ht = min(max_ht, rootdepth * 80.0 / (1 + rootdepth))
        self.leafdiam_a = leafdiam_a * adjust[shade_tol - 1]  # Adjust for 0-based indexing
        self.g = g * ss[shade_tol - 1]  # Adjust for 0-based indexing
        
        # Initialized to constants
        self.fc_fire = 0.0
        self.fc_wind = 0.0
        self.fc_degday = 0.0
        self.fc_drought = 0.0
        self.fc_flood = 0.0
    
    def light_rsp(self, al):
        """Compute available sunlight factor by tolerance class."""
        light_c1 = [1.01, 1.04, 1.11, 1.24, 1.49]
        light_c2 = [4.62, 3.44, 2.52, 1.78, 1.23]
        light_c3 = [0.05, 0.06, 0.07, 0.08, 0.09]
        
        kt = self.shade_tol - 1  # Adjust for 0-based indexing
        
        flight = light_c1[kt] * (1.0 - math.exp(-light_c2[kt] * (al - light_c3[kt])))
        flight = max(0.0, min(1.0, flight))
        
        return flight
    
    def temp_rsp(self, x):
        """Compute temperature response factor."""
        ddmin = self.deg_day_min
        ddmax = self.deg_day_max
        ddopt = self.deg_day_opt
        
        a = (ddopt - ddmin) / (ddmax - ddmin)
        b = (ddmax - ddopt) / (ddmax - ddmin)
        
        if x >= ddmax or x <= ddmin:
            ftemp = 0.0
        else:
            tmp = ((x - ddmin) / (ddopt - ddmin)) ** a
            ftemp = tmp * ((ddmax - x) / (ddmax - ddopt)) ** b
        
        self.fc_degday = ftemp
    
    def drought_rsp(self, drydays, drydays_s):
        """Compute drought response factor."""
        if self.drought_tol == 1:
            if self.conifer:
                fcdry1 = fdry(drydays_s, 1) * 0.33
            else:
                fcdry1 = fdry(drydays_s, 1) * 0.2
            
            fcdry2 = fdry(drydays, 1)
            self.fc_drought = max(fcdry1, fcdry2)
        else:
            self.fc_drought = fdry(drydays, self.drought_tol)
    
    def flood_rsp(self, floodday):
        """Compute flood response factor."""
        gama = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5]
        k = self.flood_tol - 1  # Adjust for 0-based indexing
        
        if floodday <= gama[k]:
            fflood = 1.0
        else:
            fflood = 1.0  # Original function always returns 1
        
        self.fc_flood = fflood
    
    def fire_rsp(self, fire):
        """Compute fire response factor."""
        gama = [100.0, 10.0, 1.0, 0.1, 0.01, 0.001]
        k = self.fire_tol - 1  # Adjust for 0-based indexing
        
        if fire == 1:
            resp = gama[k]
        else:
            resp = 1.0
        
        self.fc_fire = resp
    
    def poor_soil_rsp(self, n_avail):
        """Calculate poor soil nutrient response."""
        # Based on the original poor_soil_rsp function but adapted for class method
        nrc = self.lownutr_tol
        sf = max(0.0, min(1.0, n_avail))
        
        fert_c1 = [-0.6274, -0.2352, 0.2133]
        fert_c2 = [3.600, 2.771, 1.789]
        fert_c3 = [-1.994, -1.550, -1.014]
        
        # WARNING: This appears to be a bug in the original code
        nrc = 4 - nrc
        nrc = max(1, min(3, nrc))  # Clamp to valid range
        
        fpoor = fert_c1[nrc - 1] + fert_c2[nrc - 1] * sf + fert_c3[nrc - 1] * sf**2
        fpoor = max(0.0, min(1.0, fpoor))
        
        return fpoor * sf


def poor_soil_rsp(sf, nrc):
    """Compute quadratic nutrient response factors."""
    fert_c1 = [-0.6274, -0.2352, 0.2133]
    fert_c2 = [3.600, 2.771, 1.789]
    fert_c3 = [-1.994, -1.550, -1.014]
    
    # WARNING: This appears to be a bug in the original code
    nrc = 4 - nrc
    sf = min(sf, 1.0)
    
    fpoor = fert_c1[nrc - 1] + fert_c2[nrc - 1] * sf + fert_c3[nrc - 1] * sf**2
    fpoor = max(0.0, min(1.0, fpoor))
    
    return fpoor * sf


def fdry(dryday, k):
    """Compute drought response factor."""
    gama = [0.50, 0.45, 0.35, 0.25, 0.15, 0.05]
    
    tmp = max(gama[k - 1] - dryday, 0.0)  # Adjust for 0-based indexing
    return (tmp / gama[k - 1]) ** 0.5