"""
Tree module for UVAFME vegetation model.
Translated from Tree.f90
"""

import math
import numpy as np
from .constants import *
from .species import SpeciesData
from .random_utils import urand


# Constants global within this module
TC = 3.92699e-5  # Related to atomic mass of C
BETA = 1.0  # Shape parameter


class TreeData(SpeciesData):
    """Tree data structure extending SpeciesData with tree-specific attributes."""
    
    def __init__(self):
        super().__init__()
        
        # Tree-specific attributes
        self.diam_max = 0.0
        self.diam_bht = 0.0  # diameter at breast height
        self.diam_canht = 0.0  # diameter at canopy height
        self.canopy_ht = STD_HT
        self.foret_ht = 1.0
        self.forska_ht = 1.0
        self.leaf_bm = 0.0  # includes roots
        self.biomC = 0.0
        self.biomN = 0.0
        self.species_index = 0
        self.mort_marker = False
    
    def initialize_tree(self, tree_species=None, si=None):
        """Initialize tree with species data and index."""
        self.diam_bht = 0.0
        self.diam_canht = 0.0
        self.canopy_ht = STD_HT
        self.foret_ht = 1.0
        self.forska_ht = 1.0
        self.biomC = 0.0
        self.biomN = 0.0
        self.mort_marker = False
        
        if tree_species is not None:
            self.copy_species_data(tree_species)
        
        if si is not None:
            self.species_index = si
        else:
            self.species_index = 0
    
    def copy_species_data(self, species):
        """Copy species data to tree."""
        # Copy all species attributes
        for attr in dir(species):
            if not attr.startswith('_') and hasattr(self, attr):
                setattr(self, attr, getattr(species, attr))
        
        # Debug: Ensure environmental factors are copied
        if hasattr(species, 'fc_degday'):
            self.fc_degday = species.fc_degday
        if hasattr(species, 'fc_drought'):
            self.fc_drought = species.fc_drought  
        if hasattr(species, 'fc_flood'):
            self.fc_flood = species.fc_flood
    
    def copy_tree(self, tree):
        """Copy tree data from another tree."""
        self.diam_max = tree.diam_max
        self.diam_bht = tree.diam_bht
        self.diam_canht = tree.diam_canht
        self.canopy_ht = tree.canopy_ht
        self.foret_ht = tree.foret_ht
        self.forska_ht = tree.forska_ht
        self.leaf_bm = tree.leaf_bm
        self.biomC = tree.biomC
        self.biomN = tree.biomN
        self.species_index = tree.species_index
        self.mort_marker = tree.mort_marker
        self.copy_species_data(tree)
    
    def update_tree(self, tree_species):
        """Update tree with new species data."""
        self.copy_species_data(tree_species)
    
    def env_stress(self, shade):
        """Calculate environmental stress factor."""
        return self.fc_degday * self.fc_drought * self.fc_flood * shade
    
    def stem_shape(self):
        """Calculate stem shape - diameter at canopy height."""
        hc = self.canopy_ht
        h = self.forska_ht
        dbh = self.diam_bht
        
        if h <= hc or h <= STD_HT:
            dhshape = dbh
        else:
            dhshape = ((h - hc) / (h - STD_HT))**(1.0 / BETA) * dbh
        
        self.diam_canht = dhshape
    
    def forska_height(self):
        """Calculate height using Forska height-diameter relationship."""
        d = self.diam_bht
        hmax = self.max_ht
        par = self.arfa_0
        
        delta_ht = hmax - STD_HT
        self.forska_ht = STD_HT + delta_ht * (1.0 - math.exp(-(par * d / delta_ht)))
    
    def biomass_c(self):
        """Calculate total carbon biomass."""
        d = self.diam_bht
        h = self.forska_ht
        hc = self.canopy_ht
        hrt = self.rootdepth
        bulk = self.wood_bulk_dens
        
        # Create a copy for basal diameter calculation
        tree_0 = TreeData()
        tree_0.copy_tree(self)
        tree_0.canopy_ht = 0.0
        
        # Calculate stem shapes
        self.stem_shape()  # dc
        tree_0.stem_shape()  # bd
        
        # Calculate biomass components
        stembc = self.stem_biomass_c(tree_0)
        twigbc = self.twig_biomass_c()
        
        abovegr_c = stembc + twigbc
        root_c = stembc * hrt / h + twigbc / 2.0
        
        self.biomC = abovegr_c + root_c
    
    def biomass_n(self):
        """Calculate total nitrogen biomass."""
        self.biomN = self.biomC / STEM_C_N
    
    def leaf_biomass_c(self):
        """Calculate leaf biomass carbon (includes roots)."""
        self.leaf_bm = self.lai_biomass_c() * self.leafarea_c * 2.0
    
    def lai_biomass_c(self):
        """Calculate leaf area index biomass carbon."""
        dc = self.diam_canht  # diameter in cm
        leafd_a = self.leafdiam_a

        # Keep formula as-is matching Fortran
        # dc is in cm, result is used for LAI calculations
        return dc * dc * leafd_a
    
    def stem_biomass_c(self, tree_0):
        """Calculate stem biomass carbon."""
        bd = tree_0.diam_canht  # basal diameter
        h = tree_0.forska_ht
        bulk_density = tree_0.wood_bulk_dens
        
        yxd = TC * bulk_density * BETA / (BETA + 2.0)
        return yxd * bd * bd * h * 0.90
    
    def twig_biomass_c(self):
        """Calculate twig biomass carbon."""
        dc = self.diam_canht
        hc = self.canopy_ht
        h = self.forska_ht
        bulk_density = self.wood_bulk_dens
        
        yxd = TC * bulk_density * (2.0 / (BETA + 2.0) - 0.33)
        return yxd * dc * dc * (h - hc)
    
    def age_survival(self):
        """Check survival based on age tolerance."""
        check = [4.605, 6.908, 11.51]  # PER 1%, 0.1%, 0.001%
        
        k = max(0, min(len(check) - 1, self.age_tol - 1))  # Clamp to valid range
        agemax = self.max_age
        
        if agemax <= 0:
            return True  # Can't die of age if max_age is invalid
        
        rand_val = urand()
        mort_prob = check[k] / agemax
        
        # DEBUG: Print survival details for first few trees
        # if hasattr(self, '_debug_count') and self._debug_count < 3:
        #     print(f"DEBUG age_survival: age_tol={self.age_tol}, k={k}, agemax={agemax}, check[k]={check[k]}, mort_prob={mort_prob:.4f}, rand={rand_val:.4f}")
        #     self._debug_count += 1
        # elif not hasattr(self, '_debug_count'):
        #     self._debug_count = 1
        #     print(f"DEBUG age_survival: age_tol={self.age_tol}, k={k}, agemax={agemax}, check[k]={check[k]}, mort_prob={mort_prob:.4f}, rand={rand_val:.4f}")
        
        if rand_val < mort_prob:
            return False
        else:
            return True
    
    def growth_survival(self):
        """Check survival based on growth stress tolerance."""
        check = [0.31, 0.34, 0.37, 0.40, 0.43]  # 5%, 10%, 20%, 40%, 80% years
        
        k = max(0, min(len(check) - 1, self.stress_tol - 1))  # Clamp to valid range
        rand_val = urand()
        
        # DEBUG: Print survival details for first few trees
        # if hasattr(self, '_debug_growth_count') and self._debug_growth_count < 3:
        #     print(f"DEBUG growth_survival: stress_tol={self.stress_tol}, k={k}, check[k]={check[k]}, mort_marker={self.mort_marker}, rand={rand_val:.4f}")
        #     self._debug_growth_count += 1
        # elif not hasattr(self, '_debug_growth_count'):
        #     self._debug_growth_count = 1
        #     print(f"DEBUG growth_survival: stress_tol={self.stress_tol}, k={k}, check[k]={check[k]}, mort_marker={self.mort_marker}, rand={rand_val:.4f}")
        
        if self.mort_marker and rand_val < check[k]:
            return False
        else:
            return True
    
    def max_growth(self):
        """Calculate maximum growth potential."""
        dc = self.diam_canht
        d = self.diam_bht
        h = self.forska_ht
        dm = self.max_diam
        hm = self.max_ht
        g = self.g
        arfa0 = self.arfa_0
        
        beginning = arfa0 * math.exp(-arfa0 * d / (hm - STD_HT)) * d
        self.diam_max = g * d * (1.0 - d * h / dm / hm) / (2.0 * h + beginning)
    
    def get_diam_category(self):
        """Get diameter category for classification."""
        dm = self.diam_bht
        diam_category = np.zeros(NHC, dtype=int)
        
        if self.mort_marker:
            diam_category[0] = 1  # Dead trees
        
        if dm <= 8.0:
            diam_category[1] = 1
        elif dm <= 28.0:
            diam_category[2] = 1
        elif dm <= 48.0:
            diam_category[3] = 1
        elif dm <= 68.0:
            diam_category[4] = 1
        elif dm <= 88.0:
            diam_category[5] = 1
        else:
            diam_category[6] = 1
        
        return diam_category
    
    def calculate_all_metrics(self):
        """Calculate all tree metrics in proper order."""
        # Calculate height first
        self.forska_height()

        # Calculate stem shape
        self.stem_shape()

        # Calculate biomass components
        self.biomass_c()
        self.biomass_n()
        self.leaf_biomass_c()

        # Calculate growth potential
        self.max_growth()

    def is_alive(self):
        """Check if tree is alive based on survival factors."""
        if not self.age_survival():
            return False
        
        if not self.growth_survival():
            return False
        
        return True
    
    def __str__(self):
        """String representation of tree."""
        return (f"Tree({self.genus_name} {self.taxonomic_name}, "
                f"DBH={self.diam_bht:.1f}cm, H={self.forska_ht:.1f}m, "
                f"BiomC={self.biomC:.3f})")
    
    def __repr__(self):
        """Detailed representation of tree."""
        return self.__str__()