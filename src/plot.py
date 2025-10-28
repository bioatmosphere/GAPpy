"""
Plot module for UVAFME vegetation model.
Translated from Plot.f90
"""

import math
import numpy as np
from typing import List, Dict, Tuple
from .constants import *
from .species import SpeciesData
from .tree import TreeData


class PlotData:
    """Plot data structure containing trees and species information."""
    
    def __init__(self):
        self.trees: List[TreeData] = []
        self.species: List[SpeciesData] = []
        self.avail_spec = []
        self.seedbank = []
        self.seedling = []
        self.con_light = []
        self.dec_light = []
        self.nutrient = []
        self.seedling_number = 0.0
        self.numspecies = 0
        self.numtrees = 0
        self.fire = 0
        self.wind = 0
    
    def initialize_plot(self, species: List[SpeciesData], maxtrees: int, maxheight: int):
        """Initialize plot with species list and size constraints."""
        if maxtrees == 0:
            raise ValueError("Must allow at least a few trees")
        
        if maxheight == 0:
            raise ValueError("Must have a nonzero maximum height")
        
        # Initialize tree storage
        self.trees = []
        
        # Initialize height-based arrays
        self.con_light = np.zeros(maxheight)
        self.dec_light = np.zeros(maxheight)
        
        # For now, assume all species in site can potentially be in plot
        self.numspecies = len(species)
        self.seedling_number = 0.0
        
        # Initialize species-based arrays
        self.avail_spec = np.zeros(self.numspecies)
        self.nutrient = np.ones(self.numspecies)
        self.seedbank = np.zeros(self.numspecies)
        self.seedling = np.zeros(self.numspecies)
        
        # Copy species list
        self.species = species.copy()
        
        # Initialize counters and status
        self.numtrees = 0
        self.fire = 0
        self.wind = 0
        
        # Reset arrays
        self.seedbank = np.zeros(self.numspecies)
        self.seedling = np.zeros(self.numspecies)
        self.seedling_number = 0
        self.avail_spec = np.zeros(self.numspecies)
        self.dec_light = np.zeros(maxheight)
        self.con_light = np.zeros(maxheight)
        self.nutrient = np.ones(self.numspecies)
    
    def add_tree(self, tree: TreeData) -> bool:
        """Add a tree to the plot."""
        self.trees.append(tree)
        self.numtrees = len(self.trees)
        return True
    
    def remove_tree(self, index: int) -> bool:
        """Remove a tree from the plot."""
        if 0 <= index < len(self.trees):
            self.trees.pop(index)
            self.numtrees = len(self.trees)
            return True
        return False
    
    def tree_dm_cats(self, genera: List[str], field: str = "genus") -> np.ndarray:
        """Calculate diameter categories for specified genera or species."""
        numitems = len(genera)
        diam_categories = np.zeros((numitems, NHC), dtype=int)

        for ig, genus in enumerate(genera):
            # Extract the comparison value based on field type
            if field == "genus":
                # For genus, genera contains just genus names (strings)
                match_value = genus
            elif field == "species":
                # For species, genera contains tuples of (genus_name, species_id)
                # Extract just the species_id
                match_value = genus[1] if isinstance(genus, tuple) else genus
            else:
                continue

            for tree in self.trees:
                if field == "genus":
                    comp = tree.genus_name
                elif field == "species":
                    comp = tree.unique_id
                else:
                    continue

                if comp == match_value:
                    tree_diams = tree.get_diam_category()
                    diam_categories[ig, :] += tree_diams

        return diam_categories
    
    def sum_over_sg(self, genera: List[str], field: str = "genus") -> Dict[str, np.ndarray]:
        """Sum biomass and other metrics over species/genera."""
        numitems = len(genera)
        
        # Initialize output arrays
        basal_area = np.zeros(numitems)
        leaf_bm = np.zeros(numitems)
        biomC = np.zeros(numitems)
        biomN = np.zeros(numitems)
        biomC_std = np.zeros(numitems)
        biomN_std = np.zeros(numitems)
        max_ht = np.zeros(numitems)
        max_diam = np.zeros(numitems)
        n = np.zeros(numitems, dtype=int)
        
        # Sum over trees
        for is_idx, genus in enumerate(genera):
            biomC_sum_sq = 0.0
            biomN_sum_sq = 0.0

            # Extract the comparison value based on field type
            if field == "genus":
                # For genus, genera contains just genus names (strings)
                match_value = genus
            elif field == "species":
                # For species, genera contains tuples of (genus_name, species_id)
                # Extract just the species_id
                match_value = genus[1] if isinstance(genus, tuple) else genus
            else:
                continue

            for tree in self.trees:
                if field == "genus":
                    comp = tree.genus_name
                elif field == "species":
                    comp = tree.unique_id
                else:
                    continue

                if comp == match_value:
                    n[is_idx] += 1
                    
                    # Determine leaf C/N ratio
                    if tree.conifer:
                        l_cn = CON_LEAF_C_N
                    else:
                        l_cn = DEC_LEAF_C_N
                    
                    # Calculate metrics
                    basal_area[is_idx] += 0.25 * PI * tree.diam_bht**2
                    lf_biom = tree.leaf_bm
                    leaf_bm[is_idx] += lf_biom
                    tot_tree_bc = tree.biomC + lf_biom
                    tot_tree_bn = tree.biomN + lf_biom / l_cn
                    biomC[is_idx] += tot_tree_bc
                    biomN[is_idx] += tot_tree_bn
                    max_ht[is_idx] = max(max_ht[is_idx], tree.forska_ht)
                    max_diam[is_idx] = max(max_diam[is_idx], tree.diam_bht)
                    
                    # Sum of squares for standard deviation
                    biomC_sum_sq += tot_tree_bc**2
                    biomN_sum_sq += tot_tree_bn**2
            
            # Calculate standard deviations
            if n[is_idx] == 0:
                biomC_std[is_idx] = 0.0
                biomN_std[is_idx] = 0.0
            else:
                ni = 1.0 / float(n[is_idx] - 1) if n[is_idx] > 1 else 0.0
                # Protect against negative variance due to floating point errors
                var_c = max(0.0, (biomC_sum_sq - biomC[is_idx]**2 / n[is_idx]) * ni)
                var_n = max(0.0, (biomN_sum_sq - biomN[is_idx]**2 / n[is_idx]) * ni)
                biomC_std[is_idx] = math.sqrt(var_c)
                biomN_std[is_idx] = math.sqrt(var_n)
        
        return {
            'basal_area': basal_area,
            'leaf_bm': leaf_bm,
            'biomC': biomC,
            'biomN': biomN,
            'biomC_std': biomC_std,
            'biomN_std': biomN_std,
            'max_ht': max_ht,
            'max_diam': max_diam,
            'n_trees': n
        }
    
    def calculate_light_profile(self):
        """Calculate light profile through canopy."""
        # Reset light arrays
        self.con_light.fill(0.0)
        self.dec_light.fill(0.0)
        
        # Sort trees by height (tallest first)
        sorted_trees = sorted(self.trees, key=lambda t: t.forska_ht, reverse=True)
        
        # Calculate light attenuation through canopy layers
        available_light = 1.0  # Full sunlight at top
        
        for tree in sorted_trees:
            height_idx = min(int(tree.forska_ht), len(self.con_light) - 1)
            
            if tree.conifer:
                self.con_light[height_idx] = available_light
                # Conifers attenuate light more
                available_light *= 0.7
            else:
                self.dec_light[height_idx] = available_light
                # Deciduous trees attenuate less
                available_light *= 0.8
    
    def get_species_by_id(self, unique_id: str) -> SpeciesData:
        """Get species data by unique ID."""
        for species in self.species:
            if species.unique_id == unique_id:
                return species
        return None
    
    def get_canopy_shading(self, height: float) -> float:
        """Get canopy shading at a given height."""
        height_idx = min(int(height), len(self.con_light) - 1)
        # Return average of conifer and deciduous light availability
        # Convert from light availability to shading (1 - availability)
        con_light = self.con_light[height_idx] if height_idx < len(self.con_light) else 0.0
        dec_light = self.dec_light[height_idx] if height_idx < len(self.dec_light) else 0.0
        avg_light = (con_light + dec_light) / 2.0
        return 1.0 - avg_light  # Convert to shading
    
    def get_total_biomass(self) -> Tuple[float, float]:
        """Get total carbon and nitrogen biomass."""
        total_c = sum(tree.biomC for tree in self.trees)
        total_n = sum(tree.biomN for tree in self.trees)
        return total_c, total_n
    
    def get_total_basal_area(self) -> float:
        """Get total basal area."""
        return sum(0.25 * PI * tree.diam_bht**2 for tree in self.trees)
    
    def get_live_trees(self) -> List[TreeData]:
        """Get list of living trees."""
        return [tree for tree in self.trees if not tree.mort_marker]
    
    def get_dead_trees(self) -> List[TreeData]:
        """Get list of dead trees."""
        return [tree for tree in self.trees if tree.mort_marker]
    
    def remove_dead_trees(self):
        """Remove dead trees from plot."""
        self.trees = [tree for tree in self.trees if not tree.mort_marker]
        self.numtrees = len(self.trees)
    
    def get_statistics(self) -> Dict[str, float]:
        """Get basic plot statistics."""
        live_trees = self.get_live_trees()
        total_c, total_n = self.get_total_biomass()
        
        return {
            'total_trees': self.numtrees,
            'live_trees': len(live_trees),
            'dead_trees': len(self.get_dead_trees()),
            'total_biomass_c': total_c,
            'total_biomass_n': total_n,
            'total_basal_area': self.get_total_basal_area(),
            'avg_dbh': np.mean([tree.diam_bht for tree in live_trees]) if live_trees else 0.0,
            'avg_height': np.mean([tree.forska_ht for tree in live_trees]) if live_trees else 0.0,
            'max_dbh': max([tree.diam_bht for tree in live_trees]) if live_trees else 0.0,
            'max_height': max([tree.forska_ht for tree in live_trees]) if live_trees else 0.0
        }
    
    def __str__(self):
        """String representation of plot."""
        stats = self.get_statistics()
        return (f"Plot(Trees={stats['total_trees']}, "
                f"Live={stats['live_trees']}, "
                f"BiomC={stats['total_biomass_c']:.2f}, "
                f"BA={stats['total_basal_area']:.2f})")
    
    def __repr__(self):
        """Detailed representation of plot."""
        return self.__str__()