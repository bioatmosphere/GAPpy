"""
Model module for GAPpy vegetation model.
Orchestrates biogeochemical processes, forest dynamics, and ecosystem interactions.
"""

import math
import numpy as np
from typing import List, Dict, Optional
from .constants import *
from .parameters import params
from .species import SpeciesData
from .site import SiteData
from .tree import TreeData
from .plot import PlotData
from .climate import ex_rad, hargrea, cov365, cov365a
from .soil import SoilData
from .random_utils import urand, set_site_rng_seed


class ForestModel:
    """Main forest dynamics model orchestrating all processes."""
    
    def __init__(self):
        self.random_seed = None
        from .random_utils import clim_urand
        self.climate_uniform_func = clim_urand
    
    def set_site_rng_seed(self, fixed_seed: bool = False):
        """Set random number generator seed for site."""
        if fixed_seed:
            self.random_seed = 42
            set_site_rng_seed(fixed_seed=True, seed=self.random_seed)
        else:
            self.random_seed = None
            set_site_rng_seed(fixed_seed=False)
    
    def bio_geo_climate(self, site: SiteData, year: int):
        """Process biogeochemical climate interactions.

        Matches Fortran BioGeoClimate order:
        1. Prepare daily climate arrays
        2. Daily loop: soil water, decomp, accumulate degday/growdays/drydays/flooddays/PET
        3. Species climate responses using computed values
        """
        # Calculate daily climate variables
        self.calculate_daily_climate(site, year)

        # Process soil biogeochemistry (daily loop computes degday, growdays,
        # drydays, flooddays, PET — matching Fortran Model.f90 lines 161-236)
        self.process_soil_biogeochemistry(site)

        # Update species climate response factors using computed values
        # (matches Fortran Model.f90 lines 219-222)
        self.update_species_climate_responses(site)
    
    def calculate_daily_climate(self, site: SiteData, year: int):
        """Calculate daily climate variables from monthly data with stochastic variability."""
        from .random_utils import clim_nrand

        # Apply stochastic climate fluctuations to monthly values
        tmin_with_var = np.zeros(12)
        tmax_with_var = np.zeros(12)
        precip_with_var = np.zeros(12)

        # Calculate atmospheric N deposition from monthly precipitation (Fortran Model.f90:124)
        rain_n = 0.0

        for month in range(12):
            # Generate random fluctuation factors
            temp_f = clim_nrand(0.0, 1.0)
            prcp_f = clim_nrand(0.0, 1.0)

            # Clamp temperature fluctuations to [-1, 1]
            temp_f = max(-1.0, min(temp_f, 1.0))

            # Clamp precipitation fluctuations to [-0.5, 0.5]
            prcp_f = max(-0.5, min(prcp_f, 0.5))

            # Apply fluctuations scaled by standard deviations
            # If std arrays don't exist, initialize with zeros
            if not hasattr(site, 'tmin_std') or site.tmin_std is None:
                site.tmin_std = np.zeros(12)
                site.tmax_std = np.zeros(12)
                site.precip_std = np.zeros(12)

            tmin_with_var[month] = site.tmin[month] + temp_f * site.tmin_std[month]
            tmax_with_var[month] = site.tmax[month] + temp_f * site.tmax_std[month]
            precip_with_var[month] = max(site.precip[month] + prcp_f * site.precip_std[month], 0.0)

            # Accumulate atmospheric N deposition from monthly precip (matches Fortran)
            rain_n += precip_with_var[month] * PRCP_N

        # Convert monthly temperature to daily (using climate with variability)
        daily_tmin = cov365(tmin_with_var)
        daily_tmax = cov365(tmax_with_var)

        # Convert monthly precipitation to daily (with randomness and variability)
        daily_precip = cov365a(precip_with_var, self.climate_uniform_func)

        # Store daily climate for biogeochemistry processing
        site.daily_tmin = daily_tmin
        site.daily_tmax = daily_tmax
        site.daily_precip = daily_precip / 10.0  # Convert mm to cm

        # Store atmospheric N deposition (calculated from monthly precip)
        site.rain_n = rain_n

        # Store total annual rainfall as sum of monthly precip (matches Fortran Model.f90:123,231)
        site.rain = np.sum(precip_with_var)
    
    def calculate_degree_days(self, site: SiteData):
        """Calculate growing degree days."""
        # Base temperature for growing degree days
        base_temp = 5.0
        
        # Calculate from monthly data
        monthly_avg = (site.tmin + site.tmax) / 2.0
        
        # DEBUG: Print climate data for the first site
        if not hasattr(self, '_debug_climate_printed'):
            print(f"DEBUG climate: tmin={site.tmin}, tmax={site.tmax}")
            print(f"DEBUG climate: monthly_avg={monthly_avg}")
            self._debug_climate_printed = True
        
        # Sum degree days above base temperature
        degree_days = np.sum(np.maximum(monthly_avg - base_temp, 0) * 30)  # 30 days/month
        site.deg_days = degree_days
        
        # Calculate growing days (days above base temperature)
        site.grow_days = np.sum(monthly_avg > base_temp) * 30
        
        # DEBUG: Print degree days calculation
        if hasattr(self, '_debug_climate_printed') and not hasattr(self, '_debug_degdays_printed'):
            print(f"DEBUG degdays: degree_days={degree_days}, grow_days={site.grow_days}")
            self._debug_degdays_printed = True
    
    def calculate_potential_evapotranspiration(self, site: SiteData):
        """Calculate potential evapotranspiration using Hargreaves method."""
        # Use middle of year for radiation calculation
        julian_day = 180  # Mid-year approximation
        
        # Calculate extraterrestrial radiation
        erad, daylength, exradmx = ex_rad(julian_day, site.latitude)
        
        # Calculate monthly potential evapotranspiration
        monthly_pet = []
        for month in range(12):
            tmin = site.tmin[month]
            tmax = site.tmax[month]
            tavg = (tmin + tmax) / 2.0
            
            # Hargreaves equation
            pet = hargrea(tmin, tmax, tavg, erad)
            monthly_pet.append(pet)
        
        # Average daily potential evapotranspiration
        site.pot_evap_day = np.mean(monthly_pet)
    
    def update_species_climate_responses(self, site: SiteData):
        """Update species climate response factors.

        Uses degday, drydays, and flooddays computed by process_soil_biogeochemistry
        from the actual daily soil water balance (matches Fortran Model.f90:219-222).
        """
        num_species = len(site.species)
        for k in range(num_species):
            site.species[k].temp_rsp(site.deg_days)
            site.species[k].drought_rsp(site.dry_days_upper_layer,
                                        site.dry_days_base_layer)
            site.species[k].flood_rsp(site.flood_days)
    
    def process_soil_biogeochemistry(self, site: SiteData):
        """Process soil biogeochemical cycles with daily time steps.

        Matches Fortran Model.f90 BioGeoClimate daily loop (lines 161-236):
        - Uses Hargreaves PET with per-day extraterrestrial radiation
        - Accumulates degday, growdays, drydays, flooddays from soil water state
        - Normalizes drydays/flooddays by growing season length
        """
        # Transfer previous year's litter into A0 pool (Fortran Model.f90 lines 158-159)
        site.soil.A0_c0 += site.soil.C_into_A0
        site.soil.A0_n0 += site.soil.N_into_A0

        # Initialize accumulators (Fortran Model.f90 lines 142-156)
        total_C_resp = 0.0
        total_avail_N = 0.0
        total_pet = 0.0
        total_aet = 0.0
        degday = 0.0
        freeze = 0.0  # Local variable, always 0.0 (matches Fortran Model.f90:147)
        growdays = 0.0
        drydays_upper = 0.0
        drydays_base = 0.0
        flooddays = 0.0
        outwater = 0.0

        min_grow_temp = 5.0
        max_dry_parm = 1.0001
        min_flood_parm = 0.9999

        # Daily loop matching Fortran Model.f90 lines 161-206
        days_per_year = 365
        for day in range(days_per_year):
            julia = day + 1  # 1-based Julian day

            # Get daily climate values
            day_tmin = site.daily_tmin[day]
            day_tmax = site.daily_tmax[day]
            day_precip = site.daily_precip[day]  # Already in cm
            day_temp = (day_tmin + day_tmax) / 2.0

            # Extraterrestrial radiation and Hargreaves PET (Fortran Model.f90:165-167)
            erad, daylength, exradmx = ex_rad(julia, site.latitude)
            pot_ev_day = hargrea(day_tmin, day_tmax, day_temp, erad)

            # Process soil water balance with daily climate
            # Uses freeze=0.0 (local), matching Fortran Model.f90:170
            water_results = site.soil.soil_water(
                site.slope, site.leaf_area_ind, site.leaf_area_w0,
                site.sigma, freeze, day_precip, pot_ev_day
            )

            # Extract soil water factors
            act_ev_day = water_results[0]
            laiw0_scaled_by_max = water_results[1]
            laiw0_scaled_by_min = water_results[2]
            aow0_scaled_by_max = water_results[3]
            aow0_scaled_by_min = water_results[4]
            sbw0_scaled_by_max = water_results[5]
            sbw0_scaled_by_min = water_results[6]
            saw0_scaled_by_fc = water_results[7]
            saw0_scaled_by_wp = water_results[8]

            # Update canopy water content (critical!)
            site.leaf_area_w0 = water_results[9]  # Updated lai_w0

            # Process soil decomposition (Fortran Model.f90:178-180)
            avail_N, C_resp = site.soil.soil_decomp(
                0.0, 0.0, 0.0, 0.0,
                day_temp, day_precip, aow0_scaled_by_max, saw0_scaled_by_fc,
                sbw0_scaled_by_max
            )

            # Accumulate daily outputs (Fortran Model.f90:182-186)
            outwater += site.soil.runoff
            total_avail_N += max(avail_N, 0.0)
            total_pet += pot_ev_day
            total_aet += act_ev_day
            total_C_resp += C_resp

            # Accumulate degree days and growing season (Fortran Model.f90:189-192)
            if day_temp >= min_grow_temp:
                degday += (day_temp - min_grow_temp)
                growdays += 1.0

            # Accumulate dry days from soil water state (Fortran Model.f90:193-200)
            if (saw0_scaled_by_fc < max_dry_parm and
                    sbw0_scaled_by_min < max_dry_parm and
                    sbw0_scaled_by_max < max_dry_parm):
                drydays_upper += 1.0

            if saw0_scaled_by_wp < max_dry_parm:
                drydays_base += 1.0

            # Accumulate flood days (Fortran Model.f90:202-204)
            # Note: Fortran uses kron(yxd3) where yxd3 is never set (known bug),
            # so flooddays effectively stays 0. We replicate this behavior.
            # if aow0_scaled_by_min > min_flood_parm:
            #     flooddays += kron(yxd3)  # yxd3 uninitialized → 0.0 → kron(0)=0

        # Normalize by growing season (Fortran Model.f90:208-217)
        if growdays == 0:
            drydays_upper = 0.0
            drydays_base = 0.0
            flooddays = 0.0
        else:
            tmp = max(min(site.rain / total_pet, 1.0),
                      min(total_aet / total_pet, 1.0)) if total_pet > 0 else 1.0
            drydays_upper = min(drydays_upper / growdays, 1.0 - tmp)
            drydays_base = drydays_base / growdays
            flooddays = flooddays / growdays

        # Store annual totals (Fortran Model.f90:225-236)
        site.soil.avail_N = total_avail_N + site.rain_n
        site.soil.total_C_rsp = total_C_resp
        site.soil.runoff = outwater
        site.pot_evap_day = total_pet
        site.act_evap_day = total_aet
        site.grow_days = growdays
        site.deg_days = degday
        site.flood_days = flooddays
        site.dry_days_upper_layer = drydays_upper
        site.dry_days_base_layer = drydays_base
    
    def canopy(self, site: SiteData, year: int = -1):
        """Process canopy light interactions - complete Fortran translation."""
        # Light extinction parameter from Fortran
        xt = -0.40

        site.leaf_area_ind = 0.0
        num_species = len(site.species)

        # Debug for year 0
        debug_canopy = (year == 0 and hasattr(self, '_debug_year0') and self._debug_year0)

        for plot_idx, plot in enumerate(site.plots):
            ntrees = len(plot.trees)
            
            if ntrees == 0:
                # No trees - full light availability
                plot.con_light.fill(1.0)
                plot.dec_light.fill(1.0)
                plot.nutrient = [1.0] * num_species
            else:
                # Initialize leaf area density arrays
                maxheight = len(plot.con_light)
                lvd_c1 = np.zeros(maxheight)  # Deciduous LAI distribution
                lvd_c2 = np.zeros(maxheight)  # Coniferous LAI distribution  
                lvd_c3 = np.zeros(maxheight)  # Cumulative deciduous LAI
                lvd_c4 = np.zeros(maxheight)  # Cumulative coniferous LAI
                
                # Calculate LAI distribution for each tree
                for tree in plot.trees:
                    if tree.mort_marker:
                        continue
                        
                    forht = tree.forska_ht
                    canht = tree.canopy_ht
                    
                    # Calculate tree LAI using lai_biomass_c equivalent
                    tlai = tree.lai_biomass_c()
                    site.leaf_area_ind += tlai
                    
                    # Distribute LAI across height layers in canopy
                    jtmp = max(int(forht) - int(canht) + 1, 1)
                    lvd_adj = tlai / float(jtmp)
                    
                    # Add LAI to appropriate height layers
                    canht_int = max(0, min(int(canht), maxheight - 1))
                    forht_int = max(0, min(int(forht), maxheight - 1))
                    
                    if tree.conifer:
                        # Conifers contribute equally to both light arrays
                        for ih in range(canht_int, forht_int + 1):
                            if ih < maxheight:
                                lvd_c1[ih] += lvd_adj
                                lvd_c2[ih] += lvd_adj
                    else:
                        # Deciduous: full contribution to deciduous, 80% to coniferous
                        for ih in range(canht_int, forht_int + 1):
                            if ih < maxheight:
                                lvd_c1[ih] += lvd_adj
                                lvd_c2[ih] += lvd_adj * 0.8
                
                # Calculate cumulative LAI from top down
                lvd_c3[maxheight - 1] = lvd_c1[maxheight - 1]
                lvd_c4[maxheight - 1] = lvd_c2[maxheight - 1]
                
                for ih in range(1, maxheight):
                    lvd_c3[maxheight - ih - 1] = lvd_c3[maxheight - ih] + lvd_c1[maxheight - ih - 1]
                    lvd_c4[maxheight - ih - 1] = lvd_c4[maxheight - ih] + lvd_c2[maxheight - ih - 1]
                
                # Calculate light availability using Beer-Lambert law
                for ih in range(maxheight - 1):
                    if ih + 1 < maxheight:
                        plot.dec_light[ih] = math.exp(xt * lvd_c3[ih + 1] / params.plotsize)
                        plot.con_light[ih] = math.exp(xt * lvd_c4[ih + 1] / params.plotsize)

                # Debug first plot light arrays
                if debug_canopy and plot_idx == 0:
                    print(f"\n=== DEBUG Canopy Light (Plot 0, Year 0) ===")
                    print(f"  Number of trees: {ntrees}")
                    print(f"  Max height: {maxheight}")
                    print(f"  Total LAI: {site.leaf_area_ind:.4f}")
                    print(f"  dec_light[0:10]: {plot.dec_light[0:10]}")
                    print(f"  con_light[0:10]: {plot.con_light[0:10]}")
                    print(f"  lvd_c3[0:10] (cumulative dec LAI): {lvd_c3[0:10]}")
                    print(f"  lvd_c4[0:10] (cumulative con LAI): {lvd_c4[0:10]}")

        # Calculate site-level leaf area index
        site.leaf_area_ind = site.leaf_area_ind / float(site.numplots) / params.plotsize
    
    def growth(self, site: SiteData, year: int = -1):
        """Process tree growth - complete Fortran translation."""
        from .constants import HEC_TO_M2, STEM_C_N, CON_LEAF_C_N, DEC_LEAF_C_N, CON_LEAF_RATIO

        # Enable debug output for year 0
        if year == 0 and not hasattr(self, '_debug_year0'):
            self._debug_year0 = True
        elif year != 0:
            self._debug_year0 = False

        num_species = len(site.species)
        
        # Initialize soil carbon and nitrogen tracking
        site.soil.C_into_A0 = 0.0
        site.soil.N_into_A0 = 0.0
        
        N_used = 0.0
        net_prim_prodC = 0.0
        net_prim_prodN = 0.0
        biomc = 0.0
        biomn = 0.0
        
        leaf_b = 1.0 + CON_LEAF_RATIO
        growth_thresh = 0.05
        
        for plot in site.plots:
            # Reset available species tracking
            plot.avail_spec = [0.0] * num_species
            ntrees = len(plot.trees)
            N_req = 0.0
            
            if ntrees > 0:
                # Arrays to store tree data during processing
                N_stress = [0.0] * ntrees
                bleaf = [0.0] * ntrees  # Previous leaf biomass
                diam = [0.0] * ntrees   # Previous diameter
                biom_C = [0.0] * ntrees # Previous biomass C
                forska_shade = [0.0] * ntrees
                forht = [0.0] * ntrees
                khc = [0] * ntrees  # Canopy height index
                kh = [0] * ntrees   # Tree height index
                
                # FIRST TREE LOOP: Initial calculations and N requirements
                debug_first_tree = True  # Debug flag for first tree only
                for it, tree in enumerate(plot.trees):
                    if tree.mort_marker:
                        continue

                    k = tree.species_index
                    tree.update_tree(site.species[k])  # Update tree with current species data

                    # Store initial values
                    diam[it] = tree.diam_bht
                    canht = tree.canopy_ht
                    forht[it] = tree.forska_ht
                    
                    # Update available species tracking (kron: 1 if any tree exceeds threshold, else 0)
                    max_threshold = site.species[k].max_diam * growth_thresh
                    if diam[it] > max_threshold:
                        plot.avail_spec[k] = 1.0
                    
                    # Height indices for light calculations
                    khc[it] = min(int(canht), len(plot.con_light) - 1)
                    kh[it] = min(int(forht[it]), len(plot.con_light) - 1)
                    
                    # Calculate leaf biomass and max growth
                    tree.leaf_biomass_c()
                    tree.max_growth()
                    
                    # Light competition effects
                    if tree.conifer:
                        canopy_shade = site.species[k].light_rsp(plot.con_light[kh[it]])
                        forska_shade[it] = site.species[k].light_rsp(plot.con_light[khc[it]])
                    else:
                        canopy_shade = site.species[k].light_rsp(plot.dec_light[kh[it]])
                        forska_shade[it] = site.species[k].light_rsp(plot.dec_light[khc[it]])
                    
                    # Environmental stress calculation
                    N_stress[it] = tree.env_stress(canopy_shade)

                    # Debug first tree in year 0
                    if debug_first_tree and hasattr(self, '_debug_year0') and self._debug_year0:
                        print(f"\n=== DEBUG First Tree (Loop 1) ===")
                        print(f"  Species: {tree.genus_name} (k={k})")
                        print(f"  Initial DBH: {diam[it]:.2f} cm")
                        print(f"  fc_degday: {tree.fc_degday:.4f}")
                        print(f"  fc_drought: {tree.fc_drought:.4f}")
                        print(f"  fc_flood: {tree.fc_flood:.4f}")
                        print(f"  canopy_shade: {canopy_shade:.4f}")
                        print(f"  N_stress[it] = {tree.fc_degday:.4f} × {tree.fc_drought:.4f} × {tree.fc_flood:.4f} × {canopy_shade:.4f} = {N_stress[it]:.6f}")
                        print(f"  diam_max: {tree.diam_max:.4f} cm/yr")
                        debug_first_tree = False

                    # Increment diameter with growth and stress
                    tree.diam_bht = diam[it] + tree.diam_max * N_stress[it]
                    
                    # Update height for new diameter
                    tree.forska_height()
                    
                    # Save current leaf biomass, then update
                    bleaf[it] = tree.leaf_bm
                    tree.stem_shape()
                    tree.leaf_biomass_c()
                    
                    # Calculate nitrogen requirements for leaf/root growth
                    if site.species[k].conifer:
                        N_req += (leaf_b * tree.leaf_bm - bleaf[it]) / CON_LEAF_C_N
                    else:
                        N_req += tree.leaf_bm / DEC_LEAF_C_N
                    
                    # Save old biomass and calculate new
                    biom_C[it] = tree.biomC
                    tree.biomass_c()
                    tree.biomass_n()
                    
                    # Add wood growth N requirement
                    N_req += (tree.biomC - biom_C[it]) / STEM_C_N
                
                # Calculate nutrient availability
                N_req = max(N_req * HEC_TO_M2 / params.plotsize, 0.00001)
                N_supply_demand = site.soil.avail_N / N_req

                # Debug nitrogen in year 0
                if hasattr(self, '_debug_year0') and self._debug_year0:
                    print(f"\n=== DEBUG Nitrogen (between loops) ===")
                    print(f"  Available N: {site.soil.avail_N:.6f} tn/ha")
                    print(f"  Required N: {N_req:.6f} tn/ha")
                    print(f"  Supply:Demand ratio: {N_supply_demand:.4f}")

                # Update nutrient stress for all species
                for i in range(num_species):
                    plot.nutrient[i] = site.species[i].poor_soil_rsp(N_supply_demand)
                    if hasattr(self, '_debug_year0') and self._debug_year0 and i < 3:
                        print(f"  Species {i} ({site.species[i].genus_name}): nutrient = {plot.nutrient[i]:.6f}")
                
                # SECOND TREE LOOP: Final diameter increments and biomass
                debug_first_tree2 = True  # Debug flag for first tree in loop 2
                for it, tree in enumerate(plot.trees):
                    if tree.mort_marker:
                        continue

                    k = tree.species_index
                    tree.update_tree(site.species[k])  # Update tree with current species data

                    # Combined stress factor
                    fc_n = N_stress[it] * plot.nutrient[k]
                    dt = fc_n * tree.diam_max

                    # Set final diameter
                    tree.diam_bht = diam[it] + dt

                    # Check mortality threshold
                    pp = min(site.species[k].max_diam / site.species[k].max_age * 0.1, growth_thresh)

                    # Debug first tree in year 0
                    if debug_first_tree2 and hasattr(self, '_debug_year0') and self._debug_year0:
                        print(f"\n=== DEBUG First Tree (Loop 2 - Mortality Check) ===")
                        print(f"  N_stress[it]: {N_stress[it]:.6f}")
                        print(f"  plot.nutrient[k]: {plot.nutrient[k]:.6f}")
                        print(f"  fc_n = {N_stress[it]:.6f} × {plot.nutrient[k]:.6f} = {fc_n:.6f}")
                        print(f"  dt (diameter increment) = {fc_n:.6f} × {tree.diam_max:.4f} = {dt:.6f} cm")
                        print(f"  pp (minimum threshold) = min({site.species[k].max_diam:.1f}/{site.species[k].max_age:.0f} × 0.1, 0.05) = {pp:.6f} cm")
                        print(f"  Mortality check: dt <= pp? {dt:.6f} <= {pp:.6f}? {dt <= pp}")
                        print(f"  Mortality check: fc_n <= 0.05? {fc_n:.6f} <= 0.05? {fc_n <= growth_thresh}")
                        print(f"  DIES: {dt <= pp or fc_n <= growth_thresh}")
                        debug_first_tree2 = False
                    
                    # DEBUG: Show growth calculation for first few trees
                    # if hasattr(self, '_debug_growth_calc') and self._debug_growth_calc < 3:
                    #     print(f"DEBUG growth: fc_n={fc_n:.4f}, diam_max={tree.diam_max:.4f}, dt={dt:.4f}, pp={pp:.4f}, dt<=pp={dt<=pp}, fc_n<=thresh={fc_n<=growth_thresh}")
                    #     self._debug_growth_calc += 1
                    # elif not hasattr(self, '_debug_growth_calc'):
                    #     self._debug_growth_calc = 1
                    #     print(f"DEBUG growth: fc_n={fc_n:.4f}, diam_max={tree.diam_max:.4f}, dt={dt:.4f}, pp={pp:.4f}, dt<=pp={dt<=pp}, fc_n<=thresh={fc_n<=growth_thresh}")
                    
                    if dt <= pp or fc_n <= growth_thresh:
                        tree.mort_marker = True
                    else:
                        tree.mort_marker = False
                    
                    if not tree.mort_marker:
                        # Final height and shape calculations
                        tree.forska_height()
                        tree.stem_shape()
                        tree.leaf_biomass_c()
                        leafbm = tree.leaf_bm
                        
                        tree.biomass_c()
                        tree.biomass_n()
                        
                        # Calculate delta biomass and accumulate
                        d_bioC = tree.biomC - biom_C[it]
                        net_prim_prodC += d_bioC
                        N_used += d_bioC / STEM_C_N
                        
                        if tree.conifer:
                            prim_prod = leaf_b * leafbm - bleaf[it]
                            net_prim_prodC += prim_prod
                            N_used += prim_prod / CON_LEAF_C_N
                            
                            biomc += tree.biomC + leafbm
                            biomn += tree.biomN + leafbm / CON_LEAF_C_N
                        else:
                            net_prim_prodC += leafbm
                            N_used += leafbm / DEC_LEAF_C_N
                            
                            biomc += tree.biomC
                            biomn += tree.biomN
                
                # THIRD TREE LOOP: Canopy height and litter fall
                for it, tree in enumerate(plot.trees):
                    if tree.mort_marker:
                        continue
                        
                    k = tree.species_index
                    tree.update_tree(site.species[k])  # Update tree with current species data
                    forht[it] = tree.forska_ht
                    
                    # Environmental check for canopy growth
                    check = (site.species[k].fc_degday * site.species[k].fc_drought * 
                            site.species[k].fc_flood * forska_shade[it] * plot.nutrient[k])
                    
                    if check <= growth_thresh:
                        khc[it] += 1
                        if khc[it] < int(forht[it]):
                            # Increment canopy height
                            tree.canopy_ht = float(khc[it]) + 0.01
                            
                            tree.stem_shape()
                            bct = tree.biomC  # Save old biomass
                            tree.biomass_c()
                            tree.biomass_n()
                            
                            # Litter fall from canopy height change
                            d_bc = bct - tree.biomC
                            site.soil.C_into_A0 += d_bc
                            site.soil.N_into_A0 += d_bc / STEM_C_N
                            net_prim_prodC -= d_bc
                            
                            # Leaf litter
                            leafbm = tree.leaf_bm
                            tree.leaf_biomass_c()
                            d_leafb = leafbm - tree.leaf_bm
                            
                            if tree.conifer:
                                site.soil.C_into_A0 += d_leafb * leaf_b
                                site.soil.N_into_A0 += d_leafb / CON_LEAF_C_N * leaf_b
                                net_prim_prodC -= d_leafb * leaf_b
                            else:
                                site.soil.C_into_A0 += d_leafb
                                site.soil.N_into_A0 += d_leafb / DEC_LEAF_C_N
                                net_prim_prodC -= d_leafb
        
        # Final unit conversions and soil updates
        uconvert = HEC_TO_M2 / params.plotsize / site.numplots
        N_used *= uconvert
        site.soil.biomc = biomc * uconvert
        site.soil.biomn = biomn * uconvert
        site.soil.net_prim_prodC = net_prim_prodC * uconvert
        site.soil.net_C_into_A0 = site.soil.C_into_A0 * uconvert
        site.soil.N_used = N_used
        site.soil.net_prim_prodN = N_used
        site.soil.avail_N = site.soil.avail_N - site.soil.N_used
    
    def mortality(self, site: SiteData):
        """Process tree mortality - complete Fortran translation."""
        from .constants import CON_LEAF_RATIO, STEM_C_N, CON_LEAF_C_N, DEC_LEAF_C_N, HEC_TO_M2
        
        num_species = len(site.species)
        leaf_b = 1.0 + CON_LEAF_RATIO
        biomc = 0.0
        biomn = 0.0
        NPP_loss = 0.0
        NPPn_loss = 0.0
        
        for plot in site.plots:
            # Check for plot-level disturbances (fire or wind)
            fire_prob = urand()  # Random number for fire
            wind_prob = urand()  # Random number for wind
            
            if fire_prob < site.fire_prob or wind_prob < site.wind_prob:
                # PLOT-LEVEL DISTURBANCE OCCURS

                if plot.numtrees > 0:
                    # Determine disturbance type (fire takes precedence)
                    if fire_prob < site.fire_prob:
                        # FIRE DISTURBANCE
                        plot.fire = 5

                        # Fire response for each species and seedling establishment
                        for is_idx in range(num_species):
                            site.species[is_idx].fire_rsp(1)  # Fire occurred
                            plot.seedling[is_idx] = (
                                site.species[is_idx].invader * 10.0 +
                                site.species[is_idx].sprout_num * plot.avail_spec[is_idx]
                            ) * site.species[is_idx].fc_fire

                        plot.wind = 0
                    else:
                        # WIND DISTURBANCE
                        plot.wind = 3

                        # Wind disturbance seedling establishment
                        for is_idx in range(num_species):
                            plot.seedling[is_idx] = (
                                site.species[is_idx].invader +
                                plot.seedling[is_idx] +
                                site.species[is_idx].sprout_num * plot.avail_spec[is_idx]
                            )

                        plot.fire = 0

                    # Calculate total biomass from all trees before they die
                    zc = 0.0  # Total carbon
                    zn = 0.0  # Total nitrogen

                    for tree in plot.trees:
                        k = tree.species_index
                        tree.update_tree(site.species[k])  # Update with current species data

                        # Calculate leaf biomass
                        tree.leaf_biomass_c()
                        tmp = tree.leaf_bm

                        if site.species[k].conifer:
                            zc += tree.biomC + tmp * leaf_b
                            zn += tree.biomC / STEM_C_N + tmp / CON_LEAF_C_N * leaf_b
                        else:
                            zc += tree.biomC + tmp
                            zn += tree.biomC / STEM_C_N + tmp / DEC_LEAF_C_N

                    # Transfer all biomass to soil
                    site.soil.C_into_A0 += zc
                    site.soil.N_into_A0 += zn
                    # Track disturbance-killed biomass for subtraction (Fortran lines 689-693)
                    biomc += zc
                    biomn += zn
                    NPP_loss += zc
                    NPPn_loss += zn

                    # Mark that seedlings exist for renewal (Fortran line 695)
                    plot.seedling_number = 1.0

                # All trees dead now (outside numtrees guard, matching Fortran line 700)
                plot.trees.clear()
                plot.numtrees = 0

            else:
                # NO PLOT-LEVEL DISTURBANCE - Individual tree mortality
                plot.fire = 0
                plot.wind = 0
                surviving_trees = []

                for tree in plot.trees:
                    k = tree.species_index
                    tree.update_tree(site.species[k])
                    
                    # Check individual tree survival
                    if tree.growth_survival() and tree.age_survival():
                        # TREE SURVIVES
                        surviving_trees.append(tree)
                        
                        # Calculate leaf litter from surviving trees
                        tree.leaf_biomass_c()
                        leaf_bm = tree.leaf_bm
                        
                        if site.species[k].conifer:
                            # Conifers drop 30% of leaves annually
                            litter_c = leaf_bm * (leaf_b - 1.0)
                            litter_n = litter_c / CON_LEAF_C_N
                        else:
                            # Deciduous trees drop all leaves annually
                            litter_c = leaf_bm
                            litter_n = litter_c / DEC_LEAF_C_N
                        
                        site.soil.C_into_A0 += litter_c
                        site.soil.N_into_A0 += litter_n
                        NPP_loss += litter_c
                        NPPn_loss += litter_n
                    
                    else:
                        # TREE DIES
                        bmc = tree.biomC
                        tree.leaf_biomass_c()
                        leaf_bm = tree.leaf_bm
                        
                        if site.species[k].conifer:
                            # Dead conifer: all biomass to soil
                            dead_c = bmc + leaf_bm * leaf_b
                            dead_n = bmc / STEM_C_N + leaf_bm / CON_LEAF_C_N * leaf_b
                        else:
                            # Dead deciduous: all biomass to soil
                            dead_c = bmc + leaf_bm
                            dead_n = bmc / STEM_C_N + leaf_bm / DEC_LEAF_C_N
                        
                        site.soil.C_into_A0 += dead_c
                        site.soil.N_into_A0 += dead_n
                        NPP_loss += dead_c
                        NPPn_loss += dead_n
                
                # Update plot with surviving trees
                plot.trees = surviving_trees
                plot.numtrees = len(surviving_trees)
        
        # Apply unit conversions and update soil totals (Fortran Model.f90:787-792)
        uconvert = HEC_TO_M2 / params.plotsize / site.numplots
        # Subtract disturbance-killed biomass from standing totals set in growth()
        site.soil.biomc = site.soil.biomc - biomc * uconvert
        site.soil.biomn = site.soil.biomn - biomn * uconvert
        site.soil.net_prim_prodC = site.soil.net_prim_prodC - NPP_loss * uconvert
        # net_prim_prodN adjustment commented out in Fortran (Model.f90 line 792)
        # site.soil.net_prim_prodN = site.soil.net_prim_prodN - NPPn_loss * uconvert
    
    def renewal(self, site: SiteData):
        """Process tree renewal - complete Fortran translation."""
        from .constants import HEC_TO_M2, STEM_C_N, CON_LEAF_C_N, DEC_LEAF_C_N, STD_HT, CON_LEAF_RATIO
        from .random_utils import nrand

        growth_thresh = 0.05
        growth_min = 0.01
        epsilon = 1e-10
        num_species = len(site.species)
        leaf_b = 1.0 + CON_LEAF_RATIO

        # Track totals across all plots
        new_count = 0
        net_prim_prod = 0.0
        N_used = 0.0

        # Plot loop only runs if available nitrogen exists (Fortran line 823)
        if site.soil.avail_N > 0.0:
            for plot in site.plots:
                # Check if recruitment can occur (normal conditions)
                can_recruit = (plot.numtrees != 0 or
                               (plot.wind == 0 and plot.fire == 0))

                # Initialize probability variables for this plot
                prob = np.zeros(num_species)
                probsum = 0.0
                nrenew = 0

                if can_recruit:
                    # Step 1: Calculate growth capacity and regrowth
                    regrowth = np.zeros(num_species)
                    growmax = 0.0

                    for i, species in enumerate(plot.species):
                        grow_cap = (species.fc_degday * species.fc_drought *
                                    species.fc_flood * plot.nutrient[i])

                        if species.conifer:
                            light_factor = species.light_rsp(plot.con_light[0])
                        else:
                            light_factor = species.light_rsp(plot.dec_light[0])

                        regrowth[i] = grow_cap * light_factor
                        growmax = max(growmax, regrowth[i])

                        if regrowth[i] <= growth_thresh:
                            regrowth[i] = 0.0

                    # Step 2: Determine maximum recruitment number
                    max_renew = min(int(params.plotsize * growmax) - plot.numtrees,
                                    int(params.plotsize * 0.5))
                    nrenew = min(max(max_renew, 3),
                                 int(params.plotsize) - plot.numtrees)
                    nrenew = max(0, min(nrenew, params.maxtrees - plot.numtrees))

                    # Seedbank/seedling pipeline differs on first vs subsequent cycles
                    if plot.seedling_number == 0.0:
                        # First recruitment cycle
                        for i, species in enumerate(plot.species):
                            plot.seedbank[i] += (
                                species.invader +
                                species.seed_num * plot.avail_spec[i] +
                                species.sprout_num * plot.avail_spec[i]
                            )

                            if regrowth[i] >= growth_thresh:
                                plot.seedling[i] += plot.seedbank[i]
                                plot.seedbank[i] = 0.0
                            else:
                                plot.seedbank[i] *= species.seed_surv

                            species.fire_rsp(plot.fire)
                            plot.seedling[i] = (
                                plot.seedling[i] +
                                species.sprout_num * plot.avail_spec[i]
                            ) * species.fc_fire

                            # kron: 1.0 if seedling > 0, else 0.0 (before scaling)
                            plot.seedling_number = max(
                                1.0 if plot.seedling[i] > 0.0 else 0.0,
                                plot.seedling_number)

                            plot.seedling[i] *= params.plotsize

                        # Calculate probabilities after seedbank update
                        for i in range(num_species):
                            prob[i] = plot.seedling[i] * regrowth[i]
                            probsum += prob[i]

                    else:
                        # Subsequent cycles: probabilities first, then update
                        for i in range(num_species):
                            prob[i] = plot.seedling[i] * regrowth[i]
                            probsum += prob[i]

                        for i, species in enumerate(plot.species):
                            plot.seedbank[i] += (
                                species.invader +
                                species.seed_num * plot.avail_spec[i] +
                                species.sprout_num * plot.avail_spec[i]
                            )

                            if regrowth[i] >= growth_min:
                                plot.seedling[i] += plot.seedbank[i]
                                plot.seedbank[i] = 0.0
                            else:
                                plot.seedbank[i] *= species.seed_surv

                            species.fire_rsp(plot.fire)
                            plot.seedling[i] = (
                                plot.seedling[i] +
                                species.sprout_num * plot.avail_spec[i]
                            ) * species.fc_fire

                            plot.seedling[i] *= params.plotsize

                            # kron after scaling (matching Fortran order)
                            plot.seedling_number = max(
                                1.0 if plot.seedling[i] > 0.0 else 0.0,
                                plot.seedling_number)

                else:
                    # Post-disturbance (numtrees == 0 AND recent disturbance)
                    if plot.fire == 1 or plot.wind == 1:
                        for i, species in enumerate(plot.species):
                            grow_cap = (species.fc_degday * species.fc_drought *
                                        species.fc_flood)

                            prob[i] = plot.seedling[i] * grow_cap
                            probsum += prob[i]

                            plot.seedling[i] *= params.plotsize

                            plot.seedling_number = max(
                                1.0 if plot.seedling[i] > 0.0 else 0.0,
                                plot.seedling_number)

                        plot.fire = 0
                        plot.wind = 0
                    else:
                        plot.fire = max(0, plot.fire - 1)
                        plot.wind = max(0, plot.wind - 1)

                # Normalize and create cumulative probabilities
                if probsum > epsilon:
                    for i in range(num_species):
                        prob[i] /= probsum
                    for i in range(1, num_species):
                        prob[i] += prob[i - 1]
                else:
                    nrenew = 0

                # Recruit new trees
                if nrenew >= 1:
                    for _ in range(nrenew):
                        q0 = urand()
                        selected_species = 0

                        while q0 > prob[selected_species]:
                            selected_species += 1
                            if selected_species >= num_species:
                                # Fallback: random reselection (Fortran line 1015)
                                selected_species = min(
                                    int(urand(0.0, float(num_species))),
                                    num_species - 1)
                                q0 = urand()

                        new_count += 1
                        plot.seedling[selected_species] -= 1.0

                        new_tree = TreeData()
                        new_tree.initialize_tree(
                            plot.species[selected_species], selected_species)

                        z = 1.5 + nrand(0.0, 1.0)
                        z = max(0.5, min(2.5, z))

                        new_tree.diam_bht = z
                        new_tree.canopy_ht = 1.0

                        new_tree.forska_height()
                        new_tree.stem_shape()
                        new_tree.biomass_c()
                        new_tree.biomass_n()
                        new_tree.leaf_biomass_c()

                        plot.trees.append(new_tree)
                        plot.numtrees = len(plot.trees)

                        # Leaf biomass for litter/NPP (Fortran lines 1046-1072)
                        zz = new_tree.leaf_bm

                        species = plot.species[selected_species]
                        if species.conifer:
                            net_prim_prod += zz * leaf_b + new_tree.biomC
                            N_used += zz / CON_LEAF_C_N + new_tree.biomN

                            site.soil.C_into_A0 += zz * (leaf_b - 1.0)
                            site.soil.N_into_A0 += zz * (leaf_b - 1.0) / CON_LEAF_C_N
                        else:
                            net_prim_prod += new_tree.biomC + zz
                            N_used += new_tree.biomN + zz / DEC_LEAF_C_N

                            site.soil.C_into_A0 += zz
                            site.soil.N_into_A0 += zz / DEC_LEAF_C_N

                # Seedling scaling and mortality for next cycle (Fortran lines 1077-1080)
                for i in range(num_species):
                    plot.seedling[i] = (plot.seedling[i] *
                                        plot.species[i].seedling_lg) / params.plotsize

                # Reset fire indicator (Fortran line 1083)
                plot.fire = 0

        # End-of-function soil updates - ALWAYS execute (Fortran lines 1089-1111)
        uconvert = HEC_TO_M2 / params.plotsize / site.numplots
        N_used *= uconvert
        net_prim_prod *= uconvert

        # N balance
        avtt = site.soil.avail_N - N_used
        if avtt > 0.0:
            site.soil.net_N_into_A0 = avtt * min(site.soil.runoff / 1000.0, 0.1)
            site.soil.A_n0 += avtt - site.soil.net_N_into_A0
        else:
            site.soil.A_n0 += avtt
            site.soil.net_N_into_A0 = 0.0

        # Runoff leaching
        site.soil.A_n0 -= 0.00002 * site.soil.runoff
        site.soil.A_c0 -= site.soil.net_N_into_A0 * 20.0
        site.soil.BL_c0 += site.soil.net_N_into_A0 * 20.0
        site.soil.BL_n0 += site.soil.net_N_into_A0

        # Unit conversions for litter inputs
        site.soil.C_into_A0 *= uconvert
        site.soil.N_into_A0 *= uconvert
        site.soil.N_used = N_used

        # Update NPP totals
        site.soil.net_prim_prodC += net_prim_prod
        site.soil.net_prim_prodN += N_used
        site.soil.new_growth = int(new_count * uconvert)
    
    def initialize_forest(self, site: SiteData):
        """Initialize forest with starting trees."""
        for plot in site.plots:
            # Add some initial trees of different species
            for i, species in enumerate(plot.species):
                # Number of initial trees (minimum 5 per species)
                n_initial = max(1, int(2 * species.invader + 1))
                
                for _ in range(min(n_initial, params.maxtrees // len(plot.species))):
                    tree = TreeData()
                    tree.initialize_tree(species, i)
                    
                    # Initialize with random sizes
                    tree.diam_bht = urand(5.0, 30.0)  # 5-30 cm DBH
                    tree.calculate_all_metrics()
                    
                    plot.add_tree(tree)
                    
            live_count = len([t for t in plot.trees if not t.mort_marker])
            print(f"  Initialized plot with {len(plot.trees)} trees ({live_count} alive)")
    
    def run_annual_cycle(self, site: SiteData, year: int):
        """Run complete annual cycle for a site."""
        # Enable debug for year 0
        if year == 0 and not hasattr(self, '_debug_year0'):
            self._debug_year0 = True
        elif year != 0:
            self._debug_year0 = False

        # Biogeochemical processes
        self.bio_geo_climate(site, year)

        # Forest dynamics
        self.canopy(site, year)  # Pass year for debugging
        self.growth(site, year)  # Pass year for debugging
        self.mortality(site)
        self.renewal(site)
        
        # Update site-level statistics
        self.update_site_statistics(site)
    
    def update_site_statistics(self, site: SiteData):
        """Update site-level forest statistics."""
        total_biomass_c = 0.0
        total_biomass_n = 0.0
        total_basal_area = 0.0
        total_trees = 0
        
        for plot in site.plots:
            plot_stats = plot.get_statistics()
            total_biomass_c += plot_stats['total_biomass_c']
            total_biomass_n += plot_stats['total_biomass_n']
            total_basal_area += plot_stats['total_basal_area']
            total_trees += plot_stats['total_trees']
        
        # Store site-level statistics
        site.total_biomass_c = total_biomass_c
        site.total_biomass_n = total_biomass_n
        site.total_basal_area = total_basal_area
        site.total_trees = total_trees