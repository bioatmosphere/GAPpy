"""
Model module for UVAFME vegetation model.
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
        """Process biogeochemical climate interactions."""
        # Calculate daily climate variables
        self.calculate_daily_climate(site, year)
        
        # Calculate growing degree days
        self.calculate_degree_days(site)
        
        # Calculate potential evapotranspiration
        self.calculate_potential_evapotranspiration(site)
        
        # Update species climate response factors
        self.update_species_climate_responses(site)
        
        # Debug: Check if environmental factors are being set (only year 0)
        if year == 0 and site.species and len(site.species) > 0:
            s = site.species[0]  # Check first species
            print(f"  DEBUG Year {year}: Species 0 - degday={s.fc_degday:.3f}, drought={s.fc_drought:.3f}, flood={s.fc_flood:.3f}")
        
        # Process soil biogeochemistry
        self.process_soil_biogeochemistry(site)
    
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

        # Calculate daily averages for the year
        site.rain = np.mean(daily_precip) / 10.0  # Convert to cm/day

        # Calculate temperature statistics
        daily_temp_avg = (daily_tmin + daily_tmax) / 2.0
        site.freeze = np.sum(daily_temp_avg < 0) / len(daily_temp_avg)
    
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
        """Update species climate response factors."""
        for plot in site.plots:
            for i, species in enumerate(plot.species):
                # DEBUG: Print species temperature parameters for first species
                if not hasattr(self, '_debug_species_printed') and i == 0:
                    print(f"DEBUG species {i}: deg_day_min={species.deg_day_min}, deg_day_opt={species.deg_day_opt}, deg_day_max={species.deg_day_max}")
                    print(f"DEBUG temp_rsp input: site.deg_days={site.deg_days}")
                    self._debug_species_printed = True
                
                # Temperature response
                species.temp_rsp(site.deg_days)
                
                # DEBUG: Print temperature response result
                if hasattr(self, '_debug_species_printed') and not hasattr(self, '_debug_temp_rsp_printed') and i == 0:
                    print(f"DEBUG temp_rsp result: fc_degday={species.fc_degday}")
                    self._debug_temp_rsp_printed = True
                
                # Drought response (simplified)
                dry_days = max(0, 100 - site.rain * 365)  # Estimate dry days
                dry_days_surface = dry_days * 0.8  # Surface drying faster
                species.drought_rsp(dry_days, dry_days_surface)
                
                # Flood response (simplified)
                flood_days = max(0, site.rain * 365 - 1000) / 100  # Estimate flood days
                species.flood_rsp(flood_days)
                
                # Fire response (based on site probability)
                fire_occurrence = 1 if urand() < site.fire_prob else 0
                species.fire_rsp(fire_occurrence)
    
    def process_soil_biogeochemistry(self, site: SiteData):
        """Process soil biogeochemical cycles with daily time steps."""
        # Calculate annual litter inputs from existing trees
        annual_litter_c1 = 0.0  # Above-ground litter
        annual_litter_c2 = 0.0  # Below-ground litter
        annual_litter_n1 = 0.0  # Above-ground N
        annual_litter_n2 = 0.0  # Below-ground N

        for plot in site.plots:
            for tree in plot.trees:
                if tree.mort_marker:
                    # Dead tree contributes to litter
                    annual_litter_c1 += tree.biomC * 0.7  # 70% above-ground
                    annual_litter_c2 += tree.biomC * 0.3  # 30% below-ground
                    annual_litter_n1 += tree.biomN * 0.7
                    annual_litter_n2 += tree.biomN * 0.3
                else:
                    # Living tree contributes small amount (leaf fall, etc.)
                    annual_turnover = 0.1  # 10% annual turnover
                    annual_litter_c1 += tree.leaf_bm * annual_turnover
                    annual_litter_n1 += tree.leaf_bm * annual_turnover / DEC_LEAF_C_N

        # Convert annual litter to daily inputs
        daily_litter_c1 = annual_litter_c1 / 365.0
        daily_litter_c2 = annual_litter_c2 / 365.0
        daily_litter_n1 = annual_litter_n1 / 365.0
        daily_litter_n2 = annual_litter_n2 / 365.0

        # Initialize accumulators for annual totals
        total_C_resp = 0.0
        total_avail_N = 0.0
        total_act_evap = 0.0

        # Daily loop matching Fortran Model.f90 lines 161-186
        days_per_year = 365
        for day in range(days_per_year):
            # Get daily climate values
            day_tmin = site.daily_tmin[day]
            day_tmax = site.daily_tmax[day]
            day_precip = site.daily_precip[day]  # Already in cm
            day_temp = (day_tmin + day_tmax) / 2.0

            # Daily potential evapotranspiration (from daily temp)
            # Use Hamon equation like in calculate_potential_evapotranspiration
            lat_rad = site.latitude * PI / 180.0
            day_of_year = day + 1

            # Solar declination
            decl = 0.409 * math.sin(2.0 * PI * day_of_year / 365.0 - 1.39)

            # Day length calculation
            # Clamp the argument to [-1, 1] to avoid domain error at extreme latitudes
            cos_arg = -math.tan(lat_rad) * math.tan(decl)
            cos_arg = max(-1.0, min(1.0, cos_arg))
            ws = math.acos(cos_arg)
            daylength_hours = 24.0 * ws / PI

            # Potential evapotranspiration (Hamon equation)
            if day_temp > 0:
                sat_vap_density = 0.622 * 6.108 * math.exp(17.27 * day_temp / (237.3 + day_temp))
                pot_ev_day = 0.1651 * daylength_hours * sat_vap_density / (day_temp + 273.3)
            else:
                pot_ev_day = 0.0

            # Process soil water balance with daily climate
            water_results = site.soil.soil_water(
                site.slope, site.leaf_area_ind, site.leaf_area_w0,
                site.sigma, site.freeze, day_precip, pot_ev_day
            )

            # Extract soil water factors
            act_ev_day = water_results[0]
            aow0_scaled = water_results[3]  # aow0_scaled_by_max
            saw0_scaled = water_results[7]  # saw0_scaled_by_fc
            sbw0_scaled = water_results[5]  # sbw0_scaled_by_max

            # Update canopy water content (critical for flood_days calculation!)
            site.leaf_area_w0 = water_results[9]  # Updated lai_w0

            # Process soil decomposition with daily inputs
            avail_N, C_resp = site.soil.soil_decomp(
                daily_litter_c1, daily_litter_c2, daily_litter_n1, daily_litter_n2,
                day_temp, day_precip, aow0_scaled, saw0_scaled, sbw0_scaled
            )

            # Accumulate daily outputs (matches Fortran lines 183-184)
            total_C_resp += C_resp
            total_avail_N += max(avail_N, 0.0)
            total_act_evap += act_ev_day

        # Store annual totals (matches Fortran line 186, 225)
        site.soil.total_C_rsp = total_C_resp
        site.soil.avail_N = total_avail_N + site.rain_n  # Add atmospheric N from monthly precip
        site.act_evap_day = total_act_evap / days_per_year  # Average daily evap

        # Calculate final water stress indicators from last day's state
        # (water_results from last iteration)
        site.dry_days_upper_layer = max(0, 100 - water_results[3] * 100)
        site.dry_days_base_layer = max(0, 100 - water_results[5] * 100)
        site.flood_days = max(0, water_results[1] * 100 - 80)
    
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
                    
                    # Update available species tracking
                    max_threshold = site.species[k].max_diam * growth_thresh
                    if diam[it] > max_threshold:
                        plot.avail_spec[k] = max(diam[it] - max_threshold, plot.avail_spec[k])
                    
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
                
                # Determine disturbance type (fire takes precedence)
                if fire_prob < site.fire_prob:
                    # FIRE DISTURBANCE
                    plot.fire = 5
                    plot.wind = 0
                    
                    # Fire response for each species and seedling establishment
                    for is_idx in range(num_species):
                        site.species[is_idx].fire_rsp(1)  # Fire occurred
                        plot.seedling[is_idx] = (
                            site.species[is_idx].invader * 10.0 +
                            site.species[is_idx].sprout_num * plot.avail_spec[is_idx]
                        ) * site.species[is_idx].fc_fire
                else:
                    # WIND DISTURBANCE  
                    plot.wind = 3
                    plot.fire = 0
                    
                    # Wind disturbance seedling establishment
                    for is_idx in range(num_species):
                        plot.seedling[is_idx] = (
                            site.species[is_idx].invader +
                            plot.seedling[is_idx] +
                            site.species[is_idx].sprout_num * plot.avail_spec[is_idx]
                        )
                
                # All trees die in disturbance
                plot.trees.clear()
                plot.numtrees = 0
                
                # Transfer all biomass to soil
                site.soil.C_into_A0 += zc
                site.soil.N_into_A0 += zn
                NPP_loss += zc
                NPPn_loss += zn
                
            else:
                # NO PLOT-LEVEL DISTURBANCE - Individual tree mortality
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
                        
                        # Track biomass of surviving trees
                        if site.species[k].conifer:
                            biomc += tree.biomC + leaf_bm
                            biomn += tree.biomN + leaf_bm / CON_LEAF_C_N
                        else:
                            biomc += tree.biomC
                            biomn += tree.biomN
                    
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
        
        # Apply unit conversions and update soil totals
        uconvert = HEC_TO_M2 / params.plotsize / site.numplots
        site.soil.biomc = biomc * uconvert
        site.soil.biomn = biomn * uconvert
        site.soil.net_prim_prodC = site.soil.net_prim_prodC - NPP_loss * uconvert
        site.soil.net_prim_prodN = site.soil.net_prim_prodN - NPPn_loss * uconvert
    
    def renewal(self, site: SiteData):
        """Process tree renewal - complete Fortran translation."""
        from .constants import HEC_TO_M2, STEM_C_N, CON_LEAF_C_N, DEC_LEAF_C_N, STD_HT
        from .random_utils import nrand
        
        # Early return if no available nitrogen
        if site.soil.avail_N <= 0.0:
            return
        
        growth_thresh = 0.05
        epsilon = 1e-10
        num_species = len(site.species)
        
        for plot in site.plots:
            # Check if recruitment can occur
            can_recruit = (plot.numtrees != 0 or (plot.wind == 0 and plot.fire == 0))
            
            if can_recruit:
                # Step 1: Calculate growth capacity and regrowth for all species
                regrowth = np.zeros(num_species)
                growmax = 0.0
                
                for i, species in enumerate(plot.species):
                    # Growth capacity from environmental factors
                    grow_cap = (species.fc_degday * species.fc_drought * 
                               species.fc_flood * plot.nutrient[i])
                    
                    # Light response for regrowth calculation
                    if species.conifer:
                        light_factor = species.light_rsp(plot.con_light[0])
                    else:
                        light_factor = species.light_rsp(plot.dec_light[0])
                    
                    regrowth[i] = grow_cap * light_factor
                    
                    # Apply growth threshold filter
                    if regrowth[i] <= growth_thresh:
                        regrowth[i] = 0.0
                    
                    growmax = max(growmax, regrowth[i])
                
                # Step 2: Determine maximum recruitment number
                max_renew = min(int(params.plotsize * growmax) - plot.numtrees,
                               int(params.plotsize * 0.5))
                nrenew = min(max(max_renew, 3), int(params.plotsize) - plot.numtrees)
                nrenew = max(0, min(nrenew, params.maxtrees - plot.numtrees))

                # Check if this is first recruitment cycle for this plot
                # The seedling calculation differs based on seedling_number
                if plot.seedling_number == 0.0:
                    # Step 3a: First recruitment - update seedbank and seedling populations
                    for i, species in enumerate(plot.species):
                        # Update seedbank
                        plot.seedbank[i] += (
                            species.invader +
                            species.seed_num * plot.avail_spec[i] +
                            species.sprout_num * plot.avail_spec[i]
                        )

                        # Transfer from seedbank to seedling if growth conditions met
                        if regrowth[i] >= growth_thresh:
                            plot.seedling[i] += plot.seedbank[i]
                            plot.seedbank[i] = 0.0
                        else:
                            # Seedbank mortality
                            plot.seedbank[i] *= species.seed_surv

                        # Apply fire response (if fire occurred)
                        if plot.fire > 0:
                            plot.seedling[i] = (plot.seedling[i] +
                                              species.sprout_num * plot.avail_spec[i]) * species.fc_fire

                        # Update seedling_number to track maximum
                        plot.seedling_number = max(int(plot.seedling[i]), plot.seedling_number)

                        # Scale to plot size
                        plot.seedling[i] *= params.plotsize

                    # Step 4: Calculate species selection probabilities
                    prob = np.zeros(num_species)
                    probsum = 0.0

                    for i in range(num_species):
                        prob[i] = plot.seedling[i] * regrowth[i]
                        probsum += prob[i]

                else:
                    # Step 3b: Subsequent recruitments - use existing seedling populations
                    # Calculate probabilities first (before scaling)
                    prob = np.zeros(num_species)
                    probsum = 0.0

                    for i in range(num_species):
                        prob[i] = plot.seedling[i] * regrowth[i]
                        probsum += prob[i]

                    # Now update seedbank and seedling for next cycle
                    for i, species in enumerate(plot.species):
                        # Update seedbank
                        plot.seedbank[i] += (
                            species.invader +
                            species.seed_num * plot.avail_spec[i] +
                            species.sprout_num * plot.avail_spec[i]
                        )

                        # Transfer from seedbank to seedling if minimum growth met
                        growth_min = 0.0  # Minimum growth threshold
                        if regrowth[i] >= growth_min:
                            plot.seedling[i] += plot.seedbank[i]
                            plot.seedbank[i] = 0.0
                        else:
                            # Seedbank mortality
                            plot.seedbank[i] *= species.seed_surv

                        # Apply fire response
                        if plot.fire > 0:
                            plot.seedling[i] = (plot.seedling[i] +
                                              species.sprout_num * plot.avail_spec[i]) * species.fc_fire

                        # Scale to plot size
                        plot.seedling[i] *= params.plotsize

                        # Update seedling_number
                        plot.seedling_number = max(int(plot.seedling[i]), plot.seedling_number)
                
                # Proceed with recruitment if probabilities exist
                if probsum > epsilon and nrenew > 0:
                    # Normalize and create cumulative probabilities
                    for i in range(num_species):
                        prob[i] /= probsum
                    
                    for i in range(1, num_species):
                        prob[i] += prob[i-1]
                    
                    # Step 5: Recruit new trees
                    for _ in range(nrenew):
                        # Select species by probability
                        q0 = urand()
                        selected_species = 0
                        
                        while selected_species < num_species - 1 and q0 > prob[selected_species]:
                            selected_species += 1
                        
                        # Ensure valid selection
                        if plot.seedling[selected_species] > 0:
                            # Decrement seedling count
                            plot.seedling[selected_species] -= 1.0
                            
                            # Create new tree
                            new_tree = TreeData()
                            new_tree.initialize_tree(plot.species[selected_species], selected_species)
                            
                            # Initialize size with normal distribution
                            z = 1.5 + nrand(0.0, 1.0)  # Normal random with mean=0, std=1
                            z = max(0.5, min(2.5, z))   # Clamp between 0.5 and 2.5 cm
                            
                            new_tree.diam_bht = z
                            new_tree.canopy_ht = 1.0
                            new_tree.forska_ht = 1.0
                            new_tree.mort_marker = False
                            
                            # Calculate derived attributes
                            new_tree.forska_height()
                            new_tree.stem_shape()
                            new_tree.biomass_c()
                            new_tree.biomass_n()
                            new_tree.leaf_biomass_c()
                            
                            # Add to plot
                            plot.trees.append(new_tree)
                            plot.numtrees = len(plot.trees)
                            
                            # Update carbon and nitrogen pools
                            species = plot.species[selected_species]
                            if species.conifer:
                                site.soil.net_prim_prodC += new_tree.biomC + new_tree.leaf_bm * (1.0 + CON_LEAF_RATIO)
                                site.soil.net_prim_prodN += (new_tree.biomC / STEM_C_N + 
                                                           new_tree.leaf_bm / CON_LEAF_C_N * (1.0 + CON_LEAF_RATIO))
                            else:
                                site.soil.net_prim_prodC += new_tree.biomC + new_tree.leaf_bm
                                site.soil.net_prim_prodN += new_tree.biomC / STEM_C_N + new_tree.leaf_bm / DEC_LEAF_C_N
            
            # Step 6: Scale seedlings back down (reverse the plotsize multiplication)
            # and apply seedling mortality for next cycle
            for i in range(num_species):
                # Scale back to per-m² units by dividing by plotsize
                # and apply seedling survival factor
                plot.seedling[i] = (plot.seedling[i] * plot.species[i].seedling_lg) / params.plotsize

            # Reset fire indicator
            plot.fire = 0
    
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