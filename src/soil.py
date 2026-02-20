"""
Soil module for GAPpy vegetation model.
Translated from Soil.f90
"""

import numpy as np
from .constants import *

# Global constants
AO_CN_0 = 30.0
SA_CN_0 = 4.0
SB_CN_0 = 20.0
AO_RESP = 5.24e-4
SA_RESP = 1.24e-5
SB_RESP = 2.74e-7

SOIL_BASE_DEPTH = 70.0
BASE_MAX = 0.6
BASE_MIN = 0.1
AO_MIN = 0.025
AO_MAX = 0.25
LAI_MIN = 0.01
LAI_MAX = 0.15


class SoilData:
    """Soil data structure containing soil properties and state variables."""
    
    def __init__(self):
        self.A0_c0 = 0.0
        self.A_c0 = 0.0
        self.A0_n0 = 0.0
        self.A_n0 = 0.0
        self.A0_w0 = 0.0
        self.A_w0 = 0.0
        self.A_field_cap = 0.0
        self.A_perm_wp = 0.0
        self.C_into_A0 = 0.0
        self.N_into_A0 = 0.0
        self.net_C_into_A0 = 0.0
        self.net_N_into_A0 = 0.0
        self.N_used = 0.0
        self.avail_N = 0.0
        self.BL_c0 = 0.0
        self.BL_n0 = 0.0
        self.BL_w0 = 0.0
        self.base_h = 0.0
        self.biomC = 0.0
        self.biomN = 0.0
        self.net_prim_prodN = 0.0
        self.net_prim_prodC = 0.0
        self.runoff = 0.0
        self.total_C_rsp = 0.0
        self.new_growth = 0

    def soil_decomp(self, litter_c1, litter_c2, litter_n1, litter_n2, 
                    tempC, precip, aow0_scaled_by_max, saw0_scaled_by_fc, 
                    sbw0_scaled_by_max):
        """
        Soil decomposition model for computing available N and soil respiration.
        
        Parameters:
        - litter_c1/2: input C as litter (above/under ground): tc/ha/d
        - litter_n1/2: input N as litter (above/under ground): tn/ha/d
        - precip: input water as rain fall cm/d
        - tempC: daily temperature degree C
        - aow0_scaled_by_max: aow0/aowmax cm/cm
        - saw0_scaled_by_fc: saw0/safc cm/cm
        - sbw0_scaled_by_max: sbw0/sbwmax cm/cm (never used)
        
        Returns:
        - avail_N: available N for plant growth tn/ha
        - C_resp: emission to atmosphere as CO2 tc/ha/d
        - aow0_scaled_by_max: clamped value (Fortran intent(out), propagates back)
        """
        
        # Copy from object to local variables
        ao_c0 = self.A0_c0
        ao_n0 = self.A0_n0
        sa_c0 = self.A_c0
        sa_n0 = self.A_n0
        sb_c0 = self.BL_c0
        sb_n0 = self.BL_n0
        
        # Ao layer C N balance
        ao_c0 = ao_c0 + litter_c1
        ao_n0 = ao_n0 + litter_n1
        ao_cn = ao_c0 / ao_n0
        aow0_scaled_by_max = min(aow0_scaled_by_max, 0.5)
        
        aofunc = max((1.0 - (1.0 - aow0_scaled_by_max / 0.3)**2), 0.2)
        
        if tempC >= -5.0:
            tadjst = 3.0**(0.1 * (tempC - 1.0))
            tadjst1 = 2.5**(0.1 * (tempC - 1.0))
        else:
            tadjst = 0.0
            tadjst1 = 0.0
        
        resp1 = tadjst * aofunc * AO_RESP * ao_c0
        yxdn = resp1 / ao_cn
        yxdc = yxdn * AO_CN_0
        ao_c0 = ao_c0 - yxdc - resp1
        ao_n0 = ao_n0 - yxdn
        
        # Soil A layer C N balance
        sa_c0 = sa_c0 + yxdc + litter_c2
        sa_n0 = sa_n0 + yxdn + litter_n2
        sa_cn = sa_c0 / sa_n0
        
        resp2 = tadjst1 * max(1.0 - (1.0 - saw0_scaled_by_fc / 0.8)**2, 0.2) * SA_RESP * sa_c0
        
        tosb = resp2 / SB_CN_0
        avail_N = resp2 / sa_cn * max(0.5, (sa_cn - SA_CN_0) / sa_cn)
        sa_c0 = sa_c0 - resp2 - tosb
        sa_n0 = sa_n0 - avail_N
        
        # Base layer C balance
        sb_c0 = sb_c0 + tosb
        resp3 = sb_c0 * SB_RESP * tadjst1
        sb_c0 = sb_c0 - resp3
        
        # Total respiration of the 'Three' Layer
        C_resp = resp1 + resp2 + resp3
        
        # Update object state
        self.A0_c0 = ao_c0
        self.A0_n0 = ao_n0
        self.A_c0 = sa_c0
        self.A_n0 = sa_n0
        self.BL_c0 = sb_c0
        self.BL_n0 = sb_n0
        self.avail_N = avail_N

        # Return clamped aow0_scaled_by_max to match Fortran intent(out) semantics
        return avail_N, C_resp, aow0_scaled_by_max

    def soil_water(self, slope, lai, lai_w0, sigma, freeze, rain, pot_ev_day):
        """
        Soil water daily cycle model.

        Parameters:
        - slope: slope degree
        - lai: canopy leaf area index m/m
        - lai_w0: initial canopy water content
        - sigma: sigma parameter
        - freeze: freeze factor
        - rain: daily precipitation cm/day
        - pot_ev_day: daily potential evapotranspiration cm/day

        Returns (tuple of 11 values):
        - act_ev_day: daily actual evapotranspiration cm/day
        - laiw0_scaled_by_max, laiw0_scaled_by_min: canopy water scaling factors
        - aow0_scaled_by_max, aow0_scaled_by_min: A0 layer water scaling factors
        - sbw0_scaled_by_max, sbw0_scaled_by_min: base layer water scaling factors
        - saw0_scaled_by_fc, saw0_scaled_by_wp: A layer water scaling factors
        - lai_w0: updated canopy water content (must be persisted by caller!)
        - lai: clamped leaf area index (Fortran intent(inout), propagates back)
        """
        
        # Copy from object to local variables
        ao = self.A0_c0
        ao_w0 = self.A0_w0
        sa_w0 = self.A_w0
        sa_fc = self.A_field_cap
        sa_pwp = self.A_perm_wp
        sb_w0 = self.BL_w0
        runoff = self.runoff
        
        act_ev_day = 0.0
        lai = max(lai, 1.0)
        sbh = SOIL_BASE_DEPTH
        
        laiw_min = lai * LAI_MIN
        laiw_max = lai * LAI_MAX
        aow_min = ao * AO_MIN
        aow_max = ao * AO_MAX
        sbw_min = sbh * BASE_MIN
        sbw_max = sbh * BASE_MAX
        
        # Forest region water balance for underground water table
        if rain > 0.01:
            table_water = rain * sigma * freeze
            sb_w0 = min(sbw_max, sb_w0 + table_water)
        
        if pot_ev_day <= 0.0:
            laiw = min(rain + lai_w0, laiw_max)
            yxd1 = max(rain - laiw + lai_w0, 0.0)
            aow = min(ao_w0 + yxd1, aow_max)
            yxd2 = max(yxd1 - aow + ao_w0, 0.0)
            sbw = min(yxd2 + sb_w0, sbw_max)
            runoff = max(yxd2 - sbw + sb_w0, 0.0)
            lai_w0 = laiw
            ao_w0 = aow
            sb_w0 = sbw
        else:
            # Modified by BW on 2017-07-13
            lai_loss = min((laiw_max - lai_w0), rain)
            
            yxd1 = max(rain - lai_loss, 0.0)
            laiw = lai_w0 + lai_loss
            yxd = (slope / 90.0)**2
            lossslp = yxd * yxd1
            yxd2 = yxd1 - lossslp - pot_ev_day
            
            if yxd2 > 0.0:
                saw = min(sa_fc - sa_w0, yxd2) + sa_w0
                yxd3 = max(yxd2 - (saw - sa_w0), 0.0)
                aow = min(aow_max - ao_w0, yxd3) + ao_w0
                yxd4 = max(yxd3 - (aow - ao_w0), 0.0)
                sbw = min(sbw_max - sb_w0, yxd4) + sb_w0
                act_ev_day = pot_ev_day
                runoff = yxd2 + lossslp
                sa_w0 = saw
                sb_w0 = sbw
                ao_w0 = aow
            else:
                lai_w1 = min(-yxd2, lai_w0 - laiw_min)
                lai_w0 = lai_w0 - lai_w1
                act_ev_day = act_ev_day + lai_w1
                
                yxd3 = min(yxd2 + lai_w1, 0.0)
                ao_w1 = min(-yxd3, ao_w0 - aow_min)
                act_ev_day = act_ev_day + ao_w1
                ao_w0 = ao_w0 - ao_w1
                yxd4 = min(yxd3 + ao_w1, 0.0)
                sa_w1 = min(-yxd4, sa_w0 - sa_pwp)
                act_ev_day = act_ev_day + sa_w1
                sa_w0 = sa_w0 - sa_w1
                yxd5 = min(yxd4 + sa_w1, 0.0)
                sb_w1 = min(-yxd5, sb_w0 - sbw_min)
                sb_w0 = sb_w0 - sb_w1
                runoff = lossslp
        
        # Calculate scaling factors
        laiw0_scaled_by_max = lai_w0 / laiw_max
        laiw0_scaled_by_min = lai_w0 / laiw_min
        aow0_scaled_by_max = ao_w0 / aow_max
        aow0_scaled_by_min = ao_w0 / aow_min
        sbw0_scaled_by_max = sb_w0 / sbw_max
        sbw0_scaled_by_min = sb_w0 / sbw_min
        saw0_scaled_by_fc = sa_w0 / sa_fc
        saw0_scaled_by_wp = sa_w0 / sa_pwp
        
        runoff = max(runoff, 0.0)
        
        # Update object state
        self.A0_c0 = ao
        self.A0_w0 = ao_w0
        self.A_w0 = sa_w0
        self.A_field_cap = sa_fc
        self.A_perm_wp = sa_pwp
        self.BL_w0 = sb_w0
        self.base_h = sbh
        self.runoff = runoff
        
        return (act_ev_day, laiw0_scaled_by_max, laiw0_scaled_by_min,
                aow0_scaled_by_max, aow0_scaled_by_min, sbw0_scaled_by_max,
                sbw0_scaled_by_min, saw0_scaled_by_fc, saw0_scaled_by_wp,
                lai_w0, lai)