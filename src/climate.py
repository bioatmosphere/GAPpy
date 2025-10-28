"""
Climate module for UVAFME vegetation model.
Translated from Climate.f90
"""

import math
import numpy as np
from .constants import *

# Global climate variables
accumulated_tmin = 0.0
accumulated_tmax = 0.0
accumulated_precip = np.zeros(12)

def set_site_climate(same_climate, fixed_seed):
    """Set climate for a site."""
    global accumulated_tmin, accumulated_tmax, accumulated_precip
    
    # This would call set_climate_rng_seed in the original
    # For now, we'll implement basic functionality
    accumulated_tmin = 0.0
    accumulated_tmax = 0.0
    accumulated_precip = np.zeros(12)


def cov365(ta1):
    """Convert monthly state data into daily data."""
    ta = np.zeros(13)
    vta = np.zeros(381)
    
    ltmt = [16, 45, 75, 105, 136, 166, 196, 227, 258, 288, 319, 349, 381]
    
    ta[12] = ta1[0]
    ta[0:12] = ta1[0:12]
    
    for k in range(12):
        tyxd = (ta[k+1] - ta[k]) / float(ltmt[k+1] - ltmt[k])
        for md in range(ltmt[k], ltmt[k+1]):
            vta[md] = ta[k] + tyxd * float(md - ltmt[k])
    
    result = np.zeros(365)
    # Fortran: do md=16,365 (1-based) → Python: range(15, 365) (0-based, indices 15-364)
    for md in range(15, 365):
        result[md] = vta[md]

    # Fortran: do md=1,15 (1-based) → Python: range(0, 15) (0-based, indices 0-14)
    for md in range(0, 15):
        result[md] = vta[365 + md]
    
    return result


def cov365a(ta, uniform_func):
    """Convert monthly integrated data into daily data randomly."""
    ltmt = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    vta = np.zeros(365)
    
    md = 0
    for k in range(12):
        raindays = min(25.0, ta[k] / 4.0 + 1.0)
        ik = int(raindays)
        rr = ta[k] / float(ik)
        ss = raindays / ltmt[k]
        inum = ik
        
        for i in range(ltmt[k]):
            yxdr = uniform_func()
            if inum > 0:
                if yxdr <= ss:
                    vta[md] = rr
                    inum -= 1
                else:
                    vta[md] = 0.0
            else:
                vta[md] = 0.0
            md += 1
        
        if inum > 0:
            vta[md - 16] = float(inum) * rr
    
    return vta


def ex_rad(julia, latit):
    """Calculate extraterrestrial radiation."""
    rlat = DEG2RAD * latit
    
    dr = 1.0 + AC * math.cos(B * float(julia))
    dairta = AS * math.sin(B * float(julia) + PHASE)
    yxd = -math.tan(rlat) * math.tan(dairta)
    
    if yxd >= 1.0:
        omega = 0.0
    elif yxd <= -1.0:
        omega = PI
    else:
        omega = math.acos(yxd)
    
    erad = AMP * math.cos(rlat) * math.cos(dairta) * (math.sin(omega) - omega * math.cos(omega))
    daylength = DL_OMEGA * omega
    exradmx = EXRAD_COEF * dr * math.cos(rlat - dairta)
    
    return erad, daylength, exradmx


def hargrea(tmin, tmax, ta, exrad):
    """Hargreaves evaporation formulation."""
    if ta <= 0.0:
        return 0.0
    else:
        return H_COEFF * (tmax - tmin)**0.5 * (ta + H_ADDON) * exrad