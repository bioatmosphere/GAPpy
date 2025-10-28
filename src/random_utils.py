"""
Random number utilities for UVAFME vegetation model.
Translated from Random.f90

This module implements the exact random number generators from the original Fortran code,
including the Linear Congruential Generator (LCG) from Numerical Recipes for climate RNG.
"""

import random
import math
from typing import Optional


class ClimateRNG:
    """
    Climate-specific RNG using Linear Congruential Generator from Numerical Recipes.

    This implements the exact algorithm from Random.f90 clim_urand function (lines 104-170)
    to ensure reproducible results matching the Fortran implementation.
    """

    # Constants from Numerical Recipes (Fortran lines 112-114)
    M1 = 259200
    IA1 = 7141
    IC1 = 54773
    M2 = 134456
    IA2 = 8121
    IC2 = 28411
    M3 = 243000
    IA3 = 4561
    IC3 = 51349

    def __init__(self):
        self.ix1 = 0
        self.ix2 = 0
        self.ix3 = 0
        self.r = [0.0] * 97  # Shuffle table
        self.first = True
        self.reset_clim = True
        self.clim_seed = 0

    def seed(self, seed_val: int):
        """Set the seed for climate RNG."""
        self.clim_seed = seed_val
        self.reset_clim = True
        self.first = True

    def uniform(self, lb: float = 0.0, ub: float = 1.0) -> float:
        """
        Generate uniform random number using LCG with Bays-Durham shuffle.
        Exact translation of clim_urand from Random.f90 (lines 104-170).
        """
        rm1 = 1.0 / self.M1
        rm2 = 1.0 / self.M2

        # Initialize on first call or after reset (Fortran lines 140-153)
        if self.first or self.reset_clim:
            self.ix1 = abs(self.IC1 - self.clim_seed) % self.M1
            self.ix1 = (self.IA1 * self.ix1 + self.IC1) % self.M1
            self.ix2 = self.ix1 % self.M2
            self.ix1 = (self.IA1 * self.ix1 + self.IC1) % self.M1
            self.ix3 = self.ix1 % self.M3

            # Fill shuffle table (Fortran lines 146-150)
            for j in range(97):
                self.ix1 = (self.IA1 * self.ix1 + self.IC1) % self.M1
                self.ix2 = (self.IA2 * self.ix2 + self.IC2) % self.M2
                self.r[j] = (float(self.ix1) + float(self.ix2) * rm2) * rm1

            self.first = False
            self.reset_clim = False

        # Generate random number (Fortran lines 155-165)
        self.ix1 = (self.IA1 * self.ix1 + self.IC1) % self.M1
        self.ix2 = (self.IA2 * self.ix2 + self.IC2) % self.M2
        self.ix3 = (self.IA3 * self.ix3 + self.IC3) % self.M3

        j = (97 * self.ix3) // self.M3  # Fortran line 158: j=1+(97*ix3)/m3

        # Bounds check (Fortran lines 159-165)
        if j > 96 or j < 0:  # Python 0-based: 0-96 instead of Fortran 1-97
            j = 96
            print('Warning: Error in climate RNG bounds')
            ran1 = self.r[96]
        else:
            ran1 = self.r[j]
            self.r[j] = (float(self.ix1) + float(self.ix2) * rm2) * rm1

        # Scale to desired range (Fortran line 167)
        return lb + (ub - lb) * ran1


class RandomState:
    """Manages random number generator state for UVAFME."""

    def __init__(self):
        self.site_rng = random.Random()
        self.climate_rng = ClimateRNG()  # Use custom LCG for climate
        self.fixed_seed = False
        self.same_climate = True
        
    def set_site_rng_seed(self, fixed_seed: bool = False, seed: Optional[int] = None):
        """Set random number generator seed for site-specific calculations."""
        self.fixed_seed = fixed_seed
        
        if fixed_seed:
            if seed is None:
                seed = 42  # Default fixed seed
            self.site_rng.seed(seed)
        else:
            self.site_rng.seed()  # Use system time
    
    def set_climate_rng_seed(self, same_climate: bool = True, fixed_seed: bool = False, seed: Optional[int] = None):
        """
        Set random number generator seed for climate calculations.
        Translated from Random.f90 set_climate_rng_seed (lines 273-326).
        """
        self.same_climate = same_climate

        if fixed_seed:
            if seed is None:
                seed = 2345678  # Default seed from Fortran (line 278)
            self.climate_rng.seed(seed)
            if same_climate:
                self.climate_rng.reset_clim = True
        else:
            if seed is None:
                import time
                seed = int(time.time() * 1000) % 1000000  # Use milliseconds
            self.climate_rng.seed(seed)
            if same_climate:
                self.climate_rng.reset_clim = True

    def urand(self, lb: float = 0.0, ub: float = 1.0, seed: Optional[int] = None) -> float:
        """
        Generate uniform random number between lb and ub.
        Translated from Random.f90 urand (lines 20-53).
        """
        if seed is not None:
            temp_rng = random.Random(seed)
            return temp_rng.uniform(lb, ub)
        else:
            return self.site_rng.uniform(lb, ub)

    def nrand(self, mean: float = 0.0, std: float = 1.0, seed: Optional[int] = None) -> float:
        """
        Generate normal random number using Box-Muller polar method.
        Translated from Random.f90 nrand (lines 56-101).
        """
        # Box-Muller polar method (Fortran lines 79-95)
        while True:
            if seed is not None:
                x1 = self.urand(-1.0, 1.0, seed)
                x2 = self.urand(-1.0, 1.0, seed)
            else:
                x1 = self.urand(-1.0, 1.0)
                x2 = self.urand(-1.0, 1.0)

            w = x1**2 + x2**2
            if w != 0.0 and w <= 1.0:
                break

        w = math.sqrt((-2.0 * math.log(w)) / w)
        y1 = x1 * w
        # y2 = x2 * w  # Not used, but available

        return y1 * std + mean

    def clim_urand(self, lb: float = 0.0, ub: float = 1.0) -> float:
        """
        Generate uniform random number for climate calculations using LCG.
        Translated from Random.f90 clim_urand (lines 104-170).
        """
        return self.climate_rng.uniform(lb, ub)

    def clim_nrand(self, mean: float = 0.0, std: float = 1.0) -> float:
        """
        Generate normal random number for climate using Box-Muller.
        Translated from Random.f90 clim_nrand (lines 173-211).
        """
        # Box-Muller polar method (Fortran lines 195-208)
        while True:
            x1 = self.clim_urand(-1.0, 1.0)
            x2 = self.clim_urand(-1.0, 1.0)

            w = x1**2 + x2**2
            if w != 0.0 and w <= 1.0:
                break

        w = math.sqrt((-2.0 * math.log(w)) / w)
        y1 = x1 * w
        # y2 = x2 * w  # Not used

        return y1 * std + mean
    
    def get_random_seed(self) -> int:
        """Get current random seed state."""
        return self.site_rng.getstate()[1][0]
    
    def set_random_seed(self, seed: int):
        """Set random seed state."""
        self.site_rng.seed(seed)


# Global random state instance
rng_state = RandomState()

# Convenience functions that use the global state
def urand(lb: float = 0.0, ub: float = 1.0, seed: Optional[int] = None) -> float:
    """Generate uniform random number between lb and ub."""
    return rng_state.urand(lb, ub, seed)

def nrand(mean: float = 0.0, std: float = 1.0, seed: Optional[int] = None) -> float:
    """Generate normal random number with given mean and standard deviation."""
    return rng_state.nrand(mean, std, seed)

def clim_urand(lb: float = 0.0, ub: float = 1.0) -> float:
    """Generate uniform random number for climate calculations."""
    return rng_state.clim_urand(lb, ub)

def clim_nrand(mean: float = 0.0, std: float = 1.0) -> float:
    """Generate normal random number for climate calculations."""
    return rng_state.clim_nrand(mean, std)

def set_site_rng_seed(fixed_seed: bool = False, seed: Optional[int] = None):
    """Set random number generator seed for site-specific calculations."""
    rng_state.set_site_rng_seed(fixed_seed, seed)

def set_climate_rng_seed(same_climate: bool = True, fixed_seed: bool = False, seed: Optional[int] = None):
    """Set random number generator seed for climate calculations."""
    rng_state.set_climate_rng_seed(same_climate, fixed_seed, seed)

def get_random_seed() -> int:
    """Get current random seed state."""
    return rng_state.get_random_seed()

def set_random_seed(seed: int):
    """Set random seed state."""
    rng_state.set_random_seed(seed)