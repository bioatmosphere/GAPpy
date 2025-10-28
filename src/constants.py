"""
Constants module for UVAFME vegetation model.
Translated from Constants.f90
"""

import math

# Code identifiers
CODENAME = 'UVAFME'
VERSION_ID = '2012'

# Global constants
PI = 4.0 * math.atan(1.0)
DEG2RAD = 0.017453  # Temporary to match with old version

# Characters and files
MAX_NLEN = 30
MAX_FILE = 132
MAX_LINE = 256
MAX_LONG_LINE = 1000
MAX_DIR = 132
MAX_CHAR = 80
MAX_FIELDS = 100

# Standard for height measurements
STD_HT = 1.3

# Unit conversions
M_TO_CM = 100.0
HEC_TO_M2 = 10000
M2_TO_HEC = 0.0001
MM_TO_CM = 0.1

# Radiation constants
# Day length model parameters
B = 0.017214
AS = 0.409
AC = 0.033
PHASE = -1.39

# Radiation latitude dependencies
AMP = 37.58603
DL_OMEGA = 7.639437
EXRAD_COEF = 0.0820

# Hargreaves evaporation constants
H_COEFF = 0.000093876
H_ADDON = 17.8

# Climate related constants
NTEMPS = 12  # Number of temperature/precipitation values (monthly)

# Tree diameter categories
NHC = 7  # Number of height categories

# Values for invalid/missing data
RNVALID = -999.0  # Real invalid value
INVALID = -999    # Integer invalid value

MAX_DAYS_PER_YEAR = 366
DAYS_PER_YEAR = 365

# Precipitation nitrogen content
PRCP_N = 0.00002

# Global tree attributes
# Number of height categories
NHC = 7

# Conifer leaf C/N ratio
CON_LEAF_C_N = 60.0

# Deciduous leaf C/N ratio
DEC_LEAF_C_N = 40.0

# Stem C/N ratio
STEM_C_N = 450.0

# Conifer to deciduous leaf area ratio
CON_LEAF_RATIO = 0.3