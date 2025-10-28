module Constants
  ! Code identifiers
  character(len=6)      :: codename='UVAFME'
  character(len=4)      :: version_id='2012'

  ! Global constants
  real,    parameter    :: pi=4.0*atan(1.0)
  !  real,    parameter    :: deg2rad=pi/180.
  !KAH Temporary to match with old version
  real,    parameter    :: deg2rad=0.017453

  ! Characters and files
  integer, parameter    :: MAX_NLEN=30
  integer, parameter    :: MAX_FILE=132
  integer, parameter    :: MAX_LINE=256
  integer, parameter    :: MAX_LONG_LINE=1000
  integer, parameter    :: MAX_DIR =132
  integer, parameter    :: MAX_CHAR=80
  integer, parameter    :: MAX_FIELDS=100

  ! Standard for height measurements
  real, parameter       :: std_ht=1.3

  ! Unit conversions
  real, parameter       :: m_to_cm=100.0
  real, parameter       :: hec_to_m2=10000
  real, parameter       :: m2_to_hec=0.0001
  real, parameter       :: mm_to_cm=0.1

  ! Radiation constants
  ! Day length model parameters
  real, parameter       :: b=0.017214
  real, parameter       :: As=0.409
  real, parameter       :: Ac=0.033
  real, parameter       :: phase=-1.39
  ! Radiation latitude dependencies
  real, parameter       :: Amp=37.58603
  real, parameter       :: dl_omega=7.639437
  real, parameter       :: exrad_coef=0.0820

  ! Hargreaves evaporation constants
  real, parameter       :: H_coeff=0.000093876
  real, parameter       :: H_addon=17.8

  ! Climate related constants
  integer, parameter    ::  NTEMPS=12
  integer, parameter    ::  max_days_per_year=366
  integer, parameter    ::  days_per_year=365
  ! Precipitation nitrogen content
  real                  :: prcp_n=0.00002

  !  Global tree attributes
  !  Number of height categories
  integer, parameter   :: NHC=7
  !  Conifer leaf C/N ratio
  real                 :: con_leaf_c_n=60.0
  !  Deciduous leaf C/N ratio
  real                 :: dec_leaf_c_n=40.0
  !  Stem C/N ratio
  real                 :: stem_c_n=450.0
  !  Conifer to deciduous leaf area ratio ?
  real                 :: con_leaf_ratio=0.3

end module Constants
