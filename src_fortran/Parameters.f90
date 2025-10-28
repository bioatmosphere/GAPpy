module Parameters
  
  use Constants

  implicit none

  !:............................................................................:
  integer             :: clim_counter=0, rand_counter=0
  ! basic parameters
  integer             :: numyears
  integer             :: numplots
  integer             :: maxtrees, maxheight

  !Values for invalid/missing data
  real                :: rnvalid=-999.0
  integer             :: invalid=-999

  !Variables determining whether to use explicit seeds for RGNs
  logical             :: fixed_seed, same_climate
  logical             :: debug

  !Spinup
  logical             :: spinup
  integer             :: spinup_yrs
  !:............................................................................:

  ! Climate change variables
  real                ::incr_tmin_by
  real                ::incr_tmax_by  
  real                ::incr_precip_by 
  real                ::decr_tmin_by
  real                ::decr_tmax_by
  real                ::decr_precip_by
  real                ::tmin_change, tmax_change, precip_change
  integer             ::begin_change_year
  integer             ::start_gcm
  integer             ::end_gcm
  integer             ::duration_of_change
  character(len=4)    ::incr_or_decr
  integer             ::year_print_interval
  logical             ::linear_cc
  logical             ::with_clim_change,use_gcm
  logical             ::adjust_for_elev 
  logical             ::plot_level_data !Plot-level data, writes many rows
  logical             ::tree_level_data !DANGER !DANGER writes huge files

  !:............................................................................:

  ! plot parameters
  real                ::plotsize, rootdepth

  !  variables for changes to all site values
  real                ::  new_slope
  real                ::  fire_level, wind_level
  real                ::  SA_field_cap
  real                ::  A0_level_C, A0_level_N


end module Parameters
