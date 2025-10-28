module Sitelist
  use Constants
  use Site
  use GenusGroups
  use Input
  implicit none

  integer                                :: numsites

  ! At some point we may want to have a sitelist type that attaches genera
  ! and species name lists.  Right now we'll stick with just an array of sites.

contains

  subroutine initialize_sitelist(sites,species_pres)
    type(SiteData), dimension(:), allocatable, intent(inout) :: sites
    type(Groups),                              intent(out)   :: species_pres

    integer, allocatable, dimension(:)   :: site_ids

    integer                              :: j, k
    integer                              :: sx
    integer                              :: cerr, serr

    call read_sitelist(site_ids)
    numsites=size(site_ids)
    allocate(sites(numsites))
    call read_sites(site_ids,sites)
    call read_climate(sites)
    call read_climate_stds(sites)

    ! Make adjustments
    do sx=1, numsites

       if ( sites(sx)%site_id .ne. invalid ) then

          ! Initialize altitude to "no adjustment" value
          sites(sx)%altitude=rnvalid

          ! Adjust according to global parameters set by user

          if ( new_slope .ne. rnvalid)  sites(sx)%slope=new_slope
          if ( fire_level .ne. rnvalid) sites(sx)%fire_prob=fire_level
          if ( wind_level .ne. rnvalid) sites(sx)%wind_prob=wind_level
          if ( SA_field_cap .ne. rnvalid) then
               sites(sx)%soil%A_field_cap=SA_field_cap
          endif
          if ( A0_level_C .ne. rnvalid) sites(sx)%soil%A0_c0=A0_level_C
          if ( A0_level_N .ne. rnvalid) sites(sx)%soil%A0_n0=A0_level_N

          ! standard adjustments
          sites(sx)%fire_prob       = sites(sx)%fire_prob / 1000.0
          sites(sx)%wind_prob       = sites(sx)%wind_prob / 1000.0
          sites(sx)%soil%A_w0       = sites(sx)%soil%A_w0 * rootDepth
          sites(sx)%soil%A_field_cap= sites(sx)%soil%A_field_cap*rootDepth
          sites(sx)%soil%A_perm_wp  = sites(sx)%soil%A_perm_wp*rootDepth
          sites(sx)%leaf_area_ind   = 1.0

       endif
       
    end do

    ! Read altitude list
    call read_altitudes(sites)

    ! Add species list
    call initialize_spec_list(sites)

    ! Get the lists of genera and species present
    call initialize_genus_groups(species_pres,sites)

  end subroutine initialize_sitelist


  subroutine initialize_spec_list(sites)
    type(SiteData),    dimension(:), intent(inout) :: sites

    type(SpeciesData), dimension(:),  allocatable  :: species_data
    integer,           dimension(:),  allocatable  :: site_ids
    integer,           dimension(:),  allocatable  :: range_site_ids
    character(len=8),  dimension(:,:),allocatable  :: range_species_ids
    integer                                        :: num_all_species
    integer                                        :: num_range_sites
    integer                                        :: l, n, ns
    logical                                        :: use_rangelist
    integer                                        :: gcstat
    logical                                        :: is_open
    logical                                        :: found_site

    allocate(site_ids(size(sites)))
    site_ids=sites(:)%site_id

    if (allocated(species_data)) deallocate(species_data)
    call read_speciesdata(site_ids,species_data)

    numsites=size(sites)

    call read_rangelist(use_rangelist,numsites,range_site_ids,range_species_ids)

    num_range_sites=size(range_site_ids)

    if ( use_rangelist ) then
       ! Add range information to Site structures
       do n=1,numsites
          found_site=.false.
          do l=1,num_range_sites
            if (sites(n)%site_id .eq. range_site_ids(l)) then
               ! Construct the list of species present
               found_site=.true.
               call attach_species(sites(n),species_data,range_species_ids(l,:))
               exit
            endif
          enddo
          ! If site n wasn't in the rangelist, add all species
          if ( .not. found_site ) then
             call attach_species(sites(n),species_data)
          endif
       enddo

    else

       !All species present in all sites
       do n=1,numsites
          call attach_species(sites(n),species_data)
       enddo

    endif

  end subroutine initialize_spec_list


end module Sitelist
