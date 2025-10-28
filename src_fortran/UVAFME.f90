!:............................................................................:

program UVAFE

  use Constants
  use Parameters
  use Input
  use Output
  use Species
  use Site
  use Sitelist
  use GenusGroups
  use Model

  !:............................................................................:

  integer               :: sndx, year ! site and year loop indexes
  real                  :: start_time,total_time

  character(len=80)     :: filelist
  integer               :: nargs
  integer               :: print_interval

  type(SiteData),    dimension(:), allocatable :: sites
  type(Groups)                                 :: species_present

  interface
     subroutine drawBanner(numsites,species_present)
       use Parameters
       use Genusgroups
       implicit none
       type(Groups)   ::  species_present
       integer        ::  numsites
     end subroutine drawBanner

     subroutine showProgress(asite)
       use Constants
       use Parameters
       use Site
       implicit none
       type(SiteData), intent(in) :: asite
     end subroutine showProgress

  end interface

  !:............................................................................:

  ! Fortran 2003
  nargs = command_argument_count()
  if (nargs .ne. 1) then
    filelist=''
  else
    call get_command_argument(1,filelist)
  endif

  ! Prepare input files
  call initialize_inputFiles(filelist)

  ! Prepare site and species data
  call initialize_sitelist(sites,species_present)

  ! Prepare output files
  call initialize_outputFiles(species_present)

  ! write runtime vars to screen
  call drawBanner(size(sites),species_present)

  ! Start timing
  call cpu_time(start_time)

  do sndx = 1, numsites

     ! make sure the site exists
     if ( sites(sndx)%site_wmo == rnvalid ) then
          write(*,*) '             No site or climate file for site ',         &
                                   sites(sndx)%site_id
          write(*,*) '             Skipping site ',sites(sndx)%site_name
          write(*,*) 
          cycle
     endif

     ! load climate and site specific vars,then adjust for altitude if requested
     call set_site_rng_seed(fixed_seed)
     call set_site_climate(same_climate,fixed_seed)

     ! skip this site if no species are present
     if ( size(sites(sndx)%species) .eq. 0 )  then
          write(*,*) '              No species present in site ',              &
                                    sites(sndx)%site_id
          write(*,*) '              Skipping site ',sites(sndx)%site_name
          write(*,*) '             '
          cycle
     endif
     call showProgress(sites(sndx))

     ! run the model
     do year = 0, numyears
        
        call BioGeoClimate(sites(sndx),year)
        ! We will write soil/CN/clim data for last year's trees but after BioGeo
        call write_soil_data(sites(sndx),year)
        call write_site_data(sites(sndx),year)
        ! Trees !
        call Canopy(sites(sndx))
        call Growth(sites(sndx))
        call Mortality(sites(sndx))
        call Renewal(sites(sndx))

        if ( spinup ) then
           if ( year < spinup_yrs ) then
               print_interval=10*year_print_interval
           else
               print_interval=year_print_interval
           endif
        else
           print_interval=year_print_interval
        endif

        if ((mod(year,print_interval).eq.0) .or. year==numyears) then
           call write_genus_data(sites(sndx),species_present,year)
           call write_species_data(sites(sndx),species_present,year)
           if (tree_level_data) then
              call write_tree_data(sites(sndx),year)
           endif
        end if

     end do

     call cpu_time(total_time)
     write(*,*) 'Cumulative time :',total_time-start_time
     write(*,'(a80)') &
     '=============================================================================='

  end do

  call close_outputFiles

end program UVAFE

!:............................................................................:

  subroutine drawBanner(numsites,species_present)
  
    use Parameters
    use Genusgroups
    implicit none

    type(Groups), intent(in)       ::  species_present
    integer,      intent(in)       ::  numsites

    write(*,500) &
    '=============================================================================='
    write(*,500) &
    '                       UVA Forest Model Enhanced                              '
    write(*,500) &
    '               Center For Regional Environmental Studies                      '
    write(*,500) &
    '                        University of Virginia                                '
    write(*,500) &
    '                   Department of Environmental Sciences                       '
    write(*,500) &
    '=============================================================================='

    write(*,*) 'Running with parameters:'
    write(*,400) 'Number of sites:',numsites
    write(*,400) 'Number of years:',numyears
    write(*,400) 'Number of species:',species_present%numspecies
    write(*,400) 'Maximum number of trees:',maxtrees
    write(*,400) 'Maximum height of trees:',maxheight
    write(*,401) 'Plotsize:',plotsize
    write(*,401) 'Root depth:',rootdepth
    if ( with_clim_change ) then
       write(*,*) 'Running with climate change'
       write(*,400) 'Beginning at year:',begin_change_year
       write(*,400) 'Duration in years:',duration_of_change
       if ( linear_cc ) then
          if ( incr_or_decr .eq. 'incr' ) then
             write(*,401) 'Total tmin increase',incr_tmin_by
             write(*,401) 'Total tmax increase',incr_tmax_by
             write(*,401) 'Total precip increase',incr_precip_by
          else if ( incr_or_decr .eq. 'decr' ) then
             write(*,401) 'Total tmin decrease',decr_tmin_by
             write(*,401) 'Total tmax decrease',decr_tmax_by
             write(*,401) 'Total precip decrease',decr_precip_by
          endif
       else if ( use_gcm ) then
          write(*,*) 'Using GCM data:'
          write(*,400) 'GCM start year ',start_gcm
          write(*,400) 'GCM end year '  ,end_gcm
       endif
    endif

    write(*,400) 'Printing interval in years:',year_print_interval
    write(*,500) &
    '=============================================================================='
    write(*,*)


    400 format(a30,i10)
    401 format(a30,f10.3)
    402 format(a30,a)
    500 format(a80)

  end subroutine drawbanner


  subroutine showProgress(asite)
    use Parameters
    use Site
    implicit none
    type(SiteData), intent(in)     ::  asite
    integer                        ::  num_site_species

    num_site_species=size(asite%species)
    write(*,500) 'Running for site ', asite%site_id,asite%site_name
    write(*,501) 'Number of species present ',num_site_species
    if ( asite%altitude .ne. rnvalid ) then
       write(*,*) '      Site altitude adjustment ',asite%altitude
    endif
    write(*,502) asite%elevation,asite%slope,                                  &
                 asite%fire_prob,asite%wind_prob,                              &
                 asite%soil%A_field_cap,asite%soil%A0_c0,asite%soil%A0_n0

    500 format(14x,a,i10,4x,a)
    501 format(14x,a,i8)
    502 format(7x,'Site parameters: elevation ',f8.2,'   slope     ',f6.2,/    &
               23x,' fire/1000   ' f6.2,'   wind/1000 ',f6.2, /                &
               23x,' SAFC*rtdpth ',f6.2,'   A0_C      ',f6.2,'  A0_N ',f6.2)

  end subroutine showProgress


!:............................................................................:

