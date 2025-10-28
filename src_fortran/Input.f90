module Input
  use Constants
  use Parameters
  use Site 
  use IO
  use FileUtils

  implicit none

contains

  subroutine initialize_inputFiles(filelist)
    character(len=*)            ::   filelist

    call open_inputFiles(filelist)
    call initialize_parameters

  end subroutine initialize_inputFiles


  subroutine initialize_parameters

    logical          ::   plot_level_data
    logical          ::   tree_level_data

    integer          ::   nunit, ios
    logical          ::   is_open, file_exists     

    namelist /uvafme/ numyears,numplots,plotsize,rootdepth,maxheight,maxtrees, &
         with_clim_change,use_gcm,start_gcm,end_gcm,                           &
         plot_level_data,tree_level_data,                                      &
         adjust_for_elev,incr_tmin_by,incr_tmax_by,                            &
         year_print_interval, incr_precip_by,decr_tmin_by,                     &
         decr_tmax_by,decr_precip_by, begin_change_year,                       &
         duration_of_change,incr_or_decr,same_climate,fixed_seed,              &
         spinup,spinup_yrs,                                                    &
         debug

    namelist /sitevals/ new_slope,fire_level,wind_level,                       &
                        SA_field_cap,A0_level_C,A0_level_N


    ! Default values

    debug=.false.

    spinup=.false.
    spinup_yrs=500

    year_print_interval=10
    begin_change_year =300
    duration_of_change= 0
    with_clim_change  = .false.
    plot_level_data   = .false.
    tree_level_data   = .false.
    adjust_for_elev   = .false.
    same_climate = .true.
    fixed_seed   = .false.
    use_gcm      = .false.
    start_gcm = 0
    end_gcm   = 100

    numyears=1000
    numplots=200
    maxtrees=1000
    maxheight=60

    incr_tmin_by=0.00
    incr_tmax_by=0.00
    incr_precip_by=0.0000
    decr_tmin_by=0.00
    decr_tmax_by=0.00
    decr_precip_by=0.0000

    plotsize=500.
    rootdepth=0.8

    incr_or_decr='incr'

    ! Read parameters namelist.  Any variables not set in this file will
    ! take the default values above.

    inquire(rt_file,exist=file_exists,opened=is_open)

    if ( rt_file .ne. 0 .and. file_exists .and. is_open ) then
       read(rt_file,uvafme,iostat=ios) 
    else
       write(*,*) 'Unable to open runtime (parameter) file'
       write(*,*) 'Using all default values'
    endif

    if ( with_clim_change ) then
       if ( .not. use_gcm ) then
          linear_cc=.true.
          if ( duration_of_change == 0 ) then
             stop "Inconsistent input for climate change: 0 years duration"
          else if ( incr_tmin_by == 0.0 .and. decr_tmax_by==0.0 .and. &
               decr_tmin_by == 0.0 .and. decr_tmax_by==0.0 .and. &
               incr_precip_by == 0.0 .and. decr_precip_by == 0.0 ) then
             stop "Inconsistent climate change: no temp or precip changes"
          else if ( incr_or_decr == 'incr' ) then
             tmin_change=incr_tmin_by/(duration_of_change+1)
             tmax_change=incr_tmax_by/(duration_of_change+1)
             precip_change=incr_precip_by/(duration_of_change+1)
          else if ( incr_or_decr == 'decr' ) then
             if ( decr_tmin_by < 0.0 ) then
                write(*,*) "Assuming decrease intended"
                decr_tmin_by=-decr_tmin_by
             endif
             if ( decr_tmax_by < 0.0 ) then
                write(*,*) "Assuming decrease intended"
                decr_tmax_by=-decr_tmax_by
             endif
             if ( decr_precip_by < 0.0 ) then
                write(*,*) "Assuming decrease intended"
                decr_precip_by=-decr_precip_by
             endif
             tmin_change   = -decr_tmin_by/(duration_of_change+1)
             tmax_change   = -decr_tmax_by/(duration_of_change+1)
             precip_change = -decr_precip_by/(duration_of_change+1)
          endif
       else
          linear_cc=.false.
          tmin_change=0.0
          tmax_change=0.0
          precip_change=0.0
          ! For consistency, it needs to change one way or the other
          ! Note that with the GCM we will not add changes to the old data,
          ! we will read in a new array for the tmin, tmax, and precip data
          if ( duration_of_change == 0 ) then
             stop "Inconsistent input for climate change: 0 years duration"
          else
             !We start counting at 0
             end_gcm=start_gcm+duration_of_change-1
          endif
       endif

    endif

    if ( debug ) then
       fixed_seed=.true.
       same_climate=.true.
    endif


    ! Set defaults
    new_slope=rnvalid
    fire_level=rnvalid
    wind_level=rnvalid
    SA_field_cap=rnvalid
    A0_level_C=rnvalid
    A0_level_N=rnvalid

    ! Read namelist for site parameters
    if ( rt_file .ne. 0 .and. file_exists .and. is_open ) then
       read(rt_file,sitevals,iostat=ios) 
    endif

    close(nunit)

  end subroutine initialize_parameters


  subroutine read_sitelist(site_ids)
    integer,allocatable,dimension(:), intent(out) ::  site_ids

    integer,allocatable,dimension(:)              ::  temp_sites
    integer                                       ::  numsites
    integer                                       ::  ios
    integer                                       ::  lc,l
    logical                                       ::  file_exists,is_open 

    inquire(slist,exist=file_exists,opened=is_open)
    if ( slist == 0 .or. (.not. file_exists) .or. (.not. is_open) ) then
       call fatal_error('Unable to find site list file')
    endif

    ! Sitelist has no header
    numsites=count_records(slist,0)
    if ( numsites <= 0 ) then
       call fatal_error("Error in site list file")
    else
       allocate(site_ids(numsites))
    endif

    ios = 0
    lc=0
    ! do while doesn't really work quite as we'd like (we want repeat until)
    ! Also note that some compilers return nonzero iostat values for harmless
    ! errors, such as format conversions, so we can't rely on testing for ios=0
    do while ( (lc .lt. numsites ) )
       lc = lc + 1
       read(slist,*,iostat=ios) site_ids(lc)
    enddo

    if ( lc .ne. numsites ) then
       write(*,*) 'Warning : Initial Sites Scan(',numsites,') and Sites Read(',&
            lc,') do not match'
       if ( lc .gt. 0 ) then
          ! This accounts for read errors that are recoverable
          allocate(temp_sites(lc))
          temp_sites=site_ids(1:lc)
          deallocate(site_ids)
          numsites=lc
          allocate(site_ids(numsites))
          site_ids=temp_sites
       else
          stop "No site data read"
       endif
    end if

    write(logf,*) 'Site data initialized. Total read in: ', numsites

  end subroutine read_sitelist


  subroutine read_altitudes(sites)
    type(SiteData),dimension(:),intent(inout)::sites

    character(len=MAX_NLEN)                  :: sitename
    character(len=MAX_LONG_LINE)             :: header
    real                                     :: lat,long
    real                                     :: altitude
    integer                                  :: siteid
    integer                                  :: numsites
    integer                                  :: numalts
    integer                                  :: ios
    integer                                  :: n,ns
    logical                                  :: is_open, file_exists

    numsites=size(sites)

    ios = 0

    inquire(altlist,exist=file_exists,opened=is_open)
    if ( altlist .ne. 0. .and. file_exists .and. is_open ) then
       ! The altitude file has a header
       numalts=count_records(altlist,1)
  
       if ( numalts <=0 ) then
          call warning("Error in altitude file, ignoring")
          sites%altitude=rnvalid
       else

          read(altlist,*) header
          do n=1,numalts
             read(altlist,*) siteid,sitename,lat,long,altitude
             do ns=1,numsites
                if ( siteid .eq. sites(ns)%site_id ) then
                   sites(ns)%altitude=altitude
                   exit
                endif
             enddo
          enddo
       endif

    else
       sites%altitude=rnvalid
    endif

  end subroutine read_altitudes


  subroutine read_rangelist(use_rangelist,numsites,range_site_ids,             &
                                                              species_present)

    logical,                                      intent(out):: use_rangelist
    integer, dimension(:),           allocatable, intent(out):: range_site_ids
    character(len=8), dimension(:,:),allocatable, intent(out):: species_present

    character(len=MAX_LONG_LINE)                             :: line
    character(len=MAX_LINE), dimension(:),  allocatable      :: fields
    character(len=8),        dimension(:),  allocatable      :: species_ids
    integer,                 dimension(:,:),allocatable      :: species_pres
    character(len=10)                                        :: fmtstr
    real                                                     :: lat,long
    integer                                                  :: numsites
    integer                                                  :: num_range_sites
    integer                                                  :: iwmo 
    integer                                                  :: max_numspec
    integer                                                  :: lc = 0
    integer                                                  :: ios
    integer                                                  :: l,n,s,counter
    logical                                                  :: is_open
    logical                                                  :: file_exists
    logical                                                  :: file_ok=.false.

    inquire(rglist,exist=file_exists,opened=is_open)
    if ( rglist .ne. 0 .and. file_exists .and. is_open ) then
       file_ok=.true.
       use_rangelist=.true.
       ios = 0
       ! The range list file has a header
       num_range_sites=count_records(rglist,1)

       if ( num_range_sites <=0 ) then 
          call warning("Error in range file, ignoring")
          file_ok=.false.
       endif
    endif

    if (file_ok) then

       write(fmtstr,'(a2,i4,a1)') '(a',MAX_LONG_LINE,')'
       ! First line has the list of species.  Parse the line.
       read(rglist,fmtstr,iostat=ios) line
       if ( ios .eq. 0 ) then
          call split_line(adjustl(trim(line)),fields)
          max_numspec=size(fields)-3
       else
          use_rangelist=.false.
          allocate(range_site_ids(1))
          allocate(species_present(1,1))
          range_site_ids=0
          species_present='NA'
          return
       endif

       allocate(species_ids(max_numspec))
       allocate(species_pres(num_range_sites,max_numspec))

       allocate(species_present(num_range_sites,max_numspec))
       allocate(range_site_ids(num_range_sites))

       ! Set up the list of species
       do n=1,max_numspec
          call quote_strip(fields(n+3))
          species_ids(n)=fields(n+3)
       enddo

       ! Now read the species lists, per site.
       ! do while doesn't really work quite as we'd like (we want repeat until)
       lc=0
       do while ( (lc .lt. num_range_sites)  )

          lc = lc + 1
          read(rglist,*,iostat=ios) iwmo,lat, long,                            &
               (species_pres(lc,l),l=1,max_numspec)
          range_site_ids(lc)=iwmo
       end do

       do n=1,max_numspec
          do s=1,num_range_sites
             if ( species_pres(s,n) == 0 ) then
                species_present(s,n)='NP'
             else 
                species_present(s,n)=species_ids(n)
             endif
          enddo
       enddo

    else

       use_rangelist=.false.
       num_range_sites=numsites
       allocate(range_site_ids(1))
       allocate(species_present(1,1))
       range_site_ids=0
       species_present='NA'

    endif

    write(logf,*) 'Species data by site initialized: ',size(range_site_ids)

  end subroutine read_rangelist


  subroutine read_climate(sites)
    type(SiteData), dimension(:), intent(inout) :: sites

    real,           dimension(NTEMPS)   :: tmin,tmax,prcp
    integer                             :: siteid
    character(len=MAX_NLEN)             :: sitename
    real                                :: lat,long
    character(len=MAX_LONG_LINE)        :: header
    integer                             :: numsites
    integer                             :: num_all_clims
    integer                             :: ios
    integer                             :: m,na,ns
    logical                             :: is_open, file_exists

    numsites=size(sites)

    inquire(cfile,exist=file_exists,opened=is_open)
    if ( cfile == 0 .or. (.not. file_exists) .or. (.not. is_open) ) then
       stop 'Unable to find climate file'
    endif

    ! The climate file has a header line
    num_all_clims=count_records(cfile,1)

    if ( num_all_clims <=0 ) then
      call fatal_error("Fatal error in climatology file")
    endif

    ! We don't care about the header
    read(cfile,'(a)'), header

    do na=1,num_all_clims

       read(cfile,*,iostat=ios) siteid,lat,long,                               &
                                (tmin(m),m=1,NTEMPS),(tmax(m),m=1,NTEMPS),     &
                                (prcp(m),m=1,NTEMPS)

       do ns=1,numsites
          if ( sites(ns)%site_id .eq. siteid ) then
             call attach_climate(sites(ns),tmin,tmax,prcp)
             exit
          endif
       enddo

    enddo

    do ns=1,numsites
       if ( any(sites(ns)%tmin     .eq. rnvalid) .or.                          &
            sites(ns)%site_wmo .eq. rnvalid )  then
          write(logf,*) 'No climate data for site number ',sites(ns)%site_id
          sites(ns)%site_wmo=rnvalid
       endif
    enddo

  end subroutine read_climate


  subroutine read_climate_stds(sites)
    type(SiteData), dimension(:), intent(inout) :: sites

    integer                        :: siteid
    character(len=MAX_NLEN)        :: sitename
    real                           :: lat,long
    real, dimension(NTEMPS)        :: tmin,tmax,prcp
    real, dimension(NTEMPS)        :: tmin_std,tmax_std,prcp_std
    character(len=MAX_LONG_LINE)   :: header
    logical                        :: is_open, file_exists
    logical                        :: file_ok=.false.
    integer                        :: numsites, numclimstd
    integer                        :: ios
    integer                        :: m,na,ns
    ! Read or compute standard deviations

    numsites=size(sites)

    inquire(cstdfile,opened=is_open,exist=file_exists)
    if ( cstdfile .ne. 0 .and. file_exists .and.  is_open ) then
       file_ok=.true.
       ! The climate deviations file has a header
       numclimstd=count_records(cstdfile,1)
       if ( numclimstd <=0 ) then
          call warning("Error in climatology std file, ignoring")
          file_ok=.false.
       endif

    endif

    if (file_ok) then

       ! We don't care about the header
       read(cstdfile,'(a)'), header

       do na=1,numclimstd

          read(cstdfile,*,iostat=ios) siteid,lat,long,                         &
                                      (tmin_std(m),m=1,NTEMPS),                &
                                      (tmax_std(m),m=1,NTEMPS),                &
                                      (prcp_std(m),m=1,NTEMPS)

          do ns=1,numsites
             if ( sites(ns)%site_id .eq. siteid ) then
                call attach_climate_std(sites(ns),tmin_std,tmax_std,prcp_std)
                exit
             endif
          enddo

       enddo

    else

       ! Not implemented
       do ns=1,numsites
          call compute_clim_stds(sites(ns),tmin_std,tmax_std,prcp_std)
          call attach_climate_std(sites(ns),tmin_std,tmax_std,prcp_std)
       enddo

    endif

  end subroutine read_climate_stds


  subroutine read_gcm_climate(yx,start_gcm,site)
    integer,          intent(in) :: yx, start_gcm
    type(SiteData),intent(inout) :: site

    integer                      :: siteid
    character(len=MAX_NLEN)      :: sitename
    real                         :: lat,long,year
    real, dimension(NTEMPS)      :: tmin,tmax,prcp
    real, dimension(NTEMPS)      :: tmin_std,tmax_std,prcp_std
    character(len=MAX_LONG_LINE) :: header
    logical                      :: is_open, file_exists
    integer                      :: num_all_clims,numclimstd
    integer                      :: ios
    integer                      :: m,na,ns

    tmin=rnvalid
    tmax=rnvalid
    prcp=rnvalid

    if ( yx .eq. start_gcm ) then
       ! New site.  We don't assume numerical ordering in the site file.
       rewind(cgcmfile)

       ! We don't care about the header
       read(cgcmfile,'(a)'), header

    endif

    do while ( .true. )

       read(cgcmfile,*,end=10) siteid,lat,long,year,                           &
            (tmin(m),m=1,NTEMPS),(tmax(m),m=1,NTEMPS),                         &
            (prcp(m),m=1,NTEMPS)

       if ( site%site_id .eq. siteid .and. int(year) .eq. yx ) then
          call attach_climate(site,tmin,tmax,prcp)
          exit
       endif

    enddo

10  continue

    if ( any(tmin     .eq. rnvalid) .or.                                       &
         any(tmax     .eq. rnvalid) .or.                                       &
         any(prcp     .eq. rnvalid) ) then
       write(logf,*) 'No climate data for site number ',site%site_id
       site%site_wmo=rnvalid
    endif

  end subroutine read_gcm_climate


  subroutine read_sites(site_ids,sites)
    integer,        dimension(:), intent(in)    :: site_ids
    type(SiteData), dimension(:), intent(inout) :: sites

    character(MAX_LINE)                :: header
    character(len=MAX_NLEN)            :: sitename,siteregion
    integer                            :: num_all_sites,numsites
    integer                            :: siteid
    real                               :: wmo,lat,long,elevation,slope
    real                               :: fire_prob,wind_prob
    real                               :: lai,lai_w0
    real                               :: Afc,A_perm_wp,base_h
    real                               :: A0_w0,A_w0,A0_c0,A0_n0,A_c0,A_n0
    real                               :: sbase_w0,sbase_c0,sbase_n0
    real                               :: sigma
    real, dimension(NTEMPS)            :: temp_lapse,prcp_lapse
    integer                            :: ios
    integer                            :: m,na,ns
    logical                            :: file_exists,is_open

    inquire(sfile,exist=file_exists,opened=is_open)
    if ( sfile == 0 .or. (.not. file_exists) .or. (.not. is_open) ) then
       call fatal_error('Unable to find site data file')
    endif

    ! The site file has a header line
    num_all_sites=count_records(sfile,1)
    if ( num_all_sites <=0 ) then
       call fatal_error("Error in site data file")
    endif

    numsites=size(site_ids)

    sites%site_wmo=rnvalid

    ! We don't care about the header
    read(sfile,'(a)'), header

    do na=1,num_all_sites

       read(sfile,*,iostat=ios) siteid,lat,long,wmo,sitename,siteregion,       &
            elevation,slope,Afc,A_perm_wp,lai,base_h,lai_w0,                   &
            A0_w0,A_w0,sbase_w0,fire_prob,wind_prob,A0_c0,                     &
            A0_n0,A_c0,A_n0,sbase_c0,sbase_n0,                                 &
            sigma,(temp_lapse(m),m=1,NTEMPS),(prcp_lapse(m),m=1,NTEMPS)

       do ns=1,numsites

          if ( site_ids(ns) .eq. siteid ) then

             call initialize_site(sites(ns),siteid,sitename,                   &
                  siteregion,lat,long,wmo,elevation,slope,Afc,A_perm_wp,       &
                  lai,base_h,lai_w0,A0_w0,A_w0,sbase_w0,fire_prob,wind_prob,   &
                  A0_c0,A0_n0,A_c0,A_n0,sbase_c0,sbase_n0,sigma,               &
                  temp_lapse, prcp_lapse)

             exit

          endif

       enddo

    enddo

  end subroutine read_sites


  subroutine read_speciesdata(site_ids,species_data)
    integer,           dimension(:), allocatable, intent(inout):: site_ids
    type(SpeciesData), dimension(:), allocatable, intent(out)  :: species_data

    integer                                    :: l, ios, lc = 0
    real                                       :: hm1
    character(len=MAX_PATH)                    :: pathname
    character(len=MAX_LONG_LINE)               :: line
    character(len=10)                          :: fmtstr
    logical                                    :: file_exists,is_open

    integer, dimension(:), allocatable         ::  temp_sites

    character(len=MAX_NLEN)  :: genus_name
    character(len=MAX_NLEN)  :: taxonomic_name
    character(len=8)         :: unique_id
    character(len=MAX_NLEN)  :: common_name
    integer                  :: leaf_type
    integer                  :: genus_sort, genus_id, species_id
    integer                  :: shade_tol,lownutr_tol,stress_tol
    integer                  :: age_tol,dry_tol,flood_tol,fire_tol
    real                     :: max_age,max_diam,max_ht
    real                     :: wood_bulk_dens
    real                     :: leafdiam_a,leafarea_c
    real                     :: deg_day_min,deg_day_opt
    real                     :: deg_day_max
    real                     :: seed_surv,seedling_lg
    real                     :: invader
    real                     :: seed_num,sprout_num
    real                     :: arfa_0,g
    logical                  :: conifer

    integer                  :: num_all_species

    inquire(splist,exist=file_exists,opened=is_open)
    if ( splist == 0 .or. (.not. file_exists) .or. (.not. is_open) ) then
       call fatal_error('Unable to find species-list file')
    endif

    ! The species file has a header line
    num_all_species=count_records(splist,1)
    if ( num_all_species <0 ) then
       call fatal_error('Error in species-list file')
    endif
 
    allocate(species_data(num_all_species))

    ios = 0
    ! Read and discard the first line
    read(splist,'(a500)',iostat=ios) line
    if (ios .ne. 0 ) then
       write(*,*) 'Unable to read the specieslist file'
       stop
    endif

    lc=0
    ! do while doesn't really work quite as we'd like (we want repeat until)
    do while (lc .lt. num_all_species)
       lc = lc + 1
       read(splist,*,end=100)                                                  &
            genus_id,                                                          &
            genus_name,                                                        &
            genus_sort,                                                        &
            taxonomic_name,                                                    &
            common_name,                                                       &
            max_age,                                                           &
            max_diam,                                                          &
            max_ht,                                                            &
            arfa_0,                                                            &
            g,                                                                 &
            wood_bulk_dens,                                                    &
            leafdiam_a,                                                        &
            leafarea_c,                                                        &
            deg_day_min,                                                       &
            deg_day_opt,                                                       &
            deg_day_max,                                                       &
            shade_tol,                                                         &
            dry_tol,                                                           &
            flood_tol,                                                         &
            lownutr_tol,                                                       &
            fire_tol,                                                          &
            stress_tol,                                                        &
            age_tol,                                                           &
            leaf_type,                                                         &
            invader,                                                           &
            seed_num,                                                          &
            sprout_num,                                                        &
            seed_surv,                                                         &
            seedling_lg,                                                       &
            unique_id
       if (leaf_type == 0) then
          conifer=.false.
       else
          conifer=.true.
       endif

       call initialize_species(species_data(lc),lc,genus_name,taxonomic_name,  &
            unique_id,common_name,genus_id,                                    &
            shade_tol,lownutr_tol,stress_tol,age_tol,                          &
            dry_tol,flood_tol,fire_tol,max_age,max_diam,                       &
            max_ht,wood_bulk_dens,rootdepth,leafdiam_a,leafarea_c,             &
            deg_day_min,deg_day_opt,deg_day_max,                               &
            seedling_lg,invader,seed_num,sprout_num,seed_surv,                 &
            arfa_0,g,conifer)

    enddo

100 continue

    ! Adjust species count to exact number read
    if ( lc .ne. num_all_species ) then
       write(*,*) 'Warning : Initial Sites Scan(',num_all_species,             &
                                        ') and Sites Read(',lc,') do not match'
       if ( lc .gt. 0 ) then
          ! This accounts for read errors that are recoverable
          allocate(temp_sites(lc))
          temp_sites=site_ids(1:lc)
          deallocate(site_ids)
          num_all_species=lc
          allocate(site_ids(num_all_species))
          site_ids=temp_sites
       else
          stop "No species data read"
       endif
    end if

    write(logf,*) 'Species data initialized. total read in: ', num_all_species

  end subroutine read_speciesdata


end module Input
