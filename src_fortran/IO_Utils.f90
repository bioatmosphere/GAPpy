module IO

  use Constants
  use Parameters
  use GenusGroups
  use FileUtils
  use csv_file

  implicit none

contains

  !:............................................................................:
  ! Utility subroutines to open/close files and write headers

  subroutine read_file_list(filelist)

    character(len=*)               ::   filelist
    character(len=MAX_FILE)        ::   input_directory, output_directory
    character(len=MAX_DIR)         ::   climate_directory, site_directory
    character(len=MAX_DIR)         ::   sitelist_directory
    integer                        ::   nunit, ios
    logical                        ::   file_exists
    namelist /filenames/input_directory,output_directory,                      &
                        climate_directory,site_directory,sitelist_directory

    ! This reads the names of the directories for the input files. If the filelist 
    ! file is not present it sets them to default values.

    ! Initialize to default values
    input_directory='input_data'
    output_directory='output_data'
    climate_directory='input_data'
    site_directory='input_data'
    sitelist_directory='input_data'

    inquire(file=filelist,exist=file_exists)
    if ( .not. file_exists ) then
       write(*,*) 'File list file not specified, using defaults'
    else
       nunit=unit_number()
       open(nunit,file=filelist,status='unknown',iostat=ios)
       if (ios .eq. 0) then
          read(nunit,filenames,iostat=ios)
          if (ios .ne. 0) then
             write(*,*) 'Error reading ',filelist,': ', ios,' Using defaults'
          end if
       endif
       close(nunit)
    endif

    inputdir   = trim(input_directory)
    outputdir  = trim(output_directory)
    climatedir = trim(climate_directory)
    sitedir    = trim(site_directory)
    slistdir   = trim(sitelist_directory)

  end subroutine read_file_list


  subroutine open_inputFiles(filelist)
    character(len=*)           :: filelist

    integer                    :: ios
    character(len=MAX_FILE)    :: filename
    character(len=MAX_PATH)    :: pathname
    character(len=10)          :: code_id
    logical                    :: file_exists

    call get_cwd(fqcwd)

    call read_file_list(filelist)

    write(code_id,'(a,a)') codename,version_id

    ! Logfile (must be opened as an "input" file so other input can write to it)
    call build_pathname(outputdir,'log.txt',pathname)
    logf=open_file(pathname)

    ! Parameters
    call build_filename(code_id,'_runtime.txt',filename)
    call build_pathname(inputdir,filename,pathname)
    rt_file=open_file(pathname,'r')

    ! Site list
    call build_filename(code_id,'_sitelist.txt',filename)
    call build_pathname(slistdir,filename,pathname)
    slist=open_file(pathname,'r')

    ! Species list
    call build_filename(code_id,'_specieslist.csv',filename)
    call build_pathname(inputdir,filename,pathname)
    splist=open_file(pathname,'r')

    ! Site info
    call build_filename(code_id,'_site.csv',filename)
    call build_pathname(sitedir,filename,pathname)
    sfile=open_file(pathname,'r')

    ! Climate info
    call build_filename(code_id,'_climate.csv',filename)
    call build_pathname(climatedir,filename,pathname)
    cfile=open_file(pathname,'r')

    !  Optional files

    ! Climate standard-deviation info
    call build_filename(code_id,'_climate_stddev.csv',filename)
    call build_pathname(climatedir,filename,pathname)
    cstdfile=open_file(pathname,'r')

    ! Climate data from GCM
    call build_filename(code_id,'_climate_GCM.csv',filename)
    call build_pathname(climatedir,filename,pathname)
    cgcmfile=open_file(pathname,'r')

    ! Rangelist
    call build_filename(code_id,'_rangelist.csv',filename)
    call build_pathname(inputdir,filename,pathname)
    rglist=open_file(pathname,'r')
    if ( rglist < 0 ) then
        write(*,*) 'Unable to open rangelist file :', pathname
        write(*,*) 'Reverting to all species in all sites'
    endif

    ! Altitudes
    call build_filename(code_id,'_altitudes.csv',filename)
    call build_pathname(inputdir,filename,pathname)
    altlist=open_file(pathname,'r')
    if ( altlist < 0 ) then
        write(*,*) 'Unable to open altitude file :',                           &
                                              pathname(1:len_trim(pathname))
        write(*,*) 'Leaving altitudes unchanged'
    endif

  end subroutine open_inputFiles


  subroutine open_outputFiles

    character(len=MAX_PATH) :: pathname
    integer                 :: ios

    ! Open files 

    ! CARBON AND NITROGEN
    call build_pathname(outputdir,'CarbonAndNitrogen.csv',pathname)
    c_and_n=open_file(pathname)

    ! Site/Climate
    call build_pathname(outputdir,'Climate.csv',pathname)
    clim_unit=open_file(pathname)

    ! Genus data
    call build_pathname(outputdir,'Genus_Data.csv',pathname)
    biom_by_g=open_file(pathname)

    ! Species data
    call build_pathname(outputdir,'Species_Data.csv',pathname)
    biom_by_s=open_file(pathname)

    ! Plotlevel data
    if (plot_level_data) then
       ! Plotlevel output
       call build_pathname(outputdir,'Plotlevel_Genus_Data.csv',pathname)
       pl_biom_by_g=open_file(pathname)

       call build_pathname(outputdir,'Plotlevel_Species_Data.csv',pathname)
       pl_biom_by_s=open_file(pathname)
    endif

    ! Tree data
    if (tree_level_data) then
       ! Plotlevel Individual Tree Output
       call build_pathname(outputdir,'Plotlevel_Tree_Data.csv',pathname)
       pl_tree=open_file(pathname)
    endif

  end subroutine open_outputFiles


  subroutine write_headers(species_present)
    type(Groups), intent(in):: species_present

    integer                 :: i,gg,xx
    character(len=1009)     :: tmpstring
    character(len=17)       :: stdtag1
    character(len=3)        :: stdtag2 = '_SD'
    character(len=20)       :: stdtag3
    character(len=MAX_PATH) :: pathname


    ! Write file headers

    ! Carbon and nitrogen (soil)
    call csv_write(c_and_n,'siteID',.false.)
    call csv_write(c_and_n,'year',.false.)
    call csv_write(c_and_n,'a0c0',.false.)
    call csv_write(c_and_n,'ac0',.false.)
    call csv_write(c_and_n,'a0n0',.false.)
    call csv_write(c_and_n,'an0',.false.)
    call csv_write(c_and_n,'bc0',.false.)
    call csv_write(c_and_n,'bn0',.false.)
    call csv_write(c_and_n,'soilresp',.false.)
    call csv_write(c_and_n,'biomassC',.false.)
    call csv_write(c_and_n,'C_into_A0',.false.)
    call csv_write(c_and_n,'net_C_into_A0',.false.)
    call csv_write(c_and_n,'net_prim_prodC',.false.)
    call csv_write(c_and_n,'biomassN',.false.)
    call csv_write(c_and_n,'N_into_A0',.false.)
    call csv_write(c_and_n,'net_N_into_A0',.false.)
    call csv_write(c_and_n,'net_prim_prodN',.false.)
    call csv_write(c_and_n,'avail_n',.true.)

    ! Climate
    call csv_write(clim_unit,'siteID',.false.)
    call csv_write(clim_unit,'year',.false.)
    call csv_write(clim_unit,'rain',.false.)
    call csv_write(clim_unit,'pet',.false.)
    call csv_write(clim_unit,'aet',.false.)
    call csv_write(clim_unit,'grow',.false.)
    call csv_write(clim_unit,'degd',.false.)
    call csv_write(clim_unit,'dryd_upper',.false.)
    call csv_write(clim_unit,'dryd_base',.false.)
    call csv_write(clim_unit,'flood_d',.true.)

    ! Data by Genus
    call csv_write(biom_by_g,'siteID',.false.)
    call csv_write(biom_by_g,'year',.false.)
    call csv_write(biom_by_g,'genus',.false.)
    call csv_write(biom_by_g,'<0',.false.)
    call csv_write(biom_by_g,'0-8',.false.)
    call csv_write(biom_by_g,'8-28',.false.)
    call csv_write(biom_by_g,'-48',.false.)
    call csv_write(biom_by_g,'-68',.false.)
    call csv_write(biom_by_g,'-88',.false.)
    call csv_write(biom_by_g,'>88',.false.)
    call csv_write(biom_by_g,'max_diam',.false.)
    call csv_write(biom_by_g,'max_ht',.false.)
    call csv_write(biom_by_g,'leaf_area_ind',.false.)
    call csv_write(biom_by_g,'basal_area',.false.)
    call csv_write(biom_by_g,'total_biomC',.false.)
    call csv_write(biom_by_g,'pl_biomC_std',.false.)
    call csv_write(biom_by_g,'total_biomN',.false.)
    call csv_write(biom_by_g,'pl_biomN_std',.true.)

    ! Data by Species
    call csv_write(biom_by_s,'siteID',.false.)
    call csv_write(biom_by_s,'year',.false.)
    call csv_write(biom_by_s,'genus',.false.)
    call csv_write(biom_by_s,'species',.false.)
    call csv_write(biom_by_s,'<0',.false.)
    call csv_write(biom_by_s,'0-8',.false.)
    call csv_write(biom_by_s,'8-28',.false.)
    call csv_write(biom_by_s,'-48',.false.)
    call csv_write(biom_by_s,'-68',.false.)
    call csv_write(biom_by_s,'-88',.false.)
    call csv_write(biom_by_s,'>88',.false.)
    call csv_write(biom_by_s,'max_diam',.false.)
    call csv_write(biom_by_s,'max_hgt',.false.)
    call csv_write(biom_by_s,'leaf_area_ind',.false.)
    call csv_write(biom_by_s,'basal_area',.false.)
    call csv_write(biom_by_s,'total_biomC',.false.)
    call csv_write(biom_by_s,'pl_biomC_std',.false.)
    call csv_write(biom_by_s,'total_biomN',.false.)
    call csv_write(biom_by_s,'pl_biomN_std',.true.)

    if (plot_level_data) then

       ! Plotlevel Data By Genus
       call csv_write(pl_biom_by_g,'siteID',.false.)
       call csv_write(pl_biom_by_g,'year',.false.)
       call csv_write(pl_biom_by_g,'plot',.false.)
       call csv_write(pl_biom_by_g,'genus',.false.)
       call csv_write(pl_biom_by_g,'species',.false.)
       call csv_write(pl_biom_by_g,'<0',.false.)
       call csv_write(pl_biom_by_g,'0-8',.false.)
       call csv_write(pl_biom_by_g,'8-28',.false.)
       call csv_write(pl_biom_by_g,'-48',.false.)
       call csv_write(pl_biom_by_g,'-68',.false.)
       call csv_write(pl_biom_by_g,'-88',.false.)
       call csv_write(pl_biom_by_g,'>88',.false.)
       call csv_write(pl_biom_by_g,'max_diam',.false.)
       call csv_write(pl_biom_by_g,'max_hgt',.false.)
       call csv_write(pl_biom_by_g,'leaf_area_ind',.false.)
       call csv_write(pl_biom_by_g,'basal_area',.false.)
       call csv_write(pl_biom_by_g,'total biomC',.false.)
       call csv_write(pl_biom_by_g,'total biomC_std',.false.)
       call csv_write(pl_biom_by_g,'total biomN',.false.)
       call csv_write(pl_biom_by_g,'total biomN_std',.true.)

       ! Plotlevel Data By Species
       call csv_write(pl_biom_by_s,'siteID',.false.)
       call csv_write(pl_biom_by_s,'year',.false.)
       call csv_write(pl_biom_by_s,'plot',.false.)
       call csv_write(pl_biom_by_s,'genus',.false.)
       call csv_write(pl_biom_by_s,'species',.false.)
       call csv_write(pl_biom_by_s,'<0',.false.)
       call csv_write(pl_biom_by_s,'0-8',.false.)
       call csv_write(pl_biom_by_s,'8-28',.false.)
       call csv_write(pl_biom_by_s,'-48',.false.)
       call csv_write(pl_biom_by_s,'-68',.false.)
       call csv_write(pl_biom_by_s,'-88',.false.)
       call csv_write(pl_biom_by_s,'>88',.false.)
       call csv_write(pl_biom_by_s,'max_diam',.false.)
       call csv_write(pl_biom_by_s,'max_hgt',.false.)
       call csv_write(pl_biom_by_s,'leaf_area_ind',.false.)
       call csv_write(pl_biom_by_s,'basal_area',.false.)
       call csv_write(pl_biom_by_s,'total biomC',.false.)
       call csv_write(pl_biom_by_s,'total biomC_std',.false.)
       call csv_write(pl_biom_by_s,'total biomN',.false.)
       call csv_write(pl_biom_by_s,'total biomN_std',.true.)
    endif

    if (tree_level_data) then
       ! Plotlevel Individual Tree Output
       call csv_write(tld,'siteID',.false.)
       call csv_write(tld,'year',.false.)
       call csv_write(tld,'plot',.false.)
       call csv_write(tld,'genus',.false.)
       call csv_write(tld,'species',.false.)
       call csv_write(tld,'diam bh',.false.)
       call csv_write(tld,'forska_height',.false.)
       call csv_write(tld,'leaf biomass',.false.)
       call csv_write(tld,'stem biomC',.false.)
       call csv_write(tld,'stem biomN',.true.)
    endif

  end subroutine write_headers


  subroutine close_outputFiles
    logical                  :: isopen

    close(c_and_n)
    close(clim_unit)
    close(biom_by_s)
    close(biom_by_g)

    ! Optional files
    inquire(pl_biom_by_g,opened=isopen)
    if (isopen) close(pl_biom_by_g)
    inquire(pl_biom_by_s,opened=isopen)
    if (isopen) close(pl_biom_by_s)
    inquire(tld,opened=isopen)
    if (isopen) close(pl_tree)

  end subroutine close_outputFiles


end module IO
