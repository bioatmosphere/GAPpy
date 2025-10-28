module Site

  use Constants
  use Species
  use Soil
  use Plot
  use Utilities
  use lists

  implicit none

  !-----------------------------------------------------------------------------
  ! The Site module contains the definition of the SiteData type, which holds
  ! attributes (physical and geographical) of a particular site, as well as
  ! procedures directly operating on site variables
  !-----------------------------------------------------------------------------

  type SiteData
      type(SoilData)                             :: soil
      type(PlotData),   dimension(:),allocatable :: plots
      type(SpeciesData),dimension(:),allocatable :: species
      real,             dimension(:),allocatable :: fc_flood
      real,             dimension(:),allocatable :: fc_drought
      real,             dimension(:),allocatable :: fc_degday
      character(len=MAX_NLEN)                    :: region
      character(len=MAX_NLEN)                    :: site_name
      integer                                    :: site_id
      integer                                    :: numplots
      real                                       :: site_wmo
      real                                       :: latitude,longitude
      real                                       :: elevation,altitude,slope
      real                                       :: leaf_area_ind,leaf_area_w0
      real                                       :: sigma
      real, dimension(NTEMPS)                    :: temp_lapse_r, precip_lapse_r
      real, dimension(NTEMPS)                    :: tmin, tmax, precip
      real, dimension(NTEMPS)                    :: tmin_std, tmax_std 
      real, dimension(NTEMPS)                    :: precip_std
      real                                       :: rain
      real                                       :: pot_evap_day
      real                                       :: act_evap_day
      real                                       :: grow_days
      real                                       :: deg_days
      real                                       :: flood_days
      real                                       :: dry_days_upper_layer
      real                                       :: dry_days_base_layer
      real                                       :: freeze
      real                                       :: fire_prob, wind_prob
  end type SiteData


contains

  !-------------------------------------------------------------------------------
  !  Methods
  !-------------------------------------------------------------------------------

  subroutine initialize_site(self,siteid,sitename,siteregion,lat,long,wmo,     &
        elevation,slope,Afc,A_perm_wp,lai,base_h,                              &
        lai_w0,A0_w0,A_w0,sbase_w0,fire_prob,                                  &
        wind_prob,A0_c0,A0_n0,A_c0,A_n0,sbase_c0,                              &
        sbase_n0,sigma,temp_lapse, prcp_lapse)

    type(SiteData),  intent(inout) :: self
    integer                        :: siteid
    character(len=MAX_NLEN)        :: sitename,siteregion
    real                           :: wmo,lat,long,elevation,slope
    real                           :: fire_prob,wind_prob
    real                           :: lai,lai_w0
    real                           :: Afc,A_perm_wp,base_h
    real                           :: A0_w0,A_w0,A0_c0,A0_n0,A_c0,A_n0
    real                           :: sbase_w0,sbase_c0,sbase_n0
    real                           :: sigma
    real, dimension(NTEMPS)        :: temp_lapse,prcp_lapse

    integer                        :: n

    !  Initialize properties from the site file
    !  Other properties will be added by other routines

    self%site_id          = siteid
    self%site_name        = sitename
    self%region           = siteregion
    self%latitude         = lat
    self%longitude        = long
    self%site_wmo         = wmo
    self%elevation        = elevation
    self%slope            = slope
    self%sigma            = sigma
    self%temp_lapse_r     = temp_lapse
    self%precip_lapse_r   = prcp_lapse
    self%leaf_area_w0     = lai_w0
    self%leaf_area_ind    = lai
    self%fire_prob        = fire_prob
    self%wind_prob        = wind_prob

    self%soil%A_field_cap = Afc
    self%soil%A_perm_wp   = A_perm_wp
    self%soil%base_h      = base_h
    self%soil%A0_w0       = A0_w0
    self%soil%A_w0        = A_w0
    self%soil%BL_w0       = sbase_w0
    self%soil%A0_c0       = A0_c0
    self%soil%A0_n0       = A0_n0
    self%soil%A_c0        = A_c0
    self%soil%A_n0        = A_n0
    self%soil%BL_c0       = sbase_c0
    self%soil%BL_n0       = sbase_n0
    
    call adjustForAltitude(self)

  end subroutine initialize_site


  subroutine attach_climate(self,tmin,tmax,prcp)
    type(SiteData),          intent(inout):: self
    real, dimension(NTEMPS), intent(in)   :: tmin,tmax,prcp

    self%tmin    = tmin
    self%tmax    = tmax
    self%precip  = prcp*mm_to_cm

  end subroutine attach_climate


  subroutine attach_climate_std(self,tmin_std,tmax_std,prcp_std)
    type(SiteData),          intent(inout):: self
    real, dimension(NTEMPS), intent(in)   :: tmin_std,tmax_std,prcp_std

    self%tmin_std   = tmin_std
    self%tmax_std   = tmax_std
    self%precip_std = prcp_std*mm_to_cm

  end subroutine attach_climate_std


  subroutine attach_species(self,species_data,range_species_ids)
      type(SiteData),                           intent(inout) :: self
      type(SpeciesData),          dimension(:), intent(inout) :: species_data
      character(len=*), optional, dimension(:), intent(in)    :: range_species_ids

      integer                                          :: num_all_species
      integer                                          :: num_site_species
      integer                                          :: num_range_species
      integer                                          :: n, nn, ns

      ! Attach the species list for the current site
      ! Makes use of overloading = for species type.

       if ( present(range_species_ids) ) then

         num_all_species=size(species_data)
         num_range_species=size(range_species_ids)
         num_site_species=count(range_species_ids .ne. 'NP')

         if ( num_site_species == 0 ) then

            allocate(self%species(0))

         else

           do n=1,num_all_species
              do nn=1,num_range_species
                if (species_data(n)%unique_id .eq. range_species_ids(nn)) then
                   call append(self%species,species_data(n))
                endif
              enddo
            enddo

         endif

       else
          ! No range list, all species in this site

          do n=1,num_all_species
             call append(self%species,species_data(n))
          enddo

       endif

       !Now that we have the species info, initialize plots

       !Use the input parameter numplots (someday maybe this can vary?)
       self%numplots=numplots

       allocate(self%plots(self%numplots))
       do n=1,self%numplots
         call initialize_plot(self%plots(n),self%species,maxtrees,maxheight)
       enddo

   end subroutine attach_species


   subroutine adjustForAltitude(self)
     type(SiteData), intent(inout) ::  self

     real       :: tav
     integer    :: z

     tav = 0.0

     if ( self%altitude .ne. rnvalid ) then

        do  z=1,12
           self%tmax(z)=self%tmax(z)-(self%altitude-self%elevation)*         &
                self%temp_lapse_r(z)*0.01
           self%tmin(z)=self%tmin(z)-(self%altitude-self%elevation)*         &
                self%temp_lapse_r(z)*0.01
           tav = tav + 0.5 * ( self%tmax(z) + self%tmin(z) )
           self%precip(z)=(max(self%precip(z) +                              &
                (self%altitude-self%elevation)*                              &
                self%precip_lapse_r(z)*0.001,0.0))
        end do
        tav = tav / 12.0
        self%freeze = kron( -tav - 2.0)

     endif

   end subroutine adjustForAltitude


   subroutine compute_clim_stds(self,tmin_stds,tmax_stds,prcp_stds)
     type(SiteData),          intent(inout) :: self
     real, dimension(NTEMPS), intent(inout) :: tmin_stds,tmax_stds,prcp_stds
       !   Dummy for now
   end subroutine compute_clim_stds


   subroutine write_site_csv(self,site_unit)
      use csv_file
      type(SiteData), intent(in) :: self
      integer,        intent(in) :: site_unit

      call csv_write(site_unit,self%rain,.false.)
      call csv_write(site_unit,self%pot_evap_day,.false.)
      call csv_write(site_unit,self%act_evap_day,.false.)
      call csv_write(site_unit,self%grow_days,.false.)
      call csv_write(site_unit,self%deg_days,.false.)
      call csv_write(site_unit,self%dry_days_upper_layer,.false.)
      call csv_write(site_unit,self%dry_days_base_layer,.false.)
      call csv_write(site_unit,self%flood_days,.true.)

  end subroutine write_site_csv


end module Site
