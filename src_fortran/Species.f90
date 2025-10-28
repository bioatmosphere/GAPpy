module Species
   ! 
   ! Author: K. A. Holcomb
   ! Species attributes and procedures
   !
   ! Methods
   !   initialize_species
   !   copy_species overloads assignment
   !   temp_rsp(self,x)
   !      verified against original 20120727
   !   drought_rsp(self,x)
   !   light_rsp(self,x)
   !      verified against original 20120727
   !   poor_soil_rsp(self,x)
   !      verified against original 20120727
   !   fire_rsp(self,x)

   use Constants

   implicit none

   !-------------------------------------------------------------------------------
   ! The Species module contains the definition of the SpeciesData type, which
   ! holds attributes pertinent to all members of a species, and associated
   ! procedures.
   !-------------------------------------------------------------------------------

   type SpeciesData
      character(len=MAX_NLEN)  ::  genus_name
      character(len=MAX_NLEN)  ::  taxonomic_name
      character(len=8)         ::  unique_id
      character(len=MAX_NLEN)  ::  common_name
      integer                  ::  genus_id, species_id
      integer                  ::  shade_tol,lownutr_tol,stress_tol,age_tol
      integer                  ::  drought_tol,flood_tol,fire_tol
      real                     ::  max_age,max_diam,max_ht
      real                     ::  wood_bulk_dens
      real                     ::  rootdepth
      real                     ::  leafdiam_a,leafarea_c
      real                     ::  deg_day_min,deg_day_opt,deg_day_max
      real                     ::  seed_surv,seedling_lg
      real                     ::  invader
      real                     ::  seed_num,sprout_num
      real                     ::  arfa_0,g
      real                     ::  fc_fire, fc_wind
      real                     ::  fc_degday,fc_drought,fc_flood
      logical                  ::  conifer
   end type SpeciesData

   interface assignment(=)
      module procedure copy_species
   end interface

   private copy_species

contains

   !-------------------------------------------------------------------------------
   ! Methods
   !-------------------------------------------------------------------------------

   subroutine initialize_species(self,species_id,genus_name,taxonomic_name,     &
                                 unique_id,common_name,genus_id,                &
                                 shade_tol,lownutr_tol,stress_tol,age_tol,      &
                                 drought_tol,flood_tol,fire_tol,                &
                                 max_age,max_diam,max_ht,                       &
                                 wood_bulk_dens,rootdepth,leafdiam_a,           &
                                 leafarea_c,deg_day_min,deg_day_opt,deg_day_max,&
                                 seedling_lg,invader,seed_num,                  &
                                 sprout_num,seed_surv,arfa_0,g,conifer)

         class(SpeciesData),      intent(inout):: self
         character(len=MAX_NLEN), intent(in)   :: genus_name
         character(len=MAX_NLEN), intent(in)   :: taxonomic_name
         character(len=8),        intent(in)   :: unique_id
         character(len=MAX_NLEN), intent(in)   :: common_name
         integer,                 intent(in)   :: genus_id, species_id
         integer,                 intent(in)   :: shade_tol,lownutr_tol,stress_tol
         integer,                 intent(in)   :: age_tol,fire_tol
         integer,                 intent(in)   :: drought_tol,flood_tol
         real,                    intent(in)   :: max_age,max_diam,max_ht
         real,                    intent(in)   :: rootdepth,wood_bulk_dens
         real,                    intent(in)   :: leafdiam_a,leafarea_c
         real,                    intent(in)   :: deg_day_min,deg_day_opt
         real,                    intent(in)   :: deg_day_max
         real,                    intent(in)   :: seed_surv,seedling_lg
         real,                    intent(in)   :: invader
         real,                    intent(in)   :: seed_num,sprout_num
         real,                    intent(in)   :: arfa_0,g
         logical,                 intent(in)   :: conifer

         real, dimension(5)                    :: ss,adjust
         data ss/1.1,1.15,1.2,1.23,1.25/
         data adjust/1.5,1.55,1.6,1.65,1.7/


         self%species_id    =species_id
         self%genus_name    =genus_name
         self%taxonomic_name=taxonomic_name
         self%common_name   =common_name
         self%unique_id     =unique_id
         self%genus_id      =genus_id
         self%max_age       =max_age
         self%max_diam      =max_diam
         self%rootdepth     =rootdepth
         self%wood_bulk_dens=wood_bulk_dens
         self%deg_day_min   =deg_day_min
         self%deg_day_max   =deg_day_max
         self%deg_day_opt   =deg_day_opt
         self%shade_tol     =shade_tol
         self%lownutr_tol   =lownutr_tol
         self%drought_tol   =drought_tol
         self%fire_tol      =fire_tol
         self%flood_tol     =flood_tol
         self%stress_tol    =stress_tol
         self%age_tol       =age_tol
         self%conifer       =conifer
         self%invader       =invader
         self%seed_num      =seed_num
         self%sprout_num    =sprout_num
         self%seed_surv     =seed_surv
         self%seedling_lg   =seedling_lg
         self%arfa_0        =arfa_0

         ! Adjustments
         self%leafarea_c    = leafarea_c/hec_to_m2
         self%max_ht        = min(max_ht,rootdepth*80.0/(1+rootdepth))
         self%leafdiam_a    = leafdiam_a*adjust(shade_tol)
         self%g             = g*ss(shade_tol)

         ! Initialized to constants
         self%fc_fire=0.0
         self%fc_wind=0.0
         self%fc_degday=0.0
         self%fc_drought=0.0
         self%fc_flood=0.0

   end subroutine initialize_species


   subroutine copy_species(self,species_data)

      class(SpeciesData),      intent(out)  :: self
      class(SpeciesData),      intent(in)   :: species_data

      self%species_id    =species_data%species_id
      self%genus_name    =species_data%genus_name
      self%taxonomic_name=species_data%taxonomic_name
      self%common_name   =species_data%common_name
      self%unique_id     =species_data%unique_id
      self%genus_id      =species_data%genus_id
      self%max_ht        =species_data%max_ht
      self%max_age       =species_data%max_age
      self%max_diam      =species_data%max_diam
      self%rootdepth     =species_data%rootdepth
      self%wood_bulk_dens=species_data%wood_bulk_dens
      self%leafarea_c    =species_data%leafarea_c
      self%leafdiam_a    =species_data%leafdiam_a
      self%deg_day_min   =species_data%deg_day_min
      self%deg_day_max   =species_data%deg_day_max
      self%deg_day_opt   =species_data%deg_day_opt
      self%shade_tol     =species_data%shade_tol
      self%lownutr_tol   =species_data%lownutr_tol
      self%drought_tol   =species_data%drought_tol
      self%fire_tol      =species_data%fire_tol
      self%flood_tol     =species_data%flood_tol
      self%stress_tol    =species_data%stress_tol
      self%age_tol       =species_data%age_tol
      self%fc_degday     =species_data%fc_degday
      self%fc_fire       =species_data%fc_fire
      self%fc_drought    =species_data%fc_drought
      self%fc_flood      =species_data%fc_flood
      self%fc_wind       =species_data%fc_wind
      self%conifer       =species_data%conifer
      self%invader       =species_data%invader
      self%seed_num      =species_data%seed_num
      self%sprout_num    =species_data%sprout_num
      self%seed_surv     =species_data%seed_surv
      self%seedling_lg   =species_data%seedling_lg
      self%arfa_0        =species_data%arfa_0
      self%g             =species_data%g

   end subroutine copy_species


   function light_rsp(self,al)
      ! Computes available sunlight factor, by tolerance class tshade
      real                           :: light_rsp
      class(SpeciesData), intent(in) :: self
      real,               intent(in) :: al

      real                           :: flight
      integer                        :: kt

      ! Light availability
      real, dimension(5)             :: light_c1, light_c2, light_c3
      data light_c1 /1.01,1.04,1.11,1.24,1.49/
      data light_c2 /4.62,3.44,2.52,1.78,1.23/
      data light_c3 /0.05,0.06,0.07,0.08,0.09/

      kt=self%shade_tol
      
      flight = light_c1(kt) * (1.0-exp(-light_c2(kt) * (al - light_c3(kt))))
      if (flight .lt. 0.0) flight = 0.0
      if (flight .gt. 1.0) flight = 1.0
      light_rsp=flight

   end function light_rsp


   function poor_soil_rsp(sf,nrc)
      !  computes quadratic itpoorent response factors
      !  by itpoorent response class (1=intol, 3=tol) and soil fertility
      real                            :: poor_soil_rsp
      real,              intent(inout):: sf
      !     class(SpeciesData), intent(in)  :: self
      integer,           intent(inout):: nrc
      !    integer,  intent(in)            :: nrc

      !    integer                       :: nrc
      real                          :: fpoor

      ! Soil fertility
      real, dimension(3)            :: fert_c1, fert_c2, fert_c3

      data  fert_c1 /-0.6274,-0.2352,0.2133/
      data  fert_c2 /3.600,2.771,1.789/
      data  fert_c3 /-1.994,-1.550,-1.014/

      !This cannot be converted to a real "method" until it is resolved whether
      ! the following is a bug.  As it is, it changes the avail_n as well as the
      ! required nitrogen on each call.

      !!!WARNING WARNING WARNING POSSIBLE MASSIVE BUG HERE
      nrc=4-nrc
      sf = min(sf,1.0)
      fpoor=fert_c1(nrc)+fert_c2(nrc)*sf+fert_c3(nrc)*sf**2
      if (fpoor.le. 0.0) fpoor=0.0
      if (fpoor.ge. 1.0) fpoor=1.0
      poor_soil_rsp=fpoor*sf

   end function poor_soil_rsp


   subroutine fire_rsp(self, fire)
      !  ability to tolerate fire
      class(SpeciesData), intent(inout) :: self
      integer,            intent(in)    :: fire
   
      real                              :: resp
      integer                           :: k
      real, dimension(6)                :: gama
      data gama/100.0,10.0,1.00,0.1,0.01,0.001/

      k=self%fire_tol

      if (fire == 1) then
         resp=gama(k)
      else
         resp=1.0
      end if

      self%fc_fire=resp

   end subroutine fire_rsp


   subroutine temp_rsp(self,x)
      class(SpeciesData), intent(inout) :: self
      real,               intent(in)    :: x

      real                              :: ftemp
      real                              :: ddmin,ddopt,ddmax
      real                              :: a,b,tmp

      ddmin=self%deg_day_min; ddmax=self%deg_day_max; ddopt=self%deg_day_opt

      a=(ddopt-ddmin)/(ddmax-ddmin)
      b=(ddmax-ddopt)/(ddmax-ddmin)
      
      if (x.ge.ddmax .or. x .le. ddmin) then
         ftemp=0.0
      else
         tmp=((x-ddmin)/(ddopt-ddmin))**a
         ftemp=tmp*((ddmax-x)/(ddmax-ddopt))**b
      end if

      self%fc_degday=ftemp

   end subroutine temp_rsp


   subroutine drought_rsp(self,drydays,drydays_s)
      class(SpeciesData), intent(inout) :: self
      real,               intent(in)    :: drydays,drydays_s

      real                              :: fcdry1, fcdry2

         if (self%drought_tol .eq. 1) then
            if (self%conifer) then
               fcdry1=fdry(drydays_s,1)*0.33
            else
               fcdry1=fdry(drydays_s,1)*0.2
            endif

            fcdry2=fdry(drydays,1)
            self%fc_drought=max(fcdry1,fcdry2)

         else
            self%fc_drought=fdry(drydays,self%drought_tol)
         end if

   end subroutine drought_rsp


   subroutine flood_rsp(self,floodday)
      class(SpeciesData), intent(inout) :: self
      real,               intent(in)    :: floodday

      real                              :: fflood
      integer                           ::k

      real, dimension(6)   ::gama
      data gama/1.0,0.9,0.8,0.7,0.6,0.5/

      k=self%flood_tol

      ! KAH This function always returns 1
      if ( floodday .le. gama(k) ) then
         fflood=1.0
      else
         fflood=1.0
      end if

      self%fc_flood=fflood

   end subroutine flood_rsp


   function fdry(dryday,k)
      real                :: fdry
      real,    intent(in) :: dryday
      integer, intent(in) :: k

      real                :: tmp
      real, dimension(6)  :: gama
      data gama/0.50,0.45,0.35,0.25,0.15,0.05/

      tmp=max(gama(k)-dryday,0.0)
      fdry=(tmp/gama(k))**0.5

   end function fdry


   subroutine write_species_csv(self,spunit,with_tree)
      use csv_file
      use FileUtils
      class(SpeciesData), intent(in) :: self
      integer,            intent(in) :: spunit
      logical,            intent(in) :: with_tree

      call csv_write(spunit,self%genus_name,.false.)
      call csv_write(spunit,self%taxonomic_name,.false.)
      call csv_write(spunit,self%seed_surv,.false.)
      call csv_write(spunit,self%seedling_lg,.false.)
      call csv_write(spunit,self%seed_num,.false.)
      call csv_write(spunit,self%sprout_num,.false.)
      call csv_write(spunit,self%fc_fire,.false.)
      call csv_write(spunit,self%fc_wind,.false.)
      call csv_write(spunit,self%fc_degday,.false.)
      call csv_write(spunit,self%fc_drought,.false.)
      if (with_tree) then
         call csv_write(spunit,self%fc_flood,.false.)
      else
         call csv_write(spunit,self%fc_flood,.true.)
      endif

   end subroutine write_species_csv


end module Species
