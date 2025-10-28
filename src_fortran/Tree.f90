module Tree
   !
   ! Author: K. A. Holcomb
   ! 
   ! Tree attributes and procedures relevant to tree life cycle
   !
   ! Methods:
   !    initialize_tree
   !        checked 20120726
   !    copy_tree
   !        checked 20120726
   !        would like to make this overload = but it's not working yet
   !        might be that it has to be a type-bound procedure for this to work
   !    stem_shape
   !        verified against original 20120726  
   !    forska_height
   !        verified against original 20120726  
   !    stem_biomass_c
   !        verified against original 20120726  
   !    twig_biomass_c
   !        verified against original 20120726  
   !    leaf_biomass_c
   !        verified against original 20120727  
   !    lai_biomass_c
   !        verified against original 20120726  
   !    biomass_c
   !        verified against original 20120726  
   !    biomass_n
   !        verified against original 20120726  
   !    env_stress
   !        verified against original 20120727  
   !    age_survival
   !        verified against original 20120730
   !    growth_survival
   !        verified against original 20120730
   !    max_growth
   !        verified against original 20120727  

   use Constants
   use Species
   use Random

   implicit none

   ! Constants global within this module
   ! KAH
   ! Not sure what this is, it may have something to do with atomic mass of C
   ! It only appears in two functions so it's here rather than in Constants
   real,     parameter      :: tC=3.92699e-5
   ! Not sure what this is either but it can be a global in the module
   ! Unclear whether it ought to be the same everywhere, but for now it is.
   real,     parameter      :: beta = 1.0

   ! Tree object

   type, extends(SpeciesData) :: TreeData
     real                     :: diam_max
     real                     :: diam_bht
     real                     :: diam_canht
     real                     :: canopy_ht
     real                     :: foret_ht
     real                     :: forska_ht
     real                     :: leaf_bm   ! includes roots
     real                     :: biomC
     real                     :: biomN
     integer                  :: species_index
     logical                  :: mort_marker
   end type TreeData

   !   interface assignment(=)
   !     module procedure copy_tree
   !   end interface

   !   private copy_tree

contains

   !==================================================================
   ! 17 Methods
   !==================================================================

   subroutine initialize_tree(self, tree_species, si)
      class(TreeData),             intent(inout) :: self
      type(SpeciesData), optional, intent(in)    :: tree_species 
      integer,           optional, intent(in)    :: si

      ! Constructor
      self%diam_bht   =0.0
      self%diam_canht =0.0
      self%canopy_ht  =std_ht
      self%foret_ht   =1.0
      self%forska_ht  =1.0
      self%biomC      =0.0
      self%biomN      =0.0
      self%mort_marker=.false.

      if ( present(tree_species) ) then
         self%SpeciesData=tree_species
      !       call copy_species(self%SpeciesData,tree_species)
      endif

      if ( present(si) ) then
         self%species_index=si
      else
         self%species_index =0
      endif

   end subroutine initialize_tree


   subroutine copy_tree(self,tree)
      class(TreeData), intent(inout):: self
      type(TreeData),  intent(in)   :: tree

      self%diam_max      = tree%diam_max
      self%diam_bht      = tree%diam_bht
      self%diam_canht    = tree%diam_canht
      self%canopy_ht     = tree%canopy_ht
      self%foret_ht      = tree%foret_ht
      self%forska_ht     = tree%forska_ht
      self%leaf_bm       = tree%leaf_bm
      self%biomC         = tree%biomC
      self%biomN         = tree%biomN
      self%species_index = tree%species_index
      self%mort_marker   = tree%mort_marker
      self%SpeciesData   = tree%SpeciesData
      !  Overloading = for the species object works.  Overloading = for the tree
      !  object that includes an overloaded = does not (yet) work so copy_tree must
      !  be called explicitly.  Probably need type-bound procedures.
      !    In case species= doesn't work (must comment out private statement in 
      !      species module):
      !      call copy_species(self%SpeciesData,tree%SpeciesData)

   end subroutine copy_tree


   subroutine update_tree(self,tree_species)
      class(TreeData),   intent(inout)        :: self
      type(SpeciesData), intent(in)           :: tree_species 
      ! A bit of extra work required to have the more flexible object structure

      self%SpeciesData   = tree_species

   end subroutine update_tree


   function env_stress(self,shade)
      class(TreeData), intent(in)   :: self
      real                          :: env_stress
      real,              intent(in) :: shade

      env_stress=self%fc_degday*self%fc_drought*self%fc_flood*shade

   end function env_stress


   subroutine stem_shape(tree)
      !   basal diameter = dhshape(0.0,h,dbh) KAH use a tree_0 object
      !   diameter at canopy height = dhshape(hc,h,dbh)
      !   Note: mentioned above "beta" can be 1, 1.5 ,2 , 3
      !   KAH Right now beta is set to 1 for all these routines.  Is it
      !   distinct here?  If so it needs a different name.
      class(TreeData), intent(inout) :: tree

      real                           :: hc,h,dbh,dhshape

      hc=tree%canopy_ht; h=tree%forska_ht;dbh=tree%diam_bht

      if (h .le. hc .or. h .le. std_ht) then
         dhshape=dbh
      else
         dhshape=((h-hc)/(h-std_ht))**(1.0/beta)*dbh
      end if

      tree%diam_canht=dhshape

   end subroutine stem_shape


   subroutine forska_height(tree)
      !==========================================================
      !    forska_height:  height-dbh relation equation
      !    Original by yan xiaodong
      !    modified from Forska model
      !----------------------------------------------------------
      !  hmax:  input maximum tree height (m)
      !  d: input diameter (cm)
      !  arfa0: input initial slope of d-h relation
      !  Basically, most shade tolerant 1.5
      !                  second position: 1.25
      !                  middle class:  1.0
      !                  second shade intolerant 0.75
      !                  most shade intolerant 0.5       
      !  forska_height: height by Forska model computing m
      !----------------------------------------------------------
      class(TreeData), intent(inout):: tree

      real                          :: d, par, hmax
      real                          :: delta_ht

      d=tree%diam_bht
      hmax=tree%max_ht;par=tree%arfa_0
      
      delta_ht=hmax-std_ht
      tree%forska_ht= std_ht + delta_ht*( 1.0 - exp(- (par*d/delta_ht)))

   end subroutine forska_height


   subroutine biomass_c(tree)
      class(TreeData), intent(inout):: tree

      type(TreeData)                :: tree_0
      real                          :: d,h,hc,hrt,bulk
      real                          :: dc, bd
      real                          :: stembc, twigbc
      real                          :: abovegr_c, root_c

      d=tree%diam_bht;h=tree%forska_ht;hc=tree%canopy_ht;
      hrt=tree%rootdepth;bulk=tree%wood_bulk_dens

      ! This doesn't work (yet) for some reason, seems like it should.
      !    tree_0=tree
      call copy_tree(tree_0,tree)
      tree_0%canopy_ht=0.0

      call stem_shape(tree)   !dc
      call stem_shape(tree_0) !bd

      stembc=stem_biomass_c(tree_0)
      twigbc=twig_biomass_c(tree)

      abovegr_c = stembc+twigbc
      root_c    = stembc*hrt/h+twigbc/2.0

      tree%biomC=(abovegr_c+root_c)

   end subroutine biomass_c


   subroutine biomass_n(tree)
      type(TreeData), intent(inout):: tree

      tree%biomN=tree%biomC/stem_c_n

   end subroutine biomass_n


   subroutine leaf_biomass_c(tree)
      !    leaf biomass C = blai*leafa_c
      !    x2 for roots also
      type(TreeData), intent(inout) :: tree

      tree%leaf_bm = lai_biomass_c(tree)*tree%leafarea_c*2.0

   end subroutine leaf_biomass_c


   function lai_biomass_c(tree)
      real                       :: lai_biomass_c
      type(TreeData), intent(in) :: tree

      real                       :: dc, leafd_a

      dc=tree%diam_canht; leafd_a=tree%leafdiam_a
      lai_biomass_c = dc * dc * leafd_a

   end function lai_biomass_c


   function age_survival(tree)
      !===================================================================
      !    Death check by life age
      !    can reach life span (max age) as a probability 1%,0.1%,0.01%
      !    1% often seen, 0.1% some survive, 0.01% few seen
      !-------------------------------------------------------------------
      logical                    :: age_survival
      type(TreeData), intent(in) :: tree

      integer                    :: k
      real                       :: agemax
      real                       :: rand_val

      real, dimension(3)         :: check
      data check /4.605,6.908,11.51/
      !PER 1---4.605   PER 0.1---6.908   PER 0.001---11.5129

      k=tree%age_tol
      agemax=tree%max_age

      rand_val=urand()
      if (rand_val .lt. (check(k)/agemax)) then
         age_survival=.false.
      else
         age_survival=.true.
      end if

   end function age_survival


   function growth_survival(tree)
      !==================================================================
      !    Death check by growth insufficiency
      !    can survive 5,10,20,40,80 years as a probability 5%.
      !------------------------------------------------------------------
      logical                    :: growth_survival
      type(TreeData), intent(in) :: tree

      integer                    :: k
      real                       :: rand_val

      real, dimension(5)         :: check
      data check /0.31,0.34,0.37,0.40,0.43/
      ! old values  data check/0.0368,0.0722,0.1391,0.2589,0.4507/
      !     5---0.45   10,---0.26  20---0.14   40---0.07   80---0.04

      k=tree%stress_tol
      rand_val=urand()

      if (tree%mort_marker .and. rand_val .lt. check(k)) then
         growth_survival=.false.
      else
         growth_survival=.true.
      end if

   end function growth_survival
        

   subroutine max_growth(tree)

      !===========================================================
      !     max_growth: optimal growth function
      !     by Yan Xiaodong
      !-----------------------------------------------------------
      type(TreeData), intent(inout):: tree

      real                         :: d,dc,h,dm,hm,g,arfa0
      real                         :: beginning
      !-----------------------------------------------------------
      !   dc: diameter at tree's canopy height (cm) 
      !   d: input diameter (cm)
      !   h: input height (m)
      !   dm: input maximum diameter (cm)
      !   hm: input maximum height (m)
      !   g: input optimal growth parameter
      !   arfa0: input initial slope of d-h relation
      !-----------------------------------------------------------
      dc=tree%diam_canht;d=tree%diam_bht;h=tree%forska_ht
      dm=tree%max_diam;hm=tree%max_ht;g=tree%g;arfa0=tree%arfa_0

      beginning=arfa0*exp(-arfa0*d/(hm-std_ht))*d
      tree%diam_max=g*d*(1.0-d*h/dm/hm)/(2.0*h+beginning)

   end subroutine max_growth


   function stem_biomass_c(tree)
      !==========================================================
      ! stem_biomass_c(bd,h,bulk_density):  stem biomass (tC)
      !     d: diameter base height dbh (cm)
      !     h:tree height (m)
      !     dc: diameter at canopy height (cm)
      !     hc: canopy height (m)
      !     bd: basal diameter (cm)
      !     bulk_density: wood bulk_density (t/m^3)
      !     note: root stem biomass C = bstemc(bd,1.0,bulk_density)
      !----------------------------------------------------------
      ! KAH I just get fatigued trying to improve the temporary varnames

      real                       :: stem_biomass_c
      type(TreeData), intent(in) :: tree

      real                       :: bd,h,bulk_density
      real                       :: yxd

      ! BD is actually basal diameter but when we call this routine we use an
      ! object with a "canopy height" of zero.
      bd=tree%diam_canht;h=tree%forska_ht
      bulk_density=tree%wood_bulk_dens

      yxd   = tC*bulk_density*beta/(beta+2.0)
      stem_biomass_c=yxd*bd*bd*h*0.90

   end function stem_biomass_c


   function twig_biomass_c(tree)
      !==========================================================
      ! Similar to stem_biomass_c, for twigs
      !---------------------------------------------------------
      real                       :: twig_biomass_c
      type(TreeData), intent(in) :: tree

      real                       :: dc,hc,h,bulk_density
      real                       :: yxd

      dc=tree%diam_canht;hc=tree%canopy_ht;h=tree%forska_ht
      bulk_density=tree%wood_bulk_dens

      yxd=tC*bulk_density*(2.0/(beta+2.0)-0.33)
      twig_biomass_c=yxd*dc*dc*(h-hc)

   end function twig_biomass_c


   subroutine get_diam_category(self,diam_category)
      class(TreeData),         intent(in)  :: self
      integer, dimension(NHC), intent(out) :: diam_category
      real                                 :: dm

      dm=self%diam_bht
      diam_category=0

      if ( self%mort_marker ) then
         diam_category(1)=1
      endif

      if ( dm <= 8.0 ) then
         diam_category(2)=1
      else if ( dm <= 28.0 ) then
         diam_category(3)=1
      else if ( dm <= 48.0 ) then
         diam_category(4)=1
      else if ( dm <= 68.0 ) then
         diam_category(5)=1
      else if ( dm <= 88.0 ) then
         diam_category(6)=1
      else
         diam_category(7)=1
      endif

   end subroutine get_diam_category


   subroutine write_tree_csv(self,tree_unit)
      use csv_file
      class(TreeData), intent(in) :: self
      integer,         intent(in) :: tree_unit

      call csv_write(tree_unit,trim(adjustl(self%genus_name)),.false.)
      call csv_write(tree_unit,self%unique_id,.false.)
      call csv_write(tree_unit,self%diam_bht,.false.)
      call csv_write(tree_unit,self%forska_ht,.false.)
      call csv_write(tree_unit,self%leaf_bm,.false.)
      call csv_write(tree_unit,self%biomC,.false.)
      call csv_write(tree_unit,self%biomN,.true.)

   end subroutine write_tree_csv


end module Tree
