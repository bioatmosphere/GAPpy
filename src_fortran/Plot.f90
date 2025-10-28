module Plot
use Constants
use Species
use Tree

implicit none

   type PlotData
      type(TreeData),   dimension(:),allocatable :: trees
      type(SpeciesData),dimension(:),allocatable :: species
      real,             dimension(:),allocatable :: avail_spec
      real,             dimension(:),allocatable :: seedbank
      real,             dimension(:),allocatable :: seedling
      real,             dimension(:),allocatable :: con_light
      real,             dimension(:),allocatable :: dec_light
      real,             dimension(:),allocatable :: nutrient
      real                                       :: seedling_number
      integer                                    :: numspecies
      integer                                    :: numtrees
      integer                                    :: fire, wind
   end type PlotData

contains


   subroutine initialize_plot(self,species,maxtrees,maxheight)
      class(PlotData),                 intent(inout) :: self  
      type(SpeciesData), dimension(:), intent(inout) :: species
      integer,                         intent(in)    :: maxtrees,maxheight
      
      integer                                        :: n

      if (maxtrees .ne. 0) then
         allocate(self%trees(maxtrees))
      else
         stop "Must allow at least a few trees"
      endif

      if (maxheight .ne. 0) then
         allocate(self%con_light(maxheight))
         allocate(self%dec_light(maxheight))
      else 
         stop "Must have a nonzero maximum height"
      endif

      ! For now the possible species in each plot is assumed to be all for the
      ! site.  But they may not all be able to grow
      self%numspecies    =size(species)
      self%seedling_number=0.0
     
      allocate(self%seedling(self%numspecies))
      allocate(self%seedbank(self%numspecies))
      allocate(self%avail_spec(self%numspecies))
      allocate(self%nutrient(self%numspecies))
      allocate(self%species(self%numspecies))

      ! Copy the species list for the plot from the site's list. 
      do n=1,self%numspecies
         self%species(n) =species(n)
      enddo

      self%numtrees=0
      self%fire    =0
      self%wind    =0

      self%seedbank       =0
      self%seedling       =0
      self%seedling_number=0

      self%avail_spec=0.
      self%seedbank  =0.0
      self%seedling  =0.0
      self%dec_light =0.0
      self%con_light =0.0
      self%nutrient  =1.0

   end subroutine initialize_plot


   subroutine tree_dm_cats(self,genera,field,diam_categories)
   class(PlotData),                       intent(in)  :: self
   character(len=MAX_NLEN), dimension(:), intent(in)  :: genera
   character(len=8),                      intent(in)  :: field
   integer, dimension(size(genera),NHC),  intent(out) :: diam_categories

   character(len=MAX_NLEN)                            :: comp
   integer, dimension(NHC)                            :: tree_diams
   integer                                            :: numitems
   integer                                            :: ig,n

      numitems=size(genera)
   
      diam_categories=0
      do ig=1,numitems
         do n=1,self%numtrees

            if ( field .eq. 'genus' ) then
               comp=self%trees(n)%genus_name
            else if ( field .eq. 'species' ) then
               comp=self%trees(n)%unique_id
            endif

            if ( comp == genera(ig) ) then
                call get_diam_category(self%trees(n),tree_diams)
                diam_categories(ig,:)=diam_categories(ig,:)+tree_diams
            endif
         enddo
      enddo
         
   end subroutine tree_dm_cats


   subroutine sum_over_sg(self,genera,field,basal_area,leaf_bm,biomC,biomN,    &
                                          biomC_std,biomN_std,max_ht,max_diam)
   class(PlotData),                     intent(in)  :: self
   character(len=MAX_NLEN),dimension(:),intent(in)  :: genera
   character(len=8),                    intent(in)  :: field
   real, dimension(size(genera)),       intent(out) :: basal_area,leaf_bm
   real, dimension(size(genera)),       intent(out) :: biomC,biomC_std
   real, dimension(size(genera)),       intent(out) :: biomN,biomN_std
   real, dimension(size(genera)),       intent(out) :: max_ht,max_diam

   integer, dimension(size(genera))                 :: n
   character(len=MAX_NLEN)                          :: comp
   real                                             :: lf_biom, l_cn
   real                                             :: tot_tree_bc, tot_tree_bn
   real                                             :: ni
   integer                                          :: numitems
   integer                                          :: it,is

      numitems=size(genera)
      basal_area=0.0
      biomC     =0.0
      biomN     =0.0
      biomC_std =0.0
      biomN_std =0.0
      leaf_bm   =0.0
      max_ht    =0.0
      max_diam  =0.0
      n         =0

      do is=1,numitems

         do it=1,self%numtrees

           if ( field .eq. 'genus' ) then
                comp=self%trees(it)%genus_name
           else if ( field .eq. 'species' ) then
                comp=self%trees(it)%unique_id
           endif

           if ( comp == genera(is) ) then
                 n(is)=n(is)+1
                 if (self%trees(it)%conifer) then
                     l_cn=con_leaf_c_n
                 else
                     l_cn=dec_leaf_c_n
                 endif
                 basal_area(is)=basal_area(is)+                                &
                                   0.25*pi*self%trees(it)%diam_bht**2
                 lf_biom       =self%trees(it)%leaf_bm
                 leaf_bm(is)   =leaf_bm(is)+lf_biom
                 tot_tree_bc   =self%trees(it)%biomC+lf_biom
                 tot_tree_bn   =self%trees(it)%biomN+lf_biom/l_cn
                 biomC(is)     =biomC(is)+tot_tree_bc
                 biomN(is)     =biomN(is)+tot_tree_bn
                 max_ht(is)    =max(max_ht(is),self%trees(it)%forska_ht)
                 max_diam(is)  =max(max_diam(is),self%trees(it)%diam_bht)
                 ! Sum of squares
                 biomC_std(is) =biomC_std(is)+tot_tree_bc**2
                 biomN_std(is) =biomN_std(is)+tot_tree_bn**2
           endif

         enddo

      enddo

      do is=1,numitems

         if ( n(is) .eq. 0 ) then
            biomC_std(is)=0.0
            biomN_std(is)=0.0
         else 
            ni=1./float(n(is)-1)
            biomC_std(is)=sqrt((biomC_std(is)-biomC(is)**2/n(is))*ni)
            biomN_std(is)=sqrt((biomN_std(is)-biomN(is)**2/n(is))*ni)
         endif

      enddo

   end subroutine sum_over_sg


end module Plot
