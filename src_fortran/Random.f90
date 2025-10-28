!:.............................................................................:
!:                                                                             :
!:               	     Random Number Generators                                :
!:                                                                             :
!:.............................................................................:

module Random

   implicit none

   ! Variables for the climate-specific RGN
   integer   :: clim_seed=0
   integer   :: rng_seed
   logical   :: reset_clim=.true.

   private      clim_seed, rng_seed, reset_clim

contains

   function urand(lb,ub,seed)
      real                            :: urand
      real, optional, intent(in)      :: lb,ub
      real, optional, intent(in)      :: seed

      integer, dimension(1)           :: seeder
      integer                         :: iseed
      real                            :: rnd
      real                            :: lower,upper

      ! Returns a uniformly-distributed random number in the range lb to ub.
      if ( present(seed) ) then
         iseed=int(seed)
         call set_random_seed(iseed)
      endif

      if ( present(lb) ) then
         lower=lb
      else
         lower=0.0
      endif

      if ( present(ub) ) then
         upper = ub
      else
         upper = 1.0
      endif

      call random_number(rnd)

      urand = lower+(upper-lower)*rnd

      return
   end function urand


   function nrand(mean,std,seed)
      real                             :: nrand
      real, optional, intent(in)       :: mean,std
      real, optional, intent(in)       :: seed

      real                             :: rnd
      integer                          :: iseed
      real                             :: mn, st
      real                             :: x1, y1, x2, y2, w

      ! Returns a normally-distributed random number with given mean and standard dev.
      if ( present(mean) ) then
         mn = mean
      else
         mn = 0.0
      endif

      if ( present(std) ) then
         st = std
      else
         st = 1.0
      endif

      ! Box-Muller polar method
      1 continue

      if ( present(seed) ) then
         iseed=int(seed)
         x1 = urand(-1.0,1.0,seed)
         x2 = urand(-1.0,1.0,seed)
      else
         x1 = urand(-1.0,1.0)
         x2 = urand(-1.0,1.0)
      endif

      w = x1**2 + x2**2
      if ( w .eq. 0.0 .or. w .gt. 1.0 ) go to 1
      w = sqrt( (-2.0 * log( w ) ) / w )
      y1 = x1*w
      y2 = x2*w

      ! Pick one, adjust its mean and std
      nrand = y1*st + mn

      return
   end function nrand


   function clim_urand(lb,ub)
      real                             :: clim_urand
      real, optional, intent(in)       :: lb,ub

      real                             :: lower,upper
      real                             :: ran1
      real                             :: rm1, rm2
      real, dimension(97)              :: r
      integer, parameter               :: m1=259200,ia1=7141,ic1=54773
      integer, parameter               :: m2=134456,ia2=8121,ic2=28411
      integer, parameter               :: m3=243000,ia3=4561,ic3=51349
      integer                          :: ix1, ix2, ix3
      integer                          :: j
      logical                          :: first=.true.
      save

      ! Used in order to have a separate RNG for climate data.
      ! Simple-minded but decent linear congruence with a shuffle, from the 
      ! Yellow Peril (Numerical Recipes First Edition).  Passes "Runs" test
      ! but probably not DIEHARD.

      if ( present(lb) ) then
         lower=lb
      else
         lower=0.0
      endif

      if ( present(ub) ) then
         upper = ub
      else
         upper = 1.0
      endif

      rm1=1./m1
      rm2=1./m2

      if (first .or. reset_clim) then
         ix1=mod(abs(ic1-clim_seed),m1)
         ix1=mod(ia1*ix1+ic1,m1)
         ix2=mod(ix1,m2)
         ix1=mod(ia1*ix1+ic1,m1)
         ix3=mod(ix1,m3)
         do j=1,97
            ix1=mod(ia1*ix1+ic1,m1)
            ix2=mod(ia2*ix2+ic2,m2)
            r(j)=(real(ix1)+real(ix2)*rm2)*rm1
         enddo
         first=.false.
         reset_clim=.false.
      endif

      ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ia2*ix2+ic2,m2)
      ix3=mod(ia3*ix3+ic3,m3)
      j=1+(97*ix3)/m3
      if ( j .gt. 97 .or. j .lt. 1 ) then
         j=97
         write(*,*) 'Error in climate RNG'
      else
         ran1=r(j)
         r(j)=(real(ix1)+real(ix2)*rm2)*rm1
      endif

      clim_urand=lower+(upper-lower)*ran1
   
      return
   end function clim_urand


   function clim_nrand(mean,std)
      real                             :: clim_nrand
      real, optional                   :: mean,std

      real                             :: rnd

      real                             :: mn, st
      real                             :: x1, y1, x2, y2, w

      ! Returns a normally-distributed random number with given mean and standard dev.
      if ( present(mean) ) then
         mn = mean
      else
         mn = 0.0
      endif

      if ( present(std) ) then
         st = std
      else
         st = 1.0
      endif

      ! Box-Muller polar method
      1 continue

      x1=clim_urand(-1.0,1.0)
      x2=clim_urand(-1.0,1.0)

      w = x1**2 + x2**2
      if ( w .eq. 0.0 .or. w .gt. 1.0 ) go to 1

      w = sqrt( (-2.0 * log( w ) ) / w )
      y1 = x1*w
      y2 = x2*w
      ! Pick one, adjust its mean and std
      clim_nrand = y1*st+mn

      return
   end function clim_nrand


   subroutine set_site_rng_seed(fixed_seed,seed)
      logical,        intent(in)      :: fixed_seed
      real, optional, intent(in)      :: seed

      integer                         :: idate(8)
      integer                         :: iseed
      integer, parameter              :: default_seed=2345678
      logical                         :: first=.true.
      real                            :: rnd

      ! Fixed seed is for debugging.  If true, set the seed to be the default
      ! seed each it is called.  If false, get a random seed (or optionally use 
      ! a specified seed) on the first call.
      !
      ! If fixed_seed is true, reset the seed for each site.  This is done
      ! because if we ever parallelize the code, the only way to ensure 
      ! reproducibility is to start each site from the same state. It's not
      ! so important for serial runs but it doesn't take all that much time.
      ! Even for serial runs, however, it means we can run the sites in any order.

      if ( .not. fixed_seed ) then
         if ( first ) then
            if ( .not. present(seed) ) then
               call date_and_time(values=idate)
               call get_random_seed(iseed)
               ! idate(8) contains millisecond
               if ( iseed .ne. 0 ) then
                  iseed = iseed * (idate(8)) 
               else
                  iseed = default_seed * (idate(8)) 
               endif
            else
               iseed=seed
            endif

            first=.false.
            rng_seed=iseed
            call set_random_seed(rng_seed)

         endif

      else

         if ( .not. present(seed) ) then
            iseed=default_seed
         else
            iseed=seed
         endif

         rng_seed=iseed

         call set_random_seed(rng_seed)
         call random_number(rnd)

      endif

   end subroutine set_site_rng_seed


   subroutine set_climate_rng_seed(same_climate,fixed_seed,seed)
      logical,           intent(in) :: same_climate,fixed_seed
      integer, optional, intent(in) :: seed

      integer                       :: idate(8)
      integer                       :: default_seed=2345678
      integer                       :: iseed
      logical                       :: first=.true.

      ! If not using fixed seed, initialize with some randomness. 
      ! Otherwise use the default seeds.
      ! For the same climate we only set the seed on the first call, then we save
      ! it and re-initialize on later calls.
      ! As with the site-related RNG, for the fixed-seed case (debugging) we
      ! reinitialize for each new site to be independent of site ordering.

      if ( .not. fixed_seed ) then
            if ( first ) then
               if ( .not. present(seed) ) then
                  call date_and_time(values=idate)
                  call get_random_seed(iseed)
                  ! idate(7) is thousands of seconds, idate(8) contains milliseconds
                  iseed = idate(7) * (idate(8))      
               else
                  iseed=seed
               endif

               first=.false.

               clim_seed=iseed
               reset_clim=.false.

            else if ( same_climate ) then
               reset_clim=.true.

            endif

      else

            if ( .not. present(seed) ) then
               iseed=default_seed
            else
               iseed=seed
            endif

            clim_seed=iseed

            if ( same_climate ) then
               reset_clim=.true.
            endif

      endif

   end subroutine set_climate_rng_seed


   subroutine get_random_seed(seed)
      integer, intent(out)                :: seed

      integer                             :: isize
      integer, dimension(:), allocatable  :: iseed

      call random_seed(size=isize)
      if (.not. allocated(iseed)) allocate(iseed(isize))

      iseed=seed
      call random_seed(get=iseed)
      seed=iseed(1)

   end subroutine get_random_seed


   subroutine set_random_seed(seed)
      integer, intent(in)                 :: seed

      integer                             :: isize
      integer, dimension(:), allocatable  :: iseed

      call random_seed(size=isize)
      if (.not. allocated(iseed)) allocate(iseed(isize))

      iseed=seed
      call random_seed(put=iseed)

   end subroutine set_random_seed


end module Random