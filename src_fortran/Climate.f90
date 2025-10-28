module Climate

   use Constants
   use Random

   implicit none

   real                              :: accumulated_tmin, accumulated_tmax 
   real, dimension(12)               :: accumulated_precip

contains

   subroutine set_site_climate(same_climate,fixed_seed)
      logical :: same_climate,fixed_seed

      call set_climate_rng_seed(same_climate,fixed_seed)

      accumulated_tmin=0.0
      accumulated_tmax=0.0
      accumulated_precip=0.0

   end subroutine set_site_climate


   !cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
   !c
   !c	cov365: convert monthly state data into daily data
   !c	by Yan Xiaodong 1996 Ver 1.0
   !c
   !c----------------------------------------------------------------------
   subroutine cov365(ta1,vta)
      real, dimension(:), intent(in)    :: ta1
      real, dimension(:), intent(inout) :: vta               

      real, dimension(13)               :: ta
      real, dimension(381)              :: vt
      real                              :: tyxd
      integer :: ltmt(13), k, md
      data ltmt/16,45,75,105,136,166,196,227,258,288,319,349,381/

      ta(13)=ta1(1)
      do k=1,12
         ta(k)=ta1(k)
      end do
      do k=1,12
         tyxd=(ta(k+1)-ta(k))/float(ltmt(k+1)-  ltmt(k))
         do md=ltmt(k),ltmt(k+1)
            vt(md)=ta(k)+tyxd*float(md- ltmt(k))
         end do
      end do
      do md=16,365
         vta(md)=vt(md)
      end do
      do md=1,15
         vta(md)=vt(365+md)
      end do
      return
   end subroutine cov365

   !cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
   !c
   !c	cov365a: convert monthly integreted data into daily data randomly
   !c	by Yan Xiaodong 1996 Ver 1.0
   !c	number of rain days is function of rainfall amount(cm)
   !c	basic feature is: 100cm rain......25 days
   !c					  	1cm rain......1 day
   !c
   !c-----------------------------------------------------------------------	 
   subroutine cov365a(ta,vta)
      real, dimension(:), intent(in)    :: ta
      real, dimension(:), intent(inout) :: vta               

      real    :: raindays, rr, ss, yxdr
      integer :: ltmt(12), md,i, k, ik, inum
      data ltmt/31,28,31,30,31,30,31,31,30,31,30,31/
      md=0
      do k=1,12
         raindays=min1(25.0,ta(k)/4.0+1.0)
         ik=int(raindays)
         rr=ta(k)/float(ik)
         ss=raindays/ltmt(k)
         inum=ik
         do i=1, ltmt(k)
            md=md+1
            if (inum.gt.0) then
               ! old function
               !yxdr=uniform()
               yxdr=clim_urand()
               if(yxdr.le.ss) then 
                  vta(md)=rr
                  inum=inum-1
               else
                  vta(md)=0.0
               end if
            else
               vta(md)=0.0
            end if
         end do
         if (inum.gt.0) then
            vta(md-15)=float(inum)*rr
         end if
      end do
      return
   end subroutine cov365a


   !<----------------------------------------------------------------------------->
   !
   !		ex_rad: Extraterrestial Radiation
   !		By Yan Xiaodong  1996 version 1.0
   !		Input variables:	julia:	Julian Days Number
   !					latit: 	Latitude (South: "-") degree o
   !		Output variables:	erad: 	daily radiation (mj/m2/day)
   !					daylength: light dady length (hours)
   !					exradmx:noon radiation (mj/m2/min)
   !
   !<----------------------------------------------------------------------------->	
   subroutine ex_rad(julia,latit,erad,daylength,exradmx)
      integer,  intent(in) :: julia
      real,     intent(in) :: latit
      real,     intent(out):: erad,daylength,exradmx

      real    :: yxd,rlat,dr,dairta,omega
      !
      !dr=1+0.033cos(2PAI/365*julia)
      !dairta=0.409sin(2PAI/365*julia-1.39)
      !omega=arccos(-tan(latit)*tan(dairta))
      !erad=Gsc*dr*24*60/PAI*cos(latit)*cos(dairta)*(sin(omega)-omega*cos(omega))
      !Gsc=0.0820 (mj/m2/min)

      rlat=deg2rad*latit

      dr=1.0+Ac*cos(b*float(julia))
      dairta=As*sin(b*float(julia)+phase)
      yxd=-tan(rlat)*tan(dairta)

      if (yxd .ge. 1.0) then
         omega=0.0
      else if (yxd .le. -1.0) then
         omega= pi
      else
         omega=acos(yxd)
      end if

      erad=Amp*cos(rlat)*cos(dairta)*(sin(omega)-omega*cos(omega))
      daylength=dl_omega*omega
      exradmx=exrad_coef*dr*cos(rlat-dairta)

      return

   end subroutine ex_rad

   !<+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++>	

   !<----------------------------------------------------------------------------->	
   !
   !  Hargreaves Evaporation formulation
   !  by Yan Xiaodong 1996 version 1.0
   !  hargrea(tmin,tmax,exrad): function output (cm/day)
   !  tmin: minimum temperature degree oC
   !  tmax: maximum temperature degree oC
   !  exrad: extraterrestrial rediation (mj/m2/day)
   !  hargrea= 0.00023/2.45*(tmax-tmin)^0.5*(tmean+17.8)*exrad
   !         = 0.000093876*(tmax-tmin)^0.5*(tmean+17.8)*exrad
   !  I think we can model the potential evaporation more exactly by factoring air
   !  humidity (%): 
   !               100%     *1/10
   !	              90%     *3/10
   !                80%	    *5/10
   !                70%     *7/10
   !                60%	    *8/10
   !                50%     *9/10
   !                40%     *10/10
   !                30%     *11/10
   !                20%     *12/10
   !                10%     *13/10
   !                 0%     *13/10 
   !<----------------------------------------------------------------------------->	
   function hargrea(tmin,tmax,ta,exrad)
      real              :: hargrea
      real, intent(in)  :: tmin, tmax, ta, exrad

      if (ta .le. 0.0) then 
         hargrea=0.0
      else
         hargrea= H_coeff*(tmax - tmin)**0.5 * (ta + H_addon ) * exrad
      end if

   end function hargrea


end module Climate