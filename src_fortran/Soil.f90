module Soil
   ! Author: K. A. Holcomb
   !
   ! Soil attributes and procedures relevant to soil properties
   !
   ! Methods: 
   !    soil_decomp
   !         verified against original 20120724
   !    soil_water
   !         verified against original 20120724
   !

   use Constants
   use Parameters

   implicit none

   ! Global constants

   !  Initial soil values
   !  Units seem to vary

   real       :: ao_cn_0=30.0
   real       :: sa_cn_0=4.0
   real       :: sb_cn_0=20.0
   real       :: ao_resp=5.24e-4
   real       :: sa_resp=1.24e-5
   real       :: sb_resp=2.74e-7

   real       :: soil_base_depth=70.0
   real       :: base_max=0.6
   real       :: base_min=0.1
   real       :: ao_min=0.025
   real       :: ao_max=0.25
   real       :: lai_min=0.01
   real       :: lai_max=0.15

   ! Soil object

   type SoilData
      real    :: A0_c0
      real    :: A_c0
      real    :: A0_n0
      real    :: A_n0
      real    :: A0_w0
      real    :: A_w0
      real    :: A_field_cap
      real    :: A_perm_wp
      real    :: C_into_A0
      real    :: N_into_A0
      real    :: net_C_into_A0
      real    :: net_N_into_A0
      real    :: N_used
      real    :: avail_N
      real    :: BL_c0
      real    :: BL_n0
      real    :: BL_w0
      real    :: base_h
      real    :: biomC, biomN
      real    :: net_prim_prodN
      real    :: net_prim_prodC
      real    :: runoff
      real    :: total_C_rsp
      integer :: new_growth
   end type SoilData

contains

   !3 Methods

   subroutine soil_decomp(soil,litter_c1,litter_c2,litter_n1,litter_n2,         &
                           tempC,precip,aow0_ScaledByMax,saw0_ScaledByFC,        &
                           sbw0_ScaledByMax,avail_N,C_resp)      

      type(SoilData), intent(inout)  ::  soil
      real,           intent(in)     :: litter_c1,litter_c2,litter_n1,litter_n2
      real,           intent(in)     :: tempC,precip
      real,           intent(out)    :: aow0_ScaledByMax
      real,           intent(out)    :: avail_N, C_resp
      ! KAH
      ! This variable is never set or used in this routine, only in soil_water
      real,           intent(in)     :: sbw0_ScaledByMax
      ! This variable is never set or used in this routine, only in soil_water
      real,           intent(in)     :: saw0_ScaledByFC

      real                           :: ao_c0,ao_n0,sa_c0,sa_n0,sb_c0,sb_n0
      real                           :: ao_cn,sa_cn
      real                           :: resp1, resp2, resp3
      real                           :: tadjst,tadjst1,aofunc,yxdn,yxdc,tosb

      save

      !<----------------------------------------------------------------------------->
      !     soil decomposition model for computing the available N and soil
      !     respiration
      !     Original version by Yan Xiaodong
      !     Modified with improved variable names by Katherine Holcomb 2012 January
      !
      !     litter_c1/2: input C as litter(above/under ground): tc/ha/d
      !     litter_n1/2: input N as litter(above/under ground): tn/ha/d
      !     precip          : input water as rain fall cm/d
      !     ao_c0           : input/output ao layer C  tc/ha
      !     ao_n0           : input/output ao layer N  tn/ha
      !     sa_c0           : input/output soil A layer C tc/ha
      !     sa_n0           : input/output soil A layer N tn/ha
      !     sb_c0           : input/output Base Soil layer C tc/ha
      !     sb_n0           : input/output base soil layer N tn/ha
      !     tempC           : input daily temperature degree oC
      !     aow0_ScaledByMax: input aow0/aowmax    cm/cm
      !     saw0_ScaledByFC : input saw0/safc      cm/cm
      !     sbw0_ScaledByMax: input sbw0/sbwmax    cm/cm -- Never used
      !     avail_N         : output  available N for plant growth tn/ha
      !     C_resp          : emission to atmosphere as co2  tc/ha/d
      !<----------------------------------------------------------------------------->

      !Avoiding pain in converting subroutine to OOPier form
      !ao_C0=soil%A0_c0 was changed to ao_c0=soil%A0_c0 by Bin Wang on 3/23/2016
      ao_c0=soil%A0_c0;ao_n0=soil%A0_n0
      sa_c0=soil%A_c0; sa_n0=soil%A_n0
      sb_c0=soil%BL_c0;sb_n0=soil%BL_n0

      !c  this part Ao layer C N balance
      ao_c0=ao_c0+litter_c1
      ao_n0=ao_n0+litter_n1
      ao_cn=ao_c0/ao_n0
      aow0_ScaledByMax=min(aow0_ScaledByMax,0.5)

      aofunc=max((1.0-(1.0-aow0_ScaledByMax/0.3)**2),0.2)

      if (tempC .ge. -5.0) then
         tadjst =3.0**(0.1*(tempC-1.0))
         tadjst1=2.5**(0.1*(tempC-1.0))
      else
         tadjst =0.0
         tadjst1=0.0
      end if
      resp1=tadjst*aofunc*ao_resp*ao_c0
      !c     ^    *eff_c_n(aocn)
      yxdn=resp1/ao_cn
      yxdc=yxdn*ao_cn_0
      ao_c0=ao_c0-yxdc-resp1
      ao_n0=ao_n0-yxdn
      
      !!c  this part soil A layer C N balance
      sa_c0=sa_c0+yxdc+litter_c2
      sa_n0=sa_n0+yxdn+litter_n2
      sa_cn=sa_c0/sa_n0
      !KAH commented out in original
      !saw0_ScaledByFC=min(sa_w0,0.8)

      resp2=tadjst1*max(1.0-(1.0-saw0_ScaledByFC/0.8)**2,0.2)*sa_resp*sa_c0
      !^  *eff_c_n(sacn)
      !
      tosb=resp2/sb_cn_0
      avail_N=resp2/sa_cn*max(0.5,(sa_cn-sa_cn_0)/sa_cn)
      sa_c0=sa_c0-resp2-tosb
      sa_n0=sa_n0-avail_N
      
      !!Base layer C balance
      sb_c0=sb_c0+tosb
      resp3=sb_c0*sb_resp*tadjst1
      sb_c0=sb_c0-resp3
      
      !Total respiration of the 'Three' Layer
      C_resp=resp1+resp2+resp3

      soil%A0_c0=ao_c0;soil%A0_n0=ao_n0
      soil%A_c0=sa_c0;soil%A_n0=sa_n0
      soil%BL_c0=sb_c0;soil%BL_n0=sb_n0
      soil%avail_N=avail_N

   end subroutine soil_decomp


   subroutine soil_water(soil,slope,lai,lai_w0,sigma,freeze,rain,               &
                        pot_ev_day,act_ev_day,laiw0_ScaledByMax,                &
                        laiw0_ScaledByMin,aow0_ScaledByMax,aow0_ScaledByMin,    &
                        sbw0_ScaledByMax,sbw0_ScaledByMin,saw0_ScaledByFC,      &
                        saw0_ScaledByWP)
      !!KAH BUGS?
      !  There were at least two variables originally claimed to be input that were
      !  changed, lai and act_ev_day.  The latter doesn't really matter but the lai
      !  does matter, it affects a field of the soil object.
      !  Also, sbh is supposedly input but it is set to 70

      !   soil water daily cycle model
      !   by yan xiaodong 1997 Version 1.0

      type(SoilData), intent(inout) :: soil
      real,           intent(inout) :: lai,lai_w0
      real,           intent(in)    :: pot_ev_day
      real,           intent(in)    :: slope,sigma,rain,freeze
      real,           intent(out)   :: act_ev_day
      real,           intent(out)   :: laiw0_ScaledByMax,laiw0_ScaledByMin
      real,           intent(out)   :: aow0_ScaledByMax,aow0_ScaledByMin
      real,           intent(out)   :: sbw0_ScaledByMax,sbw0_ScaledByMin
      real,           intent(out)   :: saw0_ScaledByFC,saw0_ScaledByWP

      real                          :: ao_w0,ao,sa_w0,sa_fc
      real                          :: sb_w0,sbw_min,sbw_max
      real                          :: sa_pwp,sbh,petd,runoff,aetd
      real                          :: lossslp,yxd1,yxd2,yxd3,yxd4,yxd,saw,aow,sbw
      real                          :: aow_min,aow_max
      real                          :: lai_loss
      real                          :: laiw,lai_w1,laiw_min,laiw_max
      real                          :: ao_w1,sa_w1,yxd5,sb_w1
      real                          :: table_water

      save

      !<----------------------------------------------------------------------------->
      !       slope       :input slope degree o
      !       lai         :input canopy leaf area index m/m
      !       ao_w0       :input/output ao water content cm
      !       ao          :input ao dry matter amount t/ha
      !       sa_w0       :input/output soil water content cm
      !       sa_fc       :input soil field capacity cm
      !       sa_pwp      :input soil permanent wilt point cm
      !       sb_w0       :input/output base soil water content cm
      !       sbh         :input base soil depth cm
      !       pot_ev_day  :input daily potential evapotraspiration cm/day
      !       rain        :input daily precipitation cm/day
      !       act_ev_day  :output daily actual evapotraspiration cm/day
      !       runoff      :output daily runoff
      !
      !
      !         rain  \/      /\interception loss
      !       ---------------------
      !       |              |   lai 
      !       ----------------
      !              ||
      !              ||                                         aet /\
      !              ||          runoff                             ||
      !              ||       ..>>>>>>                              ||
      !       ++++++++++++++++++++ ao layer water store.......>>>>>>||
      !       --------------------                                  ||
      !                                                             ||
      !                            soil layer water store......>>>>>||
      !                                                             ||
      !       ====================                                  ||
      !                            base soil water store.......>>>>>
      !       ====================
      !
      !<----------------------------------------------------------------------------->

      ! Avoiding too much pain
      ao=soil%A0_c0;ao_w0=soil%A0_w0
      sa_w0=soil%A_w0;sa_fc=soil%A_field_cap;sa_pwp=soil%A_perm_wp
      sb_w0=soil%BL_w0
      runoff=soil%runoff

      act_ev_day=0.0
      lai=max(lai,1.0)
      sbh=soil_base_depth
      laiw_min=lai*lai_min
      laiw_max=lai*lai_max
      aow_min=ao*ao_min
      aow_max=ao*ao_max
      sbw_min=sbh*base_min
      sbw_max=sbh*base_max
      !--------------------------------
      !forest region water balance for underground water table
      if (rain.gt. 0.01) then
         table_water=rain*sigma*freeze
         sb_w0=min(sbw_max,sb_w0+table_water)
      end if
      !-------------------------------
      if (pot_ev_day.le. 0.0) then
         laiw=min(rain+lai_w0, laiw_max)
         yxd1=max(rain-laiw+lai_w0,0.0)
         aow=min(ao_w0+yxd1,aow_max)
         yxd2=max(yxd1-aow+ao_w0,0.0)
         sbw=min(yxd2+sb_w0,sbw_max)
         runoff=max(yxd2-sbw+sb_w0,0.0)
         lai_w0=laiw
         ao_w0=aow
         sb_w0=sbw
      else
         !lai_loss: the amout to be intercepted by the canopy

         !KAH Why are we multiplying rain by 0.00??? BW??
         !lai_loss=min((laiw_max-lai_w0),0.00*rain)
         
         !Modified by BW on 2017-07-13
         lai_loss=min((laiw_max - lai_w0),rain)
         
         !yxd1: the amount to be able to be on the floor
         yxd1=max(rain-lai_loss,0.0)
         laiw=lai_w0+lai_loss
         !lossslp: the amount to be runoff owing to slope
         yxd=(slope/90.0)**2
         lossslp=yxd*yxd1
         !yxd2: available rain amount for soil after slope and interception
         yxd2=yxd1-lossslp-pot_ev_day
         !
         !if available rain above 0.0, the water will be allocated to 
         !different layer by the order: 
         !    soil A layer, saw-saw0; Ao layer, aow-aow0; base soil layer,sbw-sbw0
         !
         if (yxd2.gt. 0.0) then
            saw = min(sa_fc-sa_w0, yxd2) + sa_w0
            yxd3=max(yxd2-(saw-sa_w0),0.0)
            aow =min(aow_max-ao_w0,yxd3) +ao_w0
            yxd4=max(yxd3-(aow-ao_w0),0.0)
            sbw=min(sbw_max-sb_w0,yxd4) +sb_w0
            !actual evaporation=potential evaporation
            act_ev_day=pot_ev_day
            !runoff amount= slope run off + rest part after soil fully charged
            !but why is yxd2 added???? 
            runoff=yxd2+lossslp
            sa_w0=saw
            sb_w0=sbw
            ao_w0=aow
            !once available rain under 0.0, the evaporation will extract the     &
            !water by the order:
            !     canopy,
            !     Ao layer
            !     Soil A layer
            !     Base soil layer
            !
         else
            !
            !extracted from canopy: lai_w1
            lai_w1=min(-yxd2,lai_w0-laiw_min)
            !canopy availabe water: lai_w0
            lai_w0= lai_w0-lai_w1
            act_ev_day=act_ev_day+lai_w1
            !     
            yxd3=min(yxd2+lai_w1,0.0)
            !extracted from ao layer: aow1
            ao_w1= min(-yxd3, ao_w0-aow_min)
            act_ev_day=act_ev_day+ao_w1
            !ao layer available water: aow
            ao_w0=ao_w0-ao_w1
            yxd4=min(yxd3+ao_w1,0.0)
            !extracted from soil A layer: saw1
            sa_w1=min(-yxd4,sa_w0-sa_pwp)
            act_ev_day=act_ev_day+sa_w1
            !soil A layer available water: saw
            sa_w0= sa_w0 - sa_w1
            yxd5= min(yxd4 + sa_w1,0.0)
            !extracted from base soil layer: sb_w1
            sb_w1=min(-yxd5,sb_w0-sbw_min)
            !soil B layer available water: sbw
            sb_w0=sb_w0-sb_w1
            !extracted from slope runoff: tmp
            !actual runoff: runoff
            runoff=lossslp
            !
         end if
      end if

      laiw0_ScaledByMax=lai_w0/laiw_max
      laiw0_ScaledByMin=lai_w0/laiw_min
      aow0_ScaledByMax =ao_w0/aow_max 
      aow0_ScaledByMin =ao_w0/aow_min
      sbw0_ScaledByMax =sb_w0/sbw_max
      sbw0_ScaledByMin =sb_w0/sbw_min
      saw0_ScaledByFC =sa_w0/sa_fc 
      saw0_ScaledByWP =sa_w0/sa_pwp 

      runoff=max(runoff,0.0)

      soil%A0_c0=ao;soil%A0_w0=ao_w0;soil%A_w0=sa_w0;soil%A_field_cap=sa_fc
      soil%A_perm_wp=sa_pwp;soil%BL_w0=sb_w0;soil%base_h=sbh;soil%runoff=runoff

   end subroutine soil_water


   subroutine write_soil_csv(self,soil_unit)
      use csv_file
      type(SoilData), intent(in) :: self
      integer,        intent(in) :: soil_unit

      call csv_write(soil_unit,self%A0_c0,.false.)
      call csv_write(soil_unit,self%A_c0,.false.)
      call csv_write(soil_unit,self%A0_n0,.false.)
      call csv_write(soil_unit,self%A_n0,.false.)
      call csv_write(soil_unit,self%BL_c0,.false.)
      call csv_write(soil_unit,self%BL_n0,.false.)
      call csv_write(soil_unit,self%total_C_rsp,.false.)
      call csv_write(soil_unit,self%biomC,.false.)
      call csv_write(soil_unit,self%C_into_A0,.false.)
      call csv_write(soil_unit,self%net_C_into_A0,.false.)
      call csv_write(soil_unit,self%net_prim_prodC,.false.)
      call csv_write(soil_unit,self%biomN,.false.)
      call csv_write(soil_unit,self%N_into_A0,.false.)
      call csv_write(soil_unit,self%net_N_into_A0,.false.)
      call csv_write(soil_unit,self%net_prim_prodN,.false.)
      call csv_write(soil_unit,self%avail_N,.true.)

   end subroutine write_soil_csv


end module Soil
