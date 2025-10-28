module Model

  use Parameters
  use Constants
  use Soil
  use Site
  use Species
  use Tree
  use Random
  use Climate
  use Input

  implicit none

  real, parameter                            :: grow_min=0.05
  real, parameter                            :: growth_thresh=0.05
  real                                       :: growth_min=0.01


contains


   subroutine BioGeoClimate(site,year)

      integer,                         intent(in)    :: year
      type(SiteData),                  intent(inout) :: site

      integer                        :: wmo
      integer                        :: gcm_year
      integer                        :: num_species

      real, dimension(NTEMPS)        :: t1, t2, p1
      real, dimension(NTEMPS)        :: tmptmin, tmptmax, tmpprec
      real, dimension(days_per_year) :: daytemp, daytemp_min, daytemp_max
      real, dimension(days_per_year) :: daynums, dayprecip

      real                   :: litter_c_lev1,litter_c_lev2
      real                   :: litter_n_lev1,litter_n_lev2
      real                   :: rain,rain_n,freeze
      real                   :: temp_f,prcp_f
      real                   :: total_rsp,avail_n,n_avail,C_resp,pet,aet
      real                   :: growdays,drydays_upper,drydays_base
      real                   :: flooddays,degday
      real                   :: outwater
      real                   :: exrad,daylength,exradmx
      real                   :: pot_ev_day
      real                   :: act_ev_day
      real                   :: laiw0_ScaledByMax,laiw0_ScaledByMin
      real                   :: aow0_ScaledByMax,aow0_ScaledByMin
      real                   :: sbw0_ScaledByMax,sbw0_ScaledByMin
      real                   :: saw0_ScaledByFC,saw0_ScaledByWP
      ! KAH -- POSSIBLE BUG
      real                   :: yxd3   !used but never set
      ! used to temporarily hold accumulated climate variables
      real                   :: tmpstep1,tmpstep2
      real                   :: tmp
      integer                :: i,j,k,m

      real, parameter        :: min_grow_temp =5.0
      real, parameter        :: max_dry_parm  =1.0001
      real, parameter        :: min_flood_parm=0.9999

      save   

      wmo=site%site_id
      num_species=size(site%species)

      rain  =0.0
      rain_n=0.0

      ! The user is expected to input decr_by values as positive.

      if ( linear_cc ) then
         if (year .ge. begin_change_year .AND. year .le.                         &
            (begin_change_year + duration_of_change)) then
            accumulated_tmin = accumulated_tmin + tmin_change
            accumulated_tmax = accumulated_tmax + tmax_change
            do m=1,12
               tmpstep1 = site%precip(m) + accumulated_precip(m)
               tmpstep2 = tmpstep1 * precip_change
               accumulated_precip(m) = accumulated_precip(m) + tmpstep2
            end do
         endif

      else if ( use_gcm ) then
         gcm_year=start_gcm+year-begin_change_year
         if ( gcm_year .ge. start_gcm .and. gcm_year .le. end_gcm ) then
            call read_gcm_climate(gcm_year,start_gcm,site)
         endif
      endif

      do i=1,NTEMPS

         if ( linear_cc ) then
            tmptmin(i) = site%tmin(i)   + accumulated_tmin
            tmptmax(i) = site%tmax(i)   + accumulated_tmax
            tmpprec(i) = site%precip(i) + accumulated_precip(i)
         else
            tmptmin(i) = site%tmin(i) 
            tmptmax(i) = site%tmax(i) 
            tmpprec(i) = site%precip(i) 
         endif

         ! Climate fluctuations
         temp_f=clim_nrand(0.0,1.0)
         prcp_f=clim_nrand(0.0,1.0)
         temp_f=max(-1.0,min(temp_f,1.0))

         !--------------------------------------------------------------
         !-!!!!!!!!!!!!!!!Note this part-!!!!!!!!!!!!!!!!!!!!!1
         !--------------------------------------------------------------
         ! KAH Both these can't operate. The second one overwrites the first one.
         ! This doesn't add on fluctuations that are scaled to the std.
         t1(i)=tmptmin(i)+temp_f*site%tmin_std(i)
         t2(i)=tmptmax(i)+temp_f*site%tmax_std(i)
         !t1(i)=tmptmin(i)+temp_f
         !t2(i)=tmptmax(i)+temp_f

         !forest cover can increase rainfall by maximum 15%
         prcp_f=max(-0.5,min(prcp_f,0.5))

         p1(i)=max(tmpprec(i)+prcp_f*site%precip_std(i),0.0)
         rain  =rain  +p1(i)
         rain_n=rain_n+p1(i)*prcp_n

         if ( with_clim_change ) then
         !Maybe write something here
         endif
      end do

      call cov365(t1,daytemp_min)
      call cov365(t2,daytemp_max)
      call cov365a(p1,dayprecip)

      do i=1,days_per_year
         daytemp(i)=0.5*(daytemp_min(i) + daytemp_max(i))
      end do
      !
      !Daily cycles of C, N, H2O
      !
      !Initialize accumulators
      total_rsp=0.0
      avail_n=0.0
      pet=0.0
      aet=0.0
      degday=0.0
      freeze=0.0
      growdays=0.0
      drydays_upper=0.0
      drydays_base=0.0
      flooddays=0.0
      litter_c_lev1=0.0
      litter_n_lev1=0.0
      litter_c_lev2=0.0
      litter_n_lev2=0.0
      outwater=0.0

      site%soil%A0_c0=site%soil%A0_c0+site%soil%C_into_A0
      site%soil%A0_n0=site%soil%A0_n0+site%soil%N_into_A0

      do j=1,days_per_year
         !
         !Water in each Soil layer
         !
         call ex_rad(j,site%latitude,exrad,daylength,exradmx)

         pot_ev_day=hargrea(daytemp_min(j),daytemp_max(j),daytemp(j),exrad)

         call soil_water(site%soil,site%slope,site%leaf_area_ind,              &
            site%leaf_area_w0,site%sigma,freeze,dayprecip(j),                  &
            pot_ev_day,act_ev_day,                                             &
            laiw0_ScaledByMax,laiw0_ScaledByMin,aow0_ScaledByMax,              &
            aow0_ScaledByMin,sbw0_ScaledByMax,sbw0_ScaledByMin,                &
            saw0_ScaledByFC,saw0_ScaledByWP)

         !C and N in Soil and Available N for tree growth
         !
         call soil_decomp(site%soil,litter_c_lev1,litter_c_lev2,litter_n_lev1,  &
                        litter_n_lev2,daytemp(j),dayprecip(j),aow0_ScaledByMax, &
                        saw0_ScaledByFC,sbw0_ScaledByMax,n_avail,C_resp)

         outwater = outwater+site%soil%runoff
         avail_n = max(n_avail,0.0)+avail_n
         pet = pot_ev_day+pet
         aet = act_ev_day+aet
         total_rsp = total_rsp+C_resp

         !computing degday, growing season length(growdays),dry day, &flood days
         if (daytemp(j) .ge. min_grow_temp) then
            degday=degday+(daytemp(j)- min_grow_temp)
            growdays=growdays+1.0
         end if
         if ((saw0_ScaledByFC .lt. max_dry_parm)  .and.                    &
            (sbw0_ScaledByMin .lt. max_dry_parm)  .and.                     &
            (sbw0_ScaledByMax .lt. max_dry_parm)) then
            drydays_upper=drydays_upper+1.0
         endif
         if (saw0_ScaledByWP.lt. max_dry_parm) then
            drydays_base=drydays_base+1.0
         endif
         ! POSSIBLE BUG -- yxd3 never set
         if (aow0_ScaledByMin .gt. min_flood_parm) then
            flooddays = flooddays + kron(yxd3)
         endif

      end do

      if ( growdays .eq. 0 ) then
         drydays_upper=0
         drydays_base=0
         flooddays=0
      else
         tmp = max(min(rain/pet,1.0),min(aet/pet,1.0))
         drydays_upper = min(drydays_upper/growdays,1.0-tmp)
         drydays_base = drydays_base/growdays
         flooddays = flooddays/growdays
      endif

      do k=1,num_species
         call temp_rsp(site%species(k),degday)
         call drought_rsp(site%species(k),drydays_upper,drydays_base)
         call flood_rsp(site%species(k),flooddays)
      end do

      site%soil%avail_N = avail_n+rain_n
      site%soil%total_C_rsp = total_rsp
      site%soil%runoff = outwater

      site%pot_evap_day = pet
      site%act_evap_day = aet
      site%rain = rain
      site%grow_days = growdays
      site%deg_days = degday
      site%flood_days = flooddays
      site%dry_days_upper_layer = drydays_upper
      site%dry_days_base_layer = drydays_base

   end subroutine BioGeoClimate


   subroutine Canopy(site)

      type(SiteData),              intent(inout) :: site

      real, parameter                            :: xt=-0.40

      real, dimension(maxheight)                 :: lvd_c1,lvd_c2
      real, dimension(maxheight)                 :: lvd_c3,lvd_c4
      real                                       :: forht,canht
      real                                       :: tlai, lvd_adj

      integer                                    :: num_species, ntrees
      integer                                    :: i,ip,ih,is,it
      integer                                    :: k,jtmp
      !temporary for testing
      !real :: clsum=0.0, dlsum=0.0

      site%leaf_area_ind = 0.0
      num_species = size(site%species)

      do ip=1,site%numplots

         ntrees=site%plots(ip)%numtrees

         if (ntrees .eq. 0) then
            site%plots(ip)%con_light=1.0
            site%plots(ip)%dec_light=1.0
            site%plots(ip)%nutrient =1.0

         else

            lvd_c1=0.0
            lvd_c2=0.0
            lvd_c3=0.0
            lvd_c4=0.0

            do it=1,ntrees

               k=site%plots(ip)%trees(it)%species_index

               forht=site%plots(ip)%trees(it)%forska_ht
               canht=site%plots(ip)%trees(it)%canopy_ht
               tlai=lai_biomass_c(site%plots(ip)%trees(it))
               site%leaf_area_ind=site%leaf_area_ind+tlai

               jtmp=max(int(forht)-int(canht)+1,1)
               lvd_adj=tlai/float(jtmp)

               if (site%plots(ip)%trees(it)%conifer) then
                  do ih=int(canht),int(forht)
                     lvd_c1(ih)=lvd_c1(ih)+lvd_adj
                     lvd_c2(ih)=lvd_c2(ih)+lvd_adj
                  end do
               else
                  do ih=int(canht),int(forht)
                     lvd_c1(ih)=lvd_c1(ih)+lvd_adj
                     lvd_c2(ih)=lvd_c2(ih)+lvd_adj*0.8
                  end do
               end if
            end do
            lvd_c3(maxheight)=lvd_c1(maxheight)
            lvd_c4(maxheight)=lvd_c2(maxheight)
            do ih=1,maxheight-1
               lvd_c3(maxheight-ih)=lvd_c3(maxheight-ih+1)+lvd_c1(maxheight-ih)
               lvd_c4(maxheight-ih)=lvd_c4(maxheight-ih+1)+lvd_c2(maxheight-ih)
            end do
            do ih=1,maxheight-1
               site%plots(ip)%dec_light(ih)=exp(xt*lvd_c3(ih+1)/plotsize)
               site%plots(ip)%con_light(ih)=exp(xt*lvd_c4(ih+1)/plotsize)
            end do
         end if
      end do

      site%leaf_area_ind=site%leaf_area_ind/float(site%numplots)/plotsize

   end subroutine Canopy


   subroutine Growth(site)

      type(SiteData),                intent(inout) :: site

      real,    dimension(maxtrees)                 :: N_stress
      real,    dimension(maxtrees)                 :: bleaf
      real,    dimension(maxtrees)                 :: diam, biom_C
      real,    dimension(maxtrees)                 :: forska_shade,forht
      integer, dimension(maxtrees)                 :: khc,kh
      real                                         :: canopy_shade
      real                                         :: leaf_b, leaf_bm
      real                                         :: leafbm,leaf1
      real                                         :: d_diam_max
      real                                         :: d_dd,d_leafb
      real                                         :: diameter
      real                                         :: canht
      real                                         :: fct,fc_n,z,pp,bct
      real                                         :: d_bc,d_bn,dt,d_bioC
      real                                         :: N_used,N_req
      real                                         :: prim_prod
      real                                         :: net_prim_prodC
      real                                         :: net_prim_prodN
      real                                         :: net_C_into_A0
      real                                         :: biomc, biomn
      real                                         :: check
      real                                         :: uconvert
      integer                                      :: num_species
      integer                                      :: ntrees
      integer                                      :: it,is,ip,k


      num_species=size(site%species)

      site%soil%C_into_A0=0.0
      site%soil%N_into_A0=0.0

      N_used=0.0
      net_prim_prodC=0.0
      net_prim_prodN=0.0
      biomc=0.0
      biomn=0.0

      leaf_b = 1.0 + con_leaf_ratio

      do ip=1,site%numplots

         site%plots(ip)%avail_spec=0.0
         ntrees=site%plots(ip)%numtrees
         N_req=0.0

         if (ntrees .gt. 0) then
            N_req=0.0
            do it=1,ntrees

               k=site%plots(ip)%trees(it)%species_index
               call update_tree(site%plots(ip)%trees(it),site%species(k))

               ! convenience variables to reduce table lookups
               diam(it) =site%plots(ip)%trees(it)%diam_bht
               canht    =site%plots(ip)%trees(it)%canopy_ht
               forht(it)=site%plots(ip)%trees(it)%forska_ht

               site%plots(ip)%avail_spec(k)=                                   &
                  max(kron(diam(it)-site%species(k)%max_diam*growth_thresh),   &
                           site%plots(ip)%avail_spec(k))

               khc(it)=int(canht)
               kh(it) =int(forht(it))

               call leaf_biomass_c(site%plots(ip)%trees(it))
               call max_growth(site%plots(ip)%trees(it))
               !
               !shade light effect for updating canopy height
               !
               if (site%plots(ip)%trees(it)%conifer) then
                  canopy_shade=light_rsp(site%species(k),                        &
                     site%plots(ip)%con_light(kh(it)))
                  forska_shade(it)=light_rsp(site%species(k),                    &
                     site%plots(ip)%con_light(khc(it)))
               else
                  canopy_shade=light_rsp(site%species(k),                        &
                     site%plots(ip)%dec_light(kh(it)))
                  forska_shade(it)=light_rsp(site%species(k),                    &
                     site%plots(ip)%dec_light(khc(it)))
               end if
               !
               !computing nitrogen stress
               !
               N_stress(it)=env_stress(site%plots(ip)%trees(it),canopy_shade)

               ! Increment tree diameter with growth and stress factors
               site%plots(ip)%trees(it)%diam_bht=                                &
                              site%plots(ip)%trees(it)%diam_bht+               &
                              site%plots(ip)%trees(it)%diam_max*N_stress(it)

               ! Compute forska height for this diameter
               call forska_height(site%plots(ip)%trees(it))

               ! Stem/trunk and leaves/roots
               ! Save current value of leaf biomass
               bleaf(it)=site%plots(ip)%trees(it)%leaf_bm

               ! Update with new height
               call stem_shape(site%plots(ip)%trees(it))
               call leaf_biomass_c(site%plots(ip)%trees(it))

               ! leaf and fine root growth N requirement
               ! KAH Note that order is reversed from old code -- bleaf is the
               ! previous value, not the current value

               if (site%species(k)%conifer) then
                  N_req=N_req+(leaf_b*site%plots(ip)%trees(it)%leaf_bm-          &
                                                         bleaf(it))/con_leaf_c_n
               else
                  N_req=N_req+site%plots(ip)%trees(it)%leaf_bm/dec_leaf_c_n
               end if

               ! wood growth N requirement
               !Store old value of total biomass in C
               biom_C(it)=site%plots(ip)%trees(it)%biomC
               !Compute new value
               call biomass_c(site%plots(ip)%trees(it))
               call biomass_n(site%plots(ip)%trees(it))
               !Delta adjusted for N
               N_req=N_req+(site%plots(ip)%trees(it)%biomC-biom_C(it))/stem_c_n

            end do  ! first tree loop

            N_req=max(N_req*hec_to_m2/plotsize,0.00001)
            N_req=site%soil%avail_N/N_req

            do is=1,num_species
               site%plots(ip)%nutrient(is)=                                    &
                  poor_soil_rsp(N_req,site%species(is)%lownutr_tol)
            end do

            !Actual dbh & other dimension increment, net_prim_prod, and N used
            do it=1,ntrees

               k=site%plots(ip)%trees(it)%species_index
               call update_tree(site%plots(ip)%trees(it),site%species(k))

               fc_n=N_stress(it)*site%plots(ip)%nutrient(k)
               dt=fc_n*site%plots(ip)%trees(it)%diam_max

               ! Increment old diameter, store as current tree diameter
               site%plots(ip)%trees(it)%diam_bht=diam(it)+dt

               pp=min(site%species(k)%max_diam/site%species(k)%max_age*0.1,      &
                                                            growth_thresh)

               if ( (dt .le. pp) .or. (fc_n .le. growth_thresh)) then
                  site%plots(ip)%trees(it)%mort_marker=.true.
               else
                  site%plots(ip)%trees(it)%mort_marker=.false.
               end if

               !c	computing net_prim_prod and N used

               ! compute actual height and diameter w/o intermediate adjustments
               call forska_height(site%plots(ip)%trees(it))
               call stem_shape(site%plots(ip)%trees(it))

               !Update biomass, saving value into convenience variable
               call leaf_biomass_c(site%plots(ip)%trees(it))
               leafbm=site%plots(ip)%trees(it)%leaf_bm

               call biomass_c(site%plots(ip)%trees(it))
               call biomass_n(site%plots(ip)%trees(it))

               !Delta biomass
               d_bioC=site%plots(ip)%trees(it)%biomC - biom_C(it)

               !CN used
               net_prim_prodC=net_prim_prodC+d_bioC
               N_used=N_used+d_bioC/stem_c_n 

               if (site%plots(ip)%trees(it)%conifer) then
                  prim_prod=leaf_b*leafbm-bleaf(it)
                  net_prim_prodC=net_prim_prodC+prim_prod
                  N_used=N_used+prim_prod/con_leaf_c_n

                  !Accumulate total biomasses
                  biomc=biomc+site%plots(ip)%trees(it)%biomC+leafbm
                  biomn=biomn+site%plots(ip)%trees(it)%biomN+leafbm/con_leaf_c_n

               else
               
                  net_prim_prodC=net_prim_prodC+leafbm
                  N_used=N_used+leafbm/dec_leaf_c_n

                  ! Here we do not adjust for leaf biomass, we do that only for
                  ! conifers
                  biomc=biomc+site%plots(ip)%trees(it)%biomC
                  biomn=biomn+site%plots(ip)%trees(it)%biomN

               end if

            end do  ! second tree loop

            !  computing canopy height and wood litter fall
            do it=1,ntrees

               k=site%plots(ip)%trees(it)%species_index
               call update_tree(site%plots(ip)%trees(it),site%species(k))
               forht(it)=site%plots(ip)%trees(it)%forska_ht

               check=site%species(k)%fc_degday*                                  &
                     site%species(k)%fc_drought*                                 &
                     site%species(k)%fc_flood*forska_shade(it)*                  &
                     site%plots(ip)%nutrient(k)

               if (check .le. growth_thresh) then 
                  khc(it)=khc(it)+1
                  if (khc(it) .lt.  int(forht(it))) then

                     ! Increment canopy height 
                     site%plots(ip)%trees(it)%canopy_ht=float(khc(it))+0.01

                     call stem_shape(site%plots(ip)%trees(it))
                     ! save old biomass in C
                     bct=site%plots(ip)%trees(it)%biomC
                     ! update biomass in C and N
                     call biomass_c(site%plots(ip)%trees(it))
                     call biomass_n(site%plots(ip)%trees(it))

                     d_bc=bct-site%plots(ip)%trees(it)%biomC

                     site%soil%C_into_A0=site%soil%C_into_A0+d_bc
                     site%soil%N_into_A0=site%soil%N_into_A0+d_bc/stem_c_n
                     
                     !***************
                     net_prim_prodC = net_prim_prodC - d_bc
                     !N_used = N_used - d_bc/stem_c_n
                     
                     !Previous value of leaf biomass
                     leafbm=site%plots(ip)%trees(it)%leaf_bm
                     !Update again and get difference (note: another reversal)
                     call leaf_biomass_c(site%plots(ip)%trees(it))
                     d_leafb=leafbm-site%plots(ip)%trees(it)%leaf_bm

                     if (site%species(k)%conifer) then
                        site%soil%C_into_A0=site%soil%C_into_A0+d_leafb*leaf_b
                        site%soil%N_into_A0=site%soil%N_into_A0+                 &
                                                   d_leafb/con_leaf_c_n*leaf_b
                        !!!!!!!!!!!!!!!!!
                        net_prim_prodC=net_prim_prodC-d_leafb*leaf_b
                     !  N_used=N_used-d_leafb/con_leaf_c_n*leaf_b
                     else
                        site%soil%C_into_A0=site%soil%C_into_A0+d_leafb
                        site%soil%N_into_A0=site%soil%N_into_A0+d_leafb/         &
                                                            dec_leaf_c_n
                        !!!!!!!!!!!!!!!!
                        net_prim_prodC=net_prim_prodC-d_leafb
                     ! N_used=N_used-d_leafb/con_leaf_c_n
                     end if
                  end if
               end if

            end do  !third tree loop

         end if

      end do

      uconvert=hec_to_m2/plotsize/site%numplots
      N_used=N_used*uconvert
      site%soil%biomc=biomc*uconvert
      site%soil%biomn=biomn*uconvert
      !!!!!!!!!!!!!!!!!!!!
      site%soil%net_prim_prodC=net_prim_prodC*uconvert
      
      site%soil%net_C_into_A0=site%soil%C_into_A0*uconvert
      site%soil%N_used=N_used
      site%soil%net_prim_prodN=N_used
      site%soil%avail_N=site%soil%avail_N - site%soil%N_used

   end subroutine Growth


   subroutine Mortality(site)

      type(SiteData),        intent(inout) :: site

      real                                 :: biomc,biomn
      real                                 :: leaf_b
      real                                 :: fire_prob,wind_prob
      real                                 :: leaf_bm,bmc
      real                                 :: zc,zn,dd,tmp
      real                                 :: uconvert
      real                                 :: NPP_loss,NPPn_loss
      real                                 :: net_prim_prodC,net_prim_prodN
      integer                              :: num_species
      integer                              :: jt,k          
      integer                              :: ip,it,is

      num_species=size(site%species)

      leaf_b= 1.0 + con_leaf_ratio
      biomc = 0.0
      biomn = 0.0
      !*****
      NPP_loss  = 0.0
      NPPn_loss = 0.0

      do ip=1,site%numplots

         fire_prob = urand()
         wind_prob = urand()

         if ( fire_prob < site%fire_prob .or. wind_prob < site%wind_prob ) then

            if (site%plots(ip)%numtrees > 0) then

               if ( fire_prob < site%fire_prob) then 

                  site%plots(ip)%fire=5

                  do is=1,num_species

                     call fire_rsp(site%species(is),1)
                     site%plots(ip)%seedling(is)=(site%species(is)%invader*10.0 +  &
                                                  site%species(is)%sprout_num *    &
                                                  site%plots(ip)%avail_spec(is)) * &
                                                  site%species(is)%fc_fire
                  end do

                  site%plots(ip)%wind=0

               else

                  site%plots(ip)%wind=3

                  do is=1,num_species
                     site%plots(ip)%seedling(is)=site%species(is)%invader+       &
                                    site%plots(ip)%seedling(is)+                 &
                                    site%species(is)%sprout_num*                 &
                                    site%plots(ip)%avail_spec(is)
                  end do

                  site%plots(ip)%fire=0

               end if

               zc=0.0
               zn=0.0

               do it=1,site%plots(ip)%numtrees

                  k=site%plots(ip)%trees(it)%species_index
                  call update_tree(site%plots(ip)%trees(it),site%species(k))

                  tmp=lai_biomass_c(site%plots(ip)%trees(it))*                   &
                                    site%species(k)%leafarea_c*2.0

                  if (site%species(k)%conifer) then 
                     zc=zc+site%plots(ip)%trees(it)%biomC+tmp*leaf_b
                     zn=zn+site%plots(ip)%trees(it)%biomC/stem_c_n               &
                        +tmp/con_leaf_c_n*leaf_b
                  else
                     zc=zc+site%plots(ip)%trees(it)%biomC+tmp
                     zn=zn+site%plots(ip)%trees(it)%biomC/stem_c_n               &
                                                         +tmp/dec_leaf_c_n

                  end if

               end do

               site%soil%C_into_A0 = site%soil%C_into_A0 + zc
               site%soil%N_into_A0 = site%soil%N_into_A0 + zn
               biomc = zc + biomc
               biomn = zn + biomn
               !******
               NPP_loss  = NPP_loss  + zc
               NPPn_loss = NPPn_loss + zn
               
               site%plots(ip)%seedling_number = 1.0

            end if

            !All trees dead now
            site%plots(ip)%numtrees=0

         else

            jt=0

            if (site%plots(ip)%numtrees > 0) then

               site%plots(ip)%fire=0
               site%plots(ip)%wind=0

               do it=1,site%plots(ip)%numtrees

                  k=site%plots(ip)%trees(it)%species_index
                  call update_tree(site%plots(ip)%trees(it),site%species(k))

                  ! I get really tired of trying to come up with good variables
                  ! leaf_bm is a convenience variable
                  call leaf_biomass_c(site%plots(ip)%trees(it))
                  leaf_bm=site%plots(ip)%trees(it)%leaf_bm

                  if (growth_survival(site%plots(ip)%trees(it)) .and.            &
                     age_survival(site%plots(ip)%trees(it)))   then

                     !Tree survives, copy its attributes from old list to new
                     !We start at 1 so at worst we copy 1.. to themselves.  jt
                     !should always be <=it
                     jt=jt+1
                     ! site%plots(ip)%trees(jt)=site%plots(ip)%trees(it)
                     ! assignment overloading doesn't quite work yet
                     call copy_tree(site%plots(ip)%trees(jt),&
                                    site%plots(ip)%trees(it))

                     if (site%species(k)%conifer) then
                        site%soil%C_into_A0 = leaf_bm*(leaf_b-1.0)+                &
                                                            site%soil%C_into_A0
                        site%soil%N_into_A0 = site%soil%N_into_A0+                 &
                                                leaf_bm*(leaf_b-1.0)/con_leaf_c_n
                        !!!!!!!
                        NPP_loss  = NPP_loss + leaf_bm*(leaf_b-1.0)
                        NPPn_loss = NPPn_loss+ leaf_bm*(leaf_b-1.0)/con_leaf_c_n
                        
                     else
                        site%soil%C_into_A0 = leaf_bm+site%soil%C_into_A0
                        site%soil%N_into_A0 = site%soil%N_into_A0+                 &
                                                            leaf_bm/dec_leaf_c_n
                        !!!!!
                        NPP_loss  = NPP_loss + leaf_bm
                        NPPn_loss = NPPn_loss + leaf_bm/dec_leaf_c_n
                        
                     end if

                  else
                     !Tree dies, don't copy it
                     bmc=site%plots(ip)%trees(it)%biomC

                     if (site%species(k)%conifer) then 
                        site%soil%C_into_A0 = site%soil%C_into_A0+bmc+             &
                                                            leaf_bm*leaf_b
                        site%soil%N_into_A0 = site%soil%N_into_A0+bmc/stem_c_n +   &
                                                leaf_bm/con_leaf_c_n*leaf_b
                        !!!!!!!
                        NPP_loss  = NPP_loss + leaf_bm*(leaf_b)+bmc
                        NPPn_loss = NPPn_loss + leaf_bm*(leaf_b)/con_leaf_c_n + bmc/stem_c_n                        
                        
                     else
                        site%soil%C_into_A0=site%soil%C_into_A0+bmc+leaf_bm
                        site%soil%N_into_A0=site%soil%N_into_A0+bmc/stem_c_n+    &
                                                         leaf_bm/dec_leaf_c_n
                        !!!!!!!!!!!!!                                  
                        NPP_loss  = NPP_loss + bmc+leaf_bm 
                        NPPn_loss = NPPn_loss + bmc/stem_c_n + leaf_bm/dec_leaf_c_n 
                                                      
                     end if

                  end if         ! end age/jgrow mortality check

               end do            ! end each tree loop

            endif                ! end tree loop of else clause of fire/wind loop

            site%plots(ip)%numtrees=jt

         end if                  ! end fire/wind loop

      end do                     ! end plot loop

      uconvert = hec_to_m2/plotsize/float(site%numplots)

      site%soil%biomC = site%soil%biomC-biomc*uconvert
      site%soil%biomN = site%soil%biomN-biomn*uconvert
      site%soil%net_prim_prodC = site%soil%net_prim_prodC - NPP_loss*uconvert
      !site%soil%net_prim_prodN = site%soil%net_prim_prodN - NPPn_loss*uconvert

      return

   end subroutine Mortality


   subroutine Renewal(site)

      type(SiteData),              intent(inout) :: site

      real,    dimension(size(site%species))     :: regrowth
      real,    dimension(size(site%species))     :: prob
      real                                       :: net_prim_prod,N_used
      real                                       :: leaf_b
      real                                       :: growmax,grow_cap
      real                                       :: probsum,q0
      real                                       :: z,zz,avtt
      real                                       :: uconvert
      integer                                    :: num_species
      integer                                    :: max_renew
      integer                                    :: new
      integer                                    :: ip,is,it,k
      integer                                    :: nrenew,irenew

      new=0
      net_prim_prod=0.0
      N_used=0.0
      leaf_b=1.0+con_leaf_ratio
      num_species=size(site%species)

      if (site%soil%avail_N .gt. 0.0) then 

         do ip=1,site%numplots

            if ((site%plots(ip)%numtrees .ne. 0) .or.                            &
               ((site%plots(ip)%wind .eq. 0) .and.                               &
               (site%plots(ip)%fire .eq. 0)))                                   &
               then

               growmax=0.0

               do is =1,num_species

                  grow_cap=site%species(is)%fc_degday*                           &
                           site%species(is)%fc_drought*                          &
                           site%species(is)%fc_flood*                            &
                           site%plots(ip)%nutrient(is)

                  if (site%species(is)%conifer) then
                     regrowth(is)=grow_cap*light_rsp(site%species(is),            &
                                       site%plots(ip)%con_light(1))
                  else
                     regrowth(is)=grow_cap*light_rsp(site%species(is),            &
                                       site%plots(ip)%dec_light(1))
                  end if

                  growmax=max(growmax,regrowth(is))
                  if (regrowth(is).le. growth_thresh) regrowth(is)=0.0

               end do

               !Computing the max renew number
               max_renew=min(int(plotsize*growmax)-                              &
                              site%plots(ip)%numtrees,int(plotsize*0.5))

               nrenew=min(max(max_renew,3),int(plotsize)-                        &
                              site%plots(ip)%numtrees)

               !Computing the seed bank size and seedling bank size (1/m^2)
               if (site%plots(ip)%seedling_number .eq. 0.0) then

                  do is =1, num_species
                     site%plots(ip)%seedbank(is) =                               &
                                    site%plots(ip)%seedbank(is) +                 &
                                    site%species(is)%invader +                    &
                                    site%species(is)%seed_num*                    &
                                    site%plots(ip)%avail_spec(is)+                &
                                    site%species(is)%sprout_num*                  &
                                    site%plots(ip)%avail_spec(is)

                     if (regrowth(is) .ge. growth_thresh) then
                        site%plots(ip)%seedling(is)=                             &
                                          site%plots(ip)%seedling(is)+          &
                                          site%plots(ip)%seedbank(is)
                        site%plots(ip)%seedbank(is)=0.0
                     else
                        site%plots(ip)%seedbank(is)=                             &
                                          site%plots(ip)%seedbank(is)*          &
                                          site%species(is)%seed_surv
                     endif
                     call fire_rsp(site%species(is),site%plots(ip)%fire)
                     site%plots(ip)%seedling(is)=                                &
                                       (site%plots(ip)%seedling(is)+             &
                                       site%species(is)%sprout_num*             &
                                       site%plots(ip)%avail_spec(is))*          &
                                       site%species(is)%fc_fire

                     site%plots(ip)%seedling_number=                             &
                                       max(kron(site%plots(ip)%seedling(is)),   &
                                             site%plots(ip)%seedling_number)

                     site%plots(ip)%seedling(is)=                                &
                                       site%plots(ip)%seedling(is)*plotsize
                  end do

                  probsum=0.0

                  do is =1,num_species 
                     prob(is)= site%plots(ip)%seedling(is)*regrowth(is)
                     probsum = probsum+ prob(is)
                  end do

               else

                  probsum=0.0

                  do is =1,num_species 
                     prob(is)=site%plots(ip)%seedling(is)*regrowth(is)
                     probsum = probsum + prob(is)
                  end do

                  do is=1,num_species
                     site%plots(ip)%seedbank(is)=                                &
                              site%plots(ip)%seedbank(is) +                       &
                              site%species(is)%invader+                           &
                              site%species(is)%seed_num*                          &
                              site%plots(ip)%avail_spec(is)+                      &
                              site%species(is)%sprout_num*                        &
                              site%plots(ip)%avail_spec(is)

                     if (regrowth(is) .ge. growth_min) then
                        site%plots(ip)%seedling(is)=site%plots(ip)%seedling(is)+ &
                                                   site%plots(ip)%seedbank(is)
                        site%plots(ip)%seedbank(is)=0.0

                     else
                        site%plots(ip)%seedbank(is)=                             &
                              site%plots(ip)%seedbank(is)*                       &
                              site%species(is)%seed_surv
                     endif

                     call fire_rsp(site%species(is),site%plots(ip)%fire)
                     site%plots(ip)%seedling(is)=                                &
                              (site%plots(ip)%seedling(is) +                     &
                              site%species(is)%sprout_num*                       &
                              site%plots(ip)%avail_spec(is))*                    &
                              site%species(is)%fc_fire

                     site%plots(ip)%seedling(is)=                                &
                              plotsize*site%plots(ip)%seedling(is)

                     site%plots(ip)%seedling_number=                             &
                              max(kron(site%plots(ip)%seedling(is)),             &
                                       site%plots(ip)%seedling_number)
                  end do
               end if
               !c		after disturbances
            else

               if (site%plots(ip)%fire .eq. 1 .or.                               &
                  site%plots(ip)%wind .eq. 1) then

                  probsum=0.0

                  do is =1,num_species 

                     grow_cap = site%species(is)%fc_degday*                      &
                                 site%species(is)%fc_drought*                     &
                                 site%species(is)%fc_flood

                     prob(is) = site%plots(ip)%seedling(is)*grow_cap
                     probsum  = probsum+prob(is)

                     site%plots(ip)%seedling(is)=                                &
                              site%plots(ip)%seedling(is)*plotsize

                     site%plots(ip)%seedling_number=                             &
                              max(kron(site%plots(ip)%seedling(is)),             &
                                       site%plots(ip)%seedling_number)
                  end do

                  site%plots(ip)%fire=0
                  site%plots(ip)%wind=0

               else

                  site%plots(ip)%fire=max(0,site%plots(ip)%fire-1)
                  site%plots(ip)%wind=max(0,site%plots(ip)%wind-1)

               end if
            end if

            !After disturbances									
            if (probsum.gt. epsilon(1.0)) then

               do is =1,num_species
                  prob(is)=prob(is)/probsum
               end do

               do is=2,num_species
                  prob(is)=prob(is-1)+prob(is)
               end do

            else

               nrenew=0

            end if

            it=site%plots(ip)%numtrees

            if (nrenew .ge. 1) then

               do irenew=1, nrenew

                  q0=urand()
                  is = 1

                  do while (q0 .gt. prob(is))

                     is = is + 1
                     if (is .gt. num_species) then
                        is=1+int(urand(0.0,real(num_species)))
                        q0=urand()
                     endif

                  end do

                  !Set new tree feature!

                  new=new+1
                  site%plots(ip)%seedling(is)= &
                                    site%plots(ip)%seedling(is)-1.0
                  it=it+1
                  call initialize_tree(site%plots(ip)%trees(it),                 &
                                                      site%species(is),is)

                  k=is

                  z=1.5 +nrand(0.0,1.0)

                  if (z .ge. 2.5) z=2.5
                  if (z .le. 0.5) z=0.5

                  site%plots(ip)%trees(it)%diam_bht=z
                  site%plots(ip)%trees(it)%canopy_ht=1.0

                  call forska_height(site%plots(ip)%trees(it))
                  call stem_shape(site%plots(ip)%trees(it))
                  call biomass_c(site%plots(ip)%trees(it))
                  call biomass_n(site%plots(ip)%trees(it))
                  call leaf_biomass_c(site%plots(ip)%trees(it))

                  zz=lai_biomass_c(site%plots(ip)%trees(it))*                    &
                                    site%species(k)%leafarea_c*2.0

                  if (site%species(is)%conifer) then 

                     net_prim_prod=net_prim_prod+zz*leaf_b+                      &
                                    site%plots(ip)%trees(it)%biomC

                     N_used=N_used+zz/con_leaf_c_n+                              &
                                    site%plots(ip)%trees(it)%biomN

                     site%soil%C_into_A0=site%soil%C_into_A0+zz*(leaf_b-1.0)
                     site%soil%N_into_A0=site%soil%N_into_A0+                    &
                                             zz*(leaf_b-1.0)/con_leaf_c_n

                  else

                     net_prim_prod=net_prim_prod+                                &
                                    site%plots(ip)%trees(it)%biomC + zz

                     N_used=N_used+site%plots(ip)%trees(it)%biomN        &
                                                               +zz/dec_leaf_c_n

                     site%soil%C_into_A0=site%soil%C_into_A0+zz
                     site%soil%N_into_A0=site%soil%N_into_A0+zz/dec_leaf_c_n

                  end if
               end do
            end if
            site%plots(ip)%numtrees=it

            do is =1, num_species 
               site%plots(ip)%seedling(is)=site%plots(ip)%seedling(is)*          &
                                                site%species(is)%seedling_lg/&
                                                                     plotsize
            end do

            site%plots(ip)%fire=0

         end do

      end if

      uconvert=hec_to_m2/plotsize/float(site%numplots)
      N_used=N_used*uconvert
      net_prim_prod=net_prim_prod*uconvert
      avtt=site%soil%avail_N-N_used
      if (avtt .gt. 0.0) then
         site%soil%net_N_into_A0=avtt*min(site%soil%runoff/1000.0,0.1)
         site%soil%A_n0=site%soil%A_n0+avtt-site%soil%net_N_into_A0
      else
         site%soil%A_n0=site%soil%A_n0+avtt
         site%soil%net_N_into_A0=0.0
      end if

      site%soil%A_n0=site%soil%A_n0-0.00002*site%soil%runoff
      site%soil%A_c0=site%soil%A_c0-site%soil%net_N_into_A0*20.0
      site%soil%BL_c0=site%soil%BL_c0+site%soil%net_N_into_A0*20.0
      site%soil%BL_n0=site%soil%BL_n0+site%soil%net_N_into_A0
      site%soil%C_into_A0=(site%soil%C_into_A0)*uconvert
      site%soil%N_into_A0=(site%soil%N_into_A0)*uconvert
      site%soil%N_used=N_used

      site%soil%net_prim_prodC=site%soil%net_prim_prodC+net_prim_prod
      site%soil%net_prim_prodN=site%soil%net_prim_prodN+N_used
      site%soil%new_growth=int(float(new)*uconvert)

   end subroutine Renewal

end module Model