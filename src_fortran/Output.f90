module Output

   use Constants
   use Parameters
   use Tree
   use Plot
   use GenusGroups
   use IO
   use csv_file

   implicit none

   ! Plot adjustments
   real                :: plotscale, plotadj, plotrenorm


contains


   subroutine initialize_outputFiles(species_present)
      type(Groups), intent(in)     :: species_present

      call open_outputFiles
      call write_headers(species_present)

   end subroutine initialize_outputFiles


   subroutine write_genus_data(site,species_pres,year)
      type(SiteData),   intent(in)            :: site
      type(Groups),     intent(in)            :: species_pres
      integer,          intent(in)            :: year

      integer, dimension(site%numplots,species_pres%numgenera,NHC) ::            &
                                                diam_categories
      real, dimension(site%numplots,species_pres%numgenera)        ::            &
                                                basal_area,                     &
                                                biomC, biomN, leaf_bm,          &
                                                biomC_std,biomN_std,            &
                                                max_ht,max_diam,                &
                                                plotlevel_biomc,                &
                                                plotlevel_biomn
      integer, dimension(species_pres%numgenera,NHC)             ::              &
                                                total_dm_cats
      real, dimension(species_pres%numgenera)                    ::              &
                                                total_biomc,total_biomn,        &
                                                total_basal,total_laibm,        &
                                                total_max_ht,total_max_diam,    &
                                                plot_mean_biomc,                &
                                                plot_std_biomc,                 &
                                                plot_mean_biomn,                &
                                                plot_std_biomn
      integer, dimension(species_pres%numspecies)                 ::             &
                                                spec_index
      logical, dimension(species_pres%numgenera)                  ::             &
                                                in_site

      character(len=8)                        :: field='genus'
      real                                    :: lc_n
      integer                                 :: iwmo
      integer                                 :: ip, is, ns
      logical                                 :: genus_plot_level_data
                                                

      plotscale = hec_to_m2/plotsize
      plotadj   = plotscale/float(site%numplots)
      plotrenorm= 1./plotsize/float(site%numplots)

      iwmo=site%site_id

      in_site=.false.
      do is=1,species_pres%numgenera
      do ip=1,site%numplots
         do ns=1,size(site%plots(ip)%species)
            if ( species_pres%genusgroups(is) .eq. &
                        site%plots(ip)%species(ns)%genus_name ) then
               spec_index(is)=ns
               in_site(is) = .true.
               exit
            endif
         enddo
      enddo
      enddo

      do ip=1,site%numplots
         call sum_over_sg(site%plots(ip),species_pres%genusgroups,field,         &
                        basal_area(ip,:),leaf_bm(ip,:),                          &
                        biomC(ip,:),biomN(ip,:),biomC_std(ip,:),                 &
                        biomN_std(ip,:),max_ht(ip,:),max_diam(ip,:))

         call tree_dm_cats(site%plots(ip),species_pres%genusgroups,field,        &
                                                      diam_categories(ip,:,:))

      enddo

      total_biomc=0.0
      total_biomn=0.0
      total_basal=0.0
      total_laibm=0.0
      total_max_ht  =0.0
      total_max_diam=0.0
      plotlevel_biomc=rnvalid
      plotlevel_biomn=rnvalid

      total_dm_cats=0

      do is=1,species_pres%numgenera
         do ip=1,site%numplots
            if ( in_site(is) ) then
               ns=spec_index(is)
               total_biomc(is)=total_biomc(is)+biomc(ip,is)
               total_biomn(is)=total_biomn(is)+biomn(ip,is)

               total_basal(is)=total_basal(is)+basal_area(ip,is)
               total_laibm(is)=total_laibm(is)+leaf_bm(ip,is)/                   &
                                       (site%plots(ip)%species(ns)%leafarea_c*2.0)

               total_max_ht(is)  =max(total_max_ht(is),max_ht(ip,is))
               total_max_diam(is)=max(total_max_diam(is),max_diam(ip,is))

               total_dm_cats(is,:)=total_dm_cats(is,:)+diam_categories(ip,is,:)
            endif
         enddo

         call stddev(biomc(:,is),plot_mean_biomc(is),plot_std_biomc(is),rnvalid)
         call stddev(biomn(:,is),plot_mean_biomn(is),plot_std_biomn(is),rnvalid)

      enddo

         ! In the original code the standard deviation was not rescaled to the
         ! same scale as the plotlevel_biomc.  I am doing so here but that can be
         ! changed.  The original code didn't print the std_biomn at all.
         plotlevel_biomc=plotscale*biomc
         plot_mean_biomc=plotscale*plot_mean_biomc
         plot_std_biomc=plotscale*plot_std_biomc
         plotlevel_biomn=plotscale*biomn
         plot_mean_biomn=plotscale*plot_mean_biomn
         plot_std_biomn=plotscale*plot_std_biomn

      do is=1,species_pres%numgenera
         call csv_write(biom_by_g,iwmo,.false.)
         call csv_write(biom_by_g,year,.false.)
         call csv_write(biom_by_g,trim(adjustl(species_pres%genusgroups(is))),   &
                                                                     .false.)
         if ( in_site(is) ) then

            total_basal(is)=total_basal(is)*plotrenorm
            total_laibm(is)=total_laibm(is)*plotrenorm
            total_biomc(is)=total_biomc(is)*plotadj
            total_biomn(is)=total_biomn(is)*plotrenorm*10.0
            total_dm_cats(is,:)=total_dm_cats(is,:)*plotadj

            call csv_write(biom_by_g,total_dm_cats(is,1),.false.)
            call csv_write(biom_by_g,total_dm_cats(is,2),.false.)
            call csv_write(biom_by_g,total_dm_cats(is,3),.false.)
            call csv_write(biom_by_g,total_dm_cats(is,4),.false.) 
            call csv_write(biom_by_g,total_dm_cats(is,5),.false.)
            call csv_write(biom_by_g,total_dm_cats(is,6),.false.)
            call csv_write(biom_by_g,total_dm_cats(is,7),.false.)
            call csv_write(biom_by_g,total_max_diam(is),.false.)
            call csv_write(biom_by_g,total_max_ht(is),.false.)
            call csv_write(biom_by_g,total_laibm(is),.false.)
            call csv_write(biom_by_g,total_basal(is),.false.)
            call csv_write(biom_by_g,total_biomc(is),.false.)
            call csv_write(biom_by_g,plot_std_biomc(is),.false.)
            call csv_write(biom_by_g,total_biomn(is),.false.)
            call csv_write(biom_by_g,plot_std_biomn(is),.true.)

         else

            call csv_write(biom_by_g,rnvalid,.false.)
            call csv_write(biom_by_g,rnvalid,.false.)
            call csv_write(biom_by_g,rnvalid,.false.)
            call csv_write(biom_by_g,rnvalid,.false.)
            call csv_write(biom_by_g,rnvalid,.false.)
            call csv_write(biom_by_g,rnvalid,.false.)
            call csv_write(biom_by_g,rnvalid,.false.)
            call csv_write(biom_by_g,rnvalid,.false.)
            call csv_write(biom_by_g,rnvalid,.false.)
            call csv_write(biom_by_g,rnvalid,.false.)
            call csv_write(biom_by_g,rnvalid,.false.)
            call csv_write(biom_by_g,rnvalid,.false.)
            call csv_write(biom_by_g,rnvalid,.false.)
            call csv_write(biom_by_g,rnvalid,.false.)
            call csv_write(biom_by_g,rnvalid,.true.)


         endif

      enddo

      !   For now this is internal since we don't expect to want it often
      genus_plot_level_data=.false.
      if ( genus_plot_level_data ) then

         do ip=1,site%numplots
            do is=1,species_pres%numgenera
               call csv_write(pl_biom_by_g,iwmo,.false.)
               call csv_write(pl_biom_by_g,year,.false.)
               call csv_write(pl_biom_by_g,ip,.false.)

               if ( in_site(is) ) then

                  call csv_write(pl_biom_by_g,diam_categories(ip,is,1),.false.)
                  call csv_write(pl_biom_by_g,diam_categories(ip,is,2),.false.)
                  call csv_write(pl_biom_by_g,diam_categories(ip,is,3),.false.)
                  call csv_write(pl_biom_by_g,diam_categories(ip,is,4),.false.)
                  call csv_write(pl_biom_by_g,diam_categories(ip,is,5),.false.)
                  call csv_write(pl_biom_by_g,diam_categories(ip,is,6),.false.)
                  call csv_write(pl_biom_by_g,diam_categories(ip,is,7),.false.)
                  call csv_write(pl_biom_by_g,max_diam(ip,is),.false.)
                  call csv_write(pl_biom_by_g,max_ht(ip,is),.false.)
                  call csv_write(pl_biom_by_g,leaf_bm(ip,is),.false.)
                  call csv_write(pl_biom_by_g,basal_area(ip,is),.false.)
                  call csv_write(pl_biom_by_g,biomC(ip,is),.false.)
                  call csv_write(pl_biom_by_g,biomC_std(ip,is),.false.)
                  call csv_write(pl_biom_by_g,biomN(ip,is),.false.)
                  call csv_write(pl_biom_by_g,biomN_std(ip,is),.true.)
               else
                  call csv_write(pl_biom_by_g,rnvalid,.false.)
                  call csv_write(pl_biom_by_g,rnvalid,.false.)
                  call csv_write(pl_biom_by_g,rnvalid,.false.)
                  call csv_write(pl_biom_by_g,rnvalid,.false.)
                  call csv_write(pl_biom_by_g,rnvalid,.false.)
                  call csv_write(pl_biom_by_g,rnvalid,.false.)
                  call csv_write(pl_biom_by_g,rnvalid,.false.)
                  call csv_write(pl_biom_by_g,rnvalid,.false.)
                  call csv_write(pl_biom_by_g,rnvalid,.false.)
                  call csv_write(pl_biom_by_g,rnvalid,.false.)
                  call csv_write(pl_biom_by_g,rnvalid,.false.)
                  call csv_write(pl_biom_by_g,rnvalid,.false.)
                  call csv_write(pl_biom_by_g,rnvalid,.false.)
                  call csv_write(pl_biom_by_g,rnvalid,.false.)
                  call csv_write(pl_biom_by_g,rnvalid,.true.)
               endif

            enddo
         enddo

      endif

   end subroutine write_genus_data


   subroutine write_species_data(site,species_pres,year)
      type(SiteData),   intent(in)            :: site
      type(Groups),     intent(in)            :: species_pres
      integer,          intent(in)            :: year

      integer, dimension(site%numplots,species_pres%numspecies,NHC) ::           &
                                                diam_categories
      real, dimension(site%numplots,species_pres%numspecies)      ::             &
                                                basal_area,                     &
                                                biomC, biomN, leaf_bm,          &
                                                biomC_std,biomN_std,            &
                                                max_ht,max_diam,                &
                                                plotlevel_biomc,                &
                                                plotlevel_biomn
      real, dimension(species_pres%numspecies,NHC)                ::             &
                                                total_dm_cats
      real, dimension(species_pres%numspecies)                    ::             &
                                                total_biomc,total_biomn,        &
                                                total_basal,total_laibm,        &
                                                total_max_ht,total_max_diam,    &
                                                plot_mean_biomc,                &
                                                plot_std_biomc,                 &
                                                plot_mean_biomn,                &
                                                plot_std_biomn
      integer, dimension(species_pres%numspecies)                  ::            &
                                                spec_index
      logical, dimension(species_pres%numspecies)                  ::            &
                                                in_site

      character(len=8)                        :: field='species'
      real                                    :: lc_n
      integer                                 :: iwmo
      integer                                 :: ip, is, ns
                                                

      plotscale = hec_to_m2/plotsize !converted from m2 to hecter
      plotadj   = plotscale/float(site%numplots)
      plotrenorm= 1./plotsize/float(site%numplots)

      iwmo=site%site_id

      in_site=.false.
      do is=1,species_pres%numspecies
         do ip=1,site%numplots
            do ns=1,size(site%plots(ip)%species)
               if ( species_pres%spec_names(is,2) .eq. &
                           site%plots(ip)%species(ns)%unique_id ) then
                  in_site(is) = .true.
                  spec_index(is)=ns
                  exit
               endif
            enddo
         enddo
      enddo

      do ip=1,site%numplots
         call sum_over_sg(site%plots(ip),species_pres%spec_names(:,2),field,     &
                        basal_area(ip,:),leaf_bm(ip,:),                        &
                        biomC(ip,:),biomN(ip,:),biomC_std(ip,:),               &
                        biomN_std(ip,:),max_ht(ip,:),max_diam(ip,:))

         call tree_dm_cats(site%plots(ip),species_pres%spec_names(:,2),field,    &
                                                      diam_categories(ip,:,:))
      enddo

      total_biomc=0.0
      total_biomn=0.0
      total_basal=0.0
      total_laibm=0.0
      total_max_ht  =0.0
      total_max_diam=0.0
      plotlevel_biomc=rnvalid
      plotlevel_biomn=rnvalid

      total_dm_cats=0

      do is=1,species_pres%numspecies
         if ( in_site(is) ) then
            ns=spec_index(is)
            do ip=1,site%numplots
               total_biomc(is)=total_biomc(is)+biomc(ip,is)
               total_biomn(is)=total_biomn(is)+biomn(ip,is)

               total_basal(is)=total_basal(is)+basal_area(ip,is)
               total_laibm(is)=total_laibm(is)+leaf_bm(ip,is)/                   &
                                       (site%plots(ip)%species(ns)%leafarea_c*2.0)

               total_max_ht(is)  =max(total_max_ht(is),max_ht(ip,is))
               total_max_diam(is)=max(total_max_diam(is),max_diam(ip,is))

               total_dm_cats(is,:)=total_dm_cats(is,:)+diam_categories(ip,is,:)
            enddo
         endif

         call stddev(biomc(:,is),plot_mean_biomc(is),                            &
                                          plot_std_biomc(is),rnvalid)
         call stddev(biomn(:,is),plot_mean_biomn(is),                            &
                                          plot_std_biomn(is),rnvalid)

      enddo

      ! In the original code the standard deviation was not rescaled to the
      ! same scale as the plotlevel_biomc.  I am doing so here but that can be
      ! changed.  The original code didn't print the std_biomn at all.
      plotlevel_biomc=plotscale*biomc
      plot_mean_biomc=plotscale*plot_mean_biomc
      plot_std_biomc=plotscale*plot_std_biomc
      plotlevel_biomn=plotscale*biomn
      plot_mean_biomn=plotscale*plot_mean_biomn
      plot_std_biomn=plotscale*plot_std_biomn

      do is=1,species_pres%numspecies
         call csv_write(biom_by_s,iwmo,.false.)
         call csv_write(biom_by_s,year,.false.)
         call csv_write(biom_by_s,trim(adjustl(species_pres%spec_names(is,1))),  &
                                                                     .false.)
         call csv_write(biom_by_s,trim(adjustl(species_pres%spec_names(is,2))),  &
                                                                     .false.)
         if ( in_site(is) ) then

            total_basal(is)=total_basal(is)*plotrenorm
            total_laibm(is)=total_laibm(is)*plotrenorm
            total_biomc(is)=total_biomc(is)*plotadj
            total_biomn(is)=total_biomn(is)*plotrenorm*10.0
            total_dm_cats(is,:)=total_dm_cats(is,:)*plotadj

            call csv_write(biom_by_s,total_dm_cats(is,1),.false.)
            call csv_write(biom_by_s,total_dm_cats(is,2),.false.)
            call csv_write(biom_by_s,total_dm_cats(is,3),.false.)
            call csv_write(biom_by_s,total_dm_cats(is,4),.false.) 
            call csv_write(biom_by_s,total_dm_cats(is,5),.false.)
            call csv_write(biom_by_s,total_dm_cats(is,6),.false.)
            call csv_write(biom_by_s,total_dm_cats(is,7),.false.)
            call csv_write(biom_by_s,total_max_diam(is),.false.)
            call csv_write(biom_by_s,total_max_ht(is),.false.)
            call csv_write(biom_by_s,total_laibm(is),.false.)
            call csv_write(biom_by_s,total_basal(is),.false.)
            call csv_write(biom_by_s,total_biomc(is),.false.)
            call csv_write(biom_by_s,plot_std_biomc(is),.false.)
            call csv_write(biom_by_s,total_biomn(is),.false.)
            call csv_write(biom_by_s,plot_std_biomn(is),.true.)

         else

            call csv_write(biom_by_s,rnvalid,.false.)
            call csv_write(biom_by_s,rnvalid,.false.)
            call csv_write(biom_by_s,rnvalid,.false.)
            call csv_write(biom_by_s,rnvalid,.false.)
            call csv_write(biom_by_s,rnvalid,.false.)
            call csv_write(biom_by_s,rnvalid,.false.)
            call csv_write(biom_by_s,rnvalid,.false.)
            call csv_write(biom_by_s,rnvalid,.false.)
            call csv_write(biom_by_s,rnvalid,.false.)
            call csv_write(biom_by_s,rnvalid,.false.)
            call csv_write(biom_by_s,rnvalid,.false.)
            call csv_write(biom_by_s,rnvalid,.false.)
            call csv_write(biom_by_s,rnvalid,.false.)
            call csv_write(biom_by_s,rnvalid,.false.)
            call csv_write(biom_by_s,rnvalid,.true.)

         endif

      enddo

      if ( plot_level_data ) then

         do ip=1,site%numplots
            do is=1,species_pres%numspecies
               call csv_write(pl_biom_by_s,iwmo,.false.)
               call csv_write(pl_biom_by_s,year,.false.)
               call csv_write(pl_biom_by_s,ip,.false.)

               if ( in_site(is) ) then
                  call csv_write(pl_biom_by_s,diam_categories(ip,is,1),.false.)
                  call csv_write(pl_biom_by_s,diam_categories(ip,is,2),.false.)
                  call csv_write(pl_biom_by_s,diam_categories(ip,is,3),.false.)
                  call csv_write(pl_biom_by_s,diam_categories(ip,is,4),.false.)
                  call csv_write(pl_biom_by_s,diam_categories(ip,is,5),.false.)
                  call csv_write(pl_biom_by_s,diam_categories(ip,is,6),.false.)
                  call csv_write(pl_biom_by_s,diam_categories(ip,is,7),.false.)
                  call csv_write(pl_biom_by_s,max_diam(ip,is),.false.)
                  call csv_write(pl_biom_by_s,max_ht(ip,is),.false.)
                  call csv_write(pl_biom_by_s,leaf_bm(ip,is),.false.)
                  call csv_write(pl_biom_by_s,basal_area(ip,is),.false.)
                  call csv_write(pl_biom_by_s,biomC(ip,is),.false.)
                  call csv_write(pl_biom_by_s,biomC_std(ip,is),.false.)
                  call csv_write(pl_biom_by_s,biomN(ip,is),.false.)
                  call csv_write(pl_biom_by_s,biomN_std(ip,is),.true.)
               else
                  call csv_write(pl_biom_by_s,rnvalid,.false.)
                  call csv_write(pl_biom_by_s,rnvalid,.false.)
                  call csv_write(pl_biom_by_s,rnvalid,.false.)
                  call csv_write(pl_biom_by_s,rnvalid,.false.)
                  call csv_write(pl_biom_by_s,rnvalid,.false.)
                  call csv_write(pl_biom_by_s,rnvalid,.false.)
                  call csv_write(pl_biom_by_s,rnvalid,.false.)
                  call csv_write(pl_biom_by_s,rnvalid,.false.)
                  call csv_write(pl_biom_by_s,rnvalid,.false.)
                  call csv_write(pl_biom_by_s,rnvalid,.false.)
                  call csv_write(pl_biom_by_s,rnvalid,.false.)
                  call csv_write(pl_biom_by_s,rnvalid,.false.)
                  call csv_write(pl_biom_by_s,rnvalid,.false.)
                  call csv_write(pl_biom_by_s,rnvalid,.false.)
                  call csv_write(pl_biom_by_s,rnvalid,.true.)
               endif

            enddo
         enddo

      endif

   end subroutine write_species_data


   subroutine write_site_data(site,year)
      class(SiteData), intent(in)  :: site
      integer,         intent(in)  :: year
      integer                      :: iwmo

      iwmo=site%site_id

      call csv_write(clim_unit,iwmo,.false.)
      call csv_write(clim_unit,year,.false.)
      call write_site_csv(site,clim_unit)

   end subroutine write_site_data


   subroutine write_soil_data(site,year)
      class(SiteData), intent(in)  :: site
      integer,         intent(in)  :: year
      integer                      :: iwmo

      iwmo=site%site_id

      call csv_write(c_and_n,iwmo,.false.)
      call csv_write(c_and_n,year,.false.)
      call write_soil_csv(site%soil,c_and_n)

   end subroutine write_soil_data


   subroutine write_tree_data(site,year)
      class(SiteData), intent(in)  :: site
      integer,         intent(in)  :: year
      integer                      :: iwmo
      integer                      :: ip, it

      iwmo=site%site_id

      do ip=1,site%numplots
         do it=1,site%plots(ip)%numtrees
            call csv_write(tld,iwmo,.false.)
            call csv_write(tld,year,.false.)
            call csv_write(tld,ip,.false.)
            call csv_write(tld,it,.false.)
            call write_tree_csv(site%plots(ip)%trees(it),tld)
         enddo
      enddo

   end subroutine write_tree_data


end module Output
