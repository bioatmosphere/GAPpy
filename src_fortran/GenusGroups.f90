module Genusgroups
  
  use Constants
  use Site
  use vararray

  implicit none


  type Groups
    character(len=MAX_NLEN), allocatable, dimension(:)   :: genusgroups
    character(len=MAX_NLEN), allocatable, dimension(:,:) :: spec_names
    integer                                              :: numgenera
    integer                                              :: numspecies
  end type Groups

contains

  subroutine initialize_genus_groups(group,sites)

    type(Groups),                    intent(out)   :: group
    type(SiteData),    dimension(:), intent(inout) :: sites


    type(SpeciesData), dimension(:), allocatable   :: species_data

    character(len=MAX_NLEN), dimension(:), allocatable::species_names
    integer                                        :: numsites
    integer                                        :: numspecies
    integer                                        :: num_all_species
    integer                                        :: n, ns, na, nv
    logical                                        :: already_seen

    numsites=size(sites)

    num_all_species=0

    do n=1,numsites
       num_all_species=num_all_species+size(sites(n)%species)
    enddo

    allocate(species_data(num_all_species))
    nv=1
    do ns=1,numsites
       do na=1,size(sites(ns)%species)
          species_data(nv)%unique_id=sites(ns)%species(na)%unique_id
          species_data(nv)%genus_name=sites(ns)%species(na)%genus_name
          nv=nv+1
       enddo
    enddo
    call get_unique_items(species_data(:)%unique_id,species_names)
    call get_unique_items(species_data(:)%genus_name,group%genusgroups)

    group%numspecies=size(species_names)
    group%numgenera =size(group%genusgroups)

    allocate(group%spec_names(group%numspecies,2))

    do na=1,size(species_data)
       do ns=1,group%numspecies
          if ( species_names(ns) .eq. species_data(na)%unique_id ) then
             group%spec_names(ns,1)=species_data(na)%genus_name
             group%spec_names(ns,2)=species_data(na)%unique_id
             cycle
          endif
       enddo
    enddo

  end subroutine initialize_genus_groups


  subroutine get_unique_items(array,unique_array)
     character(len=*),        dimension(:),               intent(inout):: array
     character(len=MAX_NLEN), dimension(:), allocatable,  intent(out)  ::      &
                                                                   unique_array 

     character(len=MAX_NLEN), dimension(size(array))   :: temp_array
     integer                                           :: n, ncount, nsize

     nsize=size(array)
     temp_array=array

     ! If the array has only one element then we are done
     if ( nsize .eq. 1 ) then
        allocate(unique_array(1))
        unique_array(1)=array(1)
        return
     endif

     call sort(temp_array,nsize)

     ! First count the number of unique items (only required since I want the
     ! array returned to be the correct length)
     ncount=1
     do n=2, nsize
        if ( temp_array(n-1) .ne. temp_array(n) ) then
           ncount=ncount+1
        endif
     enddo

     allocate(unique_array(ncount))

     unique_array(1)=temp_array(1)
     ncount=2
     do n=2,nsize
        if ( temp_array(n-1) .ne. temp_array(n) ) then
           unique_array(ncount)=trim(temp_array(n))
           ncount=ncount+1
        endif
     enddo

  end subroutine get_unique_items


end module GenusGroups
