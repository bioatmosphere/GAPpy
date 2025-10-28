module lists
! The module used must overload the assigment operator
use Species, list_type=>SpeciesData
implicit none

! So far I have only implemented append and delete, but that's all I need for
! most programs.

contains

   subroutine append(list,new_item)
   type(list_type), dimension(:), allocatable, intent(inout) :: list
   type(list_type),                            intent(inout) :: new_item
   type(list_type), dimension(:), allocatable                :: temp_list
   integer                                                   :: n, list_len

      if ( .not. allocated(list) ) then
          allocate(list(1))
          list(1)=new_item
      else
         list_len=size(list)
         allocate(temp_list(list_len))

         do n=1,list_len
            temp_list(n)=list(n)
         enddo

         deallocate(list)
         allocate(list(list_len+1))

         do n=1,list_len
            list(n)=temp_list(n)
         enddo
         list(list_len+1)=new_item

      endif

   end subroutine append


   subroutine delete(list,list_index)
   type(list_type), dimension(:), allocatable, intent(inout) :: list
   integer,                                    intent(in)    :: list_index
   type(list_type), dimension(:), allocatable                :: temp_list
   integer                                                   :: list_len
   integer                                                   :: n, ncount

      if ( .not. allocated(list) .or. size(list)<1 ) then
           write(*,*) "Empty list, cannot remove item"
          return
      endif

      list_len=size(list)

      if ( list_index > list_len ) then
           write(*,*) "List index not in range of list"
           return
      endif

      allocate(temp_list(list_len))

      do n=1,list_len
         temp_list(n)=list(n)
      enddo

      deallocate(list)
      allocate(list(list_len-1))

      ncount=1
      do n=1,list_len
         if ( n .ne. list_index ) then
            list(ncount)=temp_list(n)
            ncount=ncount+1
         endif
      enddo

   end subroutine delete

end module lists 
