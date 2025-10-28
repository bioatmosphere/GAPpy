module DATA_MODULE
  use Constants

  type CHAR_DATA
     character(len=MAX_NLEN)  ::  value
  end type CHAR_DATA

  type(CHAR_DATA), parameter :: empty_data=CHAR_DATA('')

end module DATA_MODULE


module vararray
    use DATA_MODULE, LIST_DATA=>CHAR_DATA,                                       &
                empty_list_data=>empty_data

    ! lists.f90
    !     Implementation of dynamically growing arrays
    !
    !     Basically this is Arjen Markus's vector.f90 module but I have
    !     renamed it from vector to lists, because most scientists don't
    !     care about fine distinctions between "linked lists" and C++-type
    !     "vectors" and this object behaves like a Python/Perl/etc. list. 
    !
    !     The module is straightforward: it defines a suitable
    !     data structure, data can be added to the list
    !     and you can retrieve data from it.
    !
    !     Note:
    !     For the function list_at() we need a parameter
    !     that represents the "empty list data" value.
    !
    !     $Id: vectors.f90,v 1.4 2008/06/12 15:12:39 relaxmike Exp $
    !
    ! Written by Arjen Markus

    type LIST
        private
        integer                                  :: not_used
        type(LIST_DATA), dimension(:), pointer :: data    => null()
    end type list

    private
    public :: LIST
    public :: LIST_DATA
    public :: list_create
    public :: list_append
    public :: list_at
    public :: list_size
    public :: list_put
    public :: list_delete_elements
    public :: list_destroy
    public :: list_insert_empty

    real, parameter :: growth_rate = 1.1

contains

    ! list_create
    !     Create a new list
    !
    ! Arguments:
    !     lst           Variable that should hold the list
    !     capacity      Initial capacity (optional)
    !
    ! Note:
    !     The fields of the list data structure are set
    !
    subroutine list_create( lst, capacity )
        type(LIST)         :: lst
        integer, optional  :: capacity

        integer            :: cap

        !
        ! Check that the list does not have any data left
        !
        if ( associated(lst%data) ) then
            call list_destroy( lst )
        endif

        if ( present(capacity) ) then
            cap = max( 1, capacity )
        else
            cap = 10
        endif
        allocate( lst%data(1:cap) )
        lst%not_used = 0
    end subroutine list_create

    ! list_destroy
    !     Destroy a list
    !
    ! Arguments:
    !     lst           list in question
    !
    subroutine list_destroy( lst )
        type(LIST)       :: lst

        !
        ! Check that the list does not have any data left
        !
        if ( associated(lst%data) ) then
            deallocate( lst%data )
        endif
        lst%not_used = 0
    end subroutine list_destroy

    ! list_size
    !     Return the number of elements in use
    !
    ! Arguments:
    !     lst           Variable that holds the list
    !
    integer function list_size( lst )
        type(LIST)       :: lst

        list_size = lst%not_used
    end function list_size

    ! list_at
    !     Get the value of the nth element of the list
    !
    ! Arguments:
    !     lst           list in question
    !     n             Index of the element whose value
    !                   should be retrieved
    !
    type(LIST_DATA) function list_at( lst, n )
        type(LIST)       :: lst
        integer          :: n

        if ( n .lt. 1 .or. n .gt. lst%not_used ) then
            list_at = empty_list_data
        else
            list_at = lst%data(n)
        endif
    end function list_at

    ! list_insert_empty
    !     Insert one or more empty elements
    !
    ! Arguments:
    !     list        list in question
    !     pos           Position to insert the empty elements
    !     number        Number of empty elements
    !
    subroutine list_insert_empty( lst, pos, number )
        type(LIST)         :: lst
        integer, intent(in)  :: pos
        integer, intent(in)  :: number

        integer              :: i

        if ( number .lt. 1 .or. pos .lt. 1 .or. pos .gt. lst%not_used ) then
            return
        endif

        if ( lst%not_used+number .ge. size(lst%data) ) then
            call list_increase_capacity( lst, lst%not_used+number )
        endif

        do i = lst%not_used,pos,-1
            lst%data(i+number) = lst%data(i)
        enddo

        do i = 1,number
            lst%data(pos+i-1) = empty_list_data
        enddo

        lst%not_used = lst%not_used + number
    end subroutine list_insert_empty

    ! list_delete_elements
    !     Delete one or more elements
    !
    ! Arguments:
    !     list        list in question
    !     pos           Position to start deletion
    !     number        Number of elements
    !
    subroutine list_delete_elements( lst, pos, number )

        type(LIST)         :: lst
        integer, intent(in)  :: pos
        integer, intent(in)  :: number

        integer              :: i

        if ( number .lt. 1 .or. pos .lt. 1 .or. pos .gt. lst%not_used ) then
            return
        endif

        do i = pos,lst%not_used-number
            lst%data(i) = lst%data(i+number)
        enddo

        lst%not_used = lst%not_used - number

    end subroutine list_delete_elements

    ! list_append
    !     Append a value to the list
    !
    ! Arguments:
    !     lst           list in question
    !     data          Data to be appended
    !
    subroutine list_append( lst, data )

        type(LIST)       :: lst
        type(LIST_DATA)  :: data

        if ( lst%not_used .ge. size(lst%data) ) then
            call list_increase_capacity( lst, lst%not_used+1 )
        endif

        lst%not_used = lst%not_used + 1
        lst%data(lst%not_used) = data

    end subroutine list_append

    ! list_put
    !     Put a value at a specific element of the list
    !     (it needs not yet exist)
    !
    ! Arguments:
    !     lst           list in question
    !     n             Index of the element
    !     data          Data to be put in the list
    !
    subroutine list_put( lst, n, data )

        type(LIST)       :: lst
        integer            :: n
        type(LIST_DATA)  :: data

        if ( n .lt. 1 ) then
            return
        endif
        if ( n .gt. size(lst%data) ) then
            call list_increase_capacity( lst, n )
        endif

        lst%not_used = max( lst%not_used, n)
        lst%data(n) = data

    end subroutine list_put

    ! list_increase_capacity
    !     Expand the array holding the data
    !
    ! Arguments:
    !     lst           list in question
    !     capacity      Minimum capacity
    !
    subroutine list_increase_capacity( lst, capacity )

        type(LIST)       :: lst
        integer            :: capacity

        integer            :: new_cap
        type(LIST_DATA), dimension(:), pointer :: new_data

        new_cap = max( capacity, nint( growth_rate * size(lst%data) ) )

        if ( new_cap .gt. size(lst%data) ) then
            allocate( new_data(1:new_cap) )
            new_data(1:lst%not_used) = lst%data(1:lst%not_used)
            new_data(lst%not_used+1:new_cap) = empty_list_data

            deallocate( lst%data )
            lst%data => new_data
        endif

    end subroutine list_increase_capacity

end module vararray
