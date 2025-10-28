module Utilities
  
  use Constants
  
  implicit none

contains

  function kron(x)
    ! Kronecker delta
    real,      intent(in)             :: x
    real                              :: kron

    if (x .gt. 0.0) then
       kron=1.0
    else
       kron=0.0
    endif

  end function kron


  function roundtoN(x,N)
    ! Round to N decimal places
    real,    intent(in)             :: x
    integer, intent(in)             :: N
    real                            :: roundtoN

    real                            :: factor

    factor=10**N
    roundtoN=float(int(x * factor  + 0.5)) / factor

  end function roundtoN


  subroutine stddev(array,arith_mean,std_deviation,missing)
    ! Compute the mean and standard deviation of an array
    real, dimension(:), intent(in)   :: array
    real,               intent(out)  :: arith_mean
    real,               intent(out)  :: std_deviation
    real, optional,     intent(in)   :: missing

    real                             :: sqr_sum
    integer                          :: nelems

       if (present(missing)) then
          nelems=count(array .ne. missing)
          arith_mean=sum(array,mask=array .ne. missing)/float(nelems)
          sqr_sum=sum(array**2,mask=array .ne. missing)
       else
          nelems=size(array)
          arith_mean=sum(array)/float(nelems)
          sqr_sum=sum(array**2)
       endif

       if ( nelems == 0 ) then
          std_deviation=0.0
          return
       endif

       !Sample standard deviation
       std_deviation=sqrt((sqr_sum-nelems*arith_mean**2)/(nelems-1))

  end subroutine stddev


  subroutine sort(array,length)
    ! Sort (because we'd like to have the genera etc. names in alphabetical
    ! order)
    character(len=*), dimension(length), intent(inout) :: array
    integer,                             intent(in)    :: length

    character(len=MAX_NLEN)                    :: temp
    integer                                    :: inc,lower,upper
    integer                                    :: i,j

    ! Sorts an array with indices running from 1 to length in place using
    ! Shell sort
    ! Shell sort is fine for small arrays (these will never be large).
    ! This is for character strings up to length NLEN.
    ! It will use lexigraphical ordering.

    lower=lbound(array,1)
    upper=ubound(array,1)

    inc=length/2

    ! We continue until the stride is the lower bound

    do while (inc > 0)
       do i=inc+lower,upper
          temp=adjustl(array(i))
          j=i
          ! Cannot use j>inc .and. array(j-inc)>temp if it doesn't short-circuit
          do while ( j > inc )
             if ( array(j-inc) > temp ) then
                array(j)=array(j-inc)
                j=j-inc
             else
                exit
             endif
          enddo
          array(j)=temp
       enddo
       if ( inc == 2 ) then
          inc=1
       else
          inc=int(real(inc)/2.2)
       endif
    enddo

  end subroutine sort

end module Utilities
