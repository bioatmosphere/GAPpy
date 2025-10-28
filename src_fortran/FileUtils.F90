module FileUtils
   use Constants
   use Parameters
   implicit none
   !*******************************************************************************
   !
   ! This module manages variables related to files and their paths.
   !
   !*******************************************************************************

   #ifdef PC
   integer, parameter      :: MAX_PATH=260
   character(len=1)        :: separator='\'
   #else
   integer, parameter      :: MAX_PATH=256
   character(len=1)        :: separator='/'
   #endif

   !Names of directories for input/output
   character(len=MAX_PATH) :: fqcwd
   character(len=MAX_PATH) :: inputdir,outputdir,climatedir,sitedir,slistdir

   integer, parameter      :: base_unit=10

   !Names of unit numbers so we can more easily tell which file is which

   ! Input
   integer                 :: rt_file=0,slist=0,sfile=0,cfile=0,cstdfile=0
   integer                 :: cgcmfile=0,splist=0,rglist=0,altlist=0

   ! Output
   integer                 :: logf
   integer                 :: c_and_n, clim_unit, tld
   integer                 :: biom_by_s, biom_by_g
   integer                 :: pl_biom_by_g, pl_biom_by_s
   integer                 :: pl_tree


contains


  subroutine get_cwd(fqcwd)
    character(len=*), intent(out) :: fqcwd
    ! May vary by compiler
    call getcwd(fqcwd)
  end subroutine get_cwd

  function unit_number()
     integer          :: unit_number
     integer          :: base_unit=10
     integer          :: iunit
     logical          :: first=.true.
     save                 

     if ( first ) then
        iunit=base_unit
        first=.false.
     else
        iunit=iunit+1
     endif
     unit_number=iunit

  end function unit_number

  function open_file(filename,mode)
     ! Opens the file filename if it can, returns a unit number for it.
     integer                                 :: open_file
     character(len=*), intent(in)            :: filename
     character(len=*), intent(in), optional  :: mode
     character(len=9)              :: fmode
     logical                       :: file_exists
     character(len=MAX_PATH)       :: fname
     character(len=MAX_CHAR)       :: message
     integer                       :: i, ios, iunit
     integer, dimension(MAX_PATH) :: farray

     if ( present(mode)) then
        if ( mode == 'r' .or. mode == 'R' ) then
           fmode='read'
        else if ( mode == 'w' .or. mode == 'W' ) then
           fmode='write'
        else if ( mode == 'rw' .or. mode == 'RW'                               &
             .or. mode == 'wr' .or. mode == 'WR' ) then
           fmode='readwrite'
        else
           fmode='readwrite'
        endif
     else
        fmode='readwrite'
     endif

     fname=trim(adjustl(filename))

     if ( fmode == 'read' .or. fmode == 'readwrite') then
     ! File needs to exist
     !This is the best I could figure out to handle uninitialized filename
     !errors (which probably contain garbage characters so they aren't blank)
        farray=0
        do i=1,len_trim(fname)
           farray(i)=ichar(fname(i:i))
        enddo
        if (any(farray > 127)) then
           call fatal_error("Invalid filename")
        endif
     endif
     
     inquire(file=fname,exist=file_exists)

     if ( file_exists .and. fmode == 'write' ) then
        write(message,'(a,a,a)') "File ",fname(1:len_trim(fname)),            &
                                  " exists. Cannot open write only."
        call warning(message)
        open_file=invalid
     else if ( .not. file_exists .and. fmode == 'read' ) then
        iunit=invalid
     else
        iunit=unit_number()
        open(iunit,file=fname,action=fmode,iostat=ios)
        if ( ios .ne. 0 ) then
           iunit=invalid
        endif
     endif

     open_file=iunit

  end function open_file

  
  function count_records(funit,nheaders)
    ! Counts the number of lines, not including the first nheaders lines.
    integer                       :: count_records
    integer, intent(in)           :: funit
    integer, intent(in)           :: nheaders

    integer                       :: ios, lpcount
    character(len=MAX_LINE)       :: line

    lpcount = 0

    do
      read(funit,*,end=10,iostat=ios) line
         if (ios==0) then
         ! success, increment counter
            lpcount = lpcount + 1
         else
            call fatal_error("Unable to read file")
         endif
    end do
   10  continue

    count_records = lpcount-nheaders

    rewind(funit)

  end function count_records


  subroutine build_pathname(subdir,filename,pathname)
    ! Returns the full path to a file
    character(len=*), intent(in)  ::  subdir
    character(len=*), intent(in)  ::  filename
    character(len=*),intent(out)  ::  pathname

    character(len=MAX_PATH)       ::  pname
    character(len=MAX_FILE)       ::  tmpsub,tmpfile
    character(len=15)             ::  fmtstr
    integer                       ::  flen, slen

    ! Reset pathname
    pathname=''

    tmpsub=subdir
    tmpfile=filename
    call ws_strip(tmpsub)
    call ws_strip(tmpfile)
    slen=len_trim(adjustl(tmpsub))
    flen=len_trim(adjustl(tmpfile))
    write(fmtstr,'(a,i2,a,a,i2,a)'),'(a',slen,',a1,','a',flen,')'

    write(pname,fmtstr) tmpsub(1:slen),separator,tmpfile(1:flen)
    call ws_strip(pname)

    pathname=pname

  end subroutine build_pathname


  subroutine build_filename(string,suffix,filename)
    ! Returns a filename from a specified string and a suffix
    character(len=*), intent(in) ::  suffix
    character(len=*), intent(in) ::  string
    character(len=*), intent(out)::  filename

    character(len=MAX_FILE/2)     ::  tmpstr,tmpsuf
    character(len=MAX_FILE)       ::  fname
    character(len=15)             ::  fmtstr
    integer                       ::  strlen,suflen

    ! Reset filename
    filename=''

    tmpstr=string
    tmpsuf=suffix
    call ws_strip(tmpstr)
    call ws_strip(tmpsuf)
    strlen=len_trim(adjustl(tmpstr))
    suflen=len_trim(adjustl(tmpsuf))
    ! Construct name of file
    write(fmtstr,'(a,i2,a,i2,a)'),'(a',strlen,',a',suflen,')'
    write(fname,fmtstr),tmpstr(1:strlen),tmpsuf(1:suflen)
    call ws_strip(fname)

    filename=fname

  end subroutine build_filename


  subroutine split_line(line,fields)
    character(len=*)                            :: line
    character(len=*),dimension(:), allocatable  ::  fields

    character(len=len_trim(line))               :: tempstrg,rem_string
    character(len=MAX_LONG_LINE),dimension(:),allocatable :: fn
    character(len=1)                            :: delimiter,ch

    integer                                     :: i, n

    integer                                     :: nfiles

    ! Parses the single character with the list of species into the
    ! separate names, and returns them as an array. Sets the parton
    ! kinds to zero.

    delimiter=','

    call ws_strip(line)
    tempstrg=adjustl(line)

    allocate(fn(MAX_FIELDS))

    n=1
    do
       call split(delimiter,tempstrg,fn(n),rem_string)
       tempstrg=adjustl(rem_string)
       if ( tempstrg .eq. '' ) then
          exit
       else
          n=n+1
       endif

    enddo

    allocate(fields(n))

    do i=1,n
       fields(i)=trim(adjustl(fn(i)))
    enddo

  end subroutine split_line


  subroutine split(delimiter,string,res_string,rem_string)
    character(len=1)                             ::  delimiter
    character(len=*)                             ::  string
    character(len=*)                             ::  res_string, rem_string

    integer                                      ::  i, ii, loc

    ! Splits a string and returns the part of the string before a delimiter, along
    ! with the remainder of the string.

    res_string=''
    rem_string=''

    string=adjustl(string)

    loc=scan(string,delimiter)

    ! Return if we don't find the delimiter in the string
    if (loc .eq. 0) then
       res_string=adjustl(string)
       rem_string=''
       return
    endif

    do i=1,loc-1
       res_string(i:i)=string(i:i)
    enddo

    do i=loc+1,len_trim(string)
       ii=i-loc
       rem_string(ii:ii)=string(i:i)
    enddo

  end subroutine split


  subroutine ws_strip(string)
    character(len=*)                             ::  string

    character(len=len(string))                   ::  tmpstring
    character(len=1)                             ::  c
    integer                                      ::  i, n, len_string

    ! Removes whitespace (spaces and tabs) from a string.

    string=adjustl(string)
    len_string=len_trim(string)

    tmpstring=''

    n=0
    do i=1, len_string
       c = string(i:i)
       if ( (c .eq. ' ') .or. (c .eq. '\t') ) then
          cycle
       else
          n=n+1
          tmpstring(n:n)=c
       endif
    enddo

    string=trim(adjustl(tmpstring))

  end subroutine ws_strip


  subroutine quote_strip(string)
    character(len=*)                             ::  string

    character(len=len(string))                   ::  tmpstring
    character(len=1)                             ::  c
    integer                                      ::  i, n, len_string

    ! Removes quotes (single or double) from a string

    string=adjustl(string)
    len_string=len_trim(string)

    tmpstring=''

    n=0
    do i=1, len_string
       c = string(i:i)
       if ( (c .eq. "'") .or. (c .eq. '"') ) then
          cycle
       else
          n=n+1
          tmpstring(n:n)=c
       endif
    enddo

    string=trim(adjustl(tmpstring))

  end subroutine quote_strip

  subroutine warning(message)
     character(len=*), intent(in)  :: message
     write(*,*) message
  end subroutine warning

  subroutine fatal_error(message)
     character(len=*), intent(in)  :: message
     write(*,*) message
     stop
  end subroutine fatal_error

end module FileUtils