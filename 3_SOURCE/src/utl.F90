
module utilities

private

public utl_casetrans,         &
       utl_int2char

public utl_duplicate_check,   &
       utl_wait_std

public utl_addquote,          &
       utl_delete_file,       &
       utl_nextunit
public utl_bomb_out,       &
       utl_writmess
       
public utl_reallocate_dblm, &
       utl_reallocate_dblv
       
public utl_get_command_line, &
       utl_parse_command_line

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                               STRING UTILITES                                    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine UTL_CASETRANS(string,hi_or_lo)

! -- Subroutine casetrans converts a string to upper or lower case.

! -- Arguments are as follows:-
!      string:    contains the string whose case must be changed
!      hi_or_lo:  must be either 'lo' or 'hi' to indicate
!                 change of case direction.

        character (len=*), intent(inout)        :: string
        character (len=*), intent(in)           :: hi_or_lo
        character                               :: alo, ahi
        integer                                 :: inc,i

        if(hi_or_lo.eq.'lo') then
          alo='A'; ahi='Z'; inc=iachar('a')-iachar('A')
        else if(hi_or_lo.eq.'hi') then
          alo='a'; ahi='z'; inc=iachar('A')-iachar('a')
        else
          return
        endif

        do i=1,len_trim(string)
          if((string(i:i).ge.alo).and.(string(i:i).le.ahi)) &
          string(i:i)=achar(iachar(string(i:i))+inc)
        enddo

        return

end subroutine UTL_CASETRANS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine UTL_INT2CHAR(value,string,nchar)

! -- Subroutine UTL_INT2CHAR converts an integer to a string.

	integer, intent(in)             :: value
	character (len=*), intent(out)  :: string
	integer, intent(in), optional   :: nchar
	character (len=12)              :: afmt
	integer                         :: llen

	string=' '
	afmt='(i    )'
	llen=min(30,len(string))
	if(present(nchar)) llen=min(llen,nchar)
	write(afmt(3:6),'(i4)') llen
	write(string(1:llen),afmt,err=100) value
	string=adjustl(string)
	if(string(1:1).eq.'*') go to 100
	return

100     string(1:llen)=repeat('#',llen)
	return

end subroutine utl_INT2CHAR



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                               GENERAL UTILITES                                   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! -- SUBROUTINE REALLOCATE_dblm expands an allocatable double precision matrix from
!    its current size to n rows by m columns.  Previous values are preserved
!    and new values are initialized to 9999999999.

subroutine utl_reallocate_dblm(A,n,m)
    double precision, pointer :: A(:,:), TMP(:,:)
    integer,       intent(in) :: m, n
    integer                   :: nold,mold
    nold=size(A,1)
    mold=size(A,2)
    allocate(TMP(nold,mold))
    TMP=A            
    deallocate(A)
    allocate(A(n,m))
    A=9999999999 ! matrix
    A(1:nold,1:mold) = TMP
    deallocate(TMP)
end subroutine utl_reallocate_dblm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! -- SUBROUTINE REALLOCATE_dblv expands an allocatable double precision vector
!    its current size to n elements.  Previous values are preserved
!    and new values are initialized to 9999999999.


subroutine utl_reallocate_dblv(A,n)
    double precision, pointer :: A(:), TMP(:)
    integer,       intent(in) :: n
    integer                   :: nold
    nold=size(A)
    allocate(TMP(nold))
    TMP=A            
    deallocate(A)
    allocate(A(n))
    A=9999999999 ! vector
    A(1:nold) = TMP
    deallocate(TMP)
end subroutine utl_reallocate_dblv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer function UTL_DUPLICATE_CHECK(ndim,aname,awork,aword)

! -- Function UTL_DUPLICATE_CHECK checks for duplicates in a list of names.
!    It is assumed that all names have already been converted to the correct
!    case, and to the same justification.

!    It is assumed that the word length in any of the character arrays is 20 characters
!    or less.

        implicit none

        integer, intent(in)              :: ndim
        character (len=*), intent(in)    :: aname(ndim)
        character (len=*), intent(inout) :: awork(ndim)
        character (len=*), intent(out)   :: aword

        character (len=20)   :: w
        integer              :: i,j,ifail

        ifail=0
        do  i=1, ndim
          awork(i) = aname(i)
        enddo

! -- Sort awork.

        do i = 2, ndim
           w = awork(i)
          do j = i,2,-1
             if ( w >= awork(j-1) ) exit
             awork(j) = awork(j-1)
          enddo
          awork(j) = w
        enddo

! -- Look for duplicates.

        do  i=2, ndim
          if(awork(i).eq.awork(i-1))then
            aword=adjustl(awork(i))
            ifail=1
            go to 9900
          endif
        enddo

9900  continue
        if(ifail.eq.0)then
          utl_duplicate_check=0
        else if(ifail.eq.1)then
          utl_duplicate_check=1
        endif
        return

end function UTL_DUPLICATE_CHECK


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine UTL_WAIT_STD(nsec)

! -- Subroutine UTL_WAIT occupies time for nsec hundredths of a second. It does this using
!    standard FORTRAN calls. This is a big consumer of resources. If the SLEEP command is
!    available, it is better to use that.

        implicit none

        integer ddate(8),iticks,iticks1,nsec

        call date_and_time(values=ddate)
        iticks=ddate(5)*360000+ddate(6)*6000+ddate(7)*100+ddate(8)/10
10      call date_and_time(values=ddate)
        iticks1=ddate(5)*360000+ddate(6)*6000+ddate(7)*100+ddate(8)/10
        if(iticks1.lt.iticks) iticks1=iticks1+8640000
        if(iticks1.lt.iticks+nsec) go to 10

        return

end subroutine UTL_WAIT_STD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                               ERROR UTILITES                                     !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine utl_bomb_out(errstruc)
    ! pulls an error message and pushes it out to standard error
    ! should create option to push to output file
    
    ! declarations    
    use error_message
    character*500            :: amessage
    integer                  :: n1, n2, i
    type (err_failure_struc) :: errstruc    

    write(6,*)
    n1=err_get_num_error(errstruc)
    do i = 1,n1
      call err_get_message(errstruc,i,amessage)
      amessage=' '//amessage
      call utl_writmess(6,amessage)
    enddo    
    
    write(6,*)
    do i = 1,n1
      call err_get_failfunc(errstruc,i,amessage)
      write(6,900) i, trim(amessage)    
    enddo    
900 format(' Message No. ',i2,':  Origin Function = ',a)
     stop 
end subroutine utl_bomb_out   



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine utl_writmess(iunit,amessage)

! -- Subroutine utl_writmess is used for writing error messages. It splits a text string
!    into lines of less than 80 characters in length, with the break made where a space
!    exists in the line.
! adapted from John Doherty's writmess code
    implicit none

	integer, intent(in)              :: iunit
    character (len=*), intent(inout) :: amessage
    integer                          :: jend,i,nblc,leadblank,itake,j
	character (len=20)               :: ablank

	ablank=' '
	itake=0
	j=0

    if(amessage.eq.' ')then
     write(iunit,*)
     return
    endif
     write(iunit,*)
	do i=1,min(20,len(amessage))
	 if(amessage(i:i).ne.' ')go to 21
20  enddo
21	leadblank=i-1
	nblc=len_trim(amessage)
5       jend=j+78-itake
	if(jend.ge.nblc) go to 100
	do i=jend,j+1,-1
	 if(amessage(i:i).eq.' ') then
	  if(itake.eq.0) then
	   write(iunit,'(a)') amessage(j+1:i)
	   itake=2+leadblank
	  else
	   write(iunit,'(a)') ablank(1:leadblank+2)//amessage(j+1:i)
	  endif
	  j=i
	  go to 5
	 endif
	enddo
	if(itake.eq.0)then
	 write(iunit,'(a)') amessage(j+1:jend)
     itake=2+leadblank
	else
     write(iunit,'(a)') ablank(1:leadblank+2)//amessage(j+1:jend)
	endif
	j=jend
	go to 5
100     jend=nblc
	if(itake.eq.0)then
	 write(iunit,'(a)') amessage(j+1:jend)
	else
	 write(iunit,'(a)') ablank(1:leadblank+2)//amessage(j+1:jend)
	endif
     return
   end subroutine utl_writmess

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                               FILE UTILITES                                      !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine UTL_ADDQUOTE(afile,aqfile)

! -- Subroutine UTL_ADDQUOTE adds quotes to a filename if it has a space in it.

        implicit none

        character (len=*), intent(in)   :: afile
        character (len=*), intent(out)  :: aqfile
        integer nbb

        if(index(trim(afile),' ').eq.0)then
          aqfile=afile
        else
          aqfile(1:1)='"'
          aqfile(2:)=trim(afile)
          nbb=len_trim(aqfile)+1
          aqfile(nbb:nbb)='"'
        endif

        return

end subroutine UTL_ADDQUOTE



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer function UTL_DELETE_FILE(afile)

! -- Function UTL_DELETE_FILE deletes a file.

        character (len=*), intent(in)      :: afile

        logical                            :: lexist
        integer                            :: iunit,ifail,ierr,jerr

        ifail=0

        inquire(file=afile,exist=lexist)
        if(lexist)then
          iunit=utl_nextunit()
          open(unit=iunit,file=afile,status='old',iostat=ierr)
          if(ierr.eq.0)then
            close(unit=iunit,status='delete',iostat=jerr)
            if(jerr.ne.0) then
              close(unit=iunit,iostat=jerr)
              ifail=1
              go to 9900
            endif
          else
            ifail=1
            go to 9900
          endif
        endif
        
9900    continue
        if(ifail.eq.1)then
          utl_delete_file=1
        else
          utl_delete_file=0
        endif

        return
        
end function UTL_DELETE_FILE

! Note that there seems to be a bug in lf95 in that if there is an error deleting
! the file there is a compiler error rather than ierr being given a nonzero value.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



integer function UTL_NEXTUNIT()

! -- Function UTL_NEXTUNIT determines the lowest unit number available for
! -- opening.

        logical::lopen

        do utl_nextunit=10,100
          inquire(unit=utl_nextunit,opened=lopen)
          if(.not.lopen) return
        enddo

        utl_nextunit=0

end function UTL_NEXTUNIT



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                            COMMAND LINE UTILITES                                 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE UTL_GET_COMMAND_LINE(COMMAND_LINE)

! -- Subroutine UTL_GET_COMMAND_LINE retreives any text supplied on the command line.
! -- This is adapted from the GET_COMMAND_LINE subroutine in the main PEST.F code from J. Doherty

      IMPLICIT NONE
      CHARACTER (LEN=*), INTENT(OUT)  :: COMMAND_LINE

      INTEGER             :: IARGC
      INTEGER             :: LLEN,NARG,IB,I,NB,IBB
      CHARACTER (LEN=100) :: ARG(4)


      COMMAND_LINE=' '

       LLEN=LEN(COMMAND_LINE)
       NARG=IARGC()
       IF(NARG.EQ.0) GO TO 100
       IB=0
       DO I=1,MIN(NARG,4)
         CALL GETARG(I,ARG(I))
         NB=LEN_TRIM(ARG(I))
         IBB=MIN(IB+NB+1,LLEN)
         COMMAND_LINE(IB+1:IBB)= ARG(I)(1:NB)
         IB=IBB
         IF(IB.GE.LLEN) GO TO 100
       enddo


100   CONTINUE
      RETURN

      END subroutine UTL_GET_COMMAND_LINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE UTL_PARSE_COMMAND_LINE(IFAIL,COMMAND_LINE,CTLFILE,RESTART)

! -- Subroutine UTL_PARSE_COMMAND_LINE parses the PEST command line.
! -- This is adapted from the PARSE_COMMAND_LINE subroutine in the main PEST.F code from J. Doherty
! -- RESTART is read in and returned, although it is not currently used (MNF 11/5/2010)
      IMPLICIT NONE

      INTEGER, INTENT(OUT)              :: IFAIL
      CHARACTER (LEN=*), INTENT(INOUT)  :: COMMAND_LINE
      CHARACTER (LEN=*), INTENT(OUT)    :: CTLFILE
      INTEGER, INTENT(OUT)              :: RESTART

      INTEGER                           :: K,IR,IS,IJ,ID,I,IO
      INTEGER                           :: IH,IM,IL,NP
      CHARACTER (LEN=1)                 :: AA

      IFAIL=0
      IF(COMMAND_LINE.EQ.' ') GO TO 9000
      AA=' '
      COMMAND_LINE=ADJUSTL(COMMAND_LINE)
      CTLFILE=COMMAND_LINE
      IF(CTLFILE(1:1).EQ.'"')THEN
        AA='"'
      ELSE IF(CTLFILE(1:1).EQ.'''')THEN
        AA=''''
      endif
      IF(AA.NE.' ') CTLFILE=CTLFILE(2:)
      I=INDEX(CTLFILE,AA)
      IF(I.LE.1) GO TO 9000
      CTLFILE=CTLFILE(1:I-1)
      RESTART=0
      IR=0
      IJ=0
      IS=0
      ID=0
      IO=0
      IH=0
      IM=0
      IR=INDEX(COMMAND_LINE,' /r ')
      IF(IR.EQ.0) IR=INDEX(COMMAND_LINE,' /R ')
      IJ=INDEX(COMMAND_LINE,' /j ')
      IF(IJ.EQ.0) IJ=INDEX(COMMAND_LINE,' /J ')
      IS=INDEX(COMMAND_LINE,' /s ')
      IF(IS.EQ.0) IS=INDEX(COMMAND_LINE,' /S ')
      ID=INDEX(COMMAND_LINE,' /d ')
      IF(ID.EQ.0) ID=INDEX(COMMAND_LINE,' /D ')
      IO=INDEX(COMMAND_LINE,' /i ')
      IF(IO.EQ.0) IO=INDEX(COMMAND_LINE,' /I ')
      IH=INDEX(COMMAND_LINE,' /h ')
      IF(IH.EQ.0) IH=INDEX(COMMAND_LINE,' /H ')
      IM=INDEX(COMMAND_LINE,' /m ')
      IF(IM.EQ.0) IM=INDEX(COMMAND_LINE,' /M ')
      IL=INDEX(COMMAND_LINE,' /l ')
      IF(IL.EQ.0) IL=INDEX(COMMAND_LINE,' /L ')

      IF(IR.NE.0)THEN
        RESTART=1
        COMMAND_LINE(IR+1:IR+2)='  '
        IR=1
      endif
      IF(IJ.NE.0)THEN
        RESTART=2
        COMMAND_LINE(IJ+1:IJ+2)='  '
        IJ=1
      endif
      IF(IS.NE.0)THEN
        RESTART=3
        COMMAND_LINE(IS+1:IS+2)='  '
        IS=1
      endif
      IF(ID.NE.0)THEN
        RESTART=4
        COMMAND_LINE(ID+1:ID+2)='  '
        ID=1
      endif
      IF(IO.NE.0)THEN
        RESTART=5
        COMMAND_LINE(IO+1:IO+2)='  '
        IO=1
      endif

      IF(IR+IJ+IS+ID+IO.GT.1) GO TO 9000
      IF(INDEX(COMMAND_LINE,' /').NE.0) GO TO 9000

! -- The following is used to handle spaces in filenames because of idiosyncracies
!    in the command line argument command.

      IF(AA.EQ.' ')THEN
        CTLFILE=COMMAND_LINE
      endif

      RETURN

9000  IFAIL=1

      RETURN
      END subroutine UTL_PARSE_COMMAND_LINE





end module utilities
