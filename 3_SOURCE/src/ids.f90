module jupiter_input_data_support

       private

! -- Size variables
       integer, parameter   :: MAXCOL=20           ! Maximum number of columns in TABLE block
       integer, parameter   :: COLWIDTH=50         ! Maximum characters in column header
       integer, parameter   :: KEYWORDWIDTH=100    ! Maximum characters in a keyword
       integer, parameter   :: KEYSTRINGWIDTH=256  ! Maximum characters in a keyword's value
       integer, parameter   :: FILEWIDTH=256       ! Maximum characters in a filename
       integer, parameter   :: LABELWIDTH=100      ! Maximum characters in a blockname
       integer, parameter   :: LINEWIDTH=2000      ! Maximum length of line in input file
       integer, parameter   :: ERRORWIDTH=2000     ! Maximum length of error message string
       
       public  COLWIDTH, ERRORWIDTH, MAXCOL, &
                KEYWORDWIDTH, KEYSTRINGWIDTH, FILEWIDTH, LABELWIDTH, LINEWIDTH

! -- Character variables

       character (len=1)    :: comma=','
       character (len=1)    :: squote=''''
       character (len=1)    :: dquote='"'
       character (len=1)    :: space=' '
       character (len=1)    :: tab=char(9)
       character (len=1)    :: equality='='

! -- General variables

       integer                    :: block_status_call=0     ! Number of calls to IDS_BLOCK_STATUS
       integer                    :: ifail=0                 ! Indicates error condition
       integer                    :: linepos                 ! Cursor position in processing of line
       integer                    :: newblockflag=0          ! Indicates that a block has just been opened
       integer                    :: endblockflag=0          ! Set to 1 if "end of block" line has already been read
       integer                    :: position=0              ! 0=between blocks; 1=in keyword block; 2=in table block
       integer                    :: tablerow=0              ! The current row when reading a table
       character (len=LABELWIDTH) :: blockname_in_last=' '   ! Blockname used on previous IDS_READ_BLOCK_TABLE call

! -- Block variables

       integer                         :: emptyfileblockflag=0   ! Flag for a FILE block that contains no files
       character (len=LABELWIDTH)      :: blocklabel=' '         ! Name of current block
       character (len=LABELWIDTH)      :: blocklabelmain=' '     ! Name of current block in main input file
       integer                         :: blockformat            ! 1=keyword; 2=table; 3=files;
       integer                         :: blockformatmain        ! 1=keyword; 2=table; 3=files;
       character (len=KEYSTRINGWIDTH)  :: nextword=' '           ! The next word (with quotes removed)
       character (len=KEYWORDWIDTH  )  :: keyword=' '            ! Current keyword
       character (len=KEYSTRINGWIDTH)  :: keystring=' '          ! String associated with keyword

! ----- File variables

       integer                    :: fileflag=1     ! 1 when reading main file; 2 when reading subsidiary file
       integer                    :: iunit          ! The unit number of main input file
       integer                    :: junit          ! The unit number of subsiduary input file
       integer                    :: kunit          ! General unit number
       integer                    :: iline=0        ! Line number in main input file
       integer                    :: jline=0        ! Line number in subsidiary input file
       integer                    :: kline          ! General line number
       character (len=10)         :: aline          ! Used for writing line numbers in error messages
       character (len=FILEWIDTH)  :: infile         ! The name of the main input file
       character (len=FILEWIDTH)  :: qinfile        ! The name of the main input file with quotes if space
       character (len=FILEWIDTH)  :: subfile        ! The name of the current subsidiary input file
       character (len=FILEWIDTH)  :: qsubfile       ! The name of the current subsidiary input file with quotes if space
       character (len=FILEWIDTH)  :: qfile          ! Used in writing error messages

! ----- Table variables

       integer                   :: ncol                     ! Number of columns in a table
       integer                   :: nrow                     ! Number of rows in a table
       integer                   :: colexport(MAXCOL)        ! Index of table columns to export
       integer                   :: maxcolexport             ! The maximum column number for exporting
       character (len=COLWIDTH)  :: columnlabel(MAXCOL)      ! Column labels
       integer                   :: endoftableblockflag = 0  ! The end of a table block has been encountered

       character (len=LINEWIDTH)       :: cline              ! Character string used for reading lines of a file
       character (len=ERRORWIDTH)      :: amessage           ! String used for error messages

! -- Visible subroutines

       public ids_block_status,              &
              ids_read_block_keywords,       &
              ids_read_block_table,          &
              ids_get_message_string

contains


       subroutine ids_block_status(ifail_out,iunit_in,junit_in,numblock,label,numrow,infile_in)

! -- Subroutine IDS_BLOCK_STATUS returns the presence/absence of a block. If the block is a TABLE
!    block it provides the number of rows in the block.

         implicit none

         integer, intent(out)           :: ifail_out
         integer, intent(in)            :: iunit_in
         integer, intent(in)            :: junit_in
         integer, intent(in)            :: numblock
         character (len=*), intent(in)  :: label(:)
         integer, intent(out)           :: numrow(:)
         character (len=*), intent(in)  :: infile_in

         logical                        :: lopened
         integer                        :: iblock,jblock
         character (len=LABELWIDTH)     :: atemp,atemp1

! -- Initialisation

         ifail_out=0
         iunit=iunit_in
         junit=junit_in
         numrow=0                 ! an array
         infile=infile_in
         call addquote(infile,qinfile)
         call initialise()

! -- The status of IUNIT is checked.

         inquire(unit=iunit,opened=lopened)
         if(.not.lopened)then
           write(amessage,3)
3          format('Programming error in call to subroutine IDS_BLOCK_STATUS: IUNIT is ', &
           'not connected to an opened file.')
           go to 9890
         endif
         inquire(unit=junit,opened=lopened)
         if(lopened)then
           write(amessage,4)
4          format('Programming error in call to subroutine IDS_BLOCK_STATUS: JUNIT is ', &
           'already connected to an opened file.')
           go to 9890
         endif

! -- A check is made that the user has not supplied duplicate block names.

         if(numblock.gt.1)then
           do iblock=1,numblock-1
             atemp=adjustl(label(iblock))
             call uppercase(atemp)
             do jblock=iblock+1,numblock
               atemp1=adjustl(label(jblock))
               call uppercase(atemp1)
               if(atemp.eq.atemp1)then
                 write(amessage,5)
5                format('Programming error - duplicate names supplied in LABEL array ',  &
                 'in call to subroutine IDS_BLOCK_STATUS.')
                 go to 9890
               endif
             enddo
           enddo
         endif

! -- We now find the requested blocks.

         do
           call find_next_block()
           if(ifail.ne.0) go to 9890
           if(linepos.eq.-999) go to 100
           if(blockformat.eq.2)then
             call read_table_info()
             if(ifail.ne.0) go to 9890
           else if((blockformat.eq.3).and.(emptyfileblockflag.eq.1))then
             nrow=0
             ncol=0
           endif

           do iblock=1,numblock
             atemp=adjustl(label(iblock))
             call uppercase(atemp)
             if(atemp.eq.blocklabel)then
               if(iblock.lt.numblock)then
                 do jblock=iblock+1,numblock
                   if(numrow(jblock).ne.0)then
                     numrow(iblock)=-999
                     go to 50
                   endif
                 enddo
               endif
               if(blockformat.eq.1)then
                 if(fileflag.eq.2)then
                   write(amessage,20) trim(blocklabel),trim(qsubfile)
20                 format('Block ',a,' in file ',a,' has a blockformat of KEYWORD and is thus ', &
                   'not permitted in a subsidiary file cited in a FILE format block in main input file.')
                   go to 9890
                 endif
                 if(numrow(iblock).ne.0)then
                   if(numrow(iblock).gt.0)then
                     write(amessage,30) trim(blocklabel),trim(qinfile)
30                   format('The name ',a,' has been ascribed to both a KEYWORD and TABLE ', &
                     'block in file ',a,'.')
                     go to 9890
                   else
                     numrow(iblock)=-999
                   endif
                 else
                   numrow(iblock)=-1
                 endif
               else
                 if(numrow(iblock).eq.-1)then
                   write(amessage,40) trim(blocklabel)
40                 format('The name ',a,' has been ascribed to both a KEYWORD and TABLE ', &
                   'block in input dataset.')
                   go to 9890
                 else
                   numrow(iblock)=numrow(iblock)+nrow
                 endif
               endif
               go to 50
             endif
           enddo
50         continue
           call find_end_of_block(.TRUE.)
           if(ifail.ne.0) go to 9890
         enddo

100      continue

         rewind(unit=iunit)
         call initialise()
         block_status_call=block_status_call+1

         return


9890     ifail_out=1
         return

       end subroutine ids_block_status





       subroutine ids_read_block_keywords(ifail_out,blockname_in,numlist,  &
       keywordname,keywordstring,keywordline)

! -- Subroutine IDS_READ_BLOCK_KEYWORDS reads keywords in a KEYWORDS block, checking them against those
!    provided in a user-supplied list.

         implicit none

         integer, intent(out)           :: ifail_out
         character (len=*), intent(in)  :: blockname_in
         integer, intent(in)            :: numlist
         character (len=*), intent(in)  :: keywordname(numlist)
         character (len=*), intent(out) :: keywordstring(numlist)
         integer, intent(out)           :: keywordline(numlist)

         integer                        :: i,ilist,jlist
         character (len=LABELWIDTH)     :: ablock
         character (len=KEYWORDWIDTH)   :: akeyword,atemp


! -- Security

         if(block_status_call.eq.0)then
           write(amessage,2)
2          format('Programming error in call to subroutine IDS_READ_BLOCK_KEYWORDS: subroutine ',  &
           'IDS_BLOCK_STATUS has not yet been called.')
           go to 9890
         endif

! -- Initialisation

         ifail_out=0
         keywordline=0         ! An array
         keywordstring=' '     ! An array

! -- We check that we have indeed finished dealing with the previous block.

         if(position.ne.0)then
           write(amessage,1) trim(blockname_in)
1          format('Programming error in call to subroutine IDS_READ_BLOCK_KEYWORDS: an attempt ',    &
           'was made to process "',a,'" block when processing of previous block was not complete. '  &
           'If the previously processed block was of TABLE format, it is possible that subroutine ', &
           'IDS_READ_BLOCK_TABLE was called fewer times than there are rows in the (collective) ',   &
           'table(s) for this block.')
           go to 9890
         endif

! -- We check that there are no duplicates in the keyword list.

         if(numlist.gt.1)then
           do ilist=1,numlist-1
             akeyword=adjustl(keywordname(ilist))
             call uppercase(akeyword)
             do jlist=ilist+1,numlist
               atemp=adjustl(keywordname(jlist))
               call uppercase(atemp)
               if(atemp.eq.akeyword)then
                 write(amessage,5)
5                format('Programming error in call to subroutine IDS_READ_BLOCK_KEYWORDS: ',  &
                 'duplicate keywords provided in KEYWORDNAME list.')
                 go to 9890
               endif
             enddo
           enddo
         endif

! -- Now we find the requested block.

         ablock=adjustl(blockname_in)
         call uppercase(ablock)
         do
           call find_next_block()
           if(ifail.ne.0) go to 9890
           if(blocklabel.eq.' ')then
             write(amessage,10) trim(ablock),trim(qinfile)
10           format('Error in call to subroutine IDS_READ_BLOCK_KEYWORDS: ',  &
             'cannot find ',a,' block through forward search in file ',a,'.')
             go to 9890
           endif
           if(blocklabel.eq.ablock) exit
           call find_end_of_block(.FALSE.)
           if(ifail.ne.0) go to 9890
         enddo
         if(blockformat.ne.1)then
           if(fileflag.eq.1)then
             qfile=qinfile
           else
             qfile=qsubfile
           endif
           write(amessage,15) trim(blocklabel),trim(qfile)
15         format('Programming error in call to subroutine IDS_READ_BLOCK_KEYWORDS: ',  &
           'block ',a,' located in file ',a,' is expected to have KEYWORD format if ',  &
           'read using this subroutine.')
           go to 9890
         endif

! -- Keywords are now read from the block.

         do
           call find_next_keyword()
           if(ifail.ne.0) go to 9890
           if(keyword.eq.' ') go to 50
           do i=1,numlist
             akeyword=adjustl(keywordname(i))
             call uppercase(akeyword)
             if(keyword.eq.akeyword)then
               if(keywordline(i).eq.0)then
                 keywordline(i)=iline
                 keywordstring(i)=keystring      ! I would prefer that this had not been converted to upper case
                 go to 40
               else
                 call BDP_INT2CHAR(iline,aline)
                 write(amessage,20) trim(keyword),trim(ablock),trim(aline),trim(qfile)
20               format(a,' keyword duplicated in block ',a,' at line ',a,' of file ',a,'.')
                 go to 9890
               endif
             endif
           enddo
40         continue
         enddo

50       continue
         return

9890     continue
         ifail_out=1

       end subroutine ids_read_block_keywords




       subroutine ids_read_block_table(ifail_out,blockname_in,numcol,columnname,columnstring,line,filename)

! -- Subroutine IDS_READ_BLOCK_TABLE returns the contents of specified columns of a table for one particular
!    line of that table.

         implicit none

         integer, intent(out)            :: ifail_out
         character (len=*), intent(in)   :: blockname_in
         integer, intent(in)             :: numcol
         character (len=*), intent(in)   :: columnname(numcol)
         character (len=*), intent(out)  :: columnstring(numcol)
         integer, intent(out)            :: line
         character (len=*), intent(out)  :: filename

         integer                         :: icol,jcol,n
         character (len=LABELWIDTH)      :: ablock
         character (len=COLWIDTH)        :: acolumn

! -- Security

         if(block_status_call.eq.0)then
           write(amessage,2)
2          format('Programming error in call to subroutine IDS_READ_BLOCK_TABLE: subroutine ',  &
           'IDS_BLOCK_STATUS has not yet been called.')
           go to 9890
         endif

! -- Initialisation

         ifail_out=0
         columnstring=' '  ! an array
         if(numcol.gt.MAXCOL)then
           write(amessage,5)
5          format('Programming error in call to subroutine IDS_READ_BLOCK_TABLE: NUMCOL argument greater ', &
           'than MAXCOL. Increase MAXCOL and re-compile program.')
           go to 9890
         endif

! -- We don't bother checking for duplicate supplied column names as this subroutine may be called
!    many times.

! -- We check where we are in the reading_of_file process.

         ablock=adjustl(blockname_in)
         call uppercase(ablock)

         if(ablock.ne.blockname_in_last) endoftableblockflag=0
         if((ablock.ne.blockname_in_last).or.(endoftableblockflag.eq.1))then

           if(position.eq.2)then
             write(amessage,7) trim(blockname_in)
7            format('Programming error in call to subroutine IDS_READ_BLOCK_TABLE: an attempt ',       &
             'was made to process "',a,'" block when processing of previous block was not complete. ', &
             'If the previously processed block was of TABLE format, it is possible that subroutine ', &
             'IDS_READ_BLOCK_TABLE was called fewer times than there are rows in the (collective) ',   &
             'table(s) for this block.')
             go to 9890
           endif
           endoftableblockflag=0

! -- We find the requested block.

10         continue
           do
             call find_next_block()
             if(ifail.ne.0) go to 9890
             if(blocklabel.eq.' ')then
               write(amessage,12) trim(ablock),trim(qinfile)
12             format('Error in call to subroutine IDS_READ_BLOCK_TABLE: cannot find ',  &
               '(all incidences of) ',a,' block through forward search in file ',a,'. It is ', &
               'likely that this subroutine has been called more times than the number of rows '  &
               'in the (collective) table(s) with this label.')
               go to 9890
             endif
             if(blocklabel.eq.ablock) then
               exit
             else
               call find_end_of_block(.FALSE.)
             endif
           enddo
           if(blockformat.ne.2)then
             if(fileflag.eq.1)then
               qfile=qinfile
             else
               qfile=qsubfile
             endif
             write(amessage,15) trim(blocklabel),trim(qfile)
15           format('Programming error in call to subroutine IDS_READ_BLOCK_TABLE: ',  &
             'block ',a,' located in file ',a,' is expected to have TABLE format if ',  &
             'read using this subroutine.')
             go to 9890
           endif
           call read_table_info()
           if(ifail.ne.0) go to 9890
           if(nrow.eq.0)then
             call find_end_of_block(.FALSE.)
             if(ifail.ne.0) go to 9890
             blockname_in_last=ablock
             go to 10
           endif

! -- It is inefficient to repeat the following code on each occasion that a line of data is read. If it
!    could be guaranteed that the user's column list would be the same on each occasion that a line is
!    read from a table, this code would need to be executed only on entry to each block.

           colexport=0                        ! an array
           maxcolexport=0
           do icol=1,ncol
             do jcol=1,numcol
               acolumn=adjustl(columnname(jcol))
               call uppercase(acolumn)
               if(acolumn.eq.columnlabel(icol))then
                 colexport(icol)=jcol
                 maxcolexport=icol
                 go to 20
               endif
             enddo
20           continue
           enddo
         endif

! -- A line is read.

         tablerow=tablerow+1
         call read_next_line()
         if(ifail.ne.0) go to 9890
         if(maxcolexport.ne.0)then
           do icol=1,maxcolexport
             call find_next_word()         ! make sure this works if we go past the end of a line
             if(ifail.ne.0) go to 9890
             n=colexport(icol)
             if(n.ne.0) columnstring(n)=nextword
           enddo
         endif
         if(fileflag.eq.1)then
           filename=infile
           line=iline
         else
           filename=subfile
           line=jline
         endif

! -- If this is the last line in this table or subtable, go to the end of the table.

         if(tablerow.eq.nrow)then
           call find_end_of_block(.FALSE.)
           endoftableblockflag=1
         endif
         blockname_in_last=ablock

         return

9890     continue
         ifail_out=1

       end subroutine ids_read_block_table



       subroutine ids_get_message_string(amessage_out)

         implicit none
         character (len=*), intent(out) :: amessage_out
         amessage_out=amessage
         return

       end subroutine ids_get_message_string




       subroutine find_next_keyword()

! -- Subroutine FIND_NEXT_KEYWORD finds the next keyword and its associated string in a keyword block.

         implicit none

         integer           :: newlineflag,n,i,nn,iflag1,iflag2,iflag3,iflag4
         character (len=4) :: atemp4

! -- Initialisation

         ifail=0
         keyword=' '
         keystring=' '
         newlineflag=0

         if(newblockflag.eq.1)then
           newblockflag=0
           call read_next_line()
           if(ifail.ne.0) return
           if(linepos.eq.-999)then
             call BDP_INT2CHAR(iline,aline)
             write(amessage,10) trim(qinfile),trim(blocklabel)
10           format('Unexpected end to file ',a,' encountered while reading ',a,' block.')
             go to 9890
           endif
           newlineflag=1
         endif

! -- Is there an "=" character on the line?

12       continue
         if(newlineflag.eq.1)then
           atemp4=cline(1:4)
           call uppercase(atemp4)
           if(atemp4.eq.'END ')then
             endblockflag=1
             call find_end_of_block(.FALSE.)
             return
           endif
           newlineflag=0

! -- We check for the location of spurious commas.

           if(index(cline,comma).ne.0)then
             iflag1=0
             iflag2=0
             iflag3=0
             iflag4=0
             do i=1,len_trim(cline)
               if(cline(i:i).eq.squote)then
                 iflag3=abs(iflag3-1)
                 iflag1=0
                 iflag2=0
               else if(cline(i:i).eq.dquote)then
                 iflag4=abs(iflag4-1)
                 iflag1=0
                 iflag2=0
               else
                 if((iflag3.eq.0).and.(iflag4.eq.0))then
                   if(cline(i:i).eq.comma) then
                     if(iflag2.ne.0)then
                       call BDP_INT2CHAR(iline,aline)
                       write(amessage,22) trim(aline),trim(qinfile)
22                     format('Comma in illegal position on line ',a,' of file ',a,'.')
                       go to 9890
                     else
                       iflag1=1
                       cycle
                     endif
                   else if(cline(i:i).eq.equality)then
                     if(iflag1.eq.1)then
                       call BDP_INT2CHAR(iline,aline)
                       write(amessage,22) trim(aline),trim(qinfile)
                       go to 9890
                     else
                       iflag2=1
                       cycle
                     endif
                   else
                     if((cline(i:i).eq.space).or.(cline(i:i).eq.tab))then
                       cycle
                     else
                       iflag1=0
                       iflag2=0
                     endif
                   endif
                 endif
               endif
             enddo
           endif

         endif
         n=index(cline(linepos+1:),equality)
         if(n.eq.0)then
           call repchar(cline(linepos+1:),comma,space)
           call repchar(cline(linepos+1:),tab,space)
           if(cline(linepos+1:).eq.' ')then
             call read_next_line()
             if(ifail.ne.0) return
             newlineflag=1
             go to 12
           endif
           call BDP_INT2CHAR(iline,aline)
           write(amessage,20) trim(aline),trim(qinfile)
20         format('Non-comment, non-keyword text found on line ',a,' of file ',a,'.')
           go to 9890
         endif

! -- The "=" symbol is deleted and we read the next two words.

         nn=linepos+n
         cline(nn:nn)=' '
         call find_next_word()
         if(ifail.ne.0) return
         if(nextword.eq.' ') then
           call BDP_INT2CHAR(iline,aline)
           write(amessage,30) trim(aline),trim(qinfile)
           go to 9890
         endif
         if(linepos.gt.nn)then
           call BDP_INT2CHAR(iline,aline)
           write(amessage,25) trim(aline),trim(qinfile)
25         format('No keyword precedes "=" symbol at line ',a,' of file ',a,'.')
           go to 9890
         endif
         keyword=nextword
         call uppercase(keyword)
         do i=linepos+1,nn
           if((cline(i:i).ne.space).and.(cline(i:i).ne.tab))then
             call BDP_INT2CHAR(iline,aline)
             write(amessage,30) trim(aline),trim(qinfile)
30           format('Improper keyword syntax at line ',a,' of file ',a,'.')
             go to 9890
           endif
         enddo
         call find_next_word()
         if(ifail.ne.0) return
         if(nextword.eq.' ')then
           call BDP_INT2CHAR(iline,aline)
           write(amessage,50) trim(keyword),trim(aline),trim(qinfile)
50         format('No value provided for ',a,' keyword at line ',a,' of file ',a,'.')
           go to 9890
         endif
         keystring=nextword
         n=len_trim(keystring)
         if(keystring(n:n).eq.equality)then
           call BDP_INT2CHAR(iline,aline)
           write(amessage,25) trim(aline),trim(qinfile)
           go to 9890
         endif

         return

9890     ifail=1
         return

       end subroutine find_next_keyword




       subroutine find_next_block()

! -- Subroutine FIND_NEXT_BLOCK finds the next block in the input file or subsidiary file.
!    If it encounters a block with a blockformat of FILE, it moves to the first block within
!    it which actually contains data.

         implicit none

         integer             :: ierr
         character (len=4)   :: atemp4

! -- Initialisation

         ifail=0
         blocklabel=' '
         emptyfileblockflag=0

5        continue
         call read_next_line()
         if(ifail.ne.0) return
         if(linepos.eq.-999)then
           if(fileflag.eq.1)return
           close(unit=junit)
           fileflag=1
           call read_next_line()
           if(ifail.ne.0) return
           if(linepos.eq.-999)then
             write(amessage,10) trim(qinfile),trim(blocklabelmain)
10           format('Unexpected end encountered to file ',a,' while reading ',a,' block.')
             go to 9890
           endif
           atemp4=cline(1:4)
           call uppercase(atemp4)
           if(atemp4.eq.'END ')then
             endblockflag=1
             call find_end_of_block(.FALSE.)
             blockformatmain=0
             blocklabelmain=' '
             go to 5
           else
             call find_next_word()
             if(ifail.ne.0) return
             subfile=nextword
             call addquote(subfile,qsubfile)
             open(unit=junit,file=subfile,status='old',iostat=ierr)
             if(ierr.ne.0)then
               call BDP_INT2CHAR(iline,aline)
               write(amessage,15) trim(qsubfile),trim(aline),trim(qinfile)
15             format('Cannot open file ',a,' cited at line ',a,' of file ',a,'.')
               go to 9890
             endif
             fileflag=2
             jline=0
             go to 5
           endif
         endif
         call find_next_word()
         call uppercase(nextword)
         if(fileflag.eq.1)then
           qfile=qinfile
           kline=iline
         else
           qfile=qsubfile
           kline=jline
         endif
         if(nextword.ne.'BEGIN')then
           call BDP_INT2CHAR(kline,aline)
           write(amessage,20) trim(aline),trim(qfile)
20         format('First word on line ',a,' of file ',a,' expected to be "BEGIN".')
           go to 9890
         endif
         call find_next_word()
         if(ifail.ne.0) return
         if(nextword.eq.' ') then
           call BDP_INT2CHAR(kline,aline)
           write(amessage,30) trim(aline),trim(qfile)
30         format('No blocklabel supplied after "BEGIN" statement at line ',a,' of file ',a,'.')
           go to 9890
         endif
         blocklabel=nextword
         call uppercase(blocklabel)
         if(fileflag.eq.1)blocklabelmain=blocklabel
         call find_next_word()
         if(ifail.ne.0) return
         if(nextword.eq.' ')then
           blockformat=1
         else if(nextword(1:7).eq.'KEYWORD')then
           blockformat=1
         else if(nextword(1:5).eq.'TABLE')then
           blockformat=2
         else if(nextword(1:4).eq.'FILE')then
           blockformat=3
         else
           call BDP_INT2CHAR(kline,aline)
           write(amessage,40) trim(aline),trim(qfile)
40         format('Unrecognised BLOCKFORMAT at line ',a,' of file ',a,'.')
           go to 9890
         endif
         if(fileflag.eq.1)then
           blockformatmain=blockformat
         else
           if(blockformat.ne.2)then
             call BDP_INT2CHAR(kline,aline)
             write(amessage,45) trim(qsubfile),trim(blocklabelmain),trim(qinfile),   &
             trim(aline),trim(qsubfile)
45           format('Only TABLE blocks are allowed in file ',a,' cited in ',a,      &
             ' block of file ',a,'. Error occurs at line ',a,' of file ',a,'.')
             go to 9890
           endif
           if(blocklabel.ne.blocklabelmain)then
             write(amessage,50) trim(blocklabel),trim(qsubfile),trim(blocklabelmain),trim(qinfile)
50           format('Block ',a,' in file ',a,' has different label from ', &
             'block ',a,' of file ',a,' in which it is cited.')
             go to 9890
           endif
         endif
         position=blockformat
         if(position.eq.2) position=0       ! This is assigned after the table header is read.
         if(blockformat.eq.3)then
           call read_next_line()
           if(ifail.ne.0) go to 9890
           atemp4=cline(1:4)
           call uppercase(atemp4)
           if(atemp4.eq.'END ')then
             emptyfileblockflag=1
             endblockflag=1
             return
           else
             call find_next_word
             if(ifail.ne.0) go to 9890
             subfile=nextword
             call addquote(subfile,qsubfile)
             open(unit=junit,file=subfile,status='old',iostat=ierr)
             if(ierr.ne.0)then
               call BDP_INT2CHAR(iline,aline)
               write(amessage,15) trim(qsubfile),trim(aline),trim(qfile)
               go to 9890
             endif
             fileflag=2
             jline=0
             go to 5
           endif
         endif

100      continue
         newblockflag=1
         go to 9999

9890     continue
         ifail=1

9999     return

       end subroutine find_next_block



       subroutine find_end_of_block(countrows)

! -- Subroutine FIND_END_OF_BLOCK finds the end of a block.

         implicit none

         logical, intent(in) :: countrows
         integer             :: icount

! -- Initialisation

         if(fileflag.eq.1)then
           qfile=qinfile
         else
           qfile=qsubfile
         endif
         icount=0

! -- The end of the block is found.

         if(endblockflag.eq.1)then
           call find_next_word()
           if(ifail.ne.0) return
           endblockflag=0
         else
           do
             call read_next_line()
             if(ifail.ne.0)return
             if(linepos.eq.-999) go to 9000
             call find_next_word()
             if(ifail.ne.0) return
             call uppercase(nextword)
             if(nextword.eq.'END') exit
             icount=icount+1
             if(countrows)then
               if(blockformat.eq.2)then
                 if(icount.gt.nrow)then
                   if(fileflag.eq.1)then
                     kline=iline
                   else
                     kline=jline
                   endif
                   call BDP_INT2CHAR(kline,aline)
                   write(amessage,10) trim(aline),trim(qfile)
10                 format('Table block END expected at line ',a,' of file ',a,'.')
                   go to 9890
                 endif
               endif
             endif
           enddo
         endif
         call find_next_word()
         if(ifail.ne.0) return
         if(nextword.eq.' ') go to 9100
         call uppercase(nextword)
         if(nextword.ne.blocklabel) go to 9100

! -- For a table block, a check is made that the number of rows is the same as stated in header.

         if(countrows)then
           if(blockformat.eq.2)then
             if(icount.ne.nrow)then
               write(amessage,50) trim(blocklabel),trim(qfile)
50             format('Number of data rows in ',a,' block of file ',a,' does not agree with NROW ',  &
               'in block table header.')
               go to 9890
             endif
           endif
         endif

! -- Variables are re-set.

         blocklabel=' '
         blockformat=0
         if(fileflag.eq.1)then
           blocklabelmain=' '
           blockformatmain=0
         endif
         position=0


         return

9000     write(amessage,9010) trim(qfile),trim(blocklabel)
9010     format('Unexpected end encountered to file ',a,' while reading ',a,' block.')
         go to 9890

9100     if(fileflag.eq.1)then
           kline=iline
         else
           kline=jline
         endif
         call BDP_INT2CHAR(kline,aline)
         write(amessage,9110) trim(aline),trim(qfile)
9110    format('Label at line ',a,' of file ',a,' is absent, or does not match label at ',  &
         'start of block.')
         go to 9890

9890     ifail=1
         return

       end subroutine find_end_of_block




       subroutine find_next_word()

! -- Subroutine FIND_NEXT_WORD finds the word following the current word. The "cursor"
!    is placed at the end of this word.

         implicit none

         integer            :: lpos,nb,i,j,k
         character (len=1)  :: aa,bb

! -- Initialisation.

         ifail=0
         nextword=' '

! -- Search for next word.

         nb=len_trim(cline)
         if(linepos.eq.nb) return
         do i=linepos+1,nb
           if((cline(i:i).ne.tab).and.(cline(i:i).ne.space)) exit
         enddo
         aa=cline(i:i)
         if(aa.eq.comma)then
           linepos=i
           return
         endif
         if((aa.eq.squote).or.(aa.eq.dquote))then
           do j=i+1,nb
             if(cline(j:j).eq.aa) go to 20
           enddo
           go to 9000
20         continue
           nextword=cline(i+1:j-1)
           if(j.eq.nb)then
             linepos=j
             return
           else
             bb=cline(j+1:j+1)
             if((bb.ne.space).and.(bb.ne.tab).and.(bb.ne.comma))then
               if(fileflag.eq.1)then
                 call BDP_INT2CHAR(iline,aline)
                 qfile=qinfile
               else
                 call BDP_INT2CHAR(jline,aline)
                 qfile=qsubfile
               endif
               write(amessage,30) trim(aline),trim(qfile)
30             format('A space, tab or comma must follow a quote: error on line ',a,' of file ',a,'.')
               go to 9890
             endif
             do k=j+1,nb
               if((cline(k:k).ne.tab).and.(cline(k:k).ne.space))then
                 if(cline(k:k).eq.comma)then
                   linepos=k
                 else
                   linepos=j
                 endif
                 return
               endif
             enddo
             linepos=j
             return
           endif
         else
           do j=i+1,nb
             if((cline(j:j).eq.space).or.(cline(j:j).eq.tab).or.(cline(j:j).eq.comma))then
               aa=cline(j:j)
               if(aa.eq.comma)then
                 if(j.eq.i+1)then
                   nextword=' '
                 else
                   nextword=cline(i:j-1)
                   if(lastcharquote()) go to 9000
                 endif
                 linepos=j
                 return
               else
                 nextword=cline(i:j-1)
                 if(lastcharquote()) go to 9000
                 do k=j+1,nb
                   if((cline(k:k).ne.tab).and.(cline(k:k).ne.space))then
                     if(cline(k:k).eq.comma)then
                       linepos=k
                     else
                       linepos=j-1
                     endif
                     return
                   endif
                 enddo
                 linepos=j-1
                 return
               endif
             endif
           enddo
           if(cline(nb:nb).eq.comma)then
             if(nb.eq.i+1)then
               nextword=' '
             else
               nextword=cline(i:nb-1)
               if(lastcharquote()) go to 9000
             endif
             linepos=nb
           else
             nextword=cline(i:nb)
             if(lastcharquote()) go to 9000
             linepos=nb
           endif
         endif

         return

9000     continue
         ifail=1
         if(fileflag.eq.1)then
           call BDP_INT2CHAR(iline,aline)
           qfile=qinfile
         else
           call BDP_INT2CHAR(jline,aline)
           qfile=qsubfile
         endif
         write(amessage,9010) trim(aline),trim(qfile)
9010     format('Unbalanced quotes at line ',a,' of file ',a,'.')
         go to 9890

9890     ifail=1
         return

       end subroutine find_next_word



       logical function lastcharquote()

! -- Subroutine lastcharquote checks for an orphaned quote at the end of a word.

         implicit none
         integer            :: n
         character (len=1)  :: aa

         lastcharquote=.false.
         if(nextword.ne.' ')then
           n=len_trim(nextword)
           aa=nextword(n:n)
           if((aa.eq.squote).or.(aa.eq.dquote))then
             lastcharquote=.true.
           endif
         endif

         return
       end function lastcharquote




       subroutine read_next_line()

! -- Subroutine READ_NEXT_LINE reads the next line of the main or subsidiary input file.
!    It keeps reading until it finds a non-empty, non-comment line. It also removes commented
!    sections from the end of a line.

         implicit none

         integer   :: i

         ifail=0
         cline=' '

         do
           if(fileflag.eq.1)then
             iline=iline+1
             kunit=iunit
           else
             jline=jline+1
             kunit=junit
           endif
           read(kunit,'(a)',err=9000,end=100) cline
10         continue
           cline=adjustl(cline)

! -- Tabs are removed.

           call repchar(cline,char(9),' ')
           if(cline(1:1).eq.'#') cycle
           i=index(cline,'#')
           if(i.ne.0) cline(i:)=' '
           if(cline.eq.' ') cycle
           exit
         enddo
         linepos=0
         return

100      linepos=-999
         return

9000     continue
         if(fileflag.eq.1)then
           call BDP_INT2CHAR(iline,aline)
           qfile=qinfile
         else
           call BDP_INT2CHAR(jline,aline)
           qfile=qsubfile
         endif
         write(amessage,9010) trim(aline),trim(qfile)
9010     format('Error encountered when reading line ',a,' of file ',a,'.')
         go to 9890

9890     ifail=1

       end subroutine read_next_line




       subroutine read_table_info()

! -- Subroutine READ_TABLE_INFO reads the header information to a table.

         implicit none

         integer     :: n,k,i,itemp,j,nb
         character (len=4)   :: acolrow
         character (len=5)   :: atemp1
         character (len=100) :: atemp

! -- Initialisation

         ifail=0
         ncol=0
         nrow=0
         columnlabel=' '        ! An array
         if(fileflag.eq.1)then
           qfile=qinfile
         else
           qfile=qsubfile
         endif

! -- The header is read

         call read_next_line()
         if(ifail.ne.0) return
         if(linepos.eq.-999) go to 9000
         call uppercase(cline)
         if(fileflag.eq.1)then
           kline=iline
         else
           kline=jline
         endif
         n=index(cline,'COLUMNLABELS')
         if(n.eq.0)then
           call BDP_INT2CHAR(kline,aline)
           write(amessage,20) trim(aline),trim(qfile)
20         format('COLUMNLABELS keyword expected on line ',a,' of file ',a,'.')
           go to 9890
         endif
         cline(n:n+11)=' '
         n=index(cline,'DATAFILES')
         if(n.ne.0)then
           call BDP_INT2CHAR(kline,aline)
           write(amessage,30) trim(aline),trim(qfile)
30         format('DATAFILES keyword at line ',a,' of file ',a,' is not supported.')
           go to 9890
         endif

! -- The number of rows and columns is read.

         do
           n=index(cline,' =')
           if(n.eq.0) go to 50
           cline(n:n+1)='= '
         enddo
50       continue
         do k=1,2
           if(k.eq.1)then
             acolrow='NROW'
           else
             acolrow='NCOL'
           endif
           n=index(cline,trim(acolrow))
           if(n.eq.0)then
             call BDP_INT2CHAR(kline,aline)
             write(amessage,60) trim(acolrow),trim(aline),trim(qfile)
60           format('Cannot find ',a,' keyword on line ',a,' of file ',a,'.')
             go to 9890
           endif
           n=index(cline,equality)
           if(n.eq.0)then
             call BDP_INT2CHAR(kline,aline)
             write(amessage,70) trim(acolrow),trim(aline),trim(qfile)
             go to 9890
           endif
           atemp1=trim(acolrow)//'='
           n=index(cline,trim(atemp1))
           if(n.eq.0)then
             call BDP_INT2CHAR(kline,aline)
             write(amessage,70) trim(acolrow),trim(aline),trim(qfile)
70           format('Improper syntax for ',a,' keyword on line ',a,' of file ',a,'.')
             go to 9890
           endif
           cline(n+4:n+4)=' '
           atemp=cline(n+5:n+100)
           if(atemp.eq.' ') go to 9100
           atemp=adjustl(atemp)
           nb=len_trim(atemp)
           do i=1,nb
             if((atemp(i:i).eq.space).or.(atemp(i:i).eq.tab).or.(atemp(i:i).eq.comma)) then
               atemp(i:)=' '
               go to 80
             endif
           enddo
           i=nb+1
80         i=i-1
           call char2int(atemp(1:i),itemp)
           if(ifail.ne.0) go to 9100
           if(itemp.lt.0)then
             call BDP_INT2CHAR(kline,aline)
             write(amessage,85) trim(acolrow),trim(aline),trim(qfile)
85           format('Negative value for ',a,' supplied on line ',a,' of file ',a,'.')
             go to 9890
           endif
           if(k.eq.1)then
             nrow=itemp
           else
             ncol=itemp
           endif
         enddo
         if(ncol.gt.MAXCOL)then
           call BDP_INT2CHAR(MAXCOL,atemp1)
           call BDP_INT2CHAR(kline,aline)
           write(amessage,90) trim(atemp1),trim(aline),trim(qfile)
90         format('Number of columns exceeds maximum allowed number of ',a,' at line ',a,   &
           ' of file 'a,'. Alter MAXCOL and re-compile program.')
           go to 9890
         endif
         if(index(cline,'GROUPNAME').ne.0)then
           call BDP_INT2CHAR(kline,aline)
           write(amessage,94) trim(aline),trim(qfile)
94         format('GROUPNAME keyword not supported on line ',a,' of file ',a,'.')
           go to 9890
         endif
         if(index(cline,'=').ne.0)then
           call BDP_INT2CHAR(kline,aline)
           write(amessage,95) trim(aline),trim(qfile)
95         format('Keywords other than NCOL and NROW are not supported on line ',a,' of file ',a,'.')
           go to 9890
         endif

! -- The column headers are read.

         call read_next_line()
         if(ifail.ne.0) return
         if(linepos.eq.-999) go to 9000
         call uppercase(cline)
         if(fileflag.eq.1)then
           kline=iline
         else
           kline=jline
         endif
         do i=1,ncol
           call find_next_word()
           if(ifail.ne.0) return
           if(nextword.eq.' ')then
             call BDP_INT2CHAR(kline,aline)
             write(amessage,100) trim(aline),trim(qfile)
100          format('Missing column label(s) at line ',a,' of file ',a,'.')
             go to 9890
           endif
           columnlabel(i)=nextword
           if(i.gt.1)then
             do j=1,i-1
               if(columnlabel(j).eq.nextword)then
                 call BDP_INT2CHAR(kline,aline)
                 write(amessage,120) trim(aline),trim(qfile)
120              format('Duplicated column label at line ',a,' of file ',a,'.')
                 go to 9890
               endif
             enddo
           endif
         enddo

         tablerow=0
         position=2
         return

9000     write(amessage,9010) trim(qfile),trim(blocklabel)
9010     format('Unexpected end to file ',a,' while reading ',a,' block.')
         go to 9890
9100     call BDP_INT2CHAR(kline,aline)
         write(amessage,9110) trim(acolrow),trim(aline),trim(qfile)
9110     format('Cannot read ',a,' specifier from line ',a,' of file ',a,'.')
         go to 9890

9890     continue
         ifail=1

         return

       end subroutine read_table_info



       subroutine BDP_INT2CHAR(ival,astring)

! -- Subroutine BDP_INT2CHAR writes an integer to a string.

         implicit none
         integer, intent(in)            :: ival
         character (len=*), intent(out) :: astring
         character (len=7)              :: afmt

         afmt='(i    )'
         write(afmt(3:6),'(i4)') len(astring)
         write(astring,afmt) ival
         astring=adjustl(astring)
         return

       end subroutine BDP_INT2CHAR



       subroutine char2int(astring,inum)

! -- Subroutine CHAR2INT reads an integer from a string.

         character (len=*), intent(in)    :: astring
         integer, intent(out)             :: inum
         character (len=7)                :: afmt

         ifail=0
         afmt='(i    )'
         write(afmt(3:6),'(i4)') len(astring)
         read(astring,afmt,err=100) inum
         return

100      ifail=1
         return

       end subroutine char2int



       subroutine initialise()

! -- Subroutine INITIALISE initialises global variables employed by the module.

         iline=0
         jline=0
         position=0
         tablerow=0
         newblockflag=0
         endblockflag=0
         emptyfileblockflag=0
         endoftableblockflag=0
         fileflag=1

         blockname_in_last=' '
         blocklabel=' '
         blocklabelmain=' '
         nextword=' '
         keyword=' '
         keystring=' '

         return
       end subroutine initialise


       subroutine repchar(astring,char1,char2)

! -- Subroutine REPCHAR replaces one character by another throughout a string.

         character (len=*), intent(inout)  :: astring
         character (len=1), intent(in)     :: char1
         character (len=1), intent(in)     :: char2

         integer :: i

         do i=1,len_trim(astring)
           if(astring(i:i).eq.char1) astring(i:i)=char2
         enddo

         return
       end subroutine repchar



       subroutine uppercase(astring)

! -- Subroutine UPPERCASE converts a string to upper case.

         character (len=*), intent(inout)        :: astring
         character                               :: alo, ahi
         integer                                 :: inc,i

         alo='a'
         ahi='z'
         inc=iachar('A')-iachar('a')

	 do i=1,len_trim(astring)
	   if((astring(i:i).ge.alo).and.(astring(i:i).le.ahi)) &
	   astring(i:i)=achar(iachar(astring(i:i))+inc)
	 enddo

	 return

       end subroutine uppercase



       subroutine addquote(afile,aqfile)

! -- Subroutine ADDQUOTE adds quotes to a filename if it has a space in it.

         character (len=*), intent(in)  :: afile
         character (len=*), intent(out) :: aqfile
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

       end subroutine addquote


end module jupiter_input_data_support



