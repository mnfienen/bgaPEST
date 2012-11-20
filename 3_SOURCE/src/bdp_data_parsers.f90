module bdp_data_parsers

        use bayes_pest_control
        use jupiter_input_data_support
        use utilities
        use error_message

contains

!********  subroutine bdp_read_cv_algorithmic(BL,cv_A,inunit,retmsg)
       subroutine bdp_read_cv_algorithmic(BL,cv_A,inunit,retmsg)
       ! SUBROUTINE to read in algorithmic control keywords
         use bayes_pest_control
         use jupiter_input_data_support
       ! DECLARATIONS
         implicit none
         type(cv_algorithmic), intent(inout)  :: cv_A
         type(tp_block),       intent(inout)  :: BL(NUM_BLOCK) 
         integer,              intent(in)     :: inunit
         character (len=ERRORWIDTH), intent(inout)  :: retmsg
         integer                              :: i  !local counter
         integer                              :: ifail
         if (BL(1)%numrows .EQ. -1) then ! check for status as KEYWORDS
          call bdp_init_keyword_vars(BL,1)
          call ids_read_block_keywords(ifail,BL(1)%label,BL(1)%numkw, &
             & BL(1)%keywords,BL(1)%keywordstring,BL(1)%keywordline)
          if (ifail.ne.0) then
           call ids_get_message_string(retmsg)
           call utl_writmess(6,retmsg)
          endif
         elseif (BL(1)%numrows .EQ. 0) then
          retmsg = 'All default values accepted for cv_algorithmic block'       
         endif
       ! Parse out the values as appropriate
         do i=1,BL(1)%numkw
          if(BL(1)%keywordline(i) .ne. 0) then
            select case (trim(BL(1)%keywords(i)))
             case ('structural_conv')
               call drealread(ifail,BL(1)%keywordstring(i), cv_A%structural_conv)
             case ('phi_conv')
               call drealread(ifail,BL(1)%keywordstring(i), cv_A%phi_conv)
             case ('it_max_structural')
               call intread(ifail,BL(1)%keywordstring(i), cv_A%it_max_structural)
             case ('it_max_phi')
               call intread(ifail,BL(1)%keywordstring(i), cv_A%it_max_phi)
             case ('it_max_bga')
               call intread(ifail,BL(1)%keywordstring(i), cv_A%it_max_bga)
             case ('bga_conv')
               call drealread(ifail,BL(1)%keywordstring(i), cv_A%bga_conv)  
             case ('linesearch')
               call intread(ifail,BL(1)%keywordstring(i), cv_A%lns_flag)
             case ('it_max_linesearch')
               call intread(ifail,BL(1)%keywordstring(i), cv_A%it_max_lns)
             case ('theta_cov_form')
               call intread(ifail,BL(1)%keywordstring(i), cv_A%theta_cov_form)
             case ('Q_compression_flag')
               call intread(ifail,BL(1)%keywordstring(i), cv_A%Q_compression_flag)
             case ('deriv_mode')
               call intread(ifail,BL(1)%keywordstring(i), cv_A%deriv_mode)
             case('posterior_cov_flag')
               call intread(ifail,BL(1)%keywordstring(i), cv_A%post_cov_flag)
             case('par_anisotropy')
               call intread(ifail,BL(1)%keywordstring(i), cv_A%par_anisotropy)
             case ('jacobian_format')
               cv_A%jacobian_format = trim(adjustl(BL(1)%keywordstring(i)))
               call UTL_CASETRANS(cv_A%jacobian_format,'lo')
             case ('jacobian_file')
               cv_A%jacfle = trim(adjustl(BL(1)%keywordstring(i)))
            end select
           endif
          enddo
          if(cv_A%deriv_mode.eq.0) then
             cv_A%jacobian_format = 'binary'
             cv_A%jacfle = 'scratch.jco'
          endif
          if (cv_A%bga_conv.eq.0.01) then
             cv_A%bga_conv = cv_A%phi_conv*10.
          endif
 end subroutine bdp_read_cv_algorithmic
     
     
!********  subroutine bdp_read_cv_prior_mean(BL,cv_PM,inunit)
       subroutine bdp_read_cv_prior_mean(BL,cv_PM,inunit,retmsg)
       ! SUBROUTING to read in algorithmic control keywords
         use bayes_pest_control
         use jupiter_input_data_support
       ! DECLARATIONS
         implicit none
         type(cv_prior_mean),  intent(inout)  :: cv_PM
         type(tp_block),       intent(inout)  :: BL(NUM_BLOCK) 
         integer,              intent(in)     :: inunit
         character (len=ERRORWIDTH), intent(inout)  :: retmsg
         integer                              :: i  !local counter
         integer                              :: ifail
         if (BL(2)%numrows .EQ. -1) then ! check for status as KEYWORDS
          call bdp_init_keyword_vars(BL,2)
          call ids_read_block_keywords(ifail,BL(2)%label,BL(2)%numkw, &
           & BL(2)%keywords,BL(2)%keywordstring,BL(2)%keywordline)
          if (ifail.ne.0) then
           call ids_get_message_string(retmsg)
           call utl_writmess(6,retmsg)
          endif
         elseif (BL(2)%numrows .EQ. 0) then
          retmsg = 'All default values accepted for cv_prior_mean block'       
         endif
       ! Parse out the values as appropriate
         do i=1,BL(2)%numkw
          if(BL(2)%keywordline(i) .ne. 0) then
           select case (trim(BL(2)%keywords(i)))
            case ('prior_betas')
              call intread(ifail,BL(2)%keywordstring(i), cv_PM%betas_flag)              
            case ('beta_cov_form')
              call intread(ifail,BL(2)%keywordstring(i), cv_PM%Qbb_form)
           end select
          endif
         enddo
         
      ! run a quick error check on prior 
end subroutine bdp_read_cv_prior_mean
     
!********  subroutine bdp_read_data_prior_mean (BL,cv_PM,d_PM,inunit,errmsg)
       subroutine bdp_read_data_prior_mean (BL,cv_PM,d_PM,cv_PAR,inunit,retmsg)
       ! SUBROUTING to read in algorithmic control keywords
         use bayes_pest_control
         use jupiter_input_data_support
       ! DECLARATIONS
         implicit none
         type(d_prior_mean),  intent(inout)        :: d_PM
         type(cv_param),       intent(inout)       :: cv_PAR
         type(cv_prior_mean),  intent(inout)       :: cv_PM
         type(tp_block),       intent(inout)       :: BL(NUM_BLOCK) 
         integer,              intent(in)          :: inunit
         character (len=ERRORWIDTH), intent(inout) :: retmsg
         character (len=COLWIDTH), pointer     :: columnname(:)
         character (len=COLWIDTH), pointer     :: columnstring(:)         
         character (len=200)                       :: filename
         character (len=10)                        :: tmp_str
         integer                                   :: numcol
         integer                                   :: i,j  !local counter
         integer                                   :: ifail, line
         
         if (BL(3)%numrows .EQ. 0) then
          retmsg = 'No input information provided in data_prior_mean block'  
          return     
         endif 
         cv_PAR%p = BL(3)%numrows !Define the number of betas
         allocate(d_PM%Partrans(cv_PAR%p))
         d_PM%Partrans   = UNINIT_INT  ! array
   
   if (cv_PM%betas_flag.eq.0) then      
      call bdp_alloc_Partrans(d_PM,columnname,columnstring,3,BL(3)%numrows)
       !**********************************************************************************
       !*************     PARAMETER TRANSFORMATION FLAG     ******************************
         do j=1,BL(3)%numrows
           call ids_read_block_table(ifail,BL(3)%label,3,columnname,columnstring,line,filename)
           tmp_str = trim(adjustl(columnstring(2)))  
            call bpc_casetrans(tmp_str,'lo')
            if (tmp_str.eq.'log') then
                d_PM%Partrans(j)=1
            elseif (tmp_str.eq.'power') then
                d_PM%Partrans(j)=2
            elseif (tmp_str.eq.'none') then
                d_PM%Partrans(j)=0
            else
              write(retmsg,60) j        
60            format('Error: Parameter transformation type not allowed.' &
              ' Check row ',i6,' of the prior mean data table. Excecution stopped.')
              call utl_writmess(6,retmsg)
              stop
            endif
            if (columnstring(3).ne.' ') call drealread(ifail, columnstring(3), d_PM%alpha(j))
         enddo
       !********************************************************************************
       !********************************************************************************
         
   else
        
        select case (cv_PM%Qbb_form)
          case (0)    ! no prior covariance provided for Qbb
             write(retmsg,61)         
61           format('Error: Cannot provide prior beta values without covariance.' &
                  ' In cv_prior_mean block, if betas_flag = 1, beta_cov_form cannot be 0. Execution stopped.  ')
             call utl_writmess(6,retmsg)
             stop        
            
          case (1)    !covariance of beta on diagonal
            ! allocate and initialize both d_PM and column variables for table reading
            call bdp_alloc_d_PM(d_PM,columnname,columnstring,5,BL(3)%numrows)
            do i=1,BL(3)%numrows
              call ids_read_block_table(ifail,BL(3)%label,5,columnname,columnstring,line,filename)
              
               !**********************************************************************************
               !*************     PARAMETER TRANSFORMATION FLAG     ******************************
                tmp_str = trim(adjustl(columnstring(2)))  
                call bpc_casetrans(tmp_str,'lo')
                  if (tmp_str.eq.'log') then
                    d_PM%Partrans(i)=1
                  elseif (tmp_str.eq.'power') then
                    d_PM%Partrans(i)=2
                  elseif (tmp_str.eq.'none') then
                    d_PM%Partrans(i)=0
                  else
                    write(retmsg,600) i        
600                  format('Error: Parameter transformation type not allowed.' &
                    ' Check row ',i6,' of the prior mean data table. Excecution stopped.')
                    call utl_writmess(6,retmsg)
                    stop
                 endif
              !********************************************************************************
              !********************************************************************************
              if (columnstring(3).ne.' ') call drealread(ifail, columnstring(3), d_PM%alpha(i))
              call drealread(ifail, columnstring(4), d_PM%beta_0(i))
              call drealread(ifail, columnstring(5), d_PM%Qbb(i,i))
            enddo
          case (2)    !full beta covariance matrix.  cv_PM%prior_betas by cv_PM%prior_betas
            numcol = cv_PAR%p + 4
            ! allocate and initialize both d_PM and column variables for table reading
            call  bdp_alloc_d_PM(d_PM,columnname,columnstring,numcol,BL(3)%numrows)
            do i=1,BL(3)%numrows
              call ids_read_block_table(ifail,BL(3)%label,numcol,columnname,columnstring,line,filename)
               
               !**********************************************************************************
               !*************     PARAMETER TRANSFORMATION FLAG     ******************************
                tmp_str = trim(adjustl(columnstring(2)))  
                call bpc_casetrans(tmp_str,'lo')
                  if (tmp_str.eq.'log') then
                    d_PM%Partrans(i)=1
                  elseif (tmp_str.eq.'power') then
                    d_PM%Partrans(i)=2
                  elseif (tmp_str.eq.'none') then
                    d_PM%Partrans(i)=0
                  else
                    write(retmsg,66) i        
66                  format('Error: Parameter transformation type not allowed.' &
                    ' Check row ',i6,' of the prior mean data table. Excecution stopped.')
                    call utl_writmess(6,retmsg)
                    stop
                 endif
              !********************************************************************************
              !********************************************************************************
              if (columnstring(3).ne.' ') call drealread(ifail, columnstring(3), d_PM%alpha(i))
              call drealread(ifail, columnstring(4), d_PM%beta_0(i))
              do j = 5,numcol
                call drealread(ifail, columnstring(j), d_PM%Qbb(i,j-4))
              enddo
            enddo
        end select
          if (associated(columnname))           deallocate(columnname)
          if (associated(columnstring))         deallocate(columnstring)
          
    endif  !(cv_PM%betas_flag.ne.0)
          
end subroutine bdp_read_data_prior_mean

    !********  subroutine bdp_read_cv_parameters(BL,cv_PAR,inunit)
           subroutine bdp_read_cv_parameters(BL,cv_PAR,inunit,retmsg)
           ! SUBROUTINE to read in parameter control keywords
             use bayes_pest_control
             use jupiter_input_data_support
           ! DECLARATIONS
             implicit none
             type(cv_param),       intent(inout)  :: cv_PAR
             type(tp_block),       intent(inout)  :: BL(NUM_BLOCK) 
             integer,              intent(in)     :: inunit
             character (len=ERRORWIDTH), intent(inout)  :: retmsg
             integer                              :: i  !local counter
             integer                              :: ifail
             
             if (BL(8)%numrows .EQ. 0) then
              retmsg = 'No input information provided in cv_parameters block'  
              return     
             endif 
              call bdp_init_keyword_vars(BL,8)
              call ids_read_block_keywords(ifail,BL(8)%label,BL(8)%numkw, &
               & BL(8)%keywords,BL(8)%keywordstring,BL(8)%keywordline)
              if (ifail.ne.0) then
               call ids_get_message_string(retmsg)
               call utl_writmess(6,retmsg)
              endif
              
           ! Parse out the values as appropriate
           do i=1,BL(8)%numkw
              if(BL(8)%keywordline(i) .ne. 0) then
               select case (trim(BL(8)%keywords(i)))
                case ('ndim')
                  call intread(ifail,BL(8)%keywordstring(i), cv_PAR%ndim)
                case default
                !!!MNF DEBUG MAKE AN ERROR EXCEPTION HERE!!!!!
               end select
              endif
             enddo
                
    end subroutine bdp_read_cv_parameters

    !********  subroutine bdp_read_cv_compression(BL,Q0_All,nbeta,inunit,retmsg)
           subroutine bdp_read_cv_compression(BL,Q0_All,nbeta,inunit,retmsg)
           ! SUBROUTINE to read in Compression block information
             use bayes_pest_control
             use jupiter_input_data_support
           ! DECLARATIONS
             implicit none
             type(Q0_compr),       intent(inout)  :: Q0_All(:)
             type(tp_block),       intent(inout)  :: BL(NUM_BLOCK) 
             integer,              intent(in)     :: inunit
             character (len=ERRORWIDTH), intent(inout)  :: retmsg
             integer                              :: i  !local counter
             integer                              :: ifail, line
             character (len=COLWIDTH), pointer    :: columnname(:)
             character (len=COLWIDTH), pointer    :: columnstring(:)         
             character (len=200)                  :: filename
             integer                              :: numcol
             integer,              intent(in)     :: nbeta                    
            
             if (bl(9)%numrows .EQ. 0) then
              retmsg = 'No input information provided in Q_compression_cv block. Q_compression_flag must be 0. Execution stopped.'
              call utl_writmess(6,retmsg)
              stop   
             endif 
              call bdp_init_keyword_vars(BL,9)
              if (nbeta.ne.bl(9)%numrows) then
                write(retmsg,10) 
10               format('Error: The number of beta associations in Q_compression_cv block differs from', &
                   ' the number of prior betas in prior_mean_data block. Excecution stopped.')
                call utl_writmess(6,retmsg)
                stop
              endif
              
        
              numcol = 5
              allocate(columnstring(numcol))
              allocate(columnname(numcol))
              columnname=(/'BetaAssoc','Toep_flag','Nrow','Ncol','Nlay'/)
              columnstring=(/' ',' '/)
              Q0_All%Nrow = UNINIT_INT
              Q0_All%Ncol = UNINIT_INT
              Q0_All%Nlay = UNINIT_INT
              do i=1,nbeta
                 call ids_read_block_table(ifail,bl(9)%label,numcol,columnname,columnstring,line,filename)
                 call intread(ifail, columnstring(1),Q0_All(i)%BetaAss)
                 call intread(ifail, columnstring(2),Q0_All(i)%Toep_flag)
                 if (Q0_All(i)%Toep_flag.ne.0) then !Assign Nrow Ncol Nlay just if toeplitz is required
                  call intread(ifail, columnstring(3),Q0_All(i)%Nrow)
                  call intread(ifail, columnstring(4),Q0_All(i)%Ncol)
                  call intread(ifail, columnstring(5),Q0_All(i)%Nlay)
                 endif
              enddo
              
    end subroutine bdp_read_cv_compression


    !********  subroutine bdp_read_parameter_groups(BL,cv_PAR,inunit,errmsg)
           subroutine bdp_read_parameter_groups(BL,cv_PAR,inunit,retmsg)
           ! SUBROUTINE to read in parameter group information
             use bayes_pest_control
             use jupiter_input_data_support
           ! DECLARATIONS
             implicit none
             type(cv_param),       intent(inout)  :: cv_PAR
             type(tp_block),       intent(inout)  :: BL(NUM_BLOCK) 
             integer,              intent(in)     :: inunit
             character (len=ERRORWIDTH), intent(inout)  :: retmsg
             integer                              :: i  !local counter
             integer                              :: ifail, line
             character (len=COLWIDTH), pointer    :: columnname(:)
             character (len=COLWIDTH), pointer    :: columnstring(:)         
             character (len=200)                  :: filename
             integer                              :: numcol
             
            
             if (bl(10)%numrows .EQ. 0) then
              retmsg = 'No input information provided in cv_parameters block'  
              return     
             endif 
              call bdp_init_keyword_vars(BL,10)
              cv_PAR%npargp=bl(10)%numrows
              allocate(cv_PAR%grp_name(cv_PAR%npargp))
              allocate(cv_PAR%grp_type(cv_PAR%npargp))     
              allocate(cv_PAR%derinc(cv_PAR%npargp))     
              numcol = 3
              allocate(columnstring(numcol))
              allocate(columnname(numcol))
              columnname=(/'groupname','grouptype','derinc'/)
              columnstring=(/' ',' ',' '/)
              
              do i=1,cv_PAR%npargp
                 call ids_read_block_table(ifail,bl(10)%label,3,columnname,columnstring,line,filename)
                 cv_PAR%grp_name(i) = trim(adjustl(columnstring(1)))
                 call intread(ifail, columnstring(2),cv_PAR%grp_type(i))
                 call DREALREAD(ifail,columnstring(3),cv_PAR%derinc(i))
              enddo
              
    end subroutine bdp_read_parameter_groups

    !********  subroutine bdp_read_data_anisotropy(BL,d_ANI,cv_PAR,inunit,retmsg)
           subroutine bdp_read_data_par_anisotropy(BL,d_ANI,cv_PAR,inunit,retmsg)
           ! SUBROUTINE TO READ ANISOTROPY INFORMATION FOR PARAMETERS
             use bayes_pest_control
             use utilities
             
           ! DECLARATIONS
             implicit none
             type (d_anisotropy)                       :: d_ANI
             type (cv_param), intent(in)               :: cv_PAR
             type(tp_block)                            :: BL(NUM_BLOCK)
             integer, intent(in)                       :: inunit
             integer                                   :: ifail, line
             character (len=ERRORWIDTH), intent(inout) :: retmsg
             character (len=COLWIDTH), pointer         :: columnname(:)
             character (len=COLWIDTH), pointer         :: columnstring(:)         
             character (len=200)                       :: filename
             integer                                   :: i, numcol
             
             if (bl(17)%numrows .EQ. 0) then
              retmsg = 'WARNING: No input information provided in data_parameter_anisotropy block'  
              call utl_writmess(6,retmsg)
              return     
             endif
             
             numcol = 4
             allocate(columnstring(numcol))
             allocate(columnname(numcol))
             allocate(d_ANI%BetaAssoc(cv_PAR%p))
             allocate(d_ANI%horiz_angle(cv_PAR%p))
             allocate(d_ANI%horiz_ratio(cv_PAR%p))
             allocate(d_ANI%vertical_ratio(cv_PAR%p))
             
             columnname=(/'BetaAssoc','horiz_angle','horiz_ratio','vertical_ratio'/)
             do i = 1,cv_PAR%p
                 call ids_read_block_table(ifail,BL(17)%label,numcol,columnname,columnstring,line,filename)
                 call intread(ifail, columnstring(1),d_ANI%BetaAssoc(i))
                 call drealread(ifail, columnstring(2),d_ANI%horiz_angle(i))
                 call drealread(ifail, columnstring(3),d_ANI%horiz_ratio(i))
                 call drealread(ifail, columnstring(4),d_ANI%vertical_ratio(i))               
             enddo
           end subroutine bdp_read_data_par_anisotropy
    
    
    !********  subroutine bdp_read_data_parameters(BL,d_PAR,npargp,inunit,retmsg)
           subroutine bdp_read_data_parameters(errstruc,BL,d_PAR,cv_PAR,cv_A,Q0_All,inunit,retmsg,miostruc)
             use bayes_pest_control
             use jupiter_input_data_support
             use error_message
             use utilities
             use model_input_output
              
           ! DECLARATIONS
             implicit none
             type (mio_struc)           :: miostruc
             type (err_failure_struc)   :: errstruc
             character (len=50)                        :: tmp_nm,tmp_str
             type(d_param),        intent(inout)       :: d_PAR
             type(cv_param),        intent(inout)      :: cv_PAR
             type(tp_block),       intent(inout)       :: BL(NUM_BLOCK)
             type (Q0_compr)                           :: Q0_All(:) 
             type (cv_algorithmic)                     :: cv_A
             integer,              intent(in)          :: inunit
             character (len=ERRORWIDTH), intent(inout) :: retmsg
             character (len=COLWIDTH), pointer     :: columnname(:)
             character (len=COLWIDTH)                  :: c1, c2
             character (len=COLWIDTH), pointer     :: columnstring(:)         
             character (len=200)                       :: filename
             integer                                   :: numcol
             integer                                   :: i,j,k,ii  !local counters
             integer                                   :: ifail, line, n1
             
             
            if (bl(11)%numrows .EQ. 0) then
              retmsg = 'No input information provided in data_parameters block'  
              return     
            endif 
            
       ! ALLOCATIONS, READING, and PARSING  
         ! parameter data  and associated allocations
            cv_PAR%npar = bl(11)%numrows
            numcol = 5 + cv_PAR%ndim
            if(mio_initialise_parameters(errstruc,miostruc,cv_PAR%npar).ne.0) then
              call utl_bomb_out(errstruc)
              n1=mio_finalise(errstruc,miostruc)
              stop 
            endif
            allocate(columnname(numcol))
            allocate(columnstring(numcol))
            allocate(d_PAR%lox(cv_PAR%npar,cv_PAR%ndim))
            allocate(d_PAR%pars(cv_PAR%npar))
            allocate(d_PAR%parnme(cv_PAR%npar))
            allocate(d_PAR%pars_old(cv_PAR%npar))
            allocate(d_PAR%SenMethod(cv_PAR%npar))
            allocate(d_PAR%BetaAssoc(cv_PAR%npar))
            allocate(d_PAR%group(cv_PAR%npar))
            allocate(d_Par%Group_type(cv_PAR%npar))
            if (cv_A%lns_flag.eq.1) then
               allocate(d_PAR%pars_lns(cv_PAR%npar))
               d_PAR%pars_lns=UNINIT_REAL
            endif
         ! initializations

                 columnname(1:5)=(/'ParamName', 'StartValue', 'GroupName', 'BetaAssoc', 'SenMethod'/)
             columnstring(1:5)=' ' ! array
             do j=1,cv_PAR%ndim
               call utl_int2char(j,tmp_str)
               columnname(j+5)='x'//tmp_str
               columnstring(j+5)=' '
             enddo
             d_PAR%lox        = UNINIT_REAL
             d_PAR%SenMethod  = UNINIT_INT  ! array
             d_PAR%BetaAssoc  = UNINIT_INT  ! array
             d_Par%Group_type = UNINIT_INT  ! array (Associate a group type based on the group name)
             d_PAR%pars       = UNINIT_REAL ! array
             d_PAR%parnme       = UNINIT_CHAR ! array
             d_PAR%pars_old   = UNINIT_REAL ! array
             d_PAR%group      = UNINIT_CHAR ! array 
          if (cv_A%Q_compression_flag.ne.0) then  !This means that we are using the compressed form of Q0; 
                                                  !we need to know the number of parameters for each beta 
                                                  !and where each different beta association start in the list 
             Q0_All%npar = 0    !Initialize the counter of the parameters for each beta     
          ! read in the data
              do j=1,cv_PAR%npar
                call ids_read_block_table(ifail,bl(11)%label,numcol,columnname,columnstring,line,filename)
                  tmp_nm = trim(adjustl(columnstring(1)))
                call bpc_casetrans(tmp_nm,'lo')
                  if(mio_put_parameter(errstruc,miostruc,j,tmp_nm).ne.0) then
                    call utl_bomb_out(errstruc)
                    n1=mio_finalise(errstruc,miostruc)
                    stop 
                  endif
                d_PAR%parnme(j) = trim(adjustl(tmp_nm))
                call drealread(ifail, columnstring(2),d_PAR%pars(j))
                d_PAR%group(j) = trim(adjustl(columnstring(3)))
                call intread(ifail, columnstring(4),d_PAR%BetaAssoc(j))
                  if (d_PAR%BetaAssoc(j).gt.cv_PAR%p) then             !*** This control avoid that a parameter
                    write(retmsg,50) d_PAR%BetaAssoc(j), cv_PAR%p, j   !*** could be associate to an undefined beta      
50                  format('Error: Beta association value ',i6, ' exceeds the number of beta values', i6, &
                    '. Check rows ',i6,' of the parameter table. Excecution stopped.')
                    call utl_writmess(6,retmsg)
                    stop
                  endif
	           
	            Q0_All(d_PAR%BetaAssoc(j))%npar = Q0_All(d_PAR%BetaAssoc(j))%npar + 1    !Counter of number of parameters
	            if (Q0_All(d_PAR%BetaAssoc(j))%npar.eq.1) Q0_All(d_PAR%BetaAssoc(j))%Beta_Start = j !Starting line of parameters with different beta
	           
	            call intread(ifail, columnstring(5),d_PAR%SenMethod(j))
                do ii = 6,numcol
                  call drealread(ifail, columnstring(ii), d_PAR%lox(j,ii-5))
                enddo             
              enddo
              
          else
           ! read in the data
              do j=1,cv_PAR%npar
                call ids_read_block_table(ifail,bl(11)%label,numcol,columnname,columnstring,line,filename)
                  tmp_nm = trim(adjustl(columnstring(1)))
                call bpc_casetrans(tmp_nm,'lo')
                  if(mio_put_parameter(errstruc,miostruc,j,tmp_nm).ne.0) then
                    call utl_bomb_out(errstruc)
                    n1=mio_finalise(errstruc,miostruc)
                    stop 
                  endif
                d_PAR%parnme(j) = trim(adjustl(tmp_nm))
                call drealread(ifail, columnstring(2),d_PAR%pars(j))
                d_PAR%group(j) = trim(adjustl(columnstring(3)))
                call intread(ifail, columnstring(4),d_PAR%BetaAssoc(j))
                  if (d_PAR%BetaAssoc(j).gt.cv_PAR%p) then             !*** This control avoid that a parameter
                    write(retmsg,500) d_PAR%BetaAssoc(j), cv_PAR%p, j   !*** could be associate to an undefined beta      
500                  format('Error: Beta association value ',i6, ' exceeds the number of beta values', i6, &
                    '. Check rows ',i6,' of the parameter table. Excecution stopped.')
                    call utl_writmess(6,retmsg)
                    stop
                  endif
	            call intread(ifail, columnstring(5),d_PAR%SenMethod(j))
                do ii = 6,numcol
                  call drealread(ifail, columnstring(ii), d_PAR%lox(j,ii-5))
                enddo             
              enddo
              
          endif   ! if compr
              
              
              do i=1, cv_PAR%npargp !Associate a group type to each parameter based on the group name
                where (d_PAR%group.eq.cv_PAR%grp_name(i))
                 d_PAR%Group_type = cv_PAR%grp_type(i)
                end where
              enddo 
                              
          if (associated(columnname))           deallocate(columnname)
          if (associated(columnstring))         deallocate(columnstring)             
    end subroutine bdp_read_data_parameters  

!********  subroutine bdp_read_cv_structural_parameters_tbl(BL,cv_S,inunit)
       subroutine bdp_read_cv_structural_parameters_tbl(BL,cv_S,nrows,inunit,retmsg)
       ! SUBROUTING to read in algorithmic control keywords
         use bayes_pest_control
         use jupiter_input_data_support
         use bpi_initializations
       ! DECLARATIONS
         implicit none
         type (err_failure_struc)                   :: errstruc
         type(cv_struct),            intent(inout)  :: cv_S
         type(tp_block),             intent(inout)  :: BL(NUM_BLOCK) 
         integer,                    intent(in)     :: inunit,nrows
         character (len=ERRORWIDTH), intent(inout)  :: retmsg
         character (len=COLWIDTH),   pointer        :: columnname(:)
         character (len=COLWIDTH),   pointer        :: columnstring(:)
         character (len=200)                        :: filename, amessage,function_name         
         integer                                    :: k,ifail,line,numcol
       
    !MD ALLOCATE AND INITIALIZE VARIABLES
         numcol = 6  
         call bdp_alloc_init_cv_S(cv_S,columnname,columnstring,numcol,nrows)
         if (BL(4)%numrows .EQ. 0) then
          retmsg = 'All default values accepted for cv_structural_parameters block' 
          return      
         endif
      ! quick error failure for number of rows based on the number of expected rows
         if (nrows .ne. BL(4)%numrows) then
            amessage = 'Number of rows for structural parameters inconsistent with the number expected based on number of beta associations'
            function_name = 'bdp_read_cv_structural_parameters_tbl'
            call err_add_error(errstruc,amessage,function_name)
            call utl_bomb_out(errstruc)
         endif 
         ! read and parse
        do k=1,nrows
          call ids_read_block_table(ifail,BL(4)%label,numcol,columnname,columnstring,line,filename)
             if (columnstring(2).ne.' ') call intread(ifail, columnstring(2),cv_S%prior_cov_mode(k))
             if (columnstring(3).ne.' ') call intread(ifail, columnstring(3),cv_S%var_type(k)) 
             if (columnstring(4).ne.' ') call intread(ifail, columnstring(4),cv_S%struct_par_opt(k))
             select case (cv_S%var_type(k))
                case(0) ! nugget variogram
                    cv_S%num_theta_type(k) = 1
                case(1) ! linear variogram
                    cv_S%num_theta_type(k) = 1
                case(2) ! exponential variogram
                    cv_S%num_theta_type(k) = 2
             end select
             if (columnstring(5).ne.' ') call intread(ifail, columnstring(5),cv_S%trans_theta(k))
             if (columnstring(6).ne.' ') call drealread(ifail, columnstring(6),cv_S%alpha_trans(k))
        enddo
       
 end subroutine bdp_read_cv_structural_parameters_tbl

!********  subroutine bdp_read_data_structural_parameters      (BL,cv_S,d_S,numzones,inunit,errmsg)
! Read the starting values of theta variables
       subroutine bdp_read_data_structural_parameters      (BL,cv_S,d_S,cv_A,numzones,inunit,retmsg)
          use bayes_pest_control
          use jupiter_input_data_support
       ! DECLARATIONS
         implicit none
         type(d_struct),             intent(inout)  :: d_S
         type(cv_struct),            intent(inout)  :: cv_S
         type(cv_algorithmic),       intent(inout)  :: cv_A
         integer,                    intent(in)     :: numzones
         type(tp_block),             intent(inout)  :: BL(NUM_BLOCK) 
         integer,                    intent(in)     :: inunit
         character (len=ERRORWIDTH), intent(inout)  :: retmsg
         character (len=COLWIDTH),   pointer        :: columnname(:)
         character (len=COLWIDTH),   pointer        :: columnstring(:)         
         character (len=200)                        :: filename
         integer                                    :: numcol
         integer                                    :: i,j  !local counter
         integer                                    :: ifail, line
                  
        if (BL(5)%numrows .EQ. 0) then
          retmsg = 'No input information provided in structural_parameters_data block'  
          return     
        endif 
   ! ALLOCATIONS, READING, and PARSING
   ! theta_0 STRUCTURAL_PARAMETERS_DATA block

      numcol = maxval(cv_S%num_theta_type) + 1 !MD numcol is the maximum number of structural parameters 
                                               !defined in num_theta_type + 1 (BetaAssoc)
      ! allocate
        call bdp_alloc_d_S(d_S, cv_A, columnname, columnstring, numcol, numzones) 
      ! read and parse
           do i=1,numzones
            call ids_read_block_table(ifail,BL(5)%label,numcol,columnname,columnstring,line,filename)
             do j=2,numcol
              call drealread(ifail,columnstring(j), d_S%theta_0(i,j-1)) !1 is slope, 2 is L
               if (d_S%theta_0(i,j-1).lt.0.) then
                  d_S%theta_0(i,j-1)=UNINIT_REAL
                  if (cv_S%num_theta_type(i).ge.j-1) then
                   write(retmsg,10) i , j
10                 format('Error: Found unexpected negative value in row',i6,' column ',i6,' of the structural paramaters data table,', &
                   ' excecution stopped.')
                   call utl_writmess(6,retmsg)
                   stop
                  endif
               else   
                  if (cv_S%num_theta_type(i).lt.j-1) then
                   d_S%theta_0(i,j-1)=UNINIT_REAL
                   write(retmsg,20) i , j
20                 format('Warning: Found unexpected positive value in row',i6,' column ',i6,' of the structural paramaters data table,', &
                   ' the value will be ignored.')
                   call utl_writmess(6,retmsg)
                   endif
               endif
             enddo
           enddo  
           d_S%theta(:,:) = d_S%theta_0  ! start out with theta_0 as the first values of theta
          if (associated(columnname))           deallocate(columnname)
          if (associated(columnstring))         deallocate(columnstring)
          end subroutine bdp_read_data_structural_parameters  
          
!********  subroutine bdp_read_structural_parameters_cov    (BL,cv_S,d_S,theta_cov_form,inunit,retmsg)
! Read the prior covariance of the structural parameter
       subroutine bdp_read_structural_parameters_cov      (BL,cv_S,d_S,theta_cov_form,inunit,retmsg)
          use bayes_pest_control
          use jupiter_input_data_support
       ! DECLARATIONS
         implicit none
         type (err_failure_struc)                  :: errstruc
         type(d_struct),       intent(inout)       :: d_S
         type(cv_struct),      intent(inout)       :: cv_S
         integer,              intent(in)          :: theta_cov_form
         type(tp_block),       intent(inout)       :: BL(NUM_BLOCK) 
         integer,              intent(in)          :: inunit
         character (len=ERRORWIDTH), intent(inout) :: retmsg
         character (len=COLWIDTH), pointer     :: columnname(:)
         character (len=COLWIDTH), pointer     :: columnstring(:)         
         character (len=200)                       :: filename, amessage,function_name
         integer                                   :: numcol, numrow
         integer                                   :: i,j  !local counter
         integer                                   :: ifail, line
         
        function_name = 'bdp_read_structural_parameters_cov'
        if (BL(6)%numrows .EQ. 0) then
          retmsg = 'No input information provided in prior_structural_parameters_data block'  
          return     
        endif 
   ! ALLOCATIONS, READING, and PARSING
   ! theta_0 and theta_cov from PRIOR_STRUCTURAL_PARAMETERS block
    
        
       th_cov_form : select case (theta_cov_form) !MD Form of theta covariance:  [0] none, [1] diag, [2] full matrix
         case(0)
           
           return
         case (1)
           numrow=sum(cv_S%num_theta_type)
           numcol=1
         case (2)
           numcol=sum(cv_S%num_theta_type)
           numrow=numcol    
        end select th_cov_form
      ! quick error failure for number of rows based on the number of expected rows
      if (numrow.ne.BL(6)%numrows) then
         amessage = 'Number of rows for theta covariance inconsistent with the number expected based on num_theta_type'
         call err_add_error(errstruc,amessage,function_name)
         call utl_bomb_out(errstruc)
      endif
      ! allocate
        call bdp_alloc_cov_S(d_S,cv_S,columnname,columnstring,numcol,numrow) 
      ! read and parse
           do i=1,numrow
            call ids_read_block_table(ifail,BL(6)%label,numcol,columnname,columnstring,line,filename)
             do j=1,numcol
              call drealread(ifail,columnstring(j), d_S%theta_cov(i,j))
             enddo
           enddo
     
          if (associated(columnname))           deallocate(columnname)
          if (associated(columnstring))         deallocate(columnstring)                     
end subroutine bdp_read_structural_parameters_cov  

!**********subroutine bdp_read_epistemic_error (BL,cv_S,d_S,theta_cov_form,inunit,retmsg)
   subroutine bdp_read_epistemic_error (BL,d_S, cv_A,cv_S,inunit,retmsg)
          use bayes_pest_control
          use jupiter_input_data_support
       ! DECLARATIONS
         implicit none
         type(d_struct),       intent(inout)       :: d_S
         type(tp_block),       intent(inout)       :: BL(NUM_BLOCK) 
         type(cv_algorithmic), intent(inout)       :: cv_A
         type (cv_struct),     intent(in)          :: cv_S
         integer,              intent(in)          :: inunit
         character (len=ERRORWIDTH), intent(inout) :: retmsg
         character (len=COLWIDTH), pointer         :: columnname(:)
         character (len=COLWIDTH), pointer         :: columnstring(:)         
         character (len=200)                       :: filename
         integer                                   :: i  !local counter
         integer                                   :: ifail, line
         
        
      ! READ IN EPISTEMIC ERROR TERM   
       
        if (BL(7)%numrows .EQ. 0) then
          retmsg = 'No input information provided in epistemic_error_term block'  
          return     
        endif 
          call bdp_init_keyword_vars(BL,7)
          call ids_read_block_keywords(ifail,BL(7)%label,BL(7)%numkw, &
           & BL(7)%keywords,BL(7)%keywordstring,BL(7)%keywordline)
          if (ifail.ne.0) then
           call ids_get_message_string(retmsg)
           call utl_writmess(6,retmsg)
          endif
         d_S%sig = UNINIT_INT  ! array
       ! Parse out the values as appropriate
         do i= 1,BL(7)%numkw
          if(BL(7)%keywordline(i) .ne. 0) then
           select case (trim(BL(7)%keywords(i)))
            case ('sig_0')
              call drealread(ifail,BL(7)%keywordstring(i), d_S%sig_0)
              d_S%sig = d_S%sig_0     ! start out with the initial value of sig_0 as the first value
            case ('sig_opt')
              call intread(ifail,BL(7)%keywordstring(i), d_S%sig_opt)
            case ('sig_p_var')
              call drealread(ifail,BL(7)%keywordstring(i), d_S%sig_p_var)
            case ('trans_sig')
              call intread(ifail,BL(7)%keywordstring(i), d_S%trans_sig)
            case ('alpha_trans')  
              call drealread(ifail,BL(7)%keywordstring(i), d_S%alpha_trans_sig)
           end select
          endif
         enddo
            
         if ((maxval(cv_S%struct_par_opt).eq.0).and.(d_S%sig_opt.eq.0).and.(cv_A%it_max_bga.gt.1)) then  !--------------------| This control avoid
           cv_A%it_max_bga =  1   !-------------------------------------------------------------------------------------------| that it_max_bga is      
           write(retmsg,10)    !----------------------------------------------------------------------------------------------| greater than 1 if 
10         format('Warning: No structural parameter optimization requested. it_max_bga in algorithmic_cv block must be 1.',& !| no structural 
                 & ' The entered value will be ignored.') !-------------------------------------------------------------------| parameter 
           call utl_writmess(6,retmsg) !--------------------------------------------------------------------------------------| optimization
         endif  !-------------------------------------------------------------------------------------------------------------| is required.
                                                                  
end subroutine bdp_read_epistemic_error



!********  subroutine bdp_read_data_model_command_line(BL,d_MOD,numcom,inunit,retmsg)
       subroutine bdp_read_data_model_command_line(BL,d_MOD,deriv_mode,inunit,retmsg)
          use bayes_pest_control
          use jupiter_input_data_support
       ! DECLARATIONS
         implicit none
         type(d_comlin),       intent(inout)       :: d_MOD
         type(tp_block),       intent(inout)       :: BL(NUM_BLOCK) 
         integer,              intent(in)          :: inunit
         integer,              intent(in)          :: deriv_mode
         character (len=ERRORWIDTH), intent(inout) :: retmsg      
         character (len=200)                       :: filename
         integer                                   :: numcol
         integer                                   :: i,j,k  !local counter
         integer                                   :: ifail, line
         
          call bdp_init_keyword_vars(BL,14)
          call ids_read_block_keywords(ifail,BL(14)%label,BL(14)%numkw, &
             & BL(14)%keywords,BL(14)%keywordstring,BL(14)%keywordline)
          if (ifail.ne.0) then
           call ids_get_message_string(retmsg)
           call utl_writmess(6,retmsg)
          endif
         ! Parse out the values as appropriate
         do i=1,BL(14)%numkw
          if(BL(14)%keywordline(i) .ne. 0) then
            select case (trim(BL(14)%keywords(i)))
             case ('Command')
               d_MOD%com = trim(adjustl(BL(14)%keywordstring(i)))
             case ('DerivCommand')
               d_MOD%dercom = trim(adjustl(BL(14)%keywordstring(i)))
            end select
           endif
          enddo
         
end subroutine bdp_read_data_model_command_line  
      
!********  subroutine bdp_read_data_model_input_output(BL,d_MIO,cv_MIO,cv_A,inunit,retmsg)
         subroutine bdp_read_data_model_input_output(BL,d_MIO,cv_MIO,cv_A,npargp,inunit,retmsg)
          use bayes_pest_control
          use jupiter_input_data_support
       ! DECLARATIONS
         implicit none
         type(d_minout),       intent(inout)       :: d_MIO
         type(cv_algorithmic), intent(in)          :: cv_A
         type (cv_minout),     intent(inout)       :: cv_MIO
         type(tp_block),       intent(inout)       :: BL(NUM_BLOCK) 
         integer,              intent(in)          :: inunit,npargp ! npargp is number of parameter groups 
         character (len=ERRORWIDTH), intent(inout) :: retmsg
         character (len=COLWIDTH), pointer     :: columnname(:)
         character (len=COLWIDTH), pointer     :: columnstring(:)         
         character (len=200)                       :: filename
         integer                                   :: numcol
         integer                                   :: i,j,k  !local counter
         integer                                   :: ifail, line
         
        nullify(d_MIO%tpl)
        nullify(d_MIO%infle)
        nullify(d_MIO%ins)
        nullify(d_MIO%outfle) 
        nullify(d_MIO%pargroup) 
        
         
        if (bl(15)%numrows .EQ. 0) then
          retmsg = 'No input information provided in model_input_files block'  
          return     
        endif 
        
   ! ALLOCATIONS, READING, and PARSING
        cv_MIO%ntplfle = bl(15)%numrows
        allocate (d_MIO%tpl(cv_MIO%ntplfle))
        allocate (d_MIO%infle(cv_MIO%ntplfle))
        allocate (d_MIO%pargroup(npargp))
        d_MIO%tpl       = UNINIT_CHAR ! array
        d_MIO%infle     = UNINIT_CHAR ! array
        d_MIO%pargroup  = UNINIT_CHAR ! array      

        if (cv_A%deriv_mode .eq. 4) then ! parameter groups used - (cv_A%deriv_mode .eq. 4)
        allocate (columnname(3))
        allocate (columnstring(3))
      ! first read input/tpl information
         columnname=(/'TemplateFile','ModInFile','groupname'/)
         columnstring=' ' ! array
         do i=1,cv_MIO%ntplfle
          call ids_read_block_table(ifail,bl(15)%label,3,columnname,columnstring,line,filename)
          d_MIO%tpl(i) = trim(adjustl(columnstring(1)))
          d_MIO%infle(i) = trim(adjustl(columnstring(2)))
          d_MIO%pargroup(i) = trim(adjustl(columnstring(3)))
         enddo   
         
        else ! no parameter groups used -  (cv_A%deriv_mode .ne. 4)
        allocate (columnname(2))
        allocate (columnstring(2))
      ! first read input/tpl information
         columnname=(/'TemplateFile','ModInFile'/)
         columnstring=' ' ! array
         do i=1,cv_MIO%ntplfle
          call ids_read_block_table(ifail,bl(15)%label,2,columnname,columnstring,line,filename)
          d_MIO%tpl(i) = trim(adjustl(columnstring(1)))
          d_MIO%infle(i) = trim(adjustl(columnstring(2)))
         enddo
        endif

        ! next read output/ins information
      if (bl(16)%numrows .EQ. 0) then
          retmsg = 'No input information provided in model_output_files block'  
          return     
        endif 
        cv_MIO%ninsfle = bl(16)%numrows
        allocate (d_MIO%ins(cv_MIO%ninsfle))
        allocate (d_MIO%outfle(cv_MIO%ninsfle))
        
         columnname=(/'InstructionFile','ModOutFile'/)
         columnstring=' '   ! array
         do i=1,cv_MIO%ninsfle
          call ids_read_block_table(ifail,bl(16)%label,2,columnname,columnstring,line,filename)
          d_MIO%ins(i) = trim(adjustl(columnstring(1)))
          d_MIO%outfle(i) = trim(adjustl(columnstring(2)))
         enddo
          if (associated(columnname))           deallocate(columnname)
          if (associated(columnstring))         deallocate(columnstring)
         
end subroutine bdp_read_data_model_input_output  

!********  subroutine bdp_read_observation_groups(BL,cv_OBS,inunit,errmsg)
       subroutine bdp_read_observation_groups(BL,cv_OBS,inunit,retmsg)
       ! SUBROUTINE to read in parameter group information
         use bayes_pest_control
         use jupiter_input_data_support
       ! DECLARATIONS
         implicit none
         type(cv_observ),       intent(inout) :: cv_OBS
         type(tp_block),       intent(inout)  :: BL(NUM_BLOCK) 
         integer,              intent(in)     :: inunit
         character (len=ERRORWIDTH), intent(inout)  :: retmsg
         integer                              :: i  !local counter
         integer                              :: ifail, line
         character (len=COLWIDTH), pointer    :: columnname(:)
         character (len=COLWIDTH), pointer    :: columnstring(:)         
         character (len=200)                  :: filename
         
         if (bl(12)%numrows .EQ. 0) then
          retmsg = 'No input information provided in observation_groups block'  
          return     
         endif 
          call bdp_init_keyword_vars(BL,12)
          cv_OBS%nobsgp=bl(12)%numrows
          allocate(cv_OBS%grp_name(cv_OBS%nobsgp))
          allocate(columnstring(cv_OBS%nobsgp))
          allocate(columnname(cv_OBS%nobsgp))
          columnname='groupname'
          columnstring=' '
          
          do i=1,cv_OBS%nobsgp
             call ids_read_block_table(ifail,bl(12)%label,1,columnname,columnstring,line,filename)
             cv_OBS%grp_name(i) = trim(adjustl(columnstring(1)))
          enddo
          
end subroutine bdp_read_observation_groups

!********  subroutine bdp_read_data_observations(BL,d_OBS,nobsgp,nobs,inunit,retmsg,miostruc)
       subroutine bdp_read_data_observations(BL,d_OBS, cv_OBS,inunit,retmsg,miostruc)
          use bayes_pest_control
          use jupiter_input_data_support
          use utilities
          use model_input_output
          use error_message
         
       ! DECLARATIONS
         implicit none
         type (mio_struc)                          :: miostruc
         type (err_failure_struc)                  :: errstruc
         character (len=50)                        :: tmp_nm
         type(d_observ),       intent(inout)       :: d_OBS
         type(cv_observ),       intent(inout)      :: cv_OBS         
         type(tp_block),       intent(inout)       :: BL(NUM_BLOCK) 
         integer,              intent(in)          :: inunit
         character (len=ERRORWIDTH), intent(inout) :: retmsg
         character (len=COLWIDTH), pointer         :: columnname(:)
         character (len=COLWIDTH), pointer         :: columnstring(:)         
         character (len=200)                       :: filename
         integer                                   :: numcol
         integer                                   :: i,j,k   !local counter
         integer                                   :: ifail, line, n1
 
         
        if (bl(13)%numrows .EQ. 0) then
          retmsg = 'No input information provided in data_observations block'  
          return     
        endif 
        
       
   ! ALLOCATIONS, READING, and PARSING
     ! observation data and associated allocations
        CV_obs%nobs = bl(13)%NUMROWS
        if(mio_initialise_observations(errstruc,miostruc,cv_OBS%nobs).ne.0) then
          call utl_bomb_out(errstruc)
          n1=mio_finalise(errstruc,miostruc)
          stop 
        endif
        k = 1
        numcol=4
        allocate(d_OBS%obs(cv_OBS%nobs))
        allocate(d_OBS%obsnme(cv_OBS%nobs))
        allocate(d_OBS%group(cv_OBS%nobs))
        allocate(d_OBS%weight(cv_OBS%nobs))
        allocate (columnname(numcol))
        allocate (columnstring(numcol))
        d_OBS%obs    = UNINIT_REAL ! array
        d_OBS%obsnme = UNINIT_CHAR ! array
        d_OBS%group  = UNINIT_CHAR ! array
        d_OBS%weight = UNINIT_REAL ! array
        
         columnname=(/'ObsName',  'ObsValue',  'GroupName',  'Weight'/)
         columnstring=' '  ! array
   
          do j=1,cv_OBS%nobs
           call ids_read_block_table(ifail,bl(13)%label,numcol,columnname,columnstring,line,filename)
           tmp_nm = trim(adjustl(columnstring(1)))
           call bpc_casetrans(tmp_nm,'lo')
           if(mio_put_observation(errstruc,miostruc,k,tmp_nm).ne.0) then
            call utl_bomb_out(errstruc)
            n1=mio_finalise(errstruc,miostruc)
            stop 
           endif
           d_OBS%obsnme(k) = tmp_nm
           call drealread(ifail, columnstring(2),d_OBS%obs(k))
           d_OBS%group(k)     =  trim(adjustl(columnstring(3)))
           call drealread(ifail, columnstring(4),d_OBS%weight(k) )
           
           k = k + 1
          enddo
          if (associated(columnname))           deallocate(columnname)
          if (associated(columnstring))         deallocate(columnstring)
end subroutine bdp_read_data_observations  

!********  subroutine bdp_alloc_d_S
       subroutine bdp_alloc_d_S(d_S, cv_A, columnname, columnstring, numcol, numzones)
         use bayes_pest_control
         use jupiter_input_data_support
       ! DECLARATIONS
         implicit none
         type(d_struct),                    intent(inout)   :: d_S
         type(cv_algorithmic),              intent(in)      :: cv_A
         character (len=COLWIDTH), pointer, intent(inout)   :: columnname(:)
         character (len=COLWIDTH), pointer, intent(inout)   :: columnstring(:) 
         character (len=3)                                  :: instr        
         integer,                           intent(in)      :: numcol,numzones
         integer                                            :: i !local counters
         
     ! *** ALLOCATIONS FOR THETA_0 ***
            allocate(columnname(numcol))
            allocate(columnstring(numcol))
            allocate(d_S%theta_0(numzones,numcol-1))
            allocate(d_S%theta(numzones,numcol-1))
            d_S%theta_0= UNINIT_REAL  !2-d array
            d_S%theta  = UNINIT_REAL  !3-d array
            columnname(1)='BetaAssoc'
            columnstring(1)=' '
            do i=2,numcol
             call int2char(i-1,instr)
             write(columnname(i),*) 'theta_0_',instr
             columnstring(i)=' '           
            enddo
 
       end subroutine bdp_alloc_d_S          

subroutine bdp_alloc_cov_S(d_S,cv_S,columnname,columnstring,numcol,numrow)
         use bayes_pest_control
         use jupiter_input_data_support
       ! DECLARATIONS
         implicit none
         type(d_struct),  intent(inout)                       :: d_S
         type(cv_struct),  intent(in)                         :: cv_S
         character (len=COLWIDTH), pointer, intent(inout) :: columnname(:)
         character (len=COLWIDTH), pointer, intent(inout) :: columnstring(:) 
         character (len=3)                                    :: instr        
         integer,                               intent(in)    :: numcol,numrow
         integer                                              :: i !local counters
         
     ! *** ALLOCATIONS FOR THETA_COV ***
            allocate(columnname(numcol))
            allocate(columnstring(numcol))
            allocate(d_S%theta_cov(numrow,numcol))
            d_S%theta_cov=UNINIT_REAL  !2-d array
             do i= 1,numcol
              call int2char(i,instr)
              write(columnname(i),*) 'theta_cov_',instr
              columnstring(i)=' '
             enddo
     
           end subroutine bdp_alloc_cov_S

!********  subroutine bdp_alloc_init_cv_S
       subroutine bdp_alloc_init_cv_S(cv_S,columnname,columnstring,numcol,numzones)
         use bayes_pest_control
         use jupiter_input_data_support
       ! DECLARATIONS
         implicit none
         type(cv_struct),                   intent(inout)     :: cv_S
         character (len=COLWIDTH), pointer, intent(inout)     :: columnname(:)
         character (len=COLWIDTH), pointer, intent(inout)     :: columnstring(:) 
         character (len=3)                                    :: instr        
         integer,                           intent(in)        :: numcol,numzones
         integer                                              :: i !local counters
         
     ! *** ALLOCATIONS FOR STRUCTURAL PARAMETERS CONTROL VARIABLES ***
            allocate(columnname(numcol))
            allocate(columnstring(numcol))
            allocate(cv_S%prior_cov_mode(numzones))
            allocate(cv_S%var_type(numzones))
            allocate(cv_S%struct_par_opt(numzones))
            allocate(cv_S%num_theta_type(numzones))
            allocate(cv_S%trans_theta(numzones))
            allocate(cv_S%alpha_trans(numzones))
            
            columnname(1:numcol)= (/ 'BetaAssoc', 'prior_cov_mode', 'var_type' , 'struct_par_opt', &
            &  'trans_theta', 'alpha_trans' /)
            columnstring(1:numcol) = ' '
            
            cv_S%prior_cov_mode =  1
            cv_S%var_type       =  1
            cv_S%struct_par_opt =  1
            cv_S%num_theta_type =  UNINIT_INT
            cv_S%trans_theta    =  0
            cv_S%alpha_trans    = 50

       end subroutine bdp_alloc_init_cv_S
       
 
       
!********  subroutine bdp_alloc_d_PM
       subroutine bdp_alloc_d_PM(d_PM,columnname,columnstring,nval,nrows)
         use bayes_pest_control
         use jupiter_input_data_support
       ! DECLARATIONS
         implicit none
         type(d_prior_mean),  intent(inout)                   :: d_PM
         character (len=COLWIDTH), pointer, intent(inout) :: columnname(:)
         character (len=COLWIDTH), pointer, intent(inout) :: columnstring(:) 
         character (len=3)                                    :: instr        
         integer,                               intent(in)    :: nval, nrows
         integer                                              :: i !local counters
         allocate (columnname(nval))
         allocate (columnstring(nval))
         columnname(1)='BetaAssoc'
         columnstring(1)=' '
         columnname(2)='Partrans'
         columnstring(2)=' '
         columnname(3)='alpha_trans'
         columnstring(3)=' '
         columnname(4)='beta_0'
         columnstring(4)=' '
         allocate (d_PM%beta_0(nrows))          
            d_PM%beta_0 = UNINIT_REAL ! an array 
         allocate (d_PM%Qbb(nrows,nrows))
         allocate (d_PM%InvQbb(nrows,nrows))
         allocate (d_PM%InvQbbB0(nrows))
            d_PM%Qbb = 0. ! an array
            d_PM%InvQbb   = UNINIT_REAL ! an array
            d_PM%InvQbbB0 = UNINIT_REAL ! an array
          do i=5,nval
           call int2char(i-4,instr)
           write(columnname(i),*) 'beta_cov_',instr
           columnstring(i)=' '
          enddo 
         allocate (d_PM%alpha(nrows))
         d_PM%alpha = 50.  
       end subroutine bdp_alloc_d_PM
       
       
       
       !********  subroutine bdp_alloc_Partrans
       subroutine bdp_alloc_Partrans(d_PM,columnname,columnstring,nval,nrows)
         use bayes_pest_control
         use jupiter_input_data_support
       ! DECLARATIONS
         implicit none
         type(d_prior_mean),  intent(inout)                   :: d_PM
         character (len=COLWIDTH), pointer, intent(inout) :: columnname(:)
         character (len=COLWIDTH), pointer, intent(inout) :: columnstring(:) 
         character (len=3)                                    :: instr        
         integer,                               intent(in)    :: nval, nrows
         integer                                              :: i !local counters
         allocate (columnname(nval))
         allocate (columnstring(nval))
         columnname(1)='BetaAssoc'
         columnstring(1)=' '
         columnname(2)='Partrans'
         columnstring(2)=' '
         columnname(3)='alpha_trans'
         columnstring(3)=' '
         allocate (d_PM%alpha(nrows))
         d_PM%alpha = 50.         
       end subroutine bdp_alloc_Partrans
       
!********  subroutine bdp_init_keyword_vars
        subroutine bdp_init_keyword_vars(BL,blnum)
       !  to allocate and initialize keywordline and keywordstring
       ! variables for reading keywords from input file 
            use bayes_pest_control
       ! ** DECLARATIONS
            implicit none
        type (tp_block)      :: BL(NUM_BLOCK)
            integer, intent(in)                              :: blnum
            integer                                          :: i
       ! ** ALLOCATIONS
            allocate(BL(blnum)%keywordline(BL(blnum)%numkw))
            allocate(BL(blnum)%keywordstring(BL(blnum)%numkw))
       ! ** INITIALIZATIONS
            BL(blnum)%keywordline   = UNINIT_INT  !this is an array
            BL(blnum)%keywordstring = UNINIT_CHAR !this is an array
        end subroutine bdp_init_keyword_vars         

!********  Subroutine INT2CHAR writes an integer to a string.
         subroutine int2char(ival,astring)
         implicit none
         integer, intent(in)            :: ival
         character (len=*), intent(out) :: astring
         character (len=7)              :: afmt

         afmt='(i    )'
         write(afmt(3:6),'(i4)') len(astring)
         write(astring,afmt) ival
         astring=adjustl(astring)
         return

       end subroutine int2char
       
end module bdp_data_parsers
