module bayes_pest_reader
    use utilities
    contains
      subroutine bpr_read(errstruc,infile,cv_A, d_A, cv_PM, d_PM, cv_S, d_S, cv_PAR, Q0_All, d_PAR, cv_OBS, &
                 &    d_OBS, d_MOD, cv_MIO, d_MIO, d_ANI, miostruc)
      
! --  Program BPR_READ reads from the input file
!     Required variables are all read or, if not found, appropriate errors are returned
!     In the case of optional variables, they are read or, if not present, default values are used.

!     a m!ke@usgs joint
!     mike fienen - mnfienen@usgs.gov  
!     V 0.0  3/5/08
!     Modified by M.D. 20/9/09

! ** DECLARATIONS **
       use jupiter_input_data_support
       use bayes_pest_control
       use bdp_data_parsers
       use bpi_initializations      
       use utilities
       use model_input_output
       use error_message
       
       implicit none
       type (mio_struc)           :: miostruc
       type (err_failure_struc)   :: errstruc
       integer                    :: ifail,inunit,i,j,n1,numlist,numcol,line,ierr,junit
       character (len=COLWIDTH)   :: columnname(MAXCOL),columnstring(MAXCOL)
       character (len=FILEWIDTH)  :: infile,outfile,filename
       character (len=ERRORWIDTH) :: errmsg
       type (cv_algorithmic), intent(inout) :: cv_A
       type (d_algorithmic),  intent(inout) :: d_A
       type (cv_prior_mean) , intent(inout) :: cv_PM
       type (d_prior_mean)  , intent(inout) :: d_PM  
       type (cv_struct)     , intent(inout) :: cv_S  
       type (d_struct)      , intent(inout) :: d_S 
       type (cv_param)      , intent(inout) :: cv_PAR
       type (Q0_compr)      , pointer       :: Q0_All(:)
       type (d_param)       , intent(inout) :: d_PAR
       type (cv_observ)     , intent(inout) :: cv_OBS
       type (d_observ)      , intent(inout) :: d_OBS
       type (d_comlin)      , intent(inout) :: d_MOD
       type (cv_minout)     , intent(inout) :: cv_MIO
       type (d_minout)      , intent(inout) :: d_MIO
       type (tp_block)                      :: BL(NUM_BLOCK)
       type (d_anisotropy)                  :: d_ANI
! ** INITIALIZATIONS **

       
! initialize all the BLOCK types
       call bpi_init_algorithmic_CVs(BL,cv_A)
       call bpi_init_prior_mean_CVs(BL,cv_PM)
       call bpi_init_prior_mean_DATA(BL,d_PM)
       call bpi_init_struct_CVs(BL,cv_S)
       call bpi_init_struct_DATA(BL,d_S)
       call bpi_init_param_CVs(BL,cv_PAR)
       call bpi_init_param_DATA(BL,d_PAR)
       call bpi_init_obs_groups(BL,cv_OBS)
       call bpi_init_obs_DATA(BL,d_OBS)
       call bpi_init_modcomlin_DATA(BL,d_MOD)
       call bpi_init_mio_CVs(BL,cv_MIO)
       call bpi_init_anisotropy_DATA(BL,d_ANI)
       if(mio_initialise(errstruc,miostruc).ne.0) then
        call utl_bomb_out(errstruc)
        n1=mio_finalise(errstruc,miostruc)
        stop 
       endif
        
! other variables
       inunit = utl_nextunit()
       errmsg=UNINIT_CHAR
! OPEN INPUT FILE
       call bpc_openfile(inunit,trim(infile),0)   ![0] at end indicates read only
       
! CHECK STATUS OF BLOCKS
       junit=utl_nextunit()   ! junit is used for reading subsidiary files
       call ids_block_status(ifail,inunit,junit,NUM_BLOCK,BL%label,BL%numrows,infile)
       if (ifail.ne.0) then
         call ids_get_message_string(errmsg)
         call utl_writmess(6,errmsg)
       endif
       
 ! READ THE VARIABLES FROM THE INPUT FILE
       
       call bdp_read_cv_algorithmic(BL,cv_A,inunit,errmsg)
       
       call bdp_read_cv_prior_mean (BL,cv_PM,inunit,errmsg)
       call bdp_read_data_prior_mean (BL,cv_PM,d_PM,cv_PAR,inunit,errmsg)
              
       call bdp_read_cv_structural_parameters_tbl (BL,cv_S,cv_PAR%p,inunit,errmsg)
       call bdp_read_data_structural_parameters (BL,cv_S,d_S,cv_A,cv_PAR%p,inunit,errmsg)
       if (cv_A%theta_cov_form.ne.0) then
          call bdp_read_structural_parameters_cov (BL,cv_S,d_S,cv_A%theta_cov_form,inunit,errmsg)
       endif
       call bdp_read_epistemic_error (BL,d_S,cv_A,cv_S,inunit,errmsg)
       call bdp_read_cv_parameters(BL,cv_PAR,inunit,errmsg)
       
       if (cv_A%Q_compression_flag.ne.0) then
          allocate(Q0_All(cv_PAR%p))
          call bdp_read_cv_compression(BL,Q0_All,cv_PAR%p,inunit,errmsg)
       else 
          allocate(Q0_All(0)) !We need this because of the external function of struct pars optimization (No memory is required)   
       endif
              
       call bdp_read_parameter_groups(BL,cv_PAR,inunit,errmsg)
       call bdp_read_data_parameters(errstruc,BL,d_PAR,cv_PAR,cv_A,Q0_All, &
                            & inunit,errmsg,miostruc)
       if (cv_PAR%npar .gt. 100000) then
         cv_A%store_Q = .FALSE.
       endif
       
       call bdp_read_observation_groups(BL,cv_OBS,inunit,errmsg)                            
       call bdp_read_data_observations(BL,d_OBS, cv_OBS, &
                            & inunit,errmsg,miostruc)
       call bpi_init_algorithmic_DATA(d_A,cv_PAR%npar,cv_OBS%nobs) !Allocate memory and initialize H matrix (nobs x npar)
       call bdp_read_data_model_command_line(BL,d_MOD,cv_A%deriv_mode,inunit,errmsg)       
       call bdp_read_data_model_input_output(BL,d_MIO,cv_MIO,cv_A,cv_PAR%npargp,inunit,errmsg)
       if (cv_A%par_anisotropy .eq. 1) then
           call bdp_read_data_par_anisotropy(BL,d_ANI,cv_PAR,inunit,errmsg)
       endif
! Deallocate BL type and close up files
       call deallocate_BL(BL)
        close (inunit)

  
! ******************        
! * SUBROUTINES
! ******************  

contains     
!********  subroutine deallocate_BL
       subroutine deallocate_BL(BL)
       ! SUBROUTINE to deallocate the pointer strings within the keywords group (BL) type
            use bayes_pest_control
       ! DECLARATIONS  
             implicit none
             type(tp_block),   intent(inout)    :: BL(NUM_BLOCK) 
             integer                            :: i
       ! deallocation
           do i=1,NUM_BLOCK
             if (associated(BL(i)%keywords))      deallocate (BL(i)%keywords)
             if (associated(BL(i)%keywordstring)) deallocate (BL(i)%keywordstring)
             if (associated(BL(i)%keywordline))   deallocate (BL(i)%keywordline)
             
             nullify(BL(i)%keywords,BL(i)%keywordstring ,BL(i)%keywordline)
           enddo
       end subroutine deallocate_BL

      end subroutine bpr_read 
end module bayes_pest_reader