program bp_main

! *********************************************
! ***     MAIN BAYES MODULE PEST PROGRAM    ***
! ***      Implementation of: Bayesian      ***
! *** Geostatistical inverse method in PEST.***
! ***                                       ***
! ***             a m!ke@usgs joint         ***
! ***          Michael N. Fienen, PhD       ***
! ***    UNITED STATES GEOLOGICAL SURVEY    ***
! ***          mnfienen@usgs.gov            ***
! ***                 AND                   ***
! ***           Marco D'Oria, PhD           ***
! ***    UNIVERSITY OF PARMA, ITALY         ***
! ***       marco.doria@unipr.it            ***
! ***           January 3, 2012             ***
! ***      Modified by M.D. 1/10/2009       ***
! ***  Further Modifications MNF/MD 11/2010 ***
! ***  Further Modifications   MD   09/2011 ***
! ***  Further Modifications MNF/MD 12/2011 ***
! ***   Version 1.0 --- November 19, 2012   ***
! *********************************************
   
   use jupiter_input_data_support
   use bayes_pest_control
   use model_input_output
   use bayes_pest_mio_setup
   use bayes_pest_model
   use bayes_pest_reader
   use bayes_pest_finalize
   use bayes_output_control
   use error_message
   use utilities
   use make_kernels
   use bayes_matrix_operations
   use param_trans
   use linesearch
   use extern_derivs
   use bayes_output_control
   use struct_param_optimization
   use posterior_cov_operations
   implicit none


!--  Main Data Arrays for OBS and PARS and ALGORITHM   
       integer                      :: n1, i, j
       integer                      :: forward_flag_der, ci95_flag
       integer                      :: s_ind, p_ind, b_ind  !Indices for structural parameters, quasi-linear and bga method loops
       type (mio_struc)             :: miostruc
       type (err_failure_struc)     :: errstruc  
       integer                      :: ifail, restart, outunit, bprunit, cparunit, cobsunit, finalparunit, postcovunit
       type (cv_algorithmic)        :: cv_A 
       type (d_algorithmic)         :: d_A
       type (cv_prior_mean)         :: cv_PM
       type (d_prior_mean)          :: d_PM  
       type (cv_struct)             :: cv_S  
       type (d_struct)              :: d_S 
       type (cv_param)              :: cv_PAR
       type (Q0_compr), pointer     :: Q0_All(:)
       type (d_param)               :: d_PAR
       type (cv_observ)             :: cv_OBS
       type (d_observ)              :: d_OBS
       type (d_comlin)              :: d_MOD
       type (cv_minout)             :: cv_MIO
       type (d_minout)              :: d_MIO
       type (kernel_XQR)            :: d_XQR
       type (d_anisotropy)          :: d_ANI
       character (len=ERRORWIDTH)   :: retmsg
       character (len=100)          :: command_line, curr_par_file, curr_resid_file,post_cov_file
       character (len=FILEWIDTH)    :: ctlfile
       character (len=FILEWIDTH)    :: casename
       character (len=FILEWIDTH)    :: atemp
       character (len=20)           :: inner_iter ! aka p_ind - this is the temporary holder for printing out the inner iteration number
       character (len=20)           :: outer_iter ! aka b_ind - this is the temporary holder for printing out the outer iteration number 
       double precision,dimension(1) :: curr_structural_conv, curr_phi_conv       !Current iteration convergence values for 
       double precision,dimension(1) :: curr_bga_conv, curr_bga_phi,prev_bga_phi  !structural parameters and quasi linear objective function 
       double precision,dimension(1) :: curr_phi !Current value for quasi linear objective function
       double precision, allocatable :: prev_struct(:) !Previous vector of theta and sigma values to be optimized for or previous objective function
       double precision, pointer    :: VV(:,:), V(:) !VV is the posterior covariance matrix, V is only the diagonal of VV 
       double precision             :: structural_conv
       double precision             :: huge_val=huge(huge_val) !Largest machine number

   nullify(Q0_All)
   nullify(VV)
   nullify(V)
! -- PRINT OUT THE BANNER INFORMATION
   call bpo_write_banner()

!-- READ AND PARSE THE COMMAND LINE TO OBTAIN THE CONTROL FILE NAME       
  call UTL_GET_COMMAND_LINE(COMMAND_LINE)
  call UTL_PARSE_COMMAND_LINE(IFAIL,COMMAND_LINE,CTLFILE,RESTART)

!-- handle the case where no control file was indicated
  IF (IFAIL.NE.0) THEN
    call bpo_write_nocmd()
    stop
  endif
  
  
! -- An extension of ".bgp" is added to the bgaPEST control file if necessary.
      i=LEN_TRIM(CTLFILE)
      IF(i.GE.5) THEN
        ATEMP=CTLFILE(I-3:I)
        CALL UTL_CASETRANS(ctlfile,'lo')
        IF(ATEMP.NE.'.bgp') THEN 
          CASENAME = CTLFILE
          CTLFILE(I+1:)='.bgp'
        ELSE
          CTLFILE = trim(CTLFILE)
          CASENAME = CTLFILE(1:I-4)
        ENDIF
      ELSE
          CASENAME = trim(CTLFILE)
          CTLFILE  = trim(CTLFILE) // '.bgp'
      ENDIF
   
!--  INITIALIZE MIO STRUCTURE SO IT CAN BE PASSED    
    if(mio_initialise(errstruc,miostruc).ne.0) then  !MD mio_initialise is an integer
      call utl_bomb_out(errstruc)
      n1=mio_finalise(errstruc,miostruc)
      stop 
    endif  


! open up the main output record file 
       bprunit = utl_nextunit()          
       call bpc_openfile(bprunit,trim(trim(casename) // '.bpr'),1) ![1] at end indicates open with write access

!--  READ INPUT FILE AND PERFORM ASSOCIATED ALLOCATIONS AND PARSING     
    call bpr_read(errstruc,CTLFILE,cv_A, d_A, cv_PM, d_PM, cv_S, d_S, cv_PAR,Q0_All, d_PAR, &
                        cv_OBS, d_OBS, d_MOD, cv_MIO, d_MIO, d_ANI,miostruc)

!--  SETUP THE MODEL INPUT AND OUTPUT INFORMATION (TPL/INPUT AND INS/OUTPUT PAIRINGS)
    call bpm_setup_mio(errstruc, cv_MIO, d_MIO,  cv_PAR%npargp, &
                        cv_OBS%nobsgp, cv_PAR%npar, cv_OBS%nobs, miostruc)
                     
!--  INITIALIZE THE INVERSE MODEL
    !Make Q0, R0, X0, and InvQbb if necessary   
    call bxq_make_X0_Q0_R0_InvQbb(d_PAR,cv_S,d_S,cv_PAR,d_XQR,cv_A,d_OBS,cv_OBS%nobs,d_PM,Q0_All,cv_PM,d_ANI)
   
    allocate(d_OBS%h(cv_OBS%nobs)) ! Allocate the current model output [y]
    
!-- IF STRUCTURAL PARAMETERS WILL BE OPTIMIZED FOR, SET UP REQUIRED INFORMATION
    if ((maxval(cv_S%struct_par_opt).eq.1).or.(d_S%sig_opt.eq.1)) then
      call bxq_theta_cov_calcs(cv_PAR,cv_S,d_S,cv_PM,cv_A)
      if (cv_A%structural_conv.ge.0.) then
        allocate(prev_struct(1))
        prev_struct = UNINIT_REAL
      else  
        allocate(prev_struct(cv_S%num_theta_opt))
        prev_struct = d_S%struct_par_opt_vec_0
      endif
    endif
    
!-- CALL THE SETUP OF EXTERNAL DERIVATIVES FILES (IF REQUIRED).  THIS HAPPENS ONLY ONCE FOR ALL BUT PARAMETERS FILE
    if ((cv_A%deriv_mode .eq. 0) .or. (cv_A%deriv_mode .eq. 4)) then
        call bxd_write_ext_PEST_files(d_MOD, cv_MIO, d_MIO, cv_OBS, cv_PAR, d_OBS,cv_A)
    endif
    
!-- WRITE THE HEADER INFORMATION TO THE BPR Run Record FILE
    call bpo_write_bpr_header(bprunit,casename,cv_PAR,cv_OBS,d_MOD, cv_A, &
                cv_MIO, d_MIO,Q0_all,cv_PM,d_PM,cv_S,d_S,d_PAR,d_ANI)

!!! --- Initialize outer loop convergence values
          curr_bga_conv = huge_val ! initialize current outer loop convergence
          prev_bga_phi  = huge_val ! initialize current outer loop convergence
          curr_bga_phi  = huge_val ! initialize current bga objective function 

    do b_ind = 1, cv_A%it_max_bga  !*********************************************************************** (more external loop)
    
    !***************************************************************************************************************************  
    !****************************** FROM HERE THE QUASI-LINEAR PARAMETER ESTIMATION LOOP ***************************************
    !***************************************************************************************************************************
    
          curr_phi_conv = huge_val !Initialize current quasi linear objective function convergence
          curr_phi      = huge_val !Initialize current quasi-linear objective function value

          do p_ind = 1, cv_A%it_max_phi !************************************************************* (first intermediate loop)
                                        !********** quasi-liner parameter estimation for given structural parameters ***********
      
             !-- RUN THE FORWARD MODEL (INCLUDES DELETING OLD OUTPUT, WRITING NEW INPUT, RUNNING MODEL, AND READING NEW OUTPUT)
             select case(cv_A%deriv_mode)
                case (0)
                    forward_flag_der = 1
                case (1)
                    forward_flag_der = 2
				case (4)
                    forward_flag_der = 4 
             end select
             call bpf_model_run(errstruc, d_MOD, cv_PAR,d_PAR, cv_OBS, cv_A,  d_OBS, d_A, forward_flag_der, miostruc)
     
            d_PAR%pars_old = d_PAR%pars   !MD At the beginning pars is the vector of the initial values of the parameters 
                                          !as read in the parameter file. Then became the best estimate. 
            
            !-- CONVERT OR NOT SENSITIVITY AND PARAMETERS IN THE ESTIMATION SPACE
            if (maxval(d_PM%Partrans).ge.1) then  !If yes, the parameters transformation is required  
               call sen_par_trans(cv_PAR, cv_OBS, d_PAR, d_A, d_PM) !Converting sensitivity and parameters in the estimation space
            endif
          
            !-- SOLVE THE BAYESIAN LINEAR SYSTEM AND CALCULATE THE OBJECTIVE FUNCTIONS           
            call  bmo_form_Qss_Qsy_HQsy(d_XQR, d_S%theta, cv_PAR, cv_OBS, cv_S, cv_A, d_A, d_PAR, Q0_All)
            call  bmo_form_Qyy(d_XQR, d_S%sig, cv_OBS, d_A)
            call  bmo_H_only_operations(d_XQR, d_A,cv_OBS,d_PAR,cv_PAR)
            call  bmo_solve_linear_system(d_XQR, d_S, d_PM, cv_PAR, cv_OBS, d_OBS, d_A, d_PAR,cv_PM)
            
            
            !-- PERFORM LINESEARCH IF REQUESTED
            if (cv_A%lns_flag.eq.1) then  !If yes, we perform the linesearch procedure  
               call lns_proc(d_XQR,d_S,cv_PAR,d_A,d_PAR,d_PM,cv_OBS,d_OBS,cv_PM,d_MOD,cv_A,p_ind,miostruc,errstruc)
            endif
            
            !-- BACK-TRANSFORM OR NOT PARAMETERS INTO PHYSICAL SPACE
            if (maxval(d_PM%Partrans).ge.1) then  !If yes, we need to back-transform the parameters in the physical space  
               call par_back_trans(cv_PAR, d_PAR, d_PM)
            endif 
            
            !Run the forward model to obtain the current modeled observation vector (consistent with the estimated parameters)  
            call bpf_model_run(errstruc, d_MOD, cv_PAR,d_PAR, cv_OBS, cv_A,  d_OBS, d_A, 0, miostruc)
            
            ! UPDATE THE OBJECTIVE FUNCTION COMPONENTS (in case of linesearch)
            if (cv_A%lns_flag.eq.1) then 
               call  cal_ob_funcs(d_XQR, d_S, d_PM, cv_PAR, cv_OBS, d_OBS,  d_A, d_PAR, cv_PM)
            endif
            
            !-- set temporary string version of iteration numbers and phi to write out
            curr_phi_conv = abs(curr_phi - d_PAR%phi_T) 
            curr_phi = d_PAR%phi_T
            
            call UTL_INT2CHAR(p_ind,inner_iter)
            call UTL_INT2CHAR(b_ind,outer_iter)  
            curr_par_file = trim(casename) // '.bpp.' // trim(outer_iter) // '_' // trim(inner_iter)
            curr_resid_file = trim(casename) // '.bre.' // trim(outer_iter) // '_' // trim(inner_iter)
            !-- Write intermediate values out to BPR record file
            call bpo_write_bpr_intermed(bprunit,p_ind,b_ind,curr_par_file,curr_resid_file,d_PAR) 
            ! --Write the intermediate parameter and residuals files
      
            cparunit = utl_nextunit()
            call bpc_openfile(cparunit,trim(curr_par_file),1) ![1] at end indicates open with write access
            call bpo_write_allpars(cv_PAR,d_PAR,cparunit)
            close(cparunit)
            cobsunit = utl_nextunit()
            call bpc_openfile(cobsunit,trim(curr_resid_file),1) ![1] at end indicates open with write access
            call bpo_write_residuals(cv_OBS,d_OBS,cobsunit)
            close(cobsunit) 
            !-- check for convergence - exit if convergence has been achieved  
            if (curr_phi_conv(1) .le. cv_A%phi_conv) then
                exit
            elseif (p_ind .ge. cv_A%it_max_phi) then
                write(retmsg,10) p_ind
10              format('Warning: Maximum number of iterations exceeded in quasi-linear parameter optimization loop during bgaPEST iteration',i4, & 
                 & '. Convergence was not achieved, but iterations will cease.')
                call utl_writmess(6,retmsg)  
            endif !- checking for convergence or exceeding maximum iterations
            
          enddo  !(first intermediate loop) quasi-linear method  --> p_ind
          curr_bga_phi = d_PAR%phi_T  ! Set the current bga outer loop convergence to equal the current value of PHI_T
    !***************************************************************************************************************************  
    !************************************* END OF QUASI-LINEAR PARAMETER ESTIMATION LOOP ***************************************
    !***************************************************************************************************************************
     
    
    !***************************************************************************************************************************  
    !********************** FROM HERE THE STRUCTURAL PARAMETER ESTIMATION LOOP  (ONLY IF REQUIRED) *****************************
    !************************************************************************** *************************************************
       if ((maxval(cv_S%struct_par_opt).eq.1).or.(d_S%sig_opt.eq.1)) then !Enter the structural pars estimation loop only if required
          if (b_ind .ge. cv_A%it_max_bga) then !-- do not re-estimate structural parameters if we have exceeded maximum number of it_max_bga
             write(retmsg,11) b_ind
11           format('Warning: Maximum number of iterations exceeded in quasi-linear parameter optimization loop during bgaPEST iteration',i4, & 
               & '. Structural parameters will not be re-calculated, so final structural parameters and posterior covariance (if requested)' &
               & ' are based on the last iteration.')
             call utl_writmess(6,retmsg) 
             cv_S%struct_par_opt = 1
               call bpo_write_bpr_final_structpar(bprunit,cv_S,d_S,cv_PAR,d_PAR)
             cv_S%struct_par_opt = 0
             exit
          else
             if ((cv_A%structural_conv.ge.0.).and.(b_ind.eq.1)) then            !In case of objective function monitoring, here the evaluation of
               call beg_str_object_fun(cv_OBS,d_OBS,d_A,cv_S,d_PM,cv_PAR,cv_PM) !the objective function with the initial parameters 
               prev_struct=cv_S%str_obj_fun                                     !(only at the first bga loop)
             endif                                                              
               
             call marginal_struct_param_optim(d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM,b_ind,cv_S%num_theta_opt)
             !Here d_S%struct_par_opt_vec is the vector of the optimized theta and sigma values
              
             if (cv_A%structural_conv.ge.0.) then
               curr_structural_conv=abs(prev_struct - cv_S%str_obj_fun) !Calculate difference between actual and previous objective function values
               structural_conv = cv_A%structural_conv                   !The structural convergenge value remain the one assigned into .bgp file
               prev_struct = cv_S%str_obj_fun                           !Assign the current value of the objective function to the previous
             else
               !curr_structural_conv = sqrt(sum((prev_struct - d_S%struct_par_opt_vec)**2)) !Calculate norm of difference between actual and previous vectors
               curr_structural_conv = sqrt(sum(((prev_struct - d_S%struct_par_opt_vec)/prev_struct)**2))!calculate the normalized root square difference 
               																							!between last iteration and current one
               structural_conv = -cv_A%structural_conv             !The structural convergenge value changes sign respect to the one assigned into .bgp file  
               prev_struct = d_S%struct_par_opt_vec                !Assign the current vector of structural parameters to the previous vector
             endif  
             
             if (curr_structural_conv(1).le.structural_conv) then !If yes, structural parameters have converged.
               call bpo_write_bpr_intermed_structpar(bprunit,cv_S,d_S,cv_PAR,d_PAR) ! write out the final structural parameters
               cv_S%struct_par_opt = 0  !Set to zero cv_S%struct_par_opt and d_S%sig_opt so the structural parameters estimation loop is no more entered.
               d_S%sig_opt = 0          !The optimized struct_par_opt_vec is used to run the quasi-linear loop that is the last one. 
               if (allocated(prev_struct)) deallocate(prev_struct) 
             endif
          endif !-- special warning if exceed maximum number of main algorithm iterations (it_max_bga) without convergence
       ! -- write the intermediate files to the BPR file
          call bpo_write_bpr_intermed_structpar(bprunit,cv_S,d_S,cv_PAR,d_PAR)
       else
         ! need to be sure that final values get written out.
         cv_S%struct_par_opt = 1
         call bpo_write_bpr_final_structpar(bprunit,cv_S,d_S,cv_PAR,d_PAR)
         cv_S%struct_par_opt = 0
         exit !If the structural pars optimization is not required or structural pars have converged (run the last quasi_linear), exit the bga_loop
       endif
    !***************************************************************************************************************************  
    !*************************** END OF STRUCTURAL PARAMETER ESTIMATION LOOP  (ONLY IF REQUIRED) *******************************
    !***************************************************************************************************************************  
    
    
    
    !*********************************************
    ! Evaluate outer bga convergence 
    !*********************************************
       curr_bga_conv = abs(prev_bga_phi - curr_bga_phi)
       if (curr_bga_conv(1) .le. cv_A%bga_conv) then
          exit
       else
         prev_bga_phi = curr_bga_phi
       endif
       
    enddo      !(more external loop) --> b_ind
    ! write out a final residuals file
    cobsunit = utl_nextunit()
    curr_resid_file = trim(casename) // '.bre.fin'
    call bpc_openfile(cobsunit,trim(curr_resid_file),1) ![1] at end indicates open with write access
    call bpo_write_residuals(cv_OBS,d_OBS,cobsunit)
    close(cobsunit) 
    
    
    !*************************************************************************************************************************
    !******** FROM HERE THE EVALUATION OF THE POSTERIOR COVARIANCE (ONLY IF REQUIRED --> cv_A%post_cov_flag = 1 **************
    !*********** The posterior covariance is the full matrix (matrix VV) in case of no compression of Q, *********************
    !********************* it is only the diagonal (vector V) in case of compression of Q ************************************
    !*************************************************************************************************************************
     if (cv_A%post_cov_flag.eq.1) then
       write(6,*) 'Calulating Posterior Covariance'
       ci95_flag = 1 ! - flag determining whether to write 95% confidence intervals or not
      call form_post_covariance(d_XQR, cv_PAR, cv_OBS, cv_S, cv_A, d_A, d_PAR,Q0_All,cv_PM,d_PM,d_S,VV,V)
      !-- write out the final parameter values and confidence intervals
        if (cv_A%Q_compression_flag .eq. 0) then  !Select if the Q0 matrix is compressed or not
          allocate(V(cv_PAR%npar))
          V = 0.D0 ! initialize the vector for V
          do i = 1,cv_PAR%npar
            V(i) = VV(i,i)
          enddo
        endif
        finalparunit = utl_nextunit()  
        curr_par_file = trim(casename) // '.bpp.fin'
        call bpc_openfile(finalparunit,trim(curr_par_file),1) ![1] at end indicates open with write access
        call bpo_write_allpars_95ci(cv_PAR,d_PAR,d_PM,V,finalparunit,ci95_flag)
        close(finalparunit)
        ! --- Also write a separate file with only the posterior covariance values
        postcovunit = utl_nextunit() 
        post_cov_file = trim(casename) // '.post.cov' 
        call bpc_openfile(postcovunit,trim(post_cov_file),1) ![1] at end indicates open with write access
        if (cv_A%Q_compression_flag.eq.0) then
          if (associated(V)) deallocate(V) 
        endif 
        call  bpo_write_posterior_covariance(cv_A%Q_compression_flag,cv_PAR,d_PAR,d_PM,V,VV,postcovunit)
        close(postcovunit)
     else ! still need to write the final parameters out
       ci95_flag = 0 ! - flag determining whether to write 95% confidence intervals or not
        allocate(V(1))
        V = 0.
        finalparunit = utl_nextunit()  
        curr_par_file = trim(casename) // '.bpp.fin'
        call bpc_openfile(finalparunit,trim(curr_par_file),1) ![1] at end indicates open with write access
        call bpo_write_allpars_95ci(cv_PAR,d_PAR,d_PM,V,finalparunit,ci95_flag)
        close(finalparunit)
        

        if (associated(V)) deallocate(V)
     endif
    !*************************************************************************************************************************
    !*********** END OF THE EVALUATION OF THE POSTERIOR COVARIANCE (ONLY IF REQUIRED --> cv_A%post_cov_flag = 1 **************
    !*************************************************************************************************************************
  
  
  write(*,*) 'Parameter estimation is complete!'
  
!-- FINALIZE and CLEANUP - deallocate all the structures
    call bpd_finalize(d_PM, d_S, cv_PAR, d_PAR, cv_OBS, d_OBS, d_MOD, d_MIO, d_XQR)


end program bp_main                   


 


