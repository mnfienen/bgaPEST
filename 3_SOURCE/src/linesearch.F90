module linesearch
 
   use jupiter_input_data_support
   use bayes_pest_control
   use model_input_output
   use bayes_pest_mio_setup
   use bayes_pest_model
   use bayes_pest_reader
   use bayes_pest_finalize
   use error_message
   use utilities
   use make_kernels
   use bayes_matrix_operations
   use param_trans
   use nelder_mead_simplex_routines
 contains
 
 subroutine lns_proc(d_XQR,d_S,cv_PAR,d_A,d_PAR,d_PM,cv_OBS,d_OBS,cv_PM,d_MOD,cv_A,it_phi,miostruc,errstruc)
 
 implicit none
        ! declarations
        type(kernel_XQR),     intent(in)    :: d_XQR
        type(d_struct),       intent(inout) :: d_S
        type(cv_param),       intent(in)    :: cv_PAR
        type (d_comlin),      intent(in)    :: d_MOD
        type(d_algorithmic),  intent(inout) :: d_A
        type(d_param),        intent(inout) :: d_PAR
        type(d_prior_mean),   intent(in)    :: d_PM
        type(cv_observ),      intent(in)    :: cv_OBS
        type(d_observ),       intent(in)    :: d_OBS
        type (cv_prior_mean), intent(in)    :: cv_PM
        type (mio_struc)                    :: miostruc
        type (err_failure_struc)            :: errstruc  
        type (cv_algorithmic),intent(inout) :: cv_A
        integer ,             intent(in)    :: it_phi       
        character (len=ERRORWIDTH)          :: retmsg
         
        integer ( kind = 4 ), parameter     :: n = 1 ! indicates the number of params to estimate. only rho here
        integer ( kind = 4 ) icount
        integer ( kind = 4 ) ifault
        integer ( kind = 4 ) konvge
        integer ( kind = 4 ) numres
        real    ( kind = 8 ) reqmin
        real    ( kind = 8 ) start(n)
        real    ( kind = 8 ) step(n)
        real    ( kind = 8 ) xmin(n)
        real    ( kind = 8 ) ynewlo
        real    ( kind = 8 ), external      :: Jmin
        
        start(1) = 0.5D+00
        step(1) = 0.5D+00
        konvge = 1
        reqmin = 1.0D-02

        call nelmin_ls ( Jmin,n,start,xmin,ynewlo,reqmin,step,konvge,cv_A%it_max_lns,icount,numres,ifault, & 
        & d_XQR,d_S,cv_PAR,cv_A,d_A,d_PAR,d_PM,cv_OBS,d_OBS,cv_PM,d_MOD,miostruc,errstruc)
        
        if (ifault.eq.2) then !Maximum number of iterations has exceeded --> Warning               
        write(retmsg,10) it_phi   
10        format('Warning: Maximum number of iterations exceeded in linesearch procedure during quasi-linear iteration',i4, & 
                 & '. The value that gives the minimum obj. funct. during the procedure will be considered.') 
          call utl_writmess(6,retmsg) 
        endif
        
        d_PAR%pars=xmin(1)*d_PAR%pars+(1-xmin(1))*d_PAR%pars_old 
        d_PAR%phi_T=ynewlo

 end subroutine lns_proc
 
end module linesearch 

 
   real ( kind = 8 )function Jmin (rho,d_XQR,d_S,cv_PAR,cv_A,d_A,d_PAR,d_PM,cv_OBS,d_OBS,cv_PM,d_MOD,miostruc,errstruc) 

   use jupiter_input_data_support
   use bayes_pest_control
   use model_input_output
   use bayes_pest_mio_setup
   use bayes_pest_model
   use bayes_pest_reader
   use bayes_pest_finalize
   use error_message
   use utilities
   use make_kernels
   use bayes_matrix_operations
   use param_trans

   implicit none
   
   double  precision rho
   type(kernel_XQR),    intent(in)     :: d_XQR
   type(d_struct),      intent(inout)  :: d_S
   type(cv_param),      intent(in)     :: cv_PAR
   type (d_comlin)                     :: d_MOD
   type(cv_algorithmic), intent(inout) :: cv_A
   type(d_algorithmic), intent(inout)  :: d_A
   type(d_param),       intent(inout)  :: d_PAR
   type(d_prior_mean),  intent(in)     :: d_PM
   type(cv_observ),     intent(in)     :: cv_OBS
   type(d_observ),      intent(in)     :: d_OBS
   type (cv_prior_mean), intent(in)    :: cv_PM
   type (mio_struc)                    :: miostruc
   type (err_failure_struc)            :: errstruc  
 
   ! -- INITIALIZATION 
   d_PAR%pars_lns=rho*d_PAR%pars+(1.-rho)*d_PAR%pars_old 
  
  !-- BACK-TRANSFORM OR NOT PARAMETERS IN THE PHYSICAL SPACE BEFORE THE FORWARD RUN
  if (maxval(d_PM%Partrans).ge.1) then  !If yes, we need to back-transform the parameters in the physical space  
     call par_back_trans_lns(cv_PAR, d_PAR, d_PM)           
  endif 
    
  !-- RUN THE FORWARD MODEL (INCLUDES DELETING OLD OUTPUT, WRITING NEW INPUT, RUNNING MODEL, AND READING NEW OUTPUT)
  call bpf_model_run(errstruc, d_MOD, cv_PAR,d_PAR, cv_OBS,  cv_A, d_OBS, d_A, 3, miostruc)
  !-- CALCULATE THE OBJECTIVE FUNCTIONS 
  call cal_ob_funcs(d_XQR, d_S, d_PM, cv_PAR, cv_OBS, d_OBS, d_A, d_PAR, cv_PM)
  
  Jmin = d_PAR%phi_T(1)
  return
end function
     