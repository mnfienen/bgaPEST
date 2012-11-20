module bpi_initializations

! Make the subroutines visible
        public  bpi_init_algorithmic_BLs,     &
                bpi_init_prior_mean_CVs,      &
                bpi_init_prior_mean_DATA,     &
                bpi_init_struct_CVs,          &
                bpi_init_struct_DATA,         &
                bpi_init_param_CVs,           &
                bpi_init_param_DATA,          &
                bpi_init_obs_CVs,             &
                bpi_init_obs_DATA,            &
                bpi_init_modcomlin_CVs,       &
                bpi_init_modcomlin_DATA,      &
                bpi_init_mio_data,            &
                bpi_init_mio_CVs
contains

!********  subroutine bpi_bpi_init_algorithmic_CVs
       subroutine bpi_init_algorithmic_CVs(BL,cv_A)
       ! SUBROUTINE to initialze algorithmic control variables
         use bayes_pest_control
       ! DECLARATIONS  
         implicit none
         type(cv_algorithmic),intent(inout)  :: cv_A
         type(tp_block),   intent(inout)     :: BL(NUM_BLOCK) 
       ! INITIALIZATIONS
            cv_A%structural_conv    = 0.001  !MD Structural parameter convergence values  
            cv_A%phi_conv           = 0.001  !MD Objective function convergence value
            cv_A%bga_conv           = 0.01   !outer loop convergence
            cv_A%it_max_structural  = 10     !MD Max number of iterations for structral parameter estimation    	    
            cv_A%it_max_phi         = 10     !MD Max number of iterations for quasi-linear estimation method  
            cv_A%it_max_bga         = 10     !MNF outer iterations number
            cv_A%lns_flag           = 0      !MD Linesearch procedure flag: [0] not perform [1] perform 
            cv_A%it_max_lns         = 4      !MD Max number of iterations for linesearch procedure
            cv_A%store_Q            = .TRUE. !flag for whether or not to store Q.  s/b false if npar is too big
            cv_A%theta_cov_form     = 0	     !Form of theta covariance:  [0] none, [1] diag, [2] full matrix
            cv_A%Q_compression_flag = 0      ! [0] none - calculate full Q0, [1] Calculate Q0 for each beta
            cv_A%deriv_mode         = 0      ! [0] use external PEST for Jacobian, [1] look for separate dercomline, [4] is parallel with Condor
            cv_A%post_cov_flag      = 0      ! [0] do not calculate posterior covariance [1] calculate post cov
            cv_A%jacobian_format    = 'binary' ! 'binary' [default] means read JCO binary, 'ascii' means read standard PEST text matrix
            cv_A%jacfle             = 'scratch.jco' ! filename for external derivatives to be read from
            cv_A%par_anisotropy     = 0      ! flag for whether anisotropy block will be read [0]=no, [1] = yes
            BL(1)%label           = 'algorithmic_cv'
            BL(1)%numrows         = UNINIT_INT
            BL(1)%numkw           = 16 
            allocate (BL(1)%keywords(BL(1)%numkw))

            BL(1)%keywords = (/'structural_conv','phi_conv','bga_conv',   &
            & 'it_max_structural','it_max_phi','it_max_bga','linesearch', &
            & 'it_max_linesearch', 'theta_cov_form',                      &
            & 'Q_compression_flag', 'store_Q', 'posterior_cov_flag',      &
            & 'deriv_mode','jacobian_format','jacobian_file','par_anisotropy'/)   
       end subroutine bpi_init_algorithmic_CVs  
        
       
!********  subroutine bpi_init_prior_mean_CVs
       subroutine bpi_init_prior_mean_CVs(BL,cv_PM)
       ! SUBROUTINE to initialze algorithmic control variables
         use bayes_pest_control
       ! DECLARATIONS  
         implicit none
         type(cv_prior_mean),intent(inout)   :: cv_PM   
         type(tp_block),   intent(inout)     :: BL(NUM_BLOCK) 
       ! INITIALIZATIONS
       ! ** N.B. ** most of these keywords/variables are pointer so they are not allocated here 
            cv_PM%betas_flag = 0    ! [0] means no prior information about beta
            cv_PM%Qbb_form   = 0	!Form of Beta covariance:  [0] none, [1] diag, [2] full matrix
            BL(2)%label           = 'prior_mean_cv'
            BL(2)%numrows         = UNINIT_INT
            BL(2)%numkw           =  2
            allocate (BL(2)%keywords(BL(2)%numkw))
            BL(2)%keywords = (/'prior_betas','beta_cov_form'/)
       end subroutine bpi_init_prior_mean_CVs  
               
               
!********  subroutine bpi_init_prior_mean_DATA
       subroutine bpi_init_prior_mean_DATA(BL,d_PM)
       ! SUBROUTINE to initialze algorithmic control variables
         use bayes_pest_control
       ! DECLARATIONS  
         implicit none
         type(d_prior_mean),intent(inout)    :: d_PM
         type(tp_block),   intent(inout)     :: BL(NUM_BLOCK) 
       ! INITIALIZATIONS
            nullify(d_PM%beta_0)
            nullify(d_PM%Qbb)
            nullify(d_PM%InvQbb)
            nullify(d_PM%InvQbbB0)
            nullify(d_PM%Partrans)
            nullify(d_PM%alpha)
            BL(3)%label           = 'prior_mean_data'
            BL(3)%numrows         = UNINIT_INT
            BL(3)%numkw           = 0 
       end subroutine bpi_init_prior_mean_DATA
       
!********  subroutine bpi_init_struct_CVs
       subroutine bpi_init_struct_CVs(BL,cv_S)
       ! SUBROUTINE to initialze structural parameter control variables
         use bayes_pest_control
       ! DECLARATIONS  
         implicit none
         type(cv_struct),intent(inout)    :: cv_S
         type(tp_block),intent(inout)     :: BL(NUM_BLOCK) 
       ! INITIALIZATIONS   
            nullify(cv_S%prior_cov_mode)
            nullify(cv_S%var_type)
            nullify(cv_S%struct_par_opt)
            nullify(cv_S%num_theta_type)
            nullify(cv_S%trans_theta)
            nullify(cv_S%alpha_trans)
            BL(4)%label         = 'structural_parameter_cv'
            BL(4)%numrows       = UNINIT_INT
            BL(4)%numkw         = 0
       end subroutine bpi_init_struct_CVs 
       
!********  subroutine bpi_init_struct_DATA
       subroutine bpi_init_struct_DATA(BL,d_S)
       ! SUBROUTINE to initialze structural parameter data
         use bayes_pest_control
       ! DECLARATIONS  
         implicit none
         type(tp_block),   intent(inout)     :: BL(NUM_BLOCK)
         type(d_struct),   intent(inout)     :: d_S 
       ! INITIALIZATIONS
            nullify(d_S%theta_0)
            nullify(d_S%struct_par_opt_vec_0)
            nullify(d_S%struct_par_opt_vec)
            nullify(d_S%theta_cov)
            nullify(d_S%invQtheta)
            nullify(d_S%theta)
            BL(5)%label           = 'structural_parameters_data'
            BL(5)%numrows         = UNINIT_INT
            BL(5)%numkw           = 0 
            BL(6)%label           = 'structural_parameters_cov'
            BL(6)%numrows         = UNINIT_INT
            BL(6)%numkw           = 0 
            
            d_S%sig_0     = UNINIT_REAL
            d_S%sig_opt   = 1
            d_S%sig_p_var = UNINIT_REAL
            d_S%trans_sig = 0
            d_S%alpha_trans_sig = 50.
            BL(7)%label           = 'epistemic_error_term'
            BL(7)%numrows         = UNINIT_INT
            BL(7)%numkw           = 5
            allocate (BL(7)%keywords(BL(7)%numkw))
            BL(7)%keywords = (/'sig_0', 'sig_opt' , 'sig_p_var', 'trans_sig', 'alpha_trans'/) 
       end subroutine bpi_init_struct_DATA 
       
!********  subroutine bpi_init_param_CVs
       subroutine bpi_init_param_CVs(BL,cv_PAR)
       ! SUBROUTINE to initialze algorithmic control variables
         use bayes_pest_control
       ! DECLARATIONS  
         implicit none
         type(cv_param),intent(inout)     :: cv_PAR
         type(tp_block),intent(inout)     :: BL(NUM_BLOCK) 
       ! INITIALIZATIONS
            nullify(cv_PAR%grp_name)
            nullify(cv_PAR%derinc)
            cv_PAR%npargp         = UNINIT_INT ! number of parameter groups
            cv_PAR%ndim           = UNINIT_INT ! number of dimensions (ok for 1, 2, or 3)         
            BL(8)%label           = 'parameter_cv'
            BL(8)%numrows         = UNINIT_INT
            BL(8)%numkw           = 1 
            allocate (BL(8)%keywords(BL(8)%numkw))
            BL(8)%keywords = (/'ndim'/)
            bl(9)%label           = 'Q_compression_cv'
            bl(9)%numrows         = UNINIT_INT
            bl(9)%numkw           = 0   
            bl(10)%label           = 'parameter_groups'
            bl(10)%numrows         = UNINIT_INT
            bl(10)%numkw           = 0   
       end subroutine bpi_init_param_CVs
       
!********  subroutine bpi_init_param_DATA
       subroutine bpi_init_param_DATA(BL,d_PAR)
       ! SUBROUTINE to initialze parameter-related block variables
         use bayes_pest_control
       ! DECLARATIONS  
         implicit none
         type(tp_block),   intent(inout)     :: BL(NUM_BLOCK)
         type(d_param),    intent(inout)       :: d_PAR 
       ! INITIALIZATIONS
            nullify(d_PAR%group)
            nullify(d_PAR%parnme)
            nullify(d_PAR%pars)
            nullify(d_PAR%pars_old)
            nullify(d_PAR%pars_lns)
            nullify(d_PAR%lox)
            nullify(d_PAR%SenMethod)
            nullify(d_PAR%BetaAssoc)
            nullify(d_PAR%Group_type)
            bl(11)%label           = 'parameter_data'
            bl(11)%numrows         = UNINIT_INT
            bl(11)%numkw           = 0 
            d_PAR%phi_T            = UNINIT_REAL
            d_PAR%phi_M            = UNINIT_REAL
            d_PAR%phi_R            = UNINIT_REAL
       end subroutine bpi_init_param_DATA 
               
                                   
!********  subroutine bpi_init_obs_CVs
       subroutine bpi_init_obs_groups(BL,cv_OBS)
       ! SUBROUTINE to initialze observation control variables
         use bayes_pest_control
       ! DECLARATIONS  
         implicit none
         type(cv_observ),intent(inout)   :: cv_OBS
         type(tp_block),intent(inout)    :: BL(NUM_BLOCK) 
       ! INITIALIZATIONS
            nullify(cv_OBS%grp_name)
            cv_OBS%nobsgp         = UNINIT_INT ! number of observation groups
            bl(12)%label           = 'observation_groups'
            bl(12)%numrows         = UNINIT_INT
            bl(12)%numkw           = 1 
            allocate (bl(12)%keywords(bl(12)%numkw))
            bl(12)%keywords = (/'nobsgp'/)
       end subroutine bpi_init_obs_groups
       
!********  subroutine bpi_init_obs_DATA
       subroutine bpi_init_obs_DATA(BL,d_OBS)
       ! SUBROUTINE to initialze algorithmic control variables
         use bayes_pest_control
       ! DECLARATIONS  
         implicit none
         type(d_observ),   intent(inout)   :: d_OBS
         type(tp_block),   intent(inout)   :: BL(NUM_BLOCK) 
       ! INITIALIZATIONS
            nullify(d_OBS%group)
            nullify(d_OBS%obs)
            nullify(d_OBS%h)
            nullify(d_OBS%weight)
            nullify(d_OBS%obsnme)
            bl(13)%label           = 'observation_data'
            bl(13)%numrows         = UNINIT_INT
            bl(13)%numkw           = 0 
       end subroutine bpi_init_obs_DATA 
               
!********  subroutine bpi_init_modcomlin_DATA
       subroutine bpi_init_modcomlin_DATA(BL,d_MOD)
       ! SUBROUTINE to initialze algorithmic control variables
         use bayes_pest_control
       ! DECLARATIONS  
         implicit none
         type(d_comlin), intent(inout)      :: d_MOD
         type(tp_block),intent(inout)    :: BL(NUM_BLOCK) 
       ! INITIALIZATIONS
            d_MOD%com = UNINIT_CHAR
            d_MOD%dercom = UNINIT_CHAR
            bl(14)%label           = 'model_command_lines'
            bl(14)%numrows         = UNINIT_INT
            bl(14)%numkw           = 2   
            allocate (bl(14)%keywords(bl(14)%numkw))
            bl(14)%keywords = (/'Command','DerivCommand'/)      
       end subroutine bpi_init_modcomlin_DATA
                                  
!********  subroutine bpi_init_mio_CVs
       subroutine bpi_init_mio_CVs(BL,cv_MIO)
       ! SUBROUTINE to initialze algorithmic control variables
         use bayes_pest_control
       ! DECLARATIONS  
         implicit none
         type(cv_minout),intent(inout)   :: cv_MIO
         type(tp_block),intent(inout)    :: BL(NUM_BLOCK) 
       ! INITIALIZATIONS
            cv_MIO%ninsfle         = UNINIT_INT ! number of instruction files
            cv_MIO%ntplfle         = UNINIT_INT ! number of template files
            bl(15)%label           = 'model_input_files'
            bl(15)%numrows         = UNINIT_INT
            bl(15)%numkw           = 0
            bl(16)%label           = 'model_output_files'
            bl(16)%numrows         = UNINIT_INT
            bl(16)%numkw           = 0  

       end subroutine bpi_init_mio_CVs
       
!********* subroutine bpi_init_anisotropy_DATA
        subroutine bpi_init_anisotropy_DATA(BL,d_ANI)
        
        use bayes_pest_control
        implicit none
        type(d_anisotropy),intent(inout)    :: d_ANI
        type(tp_block),intent(inout)        :: BL(NUM_BLOCK)    
       ! INITIALIZATIONS
            nullify(d_ANI%horiz_angle)
            nullify(d_ANI%horiz_ratio)
            nullify(d_ANI%vertical_ratio)
            nullify(d_ANI%BetaAssoc)
            BL(17)%label           = 'parameter_anisotropy'
            BL(17)%numrows         = UNINIT_INT
            BL(17)%numkw           = 0     
        end subroutine bpi_init_anisotropy_DATA
        
        
!*********  subroutine bpi_init_algorithmic_DATA(d_A,npar,nobs)
       subroutine bpi_init_algorithmic_DATA(d_A,npar,nobs)      
       !SUBROUTINE to initialize receptacle for the Jacobian (H)
       use bayes_pest_control
       !DECLARATIONS
       implicit none
       type(d_algorithmic), intent(inout)   :: d_A
       integer, intent(in)                  :: npar, nobs
       ! INITIALIZATIONS
       nullify(d_A%H)
       nullify(d_A%HX)
       nullify(d_A%HQHt)
       nullify(d_A%Hsold)
       nullify(d_A%Qsy)
       nullify(d_A%Qyy)
       nullify(d_A%beta_hat)
       nullify(d_A%ksi)
       allocate(d_A%H(nobs,npar))
       d_A%H    = UNINIT_REAL  ! matrix         
       end subroutine  bpi_init_algorithmic_DATA               
   
    end module bpi_initializations