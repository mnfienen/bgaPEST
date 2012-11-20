module bayes_pest_control
      
      
!     Module of data types to be used by the bayes pest module
!     Initializations to default values take place here as well.
!     a m!ke@usgs joint
!           
!     mike fienen - mnfienen@usgs.gov  
!     V 0.0  3/5/08  
!     Modified by M.D. 22/9/09 
!     Revised by Mike Fienen, 8/2/2011
!     Revised by M.D., 09/30/2011
!   ***************        
!   * DECLARATIONS
!   ***************
       
      implicit none
!     --general variables
      integer, parameter   	       :: UNINIT_INT=-999999 !default value for integer variables remaining uninitialized
      double precision, parameter  :: UNINIT_REAL=9.999999999d9 !default value for real variables remaining uninitialized
      character (len=12),parameter :: UNINIT_CHAR='*UNIN*'   !default value for char variables remaining uninitialized
      integer, parameter   :: NUM_BLOCK = 17                 !Is the maximum number of blocks  
      integer, parameter   :: PARNWIDTH = 12                 ! width of a parameter name  MD
      integer, parameter   :: OBSNWIDTH = 12                 ! width of a observation name MD
      integer, parameter   :: PARGROUPNMWID = 16             ! width of a parameter group name
      integer, parameter   :: OBSGROUPNMWID = 16             ! width of a parameter group name
      

!     -- types

      type  :: tp_block        !group of keywords arrays with labels and numrows vect
            character (len=50)   		    :: label
            integer              		    :: numrows
            integer              		    :: numkw
            character (len=50), pointer     :: keywords(:)
            character (len=50), pointer     :: keywordstring(:)  
            integer,            pointer     :: keywordline(:)          
      end type tp_block

   
      type  :: cv_algorithmic   ! ALGORITHMIC CONTROL VARIABLES
            double precision	::	structural_conv    !MD Structural parameter convergence values
            double precision	::	phi_conv           !MD Objective function convergence value
            double precision    ::  bga_conv           ! overall objective function convergence (outer loop)
            integer		        ::	it_max_structural  !MD Max number of iterations for struct parameters
            integer		        ::	it_max_phi         !MD Max number of iterations for objective function
            integer             ::  it_max_bga         !MD Max number of iterations for geostatistical method
            integer             ::  lns_flag           !MD Linesearch procedure flag: [0] not perform [1] perform 
            integer             ::  it_max_lns         !MD Max number of iterations for linesearch procedure
            logical             ::  store_Q            !MD TRUE --> Store Q FALSE --> Not store Q --> We need to address this option
            integer		        ::	theta_cov_form     !MD  Form of theta covariance:  [0] none, [1] diag, [2] full matrix
            integer             ::  Q_compression_flag !MD  [0] none - calculate full Q0, [1] Calculate Q0 for each beta separately 
                                                       !and if nugget store just 1, if toep_flag store just a vector instead of the matrix
            integer             ::  post_cov_flag      ! [0] means do not calculate posterior covariance [1] calculate posterior covariance   
            integer             ::  deriv_mode         ! [0] means use perturbations within PEST, [1] means separate command line will be
                                                       ! provided for running derivatives.  Typically with Adjoint State. [4] means external, parallel
                                                       ! derivatives calculation will be carried out using Condor    
            character (len=6)   ::  jacobian_format    ! 'binary' for binary Jacobian file [default] or 'ascii' for standard PEST matrix text format
                                                       ! only read if deriv_mode =1       
            character(len=100)  ::  jacfle             !jacobian file - read if cv_A%deriv_mode==1, default='scratch.jco'
            integer             ::  par_anisotropy     ! flag for whether block of  parameter anistropy values will be read [0] = no, [1] = yess
                                    
      end type cv_algorithmic
      
      type  :: d_algorithmic    ! ALGORITHMIC "GLOBALS"
            double precision, pointer  :: H(:,:)        !MD Sensitivity matrix
            double precision, pointer  :: HX(:,:)
            double precision, pointer  :: HQHt(:,:) 
            double precision, pointer  :: Hsold(:)
            double precision, pointer  :: Qsy(:,:)      !MD QHt is the cross covariance between s and y
            double precision, pointer  :: Qyy(:,:)      !MD HQHt + R (Auto-covariance matrix of the observ y) 
            double precision, pointer  :: beta_hat(:)   !MD ESTIMATED MEAN
            double precision, pointer  :: ksi(:)        
      end type d_algorithmic
            
      type  :: cv_prior_mean    ! PRIOR MEANS CONTROL VARIABLES
            integer ::  betas_flag   !MD Have or not prior informations about mean?
            integer	::	Qbb_form     !MD Form of Beta covariance:  [0] none, [1] diag, [2] full matrix
      end type cv_prior_mean
      
      type  :: d_prior_mean     ! data for prior means - pointer
            double precision,pointer	::	beta_0(:)      !MD Prior beta values
            double precision,pointer	::	Qbb(:,:)       !MD Covariance of beta
            double precision,pointer    ::  InvQbb(:,:)    !MD Inverse of covariance of beta
            double precision,pointer    ::  InvQbbB0(:)    !MD Inverse of covariance of beta * beta0
            double precision,pointer    ::  alpha(:)       !MD Exponent of power transformation
            integer, pointer            ::  Partrans(:)    !MD Vector of parameter transformation  [1] Log [0] None
      end type d_prior_mean
      
      type  :: cv_struct        ! CONTROL VARIABLES FOR STRUCTURAL PARAMETERS
            integer, pointer		    ::	prior_cov_mode(:) !MD Supplied matrix or calculated
            integer, pointer		    ::	var_type(:)  !MD Type of variogram [0] pure nugget, [1] linear, [2] exponential
            integer, pointer		    ::	struct_par_opt(:)
            integer, pointer	    	::	num_theta_type(:)
            integer, pointer	    	::	trans_theta(:)
            double precision, pointer	::	alpha_trans(:)
            double precision            ::  str_obj_fun = 0.   ! Objective function for structural parameters
            integer                     ::  num_theta_opt = 0 ! num_theta_opt is number of theta pars to optimize
      end type cv_struct
      
      type d_struct             ! DATA for structural parameters - pointer
            double precision, pointer   ::	theta_0(:,:)  !MD Initial value of theta Matrix
            double precision, pointer   ::	struct_par_opt_vec_0(:)  !Initial value of all structural parameter values to be optimized - single vector (may include sigma)
            double precision, pointer   ::	struct_par_opt_vec(:)  !Current value of all structural parameter values to be optimized - single vector
            double precision, pointer   ::	theta_cov(:,:)!MD theta covariance matrix
            double precision, pointer   ::  invQtheta(:,:)   ! prior covariance matrix for all  theta optimized pars. This may include sigma
            double precision            ::  sig_0         !MD Initial value of sigma (epistemic uncertainty parameter)
            double precision            ::  sig_p_var     !MD Variance of sigma (variance of the epistemic error) 
            integer                     ::  sig_opt       !MD Added to allow the choice to optimize or not for sig
            integer                     ::  trans_sig     !Power transform flag
            double precision            ::  alpha_trans_sig !Exponent of power transformation for sigma
            double precision, pointer   ::  theta(:,:)    !MD Structural parameters matrix
            double precision            ::  sig           !MD epistemic uncertainty parameter
      end type d_struct
      
      type  :: cv_param          ! control variables for parameters
            integer                         :: npargp  !MD Number of parameter groups 
            integer                         :: npar    !MD Total number of parameters 
            character (len=50), pointer     :: grp_name(:) !MD Name of the parameter groups
            integer, pointer                :: grp_type(:) !MD Type of groups 
            double precision, pointer       :: derinc(:) !MNF derivative increment for group
            integer                         :: ndim  !MD Spatial dimensions
            integer                         :: p     !MD Number of means
      end type cv_param
      
      type  :: Q0_compr   ! control variable type for compression of Q (Toeplitz or not) and Q0
            integer                         :: BetaAss
            integer                         :: Toep_flag
            integer                         :: Nrow
            integer                         :: Ncol
            integer                         :: Nlay
            integer                         :: Npar
            integer                         :: Beta_Start
            double precision, pointer       :: Q0_C(:,:)
      end type Q0_compr            
            
      type  :: d_param           ! data type for parameters
            character (len=PARGROUPNMWID), pointer	    :: group(:)      !MD Name of group
            character (len=PARNWIDTH), pointer     :: parnme(:)     ! parameter name
            double precision, pointer       :: pars(:)       !MD Vector of parameters
            double precision, pointer       :: pars_old(:)   !MD Sold
            double precision, pointer       :: pars_lns(:)   !MD Vector of parameters for linesearch
            double precision, pointer       :: lox(:,:)      !MD Location (coordinates)
            double precision, dimension(1)  :: phi_T         !MD Objective function Total
            double precision, dimension(1)  :: phi_M         !MD Objective function Misfit
            double precision, dimension(1)  :: phi_R         !MD Objective function Regularization
            integer, pointer                :: SenMethod(:)  !MD Sensitivity calculation method
            integer, pointer                :: BetaAssoc(:)  !MD Faces association
            integer, pointer                :: Group_type(:) !MD Vector of group type for each parameter
      end type d_param
      
      type  :: cv_observ          ! control variables for observations
            integer             :: nobsgp  !MD Number of observations groups
            integer             :: nobs    !MD Number of observations
            character (len=PARGROUPNMWID), pointer     :: grp_name(:) !MD Name of the observations groups
      end type cv_observ

      type  :: d_observ           ! data type for observations
            character (len=OBSGROUPNMWID), pointer	:: group(:) !MD Name of group
            double precision, pointer   :: obs(:)   !MD Vector of observations
            double precision, pointer   :: h(:)     !MD Current model output
            double precision, pointer   :: weight(:) !MD Weight for R matrix 
            character (len=OBSNWIDTH), pointer :: obsnme(:) ! Observation names       
      end type d_observ
        
      type  :: d_comlin        ! data type FOR COMMAND LINE ARGS
            character (len=50)	:: com !MD Command line
            character (len=50)	:: dercom !Derivatives Command line
      end type d_comlin
      
      type  :: cv_minout          ! control variables for model i/o
            integer	::	ninsfle  !MD number of instruction files
            integer	::	ntplfle  !MD number of template files
      end type cv_minout
      
      type  :: d_minout          ! data type for model i/o
            character(len=100),pointer	::	tpl(:) !MD Template file
            character(len=100),pointer	::	infle(:) !MD Input file
            character(len=100),pointer	::	ins(:) !MD Instruction file
            character(len=100),pointer	::	outfle(:) !MD Output file
            character(len=100), pointer ::  pargroup(:) ! parameter groups for parallel run setup
      end type d_minout 
      
      type  :: kernel_XQR        ! kernels of X, Q, and R
            double precision, pointer   :: X(:,:)  !MD Deterministic base functions
            double precision, pointer   :: Q0(:,:) !MD Prior covariance
            double precision, pointer   :: R0(:,:) !MD Covariance matrix of epistemic error R=sig*R0
            double precision            :: L       !MD 10 times maximum distance in Q0 matrix
      end type kernel_XQR
      
      
      type :: d_anisotropy      ! parameter anisotropy information
           double precision, pointer    :: horiz_angle(:)    ! angle, from horizontal (X) of principal variance
           double precision, pointer    :: horiz_ratio(:)    ! ratio of maximum to minimum variance
           double precision, pointer    :: vertical_ratio(:) ! ratio of horizontal to vertical variance if 3D
           integer, pointer             :: BetaAssoc(:)      ! Beta association identifier
      end type d_anisotropy
      
!   *******************************        
!   * SUBROUTINES - ALL ARE VISIBLE
!   *******************************   
      
      
contains
      subroutine bpc_openfile(unitnum,fname,RW_FLAG)
      
! -- Suroutine BPC_OPENFILE opens a file of given name and 
!       unit.  Bombs if fails
     use  utilities
        implicit none
        integer          , intent(in)   :: unitnum
        character (len=*), intent(in)   :: fname
        character (len=200)             :: errmsg
        character (len=20)              :: inout   !MD Where do you use this???
        integer                         :: ierr
        integer, intent(in)             :: RW_FLAG !MD 0 Read - 1 Write flag
        select case (RW_FLAG)
          case (0) ! assumed for input
            open(unit=unitnum,file=trim(fname),status='old',iostat=ierr)
            inout = 'input'
          case (1) ! assumed for output
            open(unit=unitnum,file=trim(fname),status='replace',action='write',iostat=ierr)
            inout = 'output'
          case default
            errmsg = 'Programming error:-> RW_FLAG can only be 0 or 1'
            call utl_writmess(6,errmsg)
            stop
          end select 
       if(ierr.ne.0)then
         write(errmsg,100) trim(inout),trim(fname)
100      format(' Cannot open ',a,' file "',a,'".')
         call utl_writmess(6,errmsg)
         stop
       endif
      end subroutine bpc_openfile
     

      
!********  Subroutine INTREAD reads an integer from a string.
        SUBROUTINE INTREAD(IFAIL,CLINE,ITEMP)
        INTEGER, intent (out)           :: IFAIL
        INTEGER, intent (out)           :: ITEMP
        CHARACTER (len=6)               :: AFMT
        CHARACTER*(*), intent(in)       :: CLINE
        IFAIL=0
        AFMT='(i   )'
        WRITE(AFMT(3:5),'(I3)') LEN(CLINE)
        READ(CLINE,AFMT,ERR=100) ITEMP
        RETURN
100     IFAIL=1
        RETURN
        end subroutine intread
        
!********  Subroutine REALREAD reads an integer from a string. MD integer or real ???? if real why f .0 ??
        SUBROUTINE SREALREAD(IFAIL,CLINE,RTEMP)
        INTEGER, intent(out)          :: IFAIL
        REAL, intent (out)            :: RTEMP
        CHARACTER (len=8)        	  :: AFMT
        CHARACTER*(*), intent(in)	  :: CLINE
        IFAIL=0
        AFMT='(F   .0)'
        WRITE(AFMT(3:5),'(I3)') LEN_TRIM(CLINE)
        READ(CLINE,AFMT,ERR=100) RTEMP
        RETURN

100     IFAIL=1
        RETURN
  end subroutine srealread
  
  !********  Subroutine DREALREAD reads an integer from a string.
        SUBROUTINE DREALREAD(IFAIL,CLINE,RTEMP)
        INTEGER, intent(out)            		:: IFAIL
        double precision, intent (out)          :: RTEMP
        CHARACTER (len=8)        			    :: AFMT
        CHARACTER*(*), intent(in)			    :: CLINE
        IFAIL=0
        AFMT='(F   .0)'
        WRITE(AFMT(3:5),'(I3)') LEN_TRIM(CLINE)
        READ(CLINE,AFMT,ERR=100) RTEMP
        RETURN

100     IFAIL=1
        RETURN
  end subroutine drealread
    
!********  Subroutine  bpc_casetrans converts a string to upper or lower case.
     ! adapted from John Doherty's mio_casetrans code
     subroutine bpc_casetrans(string,hi_or_lo)
     implicit none

	 character (len=*), intent(inout)        :: string
	 character (len=*), intent(in)           :: hi_or_lo
	 character                               :: alo, ahi
	 integer                                 :: inc,i

	 if(hi_or_lo.eq.'lo') then
	  alo='A'; ahi='Z'; inc=iachar('a')-iachar('A')
	 else if(hi_or_lo.eq.'hi') then
	  alo='a'; ahi='z'; inc=iachar('A')-iachar('a')
	 else
      write(6,*) ' *** Illegal call to subroutine CASETRANS ***'
      stop
	 endif
     do i=1,len_trim(string)
	   if((string(i:i).ge.alo).and.(string(i:i).le.ahi)) &
 	   string(i:i)=achar(iachar(string(i:i))+inc)
 	 enddo
	return
    end subroutine bpc_casetrans
    
end module bayes_pest_control