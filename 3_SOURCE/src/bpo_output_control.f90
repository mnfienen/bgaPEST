module bayes_output_control
    ! ** DECLARATIONS **
       use jupiter_input_data_support
       use bayes_pest_control
       use bdp_data_parsers
       use bpi_initializations      
       use utilities
       use model_input_output
       use error_message
    
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!            subroutine to WRITE BANNER TO THE SCREEN                      !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine bpo_write_banner()
      write(6,*)
      write(6,*) '*******************************************'
      write(6,*) '*          Welcome to bgaPEST             *'
      write(6,*) '*******************************************' 
      write(6,*) '       Version 1.0  November 19, 2012     '
      write(6,*) '           Brought to you by:             ' 
      write(6,*) '     United States Geological Survey      '
      write(6,*) '      Watermark Numerical Computing       '
      write(6,*) '            University of Parma           '
      write(6,*) '*******************************************'
      write(6,*)
      
   end subroutine bpo_write_banner
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!            subroutine to TELL USER THAT IMPROPER COMMAND LINE GIVEN      !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine bpo_write_nocmd()     
      write(6,*) 'bgaPEST is run using the command :-'
      write(6,*) 'bgaPEST bgapestfile '
      write(6,*) 'where bgapestfile is a bgaPEST control file'
   end subroutine bpo_write_nocmd
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!            subroutine to WRITE ALL PARAMETER VALUES TO A BPP FILE        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine bpo_write_allpars(cv_PAR,d_PAR,writunit)
      type (cv_param)              :: cv_PAR
      type (d_param)               :: d_PAR
      integer, intent(in)          :: writunit
      character(50)                :: outlinefmt
      character(20)                :: parwstr,pargwstr
   
    call utl_int2char(PARNWIDTH,parwstr)
    call utl_int2char(PARGROUPNMWID,pargwstr) 
    write(outlinefmt,"('(1A',A,',1A',A,',1A13,1A16)')") trim(parwstr),trim(pargwstr)
    write(writunit,trim(outlinefmt)) 'ParamName','ParamGroup','BetaAssoc','ParamVal'
    do i = 1,cv_PAR%npar
        write(outlinefmt,"('(1A',A,',1A',A,',1I13,1E16.8)')")  trim(parwstr),trim(pargwstr)
        write(writunit,trim(outlinefmt)) trim(d_PAR%parnme(i)),trim(d_PAR%group(i)), d_PAR%BetaAssoc(i),d_PAR%pars(i)
    enddo
   
   end subroutine bpo_write_allpars
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!            subroutine to WRITE ALL RESIDUAL VALUES TO A BRE FILE         !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine bpo_write_residuals(cv_OBS,d_OBS, writunit)
      type (cv_observ)                :: cv_OBS
      type (d_observ)                 :: d_OBS
      integer, intent(in)             :: writunit
      character(50)                   :: outlinefmt

    write(outlinefmt,"('(1A',I,',1A',I,',1A16,1A16)')") OBSNWIDTH,OBSGROUPNMWID    
    write(writunit,trim(outlinefmt)) 'ObsName','ObsGroup','Modeled','Measured'
    do i = 1,cv_OBS%nobs
        write(outlinefmt,"('(1A',I,',1A',I,',1E16.8,1E16.8)')") OBSNWIDTH,OBSGROUPNMWID  
        write(writunit,trim(outlinefmt)) trim(d_OBS%obsnme(i)),trim(d_OBS%group(i)),d_OBS%h(i),d_OBS%obs(i)
    enddo
           
   end subroutine bpo_write_residuals
   
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!            subroutine to WRITE POSTERIOR COVARIANCE MATRIX OR VECTOR TO OUTPUT FILE        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine bpo_write_posterior_covariance(compress_flag,cv_PAR,d_PAR,d_PM,V,VV,writunit)
      type (cv_param)                 :: cv_PAR
      type (d_param)                  :: d_PAR
      type (d_prior_mean), intent(in) :: d_PM
      integer, intent(in)             :: compress_flag
      double precision, intent(in)    :: V(:)    ! diagonal of posterior covariance in the case of compression
      double precision, intent(in)    :: VV(:,:) ! full posterior covariance matrix if no compression
      integer, intent(in)             :: writunit
      integer                         :: i,j,k
   
     select case (compress_flag)
       case (0) ! no compression - full posterior covariance matrix
          write(writunit,"(3I10)") cv_PAR%npar,cv_PAR%npar,1
            k = 0
            do i = 1,cv_PAR%npar
              do j = 1,cv_PAR%npar
                k = k+1
                write(writunit,"(1ES18.8)",advance="no") VV(i,j)
                if (MOD(k,8) .eq. 0) then
                  write(writunit,*)
                endif
              enddo
            enddo
            write(writunit,*)
       case (1) ! compressed matrix - only diagonals from posterior covariance matrix
          write(writunit,"(3I10)") cv_PAR%npar,cv_PAR%npar,-1
            do i = 1,cv_PAR%npar
              write(writunit,"(1ES18.8)"),V(i)
            enddo
     end select
     write(writunit,*) '*row and column names'
     do i = 1,cv_PAR%npar
       write(writunit,*) d_PAR%parnme(i)
     enddo
     
   end subroutine bpo_write_posterior_covariance  
   
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!            subroutine to WRITE FINAL PARAMETER VALUES with CONFIDENCE INTERVALS TO A BPP FILE        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine bpo_write_allpars_95ci(cv_PAR,d_PAR,d_PM,V,writunit,ci95_flag)
      type (cv_param)                 :: cv_PAR
      type (d_param)                  :: d_PAR
      type (d_prior_mean), intent(in) :: d_PM
      double precision, intent(in)    :: V(:)
      integer, intent(in)             :: writunit,ci95_flag
      character(50)                   :: outlinefmt
      character(20)                   :: parwstr,pargwstr
      double precision, allocatable   :: lcl(:),ucl(:),finalparvalue(:)
      integer                         :: i,j
        
        call utl_int2char(PARNWIDTH,parwstr)
        call utl_int2char(PARGROUPNMWID,pargwstr) 
        
    if (ci95_flag .eq. 1) then  ! only calculate LCL and UCL if posterior covariance was calculated
      
      allocate(finalparvalue(cv_PAR%npar))
      finalparvalue = d_PAR%pars !-- these are the optimal values, still in physical space    
      !--- handle transformations to estimation space as appropriate
      if (maxval(d_PM%Partrans).ge.1) then  !If yes, the parameters transformation is required  
        do i = 1,cv_PAR%npar 
          if (d_PM%Partrans(d_PAR%BetaAssoc(i)).eq.1) then
            finalparvalue(i) = log(d_PAR%pars(i))
          elseif (d_PM%Partrans(d_PAR%BetaAssoc(i)).eq.2) then
            finalparvalue(i)=d_PM%alpha(d_PAR%BetaAssoc(i))*((d_PAR%pars(i)**(1./d_PM%alpha(d_PAR%BetaAssoc(i))))-1)   
          endif
        enddo   
      endif
        
      allocate(lcl(cv_PAR%npar))
      allocate(ucl(cv_PAR%npar))
       
      write(outlinefmt,"('(1A',A,',1A',A,',1A13,1A16,1A16,1A16)')") trim(parwstr),trim(pargwstr)
      write(writunit,trim(outlinefmt)) 'ParamName','ParamGroup','BetaAssoc','ParamVal','95pctLCL','95pctLCL'
   
        !-- calculate LCL, and UCL
        lcl = finalparvalue
        ucl = finalparvalue
        lcl = lcl - 2*sqrt(V)
        ucl = ucl + 2*sqrt(V)
        if (allocated(finalparvalue)) deallocate(finalparvalue)     
      !-- backtransform to physical space as appropriate
      if (maxval(d_PM%Partrans).ge.1) then  !If yes, the parameters back-transformation is required  
        do i = 1,cv_PAR%npar 
         if (d_PM%Partrans(d_PAR%BetaAssoc(i)).eq.1) then
              ucl(i) = exp(ucl(i))
              lcl(i) = exp(lcl(i))
         elseif (d_PM%Partrans(d_PAR%BetaAssoc(i)).eq.2) then
              ucl(i) = ((ucl(i)/d_PM%alpha(d_PAR%BetaAssoc(i)))+1)**(d_PM%alpha(d_PAR%BetaAssoc(i)))
              lcl(i) = ((lcl(i)/d_PM%alpha(d_PAR%BetaAssoc(i)))+1)**(d_PM%alpha(d_PAR%BetaAssoc(i)))
         endif
        enddo 
      endif      
       
        !--  write parameter values, LCL, and UCL
        do j = 1,cv_PAR%npar    
            write(outlinefmt,"('(1A',A,',1A',A,',1I13,1E16.8,1E16.8,1E16.8)')")  trim(parwstr),trim(pargwstr)
            write(writunit,trim(outlinefmt)) trim(d_PAR%parnme(j)),trim(d_PAR%group(j)), d_PAR%BetaAssoc(j),d_PAR%pars(j), &
                     lcl(j),ucl(j)
        enddo
        if (allocated(lcl)) deallocate(lcl)
        if (allocated(ucl)) deallocate(ucl)
        
    else ! -- write final parameter values without confidence intervals if no posterior covariance was calculated
        write(outlinefmt,"('(1A',A,',1A',A,',1A13,1A16)')") trim(parwstr),trim(pargwstr)
        write(writunit,trim(outlinefmt)) 'ParamName','ParamGroup','BetaAssoc','ParamVal'
                    
        !--  write parameter values
        do j = 1,cv_PAR%npar    
            write(outlinefmt,"('(1A',A,',1A',A,',1I13,1E16.8)')")  trim(parwstr),trim(pargwstr)
            write(writunit,trim(outlinefmt)) trim(d_PAR%parnme(j)),trim(d_PAR%group(j)), d_PAR%BetaAssoc(j),d_PAR%pars(j)
        enddo
        
     endif ! -- ci95_flag
        
       
         
   end subroutine bpo_write_allpars_95ci
      
   
   
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!        subroutine to WRITE INITIAL VALUES TO THE RECORD (BPR) FILE       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine bpo_write_bpr_header(bprunit,casename,cv_PAR,cv_OBS, &
                d_MOD,cv_A,cv_MIO,d_MIO,Q0_all,cv_PM,d_PM,cv_S,d_S,d_PAR,d_ANI)
   
    integer, intent(in)                 :: bprunit
    integer                             :: cparunit
    type(cv_param), intent(in)          :: cv_PAR
    type(cv_observ), intent(in)         :: cv_OBS
    type (cv_struct)                    :: cv_S 
    type (d_struct)                     :: d_S  
    type(cv_prior_mean), intent(in)     :: cv_PM  
    type(d_prior_mean), intent(in)      :: d_PM
    type (Q0_compr), intent(in)         :: Q0_All(:)
    type(d_comlin), intent(in)          :: d_MOD
    type(cv_algorithmic), intent(in)    :: cv_A
    type(d_minout), intent(in)          :: d_MIO
    type(cv_minout), intent(in)         :: cv_MIO
    type(d_param), intent(in)           :: d_PAR
    type(d_anisotropy), intent(in)      :: d_ANI
    character (len=4)                   :: indent = '    '
    character (len=FILEWIDTH)           :: casename
    integer i,j ! local counters
    
    write(bprunit,*)
    write(bprunit,*)
    write(bprunit,10) casename
10  format('          bgaPEST RUN RECORD: CASE ', 1A)
!!! dimensions
    write(bprunit,*) 'Case Dimensions :-'
    write(bprunit,15) 'Number of Parameters',cv_PAR%npar   
    write(bprunit,15) 'Number of Parameter Groups',cv_PAR%npargp 
    write(bprunit,15) 'Number of Observations',cv_OBS%nobs    
    write(bprunit,15) 'Number of Observation Groups',cv_OBS%nobsgp    
15 format(A40, ':', I8)
   
!!!  Model calls
   write(bprunit,*)
   write(bprunit,*) 'Model Command Line :-'
   write(bprunit,20)  indent,d_MOD%com
   write(bprunit,*) 'Jacobian Command Line :-'
   if (cv_A%deriv_mode .eq. 1) then
      write(bprunit,20) indent,d_MOD%dercom
   else
      write(bprunit,20) indent,'Not Used'
   endif

!!! model input and output
   write(bprunit,*)
   write(bprunit,*) 'Model Interface Files:-'
   write(bprunit,20) indent,'Template files:'
   do i = 1,cv_MIO%ntplfle
     write(bprunit,25)indent,indent,d_MIO%tpl(i)
   enddo
   write(bprunit,20) indent,'for model input files:'
   do i = 1,cv_MIO%ntplfle
     write(bprunit,25)indent,indent,d_MIO%infle(i)
   enddo
   write(bprunit,20) indent,'Instruction files:'
   do i = 1,cv_MIO%ninsfle
     write(bprunit,25)indent,indent,d_MIO%ins(i)
   enddo
   write(bprunit,20) indent,'for model ouput files:'
   do i = 1,cv_MIO%ninsfle
     write(bprunit,25)indent,indent,d_MIO%outfle(i)
   enddo
   if (cv_A%deriv_mode .eq. 1) then
     write(bprunit,20) indent,'Jacobian Matrix Output File:'
     write(bprunit,25) indent,indent,cv_A%jacfle
     write(bprunit,20) indent,'Jacobian Matrix Output Format:'
     write(bprunit,25) indent,indent,cv_A%jacobian_format
     
   endif
20 format(2A )     ! single indent and str format
25 format(3A)   ! double indent and str format
30 format(1A, 1ES12.4)     ! single indent and ES format
35 format(3A,1ES12.4 )   ! double indent and ES format
40 format(2A, 1I10)     ! single indent and integer format
45 format(3A,1I10)   ! double indent and integer format
50 format(2A 1L5)     ! single indent and logical format
55 format(3A,1L5)   ! double indent and logical format
   
!!! Algorithmic control variables
    write(bprunit,*)
    write(bprunit,*) 'Algorithmic control variables:-'
    write(bprunit,20) indent,'Structural Paramter Convergence'
    write(bprunit,35) indent,indent,'structural_conv: ',cv_A%structural_conv
    write(bprunit,20) indent,'Objective Function Convergence'
    write(bprunit,35) indent,indent,'phi_conv: ',cv_A%phi_conv
    write(bprunit,20) indent,'Outermost BGA Convergence'
    write(bprunit,35) indent,indent,'bga_conv: ',cv_A%bga_conv    
    write(bprunit,20) indent,'Maximum Number of Structural Paramter Iterations'
    write(bprunit,45) indent,indent,'it_max_structural: ',cv_A%it_max_structural       
    write(bprunit,20) indent,'Maximum Number of Objective Function Iterations'
    write(bprunit,45) indent,indent,'it_max_phi: ',cv_A%it_max_phi       
    write(bprunit,20) indent,'Maximum Number of Outermost BGA Iterations'
    write(bprunit,45) indent,indent,'it_max_bga: ',cv_A%it_max_bga       
    write(bprunit,20) indent,'Linesearch Flag: [0] indicates no linesearch, [1] indicates perform linesarch'
    write(bprunit,45) indent,indent,'lns_flag: ',cv_A%lns_flag       
    write(bprunit,20) indent,'Maximum Number of Linesearch Iterations'
    write(bprunit,45) indent,indent,'it_max_lns: ',cv_A%it_max_lns
    write(bprunit,20) indent,'Logical Flag for storing Q Matrix'
    write(bprunit,55) indent,indent,'store_Q: ',cv_A%store_Q
    write(bprunit,20) indent,'Form of theta covariance: [0] none, [1], diag, [2] full matrix'
    write(bprunit,45) indent,indent,'theta_cov_form: ',cv_A%theta_cov_form
    write(bprunit,20) indent,'Compression of Q0 matrix: [0] none - compute full Q0, [1], Q0 for each beta'
    write(bprunit,45) indent,indent,'Q_ compression_flag: ',cv_A%Q_compression_flag
    write(bprunit,20) indent,'Derivatives mode: [0] External PEST Perturbations, [1] specified Jacobian command line'
    write(bprunit,45) indent,indent,'deriv_mode: ',cv_A%deriv_mode
    if (cv_A%deriv_mode .eq. 1) then
        write(bprunit,20) indent,'External derivatives calculated using file:-'
        write(bprunit,25) indent, indent, d_MOD%dercom
        write(bprunit,20) indent,'External derivatives read from file:-'
        write(bprunit,25) indent, indent, cv_A%jacfle
        write(bprunit,20) indent,'External derivatives file format:-'
        write(bprunit,25) indent, indent, cv_A%jacobian_format
    endif
    
    write(bprunit,*)        
    write(bprunit,20) indent,'Parameter Anisotropy: [0] do not read parameter_anisotropy block, [1] do read block'
    write(bprunit,45) indent,indent,'par_anisotropy: ',cv_A%par_anisotropy
    
!!! Beta associations (facies associations)
    write(bprunit,*)
    select case (cv_PAR%p)
        case (1)
            write(bprunit,56) indent,cv_PAR%p
        case default
            write(bprunit,57) indent,cv_PAR%p
    end select
56  format(1A,I4, ' Beta Association was defined:')
57  format(1A,I4, ' Beta Associations were defined:')

    if (cv_A%Q_compression_flag .ne. 0)  then
        ! write out each beta association's details
        do i = 1,cv_PAR%p
            write(bprunit,60) Q0_all(i)%BetaAss
            write(bprunit,20) indent,'Estimation space for this Beta Association: [none] physical space, [log] log transf., [power] power transf.'
            select case (d_PM%Partrans(i))
              case (0)
                write(bprunit,25) indent,indent,'none - Estimation is in the physical space'
              case (1)  
                write(bprunit,25) indent,indent,'log - Estimation is in a log transformed space'
              case (2)
                write(bprunit,74) indent,indent,'power - Estimation is in a power transformed space',' with alpha ', d_PM%alpha(i)
            end select
            write(bprunit,20) indent,'Parameter number at which this Beta Association starts'
            write(bprunit,45) indent,indent,'Beta_start: ',Q0_all(i)%Beta_start
            write(bprunit,20) indent,'Toeplitz Flag'
            write(bprunit,45) indent,indent,'Toep_flag: ',Q0_all(i)%Toep_flag
            write(bprunit,20) indent,'Number of Rows'
            write(bprunit,45) indent,indent,'Nrow: ',Q0_all(i)%Nrow
            write(bprunit,20) indent,'Number of Columns'
            write(bprunit,45) indent,indent,'Ncol: ',Q0_all(i)%Ncol
            write(bprunit,20) indent,'Number of Layers'
            write(bprunit,45) indent,indent,'Nlay: ',Q0_all(i)%Nlay
            write(bprunit,20) indent,'Number of Parameters'
            write(bprunit,45) indent,indent,'Npar: ',Q0_all(i)%Npar
        enddo
    else
        !indicate no compression specified
      do i = 1,cv_PAR%p
        write(bprunit,60) i
            write(bprunit,20) indent,'Estimation space for this Beta Association: [none] physical space, [log] log transf., [power] power transf.'
            select case (d_PM%Partrans(i))
              case (0)
                write(bprunit,25) indent,indent,'none - Estimation is in the physical space'
              case (1)  
                write(bprunit,25) indent,indent,'log - Estimation is in a log transformed space'
              case (2)
                write(bprunit,74) indent,indent,'power - Estimation is in a power transformed space',' with alpha ', d_PM%alpha(i)
            end select
       enddo
            write(bprunit,*)
            write(bprunit,65) indent,'No Compression Requested.'
    endif    
60 format(' Variables for Beta Association: ',I3)
65 format(2A)
   
!!! Derivatives Calculations
    write(bprunit,*)
    !!! placeholder here, in case we implement group-specific derivatives  

!!! Prior Means Information if Supplied
    write(bprunit,*)
    if (cv_PM%betas_flag .eq. 1) then
        write(bprunit,*) 'Prior Information on Betas:-'

        write(bprunit,20) indent,'Prior information format on prior menas (betas)'
        write(bprunit,20) indent,'[0] none, [1] diagonal only, [2] full covariance matrix'
        write(bprunit,45) indent,indent,'Qbb_form: ', cv_PM%Qbb_form    
        do i = 1,cv_PAR%p
            write(bprunit,70) indent, i
            write(bprunit,71) indent,indent,d_PM%beta_0(i)
        enddo ! i = 1,cv_PAR%p
    else
        write(bprunit,20) indent,'No Prior Information Provided for Betas'
    endif
70 format(1A, 'Starting value of Prior Mean for Beta Association: ', I4)
71 format(1A,1A, 1ES12.4)
    


!!! Structural Parameter Definitions
    write(bprunit,*)
    write(bprunit,*) 'Structural parameter control variables:-'

    
    do i = 1,cv_PAR%p
        !-- indicate the variogram type and initial structural parameter values
        write(bprunit,72) indent, i
        select case (cv_S%var_type(i))
         case (0)
            write(bprunit,73) indent, indent, 'nugget'
            write(bprunit,25) indent, indent, 'Initial structural parameter values:'
            write(bprunit,74) indent,  indent,indent, 'nugget variance', d_S%theta_0(i,1)
         case (1)
            write(bprunit,73) indent, indent, 'linear'
            write(bprunit,25) indent, indent, 'Initial structural parameter values:'
            write(bprunit,74) indent,  indent,indent, 'slope', d_S%theta_0(i,1)
         case (2)
            write(bprunit,73) indent, indent, 'exponential'
            write(bprunit,25) indent, indent, 'Initial structural parameter values:'
            write(bprunit,74) indent, indent, indent, 'variance', d_S%theta_0(i,1)
            write(bprunit,74) indent,  indent,indent, 'correlation length', d_S%theta_0(i,2)
        end select
        !-- indicate whether structural parameters are to be optimized for
            write(bprunit,25) indent,indent,'Flag determining whether structural parameters are optimized for or fixed:'
            write(bprunit,75) indent,indent,indent,'struct_par_opt: ',cv_S%struct_par_opt(i)
            if (cv_S%struct_par_opt(i) .eq. 1) then
               if (cv_S%trans_theta(i) .eq. 1) then
                   write(bprunit,25) indent, indent, 'Power transformation active'
                   write(bprunit,74)indent, indent, indent, 'alpha', cv_S%alpha_trans(i)
               endif
            endif
            
    enddo ! i = 1,cv_PAR%p
    write(bprunit,'(3A)')indent,'','Epistemic Uncertainty :-'
    write(bprunit,76) indent, indent, d_S%sig
    ! -- indicate whether epistemic uncertainty is being optimized for 
    if (d_S%sig_opt .eq. 1) then
        write(bprunit,25) indent,indent,'Flag determining whether structural parameters are optimized for or fixed:'
        write(bprunit,75) indent,indent,indent,'sig_opt: ',d_S%sig_opt
        if (d_S%trans_sig .eq. 1) then
            write(bprunit,25) indent, indent, 'Power transformation active'
            write(bprunit,74)indent, indent, indent, 'alpha', d_S%alpha_trans_sig
        endif
    endif
    
72 format(1A, 'Structural Parameters for Beta Association: ', I4)    
73 format(1A, 1A,'Variogram type: ', 1A)
74 format(4A, ' = ', 1ES13.6)
75 format(4A,1I4)
76 format(1A, 1A, 'sigma (epistemic) = ',1ES13.6)

!!! Anisotropy parameters

    if (cv_A%par_anisotropy .ne. 0) then
        write(bprunit,*)
        write(bprunit,*) 'Parameter Anisotropy variables:-' 
        do i = 1,cv_PAR%p
            !-- indicate the variogram type and initial structural parameter values
            write(bprunit,72) indent, i
            do j = 1,cv_PAR%p
                if (d_ANI%BetaAssoc(j) .eq. i) then
                ! found the right BetaAssoc, now write out the details
                write(bprunit,78) indent, indent, d_ANI%horiz_angle(i)
                write(bprunit,79) indent, indent, d_ANI%horiz_ratio(i)
                write(bprunit,80) indent, indent, d_ANI%vertical_ratio(i)
                endif
                enddo !j
        enddo !
    endif    
77 format(1A, 'Anisotropy Variables for Beta Association: ', I4) 
78 format(2A, 'horiz_angle = ', 1ES13.6)   
79 format(2A, 'horiz_ratio = ', 1ES13.6)
80 format(2A, 'vertical_ratio = ', 1ES13.6)

   
    
!!! Parameter definitions
    write(bprunit,*)
    write(bprunit,*) 'Initial Parameter Definitions:-'
    cparunit = utl_nextunit()
    !-- write out the initial parameter values as bpp.0
    call bpc_openfile(cparunit,trim(trim(casename) // '.bpp.0'),1) ![1] at end indicates open with write access
    call bpo_write_allpars(cv_PAR,d_PAR,cparunit)
    close(cparunit)
    write(bprunit,85) indent,'For initial parameter values see file :- ', trim(trim(casename) // '.bpp.0')
85  format(3A) 
    end subroutine bpo_write_bpr_header







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!        subroutine to WRITE INTERMEDIATE INFORMATION TO BOTH RECORD (BPR) FILE and TO STANDARD OUT     !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine bpo_write_bpr_intermed(bprunit,inner_iter,outer_iter,curr_par_file,curr_resid_file,d_PAR)
   
    integer, intent(in)                 :: bprunit
    integer                             :: inner_iter,outer_iter
    type(d_param), intent(in)           :: d_PAR
    character (len=4)                   :: indent = '    '
    character (len=100)                 :: curr_par_file,curr_resid_file 
    integer i,j ! local counters
    
    write(bprunit,*)
    write(bprunit,*)
    write(bprunit,*) '***********************************'
    write(bprunit,40) indent,'Outer Iteration Number: ',outer_iter
    write(bprunit,40) indent,'Inner Iteration Number: ',inner_iter
    write(bprunit,*) 
    write(bprunit,35) indent,'Total PHI : ',d_PAR%phi_T
    write(bprunit,35) indent, 'Misfit PHI : ',d_PAR%phi_M
    write(bprunit,35) indent, 'Regularization PHI : ',d_PAR%phi_R
    write(bprunit,*) 
    write(bprunit,90) indent, indent, 'For current parameter values see file :- ', curr_par_file
    write(bprunit,90) indent, indent, 'For current residual values see file :- ', curr_resid_file
30 format(1A, 1ES12.4)     ! single indent and ES format
35 format(1A,1A30,1ES12.4 )   ! double indent and ES format
36 format(3A,1ES12.4 )   ! double indent and ES format
40 format(2A, 1I10)     ! single indent and integer format
45 format(3A,1I10)   ! double indent and integer format
90  format(4A)    
    end subroutine bpo_write_bpr_intermed
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!        subroutine to WRITE INTERMEDIATE Structural Parameter INFORMATION TO BOTH RECORD (BPR) FILE and TO STANDARD OUT     !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine bpo_write_bpr_intermed_structpar(bprunit,cv_S,d_S,cv_PAR,d_PAR)
        type(d_param), intent(in)           :: d_PAR
        type (cv_param), intent(in)         :: cv_PAR
        type (cv_struct), intent(in)        :: cv_S 
        type (d_struct), intent(in)         :: d_S 
        integer, intent(in)                 :: bprunit
        integer                             :: i,j
        character (len=4)                   :: indent = '    '
    
!!! Current Structural Parameter Values
    write(bprunit,*)
    write(bprunit,200) indent, '***STRUCTURAL PARAMETERS FOR WHICH OPTIMIZATION WAS REQUESTED***'
    do i = 1,cv_PAR%p
        !-- indicate the variogram type and current structural parameter values
        
      if (cv_S%struct_par_opt(i) .eq. 1) then
        write(bprunit,172) indent, i
        select case (cv_S%var_type(i))
         case (0)
          
            write(bprunit,173) indent, indent, 'nugget'
!            write(bprunit,25) indent, indent, 'Initial structural parameter values:'
            write(bprunit,174) indent,  indent,indent, 'nugget variance', d_S%theta(i,1)
         case (1)
            write(bprunit,173) indent, indent, 'linear'
!            write(bprunit,25) indent, indent, 'Initial structural parameter values:'
            write(bprunit,174) indent,  indent,indent, 'slope', d_S%theta(i,1)
         case (2)
            write(bprunit,173) indent, indent, 'exponential'
!            write(bprunit,25) indent, indent, 'Initial structural parameter values:'
            write(bprunit,174) indent, indent, indent, 'variance', d_S%theta(i,1)
            write(bprunit,174) indent,  indent,indent, 'correlation length', d_S%theta(i,2)
        end select
      endif          
    enddo ! i = 1,cv_PAR%p
    if (d_S%sig_opt .eq. 1) then
        write(bprunit,'(3A)')indent,'','Epistemic Uncertainty :-'
        write(bprunit,176) indent, indent, d_S%sig
    endif
172 format(1A, 'Current Structural Parameters for Beta Association: ', I4)    
173 format(1A, 1A,'Variogram type: ', 1A)
174 format(4A, ' = ', 1ES13.6)
175 format(4A,1I4)
176 format(1A, 1A, 'sigma (epistemic) = ',1ES13.6)
200 format(2A )     ! single indent and str format
end subroutine bpo_write_bpr_intermed_structpar
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!        subroutine to WRITE FINAL Structural Parameter INFORMATION TO BOTH RECORD (BPR) FILE and TO STANDARD OUT     !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine bpo_write_bpr_final_structpar(bprunit,cv_S,d_S,cv_PAR,d_PAR)
        type(d_param), intent(in)           :: d_PAR
        type (cv_param), intent(in)         :: cv_PAR
        type (cv_struct), intent(in)        :: cv_S 
        type (d_struct), intent(in)         :: d_S 
        integer, intent(in)                 :: bprunit
        integer                             :: i,j
        character (len=4)                   :: indent = '    '
    
    write(bprunit,*)
    write(bprunit,300) indent, '***FINAL STRUCTURAL PARAMETERS***'
    do i = 1,cv_PAR%p
        !-- indicate the variogram type and final structural parameter values
        write(bprunit,272) indent, i
        select case (cv_S%var_type(i))
         case (0)
            write(bprunit,273) indent, indent, 'nugget'
            write(bprunit,274) indent,  indent,indent, 'nugget variance', d_S%theta(i,1)
         case (1)
            write(bprunit,273) indent, indent, 'linear'
            write(bprunit,274) indent,  indent,indent, 'slope', d_S%theta(i,1)
         case (2)
            write(bprunit,273) indent, indent, 'exponential'
            write(bprunit,274) indent, indent, indent, 'variance', d_S%theta(i,1)
            write(bprunit,274) indent,  indent,indent, 'correlation length', d_S%theta(i,2)
        end select          
    enddo ! i = 1,cv_PAR%p
      write(bprunit,'(3A)')indent,'','Epistemic Uncertainty :-'
      write(bprunit,276) indent, indent, d_S%sig
272 format(1A, 'Final Structural Parameters for Beta Association: ', I4)    
273 format(1A, 1A,'Variogram type: ', 1A)
274 format(4A, ' = ', 1ES13.6)
275 format(4A,1I4)
276 format(1A, 1A, 'Final sigma (epistemic) = ',1ES13.6)
300 format(2A )     ! single indent and str format
end subroutine bpo_write_bpr_final_structpar    
    
end module bayes_output_control