module bayes_matrix_operations

!***** Modified by M.D. 24/9/09 ********
!***** Modified MNF 8/3/2011 ***********

        use bayes_pest_control
        use utilities  
        use objective_function


contains

!*****************************************************************************************************
!***** Perform generic matrix operations ************************************************************* 
!***************************************************************************************************** 

subroutine bmo_form_Qss_Qsy_HQsy(d_XQR, theta, cv_PAR, cv_OBS, cv_S, cv_A, d_A, d_PAR,Q0_All)
        
        implicit none
        ! declarations
        type(kernel_XQR),    intent(in)     :: d_XQR
        type(cv_struct),     intent(in)     :: cv_S 
        type(cv_param),      intent(in)     :: cv_PAR
        type(cv_observ),     intent(in)     :: cv_OBS
        type(cv_algorithmic),intent(inout)  :: cv_A        
        type(d_algorithmic), intent(inout)  :: d_A
        type(d_param),       intent(inout)  :: d_PAR
        double precision,    intent(in)     :: theta(:,:)
        type(Q0_compr),      intent(in)     :: Q0_All(:)
        double precision,    allocatable    :: TMP(:,:), Qrow(:), Qss(:,:), TMP1(:,:)
        integer                             :: i, j, k, p, it, start_v, end_v
        
if (associated(d_A%Qsy))      deallocate(d_A%Qsy)
if (associated(d_A%HQHt))     deallocate(d_A%HQHt)

select case (cv_A%Q_compression_flag)  !Select if the Q0 matrix is compressed or not     
     
 case(0) !Full Q0 matrix         
    
    ! Qss is the full matrix, made up of the kernel (Q0) multiplied by the appropiate current theta values
    !*****************************************************************************************************
    !Make Qss (Q) based on Q0 and variogram type
    !*****************************************************************************************************
    select case (cv_A%store_Q)
      case (.TRUE.)
        allocate(Qss(cv_PAR%npar,cv_PAR%npar)) ! Allocation
         Qss = 0.                              ! Initialization
         do i=1, cv_PAR%npar      !Loop over all the parameters
           select case (cv_S%var_type(d_PAR%BetaAssoc(i)))
            case (0) ! means nugget ---> just multiply by theta1
              Qss(i,i)=theta(d_PAR%BetaAssoc(i),1)*d_XQR%Q0(i,i)
            case (1) ! means linear ---> we need the maximum distance and theta1
              do j=i, cv_PAR%npar
                if (d_PAR%BetaAssoc(i).eq.d_PAR%BetaAssoc(j)) then !Search in the parameters list the associated parameters  
                  Qss(i,j)=theta(d_PAR%BetaAssoc(i),1)*d_XQR%L*exp(-d_XQR%Q0(i,j)/d_XQR%L)
                  Qss(j,i)=Qss(i,j)  ! Because Qss is symmetric            
                endif   
              enddo     
            case (2) ! means exponential ---> we need theta1 and theta2
              do j=i, cv_PAR%npar
                if (d_PAR%BetaAssoc(i).eq.d_PAR%BetaAssoc(j)) then !Search in the parameters list the associated parameters  
                   Qss(i,j)=theta(d_PAR%BetaAssoc(i),1)*exp(-d_XQR%Q0(i,j)/theta(d_PAR%BetaAssoc(i),2))
                   Qss(j,i)=Qss(i,j) ! Because Qss is symmetric
                endif   
              enddo
           end select ! Variogram type
         enddo             
      case (.FALSE.) ! We need to address this option
         allocate(Qrow(cv_PAR%npar))
    end select ! store_Q
    !*****************************************************************************************************
    ! End make Qss
    !*****************************************************************************************************
 
    !*****************************************************************************************************
    ! Make Qsy which is Qss*Ht
    !*****************************************************************************************************
    allocate(d_A%Qsy(cv_PAR%npar,cv_OBS%nobs)) ! Allocation
     d_A%Qsy = UNINIT_REAL                     ! Initialization                    
    
     call dgemm('n','t',cv_PAR%npar, cv_OBS%nobs, cv_PAR%npar, &
             1.D0, Qss, cv_PAR%npar, d_A%H, cv_OBS%nobs, &
             0.D0, d_A%Qsy, cv_PAR%npar)
    !*****************************************************************************************************
    ! End make Qsy 
    !*****************************************************************************************************

 case(1) !Compressed form of Q0 matrix
 
    !*****************************************************************************************************
    !Make Qsy which is Qss*Ht based on Q0 and variogram type. Qss is calculated on fly
    !*****************************************************************************************************
    select case (cv_A%store_Q)
      case (.TRUE.)
       allocate(d_A%Qsy(cv_PAR%npar,cv_OBS%nobs)) ! Allocation
        d_A%Qsy = UNINIT_REAL                      ! Initialization  
       do p = 1, cv_PAR%p  !Loop for each beta that correspond to each different Q0_C (each beta has a separate Q0_C)
        select case (cv_S%var_type(Q0_All(p)%BetaAss)) !Here the selection of the variogram type
          case (0) ! means nugget ---> just transpose the correct portion of H and multiply by theta1 
           d_A%Qsy(Q0_All(p)%Beta_Start:Q0_All(p)%Beta_Start+Q0_All(p)%npar-1,:) =  &  ! Select the correct position of Qsy
           & theta(Q0_All(p)%BetaAss,1)* &
           & (transpose(d_A%H(:,Q0_All(p)%Beta_Start:Q0_All(p)%Beta_Start+Q0_All(p)%npar-1))) !Portion of H(p)  
          case (1) ! means linear ---> we need the maximum distance and theta1. We have 2 option: Toeplitz or not
             select case (Q0_All(p)%Toep_flag) !Selection of Toeplitz [1] or not [0]
               case(0) !Means no Toeplitz.....just compressed form. Q0(p) is the full matrix for the p-th beta
                 allocate (TMP(Q0_All(p)%npar,Q0_All(p)%npar))
                 allocate (TMP1(Q0_All(p)%npar,cv_OBS%nobs))
                 start_v = Q0_All(p)%Beta_Start
                 end_v = Q0_All(p)%Beta_Start+Q0_All(p)%npar-1
                 do it=1,Q0_All(p)%npar
                   TMP(it,:)= exp(-Q0_All(p)%Q0_C(it,:)/d_XQR%L)
                 enddo
                 call dgemm('n','t',Q0_All(p)%npar, cv_OBS%nobs, Q0_All(p)%npar, &
                  (theta(Q0_All(p)%BetaAss,1)*d_XQR%L),TMP, Q0_All(p)%npar, &
                  & d_A%H(:,start_v:end_v), cv_OBS%nobs, &
                  & 0.D0, TMP1, Q0_All(p)%npar)
                 if (allocated(TMP))      deallocate(TMP)
                 do it =1,Q0_All(p)%npar
                    d_A%Qsy(start_v+it-1,:) = TMP1(it,:)
                 enddo
                 if (allocated(TMP1))      deallocate(TMP1)
               case(1) !Means Toeplitz. Q0(p) is just a vector with the distances
                 start_v = Q0_All(p)%Beta_Start
                 end_v = Q0_All(p)%Beta_Start+Q0_All(p)%npar-1
                 call toep_mult(Q0_All, p, d_A, cv_OBS%nobs, &
                 & (theta(Q0_All(p)%BetaAss,1)),d_XQR%L,d_XQR%L, &
                 & d_A%Qsy, start_v, end_v)
             end select  !Q0_All(p)%Toep_flag)      
          case (2) ! means exponential ---> we need theta1 and theta2. We have 2 option: Toeplitz or not
             select case (Q0_All(p)%Toep_flag) !Selection of Toeplitz [1] or not [0]
               case(0) !Means no Toeplitz.....just compressed form. Q0(p) is the full matrix for the p-th beta
                 allocate (TMP(Q0_All(p)%npar,Q0_All(p)%npar))
                 allocate (TMP1(Q0_All(p)%npar,cv_OBS%nobs))
                 start_v = Q0_All(p)%Beta_Start
                 end_v = Q0_All(p)%Beta_Start+Q0_All(p)%npar-1
                 do it=1,Q0_All(p)%npar
                   TMP(it,:)= exp(-Q0_All(p)%Q0_C(it,:)/(theta(Q0_All(p)%BetaAss,2)))
                 enddo
                   call dgemm('n','t',Q0_All(p)%npar, cv_OBS%nobs, Q0_All(p)%npar, &
                  (theta(Q0_All(p)%BetaAss,1)),TMP, Q0_All(p)%npar, &
                  & d_A%H(:,start_v:end_v), cv_OBS%nobs, &
                  & 0.D0, TMP1, Q0_All(p)%npar)
                 if (allocated(TMP))      deallocate(TMP)
                 do it =1,Q0_All(p)%npar
                    d_A%Qsy(start_v+it-1,:) = TMP1(it,:)
                 enddo
                 if (allocated(TMP1))      deallocate(TMP1)
               case(1) !Means Toeplitz. Q0(p) is just a vector with the distances
                 start_v = Q0_All(p)%Beta_Start
                 end_v = Q0_All(p)%Beta_Start+Q0_All(p)%npar-1
                 call toep_mult(Q0_All, p, d_A, cv_OBS%nobs, &
                   & (theta(Q0_All(p)%BetaAss,1)),(theta(Q0_All(p)%BetaAss,2)),1.D0 , &
                   & d_A%Qsy, start_v, end_v)
             end select  !Q0_All(p)%Toep_flag)
              
          end select !(cv_S%var_type(Q0_All(p)%BetaAss))
       enddo

    case (.FALSE.) ! We need to address this option
       allocate(Qrow(cv_PAR%npar))
    end select ! store_Q
    !*****************************************************************************************************
    ! End make Qsy
    !*****************************************************************************************************

end select !(cv_A%Q_compression_flag)

!*********************************************************************************************************
!*********************************************************************************************************
!********* The next lines are valid for both the full and compressed form of Q0 cases ********************
!*********************************************************************************************************
!*********************************************************************************************************


!*****************************************************************************************************
! Make HQsy which is H*Qss*Ht
!*****************************************************************************************************
  allocate(d_A%HQHt(cv_OBS%nobs,cv_OBS%nobs)) ! Allocation
    call dgemm('n', 'n', cv_OBS%nobs, cv_OBS%nobs, cv_PAR%npar, &
            1.D0, d_A%H, cv_OBS%nobs, d_A%Qsy,  cv_PAR%npar, &
            0.D0, d_A%HQHt, cv_OBS%nobs)        
!*****************************************************************************************************
! End make HQsy 
!*****************************************************************************************************

if (allocated(Qss))      deallocate(Qss)
if (allocated(Qrow))     deallocate(Qrow)
if (allocated(TMP))      deallocate(TMP)
if (allocated(TMP1))     deallocate(TMP1)


end subroutine bmo_form_Qss_Qsy_HQsy




subroutine bmo_form_Qyy(d_XQR, sig, cv_OBS, d_A)
        
        implicit none
        ! declarations
        type(kernel_XQR),    intent(in)     :: d_XQR
        double precision,    intent(in)     :: sig
        type(cv_observ),     intent(in)     :: cv_OBS  
        type(d_algorithmic), intent(inout)  :: d_A
        integer                             :: i, j

if (associated(d_A%Qyy))     deallocate(d_A%Qyy)

!*****************************************************************************************************
! Make Qyy which is H*Qss*Ht + sig*R0 = HQsy + sig*R0
!*****************************************************************************************************
  allocate(d_A%Qyy(cv_OBS%nobs,cv_OBS%nobs)) ! Allocation
    do i = 1, cv_OBS%nobs
      do j = 1, cv_OBS%nobs
        d_A%Qyy(i,j) = d_A%HQHt(i,j) + (sig*d_XQR%R0(i,j))
      enddo
    enddo
!*****************************************************************************************************
! End make Qyy 
!*****************************************************************************************************


end subroutine bmo_form_Qyy


!*****************************************************************************************************
!***** Subroutine to perform matrix operations that only involve H and s_old *************************
!***************************************************************************************************** 

subroutine bmo_H_only_operations(d_XQR, d_A,cv_OBS,d_PAR,cv_PAR)       
        implicit none
        ! declarations
        type(kernel_XQR),    intent(in)     :: d_XQR
        type(cv_param),      intent(in)     :: cv_PAR
        type(cv_observ),     intent(in)     :: cv_OBS
        type(d_algorithmic), intent(inout)  :: d_A
        type(d_param),       intent(inout)  :: d_PAR
        
if (associated(d_A%HX))        deallocate(d_A%HX)
if (associated(d_A%Hsold))     deallocate(d_A%Hsold)

!*****************************************************************************************************
! Make H*X
!*****************************************************************************************************
  allocate(d_A%HX(cv_OBS%nobs,cv_PAR%p))  ! Allocation
    call dgemm('n', 'n', cv_OBS%nobs, cv_PAR%p, cv_PAR%npar, &
            1.D0, d_A%H, cv_OBS%nobs, d_XQR%X, cv_PAR%npar, &
            0.D0, d_A%HX, cv_OBS%nobs)
!*****************************************************************************************************
! End make H*X
!*****************************************************************************************************


!*****************************************************************************************************
! Make Hsold which is H*d_PAR%pars_old  *** Is a vector (nobs) 
!*****************************************************************************************************    
  allocate(d_A%Hsold(cv_OBS%nobs))  !Allocation
    call dgemm('n', 'n', cv_OBS%nobs, 1, cv_PAR%npar, &
            1.D0, d_A%H, cv_OBS%nobs, d_PAR%pars_old,  cv_PAR%npar, &
            0.D0, d_A%Hsold, cv_OBS%nobs)
!*****************************************************************************************************
! End make Hsold 
!*****************************************************************************************************

end subroutine bmo_H_only_operations

!*****************************************************************************************************
!****** Subroutine to create and solve the linear cokriging system ***********************************
!*****************************************************************************************************
subroutine bmo_solve_linear_system(d_XQR, d_S, d_PM, cv_PAR, cv_OBS, d_OBS, d_A, d_PAR,cv_PM)
        
        implicit none
        ! declarations
        type(kernel_XQR),    intent(in)     :: d_XQR
        type(d_struct),      intent(inout)  :: d_S
        type(d_prior_mean),  intent(in)     :: d_PM
        type(cv_param),      intent(in)     :: cv_PAR
        type(cv_observ),     intent(in)     :: cv_OBS
        type(d_algorithmic), intent(inout)  :: d_A
        type(d_param),       intent(inout)  :: d_PAR
        type(d_observ),      intent(in)     :: d_OBS
        type (cv_prior_mean), intent(in)    :: cv_PM
        double precision,    allocatable    :: LHS(:,:), RHS (:) , Soln(:), C_S(:)
        integer                             :: i, j, k

if (associated(d_A%ksi))        deallocate(d_A%ksi)
if (associated(d_A%beta_hat))   deallocate(d_A%beta_hat)

!*****************************************************************************************************
! Make RHS which is y' and -InvQbbB0 *** Is a vector (nobs + p)
!*****************************************************************************************************    
  allocate(RHS(cv_OBS%nobs+cv_PAR%p))  !Allocation
    RHS= 0.                            !Initialization (IMPORTANT THAT IS = 0.)
    RHS(1:cv_OBS%nobs) = (d_OBS%obs  - d_OBS%h) + d_A%Hsold
    if (cv_PM%betas_flag .ne. 0) then 
      RHS(cv_OBS%nobs+1:cv_OBS%nobs+cv_PAR%p) = - d_PM%InvQbbB0
    endif
!*****************************************************************************************************
! End Make RHS 
!*****************************************************************************************************    

!*****************************************************************************************************
! Make LHS which is Qyy HX and (HX)t -InvQbb *** Is a matrix (nobs + p * nobs + p) 
! Here we fill the upper side of the matrix because symmetric
!*****************************************************************************************************  
  allocate(LHS(cv_OBS%nobs+cv_PAR%p, cv_OBS%nobs+cv_PAR%p))  !Allocation
   LHS = 0.                               !Initialization (IMPORTANT THAT IS = 0.)
   LHS(1:cv_OBS%nobs,1:cv_OBS%nobs)= d_A%Qyy
   LHS(1:cv_OBS%nobs,cv_OBS%nobs+1:cv_OBS%nobs+cv_PAR%p)= d_A%HX
    if (cv_PM%betas_flag .ne. 0) then
      LHS(cv_OBS%nobs+1:cv_OBS%nobs+cv_PAR%p,cv_OBS%nobs+1:cv_OBS%nobs+cv_PAR%p)= - d_PM%InvQbb
    endif
!*****************************************************************************************************
! End Make LHS 
! Here we filled the upper side of the matrix because symmetric
!*****************************************************************************************************  


!*****************************************************************************************************
! Scale LHS which is C_S'*LHS*C_S Is a matrix (nobs+p*nobs+p) LHS is overwritten by the scaled matrix
! Scaled using the matrix C_S(ii)=abs(LHS(i,i))**(-0.5) 
! Here we calculate the upper side of the matrix because symmetric, the lower side is not correct
! but we don't use that in the solution of the system LHS * InvC_S*Soln = RHS
! RHS is C_S * RHS and is overwritten too
!*******************************************************************************************************  
  
  allocate(C_S(cv_OBS%nobs+cv_PAR%p)) !This is the vector with the diagonal values of the scaling matrix  
    C_S   = UNINIT_REAL               !Initialization
  
      do i=1,cv_OBS%nobs+cv_PAR%p
      if (LHS(i,i).eq.0.) then    !If the value on the LHS diagonal is 0 then the scaling factor will be 1
        C_S(i) = 1. 
      else
        C_S(i) = abs(LHS(i,i))**(-0.5) !We use the absolute value because Qbb is negative
      endif
      LHS(i,:) = C_S(i) * LHS(i,:)
      LHS(:,i) = C_S(i) * LHS(:,i)
      RHS(i)   = C_S(i) * RHS(i)
    enddo 

!*****************************************************************************************************
! End Scale LHS 
!*****************************************************************************************************  


!*****************************************************************************************************
! Calculate the solution of the system LHS * InvC_S*Soln = RHS 
! *** Note that here LHS and RHS are scaled using the C_S matrix *** the original was overwritten
!*****************************************************************************************************
   
   allocate(Soln(cv_OBS%nobs+cv_PAR%p))
    Soln=RHS   !Initialization, Soln will be overwritten with the solution of the system
    call SLVSSU(cv_OBS%nobs+cv_PAR%p,LHS,Soln) ! Call to solve LHS * Soln = RHS
!*** Warning: On exit LHS is overwritten by the factorization used for the solution ***
!*** The solution of this system is C_S^-1 * Soln. Soln must be rescaled.

!Rescaling the solution and assign the values to ksi and beta_hat
    allocate(d_A%ksi(cv_OBS%nobs))    !Allocation
    allocate(d_A%beta_hat(cv_PAR%p))  !Allocation   
    
    do i=1,cv_OBS%nobs                      !Rescaling the solution and assign to ksi 
      d_A%ksi(i) = C_S(i) * Soln(i)
    enddo
    do i=cv_OBS%nobs+1,cv_OBS%nobs+cv_PAR%p !Rescaling the solution and assign to beta_hat
      d_A%beta_hat(i-cv_OBS%nobs) = C_S(i) * Soln(i)
    enddo      
    
!*****************************************************************************************************
! End Calculate the solution  
!*****************************************************************************************************

if (allocated(RHS))        deallocate(RHS)
if (allocated(LHS))        deallocate(LHS)
if (allocated(Soln))       deallocate(Soln)
if (allocated(C_S))        deallocate(C_S)

!*****************************************************************************************************
! Calculate best estimate s_hat which is d_PAR%pars = X*beta_hat + Qsy * ksi
!***************************************************************************************************** 
 
  call DGEMV('n',cv_PAR%npar,cv_PAR%p,1.D0,d_XQR%X,cv_PAR%npar,  &
       d_A%beta_hat,1,0.D0,d_PAR%pars,1)
  call DGEMV('n',cv_PAR%npar,cv_OBS%nobs,1.D0,d_A%Qsy,cv_PAR%npar,  &
       d_A%ksi,1,1.D0,d_PAR%pars,1)
       
!*****************************************************************************************************
! End Calculate best estimate 
!*****************************************************************************************************   


!*****************************************************************************************************
! Calculate the objective functions (calling file calc_ob_fun)
!***************************************************************************************************** 

call cal_ob_funcs(d_XQR, d_S, d_PM, cv_PAR, cv_OBS, d_OBS, d_A, d_PAR, cv_PM)

!*****************************************************************************************************
! End Calculate the objective functions
!***************************************************************************************************** 

end subroutine bmo_solve_linear_system



!*****************************************************************************************************
!****** Subroutine to make Qsy in case of Toeplitz matrix ********************************************
!*****************************************************************************************************
subroutine toep_mult(Q0_All,ip,d_A,nobs,theta_1,theta_2,Lmax,Qsy,start_v,end_v)

type(Q0_compr),       intent(in)     :: Q0_All(:)
type(d_algorithmic),  intent(in)     :: d_A
double precision,     intent(inout)  :: Qsy(:,:) 
double precision,     intent(in)     :: theta_1,theta_2,Lmax
integer,              intent(in)     :: nobs
double precision, allocatable        :: Qtmpb(:),Qtmpg(:),Qtmpl(:),Qv(:),TMP(:)
double precision, allocatable        :: TMVSY(:)
integer                              :: ncol,nbl,nlay
integer                              :: blkg,blkl
integer                              :: i,l,k,p,it,jt,ip
integer                              :: start_v, end_v

!Note: In case of linear variogram theta_1 must be theta_1, 
!theta_2 and Lmax must be the 10 times the maximum distance in Q0_All
!In case of exponential variogram theta_1 must be theta_1,
!theta_2 must be theta_2 and Lmax must be 1

allocate (Qtmpb(Q0_All(ip)%npar))
allocate (Qtmpg(Q0_All(ip)%npar))
allocate (Qtmpl(Q0_All(ip)%npar))
allocate (Qv(Q0_All(ip)%npar))
allocate (TMP(Q0_All(ip)%npar))
allocate (TMVSY(nobs))

ncol=Q0_All(ip)%Ncol
nbl =Q0_All(ip)%Nrow
nlay=Q0_All(ip)%Nlay
Qv=0.
Qtmpb=Q0_All(ip)%Q0_C(:,1)
Qtmpg=Q0_All(ip)%Q0_C(:,1)
Qtmpl=Q0_All(ip)%Q0_C(:,1)
new_block=.true.
blkg=1
blkl=1
i=0

do p=1, Q0_All(ip)%npar !Index for all the columns of the matrix

if (i/ncol.eq.1) then
  if(blkg/nbl.eq.1)    then
     blkl = blkl+1
     Qtmpb(1:(ncol*nbl)) = Q0_All(ip)%Q0_C((ncol*nbl*(blkl-1))+1:((ncol*nbl)*(blkl-1))+(ncol*nbl),1)
     Qtmpb((ncol*nbl)+1:(ncol*nbl*nlay)) = Qtmpl(1:(ncol*nbl*nlay)-(ncol*nbl))
     Qtmpl = Qtmpb
     Qtmpg = Qtmpb
     blkg = 1
     i=1
     new_block=.true.
  else
     blkg = blkg+1
     do l= 1,(nlay*ncol*nbl),(ncol*nbl)
       Qtmpb(l:l+ncol-1)= Qtmpl((ncol*(blkg-1))+l:((ncol-1)*(blkg-1))+l+ncol)
       Qtmpb(l+ncol:ncol*nbl+(l-1)) = Qtmpg(l:(ncol*nbl)-ncol+(l-1))
     enddo
     Qtmpg=Qtmpb
     new_block=.true.
     i=1
  endif
else
i=i+1
endif
if (new_block) then
 Qv=Qtmpb
 new_block=.false.
else
 do l= 1,(nlay*ncol*nbl),(ncol*nbl)
   do k=1,nbl
     Qv(ncol*(k-1)+l)= Qtmpg(ncol*(k-1)+i+l-1)
     Qv(ncol*(k-1)+1+l:ncol*(k-1)+ncol+l-1) = Qtmpb(ncol*(k-1)+l:ncol*(k-1)+ncol-2+l)
   enddo
 enddo 
endif
Qtmpb=Qv

!******************************************************************************************************************
!********* Here, for each p, Qv is the vector that contain the value of the p-th colummn of the matrix Q **********
!********** From here Qv is available to be used in some calculation  *********************************************
!******************************************************************************************************************

!**** Here we calculate for each column of Q0 (Qv of the p-th iteration in this subroutine) H*Qt that is (Q*Ht)t **
!**** We assign the result to the p-th row (instead of column) of Qsy to obtain Q*Ht ******************************
!**** that, at the end of the loop, is the Qsy for the specified beta  ********************************************
!******************************************************************************************************************

TMP = (theta_1*Lmax*exp(-Qv/theta_2))
call DGEMV('n',nobs,Q0_All(ip)%npar,1.D0,d_A%H(:,start_v:end_v),nobs,TMP,1,0.D0,TMVSY,1)
Qsy(Q0_All(ip)%Beta_start+p-1,:)=TMVSY
enddo  !End of loop for each column of the entire Q matrix

if (allocated(Qtmpb)) deallocate(Qtmpb)
if (allocated(Qtmpg)) deallocate(Qtmpg)
if (allocated(Qtmpl)) deallocate(Qtmpl)
if (allocated(Qv))    deallocate(Qv)
if (allocated(TMP))   deallocate(TMP)
if (allocated(TMVSY)) deallocate(TMVSY)


end subroutine toep_mult
!*****************************************************************************************************
!****** End Subroutine to make Qsy in case of Toeplitz matrix ****************************************
!*****************************************************************************************************


end module bayes_matrix_operations

