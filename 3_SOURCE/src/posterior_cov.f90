module posterior_cov_operations

!***** Created by M.D. 09/27/11 ********
        use bayes_pest_control
        use utilities  
        use objective_function


contains


!*******************************************************************************************************
!****************** Subroutine to calculate posterior covariance matrix (VV) ***************************
!****** VV is the full matrix in case of no compression of Q, V is the diagonal of VV otherwise ********
!******************************************************************************************************* 

subroutine form_post_covariance(d_XQR, cv_PAR, cv_OBS, cv_S, cv_A, d_A, d_PAR,Q0_All,cv_PM,d_PM,d_S,VV,V)
        
        implicit none
        ! declarations
        type(kernel_XQR),    intent(in)     :: d_XQR
        type(cv_struct),     intent(in)     :: cv_S 
        type(cv_param),      intent(in)     :: cv_PAR
        type(cv_observ),     intent(in)     :: cv_OBS
        type(cv_algorithmic),intent(inout)  :: cv_A        
        type(d_algorithmic), intent(inout)  :: d_A
        type(d_param),       intent(inout)  :: d_PAR
        type(cv_prior_mean), intent(in)     :: cv_PM
        type (d_prior_mean), intent(in)     :: d_PM
        type (d_struct),     intent(in)     :: d_S
        type(Q0_compr),      intent(in)     :: Q0_All(:)
        double precision,    allocatable    :: TMV(:), TMP(:,:), Qrow(:), TMP1(:,:)
        double precision,    pointer        :: V(:), VV(:,:)
        integer                             :: i, j, k, p, it, start_v, end_v
      
  if (associated(d_A%Qsy))      deallocate(d_A%Qsy)
  if (associated(d_A%HQHt))     deallocate(d_A%HQHt)
  if (associated(d_A%Qyy))      deallocate(d_A%Qyy)

  
  select case (cv_A%Q_compression_flag)  !Select if the Q0 matrix is compressed or not       
      
  case(0) !Full Q0 matrix         
    
    ! Qss is the full matrix, made up of the kernel (Q0) multiplied by the appropiate current theta values
    !*****************************************************************************************************
    !Make Qss (Q) based on Q0 and variogram type (Qss is stored in VV to avoid allocation of two matrices)
    !*****************************************************************************************************
    select case (cv_A%store_Q)
      case (.TRUE.)
        allocate(VV(cv_PAR%npar,cv_PAR%npar)) ! Allocation
        VV = 0.                              ! Initialization
         do i=1, cv_PAR%npar      !Loop over all the parameters
           select case (cv_S%var_type(d_PAR%BetaAssoc(i)))
            case (0) ! means nugget ---> just multiply by theta1
              VV(i,i)=d_S%theta(d_PAR%BetaAssoc(i),1)*d_XQR%Q0(i,i)
            case (1) ! means linear ---> we need the maximum distance and theta1
              do j=i, cv_PAR%npar
                if (d_PAR%BetaAssoc(i).eq.d_PAR%BetaAssoc(j)) then !Search in the parameters list the associated parameters  
                  VV(i,j)=d_S%theta(d_PAR%BetaAssoc(i),1)*d_XQR%L*exp(-d_XQR%Q0(i,j)/d_XQR%L)
                  VV(j,i)=VV(i,j)  ! Because Qss is symmetric            
                endif   
              enddo     
            case (2) ! means exponential ---> we need theta1 and theta2
              do j=i, cv_PAR%npar
                if (d_PAR%BetaAssoc(i).eq.d_PAR%BetaAssoc(j)) then !Search in the parameters list the associated parameters  
                   VV(i,j)=d_S%theta(d_PAR%BetaAssoc(i),1)*exp(-d_XQR%Q0(i,j)/d_S%theta(d_PAR%BetaAssoc(i),2))
                   VV(j,i)=VV(i,j) ! Because Qss is symmetric
                endif   
              enddo
           end select ! Variogram type
         enddo 
         !************************************************************************************************            
         !******* if we have prior information on beta, Qss becomes Gss = Qss + XQbbXt *******************
         if (cv_PM%betas_flag.ne.0) then   !------> we have prior information about beta   
           allocate(TMP(cv_PAR%npar,cv_PAR%p)) !Calculate XQbb
           call dgemm('n', 'n', cv_PAR%npar, cv_PAR%p, cv_PAR%p, 1.D0, d_XQR%X, &
            cv_PAR%npar, d_PM%Qbb, cv_PAR%p, 0.D0,  TMP, cv_PAR%npar)
           ! now multiply XQbb (TMP) by Xt and add the result (XQbbXt) to Qss, forming Gss
           call dgemm('n', 't', cv_PAR%npar, cv_PAR%npar, cv_PAR%p, 1.D0, TMP, &
              cv_PAR%npar, d_XQR%X, cv_PAR%npar, 1.D0,  VV, cv_PAR%npar)
           if (allocated(TMP)) deallocate(TMP)
         endif
        !**************************************************************************************************
        !************** From here VV is Gss if we have prior information on beta *************************
        !**************************************************************************************************
      case (.FALSE.) ! We need to address this option
         allocate(Qrow(cv_PAR%npar))
    end select ! store_Q
    !*****************************************************************************************************
    ! End make Qss or Gss 
    !*****************************************************************************************************
 
    !*****************************************************************************************************
    ! Make Qsy which is Qss*Ht or Gss*Ht
    !*****************************************************************************************************
    allocate(d_A%Qsy(cv_PAR%npar,cv_OBS%nobs)) ! Allocation
     d_A%Qsy = UNINIT_REAL                     ! Initialization                    
    call dgemm('n','t',cv_PAR%npar, cv_OBS%nobs, cv_PAR%npar, &
             1.D0, VV, cv_PAR%npar, d_A%H, cv_OBS%nobs, &
             0.D0, d_A%Qsy, cv_PAR%npar)
    !*****************************************************************************************************
    ! End make Qsy 
    !*****************************************************************************************************
    
    !*****************************************************************************************************
    ! Make HQsy which is H*Qss*Ht or H*Gss*Ht 
    !*****************************************************************************************************
    allocate(d_A%HQHt(cv_OBS%nobs,cv_OBS%nobs)) ! Allocation
      call dgemm('n', 'n', cv_OBS%nobs, cv_OBS%nobs, cv_PAR%npar, &
            1.D0, d_A%H, cv_OBS%nobs, d_A%Qsy,  cv_PAR%npar, &
            0.D0, d_A%HQHt, cv_OBS%nobs)        
    !*****************************************************************************************************
    ! End make HQsy 
    !*****************************************************************************************************
    
    !*****************************************************************************************************
    ! Make Qyy = H*Qss*Ht + sig*R0 or H*Gss*Ht +sig*R0 = Gyy and calculate the inverse
    !*****************************************************************************************************
    allocate(d_A%Qyy(cv_OBS%nobs,cv_OBS%nobs)) ! Allocation
    do i = 1, cv_OBS%nobs
      do j = 1, cv_OBS%nobs
        d_A%Qyy(i,j) = d_A%HQHt(i,j) + (d_S%sig*d_XQR%R0(i,j))
      enddo
    enddo
    call INVGM(cv_OBS%nobs,d_A%Qyy) !Here d_A%Qyy is the inverse of Gyy
    !*****************************************************************************************************
    ! End make Qyy 
    !*****************************************************************************************************
    
    !*****************************************************************************************************
    !************* Make VV = Gss - GssHt*Gyy^-1*GssH (Gss is in VV and GssHt is in d_A%Qyy) **************
    !*****************************************************************************************************
    allocate(TMP(cv_PAR%npar,cv_OBS%nobs)) !Calculate GssHt*Gyy^-1
      call dgemm('n', 'n', cv_PAR%npar, cv_OBS%nobs, cv_OBS%nobs, 1.D0, d_A%Qsy, &
          cv_PAR%npar, d_A%Qyy, cv_OBS%nobs, 0.D0,  TMP, cv_PAR%npar)
      !now multiply GssHtGyy^-1 (TMP) by (GssHt)t, change the sign and add the result to Gss, forming VV
      call dgemm('n', 't', cv_PAR%npar, cv_PAR%npar, cv_OBS%nobs, -1.D0, TMP, &
            cv_PAR%npar, d_A%Qsy , cv_PAR%npar, 1.D0,  VV, cv_PAR%npar)
    if (allocated(TMP)) deallocate(TMP)
    !*****************************************************************************************************
    ! End make VV (The matrix VV is now the posterior covariance matrix) (FULL MATRIX)
    !*****************************************************************************************************
 case(1) !Compressed form of Q0 matrix
 
    !*****************************************************************************************************
    !Make Qsy which is Qss*Ht based on Q0 and variogram type. Qss is calculated on fly
    !*****************************************************************************************************
    select case (cv_A%store_Q)
      case (.TRUE.)
       allocate(d_A%Qsy(cv_PAR%npar,cv_OBS%nobs)) ! Allocation
        d_A%Qsy = UNINIT_REAL                      ! Initialization  
       allocate (V(cv_PAR%npar))             !Allocate the vector for the diagonal values of the posterior covariance 
        V = UNINIT_REAL                      ! Initialization  
       do p = 1, cv_PAR%p  !Loop for each beta that correspond to each different Q0_C (each beta has a separate Q0_C)
        select case (cv_S%var_type(Q0_All(p)%BetaAss)) !Here the selection of the variogram type
          case (0) ! means nugget ---> just transpose the correct portion of H and multiply by theta1 
           d_A%Qsy(Q0_All(p)%Beta_Start:Q0_All(p)%Beta_Start+Q0_All(p)%npar-1,:) =  &  ! Select the correct position of Qsy
            & d_S%theta(Q0_All(p)%BetaAss,1)* &
            & (transpose(d_A%H(:,Q0_All(p)%Beta_Start:Q0_All(p)%Beta_Start+Q0_All(p)%npar-1))) !Portion of H(p)  
           V(Q0_All(p)%Beta_Start:Q0_All(p)%Beta_Start+Q0_All(p)%npar-1)=d_S%theta(Q0_All(p)%BetaAss,1)
          case (1) ! means linear ---> we need the maximum distance and theta1. We have 2 option: Toeplitz or not
             select case (Q0_All(p)%Toep_flag) !Selection of Toeplitz [1] or not [0]
               case(0) !Means no Toeplitz.....just compressed form. Q0(p) is the full matrix for the p-th beta
                 allocate (TMP(Q0_All(p)%npar,Q0_All(p)%npar))
                 allocate (TMP1(Q0_All(p)%npar,cv_OBS%nobs))
                 start_v = Q0_All(p)%Beta_Start
                 end_v = Q0_All(p)%Beta_Start+Q0_All(p)%npar-1
                 do it=1,Q0_All(p)%npar
                   TMP(it,:)= exp(-Q0_All(p)%Q0_C(it,:)/d_XQR%L)
                   V(start_v:end_v) = TMP(it,it)*(d_S%theta(Q0_All(p)%BetaAss,1)*d_XQR%L) !Only diagonal values of Qss
                 enddo
                 call dgemm('n','t',Q0_All(p)%npar, cv_OBS%nobs, Q0_All(p)%npar, &
                  (d_S%theta(Q0_All(p)%BetaAss,1)*d_XQR%L),TMP, Q0_All(p)%npar, &
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
                 call toep_mult_post(Q0_All,p, d_A, cv_OBS%nobs, &
                 & (d_S%theta(Q0_All(p)%BetaAss,1)),d_XQR%L,d_XQR%L, d_A%Qsy , V, &
                 & start_v, end_v)
             end select  !Q0_All(p)%Toep_flag)      
          case (2) ! means exponential ---> we need theta1 and theta2. We have 2 option: Toeplitz or not
             select case (Q0_All(p)%Toep_flag) !Selection of Toeplitz [1] or not [0]
               case(0) !Means no Toeplitz.....just compressed form. Q0(p) is the full matrix for the p-th beta
                 allocate (TMP(Q0_All(p)%npar,Q0_All(p)%npar))
                 allocate (TMP1(Q0_All(p)%npar,cv_OBS%nobs))
                 start_v = Q0_All(p)%Beta_Start
                 end_v = Q0_All(p)%Beta_Start+Q0_All(p)%npar-1
                 do it=1,Q0_All(p)%npar
                   TMP(it,:)= exp(-Q0_All(p)%Q0_C(it,:)/(d_S%theta(Q0_All(p)%BetaAss,2)))
                   V(start_v:end_v) = TMP(it,it)*(d_S%theta(Q0_All(p)%BetaAss,1)) !Only diagonal values of Qss
                 enddo
                   call dgemm('n','t',Q0_All(p)%npar, cv_OBS%nobs, Q0_All(p)%npar, &
                   (d_S%theta(Q0_All(p)%BetaAss,1)),TMP, Q0_All(p)%npar, &
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
                 call toep_mult_post(Q0_All, p, d_A, cv_OBS%nobs, &
                   & (d_S%theta(Q0_All(p)%BetaAss,1)),(d_S%theta(Q0_All(p)%BetaAss,2)),1.D0 , &
                   & d_A%Qsy,V, start_v, end_v)
             end select  !Q0_All(p)%Toep_flag)
              
          end select !(cv_S%var_type(Q0_All(p)%BetaAss))
       enddo

      !************************************************************************************************            
      !********* If we have prior information on beta Qsy will be Gsy = (Qss + XQbbXt)*Ht *************
      !********* V was the diagonal of Qss and will be the diagonal of Gss = (Qss + XQbbXt) ***********
       if (cv_PM%betas_flag.ne.0) then   !------> we have prior information about beta   
         allocate (TMV(cv_OBS%nobs))
         do i = 1,cv_PAR%npar
           TMV = 0.
           do p = 1,cv_PAR%p
             start_v = Q0_All(p)%Beta_Start
             end_v   = Q0_All(p)%Beta_Start+Q0_All(p)%npar-1
             TMV = TMV + d_PM%Qbb(d_PAR%BetaAssoc(i),p)*sum(d_A%H(:,start_v:end_v), DIM=2)
           enddo
           d_A%Qsy(i,:) = d_A%Qsy(i,:) + TMV
           V(i) = V(i) + d_PM%Qbb(d_PAR%BetaAssoc(i),d_PAR%BetaAssoc(i))
         enddo    
       endif 
      !**************************************************************************************************
      !****** From here Qsy is Gsy and V the diagonal of Gss if we have prior information on beta *******
      !**************************************************************************************************
    
      !*****************************************************************************************************
      ! Make HQsy which is H*(Qss*Ht) or H*(Gss*Ht)
      !*****************************************************************************************************
       allocate(d_A%HQHt(cv_OBS%nobs,cv_OBS%nobs)) ! Allocation
         call dgemm('n', 'n', cv_OBS%nobs, cv_OBS%nobs, cv_PAR%npar, &
            1.D0, d_A%H, cv_OBS%nobs, d_A%Qsy,  cv_PAR%npar, &
            0.D0, d_A%HQHt, cv_OBS%nobs)        
      !*****************************************************************************************************
      ! End make HQsy 
      !*****************************************************************************************************
    
      !*****************************************************************************************************
      ! Make Qyy = H*Qss*Ht + sig*R0 or H*Gss*Ht +sig*R0 = Gyy and calculate the inverse
      !*****************************************************************************************************
       allocate(d_A%Qyy(cv_OBS%nobs,cv_OBS%nobs)) ! Allocation
       do i = 1, cv_OBS%nobs
        do j = 1, cv_OBS%nobs
          d_A%Qyy(i,j) = d_A%HQHt(i,j) + (d_S%sig*d_XQR%R0(i,j))
        enddo
       enddo
       call INVGM(cv_OBS%nobs,d_A%Qyy) !Here d_A%Qyy is the inverse of Gyy
      !*****************************************************************************************************
      ! End make Qyy 
      !*****************************************************************************************************
      
      !*****************************************************************************************************
      ! Make V = diag(Gss - GssHt*Gyy^-1*GssH) (Only diagonal terms in this case)
      !*****************************************************************************************************
       allocate(TMP(cv_PAR%npar,cv_OBS%nobs)) !Calculate GssHt*Gyy^-1
        call dgemm('n', 'n', cv_PAR%npar, cv_OBS%nobs, cv_OBS%nobs, 1.D0, d_A%Qsy, &
          cv_PAR%npar, d_A%Qyy, cv_OBS%nobs, 0.D0,  TMP, cv_PAR%npar)
        do i = 1,cv_PAR%npar  !Now calculate the diagonal of V = diag(Gss)-diag(GssHt*Gyy^-1*GssH)
         V(i) = V(i) - sum(TMP(i,:)*d_A%Qsy(i,:))
        enddo 
       if (allocated(TMP)) deallocate(TMP)
      !*****************************************************************************************************
      ! End make V posterior covariance (ONLY DIAGONAL TERMS)
      !*****************************************************************************************************
        
    case (.FALSE.) ! We need to address this option
       allocate(Qrow(cv_PAR%npar))
    end select ! store_Q
    !*****************************************************************************************************
    ! End make Qsy
    !*****************************************************************************************************
end select !(cv_A%Q_compression_flag)

if (allocated(Qrow))     deallocate(Qrow)

end subroutine form_post_covariance



!*****************************************************************************************************
!****** Subroutine to make Qsy in case of Toeplitz matrix ********************************************
!*****************************************************************************************************
subroutine toep_mult_post(Q0_All,ip,d_A,nobs,theta_1,theta_2,Lmax,Qsy,V,start_v,end_v)

type(Q0_compr),       intent(in)     :: Q0_All(:)
type(d_algorithmic),  intent(in)     :: d_A
double precision,     intent(inout)  :: Qsy(:,:), V(:) 
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
V(Q0_All(ip)%Beta_start+p-1) = TMP(p) !V here is the diagonal of Qss

call DGEMV('n',nobs,Q0_All(ip)%npar,1.D0,d_A%H(:,start_v:end_v),nobs,TMP,1,0.D0,TMVSY,1)
Qsy(Q0_All(ip)%Beta_start+p-1,:)=TMVSY
enddo  !End of loop for each column of the entire Q matrix

if (allocated(Qtmpb)) deallocate(Qtmpb)
if (allocated(Qtmpg)) deallocate(Qtmpg)
if (allocated(Qtmpl)) deallocate(Qtmpl)
if (allocated(Qv))    deallocate(Qv)
if (allocated(TMP))   deallocate(TMP)
if (allocated(TMVSY)) deallocate(TMVSY)

end subroutine toep_mult_post
!*****************************************************************************************************
!****** End Subroutine to make Qsy in case of Toeplitz matrix ****************************************
!*****************************************************************************************************


end module posterior_cov_operations