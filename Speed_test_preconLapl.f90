program speed_test
implicit none
        
integer:: N,M, mult

do mult=1,500,1
N=mult*32*2+2
M=mult*32
call run_resolution(N,M, mult)
enddo

contains

subroutine run_resolution (N,M, mult)
implicit none

INTEGER :: IP(N), mult
INTEGER :: N, M, I, J, K, Q
integer :: npromaw,  M_Q, M_r, point
integer :: M2(2), M3(2), MM1(2), MM2(2), iter, iterations
integer, allocatable :: JtoJM1(:,:,:), JtoJP1(:,:,:)        
REAL(Kind=4) :: start_time_all, end_time_all
DO I=1,N
  IP(I)=MOD(I+(N-2)/2-1,N-2)+1
end do

!do  npromaw=16,16 !min(M,512),4
!write(*,*) 'nproma', npromaw
!if(mod(M,npromaw)==0) then
! Allocate(JtoJM1(npromaw,M/npromaw,2))
! Allocate(JtoJP1(npromaw,M/npromaw,2))
! do K=1,M/npromaw
!   do J=1,npromaw
!      point=(k-1)*npromaw+j-1
!      JtoJM1(J,K,1)=(point+npromaw-1)/npromaw
!      JtoJM1(J,K,2)= mod(point-1,npromaw)+1
!   enddo
! enddo
! do K=1,M/npromaw
!   do J=1,npromaw
!      point=(k-1)*npromaw+j+1
!      JtoJP1(J,K,1)=(point+npromaw-1)/npromaw
!      JtoJP1(J,K,2)= mod(point-1,npromaw)+1
!   enddo
! enddo
! M_Q=npromaw
!      point=2
!      M2(1)=(point+M_Q-1)/M_Q
!      M2(2)= mod(point-1,M_Q)+1
!      point=3
!      M3(1)=(point+M_Q-1)/M_Q
!      M3(2)= mod(point-1,M_Q)+1
!
!      point=M-1
!      MM1(1)=(point+M_Q-1)/M_Q
!      MM1(2)= mod(point-1,M_Q)+1
!      point=M-2
!      MM2(1)=(point+M_Q-1)/M_Q
!      MM2(2)= mod(point-1,M_Q)+1
!
!
!call precon_ADI_Piotr_Trans_Block_ARM_RP(N, M,M/npromaw,npromaw,npromaw,JtoJM1,JtoJP1,IP,M2, M3, MM1, MM2)
!call LAPL_Trans_Block_ARM_RP(N, M,M/npromaw,npromaw,npromaw,JtoJM1,JtoJP1,IP,M2, M3, MM1, MM2)
!   write(*,*)'JtoJM1'
!   write(*,*)JtoJM1
!  ! write(*,*)'JtoJP1'
!  ! write(*,*)JtoJP1
!  else
! Allocate(JtoJM1(npromaw,M/npromaw+1,2))
! Allocate(JtoJP1(npromaw,M/npromaw+1,2))
! do K=1,M/npromaw+1
!   do J=1,npromaw
!      point=(k-1)*npromaw+j-1
!      JtoJM1(J,K,1)=(point+npromaw-1)/npromaw
!      JtoJM1(J,K,2)= mod(point-1,npromaw)+1
!   enddo
! enddo
! do K=1,M/npromaw+1
!   do J=1,npromaw
!      point=(k-1)*npromaw+j+1
!      JtoJP1(J,K,1)=(point+npromaw-1)/npromaw
!      JtoJP1(J,K,2)= mod(point-1,npromaw)+1
!      enddo
! enddo
!   write(*,*)'JtoJM1'
!   write(*,*)JtoJM1
!  ! write(*,*)'JtoJP1'
!  ! write(*,*)JtoJP1
! M_Q=npromaw
!       point=2
!      M2(1)=(point+M_Q-1)/M_Q
!      M2(2)= mod(point-1,M_Q)+1
!      point=3
!      M3(1)=(point+M_Q-1)/M_Q
!      M3(2)= mod(point-1,M_Q)+1
!
!      point=M-1
!      MM1(1)=(point+M_Q-1)/M_Q
!      MM1(2)= mod(point-1,M_Q)+1
!      point=M-2
!      MM2(1)=(point+M_Q-1)/M_Q
!      MM2(2)= mod(point-1,M_Q)+1
!call precon_ADI_Piotr_Trans_Block_ARM_RP(N, M,M/npromaw+1,npromaw,mod(M,npromaw),JtoJM1,JtoJP1,IP,M2, M3, MM1, MM2)
!call LAPL_Trans_Block_ARM_RP(N, M,M/npromaw+1,npromaw,mod(M,npromaw),JtoJM1,JtoJP1,IP,M2, M3, MM1, MM2)
!  endif
!
!  DEALLOCATE(JtoJM1)
!  DEALLOCATE(JtoJP1)
!enddo
!write(*,*) 'NO INPUT'
if(mult<60) then
        iterations=100
else
        iterations=1
endif

call cpu_time(start_time_all)

do iter=1,iterations

call precon_ADI_Piotr_Trans_ARM_RP_noInput(N,M,IP)
call LAPL_Trans_ARM_RP(N,M, IP)
end do
call cpu_time(end_time_all)

write(*,*) 'Resolution', mult, (end_time_all-start_time_all)/iterations

end subroutine

SUBROUTINE precon_ADI_Piotr_Trans_Block_ARM_RP(N, M,Q,M_Q,M_R,JtoJM1,JtoJP1,IP,M2, M3, MM1, MM2)
!SUBROUTINE precon_ADI_Piotr_Trans_Block_ARM_RP(R,QU , T_step,  A, ps, &
!&  divi,A11,A12,A21,A22,B11,B22,P0,U,V,S,S_full,N,M,IP,ID_PREC, num_of_bits,DP_Depth)
!use implicit_functions_HP

        implicit none
        integer, parameter :: HP=2
integer :: N, M, Q, K, I, J
integer :: M_Q, M_R, counter
! ws REAL(kind=2)::  ps(M_Q,0:N-2,Q), ws(M_Q,N-1,0:4,Q), qs(M_Q,0:N-1,0:2,Q)
!double precision:: QU(M_Q,N,Q), ws(M_Q,N-1,Q,0:4), aa(M_Q,Q,1:4)
Real(Kind=4):: S_full
REAL(kind=HP):: R(M_Q,N,Q) , QU(M_Q,N,Q), aa(M_Q,Q,1:4),ws(M_Q,N-1,Q,0:4)
REAL(kind=HP):: qs(M_Q,0:N-1,Q,0:2), rhs(M_Q,N,Q),  ps(M_Q,0:N-2,Q),&
        & S(M_Q,N,Q),  divi(M_Q,N,Q) , F(M_Q,N,Q), help
REAL(kind=HP) :: RPE_0,RPE_1,RPE_05,RPE_025,RPE_0125,RPE_2,RPE_3,RPE_4 , rpe_100
REAL(kind=HP) :: U_p(M_Q,N,Q),V_p(M_Q,N,Q), A11(M_Q,N,Q), A12(M_Q,N,Q),A21(M_Q,N,Q),A22(M_Q,N,Q),  &
    &      B11(M_Q,N,Q),B22(M_Q,N,Q),P0(M_Q,N,Q),U(M_Q,N,Q),V(M_Q,N,Q),  C11(M_Q,N,Q)
REAL(kind=2):: deti, det40, det41, det42, det43, det44, det3_ARM_RP,   &
                           & d11, d12, d13, d14, d21, d22, d23, d24,       &
                           & d31, d32, d33, d34, d41, d42, d43, d44,       &
                           & s1, s2, s3, s4
REAL(kind=2):: T_step, Delta_t_I, dn, dni, util, vtil, swcp, max_QX_QY
REAL(kind=4):: start_time, end_time !time_begin, time_end
integer :: point, JtoJM1(M_Q,Q,2), JtoJP1(M_Q,Q,2), M2(2), M3(2), MM1(2), MM2(2)
integer :: max_iter, iteration, il1, il2, i2, IP(N)
!REAL(Kind=HP) :: S_T(M,N), C11_T(M,N), ps_T(M,0:N-2), qs_T(M,0:N-1,0:2),divi_T(M,N)&
!        & , rhs_T(M,N)

      rpe_1 = (1.0)
    rpe_2 = (2.0)
    rpe_3 = (3.0)
    rpe_4 = (4.0)
    rpe_100 = (100.0)
    rpe_0125 = (0.125)
    rpe_025 = (0.25)
    rpe_05 = (0.5)
    rpe_0 = (0.0)
    write(*,*) 'init in subroutine'
!call rpenum_init(10)
write(*,*) 'Precision', HP 
max_iter= 2
swcp=rpe_1
!initialize the inverse of R with 0

!! into precon_prep_depth
    do K=1,Q-1,1   
      DO I=1,N
        DO J=1,M_Q
                                   
      C11(J,I,K)=0.5_HP*a11(J,I,K)
        enddo
      enddo
   enddo
   K=Q
      DO I=1,N
        DO J=1,M_R
      C11(J,I,K)=0.5_HP*a11(J,I,K)
        enddo
      enddo
Delta_t_I=1.0_HP/(0.25_HP/T_step)

!T_step_dp=T_step
!write(*,*) 'T_step,', T_step_dp
!write(*,*) 'T_step,', float(rpe_025)
!write(*,*) 'T_step,', float(rpe_025/T_step)
!
!Delta_t_I_dp=Delta_t_I
!write(*,*) 'Delta_t_I,', Delta_t_I_dp
!read(*,*)
!
!      DO J=1,M
!
!        ps(0,j)=rpe_0
!        divi(1,J)=rpe_1/(c11(2,j)+c11(N-2,j) +s(1,j)*(Delta_t_I+rpe_1))
!        ps(1,j) = c11(2,j)*divi(1,J)
!
!      enddo
!
!      DO J=1,M
!        DO I=2,N-1
!        divi(I,J)=rpe_1/(c11(i+1,j)+c11(i-1,j)*(rpe_1-ps(i-2,j))&
!              &  +s(J,I)*(Delta_t_I+rpe_1))
!        ps(J,I)=  c11(i+1,j)*divi(I,J)
!
!        enddo
!      enddo

!! end into precon_prep_depth
!Delta_t=1.0d0/T_step
!write(*,*) Delta_t




! 2) loop of following ADi iterations
do iteration=1,max_iter
! call cpu_time(start_time)
 ! 2.1 calculate new right hand side using old QU and R to get tilde{tilde(R)}
If(iteration==1) then

! 1) first iteration

    do K=1,Q-1,1   
      DO I=1,N
        DO J=1,M_Q
                                   
      rhs(J,I,K)=s(J,I,K)*( - R(J,I,K))
        enddo
      enddo
   enddo
   K=Q
      DO I=1,N
        DO J=1,M_R
      rhs(J,I,K)=s(J,I,K)*( - R(J,I,K))
        enddo
      enddo



else
  
    do K=1,Q-1,1   
      DO I=2,N-1
        DO J=1,M_Q
                                   
        U(J,I,K)=QU(J,I+1,K)-QU(J,I-1,K)
        enddo
      enddo
   enddo
   K=Q
      DO I=2,N-1
        DO J=1,M_R
        U(J,I,K)=QU(J,I+1,K)-QU(J,I-1,K)
        enddo
      enddo
     write(*,*) M2 
   K=M2(1)
      DO I=2,N-1
        DO J=M2(2),M_Q
                                   
        V(J,I,K)=QU(JtoJP1(J,K,2),I,JtoJP1(J,K,1))-QU(JtoJM1(J,K,2),I,JtoJM1(J,K,1))
        enddo
      enddo
    do K=M2(1)+1,MM1(1)-1,1   
      DO I=2,N-1
        DO J=1,M_Q
                                   
        V(J,I,K)=QU(JtoJP1(J,K,2),I,JtoJP1(J,K,1))-QU(JtoJM1(J,K,2),I,JtoJM1(J,K,1))
        enddo
      enddo
    enddo
   K=MM1(1)
      DO I=2,N-1
        DO J=1,MM1(2)
        V(J,I,K)=QU(JtoJP1(J,K,2),I,JtoJP1(J,K,1))-QU(JtoJM1(J,K,2),I,JtoJM1(J,K,1))
        enddo
      enddo
      
      do i=2,n-1

       U(1,I,1)=QU(1,I+1,1)-QU(1,I-1,1)  ! J=1 is J=1,K=1
       U(M_R,I,Q)=QU(M_R,I+1,Q)-QU(M_R,I-1,Q)      ! J=M is J=M_R, K=Q
       V(1,I,1)=QU(M2(2),i,M2(1))-QU(1,ip(I),1)  ! J=2 is J=M2(2), K=M2(1)
       V(M_R,I,Q)=QU(M_R,ip(I),Q)-QU(MM1(2),I,MM1(1))    ! J=M is J=M_R, K=Q; J=M-1 is J=MM1(2), J=MM1(1)
      enddo
     ! do i=2,n-1
     !  U(1,I)=QU(1,I+1)-QU(1,I-1)
     !  U(M,I)=QU(M,I+1)-QU(M,I-1)
     !  V(1,I)=QU(2,i)-QU(1,ip(I))
     !  V(M,I)=QU(M,ip(I))-QU(M-1,I)
     ! enddo
       CALL XBC_T_B_rp(U,N,Q,M_Q,M_R)
       CALL XBC_T_B_rp(V,N,Q,M_Q,M_R)
!V_dp(:,:)=V(:,:)
!U_dp(:,:)=U(:,:)
!    write(*,*) 'U', minval(abs(U_dp(2:n-1,:))), maxval(abs(U_dp(2:n-1,:)))
!    write(*,*) 'V', minval(abs(V_dp(2:n-1,:))), maxval(abs(V_dp(2:n-1,:)))

     do K=Q-1,1
      do I=1,n
        do J=1,m
         U_p(J,I,K)=                swcp*(a12(J,I,K)*V(J,I,K)+b11(J,I,K)*QU(J,I,K))
         ! V_p(J,I)=V(J,I)*a21(J,I)+swcp*(a22(J,I)*U(J,I)+b22(J,I)*QU(J,I))
        enddo
      enddo
    enddo
    K=Q
      do I=1,n
        do J=1,M_R
         U_p(J,I,K)=                swcp*(a12(J,I,K)*V(J,I,K)+b11(J,I,K)*QU(J,I,K))
         ! V_p(J,I)=V(J,I)*a21(J,I)+swcp*(a22(J,I)*U(J,I)+b22(J,I)*QU(J,I))
        enddo
      enddo

     do K=Q-1,1
      do I=1,n
        do J=1,m
       V_p(J,I,K)=V(J,I,K)*a21(J,I,K)+swcp*(a22(J,I,K)*U(J,I,K)+b22(J,I,K)*QU(J,I,K))
        enddo
      enddo
    enddo
    K=Q
      do I=1,n
        do J=1,M_R
       V_p(J,I,K)=V(J,I,K)*a21(J,I,K)+swcp*(a22(J,I,K)*U(J,I,K)+b22(J,I,K)*QU(J,I,K))
        enddo
      enddo


  CALL XBC_T_B_rp(U_p,N,Q,M_Q,M_R)
  CALL XBC_T_B_rp(V_p,N,Q,M_Q,M_R)
! call cpu_time(end_time)
! write(*,*) 'time inversion(rhs 1)', end_time-start_time
! call cpu_time(start_time)
!V_dp(:,:)=V(:,:)
!U_dp(:,:)=U(:,:)
!    write(*,*) 'U coeff', minval(abs(U_dp(2:n-1,:))), maxval(abs(U_dp(2:n-1,:)))
!    write(*,*) 'V coeff', minval(abs(V_dp(2:n-1,:))), maxval(abs(V_dp(2:n-1,:)))
! JtoJP1(J,K,2),I,JtoJP1(J,K,1)
   K=M2(1)
      DO I=2,N-1
        DO J=M2(2),M_Q
                                   
      F(J,I,K)= 0.5_HP*(U_p(J,I+1,K)-U_p(J,I-1,K)+V_p(JtoJP1(J,K,2),I,JtoJP1(J,K,1))-V_p(JtoJM1(J,K,2),I,JtoJM1(J,K,1)))/(S(J,I,K))
        enddo
      enddo
    do K=M2(1)+1,MM1(1)-1,1   
      DO I=2,N-1
        DO J=1,M_Q
                                   
      F(J,I,K)= 0.5_HP*(U_p(J,I+1,K)-U_p(J,I-1,K)+V_p(JtoJP1(J,K,2),I,JtoJP1(J,K,1))-V_p(JtoJM1(J,K,2),I,JtoJM1(J,K,1)))/(S(J,I,K))
        enddo
      enddo
    enddo
   K=MM1(1)
      DO I=2,N-1
        DO J=1,MM1(2)
      F(J,I,K)= 0.5_HP*(U_p(J,I+1,K)-U_p(J,I-1,K)+V_p(JtoJP1(J,K,2),I,JtoJP1(J,K,1))-V_p(JtoJM1(J,K,2),I,JtoJM1(J,K,1)))/(S(J,I,K))
        enddo
      enddo
!    DO I=2,N-1
!  DO J=2,M-1
!      F(J,I)= (U(J,I)-V_p(J-1,I))/(S(J,I))
!    end do
!  end do    
 

  DO I=2,N-1
    F(1,I,1)= 0.5_HP*(U_p(1,I+1,1)-U_p(1,I-1,1)+(V_p(M2(2),I,M2(1))+V_p(1,I,1)))/(S(1,I,1))  ! J=1 is j=1,K=1; J=2 is J=M2(2),K=M2(1)
    F(M_R,I,Q)= 0.5_HP*(U_p(M_R,I+1,Q)-U_p(M_R,I-1,Q)-(V_p(M_R,I,Q)+V_p(MM1(2),I,MM1(1))))/(S(M_R,I,Q)) ! J=M is j=M_R,K=Q; J=M1 is J=MM1(2),K=MM1(1)
  ENDDO


  CALL XBC_T_B_rp(F,N,Q, M_Q, M_R)
! call cpu_time(end_time)
! write(*,*) 'time inversion(rhs 2)', end_time-start_time
! call cpu_time(start_time)
!U_dp(:,:)=F(:,:)
!    write(*,*) 'F', minval(abs(U_dp(2:n-1,:))), maxval(abs(U_dp(2:n-1,:)))


     do K=Q-1,1
      do I=1,n
        do J=1,m
      rhs(J,I,K)=F(J,I,K) + s(J,I,K)*(Delta_t_I*QU(J,I,K)    &
                          &   - R(J,I,K))
        enddo
      enddo
    enddo
    K=Q
      do I=1,n
        do J=1,M_R
      rhs(J,I,K)=F(J,I,K) + s(J,I,K)*(Delta_t_I*QU(J,I,K)    &
                          &   - R(J,I,K))
        enddo
      enddo
 
  CALL XBC_T_B_rp(rhs,N,Q,M_Q, M_R)
!U_dp(:,:)=rhs(:,:)
!    write(*,*) 'RHS', minval(abs(U_dp(2:n-1,:))), maxval(abs(U_dp(2:n-1,:)))

end if

 
! call cpu_time(end_time)
! write(*,*) 'time inversion(rhs 3)', end_time-start_time
! call cpu_time(start_time)
     do K=Q-1,1
        do J=1,m
    ps(J,0,K)    =0.0_HP !rpe_0
    qs(J,0,K,0)  =0.0_HP !rpe_0

    divi(J,1,K)=(c11(J,2,K)+c11(J,N-2,K) +s(J,1,K)*(Delta_t_I+1.0_HP)) !rpe_1))
    divi(J,1,K)=1.0_HP/    divi(J,1,K)

    ps(J,1,K)    =c11(J,2,K)*divi(J,1,K)
    qs(J,1,K,0)  =rhs(J,1,K)*divi(J,1,K)

    qs(J,0,K,1)  =1.0_HP !rpe_1
    qs(J,1,K,1)  =0.0_HP!rpe_0

    qs(J,0,K,2)  =0.0_HP !rpe_0
    qs(J,1,K,2)  =c11(J,N-2,K)*divi(J,1,K)
        enddo
    enddo
    K=Q
        do J=1,M_R
    ps(J,0,K)    =0.0_HP !rpe_0
    qs(J,0,K,0)  =0.0_HP !rpe_0

    divi(J,1,K)=(c11(J,2,K)+c11(J,N-2,K) +s(J,1,K)*(Delta_t_I+1.0_HP)) !rpe_1))
    divi(J,1,K)=1.0_HP/    divi(J,1,K)

    ps(J,1,K)    =c11(J,2,K)*divi(J,1,K)
    qs(J,1,K,0)  =rhs(J,1,K)*divi(J,1,K)

    qs(J,0,K,1)  =1.0_HP !rpe_1
    qs(J,1,K,1)  =0.0_HP!rpe_0

    qs(J,0,K,2)  =0.0_HP !rpe_0
    qs(J,1,K,2)  =c11(J,N-2,K)*divi(J,1,K)
      enddo

  


    do K=1,Q-1,1   
      DO I=2,N-2
        DO J=1,M_Q
                                   
        divi(J,I,K)=1.0_HP/(c11(J,I+1,K)+c11(J,I-1,K)*(1.0_HP-ps(J,I-2,K))&
                  &  +s(J,I,K)*(Delta_t_I+1.0_HP))
        enddo
      enddo
   enddo
    do K=1,Q-1,1   
      DO I=2,N-2
        DO J=1,M_Q
                                   
          ps(J,I,K)=  c11(J,I+1,K)*divi(J,I,K)

       qs(J,I,K,0)= (rhs(J,I,K)+c11(J,I-1,K)*qs(J,I-2,K,0))*divi(J,I,K)
       qs(J,I,K,1)= (         c11(J,I-1,K)*qs(J,I-2,K,1))*divi(J,I,K) 
       qs(J,I,K,2)= (         c11(J,I-1,K)*qs(J,I-2,K,2))*divi(J,I,K) 
        enddo
      enddo
   enddo
   K=Q
      DO I=2,N-2
        DO J=1,M_R
        divi(J,I,K)=1.0_HP/(c11(J,I+1,K)+c11(J,I-1,K)*(1.0_HP-ps(J,I-2,K))&
                  &  +s(J,I,K)*(Delta_t_I+1.0_HP))
        enddo
      enddo
      DO I=2,N-2
        DO J=1,M_R
          ps(J,I,K)=  c11(J,I+1,K)*divi(J,I,K)

       qs(J,I,K,0)= (rhs(J,I,K)+c11(J,I-1,K)*qs(J,I-2,K,0))*divi(J,I,K)
       qs(J,I,K,1)= (         c11(J,I-1,K)*qs(J,I-2,K,1))*divi(J,I,K) 
       qs(J,I,K,2)= (         c11(J,I-1,K)*qs(J,I-2,K,2))*divi(J,I,K) 
        enddo
      enddo
!        DO I=2,N-2
!    DO J=1,M
!        divi(J,I)=(c11(J,I+1)+c11(J,I-1)*(1.0_HP-ps(J,I-2))&
!              &  +s(J,I)*(Delta_t_I+1.0_HP))
!        divi(J,I)=rpe_1/divi(J,I)
!        ps(J,I)=  c11(J,I+1)*divi(J,I)
!
!       qs(J,I,0)= (rhs(J,I)+c11(J,I-1)*qs(J,I-2,0))*divi(J,I)
!       qs(J,I,1)= (         c11(J,I-1)*qs(J,I-2,1))*divi(J,I) 
!       qs(J,I,2)= (         c11(J,I-1)*qs(J,I-2,2))*divi(J,I) 
!        enddo
!      enddo
    
! call cpu_time(end_time)
! write(*,*) 'time inversion(qs calc)', end_time-start_time

!  call cpu_time(start_time)
!      DO I=2,N-2
!      DO J=1,M
!        divi_T(J,I)=(c11_T(J,I+1)+c11_T(J,I-1)*(1.0_HP-ps_T(J,I-2))&
!              &  +s_T(J,I)*(Delta_t_I+1.0_HP))
!        divi_T(J,I)=1.0_HP/divi_T(J,I)
!        ps_T(J,I)=  c11_T(J,I+1)*divi_T(J,I)
!
!       qs_T(J,I,0)= (rhs_T(J,I)+c11_T(J,i-1)*qs_T(J,I-2,0))*divi_T(J,I)
!       qs_T(J,I,1)= (         c11_T(J,i-1)*qs_T(J,I-2,1))*divi_T(J,I) 
!       qs_T(J,I,2)= (         c11_T(J,i-1)*qs_T(J,I-2,2))*divi_T(J,I)  
!        enddo
!      enddo
!
!  call cpu_time(end_time)
!  write(*,*) 'time qs trans prepare)', end_time-start_time

  !  DO J=1,DP_Depth
  !      DO I=2,N-2
  !      divi(I,J)=(c11(i+1,j)+c11(i-1,j)*(1.0_HP-ps(i-2,j))&
  !            &  +s(i,j)*(Delta_t_I+1.0_HP))
  !      divi(I,J)=rpe_1/divi(I,J)
  !      ps(i,j)=  c11(i+1,j)*divi(I,J)

  !     qs(I,J,0)= (rhs(I,J)+c11(i-1,j)*qs(I-2,J,0))*divi(I,J)
  !     qs(I,J,1)= (         c11(i-1,j)*qs(I-2,J,1))*divi(I,J) 
  !     qs(I,J,2)= (         c11(i-1,j)*qs(I-2,J,2))*divi(I,J) 
  !      enddo
  !    enddo

  !    DO J=M+1-DP_Depth,M
  !      DO I=2,N-2
  !      divi(I,J)=(c11(i+1,j)+c11(i-1,j)*(1.0_HP-ps(i-2,j))&
  !            &  +s(i,j)*(Delta_t_I+1.0_HP))
  !      divi(I,J)=rpe_1/divi(I,J)
  !      ps(i,j)=  c11(i+1,j)*divi(I,J)

  !     qs(I,J,0)= (rhs(I,J)+c11(i-1,j)*qs(I-2,J,0))*divi(I,J)
  !     qs(I,J,1)= (         c11(i-1,j)*qs(I-2,J,1))*divi(I,J) 
  !     qs(I,J,2)= (         c11(i-1,j)*qs(I-2,J,2))*divi(I,J)  
  !      enddo
  !    enddo

!         divi_dp(:,:)=divi(:,:)
!         ps_dp(:,:)=ps(:,:)
!         qs_dp(:,:,0)=qs(:,:,0)
!         qs_dp(:,:,1)=qs(:,:,1)
!         qs_dp(:,:,2)=qs(:,:,2)
!    write(*,*) 'divi', minval(abs(divi_dp(:,:))), maxval(abs(divi_dp(:,:)))
!    write(*,*) 'ps', minval(abs(ps_dp(:,:))), maxval(abs(ps_dp(:,:)))
!    write(*,*) 'qs', minval(abs(qs_dp(:,:,0))), maxval(abs(qs_dp(:,:,0)))
!    write(*,*) 'qs', minval(abs(qs_dp(:,:,1))), maxval(abs(qs_dp(:,:,1)))
!    write(*,*) 'qs', minval(abs(qs_dp(:,:,2))), maxval(abs(qs_dp(:,:,2)))
!    read(*,*)

  !write(*,*) 'q finished'
  ! calculate the 5 linear subsystems with boundary conditions ( right boundary)
  
! call cpu_time(start_time)
     do K=Q-1,1
        do J=1,m
  ws(J,N-1,K,0)=0.0_HP !rpe_0
  ws(J,N-2,K,0)=qs(J,N-2,K,0)
  
  ws(J,N-1,K,1)=0.0_HP !rpe_0
  ws(J,N-2,K,1)=qs(J,N-2,K,1)
  
  ws(J,N-1,K,2)=1.0_HP !rpe_1
  ws(J,N-2,K,2)=0.0_HP !rpe_0
  
  ws(J,N-1,K,3)=0.0_HP !rpe_0
  ws(J,N-2,K,3)=qs(J,N-2,K,2)
  
  ws(J,N-1,K,4)=0.0_HP !rpe_0
  ws(J,N-2,K,4)=ps(J,N-2,K)
        enddo
    enddo
    K=Q
        do J=1,M_R
  ws(J,N-1,K,0)=0.0_HP !rpe_0
  ws(J,N-2,K,0)=qs(J,N-2,K,0)
  
  ws(J,N-1,K,1)=0.0_HP !rpe_0
  ws(J,N-2,K,1)=qs(J,N-2,K,1)
  
  ws(J,N-1,K,2)=1.0_HP !rpe_1
  ws(J,N-2,K,2)=0.0_HP !rpe_0
  
  ws(J,N-1,K,3)=0.0_HP !rpe_0
  ws(J,N-2,K,3)=qs(J,N-2,K,2)
  
  ws(J,N-1,K,4)=0.0_HP !rpe_0
  ws(J,N-2,K,4)=ps(J,N-2,K)
      enddo
!  DO J=1,M  
!  ws(J,N-1,0)=0.0_HP !rpe_0
!  ws(J,N-2,0)=qs(J,N-2,0)
!            
!  ws(J,N-1,1)=0.0_HP !rpe_0
!  ws(J,N-2,1)=qs(J,N-2,1)
!            
!  ws(J,N-1,2)=1.0_HP !rpe_1
!  ws(J,N-2,2)=0.0_HP !rpe_0
!            
!  ws(J,N-1,3)=0.0_HP !rpe_0
!  ws(J,N-2,3)=qs(J,N-2,2)
!            
!  ws(J,N-1,4)=0.0_HP !rpe_0
!  ws(J,N-2,4)=ps(J,N-2)
!   end do
  
!  DO J=1,M
!   DO I=N-1,2,-1
!     ws(I,J,0) = ps(I+2,J)*ws(I+2,J,0)+qs(I+2,J,0)
!     ws(I,J,1) = ps(I+2,J)*ws(I+2,J,1)+qs(I+2,J,1)
!     ws(I,J,2) = ps(I+2,J)*ws(I+2,J,2)+qs(I+2,J,2)
!     ws(I,J,3) = ps(I+2,J)*ws(I+2,J,3)+qs(I+2,J,3)
!     ws(I,J,4) = ps(I+2,J)*ws(I+2,J,4)+qs(I+2,J,4)
!    end do
!  end do 

  !! rewritten to change precision at poles as desired 
    do K=1,Q-1,1   
      DO I=N-3,1,-1
        DO J=1,M_Q
                                   
     ws(J,I,K,0) = ps(J,I,K)*ws(J,I+2,K,0)+qs(J,I,K,0)
     ws(J,I,K,1) = ps(J,I,K)*ws(J,I+2,K,1)+qs(J,I,K,1)
     ws(J,I,K,2) = ps(J,I,K)*ws(J,I+2,K,2)           
     ws(J,I,K,3) = ps(J,I,K)*ws(J,I+2,K,3)+qs(J,I,K,2)
     ws(J,I,K,4) = ps(J,I,K)*ws(J,I+2,K,4)
        enddo
      enddo
   enddo
   K=Q
      DO I=N-3,1,-1
        DO J=1,M_R
     ws(J,I,K,0) = ps(J,I,K)*ws(J,I+2,K,0)+qs(J,I,K,0)
     ws(J,I,K,1) = ps(J,I,K)*ws(J,I+2,K,1)+qs(J,I,K,1)
     ws(J,I,K,2) = ps(J,I,K)*ws(J,I+2,K,2)           
     ws(J,I,K,3) = ps(J,I,K)*ws(J,I+2,K,3)+qs(J,I,K,2)
     ws(J,I,K,4) = ps(J,I,K)*ws(J,I+2,K,4)
        enddo
      enddo


! call cpu_time(end_time)
!
! write(*,*) 'time inversion(ws prepare)', end_time-start_time

!     
!      call cpu_time(start_time)
!      DO I=N-3,1,-1
!      DO J=1,M
!                                   
!     ws(J,I,0) = ps(J,I)*ws(J,I+2,0)+qs(J,I,0)
!     ws(J,I,1) = ps(J,I)*ws(J,I+2,1)+qs(J,I,1)
!     ws(J,I,2) = ps(J,I)*ws(J,I+2,2)
!     ws(J,I,3) = ps(J,I)*ws(J,I+2,3)+qs(J,I,2)
!     ws(J,I,4) = ps(J,I)*ws(J,I+2,4)
!        enddo
!      enddo
!
! call cpu_time(end_time)
! write(*,*) 'time inversion(ws prepare) original', end_time-start_time
!      call cpu_time(start_time)
!      DO I=N-3,1,-1
!      DO J=1,M
!                                   
!     ws(J,I,2) = ps(J,I)*ws(J,I+2,2)
!        enddo
!      enddo
!        DO I=N-3,1,-1
!      DO J=1,M
!                                   
!     ws(J,I,4) = ps(J,I)*ws(J,I+2,4)
!        enddo
!      enddo
!        DO I=N-3,1,-1
!      DO J=1,M
!
!     ws(J,I,0) = ps(J,I)*ws(J,I+2,0)+qs(J,I,0)
!        enddo
!      enddo
!        DO I=N-3,1,-1
!      DO J=1,M
!                                   
!     ws(J,I,1) = ps(J,I)*ws(J,I+2,1)+qs(J,I,1)
!        enddo
!      enddo
!        DO I=N-3,1,-1
!      DO J=1,M
!
!     ws(J,I,3) = ps(J,I)*ws(J,I+2,3)+qs(J,I,2)
!        enddo
!      enddo


!         ws_dp(:,:,0)=ws(:,:,0)
!         ws_dp(:,:,1)=ws(:,:,1)
!         ws_dp(:,:,2)=ws(:,:,2)
!         ws_dp(:,:,3)=ws(:,:,3)
!         ws_dp(:,:,4)=ws(:,:,4)
!    write(*,*) 'ws', minval(abs(ws_dp(:,:,0))), maxval(abs(ws_dp(:,:,0)))
!    write(*,*) 'ws', minval(abs(ws_dp(:,:,1))), maxval(abs(ws_dp(:,:,1)))
!    write(*,*) 'ws', minval(abs(ws_dp(:,:,2))), maxval(abs(ws_dp(:,:,2)))
!    write(*,*) 'ws', minval(abs(ws_dp(:,:,3))), maxval(abs(ws_dp(:,:,3)))
!    write(*,*) 'ws', minval(abs(ws_dp(:,:,4))), maxval(abs(ws_dp(:,:,4)))
!    read(*,*)
   !! end of rewrite

  ! write(*,*) 'w finished'
  ! solve the subsystems for the coefficients a,b,g,d for each latitude
  il1=N-2
  il2=N-3
  i2 =2
! call cpu_time(end_time)
! write(*,*) 'time inversion(ws prepare) split', end_time-start_time
!
! call cpu_time(start_time)
 Do K=1,Q
  DO J=1,M_Q
    if (.not.(K==Q .AND. J>M_R)) then
      
      
       d11=ws(J,il1,K,1)-rpe_1
       d12=ws(J,il1,K,2)
       d13=ws(J,il1,K,3)
       d14=ws(J,il1,K,4)
       s1=-ws(J,il1,K,0)
       
       d21=ws(J,1  ,K,1)
       d22=ws(J,1  ,K,2)-rpe_1
       d23=ws(J,1  ,K,3)
       d24=ws(J,1  ,K,4)
       s2=-ws(J,1  ,K,0)
       
       d31=ws(J,il2,K,1)
       d32=ws(J,il2,K,2)
       d33=ws(J,il2,K,3)-rpe_1
       d34=ws(J,il2,K,4)
       s3=-ws(J,il2,K,0)
       
       d41=ws(J,i2 ,K,1)
       d42=ws(J,i2 ,K,2)
       d43=ws(J,i2 ,K,3)
       d44=ws(J,i2 ,K,4)-rpe_1
       s4=-ws(J,i2 ,K,0)

      det40=d11*& !det3_ARM_RP(d22,d23,d24,d32,d33,d34,d42,d43,d44) &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
        &(  d22*d33*d44+d23*d34*d42+d24*d32*d43    &
            &     -d42*d33*d24-d43*d34*d22-d44*d32*d23) &
       &  -d21* &!det3_ARM_RP(d12,d13,d14,d32,d33,d34,d42,d43,d44)  &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
        & ( d12*d33*d44+d13*d34*d42+d14*d32*d43    &
            &     -d42*d33*d14-d43*d34*d12-d44*d32*d32) &
       &  +d31*& !det3_ARM_RP(d12,d13,d14,d22,d23,d24,d42,d43,d44)  &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
        &  ( d12*d23*d44+d13*d24*d42+d14*d22*d43    &
            &     -d14*d23*d14-d43*d24*d12-d44*d22*d13) &
       &  -d41* & !det3_ARM_RP(d12,d13,d14,d22,d23,d24,d32,d33,d34) 
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
        & ( d12*d23*d34+d13*d24*d32+d14*d22*d33    &
            &     -d32*d23*d14-d33*d24*d12-d34*d22*d13)
     deti=rpe_1/det40
     det41=s1 *& !det3_ARM_RP(d22,d23,d24,d32,d33,d34,d42,d43,d44) &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
         & ( d22*d33*d44+d23*d34*d42+d24*d32*d43    &
         &     -d42*d33*d24-d43*d34*d22-d44*d32*d23 )&
       &  -s2 *& !det3_ARM_RP(d12,d13,d14,d32,d33,d34,d42,d43,d44)  &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
          & ( d12*d33*d44+d13*d34*d42+d14*d32*d43    &
            &     -d42*d33*d14-d43*d34*d12-d44*d32*d32) &
       &  +s3 *& !det3_ARM_RP(d12,d13,d14,d22,d23,d24,d42,d43,d44)  &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
          &( d12*d23*d44+d13*d24*d42+d14*d22*d43    &
            &     -d14*d23*d14-d43*d24*d12-d44*d22*d13)&
       &  -s4 *& !det3_ARM_RP(d12,d13,d14,d22,d23,d24,d32,d33,d34)  
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
          &( d12*d23*d34+d13*d24*d32+d14*d22*d33    &
            &     -d32*d23*d14-d33*d24*d12-d34*d22*d13)

     det42=d11*& !det3_ARM_RP( s2,d23,d24, s3,d33,d34, s4,d43,d44) &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
           &( s2*d33*d44+d23*d34*s4+d24*s3*d43    &
            &     -s4*d33*d24-d43*d34*s2-d44*s3*d23)&
       &  -d21*& !det3_ARM_RP( s1,d13,d14, s3,d33,d34, s4,d43,d44)  &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
          &(  s1*d33*d44+d13*d34*s4+d14*s3*d43    &
            &     -s4*d33*d14-d43*d34*s1-d44*s3*d13)&
       &  +d31*& !det3_ARM_RP( s1,d13,d14, s2,d23,d24, s4,d43,d44)  &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
          &(  s1*d23*d44+d13*d24*s4+d14*s2*d43    &
            &     -s4*d23*d14-d43*d24*s1-d44*s2*d13)&
       &  -d41*& !det3_ARM_RP( s1,d13,d14, s2,d23,d24, s3,d33,d34)  
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
          &(  s1*d23*d34+d13*d24*s3+d14*s2*d33    &
            &     -s3*d23*d14-d33*d24*s1-d34*s2*d13)
 
     det43=d11*& !det3_ARM_RP(d22, s2,d24,d32, s3,d34,d42, s4,d44) &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
           &( d22*s3*d44+s2*d34*d42+d24*d32*d34    &
            &     -d42*s3*d24-s4*d34*d22-d44*d32*s2) &
       &  -d21*& !det3_ARM_RP(d12, s1,d14,d32, s3,d34,d42, s4,d44)  &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
          &(  d12*s3*d44+s1*d34*d42+d14*d32*s4    &
            &     -d42*s3*d14-s4*d34*d12-d44*d32*s1)&
       &  +d31*& !det3_ARM_RP(d12, s1,d14,d22, s2,d24,d42, s4,d44)  &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
          &(  d12*s2*d44+s1*d24*d42+d14*d22*s4    &
            &     -d42*s2*d14-s4*d24*d12-d44*d22*s1)&
       &  -d41*& !det3_ARM_RP(d12, s1,d14,d22, s2,d24,d32, s3,d34)
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
         &(  d12*s2*d34+s1*d24*d32+d14*d22*s3    &
            &     -d32*s2*d14-s3*d24*d12-d34*d22*s1)

     det44=d11*& !det3_ARM_RP(d22,d23, s2,d32,d33, s3,d42,d43, s4) &
     ! function   det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
         & (  d22*d33*s4+d23*s3*d42+s2*d32*d43    &
            &     -d42*d33*s2-d43*s3*d22-s4*d32*d23) &
       &  -d21*& !det3_ARM_RP(d12,d13, s1,d32,d33, s3,d42,d43, s4)  &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
         & (   d12*d33*s4+d13*s3*d42+s1*d32*d43    &
            &     -d42*d33*s1-d43*s3*d12-s4*d32*d13) & 
       &  +d31*& !det3_ARM_RP(d12,d13, s1,d22,d23, s2,d42,d43, s4)  &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
           &( d12*d23*s4+d13*s2*d42+s1*d22*d43    &
            &     -d42*d23*s1-d43*s2*d12-s4*d22*d13 )&
       &  -d41*& !det3_ARM_RP(d12,d13, s1,d22,d23, s2,d32,d33, s3)
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
           &( d12*d23*s3+d13*s2*d32+s1*d22*d33    &
            &     -d32*d23*s1-d33*s2*d12-s3*d22*d13)
       
      aa(J,K,4)=det44*deti
      aa(J,K,3)=det43*deti
      aa(J,K,2)=det42*deti
      aa(J,K,1)=det41*deti
    end if
  enddo 
enddo
  !write(*,*) 'aa's finished'
!  DO J=1,M
!   DO I=2,N-1
!      QU(I,J)=ws(I,J,0) +aa(J,1)*ws(I,J,1) &
!                    &   +aa(J,2)*ws(I,J,2)  &
!                    &   +aa(J,3)*ws(I,J,3)  &
!                    &   +aa(J,4)*ws(I,J,4)
!  enddo
!   write(*,*) J, (QU(I,J), I=1,N)
!read(*,*)
!  enddo

  !! rewritten to change precision at poles as desired
   
! call cpu_time(end_time)
! write(*,*) 'time inversion', end_time-start_time
!
! call cpu_time(start_time)
    do K=1,Q-1,1   
      DO I=2,N-1
        DO J=1,M_Q
                                   
      QU(J,I,K)=ws(J,I,K,0) +aa(J,K,1)*ws(J,I,K,1) &
                        &   +aa(J,K,2)*ws(J,I,K,2)  &
                        &   +aa(J,K,3)*ws(J,I,K,3)  &
                        &   +aa(J,K,4)*ws(J,I,K,4)
        enddo
      enddo
   enddo
   K=Q
      DO I=2,N-1
        DO J=1,M_R
      QU(J,I,K)=ws(J,I,K,0) +aa(J,K,1)*ws(J,I,K,1) &
                        &   +aa(J,K,2)*ws(J,I,K,2)  &
                        &   +aa(J,K,3)*ws(J,I,K,3)  &
                        &   +aa(J,K,4)*ws(J,I,K,4)
        enddo
      enddo

  CALL XBC_T_B_rp(QU,N,Q,M_Q, M_R)
 !     QU_dp(:,:)=QU(:,:)
 !   write(*,*) 'QU XBC', minval(abs(QU_dp(:,:))), maxval(abs(Qu_dp(:,:)))
 !   read(*,*)
     
! call cpu_time(end_time)
! write(*,*) 'time construct QU', end_time-start_time
! 
! call cpu_time(start_time)
  call adjust_conservation_Trans_Block_ARM_RP(QU,s,S_full,n,Q,M_Q, M_R)
! call cpu_time(end_time)
! write(*,*) 'time mass fixer', end_time-start_time
  !    QU_dp(:,:)=QU(:,:)
  !    write(*,*) 'S_full',S_full
  !  write(*,*) 'QU adjust', minval(abs(QU_dp(:,:))), maxval(abs(Qu_dp(:,:)))
  !  read(*,*)
   !! end of rewrite
   
  !write(*,*) 'QU finished'






end do
  !write(*,*) 'BC QU finished'
END SUBROUTINE

SUBROUTINE precon_ADI_Piotr_Trans_ARM_RP_noInput(N,M,IP)
!use implicit_functions_HP

implicit none
integer, parameter :: HP=2
REAL(kind=4)::  S_full, start_time, end_time, R_sp(M,N), scaling
REAL(Kind=HP) :: S(M,N)
REAL(Kind=HP) :: R(M,N),QU(M,N),A11(M,N),A12(M,N),A21(M,N),A22(M,N),  &
    &      B11(M,N),B22(M,N),P0(M,N),U(M,N),V(M,N),  C11(M,N)
INTEGER :: IP(N), ID_PREC, num_of_bits,DP_Depth
INTEGER :: N, M

REAL(Kind=HP) :: max_QX_QY, F(M,N), rhs(M,N), qs(M,0:N-1,0:2), ps(M,0:N-2), &
        & ws(M,N-1,0:4), A(M,N), B(M,N), C(M,N), divi(M,N)&
        & , help
REAL(Kind=HP) :: aa(M,1:4), deti, det40, det41, det42, det43, det44, det3_ARM_RP,   &
                           & d11, d12, d13, d14, d21, d22, d23, d24,       &
                           & d31, d32, d33, d34, d41, d42, d43, d44,       &
                           & s1, s2, s3, s4
REAL(Kind=HP) ::T_step, Delta_t_I, dn, dni, util, vtil, swcp
REAL(kind=HP) :: RPE_0,RPE_1,RPE_05,RPE_025,RPE_0125,RPE_2,RPE_3,RPE_4 , rpe_100
REAL(kind=HP) :: U_p(M,N),V_p(M,N)
!REAL(Kind=4) :: 
integer :: iter, max_iter, time_scale  !! number of richardson iterations
INTEGER :: I, J, K, iteration, i2, il1, il2, npromaw

!REAL(Kind=HP) :: S_T(M,N), C11_T(M,N), ps_T(M,0:N-2), qs_T(M,0:N-1,0:2),divi_T(M,N)&
!        & , rhs_T(M,N)

DO I=1,N
  DO J=1,M
   r(J,I)=r_sp(J,I)*scaling
  end do
end do
    rpe_1 = (1.0)
    rpe_2 = (2.0)
    rpe_3 = (3.0)
    rpe_4 = (4.0)
    rpe_100 = (100.0)
    rpe_0125 = (0.125)
    rpe_025 = (0.25)
    rpe_05 = (0.5)
    rpe_0 = (0.0)
!call rpenum_init(10)
!write(*,*) 'Precision', HP 
max_iter= 2
swcp=rpe_1
!initialize the inverse of R with 0
  DO I=1,N
DO J=1,M
    QU(J,I) =rpe_0
    rhs(J,I)=rpe_0
  end do
end do

!! into precon_prep_depth
    DO I=1,N    
  DO J=1,M
 
      C11(J,I)=0.5_HP*a11(J,I)
    END DO
  END DO

Delta_t_I=1.0_HP/(0.25_HP/T_step)

!T_step_dp=T_step
!write(*,*) 'T_step,', T_step_dp
!write(*,*) 'T_step,', float(rpe_025)
!write(*,*) 'T_step,', float(rpe_025/T_step)
!
!Delta_t_I_dp=Delta_t_I
!write(*,*) 'Delta_t_I,', Delta_t_I_dp
!read(*,*)
!
!      DO J=1,M
!
!        ps(0,j)=rpe_0
!        divi(1,J)=rpe_1/(c11(2,j)+c11(N-2,j) +s(1,j)*(Delta_t_I+rpe_1))
!        ps(1,j) = c11(2,j)*divi(1,J)
!
!      enddo
!
!      DO J=1,M
!        DO I=2,N-1
!        divi(I,J)=rpe_1/(c11(i+1,j)+c11(i-1,j)*(rpe_1-ps(i-2,j))&
!              &  +s(J,I)*(Delta_t_I+rpe_1))
!        ps(J,I)=  c11(i+1,j)*divi(I,J)
!
!        enddo
!      enddo

!! end into precon_prep_depth
!Delta_t=1.0d0/T_step
!write(*,*) Delta_t




! 2) loop of following ADi iterations
do iteration=1,max_iter
!call cpu_time(start_time)
 ! 2.1 calculate new right hand side using old QU and R to get tilde{tilde(R)}
If(iteration==1) then

! 1) first iteration


        DO I=1,N
      DO J=1,M
      rhs(J,I)=s(J,I)*( - R(J,I))
        enddo
      enddo



else
  
       do I=2,n-1
      do J=2,m-1
        U(J,I)=QU(J,I+1)-QU(J,I-1)
        V(J,I)=QU(J+1,I)-QU(J-1,I)
       enddo
      enddo
      do i=2,n-1
       U(1,I)=QU(1,I+1)-QU(1,I-1)
       U(M,I)=QU(M,I+1)-QU(M,I-1)
       V(1,I)=QU(2,i)-QU(1,ip(I))
       V(M,I)=QU(M,ip(I))-QU(M-1,I)
      enddo
       CALL XBC_T_rp(U,M,N)
       CALL XBC_T_rp(V,M,N)
!V_dp(:,:)=V(:,:)
!U_dp(:,:)=U(:,:)
!    write(*,*) 'U', minval(abs(U_dp(2:n-1,:))), maxval(abs(U_dp(2:n-1,:)))
!    write(*,*) 'V', minval(abs(V_dp(2:n-1,:))), maxval(abs(V_dp(2:n-1,:)))

      do I=1,n
      do J=1,m
       U_p(J,I)=                swcp*(a12(J,I)*V(J,I)+b11(J,I)*QU(J,I))
      ! V_p(J,I)=V(J,I)*a21(J,I)+swcp*(a22(J,I)*U(J,I)+b22(J,I)*QU(J,I))
      enddo
      enddo
      do I=1,n
      do J=1,m
      ! U_p(J,I)=                swcp*(a12(J,I)*V(J,I)+b11(J,I)*QU(J,I))
       V_p(J,I)=V(J,I)*a21(J,I)+swcp*(a22(J,I)*U(J,I)+b22(J,I)*QU(J,I))
      enddo
      enddo

  CALL XBC_T_rp(U_p,M,N)
  CALL XBC_T_rp(V_p,M,N)
! call cpu_time(end_time)
! write(*,*) 'time inversion(rhs 1)', end_time-start_time
! call cpu_time(start_time)

!V_dp(:,:)=V(:,:)
!U_dp(:,:)=U(:,:)
!    write(*,*) 'U coeff', minval(abs(U_dp(2:n-1,:))), maxval(abs(U_dp(2:n-1,:)))
!    write(*,*) 'V coeff', minval(abs(V_dp(2:n-1,:))), maxval(abs(V_dp(2:n-1,:)))

    DO I=2,N-1
  DO J=2,M-1
      F(J,I)= (U_p(J,I+1)-U_p(J,I-1)+V_p(J+1,I)-V_p(J-1,I))/(S(J,I))
      !U(J,I)= U_p(J,I+1)-U_p(J,I-1)+V_p(J+1,I)
    end do
  end do    
!    DO I=2,N-1
!  DO J=2,M-1
!      F(J,I)= (U(J,I)-V_p(J-1,I))/(S(J,I))
!    end do
!  end do    
 

  DO I=2,N-1
    F(1,I)= (U_p(1,I+1)-U_p(1,I-1)+(V_p(2,I)+V_p(1,I)))/(S(1,I))
    F(M,I)= (U_p(M,I+1)-U_p(M,I-1)-(V_p(M,I)+V_p(M-1,I)))/(S(M,I))
  ENDDO


  CALL XBC_T_rp(F,M,N)
! call cpu_time(end_time)
! write(*,*) 'time inversion(rhs 2)', end_time-start_time
! call cpu_time(start_time)
!U_dp(:,:)=F(:,:)
!    write(*,*) 'F', minval(abs(U_dp(2:n-1,:))), maxval(abs(U_dp(2:n-1,:)))

      do j=1,m
      do i=1,n
       F(J,I)=0.5_HP*F(J,I) !-p(J,I)
      enddo
      enddo

  DO J=1,M
    DO I=1,N    
 
      rhs(J,I)=F(J,I) + s(J,I)*(Delta_t_I*QU(J,I)    &
                          &   - R(J,I))
    END DO
  END DO
 CALL XBC_T_rp(rhs,M,N)
!U_dp(:,:)=rhs(:,:)
!    write(*,*) 'RHS', minval(abs(U_dp(2:n-1,:))), maxval(abs(U_dp(2:n-1,:)))

end if

 
! call cpu_time(end_time)
! write(*,*) 'time inversion(rhs 3)', end_time-start_time
! call cpu_time(start_time)
  DO J=1,M

    ps(J,0)    =0.0_HP !rpe_0
    qs(J,0,0)  =0.0_HP !rpe_0

    divi(J,1)=(c11(J,2)+c11(J,N-2) +s(J,1)*(Delta_t_I+1.0_HP)) !rpe_1))
    divi(J,1)=1.0_HP/    divi(J,1)

    ps(J,1)    =c11(J,2)*divi(J,1)
    qs(J,1,0)  =rhs(J,1)*divi(J,1)

    qs(J,0,1)  =1.0_HP !rpe_1
    qs(J,1,1)  =0.0_HP!rpe_0

    qs(J,0,2)  =0.0_HP !rpe_0
    qs(J,1,2)  =c11(J,N-2)*divi(J,1)

  end do


  

!        DO I=2,N-2
!      DO J=1,M
!        divi(J,I)=(c11(J,I+1)+c11(J,I-1)*(1.0_HP-ps(J,I-2))&
!              &  +s(J,I)*(Delta_t_I+1.0_HP)) !rpe_1))
!        divi(J,I)=1.0_HP/divi(J,I)
!
!        enddo
!      enddo

! call cpu_time(end_time)
! write(*,*) 'time inversion(qs prepare)', end_time-start_time
! call cpu_time(start_time)
!        DO I=2,N-2
!      DO J=1,M
!        ps(J,I)=  c11(J,I+1)*divi(J,I)
!
!        enddo
!      enddo
!
! call cpu_time(end_time)
! write(*,*) 'time inversion(qs1 prepare)', end_time-start_time
! call cpu_time(start_time)
!      !help=qs(0,1,0)
!        DO I=2,N-2
!      DO J=1,M
!       qs(J,I,0)= (rhs(J,I)+c11(J,I-1)*2.0_HP)*divi(J,I)
!       !help=qs(I,J,0)
!        enddo
!      enddo
!     ! !help=qs(1,1,0)
!     ! DO J=1,M
!     !   DO I=3,N-2,2
!     !  qs(I,J,0)= (rhs(I,J)+c11(i-1,j)*2.0_HP)*divi(I,J)
!     !  !help=qs(I,J,0)
!     !
!     !   enddo
!     ! enddo
! call cpu_time(end_time)
! write(*,*) 'time inversion(qs1.1 prepare)', end_time-start_time
! call cpu_time(start_time)
!        DO I=2,N-2
!      DO J=1,M
!       qs(J,I,1)= (         c11(J,I-1)*2.0_HP)*divi(J,I) 
!
!        enddo
!      enddo
! call cpu_time(end_time)
! write(*,*) 'time inversion(qs1.2 prepare)', end_time-start_time

! call cpu_time(start_time)
        DO I=2,N-2
    DO J=1,M
        divi(J,I)=(c11(J,I+1)+c11(J,I-1)*(1.0_HP-ps(J,I-2))&
              &  +s(J,I)*(Delta_t_I+1.0_HP))
        divi(J,I)=rpe_1/divi(J,I)
        ps(J,I)=  c11(J,I+1)*divi(J,I)

       qs(J,I,0)= (rhs(J,I)+c11(J,I-1)*qs(J,I-2,0))*divi(J,I)
       qs(J,I,1)= (         c11(J,I-1)*qs(J,I-2,1))*divi(J,I) 
       qs(J,I,2)= (         c11(J,I-1)*qs(J,I-2,2))*divi(J,I) 
        enddo
      enddo
    
! call cpu_time(end_time)
! write(*,*) 'time inversion(qs calc)', end_time-start_time

!  call cpu_time(start_time)
!      DO I=2,N-2
!      DO J=1,M
!        divi_T(J,I)=(c11_T(J,I+1)+c11_T(J,I-1)*(1.0_HP-ps_T(J,I-2))&
!              &  +s_T(J,I)*(Delta_t_I+1.0_HP))
!        divi_T(J,I)=1.0_HP/divi_T(J,I)
!        ps_T(J,I)=  c11_T(J,I+1)*divi_T(J,I)
!
!       qs_T(J,I,0)= (rhs_T(J,I)+c11_T(J,i-1)*qs_T(J,I-2,0))*divi_T(J,I)
!       qs_T(J,I,1)= (         c11_T(J,i-1)*qs_T(J,I-2,1))*divi_T(J,I) 
!       qs_T(J,I,2)= (         c11_T(J,i-1)*qs_T(J,I-2,2))*divi_T(J,I)  
!        enddo
!      enddo
!
!  call cpu_time(end_time)
!  write(*,*) 'time qs trans prepare)', end_time-start_time

! call cpu_time(start_time)
  !  DO J=1,DP_Depth
  !      DO I=2,N-2
  !      divi(I,J)=(c11(i+1,j)+c11(i-1,j)*(1.0_HP-ps(i-2,j))&
  !            &  +s(i,j)*(Delta_t_I+1.0_HP))
  !      divi(I,J)=rpe_1/divi(I,J)
  !      ps(i,j)=  c11(i+1,j)*divi(I,J)

  !     qs(I,J,0)= (rhs(I,J)+c11(i-1,j)*qs(I-2,J,0))*divi(I,J)
  !     qs(I,J,1)= (         c11(i-1,j)*qs(I-2,J,1))*divi(I,J) 
  !     qs(I,J,2)= (         c11(i-1,j)*qs(I-2,J,2))*divi(I,J) 
  !      enddo
  !    enddo

  !    DO J=M+1-DP_Depth,M
  !      DO I=2,N-2
  !      divi(I,J)=(c11(i+1,j)+c11(i-1,j)*(1.0_HP-ps(i-2,j))&
  !            &  +s(i,j)*(Delta_t_I+1.0_HP))
  !      divi(I,J)=rpe_1/divi(I,J)
  !      ps(i,j)=  c11(i+1,j)*divi(I,J)

  !     qs(I,J,0)= (rhs(I,J)+c11(i-1,j)*qs(I-2,J,0))*divi(I,J)
  !     qs(I,J,1)= (         c11(i-1,j)*qs(I-2,J,1))*divi(I,J) 
  !     qs(I,J,2)= (         c11(i-1,j)*qs(I-2,J,2))*divi(I,J)  
  !      enddo
  !    enddo

!         divi_dp(:,:)=divi(:,:)
!         ps_dp(:,:)=ps(:,:)
!         qs_dp(:,:,0)=qs(:,:,0)
!         qs_dp(:,:,1)=qs(:,:,1)
!         qs_dp(:,:,2)=qs(:,:,2)
!    write(*,*) 'divi', minval(abs(divi_dp(:,:))), maxval(abs(divi_dp(:,:)))
!    write(*,*) 'ps', minval(abs(ps_dp(:,:))), maxval(abs(ps_dp(:,:)))
!    write(*,*) 'qs', minval(abs(qs_dp(:,:,0))), maxval(abs(qs_dp(:,:,0)))
!    write(*,*) 'qs', minval(abs(qs_dp(:,:,1))), maxval(abs(qs_dp(:,:,1)))
!    write(*,*) 'qs', minval(abs(qs_dp(:,:,2))), maxval(abs(qs_dp(:,:,2)))
!    read(*,*)

  !write(*,*) 'q finished'
  ! calculate the 5 linear subsystems with boundary conditions ( right boundary)
  
!      call cpu_time(start_time)
  DO J=1,M  
  ws(J,N-1,0)=0.0_HP !rpe_0
  ws(J,N-2,0)=qs(J,N-2,0)
            
  ws(J,N-1,1)=0.0_HP !rpe_0
  ws(J,N-2,1)=qs(J,N-2,1)
            
  ws(J,N-1,2)=1.0_HP !rpe_1
  ws(J,N-2,2)=0.0_HP !rpe_0
            
  ws(J,N-1,3)=0.0_HP !rpe_0
  ws(J,N-2,3)=qs(J,N-2,2)
            
  ws(J,N-1,4)=0.0_HP !rpe_0
  ws(J,N-2,4)=ps(J,N-2)
   end do
  
!  DO J=1,M
!   DO I=N-1,2,-1
!     ws(I,J,0) = ps(I+2,J)*ws(I+2,J,0)+qs(I+2,J,0)
!     ws(I,J,1) = ps(I+2,J)*ws(I+2,J,1)+qs(I+2,J,1)
!     ws(I,J,2) = ps(I+2,J)*ws(I+2,J,2)+qs(I+2,J,2)
!     ws(I,J,3) = ps(I+2,J)*ws(I+2,J,3)+qs(I+2,J,3)
!     ws(I,J,4) = ps(I+2,J)*ws(I+2,J,4)+qs(I+2,J,4)
!    end do
!  end do 

  !! rewritten to change precision at poles as desired 

     
!      DO I=N-3,1,-1
!      DO J=1,M
!                                   
!     ws(J,I,0) = ps(J,I)*ws(J,I+2,0)+qs(J,I,0)
!     ws(J,I,1) = ps(J,I)*ws(J,I+2,1)+qs(J,I,1)
!     ws(J,I,2) = ps(J,I)*ws(J,I+2,2)
!     ws(J,I,3) = ps(J,I)*ws(J,I+2,3)+qs(J,I,2)
!     ws(J,I,4) = ps(J,I)*ws(J,I+2,4)
!        enddo
!      enddo
!
! call cpu_time(end_time)
! write(*,*) 'time inversion(ws prepare) original', end_time-start_time
      ! call cpu_time(start_time)
      DO I=N-3,1,-1
      DO J=1,M
                                   
     ws(J,I,2) = ps(J,I)*ws(J,I+2,2)
        enddo
      enddo
        DO I=N-3,1,-1
      DO J=1,M
                                   
     ws(J,I,4) = ps(J,I)*ws(J,I+2,4)
        enddo
      enddo
        DO I=N-3,1,-1
      DO J=1,M

     ws(J,I,0) = ps(J,I)*ws(J,I+2,0)+qs(J,I,0)
        enddo
      enddo
        DO I=N-3,1,-1
      DO J=1,M
                                   
     ws(J,I,1) = ps(J,I)*ws(J,I+2,1)+qs(J,I,1)
        enddo
      enddo
        DO I=N-3,1,-1
      DO J=1,M

     ws(J,I,3) = ps(J,I)*ws(J,I+2,3)+qs(J,I,2)
        enddo
      enddo


!         ws_dp(:,:,0)=ws(:,:,0)
!         ws_dp(:,:,1)=ws(:,:,1)
!         ws_dp(:,:,2)=ws(:,:,2)
!         ws_dp(:,:,3)=ws(:,:,3)
!         ws_dp(:,:,4)=ws(:,:,4)
!    write(*,*) 'ws', minval(abs(ws_dp(:,:,0))), maxval(abs(ws_dp(:,:,0)))
!    write(*,*) 'ws', minval(abs(ws_dp(:,:,1))), maxval(abs(ws_dp(:,:,1)))
!    write(*,*) 'ws', minval(abs(ws_dp(:,:,2))), maxval(abs(ws_dp(:,:,2)))
!    write(*,*) 'ws', minval(abs(ws_dp(:,:,3))), maxval(abs(ws_dp(:,:,3)))
!    write(*,*) 'ws', minval(abs(ws_dp(:,:,4))), maxval(abs(ws_dp(:,:,4)))
!    read(*,*)
   !! end of rewrite

  ! write(*,*) 'w finished'
  ! solve the subsystems for the coefficients a,b,g,d for each latitude
  il1=N-2
  il2=N-3
  i2 =2
! call cpu_time(end_time)
! write(*,*) 'time inversion(ws prepare) split', end_time-start_time
!
! call cpu_time(start_time)

 DO J=1,M
    if (J<=DP_depth .or. J>=M+1-DP_depth) then
      
    end if
      
       d11=ws(J,il1,1)-rpe_1
       d12=ws(J,il1,2)
       d13=ws(J,il1,3)
       d14=ws(J,il1,4)
       s1=-ws(J,il1,0)
                     
       d21=ws(J,1  ,1)
       d22=ws(J,1  ,2)-rpe_1
       d23=ws(J,1  ,3)
       d24=ws(J,1  ,4)
       s2=-ws(J,1  ,0)
                     
       d31=ws(J,il2,1)
       d32=ws(J,il2,2)
       d33=ws(J,il2,3)-rpe_1
       d34=ws(J,il2,4)
       s3=-ws(J,il2,0)
                     
       d41=ws(J,i2 ,1)
       d42=ws(J,i2 ,2)
       d43=ws(J,i2 ,3)
       d44=ws(J,i2 ,4)-rpe_1
       s4=-ws(J,i2 ,0)

      det40=d11*& !det3_ARM_RP(d22,d23,d24,d32,d33,d34,d42,d43,d44) &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
        &(  d22*d33*d44+d23*d34*d42+d24*d32*d43    &
            &     -d42*d33*d24-d43*d34*d22-d44*d32*d23) &
       &  -d21* &!det3_ARM_RP(d12,d13,d14,d32,d33,d34,d42,d43,d44)  &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
        & ( d12*d33*d44+d13*d34*d42+d14*d32*d43    &
            &     -d42*d33*d14-d43*d34*d12-d44*d32*d32) &
       &  +d31*& !det3_ARM_RP(d12,d13,d14,d22,d23,d24,d42,d43,d44)  &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
        &  ( d12*d23*d44+d13*d24*d42+d14*d22*d43    &
            &     -d14*d23*d14-d43*d24*d12-d44*d22*d13) &
       &  -d41* & !det3_ARM_RP(d12,d13,d14,d22,d23,d24,d32,d33,d34) 
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
        & ( d12*d23*d34+d13*d24*d32+d14*d22*d33    &
            &     -d32*d23*d14-d33*d24*d12-d34*d22*d13)
     deti=rpe_1/det40
     det41=s1 *& !det3_ARM_RP(d22,d23,d24,d32,d33,d34,d42,d43,d44) &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
         & ( d22*d33*d44+d23*d34*d42+d24*d32*d43    &
         &     -d42*d33*d24-d43*d34*d22-d44*d32*d23 )&
       &  -s2 *& !det3_ARM_RP(d12,d13,d14,d32,d33,d34,d42,d43,d44)  &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
          & ( d12*d33*d44+d13*d34*d42+d14*d32*d43    &
            &     -d42*d33*d14-d43*d34*d12-d44*d32*d32) &
       &  +s3 *& !det3_ARM_RP(d12,d13,d14,d22,d23,d24,d42,d43,d44)  &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
          &( d12*d23*d44+d13*d24*d42+d14*d22*d43    &
            &     -d14*d23*d14-d43*d24*d12-d44*d22*d13)&
       &  -s4 *& !det3_ARM_RP(d12,d13,d14,d22,d23,d24,d32,d33,d34)  
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
          &( d12*d23*d34+d13*d24*d32+d14*d22*d33    &
            &     -d32*d23*d14-d33*d24*d12-d34*d22*d13)

     det42=d11*& !det3_ARM_RP( s2,d23,d24, s3,d33,d34, s4,d43,d44) &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
           &( s2*d33*d44+d23*d34*s4+d24*s3*d43    &
            &     -s4*d33*d24-d43*d34*s2-d44*s3*d23)&
       &  -d21*& !det3_ARM_RP( s1,d13,d14, s3,d33,d34, s4,d43,d44)  &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
          &(  s1*d33*d44+d13*d34*s4+d14*s3*d43    &
            &     -s4*d33*d14-d43*d34*s1-d44*s3*d13)&
       &  +d31*& !det3_ARM_RP( s1,d13,d14, s2,d23,d24, s4,d43,d44)  &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
          &(  s1*d23*d44+d13*d24*s4+d14*s2*d43    &
            &     -s4*d23*d14-d43*d24*s1-d44*s2*d13)&
       &  -d41*& !det3_ARM_RP( s1,d13,d14, s2,d23,d24, s3,d33,d34)  
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
          &(  s1*d23*d34+d13*d24*s3+d14*s2*d33    &
            &     -s3*d23*d14-d33*d24*s1-d34*s2*d13)
 
     det43=d11*& !det3_ARM_RP(d22, s2,d24,d32, s3,d34,d42, s4,d44) &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
           &( d22*s3*d44+s2*d34*d42+d24*d32*d34    &
            &     -d42*s3*d24-s4*d34*d22-d44*d32*s2) &
       &  -d21*& !det3_ARM_RP(d12, s1,d14,d32, s3,d34,d42, s4,d44)  &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
          &(  d12*s3*d44+s1*d34*d42+d14*d32*s4    &
            &     -d42*s3*d14-s4*d34*d12-d44*d32*s1)&
       &  +d31*& !det3_ARM_RP(d12, s1,d14,d22, s2,d24,d42, s4,d44)  &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
          &(  d12*s2*d44+s1*d24*d42+d14*d22*s4    &
            &     -d42*s2*d14-s4*d24*d12-d44*d22*s1)&
       &  -d41*& !det3_ARM_RP(d12, s1,d14,d22, s2,d24,d32, s3,d34)
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
         &(  d12*s2*d34+s1*d24*d32+d14*d22*s3    &
            &     -d32*s2*d14-s3*d24*d12-d34*d22*s1)

     det44=d11*& !det3_ARM_RP(d22,d23, s2,d32,d33, s3,d42,d43, s4) &
     ! function   det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
         & (  d22*d33*s4+d23*s3*d42+s2*d32*d43    &
            &     -d42*d33*s2-d43*s3*d22-s4*d32*d23) &
       &  -d21*& !det3_ARM_RP(d12,d13, s1,d32,d33, s3,d42,d43, s4)  &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
         & (   d12*d33*s4+d13*s3*d42+s1*d32*d43    &
            &     -d42*d33*s1-d43*s3*d12-s4*d32*d13) & 
       &  +d31*& !det3_ARM_RP(d12,d13, s1,d22,d23, s2,d42,d43, s4)  &
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
           &( d12*d23*s4+d13*s2*d42+s1*d22*d43    &
            &     -d42*d23*s1-d43*s2*d12-s4*d22*d13 )&
       &  -d41*& !det3_ARM_RP(d12,d13, s1,d22,d23, s2,d32,d33, s3)
     ! function det3_ARM_RP(r11,r12,r13,r21,r22,r23,r31,r32,r33) 
           &( d12*d23*s3+d13*s2*d32+s1*d22*d33    &
            &     -d32*d23*s1-d33*s2*d12-s3*d22*d13)
       
      aa(J,4)=det44*deti
      aa(J,3)=det43*deti
      aa(J,2)=det42*deti
      aa(J,1)=det41*deti
        !  QU_dp(1,1)=deti
       !   write(*,*) 'deti',QU_dp(1,1)
       !   read(*,*)
       !   QU_dp(1,1)=det44
       !   write(*,*) 'det44',QU_dp(1,1)
       !   read(*,*)
       !   QU_dp(1,1)=det43
       !   write(*,*) 'det43',QU_dp(1,1)
       !   read(*,*)
       !   QU_dp(1,1)=det42
       !   write(*,*) 'det42',QU_dp(1,1)
       !   read(*,*)
       !   QU_dp(1,1)=det41
       !   write(*,*) 'det41',QU_dp(1,1)
       !   read(*,*)
  end do 
  !write(*,*) 'aa's finished'
!  DO J=1,M
!   DO I=2,N-1
!      QU(I,J)=ws(I,J,0) +aa(J,1)*ws(I,J,1) &
!                    &   +aa(J,2)*ws(I,J,2)  &
!                    &   +aa(J,3)*ws(I,J,3)  &
!                    &   +aa(J,4)*ws(I,J,4)
!  enddo
!   write(*,*) J, (QU(I,J), I=1,N)
!read(*,*)
!  enddo

  !! rewritten to change precision at poles as desired
   
! call cpu_time(end_time)
! write(*,*) 'time inversion', end_time-start_time


! call cpu_time(start_time)
        DO I=1,N-1
      DO J=1,M
!          QU_dp(I,J)=QU(I,J)
          U(J,I)=ws(J,I,0) +aa(J,1)*ws(J,I,1)
!          QU_dp(I,J)=QU(I,J)
        enddo
      enddo
        DO I=1,N-1
      DO J=1,M
!          QU_dp(I,J)=QU(I,J)
          QU(J,I)=U(J,I)  &
                        &   +aa(J,2)*ws(J,I,2)
!          QU_dp(I,J)=QU(I,J)
        enddo
      enddo

        DO I=1,N-1
      DO J=1,M
!          QU_dp(I,J)=QU(I,J)
          U(J,I)=QU(J,I)  &
                        &   +aa(J,3)*ws(J,I,3)  
!          QU_dp(I,J)=QU(I,J)
        enddo
      enddo
        DO I=1,N-1
      DO J=1,M
!          QU_dp(I,J)=QU(I,J)
          QU(J,I)=U(J,I)  &
                        &   +aa(J,4)*ws(J,I,4)
!          QU_dp(I,J)=QU(I,J)
        enddo
      enddo

 !     QU_dp(:,:)=QU(:,:)
 !   write(*,*) 'QU', minval(abs(QU_dp(:,:))), maxval(abs(Qu_dp(:,:)))
 !   read(*,*)
  CALL XBC_T_rp(QU,M,N)
 !     QU_dp(:,:)=QU(:,:)
 !   write(*,*) 'QU XBC', minval(abs(QU_dp(:,:))), maxval(abs(Qu_dp(:,:)))
 !   read(*,*)
     
! call cpu_time(end_time)
! write(*,*) 'time construct QU', end_time-start_time
 
! call cpu_time(start_time)
  call adjust_conservation_Trans_ARM_RP(QU,s,S_full,n,m)
! call cpu_time(end_time)
! write(*,*) 'time mass fixer', end_time-start_time
  !    QU_dp(:,:)=QU(:,:)
  !    write(*,*) 'S_full',S_full
  !  write(*,*) 'QU adjust', minval(abs(QU_dp(:,:))), maxval(abs(Qu_dp(:,:)))
  !  read(*,*)
   !! end of rewrite
   
  !write(*,*) 'QU finished'






end do
  !write(*,*) 'BC QU finished'
END SUBROUTINE
 
SUBROUTINE XBC_T_B_rp(X,N,Q,M_Q,M_R)
implicit none

integer, parameter :: HP=2
INTEGER :: N,Q, M_Q, M_R
REAL(HP) :: X(M_Q,N,Q)

INTEGER :: J,K

do K=1,Q-1
DO J=1,M_Q
  X(J,1,K)=X(J,N-1,K)
  X(J,N,K)=X(J,2,K)
end do
enddo
K=Q
 DO J=1,M_R
  X(J,1,K)=X(J,N-1,K)
  X(J,N,K)=X(J,2,K)
 end do

END subroutine

SUBROUTINE XBC_T_sp(X,M,N)
implicit none

INTEGER :: N, M
REAL(Kind=4) :: X(M,N)

INTEGER :: J

DO J=1,M
  X(J,1)=X(J,N-1)
  X(J,N)=X(J,2)
end do
END subroutine
SUBROUTINE XBC_T_rp(X,M,N)
implicit none

integer, parameter :: HP=2
INTEGER :: N, M
REAL(HP) :: X(M,N)

INTEGER :: J

DO J=1,M
  X(J,1)=X(J,N-1)
  X(J,N)=X(J,2)
end do
END subroutine

subroutine adjust_conservation_Trans_Block_ARM_RP(p,s,S_full,n,Q,M_Q,M_R)

implicit none
integer, parameter :: HP=2
REAL(Kind=HP) :: p(M_Q,N,Q),s(M_Q,N,Q), cnst, cnst_sp(M_Q,Q)
REAL(Kind=4) ::  S_full, start_time, end_time ! cnst_sp(m)
double precision :: help1, help2
INTEGER :: N,Q, M_Q,M_R

INTEGER :: I, J,K
    !write(*,*) S_full
 !   call cpu_time(start_time)
      cnst_sp(:,:)=0.0_HP
   Do K=1,Q-1,1
     do i=2,n-1
       do j=1,M_Q
          cnst_sp(J,K)=cnst_sp(J,K)+s(J,I,K)*p(J,I,K)
        enddo
      enddo
   enddo
   K=Q
     do i=2,n-1
       do j=1,M_R
          cnst_sp(J,K)=cnst_sp(J,K)+s(J,I,K)*p(J,I,K)
        enddo
      enddo
   K=1
       do j=2,M_Q
          cnst_sp(J,1)=cnst_sp(J,1)+cnst_sp(J,K)
      enddo
   Do K=2,Q-1,1
       do j=1,M_Q
          cnst_sp(J,1)=cnst_sp(J,1)+cnst_sp(J,K)
      enddo
   enddo
   K=Q
       do j=1,M_R
          cnst_sp(J,1)=cnst_sp(J,1)+cnst_sp(J,K)
      enddo
   K=1
       do j=2,M_Q
          cnst_sp(1,1)=cnst_sp(1,1)+cnst_sp(J,1)
      enddo
     ! write(*,*) 'check S_full,', help1
      cnst=cnst_sp(1,1)/S_full
  !   call cpu_time(end_time)
  !   write(*,*) 'conservation step1', end_time-start_time
     !write(*,*) cnst

  !  call cpu_time(start_time)
   Do K=1,Q-1,1
     do i=1,n
       do j=1,M_Q
          p(J,I,K)=p(J,I,K)-cnst
        enddo
      enddo
   enddo
   K=Q
     do i=1,n
       do j=1,M_R
          p(J,I,K)=p(J,I,K)-cnst
        enddo
      enddo
  !   call cpu_time(end_time)
  !   write(*,*) 'conservation step2', end_time-start_time

end subroutine

subroutine adjust_conservation_Trans_ARM_RP(p,s,S_full,n,m)

implicit none
integer, parameter :: HP=2
REAL(Kind=HP) :: p(m,n),s(m,n), cnst, cnst_sp(m)
REAL(Kind=4) ::  S_full, start_time, end_time ! cnst_sp(m)
double precision :: help1, help2
INTEGER :: N, M

INTEGER :: I, J
    !write(*,*) S_full
   ! call cpu_time(start_time)
      cnst_sp(:)=0.0_HP
        do i=2,n-1
      do j=1,m
          !write(*,*) i, j, s(i,j), p(i,j), s(i,j)*p(i,j)
          !read(*,*)
          cnst_sp(J)=cnst_sp(J)+s(J,I)*p(J,I)
          !help1=help1+s(i,j)
          !help2=p(i,j)
          
          !write(*,*) help1, help2,help1*help2, cnst
          !read(*,*)
        enddo
      enddo
      !do j=2,m
      !    cnst_sp(1)=cnst_sp(1)+cnst_sp(J)
      !  enddo
      do j=1,m/32-1
          cnst_sp(1:32)=cnst_sp(1:32)+cnst_sp(J*32+1:J*32+32)
        enddo
      do j=2,32
          cnst_sp(1)=cnst_sp(1)+cnst_sp(J)
        enddo
     ! write(*,*) 'check S_full,', help1
      cnst=cnst_sp(1)/S_full
    ! call cpu_time(end_time)
    ! write(*,*) 'conservation step1', end_time-start_time
     !write(*,*) cnst

    ! call cpu_time(start_time)
        do i=1,n
      do j=1,m
          p(J,I)=p(J,I)-cnst
        enddo
      enddo
     ! call cpu_time(end_time)
     ! write(*,*) 'conservation step2', end_time-start_time

end subroutine

SUBROUTINE LAPL_Trans_ARM_RP(N,M,IP)

implicit none
integer, parameter :: HP=2
REAL(kind=HP):: F(M,N)
REAL(kind=4):: F_sp(M,N),P_sp(M,N)
REAL(Kind=HP) :: P(M,N),A11(M,N),A12(M,N),A21(M,N),A22(M,N),  &
    &      B11(M,N),B22(M,N),P0(M,N),U(M,N),V(M,N),S(M,N),S_L(M,N)
INTEGER :: IP(N)
INTEGER :: N, M, num_of_bits, DP_Depth

REAL(Kind=HP) :: UTIL, VTIL, scaling
INTEGER :: I, J

REAL(Kind=4) :: start_time, end_time


!DO J=1+DP_Depth,M-DP_Depth
! DO I=1,N
!  S_L(I, J)=S(I, J)
!end do
!end do


   ! call cpu_time(start_time)
  DO I=2,N-1
DO J=2,M-1
    U(J,I)=P(J,I+1)-P(J,I-1)
    V(J,I)=P(J+1,I)-P(J-1,I)
  end do
end do
  
DO I=2,N-1
  U(1,I)=P(1,I+1  )-P(1  ,I-1  )
  U(M,I)=P(M,I+1  )-P(M  ,I-1  )
  V(1,I)=P(2,I    )-P(1  ,IP(I))
  V(M,I)=P(M,IP(I))-P(M-1,I    )
ENDDO  
CALL XBC_T_rp(U,M,N)
CALL XBC_T_rp(V,M,N)

   !  call cpu_time(end_time)
   !  write(*,*) 'Lapl step1', end_time-start_time


    ! call cpu_time(start_time)

!DO J=1,M
!  DO I=1,N
!    UTIL=U(I,J)*A11(I,J)+A12(I,J)*V(I,J)+B11(I,J)*(P(I,J))  ! why only B11*P and not B11*(P-P0)
!    VTIL=V(I,J)*A21(I,J)+A22(I,J)*U(I,J)+B22(I,J)*(P(I,J))
!    U(I,J)=UTIL
!    V(I,J)=VTIL
!  ENDDO
!ENDDO

  !! rewritten to change precision at poles as desired    
        DO I=1,N
      DO J=1,M
    UTIL=U(J,I)*A11(J,I)+A12(J,I)*V(J,I)+B11(J,I)*(P(J,I))  ! why only B11*P and not B11*(P-P0)
    VTIL=V(J,I)*A21(J,I)+A22(J,I)*U(J,I)+B22(J,I)*(P(J,I))
    U(J,I)=UTIL
    V(J,I)=VTIL
        enddo
      enddo

CALL XBC_T_rp(U,M,N)
CALL XBC_T_rp(V,M,N)

   !! end of rewrite

    ! call cpu_time(end_time)
    ! write(*,*) 'Lapl step2', end_time-start_time

   !  call cpu_time(start_time)
  DO I=2,N-1
DO J=2,M-1
    F(J,I)= 0.5_HP*(U(J,I+1)-U(J,I-1)+V(J+1,I)-V(J-1,I))/S(J,I)
  end do
end do    
 
DO I=2,N-1
  F(1,I)= 0.5_HP*(U(1,I+1)-U(1,I-1)+(V(2,I)+V(1,I)))/S(1,I)
  F(M,I)= 0.5_HP*(U(M,I+1)-U(M,I-1)-(V(M,I)+V(M-1,I)))/S(M,I)
ENDDO
CALL XBC_T_rp(F,M,N)

        DO I=1,N
      DO J=1,M
      p_sp(J,I)=p(J,I)*scaling 
      enddo
      enddo
        DO I=1,N
      DO J=1,M
      F_sp(J,I)=F(J,I)*scaling 
      enddo
      enddo
        DO I=1,N
      DO J=1,M
      F_sp(J,I)=F_sp(J,I)  -p_sp(J,I) 
      enddo
      enddo
     ! call cpu_time(end_time)
     ! write(*,*) 'Lapl step3', end_time-start_time
END SUBROUTINE

SUBROUTINE LAPL_Trans_Block_ARM_RP(N, M,Q,M_Q,M_R,JtoJM1,JtoJP1,IP,M2, M3, MM1, MM2)

implicit none
integer, parameter :: HP=2
REAL(kind=HP):: F(M_Q,N,Q)
REAL(Kind=HP) :: P(M_Q,N,Q),A11(M_Q,N,Q),A12(M_Q,N,Q),A21(M_Q,N,Q),A22(M_Q,N,Q),  &
    &      B11(M_Q,N,Q),B22(M_Q,N,Q),P0(M_Q,N,Q),U(M_Q,N,Q),V(M_Q,N,Q),S(M_Q,N,Q),S_L(M_Q,N,Q), &
    &      U_p(M_Q,N,Q),V_p(M_Q,N,Q)
INTEGER :: IP(N)
INTEGER :: N, M,Q,M_Q,M_R, num_of_bits, DP_Depth

REAL(Kind=HP) :: UTIL, VTIL
INTEGER :: I, J, K
integer :: M2(2), M3(2), MM1(2), MM2(2)
integer ::JtoJM1(M_Q,Q,2), JtoJP1(M_Q,Q,2)

REAL(Kind=4) :: start_time, end_time
!DO J=1+DP_Depth,M-DP_Depth
! DO I=1,N
!  S_L(I, J)=S(I, J)
!end do
!end do

!     call cpu_time(start_time)
    do K=1,Q-1,1   
      DO I=2,N-1
        DO J=1,M_Q
                                   
        U(J,I,K)=P(J,I+1,K)-P(J,I-1,K)
        enddo
      enddo
   enddo
   K=Q
      DO I=2,N-1
        DO J=1,M_R
        U(J,I,K)=P(J,I+1,K)-P(J,I-1,K)
        enddo
      enddo
      
   K=M2(1)
      DO I=2,N-1
        DO J=M2(2),M_Q
                                   
        V(J,I,K)=P(JtoJP1(J,K,2),I,JtoJP1(J,K,1))-P(JtoJM1(J,K,2),I,JtoJM1(J,K,1))
        enddo
      enddo
    do K=M2(1)+1,MM1(1)-1,1   
      DO I=2,N-1
        DO J=1,M_Q
                                   
        V(J,I,K)=P(JtoJP1(J,K,2),I,JtoJP1(J,K,1))-P(JtoJM1(J,K,2),I,JtoJM1(J,K,1))
        enddo
      enddo
    enddo
   K=MM1(1)
      DO I=2,N-1
        DO J=1,MM1(2)
        V(J,I,K)=P(JtoJP1(J,K,2),I,JtoJP1(J,K,1))-P(JtoJM1(J,K,2),I,JtoJM1(J,K,1))
        enddo
      enddo
      
      do i=2,n-1

       U(1,I,1)=P(1,I+1,1)-P(1,I-1,1)  ! J=1 is J=1,K=1
       U(M_R,I,Q)=P(M_R,I+1,Q)-P(M_R,I-1,Q)      ! J=M is J=M_R, K=Q
       V(1,I,1)=P(M2(2),i,M2(1))-P(1,ip(I),1)  ! J=2 is J=M2(2), K=M2(1)
       V(M_R,I,Q)=P(M_R,ip(I),Q)-P(MM1(2),I,MM1(1))    ! J=M is J=M_R, K=Q; J=M-1 is J=MM1(2), J=MM1(1)
      enddo

CALL XBC_T_B_rp(U,N,Q,M_Q, M_R)
CALL XBC_T_B_rp(V,N,Q,M_Q, M_R)


!     call cpu_time(end_time)
!     write(*,*) 'Lapl step1', end_time-start_time
!
!     call cpu_time(start_time)

!DO J=1,M
!  DO I=1,N
!    UTIL=U(I,J)*A11(I,J)+A12(I,J)*V(I,J)+B11(I,J)*(P(I,J))  ! why only B11*P and not B11*(P-P0)
!    VTIL=V(I,J)*A21(I,J)+A22(I,J)*U(I,J)+B22(I,J)*(P(I,J))
!    U(I,J)=UTIL
!    V(I,J)=VTIL
!  ENDDO
!ENDDO

  !! rewritten to change precision at poles as desired    
  do k=1,Q-1,1   
    DO I=1,N
      DO J=1,M_Q
    U_p(J,I,K)=U(J,I,K)*A11(J,I,K)+A12(J,I,K)*V(J,I,K)+B11(J,I,K)*(P(J,I,K))  ! why only B11*P and not B11*(P-P0)
    V_p(J,I,K)=V(J,I,K)*A21(J,I,K)+A22(J,I,K)*U(J,I,K)+B22(J,I,K)*(P(J,I,K))
      enddo
    enddo
  enddo
  K=Q   
    DO I=1,N
      DO J=1,M_R
    U_p(J,I,K)=U(J,I,K)*A11(J,I,K)+A12(J,I,K)*V(J,I,K)+B11(J,I,K)*(P(J,I,K))  ! why only B11*P and not B11*(P-P0)
    V_p(J,I,K)=V(J,I,K)*A21(J,I,K)+A22(J,I,K)*U(J,I,K)+B22(J,I,K)*(P(J,I,K))
      enddo
    enddo

CALL XBC_T_B_rp(U_p,N,Q,M_Q, M_R)
CALL XBC_T_B_rp(V_p,N,Q,M_Q, M_R)

!     call cpu_time(end_time)
!     write(*,*) 'Lapl step2', end_time-start_time
!     call cpu_time(start_time)
   !! end of rewrite
   K=M2(1)
      DO I=2,N-1
        DO J=M2(2),M_Q
                                   
      F(J,I,K)= (U_p(J,I+1,K)-U_p(J,I-1,K)+V_p(JtoJP1(J,K,2),I,JtoJP1(J,K,1))-V_p(JtoJM1(J,K,2),I,JtoJM1(J,K,1)))/(S(J,I,K))
        enddo
      enddo
    do K=M2(1)+1,MM1(1)-1,1   
      DO I=2,N-1
        DO J=1,M_Q
                                   
      F(J,I,K)= (U_p(J,I+1,K)-U_p(J,I-1,K)+V_p(JtoJP1(J,K,2),I,JtoJP1(J,K,1))-V_p(JtoJM1(J,K,2),I,JtoJM1(J,K,1)))/(S(J,I,K))
        enddo
      enddo
    enddo
   K=MM1(1)
      DO I=2,N-1
        DO J=1,MM1(2)
      F(J,I,K)= (U_p(J,I+1,K)-U_p(J,I-1,K)+V_p(JtoJP1(J,K,2),I,JtoJP1(J,K,1))-V_p(JtoJM1(J,K,2),I,JtoJM1(J,K,1)))/(S(J,I,K))
        enddo
      enddo
!    DO I=2,N-1
!  DO J=2,M-1
!      F(J,I)= (U(J,I)-V_p(J-1,I))/(S(J,I))
!    end do
!  end do    
 

  DO I=2,N-1
    F(1,I,1)= (U_p(1,I+1,1)-U_p(1,I-1,1)+(V_p(M2(2),I,M2(1))+V_p(1,I,1)))/(S(1,I,1))  ! J=1 is j=1,K=1; J=2 is J=M2(2),K=M2(1)
    F(M_R,I,Q)= (U_p(M_R,I+1,Q)-U_p(M_R,I-1,Q)-(V_p(M_R,I,Q)+V_p(MM1(2),I,MM1(1))))/(S(M_R,I,Q)) ! J=M is j=M_R,K=Q; J=M1 is J=MM1(2),K=MM1(1)
  ENDDO

CALL XBC_T_B_rp(F,N,Q,M_Q, M_R)

!     call cpu_time(end_time)
!     write(*,*) 'Lapl step3', end_time-start_time
END SUBROUTINE

end program
