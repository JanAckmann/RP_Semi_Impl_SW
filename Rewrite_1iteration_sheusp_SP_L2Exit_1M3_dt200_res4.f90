PROGRAM IS1SPSL
use implicit_functions_SP

implicit none

INTEGER, PARAMETER :: NLON=(512)+2,NLAT=256
INTEGER, PARAMETER :: N=NLON,M=NLAT,NM=N*M
 character(len=150) :: EXP_NAME, Dp_depth_str
REAL(Kind=4) ::  U(N,M,0:1),&
               & V(N,M,0:1),  &
               & PD(N,M),     &
               & QX(N,M),     &
               & QY(N,M),     &
               & PT(N,M),     & 
               & P0(N,M),     &
               & HX(N,M),     &
               & HY(N,M),     &
               & DHX2Y(N,M),  &
               & S(N,M),      &
               & F1(N,M,0:1), &
               & F2(N,M,0:1), & 
               & E1(N,M,-1:0),&
               & E2(N,M,-1:0),&
               & UA(N,M),     &
               & VA(N,M+1),   &
               & PC(N,M),     &
               & X(N),        &
               & Y(M+1),        &
               & COR(N,M),    &
               & DIFF_U(N,M), &
               & DIFF_V(N,M), &
               & Velfx(N,M),  &
               & Velfy(N,M),  &
               & QXS(N,M),    &
               & QYS(N,M)

REAL(Kind=4) :: MGH1IHX(M),  &
               & MGH2IHY(M),  &
               & AC(M),       &
               & BC(M),       &
               & AD(M),       &
               & BD(M), A, B, C, D
REAL(Kind=4) ::U_23(N,M),&
                 & V_23(N,M),     &
                 & HX_23(N,M),     &
                 & HY_23(N,M),     &
                 & DHX2Y_23(N,M),  &
                 & S_23(N,M),      &
                 & X_23(N),        &
                 & Y_23(M+1),        &
                 & COR_23(N,M)

REAL(Kind=4) ::TIME,&
               & DX, &
               & DY, &
               & DT, &
               & mue

REAL(Kind=4) :: DX_23, &
               & DY_23, &
               & DT_23

INTEGER :: IP(N)
REAL(Kind=4) :: sum_time, sum_lp_time

!!! store old values in low prec. to calc final tend.

REAL(Kind=4) ::    PD_T(N,M), &
                  & QX_T(N,M), &
                  & QY_T(N,M)


                  
!!! the HIGH-PRECISION VARIABLES
REAL(Kind=4) :: PD_HP(N,M), &
                  & QX_HP(N,M), &
                  & QY_HP(N,M), &
                  & UA_HP(N,M), &
                  & VA_HP(N,M), &
                  & PT_HP(N,M), &
                  & P0_HP(N,M), &
                  & PD_old(N,M), &
                  & QX_old(N,M), &
                  & QY_old(N,M), &
                  & div_old(N,M), &
                  & div_new(N,M), &
                  & r_eval(N,M)
                  
                  
! CORIOLIS, GRAVITY AND EARTH RADIUS SPECIFICATION
INTEGER, PARAMETER  :: ICORIO=1
REAL(Kind=4) :: F0
REAL(Kind=4) :: G 
REAL(Kind=4) :: R

! CHARACTERISTICS OF THE FLOW
REAL(Kind=4) :: USCAL 
REAL(Kind=4) :: H00   
REAL(Kind=4) :: HMTS  

! Polar Filters
logical ::  QRelax
REAL(Kind=4) :: U0(N,M),     &
                  & V0(N,M),     &
                  & PT0(N,M),     &
                  & PD0(N,M)
REAL(Kind=4) :: Alp_REL(N,M), &
                  & atau



REAL(Kind=4) ::  F0_23
REAL(Kind=4) :: G_23
REAL(Kind=4) ::  R_23

! CHARACTERISTICS OF THE FLOW
REAL(Kind=4) :: USCAL_23 
REAL(Kind=4) ::  H00_23   
REAL(Kind=4) ::  HMTS_23, GMM_23, AMM, DETI

Integer :: IPRINT


! CONTROL TESTS: ZONAL FLOW OR ROSSBY WAVE
INTEGER:: IRHW, &
        & DP_Depth


! ! CONTROL EVALUATION OF THE PRESSURE GRADIENTS: IPS=0 CENTERED DIFFERNCING
! CHOOSING IPS=1 GIVES WEIGHTED AVERAGE OF ONE-SIDED DERIVATIVES THAT
! CONVERGES TO ONE SIDED DERIVATIVE AT THE CROSSECTION WITH THE BOTTOM
! !JA IPS=1 is deactivated
INTEGER, PARAMETER  :: IPS=0

! CREATE TAPEWR
INTEGER, PARAMETER :: IANAL=0,IRST=0,IWRITE=0
INTEGER, PARAMETER :: NFIL=50
REAL(Kind=4) :: DTFIL(NFIL),&
     & NTFIL(NFIL)

INTEGER:: NITER,  &
        & NITSM,  &
        & ICOUNT
REAL(Kind=4) :: ERROR


REAL(Kind=4) ::PI, &
     & PI2, &
     & PIH, &
     & PVEL, &
     & BETA, &
     & GI, &
     & EP
     
REAL(Kind=4) :: PI_23, &
     & PI2_23, &
     & PIH_23, &
     & PVEL_23, &
     & BETA_23, &
     & GI_23

INTEGER :: IORD, &
         & ISOR, &
         & NONOS, &
         & IDIV, &
         & ISGNSPL

INTEGER :: KMX, mpfl, liner

! grid creation
REAL(Kind=4) :: GC1, &
     & GC2, &
     & GH1, &
     & GH2

REAL(Kind=4) :: GC1_23, &
     & GC2_23, &
     & GH1_23, &
     & GH2_23
     
REAL(Kind=4) :: D0, &
     & S_full         , &
     & Exit_Cond, adv_cour 


REAL(Kind=4) :: D_Adv(N,M)
 
! the rest
REAL(Kind=4) :: GMM, SUM1, SUM0, start, finish
INTEGER :: I, J, KF, KT, NPLOT, NPRINT, NT, num_of_bits, ID_PREC
INTEGER :: stencil, ueber

logical :: codesignQ, codesignD, mountain, gcr23, save_time
mountain = .false.
 !!! RPE VARIABLES


call cpu_time(start)

  codesignQ = .FALSE.
  codesignD = .FALSE.
  gcr23     = .false.


  QRelax    =.true.
  !DP_Depth = 6
  
 !!! Experiment identifiers 
stencil=0


do ID_PREC=7,7,-5
 do IRHW = 3,3,2

  do DP_Depth=0,0,2
   write(Dp_depth_str,*) DP_Depth

  !ID_PREC=0
   EXP_NAME= 'data/sheusp_SP_1iteration_L2Exit_1M3_dt200_res4'
  ! EXP_NAME= 'data_ADI_Precon_init23'



   do num_of_bits=23,24,2
    save_time=.false.
   
    If (num_of_bits>2) then
! if .not. init23

   

  write(*,*) 'Experiment' , IRHW
  write(*,*) 'Resolution' , NLON, NLAT
  write(*,*) 'Preconditioner', ID_PREC
  write(*,*) 'Dp Zonal bands', DP_Depth
  write(*,*) 'Experiment folder ',   EXP_NAME
  write(*,*) 'Relaxation ',  QRelax


  sum_time=0.0d0
  sum_lp_time=0.0d0

 call rpenum_init(num_of_bits)


   write(*,*) 'codesign: Q | D ', codesignQ, codesignD, 'default bits: ', num_of_bits



!If (IRHW==0) then
! DATA NT,NPRINT/12096,864/
!  DT_23=100.0d0 
! ! what if it ran with the timestep of the explicit model 
! !DT_23=40.0d0   
! DT = DT_23
!else
!   DATA NT,NPRINT/2016,144/
! DT_23=600.0d0 
!if (IRHW==2) then
!NT = 4032 
!NPRINT =288
!DT_23=300.0d0 
!else
!NT = 3456 !240 
!NPRINT =216
!DT_23=400.0d0
!endif
if (IRHW==2) then
!DATA NT,NPRINT/12096,864/
NT = 2*int(2*6376) !12960  
NPRINT =int(2*797)!864
DT_23=100.0d0
KMX=4
mpfl=999999
elseif(IRHW==1) then
!DATA NT,NPRINT/12096,864/
NT = 6480 ! 6376 !int(6376*(200.0/240.0)) !12960  
NPRINT = 216 !797 !797 !int(797*(200.0/240.0)) !864
!NT = 6376 !int(6376*(200.0/240.0)) !12960  
!NPRINT =200 !797 !797 !int(797*(200.0/240.0)) !864
DT_23=200.0d0
KMX=4
atau=200.*DT_23 ! RHW4
mpfl=999999
elseif(IRHW==3) then
NT =  6480! 6376 !int(6376*(200.0/240.0)) !12960  
NPRINT =216 !797 !797 !int(797*(200.0/240.0)) !864 !DATA NT,NPRINT/12096,864/
!NT = 6376 !int(6376*(200.0/240.0)) !12960  
!NPRINT = 200 !797 !797 !797 !int(797*(200.0/240.0)) !864
DT_23=200.0d0
KMX=4
atau=2.*DT_23    !Zonal flow past Earth orography
mpfl=999999
endif
 ! what if it ran with the timestep of the explicit model 
 ! DT_23=40.0d0   
 DT = DT_23
!endif
 ! what if it ran with the timestep of the explicit model 
 ! DATA NT,NPRINT/30240,2160/
 !!! start of calculating every coefficient in HP
 
 
F0_23=1.4584E-4
F0=F0_23

 
G_23 =9.80616
G=G_23
GI_23=1.0d0/G_23
GI = GI_23

! end param test1
R_23=6371.22E+03
R=R_23


! CHARACTERISTICS OF THE FLOW
USCAL_23 = 20.!Piotrs new numbers old !5.0d0
H00_23  = 8.E3
HMTS_23  = 0.E3 

 !Piotrs new numbers old 6.E3  !! changed from HMTS  = 2.E3 in accordance to explicit
if (IRHW==2) then
!  
  mountain=.true.
  USCAL_23 = 20.0d0 
  H00_23   = 5960.0d0 
  HMTS_23  = 1.0d0 
!
elseif(IRHW==3) then
USCAL_23 = 20.!Piotrs new numbers old !5.0d0
H00_23  = 8.E3
  HMTS_23  = 1.0 
  endif

USCAL=USCAL_23
H00 = H00_23
HMTS= HMTS_23


DTFIL(:) = NFIL*40.0d0
NTFIL(:) = NFIL*2160


NITSM  = 0
ICOUNT = 0





!COMPUTATIONS:

!      CALL OPNGKS
!      CALL GSCLIP(0)
!     CALL GOPKS (6,0)        ! open GKS
!     CALL GOPWK (1,2,1)      ! open CGM workstation
!     CALL GACWK (1)          ! activate workstation
!     CALL GOPWK (2,0,8)      ! open X workstation
!     CALL GACWK (2)          ! activate workstation

!COMPUTE SOME RELEVANT CONSTANTS

PI_23=ACOS(-1.0d0)
PI = PI_23
PI2_23=2.0d0*PI_23
PI2 = PI2_23
PIH_23=0.5d0*PI_23
PIH = PIH_23
TIME=0.0d0
PVEL_23=USCAL_23/R_23
PVEL = PVEL_23
F0_23=ICORIO*F0_23
F0 = F0_23

! BETA IS AN ANGLE IN THE ZONAL FLOW TEST

BETA_23=0.0d0 
if (IRHW==2) then
BETA_23=0.0d0 
endif
BETA = BETA_23
H00_23=G_23*H00_23
H00 = H00_23


! PARAMETERS FOR MPDATA ADVECTION
IORD=2
ISOR=1
NONOS=1
IDIV=1               
ISGNSPL=0

! SMALL CONSTANT (SMALL COMPARED TO DEPTH OF THE ATMOSPHERE)
EP= 1.E-6   !! ersetzen durch tiny times factor oder was aehnliches

!COMPUTE GRID      


DX_23= PI2/FLOAT(N-2)
DX = DX_23
DY_23= PI/FLOAT(M)
DY = DY_23
! param test2


GC1_23= DT/DX                                ! DT/DX
GC1 = GC1_23
GC2_23= DT/DY                                        ! DT/DY
GC2 = GC2_23

GH1_23= .5*GC1                         ! 0.5d0*GC1
GH1 = GH1_23
GH2_23= .5*GC2                                ! 0.5d0*GC2
GH2 = GH2_23


! end param test2

DO J=1,M+1
  Y_23(J)=-PIH+(float(J)-0.5)*DY

  Y(J) = Y_23(J)

end do


DO I=2,N-1
  X_23(I)=(float(I)-1)*DX_23
  X(I) = X_23(I)
end do
X_23(1)=X_23(N-1)
X_23(N)=X_23(2)

X(1)=X(N-1)
X(N)=X(2)

DO I=1,N
  IP(I)=MOD(I+(N-2)/2-1,N-2)+1
end do
!COMPUTE METRIC TERMS FOR SPHERICAL COORDINATES

DO J=1,M
  DO I=1,N
    ! param test3


    HX_23(I,J)=R_23*COS(Y_23(J))

    HX(I,J) = HX_23(I,J)
    HY_23(I,J)=R_23
    HY(I,J) = HY_23(I,J)


    ! end param test3
    ! TEST DHX2Y(I,J)=-DT*R**2*SIN(2.*Y(J))/(HX(I,J)**2*HY(I,J))*.5
    !param test4


    S_23(I,J)=HX_23(I,J)*HY_23(I,J)
    S(I,J) = S_23(I,J)


    ! end param test4
  end do

end do

!param test5


DO I=2,N-1
  DO J=2,M-1
    DHX2Y_23(I,J)= (HX_23(I,J+1)-HX_23(I,J-1))*GC2_23/S_23(I,J)*.5
    DHX2Y(I,J) = DHX2Y_23(I,J) 
  end do

  DHX2Y_23(I,1)= (HX_23(I,  2)+HX_23(I,  1))*GC2_23/S_23(I,1)*.5
  DHX2Y(I,1) = DHX2Y_23(I,1)
  DHX2Y_23(I,M)=-(HX_23(I,  M)+HX_23(I,M-1))*GC2_23/S_23(I,M)*.5
  DHX2Y(I,M) = DHX2Y_23(I,M) 
end do


!end param test5
      
CALL XBC(DHX2Y,N,M)
CALL XBC_23(DHX2Y_23,N,M)




!CONDITIONS OF THE INITIAL STATE ***********************************

If (IRHW.EQ.2)then
  CALL TOPOGR(P0_HP,X_23,Y_23,N,M, mountain)
elseif (IRHW.EQ.3)then
  CALL EARTHTOPO(P0_HP,X_23,Y_23,N,M)
  CALL SMOOTHTOP(P0_HP,PC,IP,N,M)
else
  P0_HP(:,:)=0.0
endif

IF(IRHW.EQ.0) CALL INITZON(U_23,V_23,PT_HP,COR_23,X_23,Y_23,N,M,F0_23,BETA_23,H00_23,R_23,PVEL_23)
IF(IRHW.EQ.1) CALL INITRHW(U_23,V_23,PT_HP,COR_23,X_23,Y_23,N,M,F0_23,R_23)
IF(IRHW.EQ.2) CALL INITZON(U_23,V_23,PT_HP,COR_23,X_23,Y_23,N,M,F0_23,BETA_23,H00_23,R_23,PVEL_23)
IF(IRHW.EQ.3) CALL INITZON(U_23,V_23,PT_HP,COR_23,X_23,Y_23,N,M,F0_23,BETA_23,H00_23,R_23,PVEL_23)

If (QRelax) then
CALL POLARABS(Alp_REL,atau,Y_23,DT_23,N,M,IRHW)
else
Alp_REL(:,:)=0.0d0
endif

!! call calculateAVGdepth(PD)

IF(IRST.EQ.0) THEN
! INITIATE PRIMARY VARIABLES (PD, QX, QY)
  DO J=1,M
    DO I=1,N
      COR(I,J)=COR_23(I,J)*DT_23 
      P0_HP(I,J)=  P0_HP(I,J)*HMTS_23
      PD_HP(I,J)= max(rpe_0, PT_HP(I,J)*GI_23-P0_HP(I,J))

      QX_HP(I,J)=PD_HP(I,J)*U_23(I,J)
      QY_HP(I,J)=PD_HP(I,J)*V_23(I,J)

      QXS(I,J)=QX_HP(I,J)
      QYS(I,J)=QY_HP(I,J)

    end do

  end do

!! call calculateAVGdepth(PD)
D0=0.0d0
S_full=0.0d0

  DO J=1,M
    DO I=2,N-1
      S_full=S_full+S_23(I,J)
    end do
  end do
  write(*,*) 'Area Earth', S_full
  DO J=1,M
    DO I=2,N-1
      D0=D0+PD_HP(I,J)*S_23(I,J)/S_full
    end do
  end do
  write(*,*) 'average Depth', D0  


!! GCR parameters
If(QRelax) then

      DO J=1,M
        
       AMM=1.0d0+0.5d0*Alp_REL(1,J)
       GMM=0.5d0*COR(1,J)
       DETI=1.0d0/(AMM**2+GMM**2)
       A= AMM*DETI        !rpe_1/(Relaxation+GMM*GMM/(Relaxation))
       B= GMM*DETI        !GMM*A/(Relaxation)
       C=GH1*HY(1,J)*0.5d0
       D=GH2*HX(1,J)*0.5d0

       MGH1IHX(J)=-0.5d0*GC1/HX(1,J)
       MGH2IHY(J)=-0.5d0*GC2/HY(1,J)
       AC(J)=A*C
       BC(J)=B*C
       AD(J)=A*D
       BD(J)=B*D
      enddo

else

      DO J=1,M

       GMM=0.5d0*COR(1,J)
       A= rpe_1/(rpe_1+GMM*GMM)
       B=GMM*A
       C=GH1*HY(1,J)*0.5d0
       D=GH2*HX(1,J)*0.5d0

       MGH1IHX(J)=-0.5d0*GC1/HX(1,J)
       MGH2IHY(J)=-0.5d0*GC2/HY(1,J)
       AC(J)=A*C
       BC(J)=B*C
       AD(J)=A*D
       BD(J)=B*D
      enddo

endif
!! end gcr params
!! end gcr params
! if init23


  DO J=1,M
    DO I=1,N
      P0(I,J) = P0_HP(I,J)
      PT(I,J) = PT_HP(I,J) 
      PD(I,J) = PD_HP(I,J)
      QX(I,J) = QX_HP(I,J) 
      QY(I,J) = QY_HP(I,J) 
    end do
  end do

  DO J=1,M
    DO I=1,N
      U(I,J,0) =U_23(I,J)
      V(I,J,0) =V_23(I,J)
      U(I,J,1)=U(I,J,0)
      V(I,J,1)=V(I,J,0)
    end do
  end do
  DO J=1,M
    DO I=1,N
      PD_T(I,J)=rpe_0
      QX_T(I,J)=rpe_0
      QY_T(I,J)=rpe_0
    end do
  end do



!!! make sure that everything until here is initialized and calculated in REAL(Kind=4)
!!! before actually being downcasted
   call init_perf_markers(PD_HP(:,:)+P0(:,:),QX_HP(:,:)/PD_HP(:,:),QY_HP(:,:)/PD_HP(:,:), rpe_0, &
                   & codesignQ, codesignD, IRHW, X, Y, N, M, num_of_bits, ID_prec, EXP_NAME)

   call write_fields(PD_HP(:,:)+P0(:,:),QX_HP(:,:)/PD_HP(:,:),QY_HP(:,:)/PD_HP(:,:), rpe_0, &
& codesignQ, codesignD, IRHW, X, Y, N, M, num_of_bits, ID_prec, EXP_NAME)  !PD_HP(:,:)
  !! end initiate reduced variables
  
  
! INITIATE CORRESPONDING FORCES
  CALL PRFC0(PT,F1(:,:,0),F2(:,:,0),PD,HX,HY,IP,IPS,GH1,GH2,EP,N,M)

  DO J=1,M
    DO I=1,N 

      E1(I,J,0)=-DHX2Y(I,J)*QX(I,J)*QY(I,J)/max(PD(I,J),EP)
      E2(I,J,0)= DHX2Y(I,J)*QX(I,J)*QX(I,J)/max(PD(I,J),EP)
      E1(I,J,-1)=E1(I,J,0)
      E2(I,J,-1)=E2(I,J,0)
    end do
  end do


  DO J=1,M
    DO I=1,N 

      E1(I,J,0)=-DHX2Y(I,J)*QX(I,J)*QY(I,J)/max(PD(I,J),EP)
      E2(I,J,0)= DHX2Y(I,J)*QX(I,J)*QX(I,J)/max(PD(I,J),EP)
      E1(I,J,-1)=E1(I,J,0)
      E2(I,J,-1)=E2(I,J,0)
    end do
  end do


    DO J=1,M
      DO I=1,N

        F1(I,J,0)=F1(I,J,0)+COR(I,J)*QY(I,J)+E1(I,J,0)
        F2(I,J,0)=F2(I,J,0)-COR(I,J)*QX(I,J)+E2(I,J,0)

      end do
    end do



ELSE
  DO KF=1,NFIL
    READ(10) PD,PT,QX,QY,U,V,F1,F2,E1,E2
    IF(IWRITE.EQ.1) WRITE(9) PD,PT,QX,QY,U,V,F1,F2,E1,E2
    TIME=TIME+DTFIL(KF)*NTFIL(KF)
  end do
ENDIF

CALL DIAGNOS(QX_HP(:,:)/PD_HP(:,:),QY_HP(:,:)/PD_HP(:,:),PD_HP,&
                 & PT,HX,HY,IP,S,TIME,DX,DY,DT, SUM0,SUM1, &
                 & KT,N,M,0, NITER,NITSM,ICOUNT,ERROR, sum_time, sum_lp_time)


      ! plot not yet finished
!CALL PLOT(PT, P0, U(:,:,0), V(:,:,0),HX,HY,N,M,IRHW)


! CLOSE INITIAL CONDITIONS *************************************

! COMPUTE SOLUTION IN TIME *************************************

IF(IANAL.EQ.0) THEN
  DO KT=1,NT
    if(int(float(kt)/float(mpfl))*mpfl.eq.kt) then
      liner=1
    else
      liner=0
    endif
    !write(*,*) kt, int(float(kt)/float(mpfl))*mpfl, liner
    IPRINT=0
    IF(KT/NPRINT*NPRINT.EQ.KT) IPRINT=1

    ! COMPUTE ADVECTIVE COURANT NUMBERS
    ! COMPUTE VELOCITY PREDICTOR
    CALL VELPRD(U,V,F1,F2,PD,HX,HY,IP,N,M,GC1,GC2,EP, KMX)
      
  !IF(IPRINT.EQ.1) then
  ! call divergence (U,V) for exit condition
    Call DIVER(D_Adv(:,:),U(:,:,1),V(:,:,1),HX,HY,S,N,M,IP,1)

    Exit_Cond=maxval(ABS(D_Adv(:,:)))
    !do J=1,M
    !  write(*,*) J, maxval(ABS(D_Adv(:,J)))*D0*G_23*DT_23
    !enddo
    Exit_Cond=2.0d0*Exit_Cond*D0*G_23*DT_23
    !write(*,*) 'Exit_Cond', Exit_Cond

  !endif 

    DO J=1,M
      DO I=1,N
        U(I,J,1)=U(I,J,1)*HY(I,J)
        V(I,J,1)=V(I,J,1)*HX(I,J)
      end do
    end do
    ! COMPUTE COURANT NUMBERS AT STAGGERED TIME/SPACE POSITIONS
    DO J=1,M
      DO I=2,N-1
        UA(I,J)=(U(I,J,1)+U(I-1,J,1))*GH1
      end do
    end do

    CALL XBC(UA,N,M)
    
    DO I=1,N
      DO J=2,M
        VA(I,J)=(V(I,J,1)+V(I,J-1,1))*GH2
      end do
      VA(I,  1)=rpe_0
      VA(I,M+1)=rpe_0
    end do

   ! adv_cour=0.0d0
   !
   ! DO J=1,M
   !   DO I=1,N
   !     adv_cour=max(adv_cour, ABS(UA(I,J)/HX(I,J)/HY(I,J))+ABS(VA(I,J)/HX(I,J)/HY(I,J)))
   !   end do
   ! end do
   ! write(*,*) kt, adv_cour

! CLOSE ADVECTIVE COURANT NUMBERS

! COLLECT EXPLICIT PARTS OF CONTINUITY AND MOMENTUM EQUATIONS:

! COLLECT FROM CONTINUITY EQUATION
    DO J=1,M
      DO I=1,N
        U(I,J,1)=QX(I,J)*GC1
        V(I,J,1)=QY(I,J)*GC2
      end do
    end do

    CALL DIVER(F1(:,:,1),U(:,:,1),V(:,:,1),HX,HY,S,N,M,IP,1)


!endif
    DO J=1,M
      DO I=1,N
    ! original    F1(I,J,1)=PD(I,J)+P0(I,J)-rpe_05*F1(I,J,1)
        F1(I,J,1)=.5*F1(I,J,1)
      end do
    end do

! C--->                       ADVECTION

    DO J=1,M
      DO I=1,N
        QX(I,J)=QX(I,J)+rpe_05*F1(I,J,0)
        QY(I,J)=QY(I,J)+rpe_05*F2(I,J,0)
!       PC(I,1)=PD(I,1)
      end do
    end do


      DO J=1,M
        DO I=1,N

          PC(I,J)=PD(I,J)

        end do
      end do

    

      ! predict fluid thickness for a second order accurate implicit problem
      CALL MPDATT(UA,VA,PC,S,N,M,IORD,ISOR,NONOS,IDIV,1, IP,liner, PD_T, .FALSE.)

      CALL MPDATT(UA,VA,QX,S,N,M,IORD,ISOR,NONOS,IDIV,-1, IP,liner, QX_T, codesignQ)
      CALL MPDATT(UA,VA,QY,S,N,M,IORD,ISOR,NONOS,IDIV,-1, IP,liner, QY_T, codesignQ)
      

   IF(QRelax) then !! define reference values for polar absorbers
    IF(IRHW==1) then
           
      CALL RHWT(U0,V0,PT0,PD0,P0,COR_23,X_23,Y_23,N,M,F0_23,R_23,KT*DT_23)
       DO J=1,M
        DO I=1,N
         QXS(I,J)=U0(I,J)*PD0(I,J)
         QYS(I,J)=V0(I,J)*PD0(I,J)
        ENDDO
       enddo 
      ENDIF
    else
      !write(*,*) 'relaxation values not defined'
    endif

!--->                CORIOLIS AND METRIC FORCES
    DO J=1,M
      DO I=1,N
        UA(I,J)=QX(I,J)+rpe_05*(rpe_2*E1(I,J,0)-E1(I,J,-1)) &
                     & +rpe_05*Alp_REL(I,J)*QXS(I,J)
        !write(*,*) 'F1',I,ALP_REL(I,J), rpe_05*Alp_REL(I,J)*U0(I,J)*PD0(I,J), QX(I,J) &
        !        & ,U0(I,J),PD0(I,J)
        VA(I,J)=QY(I,J)+rpe_05*(rpe_2*E2(I,J,0)-E2(I,J,-1))&
                     & +rpe_05*Alp_REL(I,J)*QYS(I,J)
        !write(*,*) 'F1',I,ALP_REL(I,J), rpe_05*Alp_REL(I,J)*V0(I,J)*PD0(I,J), QY(I,J) &
        !        & ,V0(I,J),PD0(I,J)
        !read(*,*)
      end do   
   end do
   
   DO J=1,M
      DO I=1,N
       AMM=1.0d0+0.5d0*Alp_REL(I,J)
       GMM=0.5d0*COR(I,J)
       DETI=1.0d0/(AMM**2+GMM**2)

       QX(I,J)=(AMM*DETI*UA(I,J)+GMM*DETI*VA(I,J))*GC1
       QY(I,J)=(AMM*DETI*VA(I,J)-GMM*DETI*UA(I,J))*GC2

       
      end do
   end do


      DO J=1,M
        DO I=1,N

            QX_HP(I,J)= QX(I,J)/GC1
            QY_HP(I,J)= QY(I,J)/GC2
        end do
      end do
    
 
    
     DO J=1,M
       DO I=1,N
         E1(I,J,-1)=E1(I,J,0)
         E2(I,J,-1)=E2(I,J,0)
       end do 
     end do
    


  
    
!--->            DIVERGENCE OF ALL COLLECTED TERMS:  !! divergence of incoming mass from explicit part of new momenta QX, QY
    CALL DIVER(F2(:,:,1),QX,QY,HX,HY,S,N,M,IP,1)


    DO J=1,M
      DO I=1,N

  ! Original      F1(I,J,1)=(F1(I,J,1)-rpe_05*F2(I,J,1))*G   !! (h+h0)*g pressure minus incoming mass of 0.5 old and 0.5 new momentum
        F1(I,J,1)=(F1(I,J,1)+.5*F2(I,J,1))*G   !! (h+h0)*g pressure minus incoming mass of 0.5 old and 0.5 new momentum
        F2(I,J,1)=PT(I,J)                          !! old one without new contributions
        PC(I,J)=(PC(I,J)+P0(I,J))*G
        PD(I,J)=P0(I,J)*G
      end do
    end do

   ! CALL PRFORC(PT,E1(:,:,0),E2(:,:,0),PD,F2(:,:,1), &     !t^n estimate, O(dt/2)
   !         &   F1(:,:,0),F2(:,:,0),HX,HY,COR,N,M,IP,GC1,GC2,0,0)

    CALL PRFORC(PC,E1(:,:,0),E2(:,:,0),PD,F2(:,:,1), &      !t^{n+1}, fully O(dt^2)
            &   F1(:,:,0),F2(:,:,0),HX,HY,COR,N,M,IP,GC1,GC2,0,0)

   
  ! end gcrtest1

! COMPUTE FIRST GUESS FROM ADVECTION

  CALL  GCR_PRE(PT,F1(:,:,0),F2(:,:,0),HX,HY,S,S_full,F1(:,:,1),F2(:,:,1), &
       &      PD(:,:),E1(:,:,0),E2(:,:,0),COR,IP, &
       &      U(:,:,0),U(:,:,1),V(:,:,0),V(:,:,1),N,M,GC1,GC2,   &
           &    MGH1IHX, MGH2IHY, AC, BC, AD, BD,  &
  &      niter,nitsm,icount,error, PD_T, sum_time, sum_lp_time, ID_PREC,.FALSE., save_time,&
       &      TIME, codesignQ, IRHW, X, Y, Exit_Cond,               &
       &  EXP_NAME, iprint, num_of_bits, DP_Depth, Alp_REL)

     !! go back to REAL(Kind=4) outside of the solver
     !  read(*,*)
 If(save_time) exit  ! if iteration counter has reached 100 once, jump to next set of bits
  
   

    
   IF ( codesignD) then

      DO J=1,M
        DO I=1,N
            PT_HP(I,J)= PT_HP(I,J)+ PD_T(I,J)

        end do
      end do
      

      
      DO J=1,M
        DO I=1,N
            PT(I,J)= PT_HP(I,J)

        end do
      end do
    

     end if
       
       
  If(QRelax) then     
    CALL PRFORC_ABS(PT,F1(:,:,0),F2(:,:,0),PD,F2(:,:,1), &
       &      E1(:,:,0),E2(:,:,0),HX,HY,COR,N,M,IP,GC1,GC2,Alp_REL,1,1)

  else
    CALL PRFORC(PT,F1(:,:,0),F2(:,:,0),PD,F2(:,:,1), &
      &      E1(:,:,0),E2(:,:,0),HX,HY,COR,N,M,IP,GC1,GC2,1,1)

  endif

       
              !!! update PD, because PRFORC needs an intermediate value of PD
    IF (.not. codesignD) then
    
   
      
      DO J=1,M
        DO I=1,N
          PD(I,J)= max(EP, PT(I,J)*GI-P0(I,J))
   
        end do
      end do
   
        
      DO J=1,M
        DO I=1,N
            PD_HP(I,J)= PD(I,J)
   
        end do
      end do
    
    else
    
      
      DO J=1,M
        DO I=1,N
            PD_HP(I,J)= max(EP, PT_HP(I,J)*GI_23-P0_HP(I,J))
      !write(*,*) PD_HP(I,J)
      !read(*,*)
        end do
      end do
      

      
      DO J=1,M
        DO I=1,N
            !PT(I,J)= (PD_HP(I,J)+P0(I,J))*G
            PD(I,J)= PD_HP(I,J)
            PD_T(I,J)= rpe_0
        end do
      end do
    

     end if
     
     
! COMPUTE SOLUTION'S UPDATE

    DO J=1,M
      DO I=1,N
        QX(I,J)=(QX(I,J)+F1(I,J,0)*GI)/GC1
        QY(I,J)=(QY(I,J)+F2(I,J,0)*GI)/GC2

      end do
    end do

!! sanity check for elliptic problem
CALL DIVER(div_old(:,:),QX_old*GC1,QY_old*GC2,HX,HY,S,N,M,IP,1)
CALL DIVER(div_new(:,:),QX *GC1   ,QY *GC2   ,HX,HY,S,N,M,IP,1)


r_eval(:,:)=(0.5d0*(div_new(:,:)+div_old(:,:))*G_23-PD_old(:,:)+PT(:,:))

if (.not. (KT/NPRINT*NPRINT.NE.KT)) then
!write(*,*) div_old(1,1),div_new(1,1), PD_old(1,1), PT(1,1)
!write(*,*) 0.5d0*(div_old(1,1)+div_new(1,1))*G_23, -PD_old(1,1)+ PT(1,1)
!write(*,*) 0.5d0*(div_old(1,1)+div_new(1,1))*G_23 -PD_old(1,1)+ PT(1,1)

!write(*,*) 'r EVAL'
!  write(*,*) (maxval(abs(r_eval(:,J))), J=1,M)
!read(*,*)

endif
    

!CALL FILTRQ(QX,N,M)
!CALL FILTRQ(QY,N,M)

PD_old(:,:)=PT(:,:)
QX_old(:,:)=QX(:,:)
QY_old(:,:)=QY(:,:)

!end sanity check for elliptic problem
! COMPUTE NEW FORCES


    
! COMPUTE NEW FORCES



    CALL PRFC0(PT,F1(:,:,0),F2(:,:,0),PD,HX,HY,IP,IPS,GH1,GH2,EP,N,M)
    
    IF(IRHW.EQ.-1) THEN
      CALL SMOOTHSTATE(QX,QXS,E1(1,1,0),E2(1,1,0),IP,-1.,KT,N,M)
      CALL SMOOTHSTATE(QY,QYS,E1(1,1,0),E2(1,1,0),IP,-1.,KT,N,M)
!     CALL POLARFILT(QX,HX,HY,FZ,SCR,N,M)
!     CALL POLARFILT(QY,HX,HY,FZ,SCR,N,M)
    ENDIF

    DO J=1,M
      DO I=1,N
        E1(I,J,0)=-DHX2Y(I,J)*QX(I,J)*QY(I,J)/PD(I,J)
        E2(I,J,0)= DHX2Y(I,J)*QX(I,J)*QX(I,J)/PD(I,J)
        !write(*,*) 'E1',I,DHX2Y(I,J), E1(I,J,0), E2(I,J,0)
        !read(*,*)
      end do
    end do

    DO J=1,M
      DO I=1,N
        F1(I,J,0)=F1(I,J,0)+COR(I,J)*QY(I,J)+E1(I,J,0)     &
                   &       -ALP_REL(I,J)*(QX(I,J)-QXS(I,J))
        !if(kt>2 .and. J>=M-1) then
        !write(*,*)kt, J
        !write(*,*) 'F1',I,ALP_REL(I,J), -ALP_REL(I,J)*(QX(I,J)-QXS(I,J)), QX(I,J) &
        !        & ,U0(I,J),PD0(I,J), PD(I,J), PT(I,J)
        !endif
        F2(I,J,0)=F2(I,J,0)-COR(I,J)*QX(I,J)+E2(I,J,0)     &
                   &       -ALP_REL(I,J)*(QY(I,J)-QYS(I,J))
        !if(kt>2 .and. J>=M-1) then
        !write(*,*) 'F2',I,ALP_REL(I,J), -ALP_REL(I,J)*(QY(I,J)-QYS(I,J)), QY(I,J) &
        !        & ,V0(I,J),PD0(I,J)
        !read(*,*)
        !endif
      end do
    end do

    DO J=1,M
      DO I=1,N

        U(I,J,0)=QX(I,J)/PD(I,J)
        V(I,J,0)=QY(I,J)/PD(I,J)
      end do
    end do


!COMPUTE OUTPUTED FIELDS ****************************************

    IF(.not. (KT/NPRINT*NPRINT.NE.KT)) then
    
    
        DO J=1,M
          DO I=1,N

            QX_HP(I,J)= QX(I,J)
            QY_HP(I,J)= QY(I,J)
          end do
        end do
    

    
      IF(IWRITE.EQ.1) WRITE(9) PD,PT,QX,QY,U,V,F1,F2,E1,E2
      CALL DIAGNOS(QX_HP(:,:)/PD_HP(:,:),QY_HP(:,:)/PD_HP(:,:),PD_HP,&
                 & PT,HX,HY,IP,S,TIME,DX,DY,DT, SUM0,SUM1, &
                 & KT,N,M,1, NITER,NITSM,ICOUNT,ERROR, sum_time, sum_lp_time)
      ! plot not yet finished
      call write_perf_markers (PD_HP(:,:)+P0(:,:),QX_HP(:,:)/PD_HP(:,:),QY_HP(:,:)/PD_HP(:,:)&
                  &, TIME, codesignQ, codesignD, IRHW, X, Y, N, M, &
                  & num_of_bits,NITER,NITSM,ICOUNT, sum_time, sum_lp_time)
      call write_fields(PD_HP(:,:)+P0(:,:),QX_HP(:,:)/PD_HP(:,:),QY_HP(:,:)/PD_HP(:,:)&
                  &, TIME, codesignQ, codesignD, IRHW, X, Y, N, M, num_of_bits, ID_prec, EXP_NAME)
      !CALL PLOT(PT, P0, U(:,:,0), V(:,:,0),HX,HY,N,M,IRHW)
    end if

  end do
!!CLOSE TIME INTEGRATION

ENDIF

IF(IANAL.EQ.1) THEN
  
  TIME=0.0
  NPLOT=5
 
  DO KF=1,NFIL
    TIME=TIME+DTFIL(KF)*NTFIL(KF)/3600.0d0
    PRINT 300, TIME      
    300 FORMAT(14X,5HTIME=,F7.2,6H HOURS)
    READ(10) PD,PT,QX,QY,U,V,E1,E2,F1,F2
    IF(KF/NPLOT*NPLOT.EQ.KF) then
      ! plot not yet finished
      !CALL PLOT(PT, P0, U(:,:,0), V(:,:,0),HX,HY,N,M,IRHW)
    end if
  end do
ENDIF

endif !!spare bits from 20 to 50
end do !! end loop over all number of bits

call close_perf_markers
  end do
 end do
enddo

call cpu_time(finish)

write(*,*) 'Total runtime', finish-start
END program


subroutine topogr(h0,x,y,n,m, mountain)
use implicit_functions_SP

implicit none

INTEGER :: n, m
REAL(Kind=4) :: h0(n,m)
REAL(Kind=4) :: x(n),y(m)
REAL(Kind=4) :: hs0, Rad, x_c, y_c, dist, pi, sigma
integer :: i, j
logical :: mountain

pi = acos(-1.0d0)

hs0=1800.0d0 ! 2000.0d0
Rad=pi/9.0d0
x_c=1.0d0*pi/2.0d0
y_c=pi/8.0d0!pi/6.0d0


sigma=Rad/2.150d0

! mountain shape functions
!     dist(xln,ylt) =sqrt( (cos(ylt)*sin((xln-x0)/2))**2
!    .                      +       (sin((ylt-y0)/2))**2 )/xl
!     profm(rad)=.5*(1.+cos(pi*rad))
! special for the conical mountain on the pole ccccc

!dist(xln,ylt)=abs(ylt-y0)
!profm(rad)=amax1(0., 1.-gamm*rad/rnot)
!rnot = 2*acos(-1.)/128. * 10.
!gamm=1.
!cccccccccccccccccccccccccccccccccccccccccccccccccccc

!pi = acos(-1.)
!xl=dx * 5.
!x0 = pi
!y0 = pi*.5

do j=1,m
  do i=1,n
    h0(i,j)=0.0d0
  end do
end do

If (mountain) then
  
  do j=1,m
    do i=1,n

      dist=min(Rad, sqrt( (x(i)-x_c)**2 + (y(j)-y_c)**2) )
      h0(i,j)=hs0*(rpe_1-dist/Rad)

    end do

  end do

!  do j=1,m
!    do i=1,n
!      dist= sqrt( (x(i)-x_c)**2 + (y(j)-y_c)**2) 
!      h0(i,j)=0.4d0*hs0*(1.0d0/(sigma*sqrt(2*pi))*exp(-0.5d0*(dist/sigma)**2) )
!    end do
!  end do

end if

end subroutine

!!!!!!!!! INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INITZON(U,V,PT,COR,X,Y,N,M,F0,BETA,H00,R,Q)
use implicit_functions_SP

implicit none

REAL(Kind=4) :: PT(N,M)
REAL(Kind=4) :: U(N,M),V(N,M),COR(N,M),X(N),Y(M), F0, beta, H00, R, Q
INTEGER :: N, M


INTEGER :: I, J
 

DO J=1,M
  DO I=1,N
    COR(I,J)=F0*(-COS(X(I))*COS(Y(J))*SIN(BETA)+SIN(Y(J))*COS(BETA))
    PT(I,J)=H00-R**2.0*(F0+Q)*0.5*Q*  &
      &     (-COS(X(I))*COS(Y(J))*SIN(BETA)+SIN(Y(J))*COS(BETA))**2
  end do
end do

DO J=1,M
  DO I=1,N
    U(I,J)=Q*(COS(BETA)+TAN(Y(J))*COS(X(I))*SIN(BETA))*R*COS(Y(J))
    V(I,J)=-Q*SIN(X(I))*SIN(BETA)*R
  end do
end do
      write (6,*)  'initzon called'

END SUBROUTINE

SUBROUTINE INITRHW(U,V,PT,COR,X,Y,N,M,F0,A)
use implicit_functions_SP

implicit none
REAL(Kind=4) :: PT(N,M)

REAL(Kind=4) ::U(N,M),V(N,M),F0, A,COR(N,M),X(N),Y(M)
INTEGER :: N, M


REAL(Kind=4) ::     ATH(M), BTH(M), CTH(M), TH
REAL(Kind=4) ::     OM,K,PH0
INTEGER :: R, I, J

OM=7.848E-6
K=7.848E-6
R=4
PH0= 78.4E3
!   DATA OM,K,R,PH0/7.848E-6,7.848E-6,4,0./

DO J=1,M
  TH=Y(J)
  ATH(J)=OM*0.5*(F0+OM)*(COS(TH))**2                          &
     &  +0.25*K**2*(COS(TH))**(2*R)*( (R+1)*(COS(TH))**2      &
     &   +FLOAT(2*R**2-R-2)-2.*R**2/(COS(TH))**2 )
  BTH(J)=(F0+2.*OM)*K/FLOAT((R+1)*(R+2))*(COS(TH))**R         &
     &       *( FLOAT(R**2+2*R+2)-((R+1)*COS(TH))**2 )
  CTH(J)=0.25*K**2*(COS(TH))**(2*R)*( FLOAT(R+1)*(COS(TH))**2 &
     &       -FLOAT(R+2) )  
end do

DO J=1,M
  DO I=1,N
    COR(I,J)=F0*SIN(Y(J))
    U(I,J)=A*OM*COS(Y(J))+A*K*COS(R*X(I))                      &
       &   *(COS(Y(J)))**(R-1)*(R*(SIN(Y(J)))**2-(COS(Y(J)))**2)
      V(I,J)=-A*K*R*(COS(Y(J)))**(R-1)*SIN(Y(J))*SIN(R*X(I))
      PT(I,J)=PH0+A**2*ATH(J)+A**2*BTH(J)*COS(R*X(I))      &
       &   +A**2*CTH(J)*COS(2.*R*X(I))
  end do
end do
  write (6,*)  'initrhw called'

END SUBROUTINE


SUBROUTINE PRFC0(A,F1,F2,PD,HX,HY,IP,IPS,GH1,GH2,EP,N,M)
use implicit_functions_SP

implicit none
REAL(Kind=4) :: A(N,M),F1(N,M),F2(N,M),PD(N,M),HX(N,M),HY(N,M)
INTEGER :: IP(N)
REAL(Kind=4) :: GH1, GH2, EP
INTEGER :: IPS, N, M


REAL(Kind=4) :: GP, GN
INTEGER :: I, J

    DO I=2,N-1
       DO  J=1,M
      GP=PD(I+1,J)/(PD(I+1,J)+PD(I-1,J)+EP)*2.*IPS+(1-IPS)
      GN=PD(I-1,J)/(PD(I+1,J)+PD(I-1,J)+EP)*2.*IPS+(1-IPS)
      F1(I,J)=-GH1*PD(I,J)/HX(I,J)* &
           & ( GP*(A(I+1,J)-A(I,J))+GN*(A(I,J)-A(I-1,J)) )
      enddo
      DO J=2,M-1
      GP=PD(I,J+1)/(PD(I,J+1)+PD(I,J-1)+EP)*2.*IPS+(1-IPS)
      GN=PD(I,J-1)/(PD(I,J+1)+PD(I,J-1)+EP)*2.*IPS+(1-IPS)
      F2(I,J)=-GH2*PD(I,J)/HY(I,J)* &
          & ( GP*(A(I,J+1)-A(I,J))+GN*(A(I,J)-A(I,J-1)) )
      enddo
      GP=PD(I    ,2)/(PD(I,2)+PD(IP(I),1)+EP)*2.*IPS+(1-IPS)
      GN=PD(IP(I),1)/(PD(I,2)+PD(IP(I),1)+EP)*2.*IPS+(1-IPS)
      F2(I,1)=-GH2*PD(I,1)/HY(I,1)* &
        &  ( GP*(A(I,2)-A(I,1))+GN*(A(I,1)-A(IP(I),1)) )

      GP=PD(IP(I),M)/(PD(IP(I),M)+PD(I,M-1)+EP)*2.*IPS+(1-IPS)
      GN=PD(I  ,M-1)/(PD(IP(I),M)+PD(I,M-1)+EP)*2.*IPS+(1-IPS)
      F2(I,M)=-GH2*PD(I,M)/HY(I,M)*& 
           &  ( GP*(A(IP(I),M)-A(I,M))+GN*(A(I,M)-A(I,M-1)) )
    enddo

!DO I=2,N-1
!
!  DO J=1,M
!    F1(I,J)=-GH1*PD(I,J)/HX(I,J)*                          &
!       & (A(I+1,J)-A(I-1,J))
!  end do    
!
!  DO J=2,M-1
!    F2(I,J)=-GH2*PD(I,J)/HY(I,J)*                          &
!     &  (A(I,J+1)-A(I,J-1)) 
!  end do  
!
!  F2(I,1)=-GH2*PD(I,1)/HY(I,1)*                            &
!     &  (A(I,2)-A(IP(I),1)) 
!
!  F2(I,M)=-GH2*PD(I,M)/HY(I,M)*                            &
!     &  (A(IP(I),M)-A(I,M-1)) 
!end do
 
CALL XBC(F1,N,M)
CALL XBC(F2,N,M)

END SUBROUTINE


SUBROUTINE PRFORC( P,F1,F2,PB,P0,E1,E2,HX,HY,COR,            &
     &            N,M,IP,GC1,GC2,NOR,IRS)
use implicit_functions_SP

implicit none
REAL(Kind=4) :: P(N,M),F1(N,M),F2(N,M),PB(N,M),P0(N,M),E1(N,M),E2(N,M), &
     & HX(N,M),HY(N,M),COR(N,M)
REAL(Kind=4) :: GC1, GC2
INTEGER  :: IP(N)
INTEGER  :: N, M, NOR, IRS


INTEGER :: I, J, NM

REAL(Kind=4) :: GH1, GH2, UTILD, VTILD, GMM

GH1=rpe_05*GC1
GH2=rpe_05*GC2

DO J=2,M-1
  DO I=2,N-1
    F1(I,J)=-GH1*(P(I+1,J)-P(I-1,J))/HX(I,J)
    F2(I,J)=-GH2*(P(I,J+1)-P(I,J-1))/HY(I,J)

  end do
end do
   
DO I=2,N-1
  F1(I,1)=-GH1*(P(I+1,1)-P(I-1,1))/HX(I,1)
  F1(I,M)=-GH1*(P(I+1,M)-P(I-1,M))/HX(I,M)
  F2(I,1)=-GH2*(P(I,2)-P(IP(I),1))/HY(I,1)
  F2(I,M)=-GH2*(P(IP(I),M)-P(I,M-1))/HY(I,M)
end do

CALL XBC(F1,N,M)
CALL XBC(F2,N,M)

IF(NOR.EQ.1) THEN
!  NM=N*M
  
!  DO I=1,NM
!    UTILD=F1(I,1)*(P0(I,1)-PB(I,1))+(P(I,1)-IRS*P0(I,1))*E1(I,1)
!    VTILD=F2(I,1)*(P0(I,1)-PB(I,1))+(P(I,1)-IRS*P0(I,1))*E2(I,1)
!    GMM=.5*COR(I,1)
!    F1(I,1)=(UTILD+GMM*VTILD)/(1.+GMM**2)*GH1
!    F2(I,1)=(VTILD-GMM*UTILD)/(1.+GMM**2)*GH2
!  end do
  DO J=1,M
    DO I=1,N
      UTILD=F1(I,J)*(P0(I,J)-PB(I,J))+(P(I,J)-IRS*P0(I,J))*E1(I,J)
      VTILD=F2(I,J)*(P0(I,J)-PB(I,J))+(P(I,J)-IRS*P0(I,J))*E2(I,J)
      GMM=rpe_05*COR(I,J)
      F1(I,J)=(UTILD+GMM*VTILD)/(rpe_1+GMM**2)*GH1
      F2(I,J)=(VTILD-GMM*UTILD)/(rpe_1+GMM**2)*GH2
    end do
  end do
CALL XBC(F1,N,M)
CALL XBC(F2,N,M)
ENDIF

END SUBROUTINE


SUBROUTINE PRFORC_FIN( P,F1,F2,PB,P0,E1,E2,HX,HY,COR,       &
     &            N,M,IP,GC1,GC2,Relax_M, DT_23, NOR,IRS)
use implicit_functions_SP

implicit none
REAL(Kind=4) :: P(N,M),F1(N,M),F2(N,M),PB(N,M),P0(N,M),E1(N,M),E2(N,M), &
     & HX(N,M),HY(N,M),COR(N,M)
REAL(Kind=4) :: GC1, GC2
REAL(Kind=4) :: Relax_M(M), DT_23, Relaxation
INTEGER  :: IP(N)
INTEGER  :: N, M, NOR, IRS


INTEGER :: I, J, NM

REAL(Kind=4) :: GH1, GH2, UTILD, VTILD, GMM

GH1=rpe_05*GC1
GH2=rpe_05*GC2

DO J=2,M-1
  DO I=2,N-1
    F1(I,J)=-GH1*(P(I+1,J)-P(I-1,J))/HX(I,J)
    F2(I,J)=-GH2*(P(I,J+1)-P(I,J-1))/HY(I,J)
  end do
end do
   
DO I=2,N-1
  F1(I,1)=-GH1*(P(I+1,1)-P(I-1,1))/HX(I,1)
  F1(I,M)=-GH1*(P(I+1,M)-P(I-1,M))/HX(I,M)
  F2(I,1)=-GH2*(P(I,2)-P(IP(I),1))/HY(I,1)
  F2(I,M)=-GH2*(P(IP(I),M)-P(I,M-1))/HY(I,M)
end do

CALL XBC(F1,N,M)
CALL XBC(F2,N,M)

IF(NOR.EQ.1) THEN
!  NM=N*M
  
!  DO I=1,NM
!    UTILD=F1(I,1)*(P0(I,1)-PB(I,1))+(P(I,1)-IRS*P0(I,1))*E1(I,1)
!    VTILD=F2(I,1)*(P0(I,1)-PB(I,1))+(P(I,1)-IRS*P0(I,1))*E2(I,1)
!    GMM=.5*COR(I,1)
!    F1(I,1)=(UTILD+GMM*VTILD)/(1.+GMM**2)*GH1
!    F2(I,1)=(VTILD-GMM*UTILD)/(1.+GMM**2)*GH2
!  end do
  DO J=1,M
    DO I=1,N
      UTILD=F1(I,J)*(P0(I,J)-PB(I,J))+(P(I,J)-IRS*P0(I,J))*E1(I,J)
      VTILD=F2(I,J)*(P0(I,J)-PB(I,J))+(P(I,J)-IRS*P0(I,J))*E2(I,J)
      GMM=rpe_05*COR(I,J)
      Relaxation=1.0d0+0.5d0*DT_23*Relax_M(J)
      F1(I,J)=(UTILD+GMM/Relaxation*VTILD)/(Relaxation+GMM**2/Relaxation)*GH1
      F2(I,J)=(VTILD-GMM/Relaxation*UTILD)/(Relaxation+GMM**2/Relaxation)*GH2
    end do
  end do
CALL XBC(F1,N,M)
CALL XBC(F2,N,M)
ENDIF

END SUBROUTINE



SUBROUTINE DIVER(F, U,  V, HX,HY,S, N, M,IP,IFLG)
!!         (r,pfx,pfy,hx,hy,s,n1,n2,ip,-1)
use implicit_functions_SP
 
 implicit none
 
REAL(Kind=4) :: F(N,M),U(N,M),V(N,M),HX(N,M),HY(N,M),S(N,M)
INTEGER :: IP(N)
INTEGER :: N, M, IFLG


INTEGER:: I, J

DO J=2,M-1
  DO I=2,N-1
    F(I,J)= HY(I+1,J)*U(I+1,J)-HY(I-1,J)*U(I-1,J)        & 
        &  +HX(I,J+1)*V(I,J+1)-HX(I,J-1)*V(I,J-1)
  end do
end do    

! original
DO I=2,N-1
  F(I,1)= HY(I+1,1)*U(I+1,1)-HY(I-1,1)*U(I-1,1) &
      &  +(HX(I,2)*V(I,2)+HX(I,1)*V(I,1))
  F(I,M)= HY(I+1,M)*U(I+1,M)-HY(I-1,M)*U(I-1,M)  &
      &  -(HX(I,M)*V(I,M)+HX(I,M-1)*V(I,M-1))
end do
! end original



DO J=1,M
  DO I=2,N-1
    F(I,J)=rpe_05*IFLG*F(I,J)/S(I,J)
  end do
end do

DO J=1,M
  F(1,J)=F(N-1,J)
  F(N,J)=F(2  ,J)
end do

END SUBROUTINE


SUBROUTINE PRFORC_depth( P,F1,F2,PB,P0,E1,E2,HX,HY,COR,            &
     &            N,M,IP,GC1,GC2,NOR,IRS, num_of_bits, DP_Depth)
use implicit_functions_SP

implicit none
REAL(Kind=4) :: P(N,M),F1(N,M),F2(N,M),PB(N,M),P0(N,M),E1(N,M),E2(N,M), &
     & HX(N,M),HY(N,M),COR(N,M)
REAL(Kind=4) :: GC1, GC2
INTEGER  :: IP(N)
INTEGER  :: N, M, NOR, IRS, DP_Depth, num_of_bits


INTEGER :: I, J, NM

REAL(Kind=4) :: GH1, GH2, UTILD, VTILD, GMM

GH1=rpe_05*GC1
GH2=rpe_05*GC2

If (DP_Depth<=0) then
 DO J=2,M-1
   DO I=2,N-1
     F1(I,J)=-GH1*(P(I+1,J)-P(I-1,J))/HX(I,J)
     F2(I,J)=-GH2*(P(I,J+1)-P(I,J-1))/HY(I,J)
   end do
 end do
 
 DO I=2,N-1
   F1(I,1)=-GH1*(P(I+1,1)-P(I-1,1))/HX(I,1)
   F1(I,M)=-GH1*(P(I+1,M)-P(I-1,M))/HX(I,M)
   F2(I,1)=-GH2*(P(I,2)-P(IP(I),1))/HY(I,1)
   F2(I,M)=-GH2*(P(IP(I),M)-P(I,M-1))/HY(I,M)
 end do

CALL XBC(F1,N,M)
CALL XBC(F2,N,M)

else

      DO J=2+(DP_Depth-1),M-1-(DP_Depth-1)
        DO I=2,N-1
         F1(I,J)=-GH1*(P(I+1,J)-P(I-1,J))/HX(I,J)
         F2(I,J)=-GH2*(P(I,J+1)-P(I,J-1))/HY(I,J)
        enddo
      enddo

      DO J=2,DP_Depth
        DO I=2,N-1
         F1(I,J)=-GH1*(P(I+1,J)-P(I-1,J))/HX(I,J)
         F2(I,J)=-GH2*(P(I,J+1)-P(I,J-1))/HY(I,J)
        enddo
      enddo

      DO J=M+1-DP_Depth,M-1
        DO I=2,N-1
         F1(I,J)=-GH1*(P(I+1,J)-P(I-1,J))/HX(I,J)
         F2(I,J)=-GH2*(P(I,J+1)-P(I,J-1))/HY(I,J)
        enddo
      enddo

      DO I=2,N-1
       F1(I,1)=-GH1*(P(I+1,1)-P(I-1,1))/HX(I,1)
       F1(I,M)=-GH1*(P(I+1,M)-P(I-1,M))/HX(I,M)
       F2(I,1)=-GH2*(P(I,2)-P(IP(I),1))/HY(I,1)
       F2(I,M)=-GH2*(P(IP(I),M)-P(I,M-1))/HY(I,M)
      end do

CALL XBC(F1,N,M)
CALL XBC(F2,N,M)

endif



IF(NOR.EQ.1) THEN
!  NM=N*M
  
!  DO I=1,NM
!    UTILD=F1(I,1)*(P0(I,1)-PB(I,1))+(P(I,1)-IRS*P0(I,1))*E1(I,1)
!    VTILD=F2(I,1)*(P0(I,1)-PB(I,1))+(P(I,1)-IRS*P0(I,1))*E2(I,1)
!    GMM=.5*COR(I,1)
!    F1(I,1)=(UTILD+GMM*VTILD)/(1.+GMM**2)*GH1
!    F2(I,1)=(VTILD-GMM*UTILD)/(1.+GMM**2)*GH2
!  end do

!orig
!  DO J=1,M
!    DO I=1,N
!      UTILD=F1(I,J)*(P0(I,J)-PB(I,J))+(P(I,J)-IRS*P0(I,J))*E1(I,J)
!      VTILD=F2(I,J)*(P0(I,J)-PB(I,J))+(P(I,J)-IRS*P0(I,J))*E2(I,J)
!      GMM=rpe_05*COR(I,J)
!      F1(I,J)=(UTILD+GMM*VTILD)/(rpe_1+GMM**2)*GH1
!      F2(I,J)=(VTILD-GMM*UTILD)/(rpe_1+GMM**2)*GH2
!    end do
!  end do
!end orig

  !! rewritten to change precision at poles as desired    
      DO J=1+DP_Depth,M-DP_Depth
        DO I=1,N
      UTILD=F1(I,J)*(P0(I,J)-PB(I,J))+(P(I,J)-IRS*P0(I,J))*E1(I,J)
      VTILD=F2(I,J)*(P0(I,J)-PB(I,J))+(P(I,J)-IRS*P0(I,J))*E2(I,J)
      GMM=rpe_05*COR(I,J)
      F1(I,J)=(UTILD+GMM*VTILD)/(rpe_1+GMM**2)*GH1
      F2(I,J)=(VTILD-GMM*UTILD)/(rpe_1+GMM**2)*GH2
        enddo
      enddo

      DO J=1,DP_Depth
        DO I=1,N
      UTILD=F1(I,J)*(P0(I,J)-PB(I,J))+(P(I,J)-IRS*P0(I,J))*E1(I,J)
      VTILD=F2(I,J)*(P0(I,J)-PB(I,J))+(P(I,J)-IRS*P0(I,J))*E2(I,J)
      GMM=rpe_05*COR(I,J)
      F1(I,J)=(UTILD+GMM*VTILD)/(rpe_1+GMM**2)*GH1
      F2(I,J)=(VTILD-GMM*UTILD)/(rpe_1+GMM**2)*GH2
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=1,N
      UTILD=F1(I,J)*(P0(I,J)-PB(I,J))+(P(I,J)-IRS*P0(I,J))*E1(I,J)
      VTILD=F2(I,J)*(P0(I,J)-PB(I,J))+(P(I,J)-IRS*P0(I,J))*E2(I,J)
      GMM=rpe_05*COR(I,J)
      F1(I,J)=(UTILD+GMM*VTILD)/(rpe_1+GMM**2)*GH1
      F2(I,J)=(VTILD-GMM*UTILD)/(rpe_1+GMM**2)*GH2
        enddo
      enddo

CALL XBC(F1,N,M)
CALL XBC(F2,N,M)



ENDIF

END SUBROUTINE


SUBROUTINE DIVER_depth(F, U,  V, HX,HY,S, N, M,IP,IFLG, num_of_bits, DP_Depth)
!!         (r,pfx,pfy,hx,hy,s,n1,n2,ip,-1)
use implicit_functions_SP
 
 implicit none
 
REAL(Kind=4) :: F(N,M),U(N,M),V(N,M),HX(N,M),HY(N,M),S(N,M)

INTEGER :: IP(N)
INTEGER :: N, M, IFLG, DP_Depth, num_of_bits


INTEGER:: I, J

!DO J=2,M-1
!  DO I=2,N-1
!    F(I,J)= HY(I+1,J)*U(I+1,J)-HY(I-1,J)*U(I-1,J)        & 
!        &  +HX(I,J+1)*V(I,J+1)-HX(I,J-1)*V(I,J-1)
!  end do
!end do    

!DO I=2,N-1
!  F(I,1)= HY(I+1,1)*U(I+1,1)-HY(I-1,1)*U(I-1,1) &
!      &  +(HX(I,2)*V(I,2)+HX(I,1)*V(I,1))
!  F(I,M)= HY(I+1,M)*U(I+1,M)-HY(I-1,M)*U(I-1,M)  &
!      &  -(HX(I,M)*V(I,M)+HX(I,M-1)*V(I,M-1))
!end do

If (DP_Depth<=0) then
DO J=2,M-1
  DO I=2,N-1
    F(I,J)= HY(I+1,J)*U(I+1,J)-HY(I-1,J)*U(I-1,J)        & 
        &  +HX(I,J+1)*V(I,J+1)-HX(I,J-1)*V(I,J-1)
  end do
end do    


DO I=2,N-1
  F(I,1)= HY(I+1,1)*U(I+1,1)-HY(I-1,1)*U(I-1,1) &
      &  +(HX(I,2)*V(I,2)+HX(I,1)*V(I,1))
  F(I,M)= HY(I+1,M)*U(I+1,M)-HY(I-1,M)*U(I-1,M)  &
      &  -(HX(I,M)*V(I,M)+HX(I,M-1)*V(I,M-1))
end do

else

      DO J=2+(DP_Depth-1),M-1-(DP_Depth-1)
        DO I=2,N-1
          F(I,J)= HY(I+1,J)*U(I+1,J)-HY(I-1,J)*U(I-1,J)        & 
          &  +HX(I,J+1)*V(I,J+1)-HX(I,J-1)*V(I,J-1)
        enddo
      enddo


      DO J=2,DP_Depth
        DO I=2,N-1
          F(I,J)= HY(I+1,J)*U(I+1,J)-HY(I-1,J)*U(I-1,J)        & 
          &  +HX(I,J+1)*V(I,J+1)-HX(I,J-1)*V(I,J-1)
        enddo
      enddo

      DO J=M+1-DP_Depth,M-1
        DO I=2,N-1
          F(I,J)= HY(I+1,J)*U(I+1,J)-HY(I-1,J)*U(I-1,J)        & 
          &  +HX(I,J+1)*V(I,J+1)-HX(I,J-1)*V(I,J-1)
        enddo
      enddo

      DO I=2,N-1
        F(I,1)= HY(I+1,1)*U(I+1,1)-HY(I-1,1)*U(I-1,1) &
         &  +(HX(I,2)*V(I,2)+HX(I,1)*V(I,1))
        F(I,M)= HY(I+1,M)*U(I+1,M)-HY(I-1,M)*U(I-1,M)  &
         &  -(HX(I,M)*V(I,M)+HX(I,M-1)*V(I,M-1))
      end do

endif



!DO J=1,M
!  DO I=2,N-1
!    F(I,J)=rpe_05*IFLG*F(I,J)/S(I,J)
!  end do
!end do

  !! rewritten to change precision at poles as desired    
      DO J=1+DP_Depth,M-DP_Depth
        DO I=2,N-1
          F(I,J)=rpe_05*IFLG*F(I,J)/S(I,J)
        enddo
      enddo

      DO J=1,DP_Depth
        DO I=2,N-1
          F(I,J)=rpe_05*IFLG*F(I,J)/S(I,J)
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=2,N-1
          F(I,J)=rpe_05*IFLG*F(I,J)/S(I,J)
        enddo
      enddo
DO J=1,M
  F(1,J)=F(N-1,J)
  F(N,J)=F(2  ,J)
end do

   !! end of rewrite




END SUBROUTINE

SUBROUTINE LAP0_Piotr(A11,A12,A21,A22,B11,B22,    &
     &          PB,P0,E1,E2,HX,HY,COR,ALP_rel,N,M,GC1,GC2)
use implicit_functions_SP

implicit none
REAL(Kind=4) :: A11(N,M),A12(N,M),A21(N,M),A22(N,M),B11(N,M),B22(N,M),  &
     &      PB(N,M),P0(N,M),E1(N,M),E2(N,M),HX(N,M),HY(N,M),COR(N,M), &
     &      ALP_rel(N,M)
REAL(Kind=4) :: GC1, GC2, DETI, AMM, GMM
INTEGER :: N, M


REAL(Kind=4) :: GH1, GH2, C1, C2,  A, B, C, D
INTEGER :: I, J


      GH1=.5*GC1
      GH2=.5*GC2

DO J=1,M
  DO I=1,N
      C1=-GH1/HX(I,J)*(P0(I,J)-PB(I,J))
      C2=-GH2/HY(I,J)*(P0(I,J)-PB(I,J))
      GMM=.5*COR(I,J)
      AMM=1.+.5*ALP_rel(I,J)
      DETI=1./(AMM**2+GMM**2)
      A=AMM*DETI
      B=GMM*DETI
      A11(I,J)=-C1*A*GH1*HY(I,J)*.5
      A12(I,J)=-C2*B*GH1*HY(I,J)*.5
      A22(I,J)= C1*B*GH2*HX(I,J)*.5
      A21(I,J)=-C2*A*GH2*HX(I,J)*.5
      B11(I,J)=-   A*GH1*E1(I,J)*HY(I,J)*.5 &
      &        -   B*GH1*E2(I,J)*HY(I,J)*.5
      B22(I,J)=-   A*GH2*E2(I,J)*HX(I,J)*.5 &
      &        +   B*GH2*E1(I,J)*HX(I,J)*.5
  end do
ENDDO
end subroutine

SUBROUTINE LAP0(A11,A12,A21,A22,B11,B22,    &
     &          PB,P0,E1,E2,HX,HY,COR,N,M,GC1,GC2)
use implicit_functions_SP

implicit none
REAL(Kind=4) :: A11(N,M),A12(N,M),A21(N,M),A22(N,M),B11(N,M),B22(N,M),  &
     &      PB(N,M),P0(N,M),E1(N,M),E2(N,M),HX(N,M),HY(N,M),COR(N,M)
REAL(Kind=4) :: GC1, GC2
INTEGER :: N, M


REAL(Kind=4) :: GH1, GH2, C1, C2, GMM, A, B, C, D
INTEGER :: I, J


GH1=rpe_05*GC1
GH2=rpe_05*GC2

DO J=1,M
  DO I=1,N
    C1=-GH1/HX(I,J)*(P0(I,J)-PB(I,J))
    C2=-GH2/HY(I,J)*(P0(I,J)-PB(I,J))
    GMM=rpe_05*COR(I,J)
    A= rpe_1/(rpe_1+GMM*GMM)
    B=GMM*A
    C=GH1*HY(I,J)*rpe_05
    D=GH2*HX(I,J)*rpe_05
    A11(I,J)=-C1*A*C
    A12(I,J)=-C2*B*C
    A21(I,J)=-C2*A*D
    A22(I,J)= C1*B*D
    B11(I,J)=-   A*E1(I,J)*C                  &
       &     -   B*E2(I,J)*C
    B22(I,J)=-   A*E2(I,J)*D                  &
       &     +   B*E1(I,J)*D
  end do
ENDDO
end subroutine

SUBROUTINE LAP0_depth(A11,A12,A21,A22,B11,B22,  &
           &    PB,P0,E1,E2,HX,HY,COR,N,M,GC1,GC2, &
           &    MGH1IHX, MGH2IHY, AC, BC, AD, BD,  &
           & DP_Depth)
use implicit_functions_SP

implicit none
REAL(Kind=4) :: A11(N,M),A12(N,M),A21(N,M),A22(N,M),B11(N,M),B22(N,M),  &
     &      PB(N,M),P0(N,M),E1(N,M),E2(N,M),HX(N,M),HY(N,M),COR(N,M)

REAL(Kind=4) :: MGH1IHX(M), MGH2IHY(M), AC(M), BC(M), AD(M), BD(M)
REAL(Kind=4) ::GC1, GC2
REAL(Kind=4) :: MGH1IHX_L(M), MGH2IHY_L(M)
INTEGER :: N, M, DP_Depth


REAL(Kind=4) :: GH1, GH2, C1, C2, GMM, A, B, C, D
INTEGER :: I, J



DO J=1+DP_Depth,M-DP_Depth
  AC(J)=AC(J)
  BD(J)=BD(J)
  BC(J)=BC(J)
  AD(J)=AD(J)
  MGH1IHX_L(J)=MGH1IHX(J)
  MGH2IHY_L(J)=MGH2IHY(J)
end do

  !! rewritten to change precision at poles as desired    

      DO J=1,M
        DO I=1,N

    C1=MGH1IHX_L(J)*(P0(I,J)-PB(I,J))
    A11(I,J)=-C1*AC(J)

    A22(I,J)= C1*BD(J)

    C2=MGH2IHY_L(J)*(P0(I,J)-PB(I,J))

    A12(I,J)=-C2*BC(J)
    A21(I,J)=-C2*AD(J)

    B11(I,J)=-   E1(I,J)*AC(J)                  &
       &     -   E2(I,J)*BC(J)
    B22(I,J)=-   E2(I,J)*AD(J)                  &
       &     +   E1(I,J)*BD(J)
     ! write(*,*)I, A11(I,J),A12(I,J),A21(I,J),A22(I,J),B11(I,J),B22(I,J), E1(I,J), E2(I,J)
     ! read(*,*)
        enddo
      enddo



END SUBROUTINE


SUBROUTINE LAPL_depth(P,F,A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP, num_of_bits, DP_Depth)
use implicit_functions_SP

implicit none
REAL(Kind=4) :: P(N,M),F(N,M),A11(N,M),A12(N,M),A21(N,M),A22(N,M),  &
    &      B11(N,M),B22(N,M),P0(N,M),U(N,M),V(N,M),S(N,M),S_L(N,M)
INTEGER :: IP(N)
INTEGER :: N, M, num_of_bits, DP_Depth

REAL(Kind=4) :: UTIL, VTIL
INTEGER :: I, J



DO J=1+DP_Depth,M-DP_Depth
 DO I=1,N
  S_L(I, J)=S(I, J)
end do
end do

If (DP_Depth<=0) then

DO J=2,M-1
  DO I=2,N-1
    U(I,J)=P(I+1,J)-P(I-1,J)
    V(I,J)=P(I,J+1)-P(I,J-1)
  end do
end do
  
DO I=2,N-1
  U(I,1)=P(I+1,1)-P(I-1,1)
  U(I,M)=P(I+1,M)-P(I-1,M)
  V(I,1)=P(I,2)-P(IP(I),1)
  V(I,M)=P(IP(I),M)-P(I,M-1)
ENDDO  
CALL XBC(U,N,M)
CALL XBC(V,N,M)
else


 

      DO J=2+(DP_Depth-1),M-1-(DP_Depth-1)
        DO I=2,N-1
    U(I,J)=P(I+1,J)-P(I-1,J)
    V(I,J)=P(I,J+1)-P(I,J-1)
        enddo
      enddo


      DO J=2,DP_Depth
        DO I=2,N-1
    U(I,J)=P(I+1,J)-P(I-1,J)
    V(I,J)=P(I,J+1)-P(I,J-1)
        enddo
      enddo

      DO J=M+1-DP_Depth,M-1
        DO I=2,N-1
    U(I,J)=P(I+1,J)-P(I-1,J)
    V(I,J)=P(I,J+1)-P(I,J-1)
        enddo
      enddo

      DO I=2,N-1
  U(I,1)=P(I+1,1)-P(I-1,1)
  U(I,M)=P(I+1,M)-P(I-1,M)
  V(I,1)=P(I,2)-P(IP(I),1)
  V(I,M)=P(IP(I),M)-P(I,M-1)
      end do
CALL XBC(U,N,M)
CALL XBC(V,N,M)

endif




!DO J=1,M
!  DO I=1,N
!    UTIL=U(I,J)*A11(I,J)+A12(I,J)*V(I,J)+B11(I,J)*(P(I,J))  ! why only B11*P and not B11*(P-P0)
!    VTIL=V(I,J)*A21(I,J)+A22(I,J)*U(I,J)+B22(I,J)*(P(I,J))
!    U(I,J)=UTIL
!    V(I,J)=VTIL
!  ENDDO
!ENDDO

  !! rewritten to change precision at poles as desired    
      DO J=1+DP_Depth,M-DP_Depth
        DO I=1,N
    UTIL=U(I,J)*A11(I,J)+A12(I,J)*V(I,J)+B11(I,J)*(P(I,J))  ! why only B11*P and not B11*(P-P0)
    VTIL=V(I,J)*A21(I,J)+A22(I,J)*U(I,J)+B22(I,J)*(P(I,J))
    U(I,J)=UTIL
    V(I,J)=VTIL
        enddo
      enddo

      DO J=1,DP_Depth
        DO I=1,N
    UTIL=U(I,J)*A11(I,J)+A12(I,J)*V(I,J)+B11(I,J)*(P(I,J))  ! why only B11*P and not B11*(P-P0)
    VTIL=V(I,J)*A21(I,J)+A22(I,J)*U(I,J)+B22(I,J)*(P(I,J))
    U(I,J)=UTIL
    V(I,J)=VTIL
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=1,N
    UTIL=U(I,J)*A11(I,J)+A12(I,J)*V(I,J)+B11(I,J)*(P(I,J))  ! why only B11*P and not B11*(P-P0)
    VTIL=V(I,J)*A21(I,J)+A22(I,J)*U(I,J)+B22(I,J)*(P(I,J))
    U(I,J)=UTIL
    V(I,J)=VTIL
        enddo
      enddo
CALL XBC(U,N,M)
CALL XBC(V,N,M)

   !! end of rewrite



 



If (DP_Depth<=0) then

DO J=2,M-1
  DO I=2,N-1
    F(I,J)= (U(I+1,J)-U(I-1,J)+V(I,J+1)-V(I,J-1))/S_L(I,J)
  end do
end do    
 
DO I=2,N-1
  F(I,1)= (U(I+1,1)-U(I-1,1)+(V(I,2)+V(I,1)))/S_L(I,1)
  F(I,M)= (U(I+1,M)-U(I-1,M)-(V(I,M)+V(I,M-1)))/S_L(I,M)
ENDDO
CALL XBC(F,N,M)

else

      DO J=2+(DP_Depth-1),M-1-(DP_Depth-1)
        DO I=2,N-1
    F(I,J)= (U(I+1,J)-U(I-1,J)+V(I,J+1)-V(I,J-1))/S_L(I,J)
        enddo
      enddo




      DO J=2,DP_Depth
        DO I=2,N-1
    F(I,J)= (U(I+1,J)-U(I-1,J)+V(I,J+1)-V(I,J-1))/S(I,J)
        enddo
      enddo

      DO J=M+1-DP_Depth,M-1
        DO I=2,N-1
    F(I,J)= (U(I+1,J)-U(I-1,J)+V(I,J+1)-V(I,J-1))/S(I,J)
        enddo
      enddo

      DO I=2,N-1
  F(I,1)= (U(I+1,1)-U(I-1,1)+(V(I,2)+V(I,1)))/S(I,1)
  F(I,M)= (U(I+1,M)-U(I-1,M)-(V(I,M)+V(I,M-1)))/S(I,M)
      end do
  CALL XBC(F,N,M)

endif

! end original

!DO I=2,N-1
!  F(I,1)= (U(I+1,1)-U(I-1,1)+(V(I,2)))/S(I,1)
!  F(I,M)= (U(I+1,M)-U(I-1,M)-(V(I,M-1)))/S(I,M)
!ENDDO




END SUBROUTINE


SUBROUTINE LAPLfirst_depth(P,F,A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP, num_of_bits, DP_Depth)
use implicit_functions_SP

implicit none
REAL(Kind=4) ::  P(N,M),F(N,M),A11(N,M),A12(N,M),A21(N,M),A22(N,M),  &
    &      B11(N,M),B22(N,M),P0(N,M),U(N,M),V(N,M),S(N,M), S_L(N, M)
INTEGER :: IP(N)
INTEGER :: N, M, num_of_bits, DP_Depth

REAL(Kind=4) ::  UTIL, VTIL
INTEGER :: I, J



! truncating
DO J=1+DP_Depth,M-DP_Depth
 DO I=1,N

  S_L(I, J)=S(I, J)
end do
end do

If (DP_Depth<=0) then

DO J=2,M-1
  DO I=2,N-1
    U(I,J)=P(I+1,J)-P(I-1,J)
    V(I,J)=P(I,J+1)-P(I,J-1)
  end do
end do
  
DO I=2,N-1
  U(I,1)=P(I+1,1)-P(I-1,1)
  U(I,M)=P(I+1,M)-P(I-1,M)
  V(I,1)=P(I,2)-P(IP(I),1)
  V(I,M)=P(IP(I),M)-P(I,M-1)
ENDDO  
CALL XBC(U,N,M)
CALL XBC(V,N,M)
else

      DO J=2+(DP_Depth-1),M-1-(DP_Depth-1)
        DO I=2,N-1
    U(I,J)=P(I+1,J)-P(I-1,J)
    V(I,J)=P(I,J+1)-P(I,J-1)
        enddo
      enddo




      DO J=2,DP_Depth
        DO I=2,N-1
    U(I,J)=P(I+1,J)-P(I-1,J)
    V(I,J)=P(I,J+1)-P(I,J-1)
        enddo
      enddo

      DO J=M+1-DP_Depth,M-1
        DO I=2,N-1
    U(I,J)=P(I+1,J)-P(I-1,J)
    V(I,J)=P(I,J+1)-P(I,J-1)
        enddo
      enddo

      DO I=2,N-1
  U(I,1)=P(I+1,1)-P(I-1,1)
  U(I,M)=P(I+1,M)-P(I-1,M)
  V(I,1)=P(I,2)-P(IP(I),1)
  V(I,M)=P(IP(I),M)-P(I,M-1)
      end do
CALL XBC(U,N,M)
CALL XBC(V,N,M)


endif




!DO J=1,M
!  DO I=1,N
!    UTIL=U(I,J)*A11(I,J)+A12(I,J)*V(I,J)+B11(I,J)*(P(I,J))  ! why only B11*P and not B11*(P-P0)
!    VTIL=V(I,J)*A21(I,J)+A22(I,J)*U(I,J)+B22(I,J)*(P(I,J))
!    U(I,J)=UTIL
!    V(I,J)=VTIL
!  ENDDO
!ENDDO

  !! rewritten to change precision at poles as desired    
      DO J=1+DP_Depth,M-DP_Depth
        DO I=1,N
    UTIL=U(I,J)*A11(I,J)+A12(I,J)*V(I,J)+B11(I,J)*(P0(I,J))  ! B11*(P-P0)
    VTIL=V(I,J)*A21(I,J)+A22(I,J)*U(I,J)+B22(I,J)*(P0(I,J))
    U(I,J)=UTIL
    V(I,J)=VTIL
        enddo
      enddo



      DO J=1,DP_Depth
        DO I=1,N
    UTIL=U(I,J)*A11(I,J)+A12(I,J)*V(I,J)+B11(I,J)*(P0(I,J))  ! B11*(P-P0)
    VTIL=V(I,J)*A21(I,J)+A22(I,J)*U(I,J)+B22(I,J)*(P0(I,J))
    U(I,J)=UTIL
    V(I,J)=VTIL
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=1,N
    UTIL=U(I,J)*A11(I,J)+A12(I,J)*V(I,J)+B11(I,J)*(P0(I,J))  ! B11*(P-P0)
    VTIL=V(I,J)*A21(I,J)+A22(I,J)*U(I,J)+B22(I,J)*(P0(I,J))
    U(I,J)=UTIL
    V(I,J)=VTIL
        enddo
      enddo
CALL XBC(U,N,M)
CALL XBC(V,N,M)

   !! end of rewrite




If (DP_Depth<=0) then

DO J=2,M-1
  DO I=2,N-1
    F(I,J)= (U(I+1,J)-U(I-1,J)+V(I,J+1)-V(I,J-1))/S_L(I,J)
  end do
end do    
 
DO I=2,N-1
  F(I,1)= (U(I+1,1)-U(I-1,1)+(V(I,2)+V(I,1)))/S_L(I,1)
  F(I,M)= (U(I+1,M)-U(I-1,M)-(V(I,M)+V(I,M-1)))/S_L(I,M)
ENDDO
CALL XBC(F,N,M)

else

      DO J=2+(DP_Depth-1),M-1-(DP_Depth-1)
        DO I=2,N-1
    F(I,J)= (U(I+1,J)-U(I-1,J)+V(I,J+1)-V(I,J-1))/S_L(I,J)
        enddo
      enddo


      DO J=2,DP_Depth
        DO I=2,N-1
    F(I,J)= (U(I+1,J)-U(I-1,J)+V(I,J+1)-V(I,J-1))/S(I,J)
        enddo
      enddo

      DO J=M+1-DP_Depth,M-1
        DO I=2,N-1
    F(I,J)= (U(I+1,J)-U(I-1,J)+V(I,J+1)-V(I,J-1))/S(I,J)
        enddo
      enddo

      DO I=2,N-1
  F(I,1)= (U(I+1,1)-U(I-1,1)+(V(I,2)+V(I,1)))/S(I,1)
  F(I,M)= (U(I+1,M)-U(I-1,M)-(V(I,M)+V(I,M-1)))/S(I,M)
      end do
  CALL XBC(F,N,M)

endif

! end original

!DO I=2,N-1
!  F(I,1)= (U(I+1,1)-U(I-1,1)+(V(I,2)))/S(I,1)
!  F(I,M)= (U(I+1,M)-U(I-1,M)-(V(I,M-1)))/S(I,M)
!ENDDO




END SUBROUTINE



SUBROUTINE LAPLfirst(P,F,A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP)
use implicit_functions_SP

implicit none
REAL(Kind=4) :: P(N,M),F(N,M),A11(N,M),A12(N,M),A21(N,M),A22(N,M),  &
    &      B11(N,M),B22(N,M),P0(N,M),U(N,M),V(N,M),S(N,M)
INTEGER :: IP(N)
INTEGER :: N, M

REAL(Kind=4) :: UTIL, VTIL
INTEGER :: I, J


DO J=2,M-1
  DO I=2,N-1
    U(I,J)=P(I+1,J)-P(I-1,J)
    V(I,J)=P(I,J+1)-P(I,J-1)
  end do
end do
  
DO I=2,N-1
  U(I,1)=P(I+1,1)-P(I-1,1)
  U(I,M)=P(I+1,M)-P(I-1,M)
  V(I,1)=P(I,2)-P(IP(I),1)
  V(I,M)=P(IP(I),M)-P(I,M-1)
ENDDO   

CALL XBC(U,N,M)
CALL XBC(V,N,M)



DO J=1,M
  DO I=1,N
    UTIL=U(I,J)*A11(I,J)+A12(I,J)*V(I,J)+B11(I,J)*(P(I,J)-P0(I,J))  ! why only B11*P and not B11*(P-P0)
    VTIL=V(I,J)*A21(I,J)+A22(I,J)*U(I,J)+B22(I,J)*(P(I,J)-P0(I,J))
    U(I,J)=UTIL
    V(I,J)=VTIL
  ENDDO
ENDDO

CALL XBC(U,N,M)
CALL XBC(V,N,M)



DO J=2,M-1
  DO I=2,N-1
    F(I,J)= (U(I+1,J)-U(I-1,J)+V(I,J+1)-V(I,J-1))/S(I,J)
  end do
end do    
 

DO I=2,N-1
  F(I,1)= (U(I+1,1)-U(I-1,1)+(V(I,2)+V(I,1)))/S(I,1)
  F(I,M)= (U(I+1,M)-U(I-1,M)-(V(I,M)+V(I,M-1)))/S(I,M)
ENDDO

CALL XBC(F,N,M)


END SUBROUTINE



SUBROUTINE VELPRD(U,V,F,G,PD,HX,HY,IP,N,M,A,B,EP,KMX)
use implicit_functions_SP

implicit none
INTEGER :: N, M, KMX
REAL(Kind=4) ::  U(N,M,0:1),V(N,M,0:1),F(N,M,0:1),G(N,M,0:1),         &
     &          PD(N,M),HX(N,M),HY(N,M)
INTEGER :: IP(N)
REAL(Kind=4) :: EP, A, B

REAL(Kind=4) :: UU(N,M),VV(N,M)
REAL(Kind=4) ::  CF, C1, C2, C1H, C2H, ALFA, BETA, ALF1, BET1, ALFM, BETM
REAL(Kind=4) :: AMP, UF(N),VF(N),SCR(N)
INTEGER :: IORT, I, J

 AMP=0.0d0
 IORT=2
 CF=rpe_05
 C1=A*CF
 C2=B*CF
 C1H=C1*rpe_05
 C2H=C2*rpe_05

!  COMPUTE V+R/PHI*DT FIELD FOR LAGRANGIAN ESTIMATES
 
DO J=1,M
  DO I=1,N
    F(I,J,1)=U(I,J,0)+CF*F(I,J,0)/max(PD(I,J),EP)
    G(I,J,1)=V(I,J,0)+CF*G(I,J,0)/max(PD(I,J),EP)
  enddo
end do

! COMPUTE U AND V TO FIRST ORDER

DO J=2,M-1
  DO I=2,N-1
    ALFA=U(I,J,0)/HX(I,J)*C1
    BETA=V(I,J,0)/HY(I,J)*C2
    U(I,J,1)=F(I,J,1)-max(rpe_0,ALFA)*(F(I,J,1)-F(I-1,J,1))   &
      &              -min(rpe_0,ALFA)*(F(I+1,J,1)-F(I,J,1))   &
      &              -max(rpe_0,BETA)*(F(I,J,1)-F(I,J-1,1))   &
      &              -min(rpe_0,BETA)*(F(I,J+1,1)-F(I,J,1)) 
    V(I,J,1)=G(I,J,1)-max(rpe_0,ALFA)*(G(I,J,1)-G(I-1,J,1))   &
      &              -min(rpe_0,ALFA)*(G(I+1,J,1)-G(I,J,1))   &
      &              -max(rpe_0,BETA)*(G(I,J,1)-G(I,J-1,1))   &
      &              -min(rpe_0,BETA)*(G(I,J+1,1)-G(I,J,1))
  end do
end do

DO I=2,N-1
  ALF1=U(I,1,0)/HX(I,1)*C1
  BET1=V(I,1,0)/HY(I,1)*C2
  U(I,1,1)=F(I,1,1)-max(rpe_0,ALF1)*(F(I,1,1)-F(I-1,1,1))       & 
     &                 -min(rpe_0,ALF1)*(F(I+1,1,1)-F(I,1,1))   &
     &                 -max(rpe_0,BET1)*(F(I,1,1)-F(IP(I),1,1)) &
     &                 -min(rpe_0,BET1)*(F(I,2,1)-F(I,1,1))
  V(I,1,1)=G(I,1,1)-max(rpe_0,ALF1)*(G(I,1,1)-G(I-1,1,1))   &
     &                 -min(rpe_0,ALF1)*(G(I+1,1,1)-G(I,1,1))   & 
     &                 -max(rpe_0,BET1)*(G(I,1,1)+G(IP(I),1,1)) &
     &                 -min(rpe_0,BET1)*(G(I,2,1)-G(I,1,1))

  ALFM=U(I,M,0)/HX(I,M)*C1
  BETM=V(I,M,0)/HY(I,M)*C2
  U(I,M,1)=F(I,M,1)-max(rpe_0,ALFM)*(F(I,M,1)-F(I-1,M,1))        &
     &                 -min(rpe_0,ALFM)*(F(I+1,M,1)-F(I,M,1))    &
     &                 -max(rpe_0,BETM)*(F(I,M,1)-F(I,M-1,1))    &
     &                 -min(rpe_0,BETM)*(F(IP(I),M,1)-F(I,M,1)) 
  V(I,M,1)=G(I,M,1)-max(rpe_0,ALFM)*(G(I,M,1)-G(I-1,M,1))        &
     &                 -min(rpe_0,ALFM)*(G(I+1,M,1)-G(I,M,1))    &
     &                 -max(rpe_0,BETM)*(G(I,M,1)-G(I,M-1,1))    &
     &                 +min(rpe_0,BETM)*(G(IP(I),M,1)+G(I,M,1))
end do

CALL XBC(U(:,:,1),N,M)
CALL XBC(V(:,:,1),N,M)

IF(IORT.EQ.2) THEN
 ! COMPUTE U AND V TO SEMI-SECOND ORDER

  DO J=1,M
    DO I=1,N
      UU(I,J)=rpe_05*(U(I,J,0)+U(I,J,1))
      VV(I,J)=rpe_05*(V(I,J,0)+V(I,J,1))
    ENDDO
  ENDDO

  DO J=2,M-1
    DO I=2,N-1
      ALFA=UU(I,J)/HX(I,J)*C1
      BETA=VV(I,J)/HY(I,J)*C2
      U(I,J,1)=F(I,J,1)-max(rpe_0,ALFA)*(F(I,J,1)-F(I-1,J,1))     &
       &                 -min(rpe_0,ALFA)*(F(I+1,J,1)-F(I,J,1))   &
       &                 -max(rpe_0,BETA)*(F(I,J,1)-F(I,J-1,1))   &
       &                 -min(rpe_0,BETA)*(F(I,J+1,1)-F(I,J,1))
      V(I,J,1)=G(I,J,1)-max(rpe_0,ALFA)*(G(I,J,1)-G(I-1,J,1))     &
       &                 -min(rpe_0,ALFA)*(G(I+1,J,1)-G(I,J,1))   &
       &                 -max(rpe_0,BETA)*(G(I,J,1)-G(I,J-1,1))   &
       &                 -min(rpe_0,BETA)*(G(I,J+1,1)-G(I,J,1))
    end do
  end do

  DO I=2,N-1
    ALF1=UU(I,1)/HX(I,1)*C1
    BET1=VV(I,1)/HY(I,1)*C2
    U(I,1,1)=F(I,1,1)-max(rpe_0,ALF1)*(F(I,1,1)-F(I-1,1,1))     &
     &                 -min(rpe_0,ALF1)*(F(I+1,1,1)-F(I,1,1))   &
     &                 -max(rpe_0,BET1)*(F(I,1,1)-F(IP(I),1,1)) &
     &                 -min(rpe_0,BET1)*(F(I,2,1)-F(I,1,1)) 
    V(I,1,1)=G(I,1,1)-max(rpe_0,ALF1)*(G(I,1,1)-G(I-1,1,1))     &
      &                -min(rpe_0,ALF1)*(G(I+1,1,1)-G(I,1,1))   & 
      &                -max(rpe_0,BET1)*(G(I,1,1)+G(IP(I),1,1)) &
      &                -min(rpe_0,BET1)*(G(I,2,1)-G(I,1,1))
    ALFM=UU(I,M)/HX(I,M)*C1
    BETM=VV(I,M)/HY(I,M)*C2

    U(I,M,1)=F(I,M,1)-max(rpe_0,ALFM)*(F(I,M,1)-F(I-1,M,1))     &
     &                 -min(rpe_0,ALFM)*(F(I+1,M,1)-F(I,M,1))   &
     &                 -max(rpe_0,BETM)*(F(I,M,1)-F(I,M-1,1))   &
     &                 -min(rpe_0,BETM)*(F(IP(I),M,1)-F(I,M,1))
    V(I,M,1)=G(I,M,1)-max(rpe_0,ALFM)*(G(I,M,1)-G(I-1,M,1))     &
     &                 -min(rpe_0,ALFM)*(G(I+1,M,1)-G(I,M,1))   &
     &                 -max(rpe_0,BETM)*(G(I,M,1)-G(I,M-1,1))   &
     &                 +min(rpe_0,BETM)*(G(IP(I),M,1)+G(I,M,1))
  end do

  CALL XBC(U(:,:,1),N,M)
  CALL XBC(V(:,:,1),N,M)
ENDIF


     IF(KMX.GT.0) THEN
      DO J=1,M
       AMP=(1.-HX(1,J)/HY(1,J))**2
       DO I=1,N
        UF(I)=U(I,J,1)
        VF(I)=V(I,J,1)
       ENDDO
       CALL FILTRX(UF,SCR,AMP,N,KMX)
       CALL FILTRX(VF,SCR,AMP,N,KMX)
       DO I=1,N
        U(I,J,1)=UF(I)
        V(I,J,1)=VF(I)
       ENDDO
      enddo
     ENDIF


END SUBROUTINE


SUBROUTINE FILTRX(X,Y,AMP,N,KMX)

implicit none
REAL(Kind=4) ::  X(N),Y(N), AMP
INTEGER :: N, KMX, I, J,k, IM, IP

  DO K=1,KMX

    DO I=1,N
     IP=I+K

     IF(IP.GT.N) IP=IP-N+2
     IM=I-K
     IF(IM.LT.1) IM=IM+N-2
     Y(I)=X(I)+.25*AMP*(X(IP)-2.*X(I)+X(IM))
    ENDDO

    DO I=1,N
     X(I)=Y(I)
    ENDDO

  ENDDO   

END subroutine


SUBROUTINE FILTRQ(FLD,N,M)

implicit none

REAL(Kind=4) ::  FLD(N,M)
integer :: N, M, I


DO I=1,N
  FLD(I,1)   = FLD(I,1)  *(1.-0.5)
  FLD(I,2)   = FLD(I,2)  *(1.-0.25)
  FLD(I,3)   = FLD(I,3)  *(1.-0.125)
  FLD(I,M)   = FLD(I,M)  *(1.-0.5)
  FLD(I,M-1) = FLD(I,M-1)*(1.-0.25)
  FLD(I,M-2) = FLD(I,M-2)*(1.-0.125)
ENDDO

END subroutine








SUBROUTINE MPDATT(U1,U2,X,H,N,M,IORDs,ISOR,NONOS,IDIV,IBC, IP,liner, X_T, codes)

use implicit_functions_SP
implicit none

  INTEGER, PARAMETER  :: LW=1,MP=1-LW

  INTEGER :: N, M

  REAL(Kind=4) :: U1(N,M),U2(N,M+1),X(N,M),H(N,M), X_T(N,M)
  INTEGER :: IP(N)
  INTEGER  :: IORD, IORDs, ISOR, NONOS, IDIV, IBC 
      

  REAL(Kind=4) :: V1(N,M),V2(N,M+1),F1(N,M),F2(N,M+1)  &
      &      ,CP(N,M),CN(N,M)   
  REAL(Kind=4) :: MX(N,M),MN(N,M)
  REAL(Kind=4) :: EP, C1, C2, V1D, V2D, V2D1, V2DN
  INTEGER :: N1, N2, I, J, K, liner
  LOGICAL :: codes

  N1=N
  N2=M

  EP= 1.E-10

  IORD=IORDs
  IF(ISOR.EQ.3) IORD=MAX(IORD,3)
  if(LINER.EQ.1) IORD=1
  IF(LW.EQ.1) IORD=MIN(IORD,2)
  !write(*,*) IORD
      !! take predicted advective velocities
  DO J=1,N2
    DO I=1,N1
     V1(I,J)=U1(I,J)
    end do
  end do
  
  DO J=1,N2+1
    DO I=1,N1
     V2(I,J)=U2(I,J)
    end do
  end do

  
  IF(NONOS.EQ.1) THEN           
    DO J=2,N2-1
      DO I=2,N1-1
        MX(I,J)=max(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1))
        MN(I,J)=min(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1))
      end do
    end do
    
    DO I=2,N1-1
      MX(I,1)=max(X(I-1,1),X(I,1),X(I+1,1),IBC*X(IP(I),1),X(I,2))
      MN(I,1)=min(X(I-1,1),X(I,1),X(I+1,1),IBC*X(IP(I),1),X(I,2))
      MX(I,N2)=max(X(I-1,N2),X(I,N2),X(I+1,N2),X(I,N2-1),                &
         &                         IBC*X(IP(I),N2))
      MN(I,N2)=min(X(I-1,N2),X(I,N2),X(I+1,N2),X(I,N2-1),                &
         &                         IBC*X(IP(I),N2))
    end do  
      
    CALL XBC(MX,N1,N2)
    CALL XBC(MN,N1,N2)
  ENDIF

  
  C1=rpe_1
  C2=rpe_0
  
  DO K=1,IORD   !! k-th order correction
   
   ! COMPUTE DONOR-CELL FLUXES
    DO J=1,N2
      DO I=2,N1-1
        F1(I,J)=DONOR(C1*X(I-1,J)+C2,C1*X(I,J)+C2,V1(I,J))
      end do
    end do
    
    DO J=2,N2
      DO I=2,N1-1
        F2(I,J)=DONOR(C1*X(I,J-1)+C2,C1*X(I,J)+C2,V2(I,J))
      end do
    end do
      
    DO I=2,N1-1
      F2(I,N2+1)=DONOR(C1*X(I,N2)+C2,C1*IBC*X(IP(I),N2)+C2,V2(I,N2+1))
      F2(I,1)=DONOR(C1*IBC*X(IP(I),1)+C2,C1*X(I,1)+C2,V2(I,1))
    end do
    
    CALL XBC(F1,N1,N2)
    CALL XBC(F2,N1,N2+1)
   ! COMPUTE NEW UPWIND-SOLUTION
    
    DO J=1,N2
      DO I=2,N1-1
        X(I,J)=X(I,J)-(F1(I+1,J)-F1(I,J)+F2(I,J+1)-F2(I,J))/H(I,J)
      end do
    end do
    
    
    IF (codes) then
      DO J=1,N2
        DO I=2,N1-1
          X_T(I,J)=X_T(I,J) -(F1(I+1,J)-F1(I,J)+F2(I,J+1)-F2(I,J))/H(I,J)
        end do
      end do
      
      CALL XBC(X_T,N1,N2)
    end if
    
    
    CALL XBC(X,N1,N2)

    IF(.not. (K.EQ.IORD)) then  ! GO TO 6
    
      C1=FLOAT(MP)
      C2=FLOAT(LW)
   ! CONVERT VELOCITIES TO LOCAL STORAGE
      DO J=1,N2
        DO I=1,N1
          F1(I,J)=V1(I,J)
        end do
      end do
      
      DO J=1,N2+1
        DO I=1,N1
          F2(I,J)=V2(I,J) 
        end do
      end do 
   ! CALCULATE PSEUDO VELOCITIES      
   ! COMPUTE FIRST DIRECTION    
      
      DO J=2,N2-1
        DO I=2,N1-1
          V1(I,J)=VDYF_D(X(I-1,J),X(I,J),F1(I,J),rpe_05*(H(I-1,J)+H(I,J)),MP, EP)     &
               & +VCORR_D(F1(I,J), F2(I-1,J)+F2(I-1,J+1)+F2(I,J+1)+F2(I,J),                 &
               &   X(I-1,J-1),X(I,J-1),X(I-1,J+1),X(I,J+1),                   &
               &       rpe_05*(H(I-1,J)+H(I,J)),MP, EP)
        end do
      end do
   ! COMPUTE B.C IN Y DIRECTION
      DO I=2,N1-1
        V1(I,1)=VDYF_D(X(I-1,1),X(I,1),F1(I,1),rpe_05*(H(I-1,1)+H(I,1)),MP, EP)         &
           & +VCORR_D(F1(I,1), F2(I-1,1)+F2(I-1,2)+F2(I,2)+F2(I,1),                 &
           & IBC*X(IP(I-1),1),IBC*X(IP(I),1),X(I-1,2),X(I,2),                       &
           &           rpe_05*(H(I-1,1)+H(I,1)),MP, EP)
        V1(I,N2)=VDYF_D(X(I-1,N2),X(I,N2),F1(I,N2),rpe_05*(H(I-1,N2)+H(I,N2)),MP, EP)     &
           & +VCORR_D(F1(I,N2), F2(I-1,N2)+F2(I-1,N2+1)+F2(I,N2+1)+F2(I,N2),              &
           &   X(I-1,N2-1),X(I,N2-1),IBC*X(IP(I-1),N2),IBC*X(IP(I),N2),              &
           &           rpe_05*(H(I-1,N2)+H(I,N2)),MP, EP)        
      end do
      
      IF(IDIV.EQ.1) THEN
 
 
    ! COMPUTE FLOW-DIVERGENCE CORRECTION
        DO J=1,N2
          DO I=2,N1-1
            V1D=-VDIV1(F1(I-1,J),F1(I,J),F1(I+1,J),rpe_05*(H(I-1,J)+H(I,J)))              &
              & -VDIV2(F1(I,J),F2(I-1,J+1),F2(I,J+1),F2(I-1,J),F2(I,J),               &
              &        rpe_05*(H(I-1,J)+H(I,J)))                            
            V1(I,J)=V1(I,J)+LW*(PP(V1D)*X(I-1,J)-PN(V1D)*X(I,J))+MP*V1D
          end do
        end do
        
      ENDIF   

      
    ! COMPUTE SECOND DIRECTION
      DO J=2,N2
        DO I=2,N1-1
          V2(I,J)=VDYF_D(X(I,J-1),X(I,J),F2(I,J),rpe_05*(H(I,J-1)+H(I,J)),MP, EP)         &
             & +VCORR_D(F2(I,J), F1(I,J-1)+F1(I,J)+F1(I+1,J)+F1(I+1,J-1),             &
             &          X(I-1,J-1),X(I-1,J),X(I+1,J-1),X(I+1,J),                      &
             &         rpe_05*(H(I,J-1)+H(I,J)),MP, EP)
        end do
      end do

    ! COMPUTE B.C IN Y-DIRECTION
      DO I=2,N1-1
        V2(I,1)=VDYF_D(IBC*X(IP(I),1),X(I,1),F2(I,1),rpe_05*(H(IP(I),1)+H(I,1)),MP, EP)   &
          & +VCORR_D(F2(I,1), F1(IP(I),1)+F1(I,1)+F1(I+1,1)+F1(IP(I+1),1),            &
          &    IBC*X(IP(I-1),1),X(I-1,1),X(I+1,1),IBC*X(IP(I+1),1),                   &
          &            rpe_05*(H(IP(I),1)+H(I,1)),MP, EP)
        V2(I,N2+1)=VDYF_D(X(I,N2),IBC*X(IP(I),N2),F2(I,N2+1),                         &
          &                   rpe_05*(H(I,N2)+H(IP(I),N2)),MP, EP)                        &
          & +VCORR_D(F2(I,N2+1),F1(I,N2)+F1(IP(I),N2)+F1(IP(I+1),N2)+F1(I+1,N2),      &
          & X(I-1,N2),IBC*X(IP(I-1),N2),X(I+1,N2),IBC*X(IP(I+1),N2),                  &
          &                   rpe_05*(H(I,N2)+H(IP(I),N2)),MP, EP)
      end do
      
      IF(IDIV.EQ.1) THEN
        DO J=2,N2
          DO I=2,N1-1
            V2D=-VDIV1(F2(I,J-1),F2(I,J),F2(I,J+1),rpe_05*(H(I,J-1)+H(I,J)))              &
              & -VDIV2(F2(I,J),F1(I+1,J-1),F1(I+1,J),F1(I,J-1),F1(I,J),               &
              &   rpe_05*(H(I,J-1)+H(I,J)))
            V2(I,J)=V2(I,J)+LW*(PP(V2D)*X(I,J-1)-PN(V2D)*X(I,J))+MP*V2D
          end do
        end do
        
        DO I=2,N1-1
          V2D1=-VDIV1(-F2(IP(I),1),F2(I,1),F2(I,2),rpe_05*(H(IP(I),1)+H(I,1)))            &
            & -VDIV2( F2(I,1),F1(IP(I+1),1),F1(I+1,1),F1(IP(I),1),F1(I,1),            &
            &      rpe_05*(H(IP(I),1)+H(I,1))) 
          V2(I,1)=V2(I,1)                                                             &
            & +LW*(PP(V2D1)*IBC*X(IP(I),1)-PN(V2D1)*X(I,1))+MP*V2D1
          V2DN=-VDIV1(F2(I,N2),F2(I,N2+1),-F2(IP(I),N2+1)                             &
            &              ,rpe_05*(H(I,N2)+H(IP(I),N2)))                                 &
            &  -VDIV2(F2(I,N2+1),F1(I+1,N2),F1(IP(I+1),N2),                           &
            &     F1(I,N2),F1(IP(I),N2),rpe_05*(H(I,N2)+H(IP(I),N2))) 
          V2(I,N2+1)=V2(I,N2+1)                                                       &
            &  +LW*(PP(V2DN)*X(I,N2)-PN(V2DN)*IBC*X(IP(I),N2))+MP*V2DN
        end do
       
      ENDIF
      
!
    ! THIRD ORDER CORRECTION
      
      IF(ISOR.EQ.3) THEN
      
    ! FIRST DIRECTION
        DO J=1,N2
          DO I=3,N1-1
            V1(I,J)=V1(I,J)     +VCOR31_D(F1(I,J),                                    &
              &  X(I-2,J),X(I-1,J),X(I,J),X(I+1,J),rpe_05*(H(I-1,J)+H(I,J)),MP, EP)
          end do
   
          V1(2,J)=V1(2,J)     +VCOR31_D(F1(2,J),                                      &
            &  X(N1-2,J),X(1,J),X(2,J),X(3,J),rpe_05*(H(1,J)+H(2,J)),MP, EP)
        end do
        
        DO J=2,N2-1
          DO I=2,N1-1               
            V1(I,J)=V1(I,J)                                                           &
              &  +VCOR32_D(F1(I,J),F2(I-1,J)+F2(I-1,J+1)+F2(I,J+1)+F2(I,J),           &
              &   X(I,J-1),X(I-1,J+1),X(I-1,J-1),X(I,J+1),                            &
              &  rpe_05*(H(I-1,J)+H(I,J)),MP, EP)
          end do
        end do
    ! C B.C. FOLLOW
        DO I=2,N1-1
          V1(I,1)=V1(I,1)                                                             &
            & +VCOR32_D(F1(I,1),F2(I-1,1)+F2(I-1,2)+F2(I,2)+F2(I,1),                  &
            & IBC*X(IP(I),1),X(I-1,2),X(I,2),IBC*X(IP(I-1),1),                        &
            & rpe_05*(H(I-1,1)+H(I,1)),MP, EP)
          V1(I,N2)=V1(I,N2)                                                           &
            & +VCOR32_D(F1(I,N2),F2(I-1,N2)+F2(I-1,N2+1)+F2(I,N2+1)+F2(I,N2),         &
            & X(I,N2-1),IBC*X(IP(I-1),N2),IBC*X(IP(I),N2),X(I-1,N2-1),                &
            & rpe_05*(H(I-1,N2)+H(I,N2)),MP, EP)
        end do
      
        DO J=1,N2
          V1(1,J)=V1(N1-1,J)
          V1(N1,J)=V1(2,J)
        end do

    !  SECOND DIRECTION
        
        DO J=3,N2-1
          DO I=2,N1-1
            V2(I,J)=V2(I,J)     +VCOR31_D(F2(I,J),                                    &
              &  X(I,J-2),X(I,J-1),X(I,J),X(I,J+1),rpe_05*(H(I,J-1)+H(I,J)),MP, EP)
          end do
        end do
        
        DO I=2,N1-1
          V2(I,1)=V2(I,1)     +VCOR31_D(F2(I,1),IBC*X(IP(I),2),                       &
            &  IBC*X(IP(I),1),X(I,1),X(I,2),rpe_05*(H(IP(I),1)+H(I,1)),MP, EP)  
          V2(I,2)=V2(I,2)     +VCOR31_D(F2(I,2),                                      &
            &  IBC*X(IP(I),1),X(I,1),X(I,2),X(I,3),rpe_05*(H(I,1)+H(I,2)),MP, EP)   
          V2(I,N2)=V2(I,N2)     +VCOR31_D(F2(I,N2),X(I,N2-2),                         &
            &  X(I,N2-1),X(I,N2),IBC*X(IP(I),N2),rpe_05*(H(I,N2-1)+H(I,N2)),MP, EP)  
          V2(I,N2+1)=V2(I,N2+1) +VCOR31_D(F2(I,N2+1), X(I,N2-1),X(I,N2),              &
            &  IBC*X(IP(I),N2),IBC*X(IP(I),N2-1),rpe_05*(H(I,N2)+H(IP(I),N2)),MP, EP)  
        end do
        
        DO J=2,N2
          DO I=2,N1-1
            V2(I,J)=V2(I,J)                                                           &
              & +VCOR32_D(F2(I,J),F1(I,J-1)+F1(I+1,J-1)+F1(I+1,J)+F1(I,J),            &
              &  X(I+1,J-1),X(I-1,J),X(I-1,J-1),X(I+1,J),                             &
              &         rpe_05*(H(I,J-1)+H(I,J)),MP, EP)
          end do
        end do
        
        
      ENDIF  !! end third order correction

     ! CALL B.C IN X DIRECTION
      CALL XBC(V1,N1,N2)
      CALL XBC(V2,N1,N2+1)

      
      
      IF(.not. (NONOS.EQ.0)) then  !! if non-osc is not turned off
     ! NON-OSCILLATORY OPTION
        DO J=2,N2-1
          DO I=2,N1-1
            MX(I,J)=max(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1),MX(I,J))
            MN(I,J)=min(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1),MN(I,J))
          end do
        end do
      
        DO I=2,N1-1
          MX(I,1)=max(X(I-1,1),X(I,1),X(I+1,1),IBC*X(IP(I),1),                      &
             &                X(I,2),MX(I,1))
          MN(I,1)=min(X(I-1,1),X(I,1),X(I+1,1),IBC*X(IP(I),1),                      &
             &                X(I,2),MN(I,1))
          MX(I,N2)=max(X(I-1,N2),X(I,N2),X(I+1,N2),X(I,N2-1),                       &
             &                       IBC*X(IP(I),N2),MX(I,N2))
          MN(I,N2)=min(X(I-1,N2),X(I,N2),X(I+1,N2),X(I,N2-1),                       &
             &                      IBC*X(IP(I),N2),MN(I,N2))
        end do
        
        CALL XBC(MX,N1,N2)
        CALL XBC(MN,N1,N2)

        DO J=1,N2
          DO I=2,N1-1
            F1(I,J)=DONOR(C1*X(I-1,J)+C2,C1*X(I,J)+C2,V1(I,J))
          end do
        end do
      
        DO J=2,N2
          DO I=2,N1-1
            F2(I,J)=DONOR(C1*X(I,J-1)+C2,C1*X(I,J)+C2,V2(I,J))
          end do
        end do
      
        DO I=2,N1-1    
          F2(I,N2+1)=DONOR(C1*X(I,N2)+C2,C1*IBC*X(IP(I),N2)+C2,V2(I,N2+1))
          F2(I,1)=DONOR(C1*IBC*X(IP(I),1)+C2,C1*X(I,1)+C2,V2(I,1))
        end do
      
        CALL XBC(F1,N1,N2)
        CALL XBC(F2,N1,N2+1)

        DO J=1,N2
          DO I=2,N1-1
            CP(I,J)=(MX(I,J)-X(I,J))*H(I,J)/                                     &
               & (PN(F1(I+1,J))+PP(F1(I,J))+PN(F2(I,J+1))+PP(F2(I,J))+EP)
            CN(I,J)=(X(I,J)-MN(I,J))*H(I,J)/                                     &                        
               & (PP(F1(I+1,J))+PN(F1(I,J))+PP(F2(I,J+1))+PN(F2(I,J))+EP)
          end do
        end do
        
        CALL XBC(CP,N1,N2)
        CALL XBC(CN,N1,N2)

        IF(LW.EQ.0) THEN
          
          DO J=2,N2
            DO I=2,N1-1
              V1(I,J)=PP(V1(I,J))*                                              &
                & ( min(rpe_1,CP(I,J),CN(I-1,J))*PP(SIGN(rpe_1, X(I-1,J)))          &
                &  +min(rpe_1,CP(I-1,J),CN(I,J))*PP(SIGN(rpe_1,-X(I-1,J))) )         &
                & -PN(V1(I,J))*                                                 &
                & ( min(rpe_1,CP(I-1,J),CN(I,J))*PP(SIGN(rpe_1, X(I ,J )))          &
                &  +min(rpe_1,CP(I,J),CN(I-1,J))*PP(SIGN(rpe_1,-X(I ,J ))) )     
              V2(I,J)=PP(V2(I,J))*                                              &
                & ( min(rpe_1,CP(I,J),CN(I,J-1))*PP(SIGN(rpe_1, X(I,J-1)))          &
                &  +min(rpe_1,CP(I,J-1),CN(I,J))*PP(SIGN(rpe_1,-X(I,J-1))) )         &
                & -PN(V2(I,J))*                                                 &
                & ( min(rpe_1,CP(I,J-1),CN(I,J))*PP(SIGN(rpe_1, X(I ,J )))          &
                &  +min(rpe_1,CP(I,J),CN(I,J-1))*PP(SIGN(rpe_1,-X(I ,J ))) )     
             end do
           end do
       ! B.C. FOLLOW
       
          DO I=2,N1-1
            V2(I,1)=PP(V2(I,1))*                                              &
              & ( min(rpe_1,CP(I,1),CN(IP(I),1))*PP(SIGN(rpe_1, IBC*X(IP(I),1)))    &
              & +min(rpe_1,CP(IP(I),1),CN(I,1))*PP(SIGN(rpe_1,-IBC*X(IP(I),1))) )   &
              & -PN(V2(I,1))*                                                 &
              & ( min(rpe_1,CP(IP(I),1),CN(I,1))*PP(SIGN(rpe_1, X(I ,1 )))          &
              & +min(rpe_1,CP(I,1),CN(IP(I),1))*PP(SIGN(rpe_1,-X(I ,1 ))) )       
            V2(I,N2+1)=PP(V2(I,N2+1))*                                          &
              & ( min(rpe_1,CP(IP(I),N2),CN(I,N2))*PP(SIGN(rpe_1, X(I,N2)))         & 
              & +min(rpe_1,CP(I,N2),CN(IP(I),N2))*PP(SIGN(rpe_1,-X(I,N2))) )        &
              & -PN(V2(I,N2+1))*                                                &
              & ( min(rpe_1,CP(I,N2),CN(IP(I),N2))*PP(SIGN(rpe_1, IBC*X(IP(I),N2))) &
              & +min(rpe_1,CP(IP(I),N2),CN(I,N2))*PP(SIGN(rpe_1,-IBC*X(IP(I),N2)))) 
            V1(I,1)=PP(V1(I,1))*                                                & 
              & ( min(rpe_1,CP(I,1),CN(I-1,1))*PP(SIGN(rpe_1, X(I-1,1)))            &
              & +min(rpe_1,CP(I-1,1),CN(I,1))*PP(SIGN(rpe_1,-X(I-1,1))) )           &
              & -PN(V1(I,1))*                                                   &
              & ( min(rpe_1,CP(I-1,1),CN(I,1))*PP(SIGN(rpe_1, X(I ,1 )))            &
              & +min(rpe_1,CP(I,1),CN(I-1,1))*PP(SIGN(rpe_1,-X(I ,1 ))) )
          end do
          
        ELSE
     
          DO J=2,N2
            DO I=2,N1-1
              V1(I,J)=PP(V1(I,J))*min(rpe_1,CP(I,J),CN(I-1,J))                   &
                & -PN(V1(I,J))*min(rpe_1,CP(I-1,J),CN(I,J))
              V2(I,J)=PP(V2(I,J))*min(rpe_1,CP(I,J),CN(I,J-1))                   & 
                & -PN(V2(I,J))*min(rpe_1,CP(I,J-1),CN(I,J))
            end do
          end do
      ! B.C. FOLLOW
      
          DO I=2,N1-1
            V2(I,1)=PP(V2(I,1))*min(rpe_1,CP(I,1),CN(IP(I),1))                   &
              & -PN(V2(I,1))*min(rpe_1,CP(IP(I),1),CN(I,1))
            V2(I,N2+1)=PP(V2(I,N2+1))*min(rpe_1,CP(IP(I),N2),CN(I,N2))           &
              & -PN(V2(I,N2+1))*min(rpe_1,CP(I,N2),CN(IP(I),N2))
            V1(I,1)=PP(V1(I,1))*min(rpe_1,CP(I,1),CN(I-1,1))                     &
              & -PN(V1(I,1))*min(rpe_1,CP(I-1,1),CN(I,1))
          end do

        ENDIF
      
        CALL XBC(V1,N1,N2)
        CALL XBC(V2,N1,N2+1)
!
!         END OF NONOSCILLATORY OPTION
!
      end if
      
    end if ! IF last iteration exit loop
  
  end do
  
END SUBROUTINE   

  


SUBROUTINE XBC(X,N,M)
use implicit_functions_SP

implicit none

INTEGER :: N, M
REAL(Kind=4) :: X(N,M)

INTEGER :: J

DO J=1,M
  X(1,J)=X(N-1,J)
  X(N,J)=X(2,J)
end do

END subroutine

SUBROUTINE XBC_23(X,N,M)
use implicit_functions_SP

implicit none

INTEGER :: N, M
REAL(Kind=4) :: X(N,M)

INTEGER :: J

DO J=1,M
  X(1,J)=X(N-1,J)
  X(N,J)=X(2,J)
end do

END subroutine


SUBROUTINE DIAGNOS(U,V,PD,PT,HX,HY,IP,S,TIME,DX,DY,DT, SUM0,SUM1, &
         & KT,N,M, IFLG, NITER,NITSM,ICOUNT,ERROR, sum_time, sum_lp_time)

use implicit_functions_SP

implicit none

INTEGER :: N, M

REAL(Kind=4) ::  U(N,M),V(N,M),PD(N,M),PT(N,M),HX(N,M),HY(N,M),  &
     &        S(N,M)
REAL(Kind=4) ::  TIME,DX,DY,DT, SUM0,SUM1,ERROR
INTEGER  :: IP(N)
INTEGER  :: KT, IFLG, NITER,NITSM,ICOUNT
REAL(Kind=4) :: avg_time, sum_time, avg_lp_time, sum_lp_time, NITAV
REAL(Kind=4) ::  GC1, GC2, COUR1, COUR2, PDMX,PDMN,PDAV, SUMER, DLI
INTEGER ::  I, J



GC1=DT/DX
GC2=DT/DY

IF(IFLG.EQ.0) THEN

  KT=0
  NITER=0
  ERROR=0.0d0
  SUM0=0.0d0

  DO J=1,M
    DO I=2,N-1
      SUM0=SUM0+PD(I,J)*S(I,J)
    end do
  end do       

ENDIF

TIME=KT*DT/3600.0
      PRINT 299
  299 FORMAT(1X,1H )
      PRINT *, 'HTIME= ', TIME,'KT= ',KT
  300 FORMAT(14X,5HTIME=,F7.2,7H HOURS;,5H  KT=,I5)

  ! CHECK COURANT NUMBERS

 COUR1=0.0d0
 COUR2=0.0d0

!GC1_23= DT_23/(2.0d0*ACOS(-1.0d0)/FLOAT(N-2))                                ! DT/DX
!GC1 = GC1_23
!GC2_23= DT_23/(ACOS(-1.0d0)/FLOAT(M))                                        ! DT/DY
!GC2 = GC2_23


DO J=1,M
  DO I=2,N
    DLI=SQRT((GC1/HX(I,J))**2+(GC2/HY(I,J))**2)
    COUR1=max(COUR1,DLI*SQRT(abs(PT(I,J))))
    COUR2=max(COUR2, GC1*ABS(U(I,J)/HX(I,J))                          &
      &                  +GC2*ABS(V(I,J)/HY(I,J)) )
   ! COUR2=max(COUR2, max(DT*ABS(U(I,J)/ (HX(I,J)/FLOAT(N-2)) ),                           &
   !   &             DT*ABS(V(I,J)/ (HY(I,J)/FLOAT(M)  ) ) ) &
   !   &       )
  end do
end do

      PRINT *,'COUR1,COUR2: ', COUR1,COUR2
  301 FORMAT(4X,'COUR1,COUR2:',2E11.4)

SUM1=0.0d0
PDMX=-1.E30
PDMN= 1.E30
PDAV=0.0d0

DO J=1,M
  DO I=2,N-1
    PDMX=max(PDMX,PD(I,J))
    PDMN=min(PDMN,PD(I,J))
    PDAV=PDAV+PD(I,J)
    SUM1=SUM1+PD(I,J)*S(I,J)
  end do
end do

PDAV=PDAV/FLOAT(M*(N-2))
      write(*,*) 'SUM1, SUM0, MAX(SUM0,PDMN)', SUM1, SUM0, MAX(SUM0,PDMN)
SUMER=(SUM1-SUM0)/max(SUM0,PDMN)

      PRINT *, 'PDMX,PDMN,PDAV: ', PDMX,PDMN,PDAV,' SUMER: ',SUMER

avg_time=sum_time/MAX(ICOUNT,1)
avg_lp_time=sum_lp_time/MAX(ICOUNT,1)
NITAV=float(NITSM)/MAX(ICOUNT,1)

      write(*,*) 'ERROR:', ERROR,'NITER, NITAV (GCR): ',NITER,NITAV
!  303 FORMAT(4X,'ERROR:',E11.4,1X,'NITER, NITAV (GCR ITERATIONS):',2I4)
      write(*,*) 'Computing time per implicit solve, low precision:', avg_time, avg_lp_time

      write(*,*) 'U 1 zonal band', SUM(PD(2:N-1,1)*U(2:N-1,1))/(N-2), 'U M Zonal Band',SUM(PD(2:N-1,M)*U(2:N-1,M))/(N-2)
      write(*,*) 'U 4 zonal band', SUM(PD(2:N-1,4)*U(2:N-1,4))/(N-2), 'U M-35 Zonal Band',SUM(PD(2:N-1,M-3)*U(2:N-1,M-3))/(N-2)
END subroutine



subroutine  init_perf_markers(H_rpe,U_rpe,V_rpe, TIME_rpe, codesignQ, codesignD, &
                    &   IRHW, X_rpe, Y_rpe, N, M, bits, ID_prec, EXP_NAME)
use implicit_functions_SP

implicit none

REAL(Kind=4) ::  H_rpe(N, M), U_rpe(N, M), V_rpe(N, M), TIME_rpe, X_rpe(N), Y_rpe(M)

integer :: N, M, IRHW, bits, ID_prec 
REAL(Kind=4) :: H(N, M), U(N, M), V(N, M), X(N), Y(M), TIME
logical :: codesignQ,codesignD, itsopen

 character(len=150) path, file_name, experiment, simtime, codesQ, bit_count, Precond, EXP_NAME, codesD
INTEGER I,J


H(:, :)=H_rpe(:, :)
U(:, :)=U_rpe(:, :)
V(:, :)=V_rpe(:, :)
TIME=TIME_rpe

path ='../'//trim(adjustl(EXP_NAME))//'/Precon'
write(experiment,*) IRHW
write(Precond,*) ID_prec
write(codesQ,*) codesignQ
write(codesD,*) codesignD
write(bit_count,*) bits

inquire(unit=324, opened=itsopen)
If(.not. itsopen) then 
  file_name = trim(adjustl(Precond))//'_EXP_'//trim(adjustl(experiment)) &
      &  //'_codes_'//trim(adjustl(codesQ))//trim(adjustl(codesD))//'.txt'
  Open(unit=324, file=trim(path)//trim(file_name), status='replace', form = 'formatted')
 write(unit=324,fmt='(I5)') bits
else

 write(unit=324,fmt=*)
 write(unit=324,fmt='(I5)') bits
endif

end subroutine init_perf_markers

subroutine  close_perf_markers()
use implicit_functions_SP

 close(Unit=324)

end subroutine close_perf_markers




subroutine write_perf_markers(H_rpe,U_rpe,V_rpe, TIME_rpe,  codesignQ, codesignD, IRHW, X_rpe, Y_rpe, &
                         & N, M, bits,NITER,NITSM,ICOUNT, sum_time, sum_lp_time)
use implicit_functions_SP

implicit none

REAL(Kind=4) :: H_rpe(N, M), U_rpe(N, M), V_rpe(N, M), TIME_rpe, X_rpe(N), Y_rpe(M)

integer :: N, M, IRHW, bits,NITER,NITSM,ICOUNT, NITAV
REAL(Kind=4) :: H(N, M), U(N, M), V(N, M), X(N), Y(M), TIME, sum_time, sum_lp_time
logical ::  codesignQ, codesignD

 character(len=150) path, file_name, experiment, simtime, codesD, codesQ, bit_count
INTEGER I,J


H(:, :)=H_rpe(:, :)
U(:, :)=U_rpe(:, :)
V(:, :)=V_rpe(:, :)
TIME=TIME_rpe
X(:)=X_rpe(:)
Y(:)=Y_rpe(:)

    write(unit=324,fmt='(I5, F6.2, F8.4, F8.4)') int(time), float(NITSM)/MAX(ICOUNT,1), &
                        & sum_time/MAX(ICOUNT,1), sum_lp_time/MAX(ICOUNT,1)

sum_lp_time=0.0d0
sum_time=0.0d0
NITAV=0
NITSM=0
ICOUNT=0

end subroutine write_perf_markers

subroutine  write_residual(R_rpe, exitcond, iteration, TIME_rpe,  codesignQ, codes, &
                        &  IRHW, X_rpe, Y_rpe, N, M, bits, ID_prec,EXP_NAME)
use implicit_functions_SP




integer :: N, M, IRHW, bits, iteration
REAL(Kind=4) ::  R_rpe(N, M), TIME_rpe, X_rpe(N), Y_rpe(M)
REAL(Kind=4) :: R(N, M), X(N), Y(M), TIME, exitcond
logical :: codesignQ, codes

 character(len=150) path, file_name, experiment, simtime, codesQ, codesD, bit_count, Precond, EXP_NAME, str_iter
INTEGER I,J

R(:,:) = R_rpe(:, :)
TIME=TIME_rpe
X(:)=X_rpe(:)
Y(:)=Y_rpe(:)


write(experiment,*) IRHW
write(simtime,*) INT(TIME)
write(codesQ,*) codesignQ
write(codesD,*) codes
write(bit_count,*) bits
write(Precond,*) ID_prec
write(str_iter,*) iteration

path ='../'//trim(adjustl(EXP_NAME))//& 
                        & '/Precon'//trim(adjustl(Precond))

file_name = '_R_exp'//trim(adjustl(experiment))//'_time'//trim(adjustl(simtime))//&
     & '_iter_'//trim(adjustl(str_iter))      &
     &  //'_codes_'//trim(adjustl(codesQ))//trim(adjustl(codesD))//'_bits'//trim(adjustl(bit_count))//'.txt'
Open(unit=599, file=trim(path)//trim(file_name), status='replace', form = 'formatted')

do I=2,N-1
  do J=1,M
    write(unit=599,fmt=*) X(I), Y(J), R(I,J), exitcond
  end do
end do
 close(unit=599)


end subroutine



subroutine  write_fields(H_rpe,U_rpe,V_rpe, TIME_rpe,  codesignQ, codesignD, IRHW, X_rpe, Y_rpe, N, M, bits, ID_prec,EXP_NAME)
use implicit_functions_SP


REAL(Kind=4) :: H_rpe(N, M), U_rpe(N, M), V_rpe(N, M), TIME_rpe, X_rpe(N), Y_rpe(M)

integer :: N, M, IRHW, bits
REAL(Kind=4) :: H(N, M), U(N, M), V(N, M), X(N), Y(M), TIME
logical :: codesignQ, codesignD

 character(len=150) path, file_name, experiment, simtime, codesQ, codesD, bit_count, Precond, EXP_NAME
INTEGER I,J

H(:, :)=H_rpe(:, :)
U(:, :)=U_rpe(:, :)
V(:, :)=V_rpe(:, :)
TIME=TIME_rpe
X(:)=X_rpe(:)
Y(:)=Y_rpe(:)


write(experiment,*) IRHW
write(simtime,*) INT(TIME)
write(codesQ,*) codesignQ
write(codesD,*) codesignD
write(bit_count,*) bits
write(Precond,*) ID_prec

path ='../'//trim(adjustl(EXP_NAME))//& 
                        & '/Precon'//trim(adjustl(Precond))

file_name = '_H_exp'//trim(adjustl(experiment))//'_time'//trim(adjustl(simtime))&
     &  //'_codes_'//trim(adjustl(codesQ))//trim(adjustl(codesD))//'_bits'//trim(adjustl(bit_count))//'.txt'
Open(unit=599, file=trim(path)//trim(file_name), status='replace', form = 'formatted')

do I=2,N-1
  do J=1,M
    write(unit=599,fmt=*) X(I), Y(J), H(I,J)
  end do
end do
 close(unit=599)

file_name = '_U_exp'//trim(adjustl(experiment))//'_time'//trim(adjustl(simtime))&
     &  //'_codes_'//trim(adjustl(codesQ))//trim(adjustl(codesD))//'_bits'//trim(adjustl(bit_count))//'.txt'
Open(unit=599, file=trim(path)//trim(file_name), status='replace', form = 'formatted')

do I=2,N-1
  do J=1,M
    write(unit=599, fmt=*) X(I), Y(J), U(I,J)
  end do
end do

 close(unit=599)


file_name = '_V_exp'//trim(adjustl(experiment))//'_time'//trim(adjustl(simtime))&
     &  //'_codes_'//trim(adjustl(codesQ))//trim(adjustl(codesD))//'_bits'//trim(adjustl(bit_count))//'.txt'
Open(unit=599, file=trim(path)//trim(file_name), status='replace', form = 'formatted')

do I=2,N-1
  do J=1,M
    write(unit=599, fmt=*) X(I), Y(J), V(I,J)
  end do
end do

 close(unit=599)


end subroutine



subroutine GCR_PRE(p,pfx,pfy,hx,hy,s,S_full,b,p0,pb,e1,e2,cor,ip  &
              & ,d,q,r,ar,n1,n2,gc1,gc2, &
           &    MGH1IHX, MGH2IHY, AC, BC, AD, BD,  &
           &  niter,nitsm,icount,error, p_T, sum_time,&
           &  sum_lp_time,ID_PREC, codes, save_time , &
           & TIME, codesQ, IRHW, DX_rpe, DY_rpe, Exit_cond, EXP_NAME&
           & , iprint, num_of_bits, DP_Depth, Alp_REL)
use implicit_functions_SP


implicit none
INTEGER, parameter :: kord=4, lord=kord-1

INTEGER :: n1, n2

REAL(Kind=4) :: p(n1,n2),pfx(n1,n2),pfy(n1,n2),hx(n1,n2),hy(n1,n2),s(n1,n2), &
     &   b(n1,n2),pb(n1,n2),p0(n1,n2), S_full,                   &
     &   e1(n1,n2),e2(n1,n2),cor(n1,n2),d(n1,n2),q(n1,n2),r(n1,n2),ar(n1,n2), &
     &   p_T(n1,n2), r_HP(n1,n2), r_true(n1,n2), r0_true(n1,n2), p_true(n1,n2), &
     &   p0_true(n1, n2), b_true(n1, n2), PMB(n1, n2), PMP0(n1, n2), qr(n1,n2)
REAL(Kind=4) :: MGH1IHX(n2), MGH2IHY(n2), AC(n2), BC(n2), AD(n2), BD(n2), Alp_REL(n1,n2)
!! preconditioning
REAL(Kind=4) :: qu(n1,n2), aqu(n1,n2),  A_c(n1,n2), B_c(n1,n2), C_c(n1,n2), ps(n1+1,n2), divi(n1,n2)
INTEGER :: ID_PREC
!! end preconditioning
INTEGER :: IP(n1)
REAL(Kind=4) :: GC1, GC2,error,qrror,  max_QX_QY, epa
REAL(Kind=4) :: res_lats0(n2), res_lats(n2)
REAL(Kind=4) :: start, finish, sum_time, sum_lp_time, startLP, endLP
integer :: num_of_bits

INTEGER :: niter,nitsm,icount
INTEGER :: iprint, DP_Depth
LOGICAL :: codes, save_time


REAL(Kind=4) :: x(n1,n2,lord),ax(n1,n2,lord),ax2(lord),axaqu(lord),del(lord),  &
     & a11(n1,n2),a12(n1,n2),a21(n1,n2),a22(n1,n2),b11(n1,n2),b22(n1,n2)
REAL(Kind=4) :: a11_t(n1,n2),a12_t(n1,n2),a21_t(n1,n2),a22_t(n1,n2),b11_t(n1,n2),b22_t(n1,n2)
REAL(Kind=4) :: err0,qrr0, rax, beta, errn, qrrn, x2, y2, T_step
INTEGER :: itr, J, I, l, ll, i1, it, itmn
REAL(Kind=4) :: eps, help1, help2, quotient, lowprectime, err_true, err0_true
REAL(Kind=4) :: Exit_cond

REAL(Kind=4) ::  TIME, DX_rpe(n1), DY_rpe(n2), errnm1
 character(len=150) :: EXP_NAME
LOGICAL :: codesQ, exiting

INTEGER :: IRHW , counter

double precision :: err0_dp, err_true_dp, r_true_dp(n1,n2),pfx_dp(n1,n2),pfy_dp(n1,n2)
ps(:,:)=0.0d0
divi(:,:)=0.0d0
p_T(:,:)=0.0
!write(*,*) num_of_bits
lowprectime=0.0d0


! end if

! eps=1.e-7 tested this, did not change anything
eps=1.e-5   !! original
itr=1000
niter=0
itmn=2
exiting=.false.

epa=1.e-30


 DO J=1,n2
   DO I=1,n1
    PMB(I, J)= p(I,J)-b(I,J)
    PMP0(I,J)= p(I,J)-p0(I,J)
   enddo
 enddo

p_true(:,:)=p(:,:)
p0_true(:,:)=p0(:,:)
b_true(:,:)=b(:,:)


 DO J=1+DP_Depth,n2-DP_Depth
   DO I=1,n1
    PMB(I, J)=PMB(I, J)
    PMP0(I,J)=PMP0(I,J)
    !p(I,J)=p(I,J)
    p0(I,J)=p0(I,J)
    b(I,J)=b(I,J)

   enddo
 enddo

!call cpu_time(startLP)
DO J=1,n2
  DO I=1,n1
    r(I,J)=rpe_0
    ar(I,J)=rpe_0
    qr(I,J)=rpe_0
  enddo
enddo

do l=1,lord
  DO J=1,n2
    DO I=1,n1
      x(I,J,l)=rpe_0
      ax(I,J,l)=rpe_0
    enddo
  enddo
enddo

!call cpu_time(endLP)
lowprectime=lowprectime + endLP-startLP



  !! should be 23

!! matrix entries
call LAP0_Piotr(a11,a12,a21,a22,b11,b22,                   &
     &          pb,p0,e1,e2,hx,hy,cor,ALP_rel,n1,n2,gc1,gc2)
!call lap0_depth(a11,a12,a21,a22,b11,b22,                   &
!     &          pb,p0,e1,e2,hx,hy,cor,n1,n2,gc1,gc2, &
!           &    MGH1IHX, MGH2IHY, AC, BC, AD, BD,  &
!           & DP_Depth)

  !! should be 23
     ! Preconditioning: init
!call precon_prep_depth(T_step, A_c, B_c, C_c,a11,a12,a21,a22,b11,b22,&
!              & p0,pfx,pfy,s,n1,n2,ip,ID_PREC)
if (ID_PREC==5) then
    call precon_prep_depth(T_step, A_c, ps, divi,a11,a12,a21,a22,b11,b22,&
              & p0,pfx,pfy,s,n1,n2,ip,ID_PREC,23, DP_Depth)
elseif (ID_PREC==6) then
  call diagoc(T_step, A_c,a11,a21,s,n1,n2, DP_Depth)
elseif (ID_PREC==7) then
    call precon_prep_depth(T_step, A_c, ps, divi,a11,a12,a21,a22,b11,b22,&
              & p0,pfx,pfy,s,n1,n2,ip,ID_PREC,23, DP_Depth)
endif

  !! should be 23
!call laplfirst(p(:,:),r_HP(:,:),a11,a12,a21,a22,b11,b22, p0,   &
!     &                           pfx,pfy,s,n1,n2,ip)
! replace my laplfirst_depth
!call laplfirst_depth(p(:,:),r_HP(:,:), a11,a12,a21,a22,b11,b22,PMP0,&
!     &                     pfx,pfy,S,n1,n2,IP, 23,DP_Depth)
!with Piotr's
  CALL PRFORC_ABS(p(:,:),pfx,pfy,pb,p0, &
       &      E1(:,:),E2(:,:),HX,HY,COR,n1,n2,IP,GC1,GC2,Alp_REL,1,1)
  call diver(r,pfx,pfy,hx,hy,s,n1,n2,ip,-1)

!! calculate initial residual
!call cpu_time(startLP)
 !! should be 23
err0_dp=0.0d0
err0=0.0d0
 DO J=1,n2
   DO I=1,n1
 
    r(I,J)=rpe_05*r(I,J)-(b(I,J))
      err0=err0+r(I,J)*r(I,J)
   !   err0_dp=err0_dp+r_HP(I,J)*r_HP(I,J)
    !write(*,*) I, J, r_HP(I,J), p(I,J),b(I,J),err0 &
    !      &,E1(I,J), E2(I,j), pb(I,J), p0(I, J)

   enddo
    !read(*,*)
 enddo



if (iprint==1) then

! write(*,*) (maxval(abs(r_HP(:,J))), J=1,n2) 

endif


!write(*,*) 'wo'
!call prforc_depth(p,pfx,pfy,pb,p0,e1,e2,hx,hy,cor,n1,n2,ip,gc1,gc2,1,1, num_of_bits, DP_Depth)
!call prforc(p,pfx,pfy,pb,p0,e1,e2,hx,hy,cor,n1,n2,ip,gc1,gc2,1,1)
!write(*,*) 'prforc_depth'
!call diver_depth(r,pfx,pfy,hx,hy,s,n1,n2,ip,-1, num_of_bits, DP_Depth)
!call diver(r,pfx,pfy,hx,hy,s,n1,n2,ip,-1)
!write(*,*) 'diver_depth'

!call cpu_time(start) 

!err0=0.0d0


    err0=sqrt(err0)
    err0_dp=sqrt(err0_dp)
      ! write(*,*)err0, err0_dp, abs(err0 -err0_dp)/err0_dp, counter
      ! read(*,*)

    errnm1=err0

!! should be num_of_bits
!call cpu_time(startLP)


!err0 =maxval(ABS(r_HP(:,:)))

if (iprint==1) then
    call write_residual(r,eps*Exit_cond, niter, TIME, codesQ, codes, IRHW, DX_rpe, DY_rpe,&
                     & n1, n2, num_of_bits, 7 ,EXP_NAME)
endif

if (iprint==1) then
    call write_residual(r,eps*Exit_cond, niter, TIME, codesQ, codes, IRHW, DX_rpe, DY_rpe,&
                     & n1, n2, num_of_bits, 8 ,EXP_NAME)
endif







!call cpu_time(startLP)

call precon(r,x(:,:,1),ax(:,:,1), T_step,  A_c, ps, divi,a11,a12,a21,a22,b11,b22,p0,  &
                &   pfx,pfy,s,S_full,n1,n2,ip,ID_PREC, num_of_bits, DP_Depth)
qrr0=rpe_0
 DO J=1,n2
   DO I=1,n1
      qrr0=qrr0+x(I,J,1)*x(I,J,1)
   enddo
 enddo
qrr0=sqrt(qrr0)

  call lapl_depth(x(:,:,1),ax(:,:,1), A11,A12,A21,A22,B11,B22,P0,pfx,pfy,S,n1,n2,IP,num_of_bits,DP_Depth)




  !! rewritten to change precision at poles as desired    
  !! this sum could be problematic
 !! should be 23
      DO J=1+DP_Depth,n2-DP_Depth
        DO I=1,n1
      ax(I,J,1)=rpe_05*ax(I,J,1)-x(I,J,1)
        enddo
      enddo

      DO J=1,DP_Depth
        DO I=1,n1
      ax(I,J,1)=rpe_05*ax(I,J,1)-x(I,J,1)

        enddo
      enddo

      DO J=n2+1-DP_Depth,n2
        DO I=1,n1
      ax(I,J,1)=rpe_05*ax(I,J,1)-x(I,J,1)
        enddo
      enddo

   !! end of rewrite

!call cpu_time(endLP)
lowprectime=lowprectime + endLP-startLP

do it=1,itr



 If (it==1000) then
   save_time=.true.
   exit
 endif
  do l=1,lord
    !write(*,*) 'before ',niter

    !write(*,*) niter

     
    ax2(l)=rpe_0
    rax=rpe_0
    DO J=1,n2
      DO I=1,n1
        rax=rax+r(I,J)*ax(I,J,l)
        ax2(l)=ax2(l)+ax(I,J,l)*ax(I,J,l)
      enddo
    enddo

    ax2(l)=max(epa,ax2(l))
    beta=-rax/ax2(l)
    errn=0.0d0



!23 !!was 23 before exit criterium   
      DO J=1,n2
        DO I=1,n1
         p_T(I,J)=p_T(I,J) +beta* x(I,J,l) 
         p(I,J)=p(I,J) +beta* x(I,J,l) 
         r(I,J)  =r(I,J)   +beta*ax(I,J,l) 
        enddo
      enddo


      DO J=1,n2
        DO I=1,n1
         errn=errn+r(I,J)*r(I,J)
        enddo
      enddo




   !! end of rewrite
!    DO J=1,n2
!      DO I=1,n1
      ! If(J<3 .or. J>62) then
      ! write(*,*) I, J, 'values', r(I,J)%val, ax(I,J,l)%val, beta%val*ax(I,J,l)%val, r(I,J)%val   +beta%val*ax(I,J,l)%val
      ! end if
!         p_T(I,J)=p_T(I,J) +beta* x(I,J,l) 
        !if (.not. codes) then
        ! p(I,J)  =p(I,J)   +beta* x(I,J,l) 
!          r_HP(I,J)  =r_HP(I,J)   +beta*ax(I,J,l) 
!          r(I,J)  = r_HP(I,J)
!          errn=errn+r(I,J)*r(I,J)

         !else
         ! p_T(I,J)=p_T(I,J) +beta* x(I,J,l) 
         ! r(I,J)  =r(I,J)   +beta*ax(I,J,l)  

         ! errn=errn+r(I,J)*r(I,J)
         ! end if     
   
!      enddo
      ! If(J<3 .or. J>62) then
      !  read(*,*) 
      ! endif
!    enddo

!!! true residual

!!! true residual

r_true_dp(:,:)=0.0d0
!call laplfirst(p_true(:,:)+p_T(:,:),r_true(:,:),a11_t,a12_t,a21_t,a22_t,b11_t,b22_t, p0_true,   &
!     &                           pfx,pfy,s,n1,n2,ip)

!  CALL PRFORC_ABS(p_true(:,:)+p_T(:,:),pfx,pfy,pb,p0, &
!       &      E1(:,:),E2(:,:),HX,HY,COR,n1,n2,IP,GC1,GC2,Alp_REL,1,1)
!  call diver(r_true,pfx,pfy,hx,hy,s,n1,n2,ip,-1)
  CALL PRFORC_ABS_dp(dble(p(:,:)),pfx_dp,pfy_dp,dble(pb),dble(p0), &
       &      dble(E1(:,:)),dble(E2(:,:)),dble(HX),dble(HY),dble(COR), &
       &      n1,n2,IP,dble(GC1),dble(GC2),dble(Alp_REL),1,1)
  call diver_dp(r_true_dp,pfx_dp,pfy_dp,dble(hx),dble(hy),dble(s),n1,n2,ip,-1)
  
  DO J=1,n2
  DO I=1,n1
    r_true_dp(I,J)=0.5d0*r_true_dp(I,J)-(dble(p(I,J))-dble(p_true(I,J))+dble(b(I,J)))
 ! write(*,*), i, J, P(i,J)
  enddo
enddo
err_true_dp=0.0
      DO J=1,n2
        DO I=1,n1
         err_true_dp=err_true_dp+r_true_dp(I,J)*r_true_dp(I,J)
        enddo
      enddo
    !  write(*,*) 'errn_true model view', sqrt(err_true_dp)
   !! how well does the solver converge?
   r_true_dp(:,:)=0.0d0
!call laplfirst(p_true(:,:)+p_T(:,:),r_true(:,:),a11_t,a12_t,a21_t,a22_t,b11_t,b22_t, p0_true,   &
!     &                           pfx,pfy,s,n1,n2,ip)

!  CALL PRFORC_ABS(p_true(:,:)+p_T(:,:),pfx,pfy,pb,p0, &
!       &      E1(:,:),E2(:,:),HX,HY,COR,n1,n2,IP,GC1,GC2,Alp_REL,1,1)
!  call diver(r_true,pfx,pfy,hx,hy,s,n1,n2,ip,-1)
  CALL PRFORC_ABS_dp(dble(p_true(:,:))+dble(p_T(:,:)),pfx_dp,pfy_dp,dble(pb),dble(p0), &
       &      dble(E1(:,:)),dble(E2(:,:)),dble(HX),dble(HY),dble(COR), &
       &      n1,n2,IP,dble(GC1),dble(GC2),dble(Alp_REL),1,1)
  call diver_dp(r_true_dp,pfx_dp,pfy_dp,dble(hx),dble(hy),dble(s),n1,n2,ip,-1)
  
  DO J=1,n2
  DO I=1,n1
    r_true_dp(I,J)=0.5d0*r_true_dp(I,J)-(dble(p_T(I,J))+dble(b(I,J)))
 ! write(*,*), i, J, P(i,J)
  enddo
enddo
err_true_dp=0.0
      DO J=1,n2
        DO I=1,n1
         err_true_dp=err_true_dp+r_true_dp(I,J)*r_true_dp(I,J)
        enddo
      enddo
     ! write(*,*) 'errn_true tendency view', sqrt(err_true_dp)
      if (iprint==1) then
        !! what is the model view on convergence, what is actually added?
r_true_dp(:,:)=0.0d0
!call laplfirst(p_true(:,:)+p_T(:,:),r_true(:,:),a11_t,a12_t,a21_t,a22_t,b11_t,b22_t, p0_true,   &
!     &                           pfx,pfy,s,n1,n2,ip)

!  CALL PRFORC_ABS(p_true(:,:)+p_T(:,:),pfx,pfy,pb,p0, &
!       &      E1(:,:),E2(:,:),HX,HY,COR,n1,n2,IP,GC1,GC2,Alp_REL,1,1)
!  call diver(r_true,pfx,pfy,hx,hy,s,n1,n2,ip,-1)
  CALL PRFORC_ABS_dp(dble(p(:,:)),pfx_dp,pfy_dp,dble(pb),dble(p0), &
       &      dble(E1(:,:)),dble(E2(:,:)),dble(HX),dble(HY),dble(COR), &
       &      n1,n2,IP,dble(GC1),dble(GC2),dble(Alp_REL),1,1)
  call diver_dp(r_true_dp,pfx_dp,pfy_dp,dble(hx),dble(hy),dble(s),n1,n2,ip,-1)
  
  DO J=1,n2
  DO I=1,n1
    r_true_dp(I,J)=0.5d0*r_true_dp(I,J)-(dble(p(I,J))-dble(p_true(I,J))+dble(b(I,J)))
 ! write(*,*), i, J, P(i,J)
  enddo
enddo
err_true_dp=0.0
      DO J=1,n2
        DO I=1,n1
         err_true_dp=err_true_dp+r_true_dp(I,J)*r_true_dp(I,J)
        enddo
      enddo
      !write(*,*) 'errn_true model view', sqrt(err_true_dp)

    call write_residual(real(dble(p(:,:))-dble(p_true(:,:))),eps*Exit_cond, &
            & niter+1, TIME, codesQ, codes,&
           &  IRHW, DX_rpe, DY_rpe, n1, n2, num_of_bits, 1 ,EXP_NAME)
    call write_residual(real(r_true_dp),eps*Exit_cond, niter+1, TIME, codesQ, codes,&
           &  IRHW, DX_rpe, DY_rpe, n1, n2, num_of_bits, 7 ,EXP_NAME)

   !! how well does the solver converge?
   r_true_dp(:,:)=0.0d0
!call laplfirst(p_true(:,:)+p_T(:,:),r_true(:,:),a11_t,a12_t,a21_t,a22_t,b11_t,b22_t, p0_true,   &
!     &                           pfx,pfy,s,n1,n2,ip)

!  CALL PRFORC_ABS(p_true(:,:)+p_T(:,:),pfx,pfy,pb,p0, &
!       &      E1(:,:),E2(:,:),HX,HY,COR,n1,n2,IP,GC1,GC2,Alp_REL,1,1)
!  call diver(r_true,pfx,pfy,hx,hy,s,n1,n2,ip,-1)
  CALL PRFORC_ABS_dp(dble(p_true(:,:))+dble(p_T(:,:)),pfx_dp,pfy_dp,dble(pb),dble(p0), &
       &      dble(E1(:,:)),dble(E2(:,:)),dble(HX),dble(HY),dble(COR), &
       &      n1,n2,IP,dble(GC1),dble(GC2),dble(Alp_REL),1,1)
  call diver_dp(r_true_dp,pfx_dp,pfy_dp,dble(hx),dble(hy),dble(s),n1,n2,ip,-1)
  
  DO J=1,n2
  DO I=1,n1
    r_true_dp(I,J)=0.5d0*r_true_dp(I,J)-(dble(p_T(I,J))+dble(b(I,J)))
 ! write(*,*), i, J, P(i,J)
  enddo
enddo
err_true_dp=0.0
      DO J=1,n2
        DO I=1,n1
         err_true_dp=err_true_dp+r_true_dp(I,J)*r_true_dp(I,J)
        enddo
      enddo
     ! write(*,*) 'errn_true tendency view', sqrt(err_true_dp)
!!! end true residual

    call write_residual(real(p_T),eps*Exit_cond, &
            & niter+1, TIME, codesQ, codes,&
           &  IRHW, DX_rpe, DY_rpe, n1, n2, num_of_bits, 2 ,EXP_NAME)
    call write_residual(real(r_true_dp),eps*Exit_cond, niter+1, TIME, codesQ, codes,&
           &  IRHW, DX_rpe, DY_rpe, n1, n2, num_of_bits, 8 ,EXP_NAME)
endif


!    DO J=1,n2
!      res_lats(J)=rpe_0
!      DO I=1,n1
!   
!        res_lats(J)=res_lats(J)+r(I,J)*r(I,J)
!    
!      enddo
!    enddo
!      write(*,*) it, ( sqrt(res_lats(J)), J=1,n2)
!read(*,*)
   ! write(*,*) ( sqrt(res_lats(J)%val/res_lats0(J)%val), J=1,n2)
   ! write(*,*) 'Error: ', error, errn
   ! write(*,*) 'error .lt. eps', error .lt. eps, eps
   ! read(*,*)

    errn=sqrt(errn)
   write(*,*) niter, errn, err0
    if( niter+2 > itmn) exiting=.true.
    if(errn.ge.errnm1) exiting=.true.
    !exiting=.true.  ! one iteration
    errnm1=errn
!if(maxval(ABS(r_HP(:,:))) .lt. eps*Exit_cond) exit


! 1) get timestep length Delta_t from linear stability argument
!write(*,*) 'in precon ADI'
!max_QX_QY=(abs(A11(1,1))+abs(A21(1,1)))/(rpe_2*S(1,1))
!DO J=1,n2
!
!  DO I=1,n1
!     max_QX_QY=max(max_QX_QY,(abs(A11(I,J))+abs(A21(I,J)))/(rpe_2*S(I,J)) )     
!  ENDDO
!
!ENDDO
!write(*,*) 'precon_ADI', rpe_05/max_QX_QY



    call precon(r,qu, aqu , T_step,  A_c, ps, divi,a11,a12,a21,a22,b11,b22,p0,   &
                &   pfx,pfy,s,S_full,n1,n2,ip,ID_PREC, num_of_bits, DP_Depth)
  !    do j=1,1
  !     do i=1,n1
  !      write(*,*)  'qu',i, J, qu(i,j)
  !     read(*,*)
  !     enddo
  !    enddo
! 
 !23
  call lapl_depth(qu(:,:),aqu(:,:), A11,A12,A21,A22,B11,B22,P0,pfx,pfy,S,n1,n2,IP , num_of_bits,DP_Depth)




  !! rewritten to change precision at poles as desired    
  !! this sum could be problematic
  !! should be 23

      DO J=1+DP_Depth,n2-DP_Depth
        DO I=1,n1
      aqu(I,J)=rpe_05*aqu(I,J)-qu(I,J)
        enddo
      enddo

      DO J=1,DP_Depth
        DO I=1,n1
      aqu(I,J)=rpe_05*aqu(I,J)-qu(I,J)
        enddo
      enddo

      DO J=n2+1-DP_Depth,n2
        DO I=1,n1
      aqu(I,J)=rpe_05*aqu(I,J)-qu(I,J)
        enddo
      enddo
    niter=niter+1
   !! end of rewrite

lowprectime=lowprectime + endLP-startLP
!write(*,*) 'second precon', it



    do ll=1,l
      axaqu(ll)=rpe_0
      DO J=1,n2
        DO I=1,n1
          axaqu(ll)=axaqu(ll)+ax(I,J,ll)*aqu(I,J)
        enddo
      enddo
      del(ll)=-axaqu(ll)/ax2(ll)

 !     del(ll)=max(del(ll),rpe_05)
    enddo


    if(l.lt.lord) then
 !! should be 23   
      DO J=1,n2
        DO I=1,n1
          x(I,J,l+1)= qu(I,J)
          ax(I,J,l+1)=aqu(I,J)
        enddo
      enddo


      do ll=1,l

 !! should be 23   
      DO J=1,n2
        DO I=1,n1
            x(I,J,l+1)= x(I,J,l+1)+del(ll)* x(I,J,ll)
            ax(I,J,l+1)=ax(I,J,l+1)+del(ll)*ax(I,J,ll)
        enddo
      enddo



      enddo

    else
  !! rewritten to change precision at poles as desired 
 !! should be 23   
      DO J=1,n2
        DO I=1,n1
          x(I,J,1)= qu(I,J)+del(1)* x(I,J,1)
          ax(I,J,1)=aqu(I,J)+del(1)*ax(I,J,1)
        enddo
      enddo


   !! end of rewrite

      do ll=2,l

 !! should be 23   
      DO J=1,n2
        DO I=1,n1
            x(I,J,1 )= x(I,J,1)+del(ll)* x(I,J,ll)
            x(I,J,ll)=rpe_0
            ax(I,J,1 )=ax(I,J,1)+del(ll)*ax(I,J,ll)
            ax(I,J,ll)=rpe_0
        enddo
      enddo




      enddo

    endif


  if(exiting .eqv. .true.) exit
 
  enddo

  if(exiting .eqv. .true.) exit !! to replace the go to 200


end do
!write(*,*) niter
!  200
!niter=it
qrrn=rpe_0
 DO J=1,n2
   DO I=1,n1
      qrrn=qrrn+x(I,J,1)*x(I,J,1)
   enddo
 enddo
qrrn=sqrt(qrrn)
qrror=qrrn/qrr0

if (iprint==1) then
 write(*,*) 'Qerror', qrror 
call lap0_depth(a11,a12,a21,a22,b11,b22,                   &
     &          pb,p0,e1,e2,hx,hy,cor,n1,n2,gc1,gc2, &
           &    MGH1IHX, MGH2IHY, AC, BC, AD, BD,  &
           & 0)



! end matrices
r0_true(:,:)=0.0d0
call laplfirst(p_true(:,:),r0_true(:,:),a11,a12,a21,a22,b11,b22, p0_true,   &
     &                           pfx,pfy,s,n1,n2,ip)

DO J=1,n2
  DO I=1,n1
    r0_true(I,J)=0.5d0*r0_true(I,J)-(p_true(I,J)-b_true(I,J))
 ! write(*,*), i, J, P(i,J)
  enddo
enddo

r_true(:,:)=0.0d0
call laplfirst(p_true(:,:)+p_T(:,:),r_true(:,:),a11,a12,a21,a22,b11,b22, p0_true,   &
     &                           pfx,pfy,s,n1,n2,ip)

DO J=1,n2
  DO I=1,n1
    r_true(I,J)=0.5d0*r_true(I,J)-(p_true(I,J)+p_T(I,J)-b_true(I,J))
 ! write(*,*), i, J, P(i,J)
  enddo
enddo
err_true=0.0d0
err0_true=0.0d0

DO J=1,n2
  DO I=1,n1
    err0_true=err0_true+r0_true(I,J)*r0_true(I,J)
    err_true=err_true+r_true(I,J)*r_true(I,J)
  enddo
enddo
write(*,*) niter
 write(*,*) 'truth ACC',sqrt(err_true/err0_true),'max(rtrue)' ,&
           & maxval(ABS(r_true(:,:))),'max(r0true)', maxval(ABS(r0_true(:,:))), 'max(r)',maxval(ABS(r_HP(:,:))),&
           &'max(r0)',err0 , 'EXIT', eps*Exit_cond 

endif

!write(*,*) niter

icount=icount+1
!write(*,*) 'iterations', niter
nitsm=nitsm+niter
sum_time=sum_time+(finish-start)
sum_lp_time=sum_lp_time+lowprectime
end subroutine




subroutine precon_prep_depth(T_step,A, ps, divi, A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP,ID_PREC,num_of_bits, DP_Depth)

use implicit_functions_SP

implicit none
REAL(Kind=4) ::  R(N,M),QU(N,M),A11(N,M),A12(N,M),A21(N,M),A22(N,M),  &
    &      B11(N,M),B22(N,M),P0(N,M),U(N,M),V(N,M),S(N,M), S_L(N, M)
REAL(Kind=4) ::   A(N,M), B(N,M), C(N,M),ps(N+1,M), divi(N,M)
REAL(Kind=4) ::  AQU(N,M), T_step, Delta_t, max_QX_QY
INTEGER :: IP(N), ID_PREC, num_of_bits, DP_Depth
INTEGER :: N, M, I, J




! truncating
DO J=1+DP_Depth,M-DP_Depth
 DO I=1,N

  S_L(I, J)=S(I, J)
end do
end do



IF (ID_PREC==5) then !! ADI type preconditioner

! 1) get timestep length for preconditioner Delta_t from linear stability argument
!write(*,*) 'in precon ADI'
!max_QX_QY=(abs(A11(1,1))+abs(A21(1,1)))/(rpe_2*S(1,1))
!DO J=1,M
!  DO I=1,N
!     max_QX_QY=max(max_QX_QY,(abs(A11(I,J))+abs(A21(I,J)))/(rpe_2*S(I,J)) )     
!  ENDDO
!ENDDO
!T_step=rpe_025/max_QX_QY !0.92d0!
!Delta_t=rpe_1/T_step
!write(*,*) 'old working',Delta_t

!max_QX_QY=-1.e15
!DO J=1,M
!  DO I=1,N
!        !  write(*,*) i, j, a21(i,j), s(i,j)  , 4.0d0*max_QX_QY
!     max_QX_QY=max(max_QX_QY,(abs(A21(I,J)))/(4.0d0*S(I,J)) )   
!        !  write(*,*) i, j, a21(i,j), s(i,j)  , 4.0d0*max_QX_QY
!        !  read(*,*)
!  ENDDO
!ENDDO

max_QX_QY=-1.e15
do j=2,m-1
 do i=1,n
  max_QX_QY=amax1(max_QX_QY, 0.5*abs(a21(i,j+1)+a21(i,j-1))/s(i,j))
 enddo
enddo
do i=1,n
  max_QX_QY=amax1(max_QX_QY,0.5*abs(a21(i,2)+a21(ip(i),1))/s(i,1))
  max_QX_QY=amax1(max_QX_QY,0.5*abs(a21(ip(i),m)+a21(i,m-1))/s(i,m-1))
enddo

T_step=max_QX_QY !0.92d0!
Delta_t=rpe_1/T_step

!write(*,*) 'Delta_T',Delta_t

!max_QX_QY=(rpe_2*abs(A21(1,1)))/( ( (2.0d0*acos(-1.0d0)/Dfloat(M)) *6371.22E+03 )**2 )
!DO J=1,M
!  DO I=1,N
!     max_QX_QY=max(max_QX_QY, (rpe_2*abs(A21(I,J)))/( (2.0d0* (acos(-1.0d0)/Dfloat(M)) *6371.22E+03 )**2 ) )    
!  ENDDO
!ENDDO
!T_step=1.0d0/max_QX_QY !0.92d0!
!Delta_t=T_step
!write(*,*) Delta_t



      do J=1,M
       ps(2,J)=0.0d0
       ps(3,J)=0.0d0 
      enddo

      DO J=1+DP_Depth,M-DP_Depth
        DO I=2,N-1
     A(I,J)= -Delta_t*A11(I-1,J)/(rpe_2*S_L(I,J))
        enddo
      enddo

      DO J=1+DP_Depth,M-DP_Depth
        DO I=2,N-1
     B(I,J)= rpe_1 +Delta_t+ Delta_t* (A11(I+1,J)+A11(I-1,J))/(rpe_2*S_L(I,J))
        enddo
      enddo

      DO J=1+DP_Depth,M-DP_Depth
        DO I=2,N-1
     C(I,J)= -Delta_t*A11(I+1,J)/(rpe_2*S_L(I,J))
        enddo
      enddo

      DO J=1+DP_Depth,M-DP_Depth
        DO I=2,N-1
     divi(I,J)=A(I,J)*ps(I,J)+B(I,J)
     ps(I+2,J)= -C(I,J)/divi(I,J)
        enddo
      enddo


      DO J=1,DP_Depth
        DO I=2,N-1
     A(I,J)= -Delta_t*A11(I-1,J)/(rpe_2*S(I,J))
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=2,N-1
     A(I,J)= -Delta_t*A11(I-1,J)/(rpe_2*S(I,J))
        enddo
      enddo


      DO J=1,DP_Depth
        DO I=2,N-1
     B(I,J)= rpe_1 +Delta_t+ Delta_t* (A11(I+1,J)+A11(I-1,J))/(rpe_2*S(I,J))
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=2,N-1
     B(I,J)= rpe_1 +Delta_t+ Delta_t* (A11(I+1,J)+A11(I-1,J))/(rpe_2*S(I,J))
        enddo
      enddo


      DO J=1,DP_Depth
        DO I=2,N-1
     C(I,J)= -Delta_t*A11(I+1,J)/(rpe_2*S(I,J))
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=2,N-1
     C(I,J)= -Delta_t*A11(I+1,J)/(rpe_2*S(I,J))
        enddo
      enddo



  CALL XBC(A,N,M)
  CALL XBC(B,N,M)
  CALL XBC(C,N,M)





      DO J=1,DP_Depth
        DO I=2,N-1
     divi(I,J)=A(I,J)*ps(I,J)+B(I,J)
     ps(I+2,J)= -C(I,J)/divi(I,J)
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=2,N-1
     divi(I,J)=A(I,J)*ps(I,J)+B(I,J)
     ps(I+2,J)= -C(I,J)/divi(I,J)
        enddo
      enddo

elseIF (ID_PREC==7) then !! ADI type preconditioner

max_QX_QY=-1.e15
do j=2,m-1
 do i=1,n
  max_QX_QY=amax1(max_QX_QY, 0.5*abs(a21(i,j+1)+a21(i,j-1))/s(i,j))
 enddo
enddo
do i=1,n
  max_QX_QY=amax1(max_QX_QY,0.5*abs(a21(i,2)+a21(ip(i),1))/s(i,1))
  max_QX_QY=amax1(max_QX_QY,0.5*abs(a21(ip(i),m)+a21(i,m-1))/s(i,m-1))
enddo

T_step=max_QX_QY !0.92d0!


end if



end subroutine

subroutine precon(R,QU,AQU, T_step,A, ps, divi, A11,A12,A21,A22,B11,B22,P0,U,V,S, S_full,N,M,IP,ID_PREC, num_of_bits,DP_Depth)

use implicit_functions_SP

implicit none
REAL(Kind=4) ::  R(N,M),QU(N,M),A11(N,M),A12(N,M),A21(N,M),A22(N,M),  &
    &      B11(N,M),B22(N,M),P0(N,M),U(N,M),V(N,M),S(N,M)
REAL(Kind=4) ::  A(N,M), B(N,M), C(N,M), ps(N+1,M), divi(N,M)
REAL(Kind=4) ::  AQU(N,M), T_step, S_full
INTEGER :: IP(N), ID_PREC, num_of_bits,DP_Depth
INTEGER :: N, M, I, J

IF (ID_PREC==5) then !! ADI type preconditioner



call precon_ADI(R,QU , T_step, A, ps, divi, A11,A12,A21,A22,B11,B22,P0,U,V,S,S_full,N,M,IP,ID_PREC, num_of_bits,DP_Depth)
 ! write(*,*) 'does he get out?'
!! regardless of preconditioners L(q^(n+1)) is needed
 
  ! write(*,*) 'AQU complete'
elseif (ID_PREC==6) then
call precon_Jac(R,QU , S_full, A, ps, divi, A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP,ID_PREC, num_of_bits,DP_Depth)
elseIF (ID_PREC==7) then !! ADI type preconditioner



call precon_ADI_Piotr(R,QU , T_step, A, ps, divi, A11,A12,A21,A22,B11,B22,P0,U,V,S,S_full,N,M,IP,ID_PREC, num_of_bits,DP_Depth)
 ! write(*,*) 'does he get out?'
!! regardless of preconditioners L(q^(n+1)) is needed
 
  ! write(*,*) 'AQU complete'
elseif(ID_PREC==0) then

QU(:,:)=R(:,:)

end if



end subroutine

!   implement ADI type preconditioner based on q_n+1 =q_n + dt{ Lz(q_n+1) + Lm(q_n) + H(q_n+1) -R}
!
SUBROUTINE precon_ADI(R,QU , T_step,  A, ps, divi,A11,A12,A21,A22,B11,B22,P0,U,V,S,S_full,N,M,IP,ID_PREC, num_of_bits,DP_Depth)
use implicit_functions_SP

implicit none

REAL(Kind=4) :: R(N,M),QU(N,M),A11(N,M),A12(N,M),A21(N,M),A22(N,M),  &
    &      B11(N,M),B22(N,M),P0(N,M),U(N,M),V(N,M),S(N,M),  C11(N,M)
INTEGER :: IP(N), ID_PREC, num_of_bits,DP_Depth
INTEGER :: N, M

REAL(Kind=4) :: max_QX_QY, F(N,M), rhs(N,M), qs(N+1,M,0:4), ps(N+1,M), ws(N+1,M,0:4), A(N,M), B(N,M), C(N,M), divi(N,M)
REAL(Kind=4) :: aa(M,1:4), deti, det40, det41, det42, det43, det44, det3,   &
                           & d11, d12, d13, d14, d21, d22, d23, d24,       &
                           & d31, d32, d33, d34, d41, d42, d43, d44,       &
                           & s1, s2, s3, s4
REAL(Kind=4) ::T_step, Delta_t, dn, dni, S_full, util, vtil, swcp
!REAL(Kind=4) :: 
integer :: iter, max_iter, time_scale  !! number of richardson iterations
INTEGER :: I, J, iteration



max_iter=2
swcp=rpe_1
!initialize the inverse of R with 0
DO J=1,M
  DO I=1,N
    QU(I,J)=rpe_0
  end do
end do


Delta_t=rpe_1/T_step
!Delta_t=1.0d0/T_step
!write(*,*) Delta_t




! 2) loop of following ADi iterations
do iteration=1,max_iter
 ! 2.1 calculate new right hand side using old QU and R to get tilde{tilde(R)}
If(iteration==1) then

! 1) first iteration
  DO J=1,M
    DO I=1,N    
 
      C11(I,J)=Delta_t*0.5*a11(i,j)/s(i,J)
    END DO
  END DO

      DO J=1+DP_Depth,M-DP_Depth
        DO I=1,N
      rhs(I,J)=Delta_t*( - R(I,J))
        enddo
      enddo

      DO J=1,DP_Depth
        DO I=1,N
      rhs(I,J)=Delta_t*( - R(I,J))
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=1,N
      rhs(I,J)=Delta_t*( - R(I,J))
        enddo
      enddo


else
  
      do j=2,m-1
       do i=2,n-1
        U(i,j)=QU(i+1,j)-QU(i-1,j)
        V(i,j)=QU(i,j+1)-QU(i,j-1)
       enddo
      enddo
      do i=2,n-1
       U(i,1)=QU(i+1,1)-QU(i-1,1)
       U(i,m)=QU(i+1,m)-QU(i-1,m)
       V(i,1)=QU(i,2)-QU(ip(i),1)
       V(i,m)=QU(ip(i),m)-QU(i,m-1)
      enddo

       CALL XBC(U,N,M)
       CALL XBC(V,N,M)

      do j=1,m
      do i=1,n
       util=                swcp*(a12(i,j)*V(i,j)+b11(i,j)*QU(i,j))
       vtil=V(i,j)*a21(i,j)+swcp*(a22(i,j)*U(i,j)+b22(i,j)*QU(i,j))
       U(i,j)=util
       V(i,j)=vtil
      enddo
      enddo

  CALL XBC(U,N,M)
  CALL XBC(V,N,M)

  DO J=2,M-1
    DO I=2,N-1
      F(I,J)= (U(I+1,J)-U(I-1,J)+V(I,J+1)-V(I,J-1))/(S(I,J))
    end do
  end do    
 

  DO I=2,N-1
    F(I,1)= (U(I+1,1)-U(I-1,1)+(V(I,2)+V(I,1)))/(S(I,1))
    F(I,M)= (U(I+1,M)-U(I-1,M)-(V(I,M)+V(I,M-1)))/(S(I,M))
  ENDDO


  CALL XBC(F,N,M)
      do j=1,m
      do i=1,n
       F(i,j)=rpe_05*F(i,j) !-p(i,j)
      enddo
      enddo

  DO J=1,M
    DO I=1,N    
 
      rhs(I,J)=QU(I,J)+Delta_t*(  F(I,J)   &
                          &   - R(I,J))
    END DO
  END DO
 CALL XBC(rhs,N,M)

end if
!DO J=1,M
!
!   rhs(0,J)= rhs(N-2,J)
!
!end do 
  !write(*,*) 'right hand side finished'
 ! 2.2 solve for new QU implicitly, aka solving a tridiagonal system with periodic BC 

  ! calculate and store coefficients ai bi and ci

  !DO J=1,M
  ! DO I=2,N-1
  !   write(*,*) A(I,J)%val, B(I,J)%val, C(I,J)%val
  ! end do
  ! read(*,*)
  !end do 
  
  !DO J=1,M
  !   C(1,J)= C(N-1,J)
  !   C(0,J)= C(N-2,J)
  !end do 
  ! write(*,*) 'A B C coeffficients finished'
  ! calculate p and q's, p is the same for all, 4 q's are needed , first with (bc1, bc2) =(0,0), second with (1,0), third with (0,1)

  ! initialize the 5 linear subsystems with boundary conditions (left boundary)
 
  DO J=1,M
    ps(2,J)=rpe_0
    ps(3,J)=rpe_0

    qs(2,J,0)  =rpe_0
    qs(3,J,0)  =rpe_0

    qs(2,J,1)  =rpe_1
    qs(3,J,1)  =rpe_0

    qs(2,J,2)  =rpe_0
    qs(3,J,2)  =rpe_1

    qs(2,J,3)  =rpe_0
    qs(3,J,3)  =rpe_0

    qs(2,J,4)  =rpe_0
    qs(3,J,4)  =rpe_0
  end do


  !! rewritten to change precision at poles as desired   

!      DO J=1+DP_Depth,M-DP_Depth
!        DO I=2,N-1
!     divi(I,J)=A(I,J)*ps(I,J)+B(I,J)
!     ps(I+2,J)= -C(I,J)/divi(I,J)
!        enddo
!      enddo

!      DO J=1,DP_Depth
!        DO I=2,N-1
!     divi(I,J)=A(I,J)*ps(I,J)+B(I,J)
!     ps(I+2,J)= -C(I,J)/divi(I,J)
!        enddo
!      enddo
!
!      DO J=M+1-DP_Depth,M
!        DO I=2,N-1
!     divi(I,J)=A(I,J)*ps(I,J)+B(I,J)
!     ps(I+2,J)= -C(I,J)/divi(I,J)
!        enddo
!      enddo


   !! end of rewrite


!  DO J=1,M
!   DO I=2,N-1
!
!     qs(I+2,J,0)= (rhs(I,J)-A(I,J)*qs(I,J,0))/divi(I,J)
!     qs(I+2,J,1)= (rpe_0   -A(I,J)*qs(I,J,1))/divi(I,J)
!     qs(I+2,J,2)= (rpe_0   -A(I,J)*qs(I,J,2))/divi(I,J)
!     qs(I+2,J,3)= (rpe_0   -A(I,J)*qs(I,J,3))/divi(I,J)
!     qs(I+2,J,4)= (rpe_0   -A(I,J)*qs(I,J,4))/divi(I,J)
!
!   end do
!  end do 
  !! rewritten to change precision at poles as desired    

      DO J=1+DP_Depth,M-DP_Depth
        DO I=2,N-1
       divi(I,J)=c11(i+1,j)+c11(i-1,j)*(rpe_1-ps(i,j)) +(Delta_t+rpe_1)
        divi(I,J)=rpe_1/divi(I,J)
     ps(i+2,j)=                  c11(i+1,j)*divi(I,J) 
     qs(I+2,J,0)= (rhs(I,J)+c11(i-1,j)*qs(I,J,0))*divi(I,J)
     qs(I+2,J,1)= (c11(i-1,j)*qs(I,J,1))*divi(I,J) 
     qs(I+2,J,2)= (c11(i-1,j)*qs(I,J,2))*divi(I,J) 
     qs(I+2,J,3)= (c11(i-1,j)*qs(I,J,3))*divi(I,J) 
     qs(I+2,J,4)= (c11(i-1,j)*qs(I,J,4))*divi(I,J) 
        enddo
      enddo

      DO J=1,DP_Depth
        DO I=2,N-1
       divi(I,J)=c11(i+1,j)+c11(i-1,j)*(rpe_1-ps(i,j)) +(Delta_t+rpe_1)
        divi(I,J)=rpe_1/divi(I,J)
     ps(i+2,j)=                  c11(i+1,j)*divi(I,J) 
     qs(I+2,J,0)= (rhs(I,J)+c11(i-1,j)*qs(I,J,0))*divi(I,J)
     qs(I+2,J,1)= (c11(i-1,j)*qs(I,J,1))*divi(I,J) 
     qs(I+2,J,2)= (c11(i-1,j)*qs(I,J,2))*divi(I,J) 
     qs(I+2,J,3)= (c11(i-1,j)*qs(I,J,3))*divi(I,J) 
     qs(I+2,J,4)= (c11(i-1,j)*qs(I,J,4))*divi(I,J) 
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=2,N-1
       divi(I,J)=c11(i+1,j)+c11(i-1,j)*(rpe_1-ps(i,j)) +(Delta_t+rpe_1)
        divi(I,J)=rpe_1/divi(I,J)
     ps(i+2,j)=                  c11(i+1,j)*divi(I,J) 
     qs(I+2,J,0)= (rhs(I,J)+c11(i-1,j)*qs(I,J,0))*divi(I,J)
     qs(I+2,J,1)= (c11(i-1,j)*qs(I,J,1))*divi(I,J) 
     qs(I+2,J,2)= (c11(i-1,j)*qs(I,J,2))*divi(I,J) 
     qs(I+2,J,3)= (c11(i-1,j)*qs(I,J,3))*divi(I,J) 
     qs(I+2,J,4)= (c11(i-1,j)*qs(I,J,4))*divi(I,J) 
        enddo
      enddo



  !write(*,*) 'q finished'
  ! calculate the 5 linear subsystems with boundary conditions ( right boundary)
  
  DO J=1,M  
  ws(N  ,J,0)=rpe_0
  ws(N+1,J,0)=rpe_0

  ws(N  ,J,1)=rpe_0
  ws(N+1,J,1)=rpe_0

  ws(N  ,J,2)=rpe_0
  ws(N+1,J,2)=rpe_0

  ws(N  ,J,3)=rpe_1
  ws(N+1,J,3)=rpe_0

  ws(N  ,J,4)=rpe_0
  ws(N+1,J,4)=rpe_1
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
   
      DO J=1+DP_Depth,M-DP_Depth
        DO I=N-1,2,-1
     ws(I,J,0) = ps(I+2,J)*ws(I+2,J,0)+qs(I+2,J,0)

     ws(I,J,1) = ps(I+2,J)*ws(I+2,J,1)+qs(I+2,J,1)
     ws(I,J,2) = ps(I+2,J)*ws(I+2,J,2)+qs(I+2,J,2)
     ws(I,J,3) = ps(I+2,J)*ws(I+2,J,3)+qs(I+2,J,3)
     ws(I,J,4) = ps(I+2,J)*ws(I+2,J,4)+qs(I+2,J,4)
        enddo
      enddo

      DO J=1,DP_Depth
        DO I=N-1,2,-1
     ws(I,J,0) = ps(I+2,J)*ws(I+2,J,0)+qs(I+2,J,0)

     ws(I,J,1) = ps(I+2,J)*ws(I+2,J,1)+qs(I+2,J,1)
     ws(I,J,2) = ps(I+2,J)*ws(I+2,J,2)+qs(I+2,J,2)
     ws(I,J,3) = ps(I+2,J)*ws(I+2,J,3)+qs(I+2,J,3)
     ws(I,J,4) = ps(I+2,J)*ws(I+2,J,4)+qs(I+2,J,4)
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=N-1,2,-1
     ws(I,J,0) = ps(I+2,J)*ws(I+2,J,0)+qs(I+2,J,0)

     ws(I,J,1) = ps(I+2,J)*ws(I+2,J,1)+qs(I+2,J,1)
     ws(I,J,2) = ps(I+2,J)*ws(I+2,J,2)+qs(I+2,J,2)
     ws(I,J,3) = ps(I+2,J)*ws(I+2,J,3)+qs(I+2,J,3)
     ws(I,J,4) = ps(I+2,J)*ws(I+2,J,4)+qs(I+2,J,4)
        enddo
      enddo

   !! end of rewrite

  ! write(*,*) 'w finished'
  ! solve the subsystems for the coefficients a,b,g,d for each latitude

  DO J=1,M
    if (J<=DP_depth .or. J>=M+1-DP_depth) then
      
    end if
      
       d11=ws(N-2,J,1)-rpe_1
       d12=ws(N-2,J,2)
       d13=ws(N-2,J,3)
       d14=ws(N-2,J,4)
       s1 =-ws(N-2,J,0)

       d21=ws(N-1,J,1)
       d22=ws(N-1,J,2)-rpe_1
       d23=ws(N-1,J,3)
       d24=ws(N-1,J,4)
       s2 =-ws(N-1,J,0)

       d31=ws(2,J,1)
       d32=ws(2,J,2)
       d33=ws(2,J,3)-rpe_1
       d34=ws(2,J,4)
       s3 =-ws(2,J,0)

       d41=ws(3,J,1)
       d42=ws(3,J,2)
       d43=ws(3,J,3)
       d44=ws(3,J,4)-rpe_1
       s4 =-ws(3,J,0)

      det40=d11*det3(d22,d23,d24,d32,d33,d34,d42,d43,d44) &
        &  -d21*det3(d12,d13,d14,d32,d33,d34,d42,d43,d44)  &
        &  +d31*det3(d12,d13,d14,d22,d23,d24,d42,d43,d44)  &
        &  -d41*det3(d12,d13,d14,d22,d23,d24,d32,d33,d34) 
      deti=rpe_1/det40
      det41=s1 *det3(d22,d23,d24,d32,d33,d34,d42,d43,d44) &
        &  -s2 *det3(d12,d13,d14,d32,d33,d34,d42,d43,d44)  &
        &  +s3 *det3(d12,d13,d14,d22,d23,d24,d42,d43,d44)  &
        &  -s4 *det3(d12,d13,d14,d22,d23,d24,d32,d33,d34)  
      det42=d11*det3( s2,d23,d24, s3,d33,d34, s4,d43,d44) &
        &  -d21*det3( s1,d13,d14, s3,d33,d34, s4,d43,d44)  &
        &  +d31*det3( s1,d13,d14, s2,d23,d24, s4,d43,d44)  &
        &  -d41*det3( s1,d13,d14, s2,d23,d24, s3,d33,d34)  
      det43=d11*det3(d22, s2,d24,d32, s3,d34,d42, s4,d44) &
        &  -d21*det3(d12, s1,d14,d32, s3,d34,d42, s4,d44)  &
        &  +d31*det3(d12, s1,d14,d22, s2,d24,d42, s4,d44)  &
        &  -d41*det3(d12, s1,d14,d22, s2,d24,d32, s3,d34)
      det44=d11*det3(d22,d23, s2,d32,d33, s3,d42,d43, s4) &
        &  -d21*det3(d12,d13, s1,d32,d33, s3,d42,d43, s4)  &
        &  +d31*det3(d12,d13, s1,d22,d23, s2,d42,d43, s4)  &
        &  -d41*det3(d12,d13, s1,d22,d23, s2,d32,d33, s3)
       
      aa(J,4)=det44*deti
      aa(J,3)=det43*deti
      aa(J,2)=det42*deti
      aa(J,1)=det41*deti

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
   
      DO J=1+DP_Depth,M-DP_Depth
        DO I=2,N-1
          QU(I,J)=ws(I,J,0) +aa(J,1)*ws(I,J,1) &
                    &   +aa(J,2)*ws(I,J,2)  &
                    &   +aa(J,3)*ws(I,J,3)  &
                    &   +aa(J,4)*ws(I,J,4)
        enddo
      enddo

      DO J=1,DP_Depth
        DO I=2,N-1
          QU(I,J)=ws(I,J,0) +aa(J,1)*ws(I,J,1) &
                    &   +aa(J,2)*ws(I,J,2)  &
                    &   +aa(J,3)*ws(I,J,3)  &
                    &   +aa(J,4)*ws(I,J,4)
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=2,N-1
          QU(I,J)=ws(I,J,0) +aa(J,1)*ws(I,J,1) &
                    &   +aa(J,2)*ws(I,J,2)  &
                    &   +aa(J,3)*ws(I,J,3)  &
                    &   +aa(J,4)*ws(I,J,4)
        enddo
      enddo

  CALL XBC(QU,N,M)
 
  call adjust_conservation(QU,s,S_full,n,m)
   !! end of rewrite
   
  !write(*,*) 'QU finished'






end do
  !write(*,*) 'BC QU finished'
END SUBROUTINE

SUBROUTINE precon_ADI_Piotr(R,QU , T_step,  A, ps, &
&  divi,A11,A12,A21,A22,B11,B22,P0,U,V,S,S_full,N,M,IP,ID_PREC, num_of_bits,DP_Depth)
use implicit_functions_SP

implicit none

REAL(Kind=4) :: R(N,M),QU(N,M),A11(N,M),A12(N,M),A21(N,M),A22(N,M),  &
    &      B11(N,M),B22(N,M),P0(N,M),U(N,M),V(N,M),S(N,M),  C11(N,M)
INTEGER :: IP(N), ID_PREC, num_of_bits,DP_Depth
INTEGER :: N, M

REAL(Kind=4) :: max_QX_QY, F(N,M), rhs(N,M), qs(0:N-1,M,0:2), ps(0:N-2,M), ws(N-1,M,0:4), A(N,M), B(N,M), C(N,M), divi(N,M)
REAL(Kind=4) :: aa(M,1:4), deti, det40, det41, det42, det43, det44, det3,   &
                           & d11, d12, d13, d14, d21, d22, d23, d24,       &
                           & d31, d32, d33, d34, d41, d42, d43, d44,       &
                           & s1, s2, s3, s4
REAL(Kind=4) ::T_step, Delta_t_I, dn, dni, S_full, util, vtil, swcp
!REAL(Kind=4) :: 
integer :: iter, max_iter, time_scale  !! number of richardson iterations
INTEGER :: I, J, iteration, i2, il1, il2



max_iter=2
swcp=rpe_1
!initialize the inverse of R with 0
DO J=1,M
  DO I=1,N
    QU(I,J) =rpe_0
    rhs(I,J)=rpe_0
  end do
end do

!! into precon_prep_depth
  DO J=1,M
    DO I=1,N    
 
      C11(I,J)=rpe_05*a11(i,j)
    END DO
  END DO

Delta_t_I=rpe_1/(rpe_025/T_step)
!write(*,*) Delta_t_I, T_step


      DO J=1,M

        ps(0,j)=rpe_0
        divi(1,J)=rpe_1/(c11(2,j)+c11(N-2,j) +s(1,j)*(Delta_t_I+rpe_1))
        ps(1,j) = c11(2,j)*divi(1,J)

      enddo

      DO J=1,M
        DO I=2,N-1
        divi(I,J)=rpe_1/(c11(i+1,j)+c11(i-1,j)*(rpe_1-ps(i-2,j))&
              &  +s(i,j)*(Delta_t_I+rpe_1))
        ps(i,j)=  c11(i+1,j)*divi(I,J)

        enddo
      enddo

!! end into precon_prep_depth
!Delta_t=1.0d0/T_step
!write(*,*) Delta_t




! 2) loop of following ADi iterations
do iteration=1,max_iter
 ! 2.1 calculate new right hand side using old QU and R to get tilde{tilde(R)}
If(iteration==1) then

! 1) first iteration


      DO J=1+DP_Depth,M-DP_Depth
        DO I=1,N
      rhs(I,J)=s(i,j)*( - R(I,J))
        enddo
      enddo

      DO J=1,DP_Depth
        DO I=1,N
      rhs(I,J)=s(i,j)*( - R(I,J))
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=1,N
      rhs(I,J)=s(i,j)*( - R(I,J))
        enddo
      enddo


else
  
      do j=2,m-1
       do i=2,n-1
        U(i,j)=QU(i+1,j)-QU(i-1,j)
        V(i,j)=QU(i,j+1)-QU(i,j-1)
       enddo
      enddo
      do i=2,n-1
       U(i,1)=QU(i+1,1)-QU(i-1,1)
       U(i,m)=QU(i+1,m)-QU(i-1,m)
       V(i,1)=QU(i,2)-QU(ip(i),1)
       V(i,m)=QU(ip(i),m)-QU(i,m-1)
      enddo

       CALL XBC(U,N,M)
       CALL XBC(V,N,M)

      do j=1,m
      do i=1,n
       util=                swcp*(a12(i,j)*V(i,j)+b11(i,j)*QU(i,j))
       vtil=V(i,j)*a21(i,j)+swcp*(a22(i,j)*U(i,j)+b22(i,j)*QU(i,j))
       U(i,j)=util
       V(i,j)=vtil
      enddo
      enddo

  CALL XBC(U,N,M)
  CALL XBC(V,N,M)

  DO J=2,M-1
    DO I=2,N-1
      F(I,J)= (U(I+1,J)-U(I-1,J)+V(I,J+1)-V(I,J-1))/(S(I,J))
    end do
  end do    
 

  DO I=2,N-1
    F(I,1)= (U(I+1,1)-U(I-1,1)+(V(I,2)+V(I,1)))/(S(I,1))
    F(I,M)= (U(I+1,M)-U(I-1,M)-(V(I,M)+V(I,M-1)))/(S(I,M))
  ENDDO


  CALL XBC(F,N,M)
      do j=1,m
      do i=1,n
       F(i,j)=rpe_05*F(i,j) !-p(i,j)
      enddo
      enddo

  DO J=1,M
    DO I=1,N    
 
      rhs(I,J)=F(I,J) + s(i,j)*(Delta_t_I*QU(I,J)    &
                          &   - R(I,J))
    END DO
  END DO
 CALL XBC(rhs,N,M)

end if

 
  DO J=1,M

    ps(0,J)    =rpe_0
    qs(0,J,0)  =rpe_0

    divi(1,J)=(c11(2,j)+c11(N-2,j) +s(1,j)*(Delta_t_I+rpe_1))
    divi(1,J)=rpe_1/    divi(1,J)

    ps(1,J)    =c11(2,j)*divi(1,J)
    qs(1,J,0)  =rhs(1,j)*divi(1,J)

    qs(0,J,1)  =rpe_1
    qs(1,J,1)  =rpe_0

    qs(0,J,2)  =rpe_0
    qs(1,J,2)  =c11(N-2,j)*divi(1,J)

  end do


  

      DO J=1+DP_Depth,M-DP_Depth
        DO I=2,N-2
        divi(I,J)=(c11(i+1,j)+c11(i-1,j)*(rpe_1-ps(i-2,j))&
              &  +s(i,j)*(Delta_t_I+rpe_1))
        divi(I,J)=rpe_1/divi(I,J)
        ps(i,j)=  c11(i+1,j)*divi(I,J)

       qs(I,J,0)= (rhs(I,J)+c11(i-1,j)*qs(I-2,J,0))*divi(I,J)
       qs(I,J,1)= (         c11(i-1,j)*qs(I-2,J,1))*divi(I,J) 
       qs(I,J,2)= (         c11(i-1,j)*qs(I-2,J,2))*divi(I,J) 

        enddo
      enddo

      DO J=1,DP_Depth
        DO I=2,N-2
        divi(I,J)=(c11(i+1,j)+c11(i-1,j)*(rpe_1-ps(i-2,j))&
              &  +s(i,j)*(Delta_t_I+rpe_1))
        divi(I,J)=rpe_1/divi(I,J)
        ps(i,j)=  c11(i+1,j)*divi(I,J)

       qs(I,J,0)= (rhs(I,J)+c11(i-1,j)*qs(I-2,J,0))*divi(I,J)
       qs(I,J,1)= (         c11(i-1,j)*qs(I-2,J,1))*divi(I,J) 
       qs(I,J,2)= (         c11(i-1,j)*qs(I-2,J,2))*divi(I,J) 
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=2,N-2
        divi(I,J)=(c11(i+1,j)+c11(i-1,j)*(rpe_1-ps(i-2,j))&
              &  +s(i,j)*(Delta_t_I+rpe_1))
        divi(I,J)=rpe_1/divi(I,J)
        ps(i,j)=  c11(i+1,j)*divi(I,J)

       qs(I,J,0)= (rhs(I,J)+c11(i-1,j)*qs(I-2,J,0))*divi(I,J)
       qs(I,J,1)= (         c11(i-1,j)*qs(I-2,J,1))*divi(I,J) 
       qs(I,J,2)= (         c11(i-1,j)*qs(I-2,J,2))*divi(I,J)  
        enddo
      enddo



  !write(*,*) 'q finished'
  ! calculate the 5 linear subsystems with boundary conditions ( right boundary)
  
  DO J=1,M  
  ws(N-1  ,J,0)=rpe_0
  ws(N-2,J,0)=qs(N-2,J,0)

  ws(N-1,J,1)=rpe_0
  ws(N-2,J,1)=qs(N-2,J,1)

  ws(N-1,J,2)=rpe_1
  ws(N-2,J,2)=rpe_0

  ws(N-1,J,3)=rpe_0
  ws(N-2,J,3)=qs(N-2,J,2)

  ws(N-1,J,4)=rpe_0
  ws(N-2,J,4)=ps(N-2,j)
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
   
      DO J=1+DP_Depth,M-DP_Depth
        DO I=N-3,1,-1
     ws(I,J,0) = ps(I,J)*ws(I+2,J,0)+qs(I,J,0)

     ws(I,J,1) = ps(I,J)*ws(I+2,J,1)+qs(I,J,1)
     ws(I,J,2) = ps(I,J)*ws(I+2,J,2)
     ws(I,J,3) = ps(I,J)*ws(I+2,J,3)+qs(I,J,2)
     ws(I,J,4) = ps(I,J)*ws(I+2,J,4)
        enddo
      enddo

      DO J=1,DP_Depth
        DO I=N-3,1,-1
     ws(I,J,0) = ps(I,J)*ws(I+2,J,0)+qs(I,J,0)

     ws(I,J,1) = ps(I,J)*ws(I+2,J,1)+qs(I,J,1)
     ws(I,J,2) = ps(I,J)*ws(I+2,J,2)
     ws(I,J,3) = ps(I,J)*ws(I+2,J,3)+qs(I,J,2)
     ws(I,J,4) = ps(I,J)*ws(I+2,J,4)
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=N-3,1,-1
     ws(I,J,0) = ps(I,J)*ws(I+2,J,0)+qs(I,J,0)

     ws(I,J,1) = ps(I,J)*ws(I+2,J,1)+qs(I,J,1)
     ws(I,J,2) = ps(I,J)*ws(I+2,J,2)
     ws(I,J,3) = ps(I,J)*ws(I+2,J,3)+qs(I,J,2)
     ws(I,J,4) = ps(I,J)*ws(I+2,J,4)
        enddo
      enddo

   !! end of rewrite

  ! write(*,*) 'w finished'
  ! solve the subsystems for the coefficients a,b,g,d for each latitude
  il1=N-2
  il2=N-3
  i2 =2
  
  DO J=1,M
    if (J<=DP_depth .or. J>=M+1-DP_depth) then
      
    end if
      
       d11=ws(il1,J,1)-rpe_1
       d12=ws(il1,J,2)
       d13=ws(il1,J,3)
       d14=ws(il1,J,4)
       s1 =-ws(il1,J,0)

       d21=ws(1,J,1)
       d22=ws(1,J,2)-rpe_1
       d23=ws(1,J,3)
       d24=ws(1,J,4)
       s2 =-ws(1,J,0)

       d31=ws(il2,J,1)
       d32=ws(il2,J,2)
       d33=ws(il2,J,3)-rpe_1
       d34=ws(il2,J,4)
       s3 =-ws(il2,J,0)

       d41=ws(i2,J,1)
       d42=ws(i2,J,2)
       d43=ws(i2,J,3)
       d44=ws(i2,J,4)-rpe_1
       s4 =-ws(i2,J,0)

      det40=d11*det3(d22,d23,d24,d32,d33,d34,d42,d43,d44) &
        &  -d21*det3(d12,d13,d14,d32,d33,d34,d42,d43,d44)  &
        &  +d31*det3(d12,d13,d14,d22,d23,d24,d42,d43,d44)  &
        &  -d41*det3(d12,d13,d14,d22,d23,d24,d32,d33,d34) 
      deti=rpe_1/det40
      det41=s1 *det3(d22,d23,d24,d32,d33,d34,d42,d43,d44) &
        &  -s2 *det3(d12,d13,d14,d32,d33,d34,d42,d43,d44)  &
        &  +s3 *det3(d12,d13,d14,d22,d23,d24,d42,d43,d44)  &
        &  -s4 *det3(d12,d13,d14,d22,d23,d24,d32,d33,d34)  
      det42=d11*det3( s2,d23,d24, s3,d33,d34, s4,d43,d44) &
        &  -d21*det3( s1,d13,d14, s3,d33,d34, s4,d43,d44)  &
        &  +d31*det3( s1,d13,d14, s2,d23,d24, s4,d43,d44)  &
        &  -d41*det3( s1,d13,d14, s2,d23,d24, s3,d33,d34)  
      det43=d11*det3(d22, s2,d24,d32, s3,d34,d42, s4,d44) &
        &  -d21*det3(d12, s1,d14,d32, s3,d34,d42, s4,d44)  &
        &  +d31*det3(d12, s1,d14,d22, s2,d24,d42, s4,d44)  &
        &  -d41*det3(d12, s1,d14,d22, s2,d24,d32, s3,d34)
      det44=d11*det3(d22,d23, s2,d32,d33, s3,d42,d43, s4) &
        &  -d21*det3(d12,d13, s1,d32,d33, s3,d42,d43, s4)  &
        &  +d31*det3(d12,d13, s1,d22,d23, s2,d42,d43, s4)  &
        &  -d41*det3(d12,d13, s1,d22,d23, s2,d32,d33, s3)
       
      aa(J,4)=det44*deti
      aa(J,3)=det43*deti
      aa(J,2)=det42*deti
      aa(J,1)=det41*deti

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
   
      DO J=1+DP_Depth,M-DP_Depth
        DO I=1,N-1
          QU(I,J)=ws(I,J,0) +aa(J,1)*ws(I,J,1) &
                    &   +aa(J,2)*ws(I,J,2)  &
                    &   +aa(J,3)*ws(I,J,3)  &
                    &   +aa(J,4)*ws(I,J,4)
        enddo
      enddo

      DO J=1,DP_Depth
        DO I=1,N-1
          QU(I,J)=ws(I,J,0) +aa(J,1)*ws(I,J,1) &
                    &   +aa(J,2)*ws(I,J,2)  &
                    &   +aa(J,3)*ws(I,J,3)  &
                    &   +aa(J,4)*ws(I,J,4)
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=1,N-1
          QU(I,J)=ws(I,J,0) +aa(J,1)*ws(I,J,1) &
                    &   +aa(J,2)*ws(I,J,2)  &
                    &   +aa(J,3)*ws(I,J,3)  &
                    &   +aa(J,4)*ws(I,J,4)
        enddo
      enddo

  CALL XBC(QU,N,M)
 
  call adjust_conservation(QU,s,S_full,n,m)
   !! end of rewrite
   
  !write(*,*) 'QU finished'






end do
  !write(*,*) 'BC QU finished'
END SUBROUTINE


function det3(r11,r12,r13,r21,r22,r23,r31,r32,r33) 


implicit none
REAL(Kind=4) ::  r11,r12,r13,r21,r22,r23,r31,r32,r33, det3

        det3=    r11*r22*r33+r12*r23*r31+r13*r21*r32    &
            &     -r31*r22*r13-r32*r23*r11-r33*r21*r12

end function det3


SUBROUTINE precon_Jac(R,QU , T_step,  A, ps, divi,A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP,ID_PREC, num_of_bits,DP_Depth)
use implicit_functions_SP


implicit none

REAL(Kind=4) ::  R(N,M),QU(N,M),A11(N,M),A12(N,M),A21(N,M),A22(N,M),  &
    &      B11(N,M),B22(N,M),P0(N,M),U(N,M),V(N,M),S(N,M)
INTEGER :: IP(N), ID_PREC, num_of_bits,DP_Depth
INTEGER :: N, M

REAL(Kind=4) ::   ps(N+1,M), A(N,M), divi(N,M), rhs(N,M)

REAL(Kind=4) :: T_step, swc, betap
!REAL(Kind=4) :: 

INTEGER :: I, J, it, itr, line

itr = 0
line= 0
swc = 1.0d0


betap=2.0d0/3.0d0

do j=1,M
  do i=1,N
    QU(i,j)=-betap*R(i,j)/A(i,j)
    rhs(i,j)=0.
  enddo
enddo

call adjust_conservation(QU,s,T_step,n,m)
 ! write(*,*) 'linop in'
!call linop(QU,rhs,a11,a12,a21,a22,b11,b22,A,s,ip,swc,n,m)
  !write(*,*) 'linop out'

do it=1,itr
  !write(*,*) 'iteration', it
  do j=1,M
    do i=1,N
      QU(i,j)=betap*(rhs(i,j)-R(i,j))/A(i,j) + (1.0d0-betap)*QU(i,j)
    enddo
  enddo
  call adjust_conservation(QU,s,T_step,n,m)
     do j=1,m
      do i=1,n
        write(*,*)  'QU',i, J, QU(i,j)
        write(*,*)  'rhs',i, J, rhs(i,j)-R(i,j)
       read(*,*)
       enddo
end do
  if(it.lt.itr) call linop(QU,rhs,a11,a12,a21,a22,b11,b22,A,s,ip,swc,n,m)
enddo
!     do j=1,m
!      do i=1,n
!        write(*,*)  'QU',i, J, QU(i,j)
!        write(*,*)  'rhs',i, J, rhs(i,j)
!       read(*,*)
!       enddo
!end do
END SUBROUTINE

subroutine linop(p,r,a11,a12,a21,a22,b11,b22, &
     &               A,s,ip,swc,N,M)
use implicit_functions_SP

implicit none

REAL(Kind=4) ::  p(n,m),r(n,m), &
     &          a11(n,m),a12(n,m),a21(n,m),a22(n,m),b11(n,m),b22(n,m), &
     &          A(n,m),s(n,m)
integer :: ip(n), N, M
REAL(Kind=4) :: swc

REAL(Kind=4) ::  pfx(n,m), pfy(n,m), util, vtil
INTEGER :: I, J


      do j=2,m-1
       do i=2,n-1
        pfx(i,j)=p(i+1,j)-p(i-1,j)
        pfy(i,j)=p(i,j+1)-p(i,j-1)
       enddo
      enddo
      do i=2,n-1
       pfx(i,1)=p(i+1,1)-p(i-1,1)
       pfx(i,m)=p(i+1,m)-p(i-1,m)
       pfy(i,1)=p(i,2)-p(ip(i),1)
       pfy(i,m)=p(ip(i),m)-p(i,m-1)
      enddo
      do j=1,m
       pfx(1,j)=pfx(n-1,j)
       pfx(n,j)=pfx(2  ,j)
       pfy(1,j)=pfy(n-1,j)
       pfy(n,j)=pfy(2  ,j)
      enddo

      !do j=1,m
      !do i=1,n
      ! util=pfx(i,j)*a11(i,j)+swc*(a12(i,j)*pfy(i,j)+b11(i,j)*p(i,j))
      ! vtil=pfy(i,j)*a21(i,j)+swc*(a22(i,j)*pfx(i,j)+b22(i,j)*p(i,j))
      ! pfx(i,j)=util
      ! pfy(i,j)=vtil
      !enddo
      !enddo

      do j=1,m
      do i=1,n
       util=swc*(a12(i,j)*pfy(i,j)+b11(i,j)*p(i,j))
       vtil=swc*(a22(i,j)*pfx(i,j)+b22(i,j)*p(i,j))
       pfx(i,j)=util
       pfy(i,j)=vtil
      enddo
      enddo
    write(*,*) 'linop, middle'
      do j=2,m-1
       do i=2,n-1
        r(i,j)= ( pfx(i+1,j)   +a11(i+1,j)*p(i+2,j)       &
             &   - pfx(i-1,j)  +a11(i-1,j)*p(i-2,j)      &
             &   + pfy(i,j+1)  +a21(i,j+1)*p(i,j+2)    &
             &   - pfy(i,j-1)  +a21(i,j-1)*p(i,j-2)  )       &
             &   / s(i,j)
        write(*,*) r(i,j)
       enddo
      enddo
    write(*,*) 'linop, finish'
      do i=2,n-1
        ! r(i,1)= (pfx(i+1,1)-pfx(i-1,1)+(pfy(i,2)+pfy(i,1)))/s(i,1)
        r(i,1)= ( pfx(i+1,1)   +a11(i+1,1)*p(i+2,1)       &
             &   - pfx(i-1,1)  +a11(i-1,1)*p(i-2,1)      &
             &   + pfy(i,2)  +a21(i,2)*p(i,3)    &
             &   + pfy(i,1)  -a21(i,1)*p(ip(i),1)  )       &
             &   / s(i,1)


       !r(i,m)= (pfx(i+1,m)-pfx(i-1,m)-(pfy(i,m)+pfy(i,m-1)))/s(i,m)
        r(i,m)= ( pfx(i+1,m)   +a11(i+1,m)*p(i+2,m)       &
             &   - pfx(i-1,m)  +a11(i-1,m)*p(i-2,m)      &
             &   - pfy(i,m)  -a21(i,m)*p(ip(i),m)    &
             &   - pfy(i,m-1)  +a21(i,m-1)*p(i,m-2)  )       &
             &   / s(i,j)

      enddo
      do j=1,m
       r(1,j)=r(n-1,j)
       r(n,j)=r(2  ,j)
      enddo
      do j=1,m
      do i=1,n
       r(i,j)=0.5*r(i,j)-p(i,j)
      enddo
      enddo

end subroutine

subroutine adjust_conservation(p,s,S_full,n,m)
use implicit_functions_SP

implicit none
REAL(Kind=4) :: p(n,m),s(n,m)
REAL(Kind=4) :: S_full, cnst
INTEGER :: N, M

INTEGER :: I, J
    !write(*,*) S_full
      cnst=0.
      do j=1,m
        do i=2,n-1
          !write(*,*) i, j, s(i,j), p(i,j), s(i,j)*p(i,j)
          !read(*,*)
          cnst=cnst+s(i,j)*p(i,j)
        enddo
      enddo
      cnst=cnst/S_full
     !write(*,*) cnst

      do j=1,m
        do i=1,n
          p(i,j)=p(i,j)-cnst
        enddo
      enddo

end subroutine



subroutine diagoc(T_step,A_c,a11,a22,s, N, M, DP_depth)
use implicit_functions_SP

implicit none
      REAL(Kind=4) ::A_c(N,M),a11(N,M),a22(N,M),s(N,M)
      REAL(Kind=4) :: T_step
      integer :: itr,itmn, DP_depth, N, M

INTEGER :: I,J
T_step=0.0d0

do j=1,M
  do i=2,N-1
    T_step=T_step+s(i,j)
  enddo
enddo

do j=1,M
  do i=2,N-1
    A_c(i,j)=a11(i+1,j)+a11(i-1,j)
  enddo

  A_c(1,j)=A_c(n-1,j)  !a11(2,j)+a11(n-2,j)
  A_c(n,j)=A_c( 2 ,j)  !a11(3,j)+a11(n-1,j)
enddo

do j=2,M-1
  do i=1,N
    A_c(i,j)=A_c(i,j)+a22(i,j+1)+a22(i,j-1)
  enddo
enddo

do i=1,N
   A_c(i,1)=A_c(i,1)+a22(i,2)
   A_c(i,m)=A_c(i,m)+a22(i,m-1)
enddo

do j=1,M
  do i=1,N
    A_c(i,j)=0.5*A_c(i,j)/s(i,j)+1.0d0
  enddo
enddo

!   do i=1,n
!        write(*,*)  (A_c(i,j),  j=1,m )
!       read(*,*)
!       enddo


end subroutine

SUBROUTINE RHWT(U0,V0,PT0,PD0,P0,COR,X,Y,N,M,F0,A, TIME)
use implicit_functions_SP

implicit none
REAL(Kind=4) :: PT0(N,M), PD0(N,M), P0(N,M)

REAL(Kind=4) :: U0(N,M),V0(N,M)
REAL(Kind=4) :: COR(N,M),X(N),Y(M),F0, A, XX(N), TIME
INTEGER :: N, M


REAL(Kind=4) ::     ATH(M), BTH(M), CTH(M), TH
REAL(Kind=4) ::    OM,K,PH0

REAL(Kind=4) :: Grav, GRI, PI2, VNIU
INTEGER :: R, I, J

OM=7.848E-6
K=7.848E-6
R=4
PH0= 78.4E3
Grav =9.80616

DO J=1,M
  TH=Y(J)
  ATH(J)=OM*0.5d0*(F0+OM)*(COS(TH))**2                          &
     &  +0.25d0*K**2*(COS(TH))**(2*R)*( (R+1)*(COS(TH))**2      &
     &   +FLOAT(2*R**2-R-2)-2.0d0*R**2/(COS(TH))**2 )
  BTH(J)=(F0+2.*OM)*K/FLOAT((R+1)*(R+2))*(COS(TH))**R         &
     &       *( FLOAT(R**2+2*R+2)-((R+1)*COS(TH))**2 )
  CTH(J)=0.25d0*K**2*(COS(TH))**(2*R)*( FLOAT(R+1)*(COS(TH))**2 &
     &       -FLOAT(R+2) )  
end do

!! new location because RHW gets advective with velocity VNIU
      GRI=1./Grav
      PI2=2.*ACOS(-1.)
      VNIU = (Float(R*(3+R))*OM-F0)/Float((1+R)*(2+R))
      IF(TIME.EQ.0.) write (*,*)  'vniu:',VNIU

      DO I=1,N
       XX(I)=X(I)-VNIU*TIME
       IF(XX(I).LT.0.) XX(I)=XX(I)+PI2
       IF(XX(I).GT.PI2) XX(I)=XX(I)-PI2
      ENDDO
       XX(1)=XX(N-1)
       XX(N)=XX(2)
      !IF(TIME.EQ.0.) write (*,*)  &
      !&  'x(1),x(2),x(n-1),x(n)',x(1),x(2),x(n-1),x(n)


DO J=1,M
  DO I=1,N
      U0(I,J)=A*OM*COS(Y(J))+A*K*COS(R*XX(I))                      &
       &   *(COS(Y(J)))**(R-1)*(R*(SIN(Y(J)))**2-(COS(Y(J)))**2)
      V0(I,J)=-A*K*R*(COS(Y(J)))**(R-1)*SIN(Y(J))*SIN(R*XX(I))
      PT0(I,J)=PH0+A**2*ATH(J)+A**2*BTH(J)*COS(R*XX(I))      &
       &   +A**2*CTH(J)*COS(2.0d0*R*XX(I))
      PD0(I,J)= max(0.0d0, PT0(I,J)*GRI-P0(I,J))
  end do
end do
  !write (*,*)  'initrhw absorber called'

END SUBROUTINE

SUBROUTINE PRFORC_ABS( P,F1,F2,PB,P0,E1,E2,HX,HY,COR,       &
     &            N,M,IP,GC1,GC2,Alp_REL, NOR,IRS)
use implicit_functions_SP

implicit none
REAL(Kind=4) :: P(N,M),F1(N,M),F2(N,M),PB(N,M),P0(N,M),E1(N,M),E2(N,M), &
     & HX(N,M),HY(N,M),COR(N,M)
REAL(Kind=4) :: GC1, GC2
REAL(Kind=4) :: Alp_REL(N,M)

INTEGER  :: IP(N)
INTEGER  :: N, M, NOR, IRS


INTEGER :: I, J, NM

REAL(Kind=4) :: GH1, GH2, UTILD, VTILD, GMM, AMM, DETI, A, B

GH1=rpe_05*GC1
GH2=rpe_05*GC2

DO J=2,M-1
  DO I=2,N-1
    F1(I,J)=-GH1*(P(I+1,J)-P(I-1,J))/HX(I,J)
    F2(I,J)=-GH2*(P(I,J+1)-P(I,J-1))/HY(I,J)
  end do
end do
   
DO I=2,N-1
  F1(I,1)=-GH1*(P(I+1,1)-P(I-1,1))/HX(I,1)
  F1(I,M)=-GH1*(P(I+1,M)-P(I-1,M))/HX(I,M)
  F2(I,1)=-GH2*(P(I,2)-P(IP(I),1))/HY(I,1)
  F2(I,M)=-GH2*(P(IP(I),M)-P(I,M-1))/HY(I,M)
end do

CALL XBC(F1,N,M)
CALL XBC(F2,N,M)

IF(NOR.EQ.1) THEN
!  NM=N*M
  
!  DO I=1,NM
!    UTILD=F1(I,1)*(P0(I,1)-PB(I,1))+(P(I,1)-IRS*P0(I,1))*E1(I,1)
!    VTILD=F2(I,1)*(P0(I,1)-PB(I,1))+(P(I,1)-IRS*P0(I,1))*E2(I,1)
!    GMM=.5*COR(I,1)
!    F1(I,1)=(UTILD+GMM*VTILD)/(1.+GMM**2)*GH1
!    F2(I,1)=(VTILD-GMM*UTILD)/(1.+GMM**2)*GH2
!  end do
  DO J=1,M
    DO I=1,N
      UTILD=F1(I,J)*(P0(I,J)-PB(I,J))+(P(I,J)-IRS*P0(I,J))*E1(I,J)
      VTILD=F2(I,J)*(P0(I,J)-PB(I,J))+(P(I,J)-IRS*P0(I,J))*E2(I,J)
      
      AMM=1.+.5*Alp_REL(I,J)
      GMM=.5*COR(I,J)
      DETI=1./(AMM**2+GMM**2)
      
      A=AMM*DETI
      B=GMM*DETI

      F1(I,J)=(A*UTILD+B*VTILD)*GH1
      F2(I,J)=(A*VTILD-B*UTILD)*GH2
    end do
  end do
CALL XBC(F1,N,M)
CALL XBC(F2,N,M)
ENDIF

END SUBROUTINE


SUBROUTINE POLARABS(ALP,atau,Y,DT,N,M,IRHW)
use implicit_functions_SP

implicit none

REAL(Kind=4) :: ALP(N,M)
 
REAL(Kind=4) :: Y(M), DT
INTEGER :: N, M, IRHW

REAL(Kind=4) :: eps, pi, abswidth, atau, alpha, &
       & ymax, ymin, absy0_north, absy0_south, y0, y1, y2
INTEGER :: I, J

 DO J=1,M
 DO I=1,N
  ALP(I,J)=DT*0.
 ENDDO
 ENDDO

eps=1.e-10
pi=acos(-1.)
abswidth=pi/64.*3   !RHW4
!     atau=2.*(9.*2.*DT)
if (IRHW==1) then
      atau=200.*DT ! RHW4
elseif (IRHW==3) then
     atau=2.*DT    !Zonal flow past Earth orography
endif

alpha=1./atau

ymax = 0.5*pi
ymin =-0.5*pi
absy0_north = (ymax-abswidth)
absy0_south = (ymin+abswidth)

do j=1,m
  y0 = y(j)
  y1 = 0.
  y2 = 0.

  if(abswidth > 0. .and. atau < 1.e10) then
    y1 = max(0.0d0,-y0+absy0_south)/abswidth
    y2 = max(0.0d0, y0-absy0_north)/abswidth

!       if (y0 > 0.) y1 = exp(-(ymax-y0)/abswidth)
!       if (y0 < 0.) y2 = exp(-(y0-ymin)/abswidth)
    do i=1,n
      alp(i,j) = DT*alpha*((y1**2+y2**2)/(y1+y2+eps))
    enddo
   

  endif
enddo


end subroutine

SUBROUTINE SMOOTHSTATE(FL,FLS,SCR1,SCR2,IP,FLIP,KT,N,M)
implicit none
REAL(Kind=4) :: FL(N,M),FLS(N,M),SCR1(N,M),SCR2(N,M)
REAL(Kind=4) :: FLIP
INTEGER :: N,M, KT, IP(N)

integer :: i,j
integer :: IFLG, KITER, kit

      IFLG=1
      KITER=1

      IF(IFLG.EQ.0) THEN
      do kit=1,kiter
      do j=1,m
        do i=2,n-1
         scr1(i,j)=0.25*(fl(i+1,j)+2.*fl(i,j)+fl(i-1,j))
        enddo
         scr1(1,j)=scr1(n-1,j)
         scr1(n,j)=scr1(2,j)
      enddo
      do j=2,m-1
      do i=1,n
        scr2(i,j)=0.25*(scr1(i,j+1)+2.*scr1(i,j)+scr1(i,j-1))
      enddo
      enddo
      do i=1,n
        scr2(i,1)=0.25*(scr1(i,2)+2.*scr1(i,1)+flip*scr1(ip(i),1))
        scr2(i,m)=0.25*(flip*scr1(ip(i),m)+2.*scr1(i,m)+scr1(i,m-1))
      enddo
      enddo
      DO J=1,M
      DO I=1,N
       FLS(I,J)=SCR2(I,J)
      ENDDO
      ENDDO
      ENDIF

      IF(IFLG.EQ.1) THEN
       DO J=1,M
       DO I=1,N
        FLS(I,J)=(FLS(I,J)*KT+SCR2(I,J))/(KT+1)
       ENDDO
       ENDDO
      ENDIF

     END subroutine

     subroutine Earthtopo(h0,x,y,n,m)
implicit none
      REAL(Kind=4) :: h0(n,m),x(n),y(m)
      REAL(Kind=4) :: tp0(512,256),flon,flat, pi
      integer :: n,m, jj, ii, i ,j

      REAL(Kind=4) :: h0mx, h0mn, cycmx1, cycmxn

      pi = acos(-1.)

      do jj=1,m
        do ii=2,n-1
          read(33,*) i,j,flon,flat,tp0(i,j)

!     print*, i,j,flon,flat,tp0(i,j)
          h0(i+1,j)=tp0(i,j)
        enddo
      enddo
      call xbc(h0,n,m)

      h0mx=-1.e15
      h0mn= 1.e15
      do j=1,m
        do i=1,n
         h0mx=max(h0mx,h0(i,j))
         h0mn=min(h0mn,h0(i,j))
        enddo
      enddo
      print*, h0mx,h0mn

      cycmx1=-1.e15
      cycmxn=-1.e15
      do j=1,m
        cycmx1=max(cycmx1,abs(h0(1,j)-h0(n-1,j)))
        cycmxn=max(cycmxn,abs(h0(2,j)-h0(n,j)))
      enddo
      print*, cycmx1, cycmxn

      end subroutine

      SUBROUTINE SMOOTHTOP(FL,SCR,IP,N,M)
      implicit none
      REAL(Kind=4) :: FL(N,M),SCR(N,M)
      INTEGER :: IP(N), N, M


      REAL(Kind=4) :: h0mx, h0mn
      INTEGER :: kiter, kit, I, J

      KITER=1

      DO KIT=1,KITER
        do j=1,m
         do i=2,n-1
          scr(i,j)=0.25*(fl(i+1,j)+2.*fl(i,j)+fl(i-1,j))
         enddo
          scr(1,j)=scr(n-1,j)
          scr(n,j)=scr(2,j)
        enddo
        do j=2,m-1
         do i=1,n
          fl(i,j)=0.25*(scr(i,j+1)+2.*scr(i,j)+scr(i,j-1))
         enddo
        enddo
        do i=1,n
         fl(i,1)=0.25*(scr(i,2)+2.*scr(i,1)+scr(ip(i),1))
         fl(i,m)=0.25*(scr(ip(i),m)+2.*scr(i,m)+scr(i,m-1))
        enddo
      ENDDO

      h0mx=-1.e15
      h0mn= 1.e15
      do j=1,m
       do i=1,n
        h0mx=max(h0mx,fl(i,j))
        h0mn=min(h0mn,fl(i,j))
       enddo
      enddo
      print*, h0mx,h0mn
      
      h0mx=-1.e15
      h0mn= 1.e15
      do j=1,m
       do i=1,n
        fl(i,j)=1.289678*fl(i,j)
        h0mx=max(h0mx,fl(i,j))
        h0mn=min(h0mn,fl(i,j))
       enddo
      enddo
      print*, h0mx,h0mn

      END subroutine

SUBROUTINE PRFORC_ABS_dp( P,F1,F2,PB,P0,E1,E2,HX,HY,COR,       &
     &            N,M,IP,GC1,GC2,Alp_REL, NOR,IRS)

implicit none
double precision :: P(N,M),F1(N,M),F2(N,M),PB(N,M),P0(N,M),E1(N,M),E2(N,M), &
     & HX(N,M),HY(N,M),COR(N,M)
double precision :: GC1, GC2
double precision ::Alp_REL(N,M)

INTEGER  :: IP(N)
INTEGER  :: N, M, NOR, IRS


INTEGER :: I, J, NM

double precision :: GH1, GH2, UTILD, VTILD, GMM, AMM, DETI, A, B
GH1=0.5d0*GC1
GH2=0.5d0*GC2

DO J=2,M-1
  DO I=2,N-1
    F1(I,J)=-GH1*(P(I+1,J)-P(I-1,J))/HX(I,J)
    F2(I,J)=-GH2*(P(I,J+1)-P(I,J-1))/HY(I,J)
  end do
end do

DO I=2,N-1
  F1(I,1)=-GH1*(P(I+1,1)-P(I-1,1))/HX(I,1)
  F1(I,M)=-GH1*(P(I+1,M)-P(I-1,M))/HX(I,M)
  F2(I,1)=-GH2*(P(I,2)-P(IP(I),1))/HY(I,1)
  F2(I,M)=-GH2*(P(IP(I),M)-P(I,M-1))/HY(I,M)
end do

CALL XBC_dp(F1,N,M)
CALL XBC_dp(F2,N,M)

IF(NOR.EQ.1) THEN
!  NM=N*M

!  DO I=1,NM
!    UTILD=F1(I,1)*(P0(I,1)-PB(I,1))+(P(I,1)-IRS*P0(I,1))*E1(I,1)
!    VTILD=F2(I,1)*(P0(I,1)-PB(I,1))+(P(I,1)-IRS*P0(I,1))*E2(I,1)
!    GMM=.5*COR(I,1)
!    F1(I,1)=(UTILD+GMM*VTILD)/(1.+GMM**2)*GH1
!    F2(I,1)=(VTILD-GMM*UTILD)/(1.+GMM**2)*GH2
!  end do
  DO J=1,M
    DO I=1,N
      UTILD=F1(I,J)*(P0(I,J)-PB(I,J))+(P(I,J)-IRS*P0(I,J))*E1(I,J)
      VTILD=F2(I,J)*(P0(I,J)-PB(I,J))+(P(I,J)-IRS*P0(I,J))*E2(I,J)

      AMM=1.+.5*Alp_REL(I,J)
      GMM=.5*COR(I,J)
      DETI=1./(AMM**2+GMM**2)

      A=AMM*DETI
      B=GMM*DETI

      F1(I,J)=(A*UTILD+B*VTILD)*GH1
      F2(I,J)=(A*VTILD-B*UTILD)*GH2
    end do
  end do
CALL XBC_dp(F1,N,M)
CALL XBC_dp(F2,N,M)
ENDIF

END SUBROUTINE
SUBROUTINE DIVER_dp(F, U,  V, HX,HY,S, N, M,IP,IFLG)

 implicit none

double precision :: F(N,M),U(N,M),V(N,M),HX(N,M),HY(N,M),S(N,M)
INTEGER :: IP(N)
INTEGER :: N, M, IFLG


INTEGER:: I, J

DO J=2,M-1
  DO I=2,N-1
    F(I,J)= HY(I+1,J)*U(I+1,J)-HY(I-1,J)*U(I-1,J)        &
        &  +HX(I,J+1)*V(I,J+1)-HX(I,J-1)*V(I,J-1)
  end do
end do

! original
DO I=2,N-1
  F(I,1)= HY(I+1,1)*U(I+1,1)-HY(I-1,1)*U(I-1,1) &
      &  +(HX(I,2)*V(I,2)+HX(I,1)*V(I,1))
  F(I,M)= HY(I+1,M)*U(I+1,M)-HY(I-1,M)*U(I-1,M)  &
      &  -(HX(I,M)*V(I,M)+HX(I,M-1)*V(I,M-1))
end do
! end original



DO J=1,M
  DO I=2,N-1
    F(I,J)=0.5d0*IFLG*F(I,J)/S(I,J)
  end do
end do

DO J=1,M
  F(1,J)=F(N-1,J)
  F(N,J)=F(2  ,J)
end do

END SUBROUTINE
SUBROUTINE XBC_dp(X,N,M)


implicit none

INTEGER :: N, M
double precision :: X(N,M)

INTEGER :: J

DO J=1,M
  X(1,J)=X(N-1,J)
  X(N,J)=X(2,J)
end do

END subroutine

