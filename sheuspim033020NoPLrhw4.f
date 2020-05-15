      PROGRAM IS1SPSL
c     PARAMETER(NLON=130,NLAT=64)
c     PARAMETER(NLON=258,NLAT=128)
      PARAMETER(NLON=514,NLAT=256)
      PARAMETER(N=NLON,M=NLAT,NM=N*M)
      REAL U(N,M,0:1),V(N,M,0:1),PD(N,M),
     1     QX(N,M),QY(N,M),PT(N,M),
     2     HX(N,M),HY(N,M),DHX2Y(N,M),S(N,M),
     3     F1(N,M,0:1),F2(N,M,0:1),E1(N,M,-1:0),E2(N,M,-1:0),
     4     UA(N,M),VA(N,M+1),PC(N,M),FZ(N),SCR(N),ALP(N,M),
     5     QXS(N,M),QYS(N,M)
      COMMON/INITCON/ U0(N,M),V0(N,M),PT0(N,M),PD0(N,M),P0(N,M)
      DIMENSION X(N),XX(N),Y(M),COR(N,M)
      COMMON/IPIND/ IP(N)
      COMMON/INITVLS/ PTTRL2,PTTRLI,VLTRL2,VLTRLI,PDTRL2,PDTRLI
      DATA PTTRL2,PTTRLI/0.,-1.e-15/
      DATA VLTRL2,VLTRLI/0.,-1.e-15/
      DATA PDTRL2,PDTRLI/0.,-1.e-15/


CORIOLIS, GRAVITY AND EARTH RADIUS SPECIFICATION
      PARAMETER(ICORIO=1)
      DATA F0/1.4584E-4/
      DATA G,R/9.80616,6371.22E+03/

CHARACTERISTICS OF THE FLOW
c     DATA USCAL,H00,HMTS/20.,8.E3,2.E3/
      DATA USCAL,H00,HMTS/20.,8.E3,0.E3/
c 1 day = 2160*40s
c     DATA DT,NT,NPRINT,NPLOT/800,1594,99,398/   !max stable for 128x64
c     DATA DT,NT,NPRINT,NPLOT/400,3188,199,796/  !equivalent for 256x128
c     DATA DT,NT,NPRINT,NPLOT/200,6376,797,1594/ !equivalent for 512x256
      DATA DT,NT,NPRINT,NPLOT/200,6376,10,1594/ !equivalent for 512x256
c mountain induced Rossby wave
c     DATA DT,NT,NPRINT/400,3240,216/
c     DATA DT,NT,NPRINT/200,6480,432/  !flow over the pole
      COMMON/DISPL/ TIME,USCAL,DX,DY,H00,HMTS,DT
      COMMON/MPFIL/ LINER,MPFL
      DATA MPFL/999999/
c     DATA MPFL/1/
      INTEGER LADVEC
      DATA LADVEC/0/
      INTEGER IPRINT,IPLOT,KT0
      DATA IPRINT,IPLOT,KT0/0,0,0/

CONTROL TESTS: ZONAL FLOW OR ROSSBY WAVE
      PARAMETER(IRHW=1)
CONTROL EVALUATION OF THE PRESSURE GRADIENTS: IPS=0 CENTERED DIFFERNCING
CHOOSING IPS=1 GIVES WEIGHTED AVERAGE OF ONE-SIDED DERIVATIVES THAT
CONVERGES TO ONE SIDED DERIVATIVE AT THE CROSSECTION WITH THE BOTTOM
      PARAMETER(IPS=0)

CREATE TAPEWR
      PARAMETER(IANAL=0,IRST=0,IWRITE=1)
      PARAMETER (NFIL=4)
      DIMENSION DTFIL(NFIL),NTFIL(NFIL)
c     DATA DTFIL/NFIL*400./
c     DATA NTFIL/NFIL*797/
      DATA NTFIL/NFIL*1594/
      DATA DTFIL/NFIL*200./
      COMMON/ITERO/ ERROR,EXITCND,NITER,NITSM,ICOUNT
      DATA NITSM,ICOUNT,ERROR/0,0,1.e-10/
      COMMON/ELLIPS/ EPS,ITR,ITMN,IPRC
      common/surftot/ S_full

COMPUTATIONS:

COMPUTE SOME RELEVANT CONSTANTS
      PI=ACOS(-1.)
      PI2=2.*PI
      PIH=.5*PI
      TIME=0.
      PVEL=USCAL/R
      F0=ICORIO*F0
C BETA IS AN ANGLE IN THE ZONAL FLOW TEST
      BETA=PIH !*0.5
      H00=G*H00
      GI=1./G
C PARAMETERS FOR MPDATA ADVECTION
      IORD=2
      ISOR=1
      NONOS=1
      IDIV=1               
C SMALL CONSTANT (SMALL COMPARED TO DEPTH OF THE ATMOSPHERE)
      EP=1.E-6
C PARAMETERS FOR ELLIPTIC SOLVER 
      EPS=1.e-3 !minimum reduction of residual error
      ITR=100   !maximum # of iterations
      ITMN=1    !minimum # of itearations
! IPRC =: 0 no preconditioner, vanilla gcr; 1 diagonal preconditioner; 
!         2 diagonally preconditioned stationary  Richardson iteration
!         3 tri-diagonally preconditioned stationary  Richardson iteration
      IPRC=0 
COMPUTE GRID      
      DX=PI2/FLOAT(N-2)
      DY= PI/FLOAT(M)
      GC1=DT/DX
      GC2=DT/DY
      GH1=.5*GC1
      GH2=.5*GC2
      DO 10 J=1,M
   10 Y(J)=-PIH+(J-0.5)*DY
      DO 11 I=2,N-1
   11 X(I)=(I-1)*DX
      X(1)=X(N-1)
      X(N)=X(2)
      DO 12 I=1,N
   12 IP(I)=MOD(I+(N-2)/2-1,N-2)+1
COMPUTE METRIC TERMS FOR SPHERICAL COORDINATES
      DO 13 J=1,M
      DO 13 I=1,N
      HX(I,J)=R*COS(Y(J))
      HY(I,J)=R
CTEST DHX2Y(I,J)=-DT*R**2*SIN(2.*Y(J))/(HX(I,J)**2*HY(I,J))*.5
   13 S(I,J)=HX(I,J)*HY(I,J)
      S_full=0.0
      DO J=1,M
      DO I=2,N-1
      S_full=S_full+S(I,J)
      ENDDO
      ENDDO
      write(*,*) 'Area Earth', S_full
      DO 14 I=2,N-1
      DO 15 J=2,M-1
   15 DHX2Y(I,J)= (HX(I,J+1)-HX(I,J-1))*GC2/S(I,J)*.5
      DHX2Y(I,1)= (HX(I,  2)+HX(I,  1))*GC2/S(I,1)*.5
   14 DHX2Y(I,M)=-(HX(I,  M)+HX(I,M-1))*GC2/S(I,M)*.5
      CALL XBC(DHX2Y,N,M)

CONDITIONS OF THE INITIAL STATE ***********************************

      CALL TOPOGRAPHY(P0,X,Y,N,M)

      IF(IRHW.EQ.0) CALL INITZON(U,V,PT,COR,X,Y,N,M,F0,BETA,H00,R,PVEL)
      IF(IRHW.EQ.1) CALL INITRHW(U,V,PT,COR,X,Y,N,M,F0,R)
      CALL POLARABS(ALP,Y,DT,N,M)
      CALL DIVER(F1(1,1,1),U(1,1,0),V(1,1,0),HX,HY,S,N,M,IP,1)
      EXITCND=0.
      DO I=1,NM
      EXITCND = AMAX1(EXITCND, ABS(F1(I,1,1)))
      ENDDO
      write(*,*) 'maximal abs(div)', EXITCND

C INITIATE PRIMARY VARIABLES (PD, QX, QY)
      DO 20 I=1,NM
      COR(I,1)=COR(I,1)*DT 
      P0(I,1)=  P0(I,1)*HMTS
      PD(I,1)= AMAX1(0., PT(I,1)*GI-P0(I,1))
      QX(I,1)=PD(I,1)*U(I,1,0)
      QY(I,1)=PD(I,1)*V(I,1,0)
   20 CONTINUE

C STORE INITIAL FIELDS AND NORMS FOR LATER NORMALISATION
      DO I=1,NM
       PT0(I,1) = PT(I,1)
       PD0(I,1) = PD(I,1)
        U0(I,1) = U(I,1,0)
        V0(I,1) = V(I,1,0)
       PTTRL2 = PTTRL2 + S(I,1)*PT0(I,1)**2
       PTTRLI = AMAX1(PTTRLI,PT0(I,1))
       PDTRL2 = PDTRL2 + S(I,1)*PD0(I,1)**2
       PDTRLI = AMAX1(PDTRLI,PD0(I,1))
       VLTRL2 = VLTRL2 + S(I,1)*(U0(I,1)**2+V0(I,1))**2
       VLTRLI = AMAX1(VLTRLI, U0(I,1)**2+V0(I,1)**2)
      ENDDO

         IF(IRST.EQ.0) THEN
C INITIATE CORRESPONDING FORCES
      CALL PRFC0(PT,F1(1,1,0),F2(1,1,0),PD,HX,HY,IP,IPS,GH1,GH2,EP,N,M)
      DO 21 I=1,NM
      E1(I,1,0)=-DHX2Y(I,1)*QX(I,1)*QY(I,1)/AMAX1(PD(I,1),EP)
      E2(I,1,0)= DHX2Y(I,1)*QX(I,1)*QX(I,1)/AMAX1(PD(I,1),EP)
      E1(I,1,-1)=E1(I,1,0)
      E2(I,1,-1)=E2(I,1,0)
      F1(I,1,0)=F1(I,1,0)+COR(I,1)*QY(I,1)+E1(I,1,0)
   21 F2(I,1,0)=F2(I,1,0)-COR(I,1)*QX(I,1)+E2(I,1,0)

         ELSE
      DO 1200 KF=1,NFIL
      READ(19) PD,PT,QX,QY,U,V,F1,F2,E1,E2
      IF(IWRITE.EQ.1) WRITE(29) PD,PT,QX,QY,U,V,F1,F2,E1,E2
      TIME=TIME+DTFIL(KF)*NTFIL(KF)
 1200 CONTINUE
         ENDIF   !IRST

      CALL DIAGNOS(U,V,PD,PT,COR,HX,HY,S,X,XX,Y,IP,PI2,0,0)

CLOSE INITIAL CONDITIONS *************************************

COMPUTE SOLUTION IN TIME *************************************

      IF(IANAL.EQ.0) THEN
      DO 100 KT=1,NT
      if(kt/mpfl*mpfl.eq.kt) then
        liner=1
          else
            liner=0
              endif
       
      IPRINT=0
      IPLOT=0
      IF(KT.GE.KT0.AND.KT/NPRINT*NPRINT.EQ.KT) IPRINT=1
      IF(KT.GE.KT0.AND.KT/NPLOT*NPLOT.EQ.KT) IPLOT=1

      IF(LADVEC.eq.1 .and. KT.GT.1) go to 771
COMPUTE ADVECTIVE COURANT NUMBERS
COMPUTE VELOCITY PREDICTOR
      CALL VELPRD(U,V,F1,F2,PD,HX,HY,IP,FZ,SCR,N,M,GC1,GC2,EP)
      DO 105 I=1,NM
      U(I,1,1)=U(I,1,1)*HY(I,1)
  105 V(I,1,1)=V(I,1,1)*HX(I,1)
!! calculateAVGdepth(PD) and max abs(Div)
      D0=0.0
      DO J=1,M
      DO I=2,N-1
      D0=D0+PD(I,J)*S(I,J)/S_full
      END DO
      END DO
      if(iprint.eq.1) then 
      PRINT 299
  299 FORMAT(1X,1H )
      write(*,*) 'average Depth', D0
      endif
      CALL DIVER(F1(1,1,1),U(1,1,1),V(1,1,1),HX,HY,S,N,M,IP,1)
      EXITCND=0.
      DO I=1,NM
      EXITCND = AMAX1(EXITCND, ABS(F1(I,1,1)))
      ENDDO
      if(iprint.eq.1) write(*,*) 'maximal abs(div)', EXITCND

COMPUTE COURANT NUMBERS AT STAGGERED TIME/SPACE POSITIONS
      DO 110 J=1,M
      DO 110 I=2,N-1
  110 UA(I,J)=(U(I,J,1)+U(I-1,J,1))*GH1
      CALL XBC(UA,N,M)
      DO 111 I=1,N
      DO 112 J=2,M
  112 VA(I,J)=(V(I,J,1)+V(I,J-1,1))*GH2
      VA(I,  1)=0.
      VA(I,M+1)=0.
  111 CONTINUE
!     call compint(ua,va)

CLOSE ADVECTIVE COURANT NUMBERS

COLLECT EXPLICIT PARTS OF CONTINUITY AND MOMENTUM EQUATIONS:

COLLECT FROM CONTINUITY EQUATION
      DO 120 I=1,NM
      U(I,1,1)=QX(I,1)*GC1
  120 V(I,1,1)=QY(I,1)*GC2
      CALL DIVER(F1(1,1,1),U(1,1,1),V(1,1,1),HX,HY,S,N,M,IP,1)
      DO 121 I=1,NM
  121 F1(I,1,1)=PD(I,1)+P0(I,1)-.5*F1(I,1,1)

COLLECT FROM MOMENTUM EQUATION
C--->                       ADVECTION
      QXMC=-1.E30
      QYMC=-1.E30
      DO 130 I=1,NM
      QX(I,1)=QX(I,1)+.5*F1(I,1,0)
      QY(I,1)=QY(I,1)+.5*F2(I,1,0)
      QXMC=AMAX1(QXMC,ABS(QX(I,1)))
      QYMC=AMAX1(QYMC,ABS(QY(I,1)))
      PC(I,1)=PD(I,1)
  130 CONTINUE

  771 CONTINUE
      IF(LADVEC.eq.1) THEN 
        CALL MPDATT(UA,VA,PT,S,N,M,IORD,ISOR,NONOS,IDIV,1)
      GO TO 772
      ENDIF 
      CALL MPDATT(UA,VA,PC,S,N,M,IORD,ISOR,NONOS,IDIV,1)
      CALL MPDATT(UA,VA,QX,S,N,M,IORD,ISOR,NONOS,IDIV,-1)
      CALL MPDATT(UA,VA,QY,S,N,M,IORD,ISOR,NONOS,IDIV,-1)

C--->         CORIOLIS, METRIC & ABSORBING FORCES
      TIMSI=KT*DT
      IF(IRHW.EQ.1) THEN
       CALL RHWT(U0,V0,PT0,PD0,P0,COR,S,X,XX,Y,N,M,G,R,TIMSI)
       DO I=1,NM
        QXS(I,1)=U0(I,1)*PD0(I,1)
        QYS(I,1)=V0(I,1)*PD0(I,1)
       ENDDO
      ENDIF
      DO 139 I=1,NM
      UA(I,1)=QX(I,1)+.5*(2.*E1(I,1,0)-E1(I,1,-1))
     &               +.5*ALP(I,1)*QXS(I,1)
      VA(I,1)=QY(I,1)+.5*(2.*E2(I,1,0)-E2(I,1,-1))
     &               +.5*ALP(I,1)*QYS(I,1)
CTEST      UA(I,1)=QX(I,1)+.5*E1(I,1,0)
CTEST      VA(I,1)=QY(I,1)+.5*E2(I,1,0)
      GMM=.5*COR(I,1)
      AMM=1.+.5*ALP(I,1)
      DETI=1./(AMM**2+GMM**2)
      QX(I,1)=(AMM*DETI*UA(I,1)+GMM*DETI*VA(I,1))*GC1
      QY(I,1)=(AMM*DETI*VA(I,1)-GMM*DETI*UA(I,1))*GC2
      E1(I,1,-1)=E1(I,1,0)
  139 E2(I,1,-1)=E2(I,1,0)
C--->            DIVERGENCE OF ALL COLLECTED TERMS:
      CALL DIVER(F2(1,1,1),QX,QY,HX,HY,S,N,M,IP,1)

COMPUTE COEFFICIENTS FOR -.5*DIV(-C*PF(PT))-(PT-F1)=0 BV-PROBLEM
      DO 140 I=1,NM
      F1(I,1,1)=(F1(I,1,1)-.5*F2(I,1,1))*G
      F2(I,1,1)=PT(I,1)
      PC(I,1)=(PC(I,1)+P0(I,1))*G
  140 PD(I,1)=P0(I,1)*G
C PRFORC(...,0) computes coefficients of linearised D^{n+1}Grad H^{n+1} term 
c     CALL PRFORC(PT,E1(1,1,0),E2(1,1,0),PD,F2(1,1,1),  !t^n estimate, O(dt/2)
      CALL PRFORC(PC,E1(1,1,0),E2(1,1,0),PD,F2(1,1,1),  !t^{n+1}, fully O(dt^2)
     1            F1(1,1,0),F2(1,1,0),HX,HY,COR,ALP,N,M,IP,GC1,GC2,0,0)
COMPUTE ELLIPTIC EQUATION WITH GCR SCHEME
      CALL  GCR(PT,F1(1,1,0),F2(1,1,0),HX,HY,S,F1(1,1,1),F2(1,1,1),
     1          PD(1,1),E1(1,1,0),E2(1,1,0),COR,ALP,IP,
     2          U(1,1,0),U(1,1,1),V(1,1,0),V(1,1,1),N,M,GC1,GC2)
      CALL PRFORC(PT,F1(1,1,0),F2(1,1,0),PD,F2(1,1,1),
     1            E1(1,1,0),E2(1,1,0),HX,HY,COR,ALP,N,M,IP,GC1,GC2,1,1)
COMPUTE SOLUTION'S UPDATE
      DO 150 I=1,NM
      QX(I,1)=(QX(I,1)+F1(I,1,0)*GI)/GC1
      QY(I,1)=(QY(I,1)+F2(I,1,0)*GI)/GC2
  150 PD(I,1)= AMAX1(EP, PT(I,1)*GI-P0(I,1))

COMPUTE NEW FORCES
      CALL PRFC0(PT,F1(1,1,0),F2(1,1,0),PD,HX,HY,IP,IPS,GH1,GH2,EP,N,M)
      DO 160 I=1,NM
      E1(I,1,0)=-DHX2Y(I,1)*QX(I,1)*QY(I,1)/PD(I,1)
      E2(I,1,0)= DHX2Y(I,1)*QX(I,1)*QX(I,1)/PD(I,1)
      F1(I,1,0)=F1(I,1,0)+COR(I,1)*QY(I,1)+E1(I,1,0)
     &                   -ALP(I,1)*(QX(I,1)-QXS(I,1))
      F2(I,1,0)=F2(I,1,0)-COR(I,1)*QX(I,1)+E2(I,1,0)
     &                   -ALP(I,1)*(QY(I,1)-QYS(I,1))
      U(I,1,0)=QX(I,1)/PD(I,1)
  160 V(I,1,0)=QY(I,1)/PD(I,1)

  772 CONTINUE
COMPUTE OUTPUTED FIELDS ****************************************
      IF(KT.GE.KT0) THEN
        IF(KT/NPRINT*NPRINT.EQ.KT) 
     &    CALL DIAGNOS(U,V,PD,PT,COR,HX,HY,S,X,XX,Y,IP,PI2,KT,1)
        IF(KT/NPLOT*NPLOT.EQ.KT) THEN 
          IF(IWRITE.EQ.1) WRITE(9) PD,PT,QX,QY,U,V,F1,F2,E1,E2
        ENDIF
      ENDIF
  100 CONTINUE
CLOSE TIME INTEGRATION
      ENDIF  !IANAL

      IF(IANAL.EQ.1) THEN
      TIME=0.
      NPLOT=1
      DO 1201 KF=1,NFIL
      TIME=TIME+DTFIL(KF)*NTFIL(KF)/3600./24.
      PRINT 300, TIME      
  300 FORMAT(14X,5HTIME=,F7.2,5H DAYS)
      READ(9) PD,PT,QX,QY,U,V,F1,F2,E1,E2
 1201 CONTINUE
      ENDIF

      STOP 
      END   

CDIR$ NOLIST
      subroutine topography(h0,x,y,n,m)
      dimension h0(n,m),x(n),y(m)
      common/displ/ time,uscal,dx,dy,h00,hmts,dt
c mountain shape functions
      dist(xln,ylt) =sqrt( (cos(ylt)*sin((xln-x0)/2))**2
     .                      +       (sin((ylt-y0)/2))**2 )/xl
      profm(rad)=.5*(1.+cos(pi*rad))
c special for the conical mountain on the pole ccccc
c     dist(xln,ylt)=abs(ylt-y0)
c     profm(rad)=amax1(0., 1.-gamm*rad/rnot)
c     rnot = 2*acos(-1.)/128. * 10.
      gamm=1.
cccccccccccccccccccccccccccccccccccccccccccccccccccc

      pi = acos(-1.)
      xl=dx * 5.
      x0 = pi
      y0 = pi*.5/3.

      do 1 j=1,m
      do 1 i=1,n
      r=dist(x(i),y(j))
      if(r.le.1) then
      h0(i,j)=profm(r)
      else
      h0(i,j)=0.
      endif
    1 continue

      return
      end

      SUBROUTINE INITZON(U,V,PT,COR,X,Y,N,M,F0,BETA,H00,R,Q)
      DIMENSION U(N,M,0:1),V(N,M,0:1),PT(N,M),COR(N,M),X(N),Y(M)

      DO 20 J=1,M
      DO 20 I=1,N
      COR(I,J)=F0*(-COS(X(I))*COS(Y(J))*SIN(BETA)+SIN(Y(J))*COS(BETA))
   20 PT(I,J)=H00-R**2*(F0+Q)*0.5*Q*
     1        (-COS(X(I))*COS(Y(J))*SIN(BETA)+SIN(Y(J))*COS(BETA))**2
      DO 30 J=1,M
      DO 30 I=1,N
      U(I,J,0)=Q*(COS(BETA)+TAN(Y(J))*COS(X(I))*SIN(BETA))*R*COS(Y(J))
      V(I,J,0)=-Q*SIN(X(I))*SIN(BETA)*R
      U(I,J,1)=U(I,J,0)
   30 V(I,J,1)=V(I,J,0)
      write (6,*)  'initzon called'
      RETURN
      END

      SUBROUTINE INITRHW(U,V,PT,COR,X,Y,N,M,F0,A)
      DIMENSION U(N,M,0:1),V(N,M,0:1),PT(N,M),COR(N,M),X(N),Y(M)
      REAL K
      INTEGER R
      DATA OM,K,R,PH0/7.848E-6,7.848E-6,4,78.4E3/
      ATH(TH)=OM*0.5*(F0+OM)*(COS(TH))**2
     1   +0.25*K**2*(COS(TH))**(2*R)*( (R+1)*(COS(TH))**2
     1   +FLOAT(2*R**2-R-2)-2.*R**2/(COS(TH))**2 )
      BTH(TH)=(F0+2.*OM)*K/FLOAT((R+1)*(R+2))*(COS(TH))**R
     1       *( FLOAT(R**2+2*R+2)-((R+1)*COS(TH))**2 )
      CTH(TH)=0.25*K**2*(COS(TH))**(2*R)*( FLOAT(R+1)*(COS(TH))**2
     1       -FLOAT(R+2) )  

      DO 20 J=1,M
      DO 20 I=1,N
      COR(I,J)=F0*SIN(Y(J))
      U(I,J,0)=A*OM*COS(Y(J))+A*K*COS(R*X(I))
     1        *(COS(Y(J)))**(R-1)*(R*(SIN(Y(J)))**2-(COS(Y(J)))**2)
      V(I,J,0)=-A*K*R*(COS(Y(J)))**(R-1)*SIN(Y(J))*SIN(R*X(I))
      PT(I,J)=PH0+A**2*ATH(Y(J))+A**2*BTH(Y(J))*COS(R*X(I))
     1       +A**2*CTH(Y(J))*COS(2.*R*X(I))
      U(I,J,1)=U(I,J,0)
      V(I,J,1)=V(I,J,0)
   20 CONTINUE
      write (6,*)  'initrhw called'
      RETURN
      END

      SUBROUTINE RHWT(U0,V0,PT0,PD0,P0,COR,S,X,XX,Y,N,M,GR,A,TIME)
      DIMENSION U0(N,M),V0(N,M),PT0(N,M),PD0(N,M),P0(N,M),
     &                     COR(N,M),S(N,M),X(N),XX(N),Y(M)
      COMMON/INITVLS/ PTTRL2,PTTRLI,VLTRL2,VLTRLI,PDTRL2,PDTRLI
      REAL K,VNIU
      INTEGER R
      DATA OM,K,R,PH0/7.848E-6,7.848E-6,4,78.4E3/
      DATA F0/1.4584E-4/
      ATH(TH)=OM*0.5*(F0+OM)*(COS(TH))**2
     1   +0.25*K**2*(COS(TH))**(2*R)*( (R+1)*(COS(TH))**2
     1   +FLOAT(2*R**2-R-2)-2.*R**2/(COS(TH))**2 )
      BTH(TH)=(F0+2.*OM)*K/FLOAT((R+1)*(R+2))*(COS(TH))**R
     1       *( FLOAT(R**2+2*R+2)-((R+1)*COS(TH))**2 )
      CTH(TH)=0.25*K**2*(COS(TH))**(2*R)*( FLOAT(R+1)*(COS(TH))**2
     1       -FLOAT(R+2) )  

      OM=7.848E-6
      K =7.848E-6
      R =4
      PH0=78.4E3
      F0=1.4584E-4
      GRI=1./GR
      PI2=2.*ACOS(-1.)
      VNIU = (Float(R*(3+R))*OM-F0)/Float((1+R)*(2+R))
      IF(TIME.EQ.0.) write (6,*)  'vniu:',VNIU

      DO I=1,N
       XX(I)=X(I)-VNIU*TIME
       IF(XX(I).LT.0.) XX(I)=XX(I)+PI2
       IF(XX(I).GT.PI2) XX(I)=XX(I)-PI2
      ENDDO
       XX(1)=XX(N-1)
       XX(N)=XX(2)
      IF(TIME.EQ.0.) write (6,*) 
     &  'x(1),x(2),x(n-1),x(n)',x(1),x(2),x(n-1),x(n)

      DO 20 J=1,M
      DO 20 I=1,N
      U0(I,J)=A*OM*COS(Y(J))+A*K*COS(R*XX(I))
     1        *(COS(Y(J)))**(R-1)*(R*(SIN(Y(J)))**2-(COS(Y(J)))**2)
      V0(I,J)=-A*K*R*(COS(Y(J)))**(R-1)*SIN(Y(J))*SIN(R*XX(I))
      PT0(I,J)=PH0+A**2*ATH(Y(J))+A**2*BTH(Y(J))*COS(R*XX(I))
     1       +A**2*CTH(Y(J))*COS(2.*R*XX(I))
      PD0(I,J)= AMAX1(0., PT0(I,J)*GRI-P0(I,J))
   20 CONTINUE

      PTTRL2=0.
      PTTRLI=-1.e-15
      VLTRL2=0.
      VLTRLI=-1.e-15
      PDTRL2=0.
      PDTRLI=-1.e-15
      DO J=1,M
      DO I=1,N
       PTTRL2 = PTTRL2 + S(I,J)*PT0(I,J)**2
       PTTRLI = AMAX1(PTTRLI,PT0(I,J))
       PDTRL2 = PDTRL2 + S(I,J)*PD0(I,J)**2
       PDTRLI = AMAX1(PDTRLI,PD0(I,J))
       VLTRL2 = VLTRL2 + S(I,J)*(U0(I,J)**2+V0(I,J))**2
       VLTRLI = AMAX1(VLTRLI, U0(I,J)**2+V0(I,J)**2)
      ENDDO
      ENDDO
!     write (6,*)  'PTTRL2, PTTRLI:',PTTRL2, PTTRLI
!     write (6,*)  'VLTRL2, VLTRLI:',VLTRL2, VLTRLI
!     write (6,*)  'PDTRL2, PDTRLI:',PDTRL2, PDTRLI

      RETURN
      END

      SUBROUTINE PRFC0(A,F1,F2,PD,HX,HY,IP,IPS,GH1,GH2,EP,N,M)
      DIMENSION A(N,M),PD(N,M),HX(N,M),HY(N,M),F1(N,M),F2(N,M),IP(N)

      DO 22 I=2,N-1
      DO 23 J=1,M
      GP=PD(I+1,J)/(PD(I+1,J)+PD(I-1,J)+EP)*2.*IPS+(1-IPS)
      GN=PD(I-1,J)/(PD(I+1,J)+PD(I-1,J)+EP)*2.*IPS+(1-IPS)
   23 F1(I,J)=-GH1*PD(I,J)/HX(I,J)*
     1 ( GP*(A(I+1,J)-A(I,J))+GN*(A(I,J)-A(I-1,J)) )
      DO 24 J=2,M-1
      GP=PD(I,J+1)/(PD(I,J+1)+PD(I,J-1)+EP)*2.*IPS+(1-IPS)
      GN=PD(I,J-1)/(PD(I,J+1)+PD(I,J-1)+EP)*2.*IPS+(1-IPS)
   24 F2(I,J)=-GH2*PD(I,J)/HY(I,J)*
     1 ( GP*(A(I,J+1)-A(I,J))+GN*(A(I,J)-A(I,J-1)) )
      GP=PD(I    ,2)/(PD(I,2)+PD(IP(I),1)+EP)*2.*IPS+(1-IPS)
      GN=PD(IP(I),1)/(PD(I,2)+PD(IP(I),1)+EP)*2.*IPS+(1-IPS)
      F2(I,1)=-GH2*PD(I,1)/HY(I,1)*
     1 ( GP*(A(I,2)-A(I,1))+GN*(A(I,1)-A(IP(I),1)) )
      GP=PD(IP(I),M)/(PD(IP(I),M)+PD(I,M-1)+EP)*2.*IPS+(1-IPS)
      GN=PD(I  ,M-1)/(PD(IP(I),M)+PD(I,M-1)+EP)*2.*IPS+(1-IPS)
      F2(I,M)=-GH2*PD(I,M)/HY(I,M)*
     1 ( GP*(A(IP(I),M)-A(I,M))+GN*(A(I,M)-A(I,M-1)) )
   22 CONTINUE
      CALL XBC(F1,N,M)
      CALL XBC(F2,N,M)
      RETURN
      END

      SUBROUTINE PRFORC( P,F1,F2,PB,P0,E1,E2,HX,HY,COR,ALP,
     1                              N,M,IP,GC1,GC2,NOR,IRS)
      DIMENSION P(N,M),F1(N,M),F2(N,M),PB(N,M),P0(N,M),E1(N,M),E2(N,M),
     1                         HX(N,M),HY(N,M),COR(N,M),ALP(N,M),IP(N)

      GH1=.5*GC1
      GH2=.5*GC2

      DO 1 J=2,M-1
      DO 1 I=2,N-1
      F1(I,J)=-GH1*(P(I+1,J)-P(I-1,J))/HX(I,J)
    1 F2(I,J)=-GH2*(P(I,J+1)-P(I,J-1))/HY(I,J)
      DO 2 I=2,N-1
      F1(I,1)=-GH1*(P(I+1,1)-P(I-1,1))/HX(I,1)
      F1(I,M)=-GH1*(P(I+1,M)-P(I-1,M))/HX(I,M)
      F2(I,1)=-GH2*(P(I,2)-P(IP(I),1))/HY(I,1) 
    2 F2(I,M)=-GH2*(P(IP(I),M)-P(I,M-1))/HY(I,M)
      CALL XBC(F1,N,M)
      CALL XBC(F2,N,M)

      IF(NOR.EQ.1) THEN
      NM=N*M
      DO 3 I=1,NM
      UTILD=F1(I,1)*(P0(I,1)-PB(I,1))+(P(I,1)-IRS*P0(I,1))*E1(I,1)
      VTILD=F2(I,1)*(P0(I,1)-PB(I,1))+(P(I,1)-IRS*P0(I,1))*E2(I,1)
      GMM=.5*COR(I,1)
      AMM=1.+.5*ALP(I,1)
      DETI=1./(AMM**2+GMM**2)
      A=AMM*DETI
      B=GMM*DETI
      F1(I,1)=(A*UTILD+B*VTILD)*GH1
      F2(I,1)=(A*VTILD-B*UTILD)*GH2
    3 CONTINUE
      ENDIF

      RETURN
      END

      SUBROUTINE DIVER(F,U,V,HX,HY,S,N,M,IP,IFLG)
      DIMENSION F(N,M),U(N,M),V(N,M),HX(N,M),HY(N,M),S(N,M),IP(N)

      DO 1 J=2,M-1
      DO 1 I=2,N-1
    1 F(I,J)= HY(I+1,J)*U(I+1,J)-HY(I-1,J)*U(I-1,J) 
     .       +HX(I,J+1)*V(I,J+1)-HX(I,J-1)*V(I,J-1)
      DO 2 I=2,N-1
      F(I,1)= HY(I+1,1)*U(I+1,1)-HY(I-1,1)*U(I-1,1)
     .       +(HX(I,2)*V(I,2)+HX(I,1)*V(I,1))
    2 F(I,M)= HY(I+1,M)*U(I+1,M)-HY(I-1,M)*U(I-1,M)
     .       -(HX(I,M)*V(I,M)+HX(I,M-1)*V(I,M-1))

      DO 3 J=1,M
      DO 3 I=2,N-1
    3 F(I,J)=.5*IFLG*F(I,J)/S(I,J)

      DO 4 J=1,M
      F(1,J)=F(N-1,J)
    4 F(N,J)=F(2  ,J)
      RETURN
      END

      SUBROUTINE LAP0(A11,A12,A21,A22,B11,B22,
     .                PB,P0,E1,E2,HX,HY,COR,ALP,N,M,GC1,GC2)
      DIMENSION A11(N,M),A12(N,M),A21(N,M),A22(N,M),B11(N,M),B22(N,M),
     .PB(N,M),P0(N,M),E1(N,M),E2(N,M),HX(N,M),HY(N,M),COR(N,M),ALP(N,M)

      NM=N*M
      GH1=.5*GC1
      GH2=.5*GC2

      DO I=1,NM
      C1=-GH1/HX(I,1)*(P0(I,1)-PB(I,1))
      C2=-GH2/HY(I,1)*(P0(I,1)-PB(I,1))
      GMM=.5*COR(I,1)
      AMM=1.+.5*ALP(I,1)
      DETI=1./(AMM**2+GMM**2)
      A=AMM*DETI
      B=GMM*DETI
      A11(I,1)=-C1*A*GH1*HY(I,1)*.5
      A12(I,1)=-C2*B*GH1*HY(I,1)*.5
      A21(I,1)= C1*B*GH2*HX(I,1)*.5
      A22(I,1)=-C2*A*GH2*HX(I,1)*.5
      B11(I,1)=-   A*GH1*E1(I,1)*HY(I,1)*.5
     .         -   B*GH1*E2(I,1)*HY(I,1)*.5
      B22(I,1)=-   A*GH2*E2(I,1)*HX(I,1)*.5
     .         +   B*GH2*E1(I,1)*HX(I,1)*.5
      ENDDO

      RETURN
      END

      SUBROUTINE LAPL(P,F,A11,A12,A21,A22,B11,B22,U,V,S,N,M,IP)
      DIMENSION P(N,M),F(N,M),A11(N,M),A12(N,M),A21(N,M),A22(N,M),
     1          B11(N,M),B22(N,M),U(N,M),V(N,M),S(N,M),IP(N)

      NM=N*M

      DO 1 J=2,M-1
      DO 1 I=2,N-1
      U(I,J)=P(I+1,J)-P(I-1,J)
    1 V(I,J)=P(I,J+1)-P(I,J-1)
      DO I=2,N-1
      U(I,1)=P(I+1,1)-P(I-1,1)
      U(I,M)=P(I+1,M)-P(I-1,M)
      V(I,1)=P(I,2)-P(IP(I),1)
      V(I,M)=P(IP(I),M)-P(I,M-1)
      ENDDO
      DO J=1,M
      U(1,J)=U(N-1,J)
      U(N,J)=U(2  ,J)
      V(1,J)=V(N-1,J)
      V(N,J)=V(2  ,J)
      ENDDO

      DO I=1,NM
      UTIL=A11(I,1)*U(I,1)+A12(I,1)*V(I,1)+B11(I,1)*P(I,1)
      VTIL=A21(I,1)*U(I,1)+A22(I,1)*V(I,1)+B22(I,1)*P(I,1)
      U(I,1)=UTIL
      V(I,1)=VTIL
      ENDDO

      DO 2 J=2,M-1
      DO 2 I=2,N-1
    2 F(I,J)= (U(I+1,J)-U(I-1,J)+V(I,J+1)-V(I,J-1))/S(I,J)
      DO I=2,N-1
      F(I,1)= (U(I+1,1)-U(I-1,1)+(V(I,2)+V(I,1)))/S(I,1)
      F(I,M)= (U(I+1,M)-U(I-1,M)-(V(I,M)+V(I,M-1)))/S(I,M)
      ENDDO
      DO J=1,M
      F(1,J)=F(N-1,J)
      F(N,J)=F(2  ,J)
      ENDDO

      RETURN
      END

      subroutine gcr(p,pfx,pfy,hx,hy,s,b,p0,pb,e1,e2,cor,alp,
     1                            ip,d,qr,r,ar,n1,n2,gc1,gc2)
      dimension p(*),pfx(*),pfy(*),hx(*),hy(*),s(*),
     1          b(*),pb(*),p0(*),
     2          e1(*),e2(*),cor(*),alp(*),ip(*),
     3          d(*),qr(*),r(*),ar(*)
      common/itero/ error,exitcnd,niter,nitsm,icount
      common/ellips/ eps,itr,itmn,iprc
      common/printl/ iprint
c     PARAMETER(NLON=130,NLAT=64)
c     PARAMETER(NLON=258,NLAT=128)
      PARAMETER(NLON=514,NLAT=256)
      parameter(nm=nlon*nlat,lord=3)
      common// x(nm,lord),ax(nm,lord),ax2(lord),axar(lord),del(lord)
      dimension a11(nm),a12(nm),a21(nm),a22(nm),b11(nm),b22(nm)
c
      nml=n1*n2
      epa=1.e-38
      nlc=0

      call lap0(a11,a12,a21,a22,b11,b22,
     .          pb,p0,e1,e2,hx,hy,cor,alp,n1,n2,gc1,gc2)
      call diagoc(d,a11,a22,s)
      call precon(r,qr,ar,a11,a12,a21,a22,b11,b22,s,d,pfx,pfy,ip,iprc,0)

      do k=1,nml
        r(k)=0.
       ar(k)=0.
       qr(k)=0.
      enddo
      do l=1,lord
       do k=1,nml
         x(k,l)=0.
        ax(k,l)=0.
       enddo
      enddo
      call prforc(p,pfx,pfy,pb,p0,e1,e2,hx,hy,cor,alp,
     &                           n1,n2,ip,gc1,gc2,1,1)
      call diver(r,pfx,pfy,hx,hy,s,n1,n2,ip,-1)
      do k=1,nml
       r(k)=.5*r(k)-(p(k)-b(k))
      enddo
      call precon(r,qr,ar,a11,a12,a21,a22,b11,b22,s,d,pfx,pfy,ip,iprc,1)

      err0=0.
      do k=1,nml
       err0=err0+r(k)*r(k)
c      err0=amax1(err0,abs(r(k,1)))
      enddo
       err0=sqrt(err0)
       errnm1=err0

      do k=1,nml
       x(k,1)=qr(k)
      enddo
      call lapl(x(1,1),ax(1,1),a11,a12,a21,a22,b11,b22,
     .                              pfx,pfy,s,n1,n2,ip)
      do k=1,nml
       ax(k,1)=.5*ax(k,1)-x(k,1)
      enddo

      do 100 it=1,itr
       do l=1,lord
        ax2(l)=0.
        rax=0.
         do k=1,nml
          rax=rax+r(k)*ax(k,l)
          ax2(l)=ax2(l)+ax(k,l)*ax(k,l)
         enddo
        ax2(l)=amax1(epa,ax2(l))
        beta=-rax/ax2(l)
!       write(*,*) 'beta', beta
        errn=0.
         do k=1,nml
          p(k)=p(k)+beta* x(k,l)
          r(k)=r(k)+beta*ax(k,l)
          errn=errn+r(k)*r(k)
         enddo
c       error=errn/err0
        errn=sqrt(errn)
       write(*,*) it, l, errn, err0
       if(errn.lt.eps*err0.and.it.ge.itmn) go to 200
       if(errn.ge.errnm1) go to 200
        errnm1=errn
c      if(error.lt.eps*err0) go to 200
c      if(error.lt.eps*exitcnd) go to 200
      call precon(r,qr,ar,a11,a12,a21,a22,b11,b22,s,d,pfx,pfy,ip,iprc,1)
      call lapl(qr,ar,a11,a12,a21,a22,b11,b22,pfx,pfy,s,n1,n2,ip)
        do k=1,nml
         ar(k)=.5*ar(k)-qr(k)
        enddo
        nlc=nlc+1
         do ll=1,l
          axar(ll)=0.
           do k=1,nml
            axar(ll)=axar(ll)+ax(k,ll)*ar(k)
           enddo
          del(ll)=-axar(ll)/ax2(ll)
c         del(ll)=amax1(del(ll),0.)
         enddo
        if(l.lt.lord) then
          do k=1,nml
            x(k,l+1)=qr(k)
           ax(k,l+1)=ar(k)
          enddo
           do ll=1,l
            do k=1,nml
              x(k,l+1)= x(k,l+1)+del(ll)* x(k,ll)
             ax(k,l+1)=ax(k,l+1)+del(ll)*ax(k,ll)
            enddo
           enddo
         else
          do k=1,nml
            x(k,1)=qr(k)+del(1)* x(k,1)
           ax(k,1)=ar(k)+del(1)*ax(k,1)
          enddo
           do ll=2,l
            do k=1,nml
              x(k,1 )= x(k,1)+del(ll)* x(k,ll)
              x(k,ll)=0.
             ax(k,1 )=ax(k,1)+del(ll)*ax(k,ll)
             ax(k,ll)=0.
            enddo
           enddo
         endif
       enddo
  100 continue
  200 niter=nlc
      icount=icount+1
      nitsm=nitsm+niter
      error=errn/err0

      if(iprint.eq.1) then
c     i1=int(102.4+409.6)
c     call set(.1,.9,.3,.7,0.,2.,-0.5,0.5,1)
c     call labmod('(f4.1)','(f4.1)',4,4,2,2,20,20,0)
c     call periml(4,5,2,5)
c     x2=float(n1)
c     y2=float(n2)
c     call set(.1,.9,.3,.7,1.,x2,1.,y2,1)
c     call wtstr(cpux(i1),cpuy(200),'x/Pi',3,0,0)
c     call wtstr(cpux(22),cpuy(i1),'y/Pi',3,90,0)
c     call conrec(r,n1,n1,n2,0.,0.,0.,1,-1,-682)
c     call framE
      endiF

      return
      end

      subroutine precon(rhs,p,r,a11,a12,a21,a22,b11,b22,s,dgc,
     &                                     pfx,pfy,ip,iflg,jfl)
c     parameter(nlon=138,nlat=68)
c     PARAMETER(NLON=130,NLAT=64)
c     PARAMETER(NLON=258,NLAT=128)
      PARAMETER(NLON=514,NLAT=256)
      parameter(n=nlon,m=nlat,nm=nlon*nlat,nm1=n-1)
      dimension rhs(n,m),p(n,m),r(n,m),
     1          a11(n,m),a12(n,m),a21(n,m),a22(n,m),b11(n,m),b22(n,m),
     2          s(n,m),dgc(n,m),pfx(n,m),pfy(n,m),ip(n)
      dimension c11(n,m),po(n,m),qr(n,m)
      dimension f(0:nm1-1,m),e(0:nm1-1,m),g(0:nm1-1,m,2),q(nm1,m,4),
     1          dgh(n,m),aa(m,4)

      data betap/-1.e15/
      data itr,line/2,0/
      data swc,swcp/1.,1./
      det3(r11,r12,r13,r21,r22,r23,r31,r32,r33)=
     .           r11*r22*r33+r12*r23*r31+r13*r21*r32
     .          -r31*r22*r13-r32*r23*r11-r33*r21*r12


      IF(IFLG.EQ.0) THEN
       do j=1,m
       do i=1,n
        p(i,j)=rhs(i,j)
       enddo
       enddo
      return
      ENDIF

      IF(IFLG.EQ.1) THEN
       do j=1,m
        do i=1,n
         p(i,j)=rhs(i,j)/dgc(i,j)
        enddo
       enddo
       call adjust_conservation(p,s,n,m)
       return
      ENDIF

! The jfl.eq.0 code up to beti= may be useful for option 3. 
! Options 0 to 2 have their own definitions of beta, if required.
      if(jfl.eq.0) then
       if(line.eq.0) then
       beta=-1.e15
        if(iflg.lt.3) then
         do j=1,m
         do i=1,n
          betap=amax1(betap,(abs(a11(i,j))+abs(a22(i,j)))/(2.*s(i,j)))
         enddo
         enddo
        else
         do j=2,m-1
         do i=1,n
          betap=amax1(betap, 0.5*abs(a22(i,j+1)+a22(i,j-1))/s(i,j))
         enddo
         enddo
         do i=1,n
          betap=amax1(betap,0.5*abs(a22(i,2)+a22(ip(i),1))/s(i,1))
          betap=amax1(betap,0.5*abs(a22(ip(i),m)+a22(i,m-1))/s(i,m-1))
         enddo
        endif
       betap=0.25/betap
c        write(*,*) 'betap', betap
       else
       betap=1.
       endif
       return
      endif

      beti=1./betap !*(1-line)

      IF(IFLG.EQ.2) THEN
       betap=0.5
      do j=1,m
      do i=1,n
       p(i,j)=-betap*rhs(i,j)/dgc(i,j)
       r(i,j)=0.
      enddo
      enddo
      call adjust_conservation(p,s,n,m)
      call linop(p,r,pfx,pfy,a11,a12,a21,a22,b11,b22,dgc,s,ip,swc,iflg)

      do it=1,itr
       do j=1,m
       do i=1,n
        p(i,j)=p(i,j)+betap*(r(i,j)-rhs(i,j))/dgc(i,j)
       enddo
       enddo
       call adjust_conservation(p,s,n,m)
       if(it.lt.itr) call linop(p,r,pfx,pfy,a11,a12,a21,a22,b11,b22,
     &                                            dgc,s,ip,swc,iflg)
      enddo
      ENDIF
    
      IF(IFLG.EQ.3) THEN    
      omg=1.
      oms=1.-omg
      do j=1,m
      do i=1,n
       dgh(i,j)=0.
       po(i,j)=0.
       p(i,j)=0.
       r(i,j)=0.
       c11(i,j)=0.5*a11(i,j)
      enddo
      enddo

      if(line.eq.1) then
       do j=2,m-1
        do i=1,n
         dgh(i,j)=0.5*(a22(i,j+1)+a22(i,j-1))+1.
        enddo
       enddo
        do i=1,n
         dgh(i,1)=0.5*a22(i,2)+1.
         dgh(i,m)=0.5*a22(i,m-1)+1.
        enddo
      endif

      do 100 it=1,itr
       do j=1,m
       do i=1,n
        r(i,j)=r(i,j)+s(i,j)*(beti*p(i,j)-rhs(i,j))
     &        +dgh(j,j)*p(i,j)
       enddo
       enddo

       do j=1,m
       e(0,j)=0.
       f(0,j)=0.
       dn=c11(2,j)+c11(nm1-1,j)+s(1,j)*(beti+1.)+dgh(1,j)
       dni=1./dn
       e(1,j) = c11(2,j)*dni
       f(1,j) = r(1,j)*dni
       g(0,j,1) = 1.
       g(1,j,1) = 0.
       g(0,j,2) = 0.
       g(1,j,2) = c11(nm1-1,j)*dni
       enddo

       do i=2,nm1-1
       do j=1,m
        dn=c11(i+1,j)+c11(i-1,j)*(1.-e(i-2,j))
     &    +s(i,j)*(beti+1.)+dgh(i,j)
        dni=1./dn
        e(i,j)=                  c11(i+1,j)*dni
        f(i,j)=(c11(i-1,j)*f(i-2,j)+r(i,j))*dni
        g(i,j,1)=(c11(i-1,j)*g(i-2,j,1)   )*dni
        g(i,j,2)=(c11(i-1,j)*g(i-2,j,2)   )*dni
       enddo
       enddo

       il=nm1-1  ! remove compiler wornings for exceed array dimension
       do j=1,m
        p(il+1,j)=0.
        p(il,j)=f(il,j)
        q(il+1,j,1)=0.
        q(il,j,1)=g(il,j,1)
        q(il+1,j,2)=1.
        q(il,j,2)=0.
        q(il+1,j,3)=0.
        q(il,j,3)=g(il,j,2)
        q(il+1,j,4)=0.
        q(il,j,4)=e(il,j)
       enddo

       il=nm1-2  
       do j=1,m
       do i=il,1,-1
        p(i,j)=e(i,j)*p(i+2,j)+f(i,j)
        q(i,j,1)=e(i,j)*q(i+2,j,1)+g(i,j,1)
        q(i,j,2)=e(i,j)*q(i+2,j,2)
        q(i,j,3)=e(i,j)*q(i+2,j,3)+g(i,j,2)
        q(i,j,4)=e(i,j)*q(i+2,j,4)
       enddo
       enddo

       il1=nm1-1
       il2=nm1-2
       i2 =2
       do j=1,m
        d11= q(il1,j,1)-1.
        d12= q(il1,j,2)
        d13= q(il1,j,3)
        d14= q(il1,j,4)
        s1 =-p(il1,j)
        d21= q(1,j,1)
        d22= q(1,j,2)-1.
        d23= q(1,j,3)
        d24= q(1,j,4)
        s2 =-p(1,j)
        d31= q(il2,j,1)
        d32= q(il2,j,2)
        d33= q(il2,j,3)-1.
        d34= q(il2,j,4)
        s3 =-p(il2,j)
        d41= q(i2,j,1)
        d42= q(i2,j,2)
        d43= q(i2,j,3)
        d44= q(i2,j,4)-1.
        s4 =-p(i2,j)

        det40=d11*det3(d22,d23,d24,d32,d33,d34,d42,d43,d44)
     .       -d21*det3(d12,d13,d14,d32,d33,d34,d42,d43,d44)
     .       +d31*det3(d12,d13,d14,d22,d23,d24,d42,d43,d44)
        deti=1./det40
        det41=s1 *det3(d22,d23,d24,d32,d33,d34,d42,d43,d44)
     .       -s2 *det3(d12,d13,d14,d32,d33,d34,d42,d43,d44)
     .       +s3 *det3(d12,d13,d14,d22,d23,d24,d42,d43,d44)
     .       -s4 *det3(d12,d13,d14,d22,d23,d24,d32,d33,d34)
        det42=d11*det3( s2,d23,d24, s3,d33,d34, s4,d43,d44)
     .       -d21*det3( s1,d13,d14, s3,d33,d34, s4,d43,d44)
     .       +d31*det3( s1,d13,d14, s2,d23,d24, s4,d43,d44)
     .       -d41*det3( s1,d13,d14, s2,d23,d24, s3,d33,d34)
        det43=d11*det3(d22, s2,d24,d32, s3,d34,d42, s4,d44)
     .       -d21*det3(d12, s1,d14,d32, s3,d34,d42, s4,d44)
     .       +d31*det3(d12, s1,d14,d22, s2,d24,d42, s4,d44)
     .       -d41*det3(d12, s1,d14,d22, s2,d24,d32, s3,d34)
        det44=d11*det3(d22,d23, s2,d32,d33, s3,d42,d43, s4)
     .       -d21*det3(d12,d13, s1,d32,d33, s3,d42,d43, s4)
     .       +d31*det3(d12,d13, s1,d22,d23, s2,d42,d43, s4)
     .       -d41*det3(d12,d13, s1,d22,d23, s2,d32,d33, s3)

        aa(j,4)=det44*deti
        aa(j,3)=det43*deti
        aa(j,2)=det42*deti
        aa(j,1)=det41*deti
 
       enddo

       do j=1,m
       do i=1,nm1
       p(i,j)=p(i,j)+aa(j,1)*q(i,j,1)
     .              +aa(j,2)*q(i,j,2)
     .              +aa(j,3)*q(i,j,3)
     .              +aa(j,4)*q(i,j,4)
       enddo
       enddo
correct for round-off departures from the cyclicity in the vertical
       do j=1,m
        p(1,j)=p(nm1,j)
        p(n,j)=p(2,j)
       enddo
       call adjust_conservation(p,s,n,m)

      if(line.eq.1) then
       do j=1,m
       do i=1,n
        p(i,j)=oms*po(i,j)+omg*p(i,j)
       po(i,j)=     p(i,j)
       enddo
       enddo
      endif
      
      if(it.eq.itr) goto 101
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

      do j=1,m
      do i=1,n
       util=                  swcp*(a12(i,j)*pfy(i,j)+b11(i,j)*p(i,j))
       vtil=pfy(i,j)*a22(i,j)+swcp*(a21(i,j)*pfx(i,j)+b22(i,j)*p(i,j))
       pfx(i,j)=util
       pfy(i,j)=vtil
      enddo
      enddo

      do j=2,m-1
       do i=2,n-1
        r(i,j)= (pfx(i+1,j)-pfx(i-1,j)+pfy(i,j+1)-pfy(i,j-1))/s(i,j)
       enddo
      enddo
      do i=2,n-1
       r(i,1)= (pfx(i+1,1)-pfx(i-1,1)+(pfy(i,2)+pfy(i,1)))/s(i,1)
       r(i,m)= (pfx(i+1,m)-pfx(i-1,m)-(pfy(i,m)+pfy(i,m-1)))/s(i,m)
      enddo
      do j=1,m
       r(1,j)=r(n-1,j)
       r(n,j)=r(2  ,j)
      enddo
      do j=1,m
      do i=1,n
       r(i,j)=0.5*r(i,j) !-p(i,j)
      enddo
      enddo
  100 continue
  101 continue
      
      ENDIF

      return
      end

      subroutine adjust_conservation(p,s,n,m)
      dimension p(n,m),s(n,m)
      common/surftot/ S_full

      cnst=0.
      do j=1,m
        do i=2,n-1
          cnst=cnst+s(i,j)*p(i,j)
        enddo
      enddo
      cnst=cnst/S_full

      do j=1,m
        do i=1,n
          p(i,j)=p(i,j)-cnst
        enddo
      enddo
 
      return
      end

      subroutine diagoc(dgc,a11,a22,s)
c     parameter(nlon=138,nlat=68)
c     PARAMETER(NLON=130,NLAT=64)
c     PARAMETER(NLON=258,NLAT=128)
      PARAMETER(NLON=514,NLAT=256)
      parameter(n=nlon,m=nlat,nm=nlon*nlat)
      dimension dgc(n,m),a11(n,m),a22(n,m),s(n,m)
      common/ellips/ eps,itr,itmn,iprc

      IF(IPRC.EQ.0) THEN

      do j=1,m
       do i=1,n
        dgc(i,j)=1.
       enddo
      enddo

      ELSE

      do j=1,m
       do i=2,n-1
        dgc(i,j)=a11(i+1,j)+a11(i-1,j)
       enddo
        dgc(1,j)=dgc(n-1,j)  !a11(2,j)+a11(n-2,j)
        dgc(n,j)=dgc( 2 ,j)  !a11(3,j)+a11(n-1,j)
      enddo

      do j=2,m-1
      do i=1,n
       dgc(i,j)=dgc(i,j)+a22(i,j+1)+a22(i,j-1)
      enddo
      enddo
      do i=1,n
       dgc(i,1)=dgc(i,1)+a22(i,2)
       dgc(i,m)=dgc(i,m)+a22(i,m-1)
      enddo

      do j=1,m
       do i=1,n
        dgc(i,j)=0.5*dgc(i,j)/s(i,j)+1.
       enddo
      enddo

      ENDIF

      return
      end

      subroutine linop(p,r,pfx,pfy,a11,a12,a21,a22,b11,b22,
     &                                   dgc,s,ip,swc,iflg)
c     parameter(nlon=138,nlat=68)
c     PARAMETER(NLON=130,NLAT=64)
c     PARAMETER(NLON=258,NLAT=128)
      PARAMETER(NLON=514,NLAT=256)
      parameter(n=nlon,m=nlat,nm=nlon*nlat)
      dimension p(n,m),r(n,m),
     &          a11(n,m),a12(n,m),a21(n,m),a22(n,m),b11(n,m),b22(n,m),
     &          dgc(n,m),s(n,m),pfx(n,m),pfy(n,m),ip(n)

      IF(IFLG.EQ.3) return

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

      do j=1,m
      do i=1,n
       util=pfx(i,j)*a11(i,j)+swc*(a12(i,j)*pfy(i,j)+b11(i,j)*p(i,j))
       vtil=pfy(i,j)*a22(i,j)+swc*(a21(i,j)*pfx(i,j)+b22(i,j)*p(i,j))
       pfx(i,j)=util
       pfy(i,j)=vtil
      enddo
      enddo

      do j=2,m-1
       do i=2,n-1
        r(i,j)= (pfx(i+1,j)-pfx(i-1,j)+pfy(i,j+1)-pfy(i,j-1))/s(i,j)
       enddo
      enddo
      do i=2,n-1
       r(i,1)= (pfx(i+1,1)-pfx(i-1,1)+(pfy(i,2)+pfy(i,1)))/s(i,1)
       r(i,m)= (pfx(i+1,m)-pfx(i-1,m)-(pfy(i,m)+pfy(i,m-1)))/s(i,m)
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

      return
      end

      SUBROUTINE VELPRD(U,V,F,G,PD,HX,HY,IP,FZ,SCR,N,M,A,B,EP)
      DIMENSION U(N,M,0:1),V(N,M,0:1),F(N,M,0:1),G(N,M,0:1),
     1          PD(N,M),HX(N,M),HY(N,M),IP(N),FZ(N),SCR(N)
      COMMON/DISPL/ TIME,USCAL,DX,DY,H00,HMTS,DT
c     PARAMETER(NLON=130,NLAT=64)
c     PARAMETER(NLON=258,NLAT=128)
      PARAMETER(NLON=514,NLAT=256)
      DIMENSION UU(NLON,NLAT),VV(NLON,NLAT)
     
      AMP=0.
      IORT=2
      CF=.5
      C1=A*CF
      C2=B*CF
      C1H=C1*.5
      C2H=C2*.5
      NM=N*M
COMPUTE V+R/PHI*DT FIELD FOR LAGRANGIAN ESTIMATES
      DO 100 I=1,NM
      F(I,1,1)=U(I,1,0)+CF*F(I,1,0)/AMAX1(PD(I,1),EP)
  100 G(I,1,1)=V(I,1,0)+CF*G(I,1,0)/AMAX1(PD(I,1),EP)

COMPUTE U AND V TO FIRST ORDER
      DO 1 I=2,N-1
      DO 1 J=2,M-1
      ALFA=U(I,J,0)/HX(I,J)*C1
      BETA=V(I,J,0)/HY(I,J)*C2
      U(I,J,1)=F(I,J,1)-AMAX1(0.,ALFA)*(F(I,J,1)-F(I-1,J,1))
     1                 -AMIN1(0.,ALFA)*(F(I+1,J,1)-F(I,J,1))
     2                 -AMAX1(0.,BETA)*(F(I,J,1)-F(I,J-1,1))
     2                 -AMIN1(0.,BETA)*(F(I,J+1,1)-F(I,J,1))
      V(I,J,1)=G(I,J,1)-AMAX1(0.,ALFA)*(G(I,J,1)-G(I-1,J,1))
     1                 -AMIN1(0.,ALFA)*(G(I+1,J,1)-G(I,J,1))
     2                 -AMAX1(0.,BETA)*(G(I,J,1)-G(I,J-1,1))
     3                 -AMIN1(0.,BETA)*(G(I,J+1,1)-G(I,J,1))
    1 CONTINUE
      DO 2 I=2,N-1
      ALF1=U(I,1,0)/HX(I,1)*C1
      BET1=V(I,1,0)/HY(I,1)*C2
      U(I,1,1)=F(I,1,1)-AMAX1(0.,ALF1)*(F(I,1,1)-F(I-1,1,1))
     1                 -AMIN1(0.,ALF1)*(F(I+1,1,1)-F(I,1,1))
     2                 -AMAX1(0.,BET1)*(F(I,1,1)-F(IP(I),1,1))
     2                 -AMIN1(0.,BET1)*(F(I,2,1)-F(I,1,1))
      V(I,1,1)=G(I,1,1)-AMAX1(0.,ALF1)*(G(I,1,1)-G(I-1,1,1))
     1                 -AMIN1(0.,ALF1)*(G(I+1,1,1)-G(I,1,1))
     2                 -AMAX1(0.,BET1)*(G(I,1,1)+G(IP(I),1,1))
     3                 -AMIN1(0.,BET1)*(G(I,2,1)-G(I,1,1))
      ALFM=U(I,M,0)/HX(I,M)*C1
      BETM=V(I,M,0)/HY(I,M)*C2
      U(I,M,1)=F(I,M,1)-AMAX1(0.,ALFM)*(F(I,M,1)-F(I-1,M,1))
     1                 -AMIN1(0.,ALFM)*(F(I+1,M,1)-F(I,M,1))
     2                 -AMAX1(0.,BETM)*(F(I,M,1)-F(I,M-1,1))
     2                 -AMIN1(0.,BETM)*(F(IP(I),M,1)-F(I,M,1))
      V(I,M,1)=G(I,M,1)-AMAX1(0.,ALFM)*(G(I,M,1)-G(I-1,M,1))
     1                 -AMIN1(0.,ALFM)*(G(I+1,M,1)-G(I,M,1))
     2                 -AMAX1(0.,BETM)*(G(I,M,1)-G(I,M-1,1))
     3                 +AMIN1(0.,BETM)*(G(IP(I),M,1)+G(I,M,1))
    2 CONTINUE
      CALL XBC(U(1,1,1),N,M)
      CALL XBC(V(1,1,1),N,M)

      IF(IORT.EQ.2) THEN
COMPUTE U AND V TO SEMI-SECOND ORDER
      DO I=1,NM
        UU(I,1)=0.5*(U(I,1,0)+U(I,1,1))
        VV(I,1)=0.5*(V(I,1,0)+V(I,1,1))
      ENDDO
      DO 10 I=2,N-1
      DO 10 J=2,M-1
      ALFA=UU(I,J)/HX(I,J)*C1
      BETA=VV(I,J)/HY(I,J)*C2
      U(I,J,1)=F(I,J,1)-AMAX1(0.,ALFA)*(F(I,J,1)-F(I-1,J,1))
     1                 -AMIN1(0.,ALFA)*(F(I+1,J,1)-F(I,J,1))
     2                 -AMAX1(0.,BETA)*(F(I,J,1)-F(I,J-1,1))
     2                 -AMIN1(0.,BETA)*(F(I,J+1,1)-F(I,J,1))
      V(I,J,1)=G(I,J,1)-AMAX1(0.,ALFA)*(G(I,J,1)-G(I-1,J,1))
     1                 -AMIN1(0.,ALFA)*(G(I+1,J,1)-G(I,J,1))
     2                 -AMAX1(0.,BETA)*(G(I,J,1)-G(I,J-1,1))
     3                 -AMIN1(0.,BETA)*(G(I,J+1,1)-G(I,J,1))
   10 CONTINUE
      DO 20 I=2,N-1
      ALF1=UU(I,1)/HX(I,1)*C1
      BET1=VV(I,1)/HY(I,1)*C2
      U(I,1,1)=F(I,1,1)-AMAX1(0.,ALF1)*(F(I,1,1)-F(I-1,1,1))
     1                 -AMIN1(0.,ALF1)*(F(I+1,1,1)-F(I,1,1))
     2                 -AMAX1(0.,BET1)*(F(I,1,1)-F(IP(I),1,1))
     2                 -AMIN1(0.,BET1)*(F(I,2,1)-F(I,1,1))
      V(I,1,1)=G(I,1,1)-AMAX1(0.,ALF1)*(G(I,1,1)-G(I-1,1,1))
     1                 -AMIN1(0.,ALF1)*(G(I+1,1,1)-G(I,1,1))
     2                 -AMAX1(0.,BET1)*(G(I,1,1)+G(IP(I),1,1))
     3                 -AMIN1(0.,BET1)*(G(I,2,1)-G(I,1,1))
      ALFM=UU(I,M)/HX(I,M)*C1
      BETM=VV(I,M)/HY(I,M)*C2
      U(I,M,1)=F(I,M,1)-AMAX1(0.,ALFM)*(F(I,M,1)-F(I-1,M,1))
     1                 -AMIN1(0.,ALFM)*(F(I+1,M,1)-F(I,M,1))
     2                 -AMAX1(0.,BETM)*(F(I,M,1)-F(I,M-1,1))
     2                 -AMIN1(0.,BETM)*(F(IP(I),M,1)-F(I,M,1))
      V(I,M,1)=G(I,M,1)-AMAX1(0.,ALFM)*(G(I,M,1)-G(I-1,M,1))
     1                 -AMIN1(0.,ALFM)*(G(I+1,M,1)-G(I,M,1))
     2                 -AMAX1(0.,BETM)*(G(I,M,1)-G(I,M-1,1))
     3                 +AMIN1(0.,BETM)*(G(IP(I),M,1)+G(I,M,1))
   20 CONTINUE
      CALL XBC(U(1,1,1),N,M)
      CALL XBC(V(1,1,1),N,M)
      ENDIF

      CALL POLARFILT(U(1,1,1),HX,HY,FZ,SCR,N,M)
      CALL POLARFILT(V(1,1,1),HX,HY,FZ,SCR,N,M)
      RETURN
      END
 
      SUBROUTINE POLARFILT(FLD,HX,HY,FZ,SCR,N,M)
      DIMENSION FLD(N,M),HX(N,M),HY(N,M),FZ(N),SCR(N)

c      KMX=1   !128x64
c      KMX=2   !256x128
       KMX=4   !512x256
       IF(KMX.GT.0) THEN
       DO J=1,M
        AMP=(1.-HX(1,J)/HY(1,J))**2
        DO I=1,N
         FZ(I)=FLD(I,J)
        ENDDO
        CALL FILTRX(FZ(1),SCR(1),AMP,N,KMX)
        DO I=1,N
         FLD(I,J)=FZ(I)
        ENDDO
       ENDDO
       ENDIF

      RETURN
      END

      SUBROUTINE FILTRX(X,Y,AMP,N,KMX)
      DIMENSION X(N),Y(N)

      DO K=1,KMX

      DO I=1,N
c     IP=I+K-(I+K-1)/N*(N-2)
c     IM=I-K+(N-(I-K))/N*(N-2)
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

      RETURN
      END

      SUBROUTINE POLARABS(ALP,Y,DT,N,M)
      DIMENSION ALP(N,M),Y(M)

      DO J=1,M
      DO I=1,N
       ALP(I,J)=DT*0.
      ENDDO
      ENDDO

      eps=1.e-10
      pi=acos(-1.)
      abswidth=pi/64.*3
      atau=2.*DT
!     atau=200.*2.*DT
!      atau=2.e10
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
        y1 = amax1(0.,-y0+absy0_south)/abswidth
        y2 = amax1(0., y0-absy0_north)/abswidth
!       if (y0 > 0.) y1 = exp(-(ymax-y0)/abswidth)
!       if (y0 < 0.) y2 = exp(-(y0-ymin)/abswidth)
        do i=1,n
          alp(i,j) = DT*alpha*((y1**2+y2**2)/(y1+y2+eps))
        enddo
      endif
      enddo

      RETURN
      END

      subroutine compint(u,v)
c     PARAMETER(NLON=130,NLAT=64)
c     PARAMETER(NLON=258,NLAT=128)
      PARAMETER(NLON=514,NLAT=256)
      parameter(n=nlon,m=nlat)
      common/ipind/ ip(n)
      real u(n,m),v(n,m+1),apx(n,m,2),apy(n,m+1,2)

      ntr=4
      bt=1.

      do j=1,m
      do i=1,n
      apx(i,j,1)=u(i,j)
      enddo
      enddo
      do j=1,m+1
      do i=1,n
      apy(i,j,1)=v(i,j)
      enddo
      enddo

      DO ITR=1,NTR

      call lap1(apx,n, m ,1.)
       do j=1,m
        do i=1,n
         apx(i,j,1)=apx(i,j,1)+bt*apx(i,j,2)
       enddo
      enddo

      call lap2(apy,ip,n,m+1,1.)
       do j=1,m+1
        do i=1,n
         apy(i,j,1)=apy(i,j,1)+bt*apy(i,j,2)
        enddo
       enddo
!       do i=1,n
!        apy(i, 1 ,1)=0.
!        apy(i,m+1,1)=0.
!       enddo
      ENDDO
       do j=1,m
        do i=1,n
!        u(i,j)=apx(i,j,1)
         u(i,1)=apx(i,1,1)
         u(i,2)=apx(i,2,1)
         u(i,3)=apx(i,3,1)
         u(i,m)=apx(i,m,1)
         u(i,m-1)=apx(i,m-1,1)
         u(i,m-2)=apx(i,m-2,1)
        enddo
       enddo
!      do j=1,m+1
!       do i=1,n
!        v(i,j)=apy(i,j,1)
!       enddo
!      enddo

      return
      end

      subroutine lap1(hx,n,m,aa)
      dimension hx(n,m,2)

      do j=1,m
      do i=2,n-1
      hx(i,j,2)=0.25*(hx(i+1,j,1)-2.*aa*hx(i,j,1)+hx(i-1,j,1))
      enddo
      hx(1,j,2)=hx(n-1,j,2)
      hx(n,j,2)=hx( 2 ,j,2)
      enddo

      return
      end

      subroutine lap2(hy,ip,n,m,aa)
      dimension hy(n,m,2),ip(n)

      do i=1,n
      do j=2,m-1
      hy(i,j,2)=0.25*(hy(i,j+1,1)-2.*aa*hy(i,j,1)+hy(i,j-1,1))
      enddo
      hy(i,1,2)=0.25*(hy(i,2,1)-2.*aa*hy(i,1,1)-hy(ip(i),1,1))
      hy(i,m,2)=0.25*(-hy(ip(i),m,1)-2.*aa*hy(i,m,1)+hy(i,m-1,1))
!     hy(i,1,2)= 0.25*(hy(i,2,1)-aa*hy(i,1,1))
!     hy(i,m,2)=-0.25*(hy(i,m,1)-aa*hy(i,m-1,1))
      enddo

      return
      end

      SUBROUTINE MPDATT(U1,U2,X,H,N,M,IORD0,ISOR,NONOS,IDIV,IBC)
c     PARAMETER(N1=130,N2=64)
c     PARAMETER(N1=258,N2=128)
      PARAMETER(N1=514,N2=256)
      PARAMETER(MPG=1,MPA=1-MPG)
      DIMENSION U1(N,M),U2(N,M+1),X(N,M),H(N,M)
      DIMENSION V1(N1,N2),V2(N1,N2+1),F1(N1,N2),F2(N1,N2+1)
     *         ,CP(N1,N2),CN(N1,N2)   
      REAL MX(N1,N2),MN(N1,N2)
      COMMON/IPIND/ IP(N1)
      common/mpfil/ liner,mpfl
      DATA EP/1.E-10/
C      
      PP(Y)= AMAX1(0.,Y)
      PN(Y)=-AMIN1(0.,Y)
      DONOR(Y1,Y2,A)=PP(A)*Y1-PN(A)*Y2
              
      RAT2(Z1,Z2)=MPA*(ABS(Z2)-ABS(Z1))/(ABS(Z2)+ABS(Z1)+EP)    
c    .           *(.5+SIGN(.5,Z2*Z1))
     .           +MPG*(Z2-Z1)*.5
      RAT4(Z0,Z1,Z2,Z3)=MPA*(ABS(Z3)+ABS(Z2)-ABS(Z1)-ABS(Z0))
     .                 /(ABS(Z3)+ABS(Z2)+ABS(Z1)+ABS(Z0)+EP)
     .                 +MPG*(Z3+Z2-Z1-Z0)*.25

      VDYF(X1,X2,A,R)=(ABS(A)-A**2/R)*RAT2(X1,X2)
      VCORR(A,B,Y0,Y1,Y2,Y3,R)=-0.125*A*B/R*RAT4(Y0,Y1,Y2,Y3)
      VDIV1(A1,A2,A3,R)=0.25*A2*(A3-A1)/R
      VDIV2(A,B1,B2,B3,B4,R)=0.25*A*(B1+B2-B3-B4)/R

      VCOR31(A,X0,X1,X2,X3,R)=
     .     -(A -3.*ABS(A)*A/R+2.*A**3/R**2)/3.*RAT4(X1,X2,X0,X3)
      VCOR32(A,B,Y0,Y1,Y2,Y3,R)=
     .             0.25*B/R*(ABS(A)-2.*A**2/R)*RAT4(Y0,Y1,Y2,Y3)
C
      IORD=IORD0
      IF(ISOR.EQ.3) IORD=MAX0(IORD,3)
      if(LINER.EQ.1) IORD=1
      IF(MPG.EQ.1) IORD=MIN0(IORD,2)
CC
CC      DO 287 I=1,N1
CC  287 IP(I)=MOD(I+(N1-2)/2-1,N1-2)+1
CC
      DO 1 J=1,N2
      DO 1 I=1,N1
    1 V1(I,J)=U1(I,J)
      DO 2 J=1,N2+1
      DO 2 I=1,N1
    2 V2(I,J)=U2(I,J)

C                 
      IF(NONOS.EQ.1) THEN           
      DO 400 J=2,N2-1
      DO 400 I=2,N1-1
      MX(I,J)=AMAX1(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1))
  400 MN(I,J)=AMIN1(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1))
      DO 398 I=2,N1-1
      MX(I,1)=AMAX1(X(I-1,1),X(I,1),X(I+1,1),IBC*X(IP(I),1),X(I,2))
      MN(I,1)=AMIN1(X(I-1,1),X(I,1),X(I+1,1),IBC*X(IP(I),1),X(I,2))
      MX(I,N2)=AMAX1(X(I-1,N2),X(I,N2),X(I+1,N2),X(I,N2-1),
     1                                       IBC*X(IP(I),N2))
  398 MN(I,N2)=AMIN1(X(I-1,N2),X(I,N2),X(I+1,N2),X(I,N2-1),
     1                                  IBC*X(IP(I),N2))
      CALL XBC(MX,N1,N2)
      CALL XBC(MN,N1,N2)
      ENDIF
C    
      C1=1.
      C2=0.
                        DO 3 K=1,IORD
COMPUTE DONOR-CELL FLUXES
      DO 330 J=1,N2
      DO 330 I=2,N1-1
  330 F1(I,J)=DONOR(C1*X(I-1,J)+C2,C1*X(I,J)+C2,V1(I,J))
      DO 331 J=2,N2
      DO 331 I=2,N1-1
  331 F2(I,J)=DONOR(C1*X(I,J-1)+C2,C1*X(I,J)+C2,V2(I,J))
      DO 332 I=2,N1-1
      F2(I,N2+1)=DONOR(C1*X(I,N2)+C2,C1*IBC*X(IP(I),N2)+C2,V2(I,N2+1))
  332 F2(I,1)=DONOR(C1*IBC*X(IP(I),1)+C2,C1*X(I,1)+C2,V2(I,1))
      CALL XBC(F1,N1,N2)
      CALL XBC(F2,N1,N2+1)
COMPUTE NEW UPWIND-SOLUTION
      DO 333 J=1,N2
      DO 333 I=2,N1-1
  333 X(I,J)=X(I,J)-(F1(I+1,J)-F1(I,J)+F2(I,J+1)-F2(I,J))/H(I,J)
      CALL XBC(X,N1,N2)
C
      IF(K.EQ.IORD) GO TO 6
      C1=FLOAT(MPA)
      C2=FLOAT(MPG)
CONVERT VELOCITIES TO LOCAL STORAGE
      DO 49 J=1,N2
      DO 49 I=1,N1
   49 F1(I,J)=V1(I,J)
      DO 50 J=1,N2+1
      DO 50 I=1,N1
   50 F2(I,J)=V2(I,J) 

CALCULATE PSEUDO VELOCITIES      
COMPUTE FIRST DIRECTION    
      DO 51 J=2,N2-1
      DO 51 I=2,N1-1
   51 V1(I,J)=VDYF(X(I-1,J),X(I,J),F1(I,J),.5*(H(I-1,J)+H(I,J)))
     * +VCORR(F1(I,J), F2(I-1,J)+F2(I-1,J+1)+F2(I,J+1)+F2(I,J),
     *                 X(I-1,J-1),X(I,J-1),X(I-1,J+1),X(I,J+1),
     *                 .5*(H(I-1,J)+H(I,J)))
COMPUTE B.C IN Y DIRECTION
      DO 510 I=2,N1-1
      V1(I,1)=VDYF(X(I-1,1),X(I,1),F1(I,1),.5*(H(I-1,1)+H(I,1)))
     * +VCORR(F1(I,1), F2(I-1,1)+F2(I-1,2)+F2(I,2)+F2(I,1),
     *   IBC*X(IP(I-1),1),IBC*X(IP(I),1),X(I-1,2),X(I,2),
     *                 .5*(H(I-1,1)+H(I,1)))
  510 V1(I,N2)=VDYF(X(I-1,N2),X(I,N2),F1(I,N2),.5*(H(I-1,N2)+H(I,N2)))
     * +VCORR(F1(I,N2), F2(I-1,N2)+F2(I-1,N2+1)+F2(I,N2+1)+F2(I,N2),
     *        X(I-1,N2-1),X(I,N2-1),IBC*X(IP(I-1),N2),IBC*X(IP(I),N2),
     *                 .5*(H(I-1,N2)+H(I,N2)))        
      IF(IDIV.EQ.1) THEN
COMPUTE FLOW-DIVERGENCE CORRECTION
      DO 5101 J=1,N2
      DO 5101 I=2,N1-1
      V1D=-VDIV1(F1(I-1,J),F1(I,J),F1(I+1,J),.5*(H(I-1,J)+H(I,J)))
     *    -VDIV2(F1(I,J),F2(I-1,J+1),F2(I,J+1),F2(I-1,J),F2(I,J),
     *                 .5*(H(I-1,J)+H(I,J)))
 5101 V1(I,J)=V1(I,J)+MPG*(PP(V1D)*X(I-1,J)-PN(V1D)*X(I,J))+MPA*V1D
      ENDIF   

COMPUTE SECOND DIRECTION
      DO 52 J=2,N2
      DO 52 I=2,N1-1
   52 V2(I,J)=VDYF(X(I,J-1),X(I,J),F2(I,J),.5*(H(I,J-1)+H(I,J)))
     * +VCORR(F2(I,J), F1(I,J-1)+F1(I,J)+F1(I+1,J)+F1(I+1,J-1),
     *                 X(I-1,J-1),X(I-1,J),X(I+1,J-1),X(I+1,J),
     *                 .5*(H(I,J-1)+H(I,J)))
COMPUTE B.C IN Y-DIRECTION
      DO 520 I=2,N1-1
      V2(I,1)=VDYF(IBC*X(IP(I),1),X(I,1),F2(I,1),.5*(H(IP(I),1)+H(I,1)))
     * +VCORR(F2(I,1), F1(IP(I),1)+F1(I,1)+F1(I+1,1)+F1(IP(I+1),1),
     *         IBC*X(IP(I-1),1),X(I-1,1),X(I+1,1),IBC*X(IP(I+1),1),
     *                 .5*(H(IP(I),1)+H(I,1)))
      V2(I,N2+1)=VDYF(X(I,N2),IBC*X(IP(I),N2),F2(I,N2+1),
     1                        .5*(H(I,N2)+H(IP(I),N2)))
     *+VCORR(F2(I,N2+1),F1(I,N2)+F1(IP(I),N2)+F1(IP(I+1),N2)+F1(I+1,N2),
     *     X(I-1,N2),IBC*X(IP(I-1),N2),X(I+1,N2),IBC*X(IP(I+1),N2),
     *                        .5*(H(I,N2)+H(IP(I),N2)))
  520 CONTINUE
      IF(IDIV.EQ.1) THEN
      DO 5201 J=2,N2
      DO 5201 I=2,N1-1
      V2D=-VDIV1(F2(I,J-1),F2(I,J),F2(I,J+1),.5*(H(I,J-1)+H(I,J)))
     *    -VDIV2(F2(I,J),F1(I+1,J-1),F1(I+1,J),F1(I,J-1),F1(I,J),
     *                 .5*(H(I,J-1)+H(I,J)))
 5201 V2(I,J)=V2(I,J)+MPG*(PP(V2D)*X(I,J-1)-PN(V2D)*X(I,J))+MPA*V2D

      DO 5202 I=2,N1-1
      V2D1=-VDIV1(-F2(IP(I),1),F2(I,1),F2(I,2),.5*(H(IP(I),1)+H(I,1)))
     *     -VDIV2( F2(I,1),F1(IP(I+1),1),F1(I+1,1),F1(IP(I),1),F1(I,1),
     *                 .5*(H(IP(I),1)+H(I,1)))
      V2(I,1)=V2(I,1)
     .       +MPG*(PP(V2D1)*IBC*X(IP(I),1)-PN(V2D1)*X(I,1))+MPA*V2D1
      V2DN=-VDIV1(F2(I,N2),F2(I,N2+1),-F2(IP(I),N2+1)
     *                      ,.5*(H(I,N2)+H(IP(I),N2)))
     *     -VDIV2(F2(I,N2+1),F1(I+1,N2),F1(IP(I+1),N2),
     *            F1(I,N2),F1(IP(I),N2),.5*(H(I,N2)+H(IP(I),N2))) 
 5202 V2(I,N2+1)=V2(I,N2+1)
     .          +MPG*(PP(V2DN)*X(I,N2)-PN(V2DN)*IBC*X(IP(I),N2))+MPA*V2DN
      ENDIF
C
C THIRD ORDER CORRECTION
      IF(ISOR.EQ.3) THEN
C  FIRST DIRECTION
      DO 61 J=1,N2
      DO 610 I=3,N1-1
  610 V1(I,J)=V1(I,J)     +VCOR31(F1(I,J),
     1        X(I-2,J),X(I-1,J),X(I,J),X(I+1,J),.5*(H(I-1,J)+H(I,J)))
   61 V1(2,J)=V1(2,J)     +VCOR31(F1(2,J),
     1        X(N1-2,J),X(1,J),X(2,J),X(3,J),.5*(H(1,J)+H(2,J)))
      DO 62 J=2,N2-1
      DO 62 I=2,N1-1               
   62 V1(I,J)=V1(I,J)
     1 +VCOR32(F1(I,J),F2(I-1,J)+F2(I-1,J+1)+F2(I,J+1)+F2(I,J),
     *           X(I,J-1),X(I-1,J+1),X(I-1,J-1),X(I,J+1),
     *                   .5*(H(I-1,J)+H(I,J)))
C B.C. FOLLOW
      DO 620 I=2,N1-1
      V1(I,1)=V1(I,1)
     1 +VCOR32(F1(I,1),F2(I-1,1)+F2(I-1,2)+F2(I,2)+F2(I,1),
     *       IBC*X(IP(I),1),X(I-1,2),X(I,2),IBC*X(IP(I-1),1),
     *                   .5*(H(I-1,1)+H(I,1)))
  620 V1(I,N2)=V1(I,N2)
     1 +VCOR32(F1(I,N2),F2(I-1,N2)+F2(I-1,N2+1)+F2(I,N2+1)+F2(I,N2),
     *      X(I,N2-1),IBC*X(IP(I-1),N2),IBC*X(IP(I),N2),X(I-1,N2-1),
     *                   .5*(H(I-1,N2)+H(I,N2)))
      DO 621 J=1,N2
      V1(1,J)=V1(N1-1,J)
  621 V1(N1,J)=V1(2,J)
C
C  SECOND DIRECTION
      DO 63 J=3,N2-1
      DO 63 I=2,N1-1
   63 V2(I,J)=V2(I,J)     +VCOR31(F2(I,J),
     1        X(I,J-2),X(I,J-1),X(I,J),X(I,J+1),.5*(H(I,J-1)+H(I,J)))  
      DO 630 I=2,N1-1
      V2(I,1)=V2(I,1)     +VCOR31(F2(I,1),IBC*X(IP(I),2),
     1         IBC*X(IP(I),1),X(I,1),X(I,2),.5*(H(IP(I),1)+H(I,1)))  
      V2(I,2)=V2(I,2)     +VCOR31(F2(I,2),
     1      IBC*X(IP(I),1),X(I,1),X(I,2),X(I,3),.5*(H(I,1)+H(I,2)))  
      V2(I,N2)=V2(I,N2)     +VCOR31(F2(I,N2),X(I,N2-2),
     1    X(I,N2-1),X(I,N2),IBC*X(IP(I),N2),.5*(H(I,N2-1)+H(I,N2)))  
      V2(I,N2+1)=V2(I,N2+1) +VCOR31(F2(I,N2+1), X(I,N2-1),X(I,N2),
     1   IBC*X(IP(I),N2),IBC*X(IP(I),N2-1),.5*(H(I,N2)+H(IP(I),N2)))  
  630 CONTINUE
      DO 64 J=2,N2
      DO 64 I=2,N1-1
   64 V2(I,J)=V2(I,J)
     1 +VCOR32(F2(I,J),F1(I,J-1)+F1(I+1,J-1)+F1(I+1,J)+F1(I,J),
     *            X(I+1,J-1),X(I-1,J),X(I-1,J-1),X(I+1,J),
     *                  .5*(H(I,J-1)+H(I,J)))
      ENDIF

CALL B.C IN X DIRECTION
      CALL XBC(V1,N1,N2)
      CALL XBC(V2,N1,N2+1)
C
C     DO 23 J=1,N2
C     DO 23 I=1,N1
C  23 V1(I,J)=SIGN(1.,V1(I,J))*AMIN1(ABS(U1(I,J)),ABS(V1(I,J)))
C     DO 24 J=1,N2+1                     
C     DO 24 I=1,N1
C  24 V2(I,J)=SIGN(1.,V2(I,J))*AMIN1(ABS(U2(I,J)),ABS(V2(I,J)))
C
      IF(NONOS.EQ.0) GO TO 3
C                 NON-OSCILLATORY OPTION
      DO 401 J=2,N2-1
      DO 401 I=2,N1-1
      MX(I,J)=AMAX1(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1),MX(I,J))
  401 MN(I,J)=AMIN1(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1),MN(I,J))
      DO 3981 I=2,N1-1
      MX(I,1)=AMAX1(X(I-1,1),X(I,1),X(I+1,1),IBC*X(IP(I),1),
     1                                        X(I,2),MX(I,1))
      MN(I,1)=AMIN1(X(I-1,1),X(I,1),X(I+1,1),IBC*X(IP(I),1),
     1                                        X(I,2),MN(I,1))
      MX(I,N2)=AMAX1(X(I-1,N2),X(I,N2),X(I+1,N2),X(I,N2-1),
     1                              IBC*X(IP(I),N2),MX(I,N2))
 3981 MN(I,N2)=AMIN1(X(I-1,N2),X(I,N2),X(I+1,N2),X(I,N2-1),
     1                              IBC*X(IP(I),N2),MN(I,N2))
      CALL XBC(MX,N1,N2)
      CALL XBC(MN,N1,N2)
C
      DO 4021 J=1,N2
      DO 4021 I=2,N1-1
 4021 F1(I,J)=DONOR(C1*X(I-1,J)+C2,C1*X(I,J)+C2,V1(I,J))
      DO 403 J=2,N2
      DO 403 I=2,N1-1
  403 F2(I,J)=DONOR(C1*X(I,J-1)+C2,C1*X(I,J)+C2,V2(I,J))
      DO 410 I=2,N1-1    
      F2(I,N2+1)=DONOR(C1*X(I,N2)+C2,C1*IBC*X(IP(I),N2)+C2,V2(I,N2+1))
  410 F2(I,1)=DONOR(C1*IBC*X(IP(I),1)+C2,C1*X(I,1)+C2,V2(I,1))
      CALL XBC(F1,N1,N2)
      CALL XBC(F2,N1,N2+1)

      DO 404 J=1,N2
      DO 404 I=2,N1-1
      CP(I,J)=(MX(I,J)-X(I,J))*H(I,J)/
     1(PN(F1(I+1,J))+PP(F1(I,J))+PN(F2(I,J+1))+PP(F2(I,J))+EP)
      CN(I,J)=(X(I,J)-MN(I,J))*H(I,J)/
     1(PP(F1(I+1,J))+PN(F1(I,J))+PP(F2(I,J+1))+PN(F2(I,J))+EP)
  404 CONTINUE
      CALL XBC(CP,N1,N2)
      CALL XBC(CN,N1,N2)

      IF(MPG.EQ.0) THEN
      DO 405 J=2,N2
      DO 405 I=2,N1-1
      V1(I,J)=PP(V1(I,J))*
     1  ( AMIN1(1.,CP(I,J),CN(I-1,J))*PP(SIGN(1., X(I-1,J)))
     1   +AMIN1(1.,CP(I-1,J),CN(I,J))*PP(SIGN(1.,-X(I-1,J))) )
     2       -PN(V1(I,J))*
     2  ( AMIN1(1.,CP(I-1,J),CN(I,J))*PP(SIGN(1., X(I ,J )))
     2   +AMIN1(1.,CP(I,J),CN(I-1,J))*PP(SIGN(1.,-X(I ,J ))) )
  405 V2(I,J)=PP(V2(I,J))*
     1  ( AMIN1(1.,CP(I,J),CN(I,J-1))*PP(SIGN(1., X(I,J-1)))
     1   +AMIN1(1.,CP(I,J-1),CN(I,J))*PP(SIGN(1.,-X(I,J-1))) )
     1       -PN(V2(I,J))*
     2  ( AMIN1(1.,CP(I,J-1),CN(I,J))*PP(SIGN(1., X(I ,J )))
     2   +AMIN1(1.,CP(I,J),CN(I,J-1))*PP(SIGN(1.,-X(I ,J ))) )
C B.C. FOLLOW
      DO 4051 I=2,N1-1
      V2(I,1)=PP(V2(I,1))*
     1  ( AMIN1(1.,CP(I,1),CN(IP(I),1))*PP(SIGN(1., IBC*X(IP(I),1)))
     1   +AMIN1(1.,CP(IP(I),1),CN(I,1))*PP(SIGN(1.,-IBC*X(IP(I),1))) )
     1       -PN(V2(I,1))*
     2  ( AMIN1(1.,CP(IP(I),1),CN(I,1))*PP(SIGN(1., X(I ,1 )))
     2   +AMIN1(1.,CP(I,1),CN(IP(I),1))*PP(SIGN(1.,-X(I ,1 ))) )
      V2(I,N2+1)=PP(V2(I,N2+1))*
     1  ( AMIN1(1.,CP(IP(I),N2),CN(I,N2))*PP(SIGN(1., X(I,N2)))
     1   +AMIN1(1.,CP(I,N2),CN(IP(I),N2))*PP(SIGN(1.,-X(I,N2))) )
     1       -PN(V2(I,N2+1))*
     2( AMIN1(1.,CP(I,N2),CN(IP(I),N2))*PP(SIGN(1., IBC*X(IP(I),N2)))
     2 +AMIN1(1.,CP(IP(I),N2),CN(I,N2))*PP(SIGN(1.,-IBC*X(IP(I),N2))) )
 4051 V1(I,1)=PP(V1(I,1))*
     1  ( AMIN1(1.,CP(I,1),CN(I-1,1))*PP(SIGN(1., X(I-1,1)))
     1   +AMIN1(1.,CP(I-1,1),CN(I,1))*PP(SIGN(1.,-X(I-1,1))) )
     2       -PN(V1(I,1))*
     2  ( AMIN1(1.,CP(I-1,1),CN(I,1))*PP(SIGN(1., X(I ,1 )))
     2   +AMIN1(1.,CP(I,1),CN(I-1,1))*PP(SIGN(1.,-X(I ,1 ))) )
      ELSE
      DO 505 J=2,N2
      DO 505 I=2,N1-1
      V1(I,J)=PP(V1(I,J))*AMIN1(1.,CP(I,J),CN(I-1,J))
     2       -PN(V1(I,J))*AMIN1(1.,CP(I-1,J),CN(I,J))
  505 V2(I,J)=PP(V2(I,J))*AMIN1(1.,CP(I,J),CN(I,J-1))
     1       -PN(V2(I,J))*AMIN1(1.,CP(I,J-1),CN(I,J))
C B.C. FOLLOW
      DO 5051 I=2,N1-1
      V2(I,1)=PP(V2(I,1))*AMIN1(1.,CP(I,1),CN(IP(I),1))
     1       -PN(V2(I,1))*AMIN1(1.,CP(IP(I),1),CN(I,1))
      V2(I,N2+1)=PP(V2(I,N2+1))*AMIN1(1.,CP(IP(I),N2),CN(I,N2))
     1          -PN(V2(I,N2+1))*AMIN1(1.,CP(I,N2),CN(IP(I),N2))
 5051 V1(I,1)=PP(V1(I,1))*AMIN1(1.,CP(I,1),CN(I-1,1))
     2       -PN(V1(I,1))*AMIN1(1.,CP(I-1,1),CN(I,1))
      ENDIF
      CALL XBC(V1,N1,N2)
      CALL XBC(V2,N1,N2+1)
C
C         END OF NONOSCILLATORY OPTION
C
    3                      CONTINUE
    6 CONTINUE
      RETURN
      END   

      SUBROUTINE XBC(X,N,M)
      DIMENSION X(N,M)
      DO 1 J=1,M
      X(1,J)=X(N-1,J)
    1 X(N,J)=X(2,J)
      RETURN
      END

CDIR$ NOLIST

      SUBROUTINE DIAGNOS(U,V,PD,PT,COR,HX,HY,S,X,XX,Y,IP,PI2,KT,IFLG)
c     PARAMETER(N=130,M=64)
c     PARAMETER(N=258,M=128)
      PARAMETER(N=514,M=256)
      DIMENSION U(N,M),V(N,M),PD(N,M),PT(N,M),COR(N,M),
     1          HX(N,M),HY(N,M),S(N,M),X(N),XX(N),Y(M),IP(N)
      COMMON/INITCON/ U0(N,M),V0(N,M),PT0(N,M),PD0(N,M),P0(N,M)
      COMMON/INITVLS/ PTTRL2,PTTRLI,VLTRL2,VLTRLI,PDTRL2,PDTRLI
      COMMON/DISPL/ TIME,USCAL,DX,DY,H00,HMTS,DT
      COMMON/ITERO/ ERROR,EXITCND,NITER,NITSM,ICOUNT
      COMMON/ELLIPS/ EPS,ITR,ITMN,IPRC
      COMMON/COMP/ XLP,SUM0,SUMT
      DATA GR,ER/9.80616,6371.22E+03/
      DATA COURTMX/0./

      GC1=DT/DX
      GC2=DT/DY
      NM=N*M

      IF(IFLG.EQ.0) THEN

      COURTMX=0.
      NITER=0
      ERROR=0.
      SUM0=0.
      DO 10 J=1,M
      DO 10 I=2,N-1
   10 SUM0=SUM0+PD(I,J)*S(I,J)
      ENDIF

      TIMSI=KT*DT
      TIME=TIMSI/3600./24.
!     PRINT 299
! 299 FORMAT(1X,1H )
      PRINT 300,TIME,KT
  300 FORMAT(14X,5HTIME=,F7.2,6H DAYS;,5H  KT=,I5)

CHECK COURANT NUMBERS
      COUR1=0.
      COUR2=0.
      DO 20 I=1,NM
      DLI=SQRT((GC1/HX(I,1))**2+(GC2/HY(I,1))**2)
      COUR1=AMAX1(COUR1,DLI*SQRT(abs(PT(I,1))))
      COUR2=AMAX1(COUR2, GC1*ABS(U(I,1)/HX(I,1)) 
     1                  +GC2*ABS(V(I,1)/HY(I,1)) )
!     COUR2=AMAX1(COUR2,SQRT( (GC1*U(I,1)/HX(I,1))**2
!    1                       +(GC2*V(I,1)/HY(I,1))**2))
   20 CONTINUE
      COURTMX=amax1(COUR2,COURTMX)
      PRINT 301, COUR1,COUR2,COURTMX
  301 FORMAT(4X,'COUR1,COUR2,COURTMX:',3E11.4)

      SUMT=0.
      PDMX=-1.E30
      PDMN= 1.E30
      PDAV=0.
      DO 30 J=1,M
      DO 30 I=2,N-1
      PDMX=AMAX1(PDMX,PD(I,J))
      PDMN=AMIN1(PDMN,PD(I,J))
      PDAV=PDAV+PD(I,J)
   30 SUMT=SUMT+PD(I,J)*S(I,J)
      PDAV=PDAV/FLOAT(M*(N-2))
      SUMER=(SUMT-SUM0)/AMAX1(SUM0,1.E-30)
      PRINT 302, PDMX,PDMN,PDAV,SUMER
  302 FORMAT(4X,'PDMX,PDMN,PDAV:',3E11.4,' SUMER:',E11.4)

      NITAV=NITSM/MAX0(ICOUNT,1)
      PRINT 303, ERROR,NITER,NITAV
  303 FORMAT(4X,'ERROR:',E11.4,1X,'NITER, NITAV (GCR ITERATIONS):',2I4)

      CALL RHWT(U0,V0,PT0,PD0,P0,COR,S,X,XX,Y,N,M,GR,ER,TIMSI)

      PTERL2 = 0.
      PTERLI =-1.e15
      PDERL2 = 0.
      PDERLI =-1.e15
      QXERL2 = 0.
      QXERLI =-1.e15
      QYERL2 = 0.
      QYERLI =-1.e15
      DO I=1,NM
        PTERL2 = PTERL2 + S(I,1)*(PT(I,1)-PT0(I,1))**2
        PTERLI = AMAX1(PTERLI, ABS(PT(I,1)-PT0(I,1)))
        PDERL2 = PDERL2 + S(I,1)*(PD(I,1)-PD0(I,1))**2
        PDERLI = AMAX1(PDERLI, ABS(PD(I,1)-PD0(I,1)))
        VLERL2 =VLERL2 + S(I,1)*( (U(I,1)-U0(I,1))**2
     &                           +(V(I,1)-V0(I,1))**2 )
        VLERLI =AMAX1(VLER1, (U(I,1)-U0(I,1))**2
     &                      +(V(I,1)-V0(I,1))**2 )
      ENDDO
      PTERL2 = SQRT(PTERL2/PTTRL2)
      PTERLI = PTERLI/PTTRLI
      PDERL2 = SQRT(PDERL2/PDTRL2)
      PDERLI = PDERLI/PDTRLI
      VLERL2 = SQRT(VLERL2/VLTRL2)
      VLERLI = SQRT(VLERLI/VLTRLI)
      PRINT 401, PTERL2,PTERLI
      PRINT 402, PDERL2,PDERLI
      PRINT 403, VLERL2,VLERLI
  401 FORMAT(1X,'PTER2,PTERI;',2E11.4)
  402 FORMAT(1X,'PDER2,PDERI;',2E11.4)
  403 FORMAT(1X,'VLER2,VLERI;',2E11.4)

      RETURN
      END


CDIR$ NOLIST
      subroutine initc(xlon, ylat, rho)
c
      pi = acos(-1.)
      anot = 1.0
c      rnot = 2*pi*7./128.
c      xcent = 3*pi/2.
      rnot = 2*pi*10./128.
      xcent = 2*pi/2.
      ycent = pi/2. *1.
c     
      gamm=1.
      rho =0.
c
      dist = 2.*sqrt( (cos(ylat)*sin( (xlon-xcent)/2 ) )**2
     a     + sin((ylat-ycent)/2)**2       )
      if (dist.le.rnot) rho = anot
     a     *(1.-gamm*dist/rnot)
     b                      + rho
c
      return
      end
      subroutine dep (xar, yar, xde, yde, rvel, beta, dt)
      real xar, yar, xde, yde
c     
c     distance in radians moved per time step
      xpofs = rvel * dt
c     
c     pole of velocity rotated by angle beta
c     
      cosalf = cos(beta)
      sinalf = sin(beta)
c     
      zz = sin(yar)
      yy = sin(xar)*cos(yar)
      xx = cos(xar)*cos(yar)
c     
      yyp = yy
      xxp = xx*cosalf + zz*sinalf
      zzp = zz*cosalf - xx*sinalf
c     
      ypm = asin(zzp)
      xpm = atan2(yyp,xxp)
c     
      xpm = xpm + xpofs
c     
      zzp = sin(ypm)
      yyp = sin(xpm)*cos(ypm)
      xxp = cos(xpm)*cos(ypm)
c     
      yy = yyp
      xx = xxp*cosalf - zzp*sinalf
      zz = zzp*cosalf + xxp*sinalf
c     
      yde = asin(zz)
      xde = atan2(yy,xx)
c     
      pi = 4*atan(1.0)
c     
      if(yde.lt.-pi/2) then
         yde = -pi-yde
         xde = xde + pi/2.
      write(6,*) '-pi/2',xde,yde,xar,yar
      endif
      if(yde.gt.pi/2)  then
         yde = pi-yde
         xde = xde + pi/2.
      write(6,*) '+pi/2',xde,yde,xar,yar
      endif
      if(xde.lt.0.) then
         xde = xde + 2*pi
      endif
      if(xde.gt.2*pi) then
         xde = xde - 2*pi
      endif
c     
      return
      end   

