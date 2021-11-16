module implicit_functions_hp

!!!!!!!!!!!!!!! NEW FUNCTIONS !!!!!!!!!!!!!!!!!!!!!11
 IMPLICIT NONE                                                                                       
  SAVE                                                                                                
  INTEGER, PARAMETER :: HP=2                                                                                                    
  REAL(HP), PUBLIC :: RPE_0,RPE_1,RPE_05,RPE_025,RPE_0125,RPE_2,RPE_3,RPE_4 , rpe_100                                 
 
 contains                                                                                              
                                                                                                      
  SUBROUTINE rpenum_init(sbits)    
    
    INTEGER :: sbits 
    
    
    rpe_1 = (1.0)                                                                          
    rpe_2 = (2.0)   
    rpe_3 = (3.0)
    rpe_4 = (4.0)   
    rpe_100 = (100.0) 
    rpe_0125 = (0.125) 
    rpe_025 = (0.25)                                                                       
    rpe_05 = (0.5)                                                                         
    rpe_0 = (0.0)                                                                           
                                                                                                                                                                                                 
  END SUBROUTINE rpenum_init     
  
function DONOR(Y1,Y2,A)
implicit none
  REAL(HP) :: DONOR, Y1, Y2, A

  DONOR = max(rpe_0,A)*Y1+ min(rpe_0,A)*Y2

end function

             

function VDYF(X1,X2,A,R, EP)

implicit none

REAL(HP) :: VDYF, X1,X2,A,R, EP

  VDYF = (ABS(A)-A**2/R)*(ABS(X2)-ABS(X1))   &
     &          /(ABS(X2)+ABS(X1)+EP)
end function


function VCORR(A,B,Y1,Y2,R)

implicit none

REAL(HP) :: VCORR, A,B, Y1,Y2,R

  VCORR=-rpe_0125*A*B*Y1/(Y2*R)
end function


function VCOR31(A,X0,X1,X2,X3,R, EP)

implicit none

REAL(HP) :: VCOR31, A,X0,X1,X2,X3,R, EP

  VCOR31= -(A - rpe_3*ABS(A)*A/R+rpe_2*A**3/R**2)/rpe_3    &
     &       *(ABS(X0)+ABS(X3)-ABS(X1)-ABS(X2))     &
     &    /(ABS(X0)+ABS(X3)+ABS(X1)+ABS(X2)+EP)
end function


function VCOR32(A,B,Y1,Y2,R)

implicit none

  REAL(HP) :: VCOR32, A,B,Y1,Y2,R

  VCOR32= rpe_025*B/R*(ABS(A)-rpe_2*A**2/R)*Y1/Y2
end function


function VDIV1(A1,A2,A3,R)

IMPLICIT none
  REAL(HP) :: VDIV1, A1,A2,A3,R

  !VDIV1= rpe_025*A2*(A3-A1)/R
  VDIV1= rpe_025*(A2/R)*(A3-A1)
end function


function VDIV2(A,B1,B2,B3,B4,R)

IMPLICIT none
  REAL(HP) :: VDIV2, A,B1,B2,B3,B4,R

  !VDIV2= rpe_025*A*(B1+B2-B3-B4)/R
  VDIV2= rpe_025*(A/R)*(B1+B2-(B3+B4))
end function


function PP(Y)

implicit none
  REAL(HP) :: PP, Y

  PP = max(rpe_0,Y)
end function


function PN(Y)

implicit none
  REAL(HP) :: PN, Y

  PN = -min(rpe_0,Y) 
end function

function RAT2(Z1,Z2,MP, EP,x_scale)

implicit none
  REAL(HP) :: RAT2, Z1,Z2, EP,x_scale
  INTEGER :: MP

  RAT2=float(MP)*(ABS(Z2)-ABS(Z1))/(ABS(Z2)+ABS(Z1)+EP)*x_scale            &    
   !    .           *(.5+SIGN(.5,Z2*Z1))
    &           +(rpe_1-float(MP))*(Z2-Z1)*rpe_05
  !RAT2=float(MP)*(ABS(Z2)-ABS(Z1))/(ABS(Z2)+ABS(Z1)+EP)                  &    
  ! !    .           *(.5+SIGN(.5,Z2*Z1))
  !  &           +(rpe_1-float(MP))*(Z2-Z1)*rpe_05
end function


function RAT4(Z0,Z1,Z2,Z3,MP, EP,x_scale)

implicit none
  REAL(HP) :: RAT4, Z0,Z1,Z2,Z3, EP,x_scale
  INTEGER :: MP

  RAT4=float(MP)*(ABS(Z3)+ABS(Z2)-ABS(Z1)-ABS(Z0))                 &
    &           /(ABS(Z3)+ABS(Z2)+ABS(Z1)+ABS(Z0)+EP)*x_scale      &
    &           +(rpe_1-float(MP))*((Z3-Z1) +(Z2-Z0))*rpe_025
  !RAT4=float(MP)*(ABS(Z3)+ABS(Z2)-ABS(Z1)-ABS(Z0))                 &
  !  &           /(ABS(Z3)+ABS(Z2)+ABS(Z1)+ABS(Z0)+EP)      &
  !  &           +(rpe_1-float(MP))*(Z3+Z2-Z1-Z0)*rpe_025

end function

function VDYF_D(X1,X2,A,R,MP, EP,x_scale)

implicit none
  REAL(HP) :: VDYF_D, X1,X2,A,R, EP, AbsA, x_scale
  INTEGER :: MP
  AbsA=ABS(A)
  !VDYF_D=(ABS(A)-A**2/R)*RAT2(X1,X2,MP, EP)
  VDYF_D=AbsA*(rpe_1-AbsA/R)*RAT2(X1,X2,MP, EP,x_scale)

end function

function VDYF_D_JA(X1,X2,A,R,MP, EP,x_scale,xdiff)

implicit none
  REAL(HP) :: VDYF_D_JA, X1,X2,A,R, EP, AbsA, x_scale, xdiff
  INTEGER :: MP
  AbsA=ABS(A)
  !VDYF_D=(ABS(A)-A**2/R)*RAT2(X1,X2,MP, EP)
  VDYF_D_JA=AbsA*(rpe_1-AbsA/R)*rpe_05 *xdiff!*x_scale

end function

function VCORR_D(A,B,Y0,Y1,Y2,Y3,R,MP, EP,x_scale)

implicit none
  REAL(HP) :: VCORR_D, A,B,Y0,Y1,Y2,Y3,R, EP, x_scale
  INTEGER :: MP

  !VCORR_D=-rpe_0125*A*B/R*RAT4(Y0,Y1,Y2,Y3,MP,EP)
  VCORR_D=-rpe_0125*(A/R)*B*RAT4(Y0,Y1,Y2,Y3,MP,EP,x_scale)

end function

function VCORR_D_JA(A,B,Y0,Y1,Y2,Y3,R,MP, EP,x_scale)

implicit none
  REAL(HP) :: VCORR_D_JA, A,B,Y0,Y1,Y2,Y3,R, EP, x_scale
  INTEGER :: MP

  !VCORR_D=-rpe_0125*A*B/R*RAT4(Y0,Y1,Y2,Y3,MP,EP)
  VCORR_D_JA=-rpe_0125*(A/R)*B*((Y3+Y2) -(Y1+Y0))*rpe_025

end function

!function VCOR31_D(A,X0,X1,X2,X3,R,MP, EP)
!
!implicit none
!  REAL(HP) :: VCOR31_D, A,X0,X1,X2,X3,R, EP
!  INTEGER :: MP
!
!  VCOR31_D=        &
!  &     -(A -rpe_3*ABS(A)*A/R+rpe_2*A**3/R**2)/rpe_3*RAT4(X1,X2,X0,X3,MP, EP)
!
!end function
!
!
!function VCOR32_D(A,B,Y0,Y1,Y2,Y3,R,MP, EP)
!
!implicit none
!
!  REAL(HP) :: VCOR32_D, A,B,Y0,Y1,Y2,Y3,R, EP
!  INTEGER :: MP
!
!  VCOR32_D=                                              &
!  &     rpe_025*B/R*(ABS(A)-rpe_2*A**2/R)*RAT4(Y0,Y1,Y2,Y3,MP, EP)
!
!end function

end module implicit_functions_hp
