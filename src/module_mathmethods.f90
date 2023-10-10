!==============================================================================!
! MODULE MathMethods                                                           !
!                                                                              !
! This module contains the variables and routines related to mathematical fun- !
! ctions and numerical methods.                                                !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!
MODULE MathMethods
use Constants, only: r64, zero,one,two,half,third,pi,&
                     igfv,igfvbc
use Globals, only: gfv

implicit none

contains

subroutine GaussLaguerre(N,X,W)
!-------------------------------------------------------------------------------------!
!Purpose : Compute the zeros of Laguerre polynomial Ln(x)in the interval [0,inf],     !
!          and the corresponding weighting coefficients for Gauss-Laguerr integration.!
!                                                                                     !                                                                                                  
!Input   : n    --- Order of the Laguerre polynomial                                  !
!          X(n) --- Zeros of the Laguerre polynomial                                  !
!          W(n) --- Corresponding weighting coefficients                              !
!                                                                                     ! 
!Integral: WL  =  W * exp(X)                                                          !
!          \int_0^\infty  f(z) exp(-z) dz  =   \sum_i f(X(i)) W(i)                    !
!          \int_0^\infty  f(z) dz          =   \sum_i f(X(i)) WL(i)                    !
!-------------------------------------------------------------------------------------!
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    IMPLICIT INTEGER (I-N)
    DIMENSION X(N),W(N)
    HN=1.0D0/N
    DO 35 NR=1,N
        IF (NR.EQ.1) Z=HN
        IF (NR.GT.1) Z=X(NR-1)+HN*NR**1.27
        IT=0
10      IT=IT+1
        Z0=Z
        P=1.0D0
        DO 15 I=1,NR-1
            P=P*(Z-X(I))
15      CONTINUE
        F0=1.0D0
        F1=1.0D0-Z
        DO 20 K=2,N
            PF=((2.0D0*K-1.0D0-Z)*F1-(K-1.0D0)*F0)/K
            PD=K/Z*(PF-F1)
            F0=F1
            F1=PF
20        CONTINUE
        FD=PF/P
        Q=0.0D0
        DO 30 I=1,NR-1
            WP=1.0D0
            DO 25 J=1,NR-1
                IF (J.EQ.I) GO TO 25
                WP=WP*(Z-X(J))
25          CONTINUE
            Q=Q+WP
30      CONTINUE
        GD=(PD-Q*FD)/P
        Z=Z-FD/GD
        IF (IT.LE.40.AND.DABS((Z-Z0)/Z).GT.1.0D-15) GO TO 10
        X(NR)=Z
        W(NR)=1.0D0/(Z*PD*PD)           
35  CONTINUE
    RETURN
end subroutine GaussLaguerre

subroutine GaussHermite(N,X,W)
!-------------------------------------------------------------------------------------!
!Purpose : Compute the zeros of Hermite polynomial Ln(x) in the interval [-inf,inf],  ! 
!          and the corresponding weighting coefficients for Gauss-Hermite integration.!
!                                                                                     !
!Input   : n    --- Order of the Hermite polynomial                                   !  
!          X(n) --- Zeros of the Hermite polynomial                                   !
!          W(n) --- Corresponding weighting coefficients                              !
!                                                                                     !
!Integral: WH  =  W * exp(X**2)                                                       !
!          \int_0^\infty  f(z) exp(-z**2) dz  =   \sum_i f(X(i)) W(i)                 !
!          \int_0^\infty  f(z) dz             =   \sum_i f(X(i)) WH(i)                 !
!-------------------------------------------------------------------------------------!
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    IMPLICIT INTEGER (I-N)
    DIMENSION X(N),W(N)
    HN=1.0D0/N
    ZL=-1.1611D0+1.46D0*N**0.5
    DO 40 NR=1,N/2
        IF (NR.EQ.1) Z=ZL
        IF (NR.NE.1) Z=Z-HN*(N/2+1-NR)
        IT=0
10      IT=IT+1
        Z0=Z
        F0=1.0D0
        F1=2.0D0*Z
        DO 15 K=2,N
            HF=2.0D0*Z*F1-2.0D0*(K-1.0D0)*F0
            HD=2.0D0*K*F1
            F0=F1
            F1=HF
15      CONTINUE
        P=1.0D0
        DO 20 I=1,NR-1
            P=P*(Z-X(I))
20      CONTINUE
        FD=HF/P
        Q=0.0D0
        DO 30 I=1,NR-1
            WP=1.0D0
            DO 25 J=1,NR-1
                IF (J.EQ.I) GO TO 25
                WP=WP*(Z-X(J))
25          CONTINUE
            Q=Q+WP
30      CONTINUE
        GD=(HD-Q*FD)/P
        Z=Z-FD/GD
        IF (IT.LE.40.AND.DABS((Z-Z0)/Z).GT.1.0D-15) GO TO 10
        X(NR)=Z
        X(N+1-NR)=-Z
        R=1.0D0
        DO 35 K=1,N
            R=2.0D0*R*K
35      CONTINUE
        W(NR)=3.544907701811D0*R/(HD*HD)
        W(N+1-NR)=W(NR)
40  CONTINUE
    IF (N.NE.2*INT(N/2)) THEN
        R1=1.0D0
        R2=1.0D0
        DO 45 J=1,N
            R1=2.0D0*R1*J
            IF (J.GE.(N+1)/2) R2=R2*J
45      CONTINUE
        W(N/2+1)=0.88622692545276D0*R1/(R2*R2)
        X(N/2+1)=0.0D0
    ENDIF
    RETURN
end  subroutine GaussHermite

subroutine GaussLegendre(N,X,W)
!-------------------------------------------------------------------------------------!
!Purpose : Compute the zeros of Legendre polynomial Pn(x) in the interval [-1,1], and ! 
!          the corresponding weighting coefficients for Gauss-Legendre integration.   !
!                                                                                     !
!Input   : n    --- Order of the Legendre polynomial                                  !  
!          X(n) --- Zeros of the Legendre polynomial                                  !
!          W(n) --- Corresponding weighting coefficients                              !
!                                                                                     !
!Integral:                                                                            !
!          \int_0^{+1}     f(z) dz  =   \sum_i  f(X(i)) * W(i)                        !
!          \int_{-1}^{+1}  f(z) dz  =   \sum_i ( f(X(i)) + f(-X(i) ) * W(i)           !
!-------------------------------------------------------------------------------------!
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    IMPLICIT INTEGER (I-N)
    DIMENSION X(N),W(N)
    N0=(N+1)/2
    DO 45 NR=1,N0
        Z=DCOS(3.1415926D0*(NR-0.25D0)/N)
10      Z0=Z
        P=1.0D0
        DO 15 I=1,NR-1
            P=P*(Z-X(I))
15      CONTINUE
        F0=1.0D0
        IF (NR.EQ.N0.AND.N.NE.2*INT(N/2)) Z=0.0D0
        F1=Z
        DO 20 K=2,N
            PF=(2.0D0-1.0D0/K)*Z*F1-(1.0D0-1.0D0/K)*F0
            PD=K*(F1-Z*PF)/(1.0D0-Z*Z)
            F0=F1
            F1=PF
20      CONTINUE
        IF (Z.EQ.0.0) GO TO 40
        FD=PF/P
        Q=0.0D0
        DO 35 I=1,NR
            WP=1.0D0
            DO 30 J=1,NR
                IF (J.NE.I) WP=WP*(Z-X(J))
30          CONTINUE
            Q=Q+WP
35      CONTINUE
        GD=(PD-Q*FD)/P
        Z=Z-FD/GD
        IF (DABS(Z-Z0).GT.DABS(Z)*1.0D-15) GO TO 10
40      X(NR)=Z
        X(N+1-NR)=-Z
        W(NR)=2.0D0/((1.0D0-Z*Z)*PD*PD)
        W(N+1-NR)=W(NR)
45  CONTINUE
    RETURN
end subroutine GaussLegendre

subroutine math_gfv
!-------------------------------------------------------------------------------------!
!Purpose : Calculates sign, sqrt, factorials, etc. of integers and half int.          !
!                                                                                     !
!Input   :                                  !  
!                                                                                     !
!Variable Meaning:                                                                    !
!          iv(n)  =  (-1)**n                                                          !
!          sq(n)  =  sqrt(n)                                                          ! 
!          sqi(n) =  1/sqrt(n)                                                        !
!          sqh(n) =  sqrt(n+1/2)                                                      !
!          shi(n) =  1/sqrt(n+1/2)                                                    !
!          fak(n) =  n!                                                               !
!          ibc(m,n) = m!/(n!(m-n)!)                                                   !
!          fad(n) =  (2*n+1)!!                                                        !
!          fdi(n) =  1/(2*n+1)!!                                                      !
!          fi(n)  =  1/n!                                                             !
!          wf(n)  =  sqrt(n!)                                                         !
!          wfi(n) =  1/sqrt(n!)                                                       !
!          wfd(n) =  sqrt((2*n+1)!!)                                                  !
!          gm2(n) =  gamma(n+1/2)                                                     !
!          gmi(n) =  1/gamma(n+1/2)                                                   !
!          wg(n)  =  sqrt(gamma(n+1/2))                                               !
!          wgi(n) =  1/sqrt(gamma(n+1/2))                                             !
!-------------------------------------------------------------------------------------!
    integer :: i,m,n
    gfv%iv(0)  = +1
    gfv%sq(0)  = zero
    gfv%sqi(0) = 1.d30
    gfv%sqh(0) = sqrt(half)
    gfv%shi(0) = 1/gfv%sqh(0)
    gfv%fak(0) = one
    gfv%fad(0) = one
    gfv%fi(0)  = one
    gfv%fdi(0) = one
    gfv%wf(0)  = one
    gfv%wfi(0) = one
    gfv%wfd(0) =  one
    !gm2(0) = Gamma(1/2) = sqrt(pi)
    gfv%gm2(0) =  sqrt(pi)
    gfv%gmi(0) =  1/gfv%gm2(0)
    gfv%wg(0)  =  sqrt(gfv%gm2(0))
    gfv%wgi(0) =  1/gfv%wg(0)
    do i = 1,igfv
            gfv%iv(i)  = -gfv%iv(i-1)
            gfv%sq(i)  = dsqrt(dfloat(i))
            gfv%sqi(i) = one/gfv%sq(i)
            gfv%sqh(i) = sqrt(i+half)
            gfv%shi(i) = one/gfv%sqh(i)
            gfv%fak(i) = i*gfv%fak(i-1)
            gfv%fad(i) = (2*i+1)*gfv%fad(i-1)
            gfv%fi(i)  = one/gfv%fak(i)
            gfv%fdi(i) = one/gfv%fad(i)
            gfv%wf(i)  = gfv%sq(i)*gfv%wf(i-1)
            gfv%wfi(i) = one/gfv%wf(i)
            gfv%wfd(i) = sqrt(gfv%fad(i))
            gfv%gm2(i) = (i-half)*gfv%gm2(i-1)
            gfv%gmi(i) = one/gfv%gm2(i)
            gfv%wg(i)  = gfv%sqh(i-1)*gfv%wg(i-1)
            gfv%wgi(i) = one/gfv%wg(i)
    enddo

    gfv%ibc(0,0)= one
    do m=1,igfvbc
            do n=0,m
            gfv%ibc(m,n)=gfv%fak(m)/(gfv%fak(n)*gfv%fak(m-n)) 
            enddo   
    enddo  
    return
end subroutine math_gfv

! recursive function factorial(n) result(facto)
! !------------------------------------------------------------------------------!
! ! function factorial                                                           !
! !                                                                              !
! ! Computes the factorial: n! = n * (n-1) * ... * 1                             !
! !------------------------------------------------------------------------------!
! integer, intent(in) :: n
! real(r64) :: facto 
! if ( n <= 0 ) then 
!   facto = one  
! else
!   facto = n * factorial(n-1)
! endif
! end function



END MODULE MathMethods