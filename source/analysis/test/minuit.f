C Copyright (C) 2010, Guy Barrand. All rights reserved.
C See the file tools.license for terms.

      PROGRAM MINUIT
      IMPLICIT NONE

C////////////////////////////////////////////////////////////
C////// HBOOK ///////////////////////////////////////////////
C////////////////////////////////////////////////////////////
      INTEGER NPAWC
      PARAMETER (NPAWC = 1000000)
      INTEGER P
      COMMON /PAWC/ P(NPAWC)

      INTEGER I
      INTEGER ENTRIES

      REAL RGAUSS,RBW

C////////////////////////////////////////////////////////////
C////// MINUIT //////////////////////////////////////////////
C////////////////////////////////////////////////////////////
      INTEGER IERR
      REAL ARGS(1)
      EXTERNAL FCN

      REAL*8 ZERO
      REAL*8 HBEG,HSTP
      REAL*8 MBEG,MSTP
      REAL*8 WBEG,WSTP
C............................................................
      CALL HLIMIT(NPAWC)

      ENTRIES = 1000000

      CALL HBOOK1(10,'Gauss',100,-5.,5.,0.)
      DO I=1,ENTRIES
        CALL HFILL(10,RGAUSS(1.,2.),0.,1.4)
C        CALL HFILL(10,RGAUSS(0.,2.),0.,1.)
      ENDDO

      CALL MNINIT(5,6,7)

      ZERO = 0

      HBEG = 10000
      HSTP = 0.01
      CALL MNPARM(1,'H',HBEG,HSTP,ZERO,ZERO,IERR)
      IF(IERR.NE.0) THEN
        PRINT *,'MNPARM : error for H.'
        STOP
      ENDIF

      MBEG = 0
      MSTP = 0.01
      CALL MNPARM(2,'M',MBEG,MSTP,ZERO,ZERO,IERR)
      IF(IERR.NE.0) THEN
        PRINT *,'MNPARM : error for M.'
        STOP
      ENDIF

      WBEG = 3
      WSTP = 0.01
      CALL MNPARM(3,'W',WBEG,WSTP,ZERO,ZERO,IERR)
      IF(IERR.NE.0) THEN
        PRINT *,'MNPARM : error for W.'
        STOP
      ENDIF

      ARGS(1) = 1
      CALL MNEXCM(FCN,'SET PRI',ARGS,1,IERR,0)
      ARGS(1) = 0
      CALL MNEXCM(FCN,'MIGRAD',ARGS,0,IERR,0)

      STOP
      END


C////////////////////////////////////////////////////////////
      SUBROUTINE FCN(NPAR,GRAD,F,PARS,IFLAG,FUTIL) 
      IMPLICIT NONE
      INTEGER NPAR,IFLAG
      REAL*8 GRAD,F,PARS
      DIMENSION GRAD(*),PARS(*)
      EXTERNAL FUTIL

      REAL VBIN,XBIN

      REAL*8 CHI2,XBIND,VBIND,VAL,R
      INTEGER I,NBIN
      REAL*8 FGAUSS
      REAL HX
C............................................................
      NBIN = 100

      CHI2 = 0
      DO I=1,NBIN
        XBIN = -5 + (I-1) * 0.1 + 0.05
        VBIN = HX(10,XBIN)

        XBIND = XBIN
        VAL = FGAUSS(XBIND,PARS(1),PARS(2),PARS(3))
        VBIND = VBIN
        R = (VBIND-VAL)/0.1          
        CHI2 = CHI2 + R * R
      ENDDO

      F = CHI2
      RETURN
      END

C////////////////////////////////////////////////////////////
      REAL*8 FUNCTION FGAUSS(X,H,M,W)
      IMPLICIT NONE
      REAL*8 X,H,M,W,V
      REAL*8 EXP
C............................................................
      V = (X - M)/W
      FGAUSS = H * EXP(-0.5 * V * V)
      RETURN
      END      

C////////////////////////////////////////////////////////////
      REAL FUNCTION RGAUSS(MEAN,SIGMA)
      IMPLICIT NONE
      REAL MEAN,SIGMA

      REAL V1,V2,R,FAC
      REAL RDM
C............................................................
 10   CONTINUE
        V1 = 2.0 * RDM() - 1.0
        V2 = 2.0 * RDM() - 1.0
        R = V1 * V1 + V2 * V2
      IF(R.GT.1.0) GOTO 10
      FAC = SQRT(-2.0 * LOG(R)/R)
      RGAUSS = (V2 * FAC) * SIGMA + MEAN
      RETURN
      END      

C////////////////////////////////////////////////////////////
      REAL FUNCTION RBW(MEAN,GAMMA)
      IMPLICIT NONE
      REAL MEAN,GAMMA

      REAL R,D
      REAL RDM

      REAL M_PI_2
      DATA M_PI_2 /1.5707963/
C............................................................
      R = 2.0 * RDM() - 1.0
      D = 0.5 * GAMMA * TAN(R * M_PI_2)
      RBW = MEAN + D
      RETURN
      END      

C////////////////////////////////////////////////////////////
      REAL FUNCTION RDM()
      IMPLICIT NONE
C CERNLIB RNDM
      REAL RNDM
      REAL DUM
C The FORTRAN RAND not available if using f2c.
C      REAL RAND
C............................................................
C      RDM = RAND()
      RDM = RNDM(DUM)
      RETURN
      END      
