C Copyright (C) 2010, Guy Barrand. All rights reserved.
C See the file tools.license for terms.

      PROGRAM HBOOKNT
      IMPLICIT NONE

      INTEGER IQUEST
      COMMON/QUEST/IQUEST(100)
      INTEGER NPAWC
      PARAMETER(NPAWC=1500000)
      INTEGER P
      COMMON /PAWC/ P(NPAWC)

      INTEGER IER,ICYCLE
C............................................................
      CALL HLIMIT(NPAWC)

      IQUEST(10) = 65000        ! maximum number of records
C     option N: new, Q: to change the Quest(10) value
      CALL HROPEN(20,'ntuple','out.hbook','NQ',1024,IER)
      IF(IER.NE.0)THEN
        PRINT *,'HROPEN failed.'
        STOP
      ENDIF


C      CALL RWBOOK
C      CALL RWFILL
      CALL CWBOOK
      CALL CWFILL

      ICYCLE = 0
      CALL HROUT(0,ICYCLE,' ')
      CALL HREND('ntuple')
      CLOSE(20)

      STOP
      END

C////////////////////////////////////////////////////////////
      SUBROUTINE RWBOOK
      IMPLICIT NONE
C
      INTEGER NVAR
      PARAMETER(NVAR=3)
      CHARACTER*8 CHTAGS(NVAR)
      DATA CHTAGS/'v1','v2','v3'/
C............................................................
      CALL HBOOKN(10,'row wise',NVAR,'//RWNT',5000,CHTAGS)
      RETURN
      END

C////////////////////////////////////////////////////////////
      SUBROUTINE RWFILL
      IMPLICIT NONE
C
      INTEGER NVAR
      PARAMETER(NVAR=3)
      REAL VAR(NVAR)

      INTEGER I
      REAL RGAUSS,RBW
C............................................................
      DO I=1,1000
        VAR(1) = I
        VAR(2) = RGAUSS(1.,2.)
        VAR(3) = RBW(0.,1.)
        CALL HFN(10,VAR)
      ENDDO
C
      RETURN
      END

C////////////////////////////////////////////////////////////
      SUBROUTINE CWBOOK
      IMPLICIT NONE
C
      REAL*8 VAR1,VAR2,VAR3
      COMMON/NTUC/VAR1,VAR2,VAR3
C............................................................
      CALL HBNT(11,'column wise',' ')
      CALL HBNAME(11,'index',VAR1,'index:R*8')
      CALL HBNAME(11,'rg'   ,VAR2,'rg:R*8')
      CALL HBNAME(11,'rbw'  ,VAR3,'rbw:R*8')
      RETURN
      END

C////////////////////////////////////////////////////////////
      SUBROUTINE CWFILL
      IMPLICIT NONE
C
      REAL*8 VAR1,VAR2,VAR3
      COMMON/NTUC/VAR1,VAR2,VAR3
C
      INTEGER I
      REAL RGAUSS,RBW
C............................................................
      DO I=1,1000
        VAR1 = I
        VAR2 = RGAUSS(1.,2.)
        VAR3 = RBW(0.,1.)
        CALL HFNT(11)
      ENDDO
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
