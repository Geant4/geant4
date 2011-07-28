C Copyright (C) 2010, Guy Barrand. All rights reserved.
C See the file tools.license for terms.

C  This program produces a out.hbook file with
C an histo 1D, an histo 2D and a profile 1D.
C  With CERN-PAW, you can open the file and plot
C the histos with something as :
C   PAW > h/file 101 out.hbook
C   PAW > h/list
C   PAW > h/plot 10
C   PAW > h/plot 20
C   PAW > h/plot 30

      PROGRAM HBOOKT
      IMPLICIT NONE

      INTEGER NPAWC
      PARAMETER (NPAWC = 1000000)
      INTEGER P
      COMMON /PAWC/ P(NPAWC)

      INTEGER I
      INTEGER ENTRIES
      INTEGER PRNT

      REAL RGAUSS,RBW
      INTEGER IER,ICYCLE
C............................................................
      CALL HLIMIT(NPAWC)

C      PRNT = 1
      PRNT = 0

      ENTRIES = 100000

C Create histos directory :
      CALL HCDIR('//PAWC',' ')
      CALL HMDIR('histos',' ')
      CALL HCDIR('histos',' ')

      CALL HBOOK1(10,'Gauss',100,-5.,5.,0.)
      DO I=1,ENTRIES
        CALL HFILL(10,RGAUSS(1.,2.),0.,1.4)
      ENDDO

      CALL HBOOK2(20,'Gauss_BW',100,-5.,5.,100,-2.,2.,0.)
      DO I=1,ENTRIES
        CALL HFILL(20,RGAUSS(1.,2.),RBW(0.,1.),0.8)
      ENDDO

      CALL HBPROF(30,'Profile',100,-5.,5.,-2.,2.,' ')
      DO I=1,ENTRIES
        CALL HFILL(30,RGAUSS(1.,2.),RBW(0.,1.),1.)
      ENDDO

      CALL HROPEN(20,'histos','out.hbook','NQ',1024,ier)
      IF(IER.NE.0) THEN
        PRINT *,'main : cannot open file'
        GOTO 999
      ENDIF

      ICYCLE = 0
      CALL HROUT(0,ICYCLE,'T')
      CALL HREND('histos')
      CLOSE(20)

999   STOP
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
