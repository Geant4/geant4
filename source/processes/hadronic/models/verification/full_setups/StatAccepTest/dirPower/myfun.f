C-------------------------------------------------------------
C Last update: 20-Jan-2006 
C
C The Fortran function is a splinline, defined as followed:
C  - a polynomial of 5 degrees, if x less than a given threshold;
C  - a decreasing exponential, for other values.
C
C We have two sets of the 9 parameters (6 for the polynomial + 
C 2 for the exponential + 1 for the threshold) which have 
C been fitted from the two following distributions:
C
C  a) distributionA : visible energy for pi+ 9 GeV Pb-Sci
C     0.7 mm production cut, 20000 events, QGSP_GN, G4 8.0.
C            mean = 45.735628 +/- 0.17413024 
C
C  b) distributionB : visible energy for pi+ 9 GeV Pb-Sci
C     0.7 mm production cut, 20000 events, QGSP_GN, G4 8.0
C     but with G4LElasticB.cc .
C            mean = 52.384921 +/- 0.18443447  
C
C These distributions have a similar shape, with the
C distributionB shifted on the right w.r.t. distributionA.
C
C By changing the parameter  alpha , between 0 and 1, we
C get a smooth family of similar functions, which goes from
C the distributionA, when alpha=0.0, to distributionB, when
C alpha=1.0 .
C
C This Fortran function is used by the following kumacs:
C    check.kumac , plot.kumac , smooth.kumac 
C 
C-------------------------------------------------------------
C
       FUNCTION MYFUN(X)
C
       REAL PA(6), CA, SA, XJA
       DATA PA(1)/15.291/
       DATA PA(2)/-15.217/
       DATA PA(3)/4.3102/
       DATA PA(4)/-0.14385/ 
       DATA PA(5)/0.17296E-02/
       DATA PA(6)/-0.71385E-05/
       DATA CA/8.9573/
       DATA SA/-0.48792E-01/ 
       DATA XJA/70.0/
C
       REAL PB(6), CB, SB, XJB
       DATA PB(1)/10.942/
       DATA PB(2)/-9.4891/
       DATA PB(3)/2.2055/
       DATA PB(4)/-0.51487E-01/ 
       DATA PB(5)/0.34915E-03/
       DATA PB(6)/-0.31766E-06/
       DATA CB/9.1410/
       DATA SB/-0.46462E-01/ 
       DATA XJB/75.0/
C
       REAL A(6), C, S, XJ, VAL, alpha
C
       alpha = 1.0   ! <--- ***LOOKHERE***
C
       DO I=1,6 
          A(I) = (1.0-alpha)*PA(I) + alpha*PB(I)
          C = (1.0-alpha)*CA + alpha*CB
          S = (1.0-alpha)*SA + alpha*SB
       ENDDO
C
       XJ = (1.0-alpha)*XJA + alpha*XJB
       IF (X.LT.2.0) THEN
          VAL = 0.0
       ELSE IF (X.LT.XJ) THEN
          VAL = A(1) + A(2)*X + A(3)*X*X + A(4)*X*X*X + A(5)*X*X*X*X
     &        + A(6)*X*X*X*X*X
       ELSE
          VAL = EXP(C+S*X)
       ENDIF
C
       MYFUN=VAL
C
       END
