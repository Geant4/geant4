      PROGRAM HZECHK
      REAL AP, ZP, AT, ZT
      AP = 4.0
      ZP = 2.0
      AT = 4.0
      ZT = 2.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 4.0
      ZP = 2.0
      AT = 7.0
      ZT = 3.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 4.0
      ZP = 2.0
      AT = 12.011
      ZT = 6.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 4.0
      ZP = 2.0
      AT = 26.98154
      ZT = 13.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 4.0
      ZP = 2.0
      AT = 9.01218
      ZT = 4.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 4.0
      ZP = 2.0
      AT = 14.0067
      ZT = 7.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 4.0
      ZP = 2.0
      AT = 55.847
      ZT = 26.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 4.0
      ZP = 2.0
      AT = 180.9479
      ZT = 73.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 4.0
      ZP = 2.0
      AT = 196.9665
      ZT = 79.0
      CALL PRINTFP (AP, ZP, AT, ZT)
    
      AP = 10.0
      ZP = 5.0
      AT = 4.0
      ZT = 2.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 10.0
      ZP = 5.0
      AT = 7.0
      ZT = 3.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 10.0
      ZP = 5.0
      AT = 12.011
      ZT = 6.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 10.0
      ZP = 5.0
      AT = 26.98154
      ZT = 13.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 10.0
      ZP = 5.0
      AT = 9.01218
      ZT = 4.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 10.0
      ZP = 5.0
      AT = 14.0067
      ZT = 7.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 10.0
      ZP = 5.0
      AT = 55.847
      ZT = 26.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 10.0
      ZP = 5.0
      AT = 180.9479
      ZT = 73.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 10.0
      ZP = 5.0
      AT = 196.9665
      ZT = 79.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      
      AP = 12.0
      ZP = 6.0
      AT = 4.0
      ZT = 2.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 12.0
      ZP = 6.0
      AT = 7.0
      ZT = 3.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 12.0
      ZP = 6.0
      AT = 12.011
      ZT = 6.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 12.0
      ZP = 6.0
      AT = 26.98154
      ZT = 13.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 12.0
      ZP = 6.0
      AT = 9.01218
      ZT = 4.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 12.0
      ZP = 6.0
      AT = 14.0067
      ZT = 7.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 12.0
      ZP = 6.0
      AT = 55.847
      ZT = 26.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 12.0
      ZP = 6.0
      AT = 180.9479
      ZT = 73.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 12.0
      ZP = 6.0
      AT = 196.9665
      ZT = 79.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      
      AP = 56.0
      ZP = 26.0
      AT = 4.0
      ZT = 2.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 56.0
      ZP = 26.0
      AT = 7.0
      ZT = 3.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 56.0
      ZP = 26.0
      AT = 12.011
      ZT = 6.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 56.0
      ZP = 26.0
      AT = 26.98154
      ZT = 13.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 56.0
      ZP = 26.0
      AT = 9.01218
      ZT = 4.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 56.0
      ZP = 26.0
      AT = 14.0067
      ZT = 7.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 56.0
      ZP = 26.0
      AT = 55.847
      ZT = 26.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 56.0
      ZP = 26.0
      AT = 180.9479
      ZT = 73.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      AP = 56.0
      ZP = 26.0
      AT = 196.9665
      ZT = 79.0
      CALL PRINTFP (AP, ZP, AT, ZT)
      
      STOP
      END
C
C**************************************************************
C
      SUBROUTINE PRINTFP (AP, ZP, AT, ZT)      
      REAL AP, ZP, AT, ZT, E, RP, RT, RPT, B, F, P, EX
      E  = 100.0 * AP
      RP = 1.29*SQRT(RADIUS(AP)**2-0.84*0.84)
      RT = 1.29*SQRT(RADIUS(AT)**2-0.84*0.84)
      RPT = RP + RT
      PRINT *,""
      PRINT *, ' PROJECTILE AP = ', AP, ', ZP = ', ZP, ' RADIUS = ', RP
      PRINT *, ' TARGET     AT = ', AT, ', ZT = ', ZT, ' RADIUS = ', RT
      DO 1001, B = 0.0, RPT, RPT/10.0
        CALL BSEACH (E, AP, RP, RT, B, F, P)
        WRITE (6, '(E15.8, 5F15.7)') E, B, F, P
 1001 CONTINUE
      RETURN
      END
C
C**************************************************************
C
      SUBROUTINE BSEACH(E,AP,RP,RT,B,F,P)
C
C   PURPOSE
C   SUBROUTINE BSEACH USES A GEOMETRICAL APPROACH TO FIND ABR AND ABL
C   CALCULATES DELTAA AS A FUNCTION OF IMPACT PARAMETER
C   INCLUDES F.S.I.
C
C   DESCRIPTION OF PARAMETERS
C   E - LAB ENERGY MEV/NUCLEON
C   AP - MASS OF PROJECTILE
C   RT - RADIUS OF TARGET
C   RP - RADIUS OF PROJECTILE
C   AF - MASS OF FRAGMENT
C   B - OUTPUT; IMPACT PARAMETER
C   ABR - OUTPUT; ABRADED NUCLEONS
C   ABL - OUTPUT; ABRADED NUCLEONS
C
C   USAGE
C   CALL BSEACH(E,AP,RP,RT,AF,B,ABR,ABL)
C
      REAL RTT(5),ABRT(300,5,2),ABLT(300,5,2),BT(300,5,2)
C
c      DATA NUM,INUM/0,0/
C
C      FSI=1.
C      ES=E
C    1 CONTINUE
C      IF(NUM.EQ.0) GO TO 9006
C
C      DO 9005 INNN=1,NUM
C      IFX=INNN
C      IF(RT.EQ.RTT(INNN)) GO TO 9001
C 9005 CONTINUE
C
C 9006 INUM=INUM+1
C      IF(INUM.GT.5) INUM=1
C      NUM=NUM+1
C      IF(NUM.GT.5) NUM=5
C      GO TO 9000
C 9001 INDEX=1.01+FSI
C      IAF=AF+.51
C      B=BT(IAF,IFX,INDEX)
C      ABR=ABRT(IAF,IFX,INDEX)
C      ABL=ABLT(IAF,IFX,INDEX)
C      RETURN
C 9000 CONTINUE
      IF(AP.GT.300.) WRITE(6,6969)
      IF(AP.GT.300.) WRITE(6,6996)
C      OFSI=FSI
C      OAF=AF
C      RTT(INUM)=RT
C
C      DO 9002 IN=1,2
C      FSI=IN-1
C      UAF=AP-.5
C
C      DO 9003 AFF=.5,UAF,1.
C      AF=AFF
      ABLMIN=0.
      ABLMAX=0.
      ABRMIN=AP
      ABRMAX=0.
      BMAX=RT+RP
C      BMIN=0.
      XAVE = 16.6/(ES**0.26)
      UN=RP/BMAX
      UM=RT/RP
C      IIT=0
C   70 CONTINUE
C      IF(IIT.EQ.13)GO TO 2000
C      B=(BMAX+BMIN)/2.
      BTA=B/(RP+RT)
      IF(RT.LT.RP) GO TO 1000
      IF(B.LT.(RT-RP)) GO TO 10
      P=.125*TSQR(UM*UN)*(1./UM-2.)*((1.-BTA)/UN)**2
     +-.125*(.5*TSQR(UM*UN)*(1./UM-2.)+1.)*((1.-BTA)/UN)**3
      F=.75*TSQR(1.-UN)*((1.-BTA)/UN)**2-.125*(3.*TSQR(1.-UN)-1.)
     +*((1.-BTA)/UN)**3
      GO TO 20
   10 CONTINUE
      P=-1.
      F=1.
      GO TO 20
 1000 CONTINUE
      IF(B.LT.(RP-RT)) GO TO 1010
      P=.125*TSQR(UN*UM)*(1./UM-2.)*((1.-BTA)/UN)**2
     +-.125*(.5*TSQR(UN/UM)*(1./UM-2.)-(TSQR(1.-UM*UM)/UN-1.)
     +*TSQR((2.-UM)*UM)/UM**3)*((1.-BTA)/UN)**3
      F=.75*TSQR(1.-UN)*((1.-BTA)/UN)**2
     +-.125*(3.*TSQR(1.-UN)/UM-(1.-(1.-UM*UM)**1.5)*(1.-(1.-UM)
     +**2)**.5/UM**3)*((1.-BTA)/UN)**3
      GO TO 20
 1010 CONTINUE
      P=(TSQR(1.-UM*UM)/UN-1.)*TSQR(1.-(BTA/UN)**2)
      F=(1.-(1.-UM*UM)**1.5)*TSQR(1.-(BTA/UN)**2)
   20 ROT=ABS(B-RP)
      ROP=ABS(B-RT)
      IF(ROT.GE.RT) ROT=RT
      IF(ROP.GT.RP) ROP=RP
C
C   CLT IS LONGITUDINAL CHORD IN TARGET
C   CLP IS LONGITUDINAL CHORD IN PROJECTILE
C
      CLT=2.*TSQR(RT*RT-ROT*ROT)
      CLP=2.*TSQR(RP*RP-ROP*ROP)
      ATTEN=1.-.5*TEXP(-CLT/XAVE)-.5*TEXP(-CLP/XAVE)
C
C   ABL HERE IS DUE TO SURFACE DEFORMATION ONLY
C
      FAB=1.-F
      IF(FAB .LT. 1.E-12) FAB=0.
C   EB IS THE BINDING PER NUCLEON FOR THE ABLATION STAGE
      EB = 10.
      ABL=4.*3.1415*RP*RP*(1.+P-FAB**.6667)*.95/EB
      RO=ABS(B-RT)
C
C   FUDGE IS A SEMI-EMPRICAL CORRECTION TO DEFORMATION ENERGY
C
      FUDGE = 1 + 5*F
      IF ( RT.LT.RP.AND.B.LT.(RP-RT))FUDGE = FUDGE + 25*F*F
C      IF(RO.LT.RP)GO TO 333
C      RO=RP
C  333 AEX=1.3*CLP
C      BP=(RP*RP+B*B-RT*RT)/(2.*B+1.E-12)
C      IF (BP.LT.0.) BP=0.
C      IF(BP.GE.RP)BP=RP-1.E-12
C      CT=2.*TSQR(RP*RP-BP*BP)
C      IF(CT.LT.1.5)CT=1.5
C
C   AEX IS THE F.S.I. ENERGY CORRECTION
C
C      AEX =AEX*(1.+(CT-1.5)/3.)
C
C   USE NEXT LINE IF DON'T WANT ANY FSI
C   ABL=ABL*FUDGE+AEX*FSI*0.0
C
C      ABL=ABL*FUDGE+AEX*FSI
      EX = ABL * EB * FUDGE
      RETURN
C
C   AFP IS THE FINAL FRAGMENT MASS
C
C      AFP=AP-(F*AP+ABL)*ATTEN
C      IIT=IIT+1
C      IF(AF.GE.AFP) GO TO 21
C      BMAX=B
C      ABRMAX=AP*F*ATTEN
C      ABLMAX=ABL*ATTEN
C      GO TO 199
C   21 CONTINUE
C      BMIN=B
C      ABLMIN=ABL*ATTEN
C      ABRMIN=AP*F*ATTEN
C  199 CONTINUE
C      IF(ABS(AF-AFP).LT..0001) GO TO 2000
C      GO TO 70
C 2000 CONTINUE
C      ABL=(ABLMIN+ABLMAX)/2.
C      B=(BMAX+BMIN)/2.
C      ABR=(ABRMAX+ABRMIN)/2.
C      IAF=AF+.51
C      ABLT(IAF,INUM,IN)=ABL
C      ABRT(IAF,INUM,IN)=ABR
C      BT(IAF,INUM,IN)=B
C 9003 CONTINUE
C
C 9002 CONTINUE
C
C      FSI=OFSI
C      AF=OAF
C      GO TO 1
C      ENTRY BSEEK(E,AP,RP,RT,AF,B,ABR,ABL)
C      FSI=0.
C      ES=E
C      GO TO 1
 6996 FORMAT('YOUR VALUE OF AP IS TOO LARGE.')
 6969 FORMAT('THE BT,ABRT,ABLT ARRAYS ARE DIMENSION TO 300.')
      END
C
C**************************************************************
C
      FUNCTION TSQR(Y)
C   PURPOSE
C   TO ELIMINATE OVER/UNDER FLOW OF CPU IF SQRT IS USED
C
C   USAGE
C   RESULT = TSQR(Y)
C
      X=Y
      IF(X .LT. 1.E-37) X=1.E-37
      IF(X .GT. 1.E+37) X=1.E+37
      TSQR=SQRT(X)
      RETURN
      END
C
C**************************************************************
C
      FUNCTION TEXP(X)
C
C   PURPOSE
C   TO ELIMINATE OVER/UNDER FLOW OF CPU IF EXP IS USED
C
C   RESULT = TEXP(X)
      IF(X .LT.-80.) X=-80.
      IF(X .GT. 80.) X=80.
      TEXP =EXP(X)
      END
C
C**************************************************************
C
      FUNCTION RADIUS (A)
C
C   PURPOSE
C   GIVES RADIUS OF A NUCLEUS
C
C   DESCRIPTION OF PARAMETERS
C   A - MASS NUMBER OF A NUCLEUS
C
C   USAGE
C   RESULT = RADIUS (A)
C
      DIMENSION NA(23),RMS(23)
C
      DATA NA/1,2,3,4,6,7,9,10,11,12,13,14,15,16,17,18,19,20,22,
     +23,24,25,26/
      DATA RMS/0.85,2.095,1.976,1.671,2.57,2.41,2.519,2.45,2.42,
     +2.471,2.440,2.58,2.611,2.730,2.662,2.727,2.900,3.040,2.969,2.94,
     +3.075,3.11,3.06/
      FACT = SQRT (5./3.)
      IA = A + 0.4
      RADIUS = FACT * ( 0.84* A**(1./3.) + 0.55 )
C
      DO 1 I =1,23
      IF ( IA .EQ. NA(I)) GO TO 2
      GO TO 1
    2 RADIUS = FACT*RMS(I)
    1 CONTINUE
C
      RETURN
      END
