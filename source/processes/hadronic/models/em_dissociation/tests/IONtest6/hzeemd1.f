      PROGRAM HZEEMD1
      REAL AP, ZP, AT, ZT
      AP = 12.0
      ZP = 6.0
      AT = 12.0
      ZT = 6.0
      CALL PRINTEMD (AP, ZP, AT, ZT)
      
      AP = 12.0
      ZP = 6.0
      AT = 14.0
      ZT = 7.0
      CALL PRINTEMD (AP, ZP, AT, ZT)
      
      AP = 12.0
      ZP = 6.0
      AT = 27.0
      ZT = 13.0
      CALL PRINTEMD (AP, ZP, AT, ZT)
      
      AP = 12.0
      ZP = 6.0
      AT = 56.0
      ZT = 26.0
      CALL PRINTEMD (AP, ZP, AT, ZT)
      
      AP = 12.0
      ZP = 6.0
      AT = 181.0
      ZT = 73.0
      CALL PRINTEMD (AP, ZP, AT, ZT)
      
      AP = 12.0
      ZP = 6.0
      AT = 197.0
      ZT = 79.0
      CALL PRINTEMD (AP, ZP, AT, ZT)
      
      AP = 56.0
      ZP = 26.0
      AT = 12.0
      ZT = 6.0
      CALL PRINTEMD (AP, ZP, AT, ZT)
      
      AP = 56.0
      ZP = 26.0
      AT = 14.0
      ZT = 7.0
      CALL PRINTEMD (AP, ZP, AT, ZT)
      
      AP = 56.0
      ZP = 26.0
      AT = 27.0
      ZT = 13.0
      CALL PRINTEMD (AP, ZP, AT, ZT)
      
      AP = 56.0
      ZP = 26.0
      AT = 56.0
      ZT = 26.0
      CALL PRINTEMD (AP, ZP, AT, ZT)
      
      AP = 56.0
      ZP = 26.0
      AT = 181.0
      ZT = 73.0
      CALL PRINTEMD (AP, ZP, AT, ZT)
      
      AP = 56.0
      ZP = 26.0
      AT = 197.0
      ZT = 79.0
      CALL PRINTEMD (AP, ZP, AT, ZT)
      
      STOP
      END
C
C**************************************************************
C      
      SUBROUTINE PRINTEMD (AP, ZP, AT, ZT)
      REAL AP, ZP, AT, ZT, AF, ZF, QJ, E, MULT
      WRITE (6,'(//"PROJECTILE AP = ", F5.1, " ZP = ", F5.1)') AP, ZP
      WRITE (6,'(  "TARGET     AT = ", F5.1, " ZT = ", F5.1)') AT, ZT
      WRITE (6,'("         ENERGY  CROSS-SECTION")')
      WRITE (6,'("      [MeV/NUC]       [mbarns]")')
      MULT = 10.0**(1.0/10.0)
      DO 1010, I=0, 60, 1
        XSEC = 0.0
        E    = MULT**I
        CALL YIELDEM (AP, ZP, AT, ZT, AP-1.0,ZP, E, QJ)
        XSEC = XSEC + QJ
        CALL YIELDEM (AP, ZP, AT, ZT, AP-1.0,ZP-1.0, E, QJ)
        XSEC = XSEC + QJ
        CALL YIELDEM (AT, ZT, AP, ZP, AT-1.0,ZT, E, QJ)
        XSEC = XSEC + QJ
        CALL YIELDEM (AT, ZT, AP, ZP, AT-1.0,ZT-1.0, E, QJ)
        XSEC = XSEC + QJ
        WRITE (6,'(2F15.7)') E, XSEC
 1010 CONTINUE
      RETURN
      END
C
C**************************************************************
C
      SUBROUTINE YIELDEM(AP,ZP,AT,ZT,AF,ZF,TLAB,QJ)
C
C   PURPOSE
C   CALCULATES ELECTROMAGNETIC DISSOCIATION CROSS SECTIONS FOR ONE
C   NUCLEON REMOVAL
C
C   DESCRIPTION OF PARAMETERS
C   AP - MASS OF PROJECTILE
C   ZP - CHARGE OF PROJECTILE
C   AT - MASS OF TARGET
C   ZT - CHARGE OF TARGET
C   AF - MASS OF FRAGMENT
C   ZF - CHARGE OF FRAGMENT
C   TLAB - LAB ENERGY MEV/NUCLEON
C   QJ - ELECTROMAGNETIC CROSS SECTION FOR ONE NUCLEON REMOVAL
C
C   USAGE
C   CALL YIELDEM(AP,ZP,AT,ZT,AF,ZF,TLAB,QJ)
C
      REAL MNCSQ,II,INT,INTD,INTQ,K0,K1,NDIP,NQUAD,NP,NT,MSTAR,JAY,NU
C
      XDEE=.25
      IF(AF.NE.AP-1.) RETURN
      NT=AT-ZT
      NP=AP-ZP
      PI=3.141592653589793238D0
      FSC=0.00729735D0
      HBARC=197.32858D0
      MNCSQ=938.95D0
C
C   DIPOLE PARAMETERS
C
      JAY=36.8
      RZERO=1.18*AP**(1.D0/3.D0)
      QPRIM=17.0
      EPS=0.0768
      MSTAR=0.7*MNCSQ
      UU=3.0*JAY/(QPRIM*AP**(1.D0/3.D0))
      XFF=(1.0+EPS+3.0*UU)/(1.0+EPS+UU)
      EGDR=HBARC/SQRT(MSTAR*RZERO**2*(1.0+UU-XFF*EPS)/(8.0*JAY))
      FTRK=1.0
C
C   QUADRUPOLE PARAMETERS
C
      IF(AP.GT.100.0) FEWSR=0.9
      IF(AP.LE.100.0) FEWSR=0.6
      IF(AP.LE.40.0) FEWSR=0.3
      EGQR=63.0/AP**(1.0/3.0)
C
      IF(ZP.GE.14.0) GP=1.95*EXP(-0.075*ZP)
      IF(ZP.LT.14.0) GP=0.7
      IF(ZP.LE.8.0) GP=0.6
      IF(ZP.LT.6.0) GP=0.5
C
      GAMMA=1.0+TLAB/MNCSQ
      VEL=SQRT(1.0-1.0/GAMMA**2)
      HL=1.0/3.0
      HLL=-HL
      BMIN=1.34*(AP**HL+AT**HL-0.75*(AP**HLL+AT**HLL))
      DEEHILL=XDEE*BMIN
C
C   DEEHILL IS THE ONLY 'FUDGE' FACTOR IN THE CODE
C
      BMIN=BMIN+DEEHILL
      REDMAS=(AP*AT/(AP+AT))*MNCSQ
C
C   NOW APPLY BERTULANI LOW ENERGY CORRECTION TO BMIN
C
      BMIN=BMIN+PI*ZT*ZP*FSC*HBARC/(2.*GAMMA*REDMAS*VEL**2)
C
      SIGD=60.*NP*ZP/AP
      SIGQ=FEWSR*0.00022*ZP*AP**(2./3.)
      ECUTOF=HBARC*GAMMA*VEL/BMIN
      GD=EGDR/ECUTOF
      GQ=EGQR/ECUTOF
      CALL BESSEL(GD,K0,K1)
C
      NDIP=((2.0*ZT**2*FSC)/(EGDR*PI*VEL**2))*
     +(GD*K0*K1-0.5*VEL**2*GD**2*(K1**2-K0**2))
      CALL BESSEL(GQ,K0,K1)
      NQUAD=((2.0*ZT**2*FSC)/(EGQR*PI*VEL**4))*
     +(2.0*(1.0-VEL**2)*K1**2 + GQ*(2.0-VEL**2)**2*K0*K1
     +-0.5*VEL**4*GQ**2*(K1**2-K0**2))
      INTD=SIGD*NDIP
      INTQ=SIGQ*NQUAD*EGQR**2
      TOT=INTD+INTQ
      IF(ZP.EQ.ZF+1.) QJ=GP*TOT
      IF(ZP.EQ.ZF) QJ=(1.-GP)*TOT
C
      RETURN
  888 CONTINUE
      STOP
      END
C
C**************************************************************
C
      SUBROUTINE BESSEL(G,K0,K1)
C
C   PURPOSE
C   CALCULATES MODIFIED BESSEL FUNCTION OF SECOND KIND
C
C   DESCRIPTION OF PARAMETERS
C   G - INPUT ARGUMENT
C   K0 - OUTPUT K0(G)
C   K1 - OUTPUT K1(G)
C   A - ARRAY OF COEFFICIENTS OF APPROXIMATING POLYNOMIALS
C   B - ARRAY OF COEFFICIENTS OF APPROXIMATING POLYNOMIALS
C
C   USAGE
C   CALL BESSEL(G,K0,K1)
C
      DIMENSION A(30),B(27)
C
      REAL I0,I1,K0,K1
C
      DATA (A(I),I=1,30)/
     +3.5156229,3.0899424,1.2067492,.2659732,.0360768,.0045813,
     +.39894228,.01328592,.00225319,.00157565,.00916281,.02057706,
     +.02635537,.01647633,.00392377,.87890594,.51498869,.15084934,
     +.02658733,.00301532,.00032411,.39894228,.03988024,.00362018,
     +.00163801,.01031555,.02282967,.02895312,.01787654,.00420059/
      DATA (B(I),I=1,27)/
     +.57721566,.42278420,.23069756,.0348859,.00262698,.0001075,
     +.0000074,1.25331414,.07832358,.02189568,.01062446,.00587872,
     +.00251540,.00053208,.15443144,.67278579,.18156897,.01919402,
     +.00110404,.00004686,1.25331414,.23498619,.03655620,.01504268,
     +.00780353,.00325614,.00068245/
C
      T=G/3.75
      IF(G.LE.3.75) THEN
      I0=1.+A(1)*T**2+A(2)*T**4+A(3)*T**6+A(4)*T**8+A(5)*T**10
     ++A(6)*T**12
      I1=G*(.5+A(16)*T**2+A(17)*T**4+A(18)*T**6+A(19)*T**8
     ++A(20)*T**10+A(21)*T**12)
      ELSE
      I0=1./SQRT(G)*EXP(G)*(A(7)+A(8)/T+A(9)/T**2-A(10)/T**3
     ++A(11)/T**4-A(12)/T**5+A(13)/T**6-A(14)/T**7+A(15)/T**8)
      I1=1./SQRT(G)*EXP(G)*(A(22)-A(23)/T-A(24)/T**2+A(25)/T**3
     +-A(26)/T**4+A(27)/T**5-A(28)/T**6+A(29)/T**7-A(30)/T**8)
      END IF
      S=G/2.
      IF(G.LE.2.) THEN
      K0=-ALOG(S)*I0-B(1)+B(2)*S**2+B(3)*S**4+B(4)*S**6+B(5)*S**8
     ++B(6)*S**10+B(7)*S**12
      K1=ALOG(S)*I1+1./G*(1.+B(15)*S**2-B(16)*S**4-B(17)*S**6-
     +B(18)*S**8-B(19)*S**10-B(20)*S**12)
      ELSE
      K0=1./SQRT(G)*EXP(-G)*(B(8)-B(9)/S+B(10)/S**2-B(11)/S**3
     ++B(12)/S**4-B(13)/S**5+B(14)/S**6)
      K1=1./SQRT(G)*EXP(-G)*(B(21)+B(22)/S-B(23)/S**2+B(24)/S**3
     +-B(25)/S**4+B(26)/S**5-B(27)/S**6)
      END IF
      RETURN
      END
