C------- name of the file ----------------------------------------------
C                          DTULAP.FOR
C______________________________________________________________________
C
C        originally: JTDTU.FOR program
C        connection between DTU and JT ( hard scattering )
*
C        first parameters are taken from DTU and set to the own
C        parameter-commons of JT, then JT will be initialized
C
C______________________________________________________________________
*     revision 3.92:  adjust COMMONS,
*     caraful: DTU90 was based on a older version of DTULAP
*  ********************************************************************
      SUBROUTINE JTDTU(IOPT)
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( MAXPRO = 8 , MLINE = 1000 , MSCAHD = 250 )
      COMMON /HAPARA/ ECM,PTINI(4),Q0SQR,ALASQR,BQCD,NPD,NF,NHA,NHB
      COMMON /HAENVI/ NINDEP
      COMMON /HAOUTL/ NOUTL,NOUTER,NOUTCO
      COMMON /HAPADI/ NPDM
      COMMON /HAQQAP/ AQQAL,AQQPD,NQQAL,NQQPD
      COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
*
      CHARACTER*80 TITLE
      CHARACTER*8 PROJTY,TARGTY
C     COMMON /USER/TITLE,PROJTY,TARGTY,CMENER,ISTRUF
C    &            ,ISINGD,IDUBLD,SDFRAC,PTLAR
      COMMON /USER1/TITLE,PROJTY,TARGTY
      COMMON /USER2/CMENER,SDFRAC,PTLAR,ISTRUF,ISINGD,IDUBLD
*
      COMMON/COLLIS/S,IJPROJ,IJTAR,PTTHR,PTTHR2,IOPHRD,IJPRLU,IJTALU
C repl.      COMMON /COLLIS/ECMDTU,S,IJPROJ,IJTAR,PTTHR,IOPHRD,IJPRLU,IJTALU,PTTHR2
      COMMON /STRUFU/ISTRUM,ISTRUT
      COMMON /PTLARG/XSMAX
      COMMON /HAXSUM/XSHMX
C
C===read default values
C
      CALL HASTRT
C
C===read parameters from DTU
C
C  max. # of flavors
      NF = 4
C  partondistributions
      NPD  = ISTRUF
      NPDM = ISTRUM
C  correct scale for CTEQ PDFs
      IF((ISTRUF.GE.16).OR.(ISTRUF.LE.20)) THEN
        AQQAL = 1.D0
        AQQPD = 1.D0
      ENDIF
C  hadron a
      NHA  = IJPROJ
      IF ( IJPROJ.EQ.2 ) NHA =-1
C  hadron b
      NHB  = IJTAR
      IF ( IJTAR .EQ.2 ) NHB =-1
C  output level
      NOUTL = IOUTPA
C  cms-energy ( GeV )
C repl    ECM=ECMDTU
      ECM = CMENER
C  pt-cut ( GeV )
      PTINI(1) = PTTHR
      PTINI(2) = PTTHR2
      PTINI(3) = 0.0
      PTINI(4) = 0.0
C  maximum sum of hard x
      XSHMX    = XSMAX
C  is program called from DTU ( NINDEP=0 ) or independent ( NINDEP=1 )
      NINDEP   = IOPT
C
C===initialize JT
C
      CALL HISINI
      IF ( IOPT.EQ.0 ) CALL HARINI
      RETURN
      END
C******************************************************************************
      SUBROUTINE SELHRD(MHARD,IJPVAL,IJTVAL,PTTHRE)
C
C   select the initial parton x-fractions and flavors and the final flavors
C           for an event with mhard hard or semihard scatterings
C
C   IJPVAL,IJTVAL =0 valence quarks of projectile or target not involved
C                                                     in hard scattering
C   IJPVAL,IJTVAL =1 valence quarks of projectile or target  involved
C                                                     in hard scattering
C
C the results are in COMMON /ABRHRD/
C               XH1(I),XH2(I):     x-values of initial partons
C               IJHI1(I),IJHI2(I): flavor of initial parton
C                                  0            gluon
C                                  1,2          valence u,d quarks
C                                  11,12,13,14  sea udsc-quarks
C                                  negative     anti s or v quarks
C               IJHF1(I),IJHF2(I): flavor of final state partons
C               PHARD1(I,J),PHARD2(I,J): final part. momentum and energy
C                                J=1   PX
C                                 =2   PY
C                                 =3   PZ
C                                 =4   ENERGY  (massless partons)
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( MAXPRO = 8 , MLINE = 1000 , MSCAHD = 250 )
      CHARACTER*80 TITLE
      CHARACTER*8 PROJTY,TARGTY
C     COMMON /USER/TITLE,PROJTY,TARGTY,CMENER,ISTRUF
C    &            ,ISINGD,IDUBLD,SDFRAC,PTLAR
      COMMON /USER1/TITLE,PROJTY,TARGTY
      COMMON /USER2/CMENER,SDFRAC,PTLAR,ISTRUF,ISINGD,IDUBLD
      COMMON/COLLIS/S,IJPROJ,IJTAR,PTTHR,PTTHR2,IOPHRD,IJPRLU,IJTALU
      COMMON /ABRHRD/XH1(MSCAHD),XH2(MSCAHD),IJHI1(MSCAHD),
     *IJHI2(MSCAHD),IJHF1(MSCAHD),IJHF2(MSCAHD),PHARD1(MSCAHD,4),
     *PHARD2(MSCAHD,4)
      COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
      COMMON /HAPARA/ ECM,PTINI(4),Q0SQR,ALASQR,BQCD,NPD,NF,NHA,NHB
      COMMON /HAOUTL/ NOUTL,NOUTER,NOUTCO
      COMMON /HAEVTR/ LINE,LIN,LREC1(MLINE),LREC2(MLINE),PREC(0:3,MLINE)
      COMMON /HARSLT/ LSCAHD,LSC1HD,
     &                ETAHD(MSCAHD,2) ,PTHD(MSCAHD),
     &                XHD(MSCAHD,2)   ,VHD(MSCAHD) ,X0HD(MSCAHD,2),
     &                NINHD(MSCAHD,2) ,NOUTHD(MSCAHD,2),
     &                N0INHD(MSCAHD,2),NBRAHD(MSCAHD,2),NPROHD(MSCAHD)
      DATA X1SU/0./ , X2SU/0./
C
      IJPVAL =0
      IJTVAL =0
      IF (IOUTPA.GE.3) WRITE(6,221)
     *                  MHARD,IJPVAL,IJTVAL
  221 FORMAT (' SELHRD  ',3I10)
C  call of hard scattering routines
      CALL HAREVT(MHARD,PTTHR2)
C  select information from event-record
C   number of hard scatterings reached
      MHARD = LSCAHD
C   initial partons
      DO 10 N=1,LSCAHD
C    X-values
        XH1(N) = XHD(N,1)
        XH2(N) = XHD(N,2)
        X1SU = X1SU + XH1(N)
        X2SU = X2SU + XH2(N)
        IF( IOUTPA.GT. 6 )WRITE(6,*)N,X1SU,X2SU,XH1(N),XH2(N)
C    flavors
        III      = NINHD(N,1)
        IIIA     = ABS(III)
        IF ( IIIA.GT. 0 .AND. IIIA.LT.10 ) III    = SIGN(IIIA+10,III)
        IF ( IIIA.GE.10                  ) III    = SIGN(IIIA-10,III)
        IF ( IIIA.GE.10                  ) IJPVAL = 1
        IJHI1(N) = III
        III      = NINHD(N,2)
        IIIA     = ABS(III)
        IF ( IIIA.GT. 0 .AND. IIIA.LT.10 ) III    = SIGN(IIIA+10,III)
        IF ( IIIA.GE.10                  ) III    = SIGN(IIIA-10,III)
        IF ( IIIA.GE.10                  ) IJTVAL = 1
        IJHI2(N) = III
10      CONTINUE
C   final partons
      DO 30 N=1,LSCAHD
        I3     = 4*N-1
        I4     = 4*N
C    flavors
        IJHF1(N) = NOUTHD(N,1)
        IJHF2(N) = NOUTHD(N,2)
C    four momentum
        DO 20 J=1,3
          PHARD1(N,J) = PREC(J,I3)
20        PHARD2(N,J) = PREC(J,I4)
        PHARD1(N,4)   = PREC(0,I3)
        PHARD2(N,4)   = PREC(0,I4)
30    CONTINUE
C
C   output ( optional )
C
      IF (IOUTPA.GE.3)WRITE (6,101)
  101 FORMAT(' SELHRD OUTPUT FOR INITIAL STATE SCATTERED PARTONS')
      DO 102 I=1,LSCAHD
        IF (IOUTPA.GE.3)
     *    WRITE (6,103)I,IJPVAL,IJTVAL,IJHI1(I),IJHI2(I),XH1(I),XH2(I)
  103   FORMAT (' I,IJPVAL,IJTVAL,IJHI1,IJHI2,XH1,XH2= ',5I5,2F12.6)
  102 CONTINUE
      IF (IOUTPA.GE.3)WRITE (6,301)
  301 FORMAT(' SELHRD OUTPUT FOR FINAL STATE SCATTERED PARTONS')
      DO 302 I=1,LSCAHD
        IF (IOUTPA.GE.3)
     *    WRITE (6,303)I,IJHF1(I),IJHF2(I),(PHARD1(I,III),III=1,4)
        IF (IOUTPA.GE.3)
     *    WRITE (6,303)I,IJHF1(I),IJHF2(I),(PHARD2(I,III),III=1,4)
  303   FORMAT (' I,IJHI1,IJHI2,PHARD1 OR PHARD2 ',3I5,4F16.6)
  302 CONTINUE
      RETURN
      END
*
C______________________________________________________________________
C
C        originally:  JTWORK.FOR
C        procedures to simulate a single event
C
C
C     ( this procedures work independent to procedures in other blocks
C       if initialization was done )
C
C______________________________________________________________________
*
*  ********************************************************************
      SUBROUTINE HAREVT(MHARD,PT1IN)
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( MAXPRO = 8 , MLINE = 1000 , MSCAHD = 250 )
      COMMON /HAPARA/ ECM,PTINI(4),Q0SQR,ALASQR,BQCD,NPD,NF,NHA,NHB
      COMMON /HAENVI/ NINDEP
      COMMON /HAEVNT/ PT1,PT2,NHARD,NTRY,IHARD,ITRY,IREJEV

      PT1    = MAX(PT1IN,PTINI(1))
      PT2    = PTINI(1)
      NHARD  = MHARD
      IHARD  = 0
C     NTRY   = 5
      NTRY   = 200
C      NTRY   = 1000
      ITRY   = 0
      IREJEV = 0
      CALL HAMULT
      CALL HAOUTP
      IF ( NINDEP.EQ.1 ) CALL HISFIL
      RETURN
      END
C_______________________________________________________________________
C===========================================================================
C
C       THE FOLLOWING 5 SUBROUTINES ARE REWRITTEN BY
C        BY I.KAWRAKOW  IN ORDER TO BE ABLE 
C         TO PRODUCE GREAT NUMBER OF HARD POMERONS (>100) IN A 
C          SHORT TIME
C                                          VERSION IK.1 - 01.93
C--------------------------------------------------------------------     
      SUBROUTINE HAMULT
C--------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( MAXPRO = 8 , MLINE = 1000 , MSCAHD = 250 )
      PARAMETER ( TINY= 1.D-30, ONE=1.D0, ZSMALL=1.D-3 )
      COMMON /HAPARA/ ECM,PTINI(4),Q0SQR,ALASQR,BQCD,NPD,NF,NHA,NHB
      COMMON /HAPDCO/ NPDCOR
      COMMON /HAOUTL/ NOUTL,NOUTER,NOUTCO
      COMMON /HAEVNT/ PT1,PT2,NHARD,NTRY,IHARD,ITRY,IREJEV
      COMMON /HASCA / PTWANT,A,ALN,Z1MAX,Z1DIF,Z2MAX,Z2DIF,
     &                PT,ETAC,ETAD,X1,X2,V,U,W,W1,AXX,WEIGHT,MSPR,IREJSC
      COMMON /HAXIK / XREST,YREST,ZMAX,AXXMAX,WEMAX
      COMMON /HAEVTR/ LINE,LIN,LREC1(MLINE),LREC2(MLINE),PREC(0:3,MLINE)
      COMMON /HAXSUM/XSHMX
      INTEGER MXSECT
      COMMON /HAXSEC/ XSECTA(2,-1:MAXPRO,4),XSECT(5,-1:MAXPRO),
     &                MXSECT(0:2,-1:MAXPRO)
      COMMON /HARSLT/ LSCAHD,LSC1HD,
     &                ETAHD(MSCAHD,2) ,PTHD(MSCAHD),
     &                XHD(MSCAHD,2)   ,VHD(MSCAHD) ,X0HD(MSCAHD,2),
     &                NINHD(MSCAHD,2) ,NOUTHD(MSCAHD,2),
     &                N0INHD(MSCAHD,2),NBRAHD(MSCAHD,2),NPROHD(MSCAHD)
      ITYPE(L)      = MOD(LREC1(L),100)-50

      LINE   = 0
      LSCAHD = 0
      AA     = (2.*PT2/ECM)**2
      SA     = SQRT(AA)      
C
C  loop until event is accepted or too many attempts ( more then NTRY )
C
5     ITRY   = 0
20    ITRY   = ITRY+1
      IF(ITRY.GT.NTRY) GOTO 301
      LINE   = 0
      XREST    = XSHMX-NHARD*SA
      YREST    = XSHMX-NHARD*SA
      IF(XREST*YREST.LT.AA) THEN
        WRITE(6,*) ' ****************** HAMULT ****************** '
        WRITE(6,*) ' IT IS NOT POSSIBLE TO PRODUCE ',NHARD,' POMERONS '
        NHARD=0
        RETURN
C       STOP
      ENDIF	
      ZMAX=XREST*YREST
      AXXMAX=AA/ZMAX
      WEMAX =SQRT(1-AXXMAX)
      X1S    = 0.0
      X2S    = 0.0
      IHARD  = 0
      PTWANT = PT1
10    CONTINUE
        A        = (2.*PTWANT/ECM)**2
        SA       = SQRT(A)
        I = 5
50      I = I-1
        IF ( PT1.LT.PTINI(I) .AND. I.GT.1 ) GOTO 50
        DO 60 M=-1,MAXPRO
            XSECT(1,M) = XSECTA(1,M,I)
            XSECT(2,M) = XSECTA(2,M,I)
60     CONTINUE
       CALL HARSCA
       X1S   = X1S+X1
       X2S   = X2S+X2
       xrest=xrest-x1+sa
       yrest=yrest-x2+sa
       zmax=xrest*yrest
       IHARD=IHARD+1
       LSCAHD         = IHARD
       XHD(IHARD,1)   = X1
       XHD(IHARD,2)   = X2
       VHD(IHARD)     = V
       ETAHD(IHARD,1) = ETAC
       ETAHD(IHARD,2) = ETAD
       PTHD(IHARD)    = PT
       NPROHD(IHARD)  = MSPR
       if(zmax/a-one.lt.ZSMALL) THEN
         CALL XCHECK(X1S,X2S,LINMAX)
	 GOTO 10 
       ENDIF	 
       AXXMAX=A/ZMAX
       WEMAX=SQRT(1.-AXXMAX)
       PTWANT   = PT2
       IF(IHARD.LT.NHARD) GOTO 10
C-------------------------------------------------- NOW THE REQUIRED NUMBER
C                                              OF POMERTONS IS CREATED       
      IF ( NPDCOR.EQ.1     .AND.
     &     IHARD .GT.1     .AND.
     &     (1.-X1S)*(1.-X2S).LT.RNDM(AI)*(1.-AA*IHARD)**2 ) GOTO 5
 301  CONTINUE 
C
C  end of loop
C
C  check choice of valence quarks
      DO 120 K=1,2
        IVAL  = 0
        DO 110 I=1,IHARD
          IND = 4*(I-1)+K
          IT  = ITYPE(IND)
          IF ( ABS(IT).GT.10 .AND. IVAL.EQ.0 ) THEN
            IVAL       = 1
          ELSEIF ( ABS(IT).GT.10 .AND. IVAL.EQ.1 ) THEN
            IT         = SIGN(ABS(IT)-10,IT)
            LREC1(IND) = (LREC1(IND)/100)*100+50+IT
          ENDIF
C  fill COMMON HARSLT
          NINHD(I,K)   = IT
          NOUTHD(I,K)  = ITYPE(IND+2)
110     CONTINUE
120   CONTINUE
C
C  information if HAMULT is not able to produce the required # of scatt.
C
      IF ( IHARD.NE.NHARD .AND. NOUTER.EQ.1 ) THEN
        WRITE(6,1010) NHARD,IHARD
1010    FORMAT(' ###### HAMULT : CANNOT PRODUCE',I3,' HARD SCATT.',
     &         '; ONLY',I3,' ARE PRODUCED !!!')
      ENDIF
      RETURN
      END

C______________________________________________________________________
      SUBROUTINE RECCHK( LINMAX,X,IOPT )
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( MAXPRO = 8 , MLINE = 1000 , MSCAHD = 250 )
      COMMON /HAPARA/ ECM,PTINI(4),Q0SQR,ALASQR,BQCD,NPD,NF,NHA,NHB
      COMMON /HAPDCO/ NPDCOR
      COMMON /HAOUTL/ NOUTL,NOUTER,NOUTCO
      COMMON /HAEVNT/ PT1,PT2,NHARD,NTRY,IHARD,ITRY,IREJEV
      COMMON /HASCA / PTWANT,A,ALN,Z1MAX,Z1DIF,Z2MAX,Z2DIF,
     &                PT,ETAC,ETAD,X1,X2,V,U,W,W1,AXX,WEIGHT,MSPR,IREJSC
      COMMON /HAEVTR/ LINE,LIN,LREC1(MLINE),LREC2(MLINE),PREC(0:3,MLINE)
      INTEGER MXSECT
      COMMON /HAXSEC/ XSECTA(2,-1:MAXPRO,4),XSECT(5,-1:MAXPRO),
     &                MXSECT(0:2,-1:MAXPRO)
      COMMON /HARSLT/ LSCAHD,LSC1HD,
     &                ETAHD(MSCAHD,2) ,PTHD(MSCAHD),
     &                XHD(MSCAHD,2)   ,VHD(MSCAHD) ,X0HD(MSCAHD,2),
     &                NINHD(MSCAHD,2) ,NOUTHD(MSCAHD,2),
     &                N0INHD(MSCAHD,2),NBRAHD(MSCAHD,2),NPROHD(MSCAHD)
C
      IF( IOPT.EQ.0 ) THEN
        LSTART = LINMAX + 1
        DO 1 L = LSTART,LINE
             LP = L - 4
          PREC(1,LP) = PREC(1,L)
          PREC(2,LP) = PREC(2,L)
          PREC(3,LP) = PREC(3,L)
          PREC(0,LP) = PREC(0,L)
          LREC1( LP) = LREC1( L)
          LREC2( LP) = LREC2( L)
   1    CONTINUE
        LINE = LINE - 4
        RETURN
      ELSEIF( IOPT.EQ.1 ) THEN
        QTEST = 0.5*ECM*X
        DO 2 L=1,LINE
          PTEST = PREC(0,L)
          IF( PTEST.EQ.QTEST ) THEN
            LINMAX = L
            RETURN
          ENDIF
   2    CONTINUE
        WRITE(6,*)'  RECCHK: NO NEW LINMAX FOUND - LINMAX=',LINMAX
        RETURN
      ENDIF
      WRITE(6,*)'  RECCHK: IOPT OUT OF RANGE - 0 OR 1 - IOPT=',IOPT
      RETURN
      END
C______________________________________________________________________
      SUBROUTINE XCHECK( X1S, X2S, LINMAX )
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( MAXPRO = 8 , MLINE = 1000 , MSCAHD = 250 )
      COMMON /HAPARA/ ECM,PTINI(4),Q0SQR,ALASQR,BQCD,NPD,NF,NHA,NHB
      COMMON /HAPDCO/ NPDCOR
      COMMON /HAOUTL/ NOUTL,NOUTER,NOUTCO
      COMMON /HAEVNT/ PT1,PT2,NHARD,NTRY,IHARD,ITRY,IREJEV
      COMMON /HASCA / PTWANT,A,ALN,Z1MAX,Z1DIF,Z2MAX,Z2DIF,
     &                PT,ETAC,ETAD,X1,X2,V,U,W,W1,AXX,WEIGHT,MSPR,IREJSC
      COMMON /HAXIK / XREST,YREST,ZMAX,AXXMAX,WEMAX
      COMMON /HAEVTR/ LINE,LIN,LREC1(MLINE),LREC2(MLINE),PREC(0:3,MLINE)
      COMMON /HAXSUM/XSHMX
      INTEGER MXSECT
      COMMON /HAXSEC/ XSECTA(2,-1:MAXPRO,4),XSECT(5,-1:MAXPRO),
     &                MXSECT(0:2,-1:MAXPRO)
      COMMON /HARSLT/ LSCAHD,LSC1HD,
     &                ETAHD(MSCAHD,2) ,PTHD(MSCAHD),
     &                XHD(MSCAHD,2)   ,VHD(MSCAHD) ,X0HD(MSCAHD,2),
     &                NINHD(MSCAHD,2) ,NOUTHD(MSCAHD,2),
     &                N0INHD(MSCAHD,2),NBRAHD(MSCAHD,2),NPROHD(MSCAHD)
      PARAMETER (ONE=1D0, ZSMALL=1D-3)
C
50    CONTINUE
      IF(IHARD.LT.1) THEN
        WRITE(6,*) '  ERROR IN XCHECK : IHARD < 1 ',IHARD
	STOP
      ENDIF	
C-------------------------------------- FIND PROCESS WITH THE MAX. X
      IMAX=0
      XMAX=0.
      DO 10 I=1,IHARD
        IF(XHD(I,1).GT.XMAX) THEN
	  IMAX=I
	  XMAX=XHD(I,1)
	ENDIF
	IF(XHD(I,2).GT.XMAX) THEN
	  IMAX=I
	  XMAX=XHD(I,2)
	ENDIF
10    CONTINUE
C--------------------------------------- REJECT THIS PROCESS
      X1S=X1S-XHD(IMAX,1)
      X2S=X2S-XHD(IMAX,2)
      XREST=XREST+XHD(IMAX,1)-SQRT(A)
      YREST=YREST+XHD(IMAX,2)-SQRT(A)
      ZMAX=XREST*YREST	    
      AXXMAX=A/ZMAX
      WEMAX=SQRT(1.-AXXMAX)
      MH=0
      DO 20 I=1,IHARD
        IF(I.NE.IMAX) THEN
	  MH=MH+1
          XHD(MH,1)   = XHD(I,1)
          XHD(MH,2)   = XHD(I,2)
          VHD(MH)     = VHD(I)
          ETAHD(MH,1) = ETAHD(I,1)
          ETAHD(MH,2) = ETAHD(I,2)
          PTHD(MH)    = PTHD(I)
          NPROHD(MH)  = NPROHD(I)
	ENDIF
20    CONTINUE	  
      CALL RECCHK( 4*IMAX,XHD1,0)
      IHARD=IHARD-1
      LSCAHD=IHARD
      IF(ZMAX/A-ONE.LT.ZSMALL) GOTO 50
      RETURN
      END
C_______________________________________________________
      SUBROUTINE HAX1X2
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      COMMON /HAEVNT/ PT1,PT2,NHARD,NTRY,IHARD,ITRY,IREJEV
      COMMON /HASCA / PTWANT,A,ALN,Z1MAX,Z1DIF,Z2MAX,Z2DIF,
     &                PT,ETAC,ETAD,X1,X2,V,U,W,W1,AXX,WEIGHT,MSPR,IREJSC
C      COMMON /HAXIK / XREST,YREST,ZMAX
      COMMON /HAXIK / XREST,YREST,ZMAX,AXXMAX,WEMAX
      PARAMETER ( TINY= 1.D-30, ONE=1.D0 ,TINY6=1.D-06)

      SA=SQRT(A)
12    continue
c--------------------------------- sample z=x*y
      z=a*Exp(rndm(1.)*Log(zmax/a))
      xm=xrest
      ym=yrest
      if(xm.lt.yrest) then
        xm=yrest
        ym=xrest
      endif   
      ww=Log(xm**2/z)/Log(xm**2/a)
      if(rndm(1.1).gt.ww) goto 12
c--------------------------------- sample u=x+y
      umin=Sqrt(4.*z)
      umax=xm+z/xm
      cc=umax**2-4.*z
      if(cc.lt.0.) cc=0.
13    continue
      c=Exp(rndm(2.)*Log((umax+Sqrt(cc))/umin))
      uu=umin*(c**2+1.)/2./c
      if(uu.gt.2.*ym.and.uu.lt.ym+z/ym) goto 13
c------------------------------------- x,y from u,z
      c=uu**2-4.*z
      if(c.lt.0.)  c=0.
      c=sqrt(c)
      xtemp=(uu+c)/2.
      ytemp=(uu-c)/2.
      if(xrest.ge.yrest) then
         x=xtemp
         y=ytemp
         if(xrest.eq.yrest) then
           if(rndm(3.).gt.0.5) then
             x=ytemp
             y=xtemp
           endif
         endif
      else
         x=ytemp
         y=xtemp
      endif
      X1=X
      X2=Y
      AXX  = A/(X1*X2)
      W    = SQRT(MAX(TINY,ONE-AXX))
      W1   = AXX/(1.+W)
      RETURN
      END
C______________________________________________________________________
      SUBROUTINE HARKIN
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( MAXPRO = 8 , MLINE = 1000 , MSCAHD = 250 )
      PARAMETER ( TINY= 1.D-30, ONE=1.D0 ,TINY6=1.D-06)
      COMMON /HAPARA/ ECM,PTINI(4),Q0SQR,ALASQR,BQCD,NPD,NF,NHA,NHB
      COMMON /HAEVNT/ PT1,PT2,NHARD,NTRY,IHARD,ITRY,IREJEV
      COMMON /HASCA / PTWANT,A,ALN,Z1MAX,Z1DIF,Z2MAX,Z2DIF,
     &                PT,ETAC,ETAD,X1,X2,V,U,W,W1,AXX,WEIGHT,MSPR,IREJSC
      DIMENSION RM(-1:MAXPRO)
      COMMON /HAXIK / XREST,YREST,ZMAX,AXXMAX,WEMAX
      DATA RM / 3.31, 0.0,
     &          3.80, 0.65, 2.00, 0.65, 0.89, 0.45, 0.445, 0.89 /
      M    = MSPR
      IF     ( M.EQ.1 ) THEN
10      CALL HAX1X2
        V  =-0.5*W1/(W1+RNDM(AI)*W)
        U  =-1.-V
        R  = (1.+W)*2.25*(V*V*(3.-U*V-V/(U*U))-U)
        RMAX=RM(1)*WEMAX*(1.+WEMAX)
        WIK=R*W/RMAX
        IF(WIK.GT.1.D0) WRITE(6,*) ' HARKIN : WIK > 1 : ',M,R
C        IF ( R*W.LT.RM(1)*RNDM(AI) ) GOTO 10
        IF(WIK.LT.RNDM(AI)) GOTO 10
        IF ( RNDM(AJ).LE.0.5D0 ) V = U
      ELSEIF ( M.EQ.2 .OR. M.EQ.4 ) THEN
20      CALL HAX1X2
        WL = LOG(W1)
        V  =-EXP(-0.6931472+RNDM(AI)*WL)
        U  =-1.-V
        R  = (U*U+V*V)*((16./27.)/U-(4./3.)*V)*(WL/W)*AXX
        IF ( R*W.LT.RM(M)*RNDM(AI) ) GOTO 20
        IF ( RNDM(AJ).LE.0.5D0 ) V = U
      ELSEIF ( M.EQ.3 ) THEN
30      CALL HAX1X2
        V  =-0.5*W1/(W1+RNDM(AI)*W)
        U  =-1.-V
        R  = (1.+W)*(1.+U*U)*(1.-(4./9.)*V*V/U)
        RMAX=RM(3)*WEMAX*(1.+WEMAX)
        WIK=R*W/RMAX
        IF(WIK.GT.1.D0) WRITE(6,*) ' HARKIN : WIK > 1 : ',M,R
C        IF ( R*W.LT.RM(3)*RNDM(AI) ) GOTO 30
        IF(WIK.LT.RNDM(AI)) GOTO 30
      ELSEIF ( M.EQ.5 ) THEN
50      CALL HAX1X2
        V  =-0.5*AXX/(W1+2.*RNDM(AI)*W)
        U  =-1.-V
        R  = (4./9.)*(1.+U*U+V*V*(U*U+V*V))-(8./27.)*U*U*V
        RMAX=RM(5)*WEMAX
        WIK=R*W/RMAX
        IF(WIK.GT.1.D0) WRITE(6,*) ' HARKIN : WIK > 1 : ',M,R
C        IF ( R*W.LT.RM(5)*RNDM(AI) ) GOTO 50
        IF(WIK.LT.RNDM(AI)) GOTO 50
      ELSEIF ( M.EQ.6 ) THEN
60      CALL HAX1X2
        V  =-0.5*(1.+W)+RNDM(AI)*W
        U  =-1.-V
        R  = (4./9.)*(U*U+V*V)*AXX
        IF ( R*W.LT.RM(6)*RNDM(AI) ) GOTO 60
      ELSEIF ( M.EQ.7 ) THEN
70      CALL HAX1X2
        V  =-0.5*W1/(W1+RNDM(AI)*W)
        U  =-1.-V
        R  = (1.+W)*((2./9.)*(1.+U*U+(1.+V*V)*V*V/(U*U))-(4./27.)*V/U)
        RMAX=RM(7)*WEMAX*(1.+WEMAX)
        WIK=R*W/RMAX
        IF(WIK.GT.1.D0) WRITE(6,*) ' HARKIN : WIK > 1 : ',M,R
C        IF ( R*W.LT.RM(7)*RNDM(AI) ) GOTO 70
        IF(WIK.LT.RNDM(AI)) GOTO 70
        IF ( RNDM(AJ).LE.0.5D0 ) V = U
      ELSEIF ( M.EQ.8 ) THEN
80      CALL HAX1X2
        V  =-0.5*AXX/(W1+2.*RNDM(AI)*W)
        U  =-1.-V
        R  = (4./9.)*(1.+U*U)
        RMAX=RM(8)*WEMAX
        WIK=R*W/RMAX
        IF(WIK.GT.1.D0) WRITE(6,*) ' HARKIN : WIK > 1 : ',M,R
C        IF ( R*W.LT.RM(8)*RNDM(AI) ) GOTO 80
        IF(WIK.LT.RNDM(AI)) GOTO 80
      ELSEIF ( M.EQ.-1 ) THEN
90      CALL HAX1X2
        WL = LOG(W1)
        V  =-EXP(-0.6931472+RNDM(AI)*WL)
        U  =-1.-V
        R  = (1.+V*V)*(V/(U*U)-(4./9.))*(WL/W)*AXX
        IF ( R*W.LT.RM(-1)*RNDM(AI) ) GOTO 90
      ENDIF
C      PARAMETER ( TINY= 1.D-30, ONE=1.D0 ,TINY6 =1.D-06)
      V    = MAX(MIN(      V,-TINY6 ),-1.+TINY6 )
      U    = MAX(MIN(-1.E0-V,-TINY6 ),-1.+TINY6 )
      PT   = SQRT(U*V*X1*X2)*ECM
      ETAC = 0.5*LOG((U*X1)/(V*X2))
      ETAD = 0.5*LOG((V*X1)/(U*X2))
      RETURN
      END
C-------------------------------------------- END OF CHANGES BY IK 01.93
C===========================================================================      
C_______________________________________________________________________

      SUBROUTINE HACHEK(IOPT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      COMMON /HACUTS/ PTL,PTU,ETACL,ETACU,ETADL,ETADU
      COMMON /HASCA / PTWANT,A,ALN,Z1MAX,Z1DIF,Z2MAX,Z2DIF,
     &                PT,ETAC,ETAD,X1,X2,V,U,W,W1,AXX,WEIGHT,MSPR,IREJSC

      IOPT = 1
      IF (  PT  .LT.PTL    .OR.  PT  .GT.PTU
     & .OR. ETAC.LT.ETACL  .OR.  ETAC.GT.ETACU
     & .OR. ETAD.LT.ETADL  .OR.  ETAD.GT.ETADU ) IOPT = 0
      RETURN
      END
C______________________________________________________________________
      SUBROUTINE HAFDIS(PDS,PDA,PDB,FDISTR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( MAXPRO = 8 , MLINE = 1000 , MSCAHD = 250 )
      PARAMETER ( TINY= 1.D-30, ONE=1.D0 ,TINY6=1.D-06)
      COMMON /HAPARA/ ECM,PTINI(4),Q0SQR,ALASQR,BQCD,NPD,NF,NHA,NHB
      COMMON /HAQQAP/ AQQAL,AQQPD,NQQAL,NQQPD
      COMMON /HAEVNT/ PT1,PT2,NHARD,NTRY,IHARD,ITRY,IREJEV
      COMMON /HASCA / PTWANT,A,ALN,Z1MAX,Z1DIF,Z2MAX,Z2DIF,
     &                PT,ETAC,ETAD,X1,X2,V,U,W,W1,AXX,WEIGHT,MSPR,IREJSC
      DIMENSION PDA(-6:6),PDB(-6:6)
      INTEGER MXSECT
      COMMON /HAXSEC/ XSECTA(2,-1:MAXPRO,4),XSECT(5,-1:MAXPRO),
     &                MXSECT(0:2,-1:MAXPRO)

      FDISTR = 0.0
C  set hard scale  QQ  for alpha and partondistr.
      IF     ( NQQAL.EQ.1 ) THEN
        QQAL = AQQAL*PT*PT
      ELSEIF ( NQQAL.EQ.2 ) THEN
        QQAL = AQQAL*X1*X2*ECM*ECM
      ELSEIF ( NQQAL.EQ.3 ) THEN
        QQAL = AQQAL*X1*X2*ECM*ECM*(U*V)**(1./3.)
      ELSEIF ( NQQAL.EQ.4 ) THEN
        QQAL = AQQAL*X1*X2*ECM*ECM*U*V/(1.+V*V+U*U)
      ENDIF
      IF     ( NQQPD.EQ.1 ) THEN
        QQPD = AQQPD*PT*PT
      ELSEIF ( NQQPD.EQ.2 ) THEN
        QQPD = AQQPD*X1*X2*ECM*ECM
      ELSEIF ( NQQPD.EQ.3 ) THEN
        QQPD = AQQPD*X1*X2*ECM*ECM*(U*V)**(1./3.)
      ELSEIF ( NQQPD.EQ.4 ) THEN
        QQPD = AQQPD*X1*X2*ECM*ECM*U*V/(1.+V*V+U*U)
      ENDIF
      ALPHA = BQCD/LOG(MAX(QQAL/ALASQR,1.1*ONE))
      F  = XSECT(1,MSPR)*ALPHA**2
C calculate partondistributions
      CALL JTPDIS(X1,QQPD,NHA,MSPR,PDA)
      CALL JTPDIS(X2,QQPD,NHB,MSPR,PDB)
C calculate full distribution FDISTR
      IF ( MSPR.EQ.1  .OR.  MSPR.EQ.4 ) THEN
        PDS   = PDA(0)*PDB(0)
      ELSE
        S2    = 0.0
        S3    = 0.0
        S4    = 0.0
        S5    = 0.0
        DO 10 I=1,NF
          S2  = S2+PDA(I)*PDB(-I)+PDA(-I)*PDB( I)
          S3  = S3+PDA(I)*PDB( I)+PDA(-I)*PDB(-I)
          S4  = S4+PDA(I)+PDA(-I)
          S5  = S5+PDB(I)+PDB(-I)
10        CONTINUE
        IF     ( MSPR.EQ.2  .OR.  MSPR.EQ.5  .OR.  MSPR.EQ.6 ) THEN
          PDS = S2
        ELSEIF ( MSPR.EQ.3  .OR.  MSPR.EQ.-1 ) THEN
          PDS = PDA(0)*S5+PDB(0)*S4
        ELSEIF ( MSPR.EQ.7 ) THEN
          PDS = S3
        ELSEIF ( MSPR.EQ.8 ) THEN
          PDS = S4*S5-(S2+S3)
        ENDIF
      ENDIF
      FDISTR  = F*PDS
      RETURN
      END
C______________________________________________________________________
      SUBROUTINE HARSCA
C  HARSCA determines the type of hard subprocess, the partons taking
C  part in subprocess and the kinematic variables
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( MAXPRO = 8 , MLINE = 1000 , MSCAHD = 250 )
      COMMON /HAPARA/ ECM,PTINI(4),Q0SQR,ALASQR,BQCD,NPD,NF,NHA,NHB
      COMMON /HAOUTL/ NOUTL,NOUTER,NOUTCO
      COMMON /HAEVNT/ PT1,PT2,NHARD,NTRY,IHARD,ITRY,IREJEV
      COMMON /HASCA / PTWANT,A,ALN,Z1MAX,Z1DIF,Z2MAX,Z2DIF,
     &                PT,ETAC,ETAD,X1,X2,V,U,W,W1,AXX,WEIGHT,MSPR,IREJSC
      DIMENSION PDA(-6:6),PDB(-6:6)
      COMMON /HAEVTR/ LINE,LIN,LREC1(MLINE),LREC2(MLINE),PREC(0:3,MLINE)
      INTEGER MXSECT
      COMMON /HAXSEC/ XSECTA(2,-1:MAXPRO,4),XSECT(5,-1:MAXPRO),
     &                MXSECT(0:2,-1:MAXPRO)

      MXSECT(0,0) = 0
      XSECT(2,0)  = 0.0
      DO 15 M=-1,MAXPRO
        IF ( MXSECT(0,M).EQ.1 ) XSECT(2,0) = XSECT(2,0)+XSECT(2,M)
15    CONTINUE
C
C -------------------------------------------I
C  begin of iteration loop                   I
C                                            I
      IREJSC = 0
10    CONTINUE
      IREJSC = IREJSC+1
      IREJEV = IREJEV+1
C  find subprocess
      B      = RNDM(AI)*XSECT(2,0)
      MSPR   =-2
      SUM    = 0.0
20    MSPR   = MSPR+1
      IF ( MXSECT(0,MSPR).EQ.1 ) SUM = SUM+XSECT(2,MSPR)
      IF ( SUM.LT.B  .AND. MSPR.LT.MAXPRO ) GOTO 20
C  find kin. variables X1,X2 and V
      CALL HARKIN
C  check kin. cuts eventually given by user
      CALL HACHEK(IOPT)
      IF ( IOPT.EQ.0 ) GOTO 10
C  calculate remaining distribution
      CALL HAFDIS(PDS,PDA,PDB,F)
C  actualize counter for cross-section calculation
      IF( F .LE. 1.D-15 ) F=0.
      XSECT (3,MSPR) = XSECT (3,MSPR)+F
      XSECT (4,MSPR) = XSECT (4,MSPR)+F*F
      MXSECT(1,MSPR) = MXSECT(1,MSPR)+1
C  check F against FMAX
C
      WEIGHT = F/XSECT(2,MSPR)
      IF ( WEIGHT.LT.RNDM(AI) ) GOTO 10
C-------------------------------------------------------------------
C      IF(WEIGHT.GT.1.D0) WRITE(6,1234)F,XSECT(2,MSPR),WEIGHT
C1234  FORMAT(' HARSCA: MONTE-CARLO WEIGHT FUNCTION H/HMAX GT 1 !',/
C    * '  H = SUM OVER A,B FOR PROCESS M OF:',/
C    * '  E(M)*ALPHAS**2*XA*FA(XA,Q**2)*XB*FB(XB,Q**2)',/
C    * '  F(=H),XSECT(2,MSPR)(=HMAX), WEIGHT = H/HMAX',3E12.5)
C-------------------------------------------------------------------
C                                            I
C  end of iteration loop                     I
C -------------------------------------------I
C
C  the event is accepted now
C
C  actualize counter for accepted events
      MXSECT(2,MSPR) = MXSECT(2,MSPR)+1
      IF ( MSPR.EQ.-1 ) MSPR = 3
C  find initial partons
      SUM    = 0.0
      SCHECK = RNDM(AI)*PDS
      IF     ( MSPR.EQ.1  .OR.  MSPR.EQ.4 ) THEN
        IA = 0
        IB = 0
      ELSEIF ( MSPR.EQ.2  .OR.  MSPR.EQ.5  .OR.  MSPR.EQ.6 ) THEN
        DO 610 IA=-NF,NF
          IF ( IA.EQ.0 ) GOTO 610
          SUM  = SUM+PDA(IA)*PDB(-IA)
          IF ( SUM.GE.SCHECK ) GOTO 620
610       CONTINUE
620     IB =-IA
      ELSEIF ( MSPR.EQ.3 ) THEN
        IB     = 0
        DO 630 IA=-NF,NF
          IF ( IA.EQ.0 ) GOTO 630
          SUM  = SUM+PDA(0)*PDB(IA)
          IF ( SUM.GE.SCHECK ) GOTO 640
          SUM  = SUM+PDA(IA)*PDB(0)
          IF ( SUM.GE.SCHECK ) GOTO 650
630       CONTINUE
640     IB     = IA
        IA     = 0
650     CONTINUE
      ELSEIF ( MSPR.EQ.7 ) THEN
        DO 660 IA=-NF,NF
          IF ( IA.EQ.0 ) GOTO 660
          SUM  = SUM+PDA(IA)*PDB(IA)
          IF ( SUM.GE.SCHECK ) GOTO 670
660       CONTINUE
670     IB     = IA
      ELSEIF ( MSPR.EQ.8 ) THEN
        DO 690 IA=-NF,NF
          IF ( IA.EQ.0 ) GOTO 690
          DO 680 IB=-NF,NF
            IF ( ABS(IB).EQ.ABS(IA)  .OR.  IB.EQ.0 ) GOTO 680
            SUM = SUM+PDA(IA)*PDB(IB)
            IF ( SUM.GE.SCHECK ) GOTO 700
680         CONTINUE
690       CONTINUE
700     CONTINUE
      ENDIF
C  find final partons
      IC = IA
      ID = IB
      IF     ( MSPR.EQ.2 ) THEN
        IC = 0
        ID = 0
      ELSEIF ( MSPR.EQ.4 ) THEN
        IC = INT(FLOAT(NF+NF)*RNDM(AI))+1
        IF ( IC.GT.NF ) IC = NF-IC
        ID =-IC
      ELSEIF ( MSPR.EQ.6 ) THEN
        IC = INT(FLOAT(NF+NF-2)*RNDM(AI))+1
        IF ( IC.GT.NF-1 ) IC = NF-1-IC
        IF ( ABS(IC).EQ.ABS(IA) ) IC = SIGN(NF,IC)
        ID =-IC
      ENDIF
C
30    A1     = RNDM(AI)
      A2     = RNDM(AI)
      IF ( ((A1*A1)+(A2*A2)).GT.1.0D0 ) GOTO 30
      COSPHI = ((A1*A1)-(A2*A2))/((A1*A1)+(A2*A2))
      SINPHI = SIGN(((A1*A2)+(A1*A2))/((A1*A1)+(A2*A2)),RNDM(AI)-0.5)
C
      IF ( RNDM(AI)*PDA(IA).GT.PDA(-IA) ) IA = SIGN(ABS(IA)+10,IA)
      IF ( RNDM(AJ)*PDB(IB).GT.PDB(-IB) ) IB = SIGN(ABS(IB)+10,IB)
C  fill event record
      LINE          = LINE+1
      PREC(1,LINE)  = 0.0
      PREC(2,LINE)  = 0.0
      PREC(3,LINE)  = 0.5*ECM*X1
      PREC(0,LINE)  = PREC(3,LINE)
      LREC1(LINE)   = IA+50+100*MSPR
      LREC2(LINE)   = 01000
      LINE          = LINE+1
      PREC(1,LINE)  = 0.0
      PREC(2,LINE)  = 0.0
      PREC(3,LINE)  =-0.5*ECM*X2
      PREC(0,LINE)  =-PREC(3,LINE)
      LREC1(LINE)   = IB+50
      LREC2(LINE)   = 01000
      LINE          = LINE+1
      PREC(1,LINE)  = PT*COSPHI
      PREC(2,LINE)  = PT*SINPHI
      PREC(3,LINE)  =-0.5*ECM*(U*X1-V*X2)
      PREC(0,LINE)  =-0.5*ECM*(U*X1+V*X2)
      LREC1(LINE)   = IC+50
      LREC2(LINE)   = 11000
      LINE          = LINE+1
      PREC(1,LINE)  =-PT*COSPHI
      PREC(2,LINE)  =-PT*SINPHI
      PREC(3,LINE)  =-0.5*ECM*(V*X1-U*X2)
      PREC(0,LINE)  =-0.5*ECM*(V*X1+U*X2)
      LREC1(LINE)   = ID+50
      LREC2(LINE)   = 11000
      RETURN
      END
C_______________________________________________________________________
      SUBROUTINE HAOUTP
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( MAXPRO = 8 , MLINE = 1000 , MSCAHD = 250 )
      COMMON /HAOUTL/ NOUTL,NOUTER,NOUTCO
      COMMON /HAEVNT/ PT1,PT2,NHARD,NTRY,IHARD,ITRY,IREJEV
      COMMON /HAEVTR/ LINE,LIN,LREC1(MLINE),LREC2(MLINE),PREC(0:3,MLINE)
      COMMON /HARSLT/ LSCAHD,LSC1HD,
     &                ETAHD(MSCAHD,2) ,PTHD(MSCAHD),
     &                XHD(MSCAHD,2)   ,VHD(MSCAHD) ,X0HD(MSCAHD,2),
     &                NINHD(MSCAHD,2) ,NOUTHD(MSCAHD,2),
     &                N0INHD(MSCAHD,2),NBRAHD(MSCAHD,2),NPROHD(MSCAHD)

C  output of data for hard scattering
      IF ( NOUTL.GE.4 ) THEN
      WRITE(6,1010) NHARD,IHARD,IREJEV
1010  FORMAT(' ===HARD EVENT===    NHARD,NTRUE,REJECTIONS ',3I5,/
     &'  IA IB IC ID     XA         XB         PT       YC       YD',
     &'       PHI')
      DO 10 N=1,LSCAHD
        PHI = ATAN2(PREC(1,4*N-1),PREC(2,4*N-1))
        WRITE(6,1020) NINHD(N,1),NINHD(N,2),NOUTHD(N,1),NOUTHD(N,2),
     &             XHD(N,1),XHD(N,2),PTHD(N),ETAHD(N,1),ETAHD(N,2),PHI
1020    FORMAT(1X,4I3,2F11.7,4F9.3)
10    CONTINUE
      ENDIF
      IF ( NOUTL.GE.6 ) THEN
C  output of eventrecord
        WRITE(6,1030)
1030    FORMAT('   EVENTRECORD')
        DO 20 L=1,LINE
          WRITE(6,1040) LREC1(L),LREC2(L),(PREC(I,L),I=0,3)
20      CONTINUE
1040    FORMAT(2I12,4(1PE12.4))
        WRITE(6,1050)
1050    FORMAT(/)
      ENDIF
      RETURN
      END
C_____________________________________________________________________
C
C        original title: JTPADI.FOR
C_______________________________________________________________________
*
*  ********************************************************************
      SUBROUTINE JTPDIS(X,QQ,IHATYP,MSPR,PD)
*
*             Parton distributios
*       NPD=ISTRUF (elsewhere) 
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      COMMON /HAPARA/ ECM,PTINI(4),Q0SQR,ALASQR,BQCD,NPD,NF,NHA,NHB
      DIMENSION PD(-6:6)
      DATA ISET / 0 /
      DO 10 I=-6,6
10      PD(I) = 0.0
C  SET MAXFL - THE MAX. NUMBER OF FLAVORS TO CALCULATE
      MAXFL = NF
      IF ( MSPR.EQ.1 .OR. MSPR.EQ.4 ) MAXFL = 0
C------------------------------------------------------
C  CALL ROUTINES CALCULATING THE DISTRIBUTIONS
C  EHLQ 1/2, MRS 1/2/3, GRV, HMRS 1/2, KMRS 1/2/3/4,
C  WHERE HMRS 1/2, KMRS 1/2/3/4 CORRESP. TO HKMRS 1-6
C------------------------------------------------------
      IF     ( NPD.EQ.1 .OR.  NPD.EQ.2 ) THEN
C       CALL PDEHLQ(X,QQ,MAXFL,PD)
        WRITE(6,*) ' unsupported PDF number: ',NPD
      ELSEIF ( NPD.GE.3 .AND. NPD.LE.5 ) THEN
C       CALL PDMRS(X,QQ,MAXFL,PD)
        WRITE(6,*) ' unsupported PDF number: ',NPD
      ELSEIF(NPD.EQ.6)THEN
C       CALL PDGRV(X,QQ,PD)
        WRITE(6,*) ' unsupported PDF number: ',NPD
      ELSEIF(NPD.EQ.7)THEN
C       CALL PHKMRS(X,QQ,PD,1)
        WRITE(6,*) ' unsupported PDF number: ',NPD
      ELSEIF(NPD.EQ.8)THEN
C       CALL PHKMRS(X,QQ,PD,2)
        WRITE(6,*) ' unsupported PDF number: ',NPD
      ELSEIF(NPD.EQ.9)THEN
C       CALL PHKMRS(X,QQ,PD,3)
        WRITE(6,*) ' unsupported PDF number: ',NPD
      ELSEIF(NPD.EQ.10)THEN
C       CALL PHKMRS(X,QQ,PD,4)
        WRITE(6,*) ' unsupported PDF number: ',NPD
      ELSEIF(NPD.EQ.11)THEN
C       CALL PHKMRS(X,QQ,PD,5)
        WRITE(6,*) ' unsupported PDF number: ',NPD
      ELSEIF(NPD.EQ.12)THEN
C       CALL PHKMRS(X,QQ,PD,6)
        WRITE(6,*) ' unsupported PDF number: ',NPD
* updates of April 92, April 93
      ELSEIF((NPD.GE.13).AND.(NPD.LE.20)) THEN
C       CALL PHKMRS(X,QQ,PD,NPD-6)
        WRITE(6,*) ' unsupported PDF number: ',NPD
      ELSEIF((NPD.GE.21).AND.(NPD.LE.23)) THEN
        CALL PHKMRS(X,QQ,PD,NPD-6)
      ELSE
        WRITE(6,*) ' unsupported PDF number: ',NPD
        STOP
      ENDIF
      DO 20 I=-MAXFL,MAXFL
        IF ( PD(I).LT.1.D-15 ) PD(I) = 0.0
20      CONTINUE
C  IF ANTIPROTON CHANGE QUARK <---> ANTIQUARK
      IF ( IHATYP.EQ.-1 ) THEN
        DO 50 I=1,6
          TTTT   = PD(-I)
          PD(-I) = PD(I)
50        PD( I) = TTTT
      ENDIF
      RETURN
      END

C_______________________________________________________________________
      SUBROUTINE PHKMRS(XQ,QQ,PD,MODE)
C***************************************************************C
C                                                               C
C     ORIGINAL NAME: MRSEB( ... )                               C
C                                                               C
C     -----  VARIABLE QUARKS AND GLUONS AT SMALL X ----         C
C                                                               C
C     NEW VERSIONS !!!! JANUARY 1990  (AS DESCRIBED IN          C
C     "PARTON DISTRIBUTIONS ... " P.N. HARRIMAN, A.D. MARTIN,   C
C     R.G. ROBERTS AND W.J. STIRLING PREPRINT DTP-90-04 )       C
C                                                               C
C     NEW VERSIONS !!!! JULY 1990                               C
C     "........................ " J. KWIECINSKI, A.D. MARTIN,   C
C     R.G. ROBERTS AND W.J. STIRLING PREPRINT DTP-90-46 )       C
C      
C****************************************************************
C   All modes 1-14 dropped in dpmjet-II.4.2
C***************************************************************
C  MODE 1 CORRESPONDS TO  HARRIMAN,                             C
C  MARTIN, ROBERTS, STIRLING (EMC FIT)    WITH LAMBDA= 100 MEV  C
C  ORIGINAL: STRC27     NOW: PHMRS1                            C
C                                                               C
C  MODE 2 CORRESPONDS TO  HARRIMAN,                             C
C  MARTIN, ROBERTS, STIRLING (BCDMS FIT)  WITH LAMBDA= 190 MEV  C
C  WITH SMALL X BEHAVIOUR DETERMINED FROM FIT  "HB FIT"         C
C  ORIGINAL: STRC28     NOW: PHMRS2                            C
C                                                               C
C  MODE 3 CORRESPONDS TO  KWIECINSKI,                           C
C  MARTIN, ROBERTS, STIRLING (BCDMS FIT)  WITH LAMBDA= 190 MEV  C
C  AND XG,XQ --> CONSTANT AS X--> 0 AT Q0**2   "B0 FIT"         C
C  ORIGINAL: STRC38     NOW: PKMRS1                            C
C                                                               C
C  MODE 4 CORRESPONDS TO  KWIECINSKI,                           C
C  MARTIN, ROBERTS, STIRLING (BCDMS FIT)  WITH LAMBDA= 190 MEV  C
C  AND XG,XQ --> X**-1/2 AS X--> 0 AT Q0**2    "B- FIT"         C
C  ORIGINAL: STRC48     NOW: PKMRS2                            C
C                                                               C
C  MODE 5 CORRESPONDS TO  KWIECINSKI,                           C
C  MARTIN, ROBERTS, STIRLING (BCDMS FIT)  WITH LAMBDA= 190 MEV  C
C  AND XG,XQ --> X**-1/2 AS X--> 0 AT Q0**2    "B-(5) FIT"      C
C  I.E. WITH WEAK (R=5 GEV-1) SHADOWING INCLUDED              C
C  ORIGINAL: STRC58     NOW: PKMRS3                            C
C                                                               C
C  MODE 6 CORRESPONDS TO  KWIECINSKI,                           C
C  MARTIN, ROBERTS, STIRLING (BCDMS FIT)  WITH LAMBDA= 190 MEV  C
C  AND XG,XQ --> X**-1/2 AS X--> 0 AT Q0**2    "B-(2) FIT"      C
C  I.E. WITH STRONG  (R=2 GEV-1) SHADOWING INCLUDED           C
C  ORIGINAL: STRC68     NOW: PKMRS4                            C
C                                                               C
C****************************************************************
C   All modes 1-14 dropped in dpmjet-II.4.2
C***************************************************************
C                                                               C
C                         -*-                                   C
C                                                               C
C    (NOTE THAT X TIMES THE PARTON DISTRIBUTION FUNCTION        C
C    IS RETURNED I.E. G(X) = GLU/X ETC, AND THAT "SEA"          C
C    IS THE LIGHT QUARK SEA I.E. UBAR(X)=DBAR(X)                C
C    = SEA/X FOR A PROTON.  IF IN DOUBT, CHECK THE              C
C    MOMENTUM SUM RULE! NOTE ALSO THAT SCALE=Q**2 IN GEV**2)    C
C                                                               C
C                         -*-                                   C
C                                                               C
C***************************************************************C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
C      REAL XQ,QQ,PD
C     to get the same as outside
C     in the D IMPLICIT DOUBLE PRECISION(A-H,O-Z) world
      REAL
     & *8
     & XQ,QQ,PD
      DIMENSION PD(-6:6)
      DIMENSION PDFF(-6:2)
      X=DBLE( XQ )
      SCALE=DBLE( QQ )
C-------------------------------------------------------------------
C****************************************************************
C   All modes 1-14 dropped in dpmjet-II.4.2
C***************************************************************
C****************************************************************
C   All modes 1-14 dropped in dpmjet-II.4.2
C***************************************************************
C--------------------------------------------------------------------
C     updates Oct 98 replace DOR94LO by PO_GRV98LO (R.Engel)
      IF((MODE.EQ.15)) THEN
         SCALE2=SCALE
C        CALL DOR94LO(X,SCALE2,UV, DV, DEL, UDB, SB, GL)
         CALL PO_GRV98LO(ISET,X,SCALE2,UV,DV,US,DS,SS,GL)
         CB = 0.D0
         BB = 0.D0
         PD(-5) = BB
         PD(-4) = CB
         PD(-3) = SS
         PD(-2) = US
         PD(-1) = DS
         PD(0)  = GL
         PD(1) = DV+DS
         PD(2)  = UV+US
         PD(3)  = SS
         PD(4)  = PD(-4)
         PD(5)  = PD(-5)
         XQ= X
         QQ= SCALE
C        WRITE(6,*)' PD ',PD
         RETURN
      ENDIF
C     updates Oct 98 replace DOR94LO by PO_GRV98LO (R.Engel)
      IF((MODE.EQ.16)) THEN
         SCALE2=SCALE
C        CALL DOR94LO(X,SCALE2,UV, DV, DEL, UDB, SB, GL)
         CALL PO_GRV98LO(ISET,X,SCALE2,UV,DV,US,DS,SS,GL)
         CB = 0.D0
         BB = 0.D0
         PD(-5) = BB
         PD(-4) = CB
         PD(-3) = SS
         PD(-2) = US
         PD(-1) = DS
         PD(0)  = GL
         PD(1) = DV+DS
         PD(2)  = UV+US
         PD(3)  = SS
         PD(4)  = PD(-4)
         PD(5)  = PD(-5)
         XQ= X
         QQ= SCALE
C        WRITE(6,*)' PD ',PD
         RETURN
      ENDIF
C--------------------------------------------------------------------
C     updates Feb. 97
      IF((MODE.EQ.17)) THEN
      CALL structm(x,SCALE,Upv,Dnv,Usea,Dsea,Str,Chm,Bot,Top,Glu)
      PD(0)=  GLU
      PD(1)=  USEA+UPV
      PD(2)=  DSEA+DNV
      PD(3)=  STR
      PD(4)=  CHM
      PD(5)=  BOT
      PD(-5)= BOT
      PD(-4)= CHM
      PD(-3)= STR
      PD(-2)= DSEA
      PD(-1)= USEA 
      XQ= X
      QQ= SCALE 
        RETURN
      ENDIF
C--------------------------------------------------------------------
C--------------------------------------------------------------------
      PD(0)=  GLU
      PD(1)=  SEA+UPV
      PD(2)=  SEA+DNV
      PD(3)=  STR
      PD(4)=  CHM
      PD(5)=  BOT
      PD(-5)= BOT
      PD(-4)= CHM
      PD(-3)= STR
      PD(-2)= SEA
      PD(-1)= SEA
      XQ= X
      QQ= SCALE
C--------------------------------------------------------------------
      RETURN
      END
C
C*********************************************************************
C-----original seperate file with the name------------------------------------------------------------------
C                            DTULPTPE.FOR
C
C------------------------------------------------------------------------
C_______________________________________________________________________
C
C   PROGRAM FOR SIMULATION OF HARD SCATTERING OF HADRONIC PARTICLES
C
C     AUTHOR       K.HAHN            LEIPZIG   GDR
C_______________________________________________________________________
      SUBROUTINE LAPTAB
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( MAXPRO = 8 , MLINE = 1000 , MSCAHD = 250 )
      COMMON /HACONS/ PI,PI2,PI4,GEVTMB
      COMMON /HAPARA/ ECM,PTINI(4),Q0SQR,ALASQR,BQCD,NPD,NF,NHA,NHB
      COMMON /HAPADI/ NPDM
      COMMON /HAPDCO/ NPDCOR
      COMMON /HAQQAP/ AQQAL,AQQPD,NQQAL,NQQPD
      COMMON /HAENVI/ NINDEP
      COMMON /HAOUTL/ NOUTL,NOUTER,NOUTCO
      COMMON /HACUTS/ PTL,PTU,ETACL,ETACU,ETADL,ETADU
      COMMON /HAGAUP/ NGAUP1,NGAUP2,NGAUET,NGAUIN
      COMMON /HAEVTR/ LINE,LIN,LREC1(MLINE),LREC2(MLINE),PREC(0:3,MLINE)
      COMMON /HARSLT/ LSCAHD,LSC1HD,
     &                ETAHD(MSCAHD,2) ,PTHD(MSCAHD),
     &                XHD(MSCAHD,2)   ,VHD(MSCAHD) ,X0HD(MSCAHD,2),
     &                NINHD(MSCAHD,2) ,NOUTHD(MSCAHD,2),
     &                N0INHD(MSCAHD,2),NBRAHD(MSCAHD,2),NPROHD(MSCAHD)
      INTEGER MXSECT
      COMMON /HAXSEC/ XSECTA(2,-1:MAXPRO,4),XSECT(5,-1:MAXPRO),
     &                MXSECT(0:2,-1:MAXPRO)
      CHARACTER*8  CW
      CHARACTER*79 COMMNT
      CHARACTER*70 INP
      DIMENSION   WHAT(6)
C======================================================================
C  THE AVAILABLE CODE WORDS ARE
C
C        END     , COMMENT , ENERGYPT, PARDISTR, CUTS
C        INTPOINT, FLAVOR  , PARTICLE, OUTPUT  , INIT    ,
C        TESTINCL, TESTMC  , SUBPRON , SUBPROFF, HISOUT  ,
C        HISINI  , HARDSCAL, USER    , PARDISCO
C_______________________________________________________________________
      WRITE(6,2000)
2000  FORMAT('1***************************************************'
     & ,/,   ' MONTE-CARLO GENERATION OF HARD HADRONIC SCATTERINGS'
     & ,/,   ' ***************************************************',/)
      CALL JTDTU(1)
C  read a control card
10    CONTINUE
      READ (5,1010) INP
      IF ( INP(1:1).EQ.'-' ) GOTO 10
      WRITE(6,1011) INP
      READ(INP,1012,ERR=99) CW,WHAT
      GOTO 15
99    WRITE(6,1013)
      GOTO 10
1010  FORMAT(A70)
1011  FORMAT(' *********.* CONTROL.CARD*****.',4(9X,'.'),/,1X,A70,/)
1012  FORMAT(A8,2X,6E10.0)
1013  FORMAT('     CARD IS INCORRECT, IGNORE AND TRY NEXT CARD',/)
15    CONTINUE
C======================================================================
C======================================================================
      IF     ( CW.EQ.'END     ' ) THEN
C
C       STOPS THE PROGRAM
C_______________________________________________________________________
        WRITE(6,1030)
1030    FORMAT(' ******** END OF PROGRAM EXECUTION ********')
        RETURN
C======================================================================
      ELSEIF ( CW.EQ.'COMMENT ' ) THEN
C
C       TO INCLUDE COMMENTS IN PROGRAM OUTPUT
C  WHAT(1) - # OF CARDS FOLLOWING THE CONTROLCARD AND
C            CONTAINING THE COMMENT
C_______________________________________________________________________
        N = MAX(1,INT(WHAT(1)))
        DO 20 I=1,N
          READ(5,1040) COMMNT
20        WRITE(6,1050) COMMNT
1040    FORMAT(A79)
1050    FORMAT(1X,A79)
C======================================================================
      ELSEIF ( CW.EQ.'ENERGYPT' ) THEN
C
C       READ CMS-ENERGY AND MIN. TRANSVERSE MOMENTUM FOR JETS
C  WHAT(1) - CMS-ENERGY   ( IN GEV )       DEFAULT  540.
C  WHAT(2) - MIN. PT      ( IN GEV )       DEFAULT    2.
C  WHAT(3) - MIN. PT      ( IN GEV )       DEFAULT    2.
C  WHAT(4) - MIN. PT      ( IN GEV )       DEFAULT    2.
C  WHAT(5) - MIN. PT      ( IN GEV )       DEFAULT    2.
C
C_______________________________________________________________________
        IF ( WHAT(1).GT.0.0D0 ) ECM = WHAT(1)
        DO 22 I=1,4
          PTINI(I) = WHAT(I+1)
22      CONTINUE
C======================================================================
      ELSEIF ( CW.EQ.'PARDISTR' ) THEN
C
C       DEFINE THE PARTON DISTRIBUTION SET
C  WHAT(1) - NPD                               DEFAULT  3.
C       NPD = 1,2     : EICHTEN,HINCHLIFFE,LANE,QUIGG
C           = 3,4,5   : MARTIN,ROBERTS,STIRLING
C           = 6       : GLUECK,REYA,VOGT
C           = 7,8     : HARRIMAN,MARTIN,ROBERTS,STIRLING
C           = 9 - 12  : KWIECINSKI,MARTIN,ROBERTS,STIRLING
C
C  WHAT(2) - NPDM                              DEFAULT  0.
C       NPDM = 0,1    : CORRECTS THE PARTON DISTRIBUTION IN CASE
C                       OF NPD = 10 TO THE 1/SQRT(X) BEHAVIOUR
C                       FOR X < XMIN=1.E-05 IF NPDM = 1 AND NOT
C                       IF NPDM = 0 ( M STANDS FOR MODIFICATION )
C_______________________________________________________________________
        IPD      = INT(WHAT(1))
        IPDM     = INT(WHAT(2))
        NPD      = 3
        NPDM     = 0
*       IF ( IPD.GE.1 .AND. IPD.LE.12 ) NPD = IPD
        IF ( IPD.GE.1 .AND. IPD.LE.15 ) NPD = IPD
        IF ( IPDM.EQ.1 ) NPDM = IPDM
C======================================================================
      ELSEIF ( CW.EQ.'CUTS    ' ) THEN
C
C       set kinematic cuts on partonlevel
C  WHAT(1) - PTL     min. pt for parton               DEFAULT  0.0
C  WHAT(2) - PTU     max. pt for parton               DEFAULT  1.E+30
C  WHAT(3) - ETACL   min. rapidity for parton c       DEFAULT -1.E+30
C  WHAT(4) - ETACU   max. rapidity for parton c       DEFAULT  1.E+30
C  WHAT(5) - ETADL   min. rapidity for parton d       DEFAULT -1.E+30
C  WHAT(6) - ETADU   max. rapidity for parton d       DEFAULT  1.E+30
C
C    this cuts are additional cuts which not affect the pt-cut given
C    in the ENERGYPT card;
C    this cuts are checked during event generation and events violating
C    the cuts are rejected;
C    the defaults are such that there is really no cutting;
C_______________________________________________________________________
        PTL     = WHAT(1)
        PTU     = WHAT(2)
        ETACL   = WHAT(3)
        ETACU   = WHAT(4)
        ETADL   = WHAT(5)
        ETADU   = WHAT(6)
        IF ( PTU  .LE.PTL   ) PTU   = PTL  +1.0
        IF ( ETACU.LE.ETACL ) ETACU = ETACL+1.0
        IF ( ETADU.LE.ETADL ) ETADU = ETADL+1.0
C======================================================================
      ELSEIF ( CW.EQ.'INTPOINT' ) THEN
C
C       define number of integration points
C  WHAT(1) - NGAUP1                            DEFAULT  8.
C  WHAT(2) - NGAUP2                            DEFAULT  8.
C  WHAT(3) - NGAUET                            DEFAULT  8.
C  WHAT(4) - NGAUIN                            DEFAULT  8.
C
C    NGAUP1,NGAUP2,NGAUET : used in inclusive x-section calculation;
C                           first, second pt-integration and
C                           eta-integration respectivly
C    NGAUIN               : used in initialization ( SUBROUTINE HABINT )
C
C       max. value allowed for integration is   32 ! - For NGAUP1,P2,ET
C       max. value allowed for integration is 1000 ! - For NGAUIN.
C                                                    (4.6.91)
C_______________________________________________________________________
      IF ( WHAT(1).GE.1.D0.AND.WHAT(1).LE.32.D0) NGAUP1 = INT(WHAT(1))
      IF ( WHAT(2).GE.1.D0.AND. WHAT(2).LE.32.D0) NGAUP2 = INT(WHAT(2))
      IF ( WHAT(3).GE.1.D0.AND. WHAT(3).LE.32.D0) NGAUET = INT(WHAT(3))
      IF(WHAT(4).GE.1.D0.AND. WHAT(3).LE.1000.D0) NGAUIN = INT(WHAT(4))
C======================================================================
      ELSEIF ( CW.EQ.'FLAVOR  ' ) THEN
C       DEFINES ACTIVE FLAVORS
C  WHAT(1) - # OF FLAVORS                 DEFAULT  4
C_______________________________________________________________________
        NFF     = INT(WHAT(1))
        IF ( NFF.GE.0  .AND.  NFF .LE.6 ) NF = NFF
C======================================================================
      ELSEIF ( CW.EQ.'PARTICLE' ) THEN
C
C       TARGET AND BEAM PARTICLES
C  WHAT(1) - BEAM    ( POSITIVE Z-DIRECTION )    DEFAULT  1
C  WHAT(2) - TARGET                              DEFAULT  1
C     1 : PROTON
C    -1 : ANTIPROTON
C_______________________________________________________________________
        IHA = INT(WHAT(1))
        IF ( ABS(IHA).EQ.1 ) NHA = IHA
        IHB = INT(WHAT(2))
        IF ( ABS(IHB).EQ.1 ) NHB = IHB
C======================================================================
      ELSEIF ( CW.EQ.'OUTPUT  ' ) THEN
C
C       set output level
C  WHAT(1) - NOUTL    output level  0....                    DEFAULT  1.
C  WHAT(2) - NOUTER   error messages ( 0 - no  ; 1 - yes )   DEFAULT  1.
C_______________________________________________________________________
        IF ( WHAT(1).GE.0.D0 )                NOUTL  = INT(WHAT(1))
        IF ( WHAT(2).EQ.0.D0.OR.WHAT(2).EQ.1.D0)NOUTER = INT(WHAT(2))
        IF ( WHAT(3).GE.0.D0 )                 NOUTCO = INT(WHAT(3))
C======================================================================
      ELSEIF ( CW.EQ.'INIT    ' ) THEN
C
C       DEMANDS A NEW INITIALIZATION
C_______________________________________________________________________
        CALL HARINI
C======================================================================
      ELSEIF ( CW.EQ.'TESTINCL' ) THEN
C
C       TEST OF INCLUSIVE JET PRODUCTION
C         PLOT OF PARTONDISTRIBUTIONS USED ( GLUONS )
C         PLOT OF DIFFERENTIAL JET CROSS SECTIONS
C
C_______________________________________________________________________
        DO 35 I=1,6
          J   = INT(WHAT(I))
          IF ( J.GE.1 .AND. J.LE.4 ) CALL HATEST(J)
35        CONTINUE
C======================================================================
      ELSEIF ( CW.EQ.'TESTMC  ' ) THEN
C
C       TEST OF MONTE-CARLO JET PRODUCTION
C         PLOT OF DIFFERENTIAL JET CROSS SECTIONS
C  WHAT(1) - # OF EVENTS TO PRODUCE IN TEST       DEFAULT  100
C  WHAT(2) - # OF HARD SCATTERINGS AT HARD EVENT
C  WHAT(3) - MIN. PT FOR FIRST HARD SCATTERING
C
C_______________________________________________________________________
        NEVT = INT(WHAT(1))
        IF ( NEVT.LE.0 ) NEVT = 100
        NHARD = MAX(1,INT(WHAT(2)))
        PT1   = WHAT(3)
        CALL TIMDAT
        DO 36 I=1,NEVT
          MHARD = NHARD
          CALL HAREVT(MHARD,PT1)
36      CONTINUE
        CALL TIMDAT
C======================================================================
      ELSEIF ( CW.EQ.'SUBPRON ' ) THEN
C
C   SWITCH SUBPROCESSES ON
C
C______________________________________________________________________
        DO 40 I=1,6
          M = INT(WHAT(I))
          IF ( M.GE.1  .AND.  M.LE.MAXPRO ) MXSECT(0,M) = 1
40        CONTINUE
        MXSECT(0,-1) = MXSECT(0,3)
C======================================================================
      ELSEIF ( CW.EQ.'SUBPROFF' ) THEN
C
C   SWITCH SUBPROCESSES OFF
C
C_______________________________________________________________________
        DO 50 I=1,6
          M = INT(WHAT(I))
          IF ( M.GE.1  .AND.  M.LE.MAXPRO ) MXSECT(0,M) = 0
50        CONTINUE
        MXSECT(0,-1) = MXSECT(0,3)
C======================================================================
      ELSEIF ( CW.EQ.'HISOUT  ' ) THEN
C
C   OUTPUT OF MC RESULTS
C    WHAT()- 1   : TOTAL CROSS SECTION
C    WHAT()- 2   : PT-DISTR.
C    WHAT()- 3   : PT-DISTR.
C    WHAT()- 4   : PT-DISTR.
C    WHAT()- 5   : RAPIDITY-DISTR.
C    WHAT()- 6   : RAPIDITY-DISTR.
C
C_______________________________________________________________________
        DO 60 I=1,6
          J   = INT(WHAT(I))
          IF ( J.GE.1 .AND. J.LE.6 ) CALL HISOUT(J)
60        CONTINUE
C======================================================================
      ELSEIF ( CW.EQ.'HISINI  ' ) THEN
C
C   INITIALIZE THE HISTOGRAMS FOR MC TESTING
C_______________________________________________________________________
        CALL HISINI
C======================================================================
      ELSEIF ( CW.EQ.'HARDSCAL' ) THEN
C
C   define hard scale
C     WHAT(1) - NQQAL  : definition of hard scale in coupling constant
C     WHAT(2) - AQQAL  : factor multiplied to hard scale in coupling
C     WHAT(3) - NQQPD  : definition of hard scale in parton distr.
C     WHAT(4) - AQQPD  : factor multiplied to hard scale in parton distr.
C        the possible NQQAL(PD) are
C           1 :  QQ = AQQAL(PD) * PT**2
C           2 :  QQ = AQQAL(PD) * SP
C           3 :  QQ = AQQAL(PD) * (SP*TP*UP)**(1./3.)
C           4 :  QQ = AQQAL(PD) * (SP*TP*UP)/(SP**2+TP**2+UP**2)
C
C      DEFAULT is NQQAL = NQQPD = 1 and AQQAL = AQQPD = 1./4.
C_______________________________________________________________________
        IF ( WHAT(1).GE.1.D0.AND.WHAT(1).LE.4.D0)NQQAL = INT(WHAT(1))
        IF ( WHAT(2).GT.0.D0 )                   AQQAL =     WHAT(2)
        IF ( WHAT(3).GE.1.D0.AND.WHAT(3).LE.4.D0)NQQPD = INT(WHAT(3))
        IF ( WHAT(4).GT.0.D0 )                   AQQPD =     WHAT(4)
C    C======================================================================
C          ELSEIF ( CW.EQ.'USER    ' ) THEN
C    C
C    C_______________________________________________________________________
C            WRITE(6,9050)
C    9050    FORMAT(' -----> GIVE CONTROL TO USER ROUTINE')
C            CALL HAUSER(WHAT)
C            WRITE(6,9060)
C    9060    FORMAT(' -----> CONTROL COMES BACK FROM USER ROUTINE')
C======================================================================
      ELSEIF ( CW.EQ.'PARDISCO' ) THEN
C
C   define partoncorrelationfunction
C     WHAT(1) - NPDCOR : =0 no correlations, =1 correlations
C
C      DEFAULT is NPDCOR = 0
C_______________________________________________________________________
        IF ( WHAT(1).EQ.0.D0.OR.WHAT(1).EQ.1.D0)NPDCOR = INT(WHAT(1))
C======================================================================
      ELSE
        WRITE(6,9999)
9999    FORMAT(' ##### UNKNOWN CODEWORD;  CARD IS IGNORED ###',/)
      ENDIF
      GOTO 10
      END
C
C-----originally seperate file with the name ---------------------------
C                             JTINCL.FOR
C
C______________________________________________________________________
C
C  PROCEDURES FOR CALCULATION OF INCLUSIVE CROSS SECTIONS
C
C______________________________________________________________________
C______________________________________________________________________
      SUBROUTINE CSJ2M(PT,ETAC,ETAD,DSIGMM)
C    CALCULATION OF DIFFERENTIAL CROSS SECTION DSIG/(DETAC*DETAD*DPT)
C    FOR DIFFERENT SUBPROCESSES
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( MAXPRO = 8 , MLINE = 1000 , MSCAHD = 250 )
      PARAMETER ( TINY= 1.D-30, ONEP1=1.1D0 ,TINY6=1.D-06)
      COMMON /HACONS/ PI,PI2,PI4,GEVTMB
      COMMON /HAPARA/ ECM,PTINI(4),Q0SQR,ALASQR,BQCD,NPD,NF,NHA,NHB
      COMMON /HAQQAP/ AQQAL,AQQPD,NQQAL,NQQPD
      DOUBLE PRECISION EC,ED,XA,XB,SP,TP,UP,TT,UU,
     *       FACTOR,DSIGM(0:MAXPRO)
      DIMENSION DSIGMM(0:MAXPRO),PDA(-6:6),PDB(-6:6)

      DO 10 I=0,MAXPRO
        DSIGM(I) = 0.0
10      CONTINUE
      EC     = EXP(ETAC)
      ED     = EXP(ETAD)
C  KINETICS
      XA     = PT*(EC+ED)/ECM
      XB     = XA/(EC*ED)
      IF ( XA.GE.1.D0 .OR. XB.GE.1.D0 ) RETURN
      SP     = XA*XB*ECM*ECM
      UP     =-ECM*PT*EC*XB
      UP     = UP/SP
      TP     =-(1.+UP)
      UU     = UP*UP
      TT     = TP*TP
C  set hard scale  QQ  for alpha and partondistr.
      IF     ( NQQAL.EQ.1 ) THEN
        QQAL = AQQAL*PT*PT
      ELSEIF ( NQQAL.EQ.2 ) THEN
        QQAL = AQQAL*SP
      ELSEIF ( NQQAL.EQ.3 ) THEN
        QQAL = AQQAL*SP*(UP*TP)**(1./3.)
      ELSEIF ( NQQAL.EQ.4 ) THEN
        QQAL = AQQAL*SP*UP*TP/(1.+TT+UU)
      ENDIF
      IF     ( NQQPD.EQ.1 ) THEN
        QQPD = AQQPD*PT*PT
      ELSEIF ( NQQPD.EQ.2 ) THEN
        QQPD = AQQPD*SP
      ELSEIF ( NQQPD.EQ.3 ) THEN
        QQPD = AQQPD*SP*(UP*TP)**(1./3.)
      ELSEIF ( NQQPD.EQ.4 ) THEN
        QQPD = AQQPD*SP*UP*TP/(1.+TT+UU)
      ENDIF
C      PARAMETER ( TINY= 1.D-30, ONEP1=1.1D0 ,TINY6=1.D-06)
      ALPHA  = BQCD/LOG(MAX(QQAL/ALASQR,ONEP1))
      FACTOR = PI2*GEVTMB*PT*(ALPHA/SP)**2
C   PARTONDISTRIBUTIONS  ( MULTIPLIED BY X )
      X1 = XA
      X2 = XB
      CALL JTPDIS(X1,QQPD,NHA,0,PDA)
      CALL JTPDIS(X2,QQPD,NHB,0,PDB)
      S1    = PDA(0)*PDB(0)
      S2    = 0.0
      S3    = 0.0
      S4    = 0.0
      S5    = 0.0
      DO 20 I=1,NF
        S2  = S2+PDA(I)*PDB(-I)+PDA(-I)*PDB( I)
        S3  = S3+PDA(I)*PDB( I)+PDA(-I)*PDB(-I)
        S4  = S4+PDA(I)+PDA(-I)
        S5  = S5+PDB(I)+PDB(-I)
20    CONTINUE
C   CROSS-SECTIONS ( INCLUDING MATRIX ELEMENTS, SYMMETRY-FACTORS AND
C                                FACTORS FOR FINALSTATE-SUMMATION )
      DSIGM(1) = 2.25*(3.-((UP*TP)+UP/TT+TP/UU))
      DSIGM(6) = (4./9.)*(UU+TT)
      DSIGM(8) = (4./9.)*(1.+UU)/TT
      DSIGM(2) = (16./27.)*(UU+TT)/(UP*TP)-3.*DSIGM(6)
      DSIGM(3) = ((1.+UU)/TT)-(4./9.)*(1.+UU)/UP
      DSIGM(4) = (9./32.)*DSIGM(2)
      DSIGM(5) = DSIGM(6)+DSIGM(8)-(8./27.)*UU/TP
      DSIGM(7) = 0.5*(DSIGM(8)+(4./9.)*(1.+TT)/UU-(8./27.)/(UP*TP))

      DSIGM(1) = FACTOR*DSIGM(1)*S1
      DSIGM(2) = FACTOR*DSIGM(2)*S2
      DSIGM(3) = FACTOR*DSIGM(3)*(PDA(0)*S5+PDB(0)*S4)
      DSIGM(4) = FACTOR*DSIGM(4)*S1*NF
      DSIGM(5) = FACTOR*DSIGM(5)*S2
      DSIGM(6) = FACTOR*DSIGM(6)*S2*MAX(0,(NF-1))
      DSIGM(7) = FACTOR*DSIGM(7)*S3
      DSIGM(8) = FACTOR*DSIGM(8)*(S4*S5-(S2+S3))
C   sum over processes
      DO 40 M=1,MAXPRO
        DSIGM(0) = DSIGM(0)+DSIGM(M)
40    CONTINUE
      DO 50 M=0,MAXPRO
        DSIGMM(M) = DSIGM(M)
50    CONTINUE
      RETURN
      END
C______________________________________________________________________
      SUBROUTINE CSJ1M(PT,ETAC,DSIGM)
C  ONE JET CROSS SECTION
C   ( CSJ2M INTEGRATED OVER ETAD )
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( MAXPRO = 8 , MLINE = 1000 , MSCAHD = 250 )
C     PARAMETER ( TINY= 1.D-30 )
C     CHANGED 14.10.92 AFTER COMPARISON WITH THE ORIGINAL
C     VERSION OF DTULAP FROM 1.D-30 BACK TO THE ORIGINAL 1.D-20
      PARAMETER ( TINY= 1.D-20 )
      COMMON /HAPARA/ ECM,PTINI(4),Q0SQR,ALASQR,BQCD,NPD,NF,NHA,NHB
      COMMON /HAGAUP/ NGAUP1,NGAUP2,NGAUET,NGAUIN
      DIMENSION DSIGM(0:MAXPRO),DSIG1(0:MAXPRO)
      DIMENSION ABSZ(32),WEIG(32)

      DO 10 M=0,MAXPRO
        DSIGM(M) = 0.0
10    CONTINUE
      EC  = EXP(ETAC)
      ARG = ECM/PT
      IF  ( ARG.LE.EC .OR. ARG.LE.1./EC ) RETURN
      EDU = LOG(ARG-EC)
      EDL =-LOG(ARG-1./EC)
      NPOINT = NGAUET
      CALL GSET(EDL,EDU,NPOINT,ABSZ,WEIG)
      DO 30 I=1,NPOINT
        CALL CSJ2M(PT,ETAC,ABSZ(I),DSIG1)
        DO 20 M=0,MAXPRO
C         PCTRL= DSIG1(M)/1.E-20
          PCTRL= DSIG1(M)/TINY
          PCTRL= ABS( PCTRL )
          IF( PCTRL.GE.1.D0 ) THEN
            DSIGM(M) = DSIGM(M)+WEIG(I)*DSIG1(M)
          ENDIF
20      CONTINUE
30    CONTINUE
      RETURN
      END
C______________________________________________________________________
      SUBROUTINE CSJ1MI(PT,DSIGM)
C  ONE JET CROSS SECTION
C   ( CSJ2M INTEGRATED OVER ETAD AND ETAC )
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( MAXPRO = 8 , MLINE = 1000 , MSCAHD = 250 )
      COMMON /HAPARA/ ECM,PTINI(4),Q0SQR,ALASQR,BQCD,NPD,NF,NHA,NHB
      COMMON /HAGAUP/ NGAUP1,NGAUP2,NGAUET,NGAUIN
      DIMENSION DSIGM(0:MAXPRO),DSIG1(0:MAXPRO)
      DIMENSION ABSZ(32),WEIG(32)

      DO 10 M=0,MAXPRO
        DSIGM(M) = 0.0
10    CONTINUE
      AMT = 2.*PT/ECM
      IF ( AMT.GE.1.D0 ) RETURN
      ECU = LOG((SQRT(1.-AMT*AMT)+1.)/AMT)
      ECL =-ECU
      NPOINT = NGAUET
      CALL GSET(ECL,ECU,NPOINT,ABSZ,WEIG)
      DO 30 I=1,NPOINT
        CALL CSJ1M(PT,ABSZ(I),DSIG1)
        DO 20 M=0,MAXPRO
          DSIGM(M) = DSIGM(M)+WEIG(I)*DSIG1(M)
20      CONTINUE
30    CONTINUE
      RETURN
      END
C______________________________________________________________________
      SUBROUTINE CSHARM(DSIGM)
C
C  TOTAL HARD CROSS SECTION
C   ( CSJ2M INTEGRATED OVER PT, ETAD AND ETAC )
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( MAXPRO = 8 , MLINE = 1000 , MSCAHD = 250 )
      COMMON /HAPARA/ ECM,PTINI(4),Q0SQR,ALASQR,BQCD,NPD,NF,NHA,NHB
      COMMON /HAGAUP/ NGAUP1,NGAUP2,NGAUET,NGAUIN
      COMMON /XSECPT/ PTCUT,SIGS,DSIGH
      DIMENSION DSIGM(0:MAXPRO),DSIG1(0:MAXPRO)
      DIMENSION ABSZ(32),WEIG(32)
      DATA FAC / 3.0 /

      DO 10 M=0,MAXPRO
        DSIGM(M) = 0.0
10    CONTINUE
      IF ( PTINI(1).GE.ECM/2.D0 ) RETURN
      PTMIN  = PTINI(1)
      PTCUT  = PTMIN
      PTMAX  = MIN(FAC*PTMIN,ECM/2.)
      NPOINT = NGAUP1
      CALL CSJ1MI(PTMIN,DSIG1)
      SIG1   = DSIG1(0)
C     WRITE(6,1000) SIG1
 1000 FORMAT(1X,' d sigma/ p_t d p_t ',E12.5)
      DSIGH  = SIG1
      PTMXX  = 0.95*PTMAX
      CALL CSJ1MI(PTMXX,DSIG1)
      EX     = LOG(SIG1/(DSIG1(0)+1.E-30))/LOG(FAC)
      EX1    = 1.0-EX
      DO 50 K=1,2
        IF ( PTMIN.GE.PTMAX ) GOTO 40
        RL   = PTMIN**EX1
        RU   = PTMAX**EX1
        CALL GSET(RL,RU,NPOINT,ABSZ,WEIG)
        DO 30 I=1,NPOINT
          R  = ABSZ(I)
          PT = R**(1.0/EX1)
          CALL CSJ1MI(PT,DSIG1)
          F  = WEIG(I)*PT/(R*EX1)
          DO 20 M=0,MAXPRO
            DSIGM(M) = DSIGM(M)+F*DSIG1(M)
20        CONTINUE
30      CONTINUE
40      PTMIN  = PTMAX
        PTMAX  = ECM/2.0
        NPOINT = NGAUP2
50    CONTINUE
      RETURN
      END
C-------originally seperate file ----------------------------------------
C                        JTHIST.FOR
C
C________________________________________________________________________

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PROCEDURES TO PRINT OUT HISTOGRAMS AND INCLUSICE TESTCALCULATIONS
C     ( AND TO INITIALIZE AND FILL HISTOGRAMS ALSO )
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C______________________________________________________________________
      SUBROUTINE HISOUT(IOUT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( MAXPRO = 8 , MLINE = 1000 , MSCAHD = 250 )
      PARAMETER ( TINY= 1.D-30, ONE=1.D0 ,TINY6=1.D-06,ZERO=0.)
      COMMON /HAPARA/ ECM,PTINI(4),Q0SQR,ALASQR,BQCD,NPD,NF,NHA,NHB
      INTEGER MXSECT
      COMMON /HAXSEC/ XSECTA(2,-1:MAXPRO,4),XSECT(5,-1:MAXPRO),
     &                MXSECT(0:2,-1:MAXPRO)
      CHARACTER*18 PROC
      CHARACTER*11 PDSET,PARTIC
      COMMON /PEPROC/ PROC(0:MAXPRO),PDSET(23),PARTIC(-1:1)
C  LENGTH OF HISTO : 15000  REAL*4
      COMMON /HISTO / PT10,DPT1,ETA10,DETA1,PT20,DPT2,ETA20,DETA2,
     &                X(50,-5:5),AB(50,-5:5),HPE(50,-5:5),HEP(50,5),
     &                HPM(50,8),HEM(50,8),HP(50),HE(50),
     &                SIG(MAXPRO),STDEV(MAXPRO),
     &                FILL(12176)
C
      XSMIN =-6.
      XSSTEP= 0.1
      DO 5 J=-5,5
        DO 5 I=1,50
5         X(I,J)  =-100.
      XSECT (2,0) = XSECT (2,-1)
      MXSECT(1,0) = MXSECT(1,-1)
      MXSECT(2,0) = MXSECT(2,-1)
      DO 7 M=1,MAXPRO
        MXSECT(1,0) = MXSECT(1,0)+MXSECT(1,M)
        MXSECT(2,0) = MXSECT(2,0)+MXSECT(2,M)
7       XSECT (2,0) = XSECT(2,0)+XSECT(2,M)
      WRITE(6,1010) IOUT
1010  FORMAT(1X,20('=='),' HISTO-OUTPUT ',I2,1X,10('=='),/)
      IF ( IOUT.EQ.1 ) THEN
        WRITE(6,1040)
1040    FORMAT('      PROCESS',15X,'EVENTS',22X,'HARD CROSS SECTION',/,
     &         25X,'TOTAL  ACCEPT.',10X,'MONTE-CARLO',11X,'INCLUSIVE')
        SIGSUM = 0.0
        STDEVS = 0.0
        DO 20 M=1,MAXPRO
          SIG(M)   = 0.0
          STDEV(M) = 0.0
          IF ( MXSECT(1,M).GT.0 ) THEN
            SIG(M)   = XSECT(3,M)/MXSECT(1,M)
            STDEV(M) = SQRT(MAX(ZERO,XSECT(4,M)-XSECT(3,M)*SIG(M)))/
     &                                                     MXSECT(1,M)
          ENDIF
          IF ( M.EQ.3 .AND. MXSECT(1,-1).GT.0 ) THEN
            SIGG     = XSECT(3,-1)/MXSECT(1,-1)
            SIG(3)   = SIG(3)+SIGG
            STDEV(3) = STDEV(3)
     &       +SQRT(MAX(ZERO,XSECT(4,-1)-XSECT(3,-1)*SIGG))/MXSECT(1,-1)
          ENDIF
          SIGSUM = SIGSUM+SIG(M)
          STDEVS = STDEVS+STDEV(M)
20      CONTINUE
        MXSECT(1,3) = MXSECT(1,3)+MXSECT(1,-1)
        MXSECT(2,3) = MXSECT(2,3)+MXSECT(2,-1)
        WRITE(6,1050) PROC(0),(MXSECT(L,0),L=0,2),
     &                SIGSUM,STDEVS,XSECT(5,0)
        DO 25 M=1,MAXPRO
          IF ( MXSECT(0,M).EQ.1 ) WRITE(6,1050) PROC(M),
     &              (MXSECT(L,M),L=0,2),SIG(M),STDEV(M),XSECT(5,M)
25        CONTINUE
1050      FORMAT(A19,I3,2I8,E14.4,' +- ',E10.4,E14.4)
        MXSECT(1,3) = MXSECT(1,3)-MXSECT(1,-1)
        MXSECT(2,3) = MXSECT(2,3)-MXSECT(2,-1)
      ELSEIF ( IOUT.EQ.2 ) THEN
        FAC = XSECT(2,0)/(DPT1*MXSECT(1,0))
        DO 30 I=1,50
          AB(I,1) = PT10+(I-1)*DPT1
          IF ( HP(I).GT.1.D-35 ) X(I,1) = LOG10(FAC*HP(I))
30        CONTINUE
        WRITE(6,1060)
1060    FORMAT('  JET CROSS SECTION  PT-DISTRIBUTION',/)
        CALL PLOT(AB(1,1),X(1,1),50,1,50,PT10,DPT1,XSMIN,XSSTEP)
      ELSEIF ( IOUT.EQ.3 ) THEN
        FAC = XSECT(2,0)/(DPT1*MXSECT(1,0))
        DO 50 I=1,50
          PT = PT10+(I-1)*DPT1
          DO 40 J=1,8
            AB(I,J-6)  = PT
            IF ( HPM(I,J).GT.1.D-35 ) X(I,J-6) = LOG10(FAC*HPM(I,J))
40          CONTINUE
50        CONTINUE
        WRITE(6,1070)
1070    FORMAT('  JET CROSS SECTION  PT-DISTRIBUTION',/,
     &         '                     FOR THE DIFF. SUBPROCESSES',/)
        CALL PLOT(AB,X,400,8,50,PT10,DPT1,XSMIN,XSSTEP)
      ELSEIF ( IOUT.EQ.4 ) THEN
        FAC = XSECT(2,0)/(DPT1*DETA1*MXSECT(1,0))
        DO 70 I=1,50
          PT = PT10+(I-1)*DPT1
          DO 60 J=-5,5
            AB(I,J)  = PT
            IF ( HPE(I,J).GT.1.D-35 ) X(I,J) = LOG10(FAC*HPE(I,J))
60          CONTINUE
70        CONTINUE
        WRITE(6,1080) ETA10,-ETA10
1080    FORMAT('  JET CROSS SECTION  PT-DISTRIBUTION',/,
     &         '                     RAP.=',F5.2,'...',F4.2,/)
        CALL PLOT(AB,X,550,11,50,PT10,DPT1,XSMIN,XSSTEP)
      ELSEIF ( IOUT.EQ.5 ) THEN
        FAC = XSECT(2,0)/(DETA2*DPT2*MXSECT(1,0))
        DO 80 I=1,50
          ETA = ETA20+(I-1)*DETA2
          DO 75 J=1,5
            AB(I,J) = ETA
            IF ( HEP(I,J).GT.1.D-35 ) X(I,J) = LOG10(FAC*HEP(I,J))
75          CONTINUE
80        CONTINUE
        WRITE(6,1090) PT20,PT20+4.*DPT2
1090    FORMAT('  JET CROSS SECTION   RAP.-DISTRIBUTION',/,
     &         '                      PT=',F6.2,'...',F6.2,/)
        CALL PLOT(AB(1,1),X(1,1),250,5,50,ETA20,DETA2,XSMIN,XSSTEP)
      ELSEIF ( IOUT.EQ.6 ) THEN
        FAC = XSECT(2,0)/(DETA2*MXSECT(1,0))
        DO 100 I=1,50
          ETA = ETA20+(I-1)*DETA2
          DO 90 J=1,8
            AB(I,J-6)  = ETA
            IF ( HEM(I,J).GT.1.D-35 ) X(I,J-6) = LOG10(FAC*HEM(I,J))
90          CONTINUE
100       CONTINUE
        WRITE(6,1100)
1100    FORMAT('  JET CROSS SECTION  RAP.-DISTRIBUTION',/,
     &         '                     FOR THE DIFF. SUBPROCESSES',/)
        CALL PLOT(AB,X,400,8,50,ETA20,DETA2,XSMIN,XSSTEP)
      ENDIF
      WRITE(6,1110)
1110  FORMAT(/)
      RETURN
      END

C______________________________________________________________________
      SUBROUTINE HISINI
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( MAXPRO = 8 , MLINE = 1000 , MSCAHD = 250 )
C  LENGTH OF HISTO : 15000  REAL*4
      COMMON /HISTO / PT10,DPT1,ETA10,DETA1,PT20,DPT2,ETA20,DETA2,
     &                X(50,-5:5),AB(50,-5:5),HPE(50,-5:5),HEP(50,5),
     &                HPM(50,8),HEM(50,8),HP(50),HE(50),
     &                FILL(12192)
C
C  INITIALIZE HISTOGRAMS
      DPT1    = 1.
      DETA1   = 1.
      DPT2    = 2.
      DETA2   = 0.2
      PT10    = 1.
      ETA10   =-5.*DETA1
      PT20    = 2.
      ETA20   =-25.*DETA2
      DO 40 I=1,50
        HP(I)      = 0.0
        HE(I)      = 0.0
        DO 10 J=-5,5
10        HPE(I,J) = 0.0
        DO 20 J=1,5
20        HEP(I,J) = 0.0
        DO 30 J=1,8
          HPM(I,J) = 0.0
30        HEM(I,J) = 0.0
40      CONTINUE
      RETURN
      END
C______________________________________________________________________
      SUBROUTINE HISFIL
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( MAXPRO = 8 , MLINE = 1000 , MSCAHD = 250 )
      COMMON /HAEVTR/ LINE,LIN,LREC1(MLINE),LREC2(MLINE),PREC(0:3,MLINE)
      COMMON /HARSLT/ LSCAHD,LSC1HD,
     &                ETAHD(MSCAHD,2) ,PTHD(MSCAHD),
     &                XHD(MSCAHD,2)   ,VHD(MSCAHD) ,X0HD(MSCAHD,2),
     &                NINHD(MSCAHD,2) ,NOUTHD(MSCAHD,2),
     &                N0INHD(MSCAHD,2),NBRAHD(MSCAHD,2),NPROHD(MSCAHD)
C  LENGTH OF HISTO : 15000  REAL*4
      COMMON /HISTO / PT10,DPT1,ETA10,DETA1,PT20,DPT2,ETA20,DETA2,
     &                X(50,-5:5),AB(50,-5:5),HPE(50,-5:5),HEP(50,5),
     &                HPM(50,8),HEM(50,8),HP(50),HE(50),
     &                FILL(12192)
C
C  fill histogram
      DO 20 N=1,LSCAHD
        MSPR   = NPROHD(N)
        DO 10 K=1,2
          IPT1   = INT((PTHD(N)-PT10)/DPT1)+1
          IETA1  = INT((ETAHD(N,K)-ETA10)/DETA1+0.5)-5
          IPT2   = INT((PTHD(N)-PT20)/DPT2)+1
          IETA2  = INT((ETAHD(N,K)-ETA20)/DETA2+0.5)
          IF ( IPT1.GE. 1 .AND. IPT1.LE.50 ) THEN
            HPM(IPT1,MSPR)  = HPM(IPT1,MSPR)+1.
            HP(IPT1)        = HP(IPT1)+1.
            IF ( ABS(IETA1).LE.5 ) HPE(IPT1,IETA1) = HPE(IPT1,IETA1)+1.
          ENDIF
          IF ( IETA2.GE. 1 .AND. IETA2.LE.50 ) THEN
            HEM(IETA2,MSPR) = HEM(IETA2,MSPR)+1.
            HE(IETA2)       = HE(IETA2)+1.
            IF ( IPT2.GE.1 .AND. IPT2.LE.5 ) HEP(IETA2,IPT2) =
     &                                              HEP(IETA2,IPT2)+1.
          ENDIF
10      CONTINUE
20    CONTINUE
      RETURN
      END
C_______________________________________________________________________
      SUBROUTINE HATEST(IOUT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( MAXPRO = 8 , MLINE = 1000 , MSCAHD = 250 )
      COMMON /HAPARA/ ECM,PTINI(4),Q0SQR,ALASQR,BQCD,NPD,NF,NHA,NHB
      CHARACTER*18 PROC
      CHARACTER*11 PDSET,PARTIC
      COMMON /PEPROC/ PROC(0:MAXPRO),PDSET(23),PARTIC(-1:1)
C  LENGTH OF HISTO : 15000  REAL*4
      COMMON /HISTO / VVV(50),XS(50,6),AB(50,6),DSIG(0:MAXPRO),PD(-6:6),
     &                FILL(14328)
      IF ( IOUT.EQ.1 ) THEN
        CALL CSHARM(DSIG)
        WRITE(6,1010) ECM,PTINI(1),(PROC(M),DSIG(M),M=0,MAXPRO)
1010    FORMAT('  HARD CROSS SECTIONS FOR SINGLE PROCESSES',/,
     &   '  AT CM-ENERGY=',E8.1,' AND PTMIN=',F5.1,/,9(A25,E14.6,/))
      ELSEIF ( IOUT.EQ.2 ) THEN
C  PLOT PARTON DISTRIBUTIONS
        PDMIN  = 0.0
        PDSTEP = 0.04
        YMAX   = 5.0
        DY     = YMAX/50.
        QQ     = 10.0
        DO 10 J=1,50
          Y         = FLOAT(J-1)*DY
          VVV(J)    = 10.0**(-Y)
          DO 10 I=1,5
            XS(J,I) =-1.E+30
            AB(J,I) = Y
10          CONTINUE
        QQ   = 1.0
        DO 20 I=1,5
          QQ = QQ*10.
          DO 20 J=1,50
            CALL JTPDIS(VVV(J),QQ,1,1,PD)
            IF ( PD(0).GT.1.D-30 ) XS(J,I) = LOG10(PD(0))
20          CONTINUE
        WRITE(6,1020)
1020    FORMAT('   GLUONDISTRIBUTION OVER LOG10(X)  ( Q**2=10**I;',
     &                '   I=1...5 )')
        CALL PLOT(AB,XS,250,5,50,YMAX,-DY,PDMIN,PDSTEP)
      ELSEIF ( IOUT.EQ.3 ) THEN
        QQMIN  = 1.0
        QQSTEP = 0.1
        DO 30 I=1,50
          B         = FLOAT(I-1)*QQSTEP+QQMIN
          VVV(I)    = 10.0**B
          DO 30 J=1,5
            XS(I,J) = -1.E+30
            AB(I,J) = B
30        CONTINUE
        X   = 1.0
        DO 40 I=1,50
          X = 1.0
          DO 40 J=1,4
            X = X*0.1
            CALL JTPDIS(X,VVV(I),1,1,PD)
            IF ( PD(0).GT.1.D-30 ) XS(I,J) = LOG10(PD(0))
40          CONTINUE
        WRITE(6,1030)
1030    FORMAT('   GLUONDISTRIBUTION OVER LOG10(Q**2)  ( X=10**(-I)'
     &              ,';  I=1...4')
        CALL PLOT(AB,XS,200,4,50,QQMIN,QQSTEP,PDMIN,PDSTEP)
      ELSEIF ( IOUT.EQ.4 ) THEN
C  PLOT DIFFERENTIAL ONE JET CROSS SECTION
        XSMIN  =-6.
        XSSTEP = 0.1
        PTMIN  = 1.0
        PTSTEP = 1.0
        DO 50 I=1,50
          PT        = (I-1)*PTSTEP+PTMIN
          XS(I,1)   =-35.0
          DO 50 J=1,6
50          AB(I,J) = PT
C       DO 60 J=1,6
C         ETAC  = (J-1)*0.5
          ETAC = 0.0
          DO 60 I=1,50
            PT = AB(I,1)
            CALL CSJ1M(PT,ETAC,DSIG)
            IF ( DSIG(0).GT.1.D-30 ) XS(I,1) = LOG10(DSIG(0))
60          CONTINUE
        WRITE(6,1040)
1040    FORMAT('  DIFFERENTIAL HARD CROSS SECTION OVER PT , RAP.=0.')
        CALL PLOT(AB,XS,50,1,50,PTMIN,PTSTEP,XSMIN,XSSTEP)
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
C
C                                JTINIT.FOR
C
C---------------------------------------------------------------------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     PROCEDURES TO INITIALIZE PROGRAM
C      ( NONE OF THIS PROCEDURES IS USED DURING EVENT SIMULATION,
C        BEFORE DEMANDING EVENT SIMULATION SUBROUTINES
C        HASTRT, HARINI AND HISINI MUST BE CALLED )
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C_______________________________________________________________________
      SUBROUTINE HARINI
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( MAXPRO = 8 , MLINE = 1000 , MSCAHD = 250 )
      COMMON /HACONS/ PI,PI2,PI4,GEVTMB
      COMMON /HAPARA/ ECM,PTINI(4),Q0SQR,ALASQR,BQCD,NPD,NF,NHA,NHB
      COMMON /HAPDCO/ NPDCOR
      COMMON /HAQQAP/ AQQAL,AQQPD,NQQAL,NQQPD
      COMMON /HAOUTL/ NOUTL,NOUTER,NOUTCO
      INTEGER MXSECT
      COMMON /HAXSEC/ XSECTA(2,-1:MAXPRO,4),XSECT(5,-1:MAXPRO),
     &                MXSECT(0:2,-1:MAXPRO)
      CHARACTER*18 PROC
      CHARACTER*11 PDSET,PARTIC
      COMMON /PEPROC/ PROC(0:MAXPRO),PDSET(23),PARTIC(-1:1)
      DIMENSION DSIG(0:MAXPRO),ALAM(23),Q0S(23)
      DATA ALAM / 0.20, 0.29, 0.107, 0.250, 0.178, 0.25,
     *            0.10, 0.19, 0.190, 0.190, 0.190, 0.19,
     *            0.215,0.215,0.215,
     *            0.231,0.231,0.322, 0.247, 0.168,0.2,0.2,0.202 /
      DATA Q0S  / 5.0 , 5.0 , 5.0  , 5.0  , 5.0 , 0.2,
     *            5.0 , 5.0 , 5.0  , 5.0  , 5.0 , 5.0,
     *            5.0 , 5.0 , 5.0  , 4.0  , 4.0 , 4.0,
     *            4.0 , 4.0 ,  0.4 ,0.4 ,1.60     /
C
        IF ( NOUTL.GE.1 )CALL TIMDAT
        ALASQR = ALAM(NPD)**2
        Q0SQR  = Q0S(NPD)
        BQCD   = PI4/(11.-(2./3.)*NF)
C  check and sort the min. pt values in PTINI(1..4)
        INI = 0
        DO 30 I=1,4
          IF ( PTINI(I).LE..5D0.OR.PTINI(I).GE.ECM*.5D0)PTINI(I)=1.D+30
          IF ( PTINI(I).NE.1.D+30 ) INI = INI+1
30      CONTINUE
        DO 50 I=1,3
          DO 40 J=I+1,4
            IF ( PTINI(J).LT.PTINI(I) ) THEN
              TTT      = PTINI(J)
              PTINI(J) = PTINI(I)
              PTINI(I) = TTT
            ENDIF
40        CONTINUE
50      CONTINUE
C  calculate constant weights
        DO 10 M=-1,MAXPRO
          XSECT (3,M) = 0.0
          XSECT (4,M) = 0.0
          MXSECT(1,M) = 0
          MXSECT(2,M) = 0
10      CONTINUE
        DO 20     I = 1,4
          DO 20   M =-1,MAXPRO
            DO 20 J = 1,2
              XSECTA(J,M,I) = 0.0
20      CONTINUE
        DO 70 I=INI,1,-1
          CALL HABINT(I)
          IF ( NOUTL.GE.10 ) WRITE(6,1060) PTINI(I)
1060      FORMAT(' NORMALIZATION FOR PTMIN=',F10.4,' CALCULATED')
          CALL HAMAXI(I)
          IF ( NOUTL.GE.10 ) WRITE(6,1070) PTINI(I)
1070      FORMAT(' MAXIMA FOR PTMIN=',F10.4,' CALCULATED')
          XSECTA(1,0,I)   = PTINI(I)
          DO 60 M=-1,MAXPRO
            XSECTA(1,M,I) = XSECT(1,M)
            XSECTA(2,M,I) = XSECT(2,M)
60        CONTINUE
70      CONTINUE
C  calculate inclusive cross-sections
        CALL CSHARM(DSIG)
        DO 80 M=0,MAXPRO
          XSECT(5,M) = DSIG(M)
80      CONTINUE
C
C  print results of initialization
C
        IF ( NOUTL.GE.10 ) WRITE(6,'(/,1X,70(1H*))')
        WRITE(6,1057) PTINI(1),PDSET(NPD),SQRT(ALASQR),Q0SQR
1057    FORMAT(/,
     &         ' --- parameters of the hard scattering program ---',/,
     &         '       MIN. PT       :',F15.1,/,
     &         '       PARTON-DISTR. :',A15,/,
     &         '       LAMBDA        :',F15.3,/,
     &         '       Q0**2         :',F15.3,/)
      IF ( NOUTL.GE.1 ) THEN
        WRITE(6,1050) PARTIC(NHA),PARTIC(NHB),ECM,PTINI(1),PDSET(NPD),
     &                SQRT(ALASQR),Q0SQR,NPDCOR,NQQAL,AQQAL,NQQPD,AQQPD
1050    FORMAT(/,1X,70('*'),/,
     &         '  HARD SCATTERING PROGRAM IS INITIALIZED FOR',/,
     &         '    PROJECTILE    :',A15,/,
     &         '    TARGET        :',A15,/,
     &         '    CM-ENERGY     :',F15.1,/,
     &         '    MIN. PT       :',F15.1,/,
     &         '    PARTON-DISTR. :',A15,/,
     &         '    LAMBDA        :',F15.3,/,
     &         '    Q0**2         :',F15.3,/,
     &         '    NPDCOR        :',I15,/,
     &         '    NQQAL         :',I15,/,
     &         '    AQQAL         :',F15.3,/,
     &         '    NQQPD         :',I15,/,
     &         '    AQQPD         :',F15.3,/)
         CALL TIMDAT
      ENDIF
      RETURN
      END
C_______________________________________________________________________
      SUBROUTINE HABINT(IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( MAXPRO = 8 , MLINE = 1000 , MSCAHD = 250 )
      PARAMETER ( MXABWT = 1000 )
      PARAMETER ( ZERO=0.D0, ONE=1.D0)
      COMMON /HACONS/ PI,PI2,PI4,GEVTMB
      COMMON /HAPARA/ ECM,PTINI(4),Q0SQR,ALASQR,BQCD,NPD,NF,NHA,NHB
      COMMON /HAGAUP/ NGAUP1,NGAUP2,NGAUET,NGAUIN
      INTEGER MXSECT
      COMMON /HAXSEC/ XSECTA(2,-1:MAXPRO,4),XSECT(5,-1:MAXPRO),
     &                MXSECT(0:2,-1:MAXPRO)
      DIMENSION       ABSZ(MXABWT),WEIG(MXABWT)
      DIMENSION S(-1:MAXPRO),S1(-1:MAXPRO),S2(-1:MAXPRO),F124(-1:MAXPRO)
      DATA F124 / 1.,0.,4.,2.,2.,2.,4.,1.,4.,4. /

      A      = (2.*PTINI(IND)/ECM)**2
      ALN    = LOG(A)
      HLN    = LOG(0.5)
      NPOINT = NGAUIN
      CALL GSET(ZERO,ONE,NPOINT,ABSZ,WEIG)
      DO 10 M=-1,MAXPRO
        S1(M) = 0.0
10    CONTINUE
      DO 80 I1=1,NPOINT
        Z1   = ABSZ(I1)
        X1   = EXP(ALN*Z1)
        DO 20 M=-1,MAXPRO
          S2(M) = 0.0
20      CONTINUE
        DO 60 I2=1,NPOINT
          Z2 = (1.-Z1)*ABSZ(I2)
          X2 = EXP(ALN*Z2)
          FAXX = A/(X1*X2)
          W    = SQRT(1.-FAXX)
          W1   = FAXX/(1.+W)
          WLOG = LOG(W1)
          FWW  = FAXX*WLOG/W
          DO 30 M=-1,MAXPRO
            S(M) = 0.0
30        CONTINUE
          DO 40 I=1,NPOINT
            Z   = ABSZ(I)
            VA  =-0.5*W1/(W1+Z*W)
            UA  =-1.-VA
            VB  =-0.5*FAXX/(W1+2.*W*Z)
            UB  =-1.-VB
            VC  =-EXP(HLN+Z*WLOG)
            UC  =-1.-VC
            VE  =-0.5*(1.+W)+Z*W
            UE  =-1.-VE
            S(1)  = S(1)+(1.+W)*2.25*(VA*VA*(3.-UA*VA-VA/(UA*UA))-UA)*
     &           WEIG(I)
            S(2)  = S(2)+(VC*VC+UC*UC)*((16./27.)/UC-(4./3.)*VC)*FWW*
     &            WEIG(I)
            S(3)  = S(3)+(1.+W)*(1.+UA*UA)*(1.-(4./9.)*VA*VA/UA)*WEIG(I)
            S(5)  = S(5)+((4./9.)*(1.+UB*UB+(UB*UB+VB*VB)*VB*VB)-
     &            (8./27.)*UA*UA*VA)*WEIG(I)
            S(6)  = S(6)+(4./9.)*(UE*UE+VE*VE)*FAXX*WEIG(I)
            S(7)  = S(7)+(1.+W)*((2./9.)*(1.+UA*UA+(1.+VA*VA)*VA*VA/
     &            (UA*UA))-(4./27.)*VA/UA)*WEIG(I)
            S(8)  = S(8)+(4./9.)*(1.+UB*UB)*WEIG(I)
            S(-1) = S(-1)+(1.+VC*VC)*(VC/(UC*UC)-(4./9.))*FWW*WEIG(I)
40        CONTINUE
          S(4)    = S(2)*(9./32.)
          DO 50 M=-1,MAXPRO
            S2(M) = S2(M)+S(M)*WEIG(I2)*W
50        CONTINUE
60      CONTINUE
        DO 70 M=-1,MAXPRO
          S1(M) = S1(M)+S2(M)*(1.-Z1)*WEIG(I1)
70      CONTINUE
80    CONTINUE
      FFF    = PI*GEVTMB*ALN*ALN/(A*ECM*ECM)
      DO 90 M=-1,MAXPRO
        XSECT(1,M) = FFF*F124(M)*S1(M)
90    CONTINUE
      XSECT(1,4) = XSECT(1,4)*NF
      XSECT(1,6) = XSECT(1,6)*MAX(0,NF-1)
      RETURN
      END
C______________________________________________________________________
      SUBROUTINE HAMAXI(IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( MAXPRO = 8 , MLINE = 1000 , MSCAHD = 250 )
      PARAMETER ( NKM = 5 )
      PARAMETER ( TINY= 1.D-30 )
      COMMON /HAPARA/ ECM,PTINI(4),Q0SQR,ALASQR,BQCD,NPD,NF,NHA,NHB
      INTEGER MXSECT
      COMMON /HAXSEC/ XSECTA(2,-1:MAXPRO,4),XSECT(5,-1:MAXPRO),
     &                MXSECT(0:2,-1:MAXPRO)
      DIMENSION Z(3),D(3),FF(NKM)

      DO 40 NKON=1,NKM
        Z(1) = 0.99
        Z(2) = 0.5
        Z(3) = 0.0
        D(1) =-1.0
        D(2) = 2.0
        D(3) = 2.5
        IT   = 0
        CALL HAFDI1(NKON,Z,F2,IND)
10        IT   = IT+1
          FOLD = F2
          DO 30 I=1,3
            D(I) = D(I)/5.
            Z(I)   = Z(I)+D(I)
            CALL HAFDI1(NKON,Z,F3,IND)
            IF ( F2.GT.F3 ) Z(I) = Z(I)-D(I)
            IF ( F2.GT.F3 ) D(I) =-D(I)
20            F1   = MIN(F2,F3)
              F2   = MAX(F2,F3)
              Z(I) = Z(I)+D(I)
              CALL HAFDI1(NKON,Z,F3,IND)
              IF ( F3.GT.F2 ) GOTO 20
            ZZ     = Z(I)-D(I)
            Z(I)   = ZZ+0.5*D(I)*(F3-F1)/MAX(TINY,F2+F2-F1-F3)
            IF ( ABS(ZZ-Z(I)).GT.D(I)*0.1D0)CALL HAFDI1(NKON,Z,F1,IND)
            IF ( F1.LE.F2 ) Z(I) = ZZ
            F2     = MAX(F1,F2)
30        CONTINUE
          IF ( ABS(FOLD-F2)/F2.GT.0.002D0.OR. IT.LT.3 ) GOTO 10
        FF(NKON) = F2
40    CONTINUE
      XSECT(2,1) = FF(1)*XSECT(1,1)
      XSECT(2,2) = FF(2)*XSECT(1,2)
      XSECT(2,3) = FF(4)*XSECT(1,3)
      XSECT(2,4) = FF(1)*XSECT(1,4)
      XSECT(2,5) = FF(2)*XSECT(1,5)
      XSECT(2,6) = FF(2)*XSECT(1,6)
      XSECT(2,7) = FF(3)*XSECT(1,7)
      XSECT(2,8) = FF(5)*XSECT(1,8)
      XSECT(2,-1)= FF(4)*XSECT(1,-1)
      RETURN
      END
C______________________________________________________________________
      SUBROUTINE HAFDI1(NKON,Z,FDIS,IND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( MAXPRO = 8 , MLINE = 1000 , MSCAHD = 250 )
      PARAMETER ( NKM = 5 )
      PARAMETER ( TINY= 1.D-30, ONE=1.D0 ,TINY6=1.D-06,ZERO=0.D0)
      COMMON /HAPARA/ ECM,PTINI(4),Q0SQR,ALASQR,BQCD,NPD,NF,NHA,NHB
      COMMON /HAQQAP/ AQQAL,AQQPD,NQQAL,NQQPD
      DIMENSION F(NKM),PDA(-6:6),PDB(-6:6),Z(3)

      FDIS = 0.0
C check input values
      IF ( Z(1).LE.0.0D0  .OR.  Z(1).GE.1.0D0 ) RETURN
      IF ( Z(2).LE.0.0D0  .OR.  Z(2).GE.1.0D0 ) RETURN
      IF ( Z(3).LT.0.0D0  .OR.  Z(3).GT.1.0D0 ) RETURN
      A     = (2.*PTINI(IND)/ECM)**2
      ALN   = LOG(A)
      Y1  = EXP(ALN*Z(1))
      Y2  =-(1.-Y1)+2.*(1.-Y1)*Z(2)
      X1  = 0.5*(Y2+SQRT(Y2*Y2+4.*Y1))
      X2  = X1-Y2
      W   = SQRT(MAX(TINY,1.-A/Y1))
      V   =-0.5+W*(Z(3)-0.5)
      U   =-(1.+V)
      PT  = MAX(PTINI(IND),SQRT(U*V*Y1*ECM*ECM))
C  set hard scale  QQ  for alpha and partondistr.
      IF     ( NQQAL.EQ.1 ) THEN
        QQAL = AQQAL*PT*PT
      ELSEIF ( NQQAL.EQ.2 ) THEN
        QQAL = AQQAL*Y1*ECM*ECM
      ELSEIF ( NQQAL.EQ.3 ) THEN
        QQAL = AQQAL*Y1*ECM*ECM*(U*V)**(1./3.)
      ELSEIF ( NQQAL.EQ.4 ) THEN
        QQAL = AQQAL*Y1*ECM*ECM*U*V/(1.+V*V+U*U)
      ENDIF
      IF     ( NQQPD.EQ.1 ) THEN
        QQPD = AQQPD*PT*PT
      ELSEIF ( NQQPD.EQ.2 ) THEN
        QQPD = AQQPD*Y1*ECM*ECM
      ELSEIF ( NQQPD.EQ.3 ) THEN
        QQPD = AQQPD*Y1*ECM*ECM*(U*V)**(1./3.)
      ELSEIF ( NQQPD.EQ.4 ) THEN
        QQPD = AQQPD*Y1*ECM*ECM*U*V/(1.+V*V+U*U)
      ENDIF
      FACTOR = (BQCD/LOG(MAX(QQAL/ALASQR,1.1*ONE)))**2
C calculate partondistributions
      CALL JTPDIS(X1,QQPD,NHA,0,PDA)
      CALL JTPDIS(X2,QQPD,NHB,0,PDB)
C calculate full distribution FDIS
      DO 10 N=1,NKM
        F(N) = 0.0
10    CONTINUE
      DO 20 I=1,NF
        F(2) = F(2)+PDA(I)*PDB(-I)+PDA(-I)*PDB( I)
        F(3) = F(3)+PDA(I)*PDB( I)+PDA(-I)*PDB(-I)
        F(4) = F(4)+PDA(I)+PDA(-I)
        F(5) = F(5)+PDB(I)+PDB(-I)
20    CONTINUE
      F(1)   = PDA(0)*PDB(0)
      T      = PDA(0)*F(5)+PDB(0)*F(4)
      F(5)   = F(4)*F(5)-(F(2)+F(3))
      F(4)   = T
      FDIS   = MAX(ZERO,F(NKON)*FACTOR)
      RETURN
      END
C______________________________________________________________________
      SUBROUTINE HASTRT
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( MAXPRO = 8 , MLINE = 1000 , MSCAHD = 250 )
      COMMON /HACONS/ PI,PI2,PI4,GEVTMB
      COMMON /HAPARA/ ECM,PTINI(4),Q0SQR,ALASQR,BQCD,NPD,NF,NHA,NHB
      COMMON /HAPADI/ NPDM
      COMMON /HAPDCO/ NPDCOR
      COMMON /HAQQAP/ AQQAL,AQQPD,NQQAL,NQQPD
      COMMON /HAOUTL/ NOUTL,NOUTER,NOUTCO
      COMMON /HACUTS/ PTL,PTU,ETACL,ETACU,ETADL,ETADU
      COMMON /HAGAUP/ NGAUP1,NGAUP2,NGAUET,NGAUIN
      COMMON /HAEVTR/ LINE,LIN,LREC1(MLINE),LREC2(MLINE),PREC(0:3,MLINE)
      INTEGER MXSECT
      COMMON /HAXSEC/ XSECTA(2,-1:MAXPRO,4),XSECT(5,-1:MAXPRO),
     &                MXSECT(0:2,-1:MAXPRO)
      COMMON /HAXSUM/XSHMX
C  initialize COMMON /HAXSUM/
      XSHMX    = 1.
C  initialize COMMON /HACONS/
      PI       = 3.1415927
      PI2      = 2.*PI
      PI4      = 4.*PI
      GEVTMB   = 0.389365
C  initialize COMMON /HAPARA/
      ECM      = 1800.
      PTINI(1) = 2.0
      PTINI(2) = 0.0
      PTINI(3) = 0.0
      PTINI(4) = 0.0
      Q0SQR    = 5.0
      ALASQR   = 0.107**2
      NF       = 4
      BQCD     = PI4/(11.0-(2./3.)*NF)
      NPD      = 3
      NHA      = 1
      NHB      =-1
C  initialize COMMON /HAPADI/
      NPDM     = 0
C  initialize COMMON /HAPDCO/
      NPDCOR   = 0
C  initialize COMMON /HAQQAP/
      NQQAL    = 1
      AQQAL    = 0.25
      NQQPD    = 1
      AQQPD    = 0.25
C  initialize COMMON /HAOUTL/
      NOUTL    = 1
      NOUTER   = 1
      NOUTCO   = 0
C  initialize COMMON /HAGAUP/
      NGAUP1   = 8
      NGAUP2   = 8
      NGAUET   = 8
      NGAUIN   = 8
C  initialize COMMON /HACUTS/
      PTL      = 0.E+00
      PTU      = 1.E+30
      ETACL    =-1.E+30
      ETACU    = 1.E+30
      ETADL    =-1.E+30
      ETADU    = 1.E+30
C  initialize COMMON /HAEVTR/
      LINE       = 0
      LIN        = 0
      DO 20 L = 1,MLINE
        LREC1(L) = 0
        LREC2(L) = 0
        DO 10 I=0,3
          PREC(I,L) = 0.0
10      CONTINUE
20    CONTINUE
C  initialize COMMON /HAXSEC/
      DO 40 M=-1,MAXPRO
        DO 30 I=1,5
          XSECT(I,M) = 0.0
30      CONTINUE
        MXSECT(1,M)  = 0
        MXSECT(2,M)  = 0
        MXSECT(0,M)  = 1
40    CONTINUE
      RETURN
      END
C----------------------------------------------------------------------
      BLOCK DATA JTDATA
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( MAXPRO = 8 , MLINE = 1000 , MSCAHD = 250 )
      CHARACTER*18 PROC
      CHARACTER*11 PDSET,PARTIC
      COMMON /PEPROC/ PROC(0:MAXPRO),PDSET(23),PARTIC(-1:1)

      DATA PROC   /  'SUM OVER PROCESSES', 'G  +G  --> G  +G  ',
     &               'Q  +QB --> G  +G  ', 'G  +Q  --> G  +Q  ',
     &               'G  +G  --> Q  +QB ', 'Q  +QB --> Q  +QB ',
     &               'Q  +QB --> QS +QBS', 'Q  +Q  --> Q  +Q  ',
     &               'Q  +QS --> Q  +QS '                       /
      DATA PDSET  / ' EHLQ SET 1',' EHLQ SET 2',' MRS  SET 1',
     &              ' MRS  SET 2',' MRS  SET 3',' GRV LO    ',
     &              ' HMRS SET 1',' HMRS SET 2',' KMRS SET 1',
     &              ' KMRS SET 2',' KMRS SET 3',' KMRS SET 4',
     &              ' MRS(S0)   ',' MRS(D0)   ',' MRS(D-)   ',
     &              ' CTEQ 1M   ',' CTEQ 1MS  ',' CTEQ 1ML  ',
     &              ' CTEQ 1D   ',' CTEQ 1L   ',' GRV94LO1  ' ,
     &              ' GRV98LO   ',' CTEQ96    '/
      DATA PARTIC / ' ANTIPROTON','           ','     PROTON'   /
      END
C***********************************************************************
C
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
*    G R V  -  P R O T O N  - P A R A M E T R I Z A T I O N S     *
*                                                                 *
*                         1994 UPDATE                             *
*                                                                 *
*                 FOR A DETAILED EXPLANATION SEE                  *
*                   M. GLUECK, E.REYA, A.VOGT :                   *
*                   DO-TH 94/24  =  DESY 94-206                   *
*                    (TO APPEAR IN Z. PHYS. C)                    *
*                                                                 *
*   THE PARAMETRIZATIONS ARE FITTED TO THE EVOLVED PARTONS FOR    *
*        Q**2 / GEV**2  BETWEEN   0.4   AND  1.E6                 *
*             X         BETWEEN  1.E-5  AND   1.                  *
*   LARGE-X REGIONS, WHERE THE DISTRIBUTION UNDER CONSIDERATION   *
*   IS NEGLIGIBLY SMALL, WERE EXCLUDED FROM THE FIT.              *
*                                                                 *
*   HEAVY QUARK THRESHOLDS  Q(H) = M(H)  IN THE BETA FUNCTION :   *
*                   M(C)  =  1.5,  M(B)  =  4.5                   *
*   CORRESPONDING LAMBDA(F) VALUES IN GEV FOR  Q**2 > M(H)**2 :   *
*      LO :   LAMBDA(3)  =  0.232,   LAMBDA(4)  =  0.200,         *
*             LAMBDA(5)  =  0.153,                                *
*      NLO :  LAMBDA(3)  =  0.248,   LAMBDA(4)  =  0.200,         *
*             LAMBDA(5)  =  0.131.                                *
*   THE NUMBER OF ACTIVE QUARK FLAVOURS IS  NF = 3  EVERYWHERE    *
*   EXCEPT IN THE BETA FUNCTION, I.E. THE HEAVY QUARKS C,B,...    *
*   ARE NOT PRESENT AS PARTONS IN THE Q2-EVOLUTION.               *
*   IF NEEDED, HEAVY QUARK DENSITIES CAN BE TAKEN FROM THE 1991   *
*   GRV PARAMETRIZATION.                                          *
*                                                                 *
*   NLO DISTRIBUTIONS ARE GIVEN IN MS-BAR FACTORIZATION SCHEME    *
*   (SUBROUTINE GRV94HO) AS WELL AS IN THE DIS SCHEME (GRV94DI),  *
*   THE LEADING ORDER PARAMETRIZATION IS PROVIDED BY "GRV94LO".   *
*                                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*...INPUT PARAMETERS :
*
*    X   = MOMENTUM FRACTION
*    Q2  = SCALE Q**2 IN GEV**2
*
*...OUTPUT (ALWAYS X TIMES THE DISTRIBUTION) :
*
*    UV  = U(VAL) = U - U(BAR)
*    DV  = D(VAL) = D - D(BAR)
*    DEL = D(BAR) - U(BAR)
*    UDB = U(BAR) + D(BAR)
*    SB  = S = S(BAR)
*    GL  = GLUON
*
*...LO PARAMETRIZATION :
*
       SUBROUTINE DOR94LO (X, Q2, UV, DV, DEL, UDB, SB, GL)
       IMPLICIT DOUBLE PRECISION (A - Z)
      SAVE
       MU2  = 0.23
       LAM2 = 0.2322 * 0.2322
       S  = LOG (LOG(Q2/LAM2) / LOG(MU2/LAM2))
       DS = SQRT (S)
       S2 = S * S
       S3 = S2 * S
*...UV :
       NU  =  2.284 + 0.802 * S + 0.055 * S2
       AKU =  0.590 - 0.024 * S
       BKU =  0.131 + 0.063 * S
       AU  = -0.449 - 0.138 * S - 0.076 * S2
       BU  =  0.213 + 2.669 * S - 0.728 * S2
       CU  =  8.854 - 9.135 * S + 1.979 * S2
       DU  =  2.997 + 0.753 * S - 0.076 * S2
       UV  = DOR94FV (X, NU, AKU, BKU, AU, BU, CU, DU)
*...DV :
       ND  =  0.371 + 0.083 * S + 0.039 * S2
       AKD =  0.376
       BKD =  0.486 + 0.062 * S
       AD  = -0.509 + 3.310 * S - 1.248 * S2
       BD  =  12.41 - 10.52 * S + 2.267 * S2
       CD  =  6.373 - 6.208 * S + 1.418 * S2
       DD  =  3.691 + 0.799 * S - 0.071 * S2
       DV  = DOR94FV (X, ND, AKD, BKD, AD, BD, CD, DD)
*...DEL :
       NE  =  0.082 + 0.014 * S + 0.008 * S2
       AKE =  0.409 - 0.005 * S
       BKE =  0.799 + 0.071 * S
       AE  = -38.07 + 36.13 * S - 0.656 * S2
       BE  =  90.31 - 74.15 * S + 7.645 * S2
       CE  =  0.0
       DE  =  7.486 + 1.217 * S - 0.159 * S2
       DEL = DOR94FV (X, NE, AKE, BKE, AE, BE, CE, DE)
*...UDB :
       ALX =  1.451
       BEX =  0.271
       AKX =  0.410 - 0.232 * S
       BKX =  0.534 - 0.457 * S
       AGX =  0.890 - 0.140 * S
       BGX = -0.981
       CX  =  0.320 + 0.683 * S
       DX  =  4.752 + 1.164 * S + 0.286 * S2
       EX  =  4.119 + 1.713 * S
       ESX =  0.682 + 2.978 * S
       UDB=DOR94FW(X, S, ALX, BEX, AKX, BKX, AGX, BGX, CX, DX, EX, ESX)
*...SB :
       ALS =  0.914
       BES =  0.577
       AKS =  1.798 - 0.596 * S
       AS  = -5.548 + 3.669 * DS - 0.616 * S
       BS  =  18.92 - 16.73 * DS + 5.168 * S
       DST =  6.379 - 0.350 * S  + 0.142 * S2
       EST =  3.981 + 1.638 * S
       ESS =  6.402
       SB  = DOR94FS (X, S, ALS, BES, AKS, AS, BS, DST, EST, ESS)
*...GL :
       ALG =  0.524
       BEG =  1.088
       AKG =  1.742 - 0.930 * S
       BKG =        - 0.399 * S2
       AG  =  7.486 - 2.185 * S
       BG  =  16.69 - 22.74 * S  + 5.779 * S2
       CG  = -25.59 + 29.71 * S  - 7.296 * S2
       DG  =  2.792 + 2.215 * S  + 0.422 * S2 - 0.104 * S3
       EG  =  0.807 + 2.005 * S
       ESG =  3.841 + 0.316 * S
       GL =DOR94FW(X, S, ALG, BEG, AKG, BKG, AG, BG, CG, DG, EG, ESG)
       RETURN
       END
*
*...NLO PARAMETRIZATION (MS(BAR)) :
*
       SUBROUTINE DOR94HO (X, Q2, UV, DV, DEL, UDB, SB, GL)
       IMPLICIT DOUBLE PRECISION (A - Z)
      SAVE
       MU2  = 0.34
       LAM2 = 0.248 * 0.248
       S  = LOG (LOG(Q2/LAM2) / LOG(MU2/LAM2))
       DS = SQRT (S)
       S2 = S * S
       S3 = S2 * S
*...UV :
       NU  =  1.304 + 0.863 * S
       AKU =  0.558 - 0.020 * S
       BKU =          0.183 * S
       AU  = -0.113 + 0.283 * S - 0.321 * S2
       BU  =  6.843 - 5.089 * S + 2.647 * S2 - 0.527 * S3
       CU  =  7.771 - 10.09 * S + 2.630 * S2
       DU  =  3.315 + 1.145 * S - 0.583 * S2 + 0.154 * S3
       UV  = DOR94FV (X, NU, AKU, BKU, AU, BU, CU, DU)
*...DV :
       ND  =  0.102 - 0.017 * S + 0.005 * S2
       AKD =  0.270 - 0.019 * S
       BKD =  0.260
       AD  =  2.393 + 6.228 * S - 0.881 * S2
       BD  =  46.06 + 4.673 * S - 14.98 * S2 + 1.331 * S3
       CD  =  17.83 - 53.47 * S + 21.24 * S2
       DD  =  4.081 + 0.976 * S - 0.485 * S2 + 0.152 * S3
       DV  = DOR94FV (X, ND, AKD, BKD, AD, BD, CD, DD)
*...DEL :
       NE  =  0.070 + 0.042 * S - 0.011 * S2 + 0.004 * S3
       AKE =  0.409 - 0.007 * S
       BKE =  0.782 + 0.082 * S
       AE  = -29.65 + 26.49 * S + 5.429 * S2
       BE  =  90.20 - 74.97 * S + 4.526 * S2
       CE  =  0.0
       DE  =  8.122 + 2.120 * S - 1.088 * S2 + 0.231 * S3
       DEL = DOR94FV (X, NE, AKE, BKE, AE, BE, CE, DE)
*...UDB :
       ALX =  0.877
       BEX =  0.561
       AKX =  0.275
       BKX =  0.0
       AGX =  0.997
       BGX =  3.210 - 1.866 * S
       CX  =  7.300
       DX  =  9.010 + 0.896 * DS + 0.222 * S2
       EX  =  3.077 + 1.446 * S
       ESX =  3.173 - 2.445 * DS + 2.207 * S
       UDB=DOR94FW(X, S, ALX, BEX, AKX, BKX, AGX, BGX, CX, DX, EX, ESX)
*...SB :
       ALS =  0.756
       BES =  0.216
       AKS =  1.690 + 0.650 * DS - 0.922 * S
       AS  = -4.329 + 1.131 * S
       BS  =  9.568 - 1.744 * S
       DST =  9.377 + 1.088 * DS - 1.320 * S + 0.130 * S2
       EST =  3.031 + 1.639 * S
       ESS =  5.837 + 0.815 * S
       SB  = DOR94FS (X, S, ALS, BES, AKS, AS, BS, DST, EST, ESS)
*...GL :
       ALG =  1.014
       BEG =  1.738
       AKG =  1.724 + 0.157 * S
       BKG =  0.800 + 1.016 * S
       AG  =  7.517 - 2.547 * S
       BG  =  34.09 - 52.21 * DS + 17.47 * S
       CG  =  4.039 + 1.491 * S
       DG  =  3.404 + 0.830 * S
       EG  = -1.112 + 3.438 * S  - 0.302 * S2
       ESG =  3.256 - 0.436 * S
       GL =DOR94FW(X, S, ALG, BEG, AKG, BKG, AG, BG, CG, DG, EG, ESG)
       RETURN
       END
*
*...NLO PARAMETRIZATION (DIS) :
*
       SUBROUTINE DOR94DI (X, Q2, UV, DV, DEL, UDB, SB, GL)
       IMPLICIT DOUBLE PRECISION (A - Z)
      SAVE
       MU2  = 0.34
       LAM2 = 0.248 * 0.248
       S  = LOG (LOG(Q2/LAM2) / LOG(MU2/LAM2))
       DS = SQRT (S)
       S2 = S * S
       S3 = S2 * S
*...UV :
       NU  =  2.484 + 0.116 * S + 0.093 * S2
       AKU =  0.563 - 0.025 * S
       BKU =  0.054 + 0.154 * S
       AU  = -0.326 - 0.058 * S - 0.135 * S2
       BU  = -3.322 + 8.259 * S - 3.119 * S2 + 0.291 * S3
       CU  =  11.52 - 12.99 * S + 3.161 * S2
       DU  =  2.808 + 1.400 * S - 0.557 * S2 + 0.119 * S3
       UV  = DOR94FV (X, NU, AKU, BKU, AU, BU, CU, DU)
*...DV :
       ND  =  0.156 - 0.017 * S
       AKD =  0.299 - 0.022 * S
       BKD =  0.259 - 0.015 * S
       AD  =  3.445 + 1.278 * S + 0.326 * S2
       BD  = -6.934 + 37.45 * S - 18.95 * S2 + 1.463 * S3
       CD  =  55.45 - 69.92 * S + 20.78 * S2
       DD  =  3.577 + 1.441 * S - 0.683 * S2 + 0.179 * S3
       DV  = DOR94FV (X, ND, AKD, BKD, AD, BD, CD, DD)
*...DEL :
       NE  =  0.099 + 0.019 * S + 0.002 * S2
       AKE =  0.419 - 0.013 * S
       BKE =  1.064 - 0.038 * S
       AE  = -44.00 + 98.70 * S - 14.79 * S2
       BE  =  28.59 - 40.94 * S - 13.66 * S2 + 2.523 * S3
       CE  =  84.57 - 108.8 * S + 31.52 * S2
       DE  =  7.469 + 2.480 * S - 0.866 * S2
       DEL = DOR94FV (X, NE, AKE, BKE, AE, BE, CE, DE)
*...UDB :
       ALX =  1.215
       BEX =  0.466
       AKX =  0.326 + 0.150 * S
       BKX =  0.956 + 0.405 * S
       AGX =  0.272
       BGX =  3.794 - 2.359 * DS
       CX  =  2.014
       DX  =  7.941 + 0.534 * DS - 0.940 * S + 0.410 * S2
       EX  =  3.049 + 1.597 * S
       ESX =  4.396 - 4.594 * DS + 3.268 * S
       UDB=DOR94FW(X, S, ALX, BEX, AKX, BKX, AGX, BGX, CX, DX, EX, ESX)
*...SB :
       ALS =  0.175
       BES =  0.344
       AKS =  1.415 - 0.641 * DS
       AS  =  0.580 - 9.763 * DS + 6.795 * S  - 0.558 * S2
       BS  =  5.617 + 5.709 * DS - 3.972 * S
       DST =  13.78 - 9.581 * S  + 5.370 * S2 - 0.996 * S3
       EST =  4.546 + 0.372 * S2
       ESS =  5.053 - 1.070 * S  + 0.805 * S2
       SB  = DOR94FS (X, S, ALS, BES, AKS, AS, BS, DST, EST, ESS)
*...GL :
       ALG =  1.258
       BEG =  1.846
       AKG =  2.423
       BKG =  2.427 + 1.311 * S  - 0.153 * S2
       AG  =  25.09 - 7.935 * S
       BG  = -14.84 - 124.3 * DS + 72.18 * S
       CG  =  590.3 - 173.8 * S
       DG  =  5.196 + 1.857 * S
       EG  = -1.648 + 3.988 * S  - 0.432 * S2
       ESG =  3.232 - 0.542 * S
       GL = DOR94FW(X, S, ALG, BEG, AKG, BKG, AG, BG, CG, DG, EG, ESG)
       RETURN
       END
*
*...FUNCTIONAL FORMS OF THE PARAMETRIZATIONS :
*
       FUNCTION DOR94FV (X, N, AK, BK, A, B, C, D)
       IMPLICIT DOUBLE PRECISION (A - Z)
      SAVE
       DX = SQRT (X)
       DOR94FV=N * X**AK * (1.+ A*X**BK + X * (B + C*DX)) * (1.- X)**D
       RETURN
       END
*
       FUNCTION DOR94FW (X, S, AL, BE, AK, BK, A, B, C, D, E, ES)
       IMPLICIT DOUBLE PRECISION (A - Z)
      SAVE
       LX = LOG (1./X)
       DOR94FW = (X**AK * (A + X * (B + X*C)) * LX**BK + S**AL
     1      * DEXP (-E + SQRT (ES * S**BE * LX))) * (1.- X)**D
       RETURN
       END
*
       FUNCTION DOR94FS (X, S, AL, BE, AK, AG, B, D, E, ES)
       IMPLICIT DOUBLE PRECISION (A - Z)
      SAVE
       DX = SQRT (X)
       LX = LOG (1./X)
       DOR94FS = S**AL / LX**AK * (1.+ AG*DX + B*X) * (1.- X)**D
     1       * DEXP (-E + SQRT (ES * S**BE * LX))
       RETURN
       END
*

C---------------- end of file -----------------------------------------
