C---------------------------------------------------------------------
C                   dtcpomqj.f    (July 1993)
C-----------------------------------------------------------------------
************************************************************************
*
      SUBROUTINE SIGSHD(ECM)
*     May 1991
*     input:
*     hard cross sections (see DATA statements)
*     ECM (independent of CMENER in /USER/)
*     output:
*     bar cross sections for soft (SIGSOF) hard (SIGHAR) and diffractive
*     (SIGTRP) SCATTERING
*
*- - -  - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*
*     version determined by ISIG
C*********************************************************************
C    ISIG=1-9 dropped since dpmjet-II.4.2
C*********************************************************************
*     ISIG=1  Capella,Tran Thanh Van,Kwiecinski,PRL 58(1987)2015
*     ISIG=2  Capella,Tran Thanh Van,Kwiecinski,PRL 58(1987)2015 CHG.SIG
*     ISIG=3  ALL SUBPROCESSES HARD CROSS SECTIONS PTHR CHANGING WITH ECM
*     ISIG=4  ALL SUBPROCESSES HARD CROSS SECTIONS PTHR=3 GEV/C   MRS1
*     ISIG=5  ALL SUBPROCESSES HARD CROSS SECTIONS PTHR=2 GEV/C   MRS1
*     ISIG=6  ALL SUBPROCESSES HARD CROSS SECTIONS PTHR=1.3 GEV/C MRS1
*     ISIG=7  all subprocesses hard cross sections PTHR=2GEV/C    MRS1
*     ISIG=8  program written by Patrick and Maire
*     ISIG=9  ALL SUBPROCESSES HARD CROSS SECTIONS PTHR=1.0 GEV/C MRS1
C*********************************************************************
C    ISIG=1-9 dropped since dpmjet-II.4.2
C*********************************************************************
*     ISIG=10  Version ISIG=4 including different sets
*              of structure functions
*
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
*     CONVERSION FACTOR GEV**-2 TO MILLIBARNS
      PARAMETER (CONV=.38935D0)
      PARAMETER (PI=3.141592654D0,
     &           THREE=3.D0,
     &           TWO =2.D0,
     &           EPS=1.D-3)
      PARAMETER (THOUSA = 1000.D0)
*
* *** /OUTLEV/ controls output level for POMDI and parton X distribution
      COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
*
* *** /POMPAR/*/SIGMA/*/POMTYP/ used only in SIGMAPOM-routines (->POMDI)
      COMMON /POMTYP/IPIM,ICON,ISIG,LMAX,MMAX,NMAX,DIFEL,DIFNU
      COMMON/POMPAR/ALFA,ALFAP,A,C,AK
      COMMON /SIGMA/SIGSOF,BS,ZSOF,SIGHAR,BH,ZHAR,SIGTRP,BT,ZTRP,
     *              SIGLOO,ZLOO
*
      CHARACTER*80 TITLE
      CHARACTER*8 PROJTY,TARGTY
C     COMMON /USER/TITLE,PROJTY,TARGTY,CMENER,ISTRUF
C    &            ,ISINGD,IDUBLD,SDFRAC,PTLAR
      COMMON /USER1/TITLE,PROJTY,TARGTY
      COMMON /USER2/CMENER,SDFRAC,PTLAR,ISTRUF,ISINGD,IDUBLD
*
      COMMON /STRUFU/ISTRUM,ISTRUT
      COMMON /ALALA/ALALAM
      COMMON/COLLIS/S,IJPROJ,IJTAR,PTTHR,PTTHR2,IOPHRD,IJPRLU,IJTALU
      COMMON /HAQQAP/ AQQAL,AQQPD,NQQAL,NQQPD
*
      DIMENSION SQS(13)
      DIMENSION SQSJ(17)
      DIMENSION XSQSJ(21),XXHHJ4(21)
*
*
*
*     used in ISIG=3,5,6,7,9
      DATA XSQSJ/0.005,0.01,0.02,0.035,0.053,
     * 0.1,0.2,0.35,0.54,1.,2.,5.,
     *10.,20.,40.,100.,200.,400.,1000.,2000.,4000./
*
*     used in ISIG 10,40
      DATA SQS/1.,2.,3.,4.,5.,10.,20.,30.,40.,100.,200.,500.,1000./
*
      PI4=4.*PI
      S=ECM**2
***********************************************************************
*----------------------------------------------------------------------
*
* **  select option used
*
      GO TO (10,20,30,40,50,60,70,80,90,100),ISIG
*
   10 CONTINUE
      WRITE(6,*)' This value of ISIG no longer available ISIG=',ISIG
*----------------------------------------------------------------------
*     nach:  Capella, Tran Thanh Van, Kwiecinski, PRL 58(1987)2015
*     as used in Ranft et al. SSC 149 Eq. 4,5
*
   20 CONTINUE
      WRITE(6,*)' This value of ISIG no longer available ISIG=',ISIG
*----------------------------------------------------------------------
*     nach:  Capella,Tran Thanh Van,Kwiecinski,PRL 58(1987)2015 CHG.SIG
*
   30 CONTINUE
      WRITE(6,*)' This value of ISIG no longer available ISIG=',ISIG
*----------------------------------------------------------------------
*     nach:  all subprocesses hard cross sections LIKE 70
C            PTHR RISING WITH ECM
*            NEW SIGTRP RISING LIKE LOG(S)   A.CAPELLA 30.3.90
C            ONLY FOR USE WITH PRBLM2!!!!!!!!!!!!!!!!
   40 CONTINUE
      WRITE(6,*)' This value of ISIG no longer available ISIG=',ISIG
*----------------------------------------------------------------------
*     nach:  all subprocesses hard cross sections LIKE 70
C            PTHR RISING WITH ECM  ptthr=3gev
*            NEW SIGTRP RISING LIKE LOG(S)   A.CAPELLA 30.3.90
C            ONLY FOR USE WITH PRBLM2!!!!!!!!!!!!!!!!
   50 CONTINUE
      WRITE(6,*)' This value of ISIG no longer available ISIG=',ISIG
*----------------------------------------------------------------------
*     nach:  all subprocesses hard cross sections LIKE 70
C            PTHR RISING WITH ECM    ptthr=2 gev
*            NEW SIGTRP RISING LIKE LOG(S)   A.CAPELLA 30.3.90
C            ONLY FOR USE WITH PRBLM2!!!!!!!!!!!!!!!!
   60 CONTINUE
      WRITE(6,*)' This value of ISIG no longer available ISIG=',ISIG
*----------------------------------------------------------------------
*     nach:  all subprocesses hard cross sections LIKE 70
C            PTHR RISING WITH ECM    ptthr=1.3 GEV
*            NEW SIGTRP RISING LIKE LOG(S)   A.CAPELLA 30.3.90
C            ONLY FOR USE WITH PRBLM2!!!!!!!!!!!!!!!!
   70 CONTINUE
      WRITE(6,*)' This value of ISIG no longer available ISIG=',ISIG
*----------------------------------------------------------------------
*     nach:  all subprocesses hard cross sections
*
   80 CONTINUE
      WRITE(6,*)' This value of ISIG no longer available ISIG=',ISIG
*----------------------------------------------------------------------
*     nach:  Patrick Maires program
   90 CONTINUE
      WRITE(6,*)' This value of ISIG no longer available ISIG=',ISIG
*----------------------------------------------------------------------
*     nach:  all subprocesses hard cross sections LIKE 70
C            PTHR RISING WITH ECM    ptthr=1.0 GEV
*            NEW SIGTRP RISING LIKE LOG(S)   A.CAPELLA 30.3.90
C            ONLY FOR USE WITH PRBLM2!!!!!!!!!!!!!!!!
      RETURN
  100 CONTINUE
*----------------------------------------------------------------------
*     nach:  all subprocesses hard cross sections LIKE 70
C            PTHR RISING WITH ECM
*            NEW SIGTRP RISING LIKE LOG(S)   A.CAPELLA 30.3.90
C            ONLY FOR USE WITH PRBLM2!!!!!!!!!!!!!!!!
C                                            INTRODUCED 18.12.90
C                                            BY DIETER PERTERMANN
C     modified 11-06-92 (R.Engel)
C
C  default parameter set
      ALFA=1.076
      ALFAP=0.24
      A=40.8
      BH=3.51
      BHOO=BH
      BSOO=BH
      AK=1.5
C     ALALAM --> ALAM (see PRBLM2 and POMDI)
      ALALAM=0.0
C                       begin fixed ptthr=3GeV
C--------------------------------------------------------------------
      IF(ABS(PTTHR-THREE).LT.EPS) THEN
        WRITE(6,*)' PTTHR=3. not available in dpmjet25'
          WRITE(6,*) ' WARNING: no model parameter set available'
          WRITE(6,*) ' for this combination of PTCUT and ISTRUF'
          WRITE(6,*) ' (initialization using default values)'
          ALFA = 1.078
          A    = 42.6
          ALALAM=0.740 
          AQQAL=1.D0
          AQQPD=1.D0
          ALFAP= 0.24
          AK=2.0
      ENDIF
C                       end fixed ptthr=3GeV
C                       begin fixed ptthr=2GeV
      IF(ABS(PTTHR-TWO).LT.EPS) THEN
        WRITE(6,*)' PTTHR=2. not available in dpmjet25'
          WRITE(6,*) ' WARNING: no model parameter set available'
          WRITE(6,*) ' for this combination of PTCUT and ISTRUF'
          WRITE(6,*) ' (initialization using default values)'
          ALFA = 1.042
          A    = 64.54
          ALALAM=0.6402
          AQQAL=1.D0
          AQQPD=1.D0
          ALFAP= 0.24
          AK=2.0
      ENDIF
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C                       end fixed ptthr=2GeV
C-----------------------------------------------------------
C-----------------------------------------------------------
C-----------------------------------------------------------
C     begin  ptthr= PTTHR=2.1+0.15*(LOG10(ECM/50.))**3
C-----------------------------------------------------------
      IF(ISTRUT.EQ.1) THEN
        WRITE(6,*)' ISTRUT=1 (PTTHR=2.1+0.15*(LOG10(ECM/50.))**3)',
     *	'not available in dpmjet25'
        PTTHR=2.1+0.15*(LOG10(ECM/50.))**3
        PTTHR2=PTTHR
          WRITE(6,*) ' WARNING: no model parameter set available'
          WRITE(6,*) ' for this combination of PTCUT and ISTRUF'
          WRITE(6,*) ' (initialization using default values)'
          ALFA = 1.042
          A    = 64.54
          ALALAM=0.6402
          AQQAL=1.D0
          AQQPD=1.D0
          ALFAP= 0.24
          AK=2.0
      ENDIF
C     end   PTTHR=2.1+0.15*(LOG10(ECM/50.))**3
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      BHOO=BH
      BSOO=BH
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C-----------------------------------------------------------
C-----------------------------------------------------------
C-----------------------------------------------------------
C    begin  PTTHR=2.5+0.12*(LOG10(ECM/50.))**3
C-----------------------------------------------------------
      IF(ISTRUT.EQ.2) THEN
        PTTHR=2.5+0.12*(LOG10(ECM/50.))**3
        PTTHR2=PTTHR
        IF( ISTRUF.EQ.9 ) THEN
        WRITE(6,*)' ISTRUT=2 (PTTHR=2.5+0.12*(LOG10(ECM/50.))**3)',
     *	'and ISTRUF= ',ISTRUF ,' not available in dpmjet25'
        GO TO 778
        ELSEIF( ISTRUF.EQ.10 ) THEN
        WRITE(6,*)' ISTRUT=2 (PTTHR=2.5+0.12*(LOG10(ECM/50.))**3)',
     *	'and ISTRUF= ',ISTRUF ,' not available in dpmjet25'
        GO TO 778
        ELSEIF( ISTRUF.EQ.11 ) THEN
        WRITE(6,*)' ISTRUT=2 (PTTHR=2.5+0.12*(LOG10(ECM/50.))**3)',
     *	'and ISTRUF= ',ISTRUF ,' not available in dpmjet25'
        GO TO 778
        ELSEIF( ISTRUF.EQ.12 ) THEN
        WRITE(6,*)' ISTRUT=2 (PTTHR=2.5+0.12*(LOG10(ECM/50.))**3)',
     *	'and ISTRUF= ',ISTRUF ,' not available in dpmjet25'
        GO TO 778
        ELSEIF( ISTRUF.EQ.13 ) THEN
        WRITE(6,*)' ISTRUT=2 (PTTHR=2.5+0.12*(LOG10(ECM/50.))**3)',
     *	'and ISTRUF= ',ISTRUF ,' not available in dpmjet25'
        GO TO 778
        ELSEIF( ISTRUF.EQ.14 ) THEN
        WRITE(6,*)' ISTRUT=2 (PTTHR=2.5+0.12*(LOG10(ECM/50.))**3)',
     *	'and ISTRUF= ',ISTRUF ,' not available in dpmjet25'
        GO TO 778
        ELSEIF( ISTRUF.EQ.15 ) THEN
        WRITE(6,*)' ISTRUT=2 (PTTHR=2.5+0.12*(LOG10(ECM/50.))**3)',
     *	'and ISTRUF= ',ISTRUF ,' not available in dpmjet25'
C  CETQ PDFs with other scale
        GO TO 778
        ELSEIF( ISTRUF.EQ.16 ) THEN
        WRITE(6,*)' ISTRUT=2 (PTTHR=2.5+0.12*(LOG10(ECM/50.))**3)',
     *	'and ISTRUF= ',ISTRUF ,' not available in dpmjet25'
          AK=2.0
        GO TO 778
        ELSEIF( ISTRUF.EQ.17 ) THEN
        WRITE(6,*)' ISTRUT=2 (PTTHR=2.5+0.12*(LOG10(ECM/50.))**3)',
     *	'and ISTRUF= ',ISTRUF ,' not available in dpmjet25'
        GO TO 778
        ELSEIF( ISTRUF.EQ.18 ) THEN
        WRITE(6,*)' ISTRUT=2 (PTTHR=2.5+0.12*(LOG10(ECM/50.))**3)',
     *	'and ISTRUF= ',ISTRUF ,' not available in dpmjet25'
        GO TO 778
        ELSEIF( ISTRUF.EQ.19 ) THEN
        WRITE(6,*)' ISTRUT=2 (PTTHR=2.5+0.12*(LOG10(ECM/50.))**3)',
     *	'and ISTRUF= ',ISTRUF ,' not available in dpmjet25'
        GO TO 778
        ELSEIF( ISTRUF.EQ.20 ) THEN
        WRITE(6,*)' ISTRUT=2 (PTTHR=2.5+0.12*(LOG10(ECM/50.))**3)',
     *	'and ISTRUF= ',ISTRUF ,' not available in dpmjet25'
C              GRV94LO AK=1. only for ISTRUT = 2
        GO TO 778
	ELSEIF( ISTRUF.EQ.21 ) THEN
          ALFA = 1.0733
          ALFAP= 0.171
          A    = 47.84
          ALALAM=0.621
          BSOO=1.58
          BHOO=3.54
          AK=1.000
C              GRV94LO AK=2. only for ISTRUT = 2
	ELSEIF( ISTRUF.EQ.22 ) THEN
          ALFA = 1.0513
          ALFAP= 0.3246
          A    = 55.16
          ALALAM=0.5846
          BSOO=1.114
          BHOO=1.703
          AK=2.000
C              CTEQ96 AK=2. only for ISTRUT = 2
	ELSEIF( ISTRUF.EQ.23 ) THEN
          ALFA = 1.0448
          ALFAP= 0.372
          A    = 57.51
          ALALAM=0.566
          BSOO=0.97
          BHOO=1.47
          AK=2.000
        ELSE
 778      CONTINUE	
          WRITE(6,*) ' WARNING: no model parameter set available'
          WRITE(6,*) ' for this combination of PTCUT and ISTRUF'
          WRITE(6,*) ' (initialization using default values)'
          ALFA = 1.042
          A    = 64.54
          ALALAM=0.6402
          AQQAL=1.D0
          AQQPD=1.D0
          ALFAP= 0.24
          AK=2.0
        ENDIF
        BHOO = BHOO/CONV
        BSOO = BSOO/CONV
      ENDIF
C    end  PTTHR=2.5+0.12*(LOG10(ECM/50.))**3
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C  slopes in GeV**-2
      BS=BSOO+ALFAP*LOG(S)
      BH=BHOO
*     BS=BSOO+ALFAP*LOG(S)
C  change units to mb
      BH=BH*CONV
      BS=BS*CONV
      BT=BS
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C     BT=BS/2.
      C=40.
C                                 CHANGED 13.1.90 BY J.R.
C     C=1.8
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C     C=1.E-8
*     parametrizations of input cross sections
      SIGSOF=A*S**(ALFA-1.)
* *** hard X-section
*     read interpolation data for different
*     sets of structure functions:
*
      CALL RDXSEC(XXHHJ4)
      IF(ISTRUF.EQ.21)AK=2.
      SIGHAR=1.E-8
      IF (S.GT.2450.D0)
     *               SIGHAR=AK*0.1*(S-2450.)**0.35
      IF(ECM.GE.THOUSA*XSQSJ(2)) THEN
        DO 1031 I=1,20
          III=I+1
          IF(ECM.LT.XSQSJ(III)*THOUSA.AND.
     *       ECM.GE.THOUSA*XSQSJ(I))THEN
            DSQ=ECM-THOUSA*XSQSJ(I)
            DDSQ=THOUSA*(XSQSJ(III)-XSQSJ(I))
            DHS=(XXHHJ4(III)-XXHHJ4(I))
            SIGHAR=AK*(XXHHJ4(I)+DHS*DSQ/DDSQ)*0.5
          ENDIF
 1031   CONTINUE
      ENDIF
C
C *** trippel pomeron X-section
C
C                                  VERSION A.CAPELLA 30.3.90
      GCA=SQRT(A)
      G3CA=GCA**3
      GACA=0.42
C     BSDOCA=1.372
      BSDOCA=BSOO*CONV
C     ALSCA=.0925
      ALSCA=ALFAP*CONV
      ALNS=LOG(S)
      BSDCA=BSDOCA+2.*ALSCA*ALNS
      SIGTRP=G3CA*GACA*LOG(S/10.)/(8.*3.14*BSDCA)
      IF (SIGTRP.LT.0.D0)SIGTRP=0.01
C
      BDDCA=2.*ALSCA*ALNS
      ALO1SQ=(LOG(S/400.))**2
      ALO2SQ=(LOG(25./S))**2
      ALO3SQ=(LOG(5./20.))**2
      SIGLOO=A*GACA**2*(ALO1SQ+ALO2SQ-2.*ALO3SQ)/(32.*3.14*BDDCA)
C
      ZSOF=SIGSOF/(PI4*BS)
      ZHAR=SIGHAR/(PI4*BH)
      ZTRP=SIGTRP/(PI4*BT)
      ZLOO=SIGLOO/(PI4*BT)
C 
      WRITE(6,'(2(/1X,A))') 'SELECTED PARAMETERS:',
     &                     '===================='
      WRITE(6,'(1X,A,E12.3)')   '  ALFA   ',ALFA
      WRITE(6,'(1X,A,E12.3)')   '  ALFAP  ',ALFAP
      WRITE(6,'(1X,A,E12.3)')   '  A      ',A
      WRITE(6,'(1X,A,2E12.3)')  '  BS,BSOO',BS,BSOO*CONV
      WRITE(6,'(1X,A,2E12.3)')  '  BH,BHOO',BH,BHOO*CONV
      WRITE(6,'(1X,A,E12.3)')   '  GACA   ',GACA
      WRITE(6,'(1X,A,E12.3,/)') '  AK     ',AK
C
      RETURN
      END
*
*
************************************************************************
************************************************************************
*
      SUBROUTINE RDXSEC(XSEC)
C
C     18.12.90, Dieter Pertermann
C     15.03.93, 27.05.93 modified (R. Engel)
C
C     RDXSEC READS DATA FOR INTERPOLATION
C     OF THE TOTAL CROSS SECTIONS FOR DIFFERENT
C     SETS OF STRUCTURE FUNCTIONS. THE CHOICE
C     OF THE CORRESPONDIG DATA SET IS CONTROLED
C     BY THE OVERALL STRUCTURE FUNCTION PARAMETER
C     ISTRUF.
C
C     CMENER taken out but no action necessary as not in SUB
C
C     modified 11-05-92 (R.Engel)
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      CHARACTER*80 TITLE
      CHARACTER*8 PROJTY,TARGTY
C     COMMON /USER/TITLE,PROJTY,TARGTY,CMENER,ISTRUF
C    &            ,ISINGD,IDUBLD,SDFRAC,PTLAR
      COMMON /USER1/TITLE,PROJTY,TARGTY
      COMMON /USER2/CMENER,SDFRAC,PTLAR,ISTRUF,ISINGD,IDUBLD
*
      COMMON/COLLIS/S,IJPROJ,IJTAR,PTTHR,PTTHR2,IOPHRD,IJPRLU,IJTALU
      COMMON /STRUFU/ISTRUM,ISTRUT
C
      DIMENSION XSEC(21)
      DIMENSION XS21(189)
C
      PARAMETER(EPSIL=1.D-4,
     &          THREE=3.D0,
     &          TWO=2.D0)
C
      DATA XS21  /
     &   0.000000E+00,0.137854E-04, .02, .13, .37, 1.32,
     &   3.88, 8.02, 13.15, 24.32, 43.43, 79.69, 113.13,
     &   147.5, 180.47, 221.01, 250.37,
     &   279.4, 320.1, 349.6, 381.6,
*  total X-section ,  cteq1M  PTTHR=x GEV/C
     &       .000000E+00, .494767E-05, .02, .14, .41, 
     &    1.48, 4.17, 7.92, 11.90, 19.03, 28.59, 42.36,
     &   52.78, 62.86, 72.65, 85.61, 95.97,
     &   96.,  96.,  96.,  96.,
*  total X-section ,  cteq1MS PTTHR=x GEV/C
     &      0.000000E+00,
     &      0.517461E-05, .02, .14, .42, 1.49, 4.14, 
     &    7.87, 11.93, 19.58, 30.67, 48.39, 63.08, 
     &    78.1, 93.28, 114.33, 132.24,
     &   133.,  133.,  133.,  133.,
*  total X-section ,  cteq1ML PTTHR=x GEV/C
     &      0.000000E+00,
     &      0.717097E-05, .03, .19, .54, 1.91, 5.33, 10.11,
     &     16.16, 24.21, 36.41, 54.21, 67.92, 81.44,
     &     94.81,112.9, 127.63,
     &    128.,  128.,  128.,  128.,   
*  total X-section ,  cteq1D  PTTHR=x GEV/C
     &      0.000000E+00,
     &      0.761464E-05, .02, .17, .47, 1.56, 4.19, 
     &     7.76, 11.48, 18.11, 26.97, 39.82, 49.86, 59.35, 
     &     68.88, 81.65, 91.94,
     &    92.,  92.,  92.,  92.,
*  total X-section ,  cteq1L  PTTHR=x GEV/C
     &       .000000E+00,
     &       .620779E-05, .02, .12, .34, 1.19, 3.27, 
     &      6.16, 9.27, 14.99, 23.2, 36.85, 49.45, 
     &      64.43, 82.38, 112.06, 140.36,
     &    141.,  141.,  141.,  141.,
*  total X-section ,  GRV94LO  AK=1. PTTHR=x GEV/C
     &       .000000E+00,
     &       .620779E-05, .01, .05, .14, 0.55, 1.87, 
     &      4.29,  7.49, 14.81, 27.8, 55.99, 77.49, 
     &     105.98,138.48, 189.33, 236.37,
     &    294.,  395.,  496.,  629.,
*  total X-section ,  GRV94LO AK=2. PTTHR=x GEV/C
     &       .000000E+00,
     &       .620779E-05, .01, .10, .31, 1.16, 3.76, 
     &      8.31, 14.16, 27.11, 49.3, 90.93,129.77, 
     &     174.16,223.83, 300.20, 370.00,
     &    455.,  600.,  746.,  936.,
*  total X-section ,  CTEC96 AK=2. PTTHR=x GEV/C
     &       .000000E+00,
     &       .620779E-05, .01, .08, .27, 1.17, 4.15, 
     &      9.60, 16.75, 32.88, 61.1,125.98,169.87, 
     &     233.75,308.22, 426.95, 537.90,
     &    673.,  898., 1112., 1379./
*******************************************************************
*
      IF( ABS(PTTHR-THREE).LT.EPSIL )     THEN
          WRITE(6,*) ' ERROR RDXSEC: invalid pdf No. ',ISTRUF
          STOP
      ELSEIF( ABS(PTTHR-TWO).LT.EPSIL ) THEN
          WRITE(6,*) ' ERROR RDXSEC: invalid pdf No. ',ISTRUF
          STOP
      ELSEIF( ISTRUT.EQ.1 ) THEN
          WRITE(6,*) ' ERROR RDXSEC: invalid pdf No. ',ISTRUF
          STOP
      ELSEIF( ISTRUT.EQ.2 ) THEN
        IF( (ISTRUF.GE.9).AND.(ISTRUF.LE.20) ) THEN
          WRITE(6,*) ' ERROR RDXSEC: invalid pdf No. ',ISTRUF
          STOP
        ELSEIF( (ISTRUF.GE.21).AND.(ISTRUF.LE.23) ) THEN
          DO 311 I=1,21
            NXS = 21*(ISTRUF-15)+I
            XSEC(I)=XS21(NXS)
  311    CONTINUE
        ELSE
          WRITE(6,*) ' ERROR RDXSEC: invalid pdf No. ',ISTRUF
          STOP
        ENDIF
      ELSE
        WRITE(6,*) ' ERROR RDXSEC: PTCUT ',PTTHR,' not supported ***'
        STOP
      ENDIF
C
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE PRBLM2(ECM)
*     *     input:
*        ECM
*     output:
*        PLMN, PLMNCUmmulative, AVSOFN, AVHRDN,SIGDD/D/QEL/EL, PSOFT
*
*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*
*        Probabilities for L soft cut, M hard cut pomerons
*           and N someplace cut trippel pomerons
*
*         Aurenche Maire's  version of PRBLM
*         modified from A.M. to include calculation of x-sections
*         modified to get:
*      version with high mass diffraction (Y's and PHI's)
*              and on both sides 2 kanal low mass diffraction
* ***                                     OPTION 21.2.90 FB
*----------------------------------------------------------------------
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( ZERO=0.D0, ONE=1.D0)
      PARAMETER (CONV=0.38935D0)
      PARAMETER (PI=3.141592654D0)
      PARAMETER (MXPA25=30,MXPA26=MXPA25+1,MXPA13=13)
*     PARAMETRIZATION FOR PTMIN= 3. GEV
      PARAMETER (MXPA50=250,MXPA51=MXPA50+1)
*     PARAMETRIZATION FOR PTMIN= 2. GEV
C     PARAMETER (MXPA50=350,MXPA51=MXPA50+1)
      PARAMETER (MXPA96=96)
C     PARAMETER (MXPA96=480)
      LOGICAL LSQRT
C      PARAMETER (MXLMN=5,LSQRT=.false.)
      PARAMETER (MXLMN=5,LSQRT=.TRUE.)
      DOUBLE PRECISION DTINY
C     PARAMETER (TINY = 1.2D-38,DTINY=1.D-70,TIN=1.D-22,TINEXP = -300.D0)
C     PARAMETER (TINY = 1.2D-38,DTINY=1.D-300,TIN=1.D-22,TINEXP =
C    -700.D0)
      PARAMETER (TINY=1.2D-38,DTINY=1.D-70,TIN=1.D-22,TINEXP=-700.D0)
C      PARAMETER (TINY = 1.D-38,DTINY=tiny,TIN=1.D-22,TINEXP = -300.D0)
*     in older version used:
      PARAMETER (TINYEX  = -48.D0)
* --- -- - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - -
* *** /POLMN/ arrays having to do with cut soft and hard Pomerons
      COMMON /POLMN/PLMN(0:MXPA25,0:MXPA50,0:MXPA13),
     *              PLMNCU(0:MXPA25,0:MXPA50,0:MXPA13)
      COMMON /POLMN0/PDIFR,PHARD,PSOFT,ALFAH,BETAH,
     *              SIGTOT,SIGQEL,SIGEL,SIGINE,SIGHIN,SIGD,SIGDD
* *** /POMPAR/*/SIGMA/*/POMTYP/ used only in SIGMAPOM-routines (->POMDI)
*      (LMAX,MMAX,NMAX (max number of soft/hard/trippel pomerons)
      COMMON /POMTYP/IPIM,ICON,ISIG,LMAX,MMAX,NMAX,DIFEL,DIFNU
      COMMON /SIGMA/SIGSOF,BS,ZSOF,SIGHAR,BH,ZHAR,SIGTRP,BT,ZTRP,
     *              SIGLOO,ZLOO
      COMMON/POMPAR/ALFA,ALFAP,A,C,AK
      COMMON /SINGDI/SILMSD,SIGDI
* *** /OUTLEV/ controls output level for POMDI and parton X distribution
      COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
      COMMON /ALALA/ALALAM
      CHARACTER*80 TITLE
      CHARACTER*8 PROJTY,TARGTY
C     COMMON /USER/TITLE,PROJTY,TARGTY,CMENER,ISTRUF
C    &            ,ISINGD,IDUBLD,SDFRAC,PTLAR
      COMMON /USER1/TITLE,PROJTY,TARGTY
      COMMON /USER2/CMENER,SDFRAC,PTLAR,ISTRUF,ISINGD,IDUBLD
* --- -- - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - -
      DOUBLE PRECISION SIG,SIGP,SIGM,SIGN,SIGO
      DIMENSION SIG(0:MXPA25,0:MXPA50,0:MXPA13),
     &SIGP(0:MXPA25,0:MXPA50,0:MXPA13),SIGM(0:MXPA25,0:MXPA50,0:MXPA13),
     &SIGN(0:MXPA25,0:MXPA50,0:MXPA13),SIGO(0:MXPA25,0:MXPA50,0:MXPA13)
      DIMENSION XPNT(MXPA96),WGHT(MXPA96),
     &SSOFT(0:MXPA25),SHARD(0:MXPA50),STRPL(0:MXPA25)
C       - - required MXPA25 > NMAX - -
      DIMENSION FAK(0:MXPA13),CMBIN(0:MXPA13,0:MXPA13)
      DOUBLE PRECISION
     &       EXPSOP,EXPSOH,EXMSOP,EXMSOH,EXNSOP,EXNSOH,EXOSOP,EXOSOH,
     &       EXPHAP,EXPHAH,EXMHAP,EXMHAH,EXNHAP,EXNHAH,EXOHAP,EXOHAH,
     &       EXPTRP,EXPTRH,EXMTRP,EXMTRH,EXNTRP,EXNTRH,EXOTRP,EXOTRH,
     &       EXPLOP,EXPLOH,EXMLOP,EXMLOH,EXNLOP,EXNLOH,EXOLOP,EXOLOH,
     &       EXPEXH,EXMEXH,EXNEXH,EXOEXH,EXPEXP,EXMEXP,EXNEXP,EXOEXP
      DOUBLE PRECISION  FAPSOF,FAMSOF,FANSOF,FAOSOF,
     &                  FAPHAR,FAMHAR,FANHAR,FAOHAR,
     &                  FAPTRP,FAMTRP,FANTRP,FAOTRP,
     &                  FAPLOO,FAMLOO,FANLOO,FAOLOO
      DOUBLE PRECISION  DENOM,DENOMI,XPNTK,WGHTK,RMXLMN
     &                  ,SIGSUM,SIGINL,SIGHRI
*
*---------------------------------------------------------------------------
*
*     externe ICON option to internes NMAX=1,2,free
      IF(ICON/10.EQ.4) NMAX=2
      IF(ICON/10.EQ.5) NMAX=1
*
*     for safty
      IF( NMAX.GT.MXPA13) THEN
        WRITE(6,*)' arrays limit NMAX set to' , MXPA13
        NMAX=MXPA13
      ENDIF
      IF( MMAX.GT.MXPA50) THEN
        WRITE(6,*)' arrays limit MMAX set to' , MXPA50
        MMAX=MXPA50
      ENDIF
      IF( lMAX.GT.MXPA25) THEN
        WRITE(6,*)' arrays limit LMAX set to' , MXPA25
        LMAX=MXPA25
      ENDIF
*
      LMAXI = LMAX
      MMAXI = MMAX
      IF( NMAX.GE.3)THEN
        NMAXI = NMAX
        NNMAXI=(MXPA13-NMAXI)/(1+NMAXI)
C       aim: MXPA13 =!= NNNMAX =  NMAXI+(NMAXI+1)*NNMAXI
        NLMAXI=0
      ELSEIF( NMAX.EQ.2)THEN
        NMAXI=1
        NNMAXI=1
        NLMAXI=1
      ELSEIF( NMAX.EQ.1)THEN
        NMAXI=1
        NNMAXI=0
        NLMAXI=1
      ELSEIF( NMAX.LE.0)THEN
        NMAXI=0
        NNMAXI=0
        NLMAXI=0
      ENDIF
*
*
      LENTRY=0
*
      GOTO 111
*
*----------------------------------------------------------------------
*
      ENTRY SIGMA2(ECM)
*
*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*
      LENTRY=1
*     externe ICON option to internes NMAX=1,2,free
      IF(ICON/10.EQ.4) NMAX=2
      IF(ICON/10.EQ.5) NMAX=1
*     we drop L,M..dependent quantities, the rest is integrated for L=1,rest=0
      LMAXI = 1
      MMAXI = 0
      NMAXI = 0
      NNMAXI= 0
      NLMAXI= 0
*
*----------------------------------------------------------------------
*
* *** calculate the B-space integral for L/M soft/hard cut-pom.
* *** NPNT is the number of integration points of B-space integ.
*
 111  SIGTOT=0.0
      SIGINL=0.0
      SIGEL =0.0
      SIGELE=0.0
      SIGD  =0.0
      SIGDD =0.0
      SIGDI =0.0
      SIGDDI=0.0
      SIGINE=0.0
      SIHMDD=0.0
      SIGHMD=0.
      SIGLMD=0.
      SIGHIN=0.
      SIGIN=0.
      SIGSIN=0.0
      SIGHRI=0.0
      SIGSUM=0.0
      SIGSME=0.0
      SILMSD=0.0
      SILMDD=0.
      SLHMDD=0.
      DO 10 L=0,LMAXI
       SSOFT(L)=0.
       STRPL(L)=0.
       DO 10 M=0,MMAXI
        SHARD(M)=0.
        DO 10 N=0, MXPA13
          SIG(L,M,N)=0.
          SIGP(L,M,N)=0.
          SIGM(L,M,N)=0.
          SIGN(L,M,N)=0.
          SIGO(L,M,N)=0.
          PLMN(L,M,N)=0.
          PLMNCU(L,M,N)=0.
10    CONTINUE
*
*     get bare X-sections
      S=ECM**2
      CALL SIGSHD(ECM)
*
      IF(ALALAM.LE.1.D-2) THEN
        ALAM=0.60
      ELSE
        ALAM=ALALAM
      ENDIF
*
*     prepare integration:
      PI4 = 4.*PI
*
      IF(ECM.LT.2000.D0)THEN
        NPNT=96
        CALL GSET(ZERO,ONE,NPNT,XPNT,WGHT)
      ELSE
        NPNT=MXPA96
        CALL GSET(ZERO,ONE,NPNT,XPNT,WGHT)
      ENDIF
*
*     as low mass diffraction extra, high mass reduced:
      REDU = 1.0
      IF(IOUTPO.GE.0) WRITE (6,*) ' ALAM,REDU= ',ALAM,REDU
*
*     --here the versions started--
*     prepare factors to enter sum:
*     notation: Z(HARd/SOFt/LOOp)(Plus/Minus/N=mixed/O=mixed)
*
      ZHARP=(1.+ALAM)**2*ZHAR
      ZSOFP=(1.+ALAM)**2*ZSOF
      ZLOOP=(1.+ALAM)**2*ZLOO * REDU
      ZHARM=(1.-ALAM)**2*ZHAR
      ZSOFM=(1.-ALAM)**2*ZSOF
      ZLOOM=(1.-ALAM)**2*ZLOO * REDU
      ZHARN=(1.-ALAM**2)*ZHAR
      ZSOFN=(1.-ALAM**2)*ZSOF
      ZLOON=(1.-ALAM**2)*ZLOO * REDU
      ZHARO=(1.-ALAM**2)*ZHAR
      ZSOFO=(1.-ALAM**2)*ZSOF
      ZLOOO=(1.-ALAM**2)*ZLOO * REDU
*     one more factor at the top or at the bottom:
      ZTRPP=(1.+ALAM)**3*ZTRP * REDU
      ZTRPM=(1.-ALAM)**3*ZTRP * REDU
      ZTRPN=(1.-ALAM**2)*(1.+ALAM)*ZTRP * REDU
      ZTRPO=(1.-ALAM**2)*(1.-ALAM)*ZTRP * REDU
*
*     begin M,N,L,LL loop
*
      DO 720 L=0,LMAXI
        IF(L.EQ.0) THEN
          FAPSOF=1.
          FAMSOF=1.
          FANSOF=1.
          FAOSOF=1.
        ELSEIF(LSQRT) THEN
          FAPSOF=FAPSOF* SQRT( ZSOFP/FLOAT(L))
          FAMSOF=FAMSOF* SQRT( ZSOFM/FLOAT(L))
          FANSOF=FANSOF* SQRT( ZSOFN/FLOAT(L))
          FAOSOF=FAOSOF* SQRT( ZSOFO/FLOAT(L))
          IF (    FAPSOF .LT.DTINY )     FAPSOF=0.
          IF (    FAMSOF .LT.DTINY )     FAMSOF=0.
          IF (    FANSOF .LT.DTINY )     FANSOF=0.
          IF (    FAOSOF .LT.DTINY )     FAOSOF=0.
        ELSEIF(.NOT.LSQRT) THEN
          FAPSOF=FAPSOF*ZSOFP/FLOAT(L)
          FAMSOF=FAMSOF*ZSOFM/FLOAT(L)
          FANSOF=FANSOF*ZSOFN/FLOAT(L)
          FAOSOF=FAOSOF*ZSOFO/FLOAT(L)
          IF (FAPSOF.LT.DTINY )     FAPSOF=0.
          IF (FAMSOF.LT.DTINY )     FAMSOF=0.
          IF (FANSOF.LT.DTINY )     FANSOF=0.
          IF (FAOSOF.LT.DTINY )     FAOSOF=0.
        ENDIF
        DO 730 M=0,MMAXI
          IF(M.EQ.0) THEN
             FAPHAR=1.
             FAMHAR=1.
             FANHAR=1.
             FAOHAR=1.
          ELSEIF(LSQRT) THEN
C           WRITE(6,*)FAPHAR,ZHARP/FLOAT(M),FAMHAR,ZHARM/FLOAT(M)
            FAPHAR=FAPHAR* SQRT( ZHARP/FLOAT(M) )
            FAMHAR=FAMHAR* SQRT( ZHARM/FLOAT(M) )
            FANHAR=FANHAR* SQRT( ZHARN/FLOAT(M) )
            FAOHAR=FAOHAR* SQRT( ZHARO/FLOAT(M) )
            IF (    FAPSOF*FAPHAR .LT.DTINY )     FAPHAR=0.
            IF (    FAMSOF*FAMHAR .LT.DTINY )     FAMHAR=0.
            IF (    FANSOF*FANHAR .LT.DTINY )     FANHAR=0.
            IF (    FAOSOF*FAOHAR .LT.DTINY )     FAOHAR=0.
          ELSEIF(.NOT.LSQRT) THEN
            FAPHAR=FAPHAR*ZHARP/FLOAT(M)
            FAMHAR=FAMHAR*ZHARM/FLOAT(M)
            FANHAR=FANHAR*ZHARN/FLOAT(M)
            FAOHAR=FAOHAR*ZHARO/FLOAT(M)
            IF (FAPSOF*FAPHAR.LT.DTINY )     FAPHAR=0.
            IF (FAMSOF*FAMHAR.LT.DTINY )     FAMHAR=0.
            IF (FANSOF*FANHAR.LT.DTINY )     FANHAR=0.
            IF (FAOSOF*FAOHAR.LT.DTINY )     FAOHAR=0.
          ENDIF
          DO 740 N=0,NMAXI
            IF( N.EQ.0) THEN
               FAPTRP=1.
               FAMTRP=1.
               FANTRP=1.
               FAOTRP=1.
            ELSEIF(LSQRT) THEN
               FAPTRP=-FAPTRP* SQRT( ZTRPP/FLOAT(N) )
               FAMTRP=-FAMTRP* SQRT( ZTRPM/FLOAT(N) )
               FANTRP=-FANTRP* SQRT( ZTRPN/FLOAT(N) )
               FAOTRP=-FAOTRP* SQRT( ZTRPO/FLOAT(N) )
               IF (ABS(FAPTRP*FAPSOF*FAPHAR).LT.DTINY )     FAPTRP=0.
               IF (ABS(FAMTRP*FAMSOF*FAMHAR).LT.DTINY )     FAMTRP=0.
               IF (ABS(FANTRP*FANSOF*FANHAR).LT.DTINY )     FANTRP=0.
               IF (ABS(FAOTRP*FAOSOF*FAOHAR).LT.DTINY )     FAOTRP=0.
            ELSEIF(.NOT.LSQRT) THEN
               FAPTRP=-FAPTRP*ZTRPP/FLOAT(N)
               FAMTRP=-FAMTRP*ZTRPM/FLOAT(N)
               FANTRP=-FANTRP*ZTRPN/FLOAT(N)
               FAOTRP=-FAOTRP*ZTRPO/FLOAT(N)
               IF (ABS(FAPTRP*FAPSOF*FAPHAR).LT.DTINY )     FAPTRP=0.
               IF (ABS(FAMTRP*FAMSOF*FAMHAR).LT.DTINY )     FAMTRP=0.
               IF (ABS(FANTRP*FANSOF*FANHAR).LT.DTINY )     FANTRP=0.
               IF (ABS(FAOTRP*FAOSOF*FAOHAR).LT.DTINY )     FAOTRP=0.
            ENDIF
            DO 750 NN=0,NNMAXI
*             for compatibility no new subscript is introduced in some arrays:
              NNN=N+(NMAXI+1)*NN
*             if only first order option jump out of second order contr.:
              IF( NMAX.LE.2   .AND. N.EQ.1 .AND. NN.EQ.1 ) GO TO 750
              IF(NN.EQ.0) THEN
                FAPLOO=1.
                FAMLOO=1.
                FANLOO=1.
                FAOLOO=1.
              ELSEIF(LSQRT) THEN
                FAPLOO=-FAPLOO* SQRT( ZLOOP/FLOAT(NN))
                FAMLOO=-FAMLOO* SQRT( ZLOOM/FLOAT(NN))
                FANLOO=-FANLOO* SQRT( ZLOON/FLOAT(NN))
                FAOLOO=-FAOLOO* SQRT( ZLOOO/FLOAT(NN))
                IF(ABS(FAPLOO*FAPTRP*FAPSOF*FAPHAR).LT.DTINY )FAPLOO=0.
                IF(ABS(FAMLOO*FAMTRP*FAMSOF*FAMHAR).LT.DTINY )FAMLOO=0.
                IF(ABS(FANLOO*FANTRP*FANSOF*FANHAR).LT.DTINY )FANLOO=0.
                IF(ABS(FAOLOO*FAOTRP*FAOSOF*FAOHAR).LT.DTINY )FAOLOO=0.
              ELSEIF(.NOT.LSQRT) THEN
                FAPLOO=-FAPLOO*ZLOOP/FLOAT(NN)
                FAMLOO=-FAMLOO*ZLOOM/FLOAT(NN)
                FANLOO=-FANLOO*ZLOON/FLOAT(NN)
                FAOLOO=-FAOLOO*ZLOOO/FLOAT(NN)
                IF(ABS(FAPLOO*FAPTRP*FAPSOF*FAPHAR).LT.DTINY )FAPLOO=0.
                IF(ABS(FAMLOO*FAMTRP*FAMSOF*FAMHAR).LT.DTINY )FAMLOO=0.
                IF(ABS(FANLOO*FANTRP*FANSOF*FANHAR).LT.DTINY )FANLOO=0.
                IF(ABS(FAOLOO*FAOTRP*FAOSOF*FAOHAR).LT.DTINY )FAOLOO=0.
              ENDIF
*
*         elastic processes are not generated
          IF(L.EQ.0.AND.M.EQ.0.AND.N.EQ.0.AND.NN.EQ.0) GO TO 750
*
          DENOM=DBLE(M)/DBLE(BH)+DBLE(L)/DBLE(BS)+DBLE(N)/DBLE(BT)
     &    +DBLE(NN)/DBLE(BT)
*
          DO 735 K=1,NPNT
*
C           change intergration for large L+M+N+NN?
            IF ( (M+L+N+NN) .LE. MXLMN  ) THEN
              XPNTK=DBLE(XPNT(K))
              WGHTK=DBLE(WGHT(K))
              DENOMI=DENOM
            ELSE
              RMXLMN = DBLE(M+L+N+NN) /DBLE(MXLMN)
              XPNTK=DBLE(XPNT(K))
              WGHTK= DBLE(WGHT(K)) * XPNTK**(RMXLMN-1.)
              DENOMI= DENOM / RMXLMN
            ENDIF
*
            EXPOSP=-ZSOFP*XPNTK**(1./(DENOMI*DBLE(BS)))
            EXPOSM=-ZSOFM*XPNTK**(1./(DENOMI*DBLE(BS)))
            EXPOSN=-ZSOFN*XPNTK**(1./(DENOMI*DBLE(BS)))
            EXPOSO=-ZSOFO*XPNTK**(1./(DENOMI*DBLE(BS)))
*
            EXPOHP=-ZHARP*XPNTK**(1./(DENOMI*DBLE(BH)))
            EXPOHM=-ZHARM*XPNTK**(1./(DENOMI*DBLE(BH)))
            EXPOHN=-ZHARN*XPNTK**(1./(DENOMI*DBLE(BH)))
            EXPOHO=-ZHARO*XPNTK**(1./(DENOMI*DBLE(BH)))
*
            EXPOTP=+ZTRPP*XPNTK**(1./(DENOMI*DBLE(BT)))
            EXPOTM=+ZTRPM*XPNTK**(1./(DENOMI*DBLE(BT)))
            EXPOTN=+ZTRPN*XPNTK**(1./(DENOMI*DBLE(BT)))
            EXPOTO=+ZTRPO*XPNTK**(1./(DENOMI*DBLE(BT)))
*
            EXPOLP=+ZLOOP*XPNTK**(1./(DENOMI*DBLE(BT)))
            EXPOLM=+ZLOOM*XPNTK**(1./(DENOMI*DBLE(BT)))
            EXPOLN=+ZLOON*XPNTK**(1./(DENOMI*DBLE(BT)))
            EXPOLO=+ZLOOO*XPNTK**(1./(DENOMI*DBLE(BT)))
*
            IF(IOUTPO.GE.7) THEN
              WRITE(6,*)
     *         ' K=',K,' EXPOS/H=',EXPOSP,EXPOHP,' DENOMI/BH=',DENOMI,BH
              WRITE(6,*)
     *         ' K=',K,' EXPOS/H=',EXPOSM,EXPOHM,' DENOMI/BH=',DENOMI,BH
              WRITE(6,*)
     *         ' K=',K,' EXPOS/H=',EXPOSN,EXPOHN,' DENOMI/BH=',DENOMI,BH
              WRITE(6,*)
     *          ' K=',K,'XPNT=',XPNTK,'WGHT=',WGHTK,'DENO=',DENOMI
            ENDIF
*
*           notation:
*           EX(P=+/M=-/N=mixed/O=mixed)(EX=all/HArd,TRippel,LOop)(P/Half)
*
            IF(     EXPOSP .GT. TINEXP) THEN
              EXPSOH=EXP(0.5D00*EXPOSP)
              EXMSOH=EXP(0.5D00*EXPOSM)
              EXNSOH=EXP(0.5D00*EXPOSN)
              EXOSOH=EXP(0.5D00*EXPOSO)
            ELSE
              EXPSOH=0.
              EXMSOH=0.
              EXNSOH=0.
              EXOSOH=0.
            ENDIF
            EXPSOP=EXPSOH**2
            EXMSOP=EXMSOH**2
            EXNSOP=EXNSOH**2
            EXOSOP=EXOSOH**2
*
            IF(    EXPOHP .GT. TINEXP) THEN
              EXPHAH=EXP(0.5D00*EXPOHP)
              EXMHAH=EXP(0.5D00*EXPOHM)
              EXNHAH=EXP(0.5D00*EXPOHN)
              EXOHAH=EXP(0.5D00*EXPOHO)
            ELSE
              EXPHAH=0.
              EXMHAH=0.
              EXNHAH=0.
              EXOHAH=0.
            ENDIF
            EXPHAP=EXPHAH**2
            EXMHAP=EXMHAH**2
            EXNHAP=EXNHAH**2
            EXOHAP=EXOHAH**2
*
            IF( NMAX.GE.3) THEN
              IF( EXPOTP .GT. TINEXP) THEN
                EXPTRH=EXP(0.5D00*EXPOTP)
                EXMTRH=EXP(0.5D00*EXPOTM)
                EXNTRH=EXP(0.5D00*EXPOTN)
                EXOTRH=EXP(0.5D00*EXPOTO)
              ELSE
                EXPTRH=0.
                EXMTRH=0.
                EXNTRH=0.
                EXOTRH=0.
              ENDIF
              EXPTRP= EXPTRH**2
              EXMTRP= EXMTRH**2
              EXNTRP= EXNTRH**2
              EXOTRP= EXOTRH**2
            ELSEIF( NMAX.LE.2) THEN
                EXPTRH= 1 + 0.5*EXPOTP
                EXMTRH= 1 + 0.5*EXPOTM
                EXNTRH= 1 + 0.5*EXPOTN
                EXOTRH= 1 + 0.5*EXPOTO
                EXPTRP= 1 + EXPOTP
                EXMTRP= 1 + EXPOTM
                EXNTRP= 1 + EXPOTN
                EXOTRP= 1 + EXPOTO
            ENDIF
*
            IF( NMAX.GE.3) THEN
              IF( EXPOLP .GT. TINEXP) THEN
                EXPLOH=EXP(0.5D00*EXPOLP)
                EXMLOH=EXP(0.5D00*EXPOLM)
                EXNLOH=EXP(0.5D00*EXPOLN)
                EXOLOH=EXP(0.5D00*EXPOLO)
              ELSE
                EXPLOH=0.
                EXMLOH=0.
                EXNLOH=0.
                EXOLOH=0.
              ENDIF
              EXPLOP=EXPLOH**2
              EXMLOP=EXMLOH**2
              EXNLOP=EXNLOH**2
              EXOLOP=EXOLOH**2
            ELSEIF( NMAX.EQ.2 ) THEN
                EXPLOH= 1 + 0.5*EXPOLP
                EXMLOH= 1 + 0.5*EXPOLM
                EXNLOH= 1 + 0.5*EXPOLN
                EXOLOH= 1 + 0.5*EXPOLO
                EXPLOP= 1 + EXPOLP
                EXMLOP= 1 + EXPOLM
                EXNLOP= 1 + EXPOLN
                EXOLOP= 1 + EXPOLO
            ELSEIF( NMAX.LE.1 ) THEN
                EXPLOH= 1
                EXMLOH= 1
                EXNLOH= 1
                EXOLOH= 1
                EXPLOP= 1
                EXMLOP= 1
                EXNLOP= 1
                EXOLOP= 1
            ENDIF
*
            EXPEXH = EXPSOH *EXPHAH *EXPTRH *EXPLOH
            EXMEXH = EXMSOH *EXMHAH *EXMTRH *EXMLOH
            EXNEXH = EXNSOH *EXNHAH *EXNTRH *EXNLOH
            EXOEXH = EXOSOH *EXOHAH *EXOTRH *EXOLOH
            EXPEXP = EXPSOP *EXPHAP *EXPTRP *EXPLOP
            EXMEXP = EXMSOP *EXMHAP *EXMTRP *EXMLOP
            EXNEXP = EXNSOP *EXNHAP *EXNTRP *EXNLOP
            EXOEXP = EXOSOP *EXOHAP *EXOTRP *EXOLOP
*
            IF( ( NMAX.LE.2  .AND.  N.EQ.1 ) .OR.
     *          ( NMAX.EQ.2  .AND. NN.EQ.1 ) .OR.
     *            NMAX.EQ.0                      ) THEN
               SIGP(L,M,NNN)=SIGP(L,M,NNN)+EXPSOP *EXPHAP *WGHTK
               SIGM(L,M,NNN)=SIGM(L,M,NNN)+EXMSOP *EXMHAP *WGHTK
               SIGN(L,M,NNN)=SIGN(L,M,NNN)+EXNSOP *EXNHAP *WGHTK
               SIGO(L,M,NNN)=SIGO(L,M,NNN)+EXOSOP *EXOHAP *WGHTK
             ELSE
               SIGP(L,M,NNN)=SIGP(L,M,NNN)+EXPEXP*WGHTK
               SIGM(L,M,NNN)=SIGM(L,M,NNN)+EXMEXP*WGHTK
               SIGN(L,M,NNN)=SIGN(L,M,NNN)+EXNEXP*WGHTK
               SIGO(L,M,NNN)=SIGO(L,M,NNN)+EXOEXP*WGHTK
             ENDIF
*
*           quantities without L,M,N,NN dependence considered ones
*                              (when, chosen to get suitable weights)
            IF(L.EQ.1.AND.M.EQ.0.AND.N.EQ.0.AND.NN.EQ.0) THEN
*
              IF ( (M+L+N+NN) .GT. MXLMN  ) THEN
                WRITE(6,*)' MXLMN too low ' , MXLMN,M,L,N,NN
                RETURN
              ENDIF
              WGHFAC = WGHTK/XPNTK *PI4/DENOMI
              IF     ( NMAX.GE.3 ) THEN
                SIGELE = SIGELE + WGHFAC *
     *              0.0625*(   1.-EXPEXH   +  1.-EXMEXH
     *                        +1.-EXNEXH   +  1.-EXOEXH )**2
*               low mass diffraction:
                SILMSD = SILMSD + WGHFAC *
     *                0.125*(EXPEXH -EXMEXH)**2
                SILMDD = SILMDD + WGHFAC *
     *                0.0625*(EXPEXH+EXMEXH-EXNEXH-EXOEXH)**2
              ELSEIF( NMAX.LE.2 ) THEN
                SIGELE = SIGELE + WGHFAC *
     *           0.0625*( (   1.-EXPEXH   +  1.-EXMEXH
     *                       +1.-EXNEXH   +  1.-EXOEXH
*                            subtract second order terms in each factor:
C////CHANGED - TO + BRACKET UNNECESSARY
     *                      +(1.-EXPTRH)*(1-EXPLOH) *EXPSOH *EXPHAH
     *                      +(1.-EXMTRH)*(1-EXMLOH) *EXMSOH *EXMHAH
     *                      +(1.-EXNTRH)*(1-EXNLOH) *EXNSOH *EXNHAH
     *                      +(1.-EXOTRH)*(1-EXOLOH) *EXOSOH *EXOHAH)**2
*                            subtract second order terms of product:
     *                   - (  (2.-EXPTRH-EXPLOH) *EXPSOH *EXPHAH
     *                       +(2.-EXMTRH-EXMLOH) *EXMSOH *EXMHAH
     *                       +(2.-EXNTRH-EXNLOH) *EXNSOH *EXNHAH
     *                       +(2.-EXOTRH-EXOLOH) *EXOSOH *EXOHAH ) **2)
*               low mass diffraction:
                SILMSD = SILMSD + WGHFAC *
     *            0.125*( ( EXPEXH -EXMEXH
*                           subtract second order terms in each factor:
     *                     -(1.-EXPTRH)*(1-EXPLOH) *EXPSOH*EXPHAH
     *                     +(1.-EXMTRH)*(1-EXMLOH) *EXMSOH*EXMHAH )**2
*                           subtract second order terms of product:
     *                     -(  (2.-EXPTRH-EXPLOH) *EXPSOH *EXPHAH
     *                        -(2.-EXMTRH-EXMLOH) *EXMSOH*EXMHAH ) **2)
                SILMDD = SILMDD + WGHFAC *
     *           0.0625*( (EXPEXH+EXMEXH-EXNEXH-EXOEXH
*                         subtract second order terms in each factor:
     *                   -(1.-EXPTRH)*(1-EXPLOH) *EXPSOH *EXPHAH
     *                   -(1.-EXMTRH)*(1-EXMLOH) *EXMSOH *EXMHAH
     *                   +(1.-EXNTRH)*(1-EXNLOH) *EXNSOH *EXNHAH
     *                   +(1.-EXOTRH)*(1-EXOLOH) *EXOSOH *EXOHAH)**2
*                         subtract second order terms of product:
     *                - (  (2.-EXPTRH-EXPLOH) *EXPSOH *EXPHAH
     *                    +(2.-EXMTRH-EXMLOH) *EXMSOH *EXMHAH
     *                    -(2.-EXNTRH-EXNLOH) *EXNSOH *EXNHAH
     *                    -(2.-EXOTRH-EXOLOH) *EXOSOH *EXOHAH ) **2)
              ENDIF
              IF( NMAX.NE.2 ) THEN
                SIGTOT=SIGTOT+2.*WGHFAC*
     *              0.25*(  1.-EXPEXH  +  1.-EXMEXH +
     *                      1.-EXNEXH  +  1.-EXOEXH  )
                SIGINE = SIGINE +  WGHFAC *
     *              0.25*(  1.-EXPEXP  +  1.-EXMEXP +
     *                      1.-EXNEXP  +  1.-EXOEXP  )
*               pure-soft-inelastic (hard scatterring is included as absorbtion)
                SIGSIN=SIGSIN+ WGHFAC *
     *              0.25*(    (EXPHAP-EXPEXP)
     *                       +(EXMHAP-EXMEXP)
     *                       +(EXNHAP-EXNEXP)
     *                       +(EXOHAP-EXOEXP) )
*               hard-inelastic (soft scatterring disregarded)
                SIGHIN=SIGHIN+ WGHFAC*
     *              0.25*(  1.-EXPHAP  +  1.-EXMHAP +
     *                      1.-EXNHAP  +  1.-EXOHAP  )
              ELSEIF(  NMAX.EQ.2  ) THEN
                SIGTOT=SIGTOT+2.*WGHFAC*
     *              0.25*(  1.-EXPEXH  +  1.-EXMEXH +
     *                      1.-EXNEXH  +  1.-EXOEXH
*                            subtract second order terms in .TRH and LOH:
C/////CHANGED - TO +
     *                      +(1.-EXPTRH)*(1-EXPLOH) *EXPSOH *EXPHAH
     *                      +(1.-EXMTRH)*(1-EXMLOH) *EXMSOH *EXMHAH
     *                      +(1.-EXNTRH)*(1-EXNLOH) *EXNSOH *EXNHAH
     *                      +(1.-EXOTRH)*(1-EXOLOH) *EXOSOH *EXOHAH )
                SIGINE = SIGINE +  WGHFAC *
     *              0.25*(  1.-EXPEXP  +  1.-EXMEXP +
     *                      1.-EXNEXP  +  1.-EXOEXP
*                            subtract second order terms in .TRP and LOP:
C/////CHANGED - TO +
     *                      +(1.-EXPTRP)*(1-EXPLOP) *EXPSOP *EXPHAP
     *                      +(1.-EXMTRP)*(1-EXMLOP) *EXMSOP *EXMHAP
     *                      +(1.-EXNTRP)*(1-EXNLOP) *EXNSOP *EXNHAP
     *                      +(1.-EXOTRP)*(1-EXOLOP) *EXOSOP *EXOHAP )
*               pure-soft-inelastic (hard scatterring is included as absorbtion)
                SIGSIN=SIGSIN+ WGHFAC *
     *              0.25*(    (EXPHAP-EXPEXP)
     *                       +(EXMHAP-EXMEXP)
     *                       +(EXNHAP-EXNEXP)
     *                       +(EXOHAP-EXOEXP)
*                            subtract second order terms of 2nd column:
     *                      +(1.-EXPTRP)*(1-EXPLOP) *EXPSOP *EXPHAP
     *                      +(1.-EXMTRP)*(1-EXMLOP) *EXMSOP *EXMHAP
     *                      +(1.-EXNTRP)*(1-EXNLOP) *EXNSOP *EXNHAP
     *                      +(1.-EXOTRP)*(1-EXOLOP) *EXOSOP *EXOHAP)
*               hard-inelastic (soft scatterring disregarded)
                SIGHIN=SIGHIN+ WGHFAC*
     *              0.25*(  1.-EXPHAP  +  1.-EXMHAP +
     *                      1.-EXNHAP  +  1.-EXOHAP  )
              ENDIF
*             high mass diffraction (sep.low mass diffr. arbitrary)
*                          (naive 1/EX?TRP -> EX?TR as selected cut counts -1)
              IF( NMAX.GE.3 ) THEN
                SIGHMD=SIGHMD + WGHFAC  *
     *                       0.25*( (EXPTRP-1.)*EXPEXP
     *                             +(EXMTRP-1.)*EXMEXP
     *                             +(EXNTRP-1.)*EXNEXP
     *                             +(EXOTRP-1.)*EXOEXP)
              ELSE
                SIGHMD=SIGHMD + WGHFAC  *
     *                       0.25*( EXPOTP * EXPSOP*EXPHAP
     *                             +EXPOTM * EXMSOP*EXMHAP
     *                             +EXPOTN * EXNSOP*EXNHAP
     *                             +EXPOTO * EXOSOP*EXOHAP )
              ENDIF
              IF( NMAX.GE.3  ) THEN
                SIHMDD=SIHMDD + WGHFAC  *
     *                       0.25*( (EXPLOP-1.)*EXPEXP
     *                             +(EXMLOP-1.)*EXMEXP
     *                             +(EXNLOP-1.)*EXNEXP
     *                             +(EXOLOP-1.)*EXOEXP)
              ELSEIF (NMAX.EQ.2 ) THEN
                SIHMDD=SIHMDD + WGHFAC  *
     *                       0.25*( EXPOLP * EXPSOP*EXPHAP
     *                             +EXPOLM * EXMSOP*EXMHAP
     *                             +EXPOLN * EXNSOP*EXNHAP
     *                             +EXPOLO * EXOSOP*EXOHAP )
*              no action:
*              ELSEIF (NMAX.LE.1 ) THEN
*                SIHMDD=SIHMDD + WGHFAC  *
*     *                       0.25*(   0.   * EXPSOP*EXPHAP
*     *                             +  0.   * EXMSOP*EXMHAP
*     *                             +  0.   * EXNSOP*EXNHAP
*     *                             +  0.   * EXOSOP*EXOHAP )
              ENDIF
            ENDIF
*           ending non L,M,N,NN depending part
*
735      CONTINUE
*         ending impact-integral loop
*
          IF(ABS(FAPHAR*FAPSOF*FAPTRP*FAPLOO*SIGP(L,M,NNN)).LT.DTINY)
     &    THEN
            SIGP(L,M,NNN)=0.
          ELSEIF(LSQRT) THEN
            SIGP(L,M,NNN)=FAPHAR*FAPSOF*FAPTRP*FAPLOO*SIGP(L,M,NNN)
     *         * abs(FAPHAR*FAPSOF*FAPTRP*FAPLOO)/DENOMI*PI4
          ELSEIF(.NOT.LSQRT) THEN
            SIGP(L,M,NNN)=FAPHAR*FAPSOF*FAPTRP*FAPLOO*SIGP(L,M,NNN)
     *         /DENOMI*PI4
          ENDIF
          IF(ABS(FAMHAR*FAMSOF*FAMTRP*FAMLOO*SIGM(L,M,NNN)).LT.DTINY)
     &    THEN
            SIGM(L,M,NNN)=0.
          ELSEIF(LSQRT) THEN
            SIGM(L,M,NNN)=FAMHAR*FAMSOF*FAMTRP*FAMLOO*SIGM(L,M,NNN)
     *           * abs( FAMHAR*FAMSOF*FAMTRP*FAMLOO)/DENOMI*PI4
          ELSEIF(.NOT.LSQRT) THEN
            SIGM(L,M,NNN)=FAMHAR*FAMSOF*FAMTRP*FAMLOO*SIGM(L,M,NNN)
     *         /DENOMI*PI4
          ENDIF
          IF(ABS(FANHAR*FANSOF*FANTRP*FANLOO*SIGN(L,M,NNN)).LT.DTINY)
     &    THEN
            SIGN(L,M,NNN)=0.
          ELSEIF(LSQRT) THEN
            SIGN(L,M,NNN)=FANHAR*FANSOF*FANTRP*FANLOO*SIGN(L,M,NNN)
     *           * abs( FANHAR*FANSOF*FANTRP*FANLOO)/DENOMI*PI4
          ELSEIF(.NOT.LSQRT) THEN
            SIGN(L,M,NNN)=FANHAR*FANSOF*FANTRP*FANLOO*SIGN(L,M,NNN)
     *         /DENOMI*PI4
          ENDIF
          IF(ABS(FAOHAR*FAOSOF*FAOTRP*FAOLOO*SIGO(L,M,NNN)).LT.DTINY)
     &    THEN
           SIGO(L,M,NNN)=0.
          ELSEIF(LSQRT) THEN
            SIGO(L,M,NNN)=FAOHAR*FAOSOF*FAOTRP*FAOLOO*SIGO(L,M,NNN)
     *          * abs( FAOHAR*FAOSOF*FAOTRP*FAOLOO/DENOMI)*PI4
          ELSEIF(.NOT.LSQRT) THEN
            SIGO(L,M,NNN)=FAOHAR*FAOSOF*FAOTRP*FAOLOO*SIGO(L,M,NNN)
     *          /DENOMI*PI4
          ENDIF
*
750        CONTINUE
740      CONTINUE
730    CONTINUE
720   CONTINUE
*
* *** summing up the three contributions
*
      NNNMAX=NMAXI+(NMAXI+1)*NNMAXI
      DO 820 L=0,LMAXI
        DO 830 M=0,MMAXI
         DO 830 NNN=0,NNNMAX
          SIG(L,M,NNN)=(SIGP(L,M,NNN)+SIGM(L,M,NNN)+
     *                  SIGN(L,M,NNN)+SIGO(L,M,NNN) )/4.
830     CONTINUE
820    CONTINUE
*
*
* *** calculate summed quantities for print out
*
      DO 3 L=0,LMAXI
        DO 4 M=0,MMAXI
         DO 4 N=0,NMAXI
          DO 4 NN=0,NNMAXI
            IF( NMAX.LE.2 .AND. N.EQ.1 .AND. NN.EQ.1 ) GO TO 4
            NNN=N+(NMAXI+1)*NN
            SIGSUM=SIGSUM + SIG(L,M,NNN)
*           for options outlawing hard without soft:
            IF(M.EQ.0.OR.L.GE.1) SIGSME=SIGSME + SIG(L,M,NNN)
            SHARD(M)=SHARD(M)+SIG(L,M,NNN)
            SSOFT(L)=SSOFT(L)+SIG(L,M,NNN)
            STRPL(N)=STRPL(N)+SIG(L,M,NNN)
            SIGINL = SIGINL + SIG(L,M,NNN)
            IF(M.GE.1) SIGHRI = SIGHRI + SIG(L,M,NNN)
            IF(L.EQ.0.AND.M.EQ.0.AND.NN.EQ.0.AND.N.GE.1) THEN
              SIGDI = SIGDI  + (-1)**N*SIG(L,M,NNN)
            ELSEIF(L.EQ.0.AND.M.EQ.0.AND.N.EQ.0.AND.NN.GE.1) THEN
              SIGDDI= SIGDDI + (-1)**NN*SIG(L,M,NNN)
            ENDIF
    4   CONTINUE
    3 CONTINUE
*   elastic processes were not generated, L,M,N,NN=0 no problem
*
      SIGLMD=SILMSD+SILMDD
      SITHMD=SIGHMD+SIHMDD
      SIGD = SIGLMD + SITHMD
      SLHMDD =  SQRT(ABS(SILMDD*SIHMDD))
      SIGDD= SILMDD + SIHMDD + SLHMDD
      SIGIN=SIGINE+SIGLMD
      SIGEL=SIGTOT-SIGIN
*
* ***  print out
*
      IF(LENTRY.EQ.1.AND.IOUTPO.LE.1) RETURN
*
      WRITE(6,*)' '
      WRITE(6,*)'  --- properties of events ---'
      WRITE (6,102)
      WRITE(6,*)'  Energy=',ECM
      WRITE (6,102)
      WRITE(6,*)'  max.contributing soft/hard/diffr./doubl.diffr. cuts'
      WRITE(6,*)'                     LMAXI=  MMAXI=  NMAXI=   NNMAXI='
      WRITE(6,'(15X,4I9)')              LMAXI,MMAXI,NMAXI,NNMAXI
      WRITE(6,*)'  methode used:  '
      WRITE(6,*)'                     ISIG=   ICON=   IPIM=     '
      WRITE(6,'(15X,3I9)')                     ISIG,ICON,IPIM
      WRITE (6,102)
      WRITE(6,*)'  --- bare cross section and eikonal constants ---'
C     COMMON/POMPAR/ALFA,ALFAP,A,C,AK
C     COMMON /SIGMA/SIGSOF,BS,ZSOF,SIGHAR,BH,ZHAR,SIGTRP,BT,ZTRP,
C    *              SIGLOO,ZLOO
      WRITE(6,*)'    ALFA =',ALFA,' ALFAP =',ALFAP,' A =',A
      WRITE(6,*)'    C =',C,' AK =',AK
      WRITE(6,*)'    ALALAM =',ALALAM
      WRITE (6,102)
      WRITE(6,*)'     SIGSOF=',SIGSOF,'  BS=',BS,'  ZSOF=',ZSOF
      WRITE(6,*)'     SIGHAR=',SIGHAR,'  BH=',BH,'  ZHAR=',ZHAR
      WRITE(6,*)'     SIGTRP=',SIGTRP,'  BT=',BT,'  ZTRP=',ZTRP
      WRITE(6,*)'     SIGLOO=',SIGLOO,'  BT=',BT,'  ZLOO=',ZLOO
      WRITE (6,102)
      WRITE(6,*)'  --- observable cross sections ---'
      WRITE (6,102)
      WRITE(6,*)'     TOTAL X-SECTION         = ',SIGTOT
      WRITE(6,*)'     ELASTIC X-SECTION       = ',SIGELE
      WRITE(6,*)'     INELASTIC X-SECTION-LMD = ',SIGINE
      WRITE(6,*)'     INELASTIC X-SECTION     = ',SIGIN
      WRITE(6,*)'     HARD INEL. X-SECTION    = ',SIGHIN
      WRITE (6,102)
      WRITE(6,*)'  LOW MASS SING./DOUB.DIFFR.X-SECTION= ',SILMSD,SILMDD
      WRITE(6,*)'  => LOW MASS TOTAL DIFFRACTIV.X-SECTION=     ',SIGLMD
      WRITE(6,*)'  HIGH MASS SING./DOUB.DIFFR.X-SECTION= ',SIGDI,SIGDDI
      WRITE(6,*)'  => HIGH MASS TOTAL DIFFRACTIV.X-SECTION=    ',SITHMD
      WRITE(6,*)'  ESTIMAT.MIXED (LM+HM) DOUBL.DIFFRAC.X.SEC.= ',SLHMDD
      WRITE(6,*)'  => '
      WRITE(6,*)'     DIFFRACTIVE  X-SECTION =    ',SIGD
      WRITE(6,*)'     DOUBLY DIFFRACTIVE X-SECT. =',SIGDD
      WRITE (6,102)
*
      IF(IOUTPO.GE.0) THEN
       WRITE(6,*)'  --- observ. x-sections, altern. calculated ---'
       WRITE(6,*)'     ELASTIC X-SECTION   = ',SIGEL
       WRITE(6,*)'     INELASTIC X-SECTION-LMD = ',SIGINL
       WRITE(6,*)'     HARD INEL. X-SECTION= ',SIGHRI
       WRITE(6,*)'  HIGH MASS SING./DOUB.DIFFR.X-SECT.=',SIGHMD,SIHMDD
      WRITE(6,*)'     X-SECTION FOR (L,M,N,NN)= 1000 0100 0010 0001'
      WRITE(6,*)'            ',SIG(1,0,0),SIG(0,1,0)
     *                        ,SIG(0,0,1),SIG(0,0,2)
       WRITE (6,102)
      ENDIF
*
      IF(IOUTPO.GE.2) THEN
         WRITE (6,102)
         NNMAXP=NMAXI/2
         IF( NMAXI.LT.2)NNMAXP=1
         DO 52 N=0,NNMAXP
*          printout loops:
           DO 48 L=0,LMAXI
 48          WRITE(6,101)(SIG(L,M,N),M=0,7)
           WRITE (6,102)
           DO 50 L=0,LMAXI
 50          WRITE(6,101)(SIG(L,M,N),M=8,15)
           WRITE (6,102)
           WRITE(6,*)
     &       '  # CUT-POMERON  SSOFT X-SECT.  SHARD X-SECT.'
           DO 58 L=0,LMAXI
 58          WRITE (6,103)L,SSOFT(L),SHARD(L)
           WRITE (6,102)
*        printoutloop ends
52       CONTINUE
      ENDIF
*
* *** attribute x-sections (SIG) for CUT objects ('s and PHI's)
*                            to string configurations (PLMN) s:  **********
*
C                               CHANGED 10.1.90 BY J.R.
*     preparations:
C     just for Y and (Phi) - cuts:
      FAK(0)=1
      DO 500 I=1,NMAXI
        FAK(I)=FAK(I-1)*I
  500 CONTINUE
      DO 501 I=0,NMAXI
        DO 501 J=0,I
          CMBIN(I,J)=FAK(I)/(FAK(J)*FAK(I-J))
  501 CONTINUE
*
      TMMP=0.
      DO 5 L=0,LMAXI
      DO 5 M=0,MMAXI
       IF(ICON.EQ.44.OR.ICON.EQ.46.OR.ICON.EQ.48.
     *                             OR.ICON.EQ.54) THEN
C///test:
*        no Y or PHI cut
           PLMNTM=SIG(L,M,0)/(SIGSUM+TIN)
           PLMN(L,M,0) =  PLMNTM + PLMN(L,M,0)
           TMMP=TMMP+PLMNTM
*        Y but no PHI cut
           PLMNTM=SIG(L,M,1)/(SIGSUM+TIN)
           TMMP=TMMP+PLMNTM
           IF(L+2.LE.LMAXI) THEN
               PLMN(L+2,M,0) = (-2.)* PLMNTM + PLMN(L+2,M,0)
               PLMN(L+1,M,0) =  4.  * PLMNTM + PLMN(L+1,M,0)
           ELSE
               PLMN(LMAXI,M,0) = (-2.)* PLMNTM + PLMN(LMAXI,M,0)
               PLMN(LMAXI,M,0) =  4.  * PLMNTM + PLMN(LMAXI,M,0)
           ENDIF
           IF(L.EQ.0 .AND. M.EQ.0) THEN
             PLMN(L  ,M,1) = (-1.)* PLMNTM + PLMN(L  ,M,1)
           ELSE
             PLMN(L  ,M,0) = (-1.)* PLMNTM + PLMN(L  ,M,0)
           ENDIF
*        no Y but PHI cut
           PLMNTM=SIG(L,M,2)/(SIGSUM+TIN)
           TMMP=TMMP+PLMNTM
           IF(L+2.LE.LMAXI) THEN
               PLMN(L+2,M,0) = (-2.)* PLMNTM + PLMN(L+2,M,0)
               PLMN(L+1,M,0) =  4.  * PLMNTM + PLMN(L+1,M,0)
           ELSE
               PLMN(LMAXI,M,0) = (-2.)* PLMNTM + PLMN(LMAXI,M,0)
               PLMN(LMAXI,M,0) =  4.  * PLMNTM + PLMN(LMAXI,M,0)
           ENDIF
           IF(L.EQ.0 .AND. M.EQ.0) THEN
             PLMN(L  ,M,2) = (-1.)* PLMNTM + PLMN(L  ,M,2)
           ELSE
             PLMN(L  ,M,0) = (-1.)* PLMNTM + PLMN(L  ,M,0)
           ENDIF
C/// test end
       ELSE
        DO 51  N=0,NMAXI
         DO 51 NN=0,NNMAXI
           IF(NMAX.LE.2 .AND. N.EQ.1 .AND. NN.EQ.1) GO TO 51
           NNN=N+(NMAXI+1)*NN
*
*          to be attributed::
           PLMNTM=SIG(L,M,NNN)/(SIGSUM+TIN)
           TMMP=TMMP+PLMNTM
*
*          attribution loop for Y-cuts
           DO 511  N0CUT=0,N
           DO 511  N1CUT=0,N-N0CUT
                 N2CUT=N-N0CUT-N1CUT
*            combinatoric weight:
             CMB0=CMBIN(N,N2CUT)
             CMB1=CMBIN(N-N2CUT,N1CUT)
*
*          attribution loop for PHI-cuts
           DO 511  NN0CUT=0,NN
           DO 511  NN1CUT=0,NN-NN0CUT
                 NN2CUT=NN-NN0CUT-NN1CUT
*            combinatoric weight:
             CMBN0=CMBIN(NN,NN2CUT)
             CMBN1=CMBIN(NN-NN2CUT,NN1CUT)
*
*            attributions matrix:
*            ("L"soft,"M"hard,"N"1-diffr.,"NN"2-dif.,"NL"diffr.long in.part ):
*            obviously:
               L2STR = L
               M2STR = M
               N2STR = N0CUT
               NN2STR= NN0CUT
               NL2STR=  0
*            specialties for Y's and PHI's:
               L2STR=L2STR + N1CUT + NN1CUT +  N2CUT + NN2CUT
               IF(NMAX.LE.2)THEN
*                room to have NL2STR's:
                 NL2STR= N2CUT + NN2CUT
               ELSEIF(NMAX.GE.3)THEN
*                the extra inner piece counts here like long strings:
                 L2STR=L2STR+N2CUT+NN2CUT
               ENDIF
               IF((ICON.EQ.26.OR.ICON.EQ.36.OR.ICON.EQ.46.OR.ICON.EQ.56)
     &           .AND. (L2STR.GE.1.OR.M2STR.GE.1))THEN
                 L2STR=L2STR +  NL2STR
                 N2STR  = 0
                 NN2STR  = 0
                 NL2STR  = 0
               ENDIF
*
*            getting and checking parameter for storing:
               IF(L2STR.GT.LMAXI) L2STR=LMAXI
               IF(M2STR.GT.LMAXI) M2STR=LMAXI
               NNNSTR =N2STR +(NMAXI+1)*NN2STR
     *                       +(NNMAXI+1)*(NMAXI+1)*NL2STR
               IF(NNNSTR.GT.MXPA13) NNNSTR=MXPA13
*
*            summing contributions
               PLMN(L2STR,M2STR,NNNSTR) =  PLMNTM
     *             *CMB0*CMB1 * (-2)**N2CUT * (4)**N1CUT * (-1)**N0CUT
     *             *CMBN0*CMBN1*(-2)**NN2CUT* (4)**NN1CUT* (-1)**NN0CUT
     &           +  PLMN(L2STR,M2STR,NNNSTR)
*
  511   CONTINUE
  51    CONTINUE
      ENDIF
C///// initial general methode ends
    5 CONTINUE
       IF(ABS(TMMP-1.D0).GT..03D0)THEN
          WRITE(6,*)
     &     ' NORMALISATION ERROR SUM PLM before LMD reatribution=',TMMP
       ENDIF
*
* *** built in low mass diffraction, get averages and check normalisation:
*
*       low mass diffraction was SIGLMD=SILMSD+SILMDD
*       and mixed LM/HM diffraction SLHMDD=SQRT(SILMDD*SIHMDD)
        PLMFAC= (SIGSUM+TIN) / (SIGSUM+TIN +SIGLMD)
        PLMN(0,0,1)= PLMN(0,0,1) +
     &      ( SILMSD - SLHMDD ) / (SIGSUM+TIN)
        PLMN(0,0,2)= PLMN(0,0,2) +
     &      ( SILMDD + SLHMDD ) / (SIGSUM+TIN)
 661    CONTINUE
*     AVerage_SOft_N,AVerage_HaRD_N,SUM_over_Pl
      AVSOFN=0.
      AVHARN=0.
      AVDIFN=0.
      AVDDFN=0.
      AVDLFN=0.
      PSOFT=0.
      TEMP=0.
      TMP=0.
      TMMP=0.
      TMMP1=0.
*
*     (L,M,N,NN,NL repl. L2STR,M2... for s.,h.,Y-dif.,PHI-dif.,i.Y-dif.str.#)
      DO 6 NL=0,NLMAXI
      DO 6 NN=0,NNMAXI
      DO 6 N=0,NMAXI
        IF(NMAX.LE.2 .AND. N+NN+NL.GE.2) GO TO 6
        NNN =N +(NMAXI+1)*NN +(NNMAXI+1)*(NMAXI+1)* NL
        DO 63 M=0,MMAXI
        DO 63 L=0,LMAXI
           IF(NL.EQ.0)TMMP1  = TMMP1  + SIG(L,M,NNN)
           TMMP  = TMMP  + SIG(L,M,NNN)
           PLMN(L,M,NNN)=PLMN(L,M,NNN) * PLMFAC
           TMP  = TMP  + PLMN(L,M,NNN)
C          IF(PLMN(L,M,NNN).LT.-.000001D0)
           IF(PLMN(L,M,NNN).LT.-.000005D0)
     &             WRITE(6,*)' 0>PLMN',PLMN(L,M,NNN),L,M,N,NN,NL
           AVSOFN=AVSOFN+PLMN(L,M,NNN)*L
           AVHARN=AVHARN+PLMN(L,M,NNN)*M
           AVDIFN=AVDIFN+PLMN(L,M,NNN)*N
           AVDDFN=AVDDFN+PLMN(L,M,NNN)*NN
           AVDLFN=AVDLFN+PLMN(L,M,NNN)*NL
           IF (M.EQ.0)PSOFT=PSOFT+PLMN(L,M,NNN)
 63    CONTINUE
  6    CONTINUE
       IF(ABS(TMP-1.D0).GT..01D0)THEN
          WRITE(6,*)
     &     ' NORMALISATION ERROR SUM PLM before M reatribution=',TMP
       ENDIF
       TMMP=TMMP/SIGSUM
       TMMP1=TMMP1/SIGSUM
       IF(ABS(TMMP-1.D0).GT..01D0 .OR.ABS(TMMP1-1.D0).GT..01D0)THEN
          WRITE(6,*)
     &     ' NORMALISATION ERROR TMMP,TMMP1=',TMMP,TMMP1
       ENDIF
*
* *** reattribute purely hard scattering and get cummulant distribution
*
*     (L,M,N,NN,NL repl. L2STR,M2...
*          for soft,hard,Y-diffr.,PHI-diffr.,inner diffr.  string #)
      DO 61 NL=0,NLMAXI
      DO 61 NN=0,NNMAXI
      DO 61 N=0,NMAXI
        IF(NMAX.LE.2 .AND. N+NN+NL.GE.2) GO TO 61
        NNN =N +(NMAXI+1)*NN +(NNMAXI+1)*(NMAXI+1)* NL
        DO 612 M=0,MMAXI
        DO 611 L=0,LMAXI
* -- hard:  a pure hard scattering gets an extra soft chain, it is considered
*           a specialty of fragmentation and therefor implemented only here
*           there are are number of other options which are dropped as they
*           are hard to implement at this point
           IF (L.EQ.0.AND.M.GE.1)THEN
             PLMN(1,M,NNN)=PLMN(1,M,NNN)+PLMN(0,M,NNN)
             PLMN(0,M,NNN)=0.
           ENDIF
* -- cummulant of distribution
           TEMP  = TEMP  + PLMN(L,M,NNN)
           PLMNCU(L,M,NNN)=TEMP
  611   CONTINUE
*
          IF(IOUTPO.GE.3)WRITE (6,*)' M,(L,PLMN(L,M,N),L=0,LMAX)'
          IF(IOUTPO.GE.3)WRITE (6,106) M,(L,PLMN(L,M,N),L=0,LMAXI)
          IF(IOUTPO.GE.2)WRITE (6,*)' M,(L,PLMNCU(L,M,N),L=0,LMAX/2)'
          IF(IOUTPO.GE.2)WRITE (6,106) M,(L,PLMNCU(L,M,N),L=0,LMAXI/2)
  106     FORMAT (I3,9(I3,E11.2))
*
  612  CONTINUE
   61  CONTINUE
*
       IF(ABS(TEMP-1.D0).GT..01D0)THEN
          WRITE(6,*)' NORMALISATION ERROR SUM PLM=',TEMP
          PLMFAC=1./(TEMP+TIN)
          GO TO 661
       ENDIF
*
         IF(IOUTPO.GE.1)WRITE (6,*)
     &   '(((L,M,N,PLMN(L,M,N),N=0,2),M=0,5),L=0,7)'
         IF(IOUTPO.GE.1)WRITE (6,1106)
     &   (((L,M,N,PLMN(L,M,N),N=0,2),M=0,5),L=0,7)
         IF(IOUTPO.GE.1)WRITE (6,*)
     &   '(((L,M,N,SIG(L,M,N),N=0,2),M=0,5),L=0,7)'
         IF(IOUTPO.GE.1)WRITE (6,1106)
     &   (((L,M,N,SIG(L,M,N),N=0,2),M=0,5),L=0,7)
 1106    FORMAT (1X,3(I5,I5,I5,G12.5))
*
      PHARD=1.-PSOFT
      ALFAH=SIGHIN/(SIGINE+0.00001)
      BETAH=1.-ALFAH
      WRITE(6,116)AVSOFN,AVHARN,AVDIFN,AVDDFN,AVDLFN,
     &        PHARD,PSOFT,ALFAH,BETAH
  116 FORMAT(/'--- various averages:'/
     &       /'    AVSOFN=    AVHARN=    AVDIFN=    AVDDFN=    AVDLFN='
     &       /'   ',5F11.3
     &       /'    PHARD=     PSOFT=     ALFAH=     BETAH= '
     &       /'   ',4F11.3)
      IF(IOUTPO.GE.1)WRITE(6,*)'SIGSUM=SIGINL-LMD',SIGSUM
*
        IF(IOUTPO.GE.1)WRITE(6,610) SIGTOT,SIGINE,SIGD,SIGDD,SIGHIN
  610    FORMAT (' SIGTOT,SIGINE,SIGD,SIGDD,SIGHIN= '/' ',5E18.6)
*
101   FORMAT(' ',10E10.3)
102   FORMAT(' ')
103   FORMAT(' ',5X,I4,5X,2E15.3)
*
      RETURN
      END
*     ende problm
*
************************************************************************
*
*
      SUBROUTINE SAMPLX(L2STR,M2STR,N2STR,NN2STR,NL2STR)
*
*     input:
*        PLMNCU
*     output:
*        samples number of soft (L) and hard (M) cut pomerons from PLMNC
*                   and of (N=0/1/2) diffractive excitations (for L=M=0
*
*----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER (MXPA25=30,MXPA26=MXPA25+1,MXPA13=13)
*     PARAMETRIZATION FOR PTMIN= 3. GEV
      PARAMETER (MXPA50=250,MXPA51=MXPA50+1)
*     PARAMETRIZATION FOR PTMIN= 2. GEV
C     PARAMETER (MXPA50=350,MXPA51=MXPA50+1)
* *** /OUTLEV/ controls output level for POMDI and parton X distribution
      COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
      COMMON /POMTYP/IPIM,ICON,ISIG,LMAX,MMAX,NMAX,DIFEL,DIFNU
* *** /POLMN/ arrays having to do with cut soft and hard Pomerons
      COMMON /POLMN/PLMN(0:MXPA25,0:MXPA50,0:MXPA13),
     *              PLMNCU(0:MXPA25,0:MXPA50,0:MXPA13)
      COMMON /POLMN0/PDIFR,PHARD,PSOFT,ALFAH,BETAH,
     *              SIGTOT,SIGQEL,SIGEL,SIGINE,SIGHIN,SIGD,SIGDD
*
      PARAMETER (PI=3.141592654D0)
      DATA NPRINT/0/
*     "L"SOFT,"M"HARD
      LMAXI = LMAX
      MMAXI = MMAX
      IF(IPIM.NE.2) THEN
*       "N"diffr.,"NN"dou.dif.,("NL"long inner contr. for 2*cut Y or PHI)
        NMAXI = NMAX
        NNMAXI=0
        NLMAXI=0
      ELSEIF(IPIM.EQ.2) THEN
        IF( NMAX.GE.3)THEN
          NMAXI = NMAX
          NNMAXI=(13-NMAXI)/(1+NMAXI)
*         13 =!= NNNMAX =  NMAXI+(NMAXI+1)*NNMAXI
          NLMAXI=0
        ELSEIF( NMAX.EQ.2)THEN
          NMAXI=1
          NNMAXI=1
          NLMAXI=1
        ELSEIF( NMAX.EQ.1)THEN
          NMAXI=1
          NNMAXI=0
          NLMAXI=1
        ENDIF
      ENDIF
  111 CONTINUE
*
      X=RNDM(V)
*
      IF (X.LE.PLMNCU(0,0,0) .AND. NPRINT.LT.100)THEN
        WRITE(6,*) ' No generator of elastic events '
        WRITE(6,*) ' PLMNCU (0,0,0)  =!= 0 = ',PLMNCU(0,0,0)
        NPRINT=NPRINT+1
        GOTO 111
      ENDIF
*
      DO 5 NL=0,NLMAXI
      DO 5 NN=0,NNMAXI
      DO 5 N=0,NMAXI
        NNN =N +(NMAXI+1)*NN +(NNMAXI+1)*(NMAXI+1)* NL
        DO 6 M=0,MMAXI
        DO 7 L=0,LMAXI
*
          IF (X.LE.PLMNCU(L,M,NNN)) THEN
            L2STR=L
            M2STR=M
            N2STR=N
            NN2STR=NN
            NL2STR=NL
            RETURN
*
          ENDIF
    7    CONTINUE
    6  CONTINUE
    5 CONTINUE
*
      NPRINT=NPRINT+1
      IF(NPRINT.LT.100)  WRITE(6,*)' RAR.IN SAMPLM,PLMNCU,RND=',
     &  PLMNCU(LMAX, MMAX,NNN),X,NPRINT
      IF( PLMNCU(LMAX,MMAX,NNN) .GT. 0.1D0 ) RETURN
      IF( PLMNCU(LMAX,0,0) .GT. 0.1D0 ) RETURN
      WRITE(6,*)' RAR.IN SAMPLM- PROBLEM SEEMS BAD, DECIDE TO STOP'
      STOP
      END
*
*
************************************************************************
*
      SUBROUTINE SAMPLM(L2STR,M2STR,N2STR)
*
*     input:
*        PLMNCU
*     output:
*        samples number of soft (L) and hard (M) cut pomerons from PLMNC
*                   and of (N=0/1/2) diffractive excitations (for L=M=0
*
*----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
* *** /OUTLEV/ controls output level for POMDI and parton X distribution
      COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
      COMMON /POMTYP/IPIM,ICON,ISIG,LMAX,MMAX,NMAX,DIFEL,DIFNU
* *** /POLMN/ arrays having to do with cut soft and hard Pomerons
      PARAMETER (MXPA25=30,MXPA26=MXPA25+1,MXPA13=13)
*     PARAMETRIZATION FOR PTMIN= 3. GEV
      PARAMETER (MXPA50=250,MXPA51=MXPA50+1)
*     PARAMETRIZATION FOR PTMIN= 2. GEV
C     PARAMETER (MXPA50=350,MXPA51=MXPA50+1)
      COMMON /POLMN/PLMN(0:MXPA25,0:MXPA50,0:MXPA13),
     *              PLMNCU(0:MXPA25,0:MXPA50,0:MXPA13)
      COMMON /POLMN0/PDIFR,PHARD,PSOFT,ALFAH,BETAH,
     *              SIGTOT,SIGQEL,SIGEL,SIGINE,SIGHIN,SIGD,SIGDD
*
      PARAMETER (PI=3.141592654D0)
  111 CONTINUE
*
      X=RNDM(V)
*
      IF (X.LE.PLMNCU(0,0,0))THEN
        WRITE(6,*) ' No generator of elastic events '
        WRITE(6,*) ' PLMNCU (0,0,0)  =!= 0 = ',PLMNCU(0,0,0)
        GOTO 111
      ENDIF
*
      DO 5 N=0,NMAX
C      IF (N.GT.1) GO TO 111
       DO 6 M=0,MMAX
         DO 7 L=0,LMAX
*
          IF (X.LE.PLMNCU(L,M,N)) THEN
            L2STR=L
            M2STR=M
            N2STR=N
            RETURN
*
          ENDIF
    7    CONTINUE
    6  CONTINUE
    5 CONTINUE
*
      WRITE(6,*)' RAR.IN SAMPLM,PLMNCU,RND=',PLMNCU(LMAX,MMAX,NMAX),X
      IF( PLMNCU(LMAX,MMAX,NMAX) .GT. 0.1D0 ) RETURN
      IF( PLMNCU(LMAX,0,0) .GT. 0.1D0 ) RETURN
      WRITE(6,*)' RAR.IN SAMPLM- PROBLEM SEEMS BAD, DECIDE TO STOP'
      STOP
      END
*
*
*
*
************************************************************************
C--------------------------------------------------------------------
C
C                      dtUpom9h.for
C
C---------------------------------------------------------------------
************************************************************************
*
*         POMDI,SIGMAS,SIGSHD,PRBLM0..9,SAMPLM,SIGMA1
*
*                  routines called by code word
*                            SIGMAPOM :
*

C--------------------------------------------------------------------
C
C                      dtUpom9h.for
C
C---------------------------------------------------------------------
************************************************************************
*
*         POMDI,SIGMAS,SIGSHD,PRBLM0..9,SAMPLX,SIGMA1
*
*                  routines called by code word
*                            SIGMAPOM :
*
************************************************************************
*
*                                                 J.RANFT September 1987
      SUBROUTINE POMDI
*
*        to calculate the s-dependent  X-sections
*
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*
*     input:
*       ISIG  characterizing X sections, transmitted to SIGSHD
*       IPIM  characterizes the method to calculate SIGMA(LSOFT,mhard)
*                        IPIM=1 : integral method with low mass dif.matr
*                        IPIM=2 : int.meth.with low mass dif.matr.(2*2)
*                                 and reduced high mass diffraction
*                        IPIM=3 : integral methode
*                        IPIM=4 : integral methode with Y -cuts
*                        IPIM=5 : integral methode with Y-cuts + 2 CHANNEL EIK.
*                        rest   : not implemented or
*                                 special cases for checking
*       LMAX < MXPA25 maximal number of considered soft pomerons
*       MMAX < MXPA50 maximal number of considered hard pomerons
*       NMAX < MXPA13 maximal number of considered trippel pomerons
*                 (not used in low masss diffraction formalism)
*     output: printet and plottet in SIGMAS and
*             printed on the end of this subroutine
*
*        card XSECTION causes call to POMDI
*        card SIGMAPOM with ITEST=1 causes call to POMDI
*        POMDI calling
*           SIGMAS (  SIGMA1..3( SIGSD, (1:/3:)GSET ), PLOT)
*           PRBLM..  ( SIGSHD,  (1:)GSET )
*           (2:)SAMPLX, (1:/3:)SAMPLM
*
*----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
*
* *** /POLMN/ arrays having to do with cut soft and hard Pomerons
      PARAMETER (MXPA25=30,MXPA26=MXPA25+1,MXPA13=13)
*     PARAMETRIZATION FOR PTMIN= 3. GEV
      PARAMETER (MXPA50=250,MXPA51=MXPA50+1)
*     PARAMETRIZATION FOR PTMIN= 2. GEV
C     PARAMETER (MXPA50=350,MXPA51=MXPA50+1)
      COMMON /POLMN/PLMN(0:MXPA25,0:MXPA50,0:MXPA13),
     *              PLMNCU(0:MXPA25,0:MXPA50,0:MXPA13)
      COMMON /POLMN0/PDIFR,PHARD,PSOFT,ALFAH,BETAH,
     *              SIGTOT,SIGQEL,SIGEL,SIGINE,SIGHIN,SIGD,SIGDD
* *** /HISTOO/
      INTEGER*2 NDISLM
      COMMON /HISTOO/AS(50,9),AECM(50,9),ASIG(50,9),ALOS(50,9),
     *    ALOECM(50,9),NDISLM(0:MXPA25,0:MXPA50,0:MXPA13)
* *** /OUTLEV/ controls output level for POMDI and parton X distribution
      COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
* --- used only in SIGMAPOM-routines
* *** /POMPAR/ contains pomeron parameters used in some options
*     ALFA,ALFAP,A,BH,C,BS are soft Pomeron parameters chosen in SIGMAS
*     appearing in SIGMAS, PRBLM..,  AK is K-factor
      COMMON/POMPAR/ALFA,ALFAP,A,C,AK
* *** /SIGMA/ contains variables actually used in iteration
*     ZSOF, ZHAR, BS, BH, SIGSOF,SIGHAR ( for soft, hard, trippel-Pom.)
*     are input for X-section calculated in SIGSHD
C     (/sigma/ is out of P.A.'s CM88, containing variables used in iter.
      COMMON /SIGMA/SIGSOF,BS,ZSOF,SIGHAR,BH,ZHAR,SIGTRP,BT,ZTRP,
     *              SIGLOO,ZLOO
* *** /POMTYP/ contains parameters determining X-sections
*     IPIM,ISIG,LMAX,MMAX,NMAX  as described at "CODEWD = SIGMAPOM"
C     LMAX,MMAX,NMAX kept open to enable avoiding of numerical problem
      COMMON /POMTYP/IPIM,ICON,ISIG,LMAX,MMAX,NMAX,DIFEL,DIFNU
      COMMON /ALALA/ALALAM
      COMMON/COLLIS/SS,IJPROJ,IJTAR,PTTHR,PTTHR2,IOPHRD,IJPRLU,IJTALU
C      ECM calculated as SQRT
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      PARAMETER (PI=3.141592654D0)
      LMAXI = LMAX
      MMAXI = MMAX
      NMAXI = NMAX
      SKEEP=SS
*
*---------------------------------------------------------------------
*                                                 CALL SIGMAS
*
*  calculates and plotts out X-section at various energies
*  runs thru energies independent of the rest
*                       and contains no preperation for other work
*
          WRITE(6,'(1X/1X)')
          WRITE(6,*)' '
          WRITE(6,*)
     *' ------ testing the energy dependence of x-sections ----------'
          WRITE(6,*)' '
          IF(IOUTPO.GT.-1) WRITE(6,*)
     *'          (as function of ALAM i.e.a low mass diffr.parameter)'
          WRITE(6,*)' -----------------------------------------------'
          WRITE(6,'(1X)')
*
          DO 1007 IIJJ=1,10
          IF(IOUTPO.GT.-1 .OR. IIJJ.EQ.6)THEN
C         WITHOUT SIMGMAS FOR NORMAL USERS
C         IF(IOUTPO.GT.-1 )THEN
            ALALAM=IIJJ*0.1
            IF(IOUTPO.GT.-1) WRITE(6,1008)ALALAM
 1008       FORMAT (' ALAM= ',F10.3)
            CALL SIGMAS
          ENDIF
 1007     CONTINUE
*
*---------------------------------------------------------------------
*                                                  test sampling L,M
*
 1111   CONTINUE
*                            looping thru the energies
*
C       DO 100 III=3,9
C         S=10.**III
C         ECM=SQRT(S)
*
          S=SKEEP
          ECM=SQRT(S)
*
*                            chosing one options with a standard call
*
          IF(IPIM.EQ.2) THEN
            IF( NMAX.GE.3)THEN
              NNMAXI=(13-NMAXI)/(1+NMAXI)
*             13 =!= NNNMAX =  NMAXI+(NMAXI+1)*NNMAXI
              NLMAXI=0
            ELSEIF( NMAX.EQ.2)THEN
              NMAXI=1
              NNMAXI=1
              NLMAXI=1
            ELSEIF( NMAX.EQ.1)THEN
              NMAXI=1
              NNMAXI=0
              NLMAXI=1
            ENDIF
            CALL PRBLM2(ECM)
          ENDIF
          IF(IPIM.LT.1.AND.IPIM.GT.9)THEN
                 WRITE(6,*) 'RETURN caused by IPIM=',IPIM
                 RETURN
          ENDIF
*
*                            randomly sample events and printout
*
          WRITE (6,'(1X)' )
          WRITE (6,102)ECM,S
  102     FORMAT
     *    ('--- sample distribution for L soft and M hard inelastic'
     *   , ' pomerons (string pairs)--- '
     *    / 20X,'at ECM  = ',F10.2,' S  = ',F12.1)
          DO 31 L=0,LMAXI
            DO 32 M=0,MMAXI
            DO 32 N=0,13
              NDISLM(L,M,N)=0
   32       CONTINUE
   31     CONTINUE
*
          IF(ICON.EQ.12)GO TO 100
          DO 320 II=1,10000
            IF(IPIM.EQ.2) THEN
              CALL SAMPLX(L2STR,M2STR,N2STR,NN2STR,NL2STR)
                 NNNSTR =N2STR +(NMAXI+1)*NN2STR
     *                         +(NNMAXI+1)*(NMAXI+1)*NL2STR
              NDISLM(L2STR,M2STR,NNNSTR)=NDISLM(L2STR,M2STR,NNNSTR)+1
            ELSE
              CALL SAMPLM(L2STR,M2STR,N2STR)
              NDISLM(L2STR,M2STR,N2STR)=NDISLM(L2STR,M2STR,N2STR)+1
            ENDIF
  320     CONTINUE
*
            WRITE(6,*)
     *      '                    with no diffractive contribution'
            WRITE(6,*) ' '
            WRITE(6,*)
     *   '     ....... vertical: NSTR, horizontal MSTR .........'
            DO 3344 L=0,MIN(20,LMAXI)
 3344          WRITE(6,34)L,(NDISLM(L,M,0),M=0,20)
            WRITE(6,*) ' '
          IF(IOUTPO.GE.0)THEN
            DO 333 N=0,5
              IF(N.NE.0)THEN
               WRITE(6,*)'               WITH NSTR=',N
               DO 334 L=0,MIN(20,LMAXI)
                 WRITE(6,34)L,(NDISLM(L,M,N),M=0,20)
  334          CONTINUE
               WRITE(6,*) ' '
              ENDIF
              JMPA50 = INT(MXPA50/25)
C             WRITE(6,*) 'WIDE PLOT 0<L<25, 0<M<'
              WRITE(6,*) 'WIDE PLOT 0<L<',MXPA25,' 0<M<'
     &                   ,MXPA50,' IN STEPS OF ',JMPA50
C             DO 335 L=0,MIN(25,LMAXI)
              DO 335 L=0,MXPA25
                WRITE(6,35)L,(NDISLM(L,M,N),M=0,MXPA50,JMPA50)
  335         CONTINUE
              WRITE(6,*) ' '
  333       CONTINUE
          ENDIF
   34         FORMAT (I5,':',21I4)
   35         FORMAT (I5,26I4)
  100   CONTINUE
*     energy loop has ended
      RETURN
      END
*
*
************************************************************************
*
      SUBROUTINE SIGMAS
*
*     output:
*       energy dependence of
*       SIGMa-TOT, SIGma-INel, SIGma-Diffractive, SIGma-Hard-INelastic
*       print- and plott-out
*     using:
*       called routines  SIGSHD , SIGMA1..3 , PLOT
*
*     runs thru energies independent of the rest
*                            and contains no preperation for other parts
*
*---------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
*
* *** /OUTLEV/ controls output level for POMDI and parton X distribution
      COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
* *** /POLMN/ arrays having to do with cut soft and hard Pomerons
      PARAMETER (MXPA25=30,MXPA26=MXPA25+1,MXPA13=13)
      PARAMETER ( ZERO=0.D0, ONE=1.D0)
*     PARAMETRIZATION FOR PTMIN= 3. GEV
      PARAMETER (MXPA50=250,MXPA51=MXPA50+1)
*     PARAMETRIZATION FOR PTMIN= 2. GEV
C     PARAMETER (MXPA50=350,MXPA51=MXPA50+1)
      COMMON /POLMN/PLMN(0:MXPA25,0:MXPA50,0:MXPA13),
     *              PLMNCU(0:MXPA25,0:MXPA50,0:MXPA13)
      COMMON /POLMN0/PDIFR,PHARD,PSOFT,ALFAH,BETAH,
     *              SIGTOT,SIGQEL,SIGEL,SIGINE,SIGHIN,SIGD,SIGDD
*
* *** /POMPAR/*/SIGMA/*/POMTYP/ used only in SIGMAPOM-routines (->POMDI)
      COMMON /POMTYP/IPIM,ICON,ISIG,LMAX,MMAX,NMAX,DIFEL,DIFNU
      COMMON/POMPAR/ALFA,ALFAP,A,C,AK
      COMMON /SIGMA/SIGSOF,BS,ZSOF,SIGHAR,BH,ZHAR,SIGTRP,BT,ZTRP,
     *              SIGLOO,ZLOO
* ***
      COMMON /TOPDR/ITOPD,IDUMTP
* ***
      INTEGER*2 NDISLM
      COMMON /HISTOO/AS(50,9),AECM(50,9),ASIG(50,9),ALOS(50,9),
     *    ALOECM(50,9),NDISLM(0:MXPA25,0:MXPA50,0:MXPA13)
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      PARAMETER (PI=3.141592654D0)
*
*---------------------------------------------------------------------
*                  run thru energies
*    ------------------------------------------------------------------
*
      ISTEP=2
      IF(IOUTPO.GT.-1)ISTEP=7
      DO 100 I=1,50,ISTEP
        S=1.6**I
        ECM=SQRT(S)+3.4
        S=ECM**2
*
* *** running thru energies initializing AS,AECM,ASIG,ALOS,ALOECM *****
*
        DO 111 III=1,9
          AS(I,III)=S
          AECM(I,III)=ECM
          ALOS(I,III)=LOG10(S)
          ALOECM(I,III)=LOG10(ECM)
          ASIG(I,III)=0.
  111   CONTINUE
*
* *** calling calculation of x- section
*
        IF(IPIM.EQ.2 )THEN
          CALL SIGMA2(ECM)
          IF(I.EQ.1 .AND. IOUTPO.GE.0 ) WRITE(6,*)
     &     ' s-dep. by integr.with Y,PHI,LMD'
        ELSE
          CALL SIGMA2(ECM)
          IF(I.EQ.1 .AND. IOUTPO.GE.0 ) WRITE(6,*)
     &      ' s-dep. by integr.with Y,PHI,LMD  (DEFAULT)'
        ENDIF
*
*
C*       the preceeding coresponds to a call to SIGMA2 except extra prin
* ***   getting plot array ASIG(I,J) :
        ASIG(I,1)=SIGTOT
        ASIG(I,2)=SIGINE
        ASIG(I,3)=SIGHIN
        ASIG(I,4)=SIGSOF
        ASIG(I,5)=SIGHAR
        ASIG(I,6)=SIGTRP
        ASIG(I,7)=SIGTOT-SIGINE
        ASIG(I,8)=SIGINE-SIGHIN
        ASIG(I,9)=SIGD
        WRITE (6,1007)ECM,SIGTOT,SIGINE,SIGEL,SIGD
 1007   FORMAT (' ECM,SIGTOT,SIGINE,SIGEL,SIGD',F10.1,4E14.3)
  100 CONTINUE
*
*     energy loop ends
*---------------------------------------------------------------------
*                  print out results
*    ------------------------------------------------------------------
      WRITE (6,991)
  991 FORMAT (//' shown as line printer plott'/' with'/
*    J: drawn quantities:
     1 '  (*) SIGTOT total x-section',
     2 '  (2) SIGINE inelastic x-section'/
     3 '  (3) SIGHIN hard inelastic cross section, one or more jets',
     4 '  (4) SIGSOF   input soft x-section'/
     5 '  (5) SIGHAR   input hard x-sections',
     6 '  (6) SIGTRP   input diffractive x-section (triple pomeron)'/
     7 '  (7) SIGTOT-SIGINE  elastic x-section',
     8 '  (8) SIGINE-SIGHIN  non-hard inelastic x-section, (no jets)'/
     9 '  (9) SIGD   diffractive xross section '/
     * '  are plotted against LOG(10)of(CMENERGY)' //)
*
      CALL PLOT(ALOECM,ASIG,450,9,50,ZERO, 0.1*ONE,ZERO, 2.0*ONE)
*
*  special output on unit number 7
C  I kept it as it was
      IF (ITOPD.EQ.1) THEN
        WRITE(7,95)
   95   FORMAT(' NEW FRAME'/' SET FONT DUPLEX'/' SET SCALE X LOG'/
     *  ' SET LIMITS X FROM 1.0 TO 1E5 Y FROM 0. TO 200'/
     *  ' TITLE TOP < TOTAL,INEL. AND HARD (MINIJET) CROSS SECT.<'/
     *  ' TITLE BOTTOM <C.M.ENERGY [GEV]<'/
     *  ' TITLE < DUAL UNITARIZATION OF SOFT AND HARD CROSS SECTIONS<'/
     *  ' TITLE LEFT LINES=-1 <CROSS SECTION [MB]<'/
     *  ' TITLE 3 8.5 < SOLID = TOTAL X.S. <'/
     *  ' TITLE  < DASHED= INELASTIC X.S. <'/
     *  ' TITLE  < DOTTED= HARD X.S.<'/
     *  ' TITLE  < DOT-DASH= HARD INPUT X.S. <'/
     *  ' TITLE  < DOT-DASH= ELASTIC X.S. <')
   92   FORMAT (5F15.5)
        DO 94 IUU=1,7
          IF (IUU.EQ.4)GO TO 94
          IF (IUU.EQ.6)GO TO 94
          IF (IUU.EQ.1) WRITE(7,97)
   97           FORMAT (' SET TEXTURE SOLID')
          IF (IUU.EQ.2) WRITE(7,98)
   98           FORMAT (' SET TEXTURE DASHES')
          IF (IUU.EQ.3) WRITE(7,99)
   99           FORMAT (' SET TEXTURE DOTS')
          IF (IUU.EQ.5) WRITE(7,197)
  197           FORMAT (' SET TEXTURE DOTDASH')
          DO 93 IU=2,46
            WRITE(7,92)AECM(IU,IUU),ASIG(IU,IUU)
   93     CONTINUE
          WRITE(7,96)
   96           FORMAT (' JOIN')
   94   CONTINUE
      ENDIF
*     ending IF(ITOP=1) special output
      RETURN
      END
*
******************************************************************
*  end dtupom90
******************************************************************

