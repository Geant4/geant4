      BLOCK DATA BLKD41
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
*KEEP,PANAME.
C------------------
C
C     /PANAME/ CONTAINS PARTICLE NAMES
C        BTYPE  = LITERAL NAME OF THE PARTICLE
C
      CHARACTER*8 BTYPE
      COMMON /PANAME/ BTYPE(30)
*KEEP,NUCIMP.
      COMMON /NUCIMP/ PRMOM(5,248),TAMOM(5,248), PRMFEP,PRMFEN,TAMFEP,
     +TAMFEN, PREFEP,PREFEN,TAEFEP,TAEFEN, PREPOT(210),TAEPOT(210),
     +PREBIN,TAEBIN,FERMOD,ETACOU
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEEP,DROPPT.
      LOGICAL INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA, IPADIS,
     +ISHMAL,LPAULI
      COMMON /DROPPT/ INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA,
     +IPADIS,ISHMAL,LPAULI
*KEEP,HADTHR.
      COMMON /HADTHR/ EHADTH,INTHAD
*KEEP,NSHMAK.
      COMMON /NSHMAK/ NNSHMA,NPSHMA,NTSHMA,NSHMAC,NSHMA2
*KEEP,REJEC.
      COMMON /REJEC/ IRCO1,IRCO2,IRCO3,IRCO4,IRCO5, IRSS11,IRSS12,
     +IRSS13,IRSS14, IRSV11,IRSV12,IRSV13,IRSV14, IRVS11,IRVS12,IRVS13,
     +IRVS14, IRVV11,IRVV12,IRVV13,IRVV14
*KEEP,DSHM.
      COMMON /DSHM/ RASH,RBSH,BMAX,BSTEP,SIGSH,ROSH,GSH,
     *              BSITE(0:1,200),NSTATB,NSITEB
*KEEP,DAMP.
      COMPLEX*16 CA,CI
      COMMON /DAMP/   CA,CI,GA
*KEEP,INTMX.
      PARAMETER (INTMX=2488,INTMD=252)
*KEND.
      COMMON /DIFFRA/ ISINGD,IDIFTP,IOUDIF,IFLAGD
      COMMON /NOMIJE/ PTMIJE(10),NNMIJE(10)
      DATA PTMIJE /5.D0,7.D0,9.D0,11.D0,13.D0,15.D0,17.D0
     +,19.D0,21.D0,23.D0 /
*
*---rejection counters
      DATA IRCO1,IRCO2,IRCO3,IRCO4,IRCO5 /5*0/
      DATA IRSS11,IRSS12,IRSS13,IRSS14,IRSV11,IRSV12,IRSV13,IRSV14 /8*0/
      DATA IRVS11,IRVS12,IRVS13,IRVS14,IRVV11,IRVV12,IRVV13,IRVV14 /8*0/
C-------------------
      DATA INTHAD  /0/
*---predefinition of nuclear potentials, average binding energies, ...
      DATA PREPOT /210*0.0/
      DATA TAEPOT /210*0.0/
      DATA TAEBIN,PREBIN,FERMOD /2*0.0D0,0.6D0/
*---internal particle names
      DATA BTYPE /'PROTON  ' , 'APROTON ' , 'ELECTRON' , 'POSITRON' ,
     +'NEUTRIE ' , 'ANEUTRIE' , 'PHOTON  ' , 'NEUTRON ' , 'ANEUTRON' ,
     +'MUON+   ' , 'MUON-   ' , 'KAONLONG' , 'PION+   ' , 'PION-   ' ,
     +'KAON+   ' , 'KAON-   ' , 'LAMBDA  ' , 'ALAMBDA ' , 'KAONSHRT' ,
     +'SIGMA-  ' , 'SIGMA+  ' , 'SIGMAZER' , 'PIZERO  ' , 'KAONZERO' ,
     +'AKAONZER' , 'RESERVED' , 'BLANK   ' , 'BLANK   ' , 'BLANK   ' ,
     +'BLANK   ' /
*---print options
      DATA IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR /0, 0, 0, -1, 0,
     +0, 0,  0/
C-------------------
      DATA INTPT, FERMP, IHADSS,IHADSV,IHADVS,IHADVV, IHADA /.TRUE.,
     +.TRUE., 4*.FALSE., .TRUE./
      DATA IPADIS, ISHMAL, LPAULI /.FALSE., .FALSE., .TRUE./
C----------------------------------
      DATA NSHMAC /0/
      DATA NSHMA2 /0/
*---parameters for Glauber initialization / calculation
      DATA NSTATB, NSITEB /2000, 200/
      DATA CI /(1.0,0.0)/
*---parameters for combination of q-aq chains to color ropes
C     DATA LCOMBI /.FALSE./, NCUTOF /100/
      DATA ISINGD,IDIFTP,IOUDIF,IFLAGD /0,0,0,0/
*
      END
************************************************************************
************************************************************************
*
       SUBROUTINE DTTEST(CODEWD,WHAT,SDUM)
*
*      -- not for normal user --
*      contains input options unrecognized by DTPREP and
*      performs special initialisations or tasks for program devoloping
*
C   COMMON fully commented in DTPREP
*
*   **********************************************************************
*   *  DESCRIPTION OF THE COMMON BLOCK(S), VARIABLE(S) AND DECLARATIONS  *
*   **********************************************************************
*
*
*  /USER/ contains the parameters, expected to be modified by normal user
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      CHARACTER*80 TITLE
      CHARACTER*8 PROJTY,TARGTY
C     COMMON /USER/TITLE,PROJTY,TARGTY,CMENER,ISTRUF
C    &            ,ISINGD,IDUBLD,SDFRAC,PTLAR
      COMMON /USER1/TITLE,PROJTY,TARGTY
      COMMON /USER2/CMENER,SDFRAC,PTLAR,ISTRUF,ISINGD,IDUBLD
*
*     /COLLE/           contains the input specifying the MC. run
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      COMMON /COLLE/NEVHAD,NVERS,IHADRZ,NFILE
*
*     /COLLIS/       contains the input specifying the considered event
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      COMMON/COLLIS/S,IJPROJ,IJTAR,PTTHR,PTTHR2,IOPHRD,IJPRLU,IJTALU
*
*
*     /BOOKLT/ contains the final  particle names and PPDB-numbers
*              of 30 final particle types
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CHARACTER*8 BTYPE
      COMMON/BOOKLT/BTYPE(30),NBOOK(30)
*
*     /POLMN/           stores arrays describing probabilities of parton
*                          configurations
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       COMMON /POLMN0/PDIFR,PHARD,PSOFT,ALFAH,BETAH,
     *              SIGTOT,SIGQEL,SIGEL,SIGINE,SIGHIN,SIGD,SIGDD
*
*     /POMTYP/ contains parameters determining X-sections
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      COMMON /POMTYP/IPIM,ICON,ISIG,LMAX,MMAX,NMAX,DIFEL,DIFNU
*
*     various smaller commons
*     in alphabetical order
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      COMMON /DROPJJ/DROPJT,DROPVA
      COMMON /GLUSPL/NUGLUU,NSGLUU
      COMMON /PTLARG/XSMAX
      COMMON /PTSAMP/ ISAMPT
      COMMON /STARS/ISTAR2,ISTAR3
      COMMON /STRUFU/ISTRUM,ISTRUT
      COMMON /POPCOR/PDB,AJSDEF 
*
*  ********************************************************************
*     declarations outside of commons
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CHARACTER*8 CODEWD,SDUM
      DIMENSION WHAT(6)
*  ********************************************************************
*
*
*   **********************************************************************
*   *         Print the warning that no normal codeword                  *
*   **********************************************************************
*
      WRITE(6,9)
 9    FORMAT( ' special code word was used ')
*
*
*  The following additional CODEWD options exist at the moment:
*                   RANDOMIZ    SIGMAPOM    PARTEV      SELHARD
*       GLUSPLIT        
*                     
*   The cards marked with )+ have to be followed by data cards of
*   special format
*
*
*
*   **********************************************************************
*   *                 parse stored imput card                            *
*   **********************************************************************
*
*
*
*  *********************************************************************
*                       input card: CODEWD = RANDOMIZE
*                       Sets the SEED for the random number generators
*
*     WHAT(1)       = 1:   gets testrun      otherwise: reset
*     WHAT(2,3,4,5) = giving the SEED for the random number generators.
*           Default   ISEED1/2/3/4=12,34,56,78
*     Note.  It is advisable to use only the seeds given by the
*     program in earlier runs. Otherwise the number sequence might
*     have defects in its randomness.
C      Since 1999 RANDOMIZE initializes the RM48 generator
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*
C     IF(CODEWD.EQ.'RANDOMIZ') THEN
*
*     Reinitialize random generator
C     IF(WHAT(2).NE.0.D0) THEN
C       ISEED1=WHAT(2)
C       ISEED2=WHAT(3)
C       ISEED3=WHAT(4)
C       ISEED4=WHAT(5)
C       CALL RNDMST(ISEED1,ISEED2,ISEED3,ISEED4)
C     ENDIF
*     Test random generator and test
C     IF(WHAT(1).EQ.1) CALL RNDMTE(1)
*
C               GO TO 1
*
*
*  *********************************************************************
*                       input card: CODEWD = SIGMAPOM
*                       Defines the options for the calculation of the u
*                       tarized hard and soft multi-pomeron cross sectio
*                       and demands a testrun
*
*     WHAT(1) = ITEST testrun for ITEST = 1
*     WHAT(2) = ISIG  characterizing X sections, see SIGSHD
*                                                      default: 10
*         Only ISIG=10 kept in dpmjet-II.5
*
*     WHAT(3) = characterizes the method to calculate
*               the cut pomeron X-section SIGMA(Lsoft,Mhard,Ntrp) and
*               how to attribute them to strings
C!!!!!!!!!!!!!!!!!!!
C!!!!!!!!!!!!!!!!!!!    THE FOLLOWING DISTRIBUTIONS WITH 2 CHANNEL
C!!!!!!!!!!!!!!!!!!!    EIKONAL + H. MASS DIFFRACTION + HARD SCATTERING
C!!!!!!!!!!!!!!!!!!!
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
C
C               =482: SEE PRBLM2
C
*                                                   DEFAULT: 482
*  for all cases
*     WHAT(4) = LMAX  maximal considered number soft Pomerons
*     WHAT(5) = MMAX  maximal considered number hard Pomerons
*     WHAT(6) = NMAX  maximal considered number trippel Pomerons
*
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*
C     ELSEIF(CODEWD.EQ.'SIGMAPOM') THEN
      IF(CODEWD.EQ.'SIGMAPOM') THEN
*
                   ITEST =    WHAT(1)
                   IF(WHAT(2).NE.0.) ISIG =INT(WHAT(2))
                   IF(WHAT(3).NE.0.)  IPIM =INT(WHAT(3))
                   IF(IPIM.GT.10) THEN
                          ICON = IPIM /10
                          IPIM = IPIM-10*ICON
                   ENDIF
                   IF(IPIM.EQ.1) THEN
                                     DIFEL = WHAT(4)
                                     DIFNU = WHAT(5)
                                     LMAX =INT(WHAT(6))
                                     MMAX = LMAX
                                 ELSE
                                     LMAX =INT(WHAT(4))
                                     MMAX =INT(WHAT(5))
                                     NMAX =INT(WHAT(6))
                                 ENDIF
                   IF (ITEST.EQ.1)CALL POMDI
                   GO TO 1
*
*  ********************************************************************
*                       input card: CODEWD = GLUSPLIT
*                       Prevents the splitting of Gluons
*
*     WHAT(1)=NUGLUU                                         Default: 1.
*                    =1.  only one jet in hard gluon scattering
*     WHAT(2)=NSGLUU                                         Default: 0.
*                    =0.  two jets in soft sea gluons
*                    =1.  only one jet in soft sea gluons
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*
      ELSEIF(CODEWD.EQ.'GLUSPLIT') THEN
*
                NUGLUU  =  WHAT(1)
                NSGLUU  =  WHAT(2)
                GO TO 1
*
*  ********************************************************************
*
*  *********************************************************************
*                       input card: CODEWD = PARTEV
*                       defines the parton level collision events
*                               the X's PT's and flavors and
*                       demands a test run
*
*     WHAT(1) = 1.: testrun ;  other values: no testrun
*     WHAT(2) = number of events is NPEV                   default: 30
*     WHAT(3) = version of PARTEV is NVERS                 default: 1
*                NVERS=1  all hard partons considered to be gluons
*                         SOFT X-VALUEA BY REJECTION
*                NVERS=2  all hard partons considered to be gluons
*                         SOFT X-VALUEA BY AURENCHE -MAIRE METHOD
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*
      ELSEIF(CODEWD.EQ.'PARTEV  ') THEN
*
                ITEST = INT(WHAT(1))
                IF (WHAT(2).EQ.0.D0)NPEV=30
                IF (WHAT(2).NE.0.D0)NPEV=INT(WHAT(2))
                IF (WHAT(3).EQ.0.D0)NVERS=1
                IF (WHAT(3).NE.0.D0)NVERS=INT(WHAT(3))
*       first initialize NSOFT-NHARD event selection
*                      corresponing to choosen one options
                IF(ITEST.EQ.1) THEN
                  IF(IPIM.EQ.2)CALL PRBLM2(CMENER)
*       initialize hard scattering
                  CALL JTDTU(0)
                  CALL SAMPPT(0,PT)
C                 CALL PARTEV(NPEV)
                  CALL TIMDAT
                ENDIF
           GO TO 1
*
*  *********************************************************************
*                       input card: CODEWD = SELHARD
*                       defines the selection of X's,PT's and
*                                        flavors for hard scattering
*
*     WHAT(2) = IOPHRD selects the model                  default: 2
*   WHAT(4) = DROPJT=10:DROP FIRST HARD JET PAIR           DEFAULT: 0
*             THIS OPTION SHOULD ALLOW TO SIMULATE DRELL-YAN
*              OR W OR Z PRODUCTION EVENTS
*   WHAT(5) = PTTHR   threshold PT for hard scattering  default: 3.
*   WHAT(6) = PTTHR2 threshold PT for FIRST HARD scattering  default: 3.
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*
      ELSEIF(CODEWD.EQ.'SELHARD ') THEN
*
        IF(WHAT(2).NE.0.D0)IOPHRD=INT(WHAT(2))
        IF(WHAT(4).NE.0.D0)DROPJT=WHAT(4)
        IF(WHAT(6).NE.0.D0)PTTHR2=WHAT(6)
        IF(WHAT(5).NE.0.D0)THEN
          PTTHR=WHAT(5)
          IF(CMENER.LT.2000.0D0.AND.ISIG.EQ.3)PTTHR=WHAT(5)
          IF (CMENER.GE.2000.0D0.AND.ISIG.EQ.3)
     *                   PTTHR=0.25*LOG(CMENER/2000.)+2.
          IF(PTTHR2.LT.PTTHR)PTTHR2=PTTHR
          IF(ISTRUT.EQ.1)THEN
            PTTHR=2.1+0.15*(LOG10(CMENER/50.))**3
            PTTHR2=PTTHR
          ELSEIF(ISTRUT.EQ.2)THEN
            PTTHR=2.5+0.12*(LOG10(CMENER/50.))**3
            PTTHR2=PTTHR
          ENDIF
          WRITE(6,1244)PTTHR
 1244     FORMAT (' THRESHOLD PT FOR HARD SCATTERING PTTHR=',F12.2)
        ENDIF
      GO TO 1
*
*  *********************************************************************
*                       input card: CODEWD = XSLAPT
*                       calculates inclusive large PT cross sections
*                       and testrun of the large pt and minijet sampling
*                       in file laptabr
*
*     WHAT(I) =  not used at this level
*         special DATA CARDS are required
*         see description at beginning of LAPTABR
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ----
*
      ELSEIF(CODEWD.EQ.'XSLAPT  ') THEN
*
        CALL TIMDAT
        CALL LAPTAB
        CALL TIMDAT
      GO TO 1
*
*  *********************************************************************
*                       parameter card: CODEWD = SAMPT
*                       defines the options of soft pt sampling in
*                       subroutine SAMPPT
*     WHAT(1) = number defining the option of soft pt sampling
*               (default: 4 )
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ----
*
      ELSEIF(CODEWD.EQ.'SAMPT   ') THEN
*
        ISAMPT = INT( WHAT(1) )
        IF( ISAMPT.LT.0 .OR. ISAMPT.GT.4 ) ISAMPT=0
        GO TO 1
*
*  *********************************************************************
*         ending special parsing the code word of "input card"
*
      ENDIF
*
*               1)  not recognized cards
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ----
*     a warning will be issued in DTPREP
      RETURN
*
*               2)   recognized cards
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ----
*     action was done,
*     CODEWD is set to a value ignored in DTPREP
*
 1    CODEWD='-zzzzzzz'
      RETURN
      END
*
************************************************************************
************************************************************************
*
                      BLOCK DATA BOOKLE
*                                               
*  *********************************************************************
*                         /BOOKLT/
*
* neeeded in the following routines: DTUMAIN
*
* description
*     /BOOKLT/ contains the final  particle names and PPDB-numbers
*        BTYPE  = literal name of the particle
*        NBOOK  = the number of the particle
*                 proposed in the particle data booklet (90)
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      CHARACTER*8 BTYPE
      COMMON/BOOKLT/BTYPE(30),NBOOK(30)
*
      DATA BTYPE /'PROTON  ' , 'APROTON ' , 'ELECTRON' ,
     1            'POSITRON' , 'NEUTRIE ' , 'ANEUTRIE' ,
     2            'PHOTON  ' , 'NEUTRON ' , 'ANEUTRON' ,
     3            'MUON+   ' , 'MUON-   ' , 'KAONLONG' ,
     4            'PION+   ' , 'PION-   ' , 'KAON+   ' ,
     5            'KAON-   ' , 'LAMBDA  ' , 'ALAMBDA ' ,
     6            'KAONSHRT' , 'SIGMA-  ' , 'SIGMA+  ' ,
     7            'SIGMAZER' , 'PIZERO  ' , 'KAONZERO' ,
     9            'AKAONZER' , '        ' , '        ' ,
     Z            '        ' , '        ' , '        ' /
*
*
      DATA NBOOK / 2212      , -2212      ,  11        ,
     1            -11        ,  14        , -14        ,
     2             22        ,  2112      , -2112      ,
     3            -13        ,  13        ,  130       ,
     4             211       , -211       ,  321       ,
     5            -321       ,  3122      , -3122      ,
     6             310       ,  3114      ,  3224      ,
     7             3214      ,  111       ,  311       ,
     9            -311       ,  0         ,  0         ,
     Z              0        ,  0         ,  0         /
*
      END


C______________________________________________________________________
      SUBROUTINE SAMPPT(MODE,PT)
* pt for partons at the end of soft chains
*  this pt is sampled from the distribution
*  exp(-b*pt^2) with pt=0..ptcut
*  MODE = 0  - initialization to determine parameter b from total soft
*              and differential hard cross section
*  MODE = 1  - sample pt
*  MODE = 2  - PLOT PT DISTRIBUTION SAMPLED
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( ZERO=0.D0, ONE=1.D0)
      PARAMETER ( ALFA=0.56268D-01, BETA=0.17173D+03 )
      PARAMETER ( ACC = 0.0001D0 )
      COMMON /XSECPT/ PTCUT,SIGS,DSIGH
      COMMON /SIGMA / SIGSOF,BS,ZSOF,SIGHAR,FILL(7)
      COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
C repl. COMMON/COLLIS/ECM,S,IJPROJ,IJTAR,PTTHR,IOPHRD,IJ1LU,IJ2LU,PTTHR2
      COMMON/COLLIS/S,IJPROJ,IJTAR,PTTHR,PTTHR2,IOPHRD,IJPRLU,IJTALU
      CHARACTER*80 TITLE
      CHARACTER*8 PROJTY,TARGTY
C     COMMON /USER/TITLE,PROJTY,TARGTY,CMENER,ISTRUF
C    &            ,ISINGD,IDUBLD,SDFRAC,PTLAR
      COMMON /USER1/TITLE,PROJTY,TARGTY
      COMMON /USER2/CMENER,SDFRAC,PTLAR,ISTRUF,ISINGD,IDUBLD
*
      COMMON/PTSAMP/ ISAMPT
      DIMENSION PPTT(50),DPPTT(50)
      DATA ECM0 /0.1D0/
C     to keep identical commons for patchy ect
      ECM=CMENER
      PTTHR=2.5+0.12*(LOG10(CMENER/50.))**3
      PTCUT=PTTHR      
      CALL SIGSHD(ECM)
      IF ( MODE.EQ.0 ) THEN
        DO 201 II=1,50
          PPTT(II)=II*PTCUT/50.
          DPPTT(II)=0.
  201   CONTINUE
        SIGS  = 0.15*SIGSOF
        IF(ECM.LT.1000.)THEN
          AACUCU=0.85*(ECM-400.)/600.
          SIGS=(1.-AACUCU)*SIGSOF
        ENDIF
C*************************************************************
C
C       OPTIONS FOR SOFT PT SAMPLING
C
CWRITE(6,'(A,4E12.4)')' SAMPPT:ECM,PTCUT,SIGS,DSIGH',
C    *	ECM,PTCUT,SIGS,DSIGH
        IF(ECM0.NE.ECM)THEN
C         WRITE(6,'(A,5E12.4)')' SAMPPT:ECM,PTCUT,SIGS,DSIGH,SIGHAR',
C    *    ECM,PTCUT,SIGS,DSIGH,SIGHAR
C         WRITE(6,5559)PTCUT,SIGSOF,SIGHAR,ISAMPT
 5559     FORMAT(' SAMPPT:PTCUT,SIDSOF.SIGHATD,ISAMPT:',3E12.3,I5)
        ENDIF
        IF( ISAMPT.EQ.0 ) THEN
          C  = DSIGH/(2.*SIGS*PTCUT)
          B  = BSOFPT(ACC,C,PTCUT)
        IF( IOUTPA.GE.1 )WRITE(6,9010)MODE,ISAMPT,PTCUT,SIGS,DSIGH,B
     * ,C,SIGSOF,SIGHAR,RMIN
        ELSEIF( ISAMPT.EQ.1 ) THEN
          EB = ECM/BETA
          C = ALFA*LOG(EB)
          B  = BSOFPT(ACC,C,PTCUT)
        IF( IOUTPA.GE.1 )WRITE(6,9010)MODE,ISAMPT,PTCUT,SIGS,DSIGH,B
     * ,C,SIGSOF,SIGHAR,RMIN
        ELSEIF( ISAMPT.EQ.2 ) THEN
          B=-6.
        IF( IOUTPA.GE.1 )WRITE(6,9010)MODE,ISAMPT,PTCUT,SIGS,DSIGH,B
     * ,C,SIGSOF,SIGHAR,RMIN
        ELSEIF( ISAMPT.EQ.3 ) THEN
          B=1.E-06
        IF( IOUTPA.GE.1 )WRITE(6,9010)MODE,ISAMPT,PTCUT,SIGS,DSIGH,B
     * ,C,SIGSOF,SIGHAR,RMIN
        ELSEIF( ISAMPT.EQ.4)THEN
          AAAA=PTCUT**2*(SIGSOF+SIGHAR)
          IF (AAAA.LE.0.00001D0) THEN
            AAAA=ABS(AAAA)+0.0002
C           WRITE(6,5559)PTCUT,SIGSOF,SIGHAR
C5559       FORMAT(' SAMPPT:PTCUT,SIDSOF.SIGHATD:',3E12.3)
          ENDIF
          C=SIGHAR/AAAA
          B  = 0.5*BSOFPT(ACC,C,PTCUT)
	IF( IOUTPA.GE.1 )WRITE(6,9010)MODE,ISAMPT,PTCUT,SIGS,DSIGH,B
     * ,C,SIGSOF,SIGHAR,RMIN
        ENDIF
        ECM0=ECM
C*************************************************************
        RMIN = EXP(B*PTCUT**2)
C
C        IOUTPO=IOUTPA
C        IOUTPA=1
        IF( IOUTPA.GE.1 )WRITE(6,9010)MODE,ISAMPT,PTCUT,SIGS,DSIGH,B
     * ,C,SIGSOF,SIGHAR,RMIN
9010  FORMAT(' SAMPPT MODE,ISAMPT,PTCUT,SIGS,DSIGH,B,C,SIGSOF',
     *' SIGHAR,RMIN ',
     *         2I2,F5.2,7E13.6)
C        IOUTPA=IOUTPO
      ELSEIF ( MODE.EQ.1 ) THEN
        IF( IOUTPA.GE.1 )WRITE(6,9010)MODE,ISAMPT,PTCUT,SIGS,DSIGH,B
     * ,C,SIGSOF,SIGHAR,RMIN
        PTT   =LOG(1.0-RNDM(V)*(1.0-RMIN))/(B+0.00001D0)
        PT=SQRT(PTT)
        IIPT=PT*50./PTCUT+1.
        IIPT=MIN( IIPT,50 )
        DPPTT(IIPT)=DPPTT(IIPT)+1./(PT+0.000001D0)
C       WRITE(6,111)MODE,PTT,PT,B,RMIN
C 111   FORMAT ('SAMPPT: MODE,PTT,PT,B RMIN',I5,4E15.8)
      ELSEIF(MODE.EQ.2)THEN
        DO 202 II=1,50
          DPPTT(II)=LOG10(1.E-8+DPPTT(II))
  202   CONTINUE
        IF(IOUXEV.GE.-1)THEN
         WRITE (6,203)
  203   FORMAT(' PT DISTRIBUTION OF SOFT PARTONS AS SAMPLED IN BSOFPT')
         CALL PLOT(PPTT,DPPTT,50,1,50,ZERO,PTCUT/50.D0,ZERO,0.05D0*ONE)
        ENDIF
      ENDIF
      RETURN
      END
C****************************************************************8****
       REAL
     *      * 8
     * FUNCTION BSOFPT(ACC,CC,PPTCUT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      LOGICAL SUCCES
      COMMON /BSOFF1/C,PTCUT
      COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
      EXTERNAL BSOFC1,BSOF1
      DIMENSION X(50),Y(50)
      C=CC
      PTCUT=PPTCUT
      DO 100 I=1,50
        X(I)=-2.5+I*0.1
        Y(I)=BSOF1(X(I))
  100 CONTINUE
C     CALL PLOT (X,Y,50,1,50,-2.5D0,0.1D0,-1.D0,0.02D0)
      IF(C.LT.1.D-10) THEN
        BB=-30.
        GO TO 999
      ENDIF
C     IF (C.GT.1.) THEN
        KKKK=0
        JJJJ=0
        B1=C+3.
        B2=0.0001
        GO TO 300
  400 CONTINUE
      KKKK=KKKK+1
C     ENDIF
C     IF (C.LT.1.)THEN
        B1=-0.00001
        B2=-3.
C     ENDIF
  300 CONTINUE
      CALL ZBRAC(BSOF1,B1,B2,SUCCES)
      IF (.NOT.SUCCES)THEN
        IF (KKKK.EQ.0)GO TO 400
        JJJJ=1
      ENDIF
      IF(IOUXEV.GE.1)WRITE(6,111)B1,B2
  111 FORMAT(2F10.4)
      IF (SUCCES)THEN
        BB=RTSAFE(BSOFC1,B1,B2,ACC)
      ENDIF
      IF (JJJJ.EQ.1)BB=0.
  999 CONTINUE
      BSOFPT=BB
      RETURN
      END
      SUBROUTINE BSOFC1(B,F,DF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      COMMON /BSOFF1/C,PTCUT
      AAA=EXP(B*PTCUT**2)
      F=C*(AAA-1.)-B*AAA
      DF=C*PTCUT**2*AAA-AAA
     *                       -B*PTCUT**2*AAA
      RETURN
      END
*
*******************************************************************
*
       REAL
     *      * 8
     *  FUNCTION BSOF1(B)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      COMMON /BSOFF1/C,PTCUT
      QQQ=B*PTCUT**2
      AAA=0.
      IF(QQQ.GT.-60.) THEN
        AAA=EXP(B*PTCUT**2)
      ENDIF
      BSOF1=C*(AAA-1.)-B*AAA
C     WRITE(6,10)B,PTCUT,BSOF1,C
C  10 FORMAT (4E15.4)
      RETURN
      END
*
*******************************************************************************
       REAL
     *      * 8
     *  FUNCTION RTSAFE(FUNCD,X1,X2,XACC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER (MAXIT=200,ITEPRI=0)
      COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
      CALL FUNCD(X1,FL,DF)
C     IF(IOUXEV.GE.0.AND.ITEPRI.EQ.1) WRITE(6,9999) FL,DF
      CALL FUNCD(X2,FH,DF)
C     IF(IOUXEV.GE.0.AND.ITEPRI.EQ.1) WRITE(6,9999) FH,DF
      IF(FL*FH.GE.0.) PAUSE 'ROOT MUST BE BRACKETED'
      IF(FL.LT.0.)THEN
        XL=X1
        XH=X2
      ELSE
        XH=X1
        XL=X2
        SWAP=FL
        FL=FH
        FH=SWAP
      ENDIF
      RTSAFE=.5*(X1+X2)
      DXOLD=ABS(X2-X1)
      DX=DXOLD
      CALL FUNCD(RTSAFE,F,DF)
C     IF(IOUXEV.GE.0.AND.ITEPRI.EQ.1) WRITE(6,9998) RTSAFE,F,DF
C     IF(IOUXEV.GE.0.AND.ITEPRI.EQ.1) WRITE(6,9996)
      DO 11 J=1,MAXIT
C       IF(IOUXEV.GE.0.AND.ITEPRI.EQ.1)
C    *  WRITE(6,9997) RTSAFE,XH,XL,DXOLD,F,DF
        VR1 = VAR( RTSAFE,XH,DF,F )
        VR2 = VAR( RTSAFE,XL,DF,F )
C       IF(IOUXEV.GE.0.AND.ITEPRI.EQ.1) WRITE(6,9995) VR1,VR2
C       IF(((RTSAFE-XH)*DF-F)*((RTSAFE-XL)*DF-F).GE.0.
C    *      .OR. ABS(2.*F).GT.ABS(DXOLD*DF) ) THEN
        IF( VR1*VR2 .GE. 0.
     *      .OR. ABS(2.*F).GT.ABS(DXOLD*DF) ) THEN
          DXOLD=DX
          DX=0.5*(XH-XL)
          RTSAFE=XL+DX
          IF(XL.EQ.RTSAFE)RETURN
        ELSE
          DXOLD=DX
C         IF(IOUXEV.GE.0.AND.ITEPRI.EQ.1) WRITE(6,9999) F,DF
          DX=F/DF
          TEMP=RTSAFE
          RTSAFE=RTSAFE-DX
          IF(TEMP.EQ.RTSAFE)RETURN
        ENDIF
        IF(ABS(DX).LT.XACC) RETURN
        CALL FUNCD(RTSAFE,F,DF)
        IF(F.LT.0.) THEN
          XL=RTSAFE
          FL=F
        ELSE
          XH=RTSAFE
          FH=F
        ENDIF
11    CONTINUE
      PAUSE 'RTSAFE EXCEEDING MAXIMUM ITERATIONS'
      RETURN
9995  FORMAT('  VR1,VR2:',2E12.5)
9996  FORMAT('  RTSAFE,XH,XL,DXOLD,F,DF IN LOOP 11 J=1,MAXIT')
9997  FORMAT(3X,6E10.3)
9998  FORMAT('  RTSAFE: RTSAFE,F,DF =',3E12.5)
9999  FORMAT('  RTSAFE: F,DF =',2E12.5)
      END
*
*****************************************************************
       REAL
     *      * 8
     *  FUNCTION VAR(A,B,C,D)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER( AMBMAX = 1.0D+38, EPSI = 1.2D-38, ONE=1.D0 )
      AMB = A - B
      SIAB= SIGN(ONE,AMB)
      ABL = ABS(AMB)
      ABL = LOG10( ABL + EPSI )
      SICC= SIGN(ONE, C )
      CCL = ABS( C )
      CCL = LOG10( CCL + EPSI )
      RCHECK=ABL + CCL
      IF( RCHECK .LE. 38.D0 ) THEN
        VAR = AMB*C-D
      ELSE
        VAR = AMBMAX*SIAB*SICC - D
      ENDIF
      IF( VAR .GT. 1.0D+18 ) VAR = 1.0E+18
      IF( VAR .LT. -1.0D+18 ) VAR = -1.0E+18
      RETURN
      END
C
      SUBROUTINE ZBRAC(FUNC,X1,X2,SUCCES)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      EXTERNAL FUNC
      PARAMETER (FACTOR=1.6D0,NTRY=50)
      LOGICAL SUCCES
      IF(X1.EQ.X2)PAUSE 'You have to guess an initial range'
      F1=FUNC(X1)
      F2=FUNC(X2)
      SUCCES=.TRUE.
      DO 11 J=1,NTRY
        IF(F1*F2.LT.0.D0)RETURN
        IF(ABS(F1).LT.ABS(F2))THEN
          X1=X1+FACTOR*(X1-X2)
          F1=FUNC(X1)
        ELSE
          X2=X2+FACTOR*(X2-X1)
          F2=FUNC(X2)
        ENDIF
11    CONTINUE
      SUCCES=.FALSE.
      RETURN
      END
*
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
c     subroutine hibldr
c     COMMON /DREAC/ umo(296),plabf(296),siin(296),wk(5184),
c    *              nrk(2,268),nure(30,2)
c     read(2,1)umo,plabf,siin,wk
c   1 format (5e16.7)
c     read(2,2)nrk,nure
c   2 format (8i10)
c     return
c     end
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE CHECKF(EPN,PPN,IREJ,IORIG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
*KEEP,HKKEVT.
c     INCLUDE (HKKEVT)
      PARAMETER (AMUAMU=0.93149432D0)
      PARAMETER (NMXHKK= 89998)
c     PARAMETER (NMXHKK=25000)
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK), JMOHKK
     +(2,NMXHKK),JDAHKK(2,NMXHKK), PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK
     +(4,NMXHKK)
C
C                       WHKK(4,NMXHKK) GIVES POSITIONS AND TIMES IN
C                       PROJECTILE FRAME, THE CHAINS ARE CREATED ON
C                       THE POSITIONS OF THE PROJECTILE NUCLEONS
C                       IN THE PROJECTILE FRAME (TARGET NUCLEONS IN
C                       TARGET FRAME) BOTH POSITIONS ARE THREFORE NOT
C                       COMPLETELY CONSISTENT. THE TIMES IN THE
C                       PROJECTILE FRAME HOWEVER ARE OBTAINED BY
C                       LORENTZ TRANSFORMING FROM THE LAB SYSTEM.
C
C  Based on the proposed standard COMMON block (Sjostrand Memo 17.3,89)
C
C NMXHKK: maximum numbers of entries (partons/particles) that can be
C    stored in the commonblock.
C
C NHKK: the actual number of entries stored in current event. These are
C    found in the first NHKK positions of the respective arrays below.
C    Index IHKK, 1 <= IHKK <= NHKK, is below used to denote a given
C    entry.
C
C ISTHKK(IHKK): status code for entry IHKK, with following meanings:
C    = 0 : null entry.
C    = 1 : an existing entry, which has not decayed or fragmented.
C        This is the main class of entries which represents the
C        "final state" given by the generator.
C    = 2 : an entry which has decayed or fragmented and therefore
C        is not appearing in the final state, but is retained for
C        event history information.
C    = 3 : a documentation line, defined separately from the event
C        history. (incoming reacting
C        particles, etc.)
C    = 4 - 10 : undefined, but reserved for future standards.
C    = 11 - 20 : at the disposal of each model builder for constructs
C        specific to his program, but equivalent to a null line in the
C        context of any other program. One example is the cone defining
C        vector of HERWIG, another cluster or event axes of the JETSET
C        analysis routines.
C    = 21 - : at the disposal of users, in particular for event tracking
C        in the detector.
C
C IDHKK(IHKK) : particle identity, according to the Particle Data Group
C    standard.
C
C JMOHKK(1,IHKK) : pointer to the position where the mother is stored.
C    The value is 0 for initial entries.
C
C JMOHKK(2,IHKK) : pointer to position of second mother. Normally only
C    one mother exist, in which case the value 0 is used. In cluster
C    fragmentation models, the two mothers would correspond to the q
C    and qbar which join to form a cluster. In string fragmentation,
C    the two mothers of a particle produced in the fragmentation would
C    be the two endpoints of the string (with the range in between
C    implied).
C
C JDAHKK(1,IHKK) : pointer to the position of the first daughter. If an
C    entry has not decayed, this is 0.
C
C JDAHKK(2,IHKK) : pointer to the position of the last daughter. If an
C    entry has not decayed, this is 0. It is assumed that the daughters
C    of a particle (or cluster or string) are stored sequentially, so
C    that the whole range JDAHKK(1,IHKK) - JDAHKK(2,IHKK) contains
C    daughters. Even in cases where only one daughter is defined (e.g.
C    K0 -> K0S) both values should be defined, to make for a uniform
C    approach in terms of loop constructions.
C
C PHKK(1,IHKK) : momentum in the x direction, in GeV/c.
C
C PHKK(2,IHKK) : momentum in the y direction, in GeV/c.
C
C PHKK(3,IHKK) : momentum in the z direction, in GeV/c.
C
C PHKK(4,IHKK) : energy, in GeV.
C
C PHKK(5,IHKK) : mass, in GeV/c**2. For spacelike partons, it is allowed
C    to use a negative mass, according to PHKK(5,IHKK) = -sqrt(-m**2).
C
C VHKK(1,IHKK) : production vertex x position, in mm.
C
C VHKK(2,IHKK) : production vertex y position, in mm.
C
C VHKK(3,IHKK) : production vertex z position, in mm.
C
C VHKK(4,IHKK) : production time, in mm/c (= 3.33*10**(-12) s).
C********************************************************************
*KEEP,DELP.
      COMMON /DELP/ DELPX,DELPY,DELPZ,DELPE
*KEEP,TANUIN.
      COMMON /TANUIN/ TASUMA,TASUBI,TABI,TAMASU,TAMA,TAIMMA
*KEEP,NUCC.
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEND.
      DATA ICHECK/0/
      HELP=LOG10(EPN)
      PHELP=0.
      IF(HELP.GT.5.D0)PHELP=HELP-5.
      PTHELP=12.+PHELP*5.
      IREJ=0
      IREJJ=0
      PX=0.
      PY=0.
      PZ=0.
      PE=0.
      EEXT=0.
      EEXP=0.
      EEE1=0.
      EEEM1=0.
      EE1001=0.
      DO 10 I=1,NHKK
C                               Projectiles
        IF(ISTHKK(I).EQ.11.OR.ISTHKK(I).EQ.13) THEN
          GAM=EPN/PHKK(5,I)
          BGAM=PPN/PHKK(5,I)
          PX=PX - PHKK(1,I)
          PY=PY - PHKK(2,I)
          PZ=PZ - GAM*PHKK(3,I) - BGAM*PHKK(4,I)
          PE=PE - GAM*PHKK(4,I) - BGAM*PHKK(3,I)
        ENDIF
C                            Target
        IF(ISTHKK(I).EQ.12.OR.ISTHKK(I).EQ.14) THEN
          PX=PX - PHKK(1,I)
          PY=PY - PHKK(2,I)
          PZ=PZ - PHKK(3,I)
          PE=PE - PHKK(4,I)
        ENDIF
*                                 sum final state momenta
        IF(ISTHKK(I).EQ.1) THEN
          PX=PX + PHKK(1,I)
          PY=PY + PHKK(2,I)
          PZ=PZ + PHKK(3,I)
          PE=PE + PHKK(4,I)
        ENDIF
C                      noninteracting  Projectiles
        IF(ISTHKK(I).EQ.13.AND.JDAHKK(2,I).EQ.0) THEN
          GAM=EPN/PHKK(5,I)
          BGAM=PPN/PHKK(5,I)
          PX=PX + PHKK(1,I)
          PY=PY + PHKK(2,I)
          PZ=PZ + GAM*PHKK(3,I) + BGAM*PHKK(4,I)
          PE=PE + GAM*PHKK(4,I) + BGAM*PHKK(3,I)
        ENDIF
        IF(ISTHKK(I).EQ.14.AND.JDAHKK(2,I).EQ.0) THEN
C                      noninteracting  Targets
          PX=PX + PHKK(1,I)
          PY=PY + PHKK(2,I)
          PZ=PZ + PHKK(3,I)
          PE=PE + PHKK(4,I)
        ENDIF
        IF(ISTHKK(I).EQ.16) THEN
          IMO=JMOHKK(1,I)
          PX=PX + PHKK(1,I)
          PY=PY + PHKK(2,I)
          PZ=PZ + PHKK(3,I)
          PE=PE + PHKK(4,I)
          EEXT=EEXT + PHKK(4,I) - PHKK(4,IMO)
        ENDIF
        IF(ISTHKK(I).EQ.15) THEN
          IMO=JMOHKK(1,I)
          GAM=EPN/PHKK(5,I)
          BGAM=PPN/PHKK(5,I)
          PX=PX + PHKK(1,I)
          PY=PY + PHKK(2,I)
          PZ=PZ + GAM*PHKK(3,I) + BGAM*PHKK(4,I)
          PE=PE + GAM*PHKK(4,I) + BGAM*PHKK(3,I)
          EEXT=EEXT + PHKK(4,I) - PHKK(4,IMO)
        ENDIF
        IF(ISTHKK(I).EQ.1) THEN
	  EEE1=EEE1+PHKK(4,I)
        ENDIF
        IF(ISTHKK(I).EQ.-1) THEN
	  EEEM1=EEEM1+PHKK(4,I)
        ENDIF
        IF(ISTHKK(I).EQ.1001) THEN
	  EE1001=EE1001+PHKK(4,I)
        ENDIF
   10 CONTINUE
      EEE=EEE1+EEEM1+EE1001
      EPNTO=EPN*IP
      AEEE=EEE/EPNTO
      EEEE=EEE-EPN*IP
      AEEEE=EEE/EPN-IP
      AEEE1=EEE1/EPNTO
      AEEEM1=EEEM1/EPNTO
      AEE101=EE1001/EPNTO
      AIP=1
      AIT=IT
      AITZ=ITZ
      AIP=AIP+(AIT*AMUAMU+1.D-3*ENERGY(AIT,AITZ))/EPNTO
      DELLE=ABS(AIP-AEEE)
      ELLE=DELLE*EPNTO
      TOLE=0.030
C     TOLE=0.012
      IF(IT.LE.50)THEN
        IF(IT.EQ.IP)TOLE=0.02
C       IF(IT.EQ.IP)TOLE=0.05
      ENDIF
      IF(DELLE.GE.TOLE)IREJ=1
      IF(IREJ.EQ.1)THEN
        ICHECK=ICHECK+1
	IF(ICHECK.LE.100)THEN
          WRITE(6,'(A,I5,E10.3,5F10.4)')
     *    ' IP,EPN,AEEE,AEEEE,AEEE1,AEEEM1,AEE101:',
     *    IP,EPN,AEEE,AEEEE,AEEE1,AEEEM1,AEE101
          WRITE(6,'(A,I5,E10.3,7E12.4)')
     *    ' IP,EPN,EEE,EEEE,EEE1,EEEM1,EE1001,DELLE,ELLE:',
     *    IP,EPN,EEE,EEEE,EEE1,EEEM1,EE1001,DELLE,ELLE
        ENDIF
      ENDIF
C     PX=PX + DELPX
C     PY=PY + DELPY
C     PZ=PZ + DELPZ
C     PE=PE + DELPE
C     IF(IPRI.GT.1) THEN
C       IF (ABS(PX).GT.PTHELP.OR. ABS(PY).GT.PTHELP.OR. 
C    * 	ABS(PZ)/EPN.GT.0.025*IP.
C    +  OR. ABS(PE)/EPN.GT.0.025*IP) THEN
C         IREJJ=1
C         ICHECK=ICHECK+1
C         IF(ICHECK.LE.500.AND.IREJJ.EQ.1)THEN 
C     WRITE(6,1000) PX,PY,PZ,PE,EEXT,EEXP,DELPX,DELPY,DELPZ,DELPE,IORIG
C         ENDIF
 1000 FORMAT(' CHECKF: PX,PY,PZ,PE,EEXT,EEXP',2F7.3,2E12.3,2F7.3
     * / 8X,' DELPX/Y/Z/E',4F7.3,I10,' IORIG')
C     WRITE(6,'(8X,A,6F8.3)') ' TASUMA,TASUBI,TABI,TAMASU,TAMA,TAIMMA',
C    +TASUMA,TASUBI,TABI,TAMASU,TAMA,TAIMMA
 
      IF(IPRI.GT.1) THEN
          DO 20 I=1,NHKK
            IF(ISTHKK(I).EQ.11) THEN
              WRITE(6,1010) I,ISTHKK(I),IDHKK(I),JMOHKK(1,I),JMOHKK
     +        (2,I), JDAHKK(1,I),JDAHKK(2,I),(PHKK(KHKK,I),KHKK=1,5),
     +        (VHKK(KHKK,I),KHKK=1,4)
 
            ENDIF
            IF(ISTHKK(I).EQ.12) THEN
              WRITE(6,1010) I,ISTHKK(I),IDHKK(I),JMOHKK(1,I),JMOHKK
     +        (2,I), JDAHKK(1,I),JDAHKK(2,I),(PHKK(KHKK,I),KHKK=1,5),
     +        (VHKK(KHKK,I),KHKK=1,4)
 
            ENDIF
            IF(ISTHKK(I).EQ.1) THEN
              WRITE(6,1010) I,ISTHKK(I),IDHKK(I),JMOHKK(1,I),JMOHKK
     +        (2,I), JDAHKK(1,I),JDAHKK(2,I),(PHKK(KHKK,I),KHKK=1,5),
     +        (VHKK(KHKK,I),KHKK=1,4)
 
 1010 FORMAT (I6,I4,5I6,9(1PE10.2))
            ENDIF
            IF(ISTHKK(I).EQ.16) THEN
              IMO=JMOHKK(1,I)
              WRITE(6,1010) I,ISTHKK(I),IDHKK(I),JMOHKK(1,I),JMOHKK
     +        (2,I), JDAHKK(1,I),JDAHKK(2,I),(PHKK(KHKK,I),KHKK=1,5),
     +        (VHKK(KHKK,I),KHKK=1,4)
 
            ENDIF
   20     CONTINUE
        ENDIF
C     ENDIF
      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE CHECKN(EPN,PPN,IREJ,IORIG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
*KEEP,HKKEVT.
c     INCLUDE (HKKEVT)
      PARAMETER (AMUAMU=0.93149432D0)
      PARAMETER (NMXHKK= 89998)
c     PARAMETER (NMXHKK=25000)
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK), JMOHKK
     +(2,NMXHKK),JDAHKK(2,NMXHKK), PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK
     +(4,NMXHKK)
C
C                       WHKK(4,NMXHKK) GIVES POSITIONS AND TIMES IN
C                       PROJECTILE FRAME, THE CHAINS ARE CREATED ON
C                       THE POSITIONS OF THE PROJECTILE NUCLEONS
C                       IN THE PROJECTILE FRAME (TARGET NUCLEONS IN
C                       TARGET FRAME) BOTH POSITIONS ARE THREFORE NOT
C                       COMPLETELY CONSISTENT. THE TIMES IN THE
C                       PROJECTILE FRAME HOWEVER ARE OBTAINED BY
C                       LORENTZ TRANSFORMING FROM THE LAB SYSTEM.
C
C  Based on the proposed standard COMMON block (Sjostrand Memo 17.3,89)
C
C NMXHKK: maximum numbers of entries (partons/particles) that can be
C    stored in the commonblock.
C
C NHKK: the actual number of entries stored in current event. These are
C    found in the first NHKK positions of the respective arrays below.
C    Index IHKK, 1 <= IHKK <= NHKK, is below used to denote a given
C    entry.
C
C ISTHKK(IHKK): status code for entry IHKK, with following meanings:
C    = 0 : null entry.
C    = 1 : an existing entry, which has not decayed or fragmented.
C        This is the main class of entries which represents the
C        "final state" given by the generator.
C    = 2 : an entry which has decayed or fragmented and therefore
C        is not appearing in the final state, but is retained for
C        event history information.
C    = 3 : a documentation line, defined separately from the event
C        history. (incoming reacting
C        particles, etc.)
C    = 4 - 10 : undefined, but reserved for future standards.
C    = 11 - 20 : at the disposal of each model builder for constructs
C        specific to his program, but equivalent to a null line in the
C        context of any other program. One example is the cone defining
C        vector of HERWIG, another cluster or event axes of the JETSET
C        analysis routines.
C    = 21 - : at the disposal of users, in particular for event tracking
C        in the detector.
C
C IDHKK(IHKK) : particle identity, according to the Particle Data Group
C    standard.
C
C JMOHKK(1,IHKK) : pointer to the position where the mother is stored.
C    The value is 0 for initial entries.
C
C JMOHKK(2,IHKK) : pointer to position of second mother. Normally only
C    one mother exist, in which case the value 0 is used. In cluster
C    fragmentation models, the two mothers would correspond to the q
C    and qbar which join to form a cluster. In string fragmentation,
C    the two mothers of a particle produced in the fragmentation would
C    be the two endpoints of the string (with the range in between
C    implied).
C
C JDAHKK(1,IHKK) : pointer to the position of the first daughter. If an
C    entry has not decayed, this is 0.
C
C JDAHKK(2,IHKK) : pointer to the position of the last daughter. If an
C    entry has not decayed, this is 0. It is assumed that the daughters
C    of a particle (or cluster or string) are stored sequentially, so
C    that the whole range JDAHKK(1,IHKK) - JDAHKK(2,IHKK) contains
C    daughters. Even in cases where only one daughter is defined (e.g.
C    K0 -> K0S) both values should be defined, to make for a uniform
C    approach in terms of loop constructions.
C
C PHKK(1,IHKK) : momentum in the x direction, in GeV/c.
C
C PHKK(2,IHKK) : momentum in the y direction, in GeV/c.
C
C PHKK(3,IHKK) : momentum in the z direction, in GeV/c.
C
C PHKK(4,IHKK) : energy, in GeV.
C
C PHKK(5,IHKK) : mass, in GeV/c**2. For spacelike partons, it is allowed
C    to use a negative mass, according to PHKK(5,IHKK) = -sqrt(-m**2).
C
C VHKK(1,IHKK) : production vertex x position, in mm.
C
C VHKK(2,IHKK) : production vertex y position, in mm.
C
C VHKK(3,IHKK) : production vertex z position, in mm.
C
C VHKK(4,IHKK) : production time, in mm/c (= 3.33*10**(-12) s).
C********************************************************************
*KEEP,DELP.
      COMMON /DELP/ DELPX,DELPY,DELPZ,DELPE
*KEEP,TANUIN.
      COMMON /TANUIN/ TASUMA,TASUBI,TABI,TAMASU,TAMA,TAIMMA
*KEEP,NUCC.
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEND.
      DATA ICHECK/0/
      HELP=LOG10(EPN)
      PHELP=0.
      IF(HELP.GT.5.D0)PHELP=HELP-5.D0
      PTHELP=12.D0+PHELP*5.D0
      IREJ=0
      IREJJ=0
      PX=0.D0
      PY=0.D0
      PZ=0.D0
      PE=0.D0
      EEXT=0.D0
      EEXP=0.D0
      EEE1=0.D0
      EEEM1=0.D0
      EE1001=0.D0
      PZ1=0.D0
      PZM1=0.D0
      PZ1001=0.D0
      PX1=0.D0
      PXM1=0.0
      PX1001=0.D0
      DO 10 I=1,NHKK
C                               Projectiles
        IF(ISTHKK(I).EQ.11.OR.ISTHKK(I).EQ.13) THEN
          GAM=EPN/PHKK(5,I)
          BGAM=PPN/PHKK(5,I)
          PX=PX - PHKK(1,I)
          PY=PY - PHKK(2,I)
          PZ=PZ - GAM*PHKK(3,I) - BGAM*PHKK(4,I)
          PE=PE - GAM*PHKK(4,I) - BGAM*PHKK(3,I)
        ENDIF
C                            Target
        IF(ISTHKK(I).EQ.12.OR.ISTHKK(I).EQ.14) THEN
          PX=PX - PHKK(1,I)
          PY=PY - PHKK(2,I)
          PZ=PZ - PHKK(3,I)
          PE=PE - PHKK(4,I)
        ENDIF
*                                 sum final state momenta
        IF(ISTHKK(I).EQ.1) THEN
          PX=PX + PHKK(1,I)
          PY=PY + PHKK(2,I)
          PZ=PZ + PHKK(3,I)
          PE=PE + PHKK(4,I)
        ENDIF
C                      noninteracting  Projectiles
        IF(ISTHKK(I).EQ.13.AND.JDAHKK(2,I).EQ.0) THEN
          GAM=EPN/PHKK(5,I)
          BGAM=PPN/PHKK(5,I)
          PX=PX + PHKK(1,I)
          PY=PY + PHKK(2,I)
          PZ=PZ + GAM*PHKK(3,I) + BGAM*PHKK(4,I)
          PE=PE + GAM*PHKK(4,I) + BGAM*PHKK(3,I)
        ENDIF
        IF(ISTHKK(I).EQ.14.AND.JDAHKK(2,I).EQ.0) THEN
C                      noninteracting  Targets
          PX=PX + PHKK(1,I)
          PY=PY + PHKK(2,I)
          PZ=PZ + PHKK(3,I)
          PE=PE + PHKK(4,I)
        ENDIF
        IF(ISTHKK(I).EQ.16) THEN
          IMO=JMOHKK(1,I)
          PX=PX + PHKK(1,I)
          PY=PY + PHKK(2,I)
          PZ=PZ + PHKK(3,I)
          PE=PE + PHKK(4,I)
          EEXT=EEXT + PHKK(4,I) - PHKK(4,IMO)
        ENDIF
        IF(ISTHKK(I).EQ.15) THEN
          IMO=JMOHKK(1,I)
          GAM=EPN/PHKK(5,I)
          BGAM=PPN/PHKK(5,I)
          PX=PX + PHKK(1,I)
          PY=PY + PHKK(2,I)
          PZ=PZ + GAM*PHKK(3,I) + BGAM*PHKK(4,I)
          PE=PE + GAM*PHKK(4,I) + BGAM*PHKK(3,I)
          EEXT=EEXT + PHKK(4,I) - PHKK(4,IMO)
        ENDIF
        IF(ISTHKK(I).EQ.1) THEN
	  EEE1=EEE1+PHKK(4,I)
	  PZ1=PZ1+PHKK(3,I)
	  PX1=PX1+PHKK(1,I)
        ENDIF
        IF(ISTHKK(I).EQ.-1) THEN
	  EEEM1=EEEM1+PHKK(4,I)
	  PZM1=PZM1+PHKK(3,I)
	  PXM1=PXM1+PHKK(1,I)
        ENDIF
        IF(ISTHKK(I).EQ.1001) THEN
	  EE1001=EE1001+PHKK(4,I)
	  PZ1001=PZ1001+PHKK(3,I)
	  PX1001=PX1001+PHKK(1,I)
        ENDIF
   10 CONTINUE
      EEE=EEE1+EEEM1+EE1001
      PZPZ=PZ1+PZM1+PZ1001
      PXPX=PX1+PXM1+PX1001
C--------------------------------------------------------------
C
C                    patch to correct pz of residual nuclei
C
C--------------------------------------------------------------
      DELPZ=PPN-PZPZ
      PZPZ=PZPZ+DELPZ
      PZ1001=PZ1001+DELPZ
      EE1001=0.D0
      DO 101 I=1,NHKK
        IF(ISTHKK(I).EQ.1001) THEN
	  PHKK(3,I)=PHKK(3,I)+DELPZ
	  PHKK(4,I)=SQRT(PHKK(1,I)**2+PHKK(2,I)**2+PHKK(3,I)**2
     *                   +PHKK(5,I)**2)
	  EE1001=EE1001+PHKK(4,I)
        ENDIF
  101 CONTINUE
      EEE=EEE1+EEEM1+EE1001
C--------------------------------------------------------------
      AIP=1
      EPNTO=EPN*AIP
      EEEE=EEE-EPN*AIP
      AIP=1
      AIT=IT
      AITZ=ITZ
      BIP=EPN+(AIT*AMUAMU+1.D-3*ENERGY(AIT,AITZ))
      BMI=1.D-3*ENERGY(AIT,AITZ)
      DELLE=ABS(BIP-EEE)
C     TOLE=EPN/450000.D0
C     TOLE=EPN/2500.D0
      TOLE=0.16D0
      IF(DELLE.GE.TOLE)IREJ=1
      IF(IREJ.EQ.1)THEN
        ICHECK=ICHECK+1
	IF(ICHECK.LE.20)THEN
          WRITE(6,'(A,I5,E10.3,4F10.4)')
     *    ' IP,EPN,PXPX,PX1,PXM1,PX1001:',
     *    IP,EPN,PXPX,PX1,PXM1,PX1001
          WRITE(6,'(A,I5,E10.3,6F10.4)')
     *    ' IP,PPN,PZPZ,PZ1,PZM1,PZ1001,BIP,BMI:',
     *    IP,PPN,PZPZ,PZ1,PZM1,PZ1001,BIP,BMI
          WRITE(6,'(A,I5,E10.3,5E12.4)')
     *    ' IP,EPN,EEE,EEE1,EEEM1,EE1001,DELLE:',
     *    IP,EPN,EEE,EEE1,EEEM1,EE1001,DELLE
        ENDIF
      ENDIF
C     PX=PX + DELPX
C     PY=PY + DELPY
C     PZ=PZ + DELPZ
C     PE=PE + DELPE
C     IF(IPRI.GT.1) THEN
C       IF (ABS(PX).GT.PTHELP.OR. ABS(PY).GT.PTHELP.OR. 
C    * 	ABS(PZ)/EPN.GT.0.025D0*IP.
C    +  OR. ABS(PE)/EPN.GT.0.025D0*IP) THEN
C         IREJJ=1
C         ICHECK=ICHECK+1
C         IF(ICHECK.LE.500.AND.IREJJ.EQ.1)THEN 
C     WRITE(6,1000) PX,PY,PZ,PE,EEXT,EEXP,DELPX,DELPY,DELPZ,DELPE,IORIG
C         ENDIF
 1000 FORMAT(' CHECKF: PX,PY,PZ,PE,EEXT,EEXP',2F7.3,2E12.3,2F7.3
     * / 8X,' DELPX/Y/Z/E',4F7.3,I10,' IORIG')
C     WRITE(6,'(8X,A,6F8.3)') ' TASUMA,TASUBI,TABI,TAMASU,TAMA,TAIMMA',
C    +TASUMA,TASUBI,TABI,TAMASU,TAMA,TAIMMA
 
      IF(IPRI.GT.1) THEN
          DO 20 I=1,NHKK
            IF(ISTHKK(I).EQ.11) THEN
              WRITE(6,1010) I,ISTHKK(I),IDHKK(I),JMOHKK(1,I),JMOHKK
     +        (2,I), JDAHKK(1,I),JDAHKK(2,I),(PHKK(KHKK,I),KHKK=1,5),
     +        (VHKK(KHKK,I),KHKK=1,4)
 
            ENDIF
            IF(ISTHKK(I).EQ.12) THEN
              WRITE(6,1010) I,ISTHKK(I),IDHKK(I),JMOHKK(1,I),JMOHKK
     +        (2,I), JDAHKK(1,I),JDAHKK(2,I),(PHKK(KHKK,I),KHKK=1,5),
     +        (VHKK(KHKK,I),KHKK=1,4)
 
            ENDIF
            IF(ISTHKK(I).EQ.1) THEN
              WRITE(6,1010) I,ISTHKK(I),IDHKK(I),JMOHKK(1,I),JMOHKK
     +        (2,I), JDAHKK(1,I),JDAHKK(2,I),(PHKK(KHKK,I),KHKK=1,5),
     +        (VHKK(KHKK,I),KHKK=1,4)
 
 1010 FORMAT (I6,I4,5I6,9(1PE10.2))
            ENDIF
            IF(ISTHKK(I).EQ.16) THEN
              IMO=JMOHKK(1,I)
              WRITE(6,1010) I,ISTHKK(I),IDHKK(I),JMOHKK(1,I),JMOHKK
     +        (2,I), JDAHKK(1,I),JDAHKK(2,I),(PHKK(KHKK,I),KHKK=1,5),
     +        (VHKK(KHKK,I),KHKK=1,4)
 
            ENDIF
   20     CONTINUE
        ENDIF
C     ENDIF
      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE CHECKO(EPN,PPN,IREJ,IORIG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
*KEEP,HKKEVT.
c     INCLUDE (HKKEVT)
      PARAMETER (AMUAMU=0.93149432D0)
      PARAMETER (NMXHKK= 89998)
c     PARAMETER (NMXHKK=25000)
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK), JMOHKK
     +(2,NMXHKK),JDAHKK(2,NMXHKK), PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK
     +(4,NMXHKK)
C
C                       WHKK(4,NMXHKK) GIVES POSITIONS AND TIMES IN
C                       PROJECTILE FRAME, THE CHAINS ARE CREATED ON
C                       THE POSITIONS OF THE PROJECTILE NUCLEONS
C                       IN THE PROJECTILE FRAME (TARGET NUCLEONS IN
C                       TARGET FRAME) BOTH POSITIONS ARE THREFORE NOT
C                       COMPLETELY CONSISTENT. THE TIMES IN THE
C                       PROJECTILE FRAME HOWEVER ARE OBTAINED BY
C                       LORENTZ TRANSFORMING FROM THE LAB SYSTEM.
C
C  Based on the proposed standard COMMON block (Sjostrand Memo 17.3,89)
C
C NMXHKK: maximum numbers of entries (partons/particles) that can be
C    stored in the commonblock.
C
C NHKK: the actual number of entries stored in current event. These are
C    found in the first NHKK positions of the respective arrays below.
C    Index IHKK, 1 <= IHKK <= NHKK, is below used to denote a given
C    entry.
C
C ISTHKK(IHKK): status code for entry IHKK, with following meanings:
C    = 0 : null entry.
C    = 1 : an existing entry, which has not decayed or fragmented.
C        This is the main class of entries which represents the
C        "final state" given by the generator.
C    = 2 : an entry which has decayed or fragmented and therefore
C        is not appearing in the final state, but is retained for
C        event history information.
C    = 3 : a documentation line, defined separately from the event
C        history. (incoming reacting
C        particles, etc.)
C    = 4 - 10 : undefined, but reserved for future standards.
C    = 11 - 20 : at the disposal of each model builder for constructs
C        specific to his program, but equivalent to a null line in the
C        context of any other program. One example is the cone defining
C        vector of HERWIG, another cluster or event axes of the JETSET
C        analysis routines.
C    = 21 - : at the disposal of users, in particular for event tracking
C        in the detector.
C
C IDHKK(IHKK) : particle identity, according to the Particle Data Group
C    standard.
C
C JMOHKK(1,IHKK) : pointer to the position where the mother is stored.
C    The value is 0 for initial entries.
C
C JMOHKK(2,IHKK) : pointer to position of second mother. Normally only
C    one mother exist, in which case the value 0 is used. In cluster
C    fragmentation models, the two mothers would correspond to the q
C    and qbar which join to form a cluster. In string fragmentation,
C    the two mothers of a particle produced in the fragmentation would
C    be the two endpoints of the string (with the range in between
C    implied).
C
C JDAHKK(1,IHKK) : pointer to the position of the first daughter. If an
C    entry has not decayed, this is 0.
C
C JDAHKK(2,IHKK) : pointer to the position of the last daughter. If an
C    entry has not decayed, this is 0. It is assumed that the daughters
C    of a particle (or cluster or string) are stored sequentially, so
C    that the whole range JDAHKK(1,IHKK) - JDAHKK(2,IHKK) contains
C    daughters. Even in cases where only one daughter is defined (e.g.
C    K0 -> K0S) both values should be defined, to make for a uniform
C    approach in terms of loop constructions.
C
C PHKK(1,IHKK) : momentum in the x direction, in GeV/c.
C
C PHKK(2,IHKK) : momentum in the y direction, in GeV/c.
C
C PHKK(3,IHKK) : momentum in the z direction, in GeV/c.
C
C PHKK(4,IHKK) : energy, in GeV.
C
C PHKK(5,IHKK) : mass, in GeV/c**2. For spacelike partons, it is allowed
C    to use a negative mass, according to PHKK(5,IHKK) = -sqrt(-m**2).
C
C VHKK(1,IHKK) : production vertex x position, in mm.
C
C VHKK(2,IHKK) : production vertex y position, in mm.
C
C VHKK(3,IHKK) : production vertex z position, in mm.
C
C VHKK(4,IHKK) : production time, in mm/c (= 3.33*10**(-12) s).
C********************************************************************
      COMMON /ZENTRA/ ICENTR
*KEEP,DELP.
      COMMON /DELP/ DELPX,DELPY,DELPZ,DELPE
*KEEP,TANUIN.
      COMMON /TANUIN/ TASUMA,TASUBI,TABI,TAMASU,TAMA,TAIMMA
*KEEP,NUCC.
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEND.
      DATA ICHECK/0/
      HELP=LOG10(EPN)
      PHELP=0.
      IF(HELP.GT.5.D0)PHELP=HELP-5.
      PTHELP=12.+PHELP*5.
      IREJ=0
      IREJJ=0
      PX=0.
      PY=0.
      PZ=0.
      PE=0.
      EEXT=0.
      EEXP=0.
      EEE1=0.
      EEEM1=0.
      EE1001=0.
      DO 10 I=1,NHKK
C                               Projectiles
        IF(ISTHKK(I).EQ.11.OR.ISTHKK(I).EQ.13) THEN
          GAM=EPN/PHKK(5,I)
          BGAM=PPN/PHKK(5,I)
          PX=PX - PHKK(1,I)
          PY=PY - PHKK(2,I)
          PZ=PZ - GAM*PHKK(3,I) - BGAM*PHKK(4,I)
          PE=PE - GAM*PHKK(4,I) - BGAM*PHKK(3,I)
        ENDIF
C                            Target
        IF(ISTHKK(I).EQ.12.OR.ISTHKK(I).EQ.14) THEN
          PX=PX - PHKK(1,I)
          PY=PY - PHKK(2,I)
          PZ=PZ - PHKK(3,I)
          PE=PE - PHKK(4,I)
        ENDIF
*                                 sum final state momenta
        IF(ISTHKK(I).EQ.1) THEN
          PX=PX + PHKK(1,I)
          PY=PY + PHKK(2,I)
          PZ=PZ + PHKK(3,I)
          PE=PE + PHKK(4,I)
        ENDIF
C                      noninteracting  Projectiles
        IF(ISTHKK(I).EQ.13.AND.JDAHKK(2,I).EQ.0) THEN
          GAM=EPN/PHKK(5,I)
          BGAM=PPN/PHKK(5,I)
          PX=PX + PHKK(1,I)
          PY=PY + PHKK(2,I)
          PZ=PZ + GAM*PHKK(3,I) + BGAM*PHKK(4,I)
          PE=PE + GAM*PHKK(4,I) + BGAM*PHKK(3,I)
        ENDIF
        IF(ISTHKK(I).EQ.14.AND.JDAHKK(2,I).EQ.0) THEN
C                      noninteracting  Targets
          PX=PX + PHKK(1,I)
          PY=PY + PHKK(2,I)
          PZ=PZ + PHKK(3,I)
          PE=PE + PHKK(4,I)
        ENDIF
        IF(ISTHKK(I).EQ.16) THEN
          IMO=JMOHKK(1,I)
          PX=PX + PHKK(1,I)
          PY=PY + PHKK(2,I)
          PZ=PZ + PHKK(3,I)
          PE=PE + PHKK(4,I)
          EEXT=EEXT + PHKK(4,I) - PHKK(4,IMO)
        ENDIF
        IF(ISTHKK(I).EQ.15) THEN
          IMO=JMOHKK(1,I)
          GAM=EPN/PHKK(5,I)
          BGAM=PPN/PHKK(5,I)
          PX=PX + PHKK(1,I)
          PY=PY + PHKK(2,I)
          PZ=PZ + GAM*PHKK(3,I) + BGAM*PHKK(4,I)
          PE=PE + GAM*PHKK(4,I) + BGAM*PHKK(3,I)
          EEXT=EEXT + PHKK(4,I) - PHKK(4,IMO)
        ENDIF
        IF(ISTHKK(I).EQ.1) THEN
	  EEE1=EEE1+PHKK(3,I)
        ENDIF
   10 CONTINUE
      EEE=EEE1
      EPNTO=PPN*IP
      AEEE=EEE/EPNTO
      EEEE=EEE-PPN*IP
      AEEEE=EEE/PPN-IP
      AEEE1=EEE1/EPNTO
      AIP=1
      AIT=IT
      AITZ=ITZ
      AIP=AIP
      DELLE=ABS(AIP-AEEE)
      ELLE=DELLE*EPNTO
C     IF(DELLE.GE.0.025)IREJ=1
C     IF(IREJ.EQ.1)
C    * WRITE(6,'(A,I5,E10.3,3F10.4)')
C    *' IP,EPN,AEEE,AEEEE,AEEE1:',
C    * IP,EPN,AEEE,AEEEE,AEEE1
C     IF(IREJ.EQ.1)
C    * WRITE(6,'(A,I5,E10.3,5F10.4)')
C    *' IP,EPN,EEE,EEEE,EEE1,DELLE,ELLE:',
C    * IP,EPN,EEE,EEEE,EEE1,DELLE,ELLE
      PX=PX + DELPX
      PY=PY + DELPY
      PZ=PZ + DELPZ
      PE=PE + DELPE
      TOLE=0.025D0*IP
      IF(IP.EQ.IT.AND.IT.GT.1)TOLE=0.05D0*IP
C     IF(ICENTR.EQ.1)TOLE=TOLE*2.
      IF(EPN.LE.5.D0)TOLE=3.D0*TOLE
C     IF(IPRI.GT.1) THEN
        IF (ABS(PX).GT.PTHELP.OR. ABS(PY).GT.PTHELP.OR. 
     * 	ABS(PZ)/EPN.GT.TOLE.
     +  OR. ABS(PE)/EPN.GT.TOLE) THEN
          IREJ=1
          ICHECK=ICHECK+1
          IF(ICHECK.LE.50.AND.IREJ.EQ.1)THEN 
      WRITE(6,1000) PX,PY,PZ,PE,EEXT,EEXP,DELPX,DELPY,DELPZ,DELPE,IORIG
          ENDIF
 1000 FORMAT(' CHECKO: PX,PY,PZ,PE,EEXT,EEXP',2F7.3,2E12.3,2F7.3
     * / 8X,' DELPX/Y/Z/E',4F7.3,I10,' IORIG')
C     WRITE(6,'(8X,A,6F8.3)') ' TASUMA,TASUBI,TABI,TAMASU,TAMA,TAIMMA',
C    +TASUMA,TASUBI,TABI,TAMASU,TAMA,TAIMMA
          IF(IPRI.GE.1)THEN 
      WRITE(6,1000) PX,PY,PZ,PE,EEXT,EEXP,DELPX,DELPY,DELPZ,DELPE,IORIG
          ENDIF
      IF(IPRI.GT.1) THEN
          DO 20 I=1,NHKK
            IF(ISTHKK(I).EQ.11) THEN
              WRITE(6,1010) I,ISTHKK(I),IDHKK(I),JMOHKK(1,I),JMOHKK
     +        (2,I), JDAHKK(1,I),JDAHKK(2,I),(PHKK(KHKK,I),KHKK=1,5),
     +        (VHKK(KHKK,I),KHKK=1,4)
 
            ENDIF
            IF(ISTHKK(I).EQ.12) THEN
              WRITE(6,1010) I,ISTHKK(I),IDHKK(I),JMOHKK(1,I),JMOHKK
     +        (2,I), JDAHKK(1,I),JDAHKK(2,I),(PHKK(KHKK,I),KHKK=1,5),
     +        (VHKK(KHKK,I),KHKK=1,4)
 
            ENDIF
            IF(ISTHKK(I).EQ.1) THEN
              WRITE(6,1010) I,ISTHKK(I),IDHKK(I),JMOHKK(1,I),JMOHKK
     +        (2,I), JDAHKK(1,I),JDAHKK(2,I),(PHKK(KHKK,I),KHKK=1,5),
     +        (VHKK(KHKK,I),KHKK=1,4)
 
 1010 FORMAT (I6,I4,5I6,9(1PE10.2))
            ENDIF
            IF(ISTHKK(I).EQ.16) THEN
              IMO=JMOHKK(1,I)
              WRITE(6,1010) I,ISTHKK(I),IDHKK(I),JMOHKK(1,I),JMOHKK
     +        (2,I), JDAHKK(1,I),JDAHKK(2,I),(PHKK(KHKK,I),KHKK=1,5),
     +        (VHKK(KHKK,I),KHKK=1,4)
 
            ENDIF
   20     CONTINUE
        ENDIF
      ENDIF
          IF(IPRI.GE.1)THEN 
      WRITE(6,1000) PX,PY,PZ,PE,EEXT,EEXP,DELPX,DELPY,DELPZ,DELPE,IORIG
          ENDIF
      RETURN
      END
C
      SUBROUTINE CHECKE(EPN,PPN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
*KEEP,HKKEVT.
c     INCLUDE (HKKEVT)
      PARAMETER (NMXHKK= 89998)
c     PARAMETER (NMXHKK=25000)
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK), JMOHKK
     +(2,NMXHKK),JDAHKK(2,NMXHKK), PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK
     +(4,NMXHKK)
C
C                       WHKK(4,NMXHKK) GIVES POSITIONS AND TIMES IN
C                       PROJECTILE FRAME, THE CHAINS ARE CREATED ON
C                       THE POSITIONS OF THE PROJECTILE NUCLEONS
C                       IN THE PROJECTILE FRAME (TARGET NUCLEONS IN
C                       TARGET FRAME) BOTH POSITIONS ARE THREFORE NOT
C                       COMPLETELY CONSISTENT. THE TIMES IN THE
C                       PROJECTILE FRAME HOWEVER ARE OBTAINED BY
C                       LORENTZ TRANSFORMING FROM THE LAB SYSTEM.
C
C  Based on the proposed standard COMMON block (Sjostrand Memo 17.3,89)
C
C NMXHKK: maximum numbers of entries (partons/particles) that can be
C    stored in the commonblock.
C
C NHKK: the actual number of entries stored in current event. These are
C    found in the first NHKK positions of the respective arrays below.
C    Index IHKK, 1 <= IHKK <= NHKK, is below used to denote a given
C    entry.
C
C ISTHKK(IHKK): status code for entry IHKK, with following meanings:
C    = 0 : null entry.
C    = 1 : an existing entry, which has not decayed or fragmented.
C        This is the main class of entries which represents the
C        "final state" given by the generator.
C    = 2 : an entry which has decayed or fragmented and therefore
C        is not appearing in the final state, but is retained for
C        event history information.
C    = 3 : a documentation line, defined separately from the event
C        history. (incoming reacting
C        particles, etc.)
C    = 4 - 10 : undefined, but reserved for future standards.
C    = 11 - 20 : at the disposal of each model builder for constructs
C        specific to his program, but equivalent to a null line in the
C        context of any other program. One example is the cone defining
C        vector of HERWIG, another cluster or event axes of the JETSET
C        analysis routines.
C    = 21 - : at the disposal of users, in particular for event tracking
C        in the detector.
C
C IDHKK(IHKK) : particle identity, according to the Particle Data Group
C    standard.
C
C JMOHKK(1,IHKK) : pointer to the position where the mother is stored.
C    The value is 0 for initial entries.
C
C JMOHKK(2,IHKK) : pointer to position of second mother. Normally only
C    one mother exist, in which case the value 0 is used. In cluster
C    fragmentation models, the two mothers would correspond to the q
C    and qbar which join to form a cluster. In string fragmentation,
C    the two mothers of a particle produced in the fragmentation would
C    be the two endpoints of the string (with the range in between
C    implied).
C
C JDAHKK(1,IHKK) : pointer to the position of the first daughter. If an
C    entry has not decayed, this is 0.
C
C JDAHKK(2,IHKK) : pointer to the position of the last daughter. If an
C    entry has not decayed, this is 0. It is assumed that the daughters
C    of a particle (or cluster or string) are stored sequentially, so
C    that the whole range JDAHKK(1,IHKK) - JDAHKK(2,IHKK) contains
C    daughters. Even in cases where only one daughter is defined (e.g.
C    K0 -> K0S) both values should be defined, to make for a uniform
C    approach in terms of loop constructions.
C
C PHKK(1,IHKK) : momentum in the x direction, in GeV/c.
C
C PHKK(2,IHKK) : momentum in the y direction, in GeV/c.
C
C PHKK(3,IHKK) : momentum in the z direction, in GeV/c.
C
C PHKK(4,IHKK) : energy, in GeV.
C
C PHKK(5,IHKK) : mass, in GeV/c**2. For spacelike partons, it is allowed
C    to use a negative mass, according to PHKK(5,IHKK) = -sqrt(-m**2).
C
C VHKK(1,IHKK) : production vertex x position, in mm.
C
C VHKK(2,IHKK) : production vertex y position, in mm.
C
C VHKK(3,IHKK) : production vertex z position, in mm.
C
C VHKK(4,IHKK) : production time, in mm/c (= 3.33*10**(-12) s).
C********************************************************************
*KEEP,DELP.
      COMMON /DELP/ DELPX,DELPY,DELPZ,DELPE
*KEEP,TANUIN.
      COMMON /TANUIN/ TASUMA,TASUBI,TABI,TAMASU,TAMA,TAIMMA
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEND.
      PX=0.
      PY=0.
      PZ=0.
      PE=0.
      EEXT=0.
      EEXP=0.
      DO 10 I=1,NHKK
C                               Projectiles
        IF(ISTHKK(I).EQ.11.OR.ISTHKK(I).EQ.13) THEN
          GAM=EPN/PHKK(5,I)
          BGAM=PPN/PHKK(5,I)
          PX=PX - PHKK(1,I)
          PY=PY - PHKK(2,I)
          PZ=PZ - GAM*PHKK(3,I) - BGAM*PHKK(4,I)
          PE=PE - GAM*PHKK(4,I) - BGAM*PHKK(3,I)
        ENDIF
C                            Target
        IF(ISTHKK(I).EQ.12.OR.ISTHKK(I).EQ.14) THEN
          PX=PX - PHKK(1,I)
          PY=PY - PHKK(2,I)
          PZ=PZ - PHKK(3,I)
          PE=PE - PHKK(4,I)
        ENDIF
*                                 sum final state momenta
        IF(ISTHKK(I).EQ.1) THEN
          PX=PX + PHKK(1,I)
          PY=PY + PHKK(2,I)
          PZ=PZ + PHKK(3,I)
          PE=PE + PHKK(4,I)
        ENDIF
C                      noninteracting  Projectiles
        IF(ISTHKK(I).EQ.13.AND.JDAHKK(1,I).EQ.0) THEN
          GAM=EPN/PHKK(5,I)
          BGAM=PPN/PHKK(5,I)
          PX=PX + PHKK(1,I)
          PY=PY + PHKK(2,I)
          PZ=PZ + GAM*PHKK(3,I) + BGAM*PHKK(4,I)
          PE=PE + GAM*PHKK(4,I) + BGAM*PHKK(3,I)
        ENDIF
        IF(ISTHKK(I).EQ.14.AND.JDAHKK(1,I).EQ.0) THEN
C                      noninteracting  Targets
          PX=PX + PHKK(1,I)
          PY=PY + PHKK(2,I)
          PZ=PZ + PHKK(3,I)
          PE=PE + PHKK(4,I)
        ENDIF
        IF(ISTHKK(I).EQ.16) THEN
          IMO=JMOHKK(1,I)
          PX=PX + PHKK(1,I)
          PY=PY + PHKK(2,I)
          PZ=PZ + PHKK(3,I)
          PE=PE + PHKK(4,I)
          EEXT=EEXT + PHKK(4,I) - PHKK(4,IMO)
        ENDIF
        IF(ISTHKK(I).EQ.15) THEN
          IMO=JMOHKK(1,I)
          PX=PX + PHKK(1,I)
          PY=PY + PHKK(2,I)
          PZ=PZ + PHKK(3,I)
          PE=PE + PHKK(4,I)
          EEXT=EEXT + PHKK(4,I) - PHKK(4,IMO)
        ENDIF
   10 CONTINUE
      PX=PX + DELPX
      PY=PY + DELPY
      PZ=PZ + DELPZ
      PE=PE + DELPE
      WRITE(6,1000) PX,PY,PZ,PE,EEXT,EEXP,DELPX,DELPY,DELPZ,DELPE
 1000 FORMAT(' CHECKE: PX,PY,PZ,PE,EEXT,EEXP',6F7.3/ 8X,' DELPX/Y/Z/E',4
     +F7.3)
      WRITE(6,'(8X,A,6F8.3)') ' TASUMA,TASUBI,TABI,TAMASU,TAMA,TAIMMA',
     +TASUMA,TASUBI,TABI,TAMASU,TAMA,TAIMMA
      IF(IPRI.GT.1) THEN
        IF (ABS(PX).GT.0.004.OR. ABS(PY).GT.0.004.OR. ABS(PZ).GT.0.004.
     +  OR. ABS(PE).GT.0.004) THEN
 
 
          DO 20 I=1,NHKK
            IF(ISTHKK(I).EQ.11) THEN
              WRITE(6,1010) I,ISTHKK(I),IDHKK(I),JMOHKK(1,I),JMOHKK
     +        (2,I), JDAHKK(1,I),JDAHKK(2,I),(PHKK(KHKK,I),KHKK=1,5),
     +        (VHKK(KHKK,I),KHKK=1,4)
 
            ENDIF
            IF(ISTHKK(I).EQ.12) THEN
              WRITE(6,1010) I,ISTHKK(I),IDHKK(I),JMOHKK(1,I),JMOHKK
     +        (2,I), JDAHKK(1,I),JDAHKK(2,I),(PHKK(KHKK,I),KHKK=1,5),
     +        (VHKK(KHKK,I),KHKK=1,4)
 
            ENDIF
            IF(ISTHKK(I).EQ.1) THEN
              WRITE(6,1010) I,ISTHKK(I),IDHKK(I),JMOHKK(1,I),JMOHKK
     +        (2,I), JDAHKK(1,I),JDAHKK(2,I),(PHKK(KHKK,I),KHKK=1,5),
     +        (VHKK(KHKK,I),KHKK=1,4)
 
 1010 FORMAT (I6,I4,5I6,9(1PE10.2))
            ENDIF
            IF(ISTHKK(I).EQ.16) THEN
              IMO=JMOHKK(1,I)
              WRITE(6,1010) I,ISTHKK(I),IDHKK(I),JMOHKK(1,I),JMOHKK
     +        (2,I), JDAHKK(1,I),JDAHKK(2,I),(PHKK(KHKK,I),KHKK=1,5),
     +        (VHKK(KHKK,I),KHKK=1,4)
 
            ENDIF
   20     CONTINUE
        ENDIF
      ENDIF
      RETURN
      END
C######################################################################
C######################################################################
C######################################################################
*
C######################################################################
C######################################################################
C######################################################################
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C######################################################################
C######################################################################
C######################################################################
C######################################################################
C######################################################################
C######################################################################
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      DOUBLE PRECISION FUNCTION EBIND(IA,IZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C***
C   Binding energy for nuclei with mass number IA
C                              and atomic number IZ
C   from Shirokov & Yudin, Yad. Fizika, Nauka, Moskva 1972
C***
      DATA A1,A2,A3,A4,A5 /0.01575, 0.0178, 0.000710, 0.0237, 0.034/
C
C     WRITE (6,'(A,2I5)')' EBIND IA,IZ ',IA,IZ
      IF(IA.LE.1.OR.IZ.EQ.0)THEN
	EBIND=0
	RETURN
      ENDIF
      AA=IA
      EBIND = A1*AA - A2*AA**0. 666667- A3*IZ*IZ*AA**(-0.333333) - A4
     +*(IA-2*IZ)**2/AA
      IF (MOD(IA,2).EQ.1) THEN
        IA5=0
      ELSEIF (MOD(IZ,2).EQ.1) THEN
        IA5=1
      ELSE
        IA5=-1
      ENDIF
      EBIND=EBIND - IA5*A5*AA**(-0.75)
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE DEFAUL(EPN,PPN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
*---set default values for some parameters
*KEEP,HKKEVT.
c     INCLUDE (HKKEVT)
      PARAMETER (NMXHKK= 89998)
c     PARAMETER (NMXHKK=25000)
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK), JMOHKK
     +(2,NMXHKK),JDAHKK(2,NMXHKK), PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK
     +(4,NMXHKK)
C
C                       WHKK(4,NMXHKK) GIVES POSITIONS AND TIMES IN
C                       PROJECTILE FRAME, THE CHAINS ARE CREATED ON
C                       THE POSITIONS OF THE PROJECTILE NUCLEONS
C                       IN THE PROJECTILE FRAME (TARGET NUCLEONS IN
C                       TARGET FRAME) BOTH POSITIONS ARE THREFORE NOT
C                       COMPLETELY CONSISTENT. THE TIMES IN THE
C                       PROJECTILE FRAME HOWEVER ARE OBTAINED BY
C                       LORENTZ TRANSFORMING FROM THE LAB SYSTEM.
C
C  Based on the proposed standard COMMON block (Sjostrand Memo 17.3,89)
C
C NMXHKK: maximum numbers of entries (partons/particles) that can be
C    stored in the commonblock.
C
C NHKK: the actual number of entries stored in current event. These are
C    found in the first NHKK positions of the respective arrays below.
C    Index IHKK, 1 <= IHKK <= NHKK, is below used to denote a given
C    entry.
C
C ISTHKK(IHKK): status code for entry IHKK, with following meanings:
C    = 0 : null entry.
C    = 1 : an existing entry, which has not decayed or fragmented.
C        This is the main class of entries which represents the
C        "final state" given by the generator.
C    = 2 : an entry which has decayed or fragmented and therefore
C        is not appearing in the final state, but is retained for
C        event history information.
C    = 3 : a documentation line, defined separately from the event
C        history. (incoming reacting
C        particles, etc.)
C    = 4 - 10 : undefined, but reserved for future standards.
C    = 11 - 20 : at the disposal of each model builder for constructs
C        specific to his program, but equivalent to a null line in the
C        context of any other program. One example is the cone defining
C        vector of HERWIG, another cluster or event axes of the JETSET
C        analysis routines.
C    = 21 - : at the disposal of users, in particular for event tracking
C        in the detector.
C
C IDHKK(IHKK) : particle identity, according to the Particle Data Group
C    standard.
C
C JMOHKK(1,IHKK) : pointer to the position where the mother is stored.
C    The value is 0 for initial entries.
C
C JMOHKK(2,IHKK) : pointer to position of second mother. Normally only
C    one mother exist, in which case the value 0 is used. In cluster
C    fragmentation models, the two mothers would correspond to the q
C    and qbar which join to form a cluster. In string fragmentation,
C    the two mothers of a particle produced in the fragmentation would
C    be the two endpoints of the string (with the range in between
C    implied).
C
C JDAHKK(1,IHKK) : pointer to the position of the first daughter. If an
C    entry has not decayed, this is 0.
C
C JDAHKK(2,IHKK) : pointer to the position of the last daughter. If an
C    entry has not decayed, this is 0. It is assumed that the daughters
C    of a particle (or cluster or string) are stored sequentially, so
C    that the whole range JDAHKK(1,IHKK) - JDAHKK(2,IHKK) contains
C    daughters. Even in cases where only one daughter is defined (e.g.
C    K0 -> K0S) both values should be defined, to make for a uniform
C    approach in terms of loop constructions.
C
C PHKK(1,IHKK) : momentum in the x direction, in GeV/c.
C
C PHKK(2,IHKK) : momentum in the y direction, in GeV/c.
C
C PHKK(3,IHKK) : momentum in the z direction, in GeV/c.
C
C PHKK(4,IHKK) : energy, in GeV.
C
C PHKK(5,IHKK) : mass, in GeV/c**2. For spacelike partons, it is allowed
C    to use a negative mass, according to PHKK(5,IHKK) = -sqrt(-m**2).
C
C VHKK(1,IHKK) : production vertex x position, in mm.
C
C VHKK(2,IHKK) : production vertex y position, in mm.
C
C VHKK(3,IHKK) : production vertex z position, in mm.
C
C VHKK(4,IHKK) : production time, in mm/c (= 3.33*10**(-12) s).
C********************************************************************
*KEEP,DPAR.
C     /DPAR/   CONTAINS PARTICLE PROPERTIES
C        ANAME  = LITERAL NAME OF THE PARTICLE
C        AAM    = PARTICLE MASS IN GEV
C        GA     = DECAY WIDTH
C        TAU    = LIFE TIME OF INSTABLE PARTICLES
C        IICH   = ELECTRIC CHARGE OF THE PARTICLE
C        IIBAR  = BARYON NUMBER
C        K1,K1  = BIGIN AND END OF DECAY CHANNELS OF PARTICLE
C
      CHARACTER*8  ANAME
      COMMON /DPAR/ ANAME(210),AAM(210),GA(210),TAU(210),IICH(210),
     +IIBAR(210),K1(210),K2(210)
C------------------
*KEEP,FACTMO.
      COMMON /FACTMO/ IFACTO
*KEEP,TAUFO.
      COMMON /TAUFO/  TAUFOR,KTAUGE,ITAUVE,INCMOD
*KEEP,HADTHR.
      COMMON /HADTHR/ EHADTH,INTHAD
*KEEP,NUCC.
C     COMMON /NUCCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
C     COMMON /NUCC/   JT,JTZ,JP,JPZ,JJPROJ,JBPROJ,JJTARG,JBTARG
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
      COMMON /NUCCC/   JT,JTZ,JP,JPZ,JJPROJ,JBPROJ,JJTARG,JBTARG
*KEEP,ZENTRA.
      COMMON /ZENTRA/ ICENTR
*KEEP,CMHICO.
      COMMON /CMHICO/ CMHIS
*KEEP,RESONA.
      COMMON /RESONA/ IRESO
*KEEP,XSEADI.
      COMMON /XSEADI/ XSEACU,UNON,UNOM,UNOSEA, CVQ,CDQ,CSEA,SSMIMA,
     +SSMIMQ,VVMTHR
*KEEP,COULO.
      COMMON/COULO/ICOUL
*KEEP,EDENS.
      COMMON/EDENS/IEDEN
*KEEP,PROJK.
      COMMON /PROJK/ IPROJK
*KEEP,NUCIMP.
      COMMON /NUCIMP/ PRMOM(5,248),TAMOM(5,248), PRMFEP,PRMFEN,TAMFEP,
     +TAMFEN, PREFEP,PREFEN,TAEFEP,TAEFEN, PREPOT(210),TAEPOT(210),
     +PREBIN,TAEBIN,FERMOD,ETACOU
*KEND.
       COMMON /RECOM/IRECOM
C---------------------
*---minimum bias interactions
      ICENTR=0
*---threshold for the use of HADRIN in the primary hadron-nucleon collision
      EHADTH=5.
*---lab energy and momentum of the projectile: pion+
      EPN=200.
      IJPROJ=13
      JJPROJ=13
      PPN=SQRT((EPN-AAM(IJPROJ))*(EPN+AAM(IJPROJ)))
      IBPROJ=IIBAR(IJPROJ)
      IP=1
      IPZ=1
      JBPROJ=IIBAR(IJPROJ)
      JP=1
      JPZ=1
*---copper target
      IT=14
      ITZ=7 
      JT=14
      JTZ=7 
*---formation zone intranuclear cascade
      TAUFOR=105.
      KTAUGE=0 
      ITAUVE=1
*---inclusion of Coulomb potential for hA interactions
      ICOUL=1
      ICOULL=1
*---cascade within projectile switched off
      IPROJK=1
*---nucleus independent meson potential
      POTMES=0.002
      TAEPOT(13)=POTMES
      TAEPOT(14)=POTMES
      TAEPOT(15)=POTMES
      TAEPOT(16)=POTMES
      TAEPOT(23)=POTMES
      TAEPOT(24)=POTMES
      TAEPOT(25)=POTMES
*---definition of soft quark distributions
      XSEACU=0.05
      UNON=1.11D0
      UNOM=1.11D0
      UNOSEA=5.0D0
*---cutoff parameters for x-sampling
      CVQ=2.0D0
      CDQ=2.0D0
      CSEA=0.3D0
      SSMIMA=0.90D0
      SSMIMQ=SSMIMA**2
*---
      IRESO=0
      CMHIS=0.D0
      IEDEN=0
      IFACTO=0
C---Chain recombination
C     IRECOM=1
      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HADHAD(EPN,PPN,NHKKH1,IHTAWW,ITTA,IREJFO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C*** AT LOW ENERGIES EPN.LE.EHADTH HADRIN IS USED WITH ONE INTERACTION
C*** IN THE NUCLEUS INSTEAD OF THE DUAL PARTON MODEL (GLAUBER CASCADE)
C*** ONLY FOR HADRON-NUCLEUS COLLISIONS!
C
*KEEP,HKKEVT.
c     INCLUDE (HKKEVT)
      PARAMETER (NMXHKK= 89998)
c     PARAMETER (NMXHKK=25000)
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK), JMOHKK
     +(2,NMXHKK),JDAHKK(2,NMXHKK), PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK
     +(4,NMXHKK)
C
C                       WHKK(4,NMXHKK) GIVES POSITIONS AND TIMES IN
C                       PROJECTILE FRAME, THE CHAINS ARE CREATED ON
C                       THE POSITIONS OF THE PROJECTILE NUCLEONS
C                       IN THE PROJECTILE FRAME (TARGET NUCLEONS IN
C                       TARGET FRAME) BOTH POSITIONS ARE THREFORE NOT
C                       COMPLETELY CONSISTENT. THE TIMES IN THE
C                       PROJECTILE FRAME HOWEVER ARE OBTAINED BY
C                       LORENTZ TRANSFORMING FROM THE LAB SYSTEM.
C
C  Based on the proposed standard COMMON block (Sjostrand Memo 17.3,89)
C
C NMXHKK: maximum numbers of entries (partons/particles) that can be
C    stored in the commonblock.
C
C NHKK: the actual number of entries stored in current event. These are
C    found in the first NHKK positions of the respective arrays below.
C    Index IHKK, 1 <= IHKK <= NHKK, is below used to denote a given
C    entry.
C
C ISTHKK(IHKK): status code for entry IHKK, with following meanings:
C    = 0 : null entry.
C    = 1 : an existing entry, which has not decayed or fragmented.
C        This is the main class of entries which represents the
C        "final state" given by the generator.
C    = 2 : an entry which has decayed or fragmented and therefore
C        is not appearing in the final state, but is retained for
C        event history information.
C    = 3 : a documentation line, defined separately from the event
C        history. (incoming reacting
C        particles, etc.)
C    = 4 - 10 : undefined, but reserved for future standards.
C    = 11 - 20 : at the disposal of each model builder for constructs
C        specific to his program, but equivalent to a null line in the
C        context of any other program. One example is the cone defining
C        vector of HERWIG, another cluster or event axes of the JETSET
C        analysis routines.
C    = 21 - : at the disposal of users, in particular for event tracking
C        in the detector.
C
C IDHKK(IHKK) : particle identity, according to the Particle Data Group
C    standard.
C
C JMOHKK(1,IHKK) : pointer to the position where the mother is stored.
C    The value is 0 for initial entries.
C
C JMOHKK(2,IHKK) : pointer to position of second mother. Normally only
C    one mother exist, in which case the value 0 is used. In cluster
C    fragmentation models, the two mothers would correspond to the q
C    and qbar which join to form a cluster. In string fragmentation,
C    the two mothers of a particle produced in the fragmentation would
C    be the two endpoints of the string (with the range in between
C    implied).
C
C JDAHKK(1,IHKK) : pointer to the position of the first daughter. If an
C    entry has not decayed, this is 0.
C
C JDAHKK(2,IHKK) : pointer to the position of the last daughter. If an
C    entry has not decayed, this is 0. It is assumed that the daughters
C    of a particle (or cluster or string) are stored sequentially, so
C    that the whole range JDAHKK(1,IHKK) - JDAHKK(2,IHKK) contains
C    daughters. Even in cases where only one daughter is defined (e.g.
C    K0 -> K0S) both values should be defined, to make for a uniform
C    approach in terms of loop constructions.
C
C PHKK(1,IHKK) : momentum in the x direction, in GeV/c.
C
C PHKK(2,IHKK) : momentum in the y direction, in GeV/c.
C
C PHKK(3,IHKK) : momentum in the z direction, in GeV/c.
C
C PHKK(4,IHKK) : energy, in GeV.
C
C PHKK(5,IHKK) : mass, in GeV/c**2. For spacelike partons, it is allowed
C    to use a negative mass, according to PHKK(5,IHKK) = -sqrt(-m**2).
C
C VHKK(1,IHKK) : production vertex x position, in mm.
C
C VHKK(2,IHKK) : production vertex y position, in mm.
C
C VHKK(3,IHKK) : production vertex z position, in mm.
C
C VHKK(4,IHKK) : production time, in mm/c (= 3.33*10**(-12) s).
C********************************************************************
*KEEP,NUCC.
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
*KEEP,DPAR.
C     /DPAR/   CONTAINS PARTICLE PROPERTIES
C        ANAME  = LITERAL NAME OF THE PARTICLE
C        AAM    = PARTICLE MASS IN GEV
C        GA     = DECAY WIDTH
C        TAU    = LIFE TIME OF INSTABLE PARTICLES
C        IICH   = ELECTRIC CHARGE OF THE PARTICLE
C        IIBAR  = BARYON NUMBER
C        K1,K1  = BIGIN AND END OF DECAY CHANNELS OF PARTICLE
C
      CHARACTER*8  ANAME
      COMMON /DPAR/ ANAME(210),AAM(210),GA(210),TAU(210),IICH(210),
     +IIBAR(210),K1(210),K2(210)
C------------------
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEEP,DFINLS.
      PARAMETER (MAXFIN=10)
      COMMON /DFINLS/ ITRH(MAXFIN),CXRH(MAXFIN),CYRH(MAXFIN), CZRH
     +(MAXFIN),ELRH(MAXFIN),PLRH(MAXFIN),IRH
*KEEP,NUCIMP.
      COMMON /NUCIMP/ PRMOM(5,248),TAMOM(5,248), PRMFEP,PRMFEN,TAMFEP,
     +TAMFEN, PREFEP,PREFEN,TAEFEP,TAEFEN, PREPOT(210),TAEPOT(210),
     +PREBIN,TAEBIN,FERMOD,ETACOU
*KEND.
C------------------------------------------------------------------
C     PPN=SQRT((EPN-AAM(IJPROJ))*(EPN+AAM(IJPROJ)))
c     ipaupr=2
c     ipri=3
      IREJFO=0
      IF(IPRI.GE.2) WRITE(6,1001) IJPROJ,IJPROJ,PPN,EPN,CCCXP,CCCYP,
     +CCCZP,IHTAWW,ITTA,IELINE
 1001 FORMAT(' HADHAD 1:',
     +' IJPROJ,IJPROJ,PPN,EPN,CCCXP,CCCYP,CCCZP,IHTAWW,ITTA,IELINE'/ 2
     +I3,2E12.3,3F7.3,3I4)
      CCCXP=0.
      CCCYP=0.
      CCCZP=1.
      IELINE=0
      CALL SIHNIN(IJPROJ,ITTA,PPN,SIGHT)
      CALL SIHNEL(IJPROJ,ITTA,PPN,SIGHTE)
      SIGTOT=SIGHT + SIGHTE
      IF (SIGTOT*RNDM(BB).LE.SIGHTE)IELINE=1
      IF(IPRI.GE.2) WRITE(6,1000) IJPROJ,IJPROJ,PPN,EPN,CCCXP,CCCYP,
     +CCCZP,IHTAWW,ITTA,IELINE
 1000 FORMAT(' HADHAD 2 nach si...:',
     +' IJPROJ,IJPROJ,PPN,EPN,CCCXP,CCCYP,CCCZP,IHTAWW,ITTA,IELINE'/ 2
     +I3,2E12.3,3F7.3,3I4)
      IHADHA=0
   12 CONTINUE
      IHADHA=IHADHA+1
      IF(IPRI.GE.2) WRITE(6,1012) IJPROJ,IJPROJ,PPN,EPN,CCCXP,CCCYP,
     +CCCZP,IHTAWW,ITTA,IELINE
 1012 FORMAT(' HADHAD 12 loop:',
     +' IJPROJ,IJPROJ,PPN,EPN,CCCXP,CCCYP,CCCZP,IHTAWW,ITTA,IELINE'/ 2
     +I3,2E12.3,3F7.3,3I4)
        IF(IPRI.GE.3) THEN
          do 1212 ii=1,irh
          WRITE(6,'(I3,5(1PE12.4),I5/3X,5(1PE12.4))') II,ELRH(II),PLRH
     +    (II),CXRH(II),CYRH(II),CZRH(II),ITRH(II), (PHKK(JJJ,NHKK),JJJ
     +    =1,5)
 1212     continue
        ENDIF
*   repeated entry if Pauli blocking was active
      CALL FHAD(IJPROJ,IJPROJ,PPN,EPN,CCCXP,CCCYP,CCCZP, IHTAWW,ITTA,
     +IELINE,IREJFH)
            IF(IREJFH.EQ.1)THEN
              IREJFO=1
        IF(IPRI.GE.3) 
     + WRITE(6,'(A)')'  exit from hadhad with irejfo=1 '
              RETURN
            ENDIF
*
*   require Pauli blocking for final state nucleons
*
      IF (IHADHA.LT.3)THEN
        DO 11 II=1,IRH
          ITSEC=ITRH(II)
          IF(ITSEC.EQ.1.AND.ELRH(II).LE.TAEFEP+AAM(ITSEC))     GOTO 12
          IF(ITSEC.EQ.8.AND.ELRH(II).LE.TAEFEN+AAM(ITSEC))     GOTO 12
          IF(IIBAR(ITSEC).NE.1.AND.ELRH(II)-AAM(ITSEC)
     +                                   .LE.TAEPOT(ITSEC))  GOTO 12
   11   CONTINUE
      ENDIF
      NHKKH1=NHKK
C
      IF (IPRI.GE.2) WRITE (6,1010)IRH,NHKKH1,IHTAWW,ITTA
 1010 FORMAT (' HADHAD IRH,NHKKH1,IHTAWW,ITTA = ',4I5)
C
      IF(IPRI.GE.3) THEN
        WRITE(6,'(A/5X,A)')
     +  ' HADHAD - PARTICLE TRANSFER FROM /FINLSP/ INTO /HKKEVT/',
     +  ' II, ELRH, PLRH, CXRH, CYRH, CZRH / PHKK(1-5)'
      ENDIF
C
      ISTHKK(1)=11
      DO 10 II=1,IRH
C       IF( (ITSEC.EQ.1.AND.ELRH(II).GE.TAEFEP+AAM(ITSEC)) .OR.
C    +      (ITSEC.EQ.8.AND.ELRH(II).GE.TAEFEN+AAM(ITSEC)) )THEN
          NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)') ' HADHAD:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
          ITSEC=ITRH(II)
          IDHKK(NHKK)=MPDGHA(ITSEC)
          JMOHKK(1,NHKK)=1
          JMOHKK(2,NHKK)=IHTAWW
          JDAHKK(1,NHKK)=0
          JDAHKK(2,NHKK)=0
        PHKK(1,NHKK)=PLRH(II)*CXRH(II)
        PHKK(2,NHKK)=PLRH(II)*CYRH(II)
        PHKK(3,NHKK)=PLRH(II)*CZRH(II)
        PHKK(4,NHKK)=ELRH(II)
        IF(PHKK(4,NHKK)-AAM(ITSEC).LE.TAEPOT(ITSEC).
     +                                    AND.IIBAR(ITSEC).EQ.1)THEN
          ISTHKK(NHKK)=16
        ELSE
          ISTHKK(NHKK)=1
        ENDIF
        PHKK(5,NHKK)=AAM(ITRH(II))
C
        IF(IPRI.GE.3) THEN
          WRITE(6,'(I3,5(1PE12.4),I5/3X,5(1PE12.4),I5)') 
     +                    II,ELRH(II),PLRH
     +    (II),CXRH(II),CYRH(II),CZRH(II),ITRH(II), (PHKK(JJJ,NHKK),JJJ
     +    =1,5),IREJFO
        ENDIF
        VHKK(1,NHKK)=VHKK(1,IHTAWW)
        VHKK(2,NHKK)=VHKK(2,IHTAWW)
        VHKK(3,NHKK)=VHKK(3,IHTAWW)
        VHKK(4,NHKK)=VHKK(4,1)
C       ENDIF
   10 CONTINUE
      JDAHKK(1,1)=NHKKH1+1
      JDAHKK(2,1)=NHKK
      JDAHKK(1,IHTAWW)=NHKKH1+1
      JDAHKK(2,IHTAWW)=NHKK
c     ipaupr=0
c     ipri=0
        IF(IPRI.GE.3) 
     + WRITE(6,'(A)')'  exit from hadhad with irejfo=0 '
      RETURN
      END
      SUBROUTINE CHEBCH(IREJ,NHKKH1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
*KEEP,HKKEVT.
c     INCLUDE (HKKEVT)
      PARAMETER (NMXHKK= 89998)
c     PARAMETER (NMXHKK=25000)
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK), JMOHKK
     +(2,NMXHKK),JDAHKK(2,NMXHKK), PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK
     +(4,NMXHKK)
      COMMON /EXTEVT/ IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),
     &                IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10)
C
C                       WHKK(4,NMXHKK) GIVES POSITIONS AND TIMES IN
C                       PROJECTILE FRAME, THE CHAINS ARE CREATED ON
C                       THE POSITIONS OF THE PROJECTILE NUCLEONS
C                       IN THE PROJECTILE FRAME (TARGET NUCLEONS IN
C                       TARGET FRAME) BOTH POSITIONS ARE THREFORE NOT
C                       COMPLETELY CONSISTENT. THE TIMES IN THE
C                       PROJECTILE FRAME HOWEVER ARE OBTAINED BY
C                       LORENTZ TRANSFORMING FROM THE LAB SYSTEM.
C
C  Based on the proposed standard COMMON block (Sjostrand Memo 17.3,89)
C
C NMXHKK: maximum numbers of entries (partons/particles) that can be
C    stored in the commonblock.
C
C NHKK: the actual number of entries stored in current event. These are
C    found in the first NHKK positions of the respective arrays below.
C    Index IHKK, 1 <= IHKK <= NHKK, is below used to denote a given
C    entry.
C
C ISTHKK(IHKK): status code for entry IHKK, with following meanings:
C    = 0 : null entry.
C    = 1 : an existing entry, which has not decayed or fragmented.
C        This is the main class of entries which represents the
C        "final state" given by the generator.
C    = 2 : an entry which has decayed or fragmented and therefore
C        is not appearing in the final state, but is retained for
C        event history information.
C    = 3 : a documentation line, defined separately from the event
C        history. (incoming reacting
C        particles, etc.)
C    = 4 - 10 : undefined, but reserved for future standards.
C    = 11 - 20 : at the disposal of each model builder for constructs
C        specific to his program, but equivalent to a null line in the
C        context of any other program. One example is the cone defining
C        vector of HERWIG, another cluster or event axes of the JETSET
C        analysis routines.
C    = 21 - : at the disposal of users, in particular for event tracking
C        in the detector.
C
C IDHKK(IHKK) : particle identity, according to the Particle Data Group
C    standard.
C
C JMOHKK(1,IHKK) : pointer to the position where the mother is stored.
C    The value is 0 for initial entries.
C
C JMOHKK(2,IHKK) : pointer to position of second mother. Normally only
C    one mother exist, in which case the value 0 is used. In cluster
C    fragmentation models, the two mothers would correspond to the q
C    and qbar which join to form a cluster. In string fragmentation,
C    the two mothers of a particle produced in the fragmentation would
C    be the two endpoints of the string (with the range in between
C    implied).
C
C JDAHKK(1,IHKK) : pointer to the position of the first daughter. If an
C    entry has not decayed, this is 0.
C
C JDAHKK(2,IHKK) : pointer to the position of the last daughter. If an
C    entry has not decayed, this is 0. It is assumed that the daughters
C    of a particle (or cluster or string) are stored sequentially, so
C    that the whole range JDAHKK(1,IHKK) - JDAHKK(2,IHKK) contains
C    daughters. Even in cases where only one daughter is defined (e.g.
C    K0 -> K0S) both values should be defined, to make for a uniform
C    approach in terms of loop constructions.
C
C PHKK(1,IHKK) : momentum in the x direction, in GeV/c.
C
C PHKK(2,IHKK) : momentum in the y direction, in GeV/c.
C
C PHKK(3,IHKK) : momentum in the z direction, in GeV/c.
C
C PHKK(4,IHKK) : energy, in GeV.
C
C PHKK(5,IHKK) : mass, in GeV/c**2. For spacelike partons, it is allowed
C    to use a negative mass, according to PHKK(5,IHKK) = -sqrt(-m**2).
C
C VHKK(1,IHKK) : production vertex x position, in mm.
C
C VHKK(2,IHKK) : production vertex y position, in mm.
C
C VHKK(3,IHKK) : production vertex z position, in mm.
C
C VHKK(4,IHKK) : production time, in mm/c (= 3.33*10**(-12) s).
C********************************************************************
*KEEP,DPAR.
C     /DPAR/   CONTAINS PARTICLE PROPERTIES
C        ANAME  = LITERAL NAME OF THE PARTICLE
C        AAM    = PARTICLE MASS IN GEV
C        GA     = DECAY WIDTH
C        TAU    = LIFE TIME OF INSTABLE PARTICLES
C        IICH   = ELECTRIC CHARGE OF THE PARTICLE
C        IIBAR  = BARYON NUMBER
C        K1,K1  = BIGIN AND END OF DECAY CHANNELS OF PARTICLE
C
      CHARACTER*8  ANAME
      COMMON /DPAR/ ANAME(210),AAM(210),GA(210),TAU(210),IICH(210),
     +IIBAR(210),K1(210),K2(210)
C------------------
*KEEP,NUCC.
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEND.
      COMMON /CHABAI/CHARGI,BARNUI
      COMMON /EVAPPP/IEVAP
C-----------------------------------------------------------------------
      DATA IEVL /0/
C----------------------------------------------------------------------
      ZERO=0
      ONEONE=1
      TWOTWO=2
      NHAD=0
      NIP=IP
      AIP=IP
      IEVL=IEVL+1
      CHAEVE=0.
      BAEVE=0.
      IF(IEVAP.EQ.0)THEN
      DO 1171 I=1,NHKKH1
        IF (ISTHKK(I).EQ.13)THEN
          BAEVE=BAEVE+1
	  IF(IDHKK(I).EQ.2212)CHAEVE=CHAEVE+1.D0
        ENDIF
        IF (ISTHKK(I).EQ.14)THEN
          BAEVE=BAEVE+1
	  IF(IDHKK(I).EQ.2212)CHAEVE=CHAEVE+1.D0
        ENDIF
 1171 CONTINUE
      DO 521 I=NHKKH1,NHKK
        IF (ISTHKK(I).EQ.1.OR.ISTHKK(I).EQ.15.OR.ISTHKK(I).EQ.16)THEN
          NHAD=NHAD+1
          NRHKK=MCIHAD(IDHKK(I))
          IF (NRHKK.LE.0.OR.NRHKK.GT.410)THEN
            WRITE(6,1389)NRHKK,I,IDHKK(I),NHKKH1,NHKK
 1389       FORMAT (' distr: NRHKK ERROR ',5I10)
            NRHKK=1
          ENDIF
          ICHHKK=IICH(NRHKK)
          IBHKK=IIBAR(NRHKK)
	  CHAEVE=CHAEVE+ICHHKK
          BAEVE=BAEVE+IBHKK
        ENDIF
  521 CONTINUE
      ELSEIF(IEVAP.EQ.1)THEN
      DO 1521 I=1,NHKK
        IF (ISTHKK(I).EQ.1)THEN
          NRHKK=MCIHAD(IDHKK(I))
          ICHHKK=IICH(NRHKK)
          IBHKK=IIBAR(NRHKK)
	  CHAEVE=CHAEVE+ICHHKK
          BAEVE=BAEVE+IBHKK
	ENDIF
 1521 CONTINUE
C     WRITE(6,'(A,2F12.1)')' after isthkk=1 ',CHAEVE,BAEVE
      DO 2521 I=1,NHKK
        IF (ISTHKK(I).EQ.-1)THEN
	  IF(IDHKK(I).EQ.2112)THEN
            BAEVE=BAEVE+1
C	    WRITE(6,'(A,2F12.1)')' evap isthkk=-1',CHAEVE,BAEVE
	  ENDIF
	  IF(IDHKK(I).EQ.2212)THEN
	    CHAEVE=CHAEVE+1
            BAEVE=BAEVE+1
C	    WRITE(6,'(A,2F12.1)')' evap isthkk=-1',CHAEVE,BAEVE
	  ENDIF
	ENDIF
	IF((IDHKK(I).EQ.80000).AND.(ISTHKK(I).NE.1000))THEN
	  CHAEVE=CHAEVE+IDXRES(I)
	  BAEVE=BAEVE+IDRES(I)
C     WRITE(6,'(A,2F12.1,2I5)')' h.f.',CHAEVE,BAEVE,IDXRES(I),IDRES(I)
	ENDIF
 2521 CONTINUE
      ENDIF
      IF(IEVL.LE.10)WRITE(6,'(2A,4F10.2)')' Event charge and B-number',
     * '=',CHAEVE,BAEVE,CHARGI,BARNUI
      IF(CHAEVE-CHARGI.NE.0.D0.OR.BAEVE-BARNUI.NE.0.D0)THEN
C     DO 775 JJJ=1,200
      IF(IEVL.LE.1000)WRITE(6,'(2A,4F10.2)')'Event charge and B-numb',
     *'(violated)  =',CHAEVE,BAEVE,CHARGI,BARNUI
C 775 CONTINUE
      IREJ=1
      ENDIF
      RETURN
      END 
C****************************************************************
C
      SUBROUTINE PARPT(IFL,PT1,PT2,IPT,NEVT)
C                          Plot parton pt distribution
C 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      DIMENSION PT(50,10),YPT(50,10)
      GO TO (1,2,3),IFL
 1    CONTINUE
      DPT=0.1
      DO 10 I=1,10
        DO 10 J=1,50
          PT(J,I)=J*DPT-DPT/2.
          YPT(J,I)=1.D-50
 10   CONTINUE
      RETURN
 2    CONTINUE
      IPT1=PT1/DPT+1.
      IPT2=PT2/DPT+1.
      IF(IPT1.GT.50)IPT1=50
      IF(IPT2.GT.50)IPT2=50
      YPT(IPT1,IPT)=YPT(IPT1,IPT)+1.
      YPT(IPT2,IPT)=YPT(IPT2,IPT)+1.
      YPT(IPT1,10)=YPT(IPT1,10)+1.
      YPT(IPT2,10)=YPT(IPT2,10)+1.
      RETURN
 3    CONTINUE
      DO 30 I=1,10
        DO 30 J=1,50
          YPT(J,I)=YPT(J,I)/NEVT
          YPT(J,I)=LOG10(YPT(J,I)+1.D-18)
 30   CONTINUE
C     WRITE(6,*)' Parton pt distribution,vv=1,vsr=+2,sv=3,ss=4,zz=5,
C    *  hh=6,10=all'
C     CALL PLOT(PT,YPT,500,10,50,0.D0,DPT,-5.D0,0.1D0)
      RETURN
      END
*
*
*
*===hkkfil=============================================================*
*
      SUBROUTINE HKKFIL(IST,ID,M1,M2,PX,PY,PZ,E,NHKKAU,KORMO,ICALL)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)
      PARAMETER (TINY10=1.0D-10,TINY4=1.0D-3)

*KEEP,HKKEVT.
c     INCLUDE (HKKEVT)
      PARAMETER (NMXHKK= 89998)
c     PARAMETER (NMXHKK=25000)
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK), JMOHKK
     +(2,NMXHKK),JDAHKK(2,NMXHKK), PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK
     +(4,NMXHKK)
C
C                       WHKK(4,NMXHKK) GIVES POSITIONS AND TIMES IN
C                       PROJECTILE FRAME, THE CHAINS ARE CREATED ON
C                       THE POSITIONS OF THE PROJECTILE NUCLEONS
C                       IN THE PROJECTILE FRAME (TARGET NUCLEONS IN
C                       TARGET FRAME) BOTH POSITIONS ARE THREFORE NOT
C                       COMPLETELY CONSISTENT. THE TIMES IN THE
C                       PROJECTILE FRAME HOWEVER ARE OBTAINED BY
C                       LORENTZ TRANSFORMING FROM THE LAB SYSTEM.
C
C  Based on the proposed standard COMMON block (Sjostrand Memo 17.3,89)
C
C NMXHKK: maximum numbers of entries (partons/particles) that can be
C    stored in the commonblock.
C
C NHKK: the actual number of entries stored in current event. These are
C    found in the first NHKK positions of the respective arrays below.
C    Index IHKK, 1 <= IHKK <= NHKK, is below used to denote a given
C    entry.
C
C ISTHKK(IHKK): status code for entry IHKK, with following meanings:
C    = 0 : null entry.
C    = 1 : an existing entry, which has not decayed or fragmented.
C        This is the main class of entries which represents the
C        "final state" given by the generator.
C    = 2 : an entry which has decayed or fragmented and therefore
C        is not appearing in the final state, but is retained for
C        event history information.
C    = 3 : a documentation line, defined separately from the event
C        history. (incoming reacting
C        particles, etc.)
C    = 4 - 10 : undefined, but reserved for future standards.
C    = 11 - 20 : at the disposal of each model builder for constructs
C        specific to his program, but equivalent to a null line in the
C        context of any other program. One example is the cone defining
C        vector of HERWIG, another cluster or event axes of the JETSET
C        analysis routines.
C    = 21 - : at the disposal of users, in particular for event tracking
C        in the detector.
C
C IDHKK(IHKK) : particle identity, according to the Particle Data Group
C    standard.
C
C JMOHKK(1,IHKK) : pointer to the position where the mother is stored.
C    The value is 0 for initial entries.
C
C JMOHKK(2,IHKK) : pointer to position of second mother. Normally only
C    one mother exist, in which case the value 0 is used. In cluster
C    fragmentation models, the two mothers would correspond to the q
C    and qbar which join to form a cluster. In string fragmentation,
C    the two mothers of a particle produced in the fragmentation would
C    be the two endpoints of the string (with the range in between
C    implied).
C
C JDAHKK(1,IHKK) : pointer to the position of the first daughter. If an
C    entry has not decayed, this is 0.
C
C JDAHKK(2,IHKK) : pointer to the position of the last daughter. If an
C    entry has not decayed, this is 0. It is assumed that the daughters
C    of a particle (or cluster or string) are stored sequentially, so
C    that the whole range JDAHKK(1,IHKK) - JDAHKK(2,IHKK) contains
C    daughters. Even in cases where only one daughter is defined (e.g.
C    K0 -> K0S) both values should be defined, to make for a uniform
C    approach in terms of loop constructions.
C
C PHKK(1,IHKK) : momentum in the x direction, in GeV/c.
C
C PHKK(2,IHKK) : momentum in the y direction, in GeV/c.
C
C PHKK(3,IHKK) : momentum in the z direction, in GeV/c.
C
C PHKK(4,IHKK) : energy, in GeV.
C
C PHKK(5,IHKK) : mass, in GeV/c**2. For spacelike partons, it is allowed
C    to use a negative mass, according to PHKK(5,IHKK) = -sqrt(-m**2).
C
C VHKK(1,IHKK) : production vertex x position, in mm.
C
C VHKK(2,IHKK) : production vertex y position, in mm.
C
C VHKK(3,IHKK) : production vertex z position, in mm.
C
C VHKK(4,IHKK) : production time, in mm/c (= 3.33*10**(-12) s).
C********************************************************************
      COMMON /NNCMS/  GACMS,BGCMS,UMO,PCM,EPROJ,PPROJ
      COMMON /TRAFOP/GALAB,BGLAB,BLAB
      COMMON /PROJK/ IPROJK
      COMMON /NDON/NDONE

C     IF (MODE.GT.100) THEN
C        WRITE(LOUT,'(1X,A,I5,A,I5)')
C    &        'HKKFIL: reset NHKK = ',NHKK,' to NHKK =',NHKK-MODE+100
C        NHKK = NHKK-MODE+100
C        RETURN
C     ENDIF
      MO1  = M1
      MO2  = M2
      NHKK = NHKK+1
      IF (NHKK.GT.NMXHKK) THEN
         WRITE(LOUT,1000) NHKK
 1000    FORMAT(1X,'HKKFIL: NHKK exeeds NMXHKK = ',I7,
     &             '! program execution stopped..')
         STOP
      ENDIF
      IF (M1.LT.0) MO1 = NHKK+M1
      IF (M2.LT.0) MO2 = NHKK+M2
      ISTHKK(NHKK)   = IST
      IDHKK(NHKK)    = ID
      IF(KORMO.EQ.999)THEN
        JMOHKK(1,NHKK) = MO1
        JMOHKK(2,NHKK) = MO2
      ELSE
        JMOHKK(1,NHKK)=NHKKAU+KORMO-1
        JMOHKK(2,NHKK)=0
      ENDIF
      IF(NHKK.LE.JMOHKK(1,NHKK))THEN
C     SUBROUTINE HKKFIL(IST,ID,M1,M2,PX,PY,PZ,E,NHKKAU,KORMO,ICALL)
        WRITE(6,*)' HKKFIL(IST,ID,M1,M2,PX,PY,PZ,E,NHKKAU,KORMO)',
     *  NHKK,IST,ID,M1,M2,PX,PY,PZ,E,NHKKAU,KORMO,ICALL,JMOHKK(1,NHKK)	
      ENDIF
      JDAHKK(1,NHKK) = 0
      JDAHKK(2,NHKK) = 0
      IF (MO1.GT.0) THEN
         IF (JDAHKK(1,MO1).NE.0) THEN
            JDAHKK(2,MO1) = NHKK
         ELSE
            JDAHKK(1,MO1) = NHKK
         ENDIF
	 JDAHKK(1,MO1)=NHKKAU
      ENDIF
      IF (MO2.GT.0) THEN
         IF (JDAHKK(1,MO2).NE.0) THEN
            JDAHKK(2,MO2) = NHKK
         ELSE
            JDAHKK(1,MO2) = NHKK
         ENDIF
         JDAHKK(1,MO2) = NHKKAU
      ENDIF
      PHKK(1,NHKK) = PX
      PHKK(2,NHKK) = PY
      PHKK(3,NHKK) = PZ
      PHKK(4,NHKK) = E
      PHKK(5,NHKK) = PHKK(4,NHKK)**2-PHKK(1,NHKK)**2-
     &               PHKK(2,NHKK)**2-PHKK(3,NHKK)**2
      IF ((PHKK(5,NHKK).LT.0.0D0).AND.(ABS(PHKK(5,NHKK)).GT.TINY4))
     &   WRITE(LOUT,'(1X,A,G10.3)')
     &     'HKKFIL: negative mass**2 ',PHKK(5,NHKK)
      PHKK(5,NHKK) = SQRT(ABS(PHKK(5,NHKK)))
      IF (IST.EQ.88888.OR.IST.EQ.88887.OR.IST.EQ.88889) THEN
* special treatment for chains:
*    position of chain in Lab      = pos. of target nucleon
*    time of chain-creation in Lab = time of passage of projectile
*                                    nucleus at pos. of taget nucleus
         DO 1 I=1,3
            VHKK(I,NHKK) = VHKK(I,MO2)
    1    CONTINUE
         VHKK(4,NHKK) = VHKK(3,MO2)/BLAB-VHKK(3,MO1)/BGLAB
      ELSE
         IF(MO1.GE.1)THEN
         DO 2 I=1,4
            VHKK(I,NHKK) = VHKK(I,MO1)
            IF (IPROJK.EQ.1) THEN
              WHKK(I,NHKK) = WHKK(I,MO1)
	    ENDIF
    2    CONTINUE
         ENDIF
      ENDIF

      RETURN
      END

C*********************************************************************
C*********************************************************************

      SUBROUTINE JSPARR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      INTEGER PYCOMP

C...Purpose: to give program heading, or list an event, or particle
C...data, or current parameter values.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(4000,2),BRAT(4000),KFDP(4000,5)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/,/PYDAT3/
      CHARACTER CHAP*16,CHAN*16,CHAD(5)*16
      DIMENSION	  KBAMDP(5)

C...List parton/particle data table. Check whether to be listed.
        WRITE(MSTU(11),6800)
        MSTJ24=MSTJ(24)
        MSTJ(24)=0
        KFMAX=20883
        IF(MSTU(2).NE.0) KFMAX=MSTU(2)
C                               KF = PDG Particle number
        DO 220 KF=100,KFMAX
C                               KC = Lund particle number
        KC=PYCOMP(KF)
        IF(KC.EQ.0) GOTO 220
        IF(MSTU(14).EQ.0.AND.KF.GT.100.AND.KC.LE.100) GOTO 220
        IF(MSTU(14).GT.0.AND.KF.GT.100.AND.MAX(MOD(KF/1000,10),
     &  MOD(KF/100,10)).GT.MSTU(14)) GOTO 220
C                          BAMJET particle number
	KBAM=MCIHAD(KF)
C                          BAMJET  antiparticle number
	KABAM=MCIHAD(-KF)

C...Find  Lund particle name and mass. Print information.
        CALL PYNAME(KF,CHAP)
        IF(KF.LE.100.AND.CHAP.EQ.' '.AND.MDCY(KC,2).EQ.0) GOTO 220
C                          Lund Antiparticle Name
        CALL PYNAME(-KF,CHAN)
        PM=PYMASS(KF)
	IDC1=MDCY(KC,2)
	IDC2=MDCY(KC,2)+MDCY(KC,3)-1
        WRITE(MSTU(11),6900)KBAM,
     &	KF,KC,IDC1,IDC2,CHAP,CHAN,KCHG(KC,1),KCHG(KC,2),
     &  KCHG(KC,3),PM,PMAS(KC,2),PMAS(KC,3),PMAS(KC,4),MDCY(KC,1)
        WRITE(26,6900)KBAM,
     &	KF,KC,IDC1,IDC2,CHAP,CHAN,KCHG(KC,1),KCHG(KC,2),
     &  KCHG(KC,3),PM,PMAS(KC,2),PMAS(KC,3),PMAS(KC,4),MDCY(KC,1)

C...Particle decay: channel number, branching ration, matrix element,
C...decay products.
        IF(KF.GT.100.AND.KC.LE.100) GOTO 220
        DO 210 IDC=MDCY(KC,2),MDCY(KC,2)+MDCY(KC,3)-1
        DO 200 J=1,5
C                         Lund names of decay products
          CALL PYNAME(KFDP(IDC,J),CHAD(J))
C                         Bamjet numbers of decay products
	  KBAMDP(J)=MCIHAD(KFDP(IDC,J))
	  IF(KBAMDP(J).EQ.26)KBAMDP(J)=0
  200   CONTINUE
        WRITE(26,7001) IDC,MDME(IDC,1),MDME(IDC,2),BRAT(IDC),
     &  (KBAMDP(J),J=1,5)
  210   WRITE(MSTU(11),7000) IDC,MDME(IDC,1),MDME(IDC,2),BRAT(IDC),
     &  (CHAD(J),J=1,5)
C              The same for the antiparticle, if it exists
	IF(KABAM.NE.410)THEN
        WRITE(MSTU(11),6900)KABAM,
     &	-KF,-KC,IDC1,IDC2,CHAN,CHAP,-KCHG(KC,1),KCHG(KC,2),
     &  KCHG(KC,3),PM,PMAS(KC,2),PMAS(KC,3),PMAS(KC,4),MDCY(KC,1)
        WRITE(26,6900)KABAM,
     &	-KF,-KC,IDC1,IDC2,CHAN,CHAP,-KCHG(KC,1),KCHG(KC,2),
     &  KCHG(KC,3),PM,PMAS(KC,2),PMAS(KC,3),PMAS(KC,4),MDCY(KC,1)
        DO 211 IDC=MDCY(KC,2),MDCY(KC,2)+MDCY(KC,3)-1
        DO 201 J=1,5
C                               KC = Lund particle number
        KCDP=PYCOMP(KFDP(IDC,J))
	IF(KCDP.LE.0.OR.KCDP.GT.500)THEN
C      	WRITE(MSTU(11),'(A,I10)')' KCDP= ',KCDP
	KCDP=1
	ENDIF
C                         Bamjet numbers of decay products
	  KFDPM=-KFDP(IDC,J)
	  IF(KCHG(KCDP,3).EQ.0)KFDPM=KFDP(IDC,J)
	  KBAMDP(J)=MCIHAD(KFDPM)
	  IF(KBAMDP(J).EQ.26)KBAMDP(J)=0
C                         Lund names of decay products
          CALL PYNAME(KFDPM,CHAD(J))
  201   CONTINUE
        WRITE(26,7001) IDC,MDME(IDC,1),MDME(IDC,2),BRAT(IDC),
     &  (KBAMDP(J),J=1,5)
  211   WRITE(MSTU(11),7000) IDC,MDME(IDC,1),MDME(IDC,2),BRAT(IDC),
     &  (CHAD(J),J=1,5)
	ENDIF
  220   CONTINUE
        MSTJ(24)=MSTJ24


C...Format statements for output on unit MSTU(11) (by default 6).
 6800 FORMAT(///30X,'Particle/parton data table'//1X,'BAM',
     &1X,'ABAM',1X,'KF',1X,'KC',1X,'DCF',1X,'DCL',1X,
     &'particle',8X,'antiparticle',6X,'chg  col  anti',8X,'mass',7X,
     &'width',7X,'w-cut',5X,'lifetime',1X,'decay'/11X,'IDC',1X,'on/off',
     &1X,'ME',3X,'Br.rat.',4X,'decay products')
 6900 FORMAT(/1X,I4,I6,I4,2I5,A16,A16,3I3,1X,F12.5,2(1X,F11.5),
     &1X,F12.5,1X,I2)
 7000 FORMAT(10X,I4,2X,I3,2X,I3,2X,F8.5,4X,5A16)
 7001 FORMAT(10X,I4,2X,I3,2X,I3,2X,F8.5,4X,5I5)

      RETURN
      END

C*********************************************************************

