*-- Author :
C-------------------------------------------------------------------
C
C                 FILE DNDISTRM FORTRAN
C
C------------------------------------------------------------------
C     DEBUG SUBCHK
C     END DEBUG
      SUBROUTINE DISTR(IOP,NHKKH1,PO,IGENER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C***
C   REDUCED VERSION FOR SCORING RESULTS FROM PRIMARY AA INTERACTION
C   FROM /HKKEVT/
C   kinematics in lab (target rest system)
C      - MULTIPLICITIES
C      - number of inc generations inside the nucleus
C      - RAPIDITY DISTRIBUTION         (without cuts)
C      - pseudorapidity distribution   (without cuts)
C***
C  conventions for rapidity distributions:
C       1  Proton                       11  NEGATIVES
C       2  Neutron                      12  NBAR=9
C       3  PI+=13                       13  LAMBDA=17
C       4  PI-=14                       14  SIGMA==20,21,22
C       5  K+ =15                       15  LAMBDABAR=18
C       6  K- =16                       16  SIGMABAR=99,100,101
C       7  neutral kaons=12,19,24,25    17  THETA=97,98
C       8  pbar=2                       18  THETABAR=102,103
C       9  charged hadrons
C      10  total hadrons
C-----------------------------------------------------------------------
      DIMENSION P4P4P4(4)
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
*KEEP,NSHMAK.
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C      COMMON /NSHMAK/ NNSHMA,NPSHMA,NTSHMA,NSHMAC
      COMMON /NSHMAK/ NNSHMA,NPSHMA,NTSHMA,NSHMAC,NSHMA2
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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
C                                   modified DPMJET
       COMMON /BUFUEH/ ANNVV,ANNSS,ANNSV,ANNVS,ANNCC,
     *                 ANNDV,ANNVD,ANNDS,ANNSD,
     *                 ANNHH,ANNZZ,
     *                 PTVV,PTSS,PTSV,PTVS,PTCC,PTDV,PTVD,PTDS,PTSD,
     *                 PTHH,PTZZ,
     *                 EEVV,EESS,EESV,EEVS,EECC,EEDV,EEVD,EEDS,EESD,
     *                 EEHH,EEZZ
     *                ,ANNDI,PTDI,EEDI
     *                ,ANNZD,ANNDZ,PTZD,PTDZ,EEZD,EEDZ
C---------------------
C---------------------
      DIMENSION  XYL(50,20),YYL(50,20),YYLPS(50,20),INDX(28)
      DIMENSION DISGEN(50),XGEN(50)
      PARAMETER (NUMTYP=160)
      DIMENSION AVMULT(NUMTYP),AVE(NUMTYP),AKE(NUMTYP),AVEPT(NUMTYP)
      DATA DISGEN /50*0.0/
      DATA INDX/ 1, 8,-1,-1,-1,-1,-1, 2,12,-1,-1,7,3,4,5,6,13,15,7,14,
     *          14,14,19, 7, 7,-1,-1,-1/
C----------------------------------------------------------------------
      ONEONE=1
      ZERO=0
      TWOTWO=2
      HUNDTH=1.D-2
C     CALL DISPT(IOP,NHKKH1,PO)
C     CALL DISEVA(IOP,NHKKH1,PO,IGENER)
      GO TO (1,2,3),IOP
*  initialization call
1     CONTINUE
      ANCHSQ=0.
      DELRAP=1.
      IF(IP.EQ.1)DELRAP=0.1
      NHKKH2=NHKKH1
      IF (NHKKH1.EQ.0)NHKKH2=1
      EEO=SQRT(PO**2+AAM(NHKKH2)**2)
      WRITE(6, 1001)EEO,PO,NHKKH2,AAM(NHKKH2)
 1001 FORMAT (' EEO',F10.2,F10.2,I10,F10.2)
      DY=0.4
C
      DO 11 I=1,20
        DO 13 J=1,50
          XYL(J,I)=-2.0 + (J-1)*DY
          YYL(J,I)=1.E-18
          YYLPS(J,I)=1.E-18
13      CONTINUE
11    CONTINUE
      DO 61 I=1,NUMTYP
        AVE(I)       =1.E-18
        AVEPT(I)       =1.E-18
        AVMULT(I)    =1.E-18
   61 CONTINUE
      RETURN
C-------------------------------------------------------------------
  2   CONTINUE
C                               SEWEW secondary interactions
      CALL SEWEW(1,NHKKH1)
*  scoring of results from individual events
C                            NUMBER OF GENERATIONS INSIDE THE NUCLEUS
      IF(IGENER.LT.1.OR.IGENER.GT.50) IGENER=50
      DISGEN(IGENER)=DISGEN(IGENER) + 1.0
C
C     WRITE (18,1711)
C1711 FORMAT (' event px,py,pz,e,id ')
      NHAD=0
      AVMULC=0.
C                          HBOOK HISTOGRAMS
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     WRITE(6,'(A)')' before plomb'
      IF(IHBOOK.EQ.1)CALL PLOMB (2,P4P4P4,CCCHRG,XFXFXF,1,IJPROJ)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                Here we include the decayed resonances into the
C                        multiplicity Table
      IPRIEV=0
      DO 521 I=NHKKH1,NHKK
        IF (ISTHKK(I).EQ.2)THEN
C         WRITE(18,1712)PHKK(1,I),PHKK(2,I),PHKK(3,I),Phkk(4,I),IDHKK(I)
C         NHAD=NHAD+1
          NRHKK=MCIHAD(IDHKK(I))
          IF (NRHKK.LE.0.OR.NRHKK.GT.210)THEN
            WRITE(6,1389)NRHKK,I,IDHKK(I),NHKKH1,NHKK
            NRHKK=1
          ENDIF
          NRE=NRHKK
          NREM=NRE
          ICHHKK=IICH(NRHKK)
	  IF(NRE.GT.30)THEN
	  IF(NRE.GT.160)GO TO 521
            AVE(NREM)=AVE(NREM) + PHKK(4,I)
            AVEPT(NREM)=AVEPT(NREM) + PT
            AVMULT(NREM)=AVMULT(NREM) + 1.
	  ENDIF
        ENDIF
  521 CONTINUE
      DO 21 I=NHKKH1,NHKK
        IF (ISTHKK(I).EQ.1)THEN
C         WRITE(18,1712)PHKK(1,I),PHKK(2,I),PHKK(3,I),Phkk(4,I),IDHKK(I)
 1712     FORMAT (4E14.5,I8)
          NHAD=NHAD+1
          NRHKK=MCIHAD(IDHKK(I))
          IF (NRHKK.LE.0.OR.NRHKK.GT.210)THEN
            WRITE(6,1389)NRHKK,I,IDHKK(I),NHKKH1,NHKK
 1389       FORMAT (' DISTR: NRHKK ERROR ',5I10)
            NRHKK=1
          ENDIF
          NRE=NRHKK
          NREM=NRE
          ICHHKK=IICH(NRHKK)
*  rapidity
          PTT=PHKK(1,I)**2+PHKK(2,I)**2+0.000001
          PT=SQRT(PTT)
C         IF(PT.LT.0.5)GO TO 21
          AMT=SQRT(PTT+PHKK(5,I)**2)
C***      AMT=SQRT(PTT+AMPI**2)
C         ETOT=SQRT(AMT**2 + PHKK(3,I)**2)
C         YL=LOG((ETOT + PHKK(3,I))/AMT+1.E-18)
          YL=LOG((ABS(PHKK(3,I) + PHKK(4,I)))/AMT+1.E-18)
C         IF(YL.LT.0.75.OR.YL.GT.3.) GO TO 21  
          IF (NRE.GT.25) NRE=28
          IF (NRE.LT. 1) NRE=28
          IF(NREM.GT.NUMTYP) NREM=28
          IF(NREM.LT.1) NREM=28
          NI=INDX(NRE)
          IF (NRHKK.LE.101.AND.NRHKK.GE.99) NI=16
          IF (NRHKK.EQ.97.OR.NRHKK.EQ.98)   NI=17
          IF (NRHKK.EQ.102.OR.NRHKK.EQ.103) NI=18
          AVE(NREM)=AVE(NREM) + PHKK(4,I)
          AVEPT(NREM)=AVEPT(NREM) + PT
          AVMULT(NREM)=AVMULT(NREM) + 1.
C    TOTAL=30
          AVE(30)=AVE(30) + PHKK(4,I)
          AVEPT(30)=AVEPT(30) + PT
C    CHARGED=27
          IF (ICHHKK.NE.0)  THEN
            AVE(27)=AVE(27) + PHKK(4,I)
            AVEPT(27)=AVEPT(27) + PT
            AVMULC=AVMULC + 1.
            AVMULT(27)=AVMULT(27) + 1.
          ENDIF
          IYL=(YL+2.0)/DY + 1
          IF (IYL.LT.1)  IYL=1
          IF (IYL.GT.50) IYL=50
          IF (ICHHKK.NE.0) THEN
            YYL(IYL,9)=YYL(IYL,9)+1.
          ENDIF
          IF (ICHHKK.LT.0)    YYL(IYL,11)=YYL(IYL,11)+1.
          IF(NI.GT.0) YYL(IYL,NI)=YYL(IYL,NI)+1.
          YYL(IYL,10)=YYL(IYL,10)+1.
*  pseudorapidity
          PT=SQRT(PTT)
          PTOT=SQRT(PTT+PHKK(3,I)**2)
          YLPS=LOG((PTOT+PHKK(3,I))/PT)
          IYL=(YLPS+2.0)/DY + 1
          IF (IYL.LT.1)  IYL=1
          IF (IYL.GT.50) IYL=50
          IF (ICHHKK.NE.0) YYLPS(IYL,9)=YYLPS(IYL,9)+1.
          IF (ICHHKK.LT.0) YYLPS(IYL,11)=YYLPS(IYL,11)+1.
          IF(NI.GT.0) YYLPS(IYL,NI)=YYLPS(IYL,NI)+1.
          YYLPS(IYL,10)=YYLPS(IYL,10)+1.
C                         HBOOK HISTOGRAMS
C                          HBOOK HISTOGRAMS
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
           P4P4P4 (1) = PHKK(1,I)
           P4P4P4 (2) = PHKK(2,I)
           P4P4P4 (3) = PHKK(3,I)
           P4P4P4 (4) = PHKK(4,I)
           ITIF = NRHKK
            XFL=PHKK(4,I)/EEO	
           XFXFXF=XFL
           CCCHRG = ICHHKK
       	 IF(IHBOOK.EQ.1)CALL PLOMB (3,P4P4P4,
     *    CCCHRG,XFXFXF,ITIF,IJPROJ)
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ENDIF
21    CONTINUE
*  total multiplicity
      AVMULT(30)=AVMULT(30) + NHAD
      ANCHSQ=ANCHSQ+AVMULC**2
      ANHAD=NHAD
C                         HBOOK HISTOGRAMS
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        IF(IHBOOK.EQ.1)CALL PLOMB (4,P4P4P4,
     *                  CCCHRG,XFXFXF,1,IJPROJ)
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      RETURN
C
C--------------------------------------------------------------------
C
  3   CONTINUE
* output of final results
C      WRITE(6,'(1H1,50(1H*))')
C      WRITE(6,'(/10X,A/)')'OUTP.FR.DISTR 0.75.lt.y.lt.3., pt.gt.o.5'
C      WRITE(6,'(50(1H*))')
C                           NORMALIZE AND PRINT COUNTERS
      IF (ANNVV.GT.1.)THEN
        PTVV=PTVV/ANNVV
        EEVV=EEVV/ANNVV
      ENDIF
      IF (ANNSS.GT.1.)THEN
        PTSS=PTSS/ANNSS
        EESS=EESS/ANNSS
      ENDIF
      IF (ANNSV.GT.1.)THEN
        PTSV=PTSV/ANNSV
        EESV=EESV/ANNSV
      ENDIF
      IF (ANNVS.GT.1.)THEN
        PTVS=PTVS/ANNVS
        EEVS=EEVS/ANNVS
      ENDIF
      IF (ANNCC.GE.1.)THEN
        PTCC=PTCC/ANNCC
        EECC=EECC/ANNCC
      ENDIF
      IF (ANNVD.GT.1.)THEN
        PTVD=PTVD/ANNVD
        EEVD=EEVD/ANNVD
      ENDIF
      IF (ANNDV.GT.1.)THEN
        PTDV=PTDV/ANNDV
        EEDV=EEDV/ANNDV
      ENDIF
      IF (ANNSD.GT.1.)THEN
        PTSD=PTSD/ANNSD
        EESD=EESD/ANNSD
      ENDIF
      IF (ANNDS.GT.1.)THEN
        PTDS=PTDS/ANNDS
        EEDS=EEDS/ANNDS
      ENDIF
      IF (ANNHH.GE.1.)THEN
        PTHH=PTHH/ANNHH
        EEHH=EEHH/ANNHH
      ENDIF
      IF (ANNZZ.GE.1.)THEN
        PTZZ=PTZZ/ANNZZ
        EEZZ=EEZZ/ANNZZ
      ENDIF
      IF (NHKKH1.GT.1.)THEN
        ANNVV=ANNVV/NHKKH1
        ANNSS=ANNSS/NHKKH1
        ANNSV=ANNSV/NHKKH1
        ANNVS=ANNVS/NHKKH1
        ANNCC=ANNCC/NHKKH1
        ANNHH=ANNHH/NHKKH1
        ANNZZ=ANNZZ/NHKKH1
        ANNDV=ANNDV/NHKKH1
        ANNVD=ANNVD/NHKKH1
        ANNDS=ANNDS/NHKKH1
        ANNSD=ANNSD/NHKKH1
        ANCHSQ=ANCHSQ/NHKKH1
C       WRITE (6,7431)ANNVV,PTVV,EEVV,
C    *              ANNSS,PTSS,EESS,
C    *              ANNSV,PTSV,EESV,
C    *              ANNVS,PTVS,EEVS,
C    *              ANNDV,PTDV,EEDV,
C    *              ANNVD,PTVD,EEVD,
C    *              ANNDS,PTDS,EEDS,
C    *              ANNSD,PTSD,EESD,
C    *              ANNHH,PTHH,EEHH,
C    *              ANNZZ,PTZZ,EEZZ,
C    *              ANNCC,PTCC,EECC
C7431   FORMAT ('  VV CHAINS NN,PT ECM: ',3F12.4/
C    *          '  SS CHAINS NN,PT ECM: ',3F12.4/
C    *          '  SV CHAINS NN,PT ECM: ',3F12.4/
C    *          '  VS CHAINS NN,PT ECM: ',3F12.4/
C    *          '  DV CHAINS NN,PT ECM: ',3F12.4/
C    *          '  VD CHAINS NN,PT ECM: ',3F12.4/
C    *          '  DS CHAINS NN,PT ECM: ',3F12.4/
C    *          '  SD CHAINS NN,PT ECM: ',3F12.4/
C    *          '  HH CHAINS NN,PT ECM: ',3F12.4/
C    *          '  ZZ CHAINS NN,PT ECM: ',3F12.4/
C    *          '  CC CHAINS NN,PT ECM: ',3F12.4)
      ENDIF
C
C                            DISTRIBUTION FOR
C                            NUMBER OF GENERATIONS INSIDE THE NUCLEUS
C                            (NORMALIZED TO 1)
      AVGNOR=0.0
      AVGEN=0.0
      DO 801 IGE=1,50
        AVGNOR=AVGNOR + DISGEN(IGE)
        AVGEN=AVGEN + IGE*DISGEN(IGE)
 801  CONTINUE
      DO 802 IGE=1,50
        XGEN(IGE)=FLOAT(IGE) - 0.5
        DISGEN(IGE)=DISGEN(IGE)/AVGNOR
 802  CONTINUE
      AVGEN=AVGEN/AVGNOR
C     WRITE(*,'(A,1PE10.2//A)')
C    &      ' AVERAGE NUMBER OF GENERATIONS INSIDE THE NUCLEUS =',
C    &      AVGEN,
C    &      ' CORRESPONDING DISTRIBUTION (NORMALIZED TO 1.0)'
C     CALL PLOT(XGEN,DISGEN,50,1,50,ZERO,ONEONE,ZERO,HUNDTH)
C
      DO 310 I=1,NUMTYP
        AVMULT(I)=AVMULT(I)/NHKKH1
        AVE(I)=AVE(I)/NHKKH1
        AVEPT(I)=AVEPT(I)/(NHKKH1*AVMULT(I))
310   CONTINUE
      FFFF2=SQRT(ANCHSQ-AVMULT(27)**2)/AVMULT(27)
C     WRITE(6,7772)FFFF2,ANCHSQ,AVMULT(27)
C7772 FORMAT(' FFFF2,ANCHSQ,AVMULT(27):',3E15.5)
      WRITE(6, 64)
   64 FORMAT(' PARTICLE REF,CHAR,IBAR, MASS      AVERAGE',
     *' ENERGY, MULTIPLICITY, INELASTICITY')
      DO 62 I=1,NUMTYP
        AKE(I)=AVE(I)/EEO
        WRITE(6, 63) ANAME(I),I,IICH(I),IIBAR(I),AAM(I),
     *               AVE(I),AVMULT(I),AKE(I),AVEPT(I)
   63   FORMAT (' ',A8,3I5,F10.3,4F15.6)
   62 CONTINUE
C
      DO 33 I=1,20
      DO 35 J=1,50
        YYL(J,I)  =YYL(J,I)  /(NHKKH1*DY)
        YYLPS(J,I)=YYLPS(J,I)/(NHKKH1*DY)
35    CONTINUE
33    CONTINUE
      WRITE(6,'(1H1,11(A/))')
     &   '  Conventions for rapidity distributions: ',
     &   '     1  Proton                       11  NEGATIVES          ',
     &   '     2  Neutron                      12  NBAR=9             ',
     &   '     3  PI+=13                       13  LAMBDA=17          ',
     &   '     4  PI-=14                       14  SIGMA==20,21,22    ',
     &   '     5  K+ =15                       15  LAMBDABAR=18       ',
     &   '     6  K- =16                       16  SIGMABAR=99,100,101',
     &   '     7  neutral kaons=12,19,24,25    17  THETA=97,98        ',
     &   '     8  pbar=2                       18  THETABAR=102,103   ',
     &   '     9  charged hadrons ',
     &   '    10  total hadrons '
      WRITE(6, 66)
   66 FORMAT(' RAPIDITY DISTRIBUTION')
      WRITE(6,302)
 302  FORMAT ('   (first number gives the lower bin limit)')
      DO 36 J=1,50
        WRITE(6, 37) XYL(J,1),(YYL(J,I),I=1,10)
37      FORMAT (F10.2,10E11.3)
36    CONTINUE
      WRITE(6, 66)
      WRITE(6,302)
      DO 306 J=1,50
      WRITE(6, 37) XYL(J,1),(YYL(J,I),I=11,20)
306   CONTINUE
C
      WRITE(6, 66 )
      CALL PLOT(XYL,YYL,1000,20,50,-TWOTWO,DY,ZERO,DELRAP)
*  rescaled rapidity distribution
      CALL PLOT(XYL,YYL,1000,20,50,-TWOTWO,DY,ZERO,5.*DELRAP)
*  logarithmic rapidity distribution
      DO 38 I=1,20
      DO 38 J=1,50
        YYL(J,I)=LOG10(ABS(YYL(J,I))+1.E-8)
 38   CONTINUE
      WRITE(6, 66)
      CALL PLOT(XYL,YYL,1000,20,50,-TWOTWO,DY,-TWOTWO,5.*HUNDTH)
C
      WRITE(6,301)
  301 FORMAT ('1 PSEUDORAPIDITY DISTRIBUTION')
      WRITE(6,302)
      DO 303 J=1,50
        WRITE(6,37) XYL(J,1),(YYLPS(J,I),I=1,10)
  303 CONTINUE
      WRITE(6,301)
      WRITE(6,302)
      DO 304 J=1,50
      WRITE(6,37) XYL(J,1),(YYLPS(J,I),I=11,20)
  304 CONTINUE
      WRITE(6,301)
      CALL PLOT(XYL,YYLPS,1000,20,50,-TWOTWO,DY,ZERO,DELRAP)
      CALL PLOT(XYL,YYLPS,1000,20,50,-TWOTWO,DY,ZERO,5.*DELRAP)
      RETURN
      END
*CMZ :  1.00/02 08/01/92  17.24.06  by  H.-J. Moehring & J. Ranft
*CMZ :  1.00/01 25/11/91  15.30.00  by  H.-J. M¶hring
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE DISTRC(IOP,NHKKH1,PO,IGENER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C***
C   REDUCED VERSION FOR SCORING RESULTS FROM PRIMARY AA INTERACTION
C   FROM /HKKEVT/
C   kinematics in cms 
C      - MULTIPLICITIES
C      - number of inc generations inside the nucleus
C      - RAPIDITY DISTRIBUTION         (without cuts)
C      - pseudorapidity distribution   (without cuts)
C***
C  conventions for rapidity distributions:
C       1  Proton                       11  NEGATIVES
C       2  Neutron                      12  NBAR=9
C       3  PI+=13                       13  LAMBDA=17
C       4  PI-=14                       14  SIGMA==20,21,22
C       5  K+ =15                       15  LAMBDABAR=18
C       6  K- =16                       16  SIGMABAR=99,100,101
C       7  neutral kaons=12,19,24,25    17  THETA=97,98
C       8  pbar=2                       18  THETABAR=102,103
C       9  charged hadrons
C      10  total hadrons
C-----------------------------------------------------------------------
*KEEP,HKKEVT.
c     INCLUDE (HKKEVT)
C     PARAMETER (NMXHKK=89998)
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
*KEEP,NSHMAK.
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C      COMMON /NSHMAK/ NNSHMA,NPSHMA,NTSHMA,NSHMAC
      COMMON /NSHMAK/ NNSHMA,NPSHMA,NTSHMA,NSHMAC,NSHMA2
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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
C                                   modified DPMJET
       COMMON /BUFUEH/ ANNVV,ANNSS,ANNSV,ANNVS,ANNCC,
     *                 ANNDV,ANNVD,ANNDS,ANNSD,
     *                 ANNHH,ANNZZ,
     *                 PTVV,PTSS,PTSV,PTVS,PTCC,PTDV,PTVD,PTDS,PTSD,
     *                 PTHH,PTZZ,
     *                 EEVV,EESS,EESV,EEVS,EECC,EEDV,EEVD,EEDS,EESD,
     *                 EEHH,EEZZ
     *                ,ANNDI,PTDI,EEDI
     *                ,ANNZD,ANNDZ,PTZD,PTDZ,EEZD,EEDZ
C---------------------
      CHARACTER*80 TITLE
      CHARACTER*8 PROJTY,TARGTY
C     COMMON /USER/TITLE,PROJTY,TARGTY,CMENER,ISTRUF
C    &            ,ISINGD,IDUBLD,SDFRAC,PTLAR
      COMMON /USER1/TITLE,PROJTY,TARGTY
      COMMON /USER2/CMENER,SDFRAC,PTLAR,ISTRUF,ISINGD,IDUBLD
      DIMENSION  XYL(50,20),YYL(50,20),YYLPS(50,20),INDX(28)
      DIMENSION DISGEN(50),XGEN(50)
      PARAMETER (NUMTYP=160)
      DIMENSION AVMULT(NUMTYP),AVE(NUMTYP),AKE(NUMTYP),AVEPT(NUMTYP)
      DATA DISGEN /50*0.0/
      DATA INDX/ 1, 8,-1,-1,-1,-1,-1, 2,12,-1,-1,7,3,4,5,6,13,15,7,14,
     *          14,14,19, 7, 7,-1,-1,-1/
      DATA NIPRIE /0/
C----------------------------------------------------------------------
C     CALL DISPT(IOP,NHKKH1,PO)
      GO TO (1,2,3),IOP
*  initialization call
1     CONTINUE
C                             initialize ALICE output
C     CALL SHINIT
C
      ANCHSQ=0.
      DELRAP=1.
      IF(IP.EQ.1)DELRAP=0.1
      NHKKH2=NHKKH1
      IF (NHKKH1.EQ.0)NHKKH2=1
      EEO=SQRT(PO**2+AAM(NHKKH2)**2)
      WRITE(6, 1001)EEO,PO,NHKKH2,AAM(NHKKH2)
 1001 FORMAT (' EEO',F10.2,F10.2,I10,F10.2)
      DY=0.4
C
      DO 11 I=1,20
        DO 13 J=1,50
          XYL(J,I)=-10.0 + (J-1)*DY
          YYL(J,I)=1.E-18
          YYLPS(J,I)=1.E-18
13      CONTINUE
11    CONTINUE
      DO 61 I=1,NUMTYP
        AVE(I)       =1.E-18
        AVEPT(I)       =1.E-18
        AVMULT(I)    =1.E-18
   61 CONTINUE
      RETURN
C-------------------------------------------------------------------
  2   CONTINUE
C                               SEWEW secondary interactions
      CALL SEWEW(1,NHKKH1)
      IEVT=IEVT+1
C                  output for ALICE
C     CALL SHEVUT(IEVT)
C
*  scoring of results from individual events
C                            NUMBER OF GENERATIONS INSIDE THE NUCLEUS
      IF(IGENER.LT.1.OR.IGENER.GT.50) IGENER=50
      DISGEN(IGENER)=DISGEN(IGENER) + 1.0
C
C     WRITE (18,1711)
C1711 FORMAT (' event px,py,pz,e,id ')
      NHAD=0
      AVMULC=0.
C                Here we include the decayed resonances into the
C                        multiplicity Table
      IPRIEV=0
      DO 521 I=NHKKH1,NHKK
        IF (ISTHKK(I).EQ.2)THEN
C         WRITE(18,1712)PHKK(1,I),PHKK(2,I),PHKK(3,I),Phkk(4,I),IDHKK(I)
C         NHAD=NHAD+1
          NRHKK=MCIHAD(IDHKK(I))
          IF (NRHKK.LE.0.OR.NRHKK.GT.210)THEN
            WRITE(6,1389)NRHKK,I,IDHKK(I),NHKKH1,NHKK
            NRHKK=1
          ENDIF
          NRE=NRHKK
          NREM=NRE
          ICHHKK=IICH(NRHKK)
	  IF(NRE.GT.30)THEN
	  IF(NRE.GT.160)GO TO 521
            AVEPT(NREM)=AVEPT(NREM) + PT
            AVMULT(NREM)=AVMULT(NREM) + 1.
	  ENDIF
        ENDIF
  521 CONTINUE
      DO 21 I=NHKKH1,NHKK
        IF (ISTHKK(I).EQ.1)THEN
C         WRITE(18,1712)PHKK(1,I),PHKK(2,I),PHKK(3,I),Phkk(4,I),IDHKK(I)
 1712     FORMAT (4E14.5,I8)
          IF(PHKK(4,I).GT.CMENER/2.D0)THEN
            IPRIEV=1
            NIPRIE=NIPRIE+1
          ENDIF 
          NHAD=NHAD+1
          NRHKK=MCIHAD(IDHKK(I))
          IF (NRHKK.LE.0.OR.NRHKK.GT.210)THEN
            WRITE(6,1389)NRHKK,I,IDHKK(I),NHKKH1,NHKK
 1389       FORMAT (' DISTR: NRHKK ERROR ',5I10)
            NRHKK=1
          ENDIF
          NRE=NRHKK
          NREM=NRE
          ICHHKK=IICH(NRHKK)
*  rapidity
          PTT=PHKK(1,I)**2+PHKK(2,I)**2+0.000001
          PT=SQRT(PTT)
C         IF(PT.LT.0.5)GO TO 21
          AMT=SQRT(PTT+PHKK(5,I)**2)
C***      AMT=SQRT(PTT+AMPI**2)
C         ETOT=SQRT(AMT**2 + PHKK(3,I)**2)
C         YL=LOG((ETOT + PHKK(3,I))/AMT+1.E-18)
          YL=LOG((ABS(PHKK(3,I) + PHKK(4,I)))/AMT+1.E-18)
C         IF(YL.LT.0.75.OR.YL.GT.3.) GO TO 21  
          IF (NRE.GT.25) NRE=28
          IF (NRE.LT. 1) NRE=28
          IF(NREM.GT.NUMTYP) NREM=28
          IF(NREM.LT.1) NREM=28
          NI=INDX(NRE)
          IF (NRHKK.LE.101.AND.NRHKK.GE.99) NI=16
          IF (NRHKK.EQ.97.OR.NRHKK.EQ.98)   NI=17
          IF (NRHKK.EQ.102.OR.NRHKK.EQ.103) NI=18
          AVE(NREM)=AVE(NREM) + PHKK(4,I)
          AVEPT(NREM)=AVEPT(NREM) + PT
          AVMULT(NREM)=AVMULT(NREM) + 1.
C    TOTAL=30
          AVE(30)=AVE(30) + PHKK(4,I)
          AVEPT(30)=AVEPT(30) + PT
C    CHARGED=27
          IF (ICHHKK.NE.0)  THEN
            AVE(27)=AVE(27) + PHKK(4,I)
            AVEPT(27)=AVEPT(27) + PT
            AVMULC=AVMULC + 1.
            AVMULT(27)=AVMULT(27) + 1.
          ENDIF
          IYL=(YL+10.0)/DY + 1
          IF (IYL.LT.1)  IYL=1
          IF (IYL.GT.50) IYL=50
          IF (ICHHKK.NE.0) THEN
            YYL(IYL,9)=YYL(IYL,9)+1.
          ENDIF
          IF (ICHHKK.LT.0)    YYL(IYL,11)=YYL(IYL,11)+1.
          IF(NI.GT.0) YYL(IYL,NI)=YYL(IYL,NI)+1.
          YYL(IYL,10)=YYL(IYL,10)+1.
*  pseudorapidity
          PT=SQRT(PTT)
          PTOT=SQRT(PTT+PHKK(3,I)**2)
          YLPS=LOG((PTOT+PHKK(3,I))/PT)
          IYL=(YLPS+10.0)/DY + 1
          IF (IYL.LT.1)  IYL=1
          IF (IYL.GT.50) IYL=50
          IF (ICHHKK.NE.0) YYLPS(IYL,9)=YYLPS(IYL,9)+1.
          IF (ICHHKK.LT.0) YYLPS(IYL,11)=YYLPS(IYL,11)+1.
          IF(NI.GT.0) YYLPS(IYL,NI)=YYLPS(IYL,NI)+1.
          YYLPS(IYL,10)=YYLPS(IYL,10)+1.
        ENDIF
21    CONTINUE
      IF(IPRIEV.EQ.-1)THEN
      IF(NIPRIE.LT.20)THEN
         WRITE(6,*)' IPRIEV = 1 '
         DO 120 IHKK=1,NHKK
          WRITE(6,1050) IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),
     +    JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +    (PHKK(KHKK,IHKK),KHKK=1,5)
120   CONTINUE
 1050 FORMAT (I6,I4,5I6,5E16.8)
      ENDIF
      ENDIF
*  total multiplicity
      AVMULT(30)=AVMULT(30) + NHAD
      ANCHSQ=ANCHSQ+AVMULC**2
      ANHAD=NHAD
      RETURN
C
C--------------------------------------------------------------------
C
  3   CONTINUE
C                             output for ALICE
C     CALL SHEND
C
* output of final results
C     WRITE(6,'(1H1,50(1H*))')
C     WRITE(6,'(/10X,A/)')'OUTP.FR.DISTR 0.75.lt.y.lt.3., pt.gt.o.5'
C     WRITE(6,'(50(1H*))')
C                           NORMALIZE AND PRINT COUNTERS
      IF (ANNVV.GT.1.)THEN
        PTVV=PTVV/ANNVV
        EEVV=EEVV/ANNVV
      ENDIF
      IF (ANNSS.GT.1.)THEN
        PTSS=PTSS/ANNSS
        EESS=EESS/ANNSS
      ENDIF
      IF (ANNSV.GT.1.)THEN
        PTSV=PTSV/ANNSV
        EESV=EESV/ANNSV
      ENDIF
      IF (ANNVS.GT.1.)THEN
        PTVS=PTVS/ANNVS
        EEVS=EEVS/ANNVS
      ENDIF
      IF (ANNCC.GE.1.)THEN
        PTCC=PTCC/ANNCC
        EECC=EECC/ANNCC
      ENDIF
      IF (NHKKH1.GT.1.)THEN
        ANNVV=ANNVV/NHKKH1
        ANNSS=ANNSS/NHKKH1
        ANNSV=ANNSV/NHKKH1
        ANNVS=ANNVS/NHKKH1
        ANNCC=ANNCC/NHKKH1
        ANCHSQ=ANCHSQ/NHKKH1
        WRITE (6,7431)ANNVV,PTVV,EEVV,
     *              ANNSS,PTSS,EESS,
     *              ANNSV,PTSV,EESV,
     *              ANNVS,PTVS,EEVS,
     *              ANNCC,PTCC,EECC
 7431   FORMAT ('  VV CHAINS NN,PT ECM: ',3F12.4/
     *          '  SS CHAINS NN,PT ECM: ',3F12.4/
     *          '  SV CHAINS NN,PT ECM: ',3F12.4/
     *          '  VS CHAINS NN,PT ECM: ',3F12.4/
     *          '  CC CHAINS NN,PT ECM: ',3F12.4)
      ENDIF
C
C                            DISTRIBUTION FOR
C                            NUMBER OF GENERATIONS INSIDE THE NUCLEUS
C                            (NORMALIZED TO 1)
      AVGNOR=0.0
      AVGEN=0.0
      DO 801 IGE=1,50
        AVGNOR=AVGNOR + DISGEN(IGE)
        AVGEN=AVGEN + IGE*DISGEN(IGE)
 801  CONTINUE
      DO 802 IGE=1,50
        XGEN(IGE)=FLOAT(IGE) - 0.5
        DISGEN(IGE)=DISGEN(IGE)/AVGNOR
 802  CONTINUE
      AVGEN=AVGEN/AVGNOR
C     WRITE(*,'(A,1PE10.2//A)')
C    &      ' AVERAGE NUMBER OF GENERATIONS INSIDE THE NUCLEUS =',
C    &      AVGEN,
C    &      ' CORRESPONDING DISTRIBUTION (NORMALIZED TO 1.0)'
C     CALL PLOT(XGEN,DISGEN,50,1,50,0.0,1.0,0.0,0.01)
C
      DO 310 I=1,NUMTYP
        AVMULT(I)=AVMULT(I)/NHKKH1
        AVE(I)=AVE(I)/NHKKH1
        AVEPT(I)=AVEPT(I)/(NHKKH1*AVMULT(I))
310   CONTINUE
      FFFF2=SQRT(ANCHSQ-AVMULT(27)**2)/AVMULT(27)
      WRITE(6,7772)FFFF2,ANCHSQ,AVMULT(27)
 7772 FORMAT(' FFFF2,ANCHSQ,AVMULT(27):',3E15.5)
      WRITE(6, 64)
   64 FORMAT(' PARTICLE REF,CHAR,IBAR, MASS      AVERAGE',
     *' ENERGY, MULTIPLICITY, INELASTICITY')
      DO 62 I=1,NUMTYP
        AKE(I)=AVE(I)/CMENER
        WRITE(6, 63) ANAME(I),I,IICH(I),IIBAR(I),AAM(I),
     *               AVE(I),AVMULT(I),AKE(I),AVEPT(I)
   63   FORMAT (' ',A8,3I5,F10.3,4E15.5)
   62 CONTINUE
C
      DO 33 I=1,20
      DO 35 J=1,50
        YYL(J,I)  =YYL(J,I)  /(NHKKH1*DY)
        YYLPS(J,I)=YYLPS(J,I)/(NHKKH1*DY)
35    CONTINUE
33    CONTINUE
      WRITE(6,'(1H1,11(A/))')
     &   '  Conventions for rapidity distributions: ',
     &   '     1  Proton                       11  NEGATIVES          ',
     &   '     2  Neutron                      12  NBAR=9             ',
     &   '     3  PI+=13                       13  LAMBDA=17          ',
     &   '     4  PI-=14                       14  SIGMA==20,21,22    ',
     &   '     5  K+ =15                       15  LAMBDABAR=18       ',
     &   '     6  K- =16                       16  SIGMABAR=99,100,101',
     &   '     7  neutral kaons=12,19,24,25    17  THETA=97,98        ',
     &   '     8  pbar=2                       18  THETABAR=102,103   ',
     &   '     9  charged hadrons ',
     &   '    10  total hadrons '
      WRITE(6, 66)
   66 FORMAT(' RAPIDITY DISTRIBUTION')
      WRITE(6,302)
 302  FORMAT ('   (first number gives the lower bin limit)')
      DO 36 J=1,50
        WRITE(6, 37) XYL(J,1),(YYL(J,I),I=1,10)
37      FORMAT (F10.2,10E11.3)
36    CONTINUE
      WRITE(6, 66)
      WRITE(6,302)
      DO 306 J=1,50
      WRITE(6, 37) XYL(J,1),(YYL(J,I),I=11,20)
306   CONTINUE
C
      WRITE(6, 66 )
      CALL PLOT(XYL,YYL,1000,20,50,-10.,DY,0.,DELRAP)
*  rescaled rapidity distribution
      CALL PLOT(XYL,YYL,1000,20,50,-10.,DY,0.,5.*DELRAP)
*  logarithmic rapidity distribution
      DO 38 I=1,20
      DO 38 J=1,50
        YYL(J,I)=LOG10(ABS(YYL(J,I))+1.E-8)
 38   CONTINUE
      WRITE(6, 66)
      CALL PLOT(XYL,YYL,1000,20,50,-10.,DY,-2.0,0.05)
C
      WRITE(6,301)
  301 FORMAT ('1 PSEUDORAPIDITY DISTRIBUTION')
      WRITE(6,302)
      DO 303 J=1,50
        WRITE(6,37) XYL(J,1),(YYLPS(J,I),I=1,10)
  303 CONTINUE
      WRITE(6,301)
      WRITE(6,302)
      DO 304 J=1,50
      WRITE(6,37) XYL(J,1),(YYLPS(J,I),I=11,20)
  304 CONTINUE
      WRITE(6,301)
      CALL PLOT(XYL,YYLPS,1000,20,50,-10.,DY,0.,DELRAP)
      CALL PLOT(XYL,YYLPS,1000,20,50,-10.,DY,0.,5.*DELRAP)
      RETURN
      END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE DISTRP(IOP,NHKKH1,PO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      RETURN
      END
C--------------------------------------------Dummies----------
      SUBROUTINE DISTCO(IOP,IJPROJ,PPN,IDUMMY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      RETURN
      END
      SUBROUTINE DISTPA(IOP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      RETURN
      END
      SUBROUTINE DISRES(IOP,IJPROJ,PPN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      RETURN
      END
      SUBROUTINE EDENSI(I,J)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      RETURN
      END
      SUBROUTINE PLOMB(I,PP,CHAR,XF,ITIF,IJPROJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      RETURN
      END
   
      SUBROUTINE PLOMBC(I,PP,CHAR,XF,ITIF,IJPROJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      RETURN
      END
   
      SUBROUTINE DISPT(IOP,NHKKH1,PO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C*** 1=P,  2=N, 3=PI+, 4=PI-, 5=PIO, 6=GAM+HYP, 7=K, 8=ANUC, 9=CHARGED
C*** 10=TOT, 11=TOTHAD
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
*KEEP,NSHMAK.
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C      COMMON /NSHMAK/ NNSHMA,NPSHMA,NTSHMA,NSHMAC
      COMMON /NSHMAK/ NNSHMA,NPSHMA,NTSHMA,NSHMAC,NSHMA2
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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

C                                   modified DPMJET
       COMMON /BUFUES/ BNNVV,BNNSS,BNNSV,BNNVS,BNNCC,
     *                 BNNDV,BNNVD,BNNDS,BNNSD,
     *                 BNNHH,BNNZZ,
     *                 BPTVV,BPTSS,BPTSV,BPTVS,BPTCC,BPTDV,
     *                 BPTVD,BPTDS,BPTSD,
     *                 BPTHH,BPTZZ,
     *                 BEEVV,BEESS,BEESV,BEEVS,BEECC,BEEDV,
     *                 BEEVD,BEEDS,BEESD,
     *                 BEEHH,BEEZZ
     *                ,BNNDI,BPTDI,BEEDI
     *                ,BNNZD,BNNDZ,BPTZD,BPTDZ,BEEZD,BEEDZ
       COMMON /NCOUCS/ BCOUVV,BCOUSS,BCOUSV,BCOUVS,
     *                 BCOUZZ,BCOUHH,BCOUDS,BCOUSD,
     *                 BCOUDZ,BCOUZD,BCOUDI,
     *                 BCOUDV,BCOUVD,BCOUCC
       COMMON /BUFUEH/ ANNVV,ANNSS,ANNSV,ANNVS,ANNCC,
     *                 ANNDV,ANNVD,ANNDS,ANNSD,
     *                 ANNHH,ANNZZ,
     *                 PTVV,PTSS,PTSV,PTVS,PTCC,PTDV,PTVD,PTDS,PTSD,
     *                 PTHH,PTZZ,
     *                 EEVV,EESS,EESV,EEVS,EECC,EEDV,EEVD,EEDS,EESD,
     *                 EEHH,EEZZ
     *                ,ANNDI,PTDI,EEDI
     *                ,ANNZD,ANNDZ,PTZD,PTDZ,EEZD,EEDZ
       COMMON /NCOUCH/ ACOUVV,ACOUSS,ACOUSV,ACOUVS,
     *                 ACOUZZ,ACOUHH,ACOUDS,ACOUSD,
     *                 ACOUDZ,ACOUZD,ACOUDI,
     *                 ACOUDV,ACOUVD,ACOUCC
C---------------------
      COMMON /EVENTA/IDUMTP
C
      DIMENSION PTP(50,20),PTY(50,20),EMTY(50,20),INDX(28)
      COMMON /SIGLA/SIGLAU
      COMMON /FINAL/IFINAL
C------------------------------------------------------------------
C------------------------------------------------------------------

      DATA INDX/1,8,10,10,10,10,7,2,7,10,10,7,3,4,5,6,
     *          7,7,7,7,7,7,7,7,7,7,7,7/
C   kinematics in lab (target rest system)
C***
C  conventions for rapidity distributions:
C       1  Proton                       11  NEGATIVES
C       2  Neutron                      12  NBAR=9
C       3  PI+=13                       13  LAMBDA=17
C       4  PI-=14                       14  SIGMA==20,21,22
C       5  K+ =15                       15  LAMBDABAR=18
C       6  K- =16                       16  SIGMABAR=99,100,101
C       7  neutral kaons=12,19,24,25    17  THETA=97,98
C       8  pbar=2                       18  THETABAR=102,103
C       9  charged hadrons              19  neutral Kaons 12,19,24,25
C      10  total hadrons
C-----------------------------------------------------------------------
      DATA KPL /1/
      DATA IEVL /0/
C----------------------------------------------------------------------
      ZERO=0
      ONEONE=1
      TWOTWO=2
      NIP=IP
      AIP=IP
      DELRAP=1.
      IF(IP.EQ.1)DELRAP=0.1
C
C-----------------------------
      GO TO (1,2,3),IOP
1     CONTINUE
      DXFL=0.02
      DY=0.4
      DPT=0.10
C     DIMENSION PTP(50,20),PTY(50,20),EMTY(50,20)
      DO 102 I=1,20
	DO 102 J=1,50
	  PTP(J,I)=(J-1)*DPT
	  PTY(J,I)=1.D-12
	  EMTY(J,I)=1.D-12
 102  CONTINUE
      NHKKH2=NHKKH1
      IF (NHKKH1.EQ.0)NHKKH2=1
      EEO=SQRT(PO**2+AAM(NHKKH2)**2)
      WRITE(6, 1001)EEO,PO,NHKKH2,AAM(NHKKH2)
 1001 FORMAT (' EEO,PO ',F13.2,F13.2,I10,F10.2)
2     CONTINUE
      NHAD=0
C                          HBOOK HISTOGRAMS
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     WRITE(6,'(A)')' before plomb'
      IF(IHBOOK.EQ.1)CALL PLOMB (2,P4P4P4,CCCHRG,XFXFXF,1,IJPROJ)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO 21 I=NHKKH1,NHKK
        IF (ISTHKK(I).EQ.1)THEN
C                                  j.r. temp for s-s
          PTT=PHKK(1,I)**2+PHKK(2,I)**2+0.00001
          PT=SQRT(PTT)+0.001
C         IF(PT.LT.0.5)GO TO 21
C
          NHAD=NHAD+1
          NRHKK=MCIHAD(IDHKK(I))
           IF (NRHKK.LE.0.OR.NRHKK.GT.210)THEN
             NRHKK=1
           ENDIF
          NRE=NRHKK
          ICHHKK=IICH(NRHKK)
          IF (NRE.GT.160) NRE=28
          IF (NRE.LT. 1) NRE=28
          NREX=NRE
          IF(NREX.GE.28)NREX=28
          NI=INDX(NREX)
          NIX=NI
          IF (NRHKK.EQ.9)NIX=12
          IF (NRHKK.EQ.17.OR.NRHKK.EQ.22)NIX=13
          IF (NRHKK.LE.22.AND.NRHKK.GE.20)NIX=14
          IF (NRHKK.EQ.18.OR.NRHKK.EQ.100)NIX=15
          IF (NRHKK.LE.101.AND.NRHKK.GE.99)NIX=16
          IF (NRHKK.EQ.98)NIX=17
          IF (NRHKK.EQ.103)NIX=18
          IF (NRHKK.EQ.12.OR.NRHKK.EQ.19)NIX=19
          IF (NRHKK.EQ.24.OR.NRHKK.EQ.25)NIX=19
          IF(NRE.EQ.28) NI=7
          PTT=PHKK(1,I)**2+PHKK(2,I)**2+0.00001
          PT=SQRT(PTT)+0.001
          AMT=SQRT(PTT+PHKK(5,I)**2)
          YL=LOG((ABS(PHKK(3,I)+SQRT(PHKK(3,I)**2+AMT**2)))/AMT+1.E-18)
          YLLPS=LOG((ABS(PHKK(3,I)+SQRT(PHKK(3,I)**2+PTT)))/SQRT(PTT)
     *               +1.E-18)
          KPL=1
          IPT=PT/DPT+1.
          IF (IPT.LT.1)IPT=1
          IF (IPT.GT.50) IPT=50
          IF (ICHHKK.NE.0)PTY(IPT,9)=PTY(IPT,9)+1./PT
          IF (ICHHKK.EQ.-1)PTY(IPT,11)=PTY(IPT,11)+1./PT
          PTY(IPT,NIX)=PTY(IPT,NIX)+1./PT
          PTY(IPT,10)=PTY(IPT,10)+1./PT
          IF(YL.GT.2.3.AND.YL.LE.3.)THEN
            IAMT=AMT/DPT+1.
            IF (IAMT.LT.1)IAMT=1
            IF (IAMT.GT.48) IAMT=48
            IF (ICHHKK.NE.0)EMTY(IAMT,9)=EMTY(IAMT,9)+1./AMT
            IF (ICHHKK.EQ.-1)EMTY(IAMT,11)=EMTY(IAMT,11)+1./AMT
            EMTY(IAMT,NIX)=EMTY(IAMT,NIX)+1./AMT
            EMTY(IAMT,10)=EMTY(IAMT,10)+1./AMT
            IF(PT.GT.1.0.AND.PT.LT.2.0)THEN
              IAMT=49
              IF (ICHHKK.NE.0)EMTY(IAMT,9)=EMTY(IAMT,9)+DPT
              IF (ICHHKK.EQ.-1)EMTY(IAMT,11)=EMTY(IAMT,11)+DPT
              EMTY(IAMT,NIX)=EMTY(IAMT,NIX)+DPT
              EMTY(IAMT,10)=EMTY(IAMT,10)+DPT
            ENDIF
            IF(AMT.GT.1.72.AND.YL.LE.2.6)THEN
              IAMT=50
              IF (ICHHKK.NE.0)EMTY(IAMT,9)=EMTY(IAMT,9)+DPT
              IF (ICHHKK.EQ.-1)EMTY(IAMT,11)=EMTY(IAMT,11)+DPT
              EMTY(IAMT,NIX)=EMTY(IAMT,NIX)+DPT
              EMTY(IAMT,10)=EMTY(IAMT,10)+DPT
            ENDIF
          ENDIF
        ENDIF
21    CONTINUE
      RETURN
C
C--------------------------------------------------------------------
C
  3   CONTINUE
C          Normalization and Printing
C------------------------------------------------------------------
C------------------------------------------------------------------
      KKK=13
C             Normalization
C
C                           NORMALIZE AND PRINT COUNTERS
      IF (BNNDI.GT.1.)THEN
        BPTDI=BPTDI/BNNDI
        BEEDI=BEEDI/BNNDI
      ENDIF
      IF (BNNVV.GT.1.)THEN
        BPTVV=BPTVV/BNNVV
        BEEVV=BEEVV/BNNVV
      ENDIF
      IF (BNNSS.GT.1.)THEN
        BPTSS=BPTSS/BNNSS
        BEESS=BEESS/BNNSS
      ENDIF
      IF (BNNSV.GT.1.)THEN
        BPTSV=BPTSV/BNNSV
        BEESV=BEESV/BNNSV
      ENDIF
      IF (BNNVS.GT.1.)THEN
        BPTVS=BPTVS/BNNVS
        BEEVS=BEEVS/BNNVS
      ENDIF
      IF (BNNCC.GE.1.)THEN
        BPTCC=BPTCC/BNNCC
        BEECC=BEECC/BNNCC
      ENDIF
      IF (BNNDV.GE.1.)THEN
        BPTDV=BPTDV/BNNDV
        BEEDV=BEEDV/BNNDV
      ENDIF
      IF (BNNVD.GE.1.)THEN
        BPTVD=BPTVD/BNNVD
        BEEVD=BEEVD/BNNVD
      ENDIF
      IF (BNNDS.GE.1.)THEN
        BPTDS=BPTDS/BNNDS
        BEEDS=BEEDS/BNNDS
      ENDIF
      IF (BNNDZ.GE.1.)THEN
        BPTDZ=BPTDZ/BNNDZ
        BEEDZ=BEEDZ/BNNDZ
      ENDIF
      IF (BNNHH.GE.1.)THEN
        BPTHH=BPTHH/BNNHH
        BEEHH=BEEHH/BNNHH
      ENDIF
      IF (BNNZZ.GE.1.)THEN
        BPTZZ=BPTZZ/BNNZZ
        BEEZZ=BEEZZ/BNNZZ
      ENDIF
      IF (BNNSD.GE.1.)THEN
        BPTSD=BPTSD/BNNSD
        BEESD=BEESD/BNNSD
      ENDIF
      IF (BNNZD.GE.1.)THEN
        BPTZD=BPTZD/BNNZD
        BEEZD=BEEZD/BNNZD
      ENDIF
      IF (NHKKH1.GT.1.)THEN
        BNNVV=BNNVV/NHKKH1
        BNNDI=BNNDI/NHKKH1
        BNNSS=BNNSS/NHKKH1
        BNNSV=BNNSV/NHKKH1
        BNNVS=BNNVS/NHKKH1
        BNNCC=BNNCC/NHKKH1
        BNNVD=BNNVD/NHKKH1
        BNNDV=BNNDV/NHKKH1
        BNNDS=BNNDS/NHKKH1
        BNNSD=BNNSD/NHKKH1
        BNNDZ=BNNDZ/NHKKH1
        BNNZD=BNNZD/NHKKH1
        BNNHH=BNNHH/NHKKH1
        BNNZZ=BNNZZ/NHKKH1
        BCOUVV=BCOUVV/NHKKH1
        BCOUSS=BCOUSS/NHKKH1
        BCOUSV=BCOUSV/NHKKH1
        BCOUVS=BCOUVS/NHKKH1
        BCOUZZ=BCOUZZ/NHKKH1
        BCOUHH=BCOUHH/NHKKH1
        BCOUDS=BCOUDS/NHKKH1
        BCOUSD=BCOUSD/NHKKH1
        BCOUDZ=BCOUDZ/NHKKH1
        BCOUZD=BCOUZD/NHKKH1
        BCOUDI=BCOUDI/NHKKH1
        BCOUDV=BCOUDV/NHKKH1
        BCOUVD=BCOUVD/NHKKH1
        BCOUCC=BCOUCC/NHKKH1
        WRITE (6,7431)BNNVV,BPTVV,BEEVV,BCOUVV,
     *              BNNSS,BPTSS,BEESS,BCOUSS,
     *              BNNSV,BPTSV,BEESV,BCOUSV,
     *              BNNVS,BPTVS,BEEVS,BCOUVS,
     *              BNNCC,BPTCC,BEECC,BCOUCC,
     *              BNNDV,BPTDV,BEEDV,BCOUDV,
     *              BNNVD,BPTVD,BEEVD,BCOUVD,
     *              BNNDS,BPTDS,BEEDS,BCOUDS,
     *              BNNSD,BPTSD,BEESD,BCOUSD,
     *              BNNDZ,BPTDZ,BEEDZ,BCOUDZ,
     *              BNNZD,BPTZD,BEEZD,BCOUZD,
     *              BNNHH,BPTHH,BEEHH,BCOUHH,
     *              BNNDI,BPTDI,BEEDI,BCOUDI,
     *              BNNZZ,BPTZZ,BEEZZ,BCOUZZ
 7431   FORMAT ('  VV CHAINS NN,PT ECM: ',4F12.4/
     *          '  SS CHAINS NN,PT ECM: ',4F12.4/
     *          '  SV CHAINS NN,PT ECM: ',4F12.4/
     *          '  VS CHAINS NN,PT ECM: ',4F12.4/
     *          '  CC CHAINS NN,PT ECM: ',4F12.4/
     *          '  DV CHAINS NN,PT ECM: ',4F12.4/
     *          '  VD CHAINS NN,PT ECM: ',4F12.4/
     *          '  DS CHAINS NN,PT ECM: ',4F12.4/
     *          '  SD CHAINS NN,PT ECM: ',4F12.4/
     *          '  DZ CHAINS NN,PT ECM: ',4F12.4/
     *          '  ZD CHAINS NN,PT ECM: ',4F12.4/
     *          '  HH CHAINS NN,PT ECM: ',4F12.4/
     *          '  DI CHAINS NN,PT ECM: ',4F12.4/
     *          '  ZZ CHAINS NN,PT ECM: ',4F12.4)
      ENDIF
C
C
C
C============================================================
      DO 33 I=1,20
      DO 35 J=1,50
        DO 335 JJJJ=0,1
  335   CONTINUE
        PTY(J,I)   =PTY(J,I)   /(NHKKH1*DPT)
        EMTY(J,I)   =EMTY(J,I)   /(NHKKH1*DPT)
35    CONTINUE
33    CONTINUE
      DO 5136 J=1,50
        WRITE(6, 37)(PTY(J,I),I=1,11)
5137    FORMAT (11E11.3)
  37    FORMAT (11E11.3)
5136  CONTINUE
      WRITE(6,3654)
 3654 FORMAT(' pt-Distribution')
      DO 3614 J=1,50
        WRITE(6, 37)PTP(J,1),(PTY  (J,I),I=1,10)
 3614 CONTINUE
      DO 5914 J=1,50
        WRITE(6, 37)PTP(J,1),(PTY  (J,I),I=11,20)
 5914 CONTINUE
      WRITE(6,4654)
 4654 FORMAT(' Et-Distribution')
      WRITE(6, 37)PTP(1,1),(EMTY  (1,I),I=1,10)
      DO 4614 J=2,50
        WRITE(6, 37)PTP(J-1,1),(EMTY  (J,I),I=1,10)
        WRITE(6, 37)PTP(J,1),(EMTY  (J,I),I=1,10)
 4614 CONTINUE
      WRITE(6, 37)PTP(1,1),(EMTY  (1,I),I=11,20)
      DO 6914 J=2,50
        WRITE(6, 37)PTP(J-1,1),(EMTY  (J,I),I=11,20)
        WRITE(6, 37)PTP(J,1),(EMTY  (J,I),I=11,20)
 6914 CONTINUE
      DO 514 J=1,50
        DO 312 I=1,20
          PTY(J,I)     =LOG10(ABS(PTY(J,I))+1.E-8)
          EMTY(J,I)     =LOG10(ABS(EMTY(J,I))+1.E-8)
312     CONTINUE
        PTY(J,10)     =PTY(J,11)
        EMTY(J,10)     =EMTY(J,11)
514   CONTINUE
      WRITE(6,3654)
      CALL PLOT(PTP,PTY,1000,20,50,ZERO,DPT,-ONEONE,0.05D0)
C     CALL PLOT(PTP,PTY,1000,20,50,ZERO,DPT,ONEONE,0.05D0)
      WRITE(6,4654)

      CALL PLOT(PTP,EMTY,1000,20,50,ZERO,DPT,-6.0D0,0.10D0)
C
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HISTOG(I)
C
      RETURN
      END
C----------------------------------------------------------------
      SUBROUTINE DISEVA(IOP,NHKKH1,PO,IGENER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
C------------------------------------------------------------------
C--------------------------------evaporation-----------------------
      DIMENSION ANEVA(4),PEVAP(50,2),XPEVAP(50,2),AMEVAP(250),
     * XMEVAP(250)
      DIMENSION ANEVAP(4),AMEVPP(250)
      COMMON /FINAL/IFINAL
      COMMON /NOMIJE/ PTMIJE(10),NNMIJE(10)
      COMMON /NOMIJU/ NNMIJU(10)
       CHARACTER*8  ANAME
       COMMON /DPAR/ANAME(210),AAM(210),GA(210),TAU(210),IICH(210),
     +IIBAR(210),K1(210),K2(210)
C------------------------------------------------------------------
C------------------------------------------------------------------
      GO TO (1,2,3),IOP
C------------------------------------------------------------------
C--------------------------------evaporation-----------------------
    1 CONTINUE
      DPEVA=0.02D0
      DIT=IT/50+1
      IF(DIT.LT.1.D0)DIT=1
      IF (IFINAL.EQ.0) THEN
        DO 3112 I=1,4
          ANEVA(I)=0
          ANEVAP(I)=0
 3112   CONTINUE
        DO 3123 I=1,250
          AMEVAP(I)=0
          AMEVPP(I)=0
          XMEVAP(I)=I
 3123   CONTINUE
        DO 3113 I=1,50
          DO 3113 II=1,2
            PEVAP(I,II)=0
            XPEVAP(I,II)=I*DPEVA
 3113   CONTINUE  
      ENDIF
C        s(shower), g (grey), h(heavy) and b (black) particles
C        n- (neg prongs) sp8,sp5 (slow protons) p=0.15-0.8, 0.15-0.5
       ANSHO=0 
       ANGRE=0
       ANGRE2=0
       ANHEA=0
       ANBLA=0
       ANMIN=0
       ANSP8=0
       ANSP5=0
       RETURN
C------------------------------------------------------------------
C------------------------------------------------------------------
    2 CONTINUE
C------------------------------------------------------------------
C--------------------------------evaporation-----------------------
C     WRITE(6,'(A)')' DISTR 2 '
	DO 1121 I=NHKKH1,NHKK
          IF (ISTHKK(I).EQ.1)THEN
         NRHKK=MCIHAD(IDHKK(I))
           ICHHKK=IICH(NRHKK)
            PPTT=SQRT(PHKK(1,I)**2+PHKK(2,I)**2+0.000001)
            PTT=PPTT**2
       	  PTOT=SQRT(PPTT**2+PHKK(3,I)**2+0.000001)
           IF(ICHHKK.NE.0)THEN
            BETP=PTOT/PHKK(4,I)
	     IF(ICHHKK.EQ.-1)ANMIN=ANMIN+1
             IF(BETP.GE.0.7D0)ANSHO=ANSHO+1
             IF(BETP.LE.0.7D0)ANHEA=ANHEA+1
             IF(BETP.LE.0.2D0)ANBLA=ANBLA+1
             IF(BETP.GE.0.2D0.AND.BETP.LE.0.7D0)ANGRE=ANGRE+1
             IF(BETP.GE.0.23D0.AND.BETP.LE.0.7D0)ANGRE2=ANGRE2+1
            ENDIF
	  ENDIF
 1121   CONTINUE
      IF (IFINAL.EQ.0) THEN
        DO 3114 I=NHKKH1,NHKK
          IF(ISTHKK(I).EQ.-1) THEN
C                                      Neutrons
            IF (IDHKK(I).EQ.2112)THEN
              PPTT=(PHKK(1,I)**2+PHKK(2,I)**2+0.000001)
              PTOT=SQRT(PPTT+PHKK(3,I)**2+0.000001)
	      IF(NOBAM(I).EQ.2)THEN
                ANEVA(2)=ANEVA(2)+1.D0
                IPTOT=PTOT/DPEVA
                IF(IPTOT.LT.1)IPTOT=1
                IF(IPTOT.GT.50)IPTOT=50
                PEVAP(IPTOT,2)=PEVAP(IPTOT,2)+1.E0
              ELSEIF(NOBAM(I).EQ.1)THEN
                ANEVAP(2)=ANEVAP(2)+1.D0
              ENDIF
C                                        Protons
            ELSEIF(IDHKK(I).EQ.2212)THEN
              PPTT=(PHKK(1,I)**2+PHKK(2,I)**2+0.000001)
              PTOT=SQRT(PPTT+PHKK(3,I)**2+0.000001)
	      IF(NOBAM(I).EQ.2)THEN
                ANEVA(1)=ANEVA(1)+1.D0
                IF(PTOT.GT.0.15D0.AND.PTOT.LE.0.8D0)ANSP8=ANSP8+1
                IF(PTOT.GT.0.15D0.AND.PTOT.LE.0.5D0)ANSP5=ANSP5+1
                IPTOT=PTOT/DPEVA
                IF(IPTOT.LT.1)IPTOT=1
                IF(IPTOT.GT.50)IPTOT=50
                PEVAP(IPTOT,1)=PEVAP(IPTOT,1)+1.E0
                BETP=PTOT/PHKK(4,I)
                IF(BETP.GE.0.7D0)ANSHO=ANSHO+1
                IF(BETP.LE.0.7D0)ANHEA=ANHEA+1
                IF(BETP.LE.0.23D0)ANBLA=ANBLA+1
                IF(BETP.GE.0.2D0.AND.BETP.LE.0.7D0)ANGRE=ANGRE+1
C               IF(BETP.GE.0.23D0.AND.BETP.LE.0.7D0)
C    *                                  	ANGRE2=ANGRE2+1
                IF(PTOT.GE.0.026D0.AND.PTOT.LE.0.375D0)
     *                                  	ANGRE2=ANGRE2+1
              ELSEIF(NOBAM(I).EQ.1)THEN
                ANEVAP(1)=ANEVAP(1)+1.D0
              ENDIF
C                                        Photons
            ELSEIF(IDHKK(I).EQ.22)THEN
              PPTT=(PHKK(1,I)**2+PHKK(2,I)**2+0.000001)
              PTOT=SQRT(PPTT+PHKK(3,I)**2+0.000001)
	      IF(NOBAM(I).EQ.2)THEN
                ANEVA(3)=ANEVA(3)+1.D0
              ELSEIF(NOBAM(I).EQ.1)THEN
                ANEVAP(3)=ANEVAP(3)+1.D0
              ENDIF
C                       transform deexcitation photons into
C                       electron neutrinos (pdg=12, bamjet=5
C                       to allow separate histograms of pi-0
C                       and deexcitation photons
	      ISTHKK(I)=1
	      IDHKK(I) =12
            ENDIF
          ENDIF
          IF(IDHKK(I).EQ.80000)THEN
C                                         heavy fragments
              PPTT=(PHKK(1,I)**2+PHKK(2,I)**2+0.000001)
              PTOT=SQRT(PPTT+PHKK(3,I)**2+0.000001)
	      IF(NOBAM(I).EQ.2)THEN
                ANEVA(4)=ANEVA(4)+1.D0
                IDIT=IDRES(I)
                IF(IDIT.EQ.0)GO TO 3115
                AMEVAP(IDIT)=AMEVAP(IDIT)+1.D0
                ANHEA=ANHEA+1
                ANBLA=ANBLA+1
              ELSEIF(NOBAM(I).EQ.1)THEN
                ANEVAP(4)=ANEVAP(4)+1.D0
                IDIT=IDRES(I)
                IF(IDIT.EQ.0)GO TO 3115
                AMEVPP(IDIT)=AMEVPP(IDIT)+1.D0
              ENDIF
 3115       CONTINUE
          ENDIF
 3114   CONTINUE
      ENDIF
       RETURN
C------------------------------------------------------------------
C------------------------------------------------------------------
    3 CONTINUE
C------------------------------------------------------------------
C--------------------------------evaporation-----------------------
      IF (IFINAL.EQ.0) THEN
	WRITE(6,'(A)')' Target fragments'
        DO 3117 I=1,4
          ANEVA(I)=ANEVA(I)/NHKKH1
 3117   CONTINUE
        WRITE(6,'(A,F10.3)')' Number of evap. protons:  ',ANEVA(1)
        WRITE(6,'(A,F10.3)')' Number of evap. neutrans: ',ANEVA(2)
        WRITE(6,'(A,F10.3)')' Number of ex. gammas:     ',ANEVA(3)
        WRITE(6,'(A,F10.3)')' Number of heavy fragments:',ANEVA(4)
	WRITE(6,'(A)')' Projectile fragments'
        DO 3118 I=1,4
          ANEVAP(I)=ANEVAP(I)/NHKKH1
 3118   CONTINUE
        WRITE(6,'(A,F10.3)')
     *	' Number of evap. protons:  ',ANEVAP(1)
        WRITE(6,'(A,F10.3)')
     *	' Number of evap. neutrans: ',ANEVAP(2)
        WRITE(6,'(A,F10.3)')
     *	' Number of ex. gammas:     ',ANEVAP(3)
        WRITE(6,'(A,F10.3)')
     *	' Number of heavy fragments:',ANEVAP(4)
        ANSHO=ANSHO/NHKKH1
        ANGRE=ANGRE/NHKKH1
        ANGRE2=ANGRE2/NHKKH1
        ANHEA=ANHEA/NHKKH1
        ANBLA=ANBLA/NHKKH1
        ANMIN=ANMIN/NHKKH1
        ANSP8=ANSP8/NHKKH1
        ANSP5=ANSP5/NHKKH1
        WRITE(6,'(A,F10.3)')' Number of shower part.:',ANSHO
        WRITE(6,'(A,F10.3)')' Number of grey   part.:',ANGRE
        WRITE(6,'(A,F10.3)')' Number of grey23 part.:',ANGRE2
        WRITE(6,'(A,F10.3)')' Number of heavy  part.:',ANHEA
        WRITE(6,'(A,F10.3)')' Number of black  part.:',ANBLA
        WRITE(6,'(A,F10.3)')' Number of neg.ch.part.:',ANMIN
        WRITE(6,'(A,F10.3)')' Number of slow 8 prot.:',ANSP8
        WRITE(6,'(A,F10.3)')' Number of slow 5 prot.:',ANSP5
        DO 3128 I=1,250
          AMEVAP(I)=AMEVAP(I)/NHKKH1
          AMEVPP(I)=AMEVPP(I)/NHKKH1
	WRITE(6,'(F12.4,2E15.5)')XMEVAP(I),AMEVAP(I),AMEVPP(I)
          AMEVAP(I)=LOG10(AMEVAP(I)+1.D-9)
          AMEVPP(I)=LOG10(AMEVPP(I)+1.D-9)
 3128   CONTINUE
        WRITE(6,'(A)')' Heavy fragment spectrum'
        DIT=IT/50+1
        IF(DIT.LT.1.D0)DIT=1
        CALL PLOT(XMEVAP,AMEVAP,IT,1,IT,0.,DIT,-5.,0.1)
        DIP=IP/50+1
        IF(DIP.LT.1.D0)DIP=1
        CALL PLOT(XMEVAP,AMEVPP,IP,1,IP,0.,DIP,-5.,0.1)
        DO 3119 I=1,50
          DO 3119 II=1,2
            PEVAP(I,II)=LOG10(PEVAP(I,II)/(NHKKH1*DPEVA)+1.D-9)
 3119   CONTINUE  
      ENDIF
        WRITE(6,'(A)')' Ev. p and n momentum  spectrum'
        CALL PLOT(XPEVAP,PEVAP,100,2,50,0.,DPEVA,-5.,0.1)
C------------------------------------------------------------------
C------------------------------------------------------------------
      RETURN
      END
