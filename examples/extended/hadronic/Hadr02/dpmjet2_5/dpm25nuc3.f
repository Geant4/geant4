      SUBROUTINE SAPTRE(AM1,G1,BGX1,BGY1,BGZ1,
     &                  AM2,G2,BGX2,BGY2,BGZ2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C                   SELECT PT FOR CHAIN PAIRS, WHICH ARE RESONANCES
      B3=4.
      E1=G1*AM1
      PX1=BGX1*AM1
      PY1=BGY1*AM1
      PZ1=BGZ1*AM1
      E2=G2*AM2
      PX2=BGX2*AM2
      PY2=BGY2*AM2
      PZ2=BGZ2*AM2
C                   SAMPLE TRANSVERSE MOMENTUM LIKE IN BAMJET
C       ES DEFINED AS ES=SQRT(PT**2+AM**2)-AM
      ESMAX1=E1-AM1
      ESMAX2=E2-AM2
      ESMAX=MIN(ESMAX1,ESMAX2)
      IF(ESMAX.LE.0.05D0) RETURN
      HMA=AM1
      IF (B3*ESMAX.GT.60.D0)THEN
        EXEB=0.
      ELSE
        EXEB=EXP(-B3*ESMAX)
      ENDIF
      BEXP=HMA*(1.-EXEB)/B3
      AXEXP=(1.D0-(B3*ESMAX-1.D0)*EXEB)/B3**2
      WA=AXEXP/(BEXP+AXEXP)
      XAB=RNDM(WU)
   10 CONTINUE
      IF (XAB.LT.WA)THEN
        X=RNDM(V)
        Y=RNDM(V)
        ES=-2./(B3**2)*LOG(X*Y+1.E-7)
      ELSE
        X=RNDM(V)
        ES=ABS(-LOG(X+1.E-7)/B3)
      END IF
      IF(ES.GT.ESMAX)                                             GOTO10
      ES=ES+HMA
      HPS=SQRT((ES-HMA)*(ES+HMA))
   20 CONTINUE
      CALL DSFECF(SFE,CFE)
      SIP=SFE
      COP=CFE
      HPX=HPS*COP
      HPY=HPS*SIP
      PZ1NSQ=PZ1**2-HPS**2-2.*PX1*HPX-2.*PY1*HPY
      PZ2NSQ=PZ2**2-HPS**2+2.*PX2*HPX+2.*PY2*HPY
      IF(PZ1NSQ.LT.0.001D0.OR.PZ2NSQ.LT.0.001D0) RETURN
      PZ1=SIGN(SQRT(PZ1NSQ),PZ1)
      PZ2=SIGN(SQRT(PZ2NSQ),PZ2)
      PX1=PX1+HPX
      PY1=PY1+HPY
      PX2=PX2-HPX
      PY2=PY2-HPY
      BGX1=PX1/AM1
      BGY1=PY1/AM1
      BGZ1=PZ1/AM1
      BGX2=PX2/AM2
      BGY2=PY2/AM2
      BGZ2=PZ2/AM2
C     WRITE(6,1001) HPX,HPY
C1001 FORMAT(' HPX,HPY ',2F10.3)
      RETURN
      END
*-- Author :
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE SLTRAF(GA,BGA,EIN,PZIN,EOUT,PZOUT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PZOUT=GA*PZIN - BGA*EIN
      EOUT=GA*EIN - BGA*PZIN
      RETURN
      END
*-- Author :
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE NUCMOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C***
C   FERMI-MOMENTA FOR ALL NUCLEONS
C                 TRANSFORMED INTO NN-CMS
C   FOR INCIDENT HADRONS USE CMS MOMENTUM
C***
*KEEP,INTMX.
      PARAMETER (INTMX=2488,INTMD=252)
*KEEP,DXQX.
C     INCLUDE (XQXQ)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON /DXQX/ XPVQ(248),XPVD(248),XTVQ(248),XTVD(248), XPSQ
     +(INTMX),XPSAQ(INTMX),XTSQ(INTMX),XTSAQ(INTMX)
     *                ,XPSU(248),XTSU(248)
     *                ,XPSUT(248),XTSUT(248)
*KEEP,INTNEW.
      COMMON /INTNEW/ NVV,NSV,NVS,NSS,NDV,NVD,NDS,NSD,
     +IXPV,IXPS,IXTV,IXTS, INTVV1(248),
     +INTVV2(248),INTSV1(248),INTSV2(248), INTVS1(248),INTVS2(248),
     +INTSS1(INTMX),INTSS2(INTMX),
     +INTDV1(248),INTDV2(248),INTVD1(248),INTVD2(248),
     +INTDS1(INTMD),INTDS2(INTMD),INTSD1(INTMD),INTSD2(INTMD)
 
C  /INTNEW/
C       NVV    :   NUMBER OF INTERACTING VALENCE-VALENCE SYSTEMS
C       NSV    :   NUMBER OF INTERACTING SEA-VALENCE SYSTEMS
C       NVS    :   NUMBER OF INTERACTING VALENCE-SEA SYSTEMS
C       NSS    :   NUMBER OF INTERACTING SEA-SEA SYSTEMS
C       IXPV, IXTV : NUMBER OF GENERATED X-VALUES FOR VALENCE-QUARK
C                     SYSTEMS FROM PROJECTILE/TARGET NUCLEI
C       IXPS, IXTS : NUMBER OF GENERATED X-VALUES FOR SEA-QUARK PAIRS
C                     FROM PROJECTILE/TARGET NUCLEI
C-------------------
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,LOZUO.
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/ ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),ZUOST(INTMX),
     +INTLO(INTMX),INLOSS(INTMX)
C  /LOZUO/
C       INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C------------------
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
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
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEEP,NNCMS.
      COMMON /NNCMS/  GAMCM,BGCM,UMO,PCM,EPROJ,PPROJ
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
*KEEP,NUCIMP.
      COMMON /NUCIMP/ PRMOM(5,248),TAMOM(5,248), PRMFEP,PRMFEN,TAMFEP,
     +TAMFEN, PREFEP,PREFEN,TAEFEP,TAEFEN, PREPOT(210),TAEPOT(210),
     +PREBIN,TAEBIN,FERMOD,ETACOU
*KEEP,PROJK.
      COMMON /PROJK/ IPROJK
*KEND.
       IF(IJPROJ.EQ.5)RETURN
C
C******************************    PROJECTILE
C                            - INTERACTING PROJECTILES ISTHKK=11
      DO 10 J=1,IP
C       IF(ISTHKK(J).EQ.11) THEN
          KK=KKPROJ(J)
          PRMOM(1,J)=PHKK(1,J)
          PRMOM(2,J)=PHKK(2,J)
          GAPROJ=EPROJ/AAM(KK)
          BGPROJ=PPROJ/AAM(KK)
          CALL SLTRAF(GAPROJ,-BGPROJ, PHKK(4,J),PHKK(3,J),PRMOM4,PRMOM3)
 
          CALL SLTRAF(GAMCM,+BGCM, PRMOM4,PRMOM3,PRMOM(4,J),PRMOM(3,J))
 
          PRMOM(5,J)=SQRT( ABS((PRMOM(4,J)-AAM(KK)) *(PRMOM(4,J)+AAM(KK)
     +    )))
C       ENDIF
   10 CONTINUE
C
C------------------------------    TARGET
C                             INTERACTING TARGET NUCLEONS ISTHKK=12
      IHKK=IP
      DO 20 J=1,IT
        IHKK=IHKK + 1
C       IF(ISTHKK(IHKK).EQ.12) THEN
          KK=KKTARG(J)
          TAMOM(1,J)=PHKK(1,IHKK)
          TAMOM(2,J)=PHKK(2,IHKK)
          CALL SLTRAF(GAMCM,BGCM, PHKK(4,IHKK),PHKK(3,IHKK),TAMOM(4,J),
     +    TAMOM(3,J))
          TAMOM(5,J)=SQRT(ABS( (TAMOM(4,J)-AAM(KK)) 
     +    *(TAMOM(4,J)+AAM(KK))))
 
C       ENDIF
   20 CONTINUE
C
      IF(IPEV.GE.6) THEN
        WRITE(6,'(/A,I5/5X,A)') ' NUCMOM: IP=',IP,
     +  ' J,IPVQ(J),IPPV1(J),IPPV2(J),ISTHKK,KKPROJ,PRMOM'
        DO 30 J=1,IP
          WRITE(6,'(I4,5I3,5(1PE11.3))') J,ISTHKK(J),KKPROJ(J), IPVQ(J),
     +    IPPV1(J),IPPV2(J), (PRMOM(JJ,J),JJ=1,5)
 
   30   CONTINUE
C
        WRITE(6,'(/A,I5/5X,A)') ' NUCMOM: IT=',IT,
     +  ' J,ITVQ(J),ITTV1(J),ITTV2(J),ISTHKK,KKTARG,TAMOM'
        IHKK=IP
        DO 40 J=1,IT
          IHKK=IHKK + 1
          WRITE(6,'(I4,5I3,5(1PE11.3))') J,ISTHKK(IHKK),KKTARG(J), ITVQ
     +    (J),ITTV1(J),ITTV2(J), (TAMOM(JJ,J),JJ=1,5)
 
   40   CONTINUE
      ENDIF
C
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE FER4M(PFERM,PXT,PYT,PZT,ET,KT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C                SAMPLE FERMI MOMENTUM FROM DISTRIBUTION WITH T=0
C-----------
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
*KEEP,DROPPT.
      LOGICAL INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA, IPADIS,
     +ISHMAL,LPAULI
      COMMON /DROPPT/ INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA,
     +IPADIS,ISHMAL,LPAULI
*KEND.
C-----------
      IF (FERMP) THEN
        CALL DFERMI(PABS)
        PABS=PFERM*PABS
C                                 SAMPLE ANGLES
        CALL DPOLI(POLC,POLS)
        CALL DSFECF(SFE,CFE)
C
        CXTA=POLS*CFE
        CYTA=POLS*SFE
        CZTA=POLC
        ET=SQRT(PABS*PABS+AAM(KT)**2)
        PXT=CXTA*PABS
        PYT=CYTA*PABS
        PZT=CZTA*PABS
C
      ELSE
        ET=AAM(KT)
        PXT=0.
        PYT=0.
        PZT=0.
      ENDIF
C
      RETURN
      END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE FER4MP(IP,PFERM,PXT,PYT,PZT,ET,KT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      COMMON /FERFOR/IFERFO
C
C                SAMPLE FERMI MOMENTUM FROM DISTRIBUTION WITH T=0
C-----------
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
*KEEP,DROPPT.
      LOGICAL INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA, IPADIS,
     +ISHMAL,LPAULI
      COMMON /DROPPT/ INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA,
     +IPADIS,ISHMAL,LPAULI
*KEND.
C-----------
      IF (FERMP) THEN
        IF(IFERFO.EQ.1)THEN
	  CALL DFERMI(PABS)
          PABS=PFERM*PABS
        ENDIF
	IF(IFERFO.EQ.2)CALL DFATPR(IP,PABS)
C                                 SAMPLE ANGLES
        CALL DPOLI(POLC,POLS)
        CALL DSFECF(SFE,CFE)
C
        CXTA=POLS*CFE
        CYTA=POLS*SFE
        CZTA=POLC
        ET=SQRT(PABS*PABS+AAM(KT)**2)
        PXT=CXTA*PABS
        PYT=CYTA*PABS
        PZT=CZTA*PABS
C
      ELSE
        ET=AAM(KT)
        PXT=0.
        PYT=0.
        PZT=0.
      ENDIF
C
      RETURN
      END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE FER4MT(IT,PFERM,PXT,PYT,PZT,ET,KT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      COMMON /FERFOR/IFERFO
C
C                SAMPLE FERMI MOMENTUM FROM DISTRIBUTION WITH T=0
C-----------
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
*KEEP,DROPPT.
      LOGICAL INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA, IPADIS,
     +ISHMAL,LPAULI
      COMMON /DROPPT/ INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA,
     +IPADIS,ISHMAL,LPAULI
*KEND.
C-----------
C     WRITE(6,*)' FERMP',FERMP
      IF (FERMP) THEN
        IF(IFERFO.EQ.1)THEN
	  CALL DFERMI(PABS)
CWRITE(6,*)' PABS',PABS
          PABS=PFERM*PABS
CWRITE(6,*)' PABS',PABS
        ENDIF
	IF(IFERFO.EQ.2)CALL DFATTA(IT,PABS)
C                                 SAMPLE ANGLES
CWRITE(6,*)' PABS',PABS
        CALL DPOLI(POLC,POLS)
        CALL DSFECF(SFE,CFE)
C
        CXTA=POLS*CFE
        CYTA=POLS*SFE
        CZTA=POLC
        ET=SQRT(PABS*PABS+AAM(KT)**2)
        PXT=CXTA*PABS
        PYT=CYTA*PABS
        PZT=CZTA*PABS
C
      ELSE
        ET=AAM(KT)
        PXT=0.
        PYT=0.
        PZT=0.
      ENDIF
C
      RETURN
      END
      SUBROUTINE DFATTA(IT,PABS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C       FERMI MOMENTUM A LA C. CIOFI DEGLI ATTI ET AL PRC53(96)1689
      DIMENSION PAR10(6),PAR20(6),PAR30(6),PAR40(6),PAR50(6),
     *          PAR60(6),PAR11(6),PAR21(6),PAR31(6),PAR41(6),
     *          AIA(6),ATT(101),CATT(101),AKA(101)
      COMMON/FATTAD/DAKA(101),FATT(101)
      DATA PAR10/1.61D0,2.74D0,3.24D0,3.57D0,1.80D0,0.D0/
      DATA PAR20/2.66D0,3.33D0,3.72D0,4.97D0,4.77D0,0.D0/
      DATA PAR30/3.54D0,6.66D0,0.D0,0.D0,0.D0,0.D0/
      DATA PAR40/0.D0,0.D0,11.1D0,19.8D0,25.5D0,0.D0/
      DATA PAR50/0.D0,0.D0,0.D0,15.D0,0.D0,0.D0/
      DATA PAR60/0.D0,0.D0,0.D0,0.D0,40.3D0,0.D0/
      DATA PAR11/.426D0,.326D0,.419D0,.230D0,.275D0,0.D0/
      DATA PAR21/1.6D0,1.4D0,1.77D0,1.2D0,1.01D0,0.D0/
      DATA PAR31/.0237D0,.0263D0,.0282D0,.0286D0,.0304D0,0.D0/
      DATA PAR41/.22D0,.22D0,.22D0,.22D0,.22D0,0.D0/
      DATA AIA/12.D0,16.D0,40.D0,56.D0,208.D0,209.D0/
      DATA INIT/0/
      AIT=IT
      IF(INIT.EQ.0)THEN
C                 INITIALIZATION
C                 INTERPOLATE PARAMETERS
      DO 1 I=1,4
	IF(AIT.GE.AIA(I).AND.AIT.LT.AIA(I+1))THEN
	  DAIT=(AIT-AIA(I))/(AIA(I+1)-AIA(I))
	  DBIT=1.D0-DAIT
	  III=I
        ENDIF
   1  CONTINUE
      IF(AIT.LT.AIA(1))THEN
	DBIT=1.D0
	DAIT=0.D0
	III=1
      ENDIF
      IF(AIT.GE.AIA(5))THEN
	DBIT=1.D0
	DAIT=0.D0
	III=5
      ENDIF
      A0=DBIT*PAR10(III)+DAIT*PAR10(III+1)
      B0=DBIT*PAR20(III)+DAIT*PAR20(III+1)
      C0=DBIT*PAR30(III)+DAIT*PAR30(III+1)
      D0=DBIT*PAR40(III)+DAIT*PAR40(III+1)
      E0=DBIT*PAR50(III)+DAIT*PAR50(III+1)
      F0=DBIT*PAR60(III)+DAIT*PAR60(III+1)
      A1=DBIT*PAR11(III)+DAIT*PAR11(III+1)
      B1=DBIT*PAR21(III)+DAIT*PAR21(III+1)
      C1=DBIT*PAR31(III)+DAIT*PAR31(III+1)
      D1=DBIT*PAR41(III)+DAIT*PAR41(III+1)
      INIT=1
      DK=0.04D0
      CATT(1)=0.D0
      DO 2 I=1,101
	AI=I
	AK=(AI-1.D0)*DK
	AKA(I)=AK
	DAKA(I)=AKA(I)
	ATT(I)=AK**2*(A0*EXP(-B0*AK**2)*(1.D0+C0*AK**2+
     *         D0*AK**4+E0*AK**6+F0*AK**8)+
     *         A1*EXP(-B1*AK**2)+C1*EXP(-D1*AK**2))
	IF(I.GT.1)CATT(I)=CATT(I-1)+ATT(I)
   2  CONTINUE
      DO 3 I=1,101
	CATT(I)=CATT(I)/CATT(101)
	FATT(I)=0.D0
   3  CONTINUE
      ENDIF
C          END INITIALIZATION
      RNDFA=RNDM(V)
      DO 10 I=1,101
	IF(RNDFA.LT.CATT(I))THEN
	  IATT=I
	  GO TO 11
	ENDIF
   10 CONTINUE
   11 CONTINUE
      PABS=AKA(IATT)*0.197D0
      FATT(IATT)=FATT(IATT)+1.D0/PABS**2
      RETURN
      END
      SUBROUTINE DFATPR(IP,PABS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C       FERMI MOMENTUM A LA C. CIOFI DEGLI ATTI ET AL PRC53(96)1689
      DIMENSION PAR10(6),PAR20(6),PAR30(6),PAR40(6),PAR50(6),
     *          PAR60(6),PAR11(6),PAR21(6),PAR31(6),PAR41(6),
     *          AIA(6),ATT(101),CATT(101),AKA(101)
      DATA PAR10/1.61D0,2.74D0,3.24D0,3.57D0,1.80D0,0.D0/
      DATA PAR20/2.66D0,3.33D0,3.72D0,4.97D0,4.77D0,0.D0/
      DATA PAR30/3.54D0,6.66D0,0.D0,0.D0,0.D0,0.D0/
      DATA PAR40/0.D0,0.D0,11.1D0,19.8D0,25.5D0,0.D0/
      DATA PAR50/0.D0,0.D0,0.D0,15.D0,0.D0,0.D0/
      DATA PAR60/0.D0,0.D0,0.D0,0.D0,40.3D0,0.D0/
      DATA PAR11/.426D0,.326D0,.419D0,.230D0,.275D0,0.D0/
      DATA PAR21/1.6D0,1.4D0,1.77D0,1.2D0,1.01D0,0.D0/
      DATA PAR31/.0237D0,.0263D0,.0282D0,.0286D0,.0304D0,0.D0/
      DATA PAR41/.22D0,.22D0,.22D0,.22D0,.22D0,0.D0/
      DATA AIA/12.D0,16.D0,40.D0,56.D0,208.D0,209.D0/
      DATA INIT/0/
      AIT=IP
      IF(INIT.EQ.0)THEN
C                 INITIALIZATION
C                 INTERPOLATE PARAMETERS
      DO 1 I=1,4
	IF(AIT.GE.AIA(I).AND.AIT.LT.AIA(I+1))THEN
	  DAIT=(AIT-AIA(I))/(AIA(I+1)-AIA(I))
	  DBIT=1.D0-DAIT
	  III=I
        ENDIF
   1  CONTINUE
      IF(AIT.LT.AIA(1))THEN
	DBIT=1.D0
	DAIT=0.D0
	III=1
      ENDIF
      IF(AIT.GE.AIA(5))THEN
	DBIT=1.D0
	DAIT=0.D0
	III=5
      ENDIF
      A0=DBIT*PAR10(III)+DAIT*PAR10(III+1)
      B0=DBIT*PAR20(III)+DAIT*PAR20(III+1)
      C0=DBIT*PAR30(III)+DAIT*PAR30(III+1)
      D0=DBIT*PAR40(III)+DAIT*PAR40(III+1)
      E0=DBIT*PAR50(III)+DAIT*PAR50(III+1)
      F0=DBIT*PAR60(III)+DAIT*PAR60(III+1)
      A1=DBIT*PAR11(III)+DAIT*PAR11(III+1)
      B1=DBIT*PAR21(III)+DAIT*PAR21(III+1)
      C1=DBIT*PAR31(III)+DAIT*PAR31(III+1)
      D1=DBIT*PAR41(III)+DAIT*PAR41(III+1)
      INIT=1
      DK=0.04D0
      CATT(1)=0.D0
      DO 2 I=1,101
	AI=I
	AK=(AI-1.D0)*DK
	AKA(I)=AK
	ATT(I)=AK**2*(A0*EXP(-B0*AK**2)*(1.D0+C0*AK**2+
     *         D0*AK**4+E0*AK**6+F0*AK**8)+
     *         A1*EXP(-B1*AK**2)+C1*EXP(-D1*AK**2))
	IF(I.GT.1)CATT(I)=CATT(I-1)+ATT(I)
   2  CONTINUE
      DO 3 I=1,101
	CATT(I)=CATT(I)/CATT(101)
   3  CONTINUE
      ENDIF
C          END INITIALIZATION
      RNDFA=RNDM(V)
      DO 10 I=1,101
	IF(RNDFA.LT.CATT(I))THEN
	  IATT=I
	  GO TO 11
	ENDIF
   10 CONTINUE
   11 CONTINUE
      PABS=AKA(IATT)*0.197D0
      RETURN
      END
*-- Author :
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE FLKSAM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C                              QUARK CONTENT
C                              OF PROJECTILE AND TARGET
C                              (HADRONS / ALL NUCLEONS)
C---------------------------------------------------------------------
*KEEP,INTMX.
      PARAMETER (INTMX=2488,INTMD=252)
*KEEP,DXQX.
C     INCLUDE (XQXQ)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON /DXQX/ XPVQ(248),XPVD(248),XTVQ(248),XTVD(248), XPSQ
     +(INTMX),XPSAQ(INTMX),XTSQ(INTMX),XTSAQ(INTMX)
     *                ,XPSU(248),XTSU(248)
     *                ,XPSUT(248),XTSUT(248)
*KEEP,INTNEW.
      COMMON /INTNEW/ NVV,NSV,NVS,NSS,NDV,NVD,NDS,NSD,
     +IXPV,IXPS,IXTV,IXTS, INTVV1(248),
     +INTVV2(248),INTSV1(248),INTSV2(248), INTVS1(248),INTVS2(248),
     +INTSS1(INTMX),INTSS2(INTMX),
     +INTDV1(248),INTDV2(248),INTVD1(248),INTVD2(248),
     +INTDS1(INTMD),INTDS2(INTMD),INTSD1(INTMD),INTSD2(INTMD)
 
C  /INTNEW/
C       NVV    :   NUMBER OF INTERACTING VALENCE-VALENCE SYSTEMS
C       NSV    :   NUMBER OF INTERACTING SEA-VALENCE SYSTEMS
C       NVS    :   NUMBER OF INTERACTING VALENCE-SEA SYSTEMS
C       NSS    :   NUMBER OF INTERACTING SEA-SEA SYSTEMS
C       IXPV, IXTV : NUMBER OF GENERATED X-VALUES FOR VALENCE-QUARK
C                     SYSTEMS FROM PROJECTILE/TARGET NUCLEI
C       IXPS, IXTS : NUMBER OF GENERATED X-VALUES FOR SEA-QUARK PAIRS
C                     FROM PROJECTILE/TARGET NUCLEI
C-------------------
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,LOZUO.
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/ ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),ZUOST(INTMX),
     +INTLO(INTMX),INLOSS(INTMX)
C  /LOZUO/
C       INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C------------------
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
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
*KEEP,PROJK.
      COMMON /PROJK/ IPROJK
*KEND.
C----------
      DIMENSION IHKKQ(-6:6),IHKKQQ(-3:3,-3:3)
      DATA IHKKQ/-6,-5,-4,-3,-1,-2,0,2,1,3,4,5,6/
      DATA IHKKQQ/-3303,-3103,-3203,0,0,0,0, -3103,-1103,-2103,0,0,0,0,
     +-3203,-2103,-2203,0,0,0,0, 0,0,0,0,0,0,0, 0,0,0,0,2203,2103,3203,
     +0,0,0,0,2103,1103,3103, 0,0,0,0,3203,3103,3303/
C----------------------------------------------------------------------
C
C   FLAVORS OF VALENCE QUARKS FROM PROJECTILE HADRON/NUCLEONS
C
C-----
C
      IF(IPEV.GE.3) WRITE(6,'(A,6I4)')
     +' FLKSAM-ENTRY: IT,ITZ, IP,IPZ, IJPROJ,IBPROJ', IT,ITZ,IP,IPZ,
     +IJPROJ,IBPROJ
C
      IXPSS=IXPS
      IXTSS=IXTS
      IXPVV=IXPV
      IXTVV=IXTV
      DO 10 JP=1,IXPVV
        IFR=IFROVP(JP)
        KPROJ=KKPROJ(IFR)
        CALL FLAHAD(KPROJ,IBPROJ,IPVQ(JP),IPPV1(JP),IPPV2(JP))
C
        IF (IPEV.GE.6) WRITE (6,1000)IPVQ(JP),IPPV1(JP),IPPV2(JP)
 1000 FORMAT (' FLKSAM: IPVQ,IPPV1,IPPV2 ',3I4)
C
        JHKK=JHKKPV(JP)
        JHKKQ=JHKK - 1
        IDHKK(JHKKQ)=IHKKQ(IPVQ(JP))
        IF(IBPROJ.EQ.0) THEN
          IDHKK(JHKK)=IHKKQ(IPPV1(JP))
        ELSE
          IDHKK(JHKK)=IHKKQQ(IPPV1(JP),IPPV2(JP))
        ENDIF
C
   10 CONTINUE
C
C*********************************************************************
C
C-------------------------------SAMPLING PROJECTILE SEA FLAVORS-------
C
C*********************************************************************
C
      DO 20 N=1,IXPSS
C
        JHKKAQ=JHKKPS(N)
        JHKKQ=JHKKAQ - 1
        IDHKK(JHKKQ)=IHKKQ(IPSQ(N))
        IDHKK(JHKKAQ)=IHKKQ(IPSAQ(N))
C
   20 CONTINUE
C--------------------------------------------------------------------
C
C   FLAVORS OF VALENCE QUARKS FROM TARGET HADRON / ALL NUCLEONS
C
C-----
C
      DO 30 JT=1,IXTVV
        IFR=IFROVT(JT)
        KTARG=KKTARG(IFR)
        IBTARG=IIBAR(KTARG)
        CALL FLAHAD(KTARG,IBTARG,ITVQ(JT),ITTV1(JT),ITTV2(JT))
C
        JHKK=JHKKTV(JT)
        JHKKQ=JHKK - 1
        IDHKK(JHKKQ)=IHKKQ(ITVQ(JT))
        IF(IBTARG.EQ.0) THEN
          IDHKK(JHKK)=IHKKQ(ITTV1(JT))
        ELSE
          IDHKK(JHKK)=IHKKQQ(ITTV1(JT),ITTV2(JT))
        ENDIF
        IF (IPEV.GE.8) WRITE (6,'(A,8I4)')
     +  ' FLKSAM: KTARG,ITVQ(JT),ITTV1(JT),ITTV2(JT)', KTARG,ITVQ(JT),
     +  ITTV1(JT),ITTV2(JT), IDHKK(JHKKQ),JHKKQ,IDHKK(JHKK),JHKK
 
C
C
   30 CONTINUE
C
C*********************************************************************
C
C-------------------------------SAMPLING TARGET SEA FLAVORS-------
C
C*********************************************************************
C
      DO 40 N=1,IXTSS
C
        JHKKAQ=JHKKTS(N)
        JHKKQ=JHKKAQ - 1
        IDHKK(JHKKQ)=IHKKQ(ITSQ(N))
        IDHKK(JHKKAQ)=IHKKQ(ITSAQ(N))
C
   40 CONTINUE
C
      RETURN
      END
*-- Author :
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE FLKSAA(NN,ECM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C                              QUARK CONTENT
C                              OF PROJECTILE AND TARGET
C                              (HADRONS / ALL NUCLEONS)
C                              first run sea quark flavors
C---------------------------------------------------------------------
*KEEP,INTMX.
      PARAMETER (INTMX=2488,INTMD=252)
*KEEP,DXQX.
C     INCLUDE (XQXQ)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON /DXQX/ XPVQ(248),XPVD(248),XTVQ(248),XTVD(248), XPSQ
     +(INTMX),XPSAQ(INTMX),XTSQ(INTMX),XTSAQ(INTMX)
     *                ,XPSU(248),XTSU(248)
     *                ,XPSUT(248),XTSUT(248)
*KEEP,INTNEW.
      COMMON /INTNEW/ NVV,NSV,NVS,NSS,NDV,NVD,NDS,NSD,
     +IXPV,IXPS,IXTV,IXTS, INTVV1(248),
     +INTVV2(248),INTSV1(248),INTSV2(248), INTVS1(248),INTVS2(248),
     +INTSS1(INTMX),INTSS2(INTMX),
     +INTDV1(248),INTDV2(248),INTVD1(248),INTVD2(248),
     +INTDS1(INTMD),INTDS2(INTMD),INTSD1(INTMD),INTSD2(INTMD)
 
C  /INTNEW/
C       NVV    :   NUMBER OF INTERACTING VALENCE-VALENCE SYSTEMS
C       NSV    :   NUMBER OF INTERACTING SEA-VALENCE SYSTEMS
C       NVS    :   NUMBER OF INTERACTING VALENCE-SEA SYSTEMS
C       NSS    :   NUMBER OF INTERACTING SEA-SEA SYSTEMS
C       IXPV, IXTV : NUMBER OF GENERATED X-VALUES FOR VALENCE-QUARK
C                     SYSTEMS FROM PROJECTILE/TARGET NUCLEI
C       IXPS, IXTS : NUMBER OF GENERATED X-VALUES FOR SEA-QUARK PAIRS
C                     FROM PROJECTILE/TARGET NUCLEI
C-------------------
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,LOZUO.
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/ ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),ZUOST(INTMX),
     +INTLO(INTMX),INLOSS(INTMX)
C  /LOZUO/
C       INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C------------------
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
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
*KEEP,PROJK.
      COMMON /PROJK/ IPROJK
      COMMON /SEASU3/SEASQ 
C     COMMON /PCHARM/PCCCC
      PARAMETER (UMMM=0.3D0)
      PARAMETER (SMMM=0.5D0)
      PARAMETER (CMMM=1.3D0)
      DATA PC/0.0001D0/
*KEND.
C----------
C     DIMENSION IHKKQ(-6:6),IHKKQQ(-3:3,-3:3)
C     DATA IHKKQ/-6,-5,-4,-3,-1,-2,0,2,1,3,4,5,6/
C     DATA IHKKQQ/-3303,-3103,-3203,0,0,0,0, -3103,-1103,-2103,0,0,0,0,
C    +-3203,-2103,-2203,0,0,0,0, 0,0,0,0,0,0,0, 0,0,0,0,2203,2103,3203,
C    +0,0,0,0,2103,1103,3103, 0,0,0,0,3203,3103,3303/
C----------------------------------------------------------------------
C
C   FLAVORS OF VALENCE QUARKS FROM PROJECTILE HADRON/NUCLEONS
C
C-----
C
      DATA INICHA/0/
C----------------------------------------------------------------------
C                     Initialize Charm selection at soft chain ends
C
      IF(INICHA.EQ.0)THEN
        RX=8.
        X1=RX
        GM=2.140
        X2=UMMM
	BETOO=7.5D0
      ENDIF
      RX=8.
      X1=RX
      BETCHA=BETOO+1.3-LOG10(ECM)
      PU=DBETA(X1,X2,BETCHA)
      X2=SMMM
      PS=DBETA(X1,X2,BETCHA)
      X2=CMMM
      PC=DBETA(X1,X2,BETCHA)
C     PU1=PU/(2*PU+PS+PC)
C     PS1=PS/(2*PU+PS+PC)
      PC1=PC/(2*PU+PS+PC)
C                       changed j.r.7.12.94
C     PC=PC1/2.9
C                       changed j.r.14.12.94
C     PC=PC1/5.0
C     PC=PC1/10.0
      PC=PC1/7.0
        PU1=PU/(2*PU+PS+PC)
        PS1=PS/(2*PU+PS+PC)
      IF(INICHA.EQ.0)THEN
        INICHA=1
        WRITE(6,4567)PC,BETCHA,PU1,PS1
 4567   FORMAT(' Charm at chain ends FLKSAA: PC,BETCHA,PU,PS ',4F10.5)
      ENDIF
C----------------------------------------------------------------------
C
      IF(IPEV.GE.3) WRITE(6,'(A,6I4)')
     +' FLKSAA-ENTRY: IT,ITZ, IP,IPZ, IJPROJ,IBPROJ', IT,ITZ,IP,IPZ,
     +IJPROJ,IBPROJ
C
      IXPSS=NN
      IXTSS=NN
C
C*********************************************************************
C
C-------------------------------SAMPLING PROJECTILE SEA FLAVORS-------
C
C*********************************************************************
C
      DO 20 N=1,IXPSS
        IS=1
        RR=RNDM(V)
        IS=1.D0+RNDM(V1)*(2.D0+SEASQ)
	IF(RR.LT.PC)IS=4
        IPSQ(N)=IS
        IPSAQ(N)=-IS
        IF (IPEV.GE.8) WRITE (6,1010) N,IPSQ(N),IPSAQ(N)
 1010 FORMAT (' FLKSAA: N,IPSQ(N),IPSAQ(N) ',3I4)
C
   20 CONTINUE
C
C*********************************************************************
C
C-------------------------------SAMPLING TARGET SEA FLAVORS-------
C
C*********************************************************************
C
      DO 40 N=1,IXTSS
        IS=1
        RR=RNDM(V)
        IS=1.D0+RNDM(V1)*(2.D0+SEASQ) 
	IF(RR.LT.PC)IS=4
        ITSQ(N)=IS
        ITSAQ(N)=-IS
        IF (IPEV.GE.8) WRITE (6,1020) N,ITSQ(N),ITSAQ(N)
 1020 FORMAT (' FLKSAA: N,ITSQ(N),ITSAQ(N) ',3I4)
C
   40 CONTINUE
C
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE FLAHAD(ITYP,IBAR,IF1,IF2,IF3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C                               QUARK FLAVOR COMPOSITION FOR HADRONS
C                               ITYP : NUMBERING AS FOR BAMJET ...
C                                      LE.30 !!!!!!!!!!
C
C----------------------------------------------------------------------
C
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEND.
      DIMENSION MQUARK(3,30)
      DATA MQUARK/ 2,1,1, -2,-1,-1, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
     +2,2,1, -2,-2,-1, 0,0,0, 0,0,0, 0,0,0, 1,-2,0, 2,-1,0, 1,-3,0, 3,
     +-1,0, 1,2,3, -1,-2,-3, 0,0,0, 2,2,3, 1,1,3, 1,2,3, 1,-1,0, 2,-3,0,
     +3,-2,0, 2,-2,0, 3,-3,0, 0,0,0, 0,0,0, 0,0,0/
C----------------------------------------------------------------------
      IF(IBAR.NE.0) THEN
        IPQ1 = MQUARK(1,ITYP)
        IPQ2 = MQUARK(2,ITYP)
        IPQ3 = MQUARK(3,ITYP)
C
        IF(IPEV.GE.3) PRINT 1000, ITYP,IBAR
 1000 FORMAT(' FLAHAD: ITYP,IBAR',2I5)
C
        ISAM5=1. + 6.*RNDM(V)
        GO TO (10,20,30,40,50,60),ISAM5
   10   CONTINUE
        IF1=IPQ1
        IF2=IPQ2
        IF3=IPQ3
                                                                GO TO 70
   20   CONTINUE
        IF1=IPQ2
        IF2=IPQ3
        IF3=IPQ1
                                                                GO TO 70
   30   CONTINUE
        IF1=IPQ3
        IF2=IPQ1
        IF3=IPQ2
                                                                GO TO 70
   40   CONTINUE
        IF1=IPQ1
        IF2=IPQ3
        IF3=IPQ2
                                                                GO TO 70
   50   CONTINUE
        IF1=IPQ2
        IF2=IPQ1
        IF3=IPQ3
                                                                GO TO 70
   60   CONTINUE
        IF1=IPQ3
        IF2=IPQ2
        IF3=IPQ1
   70   CONTINUE
        IF (IPEV.GE.3) WRITE (6,1010) IF1,IF2,IF3
 1010 FORMAT (' FLAHAD: IF1,IF2,IF3 ',3I4)
      ELSE
C                         VALENCE QUARK FLAVORS FOR MESONS
        IF1=MQUARK(1,ITYP)
        IF2=MQUARK(2,ITYP)
        IF3=0
        IF(IPEV.GE.3) THEN
          WRITE(6,'(A,6I4)') ' FLAHAD (MESON): IF1,IF2,IF3', IF1,IF2,IF3
 
 
        ENDIF
      ENDIF
C
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE XKSAMP(NN,ECM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
*-----------------------------------------------------------
* SAMPLING MOMENTUM FRACTIONS OF QUARKS AND DIQUARKS
*-----------------------------------------------------------
*KEEP,DINPDA.
      COMMON /DINPDA/ IMPS(6,6),IMVE(6,6),IB08(6,21),IB10(6,21), IA08
     +(6,21),IA10(6,21), A1,B1,B2,B3,LT,LE,BET,AS,B8,AME,DIQ,ISU
*KEEP,INTMX.
      PARAMETER (INTMX=2488,INTMD=252)
      PARAMETER (AMIS=0.8D0,AMAS=2.6D0,AMIU=0.5D0,AMAU=2.6D0)
*KEEP,DXQX.
C     INCLUDE (XQXQ)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON /DXQX/ XPVQ(248),XPVD(248),XTVQ(248),XTVD(248), XPSQ
     +(INTMX),XPSAQ(INTMX),XTSQ(INTMX),XTSAQ(INTMX)
     *                ,XPSU(248),XTSU(248)
     *                ,XPSUT(248),XTSUT(248)
*KEEP,INTNEW.
      COMMON /INTNEZ/NDZ,NZD
      COMMON /INTNEW/ NVV,NSV,NVS,NSS,NDV,NVD,NDS,NSD,
     +IXPV,IXPS,IXTV,IXTS, INTVV1(248),
     +INTVV2(248),INTSV1(248),INTSV2(248), INTVS1(248),INTVS2(248),
     +INTSS1(INTMX),INTSS2(INTMX),
     +INTDV1(248),INTDV2(248),INTVD1(248),INTVD2(248),
     +INTDS1(INTMD),INTDS2(INTMD),INTSD1(INTMD),INTSD2(INTMD)
 
C  /INTNEW/
C       NVV    :   NUMBER OF INTERACTING VALENCE-VALENCE SYSTEMS
C       NSV    :   NUMBER OF INTERACTING SEA-VALENCE SYSTEMS
C       NVS    :   NUMBER OF INTERACTING VALENCE-SEA SYSTEMS
C       NSS    :   NUMBER OF INTERACTING SEA-SEA SYSTEMS
C       IXPV, IXTV : NUMBER OF GENERATED X-VALUES FOR VALENCE-QUARK
C                     SYSTEMS FROM PROJECTILE/TARGET NUCLEI
C       IXPS, IXTS : NUMBER OF GENERATED X-VALUES FOR SEA-QUARK PAIRS
C                     FROM PROJECTILE/TARGET NUCLEI
C-------------------
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,LOZUO.
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/ ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),ZUOST(INTMX),
     +INTLO(INTMX),INLOSS(INTMX)
C  /LOZUO/
C       INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C------------------
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
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
*KEEP,SHMAKL.
C     INCLUDE (SHMAKL)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON/SHMAKL/JSSH(INTMX),JTSH(INTMX),INTER1(INTMX),INTER2(INTMX)
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
*KEEP,NUCKOO.
      COMMON /NUCKOO/ PKOO(3,INTMX),TKOO(3,INTMX),PPOO(3,INTMX),
     +TPOO(3,INTMX)
*KEEP,PROJK.
      COMMON /PROJK/ IPROJK
*KEEP,XSEADI.
      COMMON /XSEADI/ XSEACU,UNON,UNOM,UNOSEA, CVQ,CDQ,CSEA,SSMIMA,
     +SSMIMQ,VVMTHR
*
*KEEP,ABRVV.
      COMMON /ABRVV/ AMCVV1(248),AMCVV2(248),GACVV1(248),GACVV2(248),
     +BGXVV1(248),BGYVV1(248),BGZVV1(248), BGXVV2(248),BGYVV2(248),
     +BGZVV2(248), NCHVV1(248),NCHVV2(248),IJCVV1(248),IJCVV2(248),
     +PQVVA1(248,4),PQVVA2(248,4), PQVVB1(248,4),PQVVB2(248,4)
*KEEP,ABRSS.
C     INCLUDE (ABRSS)
      COMMON /ABRSS/ AMCSS1(INTMX),AMCSS2(INTMX), GACSS1(INTMX),GACSS2
     +(INTMX), BGXSS1(INTMX),BGYSS1(INTMX),BGZSS1(INTMX), BGXSS2(INTMX),
     +BGYSS2(INTMX),BGZSS2(INTMX), NCHSS1(INTMX),NCHSS2(INTMX), IJCSS1
     +(INTMX),IJCSS2(INTMX), PQSSA1(INTMX,4),PQSSA2(INTMX,4), PQSSB1
     +(INTMX,4),PQSSB2(INTMX,4)
*KEEP,ABRVS.
      COMMON /ABRVS/ AMCVS1(248),AMCVS2(248),GACVS1(248),GACVS2(248),
     +BGXVS1(248),BGYVS1(248),BGZVS1(248), BGXVS2(248),BGYVS2(248),
     +BGZVS2(248), NCHVS1(248),NCHVS2(248),IJCVS1(248),IJCVS2(248),
     +PQVSA1(248,4),PQVSA2(248,4), PQVSB1(248,4),PQVSB2(248,4)
*KEEP,ABRSV.
      COMMON /ABRSV/ AMCSV1(248),AMCSV2(248),GACSV1(248),GACSV2(248),
     +BGXSV1(248),BGYSV1(248),BGZSV1(248), BGXSV2(248),BGYSV2(248),
     +BGZSV2(248), NCHSV1(248),NCHSV2(248),IJCSV1(248),IJCSV2(248),
     +PQSVA1(248,4),PQSVA2(248,4), PQSVB1(248,4),PQSVB2(248,4)
*KEND.
      LOGICAL LSEADI
      COMMON /SEADIQ/LSEADI
      COMMON/RECOM/IRECOM
      COMMON/DIQUAX/AMEDD,IDIQUA,IDIQUU
      COMMON /XSVTHR/ XSTHR,XVTHR,XDTHR,XSSTHR
      COMMON /SEAQXX/ SEAQX,SEAQXN 
      DIMENSION ISXPVQ(248),ISXPVD(248),ISXTVQ(248),ISXTVD(248)
      PARAMETER (SQMA=0.7D0)
C*******************************************************************"
C***     ACTUAL STANDARD VALUES FROM BLOCK DATA:
C
C      CSEA=1.0,  CVQ=1.,  CDQ=2.
C      UNON=2.,   UNOM=1.5,  UNOSEA=2.0
C----------------------------------
      PARAMETER (NSEA=3,NVAL=10)
      DATA ICOUN /0/
      DATA JCOUN /0/
*  NSEA: maximum number of trials to generate x's for the required number
*        of sea quark pairs for a given hadron
*    changed from 10 to 3      22/04/92
C---------------------------------------------------------------------
      JCOUN=JCOUN+1
      DO 10 I=1,IP
        JSSHS(I)=0
   10 CONTINUE
      DO 20 I=1,IT
        JTSHS(I)=0
   20 CONTINUE
      DO 30 I=1,INTMX
        ZUOSP(I)=.FALSE.
        ZUOST(I)=.FALSE.
        IF (I.GT.248)                                           GO TO 30
        ZUOVP(I)=.FALSE.
        ZUOVT(I)=.FALSE.
   30 CONTINUE
      IF(ECM.LE.1.D-3)THEN
	WRITE(6,*)' xksamp: ECM=0.D0 '
	ECM=ECM+1.D-3
      ENDIF
      XSTHR=CSEA/ECM
      IF(XSTHR.LE.1.D-12)THEN
	WRITE(6,*)' xksamp 30 : XSTHR=0.D0 ',CSEA,ECM,XSTHR
	XSTHR=XSTHR+1.D-12
      ENDIF
C-----------------------------------------------------------------
C
C                                  J.R.21.2.94
C
C----------------------------------------------------------------
C                                      j.r.12.3.97
C                                      j.r.11.4.97 part restored
      IF(IP.EQ.1) XSTHR=4./ECM**2
C                               test 28.4.97
C     IF(IP.EQ.1) XSTHR=4./ECM**2
      IF(XSTHR.LE.1.D-12)THEN
	XSTHR=XSTHR+1.D-12
      ENDIF
C----------------------------------------------------------------
C-----------------------------------------------------------------
C
C                                  J.R.16.3.95
C
C----------------------------------------------------------------
      IF(IP.GE.150.AND.IT.GE.150) XSTHR=2.5/(ECM*SQRT(ECM))
C----------------------------------------------------------------
      BSQMA=SQMA/ECM
C                       before 28.8.97
C     IF (ECM.LT.10.D0) XSTHR=((12.-ECM)/5.+1.)*CSEA/ECM
C                       28.4.97 test
      IF (ECM.LT.10.D0.AND.IP.GT.1)XSTHR=((12.-ECM)/5.+1.)*CSEA/ECM
      XVTHR=CVQ/ECM
      XDTHR=CDQ/ECM
      IF (XVTHR+XDTHR.GT.0.90D0)THEN
	XVTHR=0.95-XDTHR
	IF(XVTHR.LE.0.05D0)THEN
          WRITE (6,1000)ECM
	ENDIF
      ENDIF
      IF(ECM.LE.1.D-3)THEN
	WRITE(6,*)' xksamp: ECM=0.D0 '
	ECM=ECM+1.D-3
      ENDIF
      XSSTHR=SSMIMA/ECM
C-------------------------          20.12.91.j.r.
      IF(JCOUN.EQ.1)WRITE(6,'(A,4E15.5)')
     *' XKSAMP: XSTHR,XVTHR,XDTHR,XSSTHR ',
     * XSTHR,XVTHR,XDTHR,XSSTHR 
      IF(IPEV.GE.1)WRITE(6,'(A,4E15.5)')
     *' XKSAMP: XSTHR,XVTHR,XDTHR,XSSTHR ',
     * XSTHR,XVTHR,XDTHR,XSSTHR 
C-------------------------
C                                   TEST KINEMATICAL LIMITS
      IF (XVTHR+XDTHR.GT.0.95D0)THEN
        WRITE (6,1000)ECM
 1000   FORMAT (' PROGRAMM STOPPED IN XSAMP1 ECM = ',F6.2,' TOO SMALL')
        STOP
      ENDIF
C                  MAXIMUM NUMBER OF SEA-PAIRS ALLOWED KINEMATICALLY
C     XXSEAM=1.0 - XVTHR*(1.D0+RNDM(V1)) - XDTHR*(1.D0+RNDM(V2))
C    *             -0.01*(1.D0+5.D0*RNDM(V3))
C                                      28.4.97 test
      XXSEAM=1.0 - XVTHR*(1.D0+0.3D0*RNDM(V1)) 
     *            - XDTHR*(1.D0+0.3D0*RNDM(V2))
     *             -0.01*(1.D0+1.5D0*RNDM(V3))
C..............................................................  
C                           1/x seaquarks
      IF(SEAQXN.GE.0.75D0)THEN
	XSTHR=8.*CSEA/ECM
C                             23.5.95
	XSTHR=4.*CSEA/ECM
	XXSEAM=1.D0-XVTHR-XDTHR
C                  MAXIMUM NUMBER OF SEA-PAIRS ALLOWED KINEMATICALLY
        XXSEAM=1.0 - XVTHR*(1.D0+RNDM(V1)) - XDTHR*(1.D0+RNDM(V2))
     *             -0.01*(1.D0+5.D0*RNDM(V3))
      ENDIF
C..............................................................      
      IF(XSTHR.LE.1.D-9)THEN
	ICOUN=ICOUN+1
	IF(ICOUN.LE.50)THEN
          WRITE(6,*)' xksamp: XSTHR=0.D0 '
          WRITE(6,'(A,2E20.5,I10)')
     *            ' XXSEAM,XSTHR,NSMAX',XXSEAM,XSTHR,NSMAX
	ENDIF
	XSTHR=XSTHR+1.D-9
      ENDIF
      NSMAX=0.50*XXSEAM / XSTHR
      IF(IPEV.GE.1)WRITE(6,'(A,E15.5,I10)')
     *            ' XXSEAM,NSMAX',XXSEAM,NSMAX
*--------------------------------------------------------------------
*-------------------------------------------------------------------
C                                    Change XVTHR and XDTHR at low energies
C                                    TEST j.r. 9.2.95
      IF (XDTHR.GT.0.14D0)XDTHR=0.14D0
      IF (XVTHR.GT.0.14D0)XVTHR=0.14D0
*--------------------------------------------------------------------
C                                 PARTON X-VALUES OF INTERACTING
C                                         PROJECTILE HADRON / NUCLEONS
      IXPV=0
      IXPS=0
      UNOPRV=UNON
      IF(IJPROJ.NE.0.AND.IBPROJ.EQ.0) UNOPRV=UNOM
      IF(JCOUN.EQ.1)WRITE(6,'(A,4E15.5)')
     *' XKSAMP: XSTHR,XVTHR,XDTHR,XSSTHR ',
     * XSTHR,XVTHR,XDTHR,XSSTHR 
*  loop over projectile nucleons
      DO 100 IPP=1,IP
        IF (JSSH(IPP).NE.0) THEN
C--------------------------------------------------------------
C                          prepare diquark rejection
C--------------------------------------------------------------
          IIXPSS=IXPS
	  IIXPVV=IXPV
   99     CONTINUE
	  IXPS=IIXPSS
	  IXPV=IIXPVV
C--------------------------------------------------------------
          JIPP=JSSH(IPP)-1
          JIPP=MIN(JIPP,NSMAX)
 41       CONTINUE
          XXSEA=0.0
          IF(JIPP.GT.0) THEN
C                         j.r.11.12.97
	    XSMAX=XXSEAM - 1.5*JIPP*XSTHR
C           XSMAX=XXSEAM - 2.*JIPP*XSTHR
            IF(XSTHR.GE.XSMAX) THEN
              JIPP=JIPP-1
              GOTO 41
            ENDIF
*  x-values of sea-quark pairs
            NSCOUN=0
   40       CONTINUE
            IF(IPEV.GE.1)WRITE(6,'(A)') ' XKSAMP-40'
            XXSEA=0.0
            NSCOUN=NSCOUN+1
            IF (NSCOUN.GT.NSEA) THEN
              JIPP=JIPP-1
              NSCOUN=0
            ENDIF
            DO 70 ISQ=1,JIPP
C                                             j.r.29.4.93---
              IF(IPSQ(IXPS+1).LE.2)THEN
C..............................................................  
C                           1/sqrt(x) seaquarks
                IF(SEAQXN.LE.0.75D0)THEN
                  XPSQI=SAMPEX(XSTHR,XSMAX)
C                           1/x seaquarks
                ELSEIF(SEAQXN.GT.0.75D0)THEN
                  XPSQI=SAMPEY(XSTHR,XSMAX)
                ENDIF
C..............................................................      
		IF(IPEV.GE.1)WRITE(6,'(A,3E15.5)')
     *          'XPSQI 1:XPSQI,XSTHR,XSMAX',
     *          XPSQI,XSTHR,XSMAX
              ELSE
                IF(XSMAX.GT.XSTHR+BSQMA)THEN
                  XPSQI=SAMPXB(XSTHR+BSQMA,XSMAX,BSQMA)
	       	  IF(IPEV.GE.1)WRITE(6,'(A,4E15.5)')
     *            'XPSQI 2:XPSQI,XSTHR,XSMAX,BSQMA',
     *            XPSQI,XSTHR,XSMAX,BSQMA
                ELSE
C..............................................................  
C                           1/sqrt(x) seaquarks
                  IF(SEAQXN.LE.0.75D0)THEN
                    XPSQI=SAMPEX(XSTHR,XSMAX)
C                           1/x seaquarks
                  ELSEIF(SEAQXN.GT.0.75D0)THEN
                    XPSQI=SAMPEY(XSTHR,XSMAX)
                  ENDIF
C..............................................................      
	          IF(IPEV.GE.1)WRITE(6,'(A,3E15.5)')
     *            'XPSQI 3:XPSQI,XSTHR,XSMAX',
     *             XPSQI,XSTHR,XSMAX
                ENDIF
              ENDIF
C
              IF(IPSAQ(IXPS+1).GE.-2)THEN
C..............................................................  
C                           1/sqrt(x) seaquarks
                IF(SEAQXN.LE.0.75D0)THEN
                  XPSAQI=SAMPEX(XSTHR,XSMAX)
C                           1/x seaquarks
                ELSEIF(SEAQXN.GT.0.75D0)THEN
                  XPSAQI=SAMPEY(XSTHR,XSMAX)
                ENDIF
C..............................................................      
		IF(IPEV.GE.1)WRITE(6,'(A,3E15.5)')
     *          'XPSAQI 1:XPSAQI,XSTHR,XSMAX',
     *          XPSAQI,XSTHR,XSMAX
              ELSE
                IF(XSMAX.GT.XSTHR+BSQMA)THEN
                  XPSAQI=SAMPXB(XSTHR+BSQMA,XSMAX,BSQMA)
		  IF(IPEV.GE.1)WRITE(6,'(A,4E15.5)')
     *            'XPSAQI 2:XPSAQI,XSTHR,XSMAX,BSQMA',
     *            XPSAQI,XSTHR,XSMAX,BSQMA
                ELSE
C..............................................................  
C                           1/sqrt(x) seaquarks
                  IF(SEAQXN.LE.0.75D0)THEN
                    XPSAQI=SAMPEX(XSTHR,XSMAX)
C                           1/x seaquarks
                  ELSEIF(SEAQXN.GT.0.75D0)THEN
                    XPSAQI=SAMPEY(XSTHR,XSMAX)
                  ENDIF
C..............................................................      
		  IF(IPEV.GE.1)WRITE(6,'(A,3E15.5)')
     *            'XPSAQI 3:XPSAQI,XSTHR,XSMAX',
     *            XPSAQI,XSTHR,XSMAX
                ENDIF
              ENDIF
C                                                      ---
   50         CONTINUE
	      IF(IPEV.GE.1)
     *        WRITE(6,'(A,3E15.4)') ' XKSAMP-50: XPSQI,XSTHR,XSMAX',
     &        XPSQI,XSTHR,XSMAX
   60         CONTINUE
              IF(IPEV.GE.1)WRITE(6,'(A)') ' XKSAMP-60'
              IF(IPEV.GE.1)
     *        WRITE(6,'(A,3E15.4)') ' XKSAMP-60: XPSAQI,XSTHR,XSMAX',
     &        XPSAQI,XSTHR,XSMAX
              XXSEA=XXSEA + XPSQI + XPSAQI
              IF(XXSEA.GE.XXSEAM) THEN
                IXPS=IXPS - ISQ + 1
                GOTO 40
              ENDIF
              IXPS=IXPS+1
              IF(IPEV.GE.1)WRITE(6,'(A,I10)') ' XKSAMP-60: IXPS',IXPS
              XPSQ(IXPS)=XPSQI
              XPSAQ(IXPS)=XPSAQI
C                  Test 14.4.99	      
              XPSQ(IXPS)=XPSAQI
              XPSAQ(IXPS)=XPSQI
              IFROSP(IXPS)=IPP
              ZUOSP(IXPS)=.TRUE.
   70       CONTINUE
          ENDIF
          JSSHS(IPP)=JIPP
*  projectile valence quarks
   80     CONTINUE
          IF(XVTHR.GT.0.05D0)THEN
            IF(XVTHR.GT.1.D0-XXSEA-XDTHR)THEN
              WRITE(6,*)' xvthr,xxsea,xdthr ', XVTHR,XXSEA,XDTHR
            ENDIF
C                  TEST 15.4.99	    
C           XPVQI=BETREJ(0.5D0,UNOPRV,XVTHR,1.D0-XXSEA-XDTHR)
            XPVQI=BETREJ(0.1D0,UNOPRV,XVTHR,1.D0-XXSEA-XDTHR)
   81       CONTINUE
          ELSE
   90       CONTINUE
            IF(IPEV.GE.1)WRITE(6,'(A)') ' XKSAMP-90'
C                  TEST 15.4.99	    
C           XPVQI=DBETAR(0.5D0,UNOPRV)
            XPVQI=DBETAR(0.1D0,UNOPRV)
            IF ((XPVQI.LT.XVTHR).OR.(1.D0-XPVQI-XXSEA.LT.XDTHR))
     *      GOTO 90
          ENDIF
          XPVDI=1. - XPVQI - XXSEA
C                                    CONSISTENCY TEST
C                                    TO BE FULFILLED AUTOMATICALLY
          IF(XPVDI.LT.XDTHR) THEN
            WRITE(6,'(A/A/E12.3,4I4,3E11.3)')
     +      ' INCONSISTENT X-SAMPLING / XKSAMP / PROJECTILE',
     +      ' ECM, IP, IPP, JSSH(IPP), JIPP, XPVQI, XPVDI, XXSEA', ECM,
     +      IP, IPP, JSSH(IPP), JIPP, XPVQI, XPVDI, XXSEA
            STOP
          ENDIF
C
C--------------------------------------------------------------
C                                   diquark rejection
C                      Here we have a projectile diquark
C                      Reject it according to xd**1.5 rule
C--------------------------------------------------------------
          XTEST=XPVDI**1.5
          VV=IPP
C--------------------------------------------------------------
          IXPV=IXPV+1
          XPVQ(IXPV)=XPVQI
          XPVD(IXPV)=XPVDI
          ISXPVQ(IXPV)=0
          ISXPVD(IXPV)=0
          IFROVP(IXPV)=IPP
          ITOVP(IPP)=IXPV
          ZUOVP(IXPV)=.TRUE.
        ENDIF
  100 CONTINUE
C******************************
C                        PARTON X-VALUES OF INTERACTING TARGET NUCLEONS
       IXTV=0
       IXTS=0
       DO 170 ITT=1,IT
         IF (JTSH(ITT).NE.0) THEN
C--------------------------------------------------------------
C                          prepare diquark rejection
C--------------------------------------------------------------
           IIXTSS=IXTS
	   IIXTVV=IXTV
  169      CONTINUE
	   IXTS=IIXTSS
	   IXTV=IIXTVV
C--------------------------------------------------------------
           JITT=JTSH(ITT)-1
           JITT=MIN(JITT,NSMAX)
 111       CONTINUE
           XXSEA=0.0
           IF(JITT.GT.0) THEN
C                     j.r.11.12.97
	     XSMAX=XXSEAM -1.5*JITT*XSTHR
C            XSMAX=XXSEAM - 2.*JITT*XSTHR
             IF(XSTHR.GE.XSMAX) THEN
               JITT=JITT-1
               GOTO 111
             ENDIF
             NSCOUN=0
  110        CONTINUE
             IF(IPEV.GE.1)WRITE(6,'(A)') ' XKSAMP-110'
             XXSEA=0.0
             NSCOUN=NSCOUN+1
             IF (NSCOUN.GT.NSEA)THEN
               JITT=JITT-1
               NSCOUN=0
             ENDIF
             DO 140 ISQ=1,JITT
C                               CHANGE 23.5.90 / 13.9.90
C              IF(XSTHR.GT.0.05D0)THEN
C                                         J.R.29.4.93---
               IF(ITSQ(IXTS+1).LE.2)THEN
C..............................................................  
C                           1/sqrt(x) seaquarks
                 IF(SEAQXN.LE.0.75D0)THEN
                   XTSQI=SAMPEX(XSTHR,XSMAX)
C                           1/x seaquarks
                 ELSEIF(SEAQXN.GT.0.75D0)THEN
                   XTSQI=SAMPEY(XSTHR,XSMAX)
                 ENDIF
C..............................................................      
               ELSE
               IF(XSMAX.GT.XSTHR+BSQMA)THEN
                  XTSQI=SAMPXB(XSTHR+BSQMA,XSMAX,BSQMA)
                ELSE
C..............................................................  
C                           1/sqrt(x) seaquarks
      IF(SEAQXN.LE.0.75D0)THEN
                XTSQI=SAMPEX(XSTHR,XSMAX)
C                           1/x seaquarks
      ELSEIF(SEAQXN.GT.0.75D0)THEN
                XTSQI=SAMPEY(XSTHR,XSMAX)
      ENDIF
C..............................................................      
                ENDIF
              ENDIF
C
              IF(ITSAQ(IXTS+1).GE.-2)THEN
C..............................................................  
C                           1/sqrt(x) seaquarks
      IF(SEAQXN.LE.0.75D0)THEN
                XTSAQI=SAMPEX(XSTHR,XSMAX)
C                           1/x seaquarks
      ELSEIF(SEAQXN.GT.0.75D0)THEN
                XTSAQI=SAMPEY(XSTHR,XSMAX)
      ENDIF
C..............................................................      
              ELSE
                IF(XSMAX.GT.XSTHR+BSQMA)THEN
                  XTSAQI=SAMPXB(XSTHR+BSQMA,XSMAX,BSQMA)
                ELSE
C..............................................................  
C                           1/sqrt(x) seaquarks
      IF(SEAQXN.LE.0.75D0)THEN
                XTSAQI=SAMPEX(XSTHR,XSMAX)
C                           1/x seaquarks
      ELSEIF(SEAQXN.GT.0.75D0)THEN
                XTSAQI=SAMPEY(XSTHR,XSMAX)
      ENDIF
C..............................................................      
                ENDIF
              ENDIF
C                                                  ---
C               XTSQI=SAMPEX(XSTHR,XSMAX)
C
C               XTSAQI=SAMPEX(XSTHR,XSMAX)
C              ELSE
  120           CONTINUE
         IF(IPEV.GE.1)WRITE(6,'(A)') ' XKSAMP-120'
C               XTSQI=SAMPEX(XSTHR,XSMAX)
C               IF (XTSQI.LT.XSTHR.OR.XTSQI.GE.XSMAX)           GOTO 120
  130           CONTINUE
         IF(IPEV.GE.1)WRITE(6,'(A)') ' XKSAMP-130'
C               XTSAQI=SAMPEX(XSTHR,XSMAX)
C               IF (XTSAQI.LT.XSTHR.OR.XTSAQI.GE.XSMAX)         GOTO 130
C             ENDIF
              XXSEA=XXSEA + XTSQI + XTSAQI
              IF(XXSEA.GE.XXSEAM) THEN
                IXTS=IXTS - ISQ + 1
                                                                GOTO 110
              ENDIF
              IXTS=IXTS+1
         IF(IPEV.GE.1)WRITE(6,'(A,I10)')' XKSAMP-130: IXTS',IXTS
              XTSQ(IXTS)=XTSQI
              XTSAQ(IXTS)=XTSAQI
              IFROST(IXTS)=ITT
              ZUOST(IXTS)=.TRUE.
  140       CONTINUE
          ENDIF
          JTSHS(ITT)=JITT
C
C***                             TARGET VALENCE QUARKS
  150     CONTINUE
          IF(XVTHR.GT.0.05D0)THEN
            IF(XVTHR.GT.1.D0-XXSEA-XDTHR)THEN
              WRITE(6,*)' xvthr,xxsea,xdthr ', XVTHR,XXSEA,XDTHR
            ENDIF
C                  TEST 15.4.99	    
C           XTVQI=BETREJ(0.5D0,UNON,XVTHR,1.-XXSEA-XDTHR)
            XTVQI=BETREJ(0.1D0,UNON,XVTHR,1.-XXSEA-XDTHR)
  151     CONTINUE
          ELSE
  160       CONTINUE
         IF(IPEV.GE.1)WRITE(6,'(A)') ' XKSAMP-160'
C                  TEST 15.4.99	    
C           XTVQI=DBETAR(0.5D0,UNON)
            XTVQI=DBETAR(0.1D0,UNON)
	    XMIST=1.-XTVQI-XXSEA
	    IF(IPEV.GE.1)WRITE(6,'(A,5E15.5)')
     *      ' XTVQI,XVTHR,XXSEA,XMIST,XDTHR',
     *        XTVQI,XVTHR,XXSEA,XMIST,XDTHR
          IF((XTVQI.LT.XVTHR).OR.(1.D0-XTVQI-XXSEA.LT.XDTHR+0.0001D0))
     *                                                        GOTO 160
          ENDIF
          XTVDI=1. - XTVQI - XXSEA
C                                    CONSISTENCY TEST
C                                    TO BE FULFILLED AUTOMATICALLY
          IF(XTVDI.LT.XDTHR) THEN
            WRITE(6,'(A/A/E12.3,4I4,3E11.3)')
     +      ' INCONSISTENT X-SAMPLING / XKSAMP / TARGET',
     +      ' ECM, IT, ITT, JTSH(ITT), JITT, XTVQI, XTVDI, XXSEA', ECM,
     +      IT, ITT, JTSH(ITT), JITT, XTVQI, XTVDI, XXSEA
            STOP
          ENDIF
C
C--------------------------------------------------------------
C                                   diquark rejection
C                      Here we have a target diquark
C                      Reject it according to xd**1.5 rule
C--------------------------------------------------------------
          XTEST=XTVDI**1.5
	  VV=ITT
C	  IF(RNDM(VV).GT.XTEST)GO TO 169
C--------------------------------------------------------------
          IXTV=IXTV+1
          XTVQ(IXTV)=XTVQI
          XTVD(IXTV)=XTVDI
          ISXTVQ(IXTV)=0
          ISXTVD(IXTV)=0
          IFROVT(IXTV)=ITT
          ITOVT(ITT)=IXTV
          ZUOVT(IXTV)=.TRUE.
        ENDIF
  170 CONTINUE
C
      IF (IPEV.GE.6) THEN
        WRITE(6,1010)
 1010 FORMAT(' XKSAMP:',
     +' I,XPVQ(I),XPVD(I),IFROVP(I),ITOVP(I),ZUOVP(I),KKPROJ(I)')
        DO 180 I=1,IXPV
          WRITE(6,1020) I,XPVQ(I),XPVD(I),IFROVP(I),ITOVP(I),ZUOVP(I),
     +    KKPROJ(I)
 1020 FORMAT(I5,2E15.5,2I5,L5,I5)
  180   CONTINUE
        WRITE(6,1030)
 1030 FORMAT(' XKSAMP :  I,XPSQ(I),XPSAQ(I),IFROSP(I),ZUOSP(I)')
        DO 190 I=1,IXPS
          WRITE(6,1040) I,XPSQ(I),XPSAQ(I),IFROSP(I),ZUOSP(I)
 1040 FORMAT(I5,2E15.5,I5,L5)
  190   CONTINUE
C
        WRITE(6,1050)
 1050 FORMAT(' XKSAMP:',
     +' I,XTVQ(I),XTVD(I),IFROVT(I),ITOVT(I),ZUOVT(I),KKTARG(I)')
        DO 200 I=1,IXTV
          WRITE(6,1020) I,XTVQ(I),XTVD(I),IFROVT(I),ITOVT(I),ZUOVT(I),
     +    KKTARG(I)
  200   CONTINUE
        WRITE(6,1060)
 1060 FORMAT(' XKSAMP :  I,XTSQ(I),XTSAQ(I),IFROST(I),ZUOST(I)')
        DO 210 I=1,IXTS
          WRITE(6,1040) I,XTSQ(I),XTSAQ(I),IFROST(I),ZUOST(I)
  210   CONTINUE
      ENDIF
      IF(IPEV.GE.6) THEN
        WRITE(6,'(A)')
     +  ' XKSAMP :  I,ITOVP(I),ITOVT(I),JSSHS(I),JTSHS(I)'
        IMA=MAX(IP,IT)
        DO 220 I=1,IMA
          WRITE(6,1070) I,ITOVP(I),ITOVT(I),JSSHS(I),JTSHS(I)
 1070 FORMAT(5I5)
  220   CONTINUE
        DO 181 I=1,IXPV
	  WRITE(6,*)' I,IPVQ(I),IPPV1(I),IPPV2(I) ',
     *	  I,IPVQ(I),IPPV1(I),IPPV2(I)
  181   CONTINUE
        DO 182 I=1,IXTV
	  WRITE(6,*)' I,ITVQ(I),ITTV1(I),ITTV2(I) ',
     *	  I,ITVQ(I),ITTV1(I),ITTV2(I)
  182   CONTINUE
        DO 183 I=1,IXPS
	  WRITE(6,*)' I,IPSQ(I),IPSAQ(I) ',
     *	  I,IPSQ(I),IPSAQ(I)
  183   CONTINUE
        DO 184 I=1,IXTS
	  WRITE(6,*)' I,ITSQ(I),ITSAQ(I) ',
     *	  I,ITSQ(I),ITSAQ(I)
  184   CONTINUE
      ENDIF
C
C----------------------------------------------------------------------
C               COLLECTION OF VALENCE-VALENCE PAIRS
      NVV=0
      IF(IPEV.GE.4)WRITE(6,*)' collect v-v pairs NVV',NVV
      DO 230 I=1,NN
        INTLO(I)=.TRUE.
  230 CONTINUE
      DO 240 I=1,NN
        IIPP=INTER1(I)
        IITT=INTER2(I)
        IIPPV=ITOVP(IIPP)
        IITTV=ITOVT(IITT)
        IF(ZUOVP(IIPPV).AND.ZUOVT(IITTV)) THEN
          INTLO(I)=.FALSE.
      IF(IPEV.GE.6)WRITE(6,'(A,5I5)')
     * ' XKSAMP v-v loop IIPP,IITT,IIPPV,IITTV,NVV',IIPP,IITT,IIPPV,IITTV,NVV
          ZUOVP(IIPPV)=.FALSE.
          ZUOVT(IITTV)=.FALSE.
          NVV=NVV + 1
      IF(IPEV.GE.4)WRITE(6,*)' collect v-v pairs NVV',NVV
          NCHVV1(NVV)=0
          NCHVV2(NVV)=0
          INTVV1(NVV)=IIPPV
          INTVV2(NVV)=IITTV
C -----------------------------------------------------J.R. 6.1.92
C       AMVVP2=XTVQ(IITTV)*XPVD(IIPPV)*ECM*ECM
C       IF(AMVVP2.GT.6.D0)THEN
C                                     RESAMPLE XTVQ
C         XTVQTH=6./(XPVD(IIPPV)*ECM*ECM)
C         XTVQXX=BETREJ(0.5D0,UNOPRV,XTVQTH,XTVQ(IITTV))
C         DXTVQ=XTVQ(IITTV)-XTVQXX
C         XTVQ(IITTV)=XTVQ(IITTV)-DXTVQ
C         XTVD(IITTV)=XTVD(IITTV)+DXTVQ
C       ENDIF
C       AMVVT2=XTVD(IITTV)*XPVQ(IIPPV)*ECM*ECM
C       IF(AMVVT2.GT.6.D0)THEN
C                                     RESAMPLE XPVQ
C         XPVQTH=6./(XTVD(IITTV)*ECM*ECM)
C         XPVQXX=BETREJ(0.5D0,UNOPRV,XPVQTH,XPVQ(IIPPV))
C         DXPVQ=XPVQ(IIPPV)-XPVQXX
C         XPVQ(IIPPV)=XPVQ(IIPPV)-DXPVQ
C         XPVD(IIPPV)=XPVD(IIPPV)+DXPVQ
C       ENDIF
C--------------------------------------------------------------
        ENDIF
  240 CONTINUE
C
C                        COLLECTION OF THE SEA-VALENCE PAIRS
      NDV=0
      NSV=0
      DO 270 I=1,NN
        IF(INTLO(I)) THEN
          IIPP=INTER1(I)
          IITT=INTER2(I)
          IITTV=ITOVT(IITT)
          DO 250 J=1,IXPS
            IF(ZUOSP(J).AND.(IFROSP(J).EQ.IIPP).AND.ZUOVT(IITTV)) THEN
              ZUOSP(J)=.FALSE.
      IF(IPEV.GE.6)WRITE(6,'(A,6I5)')
     *' XKSAMP s-v loop I(NN),J(IXPS),iitt,iittv,NSV,NDV',
     +    I,J, IITT,IITTV,NSV,NDV
              ZUOVT(IITTV)=.FALSE.
              INTLO(I)=.FALSE.
C             IF(LSEADI.AND.RNDM(V).GT.AMEDD.AND.IDIQUA.EQ.1)THEN
              IF(RNDM(V).GT.AMEDD.AND.IDIQUA.EQ.1)THEN
C                               DEFINE D-V CHAINS (SEA-DIQUARK-VALENCE)
                CALL DIQSV(ECM,IITTV,J,IREJ)
                IF(IREJ.EQ.0)GO TO 260
              ENDIF
              NSV=NSV + 1
              NCHSV1(NSV)=0
              NCHSV2(NSV)=0
              INTSV1(NSV)=J
              INTSV2(NSV)=IITTV
C----------------correct sv chains to get minimum mass ------
C     IF(IP.GE.2)GO TO 5270
              AMSVQ1=XPSQ(J)*XTVD(IITTV)*ECM**2
              AMSVQ2=XPSAQ(J)*XTVQ(IITTV)*ECM**2
              JXPV=ITOVP(IIPP)
              IF(IPSQ(J).EQ.3)THEN
                IF(AMSVQ1.GT.AMAS)THEN
                  XPSQXX=(XTVD(IITTV)*ECM**2)
		  IF(XPSQXX.LE.1.D-1)XPSQXX=1.D-1
                  XPSQTH=AMAS/XPSQXX
                  XPSQXX=SAMPEX(XPSQTH,XPSQ(J))
                  DXPSQ=XPSQ(J)-XPSQXX
                  XPSQ(J)=XPSQ(J)-DXPSQ
                  XPVD(JXPV)=XPVD(JXPV)+DXPSQ
                ELSEIF(AMSVQ1.LT.AMAS)THEN
		  IF(XTVD(IITTV)*ECM**2.LE.1.D-12)THEN
		    WRITE(6,*)' xksamp: XTVD(IITTV)=0 ',IITTV
		    XTVD(IITTV)=0.1D0
		  ENDIF
                  XPSQW=AMAS/(XTVD(IITTV)*ECM**2)
                  DXPSQ=XPSQW-XPSQ(J)
                  ISXTVD(IITTV)=1
                  IF(XPVD(JXPV).GE.XDTHR+DXPSQ)THEN 
                    XPVD(JXPV)=XPVD(JXPV)-DXPSQ
                    XPSQ(J)=XPSQW
                  ENDIF
                ENDIF
                IF(AMSVQ2.GT.AMIS)THEN
                ELSEIF(AMSVQ2.LT.AMIS)THEN
		  IF(XTVQ(IITTV)*ECM**2.LE.1.D-12)THEN
		    WRITE(6,*)' xksamp: XTVQ(IITTV)=0 ',IITTV
		    XTVQ(IITTV)=0.1D0
		  ENDIF
                  XPSQW=AMIS/(XTVQ(IITTV)*ECM**2)
                  DXPSQ=XPSQW-XPSAQ(J)
                  ISXTVQ(IITTV)=1
                  IF(XPVD(JXPV).GE.XDTHR+DXPSQ)THEN 
                    XPVD(JXPV)=XPVD(JXPV)-DXPSQ
                    XPSAQ(J)=XPSQW
                  ENDIF
                ENDIF
              ELSE
                IF(AMSVQ1.GT.AMAU)THEN
		  IF(XTVD(IITTV)*ECM**2.LE.1.D-12)THEN
		    WRITE(6,*)' xksamp: XTVD(IITTV)=0 ',IITTV
		    XTVD(IITTV)=0.1D0
		  ENDIF
                  XPSQTH=AMAU/(XTVD(IITTV)*ECM**2)
                  XPSQXX=SAMPEX(XPSQTH,XPSQ(J))
                  DXPSQ=XPSQ(J)-XPSQXX
                  XPSQ(J)=XPSQ(J)-DXPSQ
                  XPVD(JXPV)=XPVD(JXPV)+DXPSQ
                ELSEIF(AMSVQ1.LT.AMAU)THEN
		  IF(XTVD(IITTV)*ECM**2.LE.1.D-12)THEN
		    WRITE(6,*)' xksamp: XTVD(IITTV)=0 ',IITTV
		    XTVD(IITTV)=0.1D0
		  ENDIF
                  XPSQW=AMAU/(XTVD(IITTV)*ECM**2)
                  DXPSQ=XPSQW-XPSQ(J)
                  ISXTVD(IITTV)=1
                  IF(XPVD(JXPV).GE.XDTHR+DXPSQ)THEN 
                    XPVD(JXPV)=XPVD(JXPV)-DXPSQ
                    XPSQ(J)=XPSQW
                  ENDIF
                ENDIF
                IF(AMSVQ2.GT.AMIU)THEN
                ELSEIF(AMSVQ2.LT.AMIU)THEN
		  IF(XTVQ(IITTV)*ECM**2.LE.1.D-12)THEN
		    WRITE(6,*)' xksamp: XTVQ(IITTV)=0 ',IITTV
		    XTVQ(IITTV)=0.1D0
		  ENDIF
                  XPSQW=AMIU/(XTVQ(IITTV)*ECM**2)
                  DXPSQ=XPSQW-XPSAQ(J)
                  ISXTVQ(IITTV)=1
                  IF(XPVD(JXPV).GE.XDTHR+DXPSQ)THEN 
                    XPVD(JXPV)=XPVD(JXPV)-DXPSQ
                    XPSAQ(J)=XPSQW
                  ENDIF
                ENDIF
              ENDIF
C5270 CONTINUE
C-----------------------------------------------------------------
                                                                GOTO 260
            ENDIF
      IF(IPEV.GE.6)WRITE(6,'(A,6I5)')
     *' XKSAMP s-v loop I(NN),J(IXPS),iitt,iittv,NSV,NDV',
     +    I,J, IITT,IITTV,NSV,NDV
  250     CONTINUE
        ENDIF
  260   CONTINUE
  270 CONTINUE
C
C                        COLLECTION OF THE VALENCE-SEA PAIRS
      NVS=0
      NVD=0
      DO 300 I=1,NN
        IF(INTLO(I)) THEN
          IIPP=INTER1(I)
          IITT=INTER2(I)
          IIPPV=ITOVP(IIPP)
          DO 280 J=1,IXTS
            IF(ZUOVP(IIPPV).AND.ZUOST(J).AND.(IFROST(J).EQ.IITT)) THEN
              ZUOST(J)=.FALSE.
              IF(IPEV.GE.6)WRITE(6,*)
     *        ' XKSAMP v-s loop IIPP,IITT,IIPPV,NVS,NVD,I,J,IXTS',
     *        IIPP,IITT,IIPPV,NVS,NVD,I,J,IXTS
              ZUOVP(IIPPV)=.FALSE.
              INTLO(I)=.FALSE.
C             IF(LSEADI.AND.RNDM(V).GT.AMEDD.AND.IDIQUA.EQ.1)THEN
              IF(RNDM(V).GT.AMEDD.AND.IDIQUA.EQ.1)THEN
C                               DEFINE V-D CHAINS (valence - sea diquark)
                CALL DIQVS(ECM,IIPPV,J,IREJ)
                IF(IPEV.GE.6)WRITE(6,*)
     *          ' XKSAMP v-s loop IIPP,IITT,IIPPV,NVS,NVD,I,J,IXTS,JXTV'
     *          ,IIPP,IITT,IIPPV,NVS,NVD,I,J,IXTS,JXTV
                IF(IREJ.EQ.0)GO TO 290
              ENDIF
              NVS=NVS + 1
              NCHVS1(NVS)=0
              NCHVS2(NVS)=0
              INTVS1(NVS)=IIPPV
              INTVS2(NVS)=J
C----------------correct vs chains to get minimum mass ------
              AMVSQ1=XPVQ(IIPPV)*XTSAQ(J)*ECM**2
              AMVSQ2=XPVD(IIPPV)*XTSQ(J)*ECM**2
              JXTV=ITOVT(IITT)
              IF(ITSQ(J).EQ.3)THEN
C               IF(AMVSQ1.GT.AMIS)THEN
                IF(AMVSQ1.LT.AMIS)THEN
		  IF(XPVQ(IIPPV)*ECM**2.LE.1.D-12)THEN
		    WRITE(6,*)' xksamp: XPVQ(IIPPV)=0 ',IIPPV
		    XPVQ(IIPPV)=0.1D0
		  ENDIF
                  XTSQW=AMIS/(XPVQ(IIPPV)*ECM**2)
                  DXTSQ=XTSQW-XTSAQ(J)
                  ISXPVQ(IIPPV)=1
                  IF(XTVD(JXTV).GE.XDTHR+DXTSQ)THEN 
                    XTVD(JXTV)=XTVD(JXTV)-DXTSQ
                    XTSAQ(J)=XTSQW
                  ENDIF
                ENDIF
                IF(AMVSQ2.GT.AMAS)THEN
		  IF(XPVD(IIPPV)*ECM**2.LE.1.D-12)THEN
		    WRITE(6,*)' xksamp: XPVD(IIPPV)=0 ',IIPPV
		    XPVD(IIPPV)=0.1D0
		  ENDIF
                  XTSQTH=AMAS/(XPVD(IIPPV)*ECM**2)
                  XTSQXX=SAMPEX(XTSQTH,XTSQ(J))
                  DXTSQ=XTSQ(J)-XTSQXX
                  XTSQ(J)=XTSQ(J)-DXTSQ
                  XTVD(JXTV)=XTVD(JXTV)+DXTSQ
                ELSEIF(AMVSQ2.LT.AMAS)THEN
		  IF(XPVD(IIPPV)*ECM**2.LE.1.D-12)THEN
		    WRITE(6,*)' xksamp: XPVD(IIPPV)=0 ',IIPPV
		    XPVD(IIPPV)=0.1D0
		  ENDIF
                  XTSQW=AMAS/(XPVD(IIPPV)*ECM**2)
                  DXTSQ=XTSQW-XTSQ(J)
                  ISXPVD(IIPPV)=1
                  IF(XTVD(JXTV).GE.XDTHR+DXTSQ)THEN 
                    XTVD(JXTV)=XTVD(JXTV)-DXTSQ
                    XTSQ(J)=XTSQW
                  ENDIF
                ENDIF
              ELSE
C               IF(AMVSQ1.GT.AMIU)THEN
                IF(AMVSQ1.LT.AMIU)THEN
		  IF(XPVQ(IIPPV)*ECM**2.LE.1.D-12)THEN
		    WRITE(6,*)' xksamp: XPVQ(IIPPV)=0 ',IIPPV
		    XPVQ(IIPPV)=0.1D0
		  ENDIF
                  XTSQW=AMIU/(XPVQ(IIPPV)*ECM**2)
                  DXTSQ=XTSQW-XTSAQ(J)
                  ISXPVQ(IIPPV)=1
                  IF(XTVD(JXTV).GE.XDTHR+DXTSQ)THEN 
                    XTVD(JXTV)=XTVD(JXTV)-DXTSQ
                    XTSAQ(J)=XTSQW
                  ENDIF
                ENDIF
                IF(AMVSQ2.GT.AMAU)THEN
		  IF(XPVD(IIPPV)*ECM**2.LE.1.D-12)THEN
		    WRITE(6,*)' xksamp: XPVD(IIPPV)=0 ',IIPPV
		    XPVD(IIPPV)=0.1D0
		  ENDIF
                  XTSQTH=AMAU/(XPVD(IIPPV)*ECM**2)
                  XTSQXX=SAMPEX(XTSQTH,XTSQ(J))
                  DXTSQ=XTSQ(J)-XTSQXX
                  XTSQ(J)=XTSQ(J)-DXTSQ
                  XTVD(JXTV)=XTVD(JXTV)+DXTSQ
                ELSEIF(AMVSQ2.LT.AMAU)THEN
		  IF(XPVD(IIPPV)*ECM**2.LE.1.D-12)THEN
		    WRITE(6,*)' xksamp: XPVD(IIPPV)=0 ',IIPPV
		    XPVD(IIPPV)=0.1D0
		  ENDIF
                  XTSQW=AMAU/(XPVD(IIPPV)*ECM**2)
                  DXTSQ=XTSQW-XTSQ(J)
                  ISXPVD(IIPPV)=1
                  IF(XTVD(JXTV).GE.XDTHR+DXTSQ)THEN 
                    XTVD(JXTV)=XTVD(JXTV)-DXTSQ
                    XTSQ(J)=XTSQW
                  ENDIF
                ENDIF
              ENDIF
C-----------------------------------------------------------------
                                                                GOTO 290
            ENDIF
  280     CONTINUE
        ENDIF
  290   CONTINUE
  300 CONTINUE
C   End loop:         COLLECTION OF THE VALENCE-SEA PAIRS
C
C                        COLLECTION OF THE SEA-SEA PAIRS
*---------------------   new version 8/03/1991  hjm
      NSS=0
      NDS=0
      NSD=0
      NDZ=0
      NZD=0
      DO 420 I=1,NN
        IF(INTLO(I)) THEN
          IIPP=INTER1(I)
          IITT=INTER2(I)
          DO 400 J=1,IXTS
            IF(ZUOST(J).AND.(IFROST(J).EQ.IITT)) THEN
              DO 390 JJ=1,IXPS
                IF(ZUOSP(JJ).AND.(IFROSP(JJ).EQ.IIPP)) THEN
                  NSS=NSS+1
                  IF(IPEV.GE.6)WRITE(6,'(A,5I5)')
     *            ' XKSAMP s-s loop IIPP,IITT,NSS',IIPP,IITT,NSS
                  NCHSS1(NSS)=0
                  NCHSS2(NSS)=0
                  IF(IPEV.GE.6)WRITE(6,*)
     *            ' XKSAMP s-s loop ,NCHSS1(NSS),NCHSS2(NSS),NSS ',
     *            NCHSS1(NSS),NCHSS2(NSS),NSS
                  INTSS1(NSS)=JJ
                  INTSS2(NSS)=J
                  INTLO(I)=.FALSE.
                  ZUOST(J)=.FALSE.
                  ZUOSP(JJ)=.FALSE.
C-------------------------------------------Mass check j.r.12/94------- 
                  SSMA1Q=XPSQ(JJ)*XTSAQ(J)*ECM**2
                  SSMA2Q=XPSAQ(JJ)*XTSQ(J)*ECM**2
                  IF(SSMA1Q.LT.1.2D0.OR.SSMA2Q.LT.1.2D0) THEN
                    ZUOST(J)=.TRUE.
                    ZUOSP(JJ)=.TRUE.
                    NSS=NSS-1
                    GO TO 410
                  ENDIF
C-------------------------------------------Mass check j.r.12/94------- 
C**********************************************************************
C**********************************************************************
C
C                                     Chain recombination option
C
                  ALLKET=(NVV+IXPS+IXTS)
		  IF(ALLKET.LE.1.D-5)THEN
		    WRITE(6,*)' xksamp ALLKET=0' , ALLKET
		    ALLKET=1.
		  ENDIf
C                 VALFRA=NVV/ALLKET
C                                       j.r.31.3.95
		  ANVVO=MIN(IXPV,IXTV)
		  ANSVO=IXTV-ANVVO
		  ANVSO=IXPV-ANVVO
		  ANSSO=(IXPV+IXPS)-ANVVO-ANSVO-ANVSO
		  IF(ANVVO+ANSSO.LE.1.D-5)THEN
		    WRITE(6,*)' xksamp (...)=0' ,ANVVO,ANSSO
		    ANSSO=1.
		  ENDIf
C		  VALFRA=1.D0
		  IF(ANVVO+ANSSO.GT.1.D-5)VALFRA=ANVVO/(ANVVO+ANSSO)
C                 IF(IRECOM.EQ.1.AND.RNDM(VALFRA).GT.VALFRA)THEN
                  IF(IRECOM.EQ.1)THEN
C--- sea-sea pair found, is there a v-v pair suitable for recombination
C    1. is there a v-v chain pair belonging to same projectile
C    2. is there a v-v chain pair belonging to same target
                    DO 4201 IVV=1,NVV
                      IF (NCHVV1(IVV).NE.99.AND.NCHVV2(IVV).NE.99)THEN
                        IXVPR=INTVV1(IVV)
                        INUCPR=IFROVP(IXVPR)
                        IXVTA=INTVV2(IVV)
                        INUCTA=IFROVT(IXVTA)
                        IF(IIPP.EQ.INUCPR.OR.IITT.EQ.INUCTA)THEN
C   suitable v-v chain pair found, calculate masses of recombined ch's
C   old chains:
C                         SSMA1Q=XPSQ(JJ)*XTSAQ(J)*ECM**2
C                         SSMA2Q=XPSAQ(JJ)*XTSQ(J)*ECM**2
C                         VVMA1Q=XPVQ(IXVPR)*XTVD(IXVTA)*ECM**2
C                         VVMA2Q=XPVD(IXVPR)*XTVQ(IXVTA)*ECM**2
C   new chains:
C                         SVMA1Q=XPSQ(JJ)*XTVD(IXVTA)*ECM**2
C                         SVMA2Q=XPSAQ(JJ)*XTVQ(IXVTA)*ECM**2
C                         VSMA1Q=XPVQ(IXVPR)*XTSAQ(J)*ECM**2
C                         VSMA2Q=XPVD(IXVPR)*XTSQ(J)*ECM**2
C
C                                  drop old v-v and s-s chains
C
                          NCHSS1(NSS)=99
                          NCHSS2(NSS)=99
                          NCHVV1(IVV)=99
                          NCHVV2(IVV)=99
                          IF(IPEV.GE.6)WRITE(6,*)
     *                    ' XKSAMP before DIQSV ,NCHSS1(NSS),',
     *                    'NCHSS2(NSS),NSS ',
     *                    NCHSS1(NSS),NCHSS2(NSS),NSS
C
C                                  assign new s-v and v-s chains
C
C             IF(LSEADI.AND.RNDM(V).GT.AMEDD.AND.IDIQUA.EQ.1)THEN
                          IF(RNDM(V).GT.AMEDD.AND.IDIQUA.EQ.1)THEN
C                           DEFINE D-V CHAINS (SEA-DIQUARK-VALENCE)
                            CALL DIQSV(ECM,IXVTA,JJ,IREJ)
                            IF(IREJ.EQ.0)GO TO 4202
                          ENDIF
                          IF(IPEV.GE.6)WRITE(6,*)
     *			  ' XKSAMP: NSS,NSV,NVS ',NSS,NSV,NVS
                          NSV=NSV+1
                          INTSV1(NSV)=JJ
                          INTSV2(NSV)=IXVTA
                          IF(IPEV.GE.6)WRITE(6,*)
     *		  	  ' XKSAMP: NSS,NSV,NVS ',NSS,NSV,NVS
C----------------correct sv chains to get minimum mass ------
                          AMSVQ1=XPSQ(JJ)*XTVD(IXVTA)*ECM**2
                          AMSVQ2=XPSAQ(JJ)*XTVQ(IXVTA)*ECM**2
                          JXPV=ITOVP(IIPP)
                          IF(IPEV.GE.6)WRITE(6,'(A,5I5)')
     *' XKSAMP s-s loop rec sv,vs IXVTA,JXPV,JJ',IXVTA,JXPV,JJ
                          IF(IPSQ(JJ).EQ.3)THEN
                            IF(AMSVQ1.GT.AMAS)THEN
		              IF(XTVD(IXVTA)*ECM**2.LE.1.D-12)THEN
		                WRITE(6,*)
     *				' xksamp: XTVD(IXVTA)=0 ',IXVTA
		                XTVD(IXVTA)=0.1D0
		              ENDIF
                              XPSQTH=AMAS/(XTVD(IXVTA)*ECM**2)
                              XPSQXX=SAMPEX(XPSQTH,XPSQ(JJ))
                              DXPSQ=XPSQ(JJ)-XPSQXX
                              XPSQ(JJ)=XPSQ(JJ)-DXPSQ
                              XPVD(JXPV)=XPVD(JXPV)+DXPSQ
                            ELSEIF(AMSVQ1.LT.AMAS)THEN
		              IF(XTVD(IXVTA)*ECM**2.LE.1.D-12)THEN
		                WRITE(6,*)
     *			        'xksamp: XTVD(IXVTA)=0 ',IXVTA
		                XTVD(IXVTA)=0.1D0
		              ENDIF
                              XPSQW=AMAS/(XTVD(IXVTA)*ECM**2)
                              DXPSQ=XPSQW-XPSQ(JJ)
                              ISXTVD(IXVTA)=1
                              IF(XPVD(JXPV).GE.XDTHR+DXPSQ)THEN 
                                XPVD(JXPV)=XPVD(JXPV)-DXPSQ
                                XPSQ(JJ)=XPSQW
                              ENDIF
                            ENDIF
C                           IF(AMSVQ2.GT.AMIS)THEN
                            IF(AMSVQ2.LT.AMIS)THEN
		              IF(XTVQ(IXVTA)*ECM**2.LE.1.D-12)THEN
		                WRITE(6,*)
     *				' xksamp: XTVQ(IXVTA)=0 ',IXVTA
		                XTVQ(IXVTA)=0.1D0
		              ENDIF
                              XPSQW=AMIS/(XTVQ(IXVTA)*ECM**2)
                              DXPSQ=XPSQW-XPSAQ(JJ)
                              ISXTVQ(IXVTA)=1
                              IF(XPVD(JXPV).GE.XDTHR+DXPSQ)THEN 
                                XPVD(JXPV)=XPVD(JXPV)-DXPSQ
C                                           s.r.xpsaq statt xpsq! 1294
                                XPSAQ(JJ)=XPSQW
                              ENDIF
                            ENDIF
                          ELSE
                            IF(AMSVQ1.GT.AMAU)THEN
		              IF(XTVD(IXVTA)*ECM**2.LE.1.D-12)THEN
		                WRITE(6,*)
     *				' xksamp: XTVD(IXVTA)=0 ',IXVTA
		                XTVD(IXVTA)=0.1D0
		              ENDIF
                              XPSQTH=AMAU/(XTVD(IXVTA)*ECM**2)
                              XPSQXX=SAMPEX(XPSQTH,XPSQ(JJ))
                              DXPSQ=XPSQ(JJ)-XPSQXX
                              XPSQ(JJ)=XPSQ(JJ)-DXPSQ
                              XPVD(JXPV)=XPVD(JXPV)+DXPSQ
                            ELSEIF(AMSVQ1.LT.AMAU)THEN
		              IF(XTVD(IXVTA)*ECM**2.LE.1.D-12)THEN
		                WRITE(6,*)
     *				' xksamp: XTVD(IXVTA)=0 ',IXVTA
		                XTVD(IXVTA)=0.1D0
		              ENDIF
                              XPSQW=AMAU/(XTVD(IXVTA)*ECM**2)
                              DXPSQ=XPSQW-XPSQ(JJ)
                              ISXTVD(IXVTA)=1
                              IF(XPVD(JXPV).GE.XDTHR+DXPSQ)THEN 
                                XPVD(JXPV)=XPVD(JXPV)-DXPSQ
                                XPSQ(JJ)=XPSQW
                              ENDIF
                            ENDIF
C                           IF(AMSVQ2.GT.AMIU)THEN
                            IF(AMSVQ2.LT.AMIU)THEN
		              IF(XTVQ(IXVTA)*ECM**2.LE.1.D-12)THEN
		                WRITE(6,*)
     *				' xksamp: XTVQ(IXVTA)=0 ',IXVTA
		                XTVQ(IXVTA)=0.1D0
		              ENDIF
                              XPSQW=AMIU/(XTVQ(IXVTA)*ECM**2)
                              DXPSQ=XPSQW-XPSAQ(JJ)
                              ISXTVQ(IXVTA)=1
                              IF(XPVD(JXPV).GE.XDTHR+DXPSQ)THEN 
                                XPVD(JXPV)=XPVD(JXPV)-DXPSQ
C                                           s.r.xpsaq statt xpsq! 1294
                                XPSAQ(JJ)=XPSQW
                              ENDIF
                            ENDIF
                          ENDIF
 4202                     CONTINUE
C-----------------------------------------------------------------
C
C                                  assign new s-v and v-s chains
C
C             IF(LSEADI.AND.RNDM(V).GT.AMEDD.AND.IDIQUA.EQ.1)THEN
                          IF(RNDM(V).GT.AMEDD.AND.IDIQUA.EQ.1)THEN
C                               DEFINE V-D CHAINS (valence - sea diquark)
                            CALL DIQVS(ECM,IXVPR,J,IREJ)
                            IF(IREJ.EQ.0)GO TO 4203
                          ENDIF
                          IF(IPEV.GE.6)WRITE(6,*)
     *			  ' XKSAMP: NSS,NSV,NVS ',NSS,NSV,NVS
                          NVS=NVS+1
                          INTVS1(NVS)=IXVPR
                          INTVS2(NVS)=J
                          IF(IPEV.GE.6)WRITE(6,*)
     *			  ' XKSAMP: NSS,NSV,NVS ',NSS,NSV,NVS
C------###-------correct vs chains to get minimum mass ------
                          AMVSQ1=XPVQ(IXVPR)*XTSAQ(J)*ECM**2
                          AMVSQ2=XPVD(IXVPR)*XTSQ(J)*ECM**2
                          JXTV=ITOVT(IITT)
                          IF(IPEV.GE.6)WRITE(6,'(A,5I5)')
     *' XKSAMP s-s loop rec vs IXVTA,JXPV,JJ',IXVTA,JXPV,JJ
                          IF(ITSQ(J).EQ.3)THEN
C                           IF(AMVSQ1.GT.AMIS)THEN
                            IF(AMVSQ1.LT.AMIS)THEN
		              IF(XPVQ(IXVPR)*ECM**2.LE.1.D-12)THEN
		                WRITE(6,*)
     *				' xksamp: XPVQ(IXVPR)=0 ',IXVPR
		                XPVQ(IXVPR)=0.1D0
		              ENDIF
                              XTSQW=AMIS/(XPVQ(IXVPR)*ECM**2)
                              DXTSQ=XTSQW-XTSAQ(J)
                              ISXPVQ(IXVPR)=1
                              IF(XTVD(JXTV).GE.XDTHR+DXTSQ)THEN 
                                XTVD(JXTV)=XTVD(JXTV)-DXTSQ
                                XTSAQ(J)=XTSQW
                              ENDIF
                            ENDIF
                            IF(AMVSQ2.GT.AMAS)THEN
		              IF(XPVD(IXVPR)*ECM**2.LE.1.D-12)THEN
		                WRITE(6,*)
     *				' xksamp: XPVD(IXVPR)=0 ',IXVPR
		                XPVD(IXVPR)=0.1D0
		              ENDIF
                              XTSQTH=AMAS/(XPVD(IXVPR)*ECM**2)
                              XTSQXX=SAMPEX(XTSQTH,XTSQ(J))
                              DXTSQ=XTSQ(J)-XTSQXX
                              XTSQ(J)=XTSQ(J)-DXTSQ
                              XTVD(JXTV)=XTVD(JXTV)+DXTSQ
                            ELSEIF(AMVSQ2.LT.AMAS)THEN
		              IF(XPVD(IXVPR)*ECM**2.LE.1.D-12)THEN
		                WRITE(6,*)
     *				' xksamp: XPVD(IXVPR)=0 ',IXVPR
		                XPVD(IXVPR)=0.1D0
		              ENDIF
                              XTSQW=AMAS/(XPVD(IXVPR)*ECM**2)
                              ISXPVD(IXVPR)=1
                              DXTSQ=XTSQW-XTSQ(J)
                              IF(XTVD(JXTV).GE.XDTHR+DXTSQ)THEN 
                                XTVD(JXTV)=XTVD(JXTV)-DXTSQ
                                XTSQ(J)=XTSQW
                              ENDIF
                            ENDIF
                          ELSE
C                           IF(AMVSQ1.GT.AMIU)THEN
                            IF(AMVSQ1.LT.AMIU)THEN
		              IF(XPVQ(IXVPR)*ECM**2.LE.1.D-12)THEN
		                WRITE(6,*)
     *				' xksamp: XPVQ(IXVPR)=0 ',IXVPR
		                XPVQ(IXVPR)=0.1D0
		              ENDIF
                              XTSQW=AMIU/(XPVQ(IXVPR)*ECM**2)
                              DXTSQ=XTSQW-XTSAQ(J)
                              ISXPVQ(IXVPR)=1
                              IF(XTVD(JXTV).GE.XDTHR+DXTSQ)THEN 
                                XTVD(JXTV)=XTVD(JXTV)-DXTSQ
                                XTSAQ(J)=XTSQW
                              ENDIF
                            ENDIF
                            IF(AMVSQ2.GT.AMAU)THEN
		              IF(XPVD(IXVPR)*ECM**2.LE.1.D-12)THEN
		                WRITE(6,*)
     *				' xksamp: XPVD(IXVPR)=0 ',IXVPR
		                XPVD(IXVPR)=0.1D0
		              ENDIF
                              XTSQTH=AMAU/(XPVD(IXVPR)*ECM**2)
                              XTSQXX=SAMPEX(XTSQTH,XTSQ(J))
                              DXTSQ=XTSQ(J)-XTSQXX
                              XTSQ(J)=XTSQ(J)-DXTSQ
                              XTVD(JXTV)=XTVD(JXTV)+DXTSQ
                            ELSEIF(AMVSQ2.LT.AMAU)THEN
		              IF(XPVD(IXVPR)*ECM**2.LE.1.D-12)THEN
		                WRITE(6,*)
     *				' xksamp: XPVD(IXVPR)=0 ',IXVPR
		                XPVD(IXVPR)=0.1D0
		              ENDIF
                              XTSQW=AMAU/(XPVD(IXVPR)*ECM**2)
                              DXTSQ=XTSQW-XTSQ(J)
                              ISXPVD(IXVPR)=1
                              IF(XTVD(JXTV).GE.XDTHR+DXTSQ)THEN 
                                XTVD(JXTV)=XTVD(JXTV)-DXTSQ
                                XTSQ(J)=XTSQW
                              ENDIF
                            ENDIF
                          ENDIF
 4203                     CONTINUE
C-----------------------------------------------------------------
C
C                                 jump out of s-s chain loop
C
                          GO TO 420
                        ENDIF
                      ENDIF
 4201               CONTINUE
                  ENDIF
C     of loop recombination  IF(IRECOM.EQ.1)THEN
C**********************************************************************
C        we continue in s-s loop
C**********************************************************************
C
C                 IF(LSEADI.AND.RNDM(V).GT.AMEDD.AND.IDIQUA.EQ.1)THEN
                  IF(RNDM(V).GT.2.D0*AMEDD-1.D0.AND.IDIQUA.EQ.1)THEN
C                               DEFINE D-S CHAINS (SEA-DIQUARK---SEA)
                    CALL DIQDSS(ECM,J,JJ,IREJ)
                    IF(IREJ.EQ.0) THEN
                      NCHSS1(NSS)=99
                      NCHSS2(NSS)=99
                      IF(IPEV.GE.6)WRITE(6,*)
     *                ' XKSAMP AFTER DIQDSS IREJ=0',
     *                ',NCHSS1(NSS),NCHSS2(NSS),NSS ',
     *                NCHSS1(NSS),NCHSS2(NSS),NSS
                      GO TO 410
                    ENDIF
                  ENDIF
C                 IF(LSEADI.AND.RNDM(V).GT.AMEDD.AND.IDIQUA.EQ.1)THEN
                  IF(RNDM(V).GT.2.D0*AMEDD-1.D0.AND.IDIQUA.EQ.1)THEN
C                               DEFINE S-D CHAINS (SEA---SEA-DIQUARK)
                    CALL DIQSSD(ECM,J,JJ,IREJ)
                    IF(IREJ.EQ.0) THEN
                      NCHSS1(NSS)=99
                      NCHSS2(NSS)=99
                      IF(IPEV.GE.6)WRITE(6,*)
     *                ' XKSAMP AFTER DIQSSD IREJ=0',
     *                ',NCHSS1(NSS),NCHSS2(NSS),NSS ',
     *                NCHSS1(NSS),NCHSS2(NSS),NSS
                      GO TO 410
                    ENDIF
                  ENDIF
                  SSMA1Q=XPSQ(JJ)*XTSAQ(J)*ECM**2
                  SSMA2Q=XPSAQ(JJ)*XTSQ(J)*ECM**2
                  IF(SSMA1Q.LT.SSMIMQ.OR.SSMA2Q.LT.SSMIMQ) THEN
                    JXPV=ITOVP(IIPP)
                    JXTV=ITOVT(IITT)
                    IF((XTVD(JXTV).GT.XDTHR+3.5D0*XSSTHR)
     *		    .AND.(XPVD(JXPV)
     +              .GT.XDTHR+3.5D0*XSSTHR)) THEN
*  maximum allowed x values for sea quarks
                      XSPMAX=1.0 - XPVQ(JXPV) - XDTHR - 1.2*XSSTHR
                      XSTMAX=1.0 - XTVQ(JXTV) - XDTHR - 1.2*XSSTHR
*  resampling of x values not possible / discard s-s interaction
                      IF((XSPMAX.LE.XSSTHR+0.05D0) .OR.(XSTMAX.LE.XSSTHR
     +                +0.05D0))                                 GOTO 380
*  resampling for projectile sea quark pair
                      ICOUS=0
  310                 CONTINUE
                      ICOUS=ICOUS + 1
                      IF(XSSTHR.GT.0.05D0) THEN
                        XPSQI=BETREJ(XSEACU,UNOSEA,XSSTHR,XSPMAX)
                        XPSAQI=BETREJ(XSEACU,UNOSEA,XSSTHR,XSPMAX)
                      ELSE
  320                   CONTINUE
                        XPSQI=DBETAR(XSEACU,UNOSEA)
                        IF(XPSQI.LT.XSSTHR.OR.XPSQI.GT.XSPMAX)  GOTO 320
  330                   CONTINUE
                        XPSAQI=DBETAR(XSEACU,UNOSEA)
                        IF(XPSAQI.LT.XSSTHR.OR.XPSAQI.GT.XSPMAX)
     +                                                          GOTO 330
                      ENDIF
*  final test of remaining x for projectile diquark
                      XPVDCO=XPVD(JXPV) - XPSQI - XPSAQI + XPSQ(JJ) +
     +                XPSAQ(JJ)
                      IF(XPVDCO.GT.XDTHR) THEN
*  projectile x sampling ok / continue with target sea
                        GOTO 340
                      ELSEIF(ICOUS.LT.5) THEN
                        GOTO 310
                      ELSE
*  too many unsuccessful attempts / discard s-s interaction
                        GOTO 380
                      ENDIF
*  resampling for target sea quark pair
  340                 CONTINUE
                      ICOUS=0
  350                 CONTINUE
                      ICOUS=ICOUS + 1
                      IF(XSSTHR.GT.0.05D0)THEN
                        XTSQI=BETREJ(XSEACU,UNOSEA,XSSTHR,XSTMAX)
                        XTSAQI=BETREJ(XSEACU,UNOSEA,XSSTHR,XSTMAX)
                      ELSE
  360                   CONTINUE
                        XTSQI=DBETAR(XSEACU,UNOSEA)
                        IF(XTSQI.LT.XSSTHR.OR.XTSQI.GT.XSTMAX)  GOTO 360
  370                   CONTINUE
                        XTSAQI=DBETAR(XSEACU,UNOSEA)
                        IF(XTSAQI.LT.XSSTHR.OR.XTSAQI.GT.XSTMAX)
     +                                                          GOTO 370
                      ENDIF
*   final test of remaining x for target diquark
                      XTVDCO=XTVD(JXTV) - XTSQI - XTSAQI + XTSQ(J) +
     +                XTSAQ(J)
                      IF(XTVDCO.LT.XDTHR) THEN
*   repeat x sampling for target sea quarks
                        IF(ICOUS.LT.5)                          GOTO 350
*   discard s-s interaction / too many unsuccessful trials
                                                                GOTO 380
                      ENDIF
*   modification of x values acceptable
                      XPVD(JXPV)=XPVDCO
                      XTVD(JXTV)=XTVDCO
                      XPSQ(JJ)=XPSQI
                      XPSAQ(JJ)=XPSAQI
                      XTSQ(J)=XTSQI
                      XTSAQ(J)=XTSAQI
                                                                GOTO 410
*   consider next s-s interaction
                    ENDIF
*   discard s-s interaction
*   resampling of x values not allowed or unsuccessful
  380               CONTINUE
                    INTLO(I)=.FALSE.
                    ZUOST(J)=.TRUE.
                    ZUOSP(JJ)=.TRUE.
                    NSS=NSS - 1
                  ENDIF
*   consider next s-s interaction
                                                                GOTO 410
                ENDIF
  390         CONTINUE
            ENDIF
  400     CONTINUE
        ENDIF
  410   CONTINUE
  420 CONTINUE
C
C                        CORRECT X-VALUES OF VALENCE QUARKS
C                        FOR NON-MATCHING SEA QUARKS
      DO 430 I=1,IXPS
        IF(ZUOSP(I)) THEN
          IIFROP=IFROSP(I)
          IITOP=ITOVP(IIFROP)
          XPVQ(IITOP)=XPVQ(IITOP) + XPSQ(I) + XPSAQ(I)
          ZUOSP(I)=.FALSE.
        ENDIF
  430 CONTINUE
      DO 440 I=1,IXTS
        IF(ZUOST(I)) THEN
          IIFROT=IFROST(I)
          IITOT=ITOVT(IIFROT)
          XTVQ(IITOT)=XTVQ(IITOT) + XTSQ(I) + XTSAQ(I)
          ZUOST(I)=.FALSE.
        ENDIF
  440 CONTINUE
C
      DO 450 I=1,IXPV
        IF(ZUOVP(I)) THEN
          IPIP=IFROVP(I)
          ISTHKK(IPIP)=13
        ENDIF
  450 CONTINUE
      DO 460 I=1,IXTV
        IF(ZUOVT(I)) THEN
          ITIT=IFROVT(I)
          ISTHKK(ITIT+IP)=14
        ENDIF
  460 CONTINUE
C
      IF(IPEV.GE.6) THEN
        WRITE(6,'(A)') ' XKSAMP: I,INTVV1,INTVV2,IFROVP,IFROVT'
        DO 470 I=1,NVV
          INUP=INTVV1(I)
          INUT=INTVV2(I)
          WRITE(6,'(5I5)') I,INUP,INUT,IFROVP(INUP),IFROVT(INUT)
  470   CONTINUE
       WRITE(6,'(A)')'XKSAMP:I(NSV),INTSV1,INTSV2,IFROSP,IFROVT'
        DO 480 I=1,NSV
          INUP=INTSV1(I)
          INUT=INTSV2(I)
          WRITE(6,'(5I5)') I,INUP,INUT,IFROSP(INUP),IFROVT(INUT)
  480   CONTINUE
        WRITE(6,'(A)') ' XKSAMP: I,INTVS1,INTVS2,IFROVP,IFROST'
        DO 490 I=1,NVS
          INUP=INTVS1(I)
          INUT=INTVS2(I)
          WRITE(6,'(5I5)') I,INUP,INUT,IFROVP(INUP),IFROST(INUT)
  490   CONTINUE
        WRITE(6,'(A)') ' XKSAMP: I,INTSS1,INTSS2,IFROSP,IFROST'
        DO 500 I=1,NSS
          INUP=INTSS1(I)
          INUT=INTSS2(I)
          WRITE(6,'(5I5)') I,INUP,INUT,IFROSP(INUP),IFROST(INUT)
  500   CONTINUE
C
        WRITE(6,'(A)')
     +  ' XKSAMP :  FINAL X-VALUES AFTER POTENTIAL CORRECTION'
        WRITE(6,1010)
        DO 510 I=1,IXPV
          WRITE(6,1020) I,XPVQ(I),XPVD(I),IFROVP(I),ITOVP(I),ZUOVP(I)
          WRITE(6,*)' I(1-IXPV),IPVQ(I),IPPV1(I),IPPV2(I) ',
     *    I,IPVQ(I),IPPV1(I),IPPV2(I)
  510   CONTINUE
        WRITE(6,1030)
        DO 520 I=1,IXPS
          WRITE(6,1040) I,XPSQ(I),XPSAQ(I),IFROSP(I),ZUOSP(I)
          WRITE(6,*)' I(1-IXPS),IPSQ(I),IPSAQ(I) ',
     *    I,IPSQ(I),IPSAQ(I)
  520   CONTINUE
        WRITE(6,1050)
        DO 530 I=1,IXTV
          WRITE(6,1020) I,XTVQ(I),XTVD(I),IFROVT(I),ITOVT(I),ZUOVT(I)
          WRITE(6,*)' I(1-IXTV),ITVQ(I),ITTV1(I),ITTV2(I) ',
     *    I,ITVQ(I),ITTV1(I),ITTV2(I)
  530   CONTINUE
        WRITE(6,1060)
        DO 540 I=1,IXTS
          WRITE(6,1040) I,XTSQ(I),XTSAQ(I),IFROST(I),ZUOST(I)
          WRITE(6,*)' I(1-IXTS),ITSQ(I),ITSAQ(I) ',
     *    I,ITSQ(I),ITSAQ(I)
  540   CONTINUE
      ENDIF
      IF(IPEV.GE.6)WRITE(6,'(A,6I5)')
     *' XKSAMP NSV,NDV,NVS,NVD',
     +    NSV,NDV,NVS,NVD
*  store properties of interacting partons into /HKKEVT/
      CALL PARHKK
      RETURN
      END
*-- Author :
*
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
      SUBROUTINE PARHKK
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C                         STORE INTERACTING PARTONS IN /HKKEVT/
C                         X-VALUES STORED IN PHKK(3,...) AND PHKK(4,...)
C                         POSITIONS OF NUCLEONS STORED IN VHKK
C                         FLAG FOR PROJECTILE VALENCE: ISTHKK=21
C                                  PROJECTILE SEA    : ISTHKK=31
C                         FLAG FOR TARGET VALENCE    : ISTHKK=22
C                                  TARGET SEA        : ISTHKK=32
*KEEP,INTMX.
      PARAMETER (INTMX=2488,INTMD=252)
*KEEP,DXQX.
C     INCLUDE (XQXQ)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON /DXQX/ XPVQ(248),XPVD(248),XTVQ(248),XTVD(248), XPSQ
     +(INTMX),XPSAQ(INTMX),XTSQ(INTMX),XTSAQ(INTMX)
     *                ,XPSU(248),XTSU(248)
     *                ,XPSUT(248),XTSUT(248)
*KEEP,INTNEW.
      COMMON /INTNEW/ NVV,NSV,NVS,NSS,NDV,NVD,NDS,NSD,
     +IXPV,IXPS,IXTV,IXTS, INTVV1(248),
     +INTVV2(248),INTSV1(248),INTSV2(248), INTVS1(248),INTVS2(248),
     +INTSS1(INTMX),INTSS2(INTMX),
     +INTDV1(248),INTDV2(248),INTVD1(248),INTVD2(248),
     +INTDS1(INTMD),INTDS2(INTMD),INTSD1(INTMD),INTSD2(INTMD)
 
C  /INTNEW/
C       NVV    :   NUMBER OF INTERACTING VALENCE-VALENCE SYSTEMS
C       NSV    :   NUMBER OF INTERACTING SEA-VALENCE SYSTEMS
C       NVS    :   NUMBER OF INTERACTING VALENCE-SEA SYSTEMS
C       NSS    :   NUMBER OF INTERACTING SEA-SEA SYSTEMS
C       IXPV, IXTV : NUMBER OF GENERATED X-VALUES FOR VALENCE-QUARK
C                     SYSTEMS FROM PROJECTILE/TARGET NUCLEI
C       IXPS, IXTS : NUMBER OF GENERATED X-VALUES FOR SEA-QUARK PAIRS
C                     FROM PROJECTILE/TARGET NUCLEI
C-------------------
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,LOZUO.
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/ ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),ZUOST(INTMX),
     +INTLO(INTMX),INLOSS(INTMX)
C  /LOZUO/
C       INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C------------------
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
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
*KEEP,SHMAKL.
C     INCLUDE (SHMAKL)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON/SHMAKL/JSSH(INTMX),JTSH(INTMX),INTER1(INTMX),INTER2(INTMX)
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEND.
C----------------------------------
      DO 10 I=1,IXPV
        NHKK=NHKK+1
        IF (NHKK.EQ.NMXHKK)THEN
          WRITE (6,'(A,2I5)') ' XKSAMP: NHKK.EQ.NMXHKK ',NHKK,NMXHKK
          RETURN
        ENDIF
        ISTHKK(NHKK)=21
        KKKHKK=IFROVP(I)
        KKK=JHKKNP(KKKHKK)
        JMOHKK(1,NHKK)=KKK
        JMOHKK(2,NHKK)=0
        JDAHKK(1,NHKK)=0
        JDAHKK(2,NHKK)=0
        PHKK(1,NHKK)=0.
        PHKK(2,NHKK)=0.
        PHKK(3,NHKK)=XPVQ(I)
        PHKK(4,NHKK)=XPVQ(I)
        PHKK(5,NHKK)=0.
C           Add here position of parton in hadron
        VHKK(1,NHKK)=VHKK(1,KKK)
        VHKK(2,NHKK)=VHKK(2,KKK)
        VHKK(3,NHKK)=VHKK(3,KKK)
        VHKK(4,NHKK)=0.
C
        IF (IPHKK.GE.3) WRITE(6,1000) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     +  JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +  (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
C
        NHKK=NHKK+1
        IF (NHKK.EQ.NMXHKK)THEN
          WRITE (6,'(A,2I5)') ' XKSAMP: NHKK.EQ.NMXHKK ',NHKK,NMXHKK
          RETURN
        ENDIF
        ISTHKK(NHKK)=21
C         KKKHKK=IFROVP(I)
        KKK=JHKKNP(KKKHKK)
        JMOHKK(1,NHKK)=KKK
        JMOHKK(2,NHKK)=0
        JDAHKK(1,NHKK)=0
        JDAHKK(2,NHKK)=0
        PHKK(1,NHKK)=0.
        PHKK(2,NHKK)=0.
        PHKK(3,NHKK)=XPVD(I)
        PHKK(4,NHKK)=XPVD(I)
        PHKK(5,NHKK)=0.
C           Add here position of parton in hadron
        VHKK(1,NHKK)=VHKK(1,KKK)
        VHKK(2,NHKK)=VHKK(2,KKK)
        VHKK(3,NHKK)=VHKK(3,KKK)
        VHKK(4,NHKK)=0.
        JHKKPV(I)=NHKK
C
        IF (IPHKK.GE.7) WRITE(6,1000) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     +  JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +  (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
 1000 FORMAT (I6,I4,5I6,9E10.2)
   10 CONTINUE
C                                             ****  PROJECTILE SEA
      DO 20 I=1,IXPS
        NHKK=NHKK+1
        IF (NHKK.EQ.NMXHKK)THEN
          WRITE (6,'(A,2I5)') ' XKSAMP: NHKK.EQ.NMXHKK ',NHKK,NMXHKK
          RETURN
        ENDIF
        ISTHKK(NHKK)=31
        KKKHKK=IFROSP(I)
        KKK=JHKKNP(KKKHKK)
        JMOHKK(1,NHKK)=KKK
        JMOHKK(2,NHKK)=0
        JDAHKK(1,NHKK)=0
        JDAHKK(2,NHKK)=0
        PHKK(1,NHKK)=0.
        PHKK(2,NHKK)=0.
        PHKK(3,NHKK)=XPSQ(I)
        PHKK(4,NHKK)=XPSQ(I)
        PHKK(5,NHKK)=0.
C           Add here position of parton in hadron
        VHKK(1,NHKK)=VHKK(1,KKK)
        VHKK(2,NHKK)=VHKK(2,KKK)
        VHKK(3,NHKK)=VHKK(3,KKK)
        VHKK(4,NHKK)=0.
C
        IF (IPHKK.GE.7) WRITE(6,1000) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     +  JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +  (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
C
        NHKK=NHKK+1
        IF (NHKK.EQ.NMXHKK)THEN
          WRITE (6,'(A,2I5)') ' XKSAMP: NHKK.EQ.NMXHKK ',NHKK,NMXHKK
          RETURN
        ENDIF
        ISTHKK(NHKK)=31
        KKKHKK=IFROSP(I)
        KKK=JHKKNP(KKKHKK)
        JMOHKK(1,NHKK)=KKK
        JMOHKK(2,NHKK)=0
        JDAHKK(1,NHKK)=0
        JDAHKK(2,NHKK)=0
        PHKK(1,NHKK)=0.
        PHKK(2,NHKK)=0.
        PHKK(3,NHKK)=XPSAQ(I)
        PHKK(4,NHKK)=XPSAQ(I)
        PHKK(5,NHKK)=0.
C           Add here position of parton in hadron
        VHKK(1,NHKK)=VHKK(1,KKK)
        VHKK(2,NHKK)=VHKK(2,KKK)
        VHKK(3,NHKK)=VHKK(3,KKK)
        VHKK(4,NHKK)=0.
        JHKKPS(I)=NHKK
C
        IF (IPHKK.GE.7) WRITE(6,1000) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     +  JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +  (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
   20 CONTINUE
C                                             *****  TARGET VALENCE
      DO 30 I=1,IXTV
        NHKK=NHKK+1
        IF (NHKK.EQ.NMXHKK)THEN
          WRITE (6,'(A,2I5)') ' XKSAMP: NHKK.EQ.NMXHKK ',NHKK,NMXHKK
          RETURN
        ENDIF
        ISTHKK(NHKK)=22
        KKKHKK=IFROVT(I)
        KKK=JHKKNT(KKKHKK)
        JMOHKK(1,NHKK)=KKK
        JMOHKK(2,NHKK)=0
        JDAHKK(1,NHKK)=0
        JDAHKK(2,NHKK)=0
        PHKK(1,NHKK)=0.
        PHKK(2,NHKK)=0.
        PHKK(3,NHKK)=XTVQ(I)
        PHKK(4,NHKK)=XTVQ(I)
        PHKK(5,NHKK)=0.
C           Add here position of parton in hadron
        VHKK(1,NHKK)=VHKK(1,KKK)
        VHKK(2,NHKK)=VHKK(2,KKK)
        VHKK(3,NHKK)=VHKK(3,KKK)
        VHKK(4,NHKK)=0.
C
        IF (IPHKK.GE.7) WRITE(6,1000) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     +  JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +  (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
C
        NHKK=NHKK+1
        IF (NHKK.EQ.NMXHKK)THEN
          WRITE (6,'(A,2I5)') ' XKSAMP: NHKK.EQ.NMXHKK ',NHKK,NMXHKK
          RETURN
        ENDIF
        ISTHKK(NHKK)=22
        KKKHKK=IFROVT(I)
        KKK=JHKKNT(KKKHKK)
        JMOHKK(1,NHKK)=KKK
        JMOHKK(2,NHKK)=0
        JDAHKK(1,NHKK)=0
        JDAHKK(2,NHKK)=0
        PHKK(1,NHKK)=0.
        PHKK(2,NHKK)=0.
        PHKK(3,NHKK)=XTVD(I)
        PHKK(4,NHKK)=XTVD(I)
        PHKK(5,NHKK)=0.
C           Add here position of parton in hadron
        VHKK(1,NHKK)=VHKK(1,KKK)
        VHKK(2,NHKK)=VHKK(2,KKK)
        VHKK(3,NHKK)=VHKK(3,KKK)
        VHKK(4,NHKK)=0.
        JHKKTV(I)=NHKK
C
        IF (IPHKK.GE.7) WRITE(6,1000) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     +  JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +  (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
   30 CONTINUE
C                                              *****  TARGET SEA
      DO 40 I=1,IXTS
        NHKK=NHKK+1
        IF (NHKK.EQ.NMXHKK)THEN
          WRITE (6,'(A,2I5)') ' XKSAMP: NHKK.EQ.NMXHKK ',NHKK,NMXHKK
          RETURN
        ENDIF
        ISTHKK(NHKK)=32
        KKKHKK=IFROST(I)
        KKK=JHKKNT(KKKHKK)
        JMOHKK(1,NHKK)=KKK
        JMOHKK(2,NHKK)=0
        JDAHKK(1,NHKK)=0
        JDAHKK(2,NHKK)=0
        PHKK(1,NHKK)=0.
        PHKK(2,NHKK)=0.
        PHKK(3,NHKK)=XTSQ(I)
        PHKK(4,NHKK)=XTSQ(I)
        PHKK(5,NHKK)=0.
C           Add here position of parton in hadron
        VHKK(1,NHKK)=VHKK(1,KKK)
        VHKK(2,NHKK)=VHKK(2,KKK)
        VHKK(3,NHKK)=VHKK(3,KKK)
        VHKK(4,NHKK)=0.
C
        IF (IPHKK.GE.7) WRITE(6,1000) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     +  JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +  (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
C
        NHKK=NHKK+1
        IF (NHKK.EQ.NMXHKK)THEN
          WRITE (6,'(A,2I5)') ' XKSAMP: NHKK.EQ.NMXHKK ',NHKK,NMXHKK
          RETURN
        ENDIF
        ISTHKK(NHKK)=32
        KKKHKK=IFROST(I)
        KKK=JHKKNT(KKKHKK)
        JMOHKK(1,NHKK)=KKK
        JMOHKK(2,NHKK)=0
        JDAHKK(1,NHKK)=0
        JDAHKK(2,NHKK)=0
        PHKK(1,NHKK)=0.
        PHKK(2,NHKK)=0.
        PHKK(3,NHKK)=XTSAQ(I)
        PHKK(4,NHKK)=XTSAQ(I)
        PHKK(5,NHKK)=0.
C           Add here position of parton in hadron
        VHKK(1,NHKK)=VHKK(1,KKK)
        VHKK(2,NHKK)=VHKK(2,KKK)
        VHKK(3,NHKK)=VHKK(3,KKK)
        VHKK(4,NHKK)=0.
        JHKKTS(I)=NHKK
C
        IF (IPHKK.GE.7) WRITE(6,1000) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     +  JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +  (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
   40 CONTINUE
      RETURN
      END
*-- Author :
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HADRKK(NHKKH1,PPN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C  HADRKK CONSTRUCTS ONE HADRONIZED EVENT FOR KK-COLLISIONS
C     ALL TYPES OF CHAINS ARE CONSIDERED
C     OPTONALLY GIVEN TYPES ARE SELECTED ACCORDING TO /DROPPT/
C
C--------------------------------------------------------------------
*KEEP,INTMX.
      PARAMETER (INTMX=2488,INTMD=252)
*KEEP,DXQX.
C     INCLUDE (XQXQ)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON /DXQX/ XPVQ(248),XPVD(248),XTVQ(248),XTVD(248), XPSQ
     +(INTMX),XPSAQ(INTMX),XTSQ(INTMX),XTSAQ(INTMX)
     *                ,XPSU(248),XTSU(248)
     *                ,XPSUT(248),XTSUT(248)
*KEEP,INTNEW.
      COMMON /INTNEW/ NVV,NSV,NVS,NSS,NDV,NVD,NDS,NSD,
     +IXPV,IXPS,IXTV,IXTS, INTVV1(248),
     +INTVV2(248),INTSV1(248),INTSV2(248), INTVS1(248),INTVS2(248),
     +INTSS1(INTMX),INTSS2(INTMX),
     +INTDV1(248),INTDV2(248),INTVD1(248),INTVD2(248),
     +INTDS1(INTMD),INTDS2(INTMD),INTSD1(INTMD),INTSD2(INTMD)
 
C  /INTNEW/
C       NVV    :   NUMBER OF INTERACTING VALENCE-VALENCE SYSTEMS
C       NSV    :   NUMBER OF INTERACTING SEA-VALENCE SYSTEMS
C       NVS    :   NUMBER OF INTERACTING VALENCE-SEA SYSTEMS
C       NSS    :   NUMBER OF INTERACTING SEA-SEA SYSTEMS
C       IXPV, IXTV : NUMBER OF GENERATED X-VALUES FOR VALENCE-QUARK
C                     SYSTEMS FROM PROJECTILE/TARGET NUCLEI
C       IXPS, IXTS : NUMBER OF GENERATED X-VALUES FOR SEA-QUARK PAIRS
C                     FROM PROJECTILE/TARGET NUCLEI
C-------------------
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,LOZUO.
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/ ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),ZUOST(INTMX),
     +INTLO(INTMX),INLOSS(INTMX)
C  /LOZUO/
C       INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C------------------
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
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
*KEEP,SHMAKL.
C     INCLUDE (SHMAKL)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON/SHMAKL/JSSH(INTMX),JTSH(INTMX),INTER1(INTMX),INTER2(INTMX)
*KEEP,NNCMS.
      COMMON /NNCMS/  GAMCM,BGCM,UMO,PCM,EPROJ,PPROJ
*KEEP,DROPPT.
      LOGICAL INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA, IPADIS,
     +ISHMAL,LPAULI
      COMMON /DROPPT/ INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA,
     +IPADIS,ISHMAL,LPAULI
*KEEP,CMHICO.
      COMMON /CMHICO/ CMHIS
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEEP,DPAR.
C     /DPAR/   CONTAINS PARTICLE PROPERTIES
C        ANAME  = LITERAL NAME OF THE PARTICLE
C        AAM    = PARTICLE MASS IN GEV
C        GA     = DECAY WIDTH
C        TAU    = LIFE TIME OF INSTABLE PARTICLES
C        IICH   = ELECTRIC CHARGE OF THE PARTICLE
C        IIBAR  = BARYON NUMBER
C        K1,K1  = BEGIN AND END OF DECAY CHANNELS OF PARTICLE
C
      CHARACTER*8  ANAME
      COMMON /DPAR/ ANAME(210),AAM(210),GA(210),TAU(210),IICH(210),
     +IIBAR(210),K1(210),K2(210)
C------------------
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
      COMMON /JSPA/PXS(40000),PYS(40000),PZS(40000),HES(40000),NNNPS
C---------------------
      COMMON /BAMCO/  NVDD
      LOGICAL LSEADI
      COMMON /SEADIQ/ LSEADI
      COMMON /MINIJ/IMINIJ,NOMJE,NOMJER,NREJEV,NOMJT,NOMJTR
      COMMON/DIQUAX/AMEDD,IDIQUA,IDIQUU
      COMMON /SECINT/ISECIN
      DATA IEVCOU/0/
C-------------
      NNNPS=0
      IEVCOU=IEVCOU+1
      NHKKH1=NHKK
      IF (IPCO.GE.1) WRITE(6,1000) NVV,NVS,NSV,NSS
 1000 FORMAT (' ENTERING HADRKK NVV,NVS,NSV,NSS '/5I5)
C----------------------------------------------------------------------
C++++++++++++++  HADRONIZE SOFT SEA-SEA CHAINS   ++++++++++++++++++++++
C---
C                                   INITIALIZE COUNTERS
      ANNVV=0.001
      ANNSS=0.001
      ANNSV=0.001
      ANNVS=0.001
      ANNCC=0.001
      ANNDV=0.001
      ANNVD=0.001
      ANNDS=0.001
      ANNSD=0.001
      ANNHH=0.001
      ANNZZ=0.001
      ANNDI=0.001
      ANNZD=0.001
      ANNDZ=0.001
      PTVV=0.
      PTSS=0.
      PTSV=0.
      PTVS=0.
      PTCC=0.
      PTDV=0.
      PTVD=0.
      PTDS=0.
      PTSD=0.
      PTHH=0.
      PTZZ=0.
      PTDI=0.
      PTZD=0.
      PTDZ=0.
      EEVV=0.
      EESS=0.
      EESV=0.
      EEVS=0.
      EECC=0.
      EEDV=0.
      EEVD=0.
      EEDS=0.
      EESD=0.
      EEHH=0.
      EEZZ=0.
      EEDI=0.
      EEZD=0.
      EEDZ=0.
C      COMMON /NCOUCH/ ACOUVV,ACOUSS,ACOUSV,ACOUVS,
C    *                 ACOUZZ,ACOUHH,ACOUDS,ACOUSD,
C    *                 ACOUDZ,ACOUZD,ACOUDI
       ACOUVV=0.
       ACOUSS=0.
       ACOUSV=0.
       ACOUVS=0.
       ACOUZZ=0.
       ACOUHH=0.
       ACOUDS=0.
       ACOUSD=0.
       ACOUDZ=0.
       ACOUZD=0.
       ACOUDI=0.
       ACOUDV=0.
       ACOUVD=0.
       ACOUCC=0.
*
      IF(IHADA.OR.IHADSS) THEN
        NVDD=0
        CALL HADRSS
      ENDIF
      IF(IHADA.OR.IHADSV) THEN
        CALL CASASV
      ENDIF
      IF(IHADA.OR.IHADVS) THEN
        CALL CASAVS
      ENDIF
        IF (IMINIJ.EQ.1) CALL HADRHH
        CALL HADRZZ
      IF(IDIQUU.EQ.1)  CALL HADRDZ
      IF(IDIQUU.EQ.1)  CALL HADRZD
C
C
C---------------------------------------------------------------
C+++++++++++++++++++    HADRONIZE sea diquark - sea CHAINS  +++++++
C
C     IF(IHADA.OR.LSEADI) THEN
      IF(IHADA) THEN
        NVDD=0
      IF(IDIQUA.EQ.1)  CALL HADRDS
      ENDIF
C
C+++++++++++++++++++    HADRONIZE sea  - sea diquark CHAINS  +++++++
C
C     IF(IHADA.OR.LSEADI) THEN
      IF(IHADA) THEN
        NVDD=0
      IF(IDIQUA.EQ.1)  CALL HADRSD
      ENDIF
C
C
C---------------------------------------------------------------
C+++++++++++++++++++    HADRONIZE SEA-VALENCE CHAINS  +++++++++++++++++
C
      IF(IHADA.OR.IHADSV) THEN
        NVDD=0
        CALL HADRSV
      ENDIF
C
C---------------------------------------------------------------
C+++++++++++++++++++    HADRONIZE sea diquark - VALENCE CHAINS  +++++++
C
C     IF(IHADA.OR.LSEADI) THEN
      IF(IHADA) THEN
        NVDD=0
      IF(IDIQUA.EQ.1)  CALL HADRDV
      ENDIF
C
C----------------------------------------------------------------------
C+++++++++++++++++++    HADRONIZE VALENCE-SEA CHAINS  +++++++++++++++++
C
      IF(IHADA.OR.IHADVS) THEN
        NVDD=0
        CALL HADRVS
      ENDIF
C
C+++++++++++++++++++    HADRONIZE valence - sea diquark CHAINS  +++++++
C
C     IF(IHADA.OR.LSEADI) THEN
      IF(IHADA) THEN
        NVDD=0
      IF(IDIQUA.EQ.1)  CALL HADRVD
      ENDIF
C
C----------------------------------------------------------------------
C                       HADRONIZE VALENCE-VALENCE CHAINS
C---
      IF(IHADA.OR.IHADVV) THEN
        NVDD=0
        CALL HADRVV
      ENDIF
C
C----------------------------------------------------------------------
C                       HADRONIZE combined (qq)-(aqaq) chains
C---
C     IF(IHADA.AND.LCOMBI) THEN
C       NVDD=15
C       CALL HADRCC
C     ENDIF
C
C---------------------------------------------------------------
C                           OPTIONAL TEST OF
C                           ENERGY-MOMENTUM CONSERVATION
C                           IN NUCLEON-NUCLEON CMS
      IF (IPCO.GE.1)THEN
        PXSU=0.
        PYSU=0.
        PZSU=0.
        ESUM=0.
        ICHSU=0
        IBASU=0
        WRITE(6,'(A)') ' HADRONS FROM HADRKK / NUCLEON-NUCLEON CMS'
        DO 10 I=NHKKH1+1,NHKK
	IF(ISTHKK(I).EQ.1)THEN
          PXSU=PXSU + PHKK(1,I)
          PYSU=PYSU + PHKK(2,I)
          PZSU=PZSU + PHKK(3,I)
          ESUM=ESUM + PHKK(4,I)
          NREF=MCIHAD(IDHKK(I))
          ICHSU=ICHSU + IICH(NREF)
          IBASU=IBASU + IIBAR(NREF)
      IF (IPCO.GE.7)
     *    WRITE(6,1010)I,(PHKK(J,I),J=1,5), IICH(NREF),IIBAR(NREF),NREF,
     +    ANAME(NREF)
 1010 FORMAT(5X,I4,5(1PE11.3),2I2,I5,A10)
	ENDIF
   10   CONTINUE
        WRITE(6,1020) PXSU,PYSU,PZSU,ESUM,ICHSU,IBASU
 1020 FORMAT(' PXSU,PYSU,PZSU,ESUM,ICHSU,IBASU'/4F10.3,2I5)
      ENDIF
C
       CALL DECHKK(NHKKH1)
C
C----------------------------------------------------------------------
C                                 LT FROM NUCLEON-NUCLEON CMS INTO LAB
C                                 PUT LAB SYSTEM PARTICLES INTO /HKKEVT/
      CMHISS=1.D0
      DO 20 I=NHKKH1+1,NHKK
        PZNN=PHKK(3,I)
        ENN =PHKK(4,I)
        IF (CMHISS.EQ.0.D0)THEN
C         PHKK(3,I) = GAMCM*PZNN + BGCM*ENN
C         PHKK(4,I) = GAMCM*ENN  + BGCM*PZNN
          PHKK(3,I) = GAMCM*PZNN + BGCM*ENN
          PHKK(4,I) = GAMCM*ENN  + BGCM*PZNN
        ENDIF
        EHECC=SQRT(PHKK(1,I)** 2+ PHKK(2,I)** 2+ PHKK(3,I)** 2+ PHKK
     +  (5,I)**2)
        IF (ABS(EHECC-PHKK(4,I)).GT.0.001D0) THEN
C            WRITE(6,'(2A/3I5,3E16.6)')
C    &         ' HADRKK: CORRECT INCONSISTENT ENERGY ',
C    *         '  IEVCOU, I,IDHKK(I), PHKK(4,I),EHECC, PHKK(5,I)',
C    *            IEVCOU, I,IDHKK(I), PHKK(4,I),EHECC, PHKK(5,I)
          PHKK(4,I)=EHECC
        ENDIF
   20 CONTINUE
C                    Secondary Interactions
      IF(ISECIN.EQ.1)CALL SEWEW(1,NHKKH1)
      KTAUAC=99
C	 IF (CMHIS.EQ.0.D0) CALL DISTR(2,NHKKH1,PPN,KTAUAC)
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                           OPTIONAL TEST OF
C                           ENERGY-MOMENTUM CONSERVATION IN LAB SYSTEM
      IF (IPCO.GE.2)THEN
        PXSU=0.
        PYSU=0.
        PZSU=0.
        ESUM=0.
        ICHSU=0
        IBASU=0
        WRITE(6,'(A)') ' HADRONS FROM HADRKK / CMS SYSTEM'
        DO 30 I=NHKKH1+1,NHKK
	  IF(ISTHKK(I).EQ.1)THEN
          PXSU=PXSU + PHKK(1,I)
          PYSU=PYSU + PHKK(2,I)
          PZSU=PZSU + PHKK(3,I)
          ESUM=ESUM + PHKK(4,I)
          NREF=MCIHAD(IDHKK(I))
          ICHSU=ICHSU + IICH(NREF)
          IBASU=IBASU + IIBAR(NREF)
      IF (IPCO.GE.7)
     *    WRITE(6,1010) I, (PHKK(J,I),J=1,5), IICH(NREF),IIBAR(NREF),
     +    NREF,ANAME(NREF)
	ENDIF
   30   CONTINUE
        WRITE(6,1020) PXSU,PYSU,PZSU,ESUM,ICHSU,IBASU
      ENDIF
C
C------------------------------------------------------------------
C
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HADRVV
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C-------------------------
C
C                       HADRONIZE VALENCE-VALENCE CHAINS
C
C                       ADD GENERATED HADRONS TO /ALLPAR/
C                          STARTING AT (NAUX + 1)
C                       AND TO /HKKEVT/ STARTING AT (NHKK + 1)
C
C-------------------------
*KEEP,INTMX.
      PARAMETER (INTMX=2488,INTMD=252)
*KEEP,DXQX.
C     INCLUDE (XQXQ)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON /DXQX/ XPVQ(248),XPVD(248),XTVQ(248),XTVD(248), XPSQ
     +(INTMX),XPSAQ(INTMX),XTSQ(INTMX),XTSAQ(INTMX)
     *                ,XPSU(248),XTSU(248)
     *                ,XPSUT(248),XTSUT(248)
*KEEP,INTNEW.
      COMMON /INTNEW/ NVV,NSV,NVS,NSS,NDV,NVD,NDS,NSD,
     +IXPV,IXPS,IXTV,IXTS, INTVV1(248),
     +INTVV2(248),INTSV1(248),INTSV2(248), INTVS1(248),INTVS2(248),
     +INTSS1(INTMX),INTSS2(INTMX),
     +INTDV1(248),INTDV2(248),INTVD1(248),INTVD2(248),
     +INTDS1(INTMD),INTDS2(INTMD),INTSD1(INTMD),INTSD2(INTMD)
 
C  /INTNEW/
C       NVV    :   NUMBER OF INTERACTING VALENCE-VALENCE SYSTEMS
C       NSV    :   NUMBER OF INTERACTING SEA-VALENCE SYSTEMS
C       NVS    :   NUMBER OF INTERACTING VALENCE-SEA SYSTEMS
C       NSS    :   NUMBER OF INTERACTING SEA-SEA SYSTEMS
C       IXPV, IXTV : NUMBER OF GENERATED X-VALUES FOR VALENCE-QUARK
C                     SYSTEMS FROM PROJECTILE/TARGET NUCLEI
C       IXPS, IXTS : NUMBER OF GENERATED X-VALUES FOR SEA-QUARK PAIRS
C                     FROM PROJECTILE/TARGET NUCLEI
C-------------------
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,LOZUO.
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/ ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),ZUOST(INTMX),
     +INTLO(INTMX),INLOSS(INTMX)
C  /LOZUO/
C       INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C------------------
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
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
*KEEP,ABRVV.
      COMMON /ABRVV/ AMCVV1(248),AMCVV2(248),GACVV1(248),GACVV2(248),
     +BGXVV1(248),BGYVV1(248),BGZVV1(248), BGXVV2(248),BGYVV2(248),
     +BGZVV2(248), NCHVV1(248),NCHVV2(248),IJCVV1(248),IJCVV2(248),
     +PQVVA1(248,4),PQVVA2(248,4), PQVVB1(248,4),PQVVB2(248,4)
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEEP,DFINPA.
      CHARACTER*8 ANF
      PARAMETER (NFIMAX=249)
      COMMON /DFINPA/ ANF(NFIMAX),PXF(NFIMAX),PYF(NFIMAX),PZF(NFIMAX),
     +HEF(NFIMAX),AMF(NFIMAX), ICHF(NFIMAX),IBARF(NFIMAX),NREF(NFIMAX)
      COMMON /DFINPZ/IORMO(NFIMAX),IDAUG1(NFIMAX),IDAUG2(NFIMAX),
     * ISTATH(NFIMAX)
*KEEP,PROJK.
      COMMON /PROJK/ IPROJK
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
       COMMON /NCOUCH/ ACOUVV,ACOUSS,ACOUSV,ACOUVS,
     *                 ACOUZZ,ACOUHH,ACOUDS,ACOUSD,
     *                 ACOUDZ,ACOUZD,ACOUDI,
     *                 ACOUDV,ACOUVD,ACOUCC
       COMMON/POPCCK/PDBCK,PDBSE,PDBSEU,
     *  IJPOCK,IREJCK,ICK4,IHAD4,ICK6,IHAD6
     *,IREJSE,ISE4,ISE6,IREJS3,ISE43,ISE63,IREJS0
     *,IHADA4,IHADA6,IREJSA,ISEA4,ISEA6,IREJA3,
     *ISEA43,ISEA63,IREJAO
C---------------------
      COMMON /ZSEA/ZSEAAV,ZSEASU,ANZSEA
C---------------------
      DIMENSION POJ(4),PAT(4)
      DATA NCALVV /0/
      IF(IPHKK.GE.6)WRITE (6,'( A)') ' hadrVV'
C-----------------------------------------------------------------
      NCALVV=NCALVV+1
      DO 50 I=1,NVV
C-----------------------drop recombined chain pairs
        IF(NCHVV1(I).EQ.99.AND.NCHVV2(I).EQ.99) GO TO 50
        IS1=INTVV1(I)
        IS2=INTVV2(I)
C
        IF (IPCO.GE.1) WRITE (6,1000) IPVQ(IS1),IPPV1(IS1),IPPV2(IS1),
     +  ITVQ(IS2),ITTV1(IS2),ITTV2(IS2), AMCVV1(I),AMCVV2(I),GACVV1(I),
     +  GACVV2(I), BGXVV1(I),BGYVV1(I),BGZVV1(I), BGXVV2(I),BGYVV2(I),
     +  BGZVV2(I), NCHVV1(I),NCHVV2(I),IJCVV1(I),IJCVV2(I), PQVVA1(I,4),
     +  PQVVA2(I,4),PQVVB1(I,4),PQVVB2(I,4)
 
 
 
 1000 FORMAT(6I5,10F9.2/10X,4I5,4F12.4)
C
C------------------------------  CHAIN 1:
C                            INCIDENT BARYONS/MESONS: QUARK-DIQUARK
C                            INCIDENT ANTIBARYONS   : AQUARK-QUARK
        IF(IBPROJ.GE.0) THEN
          IFB1=IPVQ(IS1)
          IFB2=ITTV1(IS2)
          IFB3=ITTV2(IS2)
          NOBAM=4
        ELSE
          IFB1=IPVQ(IS1)
          IFB2=ITVQ(IS2)
          IFB1=IABS(IFB1) + 6
          NOBAM=3
        ENDIF
C
        DO 10 J=1,4
          POJ(J)=PQVVA1(I,J)
          PAT(J)=PQVVA2(I,J)
   10   CONTINUE
        PT1=SQRT(POJ(1)**2+POJ(2)**2)
        PT2=SQRT(PAT(1)**2+PAT(2)**2)
        CALL PARPT(2,PT1,PT2,1,NEVT)
C------------------------------------------------------------------
C------------------------------------------------------------------
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
C               I= number of valence chain     
C               Projectile Nr ipp = IFROVP(INTVV1(I))
C               Target     Nr itt = IFROVT(INTVV2(I))
C          No of Glauber sea q at Projectile JIPP=JSSHS(IPP)
C          No of Glauber sea q at Target     JITT=JTSHS(ITT)
       IPPP = IFROVP(INTVV1(I))
       ITTT = IFROVT(INTVV2(I))
       JIPP=JSSHS(IPPP)
       JITT=JTSHS(ITTT)
C	IF(NCHVV1(I).EQ.0)THEN
C      WRITE(6,'(A,5I5)')'HADRVV: I,IPPP,ITTT,JIPP,JITT ',
C    *                     I,IPPP,ITTT,JIPP,JITT
C	ENDIF
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
      IF(IPCO.GE.1)THEN
        WRITE(6,*)' VV q-qq ,IFB1,IFB2,IFB3,',
     *	'INTVV1=IS1,INTVV2=IS2,JIPP,JITT',
     *	IFB1,IFB2,IFB3,INTVV1(I),INTVV2(I),JIPP,JITT
      ENDIF
	IF(NOBAM.EQ.3.OR.NCHVV1(I).NE.0)THEN
C       CALL HADJET(NHAD,AMCVV1(I),PAT,POJ,GACVV1(I),BGXVV1(I), BGYVV1
        CALL HADJET(NHAD,AMCVV1(I),POJ,PAT,GACVV1(I),BGXVV1(I), BGYVV1
     +  (I),BGZVV1(I),IFB1,IFB2,IFB3,IFB4, IJCVV1(I),IJCVV1(I),NOBAM,
     +  NCHVV1(I),7)
        ENDIF
	AACK=FLOAT(ICK4)/FLOAT(ICK4+IHAD4+1)
	IF((NCHVV1(I).EQ.0).AND.
     *      (NOBAM.EQ.4))THEN
        ZSEAWU=RNDM(BB)*2.D0*ZSEAAV
        RSEACK=FLOAT(JITT)*PDBSE+ ZSEAWU*PDBSEU
 	IF(IPCO.GE.1)WRITE(6,*)'HADJSE JITT,RSEACK,PDBSE 1 dpmnuc3',
     +	JITT,RSEACK,PDBSE
        IREJSS=5
        IF(RNDM(V).LE.RSEACK)THEN
	IREJSS=2
	IF(AMCVV1(I).GT.2.3D0)THEN
	  IREJSS=0
          CALL HADJSE(NHAD,AMCVV1(I),POJ,PAT,GACVV1(I),BGXVV1(I), BGYVV1
     +    (I),BGZVV1(I),IFB1,IFB2,IFB3,IFB4, IJCVV1(I),IJCVV1(I),NOBAM,
     +    NCHVV1(I),7,IREJSS,IISSQQ)
 	  IF(IPCO.GE.1)WRITE(6,'(2A,I5,F10.3,I5)')'HADJSE JITT,',
     *	  'RSEACK,IREJSS 1 dpmnuc3 ',
     +	  JITT,RSEACK,IREJSS
        ENDIF
	IF(IREJSS.GE.1)THEN
	IF(IREJSS.EQ.1)IREJSE=IREJSE+1
	IF(IREJSS.EQ.3)IREJS3=IREJS3+1
	IF(IREJSS.EQ.2)IREJS0=IREJS0+1
        CALL HADJET(NHAD,AMCVV1(I),POJ,PAT,GACVV1(I),BGXVV1(I), BGYVV1
     +  (I),BGZVV1(I),IFB1,IFB2,IFB3,IFB4, IJCVV1(I),IJCVV1(I),NOBAM,
     +  NCHVV1(I),7)
	IHAD4=IHAD4+1
	ENDIF
	IF(IREJSS.EQ.0)THEN
	  IF(IISSQQ.EQ.3)THEN
	    ISE43=ISE43+1
	  ELSE
	    ISE4=ISE4+1
	  ENDIF
	ENDIF
        ELSEIF((IJPOCK.EQ.1).AND.
     *      (AACK.LE.PDBCK))THEN
	IREJ=0
        CALL HADJCK(NHAD,AMCVV1(I),POJ,PAT,GACVV1(I),BGXVV1(I), BGYVV1
     +  (I),BGZVV1(I),IFB1,IFB2,IFB3,IFB4, IJCVV1(I),IJCVV1(I),NOBAM,
     +  NCHVV1(I),7,IREJ)
	IF(IREJ.EQ.1)THEN
	IREJCK=IREJCK+1
        CALL HADJET(NHAD,AMCVV1(I),POJ,PAT,GACVV1(I),BGXVV1(I), BGYVV1
     +  (I),BGZVV1(I),IFB1,IFB2,IFB3,IFB4, IJCVV1(I),IJCVV1(I),NOBAM,
     +  NCHVV1(I),7)
	IHAD4=IHAD4+1
	ENDIF
	IF(IREJ.EQ.0)ICK4=ICK4+1
	ELSE
        CALL HADJET(NHAD,AMCVV1(I),POJ,PAT,GACVV1(I),BGXVV1(I), BGYVV1
     +  (I),BGZVV1(I),IFB1,IFB2,IFB3,IFB4, IJCVV1(I),IJCVV1(I),NOBAM,
     +  NCHVV1(I),7)
	IHAD4=IHAD4+1
	ENDIF
	ENDIF
C------------------------------------------------------------------
C------------------------------------------------------------------
        ACOUVV=ACOUVV+1
C*** REMOVED *** 31/07/90 ***      ADD HADRONS/RESONANCES INTO
C***                               COMMON /ALLPAR/ STARTING AT NAUX
        NHKKAU=NHKK+1
        DO 20 J=1,NHAD
C
C         NHKK=NHKK+1
          IF (NHKK.EQ.NMXHKK) THEN
            WRITE (6,'(A,2I5/A)') ' HADRVV: NHKK.EQ.NMXHKK ',NHKK,NMXHKK
            RETURN
          ENDIF
C
          EHECC=SQRT(PXF(J)**2+PYF(J)**2+PZF(J)**2+AMF(J)**2)
        IF (ABS(EHECC-HEF(J)).GT.0.001D0) THEN
C           WRITE(6,'(2A/3I5,3E15.6)')
C    &            ' HADRVV / CHAIN 1 : CORRECT INCONSISTENT ENERGY ',
C    *            '  NCALVV, NHKK,NREF(J), HEF(J),EHECC, AMF(J)',
C    *            NCALVV, NHKK,NREF(J), HEF(J),EHECC, AMF(J)
          HEF(J)=EHECC
          ENDIF
          ANNVV=ANNVV+1
          EEVV=EEVV+HEF(J)
          PTVV=PTVV+SQRT(PXF(J)**2+PYF(J)**2)
C                              PUT NN-CMS HADRONS INTO /HKKEVT/
          ISTIST=1
          IF(IBARF(J).EQ.500)ISTIST=2
          CALL HKKFIL(ISTIST,MPDGHA(NREF(J)),MHKKVV(I)-3,0,
     *                PXF(J),PYF(J),PZF(J),HEF(J),NHKKAU,IORMO(J),1)
C         WRITE(6,*)' HKKFIL: NHKKAU,IORMO(J) ',ISTIST, NHKKAU,IORMO(J)
          IF(IDHKK(NHKK).EQ.99999) WRITE (6,1010) NHKK,NREF(J),IDHKK
     +    (NHKK)
 1010 FORMAT (' NHKK,NREF(J),  ',3I10)
          IMOHKK=JMOHKK(1,NHKK)
          IF(IMOHKK.LE.0.OR.IMOHKK.GT.NMXHKK)THEN
            WRITE(6,'(A,I10)')' HADRVV out of range IMOHKK= ',I10
            GO TO 2020
          ENDIF
	  IF(IREJSS.LT.0)THEN
 	  WRITE(6,*)' From HADRVV 1 first chain after HKKFIL'
          IF (IPHKK.GE.0) WRITE(6,1020) NHKK, ISTHKK(NHKK),IDHKK(NHKK),
     +    JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +    (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
          ENDIF
 1020 FORMAT (I6,I4,5I6,9E10.2)
   20   CONTINUE
C       IF(NHAD.GT.0) THEN
C         JDAHKK(1,IMOHKK)=NHKKAU
C         JDAHKK(2,IMOHKK)=NHKK
C       ENDIF
 2020   CONTINUE
C
C------------------------------    CHAIN 2
C                        INCIDENT BARYONS    :  DIQUARK-QUARK
C                        INCIDENT MESONS     :  AQUARK-QUARKC
C                        INCIDENT ANTIBARYONS:  ADIQUARK-DIQUARK
C
        IF(IBPROJ.GT.0) THEN
          IFB1=IPPV1(IS1)
          IFB2=IPPV2(IS1)
          IFB3=ITVQ(IS2)
          NOBAM=6
        ELSEIF(IBPROJ.EQ.0) THEN
          IFB1=IPPV1(IS1)
          IFB2=ITVQ(IS2)
          IFB1=IABS(IFB1) + 6
          NOBAM=3
        ELSE
          IFB1=IPPV1(IS1)
          IFB2=IPPV2(IS1)
          IFB1=IABS(IFB1) + 6
          IFB2=IABS(IFB2) + 6
          IFB3=ITTV1(IS2)
          IFB4=ITTV2(IS2)
          NOBAM=5
        ENDIF
C
        DO 30 J=1,4
          POJ(J)=PQVVB2(I,J)
          PAT(J)=PQVVB1(I,J)
   30   CONTINUE
        PT1=SQRT(POJ(1)**2+POJ(2)**2)
        PT2=SQRT(PAT(1)**2+PAT(2)**2)
        CALL PARPT(2,PT1,PT2,1,NEVT)
C***                                       POJ,PAT EXCHANGED J.R.15.2.90
C***                                       RECHANGED 19/09/90  HJM
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
C               I= number of valence chain     
C               Projectile Nr ipp = IFROVP(INTVV1(I))
C               Target     Nr itt = IFROVT(INTVV2(I))
C          No of Glauber sea q at Projectile JIPP=JSSHS(IPP)
C          No of Glauber sea q at Target     JITT=JTSHS(ITT)
       IPPP = IFROVP(INTVV1(I))
       ITTT = IFROVT(INTVV2(I))
       JIPP=JSSHS(IPPP)
       JITT=JTSHS(ITTT)
C	IF(NCHVV2(I).EQ.0)THEN
C      WRITE(6,'(A,5I5)')'HadrVV: I,IPPP,ITTT,JIPP,JITT ',
C    *                     I,IPPP,ITTT,JIPP,JITT
C	ENDIF
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
      IF(IPCO.GE.1)THEN
        WRITE(6,*)' VV qq-q ,IFB1,IFB2,IFB3,',
     *	'INTVV1=IS1,INTVV2=IS2,JIPP,JITT',
     *	IFB1,IFB2,IFB3,INTVV1(I),INTVV2(I),JIPP,JITT
      ENDIF
	IF(NOBAM.EQ.5.OR.NOBAM.EQ.3.OR.NCHVV2(I).NE.0)THEN
C       CALL HADJET(NHAD,AMCVV2(I),PAT,POJ,GACVV2(I),BGXVV2(I), BGYVV2
        CALL HADJET(NHAD,AMCVV2(I),POJ,PAT,GACVV2(I),BGXVV2(I), BGYVV2
     +  (I),BGZVV2(I),IFB1,IFB2,IFB3,IFB4, IJCVV2(I),IJCVV2(I),NOBAM,
     +  NCHVV2(I),8)
	ENDIF
	AACK=FLOAT(ICK6)/FLOAT(ICK6+IHAD6+1)
	IF((NCHVV2(I).EQ.0).AND.
     *      (NOBAM.EQ.6))THEN
        ZSEAWU=RNDM(BB)*2.D0*ZSEAAV
        RSEACK=FLOAT(JIPP)*PDBSE+ ZSEAWU*PDBSEU
 	IF(IPCO.GE.1)WRITE(6,*)'HADJSE JIPP,RSEACK,PDBSE 2 dpmnuc3',
     +	JIPP,RSEACK,PDBSE
        IREJSS=5
	IF(RNDM(V).LE.RSEACK)THEN
	IREJSS=2
	IF(AMCVV2(I).GT.2.3D0)THEN
	  IREJSS=0
          CALL HADJSE(NHAD,AMCVV2(I),POJ,PAT,GACVV2(I),BGXVV2(I), BGYVV2
     +    (I),BGZVV2(I),IFB1,IFB2,IFB3,IFB4, IJCVV2(I),IJCVV2(I),NOBAM,
     +    NCHVV2(I),8,IREJSS,IISSQQ)
 	  IF(IPCO.GE.1)WRITE(6,'(2A,I5,F10.3,I5)')'HADJSE JIPP,',
     *	  'RSEACK,IREJSS 2 dpmnux3 ',
     +	  JIPP,RSEACK,IREJSS
        ENDIF
	IF(IREJSS.GE.1)THEN
	IF(IREJSS.EQ.1)IREJSE=IREJSE+1
	IF(IREJSS.EQ.3)IREJS3=IREJS3+1
	IF(IREJSS.EQ.2)IREJS0=IREJS0+1
        CALL HADJET(NHAD,AMCVV2(I),POJ,PAT,GACVV2(I),BGXVV2(I), BGYVV2
     +  (I),BGZVV2(I),IFB1,IFB2,IFB3,IFB4, IJCVV2(I),IJCVV2(I),NOBAM,
     +  NCHVV2(I),8)
	IHAD6=IHAD6+1
	ENDIF
	IF(IREJSS.EQ.0)THEN
	  IF(IISSQQ.EQ.3)THEN
	    ISE63=ISE63+1
	  ELSE
	    ISE6=ISE6+1
	  ENDIF
	ENDIF
        ELSEIF((IJPOCK.EQ.1).AND.
     *      (AACK.LE.PDBCK))THEN
	IREJ=0
        CALL HADJCK(NHAD,AMCVV2(I),POJ,PAT,GACVV2(I),BGXVV2(I), BGYVV2
     +  (I),BGZVV2(I),IFB1,IFB2,IFB3,IFB4, IJCVV2(I),IJCVV2(I),NOBAM,
     +  NCHVV2(I),8,IREJ)
	IF(IREJ.EQ.1)THEN
	IREJCK=IREJCK+1
        CALL HADJET(NHAD,AMCVV2(I),POJ,PAT,GACVV2(I),BGXVV2(I), BGYVV2
     +  (I),BGZVV2(I),IFB1,IFB2,IFB3,IFB4, IJCVV2(I),IJCVV2(I),NOBAM,
     +  NCHVV2(I),8)
	IHAD6=IHAD6+1
	ENDIF
	IF(IREJ.EQ.0)ICK6=ICK6+1
	ELSE
        CALL HADJET(NHAD,AMCVV2(I),POJ,PAT,GACVV2(I),BGXVV2(I), BGYVV2
     +  (I),BGZVV2(I),IFB1,IFB2,IFB3,IFB4, IJCVV2(I),IJCVV2(I),NOBAM,
     +  NCHVV2(I),8)
	IHAD6=IHAD6+1
	ENDIF
	ENDIF
C                                  ADD HADRONS/RESONANCES INTO
C                                  COMMON /ALLPAR/ STARTING AT NAUX
        NHKKAU=NHKK+1
        DO 40 J=1,NHAD
C         NHKK=NHKK+1
          IF (NHKK.EQ.NMXHKK) THEN
            WRITE (6,'(A,2I5/A)') ' HADRVV: NHKK.EQ.NMXHKK ',NHKK,NMXHKK
            RETURN
          ENDIF
C
          EHECC=SQRT(PXF(J)**2+PYF(J)**2+PZF(J)**2+AMF(J)**2)
        IF (ABS(EHECC-HEF(J)).GT.0.001D0) THEN
C           WRITE(6,'(2A/3I5,3E15.6)')
C    &            ' HADRVV / CHAIN 2 : CORRECT INCONSISTENT ENERGY ',
C    *            '  NCALVV, NHKK,NREF(J), HEF(J),EHECC, AMF(J)',
C    *            NCALVV, NHKK,NREF(J), HEF(J),EHECC, AMF(J)
          HEF(J)=EHECC
          ENDIF
C                              PUT NN-CMS HADRONS INTO /HKKEVT/
          ANNVV=ANNVV+1
          EEVV=EEVV+HEF(J)
          PTVV=PTVV+SQRT(PXF(J)**2+PYF(J)**2)
          ISTIST=1
          IF(IBARF(J).EQ.500)ISTIST=2
          CALL HKKFIL(ISTIST,MPDGHA(NREF(J)),MHKKVV(I),0,
     *                PXF(J),PYF(J),PZF(J),HEF(J),NHKKAU,IORMO(J),2)
C         WRITE(6,*)' HKKFIL: NHKKAU,IORMO(J) ',ISTIST, NHKKAU,IORMO(J)
          IF(IDHKK(NHKK).EQ.99999) WRITE (6,1010)NHKK,NREF(J), IDHKK
     +    (NHKK)
          IMOHKK=JMOHKK(1,NHKK)
          IF(IMOHKK.LE.0.OR.IMOHKK.GT.NMXHKK)THEN
            WRITE(6,'(A,I10)')' HADRVV out of range IMOHKK= ',I10
            GO TO 4040
          ENDIF
	  IF(IREJSS.LT.0)THEN
 	  WRITE(6,*)' From HADRVV second chain after HKKFIL'
          IF (IPHKK.GE.0) WRITE(6,1020) NHKK, ISTHKK(NHKK),IDHKK(NHKK),
     +    JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +    (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
          ENDIF
   40   CONTINUE
C       IF(NHAD.GT.0) THEN
C         JDAHKK(1,IMOHKK)=NHKKAU
C         JDAHKK(2,IMOHKK)=NHKK
C       ENDIF
 4040   CONTINUE
   50 CONTINUE
C
C------------------------------------------------------------------
C
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HADRSV
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C-------------------------
C
C                       HADRONIZE SEA-VALENCE CHAINS
C
C                       ADD GENERATED HADRONS TO /ALLPAR/
C                          STARTING AT (NAUX + 1)
C                       AND TO /HKKEVT/ STARTING AT (NHKK + 1)
C
C---------------------------------------------------------
*KEEP,INTMX.
      PARAMETER (INTMX=2488,INTMD=252)
*KEEP,DXQX.
C     INCLUDE (XQXQ)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON /DXQX/ XPVQ(248),XPVD(248),XTVQ(248),XTVD(248), XPSQ
     +(INTMX),XPSAQ(INTMX),XTSQ(INTMX),XTSAQ(INTMX)
     *                ,XPSU(248),XTSU(248)
     *                ,XPSUT(248),XTSUT(248)
       COMMON/POPCCK/PDBCK,PDBSE,PDBSEU,
     *  IJPOCK,IREJCK,ICK4,IHAD4,ICK6,IHAD6
     *,IREJSE,ISE4,ISE6,IREJS3,ISE43,ISE63,IREJS0
     *,IHADA4,IHADA6,IREJSA,ISEA4,ISEA6,IREJA3,
     *ISEA43,ISEA63,IREJAO
*KEEP,INTNEW.
      COMMON /INTNEW/ NVV,NSV,NVS,NSS,NDV,NVD,NDS,NSD,
     +IXPV,IXPS,IXTV,IXTS, INTVV1(248),
     +INTVV2(248),INTSV1(248),INTSV2(248), INTVS1(248),INTVS2(248),
     +INTSS1(INTMX),INTSS2(INTMX),
     +INTDV1(248),INTDV2(248),INTVD1(248),INTVD2(248),
     +INTDS1(INTMD),INTDS2(INTMD),INTSD1(INTMD),INTSD2(INTMD)
 
C  /INTNEW/
C       NVV    :   NUMBER OF INTERACTING VALENCE-VALENCE SYSTEMS
C       NSV    :   NUMBER OF INTERACTING SEA-VALENCE SYSTEMS
C       NVS    :   NUMBER OF INTERACTING VALENCE-SEA SYSTEMS
C       NSS    :   NUMBER OF INTERACTING SEA-SEA SYSTEMS
C       IXPV, IXTV : NUMBER OF GENERATED X-VALUES FOR VALENCE-QUARK
C                     SYSTEMS FROM PROJECTILE/TARGET NUCLEI
C       IXPS, IXTS : NUMBER OF GENERATED X-VALUES FOR SEA-QUARK PAIRS
C                     FROM PROJECTILE/TARGET NUCLEI
C-------------------
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,LOZUO.
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/ ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),ZUOST(INTMX),
     +INTLO(INTMX),INLOSS(INTMX)
C  /LOZUO/
C       INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C------------------
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
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
*KEEP,ABRSV.
      COMMON /ABRSV/ AMCSV1(248),AMCSV2(248),GACSV1(248),GACSV2(248),
     +BGXSV1(248),BGYSV1(248),BGZSV1(248), BGXSV2(248),BGYSV2(248),
     +BGZSV2(248), NCHSV1(248),NCHSV2(248),IJCSV1(248),IJCSV2(248),
     +PQSVA1(248,4),PQSVA2(248,4), PQSVB1(248,4),PQSVB2(248,4)
*KEEP,DFINPA.
      CHARACTER*8 ANF
      PARAMETER (NFIMAX=249)
      COMMON /DFINPA/ ANF(NFIMAX),PXF(NFIMAX),PYF(NFIMAX),PZF(NFIMAX),
     +HEF(NFIMAX),AMF(NFIMAX), ICHF(NFIMAX),IBARF(NFIMAX),NREF(NFIMAX)
      COMMON /DFINPZ/IORMO(NFIMAX),IDAUG1(NFIMAX),IDAUG2(NFIMAX),
     * ISTATH(NFIMAX)
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEEP,PROJK.
      COMMON /PROJK/ IPROJK
*KEEP,NUCC.
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
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
       COMMON /NCOUCH/ ACOUVV,ACOUSS,ACOUSV,ACOUVS,
     *                 ACOUZZ,ACOUHH,ACOUDS,ACOUSD,
     *                 ACOUDZ,ACOUZD,ACOUDI,
     *                 ACOUDV,ACOUVD,ACOUCC
C---------------------
      COMMON /ZSEA/ZSEAAV,ZSEASU,ANZSEA
      COMMON /CASADI/CASAXX,ICASAD
C---------------------
      DIMENSION POJ(4),PAT(4)
      DATA NCALSV /0/
C-----------------------------------------------------------------------
      NCALSV=NCALSV+1
      DO 50 I=1,NSV
C-----------------------drop recombined chain pairs
        IF(NCHSV1(I).EQ.99.AND.NCHSV2(I).EQ.99) GO TO 50
        IS1=INTSV1(I)
        IS2=INTSV2(I)
C
        IF (IPCO.GE.6) WRITE (6,1000) IPSQ(IS1),IPSAQ(IS1),ITVQ(IS2),
     +  ITTV1(IS2),ITTV2(IS2), AMCSV1(I),AMCSV2(I),GACSV1(I),GACSV2(I),
     +  BGXSV1(I),BGYSV1(I),BGZSV1(I), BGXSV2(I),BGYSV2(I),BGZSV2(I),
     +  NCHSV1(I),NCHSV2(I),IJCSV1(I),IJCSV2(I), PQSVA1(I,4),PQSVA2
     +  (I,4),PQSVB1(I,4),PQSVB2(I,4)
 1000 FORMAT(10X,5I5,10F9.2/10X,4I5,4F12.4)
C
C++++++++++++++++++++++++++++++    CHAIN 1:  QUARK-DIQUARK   +++++++++++
        IFB1=IPSQ(IS1)
        IFB2=ITTV1(IS2)
        IFB3=ITTV2(IS2)
        DO 10 J=1,4
          POJ(J)=PQSVA1(I,J)
          PAT(J)=PQSVA2(I,J)
   10   CONTINUE
        PT1=SQRT(POJ(1)**2+POJ(2)**2)
        PT2=SQRT(PAT(1)**2+PAT(2)**2)
        CALL PARPT(2,PT1,PT2,3,NEVT)
C       IF((NCHSV1(I).NE.0.OR.NCHSV2(I).NE.0).AND.IP.NE.1)
C    &  CALL SAPTRE(AMCSV1(I),GACSV1(I),BGXSV1(I),BGYSV1(I),BGZSV1(I),
C    &              AMCSV2(I),GACSV2(I),BGXSV2(I),BGYSV2(I),BGZSV2(I))
C----------------------------------------------------------------
        IF (IPCO.GE.6)WRITE (6,1244) POJ,PAT
 1244   FORMAT ('  S-V QUARK-DIQUARK POJ,PAT ',8E12.3)
C------------------------------------------------------------------
C------------------------------------------------------------------
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
C               I= number of valence chain     
C               Target     Nr itt = IFROVT(INTSV2(I))
C          No of Glauber sea q at Target     JITT=JTSHS(ITT)
       ITTT = IFROVT(INTSV2(I))
       JITT=JTSHS(ITTT)
C	IF(NCHSV1(I).EQ.0)THEN
C      WRITE(6,'(A,3I5)')'HADRSV: I,ITTT,JITT ',
C    *                     I,ITTT,JITT
C	ENDIF
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
      IF(IPCO.GE.1)THEN
        WRITE(6,*)' SV q-qq ,IFB1,IFB2,IFB3,',
     *	'INTSV1=IS1,INTSV2=IS2,JIPPX,JITT',
     *	IFB1,IFB2,IFB3,INTSV1(I),INTSV2(I),JIPPX,JITT
       ENDIF
C------------------------------------------------------------------- 
C-------------------------------------------------------------------    
	IF((NCHSV1(I).NE.0))THEN
        CALL HADJET(NHAD,AMCSV1(I),POJ,PAT,GACSV1(I),BGXSV1(I), BGYSV1
     +  (I),BGZSV1(I),IFB1,IFB2,IFB3,IFB4, IJCSV1(I),IJCSV1(I),4,NCHSV1
     +  (I),3)
	ENDIF
	AACK=FLOAT(ICK4)/FLOAT(ICK4+IHAD4+1)
	IF((NCHSV1(I).EQ.0))THEN
        ZSEAWU=RNDM(BB)*2.D0*ZSEAAV
        RSEACK=FLOAT(JITT)*PDBSE+ ZSEAWU*PDBSEU
 	IF(IPCO.GE.1)WRITE(6,*)'HADJSE JITT,RSEACK,PDBSE 3 dpmnuc3',
     +	JITT,RSEACK,PDBSE
        IREJSS=5
        IF(RNDM(V).LE.RSEACK)THEN
	IREJSS=2
	IF(AMCSV1(I).GT.2.3D0)THEN
          IREJSS=0
          CALL HADJSE(NHAD,AMCSV1(I),POJ,PAT,GACSV1(I),BGXSV1(I),
     *	  BGYSV1
     +    (I),BGZSV1(I),IFB1,IFB2,IFB3,IFB4, IJCSV1(I),IJCSV1(I),4,
     *    NCHSV1
     +    (I),3,IREJSS,IISSQQ)
 	  IF(IPCO.GE.1)WRITE(6,'(2A,I5,F10.3,I5)')'HADJSE JITT,',
     *	  'RSEACK,IREJSS 3 dpmnuc3 ',
     +	  JITT,RSEACK,IREJSS
        ENDIF
	IF(IREJSS.GE.1)THEN
	IF(IREJSS.EQ.1)IREJSE=IREJSE+1
	IF(IREJSS.EQ.3)IREJS3=IREJS3+1
	IF(IREJSS.EQ.2)IREJS0=IREJS0+1
        CALL HADJET(NHAD,AMCSV1(I),POJ,PAT,GACSV1(I),BGXSV1(I), BGYSV1
     +  (I),BGZSV1(I),IFB1,IFB2,IFB3,IFB4, IJCSV1(I),IJCSV1(I),4,NCHSV1
     +  (I),3)
	IHAD4=IHAD4+1
	ENDIF
        IF(IREJSS.EQ.0)THEN
	  IF(IISSQQ.EQ.3)THEN
	    ISE43=ISE43+1
	  ELSE
	    ISE4=ISE4+1
	  ENDIF  
	ENDIF    
        ELSEIF((IJPOCK.EQ.1).AND.
     *      (AACK.LE.PDBCK))THEN
	IREJ=0
        CALL HADJCK(NHAD,AMCSV1(I),POJ,PAT,GACSV1(I),BGXSV1(I), BGYSV1
     +  (I),BGZSV1(I),IFB1,IFB2,IFB3,IFB4, IJCSV1(I),IJCSV1(I),4,NCHSV1
     +  (I),3,IREJ)
	IF(IREJ.EQ.1)THEN
	IREJCK=IREJCK+1
        CALL HADJET(NHAD,AMCSV1(I),POJ,PAT,GACSV1(I),BGXSV1(I), BGYSV1
     +  (I),BGZSV1(I),IFB1,IFB2,IFB3,IFB4, IJCSV1(I),IJCSV1(I),4,NCHSV1
     +  (I),3)
	IHAD4=IHAD4+1
	ENDIF
        IF(IREJ.EQ.0)ICK4=ICK4+1
	ELSE
        CALL HADJET(NHAD,AMCSV1(I),POJ,PAT,GACSV1(I),BGXSV1(I), BGYSV1
     +  (I),BGZSV1(I),IFB1,IFB2,IFB3,IFB4, IJCSV1(I),IJCSV1(I),4,NCHSV1
     +  (I),3)
	IHAD4=IHAD4+1
	ENDIF
	ENDIF
C------------------------------------------------------------------
C------------------------------------------------------------------
        ACOUSV=ACOUSV+1
C*** REMOVED 31/07/90 HJM ***      ADD HADRONS/RESONANCES INTO
C                                  COMMON /ALLPAR/ STARTING AT NAUX
        NHKKAU=NHKK+1
        PIXU=0.
        PIYU=0.
        PIZU=0.
        PIEU=0.
        DO 20 J=1,NHAD
C         NHKK=NHKK+1
          IF (NHKK.EQ.NMXHKK) THEN
            WRITE (6,'(A,2I5/A)') ' HADRSV: NHKK.EQ.NMXHKK ',NHKK,NMXHKK
            RETURN
          ENDIF
C
          EHECC=SQRT(PXF(J)**2+PYF(J)**2+PZF(J)**2+AMF(J)**2)
        IF (ABS(EHECC-HEF(J)).GT.0.001D0.AND.IBARF(J).NE.500) THEN
            WRITE(6,'(2A/3I5,3E15.6)')
     &            ' HADRSV / CHAIN 1 : CORRECT INCONSISTENT ENERGY ',
     *            '  NCALSV, NHKK,NREF(J), HEF(J),EHECC, AMF(J)',
     *            NCALSV, NHKK,NREF(J), HEF(J),EHECC, AMF(J)
          HEF(J)=EHECC
          ENDIF
          ANNSV=ANNSV+1
          EESV=EESV+HEF(J)
          PTSV=PTSV+SQRT(PXF(J)**2+PYF(J)**2)
C                               PUT NN-CMS HADRONS INTO /HKKEVT/
          ISTIST=1
          IF(IBARF(J).EQ.500)ISTIST=2
          CALL HKKFIL(ISTIST,MPDGHA(NREF(J)),MHKKSV(I)-3,0,
     *                PXF(J),PYF(J),PZF(J),HEF(J),NHKKAU,IORMO(J),3)
C         WRITE(6,*)' HKKFIL: NHKKAU,IORMO(J) ',ISTIST, NHKKAU,IORMO(J)
          IF(IDHKK(NHKK).EQ.99999) WRITE (6,1030)NHKK,NREF(J), IDHKK
     +    (NHKK)
          PIXU=PIXU+PXF(J)
          PIYU=PIYU+PYF(J)
          PIZU=PIZU+PZF(J)
          PIEU=PIEU+HEF(J)
	  IF(IREJSS.LT.0)THEN
	  WRITE(6,*)' HADRSV / CHAIN 1'
          IF (IPHKK.GE.0) WRITE(6,1010) NHKK, ISTHKK(NHKK),IDHKK(NHKK),
     +    JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +    (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
          ENDIF
   20   CONTINUE
        IF(IPCO.GE.6)WRITE(6,1644)PIXU,PIYU,PIZU,PIEU
 1644   FORMAT(' HADRSV,ch1 PIXU,PIYU,PIZU,PIEU ',4F12.5)
C       IF(NHAD.GT.0) THEN
C         JDAHKK(1,IMOHKK)=NHKKAU
C         JDAHKK(2,IMOHKK)=NHKK
C       ENDIF
C+++++++++++++++++++++++++++++   CHAIN 2:  AQUARK-QUARK  ++++++++++++++
        IFB1=IPSAQ(IS1)
        IFB2=ITVQ(IS2)
        IFB1=IABS(IFB1)+6
        DO 30 J=1,4
          POJ(J)=PQSVB2(I,J)
          PAT(J)=PQSVB1(I,J)
   30   CONTINUE
        PT1=SQRT(POJ(1)**2+POJ(2)**2)
        PT2=SQRT(PAT(1)**2+PAT(2)**2)
        CALL PARPT(2,PT1,PT2,3,NEVT)
C
      IF(IPCO.GE.1)THEN
        WRITE(6,*)' SV aq-q ,IFB1,IFB2,',
     *	'INTSV1=IS1,INTSV2=IS2,JIPPX,JITTX',
     *	IFB1,IFB2,INTSV1(I),INTSV2(I),JIPPX,JITTX
      ENDIF
C------------------------------------------------------------------- 
C-------------------------------------------------------------------    
        IF (IPCO.GE.6)WRITE (6,1244) POJ,PAT
        CALL HADJET(NHAD,AMCSV2(I),POJ,PAT,GACSV2(I),BGXSV2(I), BGYSV2
     +  (I),BGZSV2(I),IFB1,IFB2,IFB3,IFB4, IJCSV2(I),IJCSV2(I),3,NCHSV2
     +  (I),4)
C                                   ADD HADRONS/RESONANCES INTO
C                                   COMMON /ALLPAR/ STARTING AT NAUX
        NHKKAU=NHKK+1
        PIXU=0.
        PIYU=0.
        PIZU=0.
        PIEU=0.
        DO 40 J=1,NHAD
          IF (NHKK.EQ.NMXHKK) THEN
            WRITE (6,'(A,2I5/A)') ' HADRSV: NHKK.EQ.NMXHKK ', NHKK,
     +      NMXHKK
            RETURN
          ENDIF
C         NHKK=NHKK+1
C
          EHECC=SQRT(PXF(J)**2+PYF(J)**2+PZF(J)**2+AMF(J)**2)
        IF (ABS(EHECC-HEF(J)).GT.0.001D0.AND.IBARF(J).NE.500) THEN
              WRITE(6,'(2A/3I5,3E15.6)')
     &            ' HADRSV / CHAIN 2 : CORRECT INCONSISTENT ENERGY ',
     *            '  NCALSV, NHKK,NREF(J), HEF(J),EHECC, AMF(J)',
     *            NCALSV, NHKK,NREF(J), HEF(J),EHECC, AMF(J)
          HEF(J)=EHECC
          ENDIF
          ANNSV=ANNSV+1
          EESV=EESV+HEF(J)
          PTSV=PTSV+SQRT(PXF(J)**2+PYF(J)**2)
C                               PUT NN-CMS HADRONS INTO /HKKEVT/
          ISTIST=1
          IF(IBARF(J).EQ.500)ISTIST=2
          CALL HKKFIL(ISTIST,MPDGHA(NREF(J)),MHKKSV(I),0,
     *                PXF(J),PYF(J),PZF(J),HEF(J),NHKKAU,IORMO(J),4)
          IF(IDHKK(NHKK).EQ.99999) WRITE (6,1030)NHKK,NREF(J), IDHKK
     +    (NHKK)
          PIXU=PIXU+PXF(J)
          PIYU=PIYU+PYF(J)
          PIZU=PIZU+PZF(J)
          PIEU=PIEU+HEF(J)
          IF (IPHKK.GE.7) WRITE(6,1010) NHKK, ISTHKK(NHKK),IDHKK(NHKK),
     +    JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +    (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
   40   CONTINUE
        IF(IPCO.GE.6)WRITE(6,1644)PIXU,PIYU,PIZU,PIEU
C       IF(NHAD.GT.0) THEN
C         JDAHKK(1,IMOHKK)=NHKKAU
C         JDAHKK(2,IMOHKK)=NHKK
C       ENDIF
   50 CONTINUE
C----------------------------------------------------------------
C
      RETURN
 1010 FORMAT (I6,I4,5I6,9E10.2)
 1020 FORMAT (' HADRKK J.GT.NAUMAX SKIP NEXT PARTICLES ',3I10)
 1030 FORMAT (' NHKK,IDHKK(NHKK)  ',3I10)
      END
*-- Author :
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HADRSS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
*KEEP,INTMX.
      PARAMETER (INTMX=2488,INTMD=252)
*KEEP,DXQX.
C     INCLUDE (XQXQ)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON /DXQX/ XPVQ(248),XPVD(248),XTVQ(248),XTVD(248), XPSQ
     +(INTMX),XPSAQ(INTMX),XTSQ(INTMX),XTSAQ(INTMX)
     *                ,XPSU(248),XTSU(248)
     *                ,XPSUT(248),XTSUT(248)
*KEEP,INTNEW.
      COMMON /INTNEW/ NVV,NSV,NVS,NSS,NDV,NVD,NDS,NSD,
     +IXPV,IXPS,IXTV,IXTS, INTVV1(248),
     +INTVV2(248),INTSV1(248),INTSV2(248), INTVS1(248),INTVS2(248),
     +INTSS1(INTMX),INTSS2(INTMX),
     +INTDV1(248),INTDV2(248),INTVD1(248),INTVD2(248),
     +INTDS1(INTMD),INTDS2(INTMD),INTSD1(INTMD),INTSD2(INTMD)
 
C  /INTNEW/
C       NVV    :   NUMBER OF INTERACTING VALENCE-VALENCE SYSTEMS
C       NSV    :   NUMBER OF INTERACTING SEA-VALENCE SYSTEMS
C       NVS    :   NUMBER OF INTERACTING VALENCE-SEA SYSTEMS
C       NSS    :   NUMBER OF INTERACTING SEA-SEA SYSTEMS
C       IXPV, IXTV : NUMBER OF GENERATED X-VALUES FOR VALENCE-QUARK
C                     SYSTEMS FROM PROJECTILE/TARGET NUCLEI
C       IXPS, IXTS : NUMBER OF GENERATED X-VALUES FOR SEA-QUARK PAIRS
C                     FROM PROJECTILE/TARGET NUCLEI
C-------------------
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,LOZUO.
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/ ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),ZUOST(INTMX),
     +INTLO(INTMX),INLOSS(INTMX)
C  /LOZUO/
C       INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C------------------
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
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
*KEEP,ABRSS.
C     INCLUDE (ABRSS)
      COMMON /ABRSS/ AMCSS1(INTMX),AMCSS2(INTMX), GACSS1(INTMX),GACSS2
     +(INTMX), BGXSS1(INTMX),BGYSS1(INTMX),BGZSS1(INTMX), BGXSS2(INTMX),
     +BGYSS2(INTMX),BGZSS2(INTMX), NCHSS1(INTMX),NCHSS2(INTMX), IJCSS1
     +(INTMX),IJCSS2(INTMX), PQSSA1(INTMX,4),PQSSA2(INTMX,4), PQSSB1
     +(INTMX,4),PQSSB2(INTMX,4)
*KEEP,NUCC.
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
*KEEP,DFINPA.
      CHARACTER*8 ANF
      PARAMETER (NFIMAX=249)
      COMMON /DFINPA/ ANF(NFIMAX),PXF(NFIMAX),PYF(NFIMAX),PZF(NFIMAX),
     +HEF(NFIMAX),AMF(NFIMAX), ICHF(NFIMAX),IBARF(NFIMAX),NREF(NFIMAX)
      COMMON /DFINPZ/IORMO(NFIMAX),IDAUG1(NFIMAX),IDAUG2(NFIMAX),
     * ISTATH(NFIMAX)
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEEP,PROJK.
      COMMON /PROJK/ IPROJK
*KEND.
      DIMENSION POJ(4),PAT(4)
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
       COMMON /NCOUCH/ ACOUVV,ACOUSS,ACOUSV,ACOUVS,
     *                 ACOUZZ,ACOUHH,ACOUDS,ACOUSD,
     *                 ACOUDZ,ACOUZD,ACOUDI,
     *                 ACOUDV,ACOUVD,ACOUCC
C---------------------
      COMMON /PSHOW/ IPSHOW
C     COMMON /HARLUN/ IHARLU,QLUN
      COMMON /HARLUN/ QLUN,IHARLU
      COMMON /JSPAR/PXJ(1000),PYJ(1000),PZJ(1000),HEJ(1000),NNNPJ
      COMMON /JSPA/PXS(40000),PYS(40000),PZS(40000),HES(40000),NNNPS
      COMMON /NOMIJE/ PTMIJE(10),NNMIJE(10)
      COMMON /CASADI/CASAXX,ICASAD
C-----------------------------------------------------------------------
      DO 60 I=1,NSS
C-----------------------drop recombined chain pairs
        IF(NCHSS1(I).EQ.99.AND.NCHSS2(I).EQ.99) GO TO 60
        IF (INLOSS(I)) THEN
          IS1=INTSS1(I)
          IS2=INTSS2(I)
C
          IF (IPCO.GE.6) WRITE (6,1000) IPSQ(IS1),IPSAQ(IS1),ITSQ(IS2),
     +    ITSAQ(IS2), AMCSS1(I),AMCSS2(I),GACSS1(I),GACSS2(I), BGXSS1
     +    (I),BGYSS1(I),BGZSS1(I), BGXSS2(I),BGYSS2(I),BGZSS2(I), NCHSS1
     +    (I),NCHSS2(I),IJCSS1(I),IJCSS2(I), PQSSA1(I,4),PQSSA2(I,4),
     +    PQSSB1(I,4),PQSSB2(I,4)
 1000 FORMAT(10X,4I5,10F9.2/10X,4I5,4F12.4)
C
C+++++++++++++++++++++++++++++     CHAIN 1:  QUARK-AQUARK    ++++++++++
          IFB1=IPSQ(IS1)
          IFB2=ITSAQ(IS2)
          IFB2=IABS(IFB2)+6
          DO 10 J=1,4
            POJ(J)=PQSSA1(I,J)
            PAT(J)=PQSSA2(I,J)
   10     CONTINUE
        PT1=SQRT(POJ(1)**2+POJ(2)**2)
        PT2=SQRT(PAT(1)**2+PAT(2)**2)
        CALL PARPT(2,PT1,PT2,4,NEVT)
C--------------------------------------------------------------
            IHARLU=0
            QLUN=0.
            IF(IPSHOW.EQ.1)THEN
              POJPT=SQRT(POJ(2)**2+POJ(1)**2)
              PATPT=SQRT(PAT(1)**2+PAT(2)**2)
	      DO IIII=1,10
              IF(POJPT.GE.PTMIJE(IIII))NNMIJE(IIII)=
     *	      NNMIJE(IIII)+1
              IF(PATPT.GE.PTMIJE(IIII))NNMIJE(IIII)=
     *	      NNMIJE(IIII)+1
	      ENDDO
              QLUN=MIN(POJPT,PATPT)
              IF((QLUN.LT.2.5D0).OR.(AMCSS1(I).LT.5.D0))THEN
                QLUN=0.
                IHARLU=0
              ELSE 
                IHARLU=1
              ENDIF   
            ENDIF
      IF(IPCO.GE.1)THEN
        WRITE(6,*)' SS q-aq ,IFB1,IFB2,',
     *	'INTSS1=IS1,INTSS2=IS2',
     *	IFB1,IFB2,INTSS1(I),INTSS2(I)
        WRITE (6,*)' projectile sea quark IFB1=',IFB1,
     *	' from IS1=',INTSS1(I)
        WRITE(6,*)' with IPSQ(IS1),XPSQ(IS1),IFROSP(IS1)',
     *	IPSQ(IS1),XPSQ(IS1),IFROSP(IS1)
      ENDIF
        DO 798 II=1,IXPV
	  IF(IFROSP(IS1).EQ.IFROVP(II))III=II
  798   CONTINUE	
      IF(IPCO.GE.1)THEN
        WRITE (6,*)' projectile III=',III
        WRITE(6,*)' corresp. XPVQ(i),XPVD(i),IPVQ(I),IPPV1(I),IPPV2(I)',
     *   XPVQ(III),XPVD(III),IPVQ(III),IPPV1(III),IPPV2(III)
      ENDIF
C------------------------------------------------------------------- 
C                         Casado diquark option
C+++++++++++++++++++++++++++ SS    CHAIN 1:  QUARK-AQUARK    ++++++++++
C-------------------------------------------------------------------    
       IF(ICASAD.EQ.1)THEN
         IF(RNDM(VV).LE.CASAXX)THEN
	   IF(RNDM(VVV).LE.0.5D0)THEN
	     ISCASA=IPSQ(IS1) 
	     IPVCAS=IPPV1(III)
	     IPSQ(IS1)=IPVCAS
	     IPPV1(III)=ISCASA
	     IFB1=IPSQ(IS1)
      IF(IPCO.GE.1)THEN
        WRITE(6,*)' Cas SS1 q-aq 1 ,IFB1,IFB2,',
     *	'INTSS1=IS1,INTSS2=IS2,III',
     *	IFB1,IFB2,INTSS1(I),INTSS2(I),III
     *  ,'-----------------------------------------------------'
      ENDIF
	   ELSE
	     ISCASA=IPSQ(IS1) 
	     IPVCAS=IPPV2(III)
	     IPSQ(IS1)=IPVCAS
	     IPPV2(III)=ISCASA
	     IFB1=IPSQ(IS1)
      IF(IPCO.GE.1)THEN
        WRITE(6,*)' Cas SS1 q-aq 2 ,IFB1,IFB2,',
     *	'INTSS1=IS1,INTSS2=IS2,III',
     *	IFB1,IFB2,INTSS1(I),INTSS2(I),III
     *  ,'-----------------------------------------------------'
      ENDIF
	   ENDIF
	 ENDIF
       ENDIF
C------------------------------------------------------------------- 
C                         Casado diquark option
C-------------------------------------------------------------------    
          CALL HADJET(NHAD,AMCSS1(I),POJ,PAT,GACSS1(I),BGXSS1(I), BGYSS1
     +    (I),BGZSS1(I),IFB1,IFB2,IFB3,IFB4, IJCSS1(I),IJCSS1(I),3,
     +    NCHSS1(I),1)
          ACOUSS=ACOUSS+1
            IHARLU=0
            QLUN=0.
C                                  ADD HADRONS/RESONANCES INTO
C                                  COMMON /ALLPAR/ STARTING AT NAUX
          NHKKAU=NHKK+1
          DO 20 J=1,NHAD
            EHECC=SQRT(PXF(J)**2+PYF(J)**2+PZF(J)**2+AMF(J)**2)
          IF (ABS(EHECC-HEF(J)).GT.0.001D0.AND.IBARF(J).NE.500) THEN
              WRITE(6,'(A,2I5,2E16.6)')
     +        ' HADRSS: CORRECT INCONSISTENT PARTICLE ENERGY ', NHKK,
     +        NREF(J), HEF(J),EHECC
            HEF(J)=EHECC
            ENDIF
          ANNSS=ANNSS+1
          EESS=EESS+HEF(J)
          PTSS=PTSS+SQRT(PXF(J)**2+PYF(J)**2)
C                               PUT NN-CMS HADRONS INTO /HKKEVT/
C           NHKK=NHKK+1
            IF (NHKK.EQ.NMXHKK) THEN
              WRITE (6,'(A,2I5)') ' HADRSS: NHKK.EQ NMXHKK',NHKK,NMXHKK
              RETURN
            ENDIF
          ISTIST=1
          IF(IBARF(J).EQ.500)ISTIST=2
          CALL HKKFIL(ISTIST,MPDGHA(NREF(J)),MHKKSS(I)-3,0,
     *                PXF(J),PYF(J),PZF(J),HEF(J),NHKKAU,IORMO(J),5)
            IF(IDHKK(NHKK).EQ.99999) WRITE (6,1030) NHKK,NREF(J), IDHKK
     +      (NHKK)
C	    WRITE(6,*)' First chain HADRSS'
            IF (IPHKK.GE.7) WRITE(6,1010) NHKK, ISTHKK(NHKK),IDHKK
     +      (NHKK),JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK
     +      (2,NHKK),(PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK
     +      =1,4)
   20     CONTINUE
   30     CONTINUE
C       IF(NHAD.GT.0) THEN
C         JDAHKK(1,IMOHKK)=NHKKAU
C         JDAHKK(2,IMOHKK)=NHKK
C       ENDIF
          IF(NNNPJ.GE.1)THEN
            NNNPSO=NNNPS
            NNNPS=NNNPS+1
            NNNPSU=NNNPSO+NNNPJ
            DO 137 J=NNNPS,NNNPSU
              JJ=J-NNNPS+1
              IF(J.GT.40000.OR.JJ.GT.1000)THEN
C               WRITE(6,'(A,2I10)')' J.gt.40000.or.jj.gt.1000 ',J,JJ
                GO TO 137
              ENDIF
	      PXS(J)=PXJ(JJ)
	      PYS(J)=PYJ(JJ)
	      PZS(J)=PZJ(JJ)
	      HES(J)=HEJ(JJ)
  137      CONTINUE
            NNNPS=NNNPS+NNNPJ-1
          ENDIF
C
C++++++++++++++++++++++++++++++++   CHAIN 2:  AQUARK-QUARK   +++++++++
          IFB1=IPSAQ(IS1)
          IFB2=ITSQ(IS2)
          IFB1=IABS(IFB1)+6
          DO 40 J=1,4
            POJ(J)=PQSSB2(I,J)
            PAT(J)=PQSSB1(I,J)
   40     CONTINUE
        PT1=SQRT(POJ(1)**2+POJ(2)**2)
        PT2=SQRT(PAT(1)**2+PAT(2)**2)
        CALL PARPT(2,PT1,PT2,4,NEVT)
            IHARLU=0
            QLUN=0.
            IF(IPSHOW.EQ.1)THEN
              POJPT=SQRT(POJ(2)**2+POJ(1)**2)
              PATPT=SQRT(PAT(1)**2+PAT(2)**2)
	      DO IIII=1,10
              IF(POJPT.GE.PTMIJE(IIII))NNMIJE(IIII)=
     *	      NNMIJE(IIII)+1
              IF(PATPT.GE.PTMIJE(IIII))NNMIJE(IIII)=
     *	      NNMIJE(IIII)+1
	      ENDDO
              QLUN=MIN(POJPT,PATPT)
              IF((QLUN.LT.2.5D0).OR.(AMCSS2(I).LT.5.D0))THEN
                QLUN=0.
                IHARLU=0
              ELSE 
                IHARLU=1
              ENDIF   
            ENDIF
C,,
      IF(IPCO.GE.1)THEN
        WRITE(6,*)' SS aq-q ,IFB1,IFB2,',
     *	'INTSS1=IS1,INTSS2=IS2',
     *	IFB1,IFB2,INTSS1(I),INTSS2(I)
        WRITE (6,*)' target sea quark IFB2=',IFB2,
     *	' from IS2=',INTSS2(I)
        WRITE(6,*)' with ITSQ(IS2),XTSQ(IS2),IFROST(IS2)',
     *	ITSQ(IS2),XTSQ(IS2),IFROST(IS2)
      ENDIF
        DO 797 II=1,IXTV
	  IF(IFROST(IS2).EQ.IFROVT(II))III=II
  797   CONTINUE	
      IF(IPCO.GE.1)THEN
        WRITE (6,*)' projectile III=',III
        WRITE(6,*)' corresp. XTVQ(i),XTVD(i),ITVQ(I),ITTV1(I),ITTV2(I)',
     *   XTVQ(III),XTVD(III),ITVQ(III),ITTV1(III),ITTV2(III)
      ENDIF
C------------------------------------------------------------------- 
C                         Casado diquark option
C+++++++++++++++++++++++++++++ SS   CHAIN 2:  AQUARK-QUARK   +++++++++
C-------------------------------------------------------------------    
       IF(ICASAD.EQ.1)THEN
         IF(RNDM(VV).LE.CASAXX)THEN
	   IF(RNDM(VVV).LE.0.5D0)THEN
	     ISCASA=ITSQ(IS2) 
	     ITVCAS=ITTV1(III)
	     ITSQ(IS2)=ITVCAS
	     ITTV1(III)=ISCASA
	     IFB2=ITSQ(IS2)
      IF(IPCO.GE.1)THEN
        WRITE(6,*)' Cas SS2 aq-q 1 ,IFB1,IFB2,',
     *	'INTSS1=IS1,INTSS2=IS2,III',
     *	IFB1,IFB2,INTSS1(I),INTSS2(I),III
     *  ,'-----------------------------------------------------'
      ENDIF
	   ELSE
	     ISCASA=ITSQ(IS2) 
	     ITVCAS=ITTV2(III)
	     ITSQ(IS2)=ITVCAS
	     ITTV2(III)=ISCASA
	     IFB2=ITSQ(IS2)
      IF(IPCO.GE.1)THEN
        WRITE(6,*)' Cas SS2 aq-q 2 ,IFB1,IFB2,',
     *	'INTSS1=IS1,INTSS2=IS2,III',
     *	IFB1,IFB2,INTSS1(I),INTSS2(I),III
     *  ,'-----------------------------------------------------'
      ENDIF
	   ENDIF
	 ENDIF
       ENDIF
C------------------------------------------------------------------- 
C                         Casado diquark option
C-------------------------------------------------------------------    
          CALL HADJET(NHAD,AMCSS2(I),POJ,PAT,GACSS2(I),BGXSS2(I), BGYSS2
     +    (I),BGZSS2(I),IFB1,IFB2,IFB3,IFB4, IJCSS2(I),IJCSS2(I),3,
     +    NCHSS2(I),2)
            IHARLU=0
            QLUN=0.
C                                    ADD HADRONS/RESONANCES INTO
C                                    COMMON /ALLPAR/ STARTING AT NAUX
          NHKKAU=NHKK+1
          DO 50 J=1,NHAD
            EHECC=SQRT(PXF(J)**2+PYF(J)**2+PZF(J)**2+AMF(J)**2)
          IF (ABS(EHECC-HEF(J)).GT.0.001D0.AND.IBARF(J).NE.500) THEN
              WRITE(6,'(A,2I5,2E16.6)')
     +        ' HADRSS: CORRECT INCONSISTENT PARTICLE ENERGY ', NHKK,
     +        NREF(J), HEF(J),EHECC
            HEF(J)=EHECC
            ENDIF
          ANNSS=ANNSS+1
          EESS=EESS+HEF(J)
          PTSS=PTSS+SQRT(PXF(J)**2+PYF(J)**2)
C                              PUT NN-CMS HADRONS INTO /HKKEVT/
C           NHKK=NHKK+1
            IF (NHKK.EQ.NMXHKK) THEN
              WRITE (6,'(A,2I5)') ' HADRSS: NHKK.EQ NMXHKK',NHKK,NMXHKK
              RETURN
            ENDIF
          ISTIST=1
          IF(IBARF(J).EQ.500)ISTIST=2
          CALL HKKFIL(ISTIST,MPDGHA(NREF(J)),MHKKSS(I),0,
     *                PXF(J),PYF(J),PZF(J),HEF(J),NHKKAU,IORMO(J),6)
            IF(IDHKK(NHKK).EQ.99999) WRITE (6,1030) NHKK,NREF(J), IDHKK
     +      (NHKK)
C	    WRITE(6,*)' Second chain HADRSS'
            IF (IPHKK.GE.7) WRITE(6,1010) NHKK, ISTHKK(NHKK),IDHKK
     +      (NHKK),JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK
     +      (2,NHKK),(PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK
     +      =1,4)
   50     CONTINUE
C       IF(NHAD.GT.0) THEN
C         JDAHKK(1,IMOHKK)=NHKKAU
C         JDAHKK(2,IMOHKK)=NHKK
C       ENDIF
          IF(NNNPJ.GE.1)THEN
            NNNPSO=NNNPS
            NNNPS=NNNPS+1
            NNNPSU=NNNPSO+NNNPJ
            DO 187 J=NNNPS,NNNPSU
              JJ=J-NNNPS+1
              IF(J.GT.40000.OR.JJ.GT.1000)THEN
C               WRITE(6,'(A,2I10)')' J.gt.40000.or.jj.gt.1000 ',J,JJ
                GO TO 187
              ENDIF
	      PXS(J)=PXJ(JJ)
	      PYS(J)=PYJ(JJ)
	      PZS(J)=PZJ(JJ)
	      HES(J)=HEJ(JJ)
  187      CONTINUE
            NNNPS=NNNPS+NNNPJ-1
          ENDIF
        ENDIF
   60 CONTINUE
C
C--------------------------------------------------------------
C
      RETURN
 1010 FORMAT (I6,I4,5I6,9E10.2)
 1020 FORMAT (' HADRVS J.GT.NAUMAX SKIP NEXT PARTICLES ',3I10)
 1030 FORMAT (' NHKK,NREF(J),IDHKK(NHKK)  ',3I10)
 1040 FORMAT(10X,5I5,10F9.2/10X,4I5,4F12.4)
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HADRVS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C-------------------------
C
C                       HADRONIZE VALENCE-SEA CHAINS
C
C                       ADD GENERATED HADRONS TO /ALLPAR/
C                          STARTING AT (NAUX + 1)
C                       AND TO /HKKEVT/ STARTING AT (NHKK + 1)
C
C-------------------------
*KEEP,INTMX.
      PARAMETER (INTMX=2488,INTMD=252)
*KEEP,DXQX.
C     INCLUDE (XQXQ)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON /DXQX/ XPVQ(248),XPVD(248),XTVQ(248),XTVD(248), XPSQ
     +(INTMX),XPSAQ(INTMX),XTSQ(INTMX),XTSAQ(INTMX)
     *                ,XPSU(248),XTSU(248)
     *                ,XPSUT(248),XTSUT(248)
       COMMON/POPCCK/PDBCK,PDBSE,PDBSEU,
     *  IJPOCK,IREJCK,ICK4,IHAD4,ICK6,IHAD6
     *,IREJSE,ISE4,ISE6,IREJS3,ISE43,ISE63,IREJS0
     *,IHADA4,IHADA6,IREJSA,ISEA4,ISEA6,IREJA3,
     *ISEA43,ISEA63,IREJAO
*KEEP,INTNEW.
      COMMON /INTNEW/ NVV,NSV,NVS,NSS,NDV,NVD,NDS,NSD,
     +IXPV,IXPS,IXTV,IXTS, INTVV1(248),
     +INTVV2(248),INTSV1(248),INTSV2(248), INTVS1(248),INTVS2(248),
     +INTSS1(INTMX),INTSS2(INTMX),
     +INTDV1(248),INTDV2(248),INTVD1(248),INTVD2(248),
     +INTDS1(INTMD),INTDS2(INTMD),INTSD1(INTMD),INTSD2(INTMD)
 
C  /INTNEW/
C       NVV    :   NUMBER OF INTERACTING VALENCE-VALENCE SYSTEMS
C       NSV    :   NUMBER OF INTERACTING SEA-VALENCE SYSTEMS
C       NVS    :   NUMBER OF INTERACTING VALENCE-SEA SYSTEMS
C       NSS    :   NUMBER OF INTERACTING SEA-SEA SYSTEMS
C       IXPV, IXTV : NUMBER OF GENERATED X-VALUES FOR VALENCE-QUARK
C                     SYSTEMS FROM PROJECTILE/TARGET NUCLEI
C       IXPS, IXTS : NUMBER OF GENERATED X-VALUES FOR SEA-QUARK PAIRS
C                     FROM PROJECTILE/TARGET NUCLEI
C-------------------
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,LOZUO.
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/ ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),ZUOST(INTMX),
     +INTLO(INTMX),INLOSS(INTMX)
C  /LOZUO/
C       INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C------------------
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
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
*KEEP,ABRVS.
      COMMON /ABRVS/ AMCVS1(248),AMCVS2(248),GACVS1(248),GACVS2(248),
     +BGXVS1(248),BGYVS1(248),BGZVS1(248), BGXVS2(248),BGYVS2(248),
     +BGZVS2(248), NCHVS1(248),NCHVS2(248),IJCVS1(248),IJCVS2(248),
     +PQVSA1(248,4),PQVSA2(248,4), PQVSB1(248,4),PQVSB2(248,4)
*KEEP,DFINPA.
      CHARACTER*8 ANF
      PARAMETER (NFIMAX=249)
      COMMON /DFINPA/ ANF(NFIMAX),PXF(NFIMAX),PYF(NFIMAX),PZF(NFIMAX),
     +HEF(NFIMAX),AMF(NFIMAX), ICHF(NFIMAX),IBARF(NFIMAX),NREF(NFIMAX)
      COMMON /DFINPZ/IORMO(NFIMAX),IDAUG1(NFIMAX),IDAUG2(NFIMAX),
     * ISTATH(NFIMAX)
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEEP,PROJK.
      COMMON /PROJK/ IPROJK
*KEEP,NUCC.
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
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
       COMMON /NCOUCH/ ACOUVV,ACOUSS,ACOUSV,ACOUVS,
     *                 ACOUZZ,ACOUHH,ACOUDS,ACOUSD,
     *                 ACOUDZ,ACOUZD,ACOUDI,
     *                 ACOUDV,ACOUVD,ACOUCC
C---------------------
      COMMON /ZSEA/ZSEAAV,ZSEASU,ANZSEA
      COMMON /CASADI/CASAXX,ICASAD
C---------------------
      DIMENSION POJ(4),PAT(4)
C-----------------------------------------------------------------------
      DO 50 I=1,NVS
C-----------------------drop recombined chain pairs
        IF(NCHVS1(I).EQ.99.AND.NCHVS2(I).EQ.99) GO TO 50
        IS1=INTVS1(I)
        IS2=INTVS2(I)
C
        IF (IPCO.GE.6) WRITE (6,1010) IPVQ(IS1),IPPV1(IS1),IPPV2(IS1),
     +  ITSQ(IS2),ITSAQ(IS2), AMCVS1(I),AMCVS2(I),GACVS1(I),GACVS2(I),
     +  BGXVS1(I),BGYVS1(I),BGZVS1(I), BGXVS2(I),BGYVS2(I),BGZVS2(I),
     +  NCHVS1(I),NCHVS2(I),IJCVS1(I),IJCVS2(I), PQVSA1(I,4),PQVSA2
     +  (I,4),PQVSB1(I,4),PQVSB2(I,4)
C
C+++++++++++++++++++++++++++++     CHAIN 1:  QUARK-AQUARK    ++++++++++
        IFB1=IPVQ(IS1)
        IFB2=ITSAQ(IS2)
        IFB2=IABS(IFB2)+6
        DO 10 J=1,4
          POJ(J)=PQVSA1(I,J)
          PAT(J)=PQVSA2(I,J)
   10   CONTINUE
        PT1=SQRT(POJ(1)**2+POJ(2)**2)
        PT2=SQRT(PAT(1)**2+PAT(2)**2)
        CALL PARPT(2,PT1,PT2,2,NEVT)
C       IF((NCHVS1(I).NE.0.OR.NCHVS2(I).NE.0).AND.IP.NE.1)
C    &  CALL SAPTRE(AMCVS2(I),GACVS2(I),BGXVS2(I),BGYVS2(I),BGZVS2(I),
C    &              AMCVS1(I),GACVS1(I),BGXVS1(I),BGYVS1(I),BGZVS1(I))
C-----------------------------------------------------------------
C                                    POJ,PAT EXCHANGED J.R.15.2.90
C                                            RECHANGED HJM 13/2/91
      IF(IPCO.GE.1)THEN
        WRITE(6,*)' VS q-aq ,IFB1,IFB2,',
     *	'INTVS1=IS1,INTVS2=IS2,JIPPX,JITTX',
     *	IFB1,IFB2,INTVS1(I),INTVS2(I),JIPPX,JITTX
      ENDIF
        CALL HADJET(NHAD,AMCVS1(I),POJ,PAT,GACVS1(I),BGXVS1(I), BGYVS1
     +  (I),BGZVS1(I),IFB1,IFB2,IFB3,IFB4, IJCVS1(I),IJCVS1(I),3,NCHVS1
     +  (I),5)
        ACOUVS=ACOUVS+1
C                                  ADD HADRONS/RESONANCES INTO
C                                  COMMON /ALLPAR/ STARTING AT NAUX
        NHKKAU=NHKK+1
        DO 20 J=1,NHAD
          EHECC=SQRT(PXF(J)**2+PYF(J)**2+PZF(J)**2+AMF(J)**2)
        IF (ABS(EHECC-HEF(J)).GT.0.001D0.AND.IBARF(J).NE.500)THEN
            WRITE(6,'(A,2I5,2E16.6)')
     +      ' HADRVS: CORRECT INCONSISTENT PARTICLE ENERGY ', NHKK,NREF
     +      (J), HEF(J),EHECC
          HEF(J)=EHECC
          ENDIF
          ANNVS=ANNVS+1
          EEVS=EEVS+HEF(J)
          PTVS=PTVS+SQRT(PXF(J)**2+PYF(J)**2)
C                               PUT NN-CMS HADRONS INTO /HKKEVT/
C         NHKK=NHKK+1
          IF (NHKK.EQ.NMXHKK) THEN
            WRITE (6,'(A,2I5)') ' HADRVS: NHKK.EQ.NMXHKK',NHKK,NMXHKK
            RETURN
          ENDIF
          ISTIST=1
          IF(IBARF(J).EQ.500)ISTIST=2
          CALL HKKFIL(ISTIST,MPDGHA(NREF(J)),MHKKVS(I)-3,0,
     *                PXF(J),PYF(J),PZF(J),HEF(J),NHKKAU,IORMO(J),7)
          IF(IDHKK(NHKK).EQ.99999) WRITE (6,1030) NHKK, IDHKK(NHKK)
C	  WRITE(6,*)' Firt chain HADRVS'
          IF (IPHKK.GE.7) WRITE(6,1000)NHKK, ISTHKK(NHKK),IDHKK(NHKK),
     +    JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +    (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
   20   CONTINUE
C       IF(NHAD.GT.0) THEN
C         JDAHKK(1,IMOHKK)=NHKKAU
C         JDAHKK(2,IMOHKK)=NHKK
C       ENDIF
C
C++++++++++++++++++++++++++++++    CHAIN 2:  DIQUARK-QUARK   +++++++++++
        IFB1=IPPV1(IS1)
        IFB2=IPPV2(IS1)
        IFB3=ITSQ(IS2)
        DO 30 J=1,4
          POJ(J)=PQVSB2(I,J)
          PAT(J)=PQVSB1(I,J)
   30   CONTINUE
        PT1=SQRT(POJ(1)**2+POJ(2)**2)
        PT2=SQRT(PAT(1)**2+PAT(2)**2)
        CALL PARPT(2,PT1,PT2,2,NEVT)
C                                    POJ,PAT EXCHANGED J.R.15.2.90
C                                            RECHANGED HJM 21/2/91
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
C               I= number of valence chain     
C               Projectile Nr ippp= IFROVP(INTVS1(I))
C          No of Glauber sea q at Projectile JIPP=JSSHS(IPP)
       IPPP = IFROVP(INTVS1(I))
       JIPP=JSSHS(IPPP)
C	IF(NCHVS2(I).EQ.0)THEN
C      WRITE(6,'(A,3I5)')'HADRVS: I,IPPP,JIPP ',
C    *                     I,IPPP,JIPP
C	ENDIF
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
      IF(IPCO.GE.1)THEN
        WRITE(6,*)' VS qq-q ,IFB1,IFB2,IFB3,',
     *	'INTVS1=IS1,INTVS2=IS2,JIPP,JITTX',
     *	IFB1,IFB2,IFB3,INTVS1(I),INTVS2(I),JIPP,JITTX
      ENDIF
C------------------------------------------------------------------- 
C-------------------------------------------------------------------    
	IF((NCHVS2(I).NE.0))THEN
        CALL HADJET(NHAD,AMCVS2(I),POJ,PAT,GACVS2(I),BGXVS2(I), BGYVS2
     +  (I),BGZVS2(I),IFB1,IFB2,IFB3,IFB4, IJCVS2(I),IJCVS2(I),6,NCHVS2
     +  (I),6)
	ENDIF
	AACK=FLOAT(ICK6)/FLOAT(ICK6+IHAD6+1)
 	IF((NCHVS2(I).EQ.0))THEN
          ZSEAWU=RNDM(BB)*2.D0*ZSEAAV
        RSEACK=FLOAT(JIPP)*PDBSE+ ZSEAWU*PDBSEU
 	  IF(IPCO.GE.1)WRITE(6,*)'HADJSE JIPP,RSEACK,PDBSE 4 dpmnuc3',
     +    JIPP,RSEACK,PDBSE
          IREJSS=5
	  IF(RNDM(V).LE.RSEACK)THEN
	  IREJSS=2
	  IF(AMCVS2(I).GT.2.3D0)THEN
	    IREJSS=0
            CALL HADJSE(NHAD,AMCVS2(I),POJ,PAT,GACVS2(I),BGXVS2(I),
     *	    BGYVS2
     +      (I),BGZVS2(I),IFB1,IFB2,IFB3,IFB4, IJCVS2(I),IJCVS2(I),6,
     *      NCHVS2
     +      (I),6,IREJSS,IISSQQ)
 	    IF(IPCO.GE.1)WRITE(6,'(2A,I5,F10.3,I5)')'HADJSE JIPP,',
     *	    'RSEACK,IREJSS 4 dpmnuc3',
     +      JIPP,RSEACK,IREJSS
          ENDIF
	  IF(IREJSS.GE.1)THEN
	    IF(IREJSS.EQ.1)IREJSE=IREJSE+1
	    IF(IREJSS.EQ.3)IREJS3=IREJS3+1
	    IF(IREJSS.EQ.2)IREJS0=IREJS0+1
            CALL HADJET(NHAD,AMCVS2(I),POJ,PAT,GACVS2(I),
     *      BGXVS2(I), BGYVS2
     +      (I),BGZVS2(I),IFB1,IFB2,IFB3,IFB4, IJCVS2(I),
     *      IJCVS2(I),6,NCHVS2
     +      (I),6)
	    IHAD6=IHAD6+1
	  ENDIF
	  IF(IREJSS.EQ.0)THEN
	    IF(IISSQQ.EQ.3)THEN
	      ISE63=ISE63+1
	    ELSE
	      ISE6=ISE6+1
	    ENDIF
	  ENDIF  
        ELSEIF((IJPOCK.EQ.1).AND.
     *      (AACK.LE.PDBCK))THEN
	IREJ=0
        CALL HADJCK(NHAD,AMCVS2(I),POJ,PAT,GACVS2(I),BGXVS2(I), BGYVS2
     +  (I),BGZVS2(I),IFB1,IFB2,IFB3,IFB4, IJCVS2(I),IJCVS2(I),6,NCHVS2
     +  (I),6,IREJ)
	IF(IREJ.EQ.1)THEN
	IREJCK=IREJCK+1
        CALL HADJET(NHAD,AMCVS2(I),POJ,PAT,GACVS2(I),BGXVS2(I), BGYVS2
     +  (I),BGZVS2(I),IFB1,IFB2,IFB3,IFB4, IJCVS2(I),IJCVS2(I),6,NCHVS2
     +  (I),6)
	IHAD6=IHAD6+1
	ENDIF
	IF(IREJ.EQ.0)ICK6=ICK6+1
	ELSE
        CALL HADJET(NHAD,AMCVS2(I),POJ,PAT,GACVS2(I),BGXVS2(I), BGYVS2
     +  (I),BGZVS2(I),IFB1,IFB2,IFB3,IFB4, IJCVS2(I),IJCVS2(I),6,NCHVS2
     +  (I),6)
	IHAD6=IHAD6+1
	ENDIF
	ENDIF
C------------------------------------------------------------------
C------------------------------------------------------------------
C                                  ADD HADRONS/RESONANCES INTO
C                                  COMMON /ALLPAR/ STARTING AT NAUX
        NHKKAU=NHKK+1
        DO 40 J=1,NHAD
          EHECC=SQRT(PXF(J)**2+PYF(J)**2+PZF(J)**2+AMF(J)**2)
        IF (ABS(EHECC-HEF(J)).GT.0.001D0.AND.IBARF(J).NE.500) THEN
            WRITE(6,'(A,2I5,2E16.6)')
     +      ' HADRVS: CORRECT INCONSISTENT PARTICLE ENERGY ', NHKK,NREF
     +      (J), HEF(J),EHECC
          HEF(J)=EHECC
          ENDIF
          ANNVS=ANNVS+1
          EEVS=EEVS+HEF(J)
          PTVS=PTVS+SQRT(PXF(J)**2+PYF(J)**2)
C                               PUT NN-CMS HADRONS INTO /HKKEVT/
C         NHKK=NHKK+1
          IF (NHKK.EQ.NMXHKK) THEN
            WRITE (6,'(A,2I5)') ' HADRVS: NHKK.EQ.NMXHKK ',NHKK,NMXHKK
            RETURN
          ENDIF
          ISTIST=1
          IF(IBARF(J).EQ.500)ISTIST=2
          CALL HKKFIL(ISTIST,MPDGHA(NREF(J)),MHKKVS(I),0,
     *                PXF(J),PYF(J),PZF(J),HEF(J),NHKKAU,IORMO(J),8)
C         WRITE(6,*)' HKKFIL: NHKKAU,IORMO(J) ',ISTIST, NHKKAU,IORMO(J)
          IF(IDHKK(NHKK).EQ.99999) WRITE (6,1030)NHKK,NREF(J), IDHKK
     +    (NHKK)
          IF(IREJSS.LT.0)THEN 
 	  WRITE(6,*)' Second chain HADRVS'
          IF (IPHKK.GE.0) WRITE(6,1000) NHKK, ISTHKK(NHKK),IDHKK(NHKK),
     +    JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +    (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
          ENDIF
   40   CONTINUE
C       IF(NHAD.GT.0) THEN
C         JDAHKK(1,IMOHKK)=NHKKAU
C         JDAHKK(2,IMOHKK)=NHKK
C       ENDIF
   50 CONTINUE
C
      RETURN
 1000 FORMAT (I6,I4,5I6,9E10.2)
 1010 FORMAT(10X,5I5,10F9.2/10X,4I5,4F12.4)
 1020 FORMAT (' HADRVS J.GT.NAUMAX SKIP NEXT PARTICLES ',3I10)
 1030 FORMAT (' NHKK,IDHKK(NHKK)  ',3I10)
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HADJET(NHAD,AMCH,PPR,PTA,GAM,BGX,BGY,BGZ, IFB1,IFB2,
     +IFB3,IFB4,I1,I2,NOBAM,NNCH,NORIG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C  HADJET DOES ALL THE NECESSARY ROTATIONS AND LORENTZ TRANSFORMS AND
C                                               CALL CALBAM (BAMJET)
C
C         NHAD = NUMBER OF HADRONS CREATED (OUTPUT) = IHAD (CALBAM)
C******** PRODUCED PARTICLES IN COMMON /CMSRES/
C NOTE:  NOW IN /FINPAR/      HJM 21/06/90
C         AMCH = INVARIANT MASS OF JET (INPUT)
C         PPR  = 4-MOMENTUM OF FORWARD GOING PARTON (PROJECTILE)(INPUT)
C         PTA  = 4-MOMENTUM OF BACKWARD GOING PARTON (TARGET)(INPUT)
C         GAM,BGX,BGY,BGZ = LORENTZ PARAMETERS OF JET CMS RELATIVE TO
C                                              COLLISION CMS (INPUT)
C
C--------------------------------------------------------------------
C    CALBAM(NNCH,I1,I2,IFB1,IFB2,IFB3,IFB4,AMCH,NOBAM,IHAD)
C    SAMPLING OF Q-AQ, Q-QQ, QQ-AQQ CHAINS
C       USING BAMJET(IHAD,IFB1,IFB2,IFB3,IFB4,AMCH,NOBAM)-----FOR NNCH=0
C       OR    PARJET(IHAD,ICH1=I1 OR I2)------FOR NNCH=-1 OR +1
C-------------------------------------------------------------------
C   IHAD  :  NUMBER OF PRODUCED PARTICLES
C   NNCH  :  CALL BAMJET FOR NNCH=0
C            CALL PARJET FOR NNCH=+1 ICH1=I1
C                        FOR NNCH=-1 ICH1=I2
C            jet not existing for NNCH=+/-99, i.e. IHAD=0
C   PRODUCED PARTICLES IN CHAIN REST FRAME ARE IN COMMON /FINPAR/
C   AMCH  :  INVARIANT MASS OF CHAIN (BAMJET)
C
C   NOBAM :  = 3  QUARK-ANTIQUARK JET   QUARK FLAVORS : IFB1,IFB2
C                 OR ANTIQUARK-QUARK JET                  IN ANY ORDER
C
C            = 4  QUARK-DIQUARK JET, FLAVORS : QU : IFB1, DIQU :IFB2,IFB
C                 OR ANTIQUARK-ANTIDIQUARK JET
C
C
C            = 5  DIQUARK-ANTIDIQUARK JET
C                 OR ANTIDIQUARK-DIQUARK JET
C                     FLAVORS :  DIQU : IFB1,IFB2, ANTIDIQU : IFB3,IFB4
C                                IN ANY ORDER
C
C            = 6  DIQUARK-QUARK JET, FLAVORS : DIQU  : IFB1,IFB2 QU: IFB
C                 OR ANTIDIQUARK-ANTIQUARK JET
C
C   IFBI  : FLAVORS : 1,2,3,4 = U,D,S,C    7,8,9,10 = AU,AD,AS,AC
C
C   I1,I2 : NUMBER LABEL OF A HADRON CREATED BY PARJET
C
C   NORMALLY IN BAMJET THE QUARKS MOVE FORWARD  (POSITIVE Z-DIRECTION)
C                      THE QUARK FLAVORS ARE FIRST GIVEN
C   CALBAM ALLOWS EITHER THE QUARK OR ANTIQUARK (DIQU) TO MOVE FORWARD
C                      THE FORWARD GOING FLAVORS ARE GIVEN FIRST
C
C =====================================================================
*KEEP,DFINPA.
      CHARACTER*8 ANF
      PARAMETER (NFIMAX=249)
      COMMON /DFINPA/ ANF(NFIMAX),PXF(NFIMAX),PYF(NFIMAX),PZF(NFIMAX),
     +HEF(NFIMAX),AMF(NFIMAX), ICHF(NFIMAX),IBARF(NFIMAX),NREF(NFIMAX)
      COMMON /DFINPZ/IORMO(NFIMAX),IDAUG1(NFIMAX),IDAUG2(NFIMAX),
     * ISTATH(NFIMAX)
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEND.
*--------------------------------------------------------------------
*-------------------------------------------------------------------
      COMMON /JSPART/PXP(1000),PYP(1000),PZP(1000),HEPP(1000),NNNP
      COMMON /JSPAR/PXJ(1000),PYJ(1000),PZJ(1000),HEJ(1000),NNNPJ
************************************************************************
************************************************************************
      COMMON/IFRAGM/IFRAG
      DIMENSION PPR(4),PTA(4),PPRJ(4),PTAJ(4)
      LOGICAL LTESHA
      PARAMETER (TINY=1.D-10)
      DATA ICHECO/0/
      DATA ICHECA/0/
C     IPCO=6
C----------------------------------------------------------------------
C
C  CHECK HADJET ENERGY CONSERVATION PPR(4)+PTA(4)  EQ  EHAD
C                            TRANSFORM PROJECTILE INTO JET REST FRAME
C     IF(NOBAM.EQ.4.OR.NOBAM.EQ.6) IPCO=6
      IF(IPCO.GE.6) THEN
        WRITE(6,1010) GAM,BGX,BGY,BGZ,PPR,PPRJ
        WRITE(6,1000) NHAD,AMCH,PTA, IFB1,IFB2,IFB3,IFB4,I1,I2,NOBAM,
     +  NNCH,NORIG
 1000 FORMAT(10X,I10,5F10.3/10X,9I10)
      ENDIF
      IF(ABS(NNCH).EQ.99) THEN
        NHAD=0
C       IPCO=0
        RETURN
      ENDIF
C
      CALL DALTRA(GAM,-BGX,-BGY,-BGZ,PPR(1),PPR(2),PPR(3),PPR(4),PPRTOT,
     +PPRJ(1),PPRJ(2),PPRJ(3),PPRJ(4))
      CALL DALTRA(GAM,-BGX,-BGY,-BGZ,PTA(1),PTA(2),PTA(3),PTA(4),PTATOT,
     +PTAJ(1),PTAJ(2),PTAJ(3),PTAJ(4))
C
      IF(IPCO.GE.3) WRITE(6,1010)GAM,BGX,BGY,BGZ,PPR,PPRJ,PTA,PTAJ
 1010 FORMAT(' HADJET: GAM,BGX,BGY,BGZ,PPR(4),PPRJ(4) ',4F15.5/8F15.5/ 8
     +F15.5)
C                            WORK OUT COD,SID,COF,SIF OF PROJECTILE
C                                                     IN JET FRAME
      IF(PPRTOT.LT.TINY)PPRTOT=TINY
      COD= PPRJ(3)/PPRTOT
      IF(COD.GE.1.D0)COD=0.999999999999
      IF(COD.LE.-1.D0)COD=-0.999999999999
      SID=SQRT(ABS((1.D0-COD)*(1.D0+COD)))
      COF=1.
      SIF=0.
      IF((ABS(PPRJ(1)).GT.0.D0).OR.(ABS(PPRJ(2)).GT.0.D0))THEN
      IF(PPRTOT*SID.GT.1.D-9) THEN
        COF=PPRJ(1)/(SID*PPRTOT)
        SIF=PPRJ(2)/(SID*PPRTOT)
        ANORF=SQRT(ABS(COF*COF+SIF*SIF))
        COF=COF/ANORF
        SIF=SIF/ANORF
      ENDIF
      ENDIF
      IF (IPCO.GE.6) WRITE(6,1020)COD,SID,COF,SIF
 1020 FORMAT(' COD,SID,COF,SIF ',4F15.8)
C                            SAMPLE JET IN JET CMS
      CALL CALBAM(NNCH,I1,I2,IFB1,IFB2,IFB3,IFB4,AMCH,NOBAM,IHAD)
C     IF(IHAD.EQ.1)THEN
C      IPRI=1
C      IPCO=3
C     ELSE
C      IPRI=0
C      IPCO=-1
C     ENDIF
C                            NOW WE HAVE IHAD PARTICLES/RESONANCES
C                                             IN COMMON /FINPAR/
C CHECK CALBAM ENERGY CONSERVATION / jet cms
      ECAL=0.
      PXCAL=0.
      PYCAL=0.
      PZCAL=0.
      LTESHA=.FALSE.
      DO 10 I=1,IHAD
        IF(IBARF(I).EQ.500)GO TO 1011
        PXCAL=PXCAL + PXF(I)
        PYCAL=PYCAL + PYF(I)
        PZCAL=PZCAL + PZF(I)
        IF(IFRAG.EQ.0)THEN
        EHECC=SQRT(ABS(PXF(I)**2+PYF(I)**2+PZF(I)**2+AMF(I)**2))
      IF (ABS(EHECC-HEF(I)).GT.0.0001D0) THEN
        IF(IPRI.GE.1) WRITE(6,'(2A/A/3I5,3E16.6)')
     +     ' HADJET / AFTER CALBAM:',
     +     '  CORRECT INCONSISTENT ENERGY IN JET CMS',
     +            '  I, IHAD,NREF(I), HEF(I),EHECC, AMF(I)',
     *            I,IHAD,NREF(I), HEF(I),EHECC, AMF(I)
        HEF(I)=EHECC
        LTESHA=.TRUE.
        ENDIF
        ENDIF
      ECAL=ECAL + HEF(I)
 1011 CONTINUE
   10 CONTINUE
      IF (ABS(ECAL-AMCH).GT.0.005D0) LTESHA=.TRUE.
      IF (LTESHA) THEN
        ICHECA=ICHECA+1
      IF (ABS(ECAL-AMCH).GT.0.005D0) THEN
        IF(ICHECA.LE.10)WRITE(6,'(A/10I4)')
     +  ' HADJET/1:ICHECA,IFB1,IFB2,IFB3,IFB4,I1,I2,NOBAM,NNCH,NORIG',
     +     ICHECA,IFB1,
     +  IFB2,IFB3,IFB4,I1,I2,NOBAM,NNCH,NORIG
        IF(ICHECA.LE.10)WRITE(6,1030) AMCH,ECAL,PXCAL,PYCAL,PZCAL
 1030 FORMAT(' CALBAM E. CHECK (5 MeV) AMCH,ECAL,PXCAL,PYCAL,PZCAL=',
     +    /5E20.8)
      LTESHA=.FALSE.
      ENDIF
      ENDIF
      IF (IPCO.GE.3)THEN
        DO 20 I=1,IHAD
        WRITE(6,1040)I,PXF(I),PYF(I),PZF(I),HEF(I),AMF(I), ICHF(I),
     +    IBARF(I),NREF(I),ANF(I)
 1040 FORMAT(' JET SYSTEM ',I5,5F12.4,3I5,A10)
   20   CONTINUE
      ENDIF
C     CALL DECAY(IHAD,2)
C     NHAD=IHAD
C CHECK CALBAM+DECAY ENERGY CONSERVATION
C     ECAL=0.
C     PXCAL=0.
C     PYCAL=0.
C     PZCAL=0.
C     DO 1204 I=1,IHAD
C       ECAL=ECAL+HEF(I)
C       PXCAL=PXCAL+PXF(I)
C       PYCAL=PYCAL+PYF(I)
C       PZCAL=PZCAL+PZF(I)
C1204 CONTINUE
C     IOUCHA=3
C     IF (ABS(ECAL-AMCH)/AMCH.GT.0.00001D0)IOUCHA=1
C     IF (IPCO.GE.IOUCHA)WRITE(6,1203)AMCH,ECAL,PXCAL,PYCAL,PZCAL
C1203 FORMAT(' CALBAM+DECAY ENERGY CHECK  AMCH,ECAL,PXCAL,PYCAL,PZCAL '
C    *                                       /5E20.8)
C     IF (IPCO.GE.3)THEN
C       DO 143 I=1,IHAD
C         WRITE(6,144)I,PXF(I),PYF(I),PZF(I),HEF(I),AMF(I),
C    *                  ICHF(I),IBARF(I),NREF(I),ANF(I)
C 144     FORMAT(' JET SYSTEM DECAY ',I5,5F12.4,3I5,A10)
C 143   CONTINUE
C     ENDIF
C                            ROTATE JET BY COD,SID,COF,SIF
      PXCAL=0.
      PYCAL=0.
      PZCAL=0.
      LTESHA=.FALSE.
      DO 30 I=1,IHAD
      PHEC2=PXF(I)**2+PYF(I)**2+PZF(I)**2
        CALL DTRANS(PXF(I),PYF(I),PZF(I),COD,SID,COF,SIF,XX,YY,ZZ)
      PROTA2=XX**2 + YY**2 + ZZ**2
        PXF(I)=XX
        PYF(I)=YY
        PZF(I)=ZZ
      PXCAL=PXCAL + PXF(I)
      PYCAL=PYCAL + PYF(I)
      PZCAL=PZCAL + PZF(I)
C       EHECC=SQRT(ABS(PXF(I)**2+PYF(I)**2+PZF(I)**2+AMF(I)**2))
      IF(ABS(PHEC2-PROTA2).GT.0.0001D0) THEN
        WRITE(6,'(2A/3I5,3E16.6)')
     &            ' HADJET: INCONSISTENT MOMENTUM AFTER TRANS',
     *            '  I, IHAD,NREF(I), PHEC2,PROTA2, AMF(I)',
     *            I,IHAD,NREF(I), PHEC2,PROTA2, AMF(I)
C         HEF(I)=EHECC
        LTESHA=.TRUE.
        ENDIF
   30 CONTINUE
      IF (LTESHA) THEN
      WRITE(6,'(A/9I4)')
     +  ' HADJET/2: IFB1,IFB2,IFB3,IFB4,I1,I2,NOBAM,NNCH,NORIG', IFB1,
     +  IFB2,IFB3,IFB4,I1,I2,NOBAM,NNCH,NORIG
      WRITE(6,1031) PXCAL,PYCAL,PZCAL
 1031   FORMAT(' CALBAM ENERGY CHECK/2:  PXCAL,PYCAL,PZCAL='/3E20.8)
      LTESHA=.FALSE.
      ENDIF
      IF (IPCO.GE.3)THEN
        DO 40 I=1,IHAD
        WRITE(6,1050)I,PXF(I),PYF(I),PZF(I),HEF(I),AMF(I), ICHF(I),
     +    IBARF(I),NREF(I),ANF(I)
 1050 FORMAT(' ROTATET JET SYSTEM ',I5,5F12.4,3I5,A10)
   40   CONTINUE
      ENDIF
************************************************************************
************************************************************************
C                           Rotate partons
      NNNNP=NNNP
      IF(NNNP.GT.1000)NNNNP=1000
      DO 1105 I=1,NNNNP
        CALL DTRANS(PXP(I),PYP(I),PZP(I),COD,SID,COF,SIF,XX,YY,ZZ)
        PXP(I)=XX
        PYP(I)=YY
        PZP(I)=ZZ
 1105 CONTINUE
************************************************************************
************************************************************************
C                           TRANSFORM THIS JET BACK INTO CMS
C************************                  IN COMMON BLOCK/CMSRES/
C                             LORTMO USES /FINPAR/ ONLY !!!!!
C     CALL LORTRA(IHAD,1,GAM,BGX,BGY,BGZ)
      CALL LORTMO(IHAD,GAM,BGX,BGY,BGZ)
      NHAD=IHAD
C
      IF (IPCO.GE.3)THEN
        DO 50 I=1,IHAD
        WRITE(6,1060) I,PXF(I),PYF(I),PZF(I),HEF(I),AMF(I), ICHF(I),
     +    IBARF(I),NREF(I),ANF(I)
 1060 FORMAT(' CMS SYSTEM ',I5,5F12.4,3I5,A10)
   50   CONTINUE
      ENDIF
C  HADJET ENERGY CONSERVATION TEST
C  AND CONSISTENCY TEST OF PARTICLE 4-MOMENTUM
      LTESHA=.FALSE.
      EHAD=0.
      DO 60 I=1,IHAD
      IF(IBARF(I).EQ.500)GO TO 6060
      EHAD=EHAD+HEF(I)
        EHECC=SQRT(ABS(PXF(I)**2+PYF(I)**2+PZF(I)**2+AMF(I)**2))
      IF (ABS(EHECC-HEF(I)).GT.0.001D0) THEN
        IF(ABS(EHECC-HEF(I)).GT.0.1D0)WRITE(6,'(2A/4I5,3E16.6)')
     &            ' HADJET: CORRECT INCONSISTENT ENERGY AFTER LORTRA',
     *            '  NORIG, I, IHAD,NREF(I), HEF(I),EHECC, AMF(I)',
     *            NORIG, I,IHAD,NREF(I), HEF(I),EHECC, AMF(I)
        HEF(I)=EHECC
C         LTESHA=.TRUE.
        ENDIF
 6060 CONTINUE
   60 CONTINUE
************************************************************************
************************************************************************
C                           Transform partons
      NNNNP=NNNP
      IF(NNNP.GT.1000)NNNNP=1000
      NNNPJ=NNNNP
      CALL LORTRP(NNNNP,1,GAM,BGX,BGY,BGZ)
************************************************************************
************************************************************************
C     IOUCHA=3
c     EEIN=PPR(4)+PTA(4)
c     PXIN=PPR(1)+PTA(1)
c     PYIN=PPR(2)+PTA(2)
c     PZIN=PPR(3)+PTA(3)
C     PXIN=BGX*AMCH
C     PYIN=BGY*AMCH
C     PZIN=BGZ*AMCH
      EEIN=GAM*AMCH
      IF (ABS(EEIN-EHAD).GT.0.005D0) THEN
        ICHECO=ICHECO+1
      IF (ABS(EEIN-EHAD).GT.0.005D0) THEN
	IF(ICHECO.LT.10)THEN
        WRITE(6,'(A/10I5)')
     +  ' HADJET/3:ICHECO,IFB1,IFB2,IFB3,IFB4,I1,I2,NOBAM,NNCH,NORIG', 
     +    ICHECO,IFB1,
     +  IFB2,IFB3,IFB4,I1,I2,NOBAM,NNCH,NORIG
      WRITE (6,1070) EEIN,EHAD,AMCH,GAM,BGX,BGY,BGZ
 1070 FORMAT(' HADJET ENERGY CHECK (5 MeV) EEIN,EHAD,AMCH',3E20.8/
     +  20X,' GAM,BGX,BGY,BGZ ',4E20.8)
c       WRITE(6,1010)GAM,BGX,BGY,BGZ,PPR,PPRJ
	ENDIF
      ENDIF
      ENDIF
C     IPCO=0
      RETURN
      END

************************************************************************
************************************************************************
*
      SUBROUTINE LORTRP(N,NAUX,GAM,BGX,BGY,BGZ)
*
*     LORENTZ TRANSFORMATION OF  N PARTICLES IN  JSPART  TO BE
*     STORED IN  JSPAR  STARTING AT NAUX
*
*********************************************************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C     impl. mxnupa after KNO cut solved 3.92
      PARAMETER (MXNUPA=2500)
      COMMON /JSPART/PXP(1000),PYP(1000),PZP(1000),HEPP(1000),NNNP
      COMMON /JSPAR/PXJ(1000),PYJ(1000),PZJ(1000),HEJ(1000),NNNPJ
*          CHANGED JUNE 1,1987
       PARAMETER (ONE=1.D0)
       DATA IFIRST/0/
       DATA NUM/0/
       IFIRST=IFIRST+1
       NUM=NUM+1
      PXSM=0.0
      PYSM=0.0
      PZSM=0.0
      ESUM=0.
      PXSC=0.0
      PYSC=0.0
      PZSC=0.0
      ESMC=0.0
*           END OF CHANGE
      DO 1  I=1,N
      J = NAUX + I - 1
      PXSM=PXSM+PXP(I)
      PYSM=PYSM+PYP(I)
      PZSM=PZSM+PZP(I)
      ESUM=ESUM+HEPP(I)
      CALL DALTRA(GAM,BGX,BGY,BGZ,PXP(I),PYP(I),PZP(I),HEPP(I),
     *PPA,PXJ(J),PYJ(J),PZJ(J),HEJ(J))
      PXSC=PXSC+PXJ(J)
      PYSC=PYSC+PYJ(J)
      PZSC=PZSC+PZJ(J)
      ESMC=ESMC+HEJ(J)
1     CONTINUE
*     PXSM,ETC,ARE SUMS FOR BAMJET FRAGMENTS IN JET CMS
      CALL DALTRA(GAM,BGX,BGY,BGZ,PXSM,PYSM,PZSM,ESUM,
     *PPA,PXSM,PYSM,PZSM,ESUM)
*     PXSC,ETC,ARE SUMS FOR BAMJET FRAGMENTS IN PROJ,TARGET CMS
 
      PXDIF=PXSM-PXSC
      PYDIF=PYSM-PYSC
      PZDIF=PZSM-PZSC
      EDIF=ESUM-ESMC
      DIFFL=PXDIF+PYDIF+PZDIF+EDIF
      IF(DIFFL.GE.1.D-2*ONE)
     1WRITE(6,2)NUM,PXDIF,PYDIF,PZDIF,EDIF,PXSM,PXSC,
C     WRITE(6,2)NUM,PXDIF,PYDIF,PZDIF,EDIF,PXSM,PXSC,
     1PYSM,PYSC,PZSM,PZSC,ESUM,ESMC
 2    FORMAT(' ',2X,'LORTRA:NUM=',I5,2X,'PXDIF=',1PE15.6,2X,'PYDIF=',
     21PE15.6,2X,'PZDIF=',1PE15.6,2X,'EDIF=',1PE15.6/2X,'PXSM=',1PE15.6,
     32X,'PXSC=',1PE15.6,2X,'PYSM=',1PE15.6,2X,'PYSC=',1PE15.6/2X,'PZSM'
     4,1PE15.6,2X,'PZSC=',1PE15.6,2X,'ESUM=',1PE15.6,2X,'ESMC=',1PE15.6/
     52X,'LORTRA DIFFERENCES DUE TO ALTRA'/)
      RETURN
      END
*
************************************************************************
************************************************************************
*-- Author :
      SUBROUTINE COBCMA(IF1,IF2,IF3,IJNCH,NNCH,IREJ,AMCH,AMCHN,IKET)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C  REPLACE SMALL MASS BARYON CHAINS (AMCH)
C  BY OCTETT OR DECUPLETT BARYONS
C
C                       HERE ONLY THE CHAIN MASS IS CHANGED
C                      (AMCHN) BUT NO CORRECTION OF KINEMATICS!
C
C  MASS CORRECTED FOR NNCH.NE.0
C
C  IREJ=1: CHAIN GENERATION NOT ALLOWED BECAUSE OF TOO SMALL MASS
C          START FROM THE BEGINNING IN HAEVT
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
*KEEP,KETMAS.
      COMMON /KETMAS/ AM8(2),AM10(2),IB88(2),IB1010(2),AMCH1N,AMCH2N
*KEND.
C----------------
      CALL DBKLAS(IF1,IF2,IF3,IB8,IBB10)
C
      IF (IPEV.GE.2)WRITE(6,1000)IF1,IF2,IF3,IB8,IBB10
 1000 FORMAT (' COBCMA: IPQ,ITTQ1,ITTQ2,IB8,IBB10 ',5I5)
C
      AM81=AAM(IB8)
      AM101=AAM(IBB10)
      AM8(IKET)=AM81
      AM10(IKET)=AM101
      IB88(IKET)=IB8
      IB1010(IKET)=IBB10
      NNCH=0
      IJNCH=0
      IREJ=0
      AMFF1=AM101+0.3
C
C                                    j.r.10.5.93
C     IF(AMCH.LT.AMFF1) THEN
C       IREJ=1
C       RETURN
C     ENDIF
C                                     -------------
      IF(AMCH.LT.AM81) THEN
        IREJ=1
      ELSEIF (AMCH.LT.AM101)THEN
C                                PRODUCE OKTETT BARYON
C                                 CORRECT KINEMATICS
        IJNCH=IB8
        NNCH=-1
        AMCHN=AM81
      ELSEIF(AMCH.LT.AMFF1) THEN
C                                 PRODUCE DECUPLETT BARYON
C                                 CORRECT KINEMATICS
        AMCHN=AM101
        IJNCH=IBB10
        NNCH=1
      ELSE
        AMCHN=AMCH
      ENDIF
C                                 NO CORRECTIONS BUT DO  CHAIN 2
      IF(IPEV.GE.2) THEN
        WRITE(6,1010) AMCH,AMCHN,AM81,AM101
        WRITE(6,1020) IF1,IF2,IF3,IB8,IBB10,IJNCH,NNCH,IREJ
 1010 FORMAT(' COBCMA: AMCH,AMCHN,AM81,AM101', 4F13.4)
 1020 FORMAT(' COBCMA: IF1,IF2,IF3,IB8,IBB10,IJNCH,NNCH,IREJ',8I4)
      ENDIF
      RETURN
      END
*-- Author :
C
C++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE COMCMA(IFQ,IFAQ,IJNCH,NNCH,IREJ,AMCH,AMCHN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C  REPLACE SMALL MASS MESON CHAINS BY PSEUDOSCALAR OR VECTOR MESONS
C
C                       HERE ONLY THE CHAIN MASS IS CHANGED
C                      (AMCHN) BUT NO CORRECTION OF KINEMATICS!
C
C
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
*KEEP,DINPDA.
      COMMON /DINPDA/ IMPS(6,6),IMVE(6,6),IB08(6,21),IB10(6,21), IA08
     +(6,21),IA10(6,21), A1,B1,B2,B3,LT,LE,BET,AS,B8,AME,DIQ,ISU
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEND.
C------------------------
      IIFAQ=IABS(IFAQ)
      IFPS=IMPS(IIFAQ,IFQ)
      IFV=IMVE(IIFAQ,IFQ)
      IF (IPEV.GE.2)WRITE (6,1000)IIFAQ,IFQ,IFPS,IFV
 1000 FORMAT (' COMCMA',5X,' IIPPAQ,ITQ,IFPS,IFV ',4I5)
      AMPS=AAM(IFPS)
      AMV=AAM(IFV)
      NNCH=0
      IJNCH=0
      IREJ=0
      AMFF=AMV+0.3
      IF(IPEV.GE.2) WRITE(6,1010) AMCH,AMPS,AMV,IFPS,IFV
 1010 FORMAT(' AMCH,AMPS,AMV,IFPS,IFV ',3F12.4,2I10)
C                                          j.r.10.5.93
C     IF(AMCH.LT.AMFF) THEN
C       IREJ=1
C       RETURN
C     ENDIF
C                                          ------------
C
      IF(AMCH.LT.AMPS) THEN
        IREJ=1
        RETURN
      ENDIF
C
      IF (AMCH.LT.AMV) THEN
C                                PRODUCE PSEUDO SCALAR
        IJNCH=IFPS
        NNCH=-1
        AMCHN=AMPS
      ELSEIF(AMCH.LT.AMFF) THEN
C                                 PRODUCE VECTOR MESON
        IJNCH=IFV
        NNCH=1
        AMCHN=AMV
      ELSE
        AMCHN=AMCH
      ENDIF
      IF(IPEV.GE.2) THEN
        WRITE(6,1030) AMCH,AMCHN,AMPS,AMV
        WRITE(6,1020) IFQ,IFAQ,IFPS,IFV,IJNCH,NNCH,IREJ
 1030 FORMAT(' COMCMA: AMCH,AMCHN,AMPS,AMV', 4F13.4)
 1020 FORMAT(' COMCMA: IFQ,IFAQ,IFPS,IFV,IJNCH,NNCH,IREJ',8I4)
      ENDIF
C
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C   17/10/89 910191458  MEMBER NAME  MCOMCM2  (KK89.S)      F77
      SUBROUTINE COMCM2(IQ1,IQ2,IAQ1,IAQ2,NNCH,IREJ,AMCH)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
 
C--------------------------------------------------------------------
C  (QQ)-(AQ AQ) CHAIN:
C                      CHECK QUANTUM NUMBERS AND
C                      CORRECT MASS IF NECESSARY
C             REJECT IF THERE IS NO CORRESPONDING PARTICLE
C                    OR TOO LOW MASS
C--------------------------------------------------------------------
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
*KEEP,DINPDA.
      COMMON /DINPDA/ IMPS(6,6),IMVE(6,6),IB08(6,21),IB10(6,21), IA08
     +(6,21),IA10(6,21), A1,B1,B2,B3,LT,LE,BET,AS,B8,AME,DIQ,ISU
*KEND.
C--------------------------
      IREJ=0
      IIAQ1=-IAQ1
      IIAQ2=-IAQ2
      IF (IIAQ1.EQ.IQ1)                                         GO TO 10
      IF (IIAQ1.EQ.IQ2)                                         GO TO 20
      IF (IIAQ2.EQ.IQ1)                                         GO TO 30
      IF (IIAQ2.EQ.IQ2)                                         GO TO 40
C                                  REJECTION: NO CANCELLATION OF
C                                             ANY (Q-AQ) PAIR
      IREJ=1
      IF(IPEV.GE.3) THEN
        WRITE(6,'(A/5X,4I5,1PE13.5)')
     +  ' KKEVVV/COMCM2 (QU. NUMBERS): IQ1, IQ2, IAQ1, IAQ2, AMCH', IQ1,
     +  IQ2, IAQ1, IAQ2, AMCH
      ENDIF
      RETURN
C
   10 CONTINUE
C     IFPS=IMPS(IIAQ2,IQ2)
C     IFV =IMVE(IIAQ2,IQ2)
                                                                GO TO 50
   20 CONTINUE
C     IFPS=IMPS(IIAQ2,IQ1)
C     IFV =IMVE(IIAQ2,IQ1)
                                                                GO TO 50
   30 CONTINUE
C     IFPS=IMPS(IIAQ1,IQ2)
C     IFV =IMVE(IIAQ1,IQ2)
                                                                GO TO 50
   40 CONTINUE
C     IFPS=IMPS(IIAQ1,IQ1)
C     IFV =IMVE(IIAQ1,IQ1)
C
   50 CONTINUE
C     AMFPS=AAM(IFPS)
C     AMFV =AAM(IFV)
C     AMFF=AMFV+0.3
C                     EMPIRICAL DEFINITION OF AMFF
C                     TO ALLOW FOR (B-ANTIB) PAIR PRODUCTION
      AMFF=2.5
      NNCH=0
      IF (AMCH.LT.AMFF) THEN
        IREJ=1
        IF(IPEV.GE.3) THEN
          WRITE(6,'(A/5X,4I5,1PE13.5)')
     +    ' KKEVVV/COMCM2 (MASS!): IQ1, IQ2, IAQ1, IAQ2, AMCH', IQ1,
     +    IQ2, IAQ1, IAQ2, AMCH
        ENDIF
      ENDIF
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE CORMOM(AMCH1,AMCH2,AMCH1N,AMCH2N, 
     +PQ1X,PQ1Y,PQ1Z,PQ1E,PA1X,PA1Y,PA1Z,PA1E, PQ2X,PQ2Y,PQ2Z,PQ2E,PA2X,
     +PA2Y,PA2Z,PA2E, PXCH1,PYCH1,PZCH1,ECH1, PXCH2,PYCH2,PZCH2,ECH2,
     +IREJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C   CORRECT KINEMATICS IF MASS OF THE FIRST CHAIN HAS BEEN CHANGED
C                              FROM AMCH1 TO AMCH1N
C   CHAIN 1:  (XP,XTVD)
C   AMMM   :  TOTAL MASS OF TWO CHAIN SYSTEM
C   AMCH2N :  RESULTING NEW MASS FOR CHAIN 2 (OUTPUT ONLY)
C
C---                        RESCALING OF X-VALUES
C                           ACCORDING TO THE MODIFIED MASS OF CHAIN 1
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEND.
C------------------------------------
C     WRITE(6,'(A,4(7E15.5/))')' CORMOM IN',AMCH1,AMCH2,AMCH1N,AMCH2N, 
C    +PQ1X,PQ1Y,PQ1Z,PQ1E,PA1X,PA1Y,PA1Z,PA1E, PQ2X,PQ2Y,PQ2Z,PQ2E,PA2X,
C    +PA2Y,PA2Z,PA2E, PXCH1,PYCH1,PZCH1,ECH1, PXCH2,PYCH2,PZCH2,ECH2
C
      IF(AMCH1.EQ.0.D0)THEN
        IREJ=1
        WRITE(6,*) 'Error in CORMOM : AMCH1=0. Event rejected'
        RETURN
      ENDIF
      FAK=AMCH1N/AMCH1
C     WRITE(6,'(A,F10.5)')' FAK ',FAK
      AMCH1=AMCH1N
C
      PQ1XOL=PQ1X
      PQ1YOL=PQ1Y
      PQ1ZOL=PQ1Z
      PQ1EOL=PQ1E
      PA1XOL=PA1X
      PA1YOL=PA1Y
      PA1ZOL=PA1Z
      PA1EOL=PA1E
      PQ2XOL=PQ2X
      PQ2YOL=PQ2Y
      PQ2ZOL=PQ2Z
      PQ2EOL=PQ2E
      PA2XOL=PA2X
      PA2YOL=PA2Y
      PA2ZOL=PA2Z
      PA2EOL=PA2E
C
      PXCH1O=PXCH1
      PYCH1O=PYCH1
      PZCH1O=PZCH1
      ECH10=ECH1
      PXCH2O=PXCH2
      PYCH2O=PYCH2
      PZCH2O=PZCH2
      ECH20=ECH2
C
C
C---                        RESCALING OF MOMENTA FOR PARTONS OF CHAIN 1
      PQ1X=PQ1X*FAK
      PQ1Y=PQ1Y*FAK
      PQ1Z=PQ1Z*FAK
      PQ1E=PQ1E*FAK
      PA2X=PA2X*FAK
      PA2Y=PA2Y*FAK
      PA2Z=PA2Z*FAK
      PA2E=PA2E*FAK
C                          NEW MOMENTA OF PARTONS OF CHAIN 2
C                          FROM MOMENTUM CONSERVATION
      PA1X=PA1XOL+PQ1XOL-PQ1X   
      PA1Y=PA1YOL+PQ1YOL-PQ1Y   
      PA1Z=PA1ZOL+PQ1ZOL-PQ1Z   
      PA1E=PA1EOL+PQ1EOL-PQ1E   
      PQ2X=PQ2XOL+PA2XOL-PA2X   
      PQ2Y=PQ2YOL+PA2YOL-PA2Y   
      PQ2Z=PQ2ZOL+PA2ZOL-PA2Z   
      PQ2E=PQ2EOL+PA2EOL-PA2E   
C---                        NEW MOMENTUM OF CHAIN 1
      PXCH1=PQ1X+PA2X
      PYCH1=PQ1Y+PA2Y
      PZCH1=PQ1Z+PA2Z
      ECH1 =PQ1E+PA2E
C
C     WRITE(6,'(A,4(7E15.5/))')' CORMOM MOD',AMCH1,AMCH2,AMCH1N,AMCH2N, 
C    +PQ1X,PQ1Y,PQ1Z,PQ1E,PA1X,PA1Y,PA1Z,PA1E, PQ2X,PQ2Y,PQ2Z,PQ2E,PA2X,
C    +PA2Y,PA2Z,PA2E, PXCH1,PYCH1,PZCH1,ECH1, PXCH2,PYCH2,PZCH2,ECH2
C
      ROOT =(ECH1-AMCH1)*(ECH1+AMCH1)
      IF(ROOT.LT.0.D0)THEN
        IREJ=1
        WRITE(6,*)'Error in CORMOM : ROOT<0. Event rejected'
        WRITE(6,*)'ECH1=',ECH1,'   AMCH1=',AMCH1,'   ROOT=',ROOT
        RETURN
      ENDIF
      PCH1 = SQRT(ROOT) + 0.000001
C
C---                        NEW 4-MOMENTUM OF CHAIN 2
      PXCH2=PA1X+PQ2X
      PYCH2=PA1Y+PQ2Y
      PZCH2=PA1Z+PQ2Z
      ECH2 =PA1E+PQ2E
      PCH2 =SQRT(PXCH2**2+PYCH2**2+PZCH2**2)
      AMCH22=ECH2**2-PXCH2**2-PYCH2**2-PZCH2**2
C
      IF(AMCH22.LT.0.D0)THEN
        IREJ=1
        WRITE(6,*)'Error in CORMOM : AMCH22<0. Event rejected'
C        WRITE(6,*)'ECH2=',ECH2,'   PXCH2=',PXCH2,'   PYCH2=',PYCH2,
C     *            '   PZCH2=',PZCH2
        RETURN
      ENDIF
      AMCH2N=SQRT(AMCH22)
C---
      IF(IPRI.GT.1) THEN
        PXSUM=PQ1X+PA1X+PQ2X+PA2X
        PYSUM=PQ1Y+PA1Y+PQ2Y+PA2Y
        PZSUM=PQ1Z+PA1Z+PQ2Z+PA2Z
        PESUM=PQ1E+PA1E+PQ2E+PA2E
        WRITE(6,'(A)') ' CORMOM: KINEMATIC TEST FOR PARTONS'
        WRITE(6,'(A,4(1PE12.5))') ' PXSUM,PYSUM,PZSUM,PESUM', PXSUM,
     +  PYSUM,PZSUM,PESUM
      ENDIF
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE SELPT4( PTXSQ1,PTYSQ1,
     +PLQ1,EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,
     +PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1,NSELPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C  SELECT PT VALUES FOR A TWO CHAIN SYSTEM
C                            SELECT SEA QUARK AND ANTIQUARK PT-VALUES
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEEP,DROPPT.
      LOGICAL INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA, IPADIS,
     +ISHMAL,LPAULI
      COMMON /DROPPT/ INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA,
     +IPADIS,ISHMAL,LPAULI
*KEND.
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
C--------------------------------
C                                        change j.r.6.5.93
      QTXSQ1=PTXSQ1
      QTXSA1=PTXSA1
      QTXSQ2=PTXSQ2
      QTXSA2=PTXSA2
      QTYSQ1=PTYSQ1
      QTYSA1=PTYSA1
      QTYSQ2=PTYSQ2
      QTYSA2=PTYSA2
      QLQ1=PLQ1
      QLAQ1=PLAQ1
      QLQ2=PLQ2
      QLAQ2=PLAQ2
      QEQ1=EQ1
      QEAQ1=EAQ1
      QEQ2=EQ2
      QEAQ2=EAQ2
C                                        ----------------
      IANFA=0
      ITAGPT=0
C                           changed from 3.  j.r.21.8.93
      B33=16.00
      IF (IKVALA.EQ.1)B33=16.0
      ES=-2./(B33**2)*LOG(ABS(RNDM(V)*RNDM(U))+0.00000001)
      HPS=SQRT(ES*ES+2.*ES*0.94)
C............................................................
      IF (.NOT.INTPT) HPS=0.0000001
      ICOUNT=0
      IREJ=0
   10 CONTINUE
      ICOUNT=ICOUNT+1
      IF (ICOUNT.EQ.48)THEN
	HPS=0.D0
      ENDIF
      IF (ICOUNT.EQ.50)THEN
C                            REJECT EVENT
        IREJ=1
        RETURN
      ENDIF
      IF (ICOUNT.GE.1)THEN
        HPS=HPS*0.8
        CALL DSFECF(SFE,CFE)
        PTXSQ1=QTXSQ1+HPS*CFE
        PTYSQ1=QTYSQ1+HPS*SFE
        PTXSA1=QTXSA1-HPS*CFE
        PTYSA1=QTYSA1-HPS*SFE
        GO TO 111       
      ENDIF
      B33=2.*B33
      ES=-2./(B33**2)*LOG(ABS(RNDM(V)*RNDM(U))+0.00000001)
      HPS=SQRT(ES*ES+2.*ES*0.94)
C............................................................
  110 CONTINUE
      IF (.NOT.INTPT) HPS=0.0000001
C.............................................................
      CALL DSFECF(SFE,CFE)
C                                       change j.r.6.5.93
      PTXSQ1=QTXSQ1+HPS*CFE
      PTYSQ1=QTYSQ1+HPS*SFE
      PTXSA1=QTXSA1-HPS*CFE
      PTYSA1=QTYSA1-HPS*SFE
  111 CONTINUE
C                                      -----------------
      IF (IPEV.GE.7)WRITE(6,1000)PTXSQ1,PTYSQ1,PTXSA1,PTYSA1 ,PTXSQ2,
     +PTYSQ2,PTXSA2,PTYSA2
 1000 FORMAT (' PT S  ',8F12.6)
C                           KINEMATICS OF THE TWO CHAINS Q1-AQ2,AQ1-Q2
      PTTQ1=PTXSQ1**2+PTYSQ1**2
      PTTA1=PTXSA1**2+PTYSA1**2
      IF((EQ1**2.LE.PTTQ1).OR. (EAQ1**2.LE.PTTA1))            GO TO 10
C                                      
      IANFA2=0
      ITAGP2=0
      B33=16.00
      IF (IKVALA.EQ.1)B33=16.0
      ES=-2./(B33**2)*LOG(ABS(RNDM(V)*RNDM(U))+0.00000001)
      HPS=SQRT(ES*ES+2.*ES*0.94)
C............................................................
      IF (.NOT.INTPT) HPS=0.0000001
      ICOUN2=0
      IREJ=0
   12 CONTINUE
      ICOUN2=ICOUN2+1
      IF (ICOUN2.EQ.48)THEN
	HPS=0.D0
      ENDIF
      IF (ICOUN2.EQ.50)THEN
        IREJ=1
C                            REJECT EVENT
        RETURN
      ENDIF
      IF(ICOUN2.GE.1)THEN
        HPS=HPS*0.8
        CALL DSFECF(SFE,CFE)
        PTXSQ2=QTXSQ2+HPS*CFE
        PTYSQ2=QTYSQ2+HPS*SFE
        PTXSA2=QTXSA2-HPS*CFE
        PTYSA2=QTYSA2-HPS*SFE
        GO TO 113
      ENDIF
      B33=2.*B33
C
      ES=-2./(B33**2)*LOG(ABS(RNDM(V)*RNDM(U))+0.00000001)
      HPS=SQRT(ES*ES+2.*ES*0.94D0)
C............................................................
  112 CONTINUE
      IF (.NOT.INTPT) HPS=0.0000001
C.............................................................
      CALL DSFECF(SFE,CFE)
C                                       change j.r.6.5.93
      PTXSQ2=QTXSQ2+HPS*CFE
      PTYSQ2=QTYSQ2+HPS*SFE
      PTXSA2=QTXSA2-HPS*CFE
      PTYSA2=QTYSA2-HPS*SFE
  113 CONTINUE     
C                                      -----------------
C
      IF (IPEV.GE.7)WRITE(6,1000)PTXSQ1,PTYSQ1,PTXSA1,PTYSA1 ,PTXSQ2,
     +PTYSQ2,PTXSA2,PTYSA2
C                           KINEMATICS OF THE TWO CHAINS Q1-AQ2,AQ1-Q2
      PTTQ2=PTXSQ2**2+PTYSQ2**2
      PTTA2=PTXSA2**2+PTYSA2**2
      IF((EQ2**2.LE.PTTQ2).OR. (EAQ2**2.LE.PTTA2))            GO TO 12
C
C     IF(IP.GE.1)GO TO 1779
        PLQ1=SQRT(EQ1**2-PTTQ1)
        PLAQ1=SQRT(EAQ1**2-PTTA1)
        PLQ2=-SQRT(EQ2**2-PTTQ2)
        PLAQ2=-SQRT(EAQ2**2-PTTA2)
 1779 CONTINUE
C-----------
C                          CHAIN 1: Q1-AQ2     CHAIN2:  AQ1-Q2
      AMCH1Q=(EQ1+EAQ2)**2-(PTXSQ1+PTXSA2)** 2-(PTYSQ1+PTYSA2)**2-(PLQ1
     ++PLAQ2)**2
      IF (AMCH1Q.LE.0.D0)THEN
        WRITE(6,301)AMCH1Q
  301   FORMAT(' inconsistent Kinematics in SELPT AMCH1Q=',E12.3)
       WRITE(6,305) QTXSQ1,QTYSQ1,
     +QLQ1,QEQ1,QTXSA1,QTYSA1,QLAQ1,QEAQ1, QTXSQ2,QTYSQ2,QLQ2,QEQ2,
     +QTXSA2,
     +QTYSA2,QLAQ2,QEAQ2, AMCH1,AMCH2
  305 FORMAT( 'PTXSQ1,PTYSQ1,
     +PLQ1,EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,
     +PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2',5(4E15.5/))
        IREJ=1
        RETURN
      ENDIF
      AMCH1=SQRT(AMCH1Q)
C
C                          CHAIN 1: Q1-AQ2     CHAIN2:  AQ1-Q2
      AMCH2Q=(EQ2+EAQ1)**2-(PTXSQ2+PTXSA1)** 2-(PTYSQ2+PTYSA1)**2-(PLQ2
     ++PLAQ1)**2
      IF (AMCH2Q.LE.0.D0)THEN
        WRITE(6,302)AMCH2Q
  302   FORMAT(' inconsistent Kinematics in SELPT AMCH2Q=',E12.3)
       WRITE(6,305) QTXSQ1,QTYSQ1,
     +QLQ1,QEQ1,QTXSA1,QTYSA1,QLAQ1,QEAQ1, QTXSQ2,QTYSQ2,QLQ2,QEQ2,
     +QTXSA2,
     +QTYSA2,QLAQ2,QEAQ2, AMCH1,AMCH2
C        IF(ITAGPT.EQ.0)GO TO !33
        IREJ=1
        RETURN
      ENDIF
      AMCH2=SQRT(AMCH2Q)
      RETURN
      END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE SELPT( PTXSQ1,PTYSQ1,
     +PLQ1,EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,
     +PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1,
     * PTTQ2,PTTA2, NSELPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C  SELECT PT VALUES FOR A TWO CHAIN SYSTEM  DPMJET (combined
C                                           DTUNUC/TUJET method)
C                            SELECT SEA QUARK AND ANTIQUARK PT-VALUES
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEEP,DROPPT.
      LOGICAL INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA, IPADIS,
     +ISHMAL,LPAULI
      COMMON /DROPPT/ INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA,
     +IPADIS,ISHMAL,LPAULI
*KEND.
      COMMON /NNCMS/  GAMCM,BGCM,UMO,PCM,EPROJ,PPROJ
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
      data NUSEPT /0/
      data MUSEPT /0/
      IF(IPEV.GE.4)WRITE(6,6633) PTXSQ1,PTYSQ1,
     +PLQ1,EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,
     +PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1,
     *                 NSELPT
 6633 FORMAT(' selpt input: ',
     + ' PTXSQ1,PTYSQ1, +PLQ1,EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,
     + PTYSQ2,PLQ2,EQ2,PTXSA2, PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2,IREJ,
     +IKVALA, NSELPT ', 2(8E12.4/),2E12.4,2I5,2E12.4,I5)
C--------------------------------
C                                        change j.r.6.5.93
      QTXSQ1=PTXSQ1
      QTXSA1=PTXSA1
      QTXSQ2=PTXSQ2
      QTXSA2=PTXSA2
      QTYSQ1=PTYSQ1
      QTYSA1=PTYSA1
      QTYSQ2=PTYSQ2
      QTYSA2=PTYSA2
      QLQ1=PLQ1
      QLAQ1=PLAQ1
      QLQ2=PLQ2
      QLAQ2=PLAQ2
      QEQ1=EQ1
      QEAQ1=EAQ1
      QEQ2=EQ2
      QEAQ2=EAQ2
C                                        ----------------
      IANFA=0
      ITAGPT=0
      IANFA2=0
      ITAGP2=0
      IREJ=0
      ICOUNT=0
      ICOUN2=0
    1 CONTINUE
      IF ( NSELPT.EQ.0 .OR.UMO.LE.20.D0) THEN
C                           changed from 3.  j.r.21.8.93
        B33=16.0
C       IF (IKVALA.EQ.1)B33=3.7
C                                Test 12.2.96
C       IF (IKVALA.EQ.1)B33=(3.0+6.0/LOG10(UMO+10.))/2.
        IF (IKVALA.EQ.1)B33=3.0+6.0/LOG10(UMO+10.)
        IF (IKVALA.EQ.2)B33=4.0+3.0/LOG10(UMO+10.)
        ES=-2./(B33**2)*LOG(ABS(RNDM(V)*RNDM(U))+0.00000001)
        HPS=SQRT(ES*ES+2.*ES*0.94)
C............................................................
  110   CONTINUE
        IF (.NOT.INTPT) HPS=0.0000001
C.............................................................
   10   CONTINUE
        ICOUNT=ICOUNT+1
        IF (ICOUNT.EQ.48)THEN
	  HPS=0.D0
        ENDIF
        IF (ICOUNT.EQ.50)THEN
C                            REJECT EVENT
          IREJ=1
          RETURN
        ENDIF
        IF (ICOUNT.GE.2)THEN
          HPS=HPS*0.8
          PTXSQ1=PTXSQ1*0.8
          PTYSQ1=PTYSQ1*0.8
          PTXSA1=PTXSA1*0.8
          PTYSA1=PTYSA1*0.8
          CALL DSFECF(SFE,CFE)
C         PTXSQ1=QTXSQ1+HPS*CFE
C         PTYSQ1=QTYSQ1+HPS*SFE
C         PTXSA1=QTXSA1-HPS*CFE
C         PTYSA1=QTYSA1-HPS*SFE
          GO TO 111       
        ENDIF
        B33=2.*B33
        ES=-2./(B33**2)*LOG(ABS(RNDM(V)*RNDM(U))+0.00000001)
        HPS=SQRT(ES*ES+2.*ES*0.94)
C............................................................
        IF (.NOT.INTPT) HPS=0.0000001
C
      ELSEIF(NSELPT.EQ.1)THEN
          CALL SAMPPT(1,HPS)
          IF(IPEV.GE.4)WRITE(6,6638)HPS
      ELSEIF(NSELPT.EQ.2)THEN
        IF (NUSEPT.EQ.0)THEN
          CALL SAMPPT(1,HPS)
          IF(IPEV.GE.4)WRITE(6,6638)HPS
          NUSEPT=1
          USEPT=HPS
        ELSEIF(NUSEPT.EQ.1)THEN
          HPS=USEPT
        ENDIF
      ENDIF
      CALL DSFECF(SFE,CFE)
C                                       change j.r.6.5.93
      PTXSQ1=QTXSQ1+HPS*CFE
      PTYSQ1=QTYSQ1+HPS*SFE
      PTXSA1=QTXSA1-HPS*CFE
      PTYSA1=QTYSA1-HPS*SFE
      QTXSQ1=QTXSQ1*0.8
      QTYSQ1=QTYSQ1*0.8
      QTXSA1=QTXSA1*0.8
      QTYSA1=QTYSA1*0.8
  111 CONTINUE
C                                      -----------------
C
      IF (IPEV.GE.7)WRITE(6,1000)PTXSQ1,PTYSQ1,PTXSA1,PTYSA1 ,PTXSQ2,
     +PTYSQ2,PTXSA2,PTYSA2
 1000 FORMAT (' PT S  ',8F12.6)
C                           KINEMATICS OF THE TWO CHAINS Q1-AQ2,AQ1-Q2
      PTTQ1=PTXSQ1**2+PTYSQ1**2
      PTTA1=PTXSA1**2+PTYSA1**2
C                                      
      IF ( NSELPT.EQ.0.OR.UMO.LE.20.D0 ) THEN
        B33=16.0
C       IF (IKVALA.EQ.1)B33=3.7
        IF (IKVALA.EQ.1)B33=3.0+6.0/LOG10(UMO+10.)
        IF (IKVALA.EQ.2)B33=4.0+3.0/LOG10(UMO+10.)
        IREJ=0
   12   CONTINUE
        ICOUN2=ICOUN2+1
        IF (ICOUN2.EQ.48)THEN
	  HPS=0.D0
        ENDIF
        IF (ICOUN2.EQ.50)THEN
C                            REJECT EVENT
          IREJ=1
          RETURN
        ENDIF
        IF(ICOUN2.GE.2)THEN
          HPS=HPS*0.8
          PTXSQ2=PTXSQ2*0.8
          PTYSQ2=PTYSQ2*0.8
          PTXSA2=PTXSA2*0.8
          PTYSA2=PTYSA2*0.8
          CALL DSFECF(SFE,CFE)
C         PTXSQ2=QTXSQ2+HPS*CFE
C         PTYSQ2=QTYSQ2+HPS*SFE
C         PTXSA2=QTXSA2-HPS*CFE
C         PTYSA2=QTYSA2-HPS*SFE
          GO TO 113
        ENDIF
        B33=2.*B33
C
        ES=-2./(B33**2)*LOG(ABS(RNDM(V)*RNDM(U))+0.00000001)
        HPS=SQRT(ES*ES+2.*ES*0.94)
C............................................................
  112   CONTINUE
        IF (.NOT.INTPT) HPS=0.0000001
C.............................................................
      ELSEIF(NSELPT.EQ.1)THEN
        IF (MUSEPT.EQ.0)THEN
          CALL SAMPPT(1,HPS)
          IF(IPEV.GE.4)WRITE(6,6638)HPS
 6638     FORMAT (' SELPT:SAMPPT: HPS= ',E12.4)
          MUSEPT=1
          USEPTM=HPS
        ELSEIF(MUSEPT.EQ.1)THEN
          HPS=USEPTM
        ENDIF
      ENDIF
      CALL DSFECF(SFE,CFE)
C                                       change j.r.6.5.93
      PTXSQ2=QTXSQ2+HPS*CFE
      PTYSQ2=QTYSQ2+HPS*SFE
      PTXSA2=QTXSA2-HPS*CFE
      PTYSA2=QTYSA2-HPS*SFE
      QTXSQ2=QTXSQ2*0.8
      QTYSQ2=QTYSQ2*0.8
      QTXSA2=QTXSA2*0.8
      QTYSA2=QTYSA2*0.8
  113 CONTINUE     
C                                      -----------------
C
      IF (IPEV.GE.7)WRITE(6,1000)PTXSQ1,PTYSQ1,PTXSA1,PTYSA1 ,PTXSQ2,
     +PTYSQ2,PTXSA2,PTYSA2
C                           KINEMATICS OF THE TWO CHAINS Q1-AQ2,AQ1-Q2
      PTTQ1=PTXSQ1**2+PTYSQ1**2
      PTTA1=PTXSA1**2+PTYSA1**2
      PTTQ2=PTXSQ2**2+PTYSQ2**2
      PTTA2=PTXSA2**2+PTYSA2**2
      PTWQ1=SQRT(PTTQ1)
      PTWA1=SQRT(PTTA1)
      PTWQ2=SQRT(PTTQ2)
      PTWA2=SQRT(PTTA2)
      IF(PLQ1.GT.PTWQ1.AND.ABS(PLAQ2).GT.PTWQ1)THEN
	PLQ1=QLQ1-PTWQ1
	PLAQ2=QLAQ2+PTWQ1
      ELSEIF(PLQ1.GT.PTWA2.AND.ABS(PLAQ2).GT.PTWA2)THEN
	PLQ1=QLQ1-PTWA2
	PLAQ2=QLAQ2+PTWA2
      ENDIF
      IF(PLAQ1.GT.PTWA1.AND.ABS(PLQ2).GT.PTWA1)THEN
	PLAQ1=QLAQ1-PTWA1
	PLQ2=QLQ2+PTWA1
      ELSEIF(PLAQ1.GT.PTWQ2.AND.ABS(PLQ2).GT.PTWQ2)THEN
	PLAQ1=QLAQ1-PTWQ2
	PLQ2=QLQ2+PTWQ2
      ENDIF
      QLQ1=PLQ1
      QLAQ1=PLAQ1
      QLQ2=PLQ2
      QLAQ2=PLAQ2
      PTTQ1=PTTQ1+PLQ1**2
      PTTA1=PTTA1+PLAQ1**2
      PTTQ2=PTTQ2+PLQ2**2
      PTTA2=PTTA2+PLAQ2**2
        IF (INTPT) THEN
      AMTE1=0.2
      IF(AMTE1.GE.EQ1**2)AMTE1=EQ1**2/2.
      AMTE2=0.2
      IF(AMTE2.GE.EQ2**2)AMTE2=EQ2**2/2.
      AMTE3=0.2
      IF(AMTE1.GE.EAQ1**2)AMTE1=EAQ1**2/2.
      AMTE4=0.2
      IF(AMTE2.GE.EAQ2**2)AMTE2=EAQ2**2/2.
        IF((EQ1**2-AMTE1.LE.PTTQ1).OR.
     *	(EQ2**2-AMTE1.LE.PTTQ2)
     *  .OR.(EAQ1**2-AMTE3.LE.PTTA1).OR.
     *  (EAQ2**2-AMTE4.LE.PTTA2))THEN
          IF ( NSELPT.EQ.0.OR.UMO.LE.20.D0 ) THEN
            GO TO 1
          ELSE
            USEPT  = USEPT  * 0.7
            USEPTM = USEPTM * 0.7
            IF( USEPT.GT.0.01D0 .OR. USEPTM.GT.0.01D0 ) THEN
              IF (IPEV.GE.6)THEN
                WRITE(6,*)'  SELPT: JUMP AFTER REDUCTION OF USEPT'
                WRITE(6,*)'  SELPT: USEPT,USEPTM,HPS',USEPT,USEPTM,HPS
              ENDIF
              GO TO 1
            ELSE
              IREJ = 1
      IF(IPEV.GE.4)WRITE(6,6634) PTXSQ1,PTYSQ1,
     +PLQ1,EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,
     +PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1,
     *                 NSELPT
 6634 FORMAT(' selpt rejec: ',
     + ' PTXSQ1,PTYSQ1, +PLQ1,EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,
     + PTYSQ2,PLQ2,EQ2,PTXSA2, PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2,IREJ,
     +IKVALA, NSELPT ', 2(8E12.4/),2E12.4,3I5)
              RETURN
            ENDIF
          ENDIF
        ENDIF
        ENDIF
        NUSEPT=0
        MUSEPT=0
C
      IF(IP.GE.1)GO TO 1779
        QQQ1=QTXSQ1**2+QTYSQ1**2+QLQ1**2-PTTQ1
        IF(QQQ1.GT.0.D0)THEN
          PLQ1=SQRT(QQQ1)
        ELSE
          PLQ1=SQRT(EQ1**2-PTTQ1)
        ENDIF
        QQA1=QTXSA1**2+QTYSA1**2+QLAQ1**2-PTTA1
        IF(QQA1.GT.0.D0)THEN
          PLAQ1=SQRT(QQA1)
        ELSE
          PLAQ1=SQRT(EAQ1**2-PTTA1)
        ENDIF
        QQQ2=QTXSQ2**2+QTYSQ2**2+QLQ2**2-PTTQ2
        IF(QQQ2.GT.0.D0)THEN
          PLQ2=-SQRT(QQQ2)
        ELSE
          PLQ2=-SQRT(EQ2**2-PTTQ2)
        ENDIF
        QQA2=QTXSA2**2+QTYSA2**2+QLAQ2**2-PTTA2
        IF(QQA2.GT.0.D0)THEN
          PLAQ2=-SQRT(QQA2)
        ELSE
          PLAQ2=-SQRT(EAQ2**2-PTTA2)
        ENDIF
 1779 CONTINUE
C-----------
C-----------
C                          CHAIN 1: Q1-AQ2     CHAIN2:  AQ1-Q2
      AMCH1Q=(EQ1+EAQ2)**2-(PTXSQ1+PTXSA2)** 2-(PTYSQ1+PTYSA2)**2-(PLQ1
     ++PLAQ2)**2
      IF (AMCH1Q.LE.0.D0)THEN
C       IF(IANFA.EQ.0)THEN
C         IANFA=1
C         ITAGPT=1
C         GO TO 110
C       ENDIF
C       GO TO 10
        WRITE(6,301)AMCH1Q
  301   FORMAT(' inconsistent Kinematics in SELPT AMCH1Q=',E12.3)
       WRITE(6,305) QTXSQ1,QTYSQ1,
     +QLQ1,QEQ1,QTXSA1,QTYSA1,QLAQ1,QEAQ1, QTXSQ2,QTYSQ2,QLQ2,QEQ2,
     +QTXSA2,
     +QTYSA2,QLAQ2,QEAQ2, AMCH1,AMCH2
  305 FORMAT( 'PTXSQ1,PTYSQ1,
     +PLQ1,EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,
     +PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2',5(4E15.5/))
C       IF(ITAGPT.EQ.0)GO TO !33
        IREJ=1
      IF(IPEV.GE.4)WRITE(6,6635) PTXSQ1,PTYSQ1,
     +PLQ1,EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,
     +PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1,
     *                 NSELPT
 6635 FORMAT(' selpt rejec: ',
     + ' PTXSQ1,PTYSQ1, +PLQ1,EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,
     + PTYSQ2,PLQ2,EQ2,PTXSA2, PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2,IREJ,
     +IKVALA, NSELPT ', 2(8E12.4/),2E12.4,3I5)
        RETURN
      ENDIF
      AMCH1=SQRT(AMCH1Q)
C
C                          CHAIN 1: Q1-AQ2     CHAIN2:  AQ1-Q2
      AMCH2Q=(EQ2+EAQ1)**2-(PTXSQ2+PTXSA1)** 2-(PTYSQ2+PTYSA1)**2-(PLQ2
     ++PLAQ1)**2
      IF (AMCH2Q.LE.0.D0)THEN
C       IF(IANFA.EQ.0)THEN
C         IANFA=1
C         ITAGPT=1
C         GO TO 110
C       ENDIF
C       GO TO 10
        WRITE(6,302)AMCH2Q
  302   FORMAT(' inconsistent Kinematics in SELPT AMCH2Q=',E12.3)
       WRITE(6,305) QTXSQ1,QTYSQ1,
     +QLQ1,QEQ1,QTXSA1,QTYSA1,QLAQ1,QEAQ1, QTXSQ2,QTYSQ2,QLQ2,QEQ2,
     +QTXSA2,
     +QTYSA2,QLAQ2,QEAQ2, AMCH1,AMCH2
C        IF(ITAGPT.EQ.0)GO TO !33
        IREJ=1
      IF(IPEV.GE.4)WRITE(6,6636) PTXSQ1,PTYSQ1,
     +PLQ1,EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,
     +PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1,
     *                 NSELPT
 6636 FORMAT(' selpt rejec: ',
     + ' PTXSQ1,PTYSQ1, +PLQ1,EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,
     + PTYSQ2,PLQ2,EQ2,PTXSA2, PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2,IREJ,
     +IKVALA, NSELPT ', 2(8E12.4/),2E12.4,3I5)
        RETURN
      ENDIF
      AMCH2=SQRT(AMCH2Q)
      IF(IPEV.GE.4)WRITE(6,6637) PTXSQ1,PTYSQ1,
     +PLQ1,EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,
     +PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1,
     *                 NSELPT
 6637 FORMAT(' selpt exit : ',
     + ' PTXSQ1,PTYSQ1, +PLQ1,EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,
     + PTYSQ2,PLQ2,EQ2,PTXSA2, PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2,IREJ,
     +IKVALA, NSELPT ', 2(8E12.4/),2E12.4,2I5,2E12.4,I5)
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE DECHKK(NHKKH1)
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
*KEEP,DDECAC.
      PARAMETER (IDMAX9=602)
      CHARACTER*8  ZKNAME
      COMMON /DDECAC/ ZKNAME(IDMAX9),WT(IDMAX9),NZK(IDMAX9,3)
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
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
*KEND.
      DIMENSION ECMF(3),PCMF(3),CODF(3),COFF(3),SIFF(3),ITF(3)
      DIMENSION ECMFF(3),PCMFF(3)
      DIMENSION CXF(3),CYF(3),CZF(3)
      DATA ISTAB /2/
C-----------------------------------------------------
      IHKK=NHKKH1
C     IPHKK=2
      IF (IPHKK.GE.2) WRITE(6,1000) IHKK,NHKK
 1000 FORMAT(' DECHKK IHKK,NHKK= ',2I5)
C
C***                    LOOP OVER ALL PARTICLES FROM THE STACK
   10 CONTINUE
      IHKK=IHKK+1
      IF (IHKK.GT.NHKK)THEN
C       IPHKK=0
        RETURN
      ENDIF
C
      IF(ABS(ISTHKK(IHKK)).NE.1)       GOTO 10
      IQQQQ=ISTHKK(IHKK)
C
      IT=MCIHAD(IDHKK(IHKK))
C
      IF (IT.LT.1.OR.IT.GT.210) THEN
        WRITE (6,1003)IT
 1003   FORMAT (' DECHKK IT ',I10)
      ENDIF
C
C*****TEST STABLE OR UNSTABLE
C     ISTAB=2
C   ISTAB=1/2/3 MEANS  STRONG + WEAK DECAYS / ONLY STRONG DECAYS /
C   STRONG DECAYS + WEAK DECAYS FOR CHARMED PARTICLES AND TAU LEPTONS
C
C   GOTO 51 :             THERE WAS NO DECAY RETURN TO STACK
C
      IF(ISTAB.EQ.1) THEN
        IF(IT.EQ.135.OR.IT.EQ.136)                               GOTO 10
        IF(IT.GE.1.AND.IT.LE.7)                                  GOTO 10
      ELSEIF(ISTAB.EQ.2) THEN
        IF(IT.GE. 1.AND.IT.LE. 30)                               GOTO 10
        IF(IT.GE. 97.AND.IT.LE.103)                              GOTO 10
        IF(IT.GE.115.AND.IT.LE.122)                              GOTO 10
        IF(IT.GE.131.AND.IT.LE.136)                              GOTO 10
        IF(IT.EQ.109)                                            GOTO 10
        IF(IT.GE.137.AND.IT.LE.160)                              GOTO 10
      ELSEIF(ISTAB.EQ.3) THEN
        IF(IT.GE.1.AND.IT.LE.23)                                 GOTO 10
        IF(IT.GE. 97.AND.IT.LE.103)                              GOTO 10
        IF(IT.EQ.109.AND.IT.EQ.115)                              GOTO 10
        IF(IT.GE.133.AND.IT.LE.136)                              GOTO 10
      ENDIF
C***                               DECAY TO BE HANDLED
C
      PLS=SQRT(ABS(PHKK(1,IHKK)**2+PHKK(2,IHKK)**2+PHKK(3,IHKK)**2))
C
C                    Consistency check of decaying hadron
C
      AMTEST=SQRT(ABS(PHKK(4,IHKK)**2-PLS**2))
      IF(ABS(AMTEST-PHKK(5,IHKK)).GE.1.D-3)THEN
C	WRITE(6,'(A,2E15.5,I10)')' DECHKK inconsistent resonance',
C    *   AMTEST,PHKK(5,IHKK),IHKK
	PLSS=(PHKK(4,IHKK)**2-PHKK(5,IHKK))
	IF(PLSS.LE.0.D0)THEN
	  WRITE(6,'(A)')' negative momentum square!'
	  PLSS=0.D0
	ENDIF
	PLSN=SQRT(PLSS)
	AMODP=PLSN/PLS
	PHKK(1,IHKK)=PHKK(1,IHKK)*AMODP
	PHKK(2,IHKK)=PHKK(2,IHKK)*AMODP
	PHKK(3,IHKK)=PHKK(3,IHKK)*AMODP
	PLS=PLS*AMODP
      ENDIF
C
      IF(PLS.NE.0.D0) THEN
        CXS=PHKK(1,IHKK)/PLS
        CYS=PHKK(2,IHKK)/PLS
        CZS=PHKK(3,IHKK)/PLS
      ENDIF
      ELS=PHKK(4,IHKK)
      ECO=AAM(IT)
      GAM=ELS/ECO
      BGAM=PLS/ECO
C
      KZ1=K1(IT)
      VV=RNDM(V) - 1.E-17
      IIK=KZ1-1
   20 IIK=IIK+1
      IF (VV.GT.WT(IIK))                                        GO TO 20
C
C                                          IIK IS THE DECAY CHANNEL
      ITF(1)=NZK(IIK,1)
      ITF(2)=NZK(IIK,2)
C******************************   ?????????????????????
C??   IF (ITF(2)-1.LT.0) GO TO 110
C??   IF (IT2-1.LT.0) GO TO 305
      IF (ITF(2).LT.1)                                          GO TO 10
C******************************   ????????????????????
      ITF(3)=NZK(IIK,3)
C
      IF(IPHKK.GE.1) WRITE(6,1010)IT,IIK,ITF(1),ITF(2),ITF(3)
 1010 FORMAT(' DECHKK IT,IIK,IT1,IT2,IT3 ',5I5)
C  IT1,IT2, IT3 ARE THE PRODUCED PARTICLES FROM  IT
C
      IF(ITF(3).EQ.0) THEN
        NDECPR=2
        CALL DTWOPD(ECO,ECMF(1),ECMF(2),PCMF(1),PCMF(2), CODF(1),COFF
     +  (1),SIFF(1),CODF(2),COFF(2),SIFF(2), AAM(ITF(1)),AAM(ITF(2)))
        SID1=SQRT(ABS((1.-CODF(1))*(1.+CODF(1))))
        SID2=SQRT(ABS((1.-CODF(2))*(1.+CODF(2))))
        PIX1=PCMF(1)*SID1*COFF(1)
        PIY1=PCMF(1)*SID1*SIFF(1)
        PIZ1=PCMF(1)*CODF(1)
        PIX2=PCMF(2)*SID2*COFF(2)
        PIY2=PCMF(2)*SID2*SIFF(2)
        PIZ2=PCMF(2)*CODF(2)
        PIX12=PIX1+PIX2
        PIY12=PIY1+PIY2
        PIZ12=PIZ1+PIZ2
        ECM12=ECMF(1)+ECMF(2)-ECO
        IF((ABS(PIX12).GT.0.000001D0).OR.
     +     (ABS(PIY12).GT.0.000001D0).OR.
     +     (ABS(PIZ12).GT.0.000001D0).OR. 
     +     (ABS(ECM12).GT.0.000001D0))THEN
           WRITE(6,778)PIX12,PIY12,PIZ12,ECM12
  778      FORMAT(' DWOPD px,py,pz,e',4F10.6)
        ENDIF
 
      ELSE
        NDECPR=3
       CALL DTHREP(ECO,ECMF(1),ECMF(2),ECMF(3),PCMF(1),PCMF(2),PCMF(3),
     +  CODF(1),COFF(1),SIFF(1),CODF(2),COFF(2),SIFF(2), CODF(3),COFF
     +  (3),SIFF(3), AAM(ITF(1)),AAM(ITF(2)),AAM(ITF(3)))
        SID1=SQRT((1.-CODF(1))*(1.+CODF(1)))
        SID2=SQRT((1.-CODF(2))*(1.+CODF(2)))
        SID3=SQRT((1.-CODF(3))*(1.+CODF(3)))
        PIX1=PCMF(1)*SID1*COFF(1)
        PIY1=PCMF(1)*SID1*SIFF(1)
        PIZ1=PCMF(1)*CODF(1)
        PIX2=PCMF(2)*SID2*COFF(2)
        PIY2=PCMF(2)*SID2*SIFF(2)
        PIZ2=PCMF(2)*CODF(2)
        PIX3=PCMF(3)*SID3*COFF(3)
        PIY3=PCMF(3)*SID3*SIFF(3)
        PIZ3=PCMF(3)*CODF(3)
        PIX12=PIX1+PIX2+PIX3
        PIY12=PIY1+PIY2+PIY3
        PIZ12=PIZ1+PIZ2+PIZ3
        ECM12=ECMF(1)+ECMF(2)+ECMF(3)-ECO
        IF((ABS(PIX12).GT.0.000001D0).OR.
     +     (ABS(PIY12).GT.0.000001D0).OR.
     +     (ABS(PIZ12).GT.0.000001D0).OR. 
     +     (ABS(ECM12).GT.0.000001D0))THEN
           WRITE(6,779)PIX12,PIY12,PIZ12,ECM12
  779      FORMAT(' DTHEPD px,py,pz,e',4F10.6)
        ENDIF
 
      ENDIF
C
      JDAHKK(1,IHKK)=NHKK + 1
      JDAHKK(2,IHKK)=NHKK + NDECPR
      DO 30 ID=1,NDECPR
      EHECC=SQRT(ABS(PCMF(ID)** 2+ AAM(ITF(ID))**2))
        IF (ABS(EHECC-ECMF(ID)).GT.0.0001D0) THEN
             WRITE(6,'(2A/3I5,3E15.6)')
     &            ' DECHKK: CORRECT INCONSISTENT ENERGY ',
     *            '  IHKK,NHKK,ITF(ID), ECMF(ID),EHECC, AAM(ITF(ID))',
     *            IHKK,NHKK,ITF(ID), ECMF(ID),EHECC, AAM(ITF(ID))
        ENDIF
C      CALL DTRAFO(GAM,BGAM,CXS,CYS,CZS, CODF(ID),COFF(ID),
C    *SIFF(ID),PCMF
       CALL DTRAFO(GAM,BGAM,CXS,CYS,CZS, CODF(ID),COFF(ID),
     * SIFF(ID),PCMF(ID),ECMF(ID), PCMFF(ID),CXF(ID),
     *CYF(ID),CZF(ID),ECMFF(ID))
        IF (IPHKK.GE.2) WRITE(6,'(A,7E15.5/8E15.5)')' DTRAFO ',
     * GAM,BGAM,CXS,CYS,CZS, CODF(ID),COFF(ID),
     * SIFF(ID),PCMF(ID),ECMF(ID), PCMFF(ID),CXF(ID),
     *CYF(ID),CZF(ID),ECMFF(ID)
 
C*******************
C                                  THERE WAS A DECAY, DROP MOTHER IHKK
C                                  AND PUT PARTICLE INTO HKK STACK
        ISTHKK(IHKK)=2
        IF (NHKK.EQ.NMXHKK)THEN
          WRITE (6,1020)NHKK,NMXHKK
 1020 FORMAT (' NHKK.GT.NMXHKK IN DECHKK RETURN ',2I10)
C      	  IPHKK=0
          RETURN
        ENDIF
C
        NHKK=NHKK+1
        IF (NHKK.EQ.NMXHKK) THEN
          WRITE (6,'(A,2I5)') ' DECHKK: NHKK.EQ.NMXHKK ',NHKK,NMXHKK
C      	  IPHKK=0
          RETURN
        ENDIF
	IDBAM(NHKK)=ITF(ID)
        ISTHKK(NHKK)=IQQQQ
        IDHKK(NHKK)=MPDGHA(ITF(ID))
        JMOHKK(1,NHKK)=IHKK
        JMOHKK(2,NHKK)=0
        JDAHKK(1,NHKK)=0
        JDAHKK(2,NHKK)=0
        PHKK(1,NHKK)=CXF(ID)*PCMFF(ID)
        PHKK(2,NHKK)=CYF(ID)*PCMFF(ID)
        PHKK(3,NHKK)=CZF(ID)*PCMFF(ID)
      EHECC=SQRT(ABS(PCMFF(ID)** 2+ AAM(ITF(ID))**2))
        IF (ABS(EHECC-ECMFF(ID)).GT.0.003D0) THEN
             WRITE(6,'(2A/3I5,3E15.6)')
     &            ' DECHKK: CORRECT INCONSISTENT ENERGY ',
     *            '  IHKK,NHKK,ITF(ID), ECMF(ID),EHECC, AAM(ITF(ID))',
     *            IHKK,NHKK,ITF(ID), ECMFF(ID),EHECC, AAM(ITF(ID))
          ECMFF(ID)=EHECC
        ENDIF
        PHKK(4,NHKK)=ECMFF(ID)
      PHKK(5,NHKK)=AAM(ITF(ID))
        VHKK(1,NHKK)=VHKK(1,IHKK)
        VHKK(2,NHKK)=VHKK(2,IHKK)
        VHKK(3,NHKK)=VHKK(3,IHKK)
        VHKK(4,NHKK)=VHKK(4,IHKK)
C
        IF (IPHKK.GE.7) WRITE(6,1030)NHKK, ISTHKK(NHKK),IDHKK(NHKK),
     +  JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +  (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
 1030 FORMAT (I6,I4,5I6,9E10.2)
C
   30 CONTINUE
C
                                                                 GOTO 10
C
C     RETURN
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

*=== trafo ============================================================*        
*                                                                               
      SUBROUTINE DTRAFO(GAM,BGAM,CX,CY,CZ,COD,COF,SIF,P,ECM,                     
     1PL,CXL,CYL,CZL,EL)                                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C     LORENTZ TRANSFORMATION INTO THE LAB - SYSTEM                              
      SID=SQRT(1.D0-COD*COD)                                                    
      PLX=P*SID*COF                                                             
      PLY=P*SID*SIF   
      PPPT=SQRT(PLX**2+PLY**2)
      PCMZ=P*COD                                                                
      PLZ=GAM*PCMZ+BGAM*ECM                                                     
      PL=SQRT(PLX*PLX+PLY*PLY+PLZ*PLZ)                                          
      EL=GAM*ECM+BGAM*PCMZ                                                      
C     ROTATION INTO THE ORIGINAL DIRECTION                                      
      COZ=PLZ/PL                                                                
C     SIZ=SQRT((1.D0-COZ)*(1.D0+COZ))   
      SIZ=PPPT/PL
C     WRITE(6,'(A,2E25.16)')' COZ,SIZ ',COZ,SIZ
      CALL STTRAN(CX,CY,CZ,COZ,SIZ,SIF,COF,CXL,CYL,CZL)                         
      RETURN                                                                    
      END        
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                                               
      SUBROUTINE STTRAN(XO,YO,ZO,CDE,SDE,SFE,CFE,X,Y,Z)                         
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      DATA ANGLSQ/1.D-14/                                                        
C********************************************************************           
C     VERSION BY                     J. RANFT                                   
C                                    LEIPZIG                                    
C                                                                               
C     THIS IS A SUBROUTINE OF FLUKA TO GIVE NEW DIRECTION COSINES               
C                                                                               
C     INPUT VARIABLES:                                                          
C        XO,YO,ZO = ORIGINAL DIRECTION COSINES                                  
C        CDE,SDE  = COSINE AND SINE OF THE POLAR (THETA)                        
C                   ANGLE OF "SCATTERING"                                       
C        SDE      = SINE OF THE POLAR (THETA) ANGLE OF "SCATTERING"             
C        SFE,CFE  = SINE AND COSINE OF THE AZIMUTHAL (PHI) ANGLE                
C                   OF "SCATTERING"                                             
C                                                                               
C     OUTPUT VARIABLES:                                                         
C        X,Y,Z     = NEW DIRECTION COSINES                                      
C                                                                               
C     ROTATION OF COORDINATE SYSTEM (SEE CERN 64-47 )                           
C********************************************************************           
C                                                                               
*                                                                               
*  Changed by A. Ferrari                                                        
*                                                                               
*     IF (ABS(XO)-0.0001D0) 1,1,2                                               
*   1 IF (ABS(YO)-0.0001D0) 3,3,2                                               
*   3 CONTINUE                                                                  
      A = XO**2 + YO**2                                                         
      IF ( A .LT. ANGLSQ ) THEN                                                 
         X=SDE*CFE                                                              
         Y=SDE*SFE                                                              
C     Z=CDE CORRECTED AUGUST 88 PA                                              
         Z=CDE*ZO                                                               
      ELSE                                                                      
         XI=SDE*CFE                                                             
         YI=SDE*SFE                                                             
         ZI=CDE                                                                 
         A=SQRT(A)                                                              
         X=-YO*XI/A-ZO*XO*YI/A+XO*ZI                                            
         Y=XO*XI/A-ZO*YO*YI/A+YO*ZI                                             
         Z=A*YI+ZO*ZI                                                           
      END IF                                                                    
      RETURN                                                                    
      END                                                                       
*-- Author :
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE LORTMO(N,GAM,BGX,BGY,BGZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C*** LORENTZ TRANSFORMATION OF THE N PARTICLES IN  FINPAR
C
*KEEP,DFINPA.
      CHARACTER*8 ANF
      PARAMETER (NFIMAX=249)
      COMMON /DFINPA/ ANF(NFIMAX),PXF(NFIMAX),PYF(NFIMAX),PZF(NFIMAX),
     +HEF(NFIMAX),AMF(NFIMAX), ICHF(NFIMAX),IBARF(NFIMAX),NREF(NFIMAX)
      COMMON /DFINPZ/IORMO(NFIMAX),IDAUG1(NFIMAX),IDAUG2(NFIMAX),
     * ISTATH(NFIMAX)
*KEND.
C-------------------
      PARAMETER (TINY=1.D-10)
      DATA IFIRST/0/
      DATA NUM/0/
C
      IFIRST=IFIRST+1
      NUM=NUM+1
      PXSM=0.0
      PYSM=0.0
      PZSM=0.0
      ESUM=0.
      PXSC=0.0
      PYSC=0.0
      PZSC=0.0
      ESMC=0.0
C           END OF CHANGE
      DO 10 I=1,N
        PXI=PXF(I)
        PYI=PYF(I)
        PZI=PZF(I)
      EEI=HEF(I)
        PXSM=PXSM + PXI
        PYSM=PYSM + PYI
        PZSM=PZSM + PZI
        ESUM=ESUM + EEI
        CALL DALTRA(GAM,BGX,BGY,BGZ,PXI,PYI,PZI,EEI, PPA,PXF(I),PYF(I),
     +  PZF(I),HEF(I))
        PXSC=PXSC + PXF(I)
        PYSC=PYSC + PYF(I)
        PZSC=PZSC + PZF(I)
      ESMC=ESMC + HEF(I)
   10 CONTINUE
C
C     PXSM,ETC,ARE SUMS FOR BAMJET FRAGMENTS IN JET CMS
      CALL DALTRA(GAM,BGX,BGY,BGZ,PXSM,PYSM,PZSM,ESUM, PPA,PXSM,PYSM,
     +PZSM,ESUM)
C
C     PXSC,ETC,ARE SUMS FOR BAMJET FRAGMENTS IN PROJ,TARGET CMS
 
      PXDIF=PXSM-PXSC
      PYDIF=PYSM-PYSC
      PZDIF=PZSM-PZSC
      EDIF=ESUM-ESMC
      DIFFL=PXDIF+PYDIF+PZDIF+EDIF
      IF(ESUM.LT.TINY)ESUM=TINY
      DIFFL=DIFFL/ESUM
      IF(DIFFL.GE.1.D-4)WRITE(6,1000)NUM,PXDIF,PYDIF,PZDIF,EDIF,PXSM,
     +PXSC, PYSM,PYSC,PZSM,PZSC,ESUM,ESMC
 1000 FORMAT(' ',2X,'LORTRA:NUM=',I5,2X,'PXDIF=',1PE15.6,2X,'PYDIF=', 1
     +PE15.6,2X,'PZDIF=',1PE15.6,2X,'EDIF=',1PE15.6/2X,'PXSM=',1PE15.6,2
     +X,'PXSC=',1PE15.6,2X,'PYSM=',1PE15.6,2X,'PYSC=',1PE15.6/2X,'PZSM',
     +1PE15.6,2X,'PZSC=',1PE15.6,2X,'ESUM=',1PE15.6,2X,'ESMC=',1PE15.6/2
     +X,'LORTRA DIFFERENCES DUE TO ALTRA'/)
      RETURN
      END
*-- Author :
C-------------------------------------------------------------------
C
C                              FILE DMNUC3.FOR
C
C-------------------------------------------------------------------
C
      SUBROUTINE EVTEST(IREJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C            TEST OF ENERGY MOMENTUM CONSERVATION ON NIVEAU OF CHAINS
C                                      AND ON NIVEAU OF CHAIN  ENDS
C
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
*KEEP,NUCIMP.
      COMMON /NUCIMP/ PRMOM(5,248),TAMOM(5,248), PRMFEP,PRMFEN,TAMFEP,
     +TAMFEN, PREFEP,PREFEN,TAEFEP,TAEFEN, PREPOT(210),TAEPOT(210),
     +PREBIN,TAEBIN,FERMOD,ETACOU
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
*KEEP,INTMX.
      PARAMETER (INTMX=2488,INTMD=252)
*KEEP,DXQX.
C     INCLUDE (XQXQ)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON /DXQX/ XPVQ(248),XPVD(248),XTVQ(248),XTVD(248), XPSQ
     +(INTMX),XPSAQ(INTMX),XTSQ(INTMX),XTSAQ(INTMX)
     *                ,XPSU(248),XTSU(248)
     *                ,XPSUT(248),XTSUT(248)
      COMMON /INTNEZ/ NDZ,NZD
*KEEP,INTNEW.
      COMMON /INTNEW/ NVV,NSV,NVS,NSS,NDV,NVD,NDS,NSD,
     +IXPV,IXPS,IXTV,IXTS, INTVV1(248),
     +INTVV2(248),INTSV1(248),INTSV2(248), INTVS1(248),INTVS2(248),
     +INTSS1(INTMX),INTSS2(INTMX),
     +INTDV1(248),INTDV2(248),INTVD1(248),INTVD2(248),
     +INTDS1(INTMD),INTDS2(INTMD),INTSD1(INTMD),INTSD2(INTMD)
 
C  /INTNEW/
C       NVV    :   NUMBER OF INTERACTING VALENCE-VALENCE SYSTEMS
C       NSV    :   NUMBER OF INTERACTING SEA-VALENCE SYSTEMS
C       NVS    :   NUMBER OF INTERACTING VALENCE-SEA SYSTEMS
C       NSS    :   NUMBER OF INTERACTING SEA-SEA SYSTEMS
C       IXPV, IXTV : NUMBER OF GENERATED X-VALUES FOR VALENCE-QUARK
C                     SYSTEMS FROM PROJECTILE/TARGET NUCLEI
C       IXPS, IXTS : NUMBER OF GENERATED X-VALUES FOR SEA-QUARK PAIRS
C                     FROM PROJECTILE/TARGET NUCLEI
C-------------------
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,LOZUO.
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/ ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),ZUOST(INTMX),
     +INTLO(INTMX),INLOSS(INTMX)
C  /LOZUO/
C       INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C------------------
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
*KEEP,ABRSS.
C     INCLUDE (ABRSS)
      COMMON /ABRSS/ AMCSS1(INTMX),AMCSS2(INTMX), GACSS1(INTMX),GACSS2
     +(INTMX), BGXSS1(INTMX),BGYSS1(INTMX),BGZSS1(INTMX), BGXSS2(INTMX),
     +BGYSS2(INTMX),BGZSS2(INTMX), NCHSS1(INTMX),NCHSS2(INTMX), IJCSS1
     +(INTMX),IJCSS2(INTMX), PQSSA1(INTMX,4),PQSSA2(INTMX,4), PQSSB1
     +(INTMX,4),PQSSB2(INTMX,4)
*KEEP,NNCMS.
      COMMON /NNCMS/  GAMCM,BGCM,UMO,PCM,EPROJ,PPROJ
*KEEP,ABRSV.
      COMMON /ABRSV/ AMCSV1(248),AMCSV2(248),GACSV1(248),GACSV2(248),
     +BGXSV1(248),BGYSV1(248),BGZSV1(248), BGXSV2(248),BGYSV2(248),
     +BGZSV2(248), NCHSV1(248),NCHSV2(248),IJCSV1(248),IJCSV2(248),
     +PQSVA1(248,4),PQSVA2(248,4), PQSVB1(248,4),PQSVB2(248,4)
*KEEP,ABRVS.
      COMMON /ABRVS/ AMCVS1(248),AMCVS2(248),GACVS1(248),GACVS2(248),
     +BGXVS1(248),BGYVS1(248),BGZVS1(248), BGXVS2(248),BGYVS2(248),
     +BGZVS2(248), NCHVS1(248),NCHVS2(248),IJCVS1(248),IJCVS2(248),
     +PQVSA1(248,4),PQVSA2(248,4), PQVSB1(248,4),PQVSB2(248,4)
*KEEP,ABRVV.
      COMMON /ABRVV/ AMCVV1(248),AMCVV2(248),GACVV1(248),GACVV2(248),
     +BGXVV1(248),BGYVV1(248),BGZVV1(248), BGXVV2(248),BGYVV2(248),
     +BGZVV2(248), NCHVV1(248),NCHVV2(248),IJCVV1(248),IJCVV2(248),
     +PQVVA1(248,4),PQVVA2(248,4), PQVVB1(248,4),PQVVB2(248,4)
*KEEP,ABRDV.
      COMMON /ABRDV/ AMCDV1(248),AMCDV2(248),GACDV1(248),GACDV2(248),
     +BGXDV1(248),BGYDV1(248),BGZDV1(248), BGXDV2(248),BGYDV2(248),
     +BGZDV2(248), NCHDV1(248),NCHDV2(248),IJCDV1(248),IJCDV2(248),
     +PQDVA1(248,4),PQDVA2(248,4), PQDVB1(248,4),PQDVB2(248,4)
C-------------------
*KEEP,ABRVD.
      COMMON /ABRVD/ AMCVD1(248),AMCVD2(248),GACVD1(248),GACVD2(248),
     +BGXVD1(248),BGYVD1(248),BGZVD1(248), BGXVD2(248),BGYVD2(248),
     +BGZVD2(248), NCHVD1(248),NCHVD2(248),IJCVD1(248),IJCVD2(248),
     +PQVDA1(248,4),PQVDA2(248,4), PQVDB1(248,4),PQVDB2(248,4)
*KEEP,ABRDS.
      COMMON /ABRDS/ AMCDS1(248),AMCDS2(248),GACDS1(248),GACDS2(248),
     +BGXDS1(248),BGYDS1(248),BGZDS1(248), BGXDS2(248),BGYDS2(248),
     +BGZDS2(248), NCHDS1(248),NCHDS2(248),IJCDS1(248),IJCDS2(248),
     +PQDSA1(248,4),PQDSA2(248,4), PQDSB1(248,4),PQDSB2(248,4)
C-------------------
*KEEP,ABRDS.
       COMMON /ABRDZ/ AMCDZ1(INTMD),AMCDZ2(INTMD),
     +GACDZ1(INTMD),GACDZ2(INTMD),
     +BGXDZ1(INTMD),BGYDZ1(INTMD),BGZDZ1(INTMD),
     +BGXDZ2(INTMD),BGYDZ2(INTMD),
     +BGZDZ2(INTMD), NCHDZ1(INTMD),NCHDZ2(INTMD),
     +IJCDZ1(INTMD),IJCDZ2(INTMD),
     +PQDZA1(INTMD,4),PQDZA2(INTMD,4),
     +PQDZB1(INTMD,4),PQDZB2(INTMD,4),
     +IPZQ(INTMD),IPZQQ2(INTMD),ITZQ(INTMD),
     +IPZAQ(INTMD),IZAQQ2(INTMD),ITZAQ(INTMD)
     +,IDZZZ(INTMD)
C-------------------
*KEEP,ABRSD.
      COMMON /ABRSD/ AMCSD1(248),AMCSD2(248),GACSD1(248),GACSD2(248),
     +BGXSD1(248),BGYSD1(248),BGZSD1(248), BGXSD2(248),BGYSD2(248),
     +BGZSD2(248), NCHSD1(248),NCHSD2(248),IJCSD1(248),IJCSD2(248),
     +PQSDA1(248,4),PQSDA2(248,4), PQSDB1(248,4),PQSDB2(248,4)
C-------------------
*KEEP,ABRSD.
      COMMON /ABRZD/ AMCZD1(INTMD),AMCZD2(INTMD),
     +GACZD1(INTMD),GACZD2(INTMD),
     +BGXZD1(INTMD),BGYZD1(INTMD),BGZZD1(INTMD),
     +BGXZD2(INTMD),BGYZD2(INTMD),
     +BGZZD2(INTMD), NCHZD1(INTMD),NCHZD2(INTMD),
     +IJCZD1(INTMD),IJCZD2(INTMD),
     +PQZDA1(INTMD,4),PQZDA2(INTMD,4),
     +PQZDB1(INTMD,4),PQZDB2(INTMD,4),
     +IPYQ(INTMD),ITYQ(INTMD),ITYQ2(INTMD),
     +IPYAQ(INTMD),ITYAQ(INTMD),ITYAQ2(INTMD)
     +,IZDYY(INTMD)
C-------------------
*KEEP,DROPPT.
      LOGICAL INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA, IPADIS,
     +ISHMAL,LPAULI
      COMMON /DROPPT/ INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA,
     +IPADIS,ISHMAL,LPAULI
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEND.
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
      COMMON /ABRZZ/ AMCZZ1(INTMX),AMCZZ2(INTMX),
     *               GACZZ1(INTMX),GACZZ2(INTMX),
     *               BGXZZ1(INTMX),BGYZZ1(INTMX),BGZZZ1(INTMX),
     *               BGXZZ2(INTMX),BGYZZ2(INTMX),BGZZZ2(INTMX),
     *               NCHZZ1(INTMX),NCHZZ2(INTMX),
     *               IJCZZ1(INTMX),IJCZZ2(INTMX),
     *               PQZZA1(INTMX,4),PQZZA2(INTMX,4),
     *               PQZZB1(INTMX,4),PQZZB2(INTMX,4)
      COMMON /ABRHH/ AMCHH1(INTMX),AMCHH2(INTMX),
     *               GACHH1(INTMX),GACHH2(INTMX),
     *               BGXHH1(INTMX),BGYHH1(INTMX),BGZHH1(INTMX),
     *               BGXHH2(INTMX),BGYHH2(INTMX),BGZHH2(INTMX),
     *               NCHHH1(INTMX),NCHHH2(INTMX),
     *               IJCHH1(INTMX),IJCHH2(INTMX),
     *               PQHHA1(INTMX,4),PQHHA2(INTMX,4),
     *               PQHHB1(INTMX,4),PQHHB2(INTMX,4)
      COMMON /NUCJTN/NONUJ1,NONUJT,NONUS1,NONUST
      COMMON /ABRJT/XJQ1(INTMX),XJAQ1(INTMX),XJQ2(INTMX),XJAQ2(INTMX),
     *        IJJQ1(INTMX),IJJAQ1(INTMX),IJJQ2(INTMX),IJJAQ2(INTMX),
     *        AMJCH1(INTMX),AMJCH2(INTMX),GAMJH1(INTMX),GAMJH2(INTMX),
     *        BGJH1(INTMX),BGJH2(INTMX),THEJH1(INTMX),THEJH2(INTMX),
     *        BGXJH1(INTMX),BGYJH1(INTMX),BGZJH1(INTMX),
     *        BGXJH2(INTMX),BGYJH2(INTMX),BGZJH2(INTMX),
     *  PJETA1(INTMX,4),PJETA2(INTMX,4),PJETB1(INTMX,4),PJETB2(INTMX,4)
     * ,JHKKPH(INTMX),JHKKTH(INTMX),JHKKEX(INTMX),JHKKE1(INTMX)
      COMMON /ABRSOF/XSQ1(INTMX),XSAQ1(INTMX),XSQ2(INTMX),XSAQ2(INTMX),
     *        IJSQ1(INTMX),IJSAQ1(INTMX),IJSQ2(INTMX),IJSAQ2(INTMX),
     *        AMCCH1(INTMX),AMCCH2(INTMX),GAMCH1(INTMX),GAMCH2(INTMX),
     *        BGCH1(INTMX),BGCH2(INTMX),THECH1(INTMX),THECH2(INTMX),
     *        BGXCH1(INTMX),BGYCH1(INTMX),BGZCH1(INTMX),
     *        BGXCH2(INTMX),BGYCH2(INTMX),BGZCH2(INTMX),
     *        NCH1(INTMX),NCH2(INTMX),IJCH1(INTMX),IJCH2(INTMX),
     *  PSOFA1(INTMX,4),PSOFA2(INTMX,4),PSOFB1(INTMX,4),PSOFB2(INTMX,4)
     * ,JHKKPZ(INTMX),JHKKTZ(INTMX),JHKKSX(INTMX),JHKKS1(INTMX)
C------------------------
C      WRITE(6,1298)NVV,NSV,NVS,NSS,NDV,NVD,NDS,NSD,NOCC,NONUST,
C    *                                                   NONUJT
C1298 FORMAT(' NVV,NSV,NVS,NSS,NDV,NVD,NDS,NSD,NOCC,NONUST,NONUJT'
C    */11I8)
      IREJ=0
      PXBAL=0.
      PYBAL=0.
      PZBAL=0.
      PEBAL=0.
      PXSS=0.
      PYSS=0.
      PZSS=0.
      PESS=0.
      PXSV=0.
      PYSV=0.
      PZSV=0.
      PESV=0.
      PXVS=0.
      PYVS=0.
      PZVS=0.
      PEVS=0.
      PXVV=0.
      PYVV=0.
      PZVV=0.
      PEVV=0.
      PXDS=0.
      PYDS=0.
      PZDS=0.
      PEDS=0.
      PXSD=0.
      PYSD=0.
      PZSD=0.
      PESD=0.
      PXDZ=0.
      PYDZ=0.
      PZDZ=0.
      PEDZ=0.
      PXZD=0.
      PYZD=0.
      PZZD=0.
      PEZD=0.
      PXDV=0.
      PYDV=0.
      PZDV=0.
      PEDV=0.
      PXVD=0.
      PYVD=0.
      PZVD=0.
      PEVD=0.
      PXCC=0.
      PYCC=0.
      PZCC=0.
      PECC=0.
      PXZZ=0.
      PYZZ=0.
      PZZZ=0.
      PEZZ=0.
      PXHH=0.
      PYHH=0.
      PZHH=0.
      PEHH=0.
C     IF(IP.EQ.1)GAMCM=GAMCM+(NSV+NDV)*AAM(IJPROJ)/UMO
C
C     IF(IHADA.OR.IHADSS) THEN
        DO 10 N=1,NSS
          IF (INLOSS(N))THEN
            IF(ABS(NCHSS1(N)).NE.99) THEN
              PXSS=PXSS + PQSSA1(N,1) + PQSSA2(N,1)
              PYSS=PYSS + PQSSA1(N,2) + PQSSA2(N,2)
              PZSS=PZSS + PQSSA1(N,3) + PQSSA2(N,3)
              PESS=PESS + PQSSA1(N,4) + PQSSA2(N,4)
            ENDIF
            IF(ABS(NCHSS2(N)).NE.99) THEN
              PXSS=PXSS + PQSSB1(N,1) + PQSSB2(N,1)
              PYSS=PYSS + PQSSB1(N,2) + PQSSB2(N,2)
              PZSS=PZSS + PQSSB1(N,3) + PQSSB2(N,3)
              PESS=PESS + PQSSB1(N,4) + PQSSB2(N,4)
            ENDIF
          ENDIF
   10   CONTINUE
        PZBSS=GAMCM*PZSS + BGCM*PESS
        PEBSS=GAMCM*PESS + BGCM*PZSS
        PXBAL=PXSS
        PYBAL=PYSS
        PZBAL=PZBSS
        PEBAL=PEBSS
C     ENDIF
C     IF(IHADA.OR.IHADSV) THEN
        DO 20 N=1,NSV
          IF(ABS(NCHSV1(N)).NE.99) THEN
            PXSV=PXSV +PQSVA1(N,1)+PQSVA2(N,1)
            PYSV=PYSV +PQSVA1(N,2)+PQSVA2(N,2)
            PZSV=PZSV +PQSVA1(N,3)+PQSVA2(N,3)
            PESV=PESV +PQSVA1(N,4)+PQSVA2(N,4)
      IF (IPEV.GE.1)THEN
        WRITE(6,2001)
     +  PXSV,PYSV,PZSV,PESV
      ENDIF
 2001 FORMAT (
     +' SV ',4E15.5)
C
          ENDIF
          IF(ABS(NCHSV2(N)).NE.99) THEN
            PXSV=PXSV + PQSVB1(N,1)+PQSVB2(N,1)
            PYSV=PYSV + PQSVB1(N,2)+PQSVB2(N,2)
            PZSV=PZSV + PQSVB1(N,3)+PQSVB2(N,3)
            PESV=PESV + PQSVB1(N,4)+PQSVB2(N,4)
      IF (IPEV.GE.1)THEN
        WRITE(6,2001)
     +  PXSV,PYSV,PZSV,PESV
      ENDIF
          ENDIF
   20   CONTINUE
        PZBSV=GAMCM*PZSV + BGCM*PESV
        PEBSV=GAMCM*PESV + BGCM*PZSV
      IF (IPEV.GE.1)THEN
        WRITE(6,2001)
     +  PXSV,PYSV,PZBSV,PEBSV
      ENDIF
        PXBAL=PXBAL + PXSV
        PYBAL=PYBAL + PYSV
        PZBAL=PZBAL + PZBSV
        PEBAL=PEBAL + PEBSV
C     ENDIF
C     IF(IHADA.OR.IHADVS) THEN
        DO 30 N=1,NVS
          IF(ABS(NCHVS1(N)).NE.99) THEN
            PXVS=PXVS + PQVSA1(N,1) + PQVSA2(N,1)
            PYVS=PYVS + PQVSA1(N,2) + PQVSA2(N,2)
            PZVS=PZVS + PQVSA1(N,3) + PQVSA2(N,3)
            PEVS=PEVS + PQVSA1(N,4) + PQVSA2(N,4)
          ENDIF
          IF(ABS(NCHVS2(N)).NE.99) THEN
            PXVS=PXVS + PQVSB1(N,1) + PQVSB2(N,1)
            PYVS=PYVS + PQVSB1(N,2) + PQVSB2(N,2)
            PZVS=PZVS + PQVSB1(N,3) + PQVSB2(N,3)
            PEVS=PEVS + PQVSB1(N,4) + PQVSB2(N,4)
          ENDIF
   30   CONTINUE
        PZBVS=GAMCM*PZVS + BGCM*PEVS
        PEBVS=GAMCM*PEVS + BGCM*PZVS
        PXBAL=PXBAL + PXVS
        PYBAL=PYBAL + PYVS
        PZBAL=PZBAL + PZBVS
        PEBAL=PEBAL + PEBVS
        DO 31 N=1,NDS
          IF(ABS(NCHDS1(N)).NE.99) THEN
            PXDS=PXDS + PQDSA1(N,1) + PQDSA2(N,1)
            PYDS=PYDS + PQDSA1(N,2) + PQDSA2(N,2)
            PZDS=PZDS + PQDSA1(N,3) + PQDSA2(N,3)
            PEDS=PEDS + PQDSA1(N,4) + PQDSA2(N,4)
          ENDIF
          IF(ABS(NCHDS2(N)).NE.99) THEN
            PXDS=PXDS + PQDSB1(N,1) + PQDSB2(N,1)
            PYDS=PYDS + PQDSB1(N,2) + PQDSB2(N,2)
            PZDS=PZDS + PQDSB1(N,3) + PQDSB2(N,3)
            PEDS=PEDS + PQDSB1(N,4) + PQDSB2(N,4)
          ENDIF
   31   CONTINUE
        PZBDS=GAMCM*PZDS + BGCM*PEDS
        PEBDS=GAMCM*PEDS + BGCM*PZDS
        PXBAL=PXBAL + PXDS
        PYBAL=PYBAL + PYDS
        PZBAL=PZBAL + PZBDS
        PEBAL=PEBAL + PEBDS
        DO 371 N=1,NDZ
          IF(ABS(NCHDZ1(N)).NE.99) THEN
            PXDZ=PXDZ + PQDZA1(N,1) + PQDZA2(N,1)
            PYDZ=PYDZ + PQDZA1(N,2) + PQDZA2(N,2)
            PZDZ=PZDZ + PQDZA1(N,3) + PQDZA2(N,3)
            PEDZ=PEDZ + PQDZA1(N,4) + PQDZA2(N,4)
          ENDIF
          IF(ABS(NCHDZ2(N)).NE.99) THEN
            PXDZ=PXDZ + PQDZB1(N,1) + PQDZB2(N,1)
            PYDZ=PYDZ + PQDZB1(N,2) + PQDZB2(N,2)
            PZDZ=PZDZ + PQDZB1(N,3) + PQDZB2(N,3)
            PEDZ=PEDZ + PQDZB1(N,4) + PQDZB2(N,4)
          ENDIF
  371   CONTINUE
        PZBDZ=GAMCM*PZDZ + BGCM*PEDZ
        PEBDZ=GAMCM*PEDZ + BGCM*PZDZ
        PXBAL=PXBAL + PXDZ
        PYBAL=PYBAL + PYDZ
        PZBAL=PZBAL + PZBDZ
        PEBAL=PEBAL + PEBDZ
        DO 32 N=1,NSD
          IF(ABS(NCHSD1(N)).NE.99) THEN
            PXSD=PXSD + PQSDA1(N,1) + PQSDA2(N,1)
            PYSD=PYSD + PQSDA1(N,2) + PQSDA2(N,2)
            PZSD=PZSD + PQSDA1(N,3) + PQSDA2(N,3)
            PESD=PESD + PQSDA1(N,4) + PQSDA2(N,4)
          ENDIF
          IF(ABS(NCHSD2(N)).NE.99) THEN
            PXSD=PXSD + PQSDB1(N,1) + PQSDB2(N,1)
            PYSD=PYSD + PQSDB1(N,2) + PQSDB2(N,2)
            PZSD=PZSD + PQSDB1(N,3) + PQSDB2(N,3)
            PESD=PESD + PQSDB1(N,4) + PQSDB2(N,4)
          ENDIF
   32   CONTINUE
        PZBSD=GAMCM*PZSD + BGCM*PESD
        PEBSD=GAMCM*PESD + BGCM*PZSD
        PXBAL=PXBAL + PXSD
        PYBAL=PYBAL + PYSD
        PZBAL=PZBAL + PZBSD
        PEBAL=PEBAL + PEBSD
        DO 372 N=1,NZD
          IF(ABS(NCHZD1(N)).NE.99) THEN
            PXZD=PXZD + PQZDA1(N,1) + PQZDA2(N,1)
            PYZD=PYZD + PQZDA1(N,2) + PQZDA2(N,2)
            PZZD=PZZD + PQZDA1(N,3) + PQZDA2(N,3)
            PEZD=PEZD + PQZDA1(N,4) + PQZDA2(N,4)
          ENDIF
          IF(ABS(NCHZD2(N)).NE.99) THEN
            PXZD=PXZD + PQZDB1(N,1) + PQZDB2(N,1)
            PYZD=PYZD + PQZDB1(N,2) + PQZDB2(N,2)
            PZZD=PZZD + PQZDB1(N,3) + PQZDB2(N,3)
            PEZD=PEZD + PQZDB1(N,4) + PQZDB2(N,4)
          ENDIF
  372    CONTINUE
        PZBZD=GAMCM*PZZD + BGCM*PEZD
        PEBZD=GAMCM*PEZD + BGCM*PZZD
        PXBAL=PXBAL + PXZD
        PYBAL=PYBAL + PYZD
        PZBAL=PZBAL + PZBZD
        PEBAL=PEBAL + PEBZD
        DO 33 N=1,NDV
          IF(ABS(NCHDV1(N)).NE.99) THEN
            PXDV=PXDV + PQDVA1(N,1) + PQDVA2(N,1)
            PYDV=PYDV + PQDVA1(N,2) + PQDVA2(N,2)
            PZDV=PZDV + PQDVA1(N,3) + PQDVA2(N,3)
            PEDV=PEDV + PQDVA1(N,4) + PQDVA2(N,4)
          ENDIF
          IF(ABS(NCHDV2(N)).NE.99) THEN
            PXDV=PXDV + PQDVB1(N,1) + PQDVB2(N,1)
            PYDV=PYDV + PQDVB1(N,2) + PQDVB2(N,2)
            PZDV=PZDV + PQDVB1(N,3) + PQDVB2(N,3)
            PEDV=PEDV + PQDVB1(N,4) + PQDVB2(N,4)
          ENDIF
   33   CONTINUE
        PZBDV=GAMCM*PZDV + BGCM*PEDV
        PEBDV=GAMCM*PEDV + BGCM*PZDV
        PXBAL=PXBAL + PXDV
        PYBAL=PYBAL + PYDV
        PZBAL=PZBAL + PZBDV
        PEBAL=PEBAL + PEBDV
        DO 34 N=1,NVD
          IF(ABS(NCHVD1(N)).NE.99) THEN
            PXVD=PXVD + PQVDA1(N,1) + PQVDA2(N,1)
            PYVD=PYVD + PQVDA1(N,2) + PQVDA2(N,2)
            PZVD=PZVD + PQVDA1(N,3) + PQVDA2(N,3)
            PEVD=PEVD + PQVDA1(N,4) + PQVDA2(N,4)
          ENDIF
          IF(ABS(NCHVD2(N)).NE.99) THEN
            PXVD=PXVD + PQVDB1(N,1) + PQVDB2(N,1)
            PYVD=PYVD + PQVDB1(N,2) + PQVDB2(N,2)
            PZVD=PZVD + PQVDB1(N,3) + PQVDB2(N,3)
            PEVD=PEVD + PQVDB1(N,4) + PQVDB2(N,4)
          ENDIF
   34   CONTINUE
        PZBVD=GAMCM*PZVD + BGCM*PEVD
        PEBVD=GAMCM*PEVD + BGCM*PZVD
        PXBAL=PXBAL + PXVD
        PYBAL=PYBAL + PYVD
        PZBAL=PZBAL + PZBVD
        PEBAL=PEBAL + PEBVD
C     ENDIF
C     IF(IHADA.OR.IHADVV) THEN
        DO 40 N=1,NVV
          IF((NCHVV1(N).NE.99).AND.(NCHVV2(N).NE.99)) THEN
          PXVV=PXVV+PQVVA1(N,1)+PQVVA2(N,1)+PQVVB1(N,1)+PQVVB2(N,1)
          PYVV=PYVV+PQVVA1(N,2)+PQVVA2(N,2)+PQVVB1(N,2)+PQVVB2(N,2)
          PZVV=PZVV+PQVVA1(N,3)+PQVVA2(N,3)+PQVVB1(N,3)+PQVVB2(N,3)
          PEVV=PEVV+PQVVA1(N,4)+PQVVA2(N,4)+PQVVB1(N,4)+PQVVB2(N,4)
          ENDIF
   40   CONTINUE
        PZBVV=GAMCM*PZVV + BGCM*PEVV
        PEBVV=GAMCM*PEVV + BGCM*PZVV
        PXBAL=PXBAL + PXVV
        PYBAL=PYBAL + PYVV
        PZBAL=PZBAL + PZBVV
        PEBAL=PEBAL + PEBVV
C     ENDIF
C     IF(IHADA.OR.IHADSS) THEN
C	WRITE(6,*)' evtest nocc ',NOCC
C       DO 120 N=1,NOCC
C         PXCC=PXCC + POJCC(1,N) + PATCC(1,N)
C         PYCC=PYCC + POJCC(2,N) + PATCC(2,N)
C         PZCC=PZCC + POJCC(3,N) + PATCC(3,N)
C         PECC=PECC + POJCC(4,N) + PATCC(4,N)
C120    CONTINUE
C       PZBCC=GAMCM*PZCC + BGCM*PECC
C       PEBCC=GAMCM*PECC + BGCM*PZCC
C       PXBAL=PXBAL + PXCC
C       PYBAL=PYBAL + PYCC
C       PZBAL=PZBAL + PZBCC
C       PEBAL=PEBAL + PEBCC
C     ENDIF
C	WRITE(6,*)' evtest nonust ',NONUST
        DO 210 N=1,NONUST
            IF(ABS(NCHZZ1(N)).NE.99.AND.JHKKSX(N).EQ.1) THEN
            IF(ABS(NCHZZ1(N)).NE.88) THEN
              PXZZ=PXZZ + PQZZA1(N,1) + PQZZA2(N,1)
              PYZZ=PYZZ + PQZZA1(N,2) + PQZZA2(N,2)
              PZZZ=PZZZ + PQZZA1(N,3) + PQZZA2(N,3)
              PEZZ=PEZZ + PQZZA1(N,4) + PQZZA2(N,4)
            ENDIF
            ENDIF
            IF(ABS(NCHZZ2(N)).NE.99.AND.JHKKSX(N).EQ.1) THEN
            IF(ABS(NCHZZ2(N)).NE.88) THEN
              PXZZ=PXZZ + PQZZB1(N,1) + PQZZB2(N,1)
              PYZZ=PYZZ + PQZZB1(N,2) + PQZZB2(N,2)
              PZZZ=PZZZ + PQZZB1(N,3) + PQZZB2(N,3)
              PEZZ=PEZZ + PQZZB1(N,4) + PQZZB2(N,4)
            ENDIF
            ENDIF
  210   CONTINUE
        PZBZZ=GAMCM*PZZZ + BGCM*PEZZ
        PEBZZ=GAMCM*PEZZ + BGCM*PZZZ
        PXBAL=PXBAL + PXZZ
        PYBAL=PYBAL + PYZZ
        PZBAL=PZBAL + PZBZZ
        PEBAL=PEBAL + PEBZZ
C	WRITE(6,*)' evtest nonujt ',NONUJT
        DO 220 N=1,NONUJT
            IF(ABS(NCHHH1(N)).NE.99.AND.JHKKEX(N).EQ.1) THEN
              PXHH=PXHH + PQHHA1(N,1) + PQHHA2(N,1)
              PYHH=PYHH + PQHHA1(N,2) + PQHHA2(N,2)
              PZHH=PZHH + PQHHA1(N,3) + PQHHA2(N,3)
              PEHH=PEHH + PQHHA1(N,4) + PQHHA2(N,4)
            ENDIF
            IF(ABS(NCHHH2(N)).NE.99.AND.JHKKEX(N).EQ.1) THEN
              PXHH=PXHH + PQHHB1(N,1) + PQHHB2(N,1)
              PYHH=PYHH + PQHHB1(N,2) + PQHHB2(N,2)
              PZHH=PZHH + PQHHB1(N,3) + PQHHB2(N,3)
              PEHH=PEHH + PQHHB1(N,4) + PQHHB2(N,4)
            ENDIF
  220   CONTINUE
        PZBHH=GAMCM*PZHH + BGCM*PEHH
        PEBHH=GAMCM*PEHH + BGCM*PZHH
        PXBAL=PXBAL + PXHH
        PYBAL=PYBAL + PYHH
        PZBAL=PZBAL + PZBHH
        PEBAL=PEBAL + PEBHH
C
      E0000=0.D0
      P0000=0.D0
C	WRITE(6,*)' evtest ip ',IP
      DO 7767 I=1,IP
        IF(ISTHKK(I).EQ.11)E0000=E0000+PRMOM(4,I)
        IF(ISTHKK(I).EQ.11)P0000=P0000+PRMOM(3,I)
 7767 CONTINUE
C	WRITE(6,*)' evtest it ',IT
      DO 7768 II=1,IT
        I=II+IP
        IF(ISTHKK(I).EQ.12)E0000=E0000+TAMOM(4,II)
        IF(ISTHKK(I).EQ.12)P0000=P0000+TAMOM(3,II)
 7768 CONTINUE
      P000=GAMCM*P0000+BGCM*E0000
      E000=GAMCM*E0000+BGCM*P0000
      IPROJO=(PZBAL*1.001)/PPROJ
      RESIDU=ABS(E000-PEBAL)/(E000)
      IF (IPEV.GE.1)THEN
 	WRITE(6,'(A,2E15.5)')' E000,PEBAL', E000,PEBAL
        WRITE(6,1000)PXBAL,PYBAL,PZBAL,PEBAL, PXSS,PYSS,PZBSS,PEBSS,
     +  PXSV,PYSV,PZBSV,PEBSV, PXVS,PYVS,PZBVS,PEBVS, PXVV,PYVV,PZBVV,
     +  PEBVV,PXCC,PYCC,PZBCC,PEBCC,
     +  PXZZ,PYZZ,PZBZZ,PEBZZ,
     +  PXHH,PYHH,PZBHH,PEBHH,
     +  PXDS,PYDS,PZBDS,PEBDS,
     +  PXSD,PYSD,PZBSD,PEBSD,
     +  PXDZ,PYDZ,PZBDZ,PEBDZ,
     +  PXZD,PYZD,PZBZD,PEBZD,
     +  PXDV,PYDV,PZBDV,PEBDV,
     +  PXVD,PYVD,PZBVD,PEBVD
      ENDIF
      IF (RESIDU.GT.0.02D0)THEN
	IREJ=1
      ENDIF
      IF (RESIDU.GT.0.02D0.AND.IPHKK.GE.2)THEN
	IREJ=1
	WRITE(6,'(A,2E15.5)')' E000,PEBAL', E000,PEBAL
        WRITE(6,1000)PXBAL,PYBAL,PZBAL,PEBAL, PXSS,PYSS,PZBSS,PEBSS,
     +  PXSV,PYSV,PZBSV,PEBSV, PXVS,PYVS,PZBVS,PEBVS, PXVV,PYVV,PZBVV,
     +  PEBVV,PXCC,PYCC,PZBCC,PEBCC,
     +  PXZZ,PYZZ,PZBZZ,PEBZZ,
     +  PXHH,PYHH,PZBHH,PEBHH,
     +  PXDS,PYDS,PZBDS,PEBDS,
     +  PXSD,PYSD,PZBSD,PEBSD,
     +  PXDZ,PYDZ,PZBDZ,PEBDZ,
     +  PXZD,PYZD,PZBZD,PEBZD,
     +  PXDV,PYDV,PZBDV,PEBDV,
     +  PXVD,PYVD,PZBVD,PEBVD
      ENDIF
 1000 FORMAT (' 4 MOMENTUM CONS.IN EVENT LEVEL OF PARTONS',/ ' ALL',4E15
     +.5/,' SS ',4E15.5/,' SV ',4E15.5/ ' VS ',4E15.5/,' VV ',4E15.5/,
     + ' CC ',4E15.5/
     + ' ZZ ',4E15.5/
     + ' HH ',4E15.5/
     + ' DS ',4E15.5/
     + ' SD ',4E15.5/
     + ' DZ ',4E15.5/
     + ' ZD ',4E15.5/
     + ' DV ',4E15.5/
     + ' VD ',4E15.5)
C
      PXBAL=0.
      PYBAL=0.
      PZBAL=0.
      PEBAL=0.
      PXSS=0.
      PYSS=0.
      PZSS=0.
      PESS=0.
      PXSV=0.
      PYSV=0.
      PZSV=0.
      PESV=0.
      PXVS=0.
      PYVS=0.
      PZVS=0.
      PEVS=0.
      PXVV=0.
      PYVV=0.
      PZVV=0.
      PEVV=0.
      PXCC=0.
      PYCC=0.
      PZCC=0.
      PECC=0.
      PXDS=0.
      PYDS=0.
      PZDS=0.
      PEDS=0.
      PXSD=0.
      PYSD=0.
      PZSD=0.
      PESD=0.
      PXDV=0.
      PYDV=0.
      PZDV=0.
      PEDV=0.
      PXVD=0.
      PYVD=0.
      PZVD=0.
      PEVD=0.
      PXZZ=0.
      PYZZ=0.
      PZZZ=0.
      PEZZ=0.
      PXHH=0.
      PYHH=0.
      PZHH=0.
      PEHH=0.
C
C     IF(IHADA.OR.IHADSS) THEN
        DO 50 N=1,NSS
          IF (INLOSS(N))THEN
          IF(ABS(NCHSS1(N)).NE.99) THEN
            PXSS=PXSS+BGXSS1(N)*AMCSS1(N)
            PYSS=PYSS+BGYSS1(N)*AMCSS1(N)
            PZSS=PZSS+BGZSS1(N)*AMCSS1(N)
            PESS=PESS+GACSS1(N)*AMCSS1(N)
          ENDIF
          IF(ABS(NCHSS2(N)).NE.99) THEN
            PXSS=PXSS+BGXSS2(N)*AMCSS2(N)
            PYSS=PYSS+BGYSS2(N)*AMCSS2(N)
            PZSS=PZSS+BGZSS2(N)*AMCSS2(N)
            PESS=PESS+GACSS2(N)*AMCSS2(N)
          ENDIF
          ENDIF
   50   CONTINUE
        PZBSS=GAMCM*PZSS + BGCM*PESS
        PEBSS=GAMCM*PESS + BGCM*PZSS
C       DO 130 N=1,NOCC
C         PXCC=PXCC + BGXCC(N)*AMCC(N)
C         PYCC=PYCC + BGYCC(N)*AMCC(N)
C         PZCC=PZCC + BGZCC(N)*AMCC(N)
C         PECC=PECC + GACC(N)*AMCC(N)
C130    CONTINUE
        PXBAL=PXSS
        PYBAL=PYSS
        PZBAL=PZBSS
        PEBAL=PEBSS
C       PZBCC=GAMCM*PZCC + BGCM*PECC
C       PEBCC=GAMCM*PECC + BGCM*PZCC
C       PXBAL=PXBAL + PXCC
C       PYBAL=PYBAL + PYCC
C       PZBAL=PZBAL + PZBCC
C       PEBAL=PEBAL + PEBCC
C     ENDIF
C     IF(IHADA.OR.IHADSV) THEN
        DO 60 N=1,NSV
          IF(ABS(NCHSV1(N)).NE.99) THEN
          PXSV=PXSV+BGXSV1(N)*AMCSV1(N)
          PYSV=PYSV+BGYSV1(N)*AMCSV1(N)
          PZSV=PZSV+BGZSV1(N)*AMCSV1(N)
          PESV=PESV+GACSV1(N)*AMCSV1(N)
          ENDIF
          IF(ABS(NCHSV2(N)).NE.99) THEN
          PXSV=PXSV+BGXSV2(N)*AMCSV2(N)
          PYSV=PYSV+BGYSV2(N)*AMCSV2(N)
          PZSV=PZSV+BGZSV2(N)*AMCSV2(N)
          PESV=PESV+GACSV2(N)*AMCSV2(N)
          ENDIF
   60   CONTINUE
        PZBSV=GAMCM*PZSV + BGCM*PESV
        PEBSV=GAMCM*PESV + BGCM*PZSV
        PXBAL=PXBAL + PXSV
        PYBAL=PYBAL + PYSV
        PZBAL=PZBAL + PZBSV
        PEBAL=PEBAL + PEBSV
        DO 61 N=1,NDS
          IF(ABS(NCHDS1(N)).NE.99) THEN
          PXDS=PXDS+BGXDS1(N)*AMCDS1(N)
          PYDS=PYDS+BGYDS1(N)*AMCDS1(N)
          PZDS=PZDS+BGZDS1(N)*AMCDS1(N)
          PEDS=PEDS+GACDS1(N)*AMCDS1(N)
          ENDIF
          IF(ABS(NCHDS2(N)).NE.99) THEN
          PXDS=PXDS+BGXDS2(N)*AMCDS2(N)
          PYDS=PYDS+BGYDS2(N)*AMCDS2(N)
          PZDS=PZDS+BGZDS2(N)*AMCDS2(N)
          PEDS=PEDS+GACDS2(N)*AMCDS2(N)
          ENDIF
   61   CONTINUE
        PZBDS=GAMCM*PZDS + BGCM*PEDS
        PEBDS=GAMCM*PEDS + BGCM*PZDS
        PXBAL=PXBAL + PXDS
        PYBAL=PYBAL + PYDS
        PZBAL=PZBAL + PZBDS
        PEBAL=PEBAL + PEBDS
        DO 62 N=1,NSD
          IF(ABS(NCHSD1(N)).NE.99) THEN
          PXSD=PXSD+BGXSD1(N)*AMCSD1(N)
          PYSD=PYSD+BGYSD1(N)*AMCSD1(N)
          PZSD=PZSD+BGZSD1(N)*AMCSD1(N)
          PESD=PESD+GACSD1(N)*AMCSD1(N)
          ENDIF
          IF(ABS(NCHSD2(N)).NE.99) THEN
          PXSD=PXSD+BGXSD2(N)*AMCSD2(N)
          PYSD=PYSD+BGYSD2(N)*AMCSD2(N)
          PZSD=PZSD+BGZSD2(N)*AMCSD2(N)
          PESD=PESD+GACSD2(N)*AMCSD2(N)
          ENDIF
   62   CONTINUE
        PZBSD=GAMCM*PZSD + BGCM*PESD
        PEBSD=GAMCM*PESD + BGCM*PZSD
        PXBAL=PXBAL + PXSD
        PYBAL=PYBAL + PYSD
        PZBAL=PZBAL + PZBSD
        PEBAL=PEBAL + PEBSD
        DO 63 N=1,NDV
          IF(ABS(NCHDV1(N)).NE.99) THEN
          PXDV=PXDV+BGXDV1(N)*AMCDV1(N)
          PYDV=PYDV+BGYDV1(N)*AMCDV1(N)
          PZDV=PZDV+BGZDV1(N)*AMCDV1(N)
          PEDV=PEDV+GACDV1(N)*AMCDV1(N)
          ENDIF
          IF(ABS(NCHDV2(N)).NE.99) THEN
          PXDV=PXDV+BGXDV2(N)*AMCDV2(N)
          PYDV=PYDV+BGYDV2(N)*AMCDV2(N)
          PZDV=PZDV+BGZDV2(N)*AMCDV2(N)
          PEDV=PEDV+GACDV2(N)*AMCDV2(N)
          ENDIF
   63   CONTINUE
        PZBDV=GAMCM*PZDV + BGCM*PEDV
        PEBDV=GAMCM*PEDV + BGCM*PZDV
        PXBAL=PXBAL + PXDV
        PYBAL=PYBAL + PYDV
        PZBAL=PZBAL + PZBDV
        PEBAL=PEBAL + PEBDV
        DO 64 N=1,NVD
          IF(ABS(NCHVD1(N)).NE.99) THEN
          PXVD=PXVD+BGXVD1(N)*AMCVD1(N)
          PYVD=PYVD+BGYVD1(N)*AMCVD1(N)
          PZVD=PZVD+BGZVD1(N)*AMCVD1(N)
          PEVD=PEVD+GACVD1(N)*AMCVD1(N)
          ENDIF
          IF(ABS(NCHVD2(N)).NE.99) THEN
          PXVD=PXVD+BGXVD2(N)*AMCVD2(N)
          PYVD=PYVD+BGYVD2(N)*AMCVD2(N)
          PZVD=PZVD+BGZVD2(N)*AMCVD2(N)
          PEVD=PEVD+GACVD2(N)*AMCVD2(N)
          ENDIF
   64   CONTINUE
        PZBVD=GAMCM*PZVD + BGCM*PEVD
        PEBVD=GAMCM*PEVD + BGCM*PZVD
        PXBAL=PXBAL + PXVD
        PYBAL=PYBAL + PYVD
        PZBAL=PZBAL + PZBVD
        PEBAL=PEBAL + PEBVD
C     ENDIF
C     IF(IHADA.OR.IHADVS) THEN
        DO 70 N=1,NVS
          IF(ABS(NCHVS1(N)).NE.99) THEN
          PXVS=PXVS+BGXVS1(N)*AMCVS1(N)
          PYVS=PYVS+BGYVS1(N)*AMCVS1(N)
          PZVS=PZVS+BGZVS1(N)*AMCVS1(N)
          PEVS=PEVS+GACVS1(N)*AMCVS1(N)
          ENDIF
          IF(ABS(NCHVS2(N)).NE.99) THEN
          PXVS=PXVS+BGXVS2(N)*AMCVS2(N)
          PYVS=PYVS+BGYVS2(N)*AMCVS2(N)
          PZVS=PZVS+BGZVS2(N)*AMCVS2(N)
          PEVS=PEVS+GACVS2(N)*AMCVS2(N)
          ENDIF
   70   CONTINUE
        PZBVS=GAMCM*PZVS + BGCM*PEVS
        PEBVS=GAMCM*PEVS + BGCM*PZVS
        PXBAL=PXBAL + PXVS
        PYBAL=PYBAL + PYVS
        PZBAL=PZBAL + PZBVS
        PEBAL=PEBAL + PEBVS
C     ENDIF
        DO 250 N=1,NONUST
          IF(ABS(NCHZZ1(N)).NE.99.AND.JHKKSX(N).EQ.1) THEN
            PXZZ=PXZZ+BGXZZ1(N)*AMCZZ1(N)
            PYZZ=PYZZ+BGYZZ1(N)*AMCZZ1(N)
            PZZZ=PZZZ+BGZZZ1(N)*AMCZZ1(N)
            PEZZ=PEZZ+GACZZ1(N)*AMCZZ1(N)
          ENDIF
          IF(ABS(NCHZZ2(N)).NE.99.AND.JHKKSX(N).EQ.1) THEN
            PXZZ=PXZZ+BGXZZ2(N)*AMCZZ2(N)
            PYZZ=PYZZ+BGYZZ2(N)*AMCZZ2(N)
            PZZZ=PZZZ+BGZZZ2(N)*AMCZZ2(N)
            PEZZ=PEZZ+GACZZ2(N)*AMCZZ2(N)
          ENDIF
  250   CONTINUE
        PZBZZ=GAMCM*PZZZ + BGCM*PEZZ
        PEBZZ=GAMCM*PEZZ + BGCM*PZZZ
        PXBAL=PXBAL + PXZZ
        PYBAL=PYBAL + PYZZ
        PZBAL=PZBAL + PZBZZ
        PEBAL=PEBAL + PEBZZ
        DO 260 N=1,NONUJT
          IF(ABS(NCHHH1(N)).NE.99.AND.JHKKEX(N).EQ.1) THEN
            PXHH=PXHH+BGXHH1(N)*AMCHH1(N)
            PYHH=PYHH+BGYHH1(N)*AMCHH1(N)
            PZHH=PZHH+BGZHH1(N)*AMCHH1(N)
            PEHH=PEHH+GACHH1(N)*AMCHH1(N)
          ENDIF
          IF(ABS(NCHHH2(N)).NE.99.AND.JHKKEX(N).EQ.1) THEN
            PXHH=PXHH+BGXHH2(N)*AMCHH2(N)
            PYHH=PYHH+BGYHH2(N)*AMCHH2(N)
            PZHH=PZHH+BGZHH2(N)*AMCHH2(N)
            PEHH=PEHH+GACHH2(N)*AMCHH2(N)
          ENDIF
  260   CONTINUE
        PZBHH=GAMCM*PZHH + BGCM*PEHH
        PEBHH=GAMCM*PEHH + BGCM*PZHH
        PXBAL=PXBAL + PXHH
        PYBAL=PYBAL + PYHH
        PZBAL=PZBAL + PZBHH
        PEBAL=PEBAL + PEBHH
C     IF(IHADA.OR.IHADVV) THEN
        DO 80 N=1,NVV
          IF((NCHVV1(N).NE.99).AND.(NCHVV2(N).NE.99)) THEN
          PXVV=PXVV+BGXVV1(N)*AMCVV1(N)+BGXVV2(N)*AMCVV2(N)
          PYVV=PYVV+BGYVV1(N)*AMCVV1(N)+BGYVV2(N)*AMCVV2(N)
          PZVV=PZVV+BGZVV1(N)*AMCVV1(N)+BGZVV2(N)*AMCVV2(N)
          PEVV=PEVV+GACVV1(N)*AMCVV1(N)+GACVV2(N)*AMCVV2(N)
          ENDIF
   80   CONTINUE
        PZBVV=GAMCM*PZVV + BGCM*PEVV
        PEBVV=GAMCM*PEVV + BGCM*PZVV
        PXBAL=PXBAL + PXVV
        PYBAL=PYBAL + PYVV
        PZBAL=PZBAL + PZBVV
        PEBAL=PEBAL + PEBVV
C     ENDIF
C
      IF (IPEV.GE.1) WRITE(6,1010)PXBAL,PYBAL,PZBAL,
     +PEBAL, PXSS,PYSS,PZBSS,PEBSS, PXSV,PYSV,PZBSV,PEBSV, PXVS,PYVS,
     +PZBVS,PEBVS, PXVV,PYVV,PZBVV,PEBVV, PXCC,PYCC,PZBCC,PEBCC,
     +  PXDS,PYDS,PZBDS,PEBDS,
     +  PXZZ,PYZZ,PZBZZ,PEBZZ,
     +  PXHH,PYHH,PZBHH,PEBHH,
     +  PXSD,PYSD,PZBSD,PEBSD,
     +  PXDV,PYDV,PZBDV,PEBDV,
     +  PXVD,PYVD,PZBVD,PEBVD
 1010 FORMAT (' 4 MOMENTUM CONS.IN EVENT LEVEL OF CHAINS',/ ' ALL',4E15.
     +5/,' SS ',4E15.5/,' SV ',4E15.5/ ' VS ',4E15.5/,' VV ',4E15.5/,
     +  ' CC ',4E15.5/
     + ' DS ',4E15.5/
     + ' ZZ ',4E15.5/
     + ' HH ',4E15.5/
     + ' SD ',4E15.5/
     + ' DV ',4E15.5/
     + ' VD ',4E15.5)
C
      RETURN
      END
*-- Author :
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE CORVAL(AMMM,IREJ,AMCH1,AMCH2, QTX1,QTY1,QZ1,QE1,QTX2,
     +QTY2,QZ2,QE2,NORIG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C KINEMATICAL CORRECTION OF TWO-VALENCE CHAIN SYSTEM
C ACCORDING TO 2-PARTICLE KINEMATICS WITH FIXED MASSES
C
C**** WIR BRAUCHEN AUCH NOCH DIE NEUEN 4-IMPULSE DER KETTENENDEN
C
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEND.
C-----------------------------------------
      IREJ=0
      IF(AMMM.LE.AMCH1+AMCH2+0.4D0) THEN
        IREJ=1
        RETURN
      ENDIF
C
      EK1=(AMMM**2-AMCH2**2 + AMCH1**2)/(2.*AMMM)
      EK2=AMMM - EK1
      PZK1=SQRT(EK1**2 - AMCH1**2)
      PZK1=SIGN(PZK1,QZ1)
      PZK2=SQRT(EK2**2 - AMCH2**2)
      PZK2=SIGN(PZK2,QZ2)
      PXK1=0.
      PYK1=0.
      PXK2=0.
      PYK2=0.
      QTX2=PXK2
      QTY2=PYK2
      QZ2=PZK2
      QE2=EK2
      QTX1=PXK1
      QTY1=PYK1
      QZ1=PZK1
      QE1=EK1
C                      ROTATE NEW CHAIN MOMENTA
C                      INTO DIRECTION OF CHAINS BEFORE CORRECTION
C     GAM=(QE1+QE2)/AMMM
C     BGX=(QTX1+QTX2)/AMMM
C     BGY=(QTY1+QTY2)/AMMM
C     BGZ=(QZ1+QZ2)/AMMM
C
C     IF(ABS(GAM-1.D0).GT.1D-4) THEN
C       WRITE(6,'(A,I10,A/6(1PE15.5)/15X,5(1PE15.4))')
C    +  ' CORVAL: INCONSISTENT KINEMATICS OF CHAINS NORIG= ',NORIG,
C    +  ' AMMM,AMCH1,QE1,QTX1,QTY1, QZ1,AMCH2,QE2,QTX2,QTY2,QZ2',
C    +    AMMM,
C    +    AMCH1, QE1,
C    +  QTX1, QTY1, QZ1, AMCH2,QE2, QTX2, QTY2, QZ2
C       IREJ=1
C     ENDIF
C
C     CALL DALTRA(GAM,-BGX,-BGY,-BGZ,PXK1,PYK1,PZK1,EK1,PPPCH1, QTX1,
C    +QTY1,QZ1,QE1)
C     CALL DALTRA(GAM,-BGX,-BGY,-BGZ,PXK2,PYK2,PZK2,EK2,PPPCH2, QTX2,
C    +QTY2,QZ2,QE2)
C     IF(IPRI.GT.1) THEN
CC       WRITE(6,'(2A)') ' CORVAL - CORRECTION OF CHAIN MOMENTA',
C   +  ' IF MASS OF CHAIN 2 HAD TO BE CHANGED'
C     ENDIF
      RETURN
      END
 
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HADRHH
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C-------------------------
C
C                       HADRONIZE HARD CHAINS
C
C                       ADD GENERATED HADRONS TO /ALLPAR/
C                          STARTING AT (NAUX + 1)
C                       AND TO /HKKEVT/ STARTING AT (NHKK + 1)
C
C---------------------------------------------------------
      PARAMETER (INTMX=2488,INTMD=252)
      COMMON /ABRJT/XJQ1(INTMX),XJAQ1(INTMX),XJQ2(INTMX),XJAQ2(INTMX),
     *        IJJQ1(INTMX),IJJAQ1(INTMX),IJJQ2(INTMX),IJJAQ2(INTMX),
     *        AMJCH1(INTMX),AMJCH2(INTMX),GAMJH1(INTMX),GAMJH2(INTMX),
     *        BGJH1(INTMX),BGJH2(INTMX),THEJH1(INTMX),THEJH2(INTMX),
     *        BGXJH1(INTMX),BGYJH1(INTMX),BGZJH1(INTMX),
     *        BGXJH2(INTMX),BGYJH2(INTMX),BGZJH2(INTMX),
     *  PJETA1(INTMX,4),PJETA2(INTMX,4),PJETB1(INTMX,4),PJETB2(INTMX,4)
     * ,JHKKPH(INTMX),JHKKTH(INTMX),JHKKEX(INTMX),JHKKE1(INTMX)
*KEEP,INTNEW.
      COMMON /INTNEW/ NVV,NSV,NVS,NSS,NDV,NVD,NDS,NSD,
     +IXPV,IXPS,IXTV,IXTS, INTVV1(248),
     +INTVV2(248),INTSV1(248),INTSV2(248), INTVS1(248),INTVS2(248),
     +INTSS1(INTMX),INTSS2(INTMX),
     +INTDV1(248),INTDV2(248),INTVD1(248),INTVD2(248),
     +INTDS1(INTMD),INTDS2(INTMD),INTSD1(INTMD),INTSD2(INTMD)
 
C  /INTNEW/
C       NVV    :   NUMBER OF INTERACTING VALENCE-VALENCE SYSTEMS
C       NSV    :   NUMBER OF INTERACTING SEA-VALENCE SYSTEMS
C       NVS    :   NUMBER OF INTERACTING VALENCE-SEA SYSTEMS
C       NSS    :   NUMBER OF INTERACTING SEA-SEA SYSTEMS
C       IXPV, IXTV : NUMBER OF GENERATED X-VALUES FOR VALENCE-QUARK
C                     SYSTEMS FROM PROJECTILE/TARGET NUCLEI
C       IXPS, IXTS : NUMBER OF GENERATED X-VALUES FOR SEA-QUARK PAIRS
C                     FROM PROJECTILE/TARGET NUCLEI
C-------------------
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX),
     *                IFROVT(248),ITOVT(248),IFROST(INTMX),
     *             JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),JHKKNT(248),
     *                JHKKPV(INTMX),JHKKPS(INTMX),
     *                JHKKTV(INTMX),JHKKTS(INTMX),
     *                MHKKVV(INTMX),MHKKSS(INTMX),
     &                MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
C-------------------
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
C.....................................................................
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/  ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),
     *ZUOST(INTMX),
     *                INTLO(INTMX),INLOSS(INTMX)
C  /LOZUO/
C       INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C
      COMMON /ABRHH/ AMCHH1(INTMX),AMCHH2(INTMX),
     *               GACHH1(INTMX),GACHH2(INTMX),
     *               BGXHH1(INTMX),BGYHH1(INTMX),BGZHH1(INTMX),
     *               BGXHH2(INTMX),BGYHH2(INTMX),BGZHH2(INTMX),
     *               NCHHH1(INTMX),NCHHH2(INTMX),
     *               IJCHH1(INTMX),IJCHH2(INTMX),
     *               PQHHA1(INTMX,4),PQHHA2(INTMX,4),
     *               PQHHB1(INTMX,4),PQHHB2(INTMX,4)
      COMMON /HARDHA/NHARD1,NHKKHA
C
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
       COMMON /NCOUCH/ ACOUVV,ACOUSS,ACOUSV,ACOUVS,
     *                 ACOUZZ,ACOUHH,ACOUDS,ACOUSD,
     *                 ACOUDZ,ACOUZD,ACOUDI,
     *                 ACOUDV,ACOUVD,ACOUCC
C---------------------
      COMMON /PSHOW/ IPSHOW
C     COMMON /HARLUN/ IHARLU,QLUN
      COMMON /HARLUN/ QLUN,IHARLU
      COMMON /JSPAR/PXJ(1000),PYJ(1000),PZJ(1000),HEJ(1000),NNNPJ
      COMMON /JSPA/PXS(40000),PYS(40000),PZS(40000),HES(40000),NNNPS
C--------------------
*KEEP,DFINPA.
      CHARACTER*8 ANF
      PARAMETER (NFIMAX=249)
      COMMON /DFINPA/ ANF(NFIMAX),PXF(NFIMAX),PYF(NFIMAX),PZF(NFIMAX),
     +HEP(NFIMAX),AMF(NFIMAX), ICHF(NFIMAX),IBARF(NFIMAX),NREF(NFIMAX)
      COMMON /DFINPZ/IORMO(NFIMAX),IDAUG1(NFIMAX),IDAUG2(NFIMAX),
     * ISTATH(NFIMAX)
C-------------------
      PARAMETER (NMXHKK=  89998)
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &              JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),PHKK(5,NMXHKK),
     & VHKK(4,NMXHKK),WHKK(4,NMXHKK)
C
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
      COMMON /PROJK/ IPROJK
      COMMON /NUCJTN/NONUJ1,NONUJT,NONUS1,NONUST
      COMMON /GLUSPL/NUGLUU,NSGLUU
      COMMON /NOMIJE/ PTMIJE(10),NNMIJE(10)
C
      DIMENSION POJ(4),PAT(4)
      DATA NCALHH /0/
C-----------------------------------------------------------------------
        NHARD1=NHKK+1
        DO 20 I=1,NONUJT
          NCALHH=NCALHH+1
C
          IF (IPHKK.GE.2)WRITE(6,7789)NONUJT,NCALHH
 7789     FORMAT (' HADRHH NONUJT,NCALHH ',2I10)
          IF (JHKKEX(I).EQ.1)THEN
          IF (I.GT.INTMX)THEN
            WRITE (6,7744)I,INTMX
 7744       FORMAT ('  HADRHH I.GT.INTMX  ',2I10)
            RETURN
          ENDIF
C
C++++++++++++++++++++++++++++++    CHAIN 1:  QUARK-ANTIQUARK   +++++++
            IFB1=IJJQ1(I)
            IFB2=IJJAQ1(I)
            IFB2=IABS(IFB2)+6
            DO 21 J=1,4
              POJ(J)=PJETA1(I,J)
              PAT(J)=PJETA2(I,J)
  21        CONTINUE
        PT1=SQRT(POJ(1)**2+POJ(2)**2)
        PT2=SQRT(PAT(1)**2+PAT(2)**2)
        CALL PARPT(2,PT1,PT2,6,NEVT)
            IHARLU=0
            QLUN=0.
            IF(IPSHOW.EQ.1)THEN
              POJPT=SQRT(POJ(2)**2+POJ(1)**2)
              PATPT=SQRT(PAT(1)**2+PAT(2)**2)
	      DO IIII=1,10
              IF(POJPT.GE.PTMIJE(IIII))NNMIJE(IIII)=
     *	      NNMIJE(IIII)+1
              IF(PATPT.GE.PTMIJE(IIII))NNMIJE(IIII)=
     *	      NNMIJE(IIII)+1
	      ENDDO
              QLUN=MIN(POJPT,PATPT)
              IF((QLUN.LT.2.5D0).OR.(AMJCH1(I).LT.5.D0))THEN
                QLUN=0.
                IHARLU=0
              ELSE 
                IHARLU=1
              ENDIF   
            ENDIF
C----------------------------------------------------------------
            IF (GAMJH1(I).LT.0.001D0.OR.AMJCH1(I).LT.2.D0)THEN
            WRITE (6,7788)
     *      I,NHAD,AMJCH1(I),POJ,PAT,GAMJH1(I),BGXJH1(I),
     *      BGYJH1(I),BGZJH1(I),IFB1,IFB2,IFB3,IFB4,JHKKEX(I)
 7788       FORMAT (' HADRHH ',2I5,8E12.2/5E12.2,5I5)
            GO TO 9977
            ENDIF
            CALL HADJET(NHAD,AMJCH1(I),POJ,PAT,GAMJH1(I),BGXJH1(I),
     *              BGYJH1(I),BGZJH1(I),IFB1,IFB2,IFB3,IFB4,
     *              13,13,3,0,13)
            ACOUHH=ACOUHH+1
            IHARLU=0
            QLUN=0.
            NHKKAU=NHKK+1
	    IF(IPHKK.GE.3)WRITE(6,*)' HADRHH:NHKK,NHKKAU ',NHKK,NHKKAU
            IF (NHAD.GT.NFIMAX) THEN
              WRITE (6,7755)NHAD,NFIMAX
 7755         FORMAT (' NHAD.GT.NFIMAX ',2I10)
              RETURN
            ENDIF
            DO 22 J=1,NHAD
C             NHKK=NHKK+1
              IF (NHKK.EQ.NMXHKK) THEN
                WRITE (*,'(A,2I5/A)') ' HADRHH: NHKK.EQ.NMXHKK ',
     *          NHKK,NMXHKK
                RETURN
              ENDIF
C
              EHECC=SQRT(PXF(J)**2+PYF(J)**2+PZF(J)**2+AMF(J)**2)
              IF (ABS(EHECC-HEP(J)).GT.0.001D0) THEN
C               WRITE(*,'(2A/3I5,3E15.6)')
C    &            ' HADRSV / CHAIN 1 : CORRECT INCONSISTENT ENERGY ',
C    *            '  NCALHH, NHKK,NREF(J), HEP(J),EHECC, AMF(J)',
C    *            NCALHH, NHKK,NREF(J), HEP(J),EHECC, AMF(J)
                HEP(J)=EHECC
              ENDIF
              ANNHH=ANNHH+1.
              EEHH=EEHH+HEP(J)
              PTHH=SQRT(PXF(J)**2+PYF(J)**2)+PTHH
C                               PUT NN-CMS HADRONS INTO /HKKEVT/
          ISTIST=1
          IF(IBARF(J).EQ.500)ISTIST=2
          CALL HKKFIL(ISTIST,MPDGHA(NREF(J)),MHKKHH(I)-3,0,
     *                PXF(J),PYF(J),PZF(J),HEP(J),NHKKAU,IORMO(J),9)
              IF(IDHKK(NHKK).EQ.99999) WRITE (6,5009)NHKK,NREF(J),
     *                                               IDHKK(NHKK)
 	 IF(IPHKK.GE.3) WRITE(6,*)' First chain HADRHH'
              IF (IPHKK.GE.3) WRITE(6,5001) NHKK,
     *        ISTHKK(NHKK),IDHKK(NHKK),JMOHKK(1,NHKK),JMOHKK(2,NHKK),
     &        JDAHKK(1,NHKK),JDAHKK(2,NHKK),(PHKK(KHKK,NHKK),KHKK=1,5),
     &        (VHKK(KHKK,NHKK),KHKK=1,4)
   22       CONTINUE
C           JDAHKK(1,IMOHKK)=NHKKAU
C           JDAHKK(2,IMOHKK)=NHKK
          IF(NNNPJ.GE.1)THEN
            NNNPSO=NNNPS
            NNNPS=NNNPS+1
            NNNPSU=NNNPSO+NNNPJ
            DO 137 J=NNNPS,NNNPSU
              JJ=J-NNNPS+1
              IF(J.GT.40000.OR.JJ.GT.1000)THEN
C               WRITE(6,'(A,2I10)')' J.gt.40000.or.jj.gt.1000 ',J,JJ
                GO TO 137
              ENDIF
	      PXS(J)=PXJ(JJ)
	      PYS(J)=PYJ(JJ)
	      PZS(J)=PZJ(JJ)
	      HES(J)=HEJ(JJ)
  137      CONTINUE
            NNNPS=NNNPS+NNNPJ-1
          ENDIF
 9977       continue
C+++++++++++++++++++++++++++++   CHAIN 2:  AQUARK-QUARK  ++++++++++++++
            IF (NUGLUU.EQ.1) GO TO 5111
            IFB1=IJJAQ2(I)
            IFB2=IJJQ2(I)
            IFB1=IABS(IFB1)+6
            DO 23 J=1,4
              POJ(J)=PJETB1(I,J)
              PAT(J)=PJETB2(I,J)
  23        CONTINUE
        PT1=SQRT(POJ(1)**2+POJ(2)**2)
        PT2=SQRT(PAT(1)**2+PAT(2)**2)
        CALL PARPT(2,PT1,PT2,6,NEVT)
            IHARLU=0
            QLUN=0.
            IF(IPSHOW.EQ.1)THEN
              POJPT=SQRT(POJ(2)**2+POJ(1)**2)
              PATPT=SQRT(PAT(1)**2+PAT(2)**2)
	      DO IIII=1,10
              IF(POJPT.GE.PTMIJE(IIII))NNMIJE(IIII)=
     *	      NNMIJE(IIII)+1
              IF(PATPT.GE.PTMIJE(IIII))NNMIJE(IIII)=
     *	      NNMIJE(IIII)+1
	      ENDDO
              QLUN=MIN(POJPT,PATPT)
              IF((QLUN.LT.2.5D0).OR.(AMJCH2(I).LT.5.D0))THEN
                QLUN=0.
                IHARLU=0
              ELSE 
                IHARLU=1
              ENDIF   
            ENDIF
C
            CALL HADJET(NHAD,AMJCH2(I),POJ,PAT,GAMJH2(I),BGXJH2(I),
     *                BGYJH2(I),BGZJH2(I),IFB1,IFB2,IFB3,IFB4,
     *                13,13,3,0,14)
            IHARLU=0
            QLUN=0.
C                                   ADD HADRONS/RESONANCES INTO
C                                   COMMON /ALLPAR/ STARTING AT NAUX
            NHKKAU=NHKK+1
          DO 24 J=1,NHAD
C           NHKK=NHKK+1
            IF (NHKK.EQ.NMXHKK) THEN
              WRITE (*,'(A,2I5/A)') ' HADRHH: NHKK.EQ.NMXHKK ',
     &                               NHKK,NMXHKK
              RETURN
            ENDIF
C
            EHECC=SQRT(PXF(J)**2+PYF(J)**2+PZF(J)**2+AMF(J)**2)
            IF (ABS(EHECC-HEP(J)).GT.0.001D0) THEN
C             WRITE(*,'(2A/3I5,3E15.6)')
C    &            ' HADRHH / CHAIN 2 : CORRECT INCONSISTENT ENERGY ',
C    *            '  NCALHH, NHKK,NREF(J), HEP(J),EHECC, AMF(J)',
C    *            NCALHH, NHKK,NREF(J), HEP(J),EHECC, AMF(J)
              HEP(J)=EHECC
            ENDIF
            ANNHH=ANNHH+1.
            EEHH=EEHH+HEP(J)
            PTHH=SQRT(PXF(J)**2+PYF(J)**2)+PTHH
C                               PUT NN-CMS HADRONS INTO /HKKEVT/
          ISTIST=1
          IF(IBARF(J).EQ.500)ISTIST=2
          CALL HKKFIL(ISTIST,MPDGHA(NREF(J)),MHKKHH(I),0,
     *                PXF(J),PYF(J),PZF(J),HEP(J),NHKKAU,IORMO(J),10)
            IF(IDHKK(NHKK).EQ.99999) WRITE (6,5009)NHKK,NREF(J),
     *                                             IDHKK(NHKK)
C	  WRITE(6,*)' Second chain HADRHH'
            IF (IPHKK.GE.7) WRITE(6,5001) NHKK,
     *      ISTHKK(NHKK),IDHKK(NHKK),JMOHKK(1,NHKK),JMOHKK(2,NHKK),
     &      JDAHKK(1,NHKK),JDAHKK(2,NHKK),(PHKK(KHKK,NHKK),KHKK=1,5),
     &      (VHKK(KHKK,NHKK),KHKK=1,4)
   24     CONTINUE
C         JDAHKK(1,IMOHKK)=NHKKAU
C         JDAHKK(2,IMOHKK)=NHKK
          IF(NNNPJ.GE.1)THEN
            NNNPSO=NNNPS
            NNNPS=NNNPS+1
            NNNPSU=NNNPSO+NNNPJ
            DO 187 J=NNNPS,NNNPSU
              JJ=J-NNNPS+1
              IF(J.GT.40000.OR.JJ.GT.1000)THEN
C               WRITE(6,'(A,2I10)')' J.gt.40000.or.jj.gt.1000 ',J,JJ
                GO TO 187
              ENDIF
	      PXS(J)=PXJ(JJ)
	      PYS(J)=PYJ(JJ)
	      PZS(J)=PZJ(JJ)
	      HES(J)=HEJ(JJ)
  187      CONTINUE
            NNNPS=NNNPS+NNNPJ-1
          ENDIF
 5111   CONTINUE
        ENDIF
   20   CONTINUE
        CALL DECHKK(NHARD1)
        NHKKHA=NHKK
C----------------------------------------------------------------
C
      RETURN
 5001 FORMAT (I6,I4,5I6,9E10.2)
 5003 FORMAT (' HADRKK J.GT.NAUMAX SKIP NEXT PARTICLES ',3I10)
 5009 FORMAT (' NHKK,IDHKK(NHKK)  ',3I10)
      END
C
C********************************************************************
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HADRZZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C-------------------------
C
C                       HADRONIZE HARD CHAINS
C
C                       ADD GENERATED HADRONS TO /ALLPAR/
C                          STARTING AT (NAUX + 1)
C                       AND TO /HKKEVT/ STARTING AT (NHKK + 1)
C
C---------------------------------------------------------
      PARAMETER (INTMX=2488,INTMD=252)
      COMMON /ABRZZ/ AMCZZ1(INTMX),AMCZZ2(INTMX),
     *               GACZZ1(INTMX),GACZZ2(INTMX),
     *               BGXZZ1(INTMX),BGYZZ1(INTMX),BGZZZ1(INTMX),
     *               BGXZZ2(INTMX),BGYZZ2(INTMX),BGZZZ2(INTMX),
     *               NCHZZ1(INTMX),NCHZZ2(INTMX),
     *               IJCZZ1(INTMX),IJCZZ2(INTMX),
     *               PQZZA1(INTMX,4),PQZZA2(INTMX,4),
     *               PQZZB1(INTMX,4),PQZZB2(INTMX,4)
      COMMON /ABRSOF/XSQ1(INTMX),XSAQ1(INTMX),XSQ2(INTMX),XSAQ2(INTMX),
     *        IJSQ1(INTMX),IJSAQ1(INTMX),IJSQ2(INTMX),IJSAQ2(INTMX),
     *        AMCCH1(INTMX),AMCCH2(INTMX),GAMCH1(INTMX),GAMCH2(INTMX),
     *        BGCH1(INTMX),BGCH2(INTMX),THECH1(INTMX),THECH2(INTMX),
     *        BGXCH1(INTMX),BGYCH1(INTMX),BGZCH1(INTMX),
     *        BGXCH2(INTMX),BGYCH2(INTMX),BGZCH2(INTMX),
     *        NCH1(INTMX),NCH2(INTMX),IJCH1(INTMX),IJCH2(INTMX),
     *  PSOFA1(INTMX,4),PSOFA2(INTMX,4),PSOFB1(INTMX,4),PSOFB2(INTMX,4)
     * ,JHKKPZ(INTMX),JHKKTZ(INTMX),JHKKSX(INTMX),JHKKS1(INTMX)
      COMMON /NUCJTN/NONUJ1,NONUJT,NONUS1,NONUST
*KEEP,INTNEW.
      COMMON /INTNEW/ NVV,NSV,NVS,NSS,NDV,NVD,NDS,NSD,
     +IXPV,IXPS,IXTV,IXTS, INTVV1(248),
     +INTVV2(248),INTSV1(248),INTSV2(248), INTVS1(248),INTVS2(248),
     +INTSS1(INTMX),INTSS2(INTMX),
     +INTDV1(248),INTDV2(248),INTVD1(248),INTVD2(248),
     +INTDS1(INTMD),INTDS2(INTMD),INTSD1(INTMD),INTSD2(INTMD)
 
C  /INTNEW/
C       NVV    :   NUMBER OF INTERACTING VALENCE-VALENCE SYSTEMS
C       NSV    :   NUMBER OF INTERACTING SEA-VALENCE SYSTEMS
C       NVS    :   NUMBER OF INTERACTING VALENCE-SEA SYSTEMS
C       NSS    :   NUMBER OF INTERACTING SEA-SEA SYSTEMS
C       IXPV, IXTV : NUMBER OF GENERATED X-VALUES FOR VALENCE-QUARK
C                     SYSTEMS FROM PROJECTILE/TARGET NUCLEI
C       IXPS, IXTS : NUMBER OF GENERATED X-VALUES FOR SEA-QUARK PAIRS
C                     FROM PROJECTILE/TARGET NUCLEI
C-------------------
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX),
     *                IFROVT(248),ITOVT(248),IFROST(INTMX),
     *            JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),JHKKNT(248),
     *                JHKKPV(INTMX),JHKKPS(INTMX),
     *                JHKKTV(INTMX),JHKKTS(INTMX),
     *                MHKKVV(INTMX),MHKKSS(INTMX),
     &                MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
C-------------------
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
C.....................................................................
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/  ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),
     *ZUOST(INTMX),
     *                INTLO(INTMX),INLOSS(INTMX)
C  /LOZUO/
C       INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C
C
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
       COMMON /NCOUCH/ ACOUVV,ACOUSS,ACOUSV,ACOUVS,
     *                 ACOUZZ,ACOUHH,ACOUDS,ACOUSD,
     *                 ACOUDZ,ACOUZD,ACOUDI,
     *                 ACOUDV,ACOUVD,ACOUCC
C---------------------
      COMMON /PSHOW/ IPSHOW
C     COMMON /HARLUN/ IHARLU,QLUN
      COMMON /HARLUN/ QLUN,IHARLU
      COMMON /JSPAR/PXJ(1000),PYJ(1000),PZJ(1000),HEJ(1000),NNNPJ
      COMMON /JSPA/PXS(40000),PYS(40000),PZS(40000),HES(40000),NNNPS
C--------------------
*KEEP,DFINPA.
      CHARACTER*8 ANF
      PARAMETER (NFIMAX=249)
      COMMON /DFINPA/ ANF(NFIMAX),PXF(NFIMAX),PYF(NFIMAX),PZF(NFIMAX),
     +HEP(NFIMAX),AMF(NFIMAX), ICHF(NFIMAX),IBARF(NFIMAX),NREF(NFIMAX)
      COMMON /DFINPZ/IORMO(NFIMAX),IDAUG1(NFIMAX),IDAUG2(NFIMAX),
     * ISTATH(NFIMAX)
C-------------------
      PARAMETER (NMXHKK=  89998)
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &             JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),PHKK(5,NMXHKK),
     & VHKK(4,NMXHKK),WHKK(4,NMXHKK)
C
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
      COMMON /PROJK/ IPROJK
      COMMON /GLUSPL/NUGLUU,NSGLUU
      COMMON /NOMIJE/ PTMIJE(10),NNMIJE(10)
C
      DIMENSION POJ(4),PAT(4)
      DATA NCALZZ /0/
C-----------------------------------------------------------------------
        DO 20 I=1,NONUST
          IF(NCH1(I).EQ.99.OR.NCH1(I).EQ.88)GO TO 20
          IF(NCH2(I).EQ.99.OR.NCH2(I).EQ.88)GO TO 20
          NCALZZ=NCALZZ+1
C
          IF (IPHKK.GE.7)WRITE(6,7789)NONUST,NCALZZ,jhkksx(i)
 7789     FORMAT (' HADRZZ NONUST,NCALZZ,Jhkksx(i) ',3I10)
          IF (JHKKSX(I).EQ.1)THEN
          IF (I.GT.INTMX)THEN
            WRITE (6,7744)I,INTMX
 7744       FORMAT ('  HADRZZ I.GT.INTMX  ',2I10)
            RETURN
          ENDIF
C
C++++++++++++++++++++++++++++++    CHAIN 1:  QUARK-DIQUARK   +++++++++++
            IFB1=IJSQ1(I)
            IFB2=IJSAQ1(I)
            IFB2=IABS(IFB2)+6
            DO 21 J=1,4
              POJ(J)=PSOFA1(I,J)
              PAT(J)=PSOFA2(I,J)
  21        CONTINUE
        PT1=SQRT(POJ(1)**2+POJ(2)**2)
        PT2=SQRT(PAT(1)**2+PAT(2)**2)
        CALL PARPT(2,PT1,PT2,5,NEVT)
            IHARLU=0
            QLUN=0.
            IF(IPSHOW.EQ.1)THEN
              POJPT=SQRT(POJ(2)**2+POJ(1)**2)
              PATPT=SQRT(PAT(1)**2+PAT(2)**2)
	      DO IIII=1,10
              IF(POJPT.GE.PTMIJE(IIII))NNMIJE(IIII)=
     *	      NNMIJE(IIII)+1
              IF(PATPT.GE.PTMIJE(IIII))NNMIJE(IIII)=
     *	      NNMIJE(IIII)+1
	      ENDDO
              QLUN=MIN(POJPT,PATPT)
              IF((QLUN.LT.2.5D0).OR.(AMCCH1(I).LT.5.D0))THEN
                QLUN=0.
                IHARLU=0
              ELSE 
                IHARLU=1
              ENDIF   
            ENDIF
C----------------------------------------------------------------
            IF (GAMCH1(I).LT.0.001D0)WRITE (6,7788)
     *      I,NHAD,AMCCH1(I),POJ,PAT,GAMCH1(I),BGXCH1(I),
     *      BGYCH1(I),BGZCH1(I),IFB1,IFB2,IFB3,IFB4,JHKKSX(I)
 7788       FORMAT (' HADRZZ ',2I5,10E12.2/3E12.2,5I5)
            CALL HADJET(NHAD,AMCCH1(I),POJ,PAT,GAMCH1(I),BGXCH1(I),
     *              BGYCH1(I),BGZCH1(I),IFB1,IFB2,IFB3,IFB4,
     *              IJCZZ1(i),IJCZZ1(i),3,NCHZZ1(i),23)
            ACOUZZ=ACOUZZ+1
            IHARLU=0
            QLUN=0.
            NHKKAU=NHKK+1
            IF (NHAD.GT.NFIMAX) THEN
              WRITE (6,7755)NHAD,NFIMAX
 7755         FORMAT (' NHAD.GT.NFIMAX ',2I10)
              RETURN
            ENDIF
            DO 22 J=1,NHAD
C             NHKK=NHKK+1
              IF (NHKK.EQ.NMXHKK) THEN
                WRITE (*,'(A,2I5/A)') ' HADRZZ: NHKK.EQ.NMXHKK ',
     *          NHKK,NMXHKK
                RETURN
              ENDIF
C
              EHECC=SQRT(PXF(J)**2+PYF(J)**2+PZF(J)**2+AMF(J)**2)
              IF (ABS(EHECC-HEP(J)).GT.0.001D0) THEN
C               WRITE(*,'(2A/3I5,3E15.6)')
C    &            ' HADRZZ / CHAIN 1 : CORRECT INCONSISTENT ENERGY ',
C    *            '  NCALZZ, NHKK,NREF(J), HEP(J),EHECC, AMF(J)',
C    *            NCALHH, NHKK,NREF(J), HEP(J),EHECC, AMF(J)
                HEP(J)=EHECC
              ENDIF
              ANNZZ=ANNZZ+1.
              EEZZ=EEZZ+HEP(J)
              PTZZ=SQRT(PXF(J)**2+PYF(J)**2)+PTZZ
C                               PUT NN-CMS HADRONS INTO /HKKEVT/
          ISTIST=1
          IF(IBARF(J).EQ.500)ISTIST=2
          CALL HKKFIL(ISTIST,MPDGHA(NREF(J)),MHKKHH(I)-3,0,
     *                PXF(J),PYF(J),PZF(J),HEP(J),NHKKAU,IORMO(J),11)
              IF(IDHKK(NHKK).EQ.99999) WRITE (6,5009)NHKK,NREF(J),
     *                                               IDHKK(NHKK)
C	  WRITE(6,*)' First chain HADRZZ'
              IF (IPHKK.GE.7) WRITE(6,5001) NHKK,
     *        ISTHKK(NHKK),IDHKK(NHKK),JMOHKK(1,NHKK),JMOHKK(2,NHKK),
     &        JDAHKK(1,NHKK),JDAHKK(2,NHKK),(PHKK(KHKK,NHKK),KHKK=1,5),
     &        (VHKK(KHKK,NHKK),KHKK=1,4)
   22       CONTINUE
C           JDAHKK(1,IMOHKK)=NHKKAU
C           JDAHKK(2,IMOHKK)=NHKK
          IF(NNNPJ.GE.1)THEN
            NNNPSO=NNNPS
            NNNPS=NNNPS+1
            NNNPSU=NNNPSO+NNNPJ
            DO 137 J=NNNPS,NNNPSU
              JJ=J-NNNPS+1
              IF(J.GT.40000.OR.JJ.GT.1000)THEN
C               WRITE(6,'(A,2I10)')' J.gt.40000.or.jj.gt.1000 ',J,JJ
                GO TO 137
              ENDIF
	      PXS(J)=PXJ(JJ)
	      PYS(J)=PYJ(JJ)
	      PZS(J)=PZJ(JJ)
	      HES(J)=HEJ(JJ)
  137      CONTINUE
            NNNPS=NNNPS+NNNPJ-1
          ENDIF
C+++++++++++++++++++++++++++++   CHAIN 2:  AQUARK-QUARK  ++++++++++++++
            IFB1=IJSAQ2(I)
            IFB2=IJSQ2(I)
            IFB1=IABS(IFB1)+6
            DO 23 J=1,4
              POJ(J)=PSOFB1(I,J)
              PAT(J)=PSOFB2(I,J)
  23        CONTINUE
        PT1=SQRT(POJ(1)**2+POJ(2)**2)
        PT2=SQRT(PAT(1)**2+PAT(2)**2)
        CALL PARPT(2,PT1,PT2,5,NEVT)
            IHARLU=0
            QLUN=0.
            IF(IPSHOW.EQ.1)THEN
              POJPT=SQRT(POJ(2)**2+POJ(1)**2)
              PATPT=SQRT(PAT(1)**2+PAT(2)**2)
	      DO IIII=1,10
              IF(POJPT.GE.PTMIJE(IIII))NNMIJE(IIII)=
     *	      NNMIJE(IIII)+1
              IF(PATPT.GE.PTMIJE(IIII))NNMIJE(IIII)=
     *	      NNMIJE(IIII)+1
	      ENDDO
              QLUN=MIN(POJPT,PATPT)
              IF((QLUN.LT.2.5D0).OR.(AMCCH2(I).LT.5.D0))THEN
                QLUN=0.
                IHARLU=0
              ELSE 
                IHARLU=1
              ENDIF   
            ENDIF
C                                             TURN 20.8.91
            CALL HADJET(NHAD,AMCCH2(I),PAT,POJ,GAMCH2(I),BGXCH2(I),
     *                BGYCH2(I),BGZCH2(I),IFB1,IFB2,IFB3,IFB4,
     *              IJCZZ2(i),IJCZZ2(i),3,NCHZZ2(i),24)
            IHARLU=0
            QLUN=0.
C                                   ADD HADRONS/RESONANCES INTO
C                                   COMMON /ALLPAR/ STARTING AT NAUX
            NHKKAU=NHKK+1
          DO 24 J=1,NHAD
C           NHKK=NHKK+1
            IF (NHKK.EQ.NMXHKK) THEN
              WRITE (*,'(A,2I5/A)') ' HADRZZ: NHKK.EQ.NMXHKK ',
     &                               NHKK,NMXHKK
              RETURN
            ENDIF
C
            EHECC=SQRT(PXF(J)**2+PYF(J)**2+PZF(J)**2+AMF(J)**2)
            IF (ABS(EHECC-HEP(J)).GT.0.001D0) THEN
C             WRITE(*,'(2A/3I5,3E15.6)')
C    &            ' HADRZZ / CHAIN 2 : CORRECT INCONSISTENT ENERGY ',
C    *            '  NCALZZ, NHKK,NREF(J), HEP(J),EHECC, AMF(J)',
C    *            NCALZZ, NHKK,NREF(J), HEP(J),EHECC, AMF(J)
              HEP(J)=EHECC
            ENDIF
            ANNZZ=ANNZZ+1.
            EEZZ=EEZZ+HEP(J)
            PTZZ=SQRT(PXF(J)**2+PYF(J)**2)+PTZZ
C                               PUT NN-CMS HADRONS INTO /HKKEVT/
          ISTIST=1
          IF(IBARF(J).EQ.500)ISTIST=2
          CALL HKKFIL(ISTIST,MPDGHA(NREF(J)),MHKKHH(I),0,
     *                PXF(J),PYF(J),PZF(J),HEP(J),NHKKAU,IORMO(J),12)
            IF(IDHKK(NHKK).EQ.99999) WRITE (6,5009)NHKK,NREF(J),
     *                                             IDHKK(NHKK)
C	  WRITE(6,*)' Second chain HADRZZ'
            IF (IPHKK.GE.7) WRITE(6,5001) NHKK,
     *      ISTHKK(NHKK),IDHKK(NHKK),JMOHKK(1,NHKK),JMOHKK(2,NHKK),
     &      JDAHKK(1,NHKK),JDAHKK(2,NHKK),(PHKK(KHKK,NHKK),KHKK=1,5),
     &      (VHKK(KHKK,NHKK),KHKK=1,4)
   24     CONTINUE
C         JDAHKK(1,IMOHKK)=NHKKAU
C         JDAHKK(2,IMOHKK)=NHKK
          IF(NNNPJ.GE.1)THEN
            NNNPSO=NNNPS
            NNNPS=NNNPS+1
            NNNPSU=NNNPSO+NNNPJ
            DO 187 J=NNNPS,NNNPSU
              JJ=J-NNNPS+1
              IF(J.GT.40000.OR.JJ.GT.1000)THEN
C               WRITE(6,'(A,2I10)')' J.gt.40000.or.jj.gt.1000 ',J,JJ
                GO TO 187
              ENDIF
	      PXS(J)=PXJ(JJ)
	      PYS(J)=PYJ(JJ)
	      PZS(J)=PZJ(JJ)
	      HES(J)=HEJ(JJ)
  187      CONTINUE
            NNNPS=NNNPS+NNNPJ-1
          ENDIF
 5111   CONTINUE
        ENDIF
   20   CONTINUE
C----------------------------------------------------------------
C
      RETURN
 5001 FORMAT (I6,I4,5I6,9E10.2)
 5003 FORMAT (' HADRKK J.GT.NAUMAX SKIP NEXT PARTICLES ',3I10)
 5009 FORMAT (' NHKK,IDHKK(NHKK)  ',3I10)
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      SUBROUTINE QINNUC(X,Y)

C  Este programa genera una distribucion de partones de tipo gaussiana
C  centrada en el centro del hadron. Distribucion: F(b)=A*(-b**2/c).
C  La distribucion la generamos en coordenadas polares porque asi
C  tenemos primitiva.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE

      CHARACTER*80 TITLE
      CHARACTER*8 PROJTY,TARGTY
C     COMMON/USER/TITLE,PROJTY,TARGTY,CMENER,ISTRUF
C    &           ,ISINGD,IDUBLD,SDFRAC,PTLAR
      COMMON /USER1/TITLE,PROJTY,TARGTY
      COMMON /USER2/CMENER,SDFRAC,PTLAR,ISTRUF,ISINGD,IDUBLD
                
      C=4.*(0.15D-24+0.01D-24*LOG(CMENER))
  10  P=RNDM(V1)
      IF ((P). EQ .(1.D0)) THEN
         GO TO 10
      END IF
      Z=RNDM(V2)
      T=2.*3.1416*Z
      R=DSQRT(-C*DLOG(1.D00-P))
      X=R*DCOS(T)
      Y=R*DSIN(T)

      RETURN
      END

C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE CASAVS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C-------------------------
C
C                       Casado diquarks VS
C
C                       ADD GENERATED HADRONS TO /ALLPAR/
C                          STARTING AT (NAUX + 1)
C                       AND TO /HKKEVT/ STARTING AT (NHKK + 1)
C
C-------------------------
*KEEP,INTMX.
      PARAMETER (INTMX=2488,INTMD=252)
*KEEP,DXQX.
C     INCLUDE (XQXQ)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON /DXQX/ XPVQ(248),XPVD(248),XTVQ(248),XTVD(248), XPSQ
     +(INTMX),XPSAQ(INTMX),XTSQ(INTMX),XTSAQ(INTMX)
     *                ,XPSU(248),XTSU(248)
     *                ,XPSUT(248),XTSUT(248)
       COMMON/POPCCK/PDBCK,PDBSE,PDBSEU,
     *  IJPOCK,IREJCK,ICK4,IHAD4,ICK6,IHAD6
     *,IREJSE,ISE4,ISE6,IREJS3,ISE43,ISE63,IREJS0
     *,IHADA4,IHADA6,IREJSA,ISEA4,ISEA6,IREJA3,
     *ISEA43,ISEA63,IREJAO
*KEEP,INTNEW.
      COMMON /INTNEW/ NVV,NSV,NVS,NSS,NDV,NVD,NDS,NSD,
     +IXPV,IXPS,IXTV,IXTS, INTVV1(248),
     +INTVV2(248),INTSV1(248),INTSV2(248), INTVS1(248),INTVS2(248),
     +INTSS1(INTMX),INTSS2(INTMX),
     +INTDV1(248),INTDV2(248),INTVD1(248),INTVD2(248),
     +INTDS1(INTMD),INTDS2(INTMD),INTSD1(INTMD),INTSD2(INTMD)
 
C  /INTNEW/
C       NVV    :   NUMBER OF INTERACTING VALENCE-VALENCE SYSTEMS
C       NSV    :   NUMBER OF INTERACTING SEA-VALENCE SYSTEMS
C       NVS    :   NUMBER OF INTERACTING VALENCE-SEA SYSTEMS
C       NSS    :   NUMBER OF INTERACTING SEA-SEA SYSTEMS
C       IXPV, IXTV : NUMBER OF GENERATED X-VALUES FOR VALENCE-QUARK
C                     SYSTEMS FROM PROJECTILE/TARGET NUCLEI
C       IXPS, IXTS : NUMBER OF GENERATED X-VALUES FOR SEA-QUARK PAIRS
C                     FROM PROJECTILE/TARGET NUCLEI
C-------------------
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,LOZUO.
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/ ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),ZUOST(INTMX),
     +INTLO(INTMX),INLOSS(INTMX)
C  /LOZUO/
C       INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C------------------
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
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
*KEEP,ABRVS.
      COMMON /ABRVS/ AMCVS1(248),AMCVS2(248),GACVS1(248),GACVS2(248),
     +BGXVS1(248),BGYVS1(248),BGZVS1(248), BGXVS2(248),BGYVS2(248),
     +BGZVS2(248), NCHVS1(248),NCHVS2(248),IJCVS1(248),IJCVS2(248),
     +PQVSA1(248,4),PQVSA2(248,4), PQVSB1(248,4),PQVSB2(248,4)
*KEEP,DFINPA.
      CHARACTER*8 ANF
      PARAMETER (NFIMAX=249)
      COMMON /DFINPA/ ANF(NFIMAX),PXF(NFIMAX),PYF(NFIMAX),PZF(NFIMAX),
     +HEF(NFIMAX),AMF(NFIMAX), ICHF(NFIMAX),IBARF(NFIMAX),NREF(NFIMAX)
      COMMON /DFINPZ/IORMO(NFIMAX),IDAUG1(NFIMAX),IDAUG2(NFIMAX),
     * ISTATH(NFIMAX)
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEEP,PROJK.
      COMMON /PROJK/ IPROJK
*KEEP,NUCC.
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
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
       COMMON /NCOUCH/ ACOUVV,ACOUSS,ACOUSV,ACOUVS,
     *                 ACOUZZ,ACOUHH,ACOUDS,ACOUSD,
     *                 ACOUDZ,ACOUZD,ACOUDI,
     *                 ACOUDV,ACOUVD,ACOUCC
C---------------------
      COMMON /ZSEA/ZSEAAV,ZSEASU,ANZSEA
      COMMON /CASADI/CASAXX,ICASAD
C---------------------
C-----------------------------------------------------------------------
      DO 50 I=1,NVS
C-----------------------drop recombined chain pairs
        IF(NCHVS1(I).EQ.99.AND.NCHVS2(I).EQ.99) GO TO 50
        IS1=INTVS1(I)
        IS2=INTVS2(I)
C
        IF (IPCO.GE.6) WRITE (6,1010) IPVQ(IS1),IPPV1(IS1),IPPV2(IS1),
     +  ITSQ(IS2),ITSAQ(IS2), AMCVS1(I),AMCVS2(I),GACVS1(I),GACVS2(I),
     +  BGXVS1(I),BGYVS1(I),BGZVS1(I), BGXVS2(I),BGYVS2(I),BGZVS2(I),
     +  NCHVS1(I),NCHVS2(I),IJCVS1(I),IJCVS2(I), PQVSA1(I,4),PQVSA2
     +  (I,4),PQVSB1(I,4),PQVSB2(I,4)
C
C
C++++++++++++++++++++++++++++++    CHAIN 2:  DIQUARK-QUARK   +++++++++++
        IFB1=IPPV1(IS1)
        IFB2=IPPV2(IS1)
        IFB3=ITSQ(IS2)
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
C               I= number of valence chain     
C               Projectile Nr ippp= IFROVP(INTVS1(I))
C          No of Glauber sea q at Projectile JIPP=JSSHS(IPP)
       IPPP = IFROVP(INTVS1(I))
       JIPP=JSSHS(IPPP)
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
      IF(IPCO.GE.1)THEN
        WRITE(6,*)' VS qq-q ,IFB1,IFB2,IFB3,',
     *	'INTVS1=IS1,INTVS2=IS2,JIPP,JITTX',
     *	IFB1,IFB2,IFB3,INTVS1(I),INTVS2(I),JIPP,JITTX
        WRITE (6,*)' target sea quark IFB3=',IFB3,
     *	' from IS2=',INTVS2(I)
        WRITE(6,*)' with ITSQ(IS2),XTSQ(IS2),IFROST(IS2)',
     *	ITSQ(IS2),XTSQ(IS2),IFROST(IS2)
       ENDIF
        DO 797 II=1,IXTV
	  IF(IFROST(IS2).EQ.IFROVT(II))III=II
  797   CONTINUE	
      IF(IPCO.GE.1)THEN
        WRITE (6,*)' projectile III=',III
        WRITE(6,*)' corresp. XTVQ(i),XTVD(i),ITVQ(I),ITTV1(I),ITTV2(I)',
     *   XTVQ(III),XTVD(III),ITVQ(III),ITTV1(III),ITTV2(III)
       ENDIF
C------------------------------------------------------------------- 
C                         Casado diquark option
C+++++++++++++++++++++++++++++ VS   CHAIN 2:  DIQUARK-QUARK   +++++++++
C-------------------------------------------------------------------    
       IF(ICASAD.EQ.1)THEN
         IF(RNDM(VV).LE.CASAXX)THEN
	   IF(RNDM(VVV).LE.0.5D0)THEN
	     ISCASA=ITSQ(IS2) 
	     ITVCAS=ITTV1(III)
	     ITSQ(IS2)=ITVCAS
	     ITTV1(III)=ISCASA
	     IFB3=ITSQ(IS2)
      IF(IPCO.GE.1)THEN
        WRITE(6,*)' Cas VS2 qq-q 1 ,IFB1,IFB2,IFB3,',
     *	'INTVS1=IS1,INTVS2=IS2,III',
     *	IFB1,IFB2,IFB3,INTVS1(I),INTVS2(I),III
     *  ,'-----------------------------------------------------'
       ENDIF
	   ELSE
	     ISCASA=ITSQ(IS2) 
	     ITVCAS=ITTV2(III)
	     ITSQ(IS2)=ITVCAS
	     ITTV2(III)=ISCASA
	     IFB3=ITSQ(IS2)
      IF(IPCO.GE.1)THEN
        WRITE(6,*)' Cas VS2 qq-q 2 ,IFB1,IFB2,IFB3,',
     *	'INTVS1=IS1,INTVS2=IS2,III',
     *	IFB1,IFB2,IFB3,INTVS1(I),INTVS2(I),III
     *  ,'-----------------------------------------------------'
       ENDIF
	   ENDIF
	 ENDIF
       ENDIF
C------------------------------------------------------------------- 
C                         Casado diquark option
C-------------------------------------------------------------------    
   50 CONTINUE
C
      RETURN
 1010 FORMAT(10X,5I5,10F9.2/10X,4I5,4F12.4)
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE CASASV
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C-------------------------
C
C                       Casado diquarks SV
C
C                       ADD GENERATED HADRONS TO /ALLPAR/
C                          STARTING AT (NAUX + 1)
C                       AND TO /HKKEVT/ STARTING AT (NHKK + 1)
C
C---------------------------------------------------------
*KEEP,INTMX.
      PARAMETER (INTMX=2488,INTMD=252)
*KEEP,DXQX.
C     INCLUDE (XQXQ)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON /DXQX/ XPVQ(248),XPVD(248),XTVQ(248),XTVD(248), XPSQ
     +(INTMX),XPSAQ(INTMX),XTSQ(INTMX),XTSAQ(INTMX)
     *                ,XPSU(248),XTSU(248)
     *                ,XPSUT(248),XTSUT(248)
       COMMON/POPCCK/PDBCK,PDBSE,PDBSEU,
     *  IJPOCK,IREJCK,ICK4,IHAD4,ICK6,IHAD6
     *,IREJSE,ISE4,ISE6,IREJS3,ISE43,ISE63,IREJS0
     *,IHADA4,IHADA6,IREJSA,ISEA4,ISEA6,IREJA3,
     *ISEA43,ISEA63,IREJAO
*KEEP,INTNEW.
      COMMON /INTNEW/ NVV,NSV,NVS,NSS,NDV,NVD,NDS,NSD,
     +IXPV,IXPS,IXTV,IXTS, INTVV1(248),
     +INTVV2(248),INTSV1(248),INTSV2(248), INTVS1(248),INTVS2(248),
     +INTSS1(INTMX),INTSS2(INTMX),
     +INTDV1(248),INTDV2(248),INTVD1(248),INTVD2(248),
     +INTDS1(INTMD),INTDS2(INTMD),INTSD1(INTMD),INTSD2(INTMD)
 
C  /INTNEW/
C       NVV    :   NUMBER OF INTERACTING VALENCE-VALENCE SYSTEMS
C       NSV    :   NUMBER OF INTERACTING SEA-VALENCE SYSTEMS
C       NVS    :   NUMBER OF INTERACTING VALENCE-SEA SYSTEMS
C       NSS    :   NUMBER OF INTERACTING SEA-SEA SYSTEMS
C       IXPV, IXTV : NUMBER OF GENERATED X-VALUES FOR VALENCE-QUARK
C                     SYSTEMS FROM PROJECTILE/TARGET NUCLEI
C       IXPS, IXTS : NUMBER OF GENERATED X-VALUES FOR SEA-QUARK PAIRS
C                     FROM PROJECTILE/TARGET NUCLEI
C-------------------
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,LOZUO.
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/ ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),ZUOST(INTMX),
     +INTLO(INTMX),INLOSS(INTMX)
C  /LOZUO/
C       INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C------------------
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
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
*KEEP,ABRSV.
      COMMON /ABRSV/ AMCSV1(248),AMCSV2(248),GACSV1(248),GACSV2(248),
     +BGXSV1(248),BGYSV1(248),BGZSV1(248), BGXSV2(248),BGYSV2(248),
     +BGZSV2(248), NCHSV1(248),NCHSV2(248),IJCSV1(248),IJCSV2(248),
     +PQSVA1(248,4),PQSVA2(248,4), PQSVB1(248,4),PQSVB2(248,4)
*KEEP,DFINPA.
      CHARACTER*8 ANF
      PARAMETER (NFIMAX=249)
      COMMON /DFINPA/ ANF(NFIMAX),PXF(NFIMAX),PYF(NFIMAX),PZF(NFIMAX),
     +HEF(NFIMAX),AMF(NFIMAX), ICHF(NFIMAX),IBARF(NFIMAX),NREF(NFIMAX)
      COMMON /DFINPZ/IORMO(NFIMAX),IDAUG1(NFIMAX),IDAUG2(NFIMAX),
     * ISTATH(NFIMAX)
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEEP,PROJK.
      COMMON /PROJK/ IPROJK
*KEEP,NUCC.
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
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
       COMMON /NCOUCH/ ACOUVV,ACOUSS,ACOUSV,ACOUVS,
     *                 ACOUZZ,ACOUHH,ACOUDS,ACOUSD,
     *                 ACOUDZ,ACOUZD,ACOUDI,
     *                 ACOUDV,ACOUVD,ACOUCC
C---------------------
      COMMON /ZSEA/ZSEAAV,ZSEASU,ANZSEA
      COMMON /CASADI/CASAXX,ICASAD
C---------------------
      DATA NCALSV /0/
C-----------------------------------------------------------------------
      NCALSV=NCALSV+1
      DO 50 I=1,NSV
C-----------------------drop recombined chain pairs
        IF(NCHSV1(I).EQ.99.AND.NCHSV2(I).EQ.99) GO TO 50
        IS1=INTSV1(I)
        IS2=INTSV2(I)
C
        IF (IPCO.GE.6) WRITE (6,1000) IPSQ(IS1),IPSAQ(IS1),ITVQ(IS2),
     +  ITTV1(IS2),ITTV2(IS2), AMCSV1(I),AMCSV2(I),GACSV1(I),GACSV2(I),
     +  BGXSV1(I),BGYSV1(I),BGZSV1(I), BGXSV2(I),BGYSV2(I),BGZSV2(I),
     +  NCHSV1(I),NCHSV2(I),IJCSV1(I),IJCSV2(I), PQSVA1(I,4),PQSVA2
     +  (I,4),PQSVB1(I,4),PQSVB2(I,4)
 1000 FORMAT(10X,5I5,10F9.2/10X,4I5,4F12.4)
C
C++++++++++++++++++++++++++++++    CHAIN 1:  QUARK-DIQUARK   +++++++++++
        IFB1=IPSQ(IS1)
        IFB2=ITTV1(IS2)
        IFB3=ITTV2(IS2)
C------------------------------------------------------------------
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
C               I= number of valence chain     
C               Target     Nr itt = IFROVT(INTSV2(I))
C          No of Glauber sea q at Target     JITT=JTSHS(ITT)
       ITTT = IFROVT(INTSV2(I))
       JITT=JTSHS(ITTT)
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
      IF(IPCO.GE.1)THEN
        WRITE(6,*)' SV q-qq ,IFB1,IFB2,IFB3,',
     *	'INTSV1=IS1,INTSV2=IS2,JIPPX,JITT',
     *	IFB1,IFB2,IFB3,INTSV1(I),INTSV2(I),JIPPX,JITT
        WRITE (6,*)' projectile sea quark IFB1=',IFB1,
     *	' from IS1=',INTSV1(I)
        WRITE(6,*)' with IPSQ(IS1),XPSQ(IS1),IFROSP(IS1)',
     *	IPSQ(IS1),XPSQ(IS1),IFROSP(IS1)
       ENDIF
        DO 798 II=1,IXPV
	  IF(IFROSP(IS1).EQ.IFROVP(II))III=II
  798   CONTINUE	
      IF(IPCO.GE.1)THEN
        WRITE (6,*)' projectile III=',III
        WRITE(6,*)' corresp. XPVQ(i),XPVD(i),IPVQ(I),IPPV1(I),IPPV2(I)',
     *   XPVQ(III),XPVD(III),IPVQ(III),IPPV1(III),IPPV2(III)
       ENDIF
C------------------------------------------------------------------- 
C                         Casado diquark option
C++++++++++++++++++++++++++++ SV   CHAIN 1:  QUARK-DIQUARK   +++++++++++
C-------------------------------------------------------------------    
       IF(ICASAD.EQ.1)THEN
         IF(RNDM(VV).LE.CASAXX)THEN
	   IF(RNDM(VVV).LE.0.5D0)THEN
	     ISCASA=IPSQ(IS1) 
	     IPVCAS=IPPV1(III)
	     IPSQ(IS1)=IPVCAS
	     IPPV1(III)=ISCASA
	     IFB1=IPSQ(IS1)
      IF(IPCO.GE.1)THEN
        WRITE(6,*)' Cas SV1 q-qq 1 ,IFB1,IFB2,IFB3,',
     *	'INTSV1=IS1,INTSV2=IS2,JIPPX,JITT,III',
     *	IFB1,IFB2,IFB3,INTSV1(I),INTSV2(I),JIPPX,JITT,III
     *  ,'-----------------------------------------------------'
       ENDIF
	   ELSE
	     ISCASA=IPSQ(IS1) 
	     IPVCAS=IPPV2(III)
	     IPSQ(IS1)=IPVCAS
	     IPPV2(III)=ISCASA
	     IFB1=IPSQ(IS1)
      IF(IPCO.GE.1)THEN
        WRITE(6,*)' Cas SV1 q-qq 2 ,IFB1,IFB2,IFB3,',
     *	'INTSV1=IS1,INTSV2=IS2,JIPPX,JITT,III',
     *	IFB1,IFB2,IFB3,INTSV1(I),INTSV2(I),JIPPX,JITT,III
     *  ,'-----------------------------------------------------'
       ENDIF
	   ENDIF
	 ENDIF
       ENDIF
C------------------------------------------------------------------- 
C                         Casado diquark option
C-------------------------------------------------------------------    
   50 CONTINUE
C----------------------------------------------------------------
C
      RETURN
      END
