*-- Author :
C
C--------------------------------------------------------------------
      SUBROUTINE FHAD(IPRMOD,IPRO,PLAB,ELAB,CX,CY,CZ,
     *               ITHKK,ITTA,IELINE,IREJFH)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C   MODIFIED VERSION OF FERHAD / FROM FLUKA86.S (DESY)
C   01/02/90
C
C--------------------------------------------------
C***  COLLISION OF HADRON IPRO WITH NUCLEON ITHKK FROM /HKKEVT/
C         (IPRO, ITTA - CONVENTIONAL PARTICLE NUMBERING FROM FLUKA)
C     IPRO HAS LAB-ENERGY ELAB, MOMENTUM PLAB, DIRECTIONS CX,CY,CZ
C
C                      IELINI=0    INELASTIC HADRIN COLLISIONS
C                      IELINE=1      ELASTIC ELHAIN COLLISIONS
C                      IELINE= ...       ...
C***
C***  ITHKK TAKES THE FERMI-MOMENTUM FROM /HKKEVT/
C
C*** CONSERVED IS THE ENERGY, THE MOMENTUM, ELECTRIC AND BARYON. CHARGE
C*** AND STRANGENESS
C--------------------------------------------------
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
*KEEP,HADTHR.
      COMMON /HADTHR/ EHADTH,INTHAD
*KEEP,DFINLS.
      PARAMETER (MAXFIN=10)
      COMMON /DFINLS/ ITRH(MAXFIN),CXRH(MAXFIN),CYRH(MAXFIN), CZRH
     +(MAXFIN),ELRH(MAXFIN),PLRH(MAXFIN),IRH
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
*KEEP,PROJK.
      COMMON /PROJK/ IPROJK
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEND.
      ZERO=0
      IREJFH=0
C---------------------------------------------------------------
C   TEST OPTION: FORCE TYPE OF INTERACTION IF DESIRED
      IF(INTHAD.EQ.1) THEN
        IIELIN=0
      ELSEIF(INTHAD.EQ.2) THEN
        IIELIN=1
      ELSE
        IIELIN=IELINE
      ENDIF
C--------------------------------------------------
C*** COLLISION KINEMATICS
C*** LORENTZ-TRANSFORMATION INTO TARGET-NUCLEON-REST-SYSTEM
C--------------------------------------------------
      AMTAR=PHKK(5,ITHKK)
C                               RECALCULATE MOMENTA FOR KINEMATICAL
C                               CONSISTENCY OF LORENTZ-TRANSFORMATION
      PPTAR=SQRT(PHKK(1,ITHKK)**2 + PHKK(2,ITHKK)**2 + PHKK(3,ITHKK)**2)
      PPTARM=SQRT(ABS(PHKK(4,ITHKK)-AMTAR)*(PHKK(4,ITHKK)+AMTAR))
C
      BGX=PHKK(1,ITHKK)/AMTAR
      BGY=PHKK(2,ITHKK)/AMTAR
      BGZ=PHKK(3,ITHKK)/AMTAR
      IF(PPTAR.GT.ZERO) THEN
        BGX=BGX*PPTARM/PPTAR
        BGY=BGY*PPTARM/PPTAR
        BGZ=BGZ*PPTARM/PPTAR
      ENDIF
C
      PPLAB=PLAB
      GAM=PHKK(4,ITHKK)/AMTAR
      PXPRO=CX*PPLAB
      PYPRO=CY*PPLAB
      PZPRO=CZ*PPLAB
C
      ETES = ELAB + PHKK(4,ITHKK)
      PXTES= PXPRO + PHKK(1,ITHKK)
      PYTES= PYPRO + PHKK(2,ITHKK)
      PZTES= PZPRO + PHKK(3,ITHKK)
C
      CALL DALTRA(GAM,-BGX,-BGY,-BGZ, PXPRO,PYPRO,PZPRO,ELAB, PPROF,
     +PXPROF,PYPROF,PZPROF,EPROF)
      IF(IPEV.GT.3) THEN
        WRITE(6,'(2A/A,4E12.5/A,5E12.5/A,4E12.5/A,4E12.5/A,4I4)')
     &        ' FHAD:  projectile after LT into target',
     &        ' nucleon rest system',
     &        ' GAM,BGX,BGY,BGZ ', GAM,BGX,BGY,BGZ,
     &        ' PHKK(target):', (PHKK(IK,ITHKK),IK=1,5),
     &        ' PX/Y/Z/PRO, ELAB:', PXPRO,PYPRO,PZPRO,ELAB,
     &        ' Proj. after LT:', PXPROF,PYPROF,PZPROF,EPROF,
     &        ' IPRO, ITTA, ITHKK, IELINE:', IPRO,ITTA,ITHKK,IELINE
      ENDIF
      IF(EPROF.LE.AAM(IPRO)) THEN
        WRITE(6,'(2A/A,5E12.5/A,4E12.5/A,4E12.5/A,4I4)')
     &        ' FHAD: inconsistent projectile after LT into target',
     &        ' nucleon rest system',
     &        ' PHKK(target):', (PHKK(IK,ITHKK),IK=1,5),
     &        ' PX/Y/Z/PRO, ELAB:', PXPRO,PYPRO,PZPRO,ELAB,
     &        ' Proj. after LT:', PXPROF,PYPROF,PZPROF,EPROF,
     &        ' IPRO, ITTA, ITHKK, IELINE:', IPRO,ITTA,ITHKK,IELINE
        IREJFH=1
        RETURN
      ENDIF
C                                     CONSISTENCY TEST OF LOR. TRSF.
C                                     INTO TARGET REST SYSTEM
      IF(IPAUPR.GT.5) THEN
        CALL DALTRA(GAM,-BGX,-BGY,-BGZ, PHKK(1,ITHKK),PHKK(2,ITHKK),PHKK
     +  (3,ITHKK),PHKK(4,ITHKK), PPTAR,PXTARF,PYTARF,PZTARF,ETARF)
 
        WRITE(6,'(A)') 'FHAD: TARGET MOM. BEFORE/AFTER LORENTZ TRANSF.'
        WRITE(6,'(3X,A,5(1PE12.4))') ' PHKK(1-4,ITHKK)', (PHKK
     +  (JJ,ITHKK),JJ=1,4)
        WRITE(6,'(3X,A,5(1PE12.4))') ' PXTARF,PYTARF,PZTARF,ETARF',
     +  PXTARF,PYTARF,PZTARF,ETARF
        WRITE(6,'(A)') 'FHAD: PROJ. MOM. BEFORE/AFTER LORENTZ TRANSF.'
        WRITE(6,'(3X,A,5(1PE12.4))') ' PXPRO,PYPRO,PZPRO,ELAB', PXPRO,
     +  PYPRO,PZPRO,ELAB
        WRITE(6,'(3X,A,5(1PE12.4))') ' PXPROF,PYPROF,PZPROF,EPROF',
     +  PXPROF,PYPROF,PZPROF,EPROF
      ENDIF
C
C--------------------------------------------------
C***  FOR PARTICLES OF THE H-N-COLLISION, STORE THE KINEM.VARIABLES IN
C*** COMMON /FINLSP/
C--------------------------------------------------
      CXF=PXPROF/PPROF
      CYF=PYPROF/PPROF
      CZF=PZPROF/PPROF
      IREJ=0
      IF (IIELIN.EQ.0)THEN
        CALL DHADRI(IPRO,PPROF,EPROF,CXF,CYF,CZF,ITTA)
        IF(IRH.EQ.1) IREJ=1
        IF (IPAUPR.GT.2)WRITE(6,1000)IPRO,PPROF,EPROF,CXF,CYF,CZF,ITTA
 1000   FORMAT (' FHAD IPRO,PPROF,EPROF,CXF,CYF,CZF,ITTA',I5,5F10.2,I5)
      ELSEIF(IIELIN.EQ.1) THEN
        CALL ELHAIN(IPRO,PPROF,EPROF,CXF,CYF,CZF,ITTA,IREJ)
      ENDIF
      IF(IREJ.EQ.1) THEN
C                           RETURN ORIGINAL MOMENTA (SEE ELHAIN)
        IRH=2
        ITRH(1)=IPRO
        CXRH(1)=CX
        CYRH(1)=CY
        CZRH(1)=CZ
        ELRH(1)=ELAB
        PLRH(1)=PLAB
        ITRH(2)=ITTA
        CXRH(2)=PHKK(1,ITHKK)/PPTAR
        CYRH(2)=PHKK(2,ITHKK)/PPTAR
        CZRH(2)=PHKK(3,ITHKK)/PPTAR
        ELRH(2)=PHKK(4,ITHKK)
        PLRH(2)=PPTAR
        RETURN
      ENDIF
C
C--------------------------------------------------
C*** LORENTZ-TRANSFORM FROM TRS INTO LS
C--------------------------------------------------
      DO 10 III=1,IRH
      CRSUM=CXRH(III)**2 + CYRH(III)**2 + CZRH(III)**2
        IF(ABS(CRSUM-1.0).GT.1E-4) THEN
          WRITE(6,'(A,I3,1PE12.4)')
     +    ' FHAD: INCORRECT NORM. OF DIRECTION COSINES - III,CRSUM',
     +    III,CRSUM
          RCRSUM=SQRT(CRSUM)
        CXRH(III)=CXRH(III)/RCRSUM
        CYRH(III)=CYRH(III)/RCRSUM
        CZRH(III)=CZRH(III)/RCRSUM
        ENDIF
      AMI=AAM(ITRH(III))
C?      PPS=SQRT(ABS((ELRH(III)-AMI)*(ELRH(III)+AMI))+1.E-6)
      PPS=PLRH(III)
      PSX=CXRH(III)*PPS
      PSY=CYRH(III)*PPS
      PSZ=CZRH(III)*PPS
        IF(IPAUPR.GT.7) THEN
          WRITE(6,'(A,I3,6(1PE12.4))')
     +    ' FHAD: ITRH(I), AMI,ELR,PPS,PSX,PSY,PSZ', ITRH(III),AMI,ELRH
     +    (III),PPS,PSX,PSY,PSZ
        ENDIF
      CALL DALTRA(GAM,BGX,BGY,BGZ, PSX,PSY,PSZ,ELRH(III), PPPS,PPSX,
     +  PPSY,PPSZ,ELRH(III))
      IF(ELRH(III).LT.(PPPS-1D-4)) THEN
          WRITE(6,'(2A/3I3,6(1PE12.4))')
     +    ' FHAD: INCONSISTENT KINEMATICS AFTER ALTRA: ',
     +    ' IIELIN,III,ITRH(III),ELRH(III),PPPS,PSX,PSY,PSZ,AMI',IIELIN,
     +    III,ITRH(III),ELRH(III),PPPS,PSX,PSY,PSZ,AMI
        ELRH(III)=SQRT(PPPS**2 + AMI**2)
        WRITE(6,'(A,1PE12.4)') ' CORRECTED ENERGY ELRH:',ELRH(III)
      WRITE(6,'(2A/2I5,5(1PE12.4))')
     +    '      FHAD: 4-MOM. OF TARGET NUCLEON -',
     +    ' ITHKK, IDHKK , PHKK(1-4)', ITHKK, IDHKK(ITHKK), (PHKK
     +    (K,ITHKK),K=1,5)
        ENDIF
      CXRH(III)=PPSX/PPPS
      CYRH(III)=PPSY/PPPS
      CZRH(III)=PPSZ/PPPS
      PLRH(III)=PPPS
      ETES = ETES - ELRH(III)
        PXTES= PXTES- PPSX
        PYTES= PYTES- PPSY
        PZTES= PZTES- PPSZ
   10 CONTINUE
C
      IF(ABS(ETES).GT.0.001D0) THEN
C     IF(ABS(ETES).GT.0.1041D0) THEN
        IF(IPRI.GE.1) THEN
          WRITE(6,'(A,I5)') ' FHAD: TEST OF E-P CONSERVATION IELINE=',
     +    IIELIN
          WRITE(6,'(3X,A,5(1PE12.4))') ' ETES,PXTES,PYTES,PZTES', ETES,
     +    PXTES,PYTES,PZTES
          WRITE(6,1000)IPRO,PPROF,EPROF,CXF,CYF,CZF,ITTA
        ENDIF
      DO 20 III=1,IRH
        AMI=AAM(ITRH(III))
        PPS=SQRT((ELRH(III)-AMI)*(ELRH(III)+AMI))
        PSX=CXRH(III)*PPS
        PSY=CYRH(III)*PPS
        PSZ=CZRH(III)*PPS
          IF(IPRI.GE.1) THEN
            WRITE(6,'(A,I3,6(1PE12.4))')
     +      ' FHAD: ITRH(I), AMI,ELRH,PPS,PSX,PSY,PSZ', ITRH(III),AMI,
     +      ELRH(III),PPS,PSX,PSY,PSZ
          ENDIF
   20   CONTINUE
      ENDIF
C
      RETURN
      END
**sr 19-11-95: ELHAIN replaced
*
*===elhain=============================================================*
*
      SUBROUTINE ELHAIN(IP,PLA,ELAB,CX,CY,CZ,IT,IREJ)

************************************************************************
* Elastic hadron-hadron scattering.                                    *
* This is a revised version of the original.                           *
* This version dated 26.10.95 is written by S. Roesler                 *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)
      PARAMETER (TWO=2.0D0,ONE=1.0D0,OHALF=0.5D0,ZERO=0.0D0,
     &           TINY10=1.0D-10)

      PARAMETER (ENNTHR = 3.5D0)
      PARAMETER (PLOWH=0.01D0,PHIH=9.0D0,
     &           BLOWB=0.05D0,BHIB=0.2D0,
     &           BLOWM=0.1D0, BHIM=2.0D0)

      CHARACTER*8 ANAME
      COMMON /DPAR/   ANAME(210),AAM(210),GA(210),TAU(210),
     &                IICH(210),IIBAR(210),K1(210),K2(210)

      PARAMETER (MAXFIN=10)
      COMMON /DFINLS/ ITRH(MAXFIN),CXRH(MAXFIN),CYRH(MAXFIN),
     &                CZRH(MAXFIN),ELRH(MAXFIN),PLRH(MAXFIN),IRH

C     DATA TSLOPE /10.0D0/

      IREJ = 0

      PLAB = SQRT( (ELAB-AAM(IP))*(ELAB+AAM(IP)) )
      EKIN = ELAB-AAM(IP)
*   kinematical quantities in cms of the hadrons
      AMP2 = AAM(IP)**2
      AMT2 = AAM(IT)**2
      S    = AMP2+AMT2+TWO*ELAB*AAM(IT)
      ECM  = SQRT(S)
      ECMP = OHALF*ECM+(AMP2-AMT2)/(TWO*ECM)
      PCM  = SQRT( (ECMP-AAM(IP))*(ECMP+AAM(IP)) )

* nucleon-nucleon scattering at E_kin<3.5: use TSAMCS (HETC-KFA)
      IF ( ((IP.EQ.1).OR.(IP.EQ.8)).AND.
     &     ((IT.EQ.1).OR.(IT.EQ.8)).AND.(EKIN.LT.ENNTHR) ) THEN
*   TSAMCS treats pp and np only, therefore change pn into np and
*   nn into pp
         IF (IT.EQ.1) THEN
            KPROJ = IP
         ELSE
            KPROJ = 8
            IF (IP.EQ.8) KPROJ = 1
         ENDIF
         CALL TSAMCS(KPROJ,EKIN,CTCMS)
         T = TWO*PCM**2*(CTCMS-ONE)

* very crude treatment otherwise: sample t from exponential dist.
      ELSE
*   momentum transfer t
         TMAX = TWO*TWO*PCM**2
         RR = (PLAB-PLOWH)/(PHIH-PLOWH)
         IF (IIBAR(IP).NE.0) THEN
            TSLOPE = BLOWB+RR*(BHIB-BLOWB)
         ELSE
            TSLOPE = BLOWM+RR*(BHIM-BLOWM)
         ENDIF
         FMAX = EXP(-TSLOPE*TMAX)-ONE
         R = RNDM(V)
         T = LOG(ONE+R*FMAX+TINY10)/TSLOPE
         IF (T.GT.ZERO) T = LOG(ONE+R*FMAX)/TSLOPE
      ENDIF

*   target hadron in Lab after scattering
      ELRH(2) = (TWO*AMT2-T)/(TWO*AAM(IT))
      PLRH(2) = SQRT( (ELRH(2)-AAM(IT))*(ELRH(2)+AAM(IT)) )
*   projectile hadron in Lab after scattering
      ELRH(1) = ELAB+AAM(IT)-ELRH(2)
      PLRH(1) = SQRT( (ELRH(1)-AAM(IP))*(ELRH(1)+AAM(IP)) )
*   scattering angle of projectile in Lab
      CTLABP = (T-TWO*AMP2+TWO*ELAB*ELRH(1))/(TWO*PLAB*PLRH(1))
      STLABP = SQRT( (ONE-CTLABP)*(ONE+CTLABP) )
      CALL DSFECF(SPLABP,CPLABP)
*   direction cosines of projectile in Lab
**sr mod. for DPMJET: STTRAN-->DRTRAN
      CALL DRTRAN(CX,CY,CZ,CTLABP,STLABP,SPLABP,CPLABP,
     &                          CXRH(1),CYRH(1),CZRH(1))
*   scattering angle of target in Lab
      PLLABT = PLAB-CTLABP*PLRH(1)
      CTLABT = PLLABT/PLRH(2)
      STLABT = SQRT( (ONE-CTLABT)*(ONE+CTLABT) )
*   direction cosines of target in Lab
**sr mod. for DPMJET: STTRAN-->DRTRAN
      CALL DRTRAN(CX,CY,CZ,CTLABT,STLABT,-SPLABP,-CPLABP,
     &                            CXRH(2),CYRH(2),CZRH(2))
*   fill /DFINLS/
      IRH = 2
      ITRH(1) = IP
      ITRH(2) = IT

      RETURN
      END

*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C**************************************************
C
C   ELASTIC CROSS SECTION SUBROUTINES FOR HKK89
C
C   HJM 10/89
C
C******************************************************
C
      SUBROUTINE SIHNEL(IPROJ,ITAR,POO,SIEL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C********************************************************************
C     VERSION BY HJM 10/89
C     LAST CHANGES HJM 30/08/90
C                      - assignment for several particle types
C                      - low-energy parametrization for K+/-, PBAR - P
C
C     USES SIGEL (NEEDS: SIHNIN, SHPTOT, SIHAEL)
C
C     NOTE: TO BE RENEWED URGENTLY!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     INPUT VARIABLES:
C     IPROJ  = INCIENT PARTICLE TYPE
C     ITAR   = TARGET NUCLEON TYPE
C     POO    = PARTICLE MOMENTUM IN GEV/C
C
C     OUTPUT VARIABLES:
C     SIEL   = ELASTIC CROSS SECTION IN MB
C
C
C  CONVENTION: ISOSPIN INVARIANCE FOR
C                          NUCLEON-NUCLEON CROSSS SECTIONS
C                          PI+/- NUCLEON
C                          K+/-  NUCLEON
C              TO USE THE IMPROVED LOW-ENERGY CROSS SECTIONS
C              FROM SIHAEL
C********************************************************************
C
      DIMENSION SIKPP(8),SIKMP(8),SIAPP(8),P(8)
      DATA P /0.3D0,0.4D0,0.5D0,0.6D0,0.8D0,1.D0,1.5D0,2.D0/
      DATA SIKPP / 12.0D0, 12.5D0, 13.0D0, 13.0D0, 12.7D0,
     +  12.0D0, 10.2D0, 6.82D0/
      DATA SIKMP / 42.0D0, 33.0D0, 21.0D0, 16.0D0, 19.0D0,
     + 22.0D0,  9.0D0, 7.5D0 /
      DATA SIAPP / 73.0D0, 70.0D0, 62.0D0, 53.0D0, 48.0D0,
     + 43.0D0, 38.0D0, 33.0D0/
C-----------------------------------------------------------------------
      ZERO=0
      ONEONE=1
      AA=ONEONE
      PPOO=POO
      IPPR=IPROJ
C
      IF(ITAR.EQ.8) THEN
        IF(IPROJ.EQ.13) THEN
          IPPR=14
        ELSEIF(IPROJ.EQ.14) THEN
          IPPR=13
        ELSEIF(IPROJ.EQ.1) THEN
          IPPR=8
        ELSEIF(IPROJ.EQ.8) THEN
          IPPR=1
        ELSEIF(IPROJ.EQ.15) THEN
          IPPR=16
        ELSEIF(IPROJ.EQ.16) THEN
          IPPR=15
        ELSEIF(IPROJ.EQ.24) THEN
          IPPR=25
        ELSEIF(IPROJ.EQ.25) THEN
          IPPR=24
        ENDIF
      ENDIF
C
      IF(IPPR.EQ.9) THEN
        IPPR=2
      ELSEIF(IPPR.EQ.17) THEN
        IPPR=1
      ELSEIF(IPPR.EQ.18) THEN
        IPPR=2
      ELSEIF(IPPR.EQ.24) THEN
        IPPR=15
      ELSEIF(IPPR.EQ.25) THEN
        IPPR=16
      ELSEIF(IPPR.GE.20.AND.IPPR.LE.22) THEN
        IPPR=1
      ENDIF
C-----------------------------------------------------------------
C                                          K+/-, PBAR  - P
C                                          Plab < 10 GeV/c
      IF(IPPR.EQ.15.OR.IPPR.EQ.16.OR.IPPR.EQ.2) THEN
C***
        IF(PPOO.LE.2.0) THEN
C
C    CALCULATE THE MOMENTUM INDEX K FOR INTERPOLATION
C
          DO 10 JK=1,8
            IF(PPOO.LE.P(JK)) THEN
              K=JK
                                                                 GOTO 20
            ENDIF
   10     CONTINUE
          K=8
   20     CONTINUE
          KK=K-1
          IF(K.EQ.1) KK=1
C*
          IF(IPPR.EQ.15) THEN
C                                K+  - P
            S1=SIKPP(KK)
            S2=SIKPP(K)
          ELSEIF(IPPR.EQ.16) THEN
C                                K-  - P
            S1=SIKMP(KK)
            S2=SIKMP(K)
          ELSEIF(IPPR.EQ.2) THEN
C                                PBAR - P
            S1=SIAPP(KK)
            S2=SIAPP(K)
          ELSE
            WRITE(6,'(A)')
     +      ' LOGICAL ERROR IN SIHNEL - EXECUTION STOPPED'
            STOP
          ENDIF
          SIEL=S1 + (S2-S1)*(PPOO-P(KK))/(P(K)-P(KK)+1D-7)
          RETURN
C***
        ELSEIF(PPOO.LE.10.0D0) THEN
          IF(IPPR.EQ.15) THEN
C                                K+  - P
            A1=5.84
            A2=17.2
            AN=-3.06
            A3=0.206
            A4=-1.71
          ELSEIF(IPPR.EQ.16) THEN
C                                K-  - P
            A1=7.24
            A2=46.0
            AN=-4.71
            A3=0.279
            A4=-2.35
          ELSEIF(IPPR.EQ.2) THEN
C                                PBAR - P
            A1=10.6
            A2=53.1
            AN=-1.19
            A3=0.136
            A4=-1.41
          ENDIF
C
          ALP=LOG(PPOO)
          SIEL=A1 + A2*PPOO**AN + A3*ALP**2 + A4*ALP
          RETURN
C***
        ELSE
                                                                 GOTO 30
        ENDIF
      ENDIF
C-----------------------------------------------------------------
   30 CONTINUE
      CALL DSIGE(IPPR,AA,PPOO,SIEL,ZLEL)
      RETURN
      END
*-- Author :
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE DSIGE(IT,AA,POO,SEL,ZL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C********************************************************************
C     VERSION BY                     J. RANFT
C                                    LEIPZIG
C     LAST CHANGE 01. MAY 84 BY      HJM
C                                    LEIPZIG
C !!  LAST CHANGE  30/08/90    HJM:
C !!           NIZL replaced by SIHNIN, i.e. only useful
C !!                 for hadron-proton scattering in  DTUNUC
C
C
C
C     INPUT VARIABLES:
C     IT     = PARTICLE TYPE
C     AA     = ATOMIC WEIGHT OF THE NUCLEUS
C     POO    = PARTICLE MOMENTUM IN GEV/C
C
C     OUTPUT VARIABLES:
C     SI     = ELASTIC CROSS SECTION IN MB
C     ZL     = INTERACTION LENGTH IN G/CM**2
C
C     OTHER IMPORTANT VARIABLES:
C        SIG    = PROTON/NUCLEI CROSS SECTIONS
C        SEG    = PION/NUCLEI CROSS SECTIONS
C        P      = MOMENTUMS FOR WHICH THE CROSS SECTIONS ARE GIVEN IN
C                 SIG AND SEG
C        A      = NUCLEI FOR WHICH THE CROSS SECTIONS ARE GIVEN IN
C                 SIG AND SEG
C        PLAB   = MOMENTUMS FOR WHICH THE TOTAL CROSS SECTIONS ARE
C                 GIVEN IN SITO
C        SITO   = TOTAL HADRON NUCLEON CROSS SECTIONS FOR NUCLEONS,
C                 PIONS, KAONS AND ANTI-NUCLEONS.
C        ALP    =  EXPONENT OF THE PARAMETRIZATION FOR ANTI-PROTONS,
C                  RANTI-NEUTRONS AND KAONS
C        BET    =  MULTIPLIER OF PARAMETRIZATION FOR ANTI-PROTONS,
C                  ANTI-NEUTRONS AND KAONS
C
C     NOTE1: PRESENTLY CROSS SECTIONS ARE ASSUMED TO BE CONSTANT
C     ABOVE 10.0 GEV/C FOR ALL PARTICLES AND
C     BELOW 0.3 GEV/C FOR NUCLEONS AND BELOW 0.13 GEV/C FOR PIONS
C
C     NOTE2: FOR HADRONS OTHER THAN (1=PROTON,2=ANTI PROTON,8=
C     NEUTRON,9=ANTI NEUTRON,13=POSITIVE PION,14=NEGATIVE PION,15=
C     POSITIVE KAON,16=NEGATIVE KAON,24=NEUTRAL KAON,25=NEUTRAL ANTI
C     KAON) SEE TABLE ITT TO SEE THE CORRESPONDANCE
C
C     NOTE3: FOR LEPTONS AND PHOTONS PRACTICALLY ZERO CROSS SECTION
C     IS RETURNED.
C
C********************************************************************
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
*KEND.
C--------------------------------------------------------------
C--------------------------------------------------------------------
      DIMENSION SIG(13,9),SEG(16,9),P(16),A(9),ITT(39)
      DIMENSION PLAB(19),SITO(19,4),ALP(3),BET(3)
      DIMENSION REA(9,9),STOT(9)
      SAVE A, P, SIG, SEG, ITT, SITO, PLAB, ALP, BET, STOT, REA
      DATA A/9.D0,12.D0,27.D0,47.9D0,55.9D0,63.5D0,112.4D0,
     &207.2D0,238.1D0/
      DATA P/.13D0,.19D0,.25D0,.3D0,.4D0,.5D0,.6D0,.8D0,1.D0,
     &1.5D0,2.D0,3.D0,4.D0,5.D0,6.D0,10.D0/
      DATA SIG/ 485.D0,223.D0,112.D0,82.D0,66.D0,78.D0,96.D0,102.D0,
     &100.D0,98.D0,95.D0,90.D0,79.D0,
     (680.D0,348.D0,175.D0,103.D0,84.D0,87.D0,106.D0,112.D0,111.D0,
     &108.D0,107.D0,105.D0,101.D0,
     (1200.D0,738.D0,387.D0,196.D0,191.D0,200.D0,248.D0,264.D0,
     &264.D0,257.D0,252.D0,247.D0,228.D0,
     (1658.D0,1110.D0,635.D0,364.D0,332.D0,356.D0,404.D0,408.D0,
     &407.D0,404.D0,398.D0,396.D0,384.D0,
     (1730.D0,1270.D0,725.D0,400.D0,375.D0,412.D0,495.D0,505.D0,
     &495.D0,492.D0,487.D0,485.D0,475.D0,
     (1875.D0,1470.D0,835.D0,480.D0,450.D0,450.D0,535.D0,
     &580.D0,555.D0,540.D0,535.D0,530.D0,
     (525.D0,2040.D0,2160.D0,1335.D0,850.D0,740.D0,760.D0,880.D0,
     &905.D0,860.D0,840.D0,820.D0,815.D0,800.D0,
     (2340.D0,2980.D0,2270.D0,1450.D0,1230.D0,1230.D0,
     &1380.D0,1420.D0,1410.D0,1380.D0,1360.D0,
     (1350.D0,1320.D0,2680.D0,3220.D0,2530.D0,1630.D0,1420.D0,1450.D0,
     &1570.D0,1600.D0,1590.D0,1575.D0,1560.D0,1550.D0,1540.D0/
      DATA SEG/24.D0,128.D0,249.D0,256.D0,202.D0,124.D0,73.D0,60.D0,
     &64.D0,69.D0,62.D0,50.D0,44.D0,42.D0,42.D0,41.D0,21.D0,156.D0,
     (273.D0,280.D0,220.D0,212.D0,94.D0,80.D0,82.D0,85.D0,80.D0,73.D0,
     (69.D0,67.D0,66.D0,64.D0,56.D0,296.D0,560.D0,574.D0,467.D0,350.D0,
     &235.D0,210.D0,210.D0,200.D0,190.D0,183.D0,176.D0,170.D0,165.D0,
     (155.D0,100.D0,500.D0,895.D0,880.D0,690.D0,520.D0,378.D0,
     (355.D0,384.D0,373.D0,352.D0,320.D0,300.D0,288.D0,280.D0,
     &262.D0,75.D0,500.D0,965.D0,990.D0,775.D0,525.D0,410.D0,410.D0,
     (433.D0,440.D0,425.D0,395.D0,374.D0,355.D0,340.D0,303.D0,125.D0,
     (570.D0,1025.D0,1100.D0,825.D0,575.D0,418.D0,458.D0,500.D0,
     &480.D0,460.D0,440.D0,422.D0,400.D0,384.D0,355.D0,300.D0,880.D0,
     (1480.D0,1550.D0,1380.D0,940.D0,710.D0,720.D0,810.D0,760.D0,
     (740.D0,700.D0,665.D0,645.D0,620.D0,570.D0,550.D0,1475.D0,
     &2250.D0,2350.D0,1850.D0,1500.D0,
     (1120.D0,1210.D0,1480.D0,1440.D0,1400.D0,1320.D0,
     &1250.D0,1210.D0,1170.D0,1065.D0,540.D0,
     (1300.D0,2220.D0,2560.D0,1980.D0,1650.D0,1160.D0,
     &1360.D0,1600.D0,1560.D0,1510.D0,1410.D0,
     (1350.D0,1300.D0,1270.D0,1200.D0/
C     DATA ITT/1,7,0,0,0,0,0,2,8,0,0,9,3,4,6,5,1,2,9,1,1,1,3,9,10,
      DATA ITT/1,7,0,0,0,0,0,2,8,0,0,9,3,4,6,5,2,8,9,1,1,2,3,9,10,
     &         3,0,0,0,0,7,2,7,2,8,1,7,1,7/
      DATA PLAB/.3D0,.4D0,.5D0,.6D0,.7D0,.8D0,.9D0,1.D0,1.1D0,
     &1.2D0,1.3D0,1.4D0,1.5D0,2.D0,3.D0,4.D0,
     *5.D0,6.D0,10.D0/
      DATA SITO/66.8D0,63.6D0,40.35D0,31.25D0,31.1D0,
     *35.1D0,36.7D0,44.15D0,38.3D0,33.25D0,
     *29.75D0,29.3D0,29.95D0,26.55D0,24.6D0,22.95D0,
     *22.75D0,22.95D0,21.55D0,
     *12.5D0,14.1D0,13.5D0,12.75D0,12.85D0,13.9D0,15.6D0,
     *17.25D0,18.9D0,19.5D0,18.95D0,18.85D0,
     *18.45D0,18.2D0,17.5D0,17.7D0,17.5D0,17.25D0,17.4D0,
     *39.65D0,38.75D0,26.9D0,22.D0,22.D0,24.5D0,26.15D0,30.7D0,28.6D0,
     &26.4D0,24.35D0,24.1D0,24.2D0
     (,22.4D0,21.05D0,20.3D0,20.1D0,20.1D0,19.5D0,
     (280.D0,199.7D0,171.1D0,154.3D0,140.D0,130.D0,116.8D0,117.4D0,
     &111.6D0,109.D0,106.5D0,
     (102.8D0,100.D0,90.2D0,76.7D0,68.D0,62.8D0,60.7D0,56.D0/
      DATA ALP/0.823D0,0.843D0,0.630D0/
      DATA BET/1.26D0,1.31D0,0.90D0/
      DATA STOT /15.D0,20.D0,30.D0,40.D0,60.D0,80.D0,
     &100.D0,150.D0,200.D0/
      DATA REA / .20D0,.23D0,.27D0,.30D0,.35D0,.40D0,.47D0,.55D0,.60D0,
     2           .22D0,.26D0,.31D0,.35D0,.40D0,.45D0,.51D0,.59D0,.63D0,
     3           .24D0,.29D0,.36D0,.42D0,.50D0,.56D0,.60D0,.66D0,.68D0,
     4           .26D0,.32D0,.42D0,.49D0,.58D0,.63D0,.66D0,.71D0,.72D0,
     5           .27D0,.33D0,.44D0,.51D0,.61D0,.65D0,.68D0,.72D0,.74D0,
     6           .28D0,.35D0,.46D0,.53D0,.63D0,.66D0,.69D0,.73D0,.745D0,
     7           .35D0,.42D0,.53D0,.62D0,.69D0,.72D0,.74D0,.77D0,.78D0,
     8           .42D0,.51D0,.62D0,.69D0,.75D0,.77D0,.79D0,.81D0,.82D0,
     9           .44D0,.53D0,.64D0,.70D0,.76D0,.78D0,.80D0,.81D0,.82D0 /
C
C
C
C-------------------------------------------------------------
      SEL=1.0D-20
      ZL=1.0D+20
      IF(AA.LT.0.99)RETURN
      IPOL=0
      PO=POO
      IIT=ITT(IT)
      IF(IT.GT.25)IIT=0
      IF(IIT.EQ.0)RETURN
C---------------------------------------------------------
C**                          ELASTIC SCATTERING ON PROTONS
C         HJM 10/88             REASONABLE FOR P, N, PI+/-
C**
      IPR=IT
      IF(IT.EQ.23) IPR=13
      IF((AA.LT.1.5).AND. (IPR.EQ.1.OR.IPR.EQ.8.OR.IPR.EQ.13.OR.IPR.EQ.
     +14)) THEN
      EKE=SQRT(PO**2+AAM(IPR)**2) - AAM(IPR)
        CALL DSIHAE(IPR,EKE,PO,AA,SEL)
                                                                GOTO 220
      ENDIF
C**
C                            NEUTRON-NUCLEUS ELASTIC SCATTERING
C                            DATA FROM HETKFA2 FOR  EKIN .GT. 15 MEV
C                            FOR PLOTS SEE
C                               P. CLOTH ET AL.,
C                               HERMES - A MC PROGRAM SYSTEM ...
C                               JUEL-2203 (MAY 1988)
C
      IF(IT.EQ.8.AND.PO.LT.20.0D0) THEN
        IF(PO.GT.10) THEN
          IPOL=1
          PO=10.0
        ENDIF
      EKE=SQRT(PO**2+AAM(IT)**2) - AAM(IT)
        CALL DSIHAE(IT,EKE,PO,AA,SEL)
        IF(IPOL.EQ.1)                                           GOTO 240
                                                                GOTO 220
      ENDIF
C-----------------------------------------------------------
C
C
C********************************************************************
C     CALCULATE THE NEW PARTICLE NUMBER IIT:   1=P,2=N,3=PI+,4=PI-,
C     5=K-,6=K+,7=P BAR,8=N BAR,9=K ZERO ,10=K ZERO BAR
C********************************************************************
C
      IF((IIT.EQ.7).OR.(IIT.EQ.8))                              GOTO 250
      IF (PO.GT.20.D0)                                         GOTO 250
      IF(PO.LE.10.D0)                                          GOTO 10
      PO=10.
      IPOL=1
   10 CONTINUE
      IF(IIT.LE.4)                                              GO TO 40
C
C********************************************************************
C     MOMENTUM INDEX K FOR KAONS ANTI KAONS AND ANTI NUCLEONS
C********************************************************************
C
      DO 20 K=1,19
        IF(PO.LE.PLAB(K))                                       GO TO 30
   20 CONTINUE
      K=19
   30                                                           GO TO 90
C
C********************************************************************
C     CALCULATE THE MOMENTUM INDEX K FOR NUCLEONS AND PIONS
C     CALCULATE THE MASS INDEX J OF THE NUCLEUS
C********************************************************************
C
   40 CONTINUE
      DO 50 K=1,16
        IF(PO.LE.P(K))                                          GO TO 60
   50 CONTINUE
      K=16
   60 CONTINUE
      DO 80 I=2,8
        IF(AA.LE.A(I))                                          GO TO 70
                                                                GO TO 80
   70   CONTINUE
        J=I-1
                                                                GO TO 90
   80 CONTINUE
      J=8
C
C********************************************************************
C     SELECT THE FORMULEI TO BE USED FOR DIFFERENT PARTICLE TYPES
C********************************************************************
C
   90 CONTINUE
C             P , N ,PI+,PI-,K- ,K+ ,AP ,AN ,K0 ,AK0
      GO TO (100,100,110,110,120,130,140,150,130,120),IIT
C******************** PROTONS,NEUTRONS,OTHERS
  100 K=K-3
      IF(K.LT.1) K=1
      ALOGA=LOG(A(J+1)/A(J))
      AAA=AA/A(J)
      SI1=SIG(K,J)* AAA     **(LOG(SIG(K,J+1)/SIG(K,J))/ALOGA)
      IF(K.EQ.1)                                               GO TO 230
      KK=K-1
      SI2=SIG(KK,J)* AAA     **(LOG(SIG(KK,J+1)/SIG(KK,J))/ALOGA)
      K=K+3
      KK=KK+3
      SI=SI1+(PO-P(K))*(SI2-SI1)/(P(KK)-P(K))
C
      SEL=SI
C
                                                               GO TO 210
C******************** CHARGED PIONS
  110 CONTINUE
      ALOGA=LOG(A(J+1)/A(J))
      AAA=AA/A(J)
      SI1=SEG(K,J)* AAA     **(LOG(SEG(K,J+1)/SEG(K,J))/ALOGA)
      IF(K.EQ.1)                                               GO TO 230
      KK=K-1
      SI2=SEG(KK,J)* AAA     **(LOG(SEG(KK,J+1)/SEG(KK,J))/ALOGA)
      SI=SI1+(PO-P(K))*(SI2-SI1)/(P(KK)-P(K))
C
      SEL=SI
C
                                                               GO TO 210
C******************** K-,K ZERO BAR
  120 CONTINUE
      IA=1
      IS=1
                                                               GO TO 160
C******************** K+, K ZERO
  130 CONTINUE
      IA=2
      IS=2
                                                               GO TO 160
C******************** P BAR
  140 CONTINUE
C******************** N BAR
  150 CONTINUE
C
      PO=POO
                                                                GOTO 250
C
C
C********************************************************************
C     KAONS, ANTI KAONS
C********************************************************************
C
  160 KK=K-1
      IF(K.EQ.1)                                               GO TO 170
      PKK=PLAB(KK)
      SIKK=SITO(KK,IS)
      SI=(SITO(K,IS)-SIKK)*(PO-PKK)/(PLAB(K)-PKK)+SIKK
                                                               GO TO 180
  170 SI=SITO(K,IS)
  180 SI1=SI
      SI=BET(IA)*SI*AA**ALP(IA)
      IV=IT
      IF(IV.NE.24)                                             GO TO 190
      IV=15
      SI=SI*2.06
                                                               GO TO 200
  190 IF(IV.EQ.25) IV=16
C*** 151  CALL NIZL(IV,AA,PO,SINEL,ZLIN)      *** 30/08/90
  200 CALL SIHNIN(IV,1,PO,SINEL)
C
      SEL=SI-SINEL
      IF(IPOL.EQ.1)                                             GOTO 240
                                                                GOTO 210
C
C********************************************************************
C     AND NOW THE SCATTERING LENGTH IN G/CM**2
C********************************************************************
C
  210 CONTINUE
      IF(IPOL.EQ.1)                                             GOTO 240
  220 CONTINUE
C
      IF(SEL.LT.1.D-15) SEL=1.D-15
      ZL=10000.D0*AA/(6.022*SEL)
      RETURN
C
C********************************************************************
C     WE ARE IN THE LOWEST MOMENTUM BIN
C********************************************************************
C
  230 SI=SI1
      SEL=SI
                                                               GO TO 210
C***
C   ENTRY FOR SMOOTHING OF SIGEL BETWEEN 10. AND 20. GEV/C
C***
  240 CONTINUE
      PO=20.
C***
C   APPROXIMATION FOR HIGH ENERGIES
C***
  250 CONTINUE
C
      IT1=IT
      IF((IT.EQ.2).OR.(IT.EQ.9)) IT1=1
C
      STO=DSHPTO(IT1,PO)
C
C   MASS NUMBER INDEX
C***
      DO 260 IA=2,8
        IF(AA.GT.A(IA))                                         GOTO 260
        JA=IA-1
                                                                GOTO 270
  260 CONTINUE
      JA=8
  270 CONTINUE
C***
C   SIGTOT INDEX
C***
      DO 280 IS=2,8
        IF(STO.GT.STOT(IS))                                     GOTO 280
        JS=IS-1
                                                                GOTO 290
  280 CONTINUE
      JS=8
  290 CONTINUE
C
      DA1=A(JA+1)-A(JA)
      DA2=AA - A(JA)
      RR=REA(JS,JA)
      R1=RR + DA2*(REA(JS,JA+1)-RR)/DA1
      RR=REA(JS+1,JA)
      R2=RR + DA2*(REA(JS+1,JA+1)-RR)/DA1
      RACT=R1 + (STO-STOT(JS))*(R2-R1)/(STOT(JS+1)-STOT(JS))
C
C***     CALL NIZL(IT,AA,PO,SINEL1,ZLIN)          ***  30/08/90
      CALL SIHNIN(IT,1,PO,SINEL1)
      SEL1=RACT*SINEL1
      IF(IPOL.EQ.1)                                             GOTO 300
      SEL=SEL1
                                                                GOTO 220
C
  300 CONTINUE
      SEL=SEL + (SEL1-SEL)*(POO-10.)/10.
                                                                GOTO 220
C
C
C********************************************************************
C     FORMATS
C********************************************************************
C
 1000 FORMAT('          WARNING AT CALL SIGEL  ',I5)
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE DSIHAE(KPROJ,EKIN,PLAB,ANUC,SIGELA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C***
C        HJM 22/10/88
C        HJM 29/08/90 :
C                  correction of pi- proton data (Plab < 1.8 GeV/c)
C
C        CROSS SECTIONS FOR ELASTIC SCATTERING
C
C        INCLUDING - PION/NUCLEON PROTON DATA FROM BERTINI (HETKFA2)
C
C                  - ...  HIGH-ENERGY APPROXIMATION:
C                                       SIGEL/SIGTOT = CONST
C
C                  - NUCLEON-NUCLEUS DATA FROM HETKFA2
C***
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
      PARAMETER (NEN=106)
      PARAMETER (NEA=23)
      PARAMETER (NNAA=10)
      DIMENSION EKIHN(NEN),EKIHA(NEA),AMASS(NNAA)
      DIMENSION SEPIMP(NEN),SEPIPP(NEN),SEPP(NEN),SENP(NEN)
      DIMENSION SENA(NEA,NNAA),SEPA(NEA,NNAA)
      DIMENSION TSIG(2)
      DIMENSION RELTO(14)
C*************************************************************8*
*                                                                     *
*=== dblprc ==========================================================*
*                                                                     *
*---------------------------------------------------------------------*
*                                                                     *
*      Dblprc: included in any routine                                *
*                                                                     *
*  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  *
*  !!!! O N   M A C H I N E S   W H E R E   T H E   D O U B L E !!!!  *
*  !!!! P R E C I S I O N   I S   N O T   R E Q U I R E D  R E -!!!!  *
*  !!!! M O V E   T H E   D O U B L E   P R E C I S I O N       !!!!  *
*  !!!! S T A T E M E N T,  S E T   K A L G N M = 1   A N D     !!!!  *
*  !!!! C H A N G E   A L L   N U M E R I C A L   C O N S -     !!!!  *
*  !!!! T A N T S   T O   S I N G L E   P R E C I S I O N       !!!!  *
*  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  *
*                                                                     *
*         Kalgnm = real address alignment, 2 for double precision,    *
*                  1 for single precision                             *
*         Anglgb = this parameter should be set equal to the machine  *
*                  "zero" with respect to unit                        *
*         Anglsq = this parameter should be set equal to the square   *
*                  of Anglgb                                          *
*         Axcssv = this parameter should be set equal to the number   *
*                  for which unity is negligible for the machine      *
*                  accuracy                                           *
*         Andrfl = "underflow" of the machine for floating point      *
*                  operation                                          *
*         Avrflw = "overflow"  of the machine for floating point      *
*                  operation                                          *
*         Ainfnt = code "infinite"                                    *
*         Azrzrz = code "zero"                                        *
*         Einfnt = natural logarithm of the code "infinite"           *
*         Ezrzrz = natural logarithm of the code "zero"               *
*         Onemns = 1- of the machine, it is 1 - 2 x Anglgb            *
*         Onepls = 1+ of the machine, it is 1 + 2 x Anglgb            *
*         Csnnrm = maximum tolerable error on cosine normalization,   *
*                  u**2+v**2+w**2: assuming a typical anglgb relative *
*                  error on each component we would get 2xanglgb: use *
*                  4xanglgb to avoid too many normalizations          *
*         Dmxtrn = "infinite" distance for transport (cm)             *
*                                                                     *
*---------------------------------------------------------------------*
*                                                                     *
      PARAMETER ( KALGNM = 2 )
      PARAMETER ( ANGLGB = 5.0D-16 )
      PARAMETER ( ANGLSQ = 2.5D-31 )
      PARAMETER ( AXCSSV = 0.2D+16 )
      PARAMETER ( ANDRFL = 1.0D-38 )
      PARAMETER ( AVRFLW = 1.0D+38 )
      PARAMETER ( AINFNT = 1.0D+30 )
      PARAMETER ( AZRZRZ = 1.0D-30 )
      PARAMETER ( EINFNT = +69.07755278982137 D+00 )
      PARAMETER ( EZRZRZ = -69.07755278982137 D+00 )
      PARAMETER ( ONEMNS = 0.999999999999999  D+00 )
      PARAMETER ( ONEPLS = 1.000000000000001  D+00 )
      PARAMETER ( CSNNRM = 2.0D-15 )
      PARAMETER ( DMXTRN = 1.0D+08 )
*
*======================================================================*
*======================================================================*
*=========                                                   ==========*
*=========    M A T H E M A T I C A L   C O N S T A N T S    ==========*
*=========                                                   ==========*
*======================================================================*
*======================================================================*
*                                                                      *
*   Numerical constants:                                               *
*                                                                      *
*         Zerzer = 0                                                   *
*         Oneone = 1                                                   *
*         Twotwo = 2                                                   *
*         Thrthr = 3                                                   *
*         Foufou = 4                                                   *
*         Fivfiv = 5                                                   *
*         Sixsix = 6                                                   *
*         Sevsev = 7                                                   *
*         Eigeig = 8                                                   *
*         Aninen = 9                                                   *
*         Tenten = 10                                                  *
*         Hlfhlf = 1/2                                                 *
*         Onethi = 1/3                                                 *
*         Twothi = 2/3                                                 *
*         Pipipi = Circumference / diameter                            *
*         Eneper = "e", base of natural logarithm                      *
*         Sqrent = square root of "e"                                  *
*                                                                      *
*----------------------------------------------------------------------*
*
      PARAMETER ( ZERZER = 0.D+00 )
      PARAMETER ( ONEONE = 1.D+00 )
      PARAMETER ( TWOTWO = 2.D+00 )
      PARAMETER ( THRTHR = 3.D+00 )
      PARAMETER ( FOUFOU = 4.D+00 )
      PARAMETER ( FIVFIV = 5.D+00 )
      PARAMETER ( SIXSIX = 6.D+00 )
      PARAMETER ( SEVSEV = 7.D+00 )
      PARAMETER ( EIGEIG = 8.D+00 )
      PARAMETER ( ANINEN = 9.D+00 )
      PARAMETER ( TENTEN = 10.D+00 )
      PARAMETER ( HLFHLF = 0.5D+00 )
      PARAMETER ( ONETHI = ONEONE / THRTHR )
      PARAMETER ( TWOTHI = TWOTWO / THRTHR )
      PARAMETER ( PIPIPI = 3.1415926535897932270 D+00 )
      PARAMETER ( ENEPER = 2.7182818284590452354 D+00 )
      PARAMETER ( SQRENT = 1.6487212707001281468 D+00 )
*
*======================================================================*
*======================================================================*
*=========                                                   ==========*
*=========       P H Y S I C A L   C O N S T A N T S         ==========*
*=========                                                   ==========*
*======================================================================*
*======================================================================*
*                                                                      *
*   Primary constants:                                                 *
*                                                                      *
*         Clight = speed of light in cm s-1                            *
*         Avogad = Avogadro number                                     *
*         Amelgr = electron mass (g)                                   *
*         Plckbr = reduced Planck constant (erg s)                     *
*         Elccgs = elementary charge (CGS unit)                        *
*         Elcmks = elementary charge (MKS unit)                        *
*         Amugrm = Atomic mass unit (g)                                *
*         Ammumu = Muon mass (amu)                                     *
*                                                                      *
*   Derived constants:                                                 *
*                                                                      *
*         Alpfsc = Fine structure constant  = e^2/(hbar c)             *
*         Amelct = Electron mass (GeV) = 10^-16Amelgr Clight^2 / Elcmks*
*         Amugev = Atomic mass unit (GeV) = 10^-16Amelgr Clight^2      *
*                                           / Elcmks                   *
*         Ammuon = Muon mass (GeV) = Ammumu * Amugev                   *
*         Fscto2 = (Fine structure constant)^2                         *
*         Fscto3 = (Fine structure constant)^3                         *
*         Fscto4 = (Fine structure constant)^4                         *
*         Plabrc = Reduced Planck constant times the light velocity    *
*                  expressed in GeV fm                                 *
*         Rclsel = Classical electron radius (cm) = e^2 / (m_e c^2)    *
*   Conversion constants:                                              *
*         GeVMeV = from GeV to MeV                                     *
*         eMVGeV = from MeV to GeV                                     *
*         Raddeg = from radians to degrees                             *
*         Degrad = from degrees to radians                             *
*                                                                      *
*----------------------------------------------------------------------*
*
      PARAMETER ( CLIGHT = 2.99792458         D+10 )
      PARAMETER ( AVOGAD = 6.0221367          D+23 )
      PARAMETER ( AMELGR = 9.1093897          D-28 )
      PARAMETER ( PLCKBR = 1.05457266         D-27 )
      PARAMETER ( ELCCGS = 4.8032068          D-10 )
      PARAMETER ( ELCMKS = 1.60217733         D-19 )
      PARAMETER ( AMUGRM = 1.6605402          D-24 )
      PARAMETER ( AMMUMU = 0.113428913        D+00 )
*     PARAMETER ( ALPFSC = 1.D+00 / 137.035989561D+00 )
*     PARAMETER ( FSCTO2 = ALPFSC * ALPFSC )
*     PARAMETER ( FSCTO3 = FSCTO2 * ALPFSC )
*     PARAMETER ( FSCTO4 = FSCTO3 * ALPFSC )
*    It is important to set the electron mass exactly with the same
*    rounding as in the mass tables, so use the explicit expression
*     PARAMETER ( AMELCT = 1.D-16 * AMELGR * CLIGHT * CLIGHT / ELCMKS )
*    It is important to set the amu mass exactly with the same
*    rounding as in the mass tables, so use the explicit expression
*     PARAMETER ( AMUGEV = 1.D-16 * AMUGRM * CLIGHT * CLIGHT / ELCMKS )
*    It is important to set the muon mass exactly with the same
*    rounding as in the mass tables, so use the explicit expression
*     PARAMETER ( AMMUON = AMMUMU * AMUGEV ELCMKS )
*     PARAMETER ( RCLSEL = ELCCGS * ELCCGS / CLIGHT / CLIGHT / AMELGR )
      PARAMETER ( ALPFSC = 7.2973530791728595 D-03 )
      PARAMETER ( FSCTO2 = 5.3251361962113614 D-05 )
      PARAMETER ( FSCTO3 = 3.8859399018437826 D-07 )
      PARAMETER ( FSCTO4 = 2.8357075508200407 D-09 )
      PARAMETER ( PLABRC = 0.197327053        D+00 )
      PARAMETER ( AMELCT = 0.51099906         D-03 )
      PARAMETER ( AMUGEV = 0.93149432         D+00 )
      PARAMETER ( AMMUON = 0.105658389        D+00 )
      PARAMETER ( RCLSEL = 2.8179409183694872 D-13 )
      PARAMETER ( GEVMEV = 1.0                D+03 )
      PARAMETER ( EMVGEV = 1.0                D-03 )
      PARAMETER ( RADDEG = 180.D+00 / PIPIPI )
      PARAMETER ( DEGRAD = PIPIPI / 180.D+00 )
  
*$ CREATE IOUNIT.ADD
*COPY IOUNIT
*                                                                     *
*=== iounit ==========================================================*
*                                                                     *
*---------------------------------------------------------------------*
*                                                                     *
*      Iounit: included in any routine                                *
*                                                                     *
*         lunin  = standard input unit                                *
*         lunout = standard output unit                               *
*         lunerr = standard error unit                                *
*         lunber = input file for bertini nuclear data                *
*         lunech = echo file for pegs dat                             *
*         lunflu = input file for photoelectric edges and X-ray fluo- *
*                  rescence data                                      *
*         lungeo = scratch file for combinatorial geometry            *
*         lunpgs = input file for pegs material data                  *
*         lunran = output file for the final random number seed       *
*         lunxsc = input file for low energy neutron cross sections   *
*         lunrdb = unit number for reading (extra) auxiliary external *
*                  files to be closed just after reading              *
*                                                                     *
*---------------------------------------------------------------------*
*                                                                     *
      PARAMETER ( LUNIN  = 5  )
      PARAMETER ( LUNOUT = 6  )
      PARAMETER ( LUNERR = 66 )
      PARAMETER ( LUNBER = 14 )
      PARAMETER ( LUNECH = 8  )
      PARAMETER ( LUNFLU = 86 )
      PARAMETER ( LUNGEO = 16 )
      PARAMETER ( LUNPGS = 12 )
      PARAMETER ( LUNRAN = 2  )
      PARAMETER ( LUNXSC = 81 )
      PARAMETER ( LUNRDB = 1  )

*$ CREATE DIMPAR.ADD
*COPY DIMPAR
*                                                                     *
*=== dimpar ==========================================================*
*                                                                     *
*---------------------------------------------------------------------*
*                                                                     *
*      DIMPAR: included in any routine                                *
*                                                                     *
*          Mxxrgn = maximum number of regions                         *
*          Mxxmdf = maximum number of media in Fluka                  *
*          Mxxmde = maximum number of media in Emf                    *
*          Mfstck = stack dimension in Fluka                          *
*          Mestck = stack dimension in Emf                            *
*          Nallwp = number of allowed particles                       *
*          Mpdpdx = number of particle types for which EM dE/dx pro-  *
*                   cesses (ion,pair,bremss) have to be computed      *
*          Icomax = maximum number of materials for compounds (equal  *
*                   to the sum of the number of materials for every   *
*                   compound )                                        *
*          Nstbis = number of stable isotopes recorded in common iso- *
*                   top                                               *
*          Idmaxp = number of particles/resonances defined in common  *
*                   part                                              *
*                                                                     *
*---------------------------------------------------------------------*
*                                                                     *
      PARAMETER ( MXXRGN = 500  )
      PARAMETER ( MXXMDF = 56   )
      PARAMETER ( MXXMDE = 50   )
      PARAMETER ( MFSTCK = 1000 )
      PARAMETER ( MESTCK = 100  )
      PARAMETER ( NALLWP = 39   )
      PARAMETER ( MPDPDX = 8    )
      PARAMETER ( ICOMAX = 180  )
      PARAMETER ( NSTBIS = 304  )
      PARAMETER ( IDMAXP = 210  )


C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
      SAVE EKIHN, EKIHA, AMASS, SEPIMP, SEPIPP, SEPP, SENP, SENA, SEPA,
     &     TSIG , RELTO
* Statement functions: A.Ferrari 28-4-93
      SIGMLW (E) = 3.D+03 * PIPIPI
     &           / ( 1.206D+03 * E + ( -1.86D+00
     &           + 0.09415D+03 * E + 0.0001306D+06 * E**2
     &           )**2 ) + 1.D+03 * PIPIPI / ( 1.206D+03 * E
     &           + ( 0.4223D+00 + 0.13D+03 * E )**2 )
      SIGMPN (BETAPR) = 34.10D+00 / BETAPR**2 - 82.2D+00 / BETAPR
     &                + 82.2D+00
      SIGMPP (BETAPR) = 10.63D+00 / BETAPR**2 - 29.92D+00 / BETAPR
     &                + 42.9D+00
C***
C   KINETIC ENERGIES FOR TABLE LOOK-UP

      DATA EKIHN /
     &   0.00D0, 0.02D0, 0.04D0, 0.06D0, 0.08D0, 0.10D0, 0.12D0, 0.14D0,
     &   0.16D0, 0.18D0, 0.20D0, 0.22D0, 0.24D0, 0.26D0, 0.28D0, 0.30D0,
     &   0.32D0, 0.34D0, 0.36D0, 0.38D0, 0.40D0, 0.42D0, 0.44D0, 0.46D0,
     &   0.48D0, 0.50D0, 0.52D0, 0.54D0, 0.56D0, 0.58D0, 0.60D0, 0.62D0,
     &   0.64D0, 0.66D0, 0.68D0, 0.70D0, 0.72D0, 0.74D0, 0.76D0, 0.78D0,
     &   0.80D0, 0.82D0, 0.84D0, 0.86D0, 0.88D0, 0.90D0, 0.92D0, 0.94D0,
     &   0.96D0, 0.98D0, 1.00D0, 1.02D0, 1.04D0, 1.06D0, 1.08D0, 1.10D0,
     &   1.12D0, 1.14D0, 1.16D0, 1.18D0, 1.20D0, 1.22D0, 1.24D0, 1.26D0,
     &   1.28D0, 1.30D0, 1.32D0, 1.34D0, 1.36D0, 1.38D0, 1.40D0, 1.42D0,
     &   1.44D0, 1.46D0, 1.48D0, 1.50D0, 1.52D0, 1.54D0, 1.56D0, 1.58D0,
     &   1.60D0, 1.62D0, 1.64D0, 1.66D0, 1.68D0, 1.70D0, 1.72D0, 1.74D0,
     &   1.76D0, 1.78D0, 1.80D0, 1.82D0, 1.84D0, 1.86D0, 1.88D0, 1.90D0,
     &   1.92D0, 1.94D0, 1.96D0, 1.98D0, 2.00D0, 2.5D0,  3.0D0,  3.5D0,
     &   5.0D0, 10.0D0/
      DATA EKIHA /
     &   0.015D0, 0.02D0, 0.025D0, 0.03D0,  0.04D0, 0.05D0, 0.06D0,
     &   0.08D0,  0.10D0, 0.125D0, 0.15D0, 0.175D0, 0.20D0, 0.225D0,
     &   0.25D0,  0.3D0,  0.4D0,   0.6D0,  1.0D0,   2.0D0,  5.0D0,
     &   10.0D0,  22.5D0/
      DATA AMASS /
     &   4.D0, 9.D0, 12.D0, 27.D0, 47.9D0, 55.9D0, 63.5D0, 112.4D0,
     &   207.2D0, 238.1D0/
C-------------------------------------------------------------------
C
C***     PI(-)-P ELASTIC CROSS SECTION DATA
      DATA (SEPIMP(IE),IE=1,50) /
     *     1.250D+00,  1.500D+00,  1.750D+00,  2.450D+00,  3.800D+00,
     *     6.000D+00,  9.700D+00,  1.500D+01,  2.140D+01,  2.310D+01,
     *     2.295D+01,  2.070D+01,  1.795D+01,  1.550D+01,  1.360D+01,
     *     1.230D+01,  1.130D+01,  1.070D+01,  1.050D+01,  1.070D+01,
     *     1.120D+01,  1.175D+01,  1.235D+01,  1.300D+01,  1.400D+01,
     *     1.500D+01,  1.600D+01,  1.700D+01,  1.835D+01,  1.970D+01,
     *     2.050D+01,  1.915D+01,  1.770D+01,  1.650D+01,  1.570D+01,
     *     1.520D+01,  1.510D+01,  1.525D+01,  1.550D+01,  1.600D+01,
     *     1.685D+01,  1.800D+01,  2.000D+01,  2.230D+01,  2.475D+01,
     *     2.635D+01,  2.510D+01,  2.300D+01,  2.140D+01,  2.000D+01/
      DATA (SEPIMP(IE),IE=51,106) /
     *     1.870D+01,  1.750D+01,  1.670D+01,  1.585D+01,  1.505D+01,
     *     1.440D+01,  1.395D+01,  1.340D+01,  1.299D+01,  1.260D+01,
     *     1.215D+01,  1.175D+01,  1.140D+01,  1.099D+01,  1.060D+01,
     *     1.040D+01,  1.010D+01,  9.990D+00,  9.900D+00,  9.750D+00,
     *     9.600D+00,  9.550D+00,  9.450D+00,  9.350D+00,  9.250D+00,
     *     9.250D+00,  9.350D+00,  9.650D+00,  9.850D+00,  1.000D+01,
     *     1.015D+01,  1.030D+01,  1.060D+01,  1.080D+01,  1.095D+01,
     *     1.100D+01,  1.095D+01,  1.090D+01,  1.070D+01,  1.035D+01,
     *     1.000D+01,  9.600D+00,  9.050D+00,  8.550D+00,  8.200D+00,
     *     8.000D+00,  7.850D+00,  7.800D+00,  7.750D+00,  7.700D+00,
     *     7.650D+00,
     *     7.600D+00,  7.240D+00,  6.770D+00,  5.840D+00,  4.570D+00/
* *** The previous 5 points have been substituted to the erroneous
* *** ones from H.J. Mohring by A. Ferrari
C---------------------------------------------------------------------
C
C***     PI(+)-P ELASTIC CROSS SECTION DATA
      DATA (SEPIPP(IE),IE=1,50) /
     *     1.800D+00,  4.000D+00,  9.900D+00,  2.170D+01,  4.000D+01,
     *     6.580D+01,  9.680D+01,  1.392D+02,  1.800D+02,  2.000D+02,
     *     1.655D+02,  1.420D+02,  1.225D+02,  1.032D+02,  8.400D+01,
     *     6.725D+01,  5.510D+01,  4.725D+01,  4.130D+01,  3.690D+01,
     *     3.230D+01,  2.885D+01,  2.600D+01,  2.300D+01,  2.090D+01,
     *     1.875D+01,  1.675D+01,  1.500D+01,  1.340D+01,  1.200D+01,
     *     1.100D+01,  9.980D+00,  9.200D+00,  8.600D+00,  8.200D+00,
     *     8.100D+00,  8.100D+00,  8.250D+00,  8.500D+00,  8.750D+00,
     *     9.000D+00,  9.400D+00,  9.750D+00,  1.000D+01,  1.030D+01,
     *     1.075D+01,  1.130D+01,  1.200D+01,  1.275D+01,  1.330D+01/
      DATA (SEPIPP(IE),IE=51,106) /
     *     1.350D+01,  1.335D+01,  1.330D+01,  1.330D+01,  1.345D+01,
     *     1.355D+01,  1.380D+01,  1.400D+01,  1.460D+01,  1.500D+01,
     *     1.555D+01,  1.625D+01,  1.700D+01,  1.800D+01,  1.875D+01,
     *     1.920D+01,  1.925D+01,  1.890D+01,  1.830D+01,  1.790D+01,
     *     1.725D+01,  1.690D+01,  1.640D+01,  1.600D+01,  1.550D+01,
     *     1.505D+01,  1.475D+01,  1.430D+01,  1.400D+01,  1.365D+01,
     *     1.335D+01,  1.300D+01,  1.280D+01,  1.250D+01,  1.225D+01,
     *     1.205D+01,  1.195D+01,  1.175D+01,  1.150D+01,  1.135D+01,
     *     1.105D+01,  1.095D+01,  1.080D+01,  1.060D+01,  1.030D+01,
     *     1.020D+01,  1.005D+01,  9.900D+00,  9.800D+00,  9.700D+00,
     *     9.600D+00,
     *     7.350D+00,  7.200D+00,  7.000D+00,  5.800D+00,  4.800D+00/
C---------------------------------------------------------------------
C
C***     P-P ELASTIC CROSS SECTION DATA
      DATA (SEPP(IE),IE=1,50) /
     *     6.750D+02,  1.550D+02,  6.750D+01,  4.420D+01,  3.230D+01,
     *     2.800D+01,  2.520D+01,  2.370D+01,  2.300D+01,  2.275D+01,
     *     2.260D+01,  2.260D+01,  2.260D+01,  2.260D+01,  2.270D+01,
     *     2.280D+01,  2.295D+01,  2.300D+01,  2.310D+01,  2.330D+01,
     *     2.350D+01,  2.380D+01,  2.395D+01,  2.420D+01,  2.460D+01,
     *     2.485D+01,  2.500D+01,  2.530D+01,  2.565D+01,  2.600D+01,
     *     2.620D+01,  2.640D+01,  2.660D+01,  2.675D+01,  2.690D+01,
     *     2.700D+01,  2.705D+01,  2.710D+01,  2.715D+01,  2.720D+01,
     *     2.725D+01,  2.725D+01,  2.720D+01,  2.715D+01,  2.710D+01,
     *     2.700D+01,  2.695D+01,  2.680D+01,  2.670D+01,  2.660D+01/
      DATA (SEPP(IE),IE=51,106) /
     *     2.640D+01,  2.625D+01,  2.605D+01,  2.590D+01,  2.570D+01,
     *     2.545D+01,  2.525D+01,  2.500D+01,  2.480D+01,  2.470D+01,
     *     2.450D+01,  2.430D+01,  2.410D+01,  2.395D+01,  2.370D+01,
     *     2.360D+01,  2.340D+01,  2.325D+01,  2.305D+01,  2.290D+01,
     *     2.275D+01,  2.270D+01,  2.260D+01,  2.250D+01,  2.230D+01,
     *     2.225D+01,  2.210D+01,  2.200D+01,  2.195D+01,  2.190D+01,
     *     2.175D+01,  2.165D+01,  2.150D+01,  2.140D+01,  2.125D+01,
     *     2.120D+01,  2.105D+01,  2.100D+01,  2.090D+01,  2.075D+01,
     *     2.065D+01,  2.055D+01,  2.045D+01,  2.030D+01,  2.020D+01,
     *     2.005D+01,  2.000D+01,  1.995D+01,  1.980D+01,  1.975D+01,
     *     1.965D+01,
     *     17.15D+00,  14.45D+00,  13.00D+00,  11.50D+00,  10.50D+00/
C--------------------------------------------------------------------
C
C***     N-P ELASTIC CROSS SECTION DATA
      DATA (SENP(IE),IE=1,50) /
     *     1.965D+03,  4.750D+02,  2.200D+02,  1.300D+02,  9.180D+01,
     *     7.300D+01,  6.030D+01,  5.180D+01,  4.680D+01,  4.320D+01,
     *     4.080D+01,  3.910D+01,  3.760D+01,  3.650D+01,  3.550D+01,
     *     3.480D+01,  3.415D+01,  3.370D+01,  3.325D+01,  3.290D+01,
     *     3.275D+01,  3.250D+01,  3.255D+01,  3.275D+01,  3.285D+01,
     *     3.275D+01,  3.220D+01,  3.150D+01,  3.075D+01,  2.990D+01,
     *     2.875D+01,  2.775D+01,  2.695D+01,  2.630D+01,  2.590D+01,
     *     2.565D+01,  2.560D+01,  2.560D+01,  2.560D+01,  2.565D+01,
     *     2.570D+01,  2.575D+01,  2.578D+01,  2.580D+01,  2.585D+01,
     *     2.580D+01,  2.575D+01,  2.560D+01,  2.540D+01,  2.505D+01/
      DATA (SENP(IE),IE=51,106) /
     *     2.470D+01,  2.425D+01,  2.375D+01,  2.315D+01,  2.275D+01,
     *     2.230D+01,  2.200D+01,  2.175D+01,  2.155D+01,  2.145D+01,
     *     2.130D+01,  2.125D+01,  2.115D+01,  2.105D+01,  2.100D+01,
     *     2.095D+01,  2.090D+01,  2.080D+01,  2.070D+01,  2.060D+01,
     *     2.050D+01,  2.045D+01,  2.040D+01,  2.030D+01,  2.025D+01,
     *     2.020D+01,  2.015D+01,  2.010D+01,  2.005D+01,  2.002D+01,
     *     2.000D+01,  1.999D+01,  1.990D+01,  1.985D+01,  1.975D+01,
     *     1.970D+01,  1.965D+01,  1.960D+01,  1.950D+01,  1.945D+01,
     *     1.940D+01,  1.925D+01,  1.920D+01,  1.915D+01,  1.910D+01,
     *     1.900D+01,  1.898D+01,  1.895D+01,  1.890D+01,  1.880D+01,
     *     1.875D+01,
     *     17.00D+00,  14.40D+00,  12.00D+00,  11.00D+00,  10.00D+00/
C---------------------------------------------------------------------
C
C***     N-A ELASTIC CROSS SECTION DATA
C*                  NEUTRON - HELIUM
      DATA (SENA(IE,1),IE=1,NEA) /
     *     5.103D-01,  5.157D-01,  5.103D-01,  4.777D-01,  4.072D-01,
     *     3.420D-01,  2.714D-01,  1.683D-01,  6.700D-02,  6.100D-02,
     *     5.800D-02,  4.900D-02,  3.800D-02,  3.300D-02,  3.000D-02,
     *     2.400D-02,  2.300D-02,  2.900D-02,  3.600D-02,  4.100D-02,
     *     4.000D-02,  3.700D-02,  3.400D-02/
C
C*                  NEUTRON - BERYLLIUM
      DATA (SENA(IE,2),IE=1,NEA) /
     *     8.762D-01,  8.856D-01,  8.762D-01,  8.203D-01,  6.991D-01,
     *     5.873D-01,  4.661D-01,  2.890D-01,  1.401D-01,  1.305D-01,
     *     1.238D-01,  1.069D-01,  8.495D-02,  7.480D-02,  6.750D-02,
     *     5.565D-02,  5.230D-02,  6.470D-02,  7.765D-02,  8.722D-02,
     *     8.440D-02,  7.821D-02,  7.259D-02/
C
C*                  NEUTRON - CARBON
      DATA (SENA(IE,3),IE=1,NEA) /
     *     9.200D-01,  9.500D-01,  9.400D-01,  8.800D-01,  7.500D-01,
     *     6.100D-01,  5.000D-01,  3.700D-01,  1.820D-01,  1.710D-01,
     *     1.620D-01,  1.410D-01,  1.130D-01,  1.000D-01,  9.000D-02,
     *     7.500D-02,  7.000D-02,  8.600D-02,  1.020D-01,  1.140D-01,
     *     1.100D-01,  1.020D-01,  9.500D-02/
C
C*                  NEUTRON - ALUMINUM
      DATA (SENA(IE,4),IE=1,NEA) /
     *     1.090D+00,  1.180D+00,  1.240D+00,  1.280D+00,  1.260D+00,
     *     1.160D+00,  9.300D-01,  6.300D-01,  3.580D-01,  3.450D-01,
     *     3.350D-01,  2.990D-01,  2.480D-01,  2.220D-01,  2.020D-01,
     *     1.730D-01,  1.610D-01,  1.920D-01,  2.200D-01,  2.420D-01,
     *     2.370D-01,  2.220D-01,  2.060D-01/
C
C*                  NEUTRON - TITANIUM
      DATA (SENA(IE,5),IE=1,NEA) /
     *     1.029D+00,  9.469D-01,  1.091D+00,  1.284D+00,  1.591D+00,
     *     1.691D+00,  1.258D+00,  9.241D-01,  5.620D-01,  5.493D-01,
     *     5.375D-01,  4.907D-01,  4.182D-01,  3.800D-01,  3.484D-01,
     *     3.038D-01,  2.823D-01,  3.307D-01,  3.720D-01,  4.040D-01,
     *     3.959D-01,  3.743D-01,  3.517D-01/
C
C*                  NEUTRON - IRON
      DATA (SENA(IE,6),IE=1,NEA) /
     *     1.178D+00,  9.793D-01,  1.090D+00,  1.271D+00,  1.650D+00,
     *     1.799D+00,  1.339D+00,  1.009D+00,  6.223D-01,  6.132D-01,
     *     6.042D-01,  5.572D-01,  4.812D-01,  4.402D-01,  4.053D-01,
     *     3.554D-01,  3.304D-01,  3.814D-01,  4.244D-01,  4.603D-01,
     *     4.523D-01,  4.293D-01,  4.053D-01/
C
C*                  NEUTRON - COPPER
      DATA (SENA(IE,7),IE=1,NEA) /
     *     1.386D+00,  1.050D+00,  1.134D+00,  1.302D+00,  1.722D+00,
     *     1.922D+00,  1.449D+00,  1.103D+00,  6.762D-01,  6.686D-01,
     *     6.602D-01,  6.131D-01,  5.344D-01,  4.912D-01,  4.541D-01,
     *     4.004D-01,  3.728D-01,  4.273D-01,  4.725D-01,  5.103D-01,
     *     5.022D-01,  4.781D-01,  4.524D-01/
C
C*                  NEUTRON - CADMIUM
      DATA (SENA(IE,8),IE=1,NEA) /
     *     2.029D+00,  1.537D+00,  1.660D+00,  1.906D+00,  2.520D+00,
     *     2.812D+00,  2.121D+00,  1.614D+00,  1.014D+00,  1.012D+00,
     *     1.006D+00,  9.557D-01,  8.607D-01,  8.038D-01,  7.541D-01,
     *     6.775D-01,  6.334D-01,  7.080D-01,  7.669D-01,  8.156D-01,
     *     8.074D-01,  7.769D-01,  7.404D-01/
C
C*                  NEUTRON - LEAD
      DATA (SENA(IE,9),IE=1,NEA) /
     *     3.050D+00,  2.310D+00,  2.495D+00,  2.865D+00,  3.789D+00,
     *     4.228D+00,  3.188D+00,  2.426D+00,  1.536D+00,  1.538D+00,
     *     1.536D+00,  1.488D+00,  1.384D+00,  1.317D+00,  1.256D+00,
     *     1.153D+00,  1.089D+00,  1.185D+00,  1.255D+00,  1.315D+00,
     *     1.307D+00,  1.269D+00,  1.224D+00/
C
C*                  NEUTRON - URANIUM
      DATA (SENA(IE,10),IE=1,NEA) /
     *     3.346D+00,  2.535D+00,  2.738D+00,  3.143D+00,  4.157D+00,
     *     4.639D+00,  3.498D+00,  2.662D+00,  1.685D+00,  1.687D+00,
     *     1.685D+00,  1.632D+00,  1.518D+00,  1.445D+00,  1.378D+00,
     *     1.265D+00,  1.194D+00,  1.300D+00,  1.377D+00,  1.443D+00,
     *     1.434D+00,  1.392D+00,  1.343D+00/
C---------------------------------------------------------------------
C    p-A data changed by A.Ferrari: corresponding n-A data are used at
C        low energies since this is a much better approximation than
C        neglecting "tout court" the elastic scattering
C***     P-A ELASTIC CROSS SECTION DATA
C*                  PROTON - HELIUM
      DATA (SEPA(IE,1),IE=1,NEA) /
*    *   8*0.000D+00,                          6.700D-02,  6.100D-02,
     *     5.103D-01,  5.157D-01,  5.103D-01,  4.777D-01,  4.072D-01,
     *     3.420D-01,  2.714D-01,  1.683D-01,  6.700D-02,  6.100D-02,
     *     5.800D-02,  4.900D-02,  3.800D-02,  3.300D-02,  3.000D-02,
     *     2.400D-02,  2.300D-02,  2.900D-02,  3.600D-02,  4.100D-02,
     *     4.000D-02,  3.700D-02,  3.400D-02/
C
C*                  PROTON - BERYLLIUM
      DATA (SEPA(IE,2),IE=1,NEA) /
*    *   8*0.000D+00,                          1.401D-01,  1.305D-01,
     *     8.762D-01,  8.856D-01,  8.762D-01,  8.203D-01,  6.991D-01,
     *     5.873D-01,  4.661D-01,  2.890D-01,  1.401D-01,  1.305D-01,
     *     1.238D-01,  1.069D-01,  8.495D-02,  7.480D-02,  6.750D-02,
     *     5.565D-02,  5.230D-02,  6.470D-02,  7.765D-02,  8.722D-02,
     *     8.440D-02,  7.821D-02,  7.259D-02/
C
C*                  PROTON - CARBON
      DATA (SEPA(IE,3),IE=1,NEA) /
*    *   8*0.000D+00,                          1.820D-01,  1.710D-01,
     *     9.200D-01,  9.500D-01,  9.400D-01,  8.800D-01,  7.500D-01,
     *     6.100D-01,  5.000D-01,  3.700D-01,  1.820D-01,  1.710D-01,
     *     1.620D-01,  1.410D-01,  1.130D-01,  1.000D-01,  9.000D-02,
     *     7.500D-02,  7.000D-02,  8.600D-02,  1.020D-01,  1.140D-01,
     *     1.100D-01,  1.020D-01,  9.500D-02/
C
C*                  PROTON - ALUMINUM
      DATA (SEPA(IE,4),IE=1,NEA) /
*    *   8*0.000D+00,                          3.650D-01,  3.540D-01,
     *     1.090D+00,  1.180D+00,  1.240D+00,  1.280D+00,  1.260D+00,
     *     1.160D+00,  9.300D-01,  6.300D-01,  3.650D-01,  3.540D-01,
     *     3.420D-01,  3.060D-01,  2.530D-01,  2.260D-01,  2.040D-01,
     *     1.750D-01,  1.610D-01,  1.900D-01,  2.200D-01,  2.430D-01,
     *     2.370D-01,  2.220D-01,  2.070D-01/
C
C*                  PROTON - TITANIUM
      DATA (SEPA(IE,5),IE=1,NEA) /
*    *   8*0.000D+00,                          5.828D-01,  5.726D-01,
     *     1.029D+00,  9.469D-01,  1.091D+00,  1.284D+00,  1.591D+00,
     *     1.691D+00,  1.258D+00,  9.241D-01,  5.828D-01,  5.726D-01,
     *     5.594D-01,  5.100D-01,  4.310D-01,  3.897D-01,  3.561D-01,
     *     3.084D-01,  2.829D-01,  3.262D-01,  3.714D-01,  4.066D-01,
     *     3.985D-01,  3.764D-01,  3.517D-01/
C
C*                  PROTON - IRON
      DATA (SEPA(IE,6),IE=1,NEA) /
*    *   8*0.000D+00,                          6.383D-01,  6.313D-01,
     *     1.178D+00,  9.793D-01,  1.090D+00,  1.271D+00,  1.650D+00,
     *     1.799D+00,  1.339D+00,  1.009D+00,  6.383D-01,  6.313D-01,
     *     6.212D-01,  5.732D-01,  4.913D-01,  4.483D-01,  4.113D-01,
     *     3.594D-01,  3.304D-01,  3.764D-01,  4.243D-01,  4.623D-01,
     *     4.543D-01,  4.313D-01,  4.053D-01/
C
C*                  PROTON - COPPER
      DATA (SEPA(IE,7),IE=1,NEA) /
*    *   8*0.000D+00,                          6.950D-01,  6.895D-01,
     *     1.386D+00,  1.050D+00,  1.134D+00,  1.302D+00,  1.722D+00,
     *     1.922D+00,  1.449D+00,  1.103D+00,  6.950D-01,  6.895D-01,
     *     6.803D-01,  6.322D-01,  5.471D-01,  5.014D-01,  4.619D-01,
     *     4.048D-01,  3.728D-01,  4.211D-01,  4.722D-01,  5.135D-01,
     *     5.051D-01,  4.804D-01,  4.527D-01/
C
C*                  PROTON - CADMIUM
      DATA (SEPA(IE,8),IE=1,NEA) /
*    *   8*0.000D+00,                          1.045D+00,  1.043D+00,
     *     2.029D+00,  1.537D+00,  1.660D+00,  1.906D+00,  2.520D+00,
     *     2.812D+00,  2.121D+00,  1.614D+00,  1.045D+00,  1.043D+00,
     *     1.036D+00,  9.718D-01,  8.822D-01,  8.211D-01,  7.679D-01,
     *     6.828D-01,  6.325D-01,  6.951D-01,  7.647D-01,  8.232D-01,
     *     8.138D-01,  7.935D-01,  7.415D-01/
C
C*                  PROTON - LEAD
      DATA (SEPA(IE,9),IE=1,NEA) /
*    *   8*0.000D+00,                          1.589D+00,  1.584D+00,
     *     3.050D+00,  2.310D+00,  2.495D+00,  2.865D+00,  3.789D+00,
     *     4.228D+00,  3.188D+00,  2.426D+00,  1.589D+00,  1.584D+00,
     *     1.577D+00,  1.528D+00,  1.417D+00,  1.345D+00,  1.277D+00,
     *     1.159D+00,  1.086D+00,  1.159D+00,  1.252D+00,  1.331D+00,
     *     1.320D+00,  1.278D+00,  1.256D+00/
C
C*                  PROTON - URANIUM
      DATA (SEPA(IE,10),IE=1,NEA) /
*    *   8*0.000D+00,                          1.743D+00,  1.738D+00,
     *     3.346D+00,  2.535D+00,  2.738D+00,  3.143D+00,  4.157D+00,
     *     4.639D+00,  3.498D+00,  2.662D+00,  1.743D+00,  1.738D+00,
     *     1.730D+00,  1.676D+00,  1.554D+00,  1.475D+00,  1.401D+00,
     *     1.271D+00,  1.191D+00,  1.271D+00,  1.373D+00,  1.460D+00,
     *     1.448D+00,  1.402D+00,  1.378D+00/
C
      DATA RELTO / 0.175D+00, 6*0.D+00, 0.175D+00, 4*0.D+00, 0.14D+00,
     *             0.14 D+00/
C
C--------------------------------------------------------------------

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C--------------------------------------------------------------------
C
      IF(ANUC.LT.1.5D0) THEN
C                               HADRON-PROTON ELASTIC CROSS SECTIONS
        IPOL=0
        EK1=EKIN
        IF(EKIN.GT.20.0D0) THEN
          SIGELA=RELTO(KPROJ)*DSHPTO(KPROJ,PLAB)
          RETURN
        ELSEIF(EKIN.GT.10.0D0) THEN
          IPOL=1
          PO2=20.0
        EK2=SQRT(PO2**2+AAM(KPROJ)**2) - AAM(KPROJ)
          SEL2=RELTO(KPROJ)*DSHPTO(KPROJ,PO2)
          EK1=10.0
C
         ELSE IF ( EKIN .LT. 0.06D+00 ) THEN
            IF ( KPROJ .EQ. 1 ) THEN
               IF ( EKIN .LT. 0.02D+00 ) THEN
                  SIGELA = ONETHI * SIGMLW (EKIN)
               ELSE
                  BETAPR = PLAB / ( EKIN + AAM (KPROJ) )
                  SIGELA = SIGMPP (BETAPR)
               END IF
               RETURN
            ELSE IF ( KPROJ .EQ. 8 ) THEN
               IF ( EKIN .LT. 0.04D+00 ) THEN
                  SIGELA = SIGMLW (EKIN)
               ELSE
                  BETAPR = PLAB / ( EKIN + AAM (KPROJ) )
                  SIGELA = SIGMPN (BETAPR)
               END IF
               RETURN
            END IF

        ENDIF
C
        DO 10 IE=1,NEN
          IF(EK1.LT.EKIHN(IE)) THEN
            JE1=IE-1
            JE2=IE
C           WRITE (6,'(F12.8,I5,F12.8,I5,F12.8)')EK1,IE,
C    *                                           EKIHN(IE),KPROJ,EKIN
            DDEE=EKIHN(JE2) - EKIHN(JE1)
                                                                 GOTO 20
          ENDIF
   10   CONTINUE
        JE1=NEN
        JE2=NEN
        DDEE=1.
   20   CONTINUE
C****
C                                  PROTON-PROTON
        IF(KPROJ.EQ.1) THEN
          S1=SEPP(JE1)
          S2=SEPP(JE2)
C                                  NEUTRON-PROTON
        ELSEIF(KPROJ.EQ.8) THEN
          S1=SENP(JE1)
          S2=SENP(JE2)
C                                  PI(+)-PROTON
        ELSEIF(KPROJ.EQ.13) THEN
          S1=SEPIPP(JE1)
          S2=SEPIPP(JE2)
C                                  PI(-)-PROTON
        ELSEIF(KPROJ.EQ.14) THEN
          S1=SEPIMP(JE1)
          S2=SEPIMP(JE2)
C                                  UNDEFINED ENTRY CONDITIONS
        ELSE
          SIGELA=0.
          RETURN
        ENDIF
C
        SIGELA=S1 + (S2-S1)*(EK1-EKIHN(JE1))/DDEE
C
C                                  INTERPOLATION BETWEEN 10/20 GEV
        IF(IPOL.EQ.1) THEN
          SEL1=SIGELA
          SIGELA=SEL1 + (SEL2-SEL1)*(EKIN-EK1)/(EK2-EK1)
        ENDIF
C
        RETURN
C
      ENDIF
C***************************************
C                               HADRON-NUCLEUS ELASTIC CROSS SECTIONS
      DO 30 IE=1,NEA
        IF(EKIN.LT.EKIHA(IE)) THEN
          JE=IE - 1
                                                                 GOTO 40
        ENDIF
   30 CONTINUE
      IF(EKIN.EQ.EKIHA(NEA)) THEN
        JE=NEA - 1
      ELSE
        JE=-1
      ENDIF
   40 CONTINUE
C
      DO 50 IA=1,NNAA
        IF(ANUC.LT.AMASS(IA)) THEN
          JA=IA - 1
                                                                 GOTO 60
        ENDIF
   50 CONTINUE
      IF(ANUC.EQ.AMASS(NNAA)) THEN
        JA=NNAA - 1
      ELSE
        JA=-1
      ENDIF
   60 CONTINUE
C
      IF (JA) 140,110,70
   70 IF (JE) 190,150,80
   80 TEMP1=ANUC/AMASS(JA)
      TEMP2=LOG(AMASS(JA+1)/AMASS(JA))
      KE=JE
      DO 90 I=1,2
        IF(KPROJ.EQ.8) THEN
          SLOW=SENA(KE,JA)
          POWER=LOG(SENA(KE,JA+1)/SLOW)/TEMP2
        ELSE
          SLOW=SEPA(KE,JA)
          POWER=LOG(SEPA(KE,JA+1)/SLOW)/TEMP2
        ENDIF
        TSIG(I)=SLOW*TEMP1**POWER
        KE=KE+1
   90 CONTINUE
C
  100 SIGELA=TSIG(1) + (EKIN-EKIHA(JE))*(TSIG(2)-TSIG(1)) /(EKIHA(JE+1)
     +-EKIHA(JE))
 
      SIGELA=SIGELA * 1D3
      RETURN
C*
C                                          A IS LESS THAN A MIN
  110 JA=1
      TEMP1= (ANUC/AMASS(JA)) **0.66667D0
  120 IF (JE) 200,170,130
  130 IF(KPROJ.EQ.8) THEN
        TSIG(1) = SENA(JE,JA) * TEMP1
        TSIG(2) = SENA(JE+1,JA) *TEMP1
      ELSE
        TSIG(1) = SEPA(JE,JA) * TEMP1
        TSIG(2) = SEPA(JE+1,JA) *TEMP1
      ENDIF
                                                               GO TO 100
C*
C                                         A IS GREATER THAN A MAX
  140 JA=NNAA
      TEMP1= (ANUC/AMASS(JA))**.66667
                                                               GO TO 120
C*
C                                         EKIN  LT.  EMIN
  150 JE=1
  160 TEMP1=ANUC/AMASS(JA)
      TEMP2=LOG(AMASS(JA+1)/AMASS(JA))
      IF(KPROJ.EQ.8) THEN
        SLOW=SENA(JE,JA)
        POWER=LOG(SENA(JE,JA+1)/SLOW)/TEMP2
      ELSE
        SLOW=SEPA(JE,JA)
        POWER=LOG(SEPA(JE,JA+1)/SLOW)/TEMP2
      ENDIF
      SIGELA=SLOW*TEMP1**POWER
      SIGELA=SIGELA * 1D3
      RETURN
C
  170 JE=1
  180 IF(KPROJ.EQ.8) THEN
        SIGELA=SENA(JE,JA)*TEMP1
      ELSE
        SIGELA=SEPA(JE,JA)*TEMP1
      ENDIF
      SIGELA=SIGELA * 1D3
      RETURN
C*
C                                         EKIN  GT.  EMAX
  190 JE=NEA
                                                               GO TO 160
  200 JE=NEA
                                                               GO TO 180
      END
*-- Author :
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      DOUBLE PRECISION FUNCTION DSHPTO(IT,PO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------
C   TOTAL HADRON-PROTON CROSS SECTIONS
C                                        Version based on DSHNTO
C                                        in d4diff.f
C                                        used after October 1993
CC++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
      AMIT2=AAM(IT)**2
      UMO2=AMIT2 + AAM(1)**2 + 2.*AAM(1)*(PO+0.5*AMIT2/PO)
      UMO=SQRT(UMO2)
      DSHPTO=DSHNTO(IT,1,UMO)
C
      RETURN
      END
C*******************************************************************
      DOUBLE PRECISION FUNCTION XSHPTO(IT,PO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------
C   TOTAL HADRON-PROTON CROSS SECTIONS
C   PLAB.GE.10 GEV
C                                       This is the version used before
C                                       October 1993
C   PARAMETRZATION OF HJM 1982
C   CORRECTION FOR ANTIPARTICLES 16/08/90 HJM
C
C   extension below 10 GeV/c   29/08/90   HJM
C      - Plab > 4 GeV/c : Fit from REv. of Part. Prop.
C      - Plab < 4 GeV/c : SIG(el) + SIG(inel) = SIHNEL + SIHNIN
C----------------------------------------------------
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
C---------
      POO=PO
      KPRO=IT
      IF(POO.LT.4.0) THEN
        CALL SIHNEL(KPRO,1,POO,SIEL)
        CALL SIHNIN(KPRO,1,POO,SIIN)
        XSHPTO=SIEL + SIIN
        RETURN
      ENDIF
C-----------------------------------------------------
C                                     J.R.10.2.91
C                             SPECIAL TREATMENT OF P-P CROSS-SECTIONS
      IF (IT.EQ.1)THEN
        IF (POO.LT.2100.)THEN
          XSHPTO=45.6+219.*POO**(-4.23)+0.41*LOG(POO)**2-3.41*LOG(POO)
          RETURN
        ELSEIF(POO.LE.432000.) THEN
          XSHPTO=41.1+77.2*POO**(-0.68)+0.293*LOG(POO)**2-1.82*LOG(POO)
          RETURN
        ELSE
        AMIT2=AAM(1)**2
        UMO2=AMIT2 + AAM(1)**2 + 2.*AAM(1)*(PO+0.5*AMIT2/PO)
          UMO=SQRT(UMO2)
          IF (UMO.LT.1800.) THEN
            XSHPTO=61.65+(UMO-900.)*6.6/900.
            RETURN
          ELSEIF(UMO.LT.16000.)THEN
            XSHPTO=68.3+(UMO-1800.)*25./14200.
            RETURN
          ELSE
            XSHPTO=93.3+(UMO-16000.)*12.5/24000.
            RETURN
          ENDIF
        ENDIF
      ENDIF
C------------------------------------------------------
      A4=0.
      A5=0.
      A6=0.
      ITT=IT
C
      IF(POO.LT.10.0) THEN
        IPIO=0
        GOTO (10,20,260,260,260,260,260,30,40,260,260,50,60,70,80,90,
     +  100,20,50, 10,10,10,110,120,130), ITT
C***
   10   CONTINUE
C                               P - P
C                               SIGMA - P
        A1=45.6
        A2=219.
        A3=-4.23
        A4=0.410
        A5=-3.41
                                                                GOTO 140
C***
   20   CONTINUE
C                               PBAR - P
C                               LAMBDA-BAR  -  P
        A1=41.1
        A2=77.2
        A3=-0.68
        A4=0.293
        A5=-1.82
                                                                GOTO 140
C***
   30   CONTINUE
C                               N - P
        A1=47.7
        A2=100.
        A3=-4.57
        A4=0.512
        A5=-4.29
                                                                GOTO 140
C***
   40   CONTINUE
C                               NBAR - P  =  PBAR - N
        A1=41.9
        A2=96.2
        A3=-0.99
        A4=-0.154
        A5=0.0
                                                                GOTO 140
C***
   50   CONTINUE
C                               KLONG - P  =  KSHORT - P
        IF(RNDM(V).GE.0.5)                                       GOTO 80
                                                                 GOTO 90
C***
   60   CONTINUE
C                               PI+ - P
        A1=32.1
        A2=48700.
        A3=-7.85
        A4=0.540
        A5=-4.41
                                                                GOTO 140
C***
   70   CONTINUE
C                               PI-  P
        A1=33.1
        A2=15.0
        A3=-1.41
        A4=0.458
        A5=-4.06
                                                                GOTO 140
C***
   80   CONTINUE
C                               K+ - P
        A1=17.1
        A2=5.54
        A3=-2.67
        A4=0.139
        A5=-0.270
                                                                GOTO 140
C***
   90   CONTINUE
C                               K-  P
        A1=-21.1
        A2=56.2
        A3=-0.27
        A4=-0.155
        A5=6.24
                                                                GOTO 140
C***
  100   CONTINUE
C                               LAMBDA  -  P
        A1=18.0
        A2=0.121
        A3=-3.92
        A4=6.38
        A5=0.0
                                                                GOTO 140
C***
  110   CONTINUE
C                               PIZERO - P
C                               PRESENTLY PI0 = PI+
C                                could be   0.5 (PI+ & PI-) with IPIO=1
        IPIO=0
                                                                 GOTO 60
C***
  120   CONTINUE
C                               K0 - P  =  K+  - P
                                                                 GOTO 80
C                               K0 - P  =  K+  - N
C         A1=18.4
C         A2=175
C         A3=-7.85
C         A4=0.198
C         A5=-0.753
C         GOTO 200
C***
  130   CONTINUE
C                               K0BAR - P  =  K-  - P
                                                                 GOTO 90
C                               K0BAR - P  =  K-  - N
C         A1=-1040.
C         A2=1060.
C         A3=-0.03
C         A4=0.0
C         A5=27.8
C***
  140   CONTINUE
        ALP=LOG(POO)
        XSHPTO=A1 + A2*POO**A3 + A4*ALP**2 + A5*ALP
        IF(IPIO.EQ.1) THEN
          S1=XSHPTO
          IPIO=2
                                                                 GOTO 70
        ELSEIF(IPIO.EQ.2) THEN
          XSHPTO=0.5*(XSHPTO+S1)
        ENDIF
        RETURN
      ENDIF
C
      F1=1.
      AMIT2=AAM(ITT)**2
      UMO2=AMIT2 + AAM(1)**2 + 2.*AAM(1)*(PO+0.5*AMIT2/PO)
      UMO=SQRT(UMO2)
C
      GOTO (150,150,260,260,260,260,260,160,160,260,260,210,170,170,190,
     +190,240,240,210,240,240,240,250,220,230), ITT
C
  150 CONTINUE
C                               P-P
      A1=38.4
      A2=0.46
      A3=125.
      IF(ITT.EQ.1)                                              GOTO 270
C                               PBAR-P
      A5=84.1
      A6=-0.57
                                                                GOTO 270
C
  160 CONTINUE
C                              N-P = P-N
      A1=38.5
      A2=0.46
      A3=125.
      A4=15.
      IF(ITT.EQ.8)                                              GOTO 270
C                              NBAR-P = PBAR-N
      A5=77.43
      A6=-0.60
                                                                GOTO 270
C
  170 CONTINUE
      IF(UMO.LT.47.0)                                           GOTO 180
      F1=0.6667
      ITT=1
                                                                GOTO 150
  180 CONTINUE
C                             (PI-) - P
      A1=24.0
      A2=0.60
      A3=160.
      IF(ITT.EQ.14)                                             GOTO 270
C                             (PI+) - P
      A5=-7.9
      A6=-0.46
                                                                GOTO 270
C
  190 CONTINUE
      IF(UMO.LT.110.)                                           GOTO 200
      F1=0.6667
      ITT=1
                                                                GOTO 150
  200 CONTINUE
C                             (K-) - P
      A1=20.3
      A2=0.59
      A3=140.
      IF(ITT.EQ.16)                                             GOTO 270
C                             (K+) - P
      A5=-30.13
      A6=-0.58
                                                                GOTO 270
C
  210 CONTINUE
      ITT=15
      IF(RNDM(V).LT.0.5) ITT=16
                                                                GOTO 190
C
  220 CONTINUE
C***
C   K-ZERO:  SET EQUAL TO K+/PROTON
C            (SHOULD BE K+/NEUTRON)
C***
      ITT=15
                                                                GOTO 190
C
  230 CONTINUE
C***
C   K-ZERO BAR:  SET EQUAL TO K-/PROTON
C                (SHOULD BE K-/NEUTRON)
C***
      ITT=16
                                                                GOTO 190
C
  240 CONTINUE
C***
C   SIGMA +/-/0  AND  LAMBDA/LAMBDA BAR:  SET EQUAL TO P-P
C***
      ITT=1
                                                                GOTO 150
C
  250 CONTINUE
C***
C   PI0:   SET EQUAL TO PI+
C***
      ITT=13
                                                                GOTO 170
C
  260 CONTINUE
C***
C   LEPTONS AND PI0
C***
      XSHPTO=1.E-10
      RETURN
C
  270 CONTINUE
C
      XSHPTO=A1 + A2*(LOG(UMO2/A3))** 2+ A4/UMO2 + A5*UMO2**A6
 
      XSHPTO=F1*XSHPTO
C
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE SIHNIN(IPROJ,ITAR,PO,SIIN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C********************************************************************
C     VERSION OCTOBER 89 BY          H. MOEHRING
C                                    LEIPZIG
C     LAST CHANGE 30/08/90 BY HJM: modification of SIG(pi+ p) = SEGP
C
C
C     SEE: H.J. MOEHRING, HADRON-NUCLEUS INELASTIC CROSS-SECTIONS FOR
C     USE IN HADRON-CASCADE CALCULATIONS AT HIGH ENERGIES,
C     TIS DIVISION REPORT 14. OCTOBER 1983, TIS-RP/116, CERN GENEVA
C
C     INPUT VARIABLES:
C        IPROJ  = TYPE OF THE PROJECTILE
C        ITAR   = TYPE OF THE TARGET
C        PO     = PARTICLE MOMENTUM IN GEV/C
C
C     OUTPUT VARIABLES:
C        SIIN   = INTERPOLATED INELASTIC CROSS SECTION IN MILLIBARNS
C
C
C     OTHER IMPORTANT VARIABLES:
C        SIG    = PROTON/NUCLEI CROSS SECTIONS
C        SEG    = PION-/NUCLEI CROSS SECTIONS ABOVE 0.3 GEV/C
C        SEGP   = PION+/NUCLEI CROSS SECTIONS ABOVE 0.3 GEV/C
C        SIGKM  = K+ AND K0/NUCLEI CROSS SECTIONS ABOVE 3.0 GEV/C
C        SIGKP  = K+ AND K0 BAR/NUCLEI CROSS SECTIONS ABOVE 3.0 GEV/C
C        SIGAP  = ANTINUCLEON/NUCLEI CROSS SECTIONS ABOVE 3.0 GEV/C
C        SEEG   = PION/NUCLEI CROSS SECTIONS BELOW 0.3 GEV/C
C        P      = MOMENTA FOR WHICH THE CROSS SECTIONS ARE GIVEN IN
C                 SIG, SEG, SEGP, SIGKM, SIGKP AND SIGAP
C        PEE    = MOMENTA FOR WHICH THE CROSS SECTIONS ARE GIVEN IN
C                 SEEG
C        A      = NUCLEI FOR WHICH THE CROSS SECTIONS ARE GIVEN IN
C                 SIG, SEG, SEGP, SIGKM, SIGKP, SIGAP AND SEEG
C        PLAB   = MOMENTA FOR WHICH THE TOTAL CROSS SECTIONS ARE
C                 GIVEN IN TOTCRS
C        TOTCRS = TOTAL CROSS SECTIONS AS A FUNCTION OF MOMENTUM
C                 TOTCRS(K,I) WHERE K=MOMENTUM INDEX,I=REACTION TYPE
C                 I=1:NEGATIVE KAON-PROTON  = KAON ZERO BAR-NEUTRON
C                 I=2:NEGATIVE KAON-NEUTRON = KAON ZERO BAR-PROTON
C                 I=3:POSITIVE KAON-PROTON  = KAON ZERO NEUTRON
C                 I=4:POSITIVE KAON-NEUTRON = KAON ZERO-PROTON
C                 I=5:ANTI NUCLEON-NUCLEON
C
C
C     NOTE1: PRESENTLY CROSS SECTIONS ARE ASSUMED TO BE CONSTANT
C     ABOVE 10000.0 GEV/C FOR ALL PARTICLES AND
C     BELOW 0.13 GEV/C FOR PIONS AND BELOW 0.3 GEV/C FOR OTHERS
C
C     NOTE2: SEE TABLE ITT TO FIND OUT HOW DIFFERENT HADRONS
C     ARE TREATED. ALL PARTICLES WITH PARTICLE NUMBER BIGGER THAN
C     25 ARE TREATED AS PROTONS.
C
C     NOTE3: FOR LEPTONS AND PHOTONS PRACTICALLY ZERO CROSS SECTION
C     IS RETURNED.
C********************************************************************
C
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEND.
      DIMENSION SEEG(4),PEE(4),SIGKP(20),SIGKM(20),SIGAP(20)
      DIMENSION SEG(20),SIG(20),SEGP(20),P(20)
C
      DATA P/0.3D0,0.4D0,0.5D0,0.6D0,0.8D0,1.D0,
     +1.5D0,2.D0,3.D0,4.D0,5.D0,6.D0,10.D0,
     *20.D0,50.D0,100.D0,200.D0,400.D0,1000.D0,10000.D0/
C
      DATA SEEG/0.1D0,16.D0,35.D0,42.D0/
      DATA PEE /0.13D0,0.19D0,0.25D0,0.30D0/
C
      DATA (SIG(IE),IE=1,20) / 3*0.0001D0,0.1D0,1.D0,4.D0,
     +24.5D0,25.D0,27.2D0,27.8D0,
     +28.5D0,29.2D0,29.7D0, 30.5D0,31.5D0,
     +31.7D0,32.1D0,32.9D0,34.5D0,41.2D0/
 
C
      DATA (SEG(IE),IE=1,20) / 42.D0,19.D0,16.1D0,17.D0,
     +22.7D0,32.5D0,24.6D0,26.2D0,
     +25.0D0,23.7D0,23.0D0,22.5D0,22.D0,
     + 21.2D0,20.8D0,20.7D0,21.D0,21.9D0,23.8D0,28.4D0/
 
C
      DATA (SEGP(IE),IE=1,20) / 0.1D0,0.1D0,0.1D0,0.1D0,
     +11.D0,12.5D0,22.D0,19.2D0,21.5D0,
     +21.4D0,20.8D0,20.6D0,20.2D0, 19.8D0,19.9D0,
     +20.D0,20.5D0,21.5D0,23.5D0,32.2D0/
 
C
      DATA SIGAP/ 164.D0,126.D0,114.D0,98.D0,86.D0,
     +72.4D0,59.D0,57.D0,53.D0,52.D0,48.D0,45.5D0,
     +43.5D0, 40.4D0,36.5D0,35.2D0,34.5D0,34.5D0,
     +35.4D0,41.5D0/
 
C
      DATA SIGKM/ 38.D0,43.D0,23.D0,18.5D0,20.D0,
     +29.D0,25.D0,23.D0,22.5D0,21.D0,20.5D0,20.D0,
     +19.2D0, 18.5D0,17.8D0,17.8D0,18.3D0,19.2D0,
     +21.2D0,28.9D0/
 
C
      DATA SIGKP/ 0.001D0,0.001D0,0.001D0,0.001D0,
     +0.2D0,4.5D0,8.9D0,11.6D0,12.2D0,13.4D0,
     +13.6D0, 13.7D0,13.7D0,14.9D0,15.9D0,16.5D0,
     +17.4D0,18.6D0,20.9D0,28.8D0/
 
C
C---------------------------------------------------------------------
C
      SIIN=1.0D-20
      IF(ITAR.NE.1.AND.ITAR.NE.8) THEN
        IF(IPRI.GE.1) WRITE(6,'(A/A,2I5,2(1PE15.6))')
     +  ' WRONG CALL OF SIHNIN/ITAR',
     +  '    IPROJ,ITAR,PO,SIIN :', IPROJ,ITAR,PO,SIIN
 
        RETURN
      ENDIF
C
      IF(IPROJ.GE.3.AND.IPROJ.LE.7) RETURN
C
      IIPP=IPROJ
      IF(IPROJ.EQ.23) IIPP=13
      IITT=ITAR
C-----------------------------------------------------------------------
C                                      NEUTRON TARGET TO BE IMPLEMENTED!
      IF(IITT.EQ.8) IITT=1
C-----------------------------------------------------------------------
      IF(IPROJ.EQ.19.OR.IPROJ.EQ.12) THEN
        IIPP=24
        RND=RNDM(V)
        IF(RND.GT.0.5) IIPP=25
      END IF
C
C********************************************************************
C
C     CALCULATE THE MOMENTUM INDEX K
C
      DO 10 JK=1,20
        IF(PO.LE.P(JK)) THEN
          K=JK
          KK=K-1
                                                                GO TO 20
        ENDIF
   10 CONTINUE
      K=21
   20 CONTINUE
C
C*******************************************************************
C
      IF(IITT.EQ.1) THEN
C                                                PROTON TARGET
        IF(IIPP.EQ.1.OR.IIPP.EQ. 8.OR.IIPP.EQ.17.OR.(IIPP.GE.20.AND
     +  .IIPP.LE.22)) THEN
C                                        PROTON/NEUTRON
C                                        LAMBDA/ALL SIGMAS
          IF(K.EQ.1) THEN
            SIIN=SIG(K)
                                                                 GOTO 70
          ELSEIF(K.EQ.21) THEN
            SIIN=SIG(20)
                                                                 GOTO 70
          ELSE
            SI1=SIG(K)
            SI2=SIG(KK)
          ENDIF
        ELSEIF(IIPP.EQ.14) THEN
C                                        PI -
          IF(K.EQ.1) THEN
C                                        LOW ENERGY PI- (<0.3GEV/C)
            DO 30 JK=1,4
              IF(PO.LE.PEE(JK)) THEN
                kkK=JK
                                                                 GOTO 40
              ENDIF
   30       CONTINUE
            kkK=4
   40       CONTINUE
            KK=kkK-1
            IF(kkK.EQ.1) THEN
              SIIN=SEEG(kkK)
                                                                 GOTO 70
            ENDIF
            SI1=SEEG(kkK)
            SI2=SEEG(KK)
C                         INTERPOLATE LINEARLY WITH RESPECT TO MOMENTUM
C
            SIIN=SI1 + (PO-PEE(kkK))*(SI2-SI1)/(PEE(KK)-PEE(kkK))
                                                                 GOTO 70
          ELSEIF(K.EQ.21) THEN
            SIIN=SEG(20)
                                                                 GOTO 70
          ELSE
            SI1=SEG(K)
            SI2=SEG(KK)
          ENDIF
C                                                      PI +
        ELSEIF(IIPP.EQ.13) THEN
   50     CONTINUE
          IF(K.EQ.1) THEN
            SIIN=SEGP(K)
                                                                 GOTO 70
          ELSEIF(K.EQ.21) THEN
            SIIN=SEGP(20)
                                                                 GOTO 70
          ENDIF
          SI1=SEGP(K)
          SI2=SEGP(KK)
C                                                 K -  AND K0 BAR
        ELSEIF(IIPP.EQ.16.OR.IIPP.EQ.25) THEN
          IF(K.EQ.1) THEN
            SIIN=SIGKM(K)
                                                                 GOTO 70
          ELSEIF(K.EQ.21) THEN
            SIIN=SIGKM(20)
                                                                 GOTO 70
          ENDIF
          SI1=SIGKM(K)
          SI2=SIGKM(KK)
C                                                 K +  AND K0
        ELSEIF(IIPP.EQ.15.OR.IIPP.EQ.24) THEN
          IF(K.EQ.1) THEN
            SIIN=SIGKP(K)
                                                                 GOTO 70
          ELSEIF(K.EQ.21) THEN
            SIIN=SIGKP(20)
                                                                 GOTO 70
          ENDIF
          SI1=SIGKP(K)
          SI2=SIGKP(KK)
C                                                 ANTI-NUCLEONS
C                                                 ANTI-LAMBDA
        ELSEIF(IIPP.EQ.2.OR.IIPP.EQ.9.OR.IIPP.EQ.18) THEN
          IF(K.EQ.1) THEN
            SIIN=SIGAP(K)
                                                                 GOTO 70
          ELSEIF(K.EQ.21) THEN
            SIIN=SIGAP(20)
                                                                 GOTO 70
          ENDIF
          SI1=SIGAP(K)
          SI2=SIGAP(KK)
        ENDIF
      ENDIF
C
C********************************************************************
C     INTERPOLATE LINEARLY WITH RESPECT TO MOMENTUM
C
   60 SIIN=SI1 + (PO-P(K))*(SI1-SI2)/(P(K)-P(KK))
C
   70 CONTINUE
C     ZL=10000.D0*AA/(6.022D0*SIIN)
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE DHADRI(N,PLAB,ELAB,CX,CY,CZ,ITTA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C-----------------------------
C*** INPUT VARIABLES LIST:
C*** SAMPLING OF HADRON NUCLEON INTERACTION FOR (ABOUT) 0.1 LE PLAB LE 6
C*** GEV/C LABORATORY MOMENTUM REGION
C*** N    - PROJECTILE HADRON INDEX
C*** PLAB - LABORATORY MOMENTUM OF N (GEV/C)
C*** ELAB - LABORATORY ENERGY OF N (GEV)
C*** CX,CY,CZ - DIRECTION COSINES OF N IN THE LABORATORY SYSTEM
C*** ITTA - TARGET NUCLEON INDEX
C*** OUTPUT VARIABLES LIST OF PARTICLE CHARACTERISTICS IN /FINLSP/
C  IR COUNTS THE NUMBER OF PRODUCED PARTICLES
C*** ITR - PARTICLE INDEX, CXR,CYR,CZR - DIRECTION COSINES (LAB. SYST.)
C*** ELR,PLR LAB. ENERGY AND LAB. MOMENTUM OF THE SAMPLED PARTICLE
C*** RESPECT., UNITS (GEV/C AND GEV)
C----------------------------
      COMMON /DGAMRE/ REDU,AMO,AMM(15 )
      COMMON /DREDVE/ THRESH(268), IRII(17),IKII(17),IEII(17)
      COMMON/DREAC/UMO(296),PLABF(296),SIIN(296),WK(5184),
     *NRK(2,268),NURE(30,2)
      COMMON /DABLTI/ AMH(110),GAH(110),TAUH(110),ICHH(110),IBARH(110),
     +K1H(110),K2H(110)
      COMMON /DSPLI/NZK(460,3),WT(460)
      COMMON /DMETLS/ CXS(149),CYS(149), 
     +CZS(149),ELS(149),PLS(149),
     +IS,ITS(149)
      COMMON /DRUN/ RUNTES,EFTES
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
*KEEP,DFINLS.
      PARAMETER (MAXFIN=10)
      COMMON /DFINLS/ ITRH(MAXFIN),CXRH(MAXFIN),CYRH(MAXFIN), CZRH
     +(MAXFIN),ELRH(MAXFIN),PLRH(MAXFIN),IRH
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEND.
      DIMENSION ITPRF(110)
      DATA NNN/0/
      DATA UMODA/0./
      DATA ITPRF/-1,-1,5*1,-1,-1,1,1,1,-1,-1,-1,-1,6*1,-1,-1,-1,85*1/
      LOWP=0
      IF (N.LE.0.OR.N.GE.111)N=1
      IF (ITPRF( N ).GT.0 .OR. ITTA.GT.8) THEN
        GOTO 280
*       WRITE (6,1000)
*    +  ' FALSE USE OF THE PARTICLE TYPE INDEX: N, ITTA', N, ITTA
*       STOP
*1000   FORMAT (3(5H ****/),A,2I4,3(5H ****/))
*    +  45H FALSE USE OF THE PARTICLE TYPE INDEX, N,LUE ,I4,3(5H ****/))
      ENDIF
      IATMPT=0
      IF (ABS(PLAB-5.).LT.4.99999D0)                            GO TO 20
      IF(IPRI.GE.1) WRITE (6,1010) PLAB
C     STOP
 1010 FORMAT ( '  PROJECTILE HADRON MOMENTUM OUTSIDE OF THE
     + ALLOWED REGION, PLAB=',1E15.5)
 
   20 CONTINUE
      UMODAT=N*1.11111D0+ITTA*2.19291D0
      IF(UMODAT.NE.UMODA) CALL DCALUM(N,ITTA)
      UMODA=UMODAT
   30 IATMPT=0
      LOWP=LOWP+1
   40 CONTINUE
      IMACH=0
      REDU=2.
      IF (LOWP.GT.20)                                          GO TO 280
      NNN=N
      IF (NNN.EQ.N)                                             GO TO 50
      RUNTES=0.
      EFTES=0.
   50 CONTINUE
      IS=1
      IRH=0
      IST=1
      NSTAB=23
      IRE=NURE(N,1)
      IF(ITTA.GT.1) IRE=NURE(N,2)
C
C-----------------------------
C*** IE,AMT,ECM,SI DETERMINATION
C----------------------------
      CALL DSIGIN(IRE,PLAB,N,IE,AMT ,AMN,ECM,SI,ITTA)
      IANTH=-1
      IF (AMH(1).NE.0.9383D0) IANTH=1
      IF (IANTH.GE.0) SI=1.
      ECMMH=ECM
C
C-----------------------------
C    ENERGY INDEX
C  IRE CHARACTERIZES THE REACTION
C  IE IS THE ENERGY INDEX
C----------------------------
      IF (SI.LT.1.D-6)                                         GO TO 280
      IF (N.LE.NSTAB)                                           GO TO 60
      RUNTES=RUNTES+1.
      IF (RUNTES.LT.20.D0) WRITE(6,1020)N
 1020 FORMAT(3H N=,I10,30H THE PROEKTILE IS A RESONANCE )
      IF(IBARH(N).EQ.1) N=8
      IF(IBARH(N).EQ.-1)  N=9
   60 CONTINUE
      IMACH=IMACH+1
      IF (IMACH.GT.10)                                         GO TO 280
      ECM =ECMMH
      AMN2=AMN**2
      AMT2=AMT**2
      ECMN=(ECM**2+AMN2-AMT2)/(2.*ECM    )
      IF(ECMN.LE.AMN) ECMN=AMN
      PCMN=SQRT(ECMN**2-AMN2)
      GAM=(ELAB+AMT)/ECM
      BGAM=PLAB/ECM
      IF (IANTH.GE.0) ECM=2.1
C
C-----------------------------
C*** RANDOM CHOICE OF REACTION CHANNEL
C----------------------------
      IST=0
      CALL
     *DRANDM(VV)
      VV=VV-1.D-17
C
C-----------------------------
C***  PLACE REDUCED VERSION
C----------------------------
      IIEI=IEII(IRE)
      IDWK=IEII(IRE+1)-IIEI
      IIWK=IRII(IRE)
      IIKI=IKII(IRE)
C
C-----------------------------
C***  SHRINKAGE TO THE CONSIDERED ENERGY REGION FOR THE USE OF WEIGHTS
C----------------------------
      HECM=ECM
      HUMO=2.*UMO(IIEI+IDWK)-UMO(IIEI+IDWK-1)
      IF (HUMO.LT.ECM) ECM=HUMO
C
C-----------------------------
C*** INTERPOLATION PREPARATION
C----------------------------
      ECMO=UMO(IE)
      ECM1=UMO(IE-1)
      DECM=ECMO-ECM1
      DEC=ECMO-ECM
C
C-----------------------------
C*** RANDOM LOOP
C----------------------------
      IK=0
      WKK=0.
      WICOR=0.
   70 IK=IK+1
      IWK=IIWK+(IK-1)*IDWK+IE-IIEI
      WOK=WK(IWK)
      WDK=WOK-WK(IWK-1)
C
C-----------------------------
C*** TESTVARIABLE WICO/WICOR: IF CHANNEL IK HAS THE SAME WEIGHTS LIKE IK
C    GO TO NEXT CHANNEL, BECAUSE WKK((IK))-WKK((IK-1))=0, IK CAN NOT
C    CONTRIBUTE
C----------------------------
      IF (PLAB.LT.PLABF(IIEI+2)) WDK=0.
      WICO=WOK*1.23459876D0+WDK*1.735218469D0
      IF (WICO.EQ.WICOR)                                        GO TO 70
      IF (UMO(IIEI+IDWK).LT.HECM) WDK=0.
      WICOR=WICO
C
C-----------------------------
C*** INTERPOLATION IN CHANNEL WEIGHTS
C----------------------------
      EKLIM=-THRESH(IIKI+IK)
      IELIM=IEFUND(EKLIM,IRE)
      DELIM=UMO(IELIM)+EKLIM
     *+1.E-16
      DETE=(ECM-(ECMO-EKLIM)*.5)*2.
      IF (DELIM*DELIM-DETE*DETE) 90,90,80
   80 DECC=DELIM
                                                               GO TO 100
   90 DECC=DECM
  100 CONTINUE
      WKK=WOK-WDK*DEC/(DECC+1.D-9)
C
C-----------------------------
C*** RANDOM CHOICE
C----------------------------
C
      IF (VV.GT.WKK)                                            GO TO 70
C
C***IK IS THE REACTION CHANNEL
C----------------------------
      INRK=IKII(IRE)+IK
      ECM=HECM
      I1001 =0
C
  110 CONTINUE
      IT1=NRK(1,INRK)
      AM1=DAMG(IT1)
      IT2=NRK(2,INRK)
      AM2=DAMG(IT2)
      AMS=AM1+AM2
      I1001=I1001+1
      IF (I1001.GT.50)                                          GO TO 60
C
      IF (IT2*AMS.GT.IT2*ECM)                                  GO TO 110
      IT11=IT1
      IT22=IT2
      IF (IANTH.GE.0) ECM=ELAB+AMT+0.00001D0
      AM11=AM1
      AM22=AM2
      IF (IT2.GT.0)                                            GO TO 120
C
C-----------------------------
C  INCLUSION OF DIRECT RESONANCES
C  RANDOM CHOICE OF DECAY CHANNELS OF THE DIRECT RESONANCE  IT1
C------------------------
      KZ1=K1H(IT1)
      IST=IST+1
      IECO=0
      ECO=ECM
      GAM=(ELAB+AMT)/ECO
      BGAM=PLAB/ECO
      CXS(1)=CX
      CYS(1)=CY
      CZS(1)=CZ
                                                               GO TO 170
  120 CONTINUE
      CALL DRANDM(WW)
      IF(WW.LT. 0.5D0)                                         GO TO 130
      IT1=IT22
      IT2=IT11
      AM1=AM22
      AM2=AM11
  130 CONTINUE
C
C-----------------------------
C   THE FIRST PARTICLE IS DEFINED TO BE THE FORWARD GOING ONE AT SMALL T
      IBN=IBARH(N)
      IB1=IBARH(IT1)
      IT11=IT1
      IT22=IT2
      AM11=AM1
      AM22=AM2
      IF(IB1.EQ.IBN)                                           GO TO 140
      IT1=IT22
      IT2=IT11
      AM1=AM22
      AM2=AM11
  140 CONTINUE
C-----------------------------
C***IT1,IT2 ARE THE CREATED PARTICLES
C***MOMENTA AND DIRECTION COSINA IN THE CM - SYSTEM
C------------------------
      CALL DTWOPA(ECM1,ECM2,PCM1,PCM2,COD1,COD2,COF1,COF2,SIF1,SIF2,
     *IT1,IT2,ECM,ECMN,PCMN,N,AM1,AM2)
      IST=IST+1
      ITS(IST)=IT1
      AMM(IST)=AM1
C
C-----------------------------
C***TRANSFORMATION INTO LAB SYSTEM AND ROTATION
C----------------------------
      CALL DTRAFO(GAM,BGAM,CX,CY,CZ,COD1,COF1,SIF1,PCM1,ECM1,PLS(IST),
     *CXS(IST),CYS(IST),CZS(IST),ELS(IST))
      IST=IST+1
      ITS(IST)=IT2
      AMM(IST)=AM2
      CALL DTRAFO(GAM,BGAM,CX,CY,CZ,COD2,COF2,SIF2,
     *PCM2,ECM2,PLS(IST),CXS(IST),CYS(IST),CZS(IST),ELS(IST))
  150 CONTINUE
C
C-----------------------------
C***TEST   STABLE OR UNSTABLE
C----------------------------
      IF(ITS(IST).GT.NSTAB)                                    GO TO 160
      IRH=IRH+1
C
C-----------------------------
C***IRH IS THE NUMBER OF THE FINAL STABLE PARTICLE
C----------------------------
C*    IF (REDU.LT.0.D0) GO TO 1009
      ITRH(IRH)=ITS(IST)
      PLRH(IRH)=PLS(IST)
      CXRH(IRH)=CXS(IST)
      CYRH(IRH)=CYS(IST)
      CZRH(IRH)=CZS(IST)
      ELRH(IRH)=ELS(IST)
      IST=IST-1
      IF(IST.GE.1)                                             GO TO 150
                                                               GO TO 260
  160 CONTINUE
C
C  RANDOM CHOICE OF DECAY CHANNELS
C----------------------------
C
      IT=ITS(IST)
      ECO=AMM(IST)
      GAM=ELS(IST)/ECO
      BGAM=PLS(IST)/ECO
      IECO=0
      KZ1=K1H(IT)
  170 CONTINUE
      IECO=IECO+1
      CALL DRANDM(VV)
      VV=VV-1.D-17
      IIK=KZ1-1
  180 IIK=IIK+1
      IF (VV.GT.WT(IIK))                                       GO TO 180
C
C  IIK IS THE DECAY CHANNEL
C----------------------------
      IT1=NZK(IIK,1)
      I310=0
  190 CONTINUE
      I310=I310+1
      AM1=DAMG(IT1)
      IT2=NZK(IIK,2)
      AM2=DAMG(IT2)
      IF (IT2-1.LT.0)                                          GO TO 240
      IT3=NZK(IIK,3)
      AM3=DAMG(IT3)
      AMS=AM1+AM2+AM3
C
C  IF  IIK-KIN.LIM.GT.ACTUAL TOTAL CM-ENERGY, DO AGAIN RANDOM IIK-CHOICE
C----------------------------
      IF (IECO.LE.10)                                          GO TO 200
      IATMPT=IATMPT+1
      IF(IATMPT.GT.3)                                          GO TO 280
                                                                GO TO 40
  200 CONTINUE
      IF (I310.GT.50)                                          GO TO 170
      IF (AMS.GT.ECO)                                          GO TO 190
C
C  FOR THE DECAY CHANNEL
C  IT1,IT2, IT3 ARE THE PRODUCED PARTICLES FROM  IT
C----------------------------
      IF (REDU.LT.0.D0)                                        GO TO 30
      ITWTHC=0
      REDU=2.
      IF(IT3.EQ.0)                                             GO TO 220
  210 CONTINUE
      ITWTH=1
      CALL DTHREP(ECO,ECM1,ECM2,ECM3,PCM1,PCM2,PCM3,COD1,COF1,SIF1,
     *COD2,COF2,SIF2,COD3,COF3,SIF3,AM1,AM2,AM3)
                                                               GO TO 230
  220 CALL DTWOPD(ECO,ECM1,ECM2,PCM1,PCM2,COD1,COF1,SIF1,COD2,COF2,SIF2,
     +AM1,AM2)
      ITWTH=-1
      IT3=0
  230 CONTINUE
      ITWTHC=ITWTHC+1
      IF (REDU.GT.0.D0)                                        GO TO 240
      REDU=2.
      IF (ITWTHC.GT.100)                                        GO TO 30
      IF (ITWTH) 220,220,210
  240 CONTINUE
      ITS(IST  )=IT1
      IF (IT2-1.LT.0)                                          GO TO 250
      ITS(IST+1)  =IT2
      ITS(IST+2)=IT3
      RX=CXS(IST)
      RY=CYS(IST)
      RZ=CZS(IST)
      AMM(IST)=AM1
      CALL DTRAFO(GAM,BGAM,RX,RY,RZ,COD1,COF1,SIF1,PCM1,ECM1,
     *PLS(IST),CXS(IST),CYS(IST),CZS(IST),ELS(IST))
      IST=IST+1
      AMM(IST)=AM2
      CALL DTRAFO(GAM,BGAM,RX,RY,RZ,COD2,COF2,SIF2,PCM2,ECM2,
     *PLS(IST),CXS(IST),CYS(IST),CZS(IST),ELS(IST))
      IF (IT3.LE.0)                                            GO TO 250
      IST=IST+1
      AMM(IST)=AM3
      CALL DTRAFO(GAM,BGAM,RX,RY,RZ,COD3,COF3,SIF3,PCM3,ECM3,
     *PLS(IST),CXS(IST),CYS(IST),CZS(IST),ELS(IST))
  250 CONTINUE
                                                               GO TO 150
  260 CONTINUE
  270 CONTINUE
      RETURN
  280 CONTINUE
C
C----------------------------
C
C   ZERO CROSS SECTION CASE
C----------------------------
C
      IRH=1
      ITRH(1)=N
      CXRH(1)=CX
      CYRH(1)=CY
      CZRH(1)=CZ
      ELRH(1)=ELAB
      PLRH(1)=PLAB
      RETURN
      END
*-- Author :
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      BLOCK DATA  RUNTT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/DRUN/RUNT(2)
      DATA RUNT/100.D0,100.D0/
      END
*-- Author :
      DOUBLE PRECISION FUNCTION REXP(W)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C*****EXPONENTIELL VERTEILTE ZUFALLSZAHL, VERTEILUNG F=EXP(-W*X)
      A=-1.
   10 CALL DRANDM(XO)
      A=A+1.
      B=0.
      N=0
      CALL DRANDM(V)
   20 B=B+V
      N=N+1
      IF(B-XO) 20,20,30
   30 NN=N/2
      IF(N.EQ.2*NN)                                             GO TO 10
      REXP=(A+XO)/W
      RETURN
      END
*-- Author :
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      BLOCK DATA NONAME
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C*****  BLOCK DATA WITHOUT NAMES
C     INTEGER * 2 IEII,IKII
      COMMON /DSLOPE/ AMBMBB(75)
      COMMON /DREDVE/ THRESH(268), IRII(17),IKII(17),IEII(17)
C     DATAS     DATAS    DATAS      DATAS     DATAS
C******          *********
      DATA IKII/ 0, 15, 41, 67, 82, 93, 110, 133, 148, 159, 172, 183,
     &           207, 224, 241, 252, 268 /
      DATA IEII/ 0, 21, 46, 71, 92, 109, 126, 143, 160, 173, 186, 199,
     &           220, 241, 262, 279, 296 /
      DATA IRII/ 0, 315, 965, 1615, 1930, 2117, 2406, 2797, 3052, 3195,
     &           3364, 3507, 4011, 4368, 4725, 4912, 5184/
 
C
C     MASSES FOR THE SLOPE B(M) IN GEV
C     SLOPE B(M) FOR AN MESONIC SYSTEM
C     SLOPE B(M) FOR A BARYONIC SYSTEM

*
      DATA AMBMBB/  0.8D0, 0.85D0,  0.9D0, 0.95D0, 1.D0,
     &     1.05D0,  1.1D0, 1.15D0,  1.2D0, 1.25D0,
     &      1.3D0,  1.35D0, 1.4D0,  1.45D0,  1.5D0,
     &     1.55D0,  1.6D0,  1.65D0, 1.7D0,   1.75D0,
     &      1.8D0,  1.85D0, 1.9D0,  1.95D0,  2.D0,
     &     15.6D0, 14.95D0, 14.3D0, 13.65D0, 13.D0,
     &    12.35D0, 11.7D0, 10.85D0, 10.D0,  9.15D0,
     &      8.3D0,  7.8D0,  7.3D0,  7.25D0,  7.2D0,
     &     6.95D0,  6.7D0,  6.6D0,  6.5D0,   6.3D0,
     &      6.1D0,  5.85D0, 5.6D0,  5.35D0,  5.1D0,
     &      15.D0,   15.D0, 15.D0,  15.D0,   15.D0, 15.D0, 15.D0,
     &     14.2D0,  13.4D0, 12.6D0,
     &     11.8D0, 11.2D0, 10.6D0,  9.8D0,    9.D0,
     &     8.25D0,  7.5D0, 6.25D0,  5.D0,    4.5D0, 5*4.D0 /
*
      END
*-- Author :
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DOUBLE PRECISION FUNCTION DAMG(IT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C*** RANDOM SELECTION OF MASSES OF DECAYING PARTICLES
C     INTEGER * 2 ICH,IBAR,K1,K2,NZK,NRK
      COMMON/DABLTI/AMM(110),GA(110),TAU(110),ICH(110)
     *,IBAR(110),K1(110),K2(110)
      COMMON /DGAMRE/ REDU,AMO,AM (15 )
      DIMENSION GASUNI(14)
      DATA GASUNI/
     *-1.D0,-.98D0,-.95D0,-.87D0,-.72D0,-.48D0,
     *-.17D0,.17D0,.48D0,.72D0,.87D0,.95D0,.98D0,1.D0/
      DATA GAUNO/2.352D0/
      DATA GAUNON/2.4D0/
      DATA IO/14/
      DATA NSTAB/23/
      I=1
      IF (IT.LE.0)                                              GO TO 30
      IF (IT.LE.NSTAB)                                          GO TO 20
      DGAUNI=GAUNO*GAUNON/(IO-1.)
      CALL DRANDM(VV)
      VV=VV*2.-1.+1.D-16
   10 CONTINUE
      VO=GASUNI(I)
      I=I+1
      V1=GASUNI(I)
      IF (VV.GT.V1)                                             GO TO 10
      UNIGA=DGAUNI*(I-2.+(VV-VO+1.E-16)/(V1-VO)-(IO-1.)*.5)
      DAM=GA(IT)*UNIGA/GAUNO
      AAM=AMM(IT)+DAM
      DAMG=AAM
      RETURN
   20 CONTINUE
      DAMG=AMM(IT)
      RETURN
   30 CONTINUE
      DAMG=0.
      RETURN
      END
*-- Author :
      SUBROUTINE DCALUM(N,ITTA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C*** C.M.S.-ENERGY AND REACTION CHANNEL THRESHOLD CALCULATION
C     INTEGER * 2 ICH,IBAR,K1,K2,NZK,NRK
C    *,IEII,IKII,NURE
      COMMON/DABLTI/AM(110),GA(110),TAU(110),ICH(110)
     *,IBAR(110),K1(110),K2(110)
      COMMON /DREDVE/ THRESH(268), IRII(17),IKII(17),IEII(17)
      COMMON/DSPLI/NZK(460,3),WT(460)
      COMMON /DREAC/UMO( 296),PLABF( 296),SIIN( 296),WK( 5184),
     *NRK(2, 268),NURE(30,2)
      IRE=NURE(N,ITTA/8+1)
      IEO=IEII(IRE)+1
      IEE=IEII(IRE +1)
      AM1=AM(N   )
      AM12=AM1**2
      AM2=AM(ITTA)
      AM22=AM2**2
      DO 10 IE=IEO,IEE
        PLAB2=PLABF(IE)**2
        ELAB=SQRT(AM12+AM22+2.*SQRT(PLAB2+AM12)*AM2)
        UMO(IE)=ELAB
   10 CONTINUE
      IKO=IKII(IRE)+1
      IKE=IKII(IRE +1)
      UMOO=UMO(IEO)
      DO 30 IK=IKO,IKE
        IF(NRK(2,IK).GT.0)                                      GO TO 30
        IKI=NRK(1,IK)
        AMSS=5.
        K11=K1(IKI)
        K22=K2(IKI)
        DO 20 IK1=K11,K22
          IN=NZK(IK1,1)
          AMS=AM(IN)
          IN=NZK(IK1,2)
          IF(IN.GT.0)AMS=AMS+AM(IN)
          IN=NZK(IK1,3)
          IF(IN.GT.0) AMS=AMS+AM(IN)
          IF (AMS.LT.AMSS) AMSS=AMS
   20   CONTINUE
        IF(UMOO.LT.AMSS) UMOO=AMSS
        THRESH(IK)=UMOO
   30 CONTINUE
      RETURN
      END
*-- Author :
      SUBROUTINE DCHANH
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/DABLTI/AM(110),GA(110),TAU(110),ICH(110)
     *,IBAR(110),K1(110),K2(110)
      COMMON/DSPLI/ NZK(460,3),WT(460)
      DIMENSION HWT(460)
      COMMON /DREDVE/ THRESH(268), IRII(17),IKII(17),IEII(17)
      DIMENSION HWK(40)
      COMMON /DREAC/UMO( 296),PLABF( 296),SIIN( 296),WK( 5184),
     *NRK(2, 268),NURE(30,2)
      DIMENSION SI(5184)
      EQUIVALENCE (WK(1),SI(1))
C--------------------
C*** USE ONLY FOR DATAPREPARATION OF PURE HADRIN
C*** CALCULATION OF REACTION- AND DECAY-CHANNEL-WEIGHTS,
C*** THRESHOLD ENERGIES+MOMENTA OF REACTION CHNLS.
C*** CHANGE OF WT- AND WK-INPUTDATA INTO WEIGHTS FOR THE M.-C.-PROCEDURE
C*** (ADDED ONE TO EACH OTHER FOR CORRESPONDING CHANNELS)
C--------------------------
      IREG=16
      DO 90 IRE=1,IREG
        IWKO=IRII(IRE)
        IEE=IEII(IRE+1)-IEII(IRE)
        IKE=IKII(IRE+1)-IKII(IRE)
        IEO=IEII(IRE)+1
        IIKA=IKII(IRE)
*   modifications to suppress elestic scattering  24/07/91
        DO 80 IE=1,IEE
          SIS=1.E-14
          SINORC=0.0
          DO 10 IK=1,IKE
            IWK=IWKO+IEE*(IK-1)+IE
            IF(NRK(2,IIKA+IK).EQ.0) SINORC=1.0
            SIS=SIS+SI(IWK)*SINORC
   10     CONTINUE
          SIIN(IEO+IE-1)=SIS
          SIO=0.
          IF (SIS.GE.1.D-12)                                    GO TO 20
          SIS=1.
          SIO=1.
   20     CONTINUE
          SINORC=0.0
          DO 30 IK=1,IKE
            IWK=IWKO+IEE*(IK-1)+IE
            IF(NRK(2,IIKA+IK).EQ.0) SINORC=1.0
            SIO=SIO+SI(IWK)*SINORC/SIS
            HWK(IK)=SIO
   30     CONTINUE
          DO 40 IK=1,IKE
            IWK=IWKO+IEE*(IK-1)+IE
   40     WK(IWK)=HWK(IK)
          IIKI=IKII(IRE)
          DO 70 IK=1,IKE
            AM111=0.
            INRK1=NRK(1,IIKI+IK)
            IF (INRK1.GT.0) AM111=AM(INRK1)
            AM222=0.
            INRK2=NRK(2,IIKI+IK)
            IF (INRK2.GT.0) AM222=AM(INRK2)
            THRESH(IIKI+IK)=AM111 +AM222
            IF (INRK2-1.GE.0)                                   GO TO 60
            INRKK=K1(INRK1)
            AMSS=5.
            INRKO=K2(INRK1)
            DO 50 INRK1=INRKK,INRKO
              INZK1=NZK(INRK1,1)
              INZK2=NZK(INRK1,2)
              INZK3=NZK(INRK1,3)
              IF (INZK1.LE.0.OR.INZK1.GT.110)                   GO TO 50
              IF (INZK2.LE.0.OR.INZK2.GT.110)                   GO TO 50
              IF (INZK3.LE.0.OR.INZK3.GT.110)                   GO TO 50
C     WRITE (6,310)INRK1,INZK1,INZK2,INZK3
 1000 FORMAT (4I10)
              AMS=AM(INZK1)+AM(INZK2)
              IF (INZK3-1.GE.0) AMS=AMS+AM(INZK3)
              IF (AMSS.GT.AMS) AMSS=AMS
   50       CONTINUE
            AMS=AMSS
            IF (AMS.LT.UMO(IEO)) AMS=UMO(IEO)
            THRESH(IIKI+IK)=AMS
   60       CONTINUE
   70     CONTINUE
   80   CONTINUE
   90 CONTINUE
      DO 100 J=1,460
  100 HWT(J)=0.
      DO 120 I=1,110
        IK1=K1(I)
        IK2=K2(I)
        HV=0.
        IF (IK2.GT.460)IK2=460
        IF (IK1.LE.0)IK1=1
        DO 110 J=IK1,IK2
          HV=HV+WT(J)
          HWT(J)=HV
          JI=J
  110   CONTINUE
C       IF (ABS(HV-1.).GT.1.E-4) WRITE(6,1010)I,JI,HV
C1010 FORMAT (35H ERROR IN HWT, FALSE USE OF CHANWH ,2I6,F10.2)
  120 CONTINUE
      DO 130 J=1,460
  130 WT(J)=HWT(J)
      RETURN
      END
*-- Author :
      SUBROUTINE DHADDE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 ANAMC,ZKNAMC,ANAME,ZKNAME,RKNAME,ANAMZ,ZKNAMZ
C
      COMMON/DPAR/ANAMC(210),AMC(210),GAC(210),TAUC(210),ICHC(210),
     *IBARC(210),K1C(210),K2C(210)
C      COMMON /DDECAC/ ZKNAMC(540),NZKC(540,3),WTC(540)
  
      PARAMETER (IDMAX9=602)
C      CHARACTER*8 ZKNAME
      COMMON/DDECAC/ ZKNAMC(IDMAX9),WTC(IDMAX9),NZKC(IDMAX9,3)


      COMMON /DNAMS/ ANAME(110),ZKNAME(460),RKNAME(268)
      COMMON /DABLTI/AM(110),GA(110),TAU(110),ICH(110)
     *,IBAR(110),K1(110),K2(110)
      COMMON /DSPLI/NZK(460,3),WT(460)
      COMMON /DADDHP/ AMZ(16),GAZ(16),TAUZ(16),ICHZ(16),IBARZ(16),K1Z
     +(16),K2Z(16),WTZ(153),II22, NZKZ(153,3)
 
      COMMON/DADDHN/ANAMZ(16),ZKNAMZ(153)
      DATA IRETUR/0/
      IRETUR=IRETUR+1
      AM(31)=0.48
      IF (IRETUR.GT.1) RETURN
      DO 10 I=1,94
        ANAME(I)=ANAM C(I)
        AM (I)=AMC(I)
        GA( I)=GA C(I)
        TAU( I)=TAU C(I)
        ICH( I)=ICH C(I)
        IBAR( I)=IBARC (I)
        K1( I)=K1C(I)
        K2( I)=K2 C(I)
   10 CONTINUE
      AM(1)=0.9383D0
      AM(2)=AM(1)
      DO 20 I=26,30
        K1(I)=452
        K2(I)=452
   20 CONTINUE
      DO 30 I=1,307
        ZKNAME(I)=ZKNAMC(I)
        WT( I)=WT C(I)
        NZK( I,1)=NZK C(I, 1)
        NZK( I,2)=NZK C(I, 2)
        NZK( I,3)=NZK C(I, 3)
   30 CONTINUE
      DO 40 I=1,16
        L=I+94
        ANAME(L)=ANAMZ (I)
        AM (L)=AMZ(I)
        GA( L)=GA Z(I)
        TAU( L)=TAU Z(I)
        ICH( L)=ICH Z(I)
        IBAR( L)=IBARZ (I)
        K1( L)=K1Z(I)
        K2( L)=K2 Z(I)
   40 CONTINUE
      DO 50 I=1,153
        L=I+307
        ZKNAME(L)=ZKNAMZ(I)
        WT( L)=WT Z(I)
        NZK( L,3)=NZK Z(I, 3)
        NZK( L,2)=NZK Z(I, 2)
        NZK( L,1)=NZK Z(I, 1)
   50 CONTINUE
      RETURN
      END
*-- Author :
      FUNCTION IEFUND(PL,IRE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C*****IEFUN CALCULATES A MOMENTUM INDEX
C     INTEGER * 2 ICH,IBAR,K1,K2,NZK,NRK
C    *,IEII,IKII,NURE
      COMMON /DRUN/ RUNTES,EFTES
      COMMON /DREDVE/ THRESH(268), IRII(17),IKII(17),IEII(17)
      COMMON /DREAC/UMO( 296),PLABF( 296),SIIN( 296),WK( 5184),
     *NRK(2, 268),NURE(30,2)
      IPLA=IEII(IRE)+1
     *+1
      IPLE=IEII(IRE+1)
      IF (PL.LT.0.)                                             GO TO 30
      DO 10 I=IPLA,IPLE
        J=I-IPLA+1
        IF (PL.LE.PLABF(I))                                     GO TO 60
   10 CONTINUE
      I=IPLE
      IF ( EFTES.GT.40.D0)                                      GO TO 20
      EFTES=EFTES+1.
      WRITE(6,1000)PL,J
   20 CONTINUE
                                                                GO TO 70
   30 CONTINUE
      DO 40 I=IPLA,IPLE
        J=I-IPLA+1
        IF (-PL.LE.UMO(I))                                      GO TO 60
   40 CONTINUE
      I=IPLE
      IF ( EFTES.GT.40.D0)                                      GO TO 50
      EFTES=EFTES+1.
      WRITE(6,1000)PL,I
   50 CONTINUE
   60 CONTINUE
   70 CONTINUE
      IEFUND=I
      RETURN
 1000 FORMAT(14H PLAB OR -ECM=,E12.4,27H IS OUT OF CONSIDERED RANGE ,
     +7H IEFUN=,I5)
      END
*-- Author :
      SUBROUTINE DSIGIN(IRE ,PLAB,N,IE ,AMT ,AMN,ECM ,SI ,ITAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/DABLTI/AM(110),GA(110),TAU(110),ICH(110)
     *,IBAR(110),K1(110),K2(110)
      COMMON /DREDVE/ THRESH(268), IRII(17),IKII(17),IEII(17)
      COMMON /DREAC/UMO( 296),PLABF( 296),SIIN( 296),WK( 5184),
     *NRK(2, 268),NURE(30,2)
      IE=IEFUND(PLAB,IRE)
      IF (IE.LE.IEII(IRE)) IE=IE+1
      AMT=AM(ITAR)
      AMN=AM(N)
      AMN2=AMN*AMN
      AMT2=AMT*AMT
      ECM=SQRT(AMN2+AMT2+2.*AMT*SQRT(AMN2+PLAB**2))
C*** INTERPOLATION PREPARATION
      ECMO=UMO(IE)
      ECM1=UMO(IE-1)
      DECM=ECMO-ECM1
      DEC=ECMO-ECM
      IIKI=IKII(IRE)+1
      EKLIM=-THRESH(IIKI)
      WOK=SIIN(IE)
      WDK=WOK-SIIN(IE-1)
      IF (ECM.GT.ECMO) WDK=0.
C*** INTERPOLATION IN CHANNEL WEIGHTS
      IELIM=IEFUND(EKLIM,IRE)
      DELIM=UMO(IELIM)+EKLIM
     *+1.E-16
      DETE=(ECM-(ECMO-EKLIM)*.5)*2.
      IF (DELIM*DELIM-DETE*DETE) 20,20,10
   10 DECC=DELIM
                                                                GO TO 30
   20 DECC=DECM
   30 CONTINUE
      WKK=WOK-WDK*DEC/(DECC+1.D-9)
      IF (WKK.LT.0.) WKK=0.
      SI=WKK+1.D-12
      IF (-EKLIM.GT.ECM) SI=1.D-14
      RETURN
      END
*-- Author :
      SUBROUTINE DTCHOI(T,P,PP,E,EE,I,II,N,AM1,AM2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ****************************
C     TCHOIC CALCULATES A RANDOM VALUE
C     FOR THE FOUR-MOMENTUM-TRANSFER T
C     ****************************
      COMMON /DABLTI/ AM(110),GA(110),TAU(110),ICH(110), IBAR(110),K1
     +(110),K2(110)
      COMMON /DSLOPE/ SM(25),BBM(25),BBB(25)
      AMA=AM1
      AMB=AM2
      IF (I.GT.30.AND.II.GT.30)                                 GO TO 20
      III=II
      AM3=AM2
      IF (I.LE.30)                                              GO TO 10
      III=I
      AM3=AM1
   10 CONTINUE
                                                                GO TO 30
   20 CONTINUE
      III=II
      AM3=AM2
      IF (AMA.LE.AMB)                                           GO TO 30
      III=I
      AM3=AM1
   30 CONTINUE
      IB=IBAR(III)
      AMA=AM3
      K=(AMA-0.75)/0.05
      IF (K-2.LT.0) K=1
      IF (K-26.GE.0) K=25
      IF (IB)50,40,50
   40 BM=BBM(K)
                                                                GO TO 60
   50 BM=BBB(K)
   60 CONTINUE
C     NORMALIZATION
      TMIN=-2.*(E*EE-P*PP)+AM(N)**2+AM1  **2
      TMAX=-2.*(E*EE+P*PP)+AM(N)**2+AM1  **2
      CALL DRANDM(VB)
**sr 19-11-95  
C     IF (VB.LT.0.2D0) BM=BM*0.1
C    **0.5
      BM = BM*5.0D0
**
      TMI=BM*TMIN
      TMA=BM*TMAX
      ETMA=0.
      IF (ABS(TMA).GT.120.D0)                                   GO TO 70
      ETMA=EXP(TMA)
   70 CONTINUE
      AN=(1./BM)*(EXP(TMI)-ETMA)
C*** RANDOM CHOICE OF THE T - VALUE
      CALL DRANDM(R)
      T=(1./BM)*LOG(ETMA+R*AN*BM)
      RETURN
      END
*-- Author :
      SUBROUTINE DTWOPA(E1,E2,P1,P2,COD1,COD2,COF1,COF2,SIF1,SIF2,
     1IT1,IT2,UMOO,ECM,P,N,AM1,AM2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ******************************************************
C     QUASI TWO PARTICLE PRODUCTION
C     TWOPAR CALCULATES THE ENERGYS AND THE MOMENTA
C     FOR THE CREATED PARTICLES OR RESONANCES IT1 AND IT2
C     IN THE CM - SYSTEM
C     COD1,COD2,COF1,COF2,SIF1,SIF2 ARE THE ANGLES FOR
C     SPHERICAL COORDINATES
C     ******************************************************
      COMMON /DABLTI/ AM(110),GA(110),TAU(110),ICH(110), IBAR(110),K1
     +(110),K2(110)
C
      AMA=AM1
      AMB=AM2
      AMA2=AMA*AMA
      E1=((UMOO-AMB)*(UMOO+AMB) + AMA2)/(2.*UMOO)
      E2=UMOO - E1
      IF (E1.LT.AMA*1.00001D0) E1=AMA*1.00001D0
      AMTE=(E1-AMA)*(E1+AMA)
      AMTE=AMTE+1.D-18
      P1=SQRT(AMTE)
      P2=P1
C     / P2 / = / P1 /  BUT OPPOSITE DIRECTIONS
C     DETERMINATION  OF  THE ANGLES
C     COS(THETA1)=COD1      COS(THETA2)=COD2
C     SIN(PHI1)=SIF1        SIN(PHI2)=SIF2
C     COS(PHI1)=COF1        COS(PHI2)=COF2
C     PHI IS UNIFORMLY DISTRIBUTED IN ( 0,2*PI )
      CALL DCOSI(COF1,SIF1)
      COF2=-COF1
      SIF2=-SIF1
C     CALCULATION OF THETA1
      CALL DTCHOI(TR,P,P1,ECM,E1,IT1,IT2,N,AM1,AM2)
      COD1=(TR-AMA2-AM(N)*AM(N)+2.*ECM*E1)/(2.*P*P1+1.E-18)
      IF (COD1.GT.0.9999999D0) COD1=0.9999999D0
      COD2=-COD1
      RETURN
      END
*-- Author :
      SUBROUTINE DCOSI(SFE,CFE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
   10 X=RNDM(A)
      Y=RNDM(B)
      XX=X*X
      YY=Y*Y
      XY=XX+YY
      IF(XY.GT.1.)                                              GO TO 10
      CFE=(XX-YY)/XY
      SFE=2.*X*Y/XY
      IF(RNDM(C).LT.0.5D0)                                      GO TO 20
      RETURN
   20 SFE=-SFE
      RETURN
      END
*-- Author :
      SUBROUTINE DGAUSS(X,A,S)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA IS/0/
      IF(IS.NE.0)                                               GO TO 10
      CALL DRANDM(X)
      RO=SQRT(ABS(2.*LOG(X)))
      CALL DCOSI(SFE,CFE)
      X=RO*SFE
      IS=1
                                                                GO TO 20
   10 X=RO*CFE
      IS=0
   20 X=A+X*S
      RETURN
      END
*-- Author :
      SUBROUTINE DRANDM(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      X=RNDM(V)
      RETURN
      END
*-- Author :
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

*=== blkdt5 ===========================================================*
*==                                                                    *
      BLOCK DATA ZK
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 ZKNAM5,ZKNAM6,ANAMZ,ZKNAM4
      COMMON /DADDHP/ AMZ(16),GAZ(16),TAUZ(16),ICHZ(16),IBARZ(16),
     &       K1Z(16),K2Z(16),WTZ(153),II22,NZK1(153),NZK2(153),NZK3(153)
      COMMON /DADDHN/ ANAMZ(16),ZKNAM4(9),ZKNAM5(90),ZKNAM6(54)
*     Particle masses in GeV                                           *
      DATA AMZ/ 3*2.2D0, 0.9576D0, 3*1.887D0, 2.4D0, 2.03D0, 2*1.44D0,
     &          2*1.7D0, 3*0.D0/
*     Resonance width Gamma in GeV                                     *
      DATA GAZ/ 3*.2D0, .1D0, 4*.2D0, .18D0, 2*.2D0, 2*.15D0, 3*0.D0 /
*     Mean life time in seconds                                        *
      DATA TAUZ / 16*0.D0 /
*     Charge of particles and resonances                               *
      DATA ICHZ/ 0, 1, 3*0, 1, -1, 0, 1, -1, 0, 0, 1 , 3*0 /
*     Baryonic charge                                                  *
      DATA IBARZ/ 2, 7*0, 1, -1, -1, 1, 1, 3*0 /
*     First number of decay channels used for resonances               *
*     and decaying particles                                           *
      DATA K1Z/ 308,310,313,317,322,365,393,421,425,434,440,446,449,
     &          3*460/
*     Last number of decay channels used for resonances                *
*     and decaying particles                                           *
      DATA K2Z/ 309,312,316,321,364,392,420,424,433,439,445,448,451,
     &          3*460/
*     Weight of decay channel                                          *
      DATA WTZ/ .17D0, .83D0, 2*.33D0, .34D0, .17D0, 2*.33D0, .17D0,
     & .01D0, .13D0, .36D0, .27D0, .23D0, .0014D0, .0029D0, .0014D0,
     & .0029D0, 4*.0007D0, .0517D0, .0718D0, .0144D0, .0431D0, .0359D0,
     & .0718D0, .0014D0, .0273D0, .0014D0, .0431D0, 2*.0129D0, .0259D0,
     & .0517D0, .0359D0, .0014D0, 2*.0144D0, .0129D0, .0014D0, .0259D0,
     & .0359D0, .0072D0, .0474D0, .0948D0, .0259D0, .0072D0, .0144D0,
     & .0287D0, .0431D0, .0144D0, .0287D0, .0474D0, .0144D0, .0075D0,
     & .0057D0, .0019D0, .0038D0, .0095D0, 2*.0014D0, .0191D0, .0572D0,
     & .1430D0, 2*.0029D0, 5*.0477D0, .0019D0, .0191D0, .0686D0,.0172D0,
     & .0095D0, .1888D0, .0172D0, .0191D0, .0381D0, 2*.0571D0, .0190D0,
     & .0057D0, .0019D0, .0038D0, .0095D0, .0014D0, .0014D0, .0191D0,
     & .0572D0, .1430D0, 2*.0029D0, 5*.0477D0, .0019D0, .0191D0,.0686D0,
     & .0172D0, .0095D0, .1888D0, .0172D0, .0191D0, .0381D0, 2*.0571D0,
     & .0190D0, 4*.25D0, 2*.2D0, .12D0, .1D0, .07D0, .07D0, .14D0,
     & 2*.05D0, .0D0, .3334D0, .2083D0, 2*.125D0, .2083D0, .0D0, .125D0,
     & .2083D0, .3334D0, .2083D0, .125D0, .3D0, .05D0, .65D0, .3D0,
     & .05D0, .65D0, 9*1.D0 /
*     Particle numbers in decay channel                                *
      DATA NZK1/ 8, 1, 2, 9, 1, 2, 9, 2, 9, 7, 13, 31, 15, 24, 23, 13,
     & 23, 13, 2*23, 14, 13, 23, 31, 98, 2*33, 32, 23, 14, 13, 35, 2*23,
     & 14, 13, 33, 23, 98, 31, 23, 14, 13, 35, 2*33, 32, 23, 35, 33, 32,
     & 98, 5*35, 4*13, 23, 13, 98, 32, 33, 23, 13, 23, 13, 14, 13, 32,
     & 13, 98, 23, 13, 2*32, 13, 33, 32, 98, 2*35, 4*14, 23, 14, 98,
     & 2*34, 23, 14, 23, 2*14, 13, 34, 14, 98, 23, 14, 2*34, 14, 33, 32,
     & 98, 2*35, 104, 61, 105, 62, 1, 17, 21, 17, 22, 2*21, 22, 21, 2,
     & 67, 68, 69, 2, 2*9, 68, 69, 70, 2, 9, 2*24, 15, 2*25, 16, 9*0/
      DATA NZK2/ 2*8, 1, 8, 9, 2*8, 2*1, 7, 14, 13, 16, 25, 23, 14, 23,
     & 14, 31, 33, 32, 34, 35, 31, 23, 31, 33, 34, 31, 32, 34, 31, 33,
     & 32, 2*33, 35, 31, 33, 31, 33, 32, 34, 35, 31, 33, 34, 35, 31,
     & 4*33, 32, 3*35,  2*23, 13, 31, 32, 33, 13, 31, 32, 2*31, 32, 33,
     & 32, 32, 35, 31, 2*32, 33, 31, 33, 35, 33, 3*32, 35, 2*23, 14,
     & 31, 34, 33, 14, 31, 33, 2*31, 34, 32, 33, 34, 35, 31, 2*34, 33,
     & 31, 33, 35, 33, 2*34, 33, 35, 1, 2, 8, 9, 25, 13, 35, 2*32, 33,
     & 31, 13, 23, 31, 13, 23, 14, 79, 80, 31, 13, 23, 14, 78, 79, 8,
     & 1, 8, 1, 8, 1, 9*0 /
      DATA NZK3/ 23, 14, 2*13, 23, 13, 2*23, 14, 0, 7, 14, 4*0, 2*23,
     & 10*0, 33, 2*31, 0, 33, 34, 32, 34, 0, 35, 0, 31, 3*35, 0, 3*31,
     & 35, 31, 33, 34, 31, 33, 34, 31, 33, 35, 0, 23, 14, 6*0, 32, 3*33,
     & 32, 34, 0, 35, 0, 2*35, 2*31, 35, 32, 34, 31, 33, 32, 0, 23, 13,
     & 6*0, 34, 2*33, 34, 33, 34, 0, 35, 0,2*35, 2*31, 35, 2*34, 31,
     & 2*34, 25*0, 23, 2*14, 23, 2*13, 9*0 /
*     Particle  names                                                  *
      DATA ANAMZ / 'NNPI', 'ANPPI', 'ANNPI', ' ETS  ',' PAP  ',' PAN  ',
     & 'APN', 'DEO   ', 'S+2030', 'AN*-14', 'AN*014','KONPI ','AKOPPI',
     & 3*'BLANK' /
*     Name of decay channel                                            *
      DATA ZKNAM4/'NNPI0','PNPI-','APPPI+','ANNPI+','ANPPI0','APNPI+',
     & 'ANNPI0','APPPI0','ANPPI-'/
      DATA ZKNAM5/' GAGA ','P+P-GA','ETP+P-','K+K-  ','K0AK0 ',
     & ' POPO ',' P+P- ','POPOPO','P+P0P-','P0ET  ','&0R0  ','P-R+  ',
     & 'P+R-  ','POOM  ',' ETET ','ETSP0 ','R0ET  ',' R0R0 ','R+R-  ',
     & 'P0ETR0','P-ETR+','P+ETR-',' OMET ','P0R0R0','P0R+R-','P-R+R0',
     & 'P+R-R0','R0OM  ','P0ETOM','ETSR0 ','ETETET','P0R0OM','P-R+OM',
     & 'P+R-OM','OMOM  ','R0ETET','R0R0ET','R+R-ET','P0OMOM','OMETET',
     & 'R0R0R0','R+R0R-','ETSRET','OMR0R0','OMR+R-','OMOMET','OMOMR0',
     & 'OMOMOM',
     & ' P+PO ','P+POPO','P+P+P-','P+ET  ','P0R+  ','P+R0  ','ETSP+ ',
     & 'R+ET  ',' R0R+ ','POETR+','P+ETR0','POR+R-','P+R0R0','P-R+R+',
     & 'P+R-R+','R+OM  ','P+ETOM','ETSR+ ','POR+OM','P+R0OM','R+ETET',
     & 'R+R0ET','P+OMOM','R0R0R+','R+R+R-','ETSR+E','OMR+R0','OMOMR+',
     & 'P-PO  ','P-POPO','P-P-P+','P-ET  ','POR-  ','P-R0  ','ETSP- ',
     & 'R-ET  ','R-R0  ','POETR-','P-ETR0','POR-R0','P-R+R-','P-R0R0'/
      DATA ZKNAM6/'P+R-R-','R-OM  ','P-ETOM','ETSR- ','POR-OM','P-R0OM',
     & 'R-ETET','R-R0ET','P-OMOM','R0R0R-','R+R-R-','ETSR-E','OMR0R-',
     & 'OMOMR-', 'PAN-14','APN+14','NAN014','ANN014','PAKO  ','LPI+  ',
     & 'SI+OM','LAMRO+','SI0RO+','SI+RO0','SI+ETA','SI0PI+','SI+PI0',
     & 'APETA ','AN=P+ ','AN-PO ','ANOPO ','APRHOO','ANRHO-','ANETA ',
     & 'AN-P+ ','AN0PO ','AN+P- ','APRHO+','ANRHO0',
     & 'KONPIO','KOPPI-','K+NPI-','AKOPPO','AKONP+','K-PPI+',
     & 9*'BLANK'/
*=                                               end*block.zk      *
      END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*-- Author :
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

*$ CREATE BLKDT3.FOR
*COPY BLKDT3
*
*=== blkdt3 ===========================================================*
*==                                                                    *
      BLOCK DATA BLKD43
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      INCLUDE '(DBLPRC)'
C      INCLUDE '(DIMPAR)'
C      INCLUDE '(IOUNIT)'
* * * reaction channel cross section data                              *
C      INCLUDE '(REAC)'

*$ CREATE REAC.ADD
*COPY REAC
*
*=== reac =============================================================*
*
*----------------------------------------------------------------------*
*                                                                      *
*     Created on 10 december 1991  by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 10-dec-91     by    Alfredo Ferrari               *
*                                                                      *
*     This is the original common reac of Hadrin                       *
*                                                                      *
*----------------------------------------------------------------------*
*
      COMMON /DREAC/ UMO   ( 296), PLABF ( 296), SIIN  ( 296),
     +              WK    (5184), NRK  (2,268), NURE  (30,2)

      DIMENSION
     & UMOPI(92), UMOKC(68), UMOP(39), UMON(63), UMOK0(34),
     & PLAPI(92), PLAKC(68), PLAP(39), PLAN(63), PLAK0(34),
     & SPIKP1(315), SPIKPU(278), SPIKPV(372),
     & SPIKPW(278), SPIKPX(372), SPIKP4(315),
     & SPIKP5(187), SPIKP6(289),
     & SKMPEL(102), SPIKP7(289), SKMNEL(68), SPIKP8(187),
     & SPIKP9(143), SPIKP0(169), SPKPV(143),
     & SAPPEL(105), SPIKPE(399), SAPNEL(84), SPIKPZ(273),
     & SANPEL(84) , SPIKPF(273),
     & SPKP15(187), SPKP16(272),
     & NRKPI(164), NRKKC(132), NRKP(70), NRKN(116), NRKK0(54),
     & NURELN(60)
*
       DIMENSION NRKLIN(532)
       EQUIVALENCE (NRK(1,1), NRKLIN(1))
       EQUIVALENCE (   UMO(  1),  UMOPI(1)), (   UMO( 93),  UMOKC(1))
       EQUIVALENCE (   UMO(161),   UMOP(1)), (   UMO(200),   UMON(1))
       EQUIVALENCE (   UMO(263),  UMOK0(1))
       EQUIVALENCE ( PLABF(  1),  PLAPI(1)), ( PLABF( 93),  PLAKC(1))
       EQUIVALENCE ( PLABF(161),   PLAP(1)), ( PLABF(200),   PLAN(1))
       EQUIVALENCE ( PLABF(263),  PLAK0(1))
       EQUIVALENCE (   WK(   1), SPIKP1(1)), (   WK( 316), SPIKPU(1))
       EQUIVALENCE (   WK( 594), SPIKPV(1)), (   WK( 966), SPIKPW(1))
       EQUIVALENCE (   WK(1244), SPIKPX(1)), (   WK(1616), SPIKP4(1))
       EQUIVALENCE (   WK(1931), SPIKP5(1)), (   WK(2118), SPIKP6(1))
       EQUIVALENCE (   WK(2407), SKMPEL(1)), (   WK(2509), SPIKP7(1))
       EQUIVALENCE (   WK(2798), SKMNEL(1)), (   WK(2866), SPIKP8(1))
       EQUIVALENCE (   WK(3053), SPIKP9(1)), (   WK(3196), SPIKP0(1))
       EQUIVALENCE (   WK(3365),  SPKPV(1)), (   WK(3508), SAPPEL(1))
       EQUIVALENCE (   WK(3613), SPIKPE(1)), (   WK(4012), SAPNEL(1))
       EQUIVALENCE (   WK(4096), SPIKPZ(1)), (   WK(4369), SANPEL(1))
       EQUIVALENCE (   WK(4453), SPIKPF(1)), (   WK(4726), SPKP15(1))
       EQUIVALENCE (   WK(4913), SPKP16(1))
       EQUIVALENCE (NRK(1,1), NRKLIN(1))
       EQUIVALENCE (NRKLIN(   1), NRKPI(1)), (NRKLIN( 165), NRKKC(1))
       EQUIVALENCE (NRKLIN( 297),  NRKP(1)), (NRKLIN( 367),  NRKN(1))
       EQUIVALENCE (NRKLIN( 483), NRKK0(1))
       EQUIVALENCE (NURE(1,1), NURELN(1))
*
**** pi- p data                                                        *
**** pi+ n data                                                        *
      DATA PLAPI / 0.D0, .3D0, .5D0, .6D0, .7D0, .8D0, .9D0, .95D0,1.D0,
     & 1.15D0, 1.3D0, 1.5D0, 1.6D0, 1.8D0, 2.D0, 2.3D0, 2.5D0, 2.8D0,
     & 3.D0, 3.5D0, 4.D0, 0.D0, .285D0, .4D0, .45D0, .5D0, .6D0, .7D0,
     & .75D0, .8D0, .85D0, .9D0, 1.D0, 1.15D0, 1.3D0, 1.5D0, 1.6D0,
     & 1.8D0, 2.D0, 2.3D0, 2.5D0, 2.8D0, 3.D0, 3.5D0, 4.D0, 4.5D0, 0.D0,
     & .285D0, .4D0, .45D0, .5D0, .6D0, .7D0, .75D0, .8D0, .85D0, .9D0,
     & 1.D0, 1.15D0, 1.3D0, 1.5D0, 1.6D0, 1.8D0, 2.D0, 2.3D0, 2.5D0,
     & 2.8D0, 3.D0, 3.5D0, 4.D0, 4.5D0, 0.D0, .3D0, .5D0, .6D0, .7D0,
     & .8D0, .9D0, .95D0, 1.D0, 1.15D0, 1.3D0, 1.5D0, 1.6D0, 1.8D0,
     & 2.D0, 2.3D0, 2.5D0, 2.8D0, 3.D0, 3.5D0, 4.D0 /
      DATA PLAKC /
     &   0.D0,  .58D0,   .8D0, 1.01D0, 1.23D0, 1.45D0, 1.68D0, 1.94D0,
     & 2.18D0, 2.42D0, 2.68D0, 2.96D0, 3.24D0,
     & 3.51D0, 3.84D0, 4.16D0, 4.49D0,
     &   0.D0,  .58D0,   .8D0, 1.01D0, 1.23D0, 1.45D0, 1.68D0, 1.94D0,
     & 2.18D0, 2.42D0, 2.68D0, 2.96D0, 3.24D0,
     & 3.51D0, 3.84D0, 4.16D0, 4.49D0,
     &   0.D0,  .58D0,   .8D0, 1.01D0, 1.23D0, 1.45D0, 1.68D0, 1.94D0,
     & 2.18D0, 2.42D0, 2.68D0, 2.96D0, 3.24D0,
     & 3.51D0, 3.84D0, 4.16D0, 4.49D0,
     &   0.D0,  .58D0,   .8D0, 1.01D0, 1.23D0, 1.45D0, 1.68D0, 1.94D0,
     & 2.18D0, 2.42D0, 2.68D0, 2.96D0, 3.24D0,
     & 3.51D0, 3.84D0, 4.16D0, 4.49D0/
      DATA PLAK0 /
     &   0.D0,  .58D0,   .8D0, 1.01D0, 1.23D0, 1.45D0, 1.68D0, 1.94D0,
     & 2.18D0, 2.42D0, 2.68D0, 2.96D0, 3.24D0,
     & 3.51D0, 3.84D0, 4.16D0, 4.49D0,
     &   0.D0,  .58D0,   .8D0, 1.01D0, 1.23D0, 1.45D0, 1.68D0, 1.94D0,
     & 2.18D0, 2.42D0, 2.68D0, 2.96D0, 3.24D0,
     & 3.51D0, 3.84D0, 4.16D0, 4.49D0/
*                 pp   pn   np   nn                                    *
      DATA PLAP /
     &   0.D0, 1.06D0, 1.34D0, 1.63D0, 1.92D0, 2.2D0, 2.5D0,2.8D0,3.1D0,
     & 3.43D0, 3.75D0, 4.07D0, 4.43D0,
     &   0.D0, 1.06D0, 1.34D0, 1.63D0, 1.92D0, 2.2D0, 2.5D0,2.8D0,3.1D0,
     & 3.43D0, 3.75D0, 4.07D0, 4.43D0,
     &   0.D0, 1.06D0, 1.34D0, 1.63D0, 1.92D0, 2.2D0, 2.5D0,2.8D0,3.1D0,
     & 3.43D0, 3.75D0, 4.07D0, 4.43D0 /
*    app   apn   anp   ann                                             *
      DATA PLAN /
     &  0.D0,   1.D-3,   .1D0,   .2D0,   .3D0,  .4D0,  .5D0, .6D0,
     & .74D0,  1.06D0, 1.34D0, 1.63D0, 1.92D0, 2.2D0, 2.5D0,2.8D0,3.1D0,
     & 3.43D0, 3.75D0, 4.07D0, 4.43D0,
     &  0.D0,   1.D-3,   .1D0,   .2D0,   .3D0,  .4D0,  .5D0, .6D0,
     & .74D0,  1.06D0, 1.34D0, 1.63D0, 1.92D0, 2.2D0, 2.5D0,2.8D0,3.1D0,
     & 3.43D0, 3.75D0, 4.07D0, 4.43D0,
     &  0.D0,   1.D-3,   .1D0,   .2D0,   .3D0,  .4D0,  .5D0, .6D0,
     & .74D0,  1.06D0, 1.34D0, 1.63D0, 1.92D0, 2.2D0, 2.5D0,2.8D0,3.1D0,
     & 3.43D0, 3.75D0, 4.07D0, 4.43D0  /
      DATA SIIN / 296*0.D0 /
      DATA UMOPI/ 1.08D0,1.233D0,1.302D0,1.369D0,1.496D0,
     & 1.557D0,1.615D0,1.6435D0,
     & 1.672D0,1.753D0,1.831D0,1.930D0,1.978D0,2.071D0,2.159D0,
     & 2.286D0,2.366D0,2.482D0,2.56D0,
     & 2.735D0,2.90D0,
     &             1.08D0,1.222D0,1.302D0,1.3365D0,1.369D0,1.434D0,
     & 1.496D0,1.527D0,1.557D0,
     & 1.586D0,1.615D0,1.672D0,1.753D0,1.831D0,1.930D0,1.978D0,
     & 2.071D0,2.159D0,2.286D0,2.366D0,
     & 2.482D0,2.560D0,2.735D0,2.90D0,3.06D0,
     &             1.08D0,1.222D0,1.302D0,1.3365D0,1.369D0,1.434D0,
     & 1.496D0,1.527D0,1.557D0,
     & 1.586D0,1.615D0,1.672D0,1.753D0,1.831D0,1.930D0,1.978D0,
     & 2.071D0,2.159D0,2.286D0,2.366D0,
     & 2.482D0,2.560D0,2.735D0,2.90D0,3.06D0,
     &                   1.08D0,1.233D0,1.302D0,1.369D0,1.496D0,
     & 1.557D0,1.615D0,1.6435D0,
     & 1.672D0,1.753D0,1.831D0,1.930D0,1.978D0,2.071D0,2.159D0,
     & 2.286D0,2.366D0,2.482D0,2.56D0,
     &  2.735D0, 2.90D0/
      DATA UMOKC/ 1.44D0,
     &  1.598D0,1.7D0,1.8D0,1.9D0,2.0D0,2.1D0,2.2D0,2.3D0,2.4D0,2.5D0,
     & 2.6D0,2.7D0,2.8D0,2.9D0,3.0D0,
     & 3.1D0,1.44D0,
     &  1.598D0,1.7D0,1.8D0,1.9D0,2.0D0,2.1D0,2.2D0,2.3D0,2.4D0,2.5D0,
     & 2.6D0,2.7D0,2.8D0,2.9D0,3.0D0,
     & 3.1D0,1.44D0,
     &  1.598D0,1.7D0,1.8D0,1.9D0,2.0D0,2.1D0,2.2D0,2.3D0,2.4D0,2.5D0,
     & 2.6D0,2.7D0,2.8D0,2.9D0,3.0D0,
     & 3.1D0,1.44D0,
     &  1.598D0,1.7D0,1.8D0,1.9D0,2.0D0,2.1D0,2.2D0,2.3D0,2.4D0,2.5D0,
     & 2.6D0,2.7D0,2.8D0,2.9D0,3.0D0,
     &  3.1D0/
      DATA UMOK0/ 1.44D0,
     &  1.598D0,1.7D0,1.8D0,1.9D0,2.0D0,2.1D0,2.2D0,2.3D0,2.4D0,2.5D0,
     & 2.6D0,2.7D0,2.8D0,2.9D0,3.0D0,
     & 3.1D0,1.44D0,
     &  1.598D0,1.7D0,1.8D0,1.9D0,2.0D0,2.1D0,2.2D0,2.3D0,2.4D0,2.5D0,
     & 2.6D0,2.7D0,2.8D0,2.9D0,3.0D0,
     &  3.1D0/
*                 pp   pn   np   nn                                    *
      DATA UMOP/
     & 1.88D0,2.102D0,2.2D0,2.3D0,2.4D0,2.5D0,2.6D0,2.7D0,2.8D0,2.9D0,
     & 3.D0,3.1D0,3.2D0,
     & 1.88D0,2.102D0,2.2D0,2.3D0,2.4D0,2.5D0,2.6D0,2.7D0,2.8D0,2.9D0,
     & 3.D0,3.1D0,3.2D0,
     & 1.88D0,2.102D0,2.2D0,2.3D0,2.4D0,2.5D0,2.6D0,2.7D0,2.8D0,2.9D0,
     & 3.D0,3.1D0,3.2D0/
*    app   apn   anp   ann                                             *
      DATA UMON /
     & 1.877D0,1.87701D0,1.879D0,1.887D0,1.9D0,1.917D0,1.938D0,1.962D0,
     & 2.D0,2.102D0,2.2D0,2.3D0,2.4D0,2.5D0,2.6D0,2.7D0,2.8D0,2.9D0,
     & 3.D0,3.1D0,3.2D0,
     & 1.877D0,1.87701D0,1.879D0,1.887D0,1.9D0,1.917D0,1.938D0,1.962D0,
     & 2.D0,2.102D0,2.2D0,2.3D0,2.4D0,2.5D0,2.6D0,2.7D0,2.8D0,2.9D0,
     & 3.D0,3.1D0,3.2D0,
     & 1.877D0,1.87701D0,1.879D0,1.887D0,1.9D0,1.917D0,1.938D0,1.962D0,
     & 2.D0,2.102D0,2.2D0,2.3D0,2.4D0,2.5D0,2.6D0,2.7D0,2.8D0,2.9D0,
     &  3.D0,3.1D0,3.2D0/
**** reaction channel state particles                                  *
      DATA NRKPI / 13, 1, 15, 21, 81, 0, 13, 54, 23, 53, 13, 63, 13, 58,
     & 23, 57, 13, 65, 1, 32, 53, 31, 54, 32, 53, 33, 53, 35, 63, 32,
     & 13, 8, 23, 1, 17, 15, 21, 24, 22, 15, 82, 0, 61, 0, 13, 55, 23,
     & 54, 14, 53, 13, 64, 23, 63, 13, 59, 23, 58, 14, 57, 13, 66, 23,
     & 65, 1, 31, 8, 32, 1, 33, 1, 35, 54, 31, 55, 32, 54, 33, 53, 34,
     & 54, 35, 14, 1, 23, 8, 17, 24, 20, 15, 22, 24, 83, 0, 62, 0, 14,
     & 54, 23, 55, 13, 56, 14, 63, 23, 64, 14, 58, 23, 59, 13, 60, 14,
     & 65, 23, 66, 8, 31, 1, 34, 8, 33, 8, 35, 55, 31, 54, 34, 55, 33,
     & 56, 32, 55, 35, 14, 8, 24, 20, 84, 0, 14, 55, 23, 56, 14, 64, 14,
     & 59, 23, 60, 14, 66, 8, 34, 56, 31, 55, 34, 56, 33, 56, 35, 64,34/
      DATA NRKKC/ 15, 1, 89, 0, 24, 53, 15, 54, 1, 36, 1, 40, 1, 44, 36,
     & 63, 15, 63, 45, 53, 44, 54, 15, 8, 24, 1, 91, 0, 24, 54, 15, 55,
     & 8, 36, 1, 37, 8, 40, 1, 41, 8, 44, 1, 45, 36, 64, 37, 63, 15, 64,
     & 24, 63, 45, 54, 44, 55, 16, 1, 25, 8, 17, 23, 21, 14, 20,
     & 13, 22, 23, 90, 0, 38, 1, 39, 8, 16, 54, 25, 55, 1, 42, 8, 43,
     & 16, 63, 25, 64, 39, 64, 38, 63, 46, 54, 47, 55, 8, 47, 1, 46, 52,
     & 0, 51, 0, 16, 8, 17, 14, 20, 23, 22, 14, 92, 0, 8, 38, 16, 55,
     & 25, 56, 8, 42, 16, 64, 38, 64, 46, 55, 47, 56, 8, 46, 94, 0 /
*                                                                      *
*   k0 p   k0 n   ak0 p   ak/ n                                        *
*                                                                      *
      DATA NRKK0 / 24, 8, 106, 0, 15, 56, 24, 55, 37, 8, 41, 8, 45, 8,
     & 37, 64, 24, 64, 44, 56, 45, 55, 25, 1, 17, 13,   22, 13, 21, 23,
     & 107, 0, 39, 1, 25, 54, 16, 53, 43, 1, 25, 63, 39, 63, 47, 54, 46,
     & 53, 47, 1, 103, 0, 93, 0/
*   pp  pn   np   nn                                                   *
      DATA NRKP / 1, 1, 85, 0, 8, 53, 1, 54, 1, 63, 8, 57, 1, 58, 2*54,
     & 53, 55, 63, 54, 64, 53, 1, 8, 86, 0, 8, 54, 1, 55, 8, 63, 1, 64,
     & 8, 58, 1, 59, 64, 54, 63, 55, 54, 55, 53, 56, 77, 0, 2*8, 95, 0,
     & 8, 55, 1, 56, 8, 64, 8, 59, 1, 60, 2*55, 54, 56, 64, 55, 63, 56 /
*     app   apn   anp   ann                                            *
      DATA NRKN/ 1, 2, 17, 18, 15, 16, 8, 9, 13, 14, 99, 0, 87, 0, 1,
     & 68, 8, 69, 2, 54, 9, 55, 102, 0, 2, 63, 9, 64, 1, 75, 8, 76, 53,
     & 67, 54, 68, 55, 69, 56, 70, 63, 68, 64, 69, 75, 54, 76, 55, 2, 8,
     & 18, 20, 16, 24, 14, 23, 101, 0, 88, 0, 2, 55, 9, 56, 1, 67, 8,
     & 68, 2, 64, 8, 75, 2, 59, 8, 72, 68, 55, 67, 54, 69, 56, 1, 9, 18,
     & 21, 15, 25, 13, 23, 100, 0, 96, 0, 2, 53, 9, 54, 1, 69, 8, 70, 1,
     & 76, 9, 63, 1, 73, 9, 58, 55, 70, 53, 68, 54, 69 /
**** channel cross section                                             *
      DATA SPIKP1/ 0.D0, 300.D0, 40.D0, 20.D0, 13.D0,8.5D0,8.D0, 9.5D0,
     & 12.D0,14.D0,15.5D0,20.D0,17.D0,13.D0,10.D0,9.D0,8.5D0,8.D0,7.8D0,
     & 7.3D0, 6.7D0, 9*0.D0,.23D0,.35D0,.7D0,.52D0,.4D0,.3D0,.2D0,.15D0,
     & .13D0, .11D0, .09D0, .07D0, 0.D0, .033D0,.8D0,1.35D0,1.35D0,.5D0,
     & 15*0.D0, 3*0.D0,.00D0,0.80D0,2.2D0,3.6D0,4.6D0,4.7D0,3.5D0,2.4D0,
     &1.8D0,1.4D0,.75D0,.47D0,.25D0,.13D0,.08D0,6*0.D0,0.D0,1.2D0,3.3D0,
     & 5.4D0,6.9D0,7.3D0,5.3D0,3.6D0,2.7D0,2.2D0,1.1D0,.73D0,.4D0,.22D0,
     & .12D0,9*0.D0,.0D0,0.D0,2.0D0,4.4D0,6.8D0,9.9D0,7.9D0,6.0D0,3.8D0,
     &2.5D0,2.D0,1.4D0,1.D0,.6D0,.35D0,10*0.D0,.25D0,.55D0,.75D0,1.25D0,
     & 1.9D0,2.D0,1.8D0,1.5D0,1.25D0,1.D0,.8D0,6*0.D0,4*0.D0,.4D0,.85D0,
     & 1.1D0, 1.85D0, 2.8D0, 3.D0,2.7D0,2.2D0,1.85D0,1.5D0,1.2D0,6*0.D0,
     & 6*0.D0, .5D0, 1.2D0, 1.7D0, 3.4D0, 5.2D0, 6.4D0, 6.1D0, 5.6D0,
     & 5.2D0, 6*0.D0, 2*0.D0, .0D0, 1.D0, 3.3D0, 5.2D0, 4.45D0, 3.6D0,
     & 2.75D0, 1.9D0, 1.65D0, 1.3D0, .95D0, .6D0, .45D0, 6*0.D0, 3*0.D0,
     & .0D0, .45D0, 1.4D0, 1.5D0, 1.1D0, .85D0, .5D0, .3D0, .2D0, .15D0,
     & 8*0.D0, 5*0.D0, .0D0, .0D0, .6D0, .8D0, .95D0, .8D0, .7D0, .6D0,
     & .5D0, .4D0, 6*0.D0, 5*0.D0, .0D0, .00D0, .85D0, 1.2D0, 1.4D0,
     & 1.2D0, 1.05D0, .9D0, .7D0, .55D0, 6*0.D0, 5*0.D0, .0D0, .00D0,
     & 1.D0, 1.5D0, 3.5D0, 4.15D0, 3.7D0, 2.7D0, 2.3D0, 1.75D0, 6*0.D0,
     & 10*0.D0, .5D0, 2.0D0, 3.3D0, 5.4D0, 7.D0 /
**** pi+ n data                                                        *
      DATA SPIKPU/   0.D0, 25.D0, 13.D0,  11.D0, 10.5D0, 14.D0,  20.D0,
     & 20.D0, 16.D0, 14.D0, 19.D0, 28.D0, 17.5D0, 13.5D0, 12.D0, 10.5D0,
     & 10.D0, 10.D0, 9.5D0,  9.D0, 8.D0, 7.5D0, 7.D0, 6.5D0, 6.D0, 0.D0,
     & 48.D0, 19.D0, 15.D0, 11.5D0, 10.D0, 8.D0, 6.5D0,   5.5D0,  4.8D0,
     & 4.2D0, 7.5D0, 3.4D0,  2.5D0, 2.5D0, 2.1D0, 1.4D0,   1.D0,   .8D0,
     &  .6D0, .46D0,  .3D0, .2D0, .15D0, .13D0, 11*0.D0,  .95D0,  .65D0,
     & .48D0, .35D0,  .2D0, .18D0, .17D0, .16D0,  .15D0,   .1D0,  .09D0,
     & .065D0, .05D0, .04D0, 12*0.D0, .2D0, .25D0, .25D0,  .2D0,   .1D0,
     & .08D0, .06D0, .045D0,   .03D0, .02D0, .01D0,      .005D0, .003D0,
     & 12*0.D0, .3D0, .24D0,   .18D0, .15D0, .13D0,  .12D0, .11D0, .1D0,
     & .09D0,  .08D0, .05D0,   .04D0, .03D0,  0.D0, 0.16D0, .7D0, 1.3D0,
     & 3.1D0,  4.5D0,  2.D0, 18*0.D0, 3*.0D0,  0.D0, 0.D0, 4.0D0, 11.D0,
     & 11.4D0, 10.3D0, 7.5D0, 6.8D0, 4.75D0, 2.5D0,  1.5D0, .9D0, .55D0,
     &  .35D0, 13*0.D0, .1D0, .34D0, .5D0, .8D0, 1.1D0,   2.25D0, 3.3D0,
     & 2.3D0, 1.6D0, .95D0, .45D0, .28D0, .15D0, 10*0.D0, 2*0.D0, .17D0,
     & .64D0,  1.D0, 1.5D0, 2.1D0, 4.25D0, 6.2D0,  4.4D0,   3.D0, 1.8D0,
     &  .9D0, .53D0, .28D0,      10*0.D0, 2*0.D0,  .25D0,  .82D0,
     & 1.3D0, 1.9D0, 2.8D0, 5.5D0 , 8.D0,  5.7D0, 3.9D0, 2.35D0, 1.15D0,
     & .69D0, .37D0, 10*0.D0,     7*0.D0,   .0D0, .34D0,  1.5D0, 3.47D0,
     & 5.87D0, 6.23D0, 4.27D0, 2.6D0, 1.D0, .6D0,  .3D0,  .15D0, 6*0.D0/
*
      DATA SPIKPV/ 7*0.D0, .00D0, .16D0, .75D0, 1.73D0, 2.93D0, 3.12D0,
     & 2.13D0, 1.3D0, .5D0, .3D0, .15D0, .08D0, 6*0.D0, 10*0.D0, .2D0,
     & .6D0, .92D0, 2.4D0, 4.9D0, 6.25D0, 5.25D0, 3.5D0, 2.15D0, 1.4D0,
     & 1.D0, .7D0, 13*0.D0, .13D0, .4D0, .62D0, 1.6D0, 3.27D0, 4.17D0,
     & 3.5D0, 2.33D0, 1.43D0, .93D0, .66D0, .47D0, 13*0.D0, .07D0, .2D0,
     & .31D0, .8D0, 1.63D0, 2.08D0, 1.75D0, 1.17D0, .72D0, .47D0, .34D0,
     & .23D0, 17*0.D0, .33D0, 1.D0, 1.8D0, 2.67D0, 5.33D0, 6.D0, 5.53D0,
     & 5.D0, 17*0.D0, .17D0, .5D0, .9D0, 1.83D0, 2.67D0, 3.0D0, 2.77D0,
     & 2.5D0, 3*0.D0, 3*0.D0, 1.D0, 3.3D0, 2.8D0, 2.5D0, 2.3D0, 1.8D0,
     & 1.5D0, 1.1D0, .8D0, .7D0, .55D0, .3D0, 10*0.D0, 9*0.D0, .1D0,
     & .4D0, 1.D0, 1.4D0, 2.2D0, 2.5D0, 2.2D0, 1.65D0, 1.35D0, 1.1D0,
     & .8D0, .6D0, .4D0, 12*0.D0, .15D0, .6D0, 1.5D0, 2.1D0, 3.3D0,
     & 3.8D0, 3.3D0, 2.45D0, 2.05D0, 1.65D0, 1.2D0, .9D0, .6D0, 3*0.D0,
     & 9*0.D0, .10D0, .2D0, .5D0, .7D0, 1.3D0, 1.55D0, 1.9D0, 1.8D0,
     & 1.55D0, 1.35D0, 1.15D0, .95D0, .7D0, 13*0.D0, .2D0, .5D0, .7D0,
     & 1.3D0, 1.55D0, 1.9D0, 1.8D0, 1.55D0, 1.35D0, 1.15D0, .95D0, .7D0,
     & 17*0.D0, .2D0, .5D0, .85D0, 2.D0, 2.15D0, 2.05D0, 1.75D0, 1.D0,
     & 17*0.D0, .13D0, .33D0, .57D0, 1.33D0, 1.43D0, 1.36D0, 1.17D0,
     & .67D0, 17*0.D0, .07D0, .17D0, .28D0, .67D0, .72D0, .69D0, .58D0,
     & .33D0,17*0.D0,.4D0, .7D0, 1.D0, 1.6D0, 1.8D0, 2.3D0,1.9D0,1.7D0 /
**** pi- p data                                                        *
      DATA SPIKPW/ 0.D0, 25.D0, 13.D0, 11.D0, 10.5D0, 14.D0, 2*20.D0,
     & 16.D0, 14.D0, 19.D0, 28.D0, 17.5D0, 13.5D0, 12.D0, 10.5D0,
     & 2*10.D0, 9.5D0, 9.D0, 8.D0, 7.5D0, 7.D0, 6.5D0, 6.D0, 0.D0,
     & 48.D0, 19.D0, 15.D0, 11.5D0, 10.D0, 8.D0, 6.5D0, 5.5D0, 4.8D0,
     & 4.2D0, 7.5D0, 3.4D0, 2*2.5D0, 2.1D0, 1.4D0, 1.D0, .8D0, .6D0,
     & .46D0, .3D0, .2D0, .15D0, .13D0, 11*0.D0, .95D0, .65D0, .48D0,
     & .35D0, .2D0, .18D0, .17D0, .16D0, .15D0, .1D0, .09D0, .065D0,
     & .05D0, .04D0, 12*0.D0, .2D0, 2*.25D0, .2D0, .1D0, .08D0, .06D0,
     & .045D0, .03D0, .02D0, .01D0, .005D0, .003D0, 12*0.D0, .3D0,
     & .24D0, .18D0, .15D0, .13D0, .12D0, .11D0, .1D0, .09D0, .08D0,
     & .05D0, .04D0, .03D0, 0.D0, 0.16D0, .7D0, 1.3D0, 3.1D0, 4.5D0,
     & 2.D0, 23*0.D0, 4.0D0, 11.D0, 11.4D0, 10.3D0, 7.5D0, 6.8D0,
     & 4.75D0, 2.5D0, 1.5D0, .9D0, .55D0, .35D0, 13*0.D0, .1D0, .34D0,
     & .5D0, .8D0, 1.1D0, 2.25D0, 3.3D0, 2.3D0, 1.6D0, .95D0, .45D0,
     & .28D0, .15D0, 12*0.D0, .17D0, .64D0, 1.D0, 1.5D0, 2.1D0, 4.25D0,
     & 6.2D0, 4.4D0, 3.D0, 1.8D0, .9D0, .53D0, .28D0, 12*0.D0, .25D0,
     & .82D0, 1.3D0, 1.9D0, 2.8D0, 5.5D0, 8.D0, 5.7D0, 3.9D0, 2.35D0,
     & 1.15D0, .69D0, .37D0, 18*0.D0, .34D0, 1.5D0, 3.47D0, 5.87D0,
     & 6.23D0, 4.27D0, 2.6D0, 1.D0, .6D0, .3D0, .15D0, 6*0.D0/
*
      DATA SPIKPX/ 8*0.D0, .16D0, .75D0, 1.73D0, 2.93D0, 3.12D0,
     & 2.13D0, 1.3D0, .5D0, .3D0, .15D0, .08D0, 16*0.D0, .2D0, .6D0,
     & .92D0, 2.4D0, 4.9D0, 6.25D0, 5.25D0, 3.5D0, 2.15D0, 1.4D0, 1.D0,
     & .7D0, 13*0.D0, .13D0, .4D0, .62D0, 1.6D0, 3.27D0, 4.17D0, 3.5D0,
     & 2.33D0, 1.43D0, .93D0, .66D0, .47D0, 13*0.D0, .07D0, .2D0, .31D0,
     & .8D0, 1.63D0, 2.08D0, 1.75D0, 1.17D0, .72D0, .47D0, .34D0, .23D0,
     & 17*0.D0, .33D0, 1.D0, 1.8D0, 2.67D0, 5.33D0, 6.D0, 5.53D0, 5.D0,
     & 17*0.D0, .17D0, .5D0, .9D0, 1.83D0, 2.67D0, 3.0D0, 2.77D0, 2.5D0,
     & 6*0.D0, 1.D0, 3.3D0, 2.8D0, 2.5D0, 2.3D0, 1.8D0, 1.5D0, 1.1D0,
     & .8D0, .7D0, .55D0, .3D0, 19*0.D0, .1D0, .4D0, 1.D0, 1.4D0, 2.2D0,
     & 2.5D0, 2.2D0, 1.65D0, 1.35D0, 1.1D0, .8D0, .6D0, .4D0, 12*0.D0,
     & .15D0, .6D0, 1.5D0, 2.1D0, 3.3D0, 3.8D0, 3.3D0, 2.45D0, 2.05D0,
     & 1.65D0, 1.2D0, .9D0, .6D0, 12*0.D0, .10D0, .2D0, .5D0, .7D0,
     & 1.3D0, 1.55D0, 1.9D0, 1.8D0, 1.55D0, 1.35D0, 1.15D0, .95D0, .7D0,
     & 13*0.D0, .2D0, .5D0, .7D0, 1.3D0, 1.55D0, 1.9D0, 1.8D0, 1.55D0,
     & 1.35D0, 1.15D0, .95D0, .7D0, 17*0.D0, .2D0, .5D0, .85D0, 2.D0,
     & 2.15D0, 2.05D0, 1.75D0, 1.D0, 17*0.D0, .13D0, .33D0, .57D0,
     & 1.33D0, 1.43D0, 1.36D0, 1.17D0, .67D0, 17*0.D0, .07D0, .17D0,
     & .28D0, .67D0, .72D0, .69D0, .58D0, .33D0, 17*0.D0, .4D0, .7D0,
     & 1.D0, 1.6D0, 1.8D0, 2.3D0, 1.9D0, 1.7D0 /
**** pi- n data                                                        *
      DATA SPIKP4 / 0.D0, 300.D0, 40.D0, 20.D0, 13.D0, 8.5D0, 8.D0,
     & 9.5D0, 12.D0, 14.D0, 15.5D0, 20.D0, 17.D0, 13.D0, 10.D0, 9.D0,
     & 8.5D0, 8.D0, 7.8D0, 7.3D0, 6.7D0, 9*0.D0, .23D0, .35D0, .7D0,
     & .52D0, .4D0, .3D0, .2D0, .15D0, .13D0, .11D0, .09D0, .07D0, 0.D0,
     & .033D0, .8D0, 2*1.35D0, .5D0, 19*0.D0, 0.8D0, 2.2D0, 3.6D0,
     & 4.6D0, 4.7D0, 3.5D0, 2.4D0, 1.8D0, 1.4D0, .75D0, .47D0, .25D0,
     & .13D0, .08D0, 7*0.D0, 1.2D0, 3.3D0, 5.4D0, 6.9D0, 7.3D0, 5.3D0,
     & 3.6D0, 2.7D0, 2.2D0, 1.1D0, .73D0, .4D0, .22D0, .12D0, 11*0.D0,
     & 2.0D0, 4.4D0, 6.8D0, 9.9D0, 7.9D0, 6.0D0, 3.8D0, 2.5D0, 2.D0,
     & 1.4D0, 1.D0, .6D0, .35D0, 10*0.D0, .25D0, .55D0, .75D0, 1.25D0,
     & 1.9D0, 2.D0, 1.8D0, 1.5D0, 1.25D0, 1.D0, .8D0, 10*0.D0, .4D0,
     & .85D0, 1.1D0, 1.85D0, 2.8D0, 3.D0, 2.7D0, 2.2D0, 1.85D0, 1.5D0,
     & 1.2D0, 12*0.D0, .5D0, 1.2D0, 1.7D0, 3.4D0, 5.2D0, 6.4D0, 6.1D0,
     & 5.6D0, 5.2D0, 9*0.D0, 1.D0, 3.3D0, 5.2D0, 4.45D0, 3.6D0, 2.75D0,
     & 1.9D0, 1.65D0, 1.3D0, .95D0, .6D0, .45D0, 10*0.D0, .45D0, 1.4D0,
     & 1.5D0, 1.1D0, .85D0, .5D0, .3D0, .2D0, .15D0, 15*0.D0, .6D0,
     & .8D0, .95D0, .8D0, .7D0, .6D0, .5D0, .4D0, 13*0.D0, .85D0, 1.2D0,
     & 1.4D0, 1.2D0, 1.05D0, .9D0, .7D0, .55D0, 13*0.D0, 1.D0, 1.5D0,
     & 3.5D0, 4.15D0, 3.7D0, 2.7D0, 2.3D0, 1.75D0, 16*0.D0, .5D0, 2.0D0,
     & 3.3D0, 5.4D0, 7.D0 /
**** k+  p data                                                        *
      DATA SPIKP5/ 0.D0, 20.D0, 14.D0, 12.D0, 11.5D0, 10.D0, 8.D0,
     & 7.D0, 6.D0, 5.5D0, 5.3D0, 5.D0, 4.5D0, 4.4D0, 3.8D0, 3.D0, 2.8D0,
     & 0.D0, .5D0, 1.15D0, 2.D0, 1.3D0, .8D0, .45D0, 13*0.D0, 0.9D0,
     & 2.5D0, 3.D0, 2.5D0, 2.3D0, 2.D0, 1.7D0, 1.5D0, 1.2D0, .9D0, .6D0,
     & .45D0, .21D0, .2D0, 3*0.D0, .9D0, 2.5D0, 3.D0, 2.5D0, 2.3D0,
     & 2.D0, 1.7D0, 1.5D0, 1.2D0, .9D0, .6D0, .45D0, .21D0, .2D0,
     & 4*0.D0, 1.D0, 2.1D0, 2.6D0, 2.3D0, 2.1D0, 1.8D0, 1.7D0, 1.4D0,
     & 1.2D0, 1.05D0, .9D0, .66D0, .5D0, 7*0.D0, .3D0, 2*1.D0, .9D0,
     & .7D0, .4D0, .3D0, .2D0, 11*0.D0, .1D0, 1.D0, 2.2D0, 3.5D0, 4.2D0,
     & 4.55D0, 4.85D0, 4.9D0, 10*0.D0, .2D0, .7D0, 1.6D0, 2.5D0, 2.2D0,
     & 1.71D0, 1.6D0, 6*0.D0, 1.4D0, 3.8D0, 5.D0, 4.7D0, 4.4D0, 4.D0,
     & 3.5D0, 2.85D0, 2.35D0, 2.01D0, 1.8D0, 12*0.D0, .1D0, .8D0,2.05D0,
     & 3.31D0, 3.5D0, 12*0.D0, .034D0, .2D0, .75D0, 1.04D0, 1.24D0 /
**** k+  n data                                                        *
      DATA SPIKP6/ 0.D0, 6.D0, 11.D0, 13.D0, 6.D0, 5.D0, 3.D0, 2.2D0,
     & 1.5D0, 1.2D0, 1.D0, .7D0, .6D0, .5D0, .45D0, .35D0, .3D0, 0.D0,
     & 6.D0, 11.D0, 13.D0, 6.D0, 5.D0, 3.D0, 2.2D0, 1.5D0, 1.2D0, 1.D0,
     & .7D0, .6D0, .5D0, .45D0, .35D0, .3D0, 0.D0, .5D0, 1.3D0, 2.8D0,
     & 2.3D0, 1.6D0, .9D0, 13*0.D0, 0.9D0, 2.5D0, 3.D0, 2.5D0, 2.3D0,
     & 2.D0, 1.7D0, 1.5D0,1.2D0,.9D0,.6D0,.45D0,.21D0,.2D0,3*0.D0,0.9D0,
     & 2.5D0, 3.D0, 2.5D0, 2.3D0,2.D0,1.7D0,1.5D0,1.2D0,.9D0,.6D0,.45D0,
     & .21D0, .2D0,4*0.D0,1.D0,2.1D0,2.6D0,2.3D0,2.D0,1.8D0,1.7D0,1.4D0,
     & 1.2D0,1.15D0,.9D0,.66D0,.5D0,4*0.D0,1.D0,2.1D0,2.6D0,2.3D0,2.1D0,
     & 1.8D0,1.7D0,1.4D0,1.2D0, 1.15D0, .9D0, .66D0, .5D0, 7*0.D0, .3D0,
     & 2*1.D0, .9D0, .7D0, .4D0, .35D0, .2D0, 9*0.D0, .3D0, 2*1.D0,.9D0,
     & .7D0, .4D0, .35D0, .2D0, 11*0.D0, .1D0, 1.D0, 2.4D0,3.5D0,4.25D0,
     & 4.55D0, 4.85D0, 4.9D0, 9*0.D0, .1D0, 1.D0, 2.4D0, 3.5D0, 4.25D0,
     & 4.55D0, 4.85D0, 4.9D0, 10*0.D0, .2D0, .7D0, 1.6D0, 2.5D0, 2.2D0,
     & 1.71D0, 1.6D0, 10*0.D0, .2D0, .7D0, 1.6D0, 2.5D0, 2.2D0, 1.71D0,
     & 1.6D0, 6*0.D0, 1.4D0, 3.8D0, 5.D0, 4.7D0,4.4D0,4.D0,3.5D0,2.85D0,
     & 2.35D0, 2.01D0, 1.8D0, 6*0.D0, 1.4D0,3.8D0,5.D0,4.7D0,4.4D0,4.D0,
     & 3.5D0,2.85D0,2.35D0,2.01D0,1.8D0,12*0.D0,.1D0,.8D0,2.05D0,3.31D0,
     & 3.5D0, 12*0.D0, .034D0,.2D0,.75D0,1.04D0,1.24D0 /
**** k-  p data                                                        *
      DATA SKMPEL/ 0.D0, 35.D0, 22.D0, 25.D0, 17.D0, 9.D0, 9.5D0, 8.D0,
     &     7.D0, 6.5D0, 6.1D0, 5.D0, 4.8D0, 4.6D0, 4.45D0, 4.3D0, 4.2D0,
     &    0.D0, 8.D0, 3.5D0, 8.D0, 3.D0, 1.9D0, 1.7D0, 1.D0, .9D0, .8D0,
     &    .75D0, .5D0, .42D0, .38D0, .34D0, .25D0, .2D0,
     &    0.D0, 3.D0, 3.2D0, 3.5D0, 1.5D0, 1.4D0, 1.1D0, .6D0, .5D0,
     &    .35D0, .28D0, .25D0, .18D0, .12D0, .1D0, .08D0, .04D0,
     &    0.D0, 8.5D0, 2.4D0, 1.7D0, 1.3D0, 1.3D0, 1.1D0, .5D0,
     &    .4D0, .4D0, .35D0, .3D0, .28D0, .2D0, .16D0, .13D0, .11D0,
     &    0.D0, 7.D0, 4.8D0, 1.4D0, 1.9D0, .9D0, .4D0, .2D0, .13D0,
     &    .1D0, .08D0, .06D0, .04D0, .02D0, .015D0, .01D0, .01D0,
     &    0.D0, 5.5D0, 1.D0, .8D0, .75D0, .32D0, .2D0, .1D0, .09D0,
     &    .08D0, .065D0, .05D0, .04D0, .022D0, .017D0, 2*.01D0/
      DATA SPIKP7 / 0.D0, .56D0, 1.46D0, 3.16D0, 2.01D0, 1.28D0, .74D0,
     & 14*0.D0, 1.13D0, 2.61D0, 2.91D0, 2.58D0, 2.35D0, 2.02D0,
     & 1.91D0, 1.57D0, 1.35D0, 1.29D0, 1.01D0, .74D0, .65D0, 4*0.D0,
     & 1.13D0, 2.61D0, 2.91D0, 2.58D0, 2.35D0, 2.02D0, 1.91D0, 1.57D0,
     & 1.35D0, 1.29D0, 1.01D0, .74D0, .65D0,  3*0.D0, 1.0D0, 3.03D0,
     & 3.36D0, 2.8D0, 2.58D0, 2.24D0, 1.91D0, 1.68D0, 1.35D0, 1.01D0,
     & .67D0, .5D0, .24D0, .23D0, 3*0.D0, 1.0D0, 3.03D0, 3.36D0, 2.8D0,
     & 2.58D0, 2.24D0, 1.91D0, 1.68D0, 1.35D0, 1.01D0, .67D0, .5D0,
     & .24D0, .23D0, 7*0.D0, .34D0, 1.12D0, 1.12D0, 1.01D0, .78D0,
     & .45D0, .39D0, .22D0, .07D0, 0.D0, 7*0.D0, .34D0, 1.12D0, 1.12D0,
     & 1.01D0, .78D0, .45D0, .39D0, .22D0, .07D0, 0.D0, 6*0.D0, 1.71D0,
     & 4.26D0, 5.6D0, 5.57D0, 4.93D0, 4.48D0, 3.92D0, 3.19D0, 2.63D0,
     & 2.25D0, 2.D0, 6*0.D0, 1.71D0, 4.26D0, 5.6D0, 5.57D0, 4.93D0,
     & 4.48D0, 3.92D0, 3.19D0, 2.63D0, 2.25D0, 2.D0, 10*0.D0, .22D0,
     & .8D0, .75D0, 1.D0, 1.3D0, 1.5D0, 1.3D0, 10*0.D0, .22D0, .8D0,
     & .75D0, 1.D0, 1.3D0, 1.5D0, 1.3D0, 13*0.D0, .1D0, .3D0, .7D0,1.D0,
     & 13*0.D0, .1D0, .3D0, .7D0, 1.D0, 9*0.D0, .11D0, 1.72D0, 2.69D0,
     & 3.92D0, 4.76D0, 5.10D0, 5.44D0, 5.3D0, 9*0.D0, .11D0, 1.72D0,
     & 2.69D0, 3.92D0, 4.76D0, 5.1D0, 5.44D0, 5.3D0, 5*0.D0,9.2D0,4.7D0,
     & 1.9D0, 10*0.D0, 2.5D0, 15.D0, 21.5D0, 15.3D0, 3.D0, 1.5D0,
     & 10*0.D0/
***** k- n data                                                        *
      DATA SKMNEL/0.D0, 4.D0, 9.5D0, 20.D0, 13.D0, 9.5D0, 6.D0, 4.4D0,
     &        3.D0, 2.4D0, 2.D0, 1.4D0, 1.2D0, 1.D0, .9D0, .7D0, .6D0,
     &        0.D0, 4.5D0, 6.D0, 5.D0, 2.5D0, 2.D0, 1.7D0, 2.1D0,
     &        1.9D0, .9D0, .5D0, .3D0, .24D0, .2D0, .18D0, .1D0, .09D0,
     &        0.D0, 1.8D0, 2.D0, 1.1D0, .9D0, .5D0, .5D0, .4D0, .4D0,
     &        .2D0, .1D0, .06D0, .05D0, .04D0, .03D0, .02D0, .02D0,
     &        0.D0, 1.5D0, 2.D0, .9D0, 1.1D0, .4D0, .6D0, .7D0, .65D0,
     &       .3D0, .17D0, .1D0, .08D0, .07D0, .06D0, .04D0, .03D0/
      DATA SPIKP8/0.D0, .56D0, 1.29D0, 2.26D0, 1.01D0, .64D0, .37D0,
     &  14*0.D0, 1.13D0, 2.61D0, 2.91D0, 2.58D0, 2.35D0, 2.02D0,
     &  1.91D0, 1.57D0, 1.35D0, 1.29D0, 1.01D0, .74D0, .65D0,
     &  3*0.D0, 1.D0, 3.03D0, 3.36D0, 2.8D0, 2.58D0, 2.24D0,
     &  1.91D0, 1.68D0, 1.35D0, 1.01D0, .67D0, .5D0, .24D0, .23D0,
     &  3*0.D0, 1.D0, 3.03D0, 3.36D0, 2.8D0, 2.58D0, 2.24D0,
     &  1.91D0, 1.68D0, 1.35D0, 1.01D0, .67D0, .5D0, .24D0, .23D0,
     &  7*0.D0, .34D0, 1.12D0, 1.12D0, 1.01D0, .78D0, .45D0,
     &  .39D0, .22D0, .07D0, 0.D0,
     &  6*0.D0, 1.71D0, 4.26D0, 5.6D0, 5.57D0, 4.93D0,
     &  4.48D0, 3.92D0, 3.19D0, 2.63D0, 2.25D0, 2.D0,
     &  10*0.D0, .22D0, .8D0, .75D0, 1.D0, 1.3D0, 1.5D0, 1.3D0,
     &  13*0.D0, .1D0, .3D0, .7D0, 1.D0,
     &  13*0.D0, .1D0, .3D0, .7D0, 1.D0,
     &  9*0.D0, .11D0, 1.72D0, 2.69D0, 3.92D0, 4.76D0,
     &  5.10D0, 5.44D0, 5.3D0,
     &  4*0.D0, 0.00D0, 9.2D0, 4.7D0, 1.9D0, 9*0.D0/
*****  p p data                                                        *
      DATA SPIKP9/ 0.D0, 24.D0, 25.D0, 27.D0, 23.D0, 21.D0, 20.D0,
     &              19.D0, 17.D0, 15.5D0, 14.D0, 13.5D0, 13.D0,
     &              0.D0, 3.6D0, 1.7D0, 10*0.D0,
     &              .0D0, 0.D0, 8.7D0, 17.7D0, 18.8D0, 15.9D0,
     &              11.7D0, 8.D0, 6.D0, 5.3D0, 4.5D0, 3.9D0, 3.5D0,
     &              .0D0, .0D0, 2.8D0, 5.8D0, 6.2D0, 5.1D0, 3.8D0,
     &              2.7D0, 2.1D0, 1.8D0, 1.5D0, 1.3D0, 1.1D0,
     &              5*0.D0, 4.6D0, 10.2D0, 15.1D0,
     &              16.9D0, 16.5D0, 11.D0, 5.5D0, 3.5D0,
     &              10*0.D0, 4.3D0, 7.6D0, 9.D0,
     &              10*0.D0, 1.7D0, 2.6D0, 3.D0,
     &              6*0.D0, .3D0, .6D0, 1.D0, 1.6D0, 1.3D0, .8D0, .6D0,
     &              6*0.D0, .7D0, 1.2D0, 1.8D0, 2.5D0, 1.8D0, 1.3D0,
     &              1.2D0, 10*0.D0, .6D0, 1.4D0, 1.7D0,
     &              10*0.D0, 1.9D0, 4.1D0, 5.2D0/
*****  p n data                                                        *
      DATA SPIKP0/ 0.D0, 24.D0, 25.D0, 27.D0, 23.D0, 21.D0, 20.D0,
     &              19.D0, 17.D0, 15.5D0, 14.D0, 13.5D0, 13.D0,
     &              0.D0, 1.8D0, .2D0,  12*0.D0,
     &              3.2D0, 6.05D0, 9.9D0, 5.1D0,
     &              3.8D0, 2.7D0, 1.9D0, 1.5D0, 1.4D0, 1.3D0, 1.1D0,
     &              2*.0D0, 3.2D0, 6.05D0, 9.9D0, 5.1D0,
     &              3.8D0, 2.7D0, 1.9D0, 1.5D0, 1.4D0, 1.3D0, 1.1D0,
     &              5*0.D0, 4.6D0, 10.2D0, 15.1D0,
     &              16.4D0, 15.2D0, 11.D0, 5.4D0, 3.5D0,
     &              5*0.D0, 4.6D0, 10.2D0, 15.1D0,
     &              16.4D0, 15.2D0, 11.D0, 5.4D0, 3.5D0,
     &              10*0.D0, .7D0, 5.1D0, 8.D0,
     &              10*0.D0, .7D0, 5.1D0, 8.D0,
     &              10*.0D0, .3D0, 2.8D0, 4.7D0,
     &              10*.0D0, .3D0, 2.8D0, 4.7D0,
     &              7*0.D0, 1.2D0, 2.5D0, 3.5D0, 6.D0, 5.3D0, 2.9D0,
     &              7*0.D0, 1.7D0, 3.6D0, 5.4D0, 9.D0, 7.6D0, 4.2D0,
     &              5*0.D0, 7.7D0, 6.1D0, 2.9D0, 5*0.D0/
*   nn - data                                                          *
*                                                                      *
      DATA SPKP V/ 0.D0, 24.D0, 25.D0, 27.D0, 23.D0, 21.D0, 20.D0,
     &              19.D0, 17.D0, 15.5D0, 14.D0, 13.5D0, 13.D0,
     &              0.D0, 3.6D0, 1.7D0, 12*0.D0,
     &              8.7D0, 17.7D0, 18.8D0, 15.9D0,
     &              11.7D0, 8.D0, 6.D0, 5.3D0, 4.5D0, 3.9D0, 3.5D0,
     &              .0D0, .0D0, 2.8D0, 5.8D0, 6.2D0, 5.1D0, 3.8D0,
     &              2.7D0, 2.1D0, 1.8D0, 1.5D0, 1.3D0, 1.1D0,
     &              5*0.D0, 4.6D0, 10.2D0, 15.1D0, 16.9D0, 16.5D0,
     &              11.D0, 5.5D0, 3.5D0,
     &              10*0.D0, 4.3D0, 7.6D0, 9.D0,
     &              10*0.D0, 1.7D0, 2.6D0, 3.D0,
     &              6*0.D0, .3D0, .6D0, 1.D0, 1.6D0, 1.3D0, .8D0, .6D0,
     &              6*0.D0, .7D0, 1.2D0, 1.8D0, 2.5D0, 1.8D0, 1.3D0,
     &              1.2D0, 10*0.D0, .6D0, 1.4D0, 1.7D0,
     &              10*0.D0, 1.9D0, 4.1D0, 5.2D0/
****************   ap - p - data                                       *
      DATA SAPPEL/ 0.D0,  176.D0, 160.D0, 105.D0, 75.D0, 68.D0, 65.D0,
     &  50.D0,  50.D0, 43.D0, 42.D0, 40.5D0, 35.D0, 30.D0, 28.D0,
     &  25.D0,  22.D0, 21.D0, 20.D0, 18.D0, 17.D0,  11*0.D0,
     &  .05D0,  .15D0, .18D0, .2D0, .2D0, .3D0, .4D0, .6D0, .7D0, .85D0,
     &  0.D0,  1.D0, .9D0, .46D0, .3D0, .23D0, .18D0, .16D0, .14D0,
     &  .1D0,  .08D0, .05D0, .02D0, .015D0, 4*.011D0, 3*.005D0,
     &  0.D0,  55.D0, 50.D0, 25.D0, 15.D0, 15.D0, 14.D0, 12.D0,
     &  10.D0,  7.D0, 6.D0, 4.D0, 3.3D0, 2.8D0, 2.4D0, 2.D0, 1.8D0,
     &  1.55D0,  1.3D0, .95D0, .75D0,
     &  0.D0,  3.3D0, 3.D0, 1.5D0, 1.D0, .7D0, .4D0, .35D0, .4D0,
     &  .25D0,  .18D0, .08D0, .04D0, .03D0, .023D0, .016D0, .014D0,
     & .01D0,  .008D0, .006D0, .005D0/
      DATA SPIKPE/0.D0, 215.D0, 193.D0, 170.D0, 148.D0, 113.D0, 97.D0,
     & 84.D0, 78.D0, 68.D0, 64.D0, 61.D0, 46.D0, 36.D0, 31.3D0, 28.5D0,
     & 25.7D0, 22.6D0, 21.4D0, 20.7D0, 19.9D0,
     & 9*0.D0, 2.D0, 2.5D0, .2D0, 19*0.D0, .3D0, 1.4D0, 2.2D0, 1.2D0,
     & 1.1D0, 1.D0, .8D0, .6D0, .5D0, .4D0, .3D0, 10*0.D0, .3D0, 1.4D0,
     & 2.2D0, 1.2D0, 1.1D0, 1.D0, .8D0, .6D0, .5D0, .4D0, .3D0, 10*0.D0,
     & .3D0, 1.4D0, 2.2D0, 1.2D0, 1.1D0, 1.D0, .8D0, .6D0, .5D0, .4D0,
     & .3D0, 10*0.D0, .3D0, 1.4D0, 2.2D0, 1.2D0, 1.1D0, 1.D0, .8D0,
     & .6D0, .5D0, .4D0, .3D0, 9*0.D0, .6D0, 2.5D0, 5.D0, 5.2D0, 5.1D0,
     & 5.4D0, 5.8D0, 2.8D0, 2.1D0, 1.8D0, 1.6D0, 1.2D0, 13*0.D0, 1.3D0,
     & 1.5D0, 2.D0, 2.5D0, 2.5D0, 2.3D0, 1.8D0, 1.4D0, 13*0.D0, 1.3D0,
     & 1.5D0, 2.D0, 2.5D0, 2.5D0, 2.3D0, 1.8D0, 1.4D0, 13*0.D0, 1.3D0,
     & 1.5D0, 2.D0, 2.5D0, 2.5D0, 2.3D0, 1.8D0, 1.4D0, 13*0.D0, 1.3D0,
     & 1.5D0, 2.D0, 2.5D0, 2.5D0, 2.3D0, 1.8D0, 1.4D0, 14*0.D0, .2D0,
     & .5D0, 1.1D0, 1.6D0, 1.4D0, 1.1D0, .9D0, 14*0.D0, .2D0, .5D0,
     & 1.1D0, 1.6D0, 1.4D0, 1.1D0, .9D0, 14*0.D0, .2D0, .5D0, 1.1D0,
     & 1.6D0, 1.4D0, 1.1D0, .9D0, 14*0.D0, .2D0, .5D0, 1.1D0, 1.6D0,
     & 1.4D0, 1.1D0, .9D0, 17*0.D0, .3D0, 1.6D0, 2.6D0, 3.6D0, 17*0.D0,
     & .3D0, 1.6D0, 2.6D0, 3.6D0, 17*0.D0, .3D0, 1.6D0, 2.6D0,
     & 3.6D0, 17*0.D0, .3D0, 1.6D0, 2.6D0, 3.6D0 /
****************   ap - n - data                                       *
      DATA SAPNEL/
     & 0.D0,  176.D0, 160.D0, 105.D0, 75.D0,  68.D0, 65.D0,
     & 50.D0, 50.D0,  43.D0,  42.D0,  40.5D0, 35.D0, 30.D0,  28.D0,
     & 25.D0, 22.D0,  21.D0,  20.D0,  18.D0,  17.D0, 11*0.D0,
     & .05D0, .15D0, .18D0,  .2D0,    .2D0,  .3D0,  .4D0,   .6D0,  .7D0,
     & .85D0,  0.D0,  1.D0,  .9D0,    .46D0, .3D0,  .23D0, .18D0, .16D0,
     & .14D0,  .1D0, .08D0, .05D0,    .02D0, .015D0, 4*.011D0, 3*.005D0,
     & 0.D0,  3.3D0,  3.D0, 1.5D0,     1.D0, .7D0,  .4D0,  .35D0, .4D0,
     & .25D0, .18D0, .08D0, .04D0,    .03D0, .023D0, .016D0, .014D0,
     & .01D0, .008D0, .006D0, .005D0 /
       DATA SPIKPZ/ 0.D0, 215.D0, 193.D0, 170.D0, 148.D0, 113.D0, 97.D0,
     &  84.D0, 78.D0, 68.D0, 64.D0, 61.D0, 46.D0, 36.D0, 31.3D0, 28.5D0,
     & 25.7D0, 22.6D0, 21.4D0, 20.7D0, 19.9D0, 9*0.D0, 2.4D0, .2D0,
     & 20*0.D0, 1.8D0, 2.8D0, 3.6D0, 2.3D0, 1.8D0, 1.5D0, 1.3D0, 1.D0,
     & .7D0, .5D0, .3D0, 10*0.D0, 1.8D0, 2.8D0, 3.6D0, 2.3D0, 1.8D0,
     & 1.5D0, 1.3D0, 1.D0, .7D0, .5D0, .3D0, 10*0.D0, 1.8D0, 2.8D0,
     & 3.6D0, 2.3D0, 1.8D0, 1.5D0, 1.3D0, 1.D0, .7D0, .5D0, .3D0,
     & 10*0.D0, 1.8D0, 2.8D0, 3.6D0, 2.3D0, 1.8D0, 1.5D0, 1.3D0, 1.D0,
     & .7D0, .5D0, .3D0, 13*0.D0, 5.2D0, 8.7D0, 11.4D0, 14.D0, 11.9D0,
     & 7.6D0, 6.D0, 5.D0, 13*0.D0, 5.2D0, 8.7D0, 11.4D0, 14.D0, 11.9D0,
     & 7.6D0, 6.D0, 5.D0, 18*0.D0, 1.D0, 4.9D0, 8.5D0, 18*0.D0, 1.D0,
     & 4.9D0, 8.5D0,  15*0.D0, 1.9D0, 2.3D0, 4.D0, 6.5D0, 5.2D0, 3.4D0,
     & 15*0.D0, 1.9D0, 2.3D0, 4.D0, 6.5D0, 5.2D0, 3.4D0, 15*0.D0, 1.9D0,
     & 2.3D0, 4.D0, 6.5D0, 5.2D0, 3.4D0 /
*                                                                      *
*                                                                      *
****************   an - p - data                                       *
*                                                                      *
      DATA SANPEL/
     & 0.D0,  176.D0, 160.D0, 105.D0, 75.D0, 68.D0, 65.D0, 50.D0,
     & 50.D0, 43.D0,  42.D0,  40.5D0, 35.D0, 30.D0, 28.D0,
     & 25.D0, 22.D0,  21.D0,  20.D0,  18.D0, 17.D0, 11*0.D0, .05D0,
     & .15D0, .18D0,   .2D0,   .2D0,   .3D0,  .4D0, .6D0,   .7D0, .85D0,
     & 0.D0,   1.D0,   .9D0,  .46D0,  .3D0,  .23D0, .18D0, .16D0, .14D0,
     & .1D0,  .08D0,  .05D0,  .02D0, .015D0, 4*.011D0, 3*.005D0,
     & 0.D0,  3.3D0,  3.D0, 1.5D0, 1.D0, .7D0, .4D0, .35D0, .4D0, .25D0,
     & .18D0, .08D0, .04D0, .03D0, .023D0, .016D0, .014D0,
     & .01D0, .008D0, .006D0, .005D0 /
      DATA SPIKPF/ 0.D0, 215.D0, 193.D0, 170.D0, 148.D0, 113.D0, 97.D0,
     & 84.D0, 78.D0, 68.D0, 64.D0, 61.D0, 46.D0, 36.D0, 31.3D0, 28.5D0,
     & 25.7D0, 22.6D0, 21.4D0, 20.7D0, 19.9D0, 9*0.D0, 2.4D0, .2D0,
     & 20*0.D0, 1.8D0, 2.8D0, 3.6D0, 2.3D0, 1.8D0, 1.5D0, 1.3D0, 1.D0,
     & .7D0, .5D0, .3D0, 10*0.D0, 1.8D0, 2.8D0, 3.6D0, 2.3D0, 1.8D0,
     & 1.5D0, 1.3D0, 1.D0, .7D0, .5D0, .3D0, 10*0.D0, 1.8D0, 2.8D0,
     & 3.6D0, 2.3D0, 1.8D0, 1.5D0, 1.3D0, 1.D0, .7D0, .5D0, .3D0,
     & 10*0.D0, 1.8D0, 2.8D0, 3.6D0, 2.3D0, 1.8D0, 1.5D0, 1.3D0, 1.D0,
     & .7D0, .5D0, .3D0, 13*0.D0, 5.2D0, 8.7D0, 11.4D0, 14.D0, 11.9D0,
     & 7.6D0, 6.D0, 5.D0, 13*0.D0, 5.2D0, 8.7D0, 11.4D0, 14.D0, 11.9D0,
     & 7.6D0, 6.D0, 5.D0, 18*0.D0, 1.D0, 4.9D0, 8.5D0, 18*0.D0, 1.D0,
     & 4.9D0, 8.5D0, 15*0.D0, 1.9D0, 2.3D0, 4.D0, 6.5D0, 5.2D0, 3.4D0,
     & 15*0.D0, 1.9D0, 2.3D0, 4.D0, 6.5D0, 5.2D0, 3.4D0, 15*0.D0, 1.9D0,
     & 2.3D0, 4.D0, 6.5D0, 5.2D0, 3.4D0 /
****  ko - n - data                                                    *
      DATA SPKP15/0.D0, 20.D0, 14.D0, 12.D0, 11.5D0, 10.D0, 8.D0, 7.D0,
     &      6.D0, 5.5D0, 5.3D0, 5.D0, 4.5D0, 4.4D0, 3.8D0, 3.D0, 2.8D0,
     &      0.D0, .5D0, 1.15D0, 2.D0, 1.3D0, .8D0, .45D0, 10*0.D0,
     &    3*0.D0, 0.9D0, 2.5D0, 3.D0, 2.5D0, 2.3D0, 2.D0, 1.7D0,
     &     1.5D0, 1.2D0, .9D0, .6D0, .45D0, .21D0, .2D0,
     &    3*0.D0, 0.9D0, 2.5D0, 3.D0, 2.5D0, 2.3D0, 2.D0, 1.7D0,
     &     1.5D0, 1.2D0, .9D0, .6D0, .45D0, .21D0, .2D0,
     &    4*0.D0, 1.D0, 2.1D0, 2.6D0, 2.3D0, 2.1D0, 1.8D0, 1.7D0,
     &     1.4D0, 1.2D0, 1.05D0, .9D0, .66D0,  .5D0,
     &    7*0.D0, .3D0, 1.D0, 1.D0, .9D0, .7D0, .4D0, .30D0, .2D0,
     &   11*0.D0, .1D0, 1.D0, 2.2D0, 3.5D0, 4.20D0, 4.55D0,
     &    4.85D0, 4.9D0,
     &   10*0.D0, .2D0, .7D0, 1.6D0, 2.5D0, 2.2D0, 1.71D0, 1.6D0,
     &    6*0.D0, 1.4D0, 3.8D0, 5.D0, 4.7D0, 4.4D0, 4.D0, 3.5D0,
     &    2.85D0, 2.35D0, 2.01D0, 1.8D0,
     &   12*0.D0, .1D0, .8D0, 2.05D0, 3.31D0, 3.5D0,
     &   12*0.D0, .034D0, .20D0, .75D0, 1.04D0, 1.24D0  /
**** ako - p - data                                                    *
      DATA SPKP16/ 0.D0, 4.D0, 9.5D0, 20.D0, 13.D0, 9.5D0, 6.D0, 4.4D0,
     & 3.D0, 2.4D0, 2.D0, 1.4D0, 1.2D0, 1.D0, .9D0, .7D0, .6D0, 0.D0,
     & 4.5D0, 6.D0, 5.D0, 2.5D0, 2.D0, 1.7D0, 2.1D0, 1.9D0, .9D0, .5D0,
     & .3D0, .24D0, .2D0, .18D0, .1D0, .09D0, 0.D0, 1.8D0, 2.D0, 1.1D0,
     & .9D0, .5D0, .5D0, .4D0, .4D0, .2D0, .1D0, .06D0, .05D0, .04D0,
     & .03D0, .02D0, .02D0, 0.D0, 1.5D0, 2.D0, .9D0, 1.1D0, .4D0, .6D0,
     & .7D0, .65D0, .3D0, .17D0, .1D0, .08D0, .07D0, .06D0, .04D0,
     & .03D0, 0.D0, .56D0, 1.29D0, 2.26D0, 1.01D0, .64D0, .37D0,
     & 14*0.D0, 1.13D0, 2.61D0, 2.91D0, 2.58D0, 2.35D0, 2.02D0, 1.91D0,
     & 1.57D0, 1.35D0, 1.29D0, 1.01D0, .74D0, .65D0, 3*0.D0, 1.0D0,
     & 3.03D0, 3.36D0, 2.8D0, 2.58D0, 2.24D0, 1.91D0, 1.68D0, 1.35D0,
     & 1.01D0, .67D0, .5D0, .24D0, .23D0, 3*0.D0, 1.0D0, 3.03D0, 3.36D0,
     & 2.8D0, 2.58D0, 2.24D0, 1.91D0, 1.68D0, 1.35D0, 1.01D0, .67D0,
     & .5D0, .24D0, .23D0, 7*0.D0, .34D0, 1.12D0, 1.12D0, 1.01D0, .78D0,
     & .45D0, .39D0, .22D0, .07D0, 7*0.D0, 1.71D0, 4.26D0, 5.6D0,5.57D0,
     & 4.93D0, 4.48D0, 3.92D0, 3.19D0, 2.63D0, 2.25D0, 2.D0, 10*0.D0,
     & .22D0, .8D0, .75D0, 1.D0, 1.3D0, 1.5D0, 1.3D0, 13*0.D0, .1D0,
     & .3D0, .7D0, 1.D0, 13*0.D0, .1D0, .3D0, .7D0, 1.D0, 9*0.D0, .11D0,
     & 1.72D0, 2.69D0, 3.92D0, 4.76D0, 5.10D0, 5.44D0, 5.3D0, 5*0.D0,
     & 9.2D0, 4.7D0, 1.9D0, 9*0.D0, .0D0,2.5D0,15.D0,
     & 21.5D0, 15.3D0, 3.D0, 1.5D0, 10*0.D0 /
      DATA NURELN/9, 12, 5*0, 10, 14, 3*0, 1, 3, 5, 7, 6*0, 2, 6, 16,
     & 5*0, 10, 13, 5*0, 11, 12, 3*0, 2, 4, 6, 8, 6*0, 3, 15, 7, 5*0 /
*=                                               end*block.blkdt3      *
      END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      BLOCK DATA REACCH
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C*** REACTION CHANNEL CROSS SECTION DATA
C     INTEGER * 2
C    *                            NRKPI,NRKKC,NRKK0,NRKP,NRKN
C    *          ,NURE
      COMMON /DREACC/
     *UMOPI(  92),UMOKC(  68),UMOP( 39),UMON( 63),UMOK0( 34),
     *PLAPI(  92),PLAKC(  68),PLAP( 39),PLAN( 63),PLAK0( 34),
     *SIIN(296),
     *SPIKP1( 315),SPIKP U(278),SPIKP V(372),
     *SPIKP W(278),SPIKP X(372),SPIKP 4(315),
     *SPIKP 5(187),SPIKP 6(306),
     *S KMPEL(102),SPIKP 7(289),S KMNEL( 68),SPIKP 8(187),
     *SPIKP 9(143),SPIKP 0(169),SPKP V(143),
     *S APPEL(105),SPIKP E(399),S APNEL( 84),SPIKP Z(273),
     *S ANPEL( 84),SPIKP F(273),
     *SPKP15(187),SPKP16(255),
     *NRKPI( 164),NRKKC( 134),NRKP( 70),NRKN( 116),NRKK0( 52),
     *NURE(60)
C111111111111111111111111111111111111111111111111111111111111111111
*
**** pi- p data                                                        *
**** pi+ n data                                                        *
      DATA PLAPI / 0.D0, .3D0, .5D0, .6D0, .7D0, .8D0, .9D0, .95D0,1.D0,
     & 1.15D0, 1.3D0, 1.5D0, 1.6D0, 1.8D0, 2.D0, 2.3D0, 2.5D0, 2.8D0,
     & 3.D0, 3.5D0, 4.D0, 0.D0, .285D0, .4D0, .45D0, .5D0, .6D0, .7D0,
     & .75D0, .8D0, .85D0, .9D0, 1.D0, 1.15D0, 1.3D0, 1.5D0, 1.6D0,
     & 1.8D0, 2.D0, 2.3D0, 2.5D0, 2.8D0, 3.D0, 3.5D0, 4.D0, 4.5D0, 0.D0,
     & .285D0, .4D0, .45D0, .5D0, .6D0, .7D0, .75D0, .8D0, .85D0, .9D0,
     & 1.D0, 1.15D0, 1.3D0, 1.5D0, 1.6D0, 1.8D0, 2.D0, 2.3D0, 2.5D0,
     & 2.8D0, 3.D0, 3.5D0, 4.D0, 4.5D0, 0.D0, .3D0, .5D0, .6D0, .7D0,
     & .8D0, .9D0, .95D0, 1.D0, 1.15D0, 1.3D0, 1.5D0, 1.6D0, 1.8D0,
     & 2.D0, 2.3D0, 2.5D0, 2.8D0, 3.D0, 3.5D0, 4.D0 /
      DATA PLAKC /
     &   0.D0,  .58D0,   .8D0, 1.01D0, 1.23D0, 1.45D0, 1.68D0, 1.94D0,
     & 2.18D0, 2.42D0, 2.68D0, 2.96D0, 3.24D0,
     & 3.51D0, 3.84D0, 4.16D0, 4.49D0,
     &   0.D0,  .58D0,   .8D0, 1.01D0, 1.23D0, 1.45D0, 1.68D0, 1.94D0,
     & 2.18D0, 2.42D0, 2.68D0, 2.96D0, 3.24D0,
     & 3.51D0, 3.84D0, 4.16D0, 4.49D0,
     &   0.D0,  .58D0,   .8D0, 1.01D0, 1.23D0, 1.45D0, 1.68D0, 1.94D0,
     & 2.18D0, 2.42D0, 2.68D0, 2.96D0, 3.24D0,
     & 3.51D0, 3.84D0, 4.16D0, 4.49D0,
     &   0.D0,  .58D0,   .8D0, 1.01D0, 1.23D0, 1.45D0, 1.68D0, 1.94D0,
     & 2.18D0, 2.42D0, 2.68D0, 2.96D0, 3.24D0,
     & 3.51D0, 3.84D0, 4.16D0, 4.49D0/
      DATA PLAK0 /
     &   0.D0,  .58D0,   .8D0, 1.01D0, 1.23D0, 1.45D0, 1.68D0, 1.94D0,
     & 2.18D0, 2.42D0, 2.68D0, 2.96D0, 3.24D0,
     & 3.51D0, 3.84D0, 4.16D0, 4.49D0,
     &   0.D0,  .58D0,   .8D0, 1.01D0, 1.23D0, 1.45D0, 1.68D0, 1.94D0,
     & 2.18D0, 2.42D0, 2.68D0, 2.96D0, 3.24D0,
     & 3.51D0, 3.84D0, 4.16D0, 4.49D0/
*                 pp   pn   np   nn                                    *
      DATA PLAP /
     &   0.D0, 1.06D0, 1.34D0, 1.63D0, 1.92D0, 2.2D0, 2.5D0,2.8D0,3.1D0,
     & 3.43D0, 3.75D0, 4.07D0, 4.43D0,
     &   0.D0, 1.06D0, 1.34D0, 1.63D0, 1.92D0, 2.2D0, 2.5D0,2.8D0,3.1D0,
     & 3.43D0, 3.75D0, 4.07D0, 4.43D0,
     &   0.D0, 1.06D0, 1.34D0, 1.63D0, 1.92D0, 2.2D0, 2.5D0,2.8D0,3.1D0,
     & 3.43D0, 3.75D0, 4.07D0, 4.43D0 /
*    app   apn   anp   ann                                             *
      DATA PLAN /
     &  0.D0,   1.D-3,   .1D0,   .2D0,   .3D0,  .4D0,  .5D0, .6D0,
     & .74D0,  1.06D0, 1.34D0, 1.63D0, 1.92D0, 2.2D0, 2.5D0,2.8D0,3.1D0,
     & 3.43D0, 3.75D0, 4.07D0, 4.43D0,
     &  0.D0,   1.D-3,   .1D0,   .2D0,   .3D0,  .4D0,  .5D0, .6D0,
     & .74D0,  1.06D0, 1.34D0, 1.63D0, 1.92D0, 2.2D0, 2.5D0,2.8D0,3.1D0,
     & 3.43D0, 3.75D0, 4.07D0, 4.43D0,
     &  0.D0,   1.D-3,   .1D0,   .2D0,   .3D0,  .4D0,  .5D0, .6D0,
     & .74D0,  1.06D0, 1.34D0, 1.63D0, 1.92D0, 2.2D0, 2.5D0,2.8D0,3.1D0,
     & 3.43D0, 3.75D0, 4.07D0, 4.43D0  /
      DATA SIIN / 296*0.D0 /
      DATA UMOPI/ 1.08D0,1.233D0,1.302D0,1.369D0,1.496D0,
     & 1.557D0,1.615D0,1.6435D0,
     & 1.672D0,1.753D0,1.831D0,1.930D0,1.978D0,2.071D0,2.159D0,
     & 2.286D0,2.366D0,2.482D0,2.56D0,
     & 2.735D0,2.90D0,
     &             1.08D0,1.222D0,1.302D0,1.3365D0,1.369D0,1.434D0,
     & 1.496D0,1.527D0,1.557D0,
     & 1.586D0,1.615D0,1.672D0,1.753D0,1.831D0,1.930D0,1.978D0,
     & 2.071D0,2.159D0,2.286D0,2.366D0,
     & 2.482D0,2.560D0,2.735D0,2.90D0,3.06D0,
     &             1.08D0,1.222D0,1.302D0,1.3365D0,1.369D0,1.434D0,
     & 1.496D0,1.527D0,1.557D0,
     & 1.586D0,1.615D0,1.672D0,1.753D0,1.831D0,1.930D0,1.978D0,
     & 2.071D0,2.159D0,2.286D0,2.366D0,
     & 2.482D0,2.560D0,2.735D0,2.90D0,3.06D0,
     &                   1.08D0,1.233D0,1.302D0,1.369D0,1.496D0,
     & 1.557D0,1.615D0,1.6435D0,
     & 1.672D0,1.753D0,1.831D0,1.930D0,1.978D0,2.071D0,2.159D0,
     & 2.286D0,2.366D0,2.482D0,2.56D0,
     &  2.735D0, 2.90D0/
      DATA UMOKC/ 1.44D0,
     &  1.598D0,1.7D0,1.8D0,1.9D0,2.0D0,2.1D0,2.2D0,2.3D0,2.4D0,2.5D0,
     & 2.6D0,2.7D0,2.8D0,2.9D0,3.0D0,
     & 3.1D0,1.44D0,
     &  1.598D0,1.7D0,1.8D0,1.9D0,2.0D0,2.1D0,2.2D0,2.3D0,2.4D0,2.5D0,
     & 2.6D0,2.7D0,2.8D0,2.9D0,3.0D0,
     & 3.1D0,1.44D0,
     &  1.598D0,1.7D0,1.8D0,1.9D0,2.0D0,2.1D0,2.2D0,2.3D0,2.4D0,2.5D0,
     & 2.6D0,2.7D0,2.8D0,2.9D0,3.0D0,
     & 3.1D0,1.44D0,
     &  1.598D0,1.7D0,1.8D0,1.9D0,2.0D0,2.1D0,2.2D0,2.3D0,2.4D0,2.5D0,
     & 2.6D0,2.7D0,2.8D0,2.9D0,3.0D0,
     &  3.1D0/
      DATA UMOK0/ 1.44D0,
     &  1.598D0,1.7D0,1.8D0,1.9D0,2.0D0,2.1D0,2.2D0,2.3D0,2.4D0,2.5D0,
     & 2.6D0,2.7D0,2.8D0,2.9D0,3.0D0,
     & 3.1D0,1.44D0,
     &  1.598D0,1.7D0,1.8D0,1.9D0,2.0D0,2.1D0,2.2D0,2.3D0,2.4D0,2.5D0,
     & 2.6D0,2.7D0,2.8D0,2.9D0,3.0D0,
     &  3.1D0/
*                 pp   pn   np   nn                                    *
      DATA UMOP/
     & 1.88D0,2.102D0,2.2D0,2.3D0,2.4D0,2.5D0,2.6D0,2.7D0,2.8D0,2.9D0,
     & 3.D0,3.1D0,3.2D0,
     & 1.88D0,2.102D0,2.2D0,2.3D0,2.4D0,2.5D0,2.6D0,2.7D0,2.8D0,2.9D0,
     & 3.D0,3.1D0,3.2D0,
     & 1.88D0,2.102D0,2.2D0,2.3D0,2.4D0,2.5D0,2.6D0,2.7D0,2.8D0,2.9D0,
     & 3.D0,3.1D0,3.2D0/
*    app   apn   anp   ann                                             *
      DATA UMON /
     & 1.877D0,1.87701D0,1.879D0,1.887D0,1.9D0,1.917D0,1.938D0,1.962D0,
     & 2.D0,2.102D0,2.2D0,2.3D0,2.4D0,2.5D0,2.6D0,2.7D0,2.8D0,2.9D0,
     & 3.D0,3.1D0,3.2D0,
     & 1.877D0,1.87701D0,1.879D0,1.887D0,1.9D0,1.917D0,1.938D0,1.962D0,
     & 2.D0,2.102D0,2.2D0,2.3D0,2.4D0,2.5D0,2.6D0,2.7D0,2.8D0,2.9D0,
     & 3.D0,3.1D0,3.2D0,
     & 1.877D0,1.87701D0,1.879D0,1.887D0,1.9D0,1.917D0,1.938D0,1.962D0,
     & 2.D0,2.102D0,2.2D0,2.3D0,2.4D0,2.5D0,2.6D0,2.7D0,2.8D0,2.9D0,
     &  3.D0,3.1D0,3.2D0/

C111111111111111111111111111111111111111111111111111111111111111111
C111111111111111111111111111111111111111111111111111111111111111111
C*** REACTION CAHNNEL STATE PARTICLES
      DATA NRK PI/
     *13,1,15,21,81,0,
     *13,54,23,53,13,63,13,58,23,57,13,65,1,32,53,31,54,32,53,33,53,35,
     *63,32,
     *13,8,23,1,17,15,21,24,22,15,82,0,
     *61,0,13,55,23,54,14,53,13,64,
     *23,63,13,59,23,58,14,57,13,66,23,65,1,31,8,32,1,33,1,35,54,31,55,
     *32,54,33,53,34,54,35,
     *14,1,23,8,17,24,20,15,22,24,83,0,
     *62,0,14,54,23,55,13,56,14,63,
     *23,64,14,58,23,59,13,60,14,65,23,66,8,31,1,34,8,33,8,35,55,31,54,
     *34,55,33,56,32,55,35,
     *14,8,24,20,84,0,
     *14,55,23,56,14,64,14,59,23,60,14,66,8,34,56,31,55,34,56,33,56,35,
     *64,34
     F/
      DATA NRK KC/
     *15,1,89,0,
     *24,53,15,54,1,36,1,40,1,44,36,63,15,63,45,53,44,54,
     *15,8,24,1,91,0,
     *24,54,15,55,8,36,1,37,8,40,1,41,8,44,1,45,36,64,37,63,15,64,24,63,
     *45,54,44,55,93,0,
     *16,1,25,8,  17,23,21,14,  20,13,22,23,  90,0,
     *38,1,39,8,16,54,25,55,1,42,8,43,16,63,25,64,39,64,38,63,46,54,
     *47,55,8,47,1,46,52,0,51,0,
     *16,8,17,14,20,23,22,14,92,0,
     *8,38,16,55,25,56,8,42,16,64,38,64,46,55,47,56,8,46,94,0
     */
C
C   K0 P   K0 N   AK0 P   AK/ N
C
      DATA NRK K0/
     *24,8,
     *106,0,15,56,24,55,37,8,41,8,45,8,37,64,24,64,44,56,45,55,
     *25,1,17,13,  22,13,21,23,
     *107,0,39,1,25,54,16,53,43,1,25,63,39,63,47,54,46,53,47,1,103,0/
C   PP  PN   NP   NN
      DATA NRK P/1,1,85,0,
     F8,53,1,54,1,63,8,57,1,58,54,54,53,55,63,54,64,53,
     *1,8,86,0,
     F8,54,1,55,8,63,1,64,8,58,1,59,64,54,63,55,54,55,53,56,77,0,
     *8,8,
     *95,0,8,55,1,56,8,64,8,59,1,60,55,55,54,56,64,55,63,56/
C     APP   APN   ANP   ANN
      DATA NRK N/ 1,2,17,18,15,16,8,9,13,14,99,0,87,0, 1,68,8,69,2,54,9,
     +55,102,0, 2,63,9,64,1,75,8,76,53,67,54,68, 55,69,56,70,63,68,64,
     +69,75,54,76,55, 2,8,18,20, 16,24,14,23, 101,0,88,0, 2,55,9,56,1,
     +67,8,68,2,64,8,75,2,59,8,72,68,55,67,54,69,56, 1,9,18,21,15,25,13,
     +23,100,0, 96,0,2,53,9,54,1,69,8,70,1,76,9,63,1,73,9,58,55,70,53,
     +68,54,69/
C222222222222222222222222222222222222222222222222222222222222222222222 
**** channel cross section                                             *
      DATA SPIKP1/ 0.D0, 300.D0, 40.D0, 20.D0, 13.D0,8.5D0,8.D0, 9.5D0,
     & 12.D0,14.D0,15.5D0,20.D0,17.D0,13.D0,10.D0,9.D0,8.5D0,8.D0,7.8D0,
     & 7.3D0, 6.7D0, 9*0.D0,.23D0,.35D0,.7D0,.52D0,.4D0,.3D0,.2D0,.15D0,
     & .13D0, .11D0, .09D0, .07D0, 0.D0, .033D0,.8D0,1.35D0,1.35D0,.5D0,
     & 15*0.D0, 3*0.D0,.00D0,0.80D0,2.2D0,3.6D0,4.6D0,4.7D0,3.5D0,2.4D0,
     &1.8D0,1.4D0,.75D0,.47D0,.25D0,.13D0,.08D0,6*0.D0,0.D0,1.2D0,3.3D0,
     & 5.4D0,6.9D0,7.3D0,5.3D0,3.6D0,2.7D0,2.2D0,1.1D0,.73D0,.4D0,.22D0,
     & .12D0,9*0.D0,.0D0,0.D0,2.0D0,4.4D0,6.8D0,9.9D0,7.9D0,6.0D0,3.8D0,
     &2.5D0,2.D0,1.4D0,1.D0,.6D0,.35D0,10*0.D0,.25D0,.55D0,.75D0,1.25D0,
     & 1.9D0,2.D0,1.8D0,1.5D0,1.25D0,1.D0,.8D0,6*0.D0,4*0.D0,.4D0,.85D0,
     & 1.1D0, 1.85D0, 2.8D0, 3.D0,2.7D0,2.2D0,1.85D0,1.5D0,1.2D0,6*0.D0,
     & 6*0.D0, .5D0, 1.2D0, 1.7D0, 3.4D0, 5.2D0, 6.4D0, 6.1D0, 5.6D0,
     & 5.2D0, 6*0.D0, 2*0.D0, .0D0, 1.D0, 3.3D0, 5.2D0, 4.45D0, 3.6D0,
     & 2.75D0, 1.9D0, 1.65D0, 1.3D0, .95D0, .6D0, .45D0, 6*0.D0, 3*0.D0,
     & .0D0, .45D0, 1.4D0, 1.5D0, 1.1D0, .85D0, .5D0, .3D0, .2D0, .15D0,
     & 8*0.D0, 5*0.D0, .0D0, .0D0, .6D0, .8D0, .95D0, .8D0, .7D0, .6D0,
     & .5D0, .4D0, 6*0.D0, 5*0.D0, .0D0, .00D0, .85D0, 1.2D0, 1.4D0,
     & 1.2D0, 1.05D0, .9D0, .7D0, .55D0, 6*0.D0, 5*0.D0, .0D0, .00D0,
     & 1.D0, 1.5D0, 3.5D0, 4.15D0, 3.7D0, 2.7D0, 2.3D0, 1.75D0, 6*0.D0,
     & 10*0.D0, .5D0, 2.0D0, 3.3D0, 5.4D0, 7.D0 /
**** pi+ n data                                                        *
      DATA SPIKPU/   0.D0, 25.D0, 13.D0,  11.D0, 10.5D0, 14.D0,  20.D0,
     & 20.D0, 16.D0, 14.D0, 19.D0, 28.D0, 17.5D0, 13.5D0, 12.D0, 10.5D0,
     & 10.D0, 10.D0, 9.5D0,  9.D0, 8.D0, 7.5D0, 7.D0, 6.5D0, 6.D0, 0.D0,
     & 48.D0, 19.D0, 15.D0, 11.5D0, 10.D0, 8.D0, 6.5D0,   5.5D0,  4.8D0,
     & 4.2D0, 7.5D0, 3.4D0,  2.5D0, 2.5D0, 2.1D0, 1.4D0,   1.D0,   .8D0,
     &  .6D0, .46D0,  .3D0, .2D0, .15D0, .13D0, 11*0.D0,  .95D0,  .65D0,
     & .48D0, .35D0,  .2D0, .18D0, .17D0, .16D0,  .15D0,   .1D0,  .09D0,
     & .065D0, .05D0, .04D0, 12*0.D0, .2D0, .25D0, .25D0,  .2D0,   .1D0,
     & .08D0, .06D0, .045D0,   .03D0, .02D0, .01D0,      .005D0, .003D0,
     & 12*0.D0, .3D0, .24D0,   .18D0, .15D0, .13D0,  .12D0, .11D0, .1D0,
     & .09D0,  .08D0, .05D0,   .04D0, .03D0,  0.D0, 0.16D0, .7D0, 1.3D0,
     & 3.1D0,  4.5D0,  2.D0, 18*0.D0, 3*.0D0,  0.D0, 0.D0, 4.0D0, 11.D0,
     & 11.4D0, 10.3D0, 7.5D0, 6.8D0, 4.75D0, 2.5D0,  1.5D0, .9D0, .55D0,
     &  .35D0, 13*0.D0, .1D0, .34D0, .5D0, .8D0, 1.1D0,   2.25D0, 3.3D0,
     & 2.3D0, 1.6D0, .95D0, .45D0, .28D0, .15D0, 10*0.D0, 2*0.D0, .17D0,
     & .64D0,  1.D0, 1.5D0, 2.1D0, 4.25D0, 6.2D0,  4.4D0,   3.D0, 1.8D0,
     &  .9D0, .53D0, .28D0,      10*0.D0, 2*0.D0,  .25D0,  .82D0,
     & 1.3D0, 1.9D0, 2.8D0, 5.5D0 , 8.D0,  5.7D0, 3.9D0, 2.35D0, 1.15D0,
     & .69D0, .37D0, 10*0.D0,     7*0.D0,   .0D0, .34D0,  1.5D0, 3.47D0,
     & 5.87D0, 6.23D0, 4.27D0, 2.6D0, 1.D0, .6D0,  .3D0,  .15D0, 6*0.D0/
*
      DATA SPIKPV/ 7*0.D0, .00D0, .16D0, .75D0, 1.73D0, 2.93D0, 3.12D0,
     & 2.13D0, 1.3D0, .5D0, .3D0, .15D0, .08D0, 6*0.D0, 10*0.D0, .2D0,
     & .6D0, .92D0, 2.4D0, 4.9D0, 6.25D0, 5.25D0, 3.5D0, 2.15D0, 1.4D0,
     & 1.D0, .7D0, 13*0.D0, .13D0, .4D0, .62D0, 1.6D0, 3.27D0, 4.17D0,
     & 3.5D0, 2.33D0, 1.43D0, .93D0, .66D0, .47D0, 13*0.D0, .07D0, .2D0,
     & .31D0, .8D0, 1.63D0, 2.08D0, 1.75D0, 1.17D0, .72D0, .47D0, .34D0,
     & .23D0, 17*0.D0, .33D0, 1.D0, 1.8D0, 2.67D0, 5.33D0, 6.D0, 5.53D0,
     & 5.D0, 17*0.D0, .17D0, .5D0, .9D0, 1.83D0, 2.67D0, 3.0D0, 2.77D0,
     & 2.5D0, 3*0.D0, 3*0.D0, 1.D0, 3.3D0, 2.8D0, 2.5D0, 2.3D0, 1.8D0,
     & 1.5D0, 1.1D0, .8D0, .7D0, .55D0, .3D0, 10*0.D0, 9*0.D0, .1D0,
     & .4D0, 1.D0, 1.4D0, 2.2D0, 2.5D0, 2.2D0, 1.65D0, 1.35D0, 1.1D0,
     & .8D0, .6D0, .4D0, 12*0.D0, .15D0, .6D0, 1.5D0, 2.1D0, 3.3D0,
     & 3.8D0, 3.3D0, 2.45D0, 2.05D0, 1.65D0, 1.2D0, .9D0, .6D0, 3*0.D0,
     & 9*0.D0, .10D0, .2D0, .5D0, .7D0, 1.3D0, 1.55D0, 1.9D0, 1.8D0,
     & 1.55D0, 1.35D0, 1.15D0, .95D0, .7D0, 13*0.D0, .2D0, .5D0, .7D0,
     & 1.3D0, 1.55D0, 1.9D0, 1.8D0, 1.55D0, 1.35D0, 1.15D0, .95D0, .7D0,
     & 17*0.D0, .2D0, .5D0, .85D0, 2.D0, 2.15D0, 2.05D0, 1.75D0, 1.D0,
     & 17*0.D0, .13D0, .33D0, .57D0, 1.33D0, 1.43D0, 1.36D0, 1.17D0,
     & .67D0, 17*0.D0, .07D0, .17D0, .28D0, .67D0, .72D0, .69D0, .58D0,
     & .33D0,17*0.D0,.4D0, .7D0, 1.D0, 1.6D0, 1.8D0, 2.3D0,1.9D0,1.7D0 /
**** pi- p data                                                        *
      DATA SPIKPW/ 0.D0, 25.D0, 13.D0, 11.D0, 10.5D0, 14.D0, 2*20.D0,
     & 16.D0, 14.D0, 19.D0, 28.D0, 17.5D0, 13.5D0, 12.D0, 10.5D0,
     & 2*10.D0, 9.5D0, 9.D0, 8.D0, 7.5D0, 7.D0, 6.5D0, 6.D0, 0.D0,
     & 48.D0, 19.D0, 15.D0, 11.5D0, 10.D0, 8.D0, 6.5D0, 5.5D0, 4.8D0,
     & 4.2D0, 7.5D0, 3.4D0, 2*2.5D0, 2.1D0, 1.4D0, 1.D0, .8D0, .6D0,
     & .46D0, .3D0, .2D0, .15D0, .13D0, 11*0.D0, .95D0, .65D0, .48D0,
     & .35D0, .2D0, .18D0, .17D0, .16D0, .15D0, .1D0, .09D0, .065D0,
     & .05D0, .04D0, 12*0.D0, .2D0, 2*.25D0, .2D0, .1D0, .08D0, .06D0,
     & .045D0, .03D0, .02D0, .01D0, .005D0, .003D0, 12*0.D0, .3D0,
     & .24D0, .18D0, .15D0, .13D0, .12D0, .11D0, .1D0, .09D0, .08D0,
     & .05D0, .04D0, .03D0, 0.D0, 0.16D0, .7D0, 1.3D0, 3.1D0, 4.5D0,
     & 2.D0, 23*0.D0, 4.0D0, 11.D0, 11.4D0, 10.3D0, 7.5D0, 6.8D0,
     & 4.75D0, 2.5D0, 1.5D0, .9D0, .55D0, .35D0, 13*0.D0, .1D0, .34D0,
     & .5D0, .8D0, 1.1D0, 2.25D0, 3.3D0, 2.3D0, 1.6D0, .95D0, .45D0,
     & .28D0, .15D0, 12*0.D0, .17D0, .64D0, 1.D0, 1.5D0, 2.1D0, 4.25D0,
     & 6.2D0, 4.4D0, 3.D0, 1.8D0, .9D0, .53D0, .28D0, 12*0.D0, .25D0,
     & .82D0, 1.3D0, 1.9D0, 2.8D0, 5.5D0, 8.D0, 5.7D0, 3.9D0, 2.35D0,
     & 1.15D0, .69D0, .37D0, 18*0.D0, .34D0, 1.5D0, 3.47D0, 5.87D0,
     & 6.23D0, 4.27D0, 2.6D0, 1.D0, .6D0, .3D0, .15D0, 6*0.D0/
*
      DATA SPIKPX/ 8*0.D0, .16D0, .75D0, 1.73D0, 2.93D0, 3.12D0,
     & 2.13D0, 1.3D0, .5D0, .3D0, .15D0, .08D0, 16*0.D0, .2D0, .6D0,
     & .92D0, 2.4D0, 4.9D0, 6.25D0, 5.25D0, 3.5D0, 2.15D0, 1.4D0, 1.D0,
     & .7D0, 13*0.D0, .13D0, .4D0, .62D0, 1.6D0, 3.27D0, 4.17D0, 3.5D0,
     & 2.33D0, 1.43D0, .93D0, .66D0, .47D0, 13*0.D0, .07D0, .2D0, .31D0,
     & .8D0, 1.63D0, 2.08D0, 1.75D0, 1.17D0, .72D0, .47D0, .34D0, .23D0,
     & 17*0.D0, .33D0, 1.D0, 1.8D0, 2.67D0, 5.33D0, 6.D0, 5.53D0, 5.D0,
     & 17*0.D0, .17D0, .5D0, .9D0, 1.83D0, 2.67D0, 3.0D0, 2.77D0, 2.5D0,
     & 6*0.D0, 1.D0, 3.3D0, 2.8D0, 2.5D0, 2.3D0, 1.8D0, 1.5D0, 1.1D0,
     & .8D0, .7D0, .55D0, .3D0, 19*0.D0, .1D0, .4D0, 1.D0, 1.4D0, 2.2D0,
     & 2.5D0, 2.2D0, 1.65D0, 1.35D0, 1.1D0, .8D0, .6D0, .4D0, 12*0.D0,
     & .15D0, .6D0, 1.5D0, 2.1D0, 3.3D0, 3.8D0, 3.3D0, 2.45D0, 2.05D0,
     & 1.65D0, 1.2D0, .9D0, .6D0, 12*0.D0, .10D0, .2D0, .5D0, .7D0,
     & 1.3D0, 1.55D0, 1.9D0, 1.8D0, 1.55D0, 1.35D0, 1.15D0, .95D0, .7D0,
     & 13*0.D0, .2D0, .5D0, .7D0, 1.3D0, 1.55D0, 1.9D0, 1.8D0, 1.55D0,
     & 1.35D0, 1.15D0, .95D0, .7D0, 17*0.D0, .2D0, .5D0, .85D0, 2.D0,
     & 2.15D0, 2.05D0, 1.75D0, 1.D0, 17*0.D0, .13D0, .33D0, .57D0,
     & 1.33D0, 1.43D0, 1.36D0, 1.17D0, .67D0, 17*0.D0, .07D0, .17D0,
     & .28D0, .67D0, .72D0, .69D0, .58D0, .33D0, 17*0.D0, .4D0, .7D0,
     & 1.D0, 1.6D0, 1.8D0, 2.3D0, 1.9D0, 1.7D0 /
**** pi- n data                                                        *
      DATA SPIKP4 / 0.D0, 300.D0, 40.D0, 20.D0, 13.D0, 8.5D0, 8.D0,
     & 9.5D0, 12.D0, 14.D0, 15.5D0, 20.D0, 17.D0, 13.D0, 10.D0, 9.D0,
     & 8.5D0, 8.D0, 7.8D0, 7.3D0, 6.7D0, 9*0.D0, .23D0, .35D0, .7D0,
     & .52D0, .4D0, .3D0, .2D0, .15D0, .13D0, .11D0, .09D0, .07D0, 0.D0,
     & .033D0, .8D0, 2*1.35D0, .5D0, 19*0.D0, 0.8D0, 2.2D0, 3.6D0,
     & 4.6D0, 4.7D0, 3.5D0, 2.4D0, 1.8D0, 1.4D0, .75D0, .47D0, .25D0,
     & .13D0, .08D0, 7*0.D0, 1.2D0, 3.3D0, 5.4D0, 6.9D0, 7.3D0, 5.3D0,
     & 3.6D0, 2.7D0, 2.2D0, 1.1D0, .73D0, .4D0, .22D0, .12D0, 11*0.D0,
     & 2.0D0, 4.4D0, 6.8D0, 9.9D0, 7.9D0, 6.0D0, 3.8D0, 2.5D0, 2.D0,
     & 1.4D0, 1.D0, .6D0, .35D0, 10*0.D0, .25D0, .55D0, .75D0, 1.25D0,
     & 1.9D0, 2.D0, 1.8D0, 1.5D0, 1.25D0, 1.D0, .8D0, 10*0.D0, .4D0,
     & .85D0, 1.1D0, 1.85D0, 2.8D0, 3.D0, 2.7D0, 2.2D0, 1.85D0, 1.5D0,
     & 1.2D0, 12*0.D0, .5D0, 1.2D0, 1.7D0, 3.4D0, 5.2D0, 6.4D0, 6.1D0,
     & 5.6D0, 5.2D0, 9*0.D0, 1.D0, 3.3D0, 5.2D0, 4.45D0, 3.6D0, 2.75D0,
     & 1.9D0, 1.65D0, 1.3D0, .95D0, .6D0, .45D0, 10*0.D0, .45D0, 1.4D0,
     & 1.5D0, 1.1D0, .85D0, .5D0, .3D0, .2D0, .15D0, 15*0.D0, .6D0,
     & .8D0, .95D0, .8D0, .7D0, .6D0, .5D0, .4D0, 13*0.D0, .85D0, 1.2D0,
     & 1.4D0, 1.2D0, 1.05D0, .9D0, .7D0, .55D0, 13*0.D0, 1.D0, 1.5D0,
     & 3.5D0, 4.15D0, 3.7D0, 2.7D0, 2.3D0, 1.75D0, 16*0.D0, .5D0, 2.0D0,
     & 3.3D0, 5.4D0, 7.D0 /
**** k+  p data                                                        *
      DATA SPIKP5/ 0.D0, 20.D0, 14.D0, 12.D0, 11.5D0, 10.D0, 8.D0,
     & 7.D0, 6.D0, 5.5D0, 5.3D0, 5.D0, 4.5D0, 4.4D0, 3.8D0, 3.D0, 2.8D0,
     & 0.D0, .5D0, 1.15D0, 2.D0, 1.3D0, .8D0, .45D0, 13*0.D0, 0.9D0,
     & 2.5D0, 3.D0, 2.5D0, 2.3D0, 2.D0, 1.7D0, 1.5D0, 1.2D0, .9D0, .6D0,
     & .45D0, .21D0, .2D0, 3*0.D0, .9D0, 2.5D0, 3.D0, 2.5D0, 2.3D0,
     & 2.D0, 1.7D0, 1.5D0, 1.2D0, .9D0, .6D0, .45D0, .21D0, .2D0,
     & 4*0.D0, 1.D0, 2.1D0, 2.6D0, 2.3D0, 2.1D0, 1.8D0, 1.7D0, 1.4D0,
     & 1.2D0, 1.05D0, .9D0, .66D0, .5D0, 7*0.D0, .3D0, 2*1.D0, .9D0,
     & .7D0, .4D0, .3D0, .2D0, 11*0.D0, .1D0, 1.D0, 2.2D0, 3.5D0, 4.2D0,
     & 4.55D0, 4.85D0, 4.9D0, 10*0.D0, .2D0, .7D0, 1.6D0, 2.5D0, 2.2D0,
     & 1.71D0, 1.6D0, 6*0.D0, 1.4D0, 3.8D0, 5.D0, 4.7D0, 4.4D0, 4.D0,
     & 3.5D0, 2.85D0, 2.35D0, 2.01D0, 1.8D0, 12*0.D0, .1D0, .8D0,2.05D0,
     & 3.31D0, 3.5D0, 12*0.D0, .034D0, .2D0, .75D0, 1.04D0, 1.24D0 /
 
C222222222222222222222222222222222222222222222222222222222222222222222 
C2222222222222222222222222222222222222222222222222222222222222222222222 
 
 
 
C**** K+ N DATA
      DATA S PIKP6/ 0.,6.,11.,13.,6.,5.,3.,2.2,1.5,1.2,1.,.7,.6,.5,.45,.
     +35,.3, 0.,6.,11.,13.,6.,5.,3.,2.2,1.5,1.2,1.,.7,.6,.5,.45,.35,.3,
     +0.,.5,1.3,2.8,2.3,1.6,.9,10*0., 3*0.,0.9,2.5,3.,2.5,2.3,2.,1.7,
     +1.5,1.2,.9,.6,.45,.21,.2, 3*0.,0.9,2.5,3.,2.5,2.3,2.,1.7,1.5,1.2,.
     +9,.6,.45,.21,.2, 4*0.,1.0,2.1,2.6,2.3,2.0,1.8,1.7,1.4,1.2,1.15,.9,
     +.66, .5, 4*0.,1.0,2.1,2.6,2.3,2.1,1.8,1.7,1.4,1.2,1.15,.9,.66, .5,
     +7*0.,.3,1.,1.,.9,.7,.4,.35,.2,.00,0., 7*0.,.3,1.,1.,.9,.7,.4,.35,.
     +2,.00,0., 9*0.,.1,1.,2.4,3.5,4.25,4.55,4.85,4.9, 9*0.,.1,1.,2.4,
     +3.5,4.25,4.55,4.85,4.9, 10*0.,.2,.7,1.6,2.5,2.2,1.71,1.6, 10*0.,.
     +2,.7,1.6,2.5,2.2,1.71,1.6, 6*0.,1.4,3.8,5.,4.7,4.4,4.,3.5,2.85,
     +2.35,2.01,1.8, 6*0.,1.4,3.8,5.,4.7,4.4,4.,3.5,2.85,2.35,2.01,1.8,
     +12*0.,.1,.8,2.05,3.31,3.5, 12*0.,.034,.20,.75,1.04,1.24, .0,2.5,
     +15.,21.5,15.3,3.,1.5,10*0./
 
C333333333333333333333333333333333333333333333333333333333333333333333

**** k-  p data                                                        *
      DATA SKMPEL/ 0.D0, 35.D0, 22.D0, 25.D0, 17.D0, 9.D0, 9.5D0, 8.D0,
     &     7.D0, 6.5D0, 6.1D0, 5.D0, 4.8D0, 4.6D0, 4.45D0, 4.3D0, 4.2D0,
     &    0.D0, 8.D0, 3.5D0, 8.D0, 3.D0, 1.9D0, 1.7D0, 1.D0, .9D0, .8D0,
     &    .75D0, .5D0, .42D0, .38D0, .34D0, .25D0, .2D0,
     &    0.D0, 3.D0, 3.2D0, 3.5D0, 1.5D0, 1.4D0, 1.1D0, .6D0, .5D0,
     &    .35D0, .28D0, .25D0, .18D0, .12D0, .1D0, .08D0, .04D0,
     &    0.D0, 8.5D0, 2.4D0, 1.7D0, 1.3D0, 1.3D0, 1.1D0, .5D0,
     &    .4D0, .4D0, .35D0, .3D0, .28D0, .2D0, .16D0, .13D0, .11D0,
     &    0.D0, 7.D0, 4.8D0, 1.4D0, 1.9D0, .9D0, .4D0, .2D0, .13D0,
     &    .1D0, .08D0, .06D0, .04D0, .02D0, .015D0, .01D0, .01D0,
     &    0.D0, 5.5D0, 1.D0, .8D0, .75D0, .32D0, .2D0, .1D0, .09D0,
     &    .08D0, .065D0, .05D0, .04D0, .022D0, .017D0, 2*.01D0/
      DATA SPIKP7 / 0.D0, .56D0, 1.46D0, 3.16D0, 2.01D0, 1.28D0, .74D0,
     & 14*0.D0, 1.13D0, 2.61D0, 2.91D0, 2.58D0, 2.35D0, 2.02D0,
     & 1.91D0, 1.57D0, 1.35D0, 1.29D0, 1.01D0, .74D0, .65D0, 4*0.D0,
     & 1.13D0, 2.61D0, 2.91D0, 2.58D0, 2.35D0, 2.02D0, 1.91D0, 1.57D0,
     & 1.35D0, 1.29D0, 1.01D0, .74D0, .65D0,  3*0.D0, 1.0D0, 3.03D0,
     & 3.36D0, 2.8D0, 2.58D0, 2.24D0, 1.91D0, 1.68D0, 1.35D0, 1.01D0,
     & .67D0, .5D0, .24D0, .23D0, 3*0.D0, 1.0D0, 3.03D0, 3.36D0, 2.8D0,
     & 2.58D0, 2.24D0, 1.91D0, 1.68D0, 1.35D0, 1.01D0, .67D0, .5D0,
     & .24D0, .23D0, 7*0.D0, .34D0, 1.12D0, 1.12D0, 1.01D0, .78D0,
     & .45D0, .39D0, .22D0, .07D0, 0.D0, 7*0.D0, .34D0, 1.12D0, 1.12D0,
     & 1.01D0, .78D0, .45D0, .39D0, .22D0, .07D0, 0.D0, 6*0.D0, 1.71D0,
     & 4.26D0, 5.6D0, 5.57D0, 4.93D0, 4.48D0, 3.92D0, 3.19D0, 2.63D0,
     & 2.25D0, 2.D0, 6*0.D0, 1.71D0, 4.26D0, 5.6D0, 5.57D0, 4.93D0,
     & 4.48D0, 3.92D0, 3.19D0, 2.63D0, 2.25D0, 2.D0, 10*0.D0, .22D0,
     & .8D0, .75D0, 1.D0, 1.3D0, 1.5D0, 1.3D0, 10*0.D0, .22D0, .8D0,
     & .75D0, 1.D0, 1.3D0, 1.5D0, 1.3D0, 13*0.D0, .1D0, .3D0, .7D0,1.D0,
     & 13*0.D0, .1D0, .3D0, .7D0, 1.D0, 9*0.D0, .11D0, 1.72D0, 2.69D0,
     & 3.92D0, 4.76D0, 5.10D0, 5.44D0, 5.3D0, 9*0.D0, .11D0, 1.72D0,
     & 2.69D0, 3.92D0, 4.76D0, 5.1D0, 5.44D0, 5.3D0, 5*0.D0,9.2D0,4.7D0,
     & 1.9D0, 10*0.D0, 2.5D0, 15.D0, 21.5D0, 15.3D0, 3.D0, 1.5D0,
     & 10*0.D0/
C333333333333333333333333333333333333333333333333333333333333333333333 
C3333333333333333333333333333333333333333333333333333333333333333333333
 
 
 
 
C**** K- N DATA
      DATA SKMNEL/
     *0.,4.,9.5,20.,13.,9.5,6.,4.4,3.,2.4,2.,1.4,1.2,1.,.9,.7,.6,
     *0.,4.5,6.,5.,2.5,2.,1.7,2.1,1.9,.9,.5,.3,.24,.2,.18,.1,.09,
     *0.,1.8,2.,1.1,.9,.5,.5,.4,.4,.2,.1,.06,.05,.04,.03,.02,.02,
     *0.,1.5,2.,.9,1.1,.4,.6,.7,.65,.3,.17,.1,.08,.07,.06,.04,.03/
      DATA S PIKP8/ 0.,.56,1.29,2.26,1.01,.64,.37,10*0., 4*0.,1.13,2.61,
     +2.91,2.58,2.35,2.02,1.91,1.57,1.35,1.29,1.01,.74, .65, 3*0.,1.00,
     +3.03,3.36,2.8,2.58,2.24,1.91,1.68,1.35,1.01,.67,.5,.24, .23, 3*0.,
     +1.00,3.03,3.36,2.8,2.58,2.24,1.91,1.68,1.35,1.01,.67,.5,.24, .23,
     +7*0.,.34,1.12,1.12,1.01,.78,.45,.39,.22,.07,0., 6*0.,1.71,4.26,
     +5.6,5.57,4.93,4.48,3.92,3.19,2.63,2.25,2., 10*0.,.22,.8,.75,1.,
     +1.3,1.5,1.3, 13*0.,.1,.3,.7,1., 13*0.,.1,.3,.7,1., 9*0.,.11,1.72,
     +2.69,3.92,4.76,5.10,5.44,5.3, 4*0.,0.00,9.2,4.7,1.9,9*0. /
 
 
 
 
 
 
 
 
C****  P P DATA
      DATA S PIKP9/ 0.,24.,25.,27.,23.,21.,20.,19.,17.,15.5,14.,13.5,
     +13., 0.,3.6,1.7, 10*0., .0,0.0,8.7,17.7,18.8,15.9,11.7,8.,6.,5.3,
     +4.5,3.9,3.5, .0,.0,2.8,5.8,6.2,5.1,3.8,2.7,2.1,1.8,1.5,1.3,1.1, 4
     +*0.,0.0,4.6,10.2,15.1,16.9,16.5,11.,5.5,3.5, 10*0.,4.3,7.6,9., 10
     +*0.,1.7,2.6,3., 6*0.,.3,.6,1.,1.6,1.3,.8,.6, 6*0.,.7,1.2,1.8,2.5,
     +1.8,1.3,1.2, 10*0.,.6,1.4,1.7, 10*0.,1.9,4.1,5.2/
 
 
C444444444444444444444444444444444444444444444444444444444444444444444
*****  p n data                                                        *
      DATA SPIKP0/ 0.D0, 24.D0, 25.D0, 27.D0, 23.D0, 21.D0, 20.D0,
     &              19.D0, 17.D0, 15.5D0, 14.D0, 13.5D0, 13.D0,
     &              0.D0, 1.8D0, .2D0,  12*0.D0,
     &              3.2D0, 6.05D0, 9.9D0, 5.1D0,
     &              3.8D0, 2.7D0, 1.9D0, 1.5D0, 1.4D0, 1.3D0, 1.1D0,
     &              2*.0D0, 3.2D0, 6.05D0, 9.9D0, 5.1D0,
     &              3.8D0, 2.7D0, 1.9D0, 1.5D0, 1.4D0, 1.3D0, 1.1D0,
     &              5*0.D0, 4.6D0, 10.2D0, 15.1D0,
     &              16.4D0, 15.2D0, 11.D0, 5.4D0, 3.5D0,
     &              5*0.D0, 4.6D0, 10.2D0, 15.1D0,
     &              16.4D0, 15.2D0, 11.D0, 5.4D0, 3.5D0,
     &              10*0.D0, .7D0, 5.1D0, 8.D0,
     &              10*0.D0, .7D0, 5.1D0, 8.D0,
     &              10*.0D0, .3D0, 2.8D0, 4.7D0,
     &              10*.0D0, .3D0, 2.8D0, 4.7D0,
     &              7*0.D0, 1.2D0, 2.5D0, 3.5D0, 6.D0, 5.3D0, 2.9D0,
     &              7*0.D0, 1.7D0, 3.6D0, 5.4D0, 9.D0, 7.6D0, 4.2D0,
     &              5*0.D0, 7.7D0, 6.1D0, 2.9D0, 5*0.D0/
*   nn - data                                                          *
*                                                                      *
      DATA SPKP V/ 0.D0, 24.D0, 25.D0, 27.D0, 23.D0, 21.D0, 20.D0,
     &              19.D0, 17.D0, 15.5D0, 14.D0, 13.5D0, 13.D0,
     &              0.D0, 3.6D0, 1.7D0, 12*0.D0,
     &              8.7D0, 17.7D0, 18.8D0, 15.9D0,
     &              11.7D0, 8.D0, 6.D0, 5.3D0, 4.5D0, 3.9D0, 3.5D0,
     &              .0D0, .0D0, 2.8D0, 5.8D0, 6.2D0, 5.1D0, 3.8D0,
     &              2.7D0, 2.1D0, 1.8D0, 1.5D0, 1.3D0, 1.1D0,
     &              5*0.D0, 4.6D0, 10.2D0, 15.1D0, 16.9D0, 16.5D0,
     &              11.D0, 5.5D0, 3.5D0,
     &              10*0.D0, 4.3D0, 7.6D0, 9.D0,
     &              10*0.D0, 1.7D0, 2.6D0, 3.D0,
     &              6*0.D0, .3D0, .6D0, 1.D0, 1.6D0, 1.3D0, .8D0, .6D0,
     &              6*0.D0, .7D0, 1.2D0, 1.8D0, 2.5D0, 1.8D0, 1.3D0,
     &              1.2D0, 10*0.D0, .6D0, 1.4D0, 1.7D0,
     &              10*0.D0, 1.9D0, 4.1D0, 5.2D0/
****************   ap - p - data                                       *
      DATA SAPPEL/ 0.D0,  176.D0, 160.D0, 105.D0, 75.D0, 68.D0, 65.D0,
     &  50.D0,  50.D0, 43.D0, 42.D0, 40.5D0, 35.D0, 30.D0, 28.D0,
     &  25.D0,  22.D0, 21.D0, 20.D0, 18.D0, 17.D0,  11*0.D0,
     &  .05D0,  .15D0, .18D0, .2D0, .2D0, .3D0, .4D0, .6D0, .7D0, .85D0,
     &  0.D0,  1.D0, .9D0, .46D0, .3D0, .23D0, .18D0, .16D0, .14D0,
     &  .1D0,  .08D0, .05D0, .02D0, .015D0, 4*.011D0, 3*.005D0,
     &  0.D0,  55.D0, 50.D0, 25.D0, 15.D0, 15.D0, 14.D0, 12.D0,
     &  10.D0,  7.D0, 6.D0, 4.D0, 3.3D0, 2.8D0, 2.4D0, 2.D0, 1.8D0,
     &  1.55D0,  1.3D0, .95D0, .75D0,
     &  0.D0,  3.3D0, 3.D0, 1.5D0, 1.D0, .7D0, .4D0, .35D0, .4D0,
     &  .25D0,  .18D0, .08D0, .04D0, .03D0, .023D0, .016D0, .014D0,
     & .01D0,  .008D0, .006D0, .005D0/
      DATA SPIKPE/0.D0, 215.D0, 193.D0, 170.D0, 148.D0, 113.D0, 97.D0,
     & 84.D0, 78.D0, 68.D0, 64.D0, 61.D0, 46.D0, 36.D0, 31.3D0, 28.5D0,
     & 25.7D0, 22.6D0, 21.4D0, 20.7D0, 19.9D0,
     & 9*0.D0, 2.D0, 2.5D0, .2D0, 19*0.D0, .3D0, 1.4D0, 2.2D0, 1.2D0,
     & 1.1D0, 1.D0, .8D0, .6D0, .5D0, .4D0, .3D0, 10*0.D0, .3D0, 1.4D0,
     & 2.2D0, 1.2D0, 1.1D0, 1.D0, .8D0, .6D0, .5D0, .4D0, .3D0, 10*0.D0,
     & .3D0, 1.4D0, 2.2D0, 1.2D0, 1.1D0, 1.D0, .8D0, .6D0, .5D0, .4D0,
     & .3D0, 10*0.D0, .3D0, 1.4D0, 2.2D0, 1.2D0, 1.1D0, 1.D0, .8D0,
     & .6D0, .5D0, .4D0, .3D0, 9*0.D0, .6D0, 2.5D0, 5.D0, 5.2D0, 5.1D0,
     & 5.4D0, 5.8D0, 2.8D0, 2.1D0, 1.8D0, 1.6D0, 1.2D0, 13*0.D0, 1.3D0,
     & 1.5D0, 2.D0, 2.5D0, 2.5D0, 2.3D0, 1.8D0, 1.4D0, 13*0.D0, 1.3D0,
     & 1.5D0, 2.D0, 2.5D0, 2.5D0, 2.3D0, 1.8D0, 1.4D0, 13*0.D0, 1.3D0,
     & 1.5D0, 2.D0, 2.5D0, 2.5D0, 2.3D0, 1.8D0, 1.4D0, 13*0.D0, 1.3D0,
     & 1.5D0, 2.D0, 2.5D0, 2.5D0, 2.3D0, 1.8D0, 1.4D0, 14*0.D0, .2D0,
     & .5D0, 1.1D0, 1.6D0, 1.4D0, 1.1D0, .9D0, 14*0.D0, .2D0, .5D0,
     & 1.1D0, 1.6D0, 1.4D0, 1.1D0, .9D0, 14*0.D0, .2D0, .5D0, 1.1D0,
     & 1.6D0, 1.4D0, 1.1D0, .9D0, 14*0.D0, .2D0, .5D0, 1.1D0, 1.6D0,
     & 1.4D0, 1.1D0, .9D0, 17*0.D0, .3D0, 1.6D0, 2.6D0, 3.6D0, 17*0.D0,
     & .3D0, 1.6D0, 2.6D0, 3.6D0, 17*0.D0, .3D0, 1.6D0, 2.6D0,
     & 3.6D0, 17*0.D0, .3D0, 1.6D0, 2.6D0, 3.6D0 /
****************   ap - n - data                                       *
      DATA SAPNEL/
     & 0.D0,  176.D0, 160.D0, 105.D0, 75.D0,  68.D0, 65.D0,
     & 50.D0, 50.D0,  43.D0,  42.D0,  40.5D0, 35.D0, 30.D0,  28.D0,
     & 25.D0, 22.D0,  21.D0,  20.D0,  18.D0,  17.D0, 11*0.D0,
     & .05D0, .15D0, .18D0,  .2D0,    .2D0,  .3D0,  .4D0,   .6D0,  .7D0,
     & .85D0,  0.D0,  1.D0,  .9D0,    .46D0, .3D0,  .23D0, .18D0, .16D0,
     & .14D0,  .1D0, .08D0, .05D0,    .02D0, .015D0, 4*.011D0, 3*.005D0,
     & 0.D0,  3.3D0,  3.D0, 1.5D0,     1.D0, .7D0,  .4D0,  .35D0, .4D0,
     & .25D0, .18D0, .08D0, .04D0,    .03D0, .023D0, .016D0, .014D0,
     & .01D0, .008D0, .006D0, .005D0 /
       DATA SPIKPZ/ 0.D0, 215.D0, 193.D0, 170.D0, 148.D0, 113.D0, 97.D0,
     &  84.D0, 78.D0, 68.D0, 64.D0, 61.D0, 46.D0, 36.D0, 31.3D0, 28.5D0,
     & 25.7D0, 22.6D0, 21.4D0, 20.7D0, 19.9D0, 9*0.D0, 2.4D0, .2D0,
     & 20*0.D0, 1.8D0, 2.8D0, 3.6D0, 2.3D0, 1.8D0, 1.5D0, 1.3D0, 1.D0,
     & .7D0, .5D0, .3D0, 10*0.D0, 1.8D0, 2.8D0, 3.6D0, 2.3D0, 1.8D0,
     & 1.5D0, 1.3D0, 1.D0, .7D0, .5D0, .3D0, 10*0.D0, 1.8D0, 2.8D0,
     & 3.6D0, 2.3D0, 1.8D0, 1.5D0, 1.3D0, 1.D0, .7D0, .5D0, .3D0,
     & 10*0.D0, 1.8D0, 2.8D0, 3.6D0, 2.3D0, 1.8D0, 1.5D0, 1.3D0, 1.D0,
     & .7D0, .5D0, .3D0, 13*0.D0, 5.2D0, 8.7D0, 11.4D0, 14.D0, 11.9D0,
     & 7.6D0, 6.D0, 5.D0, 13*0.D0, 5.2D0, 8.7D0, 11.4D0, 14.D0, 11.9D0,
     & 7.6D0, 6.D0, 5.D0, 18*0.D0, 1.D0, 4.9D0, 8.5D0, 18*0.D0, 1.D0,
     & 4.9D0, 8.5D0,  15*0.D0, 1.9D0, 2.3D0, 4.D0, 6.5D0, 5.2D0, 3.4D0,
     & 15*0.D0, 1.9D0, 2.3D0, 4.D0, 6.5D0, 5.2D0, 3.4D0, 15*0.D0, 1.9D0,
     & 2.3D0, 4.D0, 6.5D0, 5.2D0, 3.4D0 /
*                                                                      *
*                                                                      *
****************   an - p - data                                       *
*                                                                      *
      DATA SANPEL/
     & 0.D0,  176.D0, 160.D0, 105.D0, 75.D0, 68.D0, 65.D0, 50.D0,
     & 50.D0, 43.D0,  42.D0,  40.5D0, 35.D0, 30.D0, 28.D0,
     & 25.D0, 22.D0,  21.D0,  20.D0,  18.D0, 17.D0, 11*0.D0, .05D0,
     & .15D0, .18D0,   .2D0,   .2D0,   .3D0,  .4D0, .6D0,   .7D0, .85D0,
     & 0.D0,   1.D0,   .9D0,  .46D0,  .3D0,  .23D0, .18D0, .16D0, .14D0,
     & .1D0,  .08D0,  .05D0,  .02D0, .015D0, 4*.011D0, 3*.005D0,
     & 0.D0,  3.3D0,  3.D0, 1.5D0, 1.D0, .7D0, .4D0, .35D0, .4D0, .25D0,
     & .18D0, .08D0, .04D0, .03D0, .023D0, .016D0, .014D0,
     & .01D0, .008D0, .006D0, .005D0 /
      DATA SPIKPF/ 0.D0, 215.D0, 193.D0, 170.D0, 148.D0, 113.D0, 97.D0,
     & 84.D0, 78.D0, 68.D0, 64.D0, 61.D0, 46.D0, 36.D0, 31.3D0, 28.5D0,
     & 25.7D0, 22.6D0, 21.4D0, 20.7D0, 19.9D0, 9*0.D0, 2.4D0, .2D0,
     & 20*0.D0, 1.8D0, 2.8D0, 3.6D0, 2.3D0, 1.8D0, 1.5D0, 1.3D0, 1.D0,
     & .7D0, .5D0, .3D0, 10*0.D0, 1.8D0, 2.8D0, 3.6D0, 2.3D0, 1.8D0,
     & 1.5D0, 1.3D0, 1.D0, .7D0, .5D0, .3D0, 10*0.D0, 1.8D0, 2.8D0,
     & 3.6D0, 2.3D0, 1.8D0, 1.5D0, 1.3D0, 1.D0, .7D0, .5D0, .3D0,
     & 10*0.D0, 1.8D0, 2.8D0, 3.6D0, 2.3D0, 1.8D0, 1.5D0, 1.3D0, 1.D0,
     & .7D0, .5D0, .3D0, 13*0.D0, 5.2D0, 8.7D0, 11.4D0, 14.D0, 11.9D0,
     & 7.6D0, 6.D0, 5.D0, 13*0.D0, 5.2D0, 8.7D0, 11.4D0, 14.D0, 11.9D0,
     & 7.6D0, 6.D0, 5.D0, 18*0.D0, 1.D0, 4.9D0, 8.5D0, 18*0.D0, 1.D0,
     & 4.9D0, 8.5D0, 15*0.D0, 1.9D0, 2.3D0, 4.D0, 6.5D0, 5.2D0, 3.4D0,
     & 15*0.D0, 1.9D0, 2.3D0, 4.D0, 6.5D0, 5.2D0, 3.4D0, 15*0.D0, 1.9D0,
     & 2.3D0, 4.D0, 6.5D0, 5.2D0, 3.4D0 /
****  ko - n - data                                                    *
      DATA SPKP15/0.D0, 20.D0, 14.D0, 12.D0, 11.5D0, 10.D0, 8.D0, 7.D0,
     &      6.D0, 5.5D0, 5.3D0, 5.D0, 4.5D0, 4.4D0, 3.8D0, 3.D0, 2.8D0,
     &      0.D0, .5D0, 1.15D0, 2.D0, 1.3D0, .8D0, .45D0, 10*0.D0,
     &    3*0.D0, 0.9D0, 2.5D0, 3.D0, 2.5D0, 2.3D0, 2.D0, 1.7D0,
     &     1.5D0, 1.2D0, .9D0, .6D0, .45D0, .21D0, .2D0,
     &    3*0.D0, 0.9D0, 2.5D0, 3.D0, 2.5D0, 2.3D0, 2.D0, 1.7D0,
     &     1.5D0, 1.2D0, .9D0, .6D0, .45D0, .21D0, .2D0,
     &    4*0.D0, 1.D0, 2.1D0, 2.6D0, 2.3D0, 2.1D0, 1.8D0, 1.7D0,
     &     1.4D0, 1.2D0, 1.05D0, .9D0, .66D0,  .5D0,
     &    7*0.D0, .3D0, 1.D0, 1.D0, .9D0, .7D0, .4D0, .30D0, .2D0,
     &   11*0.D0, .1D0, 1.D0, 2.2D0, 3.5D0, 4.20D0, 4.55D0,
     &    4.85D0, 4.9D0,
     &   10*0.D0, .2D0, .7D0, 1.6D0, 2.5D0, 2.2D0, 1.71D0, 1.6D0,
     &    6*0.D0, 1.4D0, 3.8D0, 5.D0, 4.7D0, 4.4D0, 4.D0, 3.5D0,
     &    2.85D0, 2.35D0, 2.01D0, 1.8D0,
     &   12*0.D0, .1D0, .8D0, 2.05D0, 3.31D0, 3.5D0,
     &   12*0.D0, .034D0, .20D0, .75D0, 1.04D0, 1.24D0  /

C444444444444444444444444444444444444444444444444444444444444444444444 
C44444444444444444444444444444444444444444444444444444444444444444
C*** AKO - P - DATA
      DATA SPKP16/
     *0.,4.,9.5,20.,13.,9.5,6.,4.4,3.,2.4,2.,1.4,1.2,1.,.9,.7,.6,
     *0.,4.5,6.,5.,2.5,2.,1.7,2.1,1.9,.9,.5,.3,.24,.2,.18,.1,.09,
     *0.,1.8,2.,1.1,.9,.5,.5,.4,.4,.2,.1,.06,.05,.04,.03,.02,.02,
     *0.,1.5,2.,.9,1.1,.4,.6,.7,.65,.3,.17,.1,.08,.07,.06,.04,.03,
     *0.,.56,1.29,2.26,1.01,.64,.37,10*0.,
     *4*0.,1.13,2.61,2.91,2.58,2.35,2.02,1.91,1.57,1.35,1.29,1.01,.74,
     *.65,
     *3*0.,1.00,3.03,3.36,2.8,2.58,2.24,1.91,1.68,1.35,1.01,.67,.5,.24,
     *.23,
     *3*0.,1.00,3.03,3.36,2.8,2.58,2.24,1.91,1.68,1.35,1.01,.67,.5,.24,
     *.23,
     *7*0.,.34,1.12,1.12,1.01,.78,.45,.39,.22,.07,0.,
     *6*0.,1.71,4.26,5.6,5.57,4.93,4.48,3.92,3.19,2.63,2.25,2.,
     *10*0.,.22,.8,.75,1.,1.3,1.5,1.3,
     *13*0.,.1,.3,.7,1.,    13*0.,.1,.3,.7,1.,
     *9*0.,.11,1.72,2.69,3.92,4.76,5.10,5.44,5.3,
     *4*0.,0.00,9.2,4.7,1.9,9*0.
     */
      DATA NURE/9,12,5*0,10,14,3*0,1,3,5,7,6*0,2,6,16,5*0,
     *10,13,5*0,11,12,3*0,2,4,6,8,6*0,3,15,7,5*0/
      END
**sr 19-11-95: DNUPRE removed
*-- Author :
C
C*******************************************************************
C
      SUBROUTINE TSAMCS(KPROJ,EKIN,CST)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  HJM 24/10/88
C
C                 SAMPLING OF COS(THETA)
C                 FOR NUCLEON-PROTON ELASTIC SCATTERING
C                 ACCORDING TO HETKFA2/BERTINI PARAMETRIZATION
C
C-------------------------------------------------------------  ------
C
      DIMENSION DCLIN(195),DCHN(143),DCHNA(36),DCHNB(60)
      DIMENSION PDCI(60),PDCH(55)
C
C--------------------------------------------------------------------
C
      DATA (DCLIN(I),I=1,80) /
C***     DCLN ARRAY
     *     5.000D-01,  1.000D+00,  0.000D+00,  1.000D+00,  0.000D+00,
     *     4.993D-01,  9.881D-01,  5.963D-02,  9.851D-01,  5.945D-02,
     *     4.936D-01,  8.955D-01,  5.224D-01,  8.727D-01,  5.091D-01,
     *     4.889D-01,  8.228D-01,  8.859D-01,  7.871D-01,  8.518D-01,
     *     4.874D-01,  7.580D-01,  1.210D+00,  7.207D-01,  1.117D+00,
     *     4.912D-01,  6.969D-01,  1.516D+00,  6.728D-01,  1.309D+00,
     *     5.075D-01,  6.471D-01,  1.765D+00,  6.667D-01,  1.333D+00,
     *     5.383D-01,  6.054D-01,  1.973D+00,  7.059D-01,  1.176D+00,
     *     5.397D-01,  5.990D-01,  2.005D+00,  7.023D-01,  1.191D+00,
     *     5.336D-01,  6.083D-01,  1.958D+00,  6.959D-01,  1.216D+00,
     *     5.317D-01,  6.075D-01,  1.962D+00,  6.897D-01,  1.241D+00,
     *     5.300D-01,  6.016D-01,  1.992D+00,  6.786D-01,  1.286D+00,
     *     5.281D-01,  6.063D-01,  1.969D+00,  6.786D-01,  1.286D+00,
     *     5.280D-01,  5.960D-01,  2.020D+00,  6.667D-01,  1.333D+00,
     *     5.273D-01,  5.920D-01,  2.040D+00,  6.604D-01,  1.358D+00,
     *     5.273D-01,  5.862D-01,  2.069D+00,  6.538D-01,  1.385D+00/
      DATA (DCLIN(I),I=81,160) /
C***     DCIN ARRAY
     *     5.223D-01,  5.980D-01,  2.814D+00,  6.538D-01,  1.385D+00,
     *     5.202D-01,  5.969D-01,  2.822D+00,  6.471D-01,  1.412D+00,
     *     5.183D-01,  5.881D-01,  2.883D+00,  6.327D-01,  1.469D+00,
     *     5.159D-01,  5.866D-01,  2.894D+00,  6.250D-01,  1.500D+00,
     *     5.133D-01,  5.850D-01,  2.905D+00,  6.170D-01,  1.532D+00,
     *     5.106D-01,  5.833D-01,  2.917D+00,  6.087D-01,  1.565D+00,
     *     5.084D-01,  5.801D-01,  2.939D+00,  6.000D-01,  1.600D+00,
     *     5.063D-01,  5.763D-01,  2.966D+00,  5.909D-01,  1.636D+00,
     *     5.036D-01,  5.730D-01,  2.989D+00,  5.814D-01,  1.674D+00,
     *     5.014D-01,  5.683D-01,  3.022D+00,  5.714D-01,  1.714D+00,
     *     4.986D-01,  5.641D-01,  3.051D+00,  5.610D-01,  1.756D+00,
     *     4.964D-01,  5.580D-01,  3.094D+00,  5.500D-01,  1.800D+00,
     *     4.936D-01,  5.573D-01,  3.099D+00,  5.431D-01,  1.827D+00,
     *     4.909D-01,  5.509D-01,  3.144D+00,  5.313D-01,  1.875D+00,
     *     4.885D-01,  5.512D-01,  3.142D+00,  5.263D-01,  1.895D+00,
     *     4.857D-01,  5.437D-01,  3.194D+00,  5.135D-01,  1.946D+00/
      DATA (DCLIN(I),I=161,195) /
     *     4.830D-01,  5.353D-01,  3.253D+00,  5.000D-01,  2.000D+00,
     *     4.801D-01,  5.323D-01,  3.274D+00,  4.915D-01,  2.034D+00,
     *     4.770D-01,  5.228D-01,  3.341D+00,  4.767D-01,  2.093D+00,
     *     4.738D-01,  5.156D-01,  3.391D+00,  4.643D-01,  2.143D+00,
     *     4.701D-01,  5.010D-01,  3.493D+00,  4.444D-01,  2.222D+00,
     *     4.672D-01,  4.990D-01,  3.507D+00,  4.375D-01,  2.250D+00,
     *     4.634D-01,  4.856D-01,  3.601D+00,  4.194D-01,  2.323D+00/
C
      DATA PDCI /
     *     4.400D+02,  1.896D-01,  1.931D-01,  1.982D-01,  1.015D-01,
     *     1.029D-01,  4.180D-02,  4.228D-02,  4.282D-02,  4.350D-02,
     *     2.204D-02,  2.236D-02,  5.900D+02,  1.433D-01,  1.555D-01,
     *     1.774D-01,  1.000D-01,  1.128D-01,  5.132D-02,  5.600D-02,
     *     6.158D-02,  6.796D-02,  3.660D-02,  3.820D-02,  6.500D+02,
     *     1.192D-01,  1.334D-01,  1.620D-01,  9.527D-02,  1.141D-01,
     *     5.283D-02,  5.952D-02,  6.765D-02,  7.878D-02,  4.796D-02,
     *     6.957D-02,  8.000D+02,  4.872D-02,  6.694D-02,  1.152D-01,
     *     9.348D-02,  1.368D-01,  6.912D-02,  7.953D-02,  9.577D-02,
     *     1.222D-01,  7.755D-02,  9.525D-02,  1.000D+03,  3.997D-02,
     *     5.456D-02,  9.804D-02,  8.084D-02,  1.208D-01,  6.520D-02,
     *     8.233D-02,  1.084D-01,  1.474D-01,  9.328D-02,  1.093D-01/
C
      DATA PDCH /
     *     1.000D+03,  9.453D-02,  9.804D-02,  8.084D-02,  1.208D-01,
     *     6.520D-02,  8.233D-02,  1.084D-01,  1.474D-01,  9.328D-02,
     *     1.093D-01,  1.400D+03,  1.072D-01,  7.450D-02,  6.645D-02,
     *     1.136D-01,  6.750D-02,  8.580D-02,  1.110D-01,  1.530D-01,
     *     1.010D-01,  1.350D-01,  2.170D+03,  4.004D-02,  3.013D-02,
     *     2.664D-02,  5.511D-02,  4.240D-02,  7.660D-02,  1.364D-01,
     *     2.300D-01,  1.670D-01,  2.010D-01,  2.900D+03,  1.870D-02,
     *     1.804D-02,  1.320D-02,  2.970D-02,  2.860D-02,  5.160D-02,
     *     1.020D-01,  2.400D-01,  2.250D-01,  3.370D-01,  4.400D+03,
     *     1.196D-03,  8.784D-03,  1.517D-02,  2.874D-02,  2.488D-02,
     *     4.464D-02,  8.330D-02,  2.008D-01,  2.360D-01,  3.567D-01/
C
      DATA (DCHN(I),I=1,90) /
     *     4.770D-01,  4.750D-01,  4.715D-01,  4.685D-01,  4.650D-01,
     *     4.610D-01,  4.570D-01,  4.550D-01,  4.500D-01,  4.450D-01,
     *     4.405D-01,  4.350D-01,  4.300D-01,  4.250D-01,  4.200D-01,
     *     4.130D-01,  4.060D-01,  4.000D-01,  3.915D-01,  3.840D-01,
     *     3.760D-01,  3.675D-01,  3.580D-01,  3.500D-01,  3.400D-01,
     *     3.300D-01,  3.200D-01,  3.100D-01,  3.000D-01,  2.900D-01,
     *     2.800D-01,  2.700D-01,  2.600D-01,  2.500D-01,  2.400D-01,
     *     2.315D-01,  2.240D-01,  2.150D-01,  2.060D-01,  2.000D-01,
     *     1.915D-01,  1.850D-01,  1.780D-01,  1.720D-01,  1.660D-01,
     *     1.600D-01,  1.550D-01,  1.500D-01,  1.450D-01,  1.400D-01,
     *     1.360D-01,  1.320D-01,  1.280D-01,  1.250D-01,  1.210D-01,
     *     1.180D-01,  1.150D-01,  1.120D-01,  1.100D-01,  1.070D-01,
     *     1.050D-01,  1.030D-01,  1.010D-01,  9.900D-02,  9.700D-02,
     *     9.550D-02,  9.480D-02,  9.400D-02,  9.200D-02,  9.150D-02,
     *     9.100D-02,  9.000D-02,  8.990D-02,  8.900D-02,  8.850D-02,
     *     8.750D-02,  8.700D-02,  8.650D-02,  8.550D-02,  8.500D-02,
     *     8.499D-02,  8.450D-02,  8.350D-02,  8.300D-02,  8.250D-02,
     *     8.150D-02,  8.100D-02,  8.030D-02,  8.000D-02,  7.990D-02/
      DATA (DCHN(I),I=91,143) /
     *     7.980D-02,  7.950D-02,  7.900D-02,  7.860D-02,  7.800D-02,
     *     7.750D-02,  7.650D-02,  7.620D-02,  7.600D-02,  7.550D-02,
     *     7.530D-02,  7.500D-02,  7.499D-02,  7.498D-02,  7.480D-02,
     *     7.450D-02,  7.400D-02,  7.350D-02,  7.300D-02,  7.250D-02,
     *     7.230D-02,  7.200D-02,  7.100D-02,  7.050D-02,  7.020D-02,
     *     7.000D-02,  6.999D-02,  6.995D-02,  6.993D-02,  6.991D-02,
     *     6.990D-02,  6.870D-02,  6.850D-02,  6.800D-02,  6.780D-02,
     *     6.750D-02,  6.700D-02,  6.650D-02,  6.630D-02,  6.600D-02,
     *     6.550D-02,  6.525D-02,  6.510D-02,  6.500D-02,  6.499D-02,
     *     6.498D-02,  6.496D-02,  6.494D-02,  6.493D-02,  6.490D-02,
     *     6.488D-02,  6.485D-02,  6.480D-02/
C
      DATA DCHNA /
     *     6.300D+02,  7.810D-02,  1.421D-01,  1.979D-01,  2.479D-01,
     *     3.360D-01,  5.400D-01,  7.236D-01,  1.000D+00,  1.540D+03,
     *     2.225D-01,  3.950D-01,  5.279D-01,  6.298D-01,  7.718D-01,
     *     9.405D-01,  9.835D-01,  1.000D+00,  2.560D+03,  2.625D-01,
     *     4.550D-01,  5.963D-01,  7.020D-01,  8.380D-01,  9.603D-01,
     *     9.903D-01,  1.000D+00,  3.520D+03,  4.250D-01,  6.875D-01,
     *     8.363D-01,  9.163D-01,  9.828D-01,  1.000D+00,  1.000D+00,
     *     1.000D+00/
C
      DATA DCHNB /
     *     6.300D+02,  3.800D-02,  7.164D-02,  1.275D-01,  2.171D-01,
     *     3.227D-01,  4.091D-01,  5.051D-01,  6.061D-01,  7.074D-01,
     *     8.434D-01,  1.000D+00,  2.040D+03,  1.200D-01,  2.115D-01,
     *     3.395D-01,  5.295D-01,  7.251D-01,  8.511D-01,  9.487D-01,
     *     9.987D-01,  1.000D+00,  1.000D+00,  1.000D+00,  2.200D+03,
     *     1.344D-01,  2.324D-01,  3.754D-01,  5.674D-01,  7.624D-01,
     *     8.896D-01,  9.808D-01,  1.000D+00,  1.000D+00,  1.000D+00,
     *     1.000D+00,  2.850D+03,  2.330D-01,  4.130D-01,  6.610D-01,
     *     9.010D-01,  9.970D-01,  1.000D+00,  1.000D+00,  1.000D+00,
     *     1.000D+00,  1.000D+00,  1.000D+00,  3.500D+03,  3.300D-01,
     *     5.450D-01,  7.950D-01,  1.000D+00,  1.000D+00,  1.000D+00,
     *     1.000D+00,  1.000D+00,  1.000D+00,  1.000D+00,  1.000D+00/
C
C---------------------------------------------------------------
      CST=1D0
C
C*                      IS THE KINETIC ENERGY GREATER THAN LIMIT ?
C
      IF (EKIN.GT.3.5D0) RETURN
C
      IF(KPROJ.EQ.8) GOTO 101
      IF(KPROJ.EQ.1) GOTO 102
C*                                             INVALID REACTION
      WRITE(6,'(A,I5/A)')
     &        ' INVALID PARTICLE TYPE IN DNUPRE - KPROJ=',KPROJ,
     &        ' COS(THETA) = 1D0 RETURNED'
      RETURN
C-------------------------------- NP ELASTIC SCATTERING----------
101   CONTINUE
      IF (EKIN.GT.0.740D0)GOTO 1000
      IF (EKIN.LT.0.300D0)THEN
C                                 EKIN .LT. 300 MEV
         IDAT=1
      ELSE
C                                 300 MEV < EKIN < 740 MEV
         IDAT=6
      END IF
C
      ENER=EKIN
      IE=ABS(ENER/0.020D0)
      UNIV=(ENER-FLOAT(IE)*0.020D0)/0.020D0
C                                            FORWARD/BACKWARD DECISION
      K=IDAT+5*IE
      BWFW=(DCLIN(K+5)-DCLIN(K))*UNIV + DCLIN(K)
      IF (RNDM(V).LT.BWFW)THEN
         VALUE2=-1D0
         K=K+1
      ELSE
         VALUE2=1D0
         K=K+3
      END IF
C
      COEF=(DCLIN(K+5)-DCLIN(K))*UNIV + DCLIN(K)
      RND=RNDM(V)
C
      IF(RND.LT.COEF)THEN
         CST=RNDM(V)
         CST=CST*VALUE2
      ELSE
         R1=RNDM(V)
         R2=RNDM(V)
         R3=RNDM(V)
         R4=RNDM(V)
C
         IF(VALUE2.GT.0.0)THEN
            CST=MAX(R1,R2,R3,R4)
            GOTO 1500
         ELSE
            R5=RNDM(V)
C
            IF (IDAT.EQ.1)THEN
               CST=-MAX(R1,R2,R3,R4,R5)
            ELSE
               R6=RNDM(V)
               R7=RNDM(V)
               CST=-MAX(R1,R2,R3,R4,R5,R6,R7)
            END IF
C
         END IF
C
      END IF
C
      GOTO 1500
C
C********                                EKIN  .GT.  0.74 GEV
C
1000  ENER=EKIN - 0.66
C     IE=ABS(ENER/0.02)
      IE=ENER/0.02
      EMEV=EKIN*1D3
C
      UNIV=(ENER-FLOAT(IE)*0.020D0)/0.020D0
      K=IE
      BWFW=(DCHN(K+1)-DCHN(K))*UNIV + DCHN(K)
      RND=RNDM(V)
C                                        FORWARD NEUTRON
      IF (RND.GE.BWFW)THEN
         DO 1200 K=10,36,9
           IF (DCHNA(K).GT.EMEV) THEN
              UNIVE=(EMEV-DCHNA(K-9))/(DCHNA(K)-DCHNA(K-9))
              UNIV=RNDM(V)
              DO 1100 I=1,8
                 II=K+I
                 P=(DCHNA(II)-DCHNA(II-9))*UNIVE + DCHNA(II-9)
C
                 IF (P.GT.UNIV)THEN
                    UNIV=RNDM(V)
                    FLTI=FLOAT(I)-UNIV
                    GOTO(290,290,290,290,330,340,350,360) I
                 END IF
 1100         CONTINUE
           END IF
 1200    CONTINUE
C
      ELSE
C                                        BACKWARD NEUTRON
         DO 1400 K=13,60,12
            IF (DCHNB(K).GT.EMEV) THEN
               UNIVE=(EMEV-DCHNB(K-12))/(DCHNB(K)-DCHNB(K-12))
               UNIV=RNDM(V)
               DO 1300 I=1,11
                 II=K+I
                 P=(DCHNB(II)-DCHNB(II-12))*UNIVE + DCHNB(II-12)
C
                 IF (P.GT.UNIV)THEN
                   UNIV=RNDM(V)
                   FLTI=FLOAT(I)-UNIV
                   GOTO(120,120,140,150,160,160,180,190,200,210,220) I
                 END IF
 1300          CONTINUE
            END IF
 1400    CONTINUE
      END IF
C
120   CST=1.0D-2*FLTI-1.0D0
      GOTO 1500
140   CST=2.0D-2*UNIV-0.98D0
      GOTO 1500
150   CST=4.0D-2*UNIV-0.96D0
      GOTO 1500
160   CST=6.0D-2*FLTI-1.16D0
      GOTO 1500
180   CST=8.0D-2*UNIV-0.80D0
      GOTO 1500
190   CST=1.0D-1*UNIV-0.72D0
      GOTO 1500
200   CST=1.2D-1*UNIV-0.62D0
      GOTO 1500
210   CST=2.0D-1*UNIV-0.50D0
      GOTO 1500
220   CST=3.0D-1*(UNIV-1.0D0)
      GOTO 1500
C
290   CST=1.0D0-2.5d-2*FLTI
      GOTO 1500
330   CST=0.85D0+0.5D-1*UNIV
      GOTO 1500
340   CST=0.70D0+1.5D-1*UNIV
      GOTO 1500
350   CST=0.50D0+2.0D-1*UNIV
      GOTO 1500
360   CST=0.50D0*UNIV
C
1500  RETURN
C
C-----------------------------------  PP ELASTIC SCATTERING -------
C
 102  CONTINUE
      EMEV=EKIN*1D3
C
      IF (EKIN.LE.0.500D0) THEN
         RND=RNDM(V)
         CST=2.0D0*RND-1.0D0
         RETURN
C
      ELSEIF (EKIN.LT.1.0D0) THEN
         DO 2200 K=13,60,12
            IF (PDCI(K).GT.EMEV) THEN
               UNIVE=(EMEV-PDCI(K-12))/(PDCI(K)-PDCI(K-12))
               UNIV=RNDM(V)
               SUM=0
               DO 2100 I=1,11
                 II=K+I
                 SUM=SUM + (PDCI(II)-PDCI(II-12))*UNIVE + PDCI(II-12)
C
                 IF (UNIV.LT.SUM)THEN
                   UNIV=RNDM(V)
                   FLTI=FLOAT(I)-UNIV
                   GOTO(55,55,55,60,60,65,65,65,65,70,70) I
                 END IF
 2100          CONTINUE
            END IF
 2200    CONTINUE
      ELSE
         DO 2400 K=12,55,11
            IF (PDCH(K).GT.EMEV) THEN
              UNIVE=(EMEV-PDCH(K-11))/(PDCH(K)-PDCH(K-11))
              UNIV=RNDM(V)
              SUM=0.0
              DO 2300 I=1,10
                II=K+I
                SUM=SUM + (PDCH(II)-PDCH(II-11))*UNIVE + PDCH(II-11)
C
                IF (UNIV.LT.SUM)THEN
                  UNIV=RNDM(V)
                  FLTI=UNIV+FLOAT(I)
                  GOTO(50,55,60,60,65,65,65,65,70,70) I
                END IF
 2300         CONTINUE
            END IF
 2400    CONTINUE
      END IF
C
50    CST=0.4D0*UNIV
      GOTO 2500
55    CST=0.2D0*FLTI
      GOTO 2500
60    CST=0.3D0+0.1D0*FLTI
      GOTO 2500
65    CST=0.6D0+0.04D0*FLTI
      GOTO 2500
70    CST=0.78D0+0.02D0*FLTI
C
2500  CONTINUE
      IF (RNDM(V).GT.0.5D0) CST=-CST
C
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
c/* *************************************************************
c/*    27/10/89
c/*    TEST OF CROSS SECTION routines
C      actually applied in dtunuc90:
c/*            SHPTOT (called from SHMAKI)
c/*            SIHNIN, SIHNEL (called from FOZOKL)
c/* *************************************************************
C
      SUBROUTINE SIGTES
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
*KEEP,PANAME.
C------------------
C
C     /PANAME/ CONTAINS PARTICLE NAMES
C        BTYPE  = LITERAL NAME OF THE PARTICLE
C
      CHARACTER*8 BTYPE
      COMMON /PANAME/ BTYPE(30)
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
      PARAMETER (NPA=45,IPMAX=16)
      DIMENSION KPROJ(IPMAX)
      DIMENSION P(NPA),EK(NPA),PL(180)
      DIMENSION SIGMA(NPA,4)
      DATA P /0.13, 0.19, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
     +1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 10.0, 20., 50.,
     +100., 200., 400., 1000., 5000., 10000.,
     +20000.,50000.,100000.,200000.,500000.,1000000.,2000000.,
     +5000000.,10000000.,20000000.,50000000.,100000000.,
     +200000000.,500000000.,1000000000./
      DATA KPROJ / 1, 8, 13, 14, 15, 16, 24, 25, 2, 9, 17, 18, 20, 21,
     +22, 23/
C********************************************************************
      WRITE(6,'(1H1)')
      WRITE(6,'(2A)') ' TEST OF CROSS SECTIONS ROUTINES APPLIED:',
     +' SHPTOT, SIHNIN, SIHNEL(SIGEL)'
      WRITE(6,'(A,5X,A)') ' PRINTED FOR EACH MOMENTUM/KINETIC ENERGY:',
     +' SIHNEL , SIHNIN / (SIHNEL+SIHNIN), SHPTOT'
 
C                                   LOOP OVER PROJECTILES
C     AA=1.0
      DO 50 IPR=1,IPMAX
        KPRO=KPROJ(IPR)
        DO 20 IE=1,NPA
          EK(IE)=SQRT(P(IE)**2+AAM(KPRO)**2) - AAM(KPRO)
C         CALL SIGEL(KPRO,AA,P(IE),SIEL1,ZLEL)
          CALL SIHNEL(KPRO,1,P(IE),SIEL)
C         CALL NIZL(KPRO,AA,P(IE),SIIN1,ZLIN)
          CALL SIHNIN(KPRO,1,P(IE),SIIN)
          SIGMA(IE,2)=SIIN
          SIGMA(IE,1)=SIEL
          SIGMA(IE,4)=SIIN + SIEL
          SIGMA(IE,3)=DSHPTO(KPRO,P(IE))
          DO 10 IM=0,3
            IIE=IM*NPA +IE
            PL(IIE)=LOG10(P(IE))
   10     CONTINUE
   20   CONTINUE
C
        WRITE(6,'(1H1,A,I3,2A)') ' PROJECTILE IPRO=', KPRO, '   - ',
     +  BTYPE(KPRO)
        WRITE(6,'(//A/2A/)')
     +  ' TABLE OF CROSS SECTION VALUES VERSUS PLAB / EKIN:',
     +  ' PLAB, EKIN  /  SIG(EL), SIG(INEL)  /',
     +  '  SIG(TOT)-USED,SIG(EL)+SIG(INEL) - TEST PRINT'
C
        DO 30 IE=1,NPA
          WRITE(6,1000) P(IE),EK(IE), (SIGMA(IE,IS),IS=1,4)
   30   CONTINUE
C
        DO 40 IM=1,3
          DO 40 IE=1,NPA
            SIGMA(IE,IM)=LOG10(SIGMA(IE,IM))
   40   CONTINUE
        WRITE(6,'(///A//)')
     +  ' DOUBLE-LOG PLOT OF SIGMA(EL/INEL/TOT/SUM EL+INEL) VS PLAB'
        NPOI=4*NPA
        CALL PLOT(PL,SIGMA,NPOI,4,NPA,-1.0D0,0.2D0,0.0D0,0.050D0)
C
   50 CONTINUE
C*****************
 1000 FORMAT(2(1PE10.2),5X,4(1PE12.3))
C
      END
