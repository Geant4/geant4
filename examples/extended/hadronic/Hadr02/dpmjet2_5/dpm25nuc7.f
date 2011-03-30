C-------------------------------------------------------------
C                    was dpmlund.f
C-------------------------------------------------------------
C     SUBROUTINE TESLUN
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     SAVE
C     COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
C     COMMON/POPCOR/PDB,AJSDEF
C     COMMON/HARLUN/IHARLU,QLUN
C     COMMON /HARLUN/ QLUN,IHARLU
C--------------------------------------------------
C                       IHARLU=0  soft jets fragmenting
C                       IHARLU=1  hard jets fragmenting
C                                 with final state evolution
C                                 in JETSET
C                       QLUN = Mass of hard partons at end of
C                              chain
C-----------------------------------------------------------
C     AJSDEF=0.
C     PDB=0.5
C     IPRI=2
c     CALL LUNDIN
C     IHARLU=0
C     QLUN=0.
C     IF(IPRI.GE.2) WRITE(6,111)QLUN
c 111 FORMAT(' QLUN= ',F10.2)
C     CALL BAMLUN(IHAD,1,7,0,0,25.,3,IREJ)
C     IF(IPRI.GE.2) WRITE(6,111)QLUN
C     CALL BAMLUN(IHAD,2,8,0,0,25.,3,IREJ)
C     IF(IPRI.GE.2) WRITE(6,111)QLUN
C     CALL BAMLUN(IHAD,3,9,0,0,25.,3,IREJ)
C     IHARLU=1
C     QLUN=8.
C     IF(IPRI.GE.2) WRITE(6,111)QLUN
C     CALL BAMLUN(IHAD,1,7,0,0,25.,3,IREJ)
C     IF(IPRI.GE.2) WRITE(6,111)QLUN
C     CALL BAMLUN(IHAD,2,8,0,0,25.,3,IREJ)
C     IF(IPRI.GE.2) WRITE(6,111)QLUN
C     CALL BAMLUN(IHAD,3,9,0,0,25.,3,IREJ)
c     QLUN=6.
C     IF(IPRI.GE.2) WRITE(6,111)QLUN
C     CALL BAMLUN(IHAD,1,7,0,0,25.,3,IREJ)
C     IF(IPRI.GE.2) WRITE(6,111)QLUN
C     CALL BAMLUN(IHAD,2,8,0,0,25.,3,IREJ)
C     IF(IPRI.GE.2) WRITE(6,111)QLUN
C 
C     CALL BAMLUN(IHAD,2,8,0,0,25.,3,IREJ)
C     IF(IPRI.GE.2) WRITE(6,111)QLUN
C     CALL BAMLUN(IHAD,3,9,0,0,25.,3,IREJ)
C     QLUN=2.
C     IF(IPRI.GE.2) WRITE(6,111)QLUN
C     CALL BAMLUN(IHAD,1,7,0,0,25.,3,IREJ)
C     IF(IPRI.GE.2) WRITE(6,111)QLUN
C     CALL BAMLUN(IHAD,2,8,0,0,25.,3,IREJ)
C     IF(IPRI.GE.2) WRITE(6,111)QLUN
C     CALL BAMLUN(IHAD,3,9,0,0,25.,3,IREJ)
C     IPRI=0
C     IHARLU=0
C     QLUN=0.
C     RETURN
C     END
      SUBROUTINE LUNDIN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      INTEGER PYCOMP
C
C     INITIALIZATION FOR JETSET-7.3 CALL IN DTUNUC 1.04 (J.R. 6/93)
C
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(4000,2),BRAT(4000),KFDP(4000,5)
      COMMON/POPCOR/PDB,AJSDEF
      COMMON/PROMU/IPROMU
      COMMON/IFRAGM/IFRAG
C     COMMON/HARLUN/IHARLU,QLUN
      COMMON /HARLUN/ QLUN,IHARLU
      DATA IPROMM/0/
C--------------------------------------------------
C                       IHARLU=0  soft jets fragmenting
C                       IHARLU=1  hard jets fragmenting
C                                 with final state evolution
C                                 in JEYTSET
C                       QLUN = Mass of hard partons at end of
C                              chain
C-----------------------------------------------------------
C                                defaults for parton showering
C                                QCD type branchings
      MSTJ(41)=1
C                                coherent branching, angular
C                                ordering
      MSTJ(42)=2
      MSTJ(43)=4
      MSTJ(44)=2
C                                 only u,d,s,c quarks
      MSTJ(45)=4
      MSTJ(46)=0
C                                  no lowest order corrections
      MSTJ(47)=0
      MSTJ(48)=0
      MSTJ(49)=0
C                                  Lambda in running alpha(s)
      PARJ(81)=0.40D0
C                                  M(min) cut-off
      PARJ(82)=1.D0
      Parj(83)=1.D0
C------------------------------------------------------------------
      IHARLU=0
      QLUN=0.D0
C-----------------------------------------------------------------
      IF(IPROMU.NE.0)IPROMM=IPROMU
C     IPROMM=1
C----------------------------------------------------------------
C        if AJSDEF=1 set default values for all jetset parameters
      IF(AJSDEF.EQ.1.D0) GO TO 100
C----------------------------------------------------------------
C                           switch off popcorn fragmentation
C                            MSTJ(12) default 2 : popcorn allowed
      IF(PDB.EQ.0.D0)MSTJ(12)=1 
C                           TEST proton x distribution 19.11.97
      MSTJ(12)=3
      PARJ(19)=0.1
C----------------------------------------------------------------
C                           MODIFY POPCORM mechanism
C                           PARJ(5), default:0.5
      IF(PDB.GT.0.D0)PARJ(5)=PDB
C
C----------------------------------------------------------------
C
C                   DEFINE Lund parameters
C                       IFRAG=1  
C
C----------------------------------------------------------------
C
      IF(IFRAG.EQ.1)THEN   
C                 Probability for meson with spin 1 (d=0.5) 
      PARJ(11)=0.6D0
C----------------------------------------------------------------
C                           Lund b-parameter (default=0.9)
      PARJ(42)=0.5D0
C----------------------------------------------------------------
C                           Lund a-parameter (default=0.5)
      PARJ(41)=0.2D0 
C----------------------------------------------------------------
C                d=1 Lund Fragmentation 2: F-F fragmentation   
C     MSTJ(11)=1
C                           F.-F. c-parameter (default=0.77)
C     PARJ(51)=0.99
C     PARJ(52)=0.99
C     PARJ(53)=0.99
C----------------------------------------------------------------
C                           Lund sigma parameter in pt distribution
C                                             (default=0.35)
      PARJ(21)=0.42D0
C----------------------------------------------------------------
C        Some of these parmeters are changed in BAMLUN
C!!!!!!!!

C----------------------------------------------------------------
C                           Diquark supression (default 0.1)
      PARJ(1)=0.07D0
C     PARJ(1)=0.06D0
C----------------------------------------------------------------
C                           STRANGENESS supression (default 0.3)
C     PARJ(2)=0.27D0
C                                        6.1.95
      PARJ(2)=0.25D0
C----------------------------------------------------------------
C                          Extra supression of strange diquarks
C                           Default 0.4
C     PARJ(3)=0.4D0
C                                        6.1.95
      PARJ(3)=0.3D0
C     PARJ(3)=0.5D0
C----------------------------------------------------------------
C                          Extra supression of Strangeness in
C                          B-M-B Situation (defaults 0.5)
C     PARJ(6)=0.75D0
C     PARJ(7)=0.75D0
C                                        6.1.95
      PARJ(6)=0.50D0
      PARJ(7)=0.50D0
C     PARJ(6)=0.75D0
C     PARJ(7)=0.75D0
C----------------------------------------------------------------
C
C----------------------------------------------------------------
C
C                   DEFINE Lund parameters
C                       IFRAG=10 
C                low energy tests 
C
C----------------------------------------------------------------
      ELSEIF(IFRAG.EQ.10)THEN
C                 Probability for meson with spin 1 (d=0.5) 
      PARJ(11)=0.6D0
C----------------------------------------------------------------
C                           Lund b-parameter (default=0.9)
      PARJ(42)=0.5
C----------------------------------------------------------------
C                           Lund a-parameter (default=0.5)
C     PARJ(41)=0.51
C                                      14.3.95 like 1
      PARJ(41)=0.20D0
C----------------------------------------------------------------
C                d=1 Lund Fragmentation 2: F-F fragmentation   
C     MSTJ(11)=1
C                           F.-F. c-parameter (default=0.77)
C     PARJ(51)=0.99
C     PARJ(52)=0.99
C     PARJ(53)=0.99
C----------------------------------------------------------------
C                           Lund sigma parameter in pt distribution
C                                             (default=0.35)
      PARJ(21)=0.42D0
C----------------------------------------------------------------
C        Some of these parmeters are changed in BAMLUN
C!!!!!!!!

C----------------------------------------------------------------
C                           Diquark supression (default 0.1)
      PARJ(1)=0.07D0
C----------------------------------------------------------------
C                           STRANGENESS supression (default 0.3)
C     PARJ(2)=0.27D0
C                                        6.1.95
      PARJ(2)=0.25D0
C----------------------------------------------------------------
C                          Extra supression of strange diquarks
C                           Default 0.4
C     PARJ(3)=0.4D0
C                                        6.1.95
      PARJ(3)=0.3D0
C----------------------------------------------------------------
C                          Extra supression of Strangeness in
C                          B-M-B Situation (defaults 0.5)
C     PARJ(6)=0.75D0
C     PARJ(7)=0.75D0
C                                        6.1.95
      PARJ(6)=0.50D0
      PARJ(7)=0.50D0
C----------------------------------------------------------------
C
C----------------------------------------------------------------
C----------------------------------------------------------------
      ENDIF
C----------------------------------------------------------------
  100 CONTINUE
      WRITE (6,2355)PDB,AJSDEF,MSTJ(12),PARJ(42),PARJ(21)
 2355 FORMAT( ' LUNDIN initialization PDB,AJSDEF= ',2F10.3/
     + ' MSTJ(12) popcorn default=2  : ',I10/
     + ' PARJ(42 )Lund b,default=0.9 : ',F10.3/
     + ' PARJ(21) sigma in pt distr,default=0.35 : ',F10.3)
C----------------------------------------------------------------
C                           Prevent particles dacaying
C                 KOS
      KC=PYCOMP(310)
      MDCY(KC,1)=0
C                 PIO
      KC=PYCOMP(111)
      MDCY(KC,1)=0
C                 LAMBDA
      KC=PYCOMP(3122)
      MDCY(KC,1)=0
C                 ALAMBDA
      KC=PYCOMP(-3122)
      MDCY(KC,1)=0
C                 SIG+
      KC=PYCOMP(3222)
      MDCY(KC,1)=0
C                 ASIG+
      KC=PYCOMP(-3222)
      MDCY(KC,1)=0
C                 SIG-
      KC=PYCOMP(3112)
      MDCY(KC,1)=0
C                 ASIG-
      KC=PYCOMP(-3112)
      MDCY(KC,1)=0
C                 SIG0
C     KC=PYCOMP(3212)
C     MDCY(KC,1)=0
C                 ASIG0
C     KC=PYCOMP(-3212)
C     MDCY(KC,1)=0
C                 TET0
      KC=PYCOMP(3322)
      MDCY(KC,1)=0
C                 ATET0
      KC=PYCOMP(-3322)
      MDCY(KC,1)=0
C                 TET-
      KC=PYCOMP(3312)
      MDCY(KC,1)=0
C                 ATET-
      KC=PYCOMP(-3312)
      MDCY(KC,1)=0
C                 OMEGA-
      KC=PYCOMP(3334)
      MDCY(KC,1)=0
C                 AOMEGA-
      KC=PYCOMP(-3334)
      MDCY(KC,1)=0
C                 TAU
      KC=PYCOMP(15)
      MDCY(KC,1)=0
C
C               HAVE CHARMED MESONS DECAYING for IPROMM=1
C
      IF(IPROMM.EQ.1)GO TO 199
C                 D+
      KC=PYCOMP(411)
      MDCY(KC,1)=0
C                 D-
      KC=PYCOMP(-411)
      MDCY(KC,1)=0
C                 D0
      KC=PYCOMP(421)
      MDCY(KC,1)=0
C                 A-D0
      KC=PYCOMP(-421)
      MDCY(KC,1)=0
C                 DS+
      KC=PYCOMP(431)
      MDCY(KC,1)=0
C                 A-DS+
      KC=PYCOMP(-431)
      MDCY(KC,1)=0
C                ETAC 
      KC=PYCOMP(441)
      MDCY(KC,1)=0
  199 CONTINUE
C                LAMBDAC+ 
      KC=PYCOMP(4122)
      MDCY(KC,1)=0
C                A-LAMBDAC+ 
      KC=PYCOMP(-4122)
      MDCY(KC,1)=0
C                SIGMAC++ 
      KC=PYCOMP(4222)
      MDCY(KC,1)=0
C                SIGMAC+ 
      KC=PYCOMP(4212)
      MDCY(KC,1)=0
C                SIGMAC0 
      KC=PYCOMP(4112)
      MDCY(KC,1)=0
C                A-SIGMAC++ 
      KC=PYCOMP(-4222)
      MDCY(KC,1)=0
C                A-SIGMAC+ 
      KC=PYCOMP(-4212)
      MDCY(KC,1)=0
C                A-SIGMAC0 
      KC=PYCOMP(-4112)
      MDCY(KC,1)=0
C                KSIC+ 
      KC=PYCOMP(4232)
      MDCY(KC,1)=0
C                KSIC0 
      KC=PYCOMP(4132)
      MDCY(KC,1)=0
C                A-KSIC+ 
      KC=PYCOMP(-4232)
      MDCY(KC,1)=0
C                A-KSIC0 
      KC=PYCOMP(-4132)
      MDCY(KC,1)=0
      RETURN
      END
      SUBROUTINE BAMLUN(IHAD,KFA1,KFA2,KFA3,KFA4,AEO,IOPT,IREJ)
C
C     IOPT = 3  q -- aq Jet
C     IOPT = 4  q -- qq Jet
C     IOPT = 5  qq -- aqaq Jet
C     IOPT = 10 q -- q -- q Jet Capella-Kopeliovich
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C     INTERFACE FOR JETSET-7.3 CALL IN DTUNUC 1.04 (J.R. 6/93)
C
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      CHARACTER*8 ANF
      PARAMETER (NFIMAX=249)
      COMMON /DFINPA/ ANF(NFIMAX),PXF(NFIMAX),PYF(NFIMAX),PZF(NFIMAX),
     +HEF(NFIMAX),AMF(NFIMAX), ICHF(NFIMAX),IBARF(NFIMAX),NREF(NFIMAX)
      COMMON /DFINPZ/IORMO(NFIMAX),IDAUG1(NFIMAX),IDAUG2(NFIMAX),
     * ISTATH(NFIMAX)
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
      CHARACTER*8 ANAME
      COMMON/DPAR/ANAME(210),AM(210),GA(210),TAU(210),ICH(210),
     *IBAR(210),K1(210),K2(210)
      COMMON/DIFFRA/ISINGD,IDIFTP,IOUDIF,IFLAGD
      COMMON/CAPKOP/XXX1,XXX3
      DIMENSION NBAMLU(10)
      DIMENSION KORIII(4000)
      COMMON/PROMU/IPROMU
      COMMON/POPCOR/PDB,AJSDEF
      COMMON/IFRAGM/IFRAG
C     COMMON/HARLUN/IHARLU,QLUN
      COMMON /HARLUN/ QLUN,IHARLU
      COMMON /NNCMS/ GAMCM,BGCM,UMO,PCM,EPROJ,PPROJ
      COMMON /JSPART/PXP(1000),PYP(1000),PZP(1000),HEP(1000),NNNP
C--------------------------------------------------
C                       IHARLU=0  soft jets fragmenting
C                       IHARLU=1  hard jets fragmenting
C                                 with final state evolution
C                                 in JEYTSET
C                       QLUN = Mass of hard partons at end of
C                              chain
C-----------------------------------------------------------
      DIMENSION IJOIN(3)
      COMMON /JNI/JNI
      COMMON /NDON/NDONE
      DATA IJOIN/1,2,3/
      DATA IPROMM/0/
      DATA NBAMLU /2,1,3,4,5,6,-2,-1,-3,-4/
      DATA ICOUN/0/
      DATA IWARN/0/
      DATA NBAML/0/
      IF(IPROMU.NE.0)IPROMM=IPROMU
C
C---------------------------------------------------------
C
C               Change Lund parameters depending on energy
C
C---------------------------------------------------------
C
C
C----------------------------------------------------------------
C
C                   DEFINE Lund parameters
C                       IFRAG=1  
C
C----------------------------------------------------------------
C
      IF(IFRAG.EQ.1)THEN   
C--------------------------------------------------------------
C--------------------------------------------------------------
C                                Test LUND default 2.2.99
C--------------------------------------------------------------
C--------------------------------------------------------------
      IF(AEO.LE.15.D0)THEN
C                           Lund b-parameter (default=0.9)
C                     2.2.99
C       PARJ(42)=1.55D0
C                    22.4.99
        PARJ(42)=1.20D0
C                           Lund a-parameter (default=0.5)
C                     2.2.99
c       PARJ(41)=0.14D0
C                    22.4.99
        PARJ(41)=0.30D0
C                           Lund sigma parameter in pt distribution
C                                             (default=0.35)
         PARJ(21)=0.35D0
         PARJ(23)=0.2D0
         PARJ(24)=2.D0
      ELSEIF(AEO.GE.30.D0)THEN
C                           Lund b-parameter (default=0.9)
C                     2.2.99
        PARJ(42)=1.1D0
C                           Lund a-parameter (default=0.5)
C                     2.2.99
        PARJ(41)=0.40D0
C                           Lund sigma parameter in pt distribution
C                                             (default=0.35)
         PARJ(21)=0.38D0
         PARJ(23)=0.2D0
         PARJ(24)=2.D0
      ELSE
        AINTER=(AEO-15.D0)/15.D0
C                           Lund b-parameter (default=0.9)
C                     2.2.99    (das war wohl falsch)
C       PARJ(42)=1.15D0-AINTER*0.45D0
C                    22.4.99
        PARJ(42)=1.20D0-AINTER*0.10D0
C                           Lund a-parameter (default=0.5)
C                     2.2.99
C       PARJ(41)=0.14D0+AINTER*0.15D0 
C                    22.4.99
        PARJ(41)=0.30D0+AINTER*0.10D0 
C                           Lund sigma parameter in pt distribution
C                                             (default=0.35)
         PARJ(21)=0.35D0+AINTER*0.03D0
         PARJ(23)=0.2D0
         PARJ(24)=2.D0
      ENDIF
C                 Probability for meson with spin 1 (d=0.5) 
	PARJ(11)=0.4D0
C                Extra supression for spin 3/2 baryons (d=1.)
	PARJ(18)=0.3D0
C                Extra supression for eta production (d=1.)
C                     taken out since no effect on eta production 5.3.97
C	PARJ(25)=1.3D0
C
C
C----------------------------------------------------------------
C
C                   DEFINE Lund parameters
C                       IFRAG=10  
C                 Low energy test
C
C----------------------------------------------------------------
      ELSEIF(IFRAG.EQ.10)THEN
        IF(AEO.LT.4.D0)THEN
C                          Supress S=3/2 Baryons at low energy
C          PARJ(18)=0.D0
C                           Lund b-parameter (default=0.9)
           IF(PPROJ.LE.30.D0)THEN
             PARJ(42)= 0.6D0
           ELSEIF(PPROJ.GE.100.D0)THEN
             PARJ(42)= 6.D0
           ELSE
             DUPAR=(PPROJ-30.D0)/70.D0
             PARJ(42)=0.6D0+DUPAR*5.4D0
           ENDIF
C                           Popcorn
           MSTJ(12)=2
           PARJ(5)=PDB
C                           Lund sigma parameter in pt distribution
C                                             (default=0.35)
           PARJ(21)=0.65D0
           PARJ(23)=0.2D0
           PARJ(24)=2.D0
C                                     extension 16.12--------------------
        ELSEIF(AEO.LT.7.D0)THEN
C                                    4--7  GeV
           AINTER=-(AEO- 7.D0)/3.D0
C                           Lund b-parameter (default=0.9)
           IF(PPROJ.LE.30.D0)THEN
             DD42=0.25D0
           ELSEIF(PPROJ.GE.100.D0)THEN
             DD42=5.65D0
           ELSE
             DUPAR=(PPROJ-30.D0)/70.D0
             DD42=0.25D0+DUPAR*5.4D0
           ENDIF
           PARJ(42)=0.35D0+AINTER*DD42
C                           Lund sigma parameter in pt distribution
C                                             (default=0.35)
           PARJ(21)=0.65D0
           PARJ(23)=0.2D0
           PARJ(24)=2.D0
C                           Popcorn
           MSTJ(12)=2
           PARJ(5)=PDB
C                                     extension 16.12--------------------
      ELSEIF(AEO.LT.10.D0)THEN
C-----------------------------------end test-15.12-----------------------
C                               Test unsuccessful p-p multiplicities too low
C     IF(AEO.LT.10.D0)THEN
C                           Lund b-parameter (default=0.9)
C        PARJ(42)=6.0         before 4.2.94
C        PARJ(42)=0.5         before 21.3.94
         PARJ(42)=0.35D0
C                           Lund sigma parameter in pt distribution
C                                             (default=0.35)
         PARJ(21)=0.55D0
         PARJ(21)=0.65D0
         PARJ(23)=0.2D0
         PARJ(24)=2.D0
C                           Popcorn
           MSTJ(12)=2
           PARJ(5)=PDB
      ELSEIF(AEO.GE.30.D0)THEN
C                           Lund b-parameter (default=0.9)
C        PARJ(42)=1.1 D0      before 21.3.94 
         PARJ(42)=0.5D0
         PARJ(42)=0.35D0
C                           Lund sigma parameter in pt distribution
C                                             (default=0.35)
         PARJ(21)=0.42D0
         PARJ(21)=0.52D0
         PARJ(23)=0.2D0
         PARJ(24)=2.D0
C                           Popcorn
           MSTJ(12)=2
           PARJ(5)=PDB
      ELSE
         AINTER=-(AEO-30.D0)/20.D0
C                           Lund b-parameter (default=0.9)
C        PARJ(42)=1.4D0+AINTER*4.6D0     before 4.2.94
c        PARJ(42)=1.1D0-AINTER*0.6D0     before 21.3.94
         PARJ(42)=0.5D0-AINTER*0.15D0
         PARJ(42)=0.35D0-AINTER*0.00D0
C                           Lund sigma parameter in pt distribution
C                                             (default=0.35)
         PARJ(21)=0.42D0+AINTER*0.13D0
         PARJ(21)=0.52D0+AINTER*0.13D0
         PARJ(23)=0.2D0
         PARJ(24)=2.D0
C                           Popcorn
           MSTJ(12)=2
           PARJ(5)=PDB
      ENDIF
C
      ENDIF
C     IF(IOPT.EQ.5)THEN
C       PARJ(21)=0.46
C     ELSEIF(IOPT.EQ.3)THEN
C       PARJ(21)=0.46
C     ELSEIF(IOPT.EQ.4)THEN
C       PARJ(21)=0.46
C     ENDIF
C
C---------------------------------------------------------
C
C               Change Lund parameters depending on energy
C                           SINGLE DIFFRACTIVE EVENTS
C
C---------------------------------------------------------
C
      IF(IFLAGD.EQ.1)THEN
C                           Lund b-parameter (default=0.9)
           PARJ(42)=0.23D0
C                           Lund sigma parameter in pt distribution
C                                             (default=0.35)
           PARJ(21)=0.55D0
C        IF(AEO.LE.6.D0)THEN
C                           Lund b-parameter (default=0.9)
C   
C          PARJ(42)=0.15
C        ELSEIF(AEO.GE.30.D0)THEN
C                           Lund b-parameter (default=0.9)
C          PARJ(42)=0.4
C        ELSE
C          AINTER=-(AEO-30.)/24.
C                           Lund b-parameter (default=0.9)
C          PARJ(42)=0.4+AINTER*0.25
C        ENDIF 
C        IF(AEO.LE.6.D0)THEN
C                           Lund sigma parameter in pt distribution
C                                             (default=0.35)
C          PARJ(21)=0.85
C        ELSEIF(AEO.LE.10.D0)THEN
C          AINTER=-(AEO-7.)/4.
C                           Lund sigma parameter in pt distribution
C                                             (default=0.35)
C          PARJ(21)=0.55+AINTER*0.30
C        ELSEIF(AEO.GE.30.D0)THEN
C                           Lund sigma parameter in pt distribution
C                                             (default=0.35)
C          PARJ(21)=0.40
C        ELSE
C          AINTER=-(AEO-30.)/20.
C                           Lund sigma parameter in pt distribution
C                                             (default=0.35)
C          PARJ(21)=0.40+AINTER*0.15
C       ENDIF
      ENDIF
C
C--------------------------------------------------------
C
C--------------------------------------------------------
C
C       Check and modify flavors at chain ends depending on energy
C       and define Lund flavor parameters
C
C--------------------------------------------------------
C
C        IOPT = 3              Quark-Antiquark
C
C--------------------------------------------------------
C
  111 CONTINUE
      IF(IOPT.EQ.3)THEN
	IF((KFA1.LE.0.OR.KFA1.GT.10).OR.
     *     (KFA2.LE.0.OR.KFA2.GT.10))THEN
	  KFA11=KFA1
	  KFA22=KFA2
	  KFA1=1
	  KFA2=7
	  IWARN=IWARN+1
	  IF(IWARN.LE.20)WRITE(6,*)
     *	  ' BAMLUN KFA1,KFA2:',KFA11,KFA22,
     *                            ' changed into :',KFA1,KFA2	  
	ENDIF
	IF(AEO.LE.4.2D0)THEN
	IF((KFA1.EQ.4).AND.
     *     (KFA2.EQ.10))THEN
          KFA11=KFA1
          KFA22=KFA2
	  KFA1=1
	  IWARN=IWARN+1
	  IF(IWARN.LE.20)WRITE(6,*)
     *	  ' BAMLUN KFA1,KFA2 at energy AEO:',KFA11,KFA22,AEO,
     *    ' changed into:',KFA1,KFA2
	ELSEIF((KFA1.EQ.10).AND.
     *     (KFA2.EQ.4))THEN
          KFA11=KFA1
          KFA22=KFA2
	  KFA2=1
	  IWARN=IWARN+1
	  IF(IWARN.LE.20)WRITE(6,*)
     *	  ' BAMLUN KFA1,KFA2 at energy AEO:',KFA11,KFA22,AEO,
     *    ' changed into:',KFA1,KFA2
        ENDIF
        ENDIF
	IF(AEO.LE.2.5D0)THEN
	IF(KFA1.EQ.4
     *     )THEN
          KFA11=KFA1
          KFA22=KFA2
	  KFA1=1
	  IWARN=IWARN+1
	  IF(IWARN.LE.20)WRITE(6,*)
     *	  ' BAMLUN KFA1,KFA2 at energy AEO:',KFA11,KFA22,AEO,
     *    ' changed into:',KFA1,KFA2
	ELSEIF(KFA1.EQ.10
     *     )THEN
          KFA11=KFA1
          KFA22=KFA2
	  KFA1=7
	  IWARN=IWARN+1
	  IF(IWARN.LE.20)WRITE(6,*)
     *	  ' BAMLUN KFA1,KFA2 at energy AEO:',KFA11,KFA22,AEO,
     *    ' changed into:',KFA1,KFA2
	ELSEIF(
     *     KFA2.EQ.4)THEN
          KFA11=KFA1
          KFA22=KFA2
	  KFA2=1
	  IWARN=IWARN+1
	  IF(IWARN.LE.20)WRITE(6,*)
     *	  ' BAMLUN KFA1,KFA2 at energy AEO:',KFA11,KFA22,AEO,
     *    ' changed into:',KFA1,KFA2
	ELSEIF(
     *     KFA2.EQ.10)THEN
          KFA11=KFA1
          KFA22=KFA2
	  KFA2=7
	  IWARN=IWARN+1
	  IF(IWARN.LE.20)WRITE(6,*)
     *	  ' BAMLUN KFA1,KFA2 at energy AEO:',KFA11,KFA22,AEO,
     *    ' changed into:',KFA1,KFA2
	ENDIF
	ENDIF
	IF(((KFA1.EQ.3.OR.KFA1.EQ.9).AND.
     *     (KFA2.EQ.3.OR.KFA2.EQ.9)).AND.AEO.LE.1.5D0)THEN
          KFA11=KFA1
          KFA22=KFA2
	  KFA1=2
	  IWARN=IWARN+1
	  IF(IWARN.LE.20)WRITE(6,*)
     *	  ' BAMLUN KFA1,KFA2 at energy AEO:',KFA11,KFA22,AEO,
     *    ' changed into:',KFA1,KFA2
        ENDIF
	IF(((KFA1.EQ.3.OR.KFA1.EQ.9).OR.
     *     (KFA2.EQ.3.OR.KFA2.EQ.9)).AND.AEO.LE.1.0D0)THEN
          KFA11=KFA1
          KFA22=KFA2
	  KFA1=2
	  KFA2=8
	  IWARN=IWARN+1
	  IF(IWARN.LE.20)WRITE(6,*)
     *	  ' BAMLUN KFA1,KFA2 at energy AEO:',KFA11,KFA22,AEO,
     *    ' changed into:',KFA1,KFA2
	ENDIF
        IF(AEO.LT.0.8D0)THEN
          IREJ=1
	  NBAML=NBAML+1
          IF(NBAML.LT.20)WRITE(6,*)' REJ. IN BAMLUN q-aq A < 0.8',AEO
          RETURN
        ENDIF
        IFLA1=NBAMLU(KFA1)
        IFLA2=NBAMLU(KFA2)
C
C--------------------------------------------------------
C
C            IOPT = 4          Quark-Diquark
C
C--------------------------------------------------------
C
      ELSEIF(IOPT.EQ.4.OR.IOPT.EQ.6)THEN
        IF(AEO.LT.1.3D0)THEN
          IREJ=1
	  NBAML=NBAML+1
          IF(NBAML.LT.20)WRITE(6,*)' REJ. IN BAMLUN q-qq E< 1.5 ',AEO
          RETURN
        ENDIF
	IF((KFA1.LE.0.OR.KFA1.GT.10).OR.
     *     (KFA2.LE.0.OR.KFA2.GT.10).OR. 
     *     (KFA2.LE.0.OR.KFA3.GT.10))THEN
          KFA11=KFA1
          KFA22=KFA2
          KFA33=KFA3
	  KFA1=1
	  KFA2=2
	  KFA3=1
	  IWARN=IWARN+1
	  IF(IWARN.LE.20)WRITE(6,*)
     *	  ' BAMLUN IOPT KFA1,KFA2,KFA3:',
     *	  IOPT,KFA11,KFA22,
     *                 KFA33,' changed into :',KFA1,KFA2,KFA3	  
     *    ,IOPT
	ENDIF
	IF(AEO.LT.3.5D0)THEN
        IF(KFA1.EQ.4
     *    )THEN
          KFA11=KFA1
          KFA22=KFA2
          KFA33=KFA3
          KFA1=1
	  IF(IWARN.LE.20)WRITE(6,*)
     *	  ' BAMLUN KFA1,KFA2,KFA3 at energy AEO:',KFA11,KFA22,KFA33,AEO,
     *    ' changed into:',KFA1,KFA2,KFA3
     *    ,IOPT
        ELSEIF(
     *    KFA2.EQ.4
     *    )THEN
          KFA11=KFA1
          KFA22=KFA2
          KFA33=KFA3
          KFA2=1
	  IF(IWARN.LE.20)WRITE(6,*)
     *	  ' BAMLUN KFA1,KFA2,KFA3 at energy AEO:',KFA11,KFA22,KFA33,AEO,
     *    ' changed into:',KFA1,KFA2,KFA3
     *    ,IOPT
        ELSEIF(
     *    KFA3.EQ.4)THEN
          KFA11=KFA1
          KFA33=KFA3
          KFA22=KFA2
          KFA3=1
	  IWARN=IWARN+1
	  IF(IWARN.LE.20)WRITE(6,*)
     *	  ' BAMLUN KFA1,KFA2,KFA3 at energy AEO:',KFA11,KFA22,KFA33,AEO,
     *    ' changed into:',KFA1,KFA2,KFA3
     *    ,IOPT
        ELSEIF(KFA1.EQ.10
     *    )THEN
          KFA11=KFA1
          KFA22=KFA2
          KFA33=KFA3
          KFA1=7
	  IF(IWARN.LE.20)WRITE(6,*)
     *	  ' BAMLUN KFA1,KFA2,KFA3 at energy AEO:',KFA11,KFA22,KFA33,AEO,
     *    ' changed into:',KFA1,KFA2,KFA3
     *    ,IOPT
        ELSEIF(
     *    KFA2.EQ.10
     *    )THEN
          KFA11=KFA1
          KFA22=KFA2
          KFA33=KFA3
          KFA2=7
	  IF(IWARN.LE.20)WRITE(6,*)
     *	  ' BAMLUN KFA1,KFA2,KFA3 at energy AEO:',KFA11,KFA22,KFA33,AEO,
     *    ' changed into:',KFA1,KFA2,KFA3
     *    ,IOPT
        ELSEIF(
     *    KFA3.EQ.10)THEN
          KFA11=KFA1
          KFA22=KFA2
          KFA33=KFA3
          KFA3=7
	  IWARN=IWARN+1
	  IF(IWARN.LE.20)WRITE(6,*)
     *	  ' BAMLUN KFA1,KFA2,KFA3 at energy AEO:',KFA11,KFA22,KFA33,AEO,
     *    ' changed into:',KFA1,KFA2,KFA3
     *    ,IOPT
        ENDIF
        ENDIF
	IF(AEO.LT.4.2D0)THEN
        IF(KFA1.EQ.4
     *    )THEN
          KFA11=KFA1
          KFA22=KFA2
          KFA33=KFA3
          KFA1=1
	  IF(IWARN.LE.20)WRITE(6,*)
     *	  ' BAMLUN KFA1,KFA2,KFA3 at energy AEO:',KFA11,KFA22,KFA33,AEO,
     *    ' changed into:',KFA1,KFA2,KFA3
     *    ,IOPT
        ELSEIF(
     *    KFA2.EQ.4
     *    )THEN
          KFA11=KFA1
          KFA22=KFA2
          KFA33=KFA3
          KFA2=1
	  IF(IWARN.LE.20)WRITE(6,*)
     *	  ' BAMLUN KFA1,KFA2,KFA3 at energy AEO:',KFA11,KFA22,KFA33,AEO,
     *    ' changed into:',KFA1,KFA2,KFA3
     *    ,IOPT
        ELSEIF(
     *    KFA3.EQ.4)THEN
          KFA11=KFA1
          KFA33=KFA3
          KFA22=KFA2
          KFA3=1
	  IWARN=IWARN+1
	  IF(IWARN.LE.20)WRITE(6,*)
     *	  ' BAMLUN KFA1,KFA2,KFA3 at energy AEO:',KFA11,KFA22,KFA33,AEO,
     *    ' changed into:',KFA1,KFA2,KFA3
     *    ,IOPT
        ELSEIF(KFA1.EQ.10
     *    )THEN
          KFA11=KFA1
          KFA22=KFA2
          KFA33=KFA3
          KFA1=7
	  IF(IWARN.LE.20)WRITE(6,*)
     *	  ' BAMLUN KFA1,KFA2,KFA3 at energy AEO:',KFA11,KFA22,KFA33,AEO,
     *    ' changed into:',KFA1,KFA2,KFA3
     *    ,IOPT
        ELSEIF(
     *    KFA2.EQ.10
     *    )THEN
          KFA11=KFA1
          KFA22=KFA2
          KFA33=KFA3
          KFA2=7
	  IF(IWARN.LE.20)WRITE(6,*)
     *	  ' BAMLUN KFA1,KFA2,KFA3 at energy AEO:',KFA11,KFA22,KFA33,AEO,
     *    ' changed into:',KFA1,KFA2,KFA3
     *    ,IOPT
        ELSEIF(
     *    KFA3.EQ.10)THEN
          KFA11=KFA1
          KFA22=KFA2
          KFA33=KFA3
          KFA3=7
	  IWARN=IWARN+1
	  IF(IWARN.LE.20)WRITE(6,*)
     *	  ' BAMLUN KFA1,KFA2,KFA3 at energy AEO:',KFA11,KFA22,KFA33,AEO,
     *    ' changed into:',KFA1,KFA2,KFA3
     *    ,IOPT
        ENDIF
        ENDIF
        IF(AEO.LT.1.25D0)THEN
          KFA11=KFA1
          KFA22=KFA2
          KFA33=KFA3
          KFA1=1
          KFA2=7
          KFA3=0
	  IWARN=IWARN+1
	  IF(IWARN.LE.20)WRITE(6,*)
     *	  ' BAMLUN KFA1,KFA2,KFA3 at energy AEO:',KFA11,KFA22,KFA33,AEO,
     *    ' changed into:',KFA1,KFA2,KFA3
     *    ,IOPT
          IOPT=3
          GO TO 111
        ENDIF
        IFLA1=NBAMLU(KFA1)
        IFL2=NBAMLU(KFA2)
        IFL3=NBAMLU(KFA3)
        IF(ABS(IFL3).GT.ABS(IFL2))THEN
          IFL22=IFL2
          IFL2=IFL3
          IFL3=IFL22
        ENDIF
        IFLA2=1000*ABS(IFL2)+100*ABS(IFL3)
        IFLZ=3
        IF(ABS(IFL3).LT.ABS(IFL2).AND.RNDM(U).LE.0.25D0)IFLZ=1
        IFLA2=IFLA2+IFLZ
        IF(IFL2.LT.0)IFLA2=-IFLA2
C       IF(IFLA1.EQ.-1.AND.IFLA2.EQ.-3303)IPRI=2
C
C--------------------------------------------------------
C
C            IOPT = 5         DiQuark-AntiDiquark
C
C--------------------------------------------------------
C
      ELSEIF(IOPT.EQ.5)THEN
        IFL1=NBAMLU(KFA1)
        IFL2=NBAMLU(KFA2)
        IF(ABS(IFL2).GT.ABS(IFL1))THEN
          IFL11=IFL1
          IFL1=IFL2
          IFL2=IFL11
        ENDIF
        IFL3=NBAMLU(KFA3)
        IFL4=NBAMLU(KFA4)
        IF(ABS(IFL4).GT.ABS(IFL3))THEN
          IFL33=IFL3
          IFL3=IFL4
          IFL4=IFL33
        ENDIF
        IFLA1=1000*ABS(IFL1)+100*ABS(IFL2)
        IFLA2=1000*ABS(IFL3)+100*ABS(IFL4)
        IFLZ=3
        IF(ABS(IFL2).LT.ABS(IFL1).AND.RNDM(U).LE.0.25D0)IFLZ=1
        IFLA1=IFLA1+IFLZ
        IF(IFL1.LT.0)IFLA1=-IFLA1
        IFLZ=3
        IF(ABS(IFL4).LT.ABS(IFL3).AND.RNDM(U).LE.0.25D0)IFLZ=1
        IFLA2=IFLA2+IFLZ
        IF(IFL3.LT.0)IFLA2=-IFLA2
C
C--------------------------------------------------------
C
C            IOPT = 10        DiQuark-AntiDiquark
C
C--------------------------------------------------------
C
      ELSEIF(IOPT.EQ.10)THEN
        IFLA1=NBAMLU(KFA1)
        IFLA2=NBAMLU(KFA2)
        IFLA3=NBAMLU(KFA3)
      ENDIF
      IF(IPRI.GE.1)WRITE(6,*)IFLA1,IFLA2,AEO,IFLA3
C     WRITE(6,103)IFLA1,IFLA2,AEO,' BAMLUN'
  103 FORMAT(' BAMLUN',2I10,F10.2,I10)
C
C---------------------------------------------------------------------------
C
C               JETSET Call
C
C---------------------------------------------------------------------------
C
C--------------------------------------------------
C                       IHARLU=0  soft jets fragmenting
C                       IHARLU=1  hard jets fragmenting
C                                 with final state evolution
C                                 in JETSET
C                       QLUN = Mass of hard partons at end of
C                              chain
C-----------------------------------------------------------
C     CALL PY2ENT(0,IFLA1,IFLA2,AEO)
C----------------------------------------------------------
      IF(IOPT.EQ.10)THEN
	XX1=XXX1
	XX3=XXX3
	XX2=1.D0-XX1
	ICOU=0
 1234   CONTINUE
	ICOU=ICOU+1
	IF(ICOU.GE.100)THEN
C	  WRITE(6,*)" Stop in dpmlund IOPT=10"
	  STOP
        ENDIF
C                       Select sea q--aq pair
	CALL XSEAPA(AEO,XXXX,ISQ,ISAQ,XSQ,XSAQ,IREJ)
	WRITE(6,*)ISQ,ISAQ,XSQ,XSAQ,XX1
C	IF(XSQ+XSAQ.GE.XX1/2.D0)GO TO 1234
	IF(XSQ.GE.XX1/2.D0.OR.XSAQ.GE.XX3/2.D0)GO TO 1234
C	IF(XSQ.LE.XX2)GO TO 1234
	XX1=XX1-XSQ
	XX3=XX3-XSAQ
        IFLASQ=NBAMLU(ISQ)
        IFLASA=NBAMLU(ISAQ)

	IF(IFLA1.GT.0)THEN
C                      Form diquark out of IFLA2 and IFLASQ
          IF(ABS(IFLASQ).GT.ABS(IFLA2))THEN
            IFl11=IFLA2
            IFLA2=IFLASQ
            IFLASQ=IFL11
          ENDIF
          IFLAD=1000*ABS(IFLA2)+100*ABS(IFLASQ)
          IFLZ=3
          IF(ABS(IFLASQ).LT.ABS(IFLA2).AND.RNDM(U).LE.0.25D0)IFLZ=1
          IFLAD=IFLAD+IFLZ
        	WRITE(6,'(4I10)')IFLA1,IFLAD,IFLASA,IFLA3
	CALL PY4ENT(1,IFLA1,IFLAD,IFLASA,IFLA3,AEO,
     *	XX1,XX2+XSQ,XX3,0.D0,4.D0*XX1*XX3)
C    *	XX1,XX2+XSQ,XX3,4.D0*XX1*(XX2+XSQ),4.D0*XX1*XX3)
C     *	XX1,XX2+XSQ,XX3,(XX1-XX2-XSQ)**2,4.D0*XX1*XX3)
          IJOIN(1)=1 
	IJOIN(2)=3 
	CALL PYJOIN(2,IJOIN) 
        IJOIN(1)=2 
	IJOIN(2)=4 
	CALL PYJOIN(2,IJOIN) 
	CALL PYEXEC
	ELSEIF(IFLA1.LT.0)THEN
C                      Form anti-diquark out of IFLA2 and IFLASA
        IF(ABS(IFLASA).GT.ABS(IFLA2))THEN
          IFl11=IFLA2
          IFLA2=IFLASA
          IFLASA=IFL11
        ENDIF
        IFLAD=1000*ABS(IFLA2)+100*ABS(IFLASA)
        IFLZ=3
        IF(ABS(IFLASA).LT.ABS(IFLA2).AND.RNDM(U).LE.0.25D0)IFLZ=1
        IFLAD=IFLAD+IFLZ
        IF(IFLA2.LT.0)IFLAD=-IFLAD
	WRITE(6,'(4I10)')IFLA1,IFLAD,IFLASQ,IFLA3
	CALL PY4ENT(1,IFLA1,IFLAD,IFLASQ,IFLA3,AEO,
     *  XX1,XX2+XSAQ,XX3,0.D0,XX1*XX3)
        IJOIN(1)=1 
	IJOIN(2)=3 
	CALL PYJOIN(2,IJOIN) 
        IJOIN(1)=2 
	IJOIN(2)=4 
	CALL PYJOIN(2,IJOIN) 
	CALL PYEXEC
	ENDIF
      ELSE
        IF(IHARLU.EQ.0)THEN
          CALL PY2ENT(0,IFLA1,IFLA2,AEO)
          PXP(1)=0.
          PYP(1)=0.
          PZP(1)=AEO/2.D0
          HEP(1)=AEO/2.D0
          PXP(2)=0.D0
          PYP(2)=0.D0
          PZP(2)=AEO/2.D0
          HEP(2)=AEO/2.D0
          NNNP=2
        ELSEIF(IHARLU.EQ.1)THEN
          CALL PY2ENT(-1,IFLA1,IFLA2,AEO)
          CALL PYSHOW(1,2,QLUN)
          DO 201 I=1,N
            IF(K(I,1).EQ.3)THEN
              NNNP=NNNP+1
              IF(NNNP.LT.1000)THEN
              PXP(NNNP)=P(I,1)
              PYP(NNNP)=P(I,2)
              PZP(NNNP)=P(I,3)
              HEP(NNNP)=P(I,4)
              ENDIF
            ENDIF
  201     CONTINUE
          IF(IPRI.GE.2)CALL PYLIST(1)
          CALL PYEXEC
C         CALL PYLIST(1)
        ENDIF
      ENDIF
C                                     Force all particles wanted to decay
      DO 1111 IIII=1,N
        IF(K(IIII,1).EQ.4)K(IIII,1)=5
 1111 CONTINUE
      CALL PYEXEC
C
      IF(IPRI.GE.2)WRITE(6,*)' After PYEXEC'
      IF(IPRI.GE.2)CALL PYLIST(1)
C         CALL PYLIST(1)
C         IF(IPRI.GE.2)CALL PYLIST(1)
C
C---------------------------------------------------------------------------
C
C               Edit JETSET event
C
C---------------------------------------------------------------------------
C
      CALL PYEDIT(12)
      ICOUN=ICOUN+1
      IF(IPRI.GE.2)WRITE(6,*)' After PYEDIT'
C         CALL PYLIST(1)
      IF(IPRI.GE.2)CALL PYLIST(1)
      IF(IHARLU.EQ.1.AND.NDONE.EQ.-107801)THEN
        WRITE(6,*)'NDONE ',NDONE
        CALL PYLIST(1)
      ENDIF
C
C---------------------------------------------------------------------------
C
C               Move JETSET event into BAMJET COMMON
C
C---------------------------------------------------------------------------
C
      IHAD=0
C                      EXISTING PARTICLES
         IF(IPRI.GE.2) WRITE(6,*)' DPMJET COMMON  particles'
      IORII=0	  
      KORJJJ=0
      KORKKK=-2
      DO 101 I=1,N
C                            fragmented string      
        IF((K(I,1).EQ.11).AND.(K(I,2).EQ.92))THEN
C                              J.R.    11.2.2000	
C---------------------------------------------------
C         IF(KORJJJ.NE.0)GO TO 101
C---------------------------------------------------
	  IORII=999
	  KORII=I
C         KORJJJ is the first string	  
	  IF(KORJJJ.EQ.0)KORJJJ=I
	ENDIF  
C                            fragmented cluster      
        IF((K(I,1).EQ.11).AND.(K(I,2).EQ.91))THEN
	  IORII=999
	  KORII=I
C         KORJJJ is the first string	  
	  IF(KORJJJ.EQ.0)KORJJJ=I
	ENDIF  
C       KORIII(I) is the string from which particle I comes	
	KORIII(I)=KORII
C       undecayed particle, decayed particle        
	IF((K(I,1).EQ.1).OR.(K(I,1).EQ.4).OR.(K(I,1).EQ.11)
     *	                .OR.(K(I,1).EQ.15))THEN
C                    strings and clusters     
          IF((K(I,2).EQ.91).OR.(K(I,2).EQ.92).OR.(K(I,2).EQ.94))THEN
            KORKKK=KORKKK+1
	    GO TO 101
          ENDIF
C                        quarks	  
          IF(ABS(K(I,2)).LE.6) GO TO 101
          IF(ABS(K(I,2)).GE.1000) THEN
	    KOO=ABS(K(I,2))
	    KAA=KOO/100
	    KLL=KOO-KAA*100
C                       diquarks	    
	    IF(KLL.EQ.1.OR.KLL.EQ.3)GO TO 101
	  ENDIF
          IHAD=IHAD+1
          IF(IHAD.GT.NFIMAX)THEN
            WRITE(6,1112)IHAD
 1112       FORMAT(' BAMLUN: IHAD.GT.NFIMAX INCREASE NFIMAX IHAD=',I10)
            GO TO 101 
          ENDIF
	  IF(K(I,3).EQ.KORIII(I))THEN
C           Hadrons from string
            IORMO(IHAD)=IORII
	  ELSE
C           Hadrons from decaying particle	  
	    IF(KORJJJ.EQ.KORIII(K(I,3)))THEN
C             IORMO(IHAD)=K(I,3)-KORJJJ
C                                     J.R.15.2.2000
              IORMO(IHAD)=K(I,3)-KORJJJ
	    ELSE
C             IORMO(IHAD)=K(I,3)-KORJJJ-1
C                                     J.R.15.2.2000
              IORMO(IHAD)=K(I,3)-KORJJJ-1-KORKKK
	    ENDIF  
	  ENDIF
          PXF(IHAD)=P(I,1)
          PYF(IHAD)=P(I,2)
          PZF(IHAD)=P(I,3)
          HEF(IHAD)=P(I,4)
          NREFB=MCIHAD(K(I,2))
	  AMASSS=PYMASS(K(I,2))
          IF(IPRI.GE.3)WRITE(6,'(4I10)')I,K(I,1),K(I,2),NREFB
          IF(NREFB.LT.1.OR.NREFB.GT.183)NREFB=4
          NREF(IHAD)=NREFB
          ANF(IHAD)=ANAME(NREFB)
C         AMF(IHAD)=AM(NREFB)
          AMF(IHAD)=AMASSS
C	  WRITE(6,*)' IHAD,AMF(IHAD),AMASSS ',IHAD,AMF(IHAD),AMASSS
          ICHF(IHAD)=ICH(NREFB)
          IBARF(IHAD)=IBAR(NREFB)
	  ISTATH(IHAD)=1
          IF(K(I,1).EQ.11.OR.K(I,1).EQ.15)THEN
	    IBARF(IHAD)=500
	    ISTATH(IHAD)=2
	  ENDIF  
C
C         IF(NREFB.EQ.31.OR.NREFB.EQ.95)THEN
C	  WRITE(6,*)'K(I,1),IBARF(IHAD),ISTATH(IHAD)',
C    *	  K(I,1),IBARF(IHAD),ISTATH(IHAD)
C	  WRITE(6,*)ICOUN,IHAD,K(I,2),NREFB,K(I,3),
C    *	  ISTATH(IHAD),IORMO(IHAD),ANF(IHAD),
C    *                PXF(IHAD),PYF(IHAD),
C    *                PZF(IHAD),HEF(IHAD),AMF(IHAD),' BAMLUN'
C         ENDIF
          IF(NDONE.EQ.-107801.AND.IHARLU.EQ.1)THEN
	  WRITE(6,*)I,ICOUN,IHAD,K(I,2),NREFB,K(I,3),
     *	  ISTATH(IHAD),IBARF(IHAD),IORMO(IHAD),ANF(IHAD),
     *                PXF(IHAD),PYF(IHAD),
     *                PZF(IHAD),HEF(IHAD),AMF(IHAD),KORIII(I),
     *     KORJJJ,' BAMLUN'
          ENDIF
  102     FORMAT(' BAMLUN',2I5,5I10,A8,5E12.3)
C
        ENDIF
  101 CONTINUE
C                            DECAYED RESONAMCES TO BE PUT INTO EVENT RECORD
C                            NOTE: PARTICLES with IBARF=500 are decayed
      IF(JNI.NE.7)GO TO 777
      IF(JNI.NE.7)THEN
C         WRITE(6,*)' DPMJET COMMON decayed particles'
      IORII=0	  
      DO 105 I=1,N
        IF(K(I,1).EQ.11.AND.K(I,2).EQ.92)THEN
	  IORII=999
	  KORII=I
	ENDIF  
        IF(K(I,1).EQ.11.OR.K(I,1).EQ.15)THEN
          IF(K(I,2).EQ.92) GO TO 105
          IF(K(I,2).EQ.94) GO TO 105
          IF(ABS(K(I,2)).LE.6) GO TO 105
          IHAD=IHAD+1
          IF(IHAD.GT.NFIMAX)THEN
            WRITE(6,1112)IHAD
            GO TO 107 
          ENDIF
	  IF(K(I,3).EQ.KORII)THEN
C           Hadrons from string
            IORMO(IHAD)=IORII
	  ELSE
	    IORMO(IHAD)=K(I,3)-KORII
	  ENDIF
          PXF(IHAD)=P(I,1)
          PYF(IHAD)=P(I,2)
          PZF(IHAD)=P(I,3)
          HEF(IHAD)=P(I,4)
          NREFB=MCIHAD(K(I,2))
	  AMASSS=PYMASS(K(I,2))
          IF(IPRI.GE.3)WRITE(6,'(5I10)')I,K(I,1),K(I,2),K(I,3),NREFB
          IF(NREFB.LT.1.OR.NREFB.GT.183)NREFB=4
          NREF(IHAD)=NREFB
          ANF(IHAD)=ANAME(NREFB)
C         AMF(IHAD)=AM(NREFB)
          AMF(IHAD)=AMASSS
          ICHF(IHAD)=ICH(NREFB)
          IBARF(IHAD)=500
C
          IF(IPRI.GE.3)WRITE(6,102)ICOUN,IHAD,K(I,2),NREFB,K(I,3),
     *	  IORMO(IHAD),ANF(IHAD),
     *                PXF(IHAD),PYF(IHAD),
     *                PZF(IHAD),HEF(IHAD),AMF(IHAD)
C
        ENDIF
  105 CONTINUE
      ENDIF
  777 CONTINUE
  107 CONTINUE
C
C---------------------------------------------------------------------------
C
      IF(IFLA1.EQ.-1.AND.IFLA2.EQ.-3303)IPRI=0
      RETURN
      END
      SUBROUTINE XSEAPA(ECM,XXXX,IPSQ1,IPSAQ1,XPSQ1,XPSAQ1,IREJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      COMMON /SEASU3/SEASQ
          IREJ=0
          IF(ECM.LT.20.D0)THEN
            XPTHRO=1.5D0*LOG10(ECM/200.D0)+3.5D0
          ELSEIF(ECM.GE.20.D0)THEN
            XPTHRO=5.0D0
          ENDIF
          XPTHR=1.5D0*XPTHRO/(ECM**1.5D0*14.D0)
C         XSTHR2=1.5D0*XPTHR2/(ECM**1.5D0*14.D0)
C         WRITE(6,*)' XSEAPA:XPTHR ',XPTHR
	  I=1
          AI=I-1
          XPTHRX=XPTHR-0.5D0*AI/ECM**2
          IF (XPTHRX.LT.4.D0/ECM**2)XPTHRX=4.D0/ECM**2
C---------------------------------------------------------------
C                              SPLIT GLUON INTO TWO SEA QUARKS
C                              FLAVORS OF SEA QUARKS
C         INCREASE s-quark fraction to acount for larger rejection
          SEASQQ=2.D0*SEASQ
          AI=I
          BI=I+I
          IPSQ1=1.+RNDM(AI)*(2.D0+SEASQQ)
          IPSAQ1=-IPSQ1
          IPSAQ1=6+IPSQ1
C         WRITE(6,*)' XSEAPA SEASQQ,IPSQ1,IPSAQ1',SEASQQ,IPSQ1,IPSAQ1
C                              X-FRAXTIONS OF SEA QUARKS
C------------------------------------------------------j.r.29.4.93
	          SOX1=0.3D0
              NCOU=0
  500         CONTINUE
              NCOU=NCOU+1
	      IF(IPSQ1.EQ.3)THEN
                IF(NCOU.GE.200)THEN
	 	  IREJ=3
	          RETURN
		ENDIF
              ELSE		
                IF(NCOU.GE.50)THEN
	          IREJ=1
	          RETURN
		ENDIF  
	      ENDIF
              XGLU1=SAMPEX(XPTHRX,SOX1)
          IF(IPSQ1.LE.2)THEN
            XPSQ1=(0.2D0+(0.36D0*RNDM(AI))**0.50D0)*XGLU1
            XPSAQ1=XGLU1-XPSQ1
          ELSEIF(IPSQ1.EQ.3)THEN
C           CHANGE J.R. 18.2.99
            IF (XPSQ1.LE.0.3D0/ECM)GO TO 500
            IF (XPSAQ1.LE.0.3D0/ECM)GO TO 500
          ENDIF
          IF(XPSAQ1.GE.XXXX)GO TO 500
	  RETURN
	  END
C------------------------------------------------------------
C             was dpmtcbsh.f
C------------------------------------------------------------
C
C     file with changes of r.e. 20.12.93 old file: d4tcbshm.fold
*-- Author :
C        PROGRAM SHMAKOV    11/1986
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  THE PROGRAM SIMULATES THE STARS OF THE SECONDARY PARTICLES IN NUCLEUS
C  -NUCLEUS COLLISIONS AT HIGH ENERGIES. THE SIMULATION IS BASED ON
C  GLAUBER APPROACH. EACH STAR IS MADE UP AS THE SUM OVER THE ORDERED N-
C  N STARS. MESON PRODUCTION IS ONLY TAKEN INTO ACCOUNT.
C       LET'S DISCRIBE THE MAIN PARAMETERS OF THE PROGRAM IN FORM:
C  FIRST LIST OF THE PARAMETERS - MEANING - (BOUNDS); AND SO ON.
C  NA,NB - AMOUNT OF NUCLONS IN NUCLEI A AND B - (GT.1,LT.41,
C  PREFERENTIALLY NA LT. NB)
C  NCA,NCB - AMOUNT OF PROTONS IN THE NUCLEI
C  RA,RB - THE RADII TO USE INTO FORMULARS FOR THE NUCLEAR DENSITIES -(F
C  SIG,G,RO - THE PARAMETERS OF THE N-N ELASTIC AMPLITUDE,
C  AMLITUDE(B)=SIG*G*(1-I*RO)*EXP(-G*B**2)/2PI - (FM**2,FM**(-2) AND
C  DIMENSIONLESS ACCORDIANLY)
C  SE2 - QUADRATE OF THE TOTAL 4-MOMENTUM BOTH NUCLEI - (GEV**2,
C  GT. (NA+NB)**2,PREFERENTIALLY GT. 4*(NA+NB)**2)
C  TF01,TF02 - FERMI ENERGIES IN THE CENTERS OF THE NUCLEI - (GEV)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE SHMAKI(NA,NCA,NB,NCB,RPROJ,RTARG,PPN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C*** SHMAKOV INITIALIZATION
*KEEP,DSHM.
      COMMON /DSHM/ RASH,RBSH,BMAX,BSTEP,SIGSH,ROSH,GSH,
     *              BSITE(0:1,200),NSTATB,NSITEB
      COMMON /DSHMS/ SIGSHS
*KEEP,RTAR.
      COMMON /RTAR/ RTARNU
*KEEP,NUCC.
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
*KEND.
      COMMON /ZENTRA/ ICENTR
      COMMON /SIGLA/SIGLAU
      COMMON /KGLAUB/JGLAUB
C--------------------------------------
      RNA=NA
      RNB=NB
      RASH=1.12D0*RNA**0.33D0
      RBSH=1.12D0*RNB**0.33D0
      IF(JGLAUB.EQ.1)THEN
        IF(NA.EQ.9)RASH=2.52D0
        IF(NA.EQ.10)RASH=2.45D0
        IF(NA.EQ.11)RASH=2.37D0
        IF(NA.EQ.12)RASH=2.45D0
        IF(NA.EQ.13)RASH=2.44D0
        IF(NA.EQ.14)RASH=2.55D0
        IF(NA.EQ.15)RASH=2.58D0
        IF(NA.EQ.16)RASH=2.71D0
        IF(NA.EQ.17)RASH=2.66D0
        IF(NA.EQ.18)RASH=2.71D0
        IF(NB.EQ.9)RBSH=2.52D0
        IF(NB.EQ.10)RBSH=2.45D0
        IF(NB.EQ.11)RBSH=2.37D0
        IF(NB.EQ.12)RBSH=2.45D0
        IF(NB.EQ.13)RBSH=2.44D0
        IF(NB.EQ.14)RBSH=2.55D0
        IF(NB.EQ.15)RBSH=2.58D0
        IF(NB.EQ.16)RBSH=2.71D0
        IF(NB.EQ.17)RBSH=2.66D0
        IF(NB.EQ.18)RBSH=2.71D0
      ENDIF
      WRITE(6,*)' SHMAKI: RASH, RBSH = ',RASH,RBSH
      RPROJ=RASH
      RTARG=RBSH
      RTARNU=RBSH
      NSTATB=2000
      IF((ICENTR.EQ.1).AND.(NA.GE.200).AND.(NB.GE.200))THEN
	NSTATB=1000
      ENDIF
      NSITEB=200
      WRITE(6, 1010)NSTATB,NSITEB
 1000 FORMAT(2I10)
 1010 FORMAT(' STATISTIC ON POINT ON B AND NUMBER OF POINTS',2I10)
C*** PARAMETERS IN NN-AMP.
      SIGSH=4.3
      PPNN=PPN
      WRITE(6,*)' SMAKI, PPN = ',PPN
      SIGSH=DSHPTO(IJPROJ,PPN)/10.D0
C     IF(NA.GT.0.D0)
C    *SIGSH=(DSHPTO(IJPROJ,PPNN)-SIGSDS(PPNN*2.D0))/10.D0
      IF(JGLAUB.EQ.2)THEN
        SIGSHS=(DSHPTO(IJPROJ,PPNN)-SIGSDS(PPNN*2.D0))/10.D0
      ELSEIF(JGLAUB.EQ.1)THEN
        SIGSHS=(DSHPTO(IJPROJ,PPNN))/10.D0
      ENDIF
      IF(IJPROJ.EQ.5)THEN
      SIGSHS=(DSHPTO(1,PPNN)-SIGSDS(PPNN*2.D0))/10.D0
      ENDIF
      SIGSH=SIGSHS
      GSH=1.6
      BSLOPE=10.D0
      SSS=2.*PPN
      ECM=SQRT(SSS)
      IF (IJPROJ.LE.12)BSLOPE=8.5D0*(1.D0+0.065D0*LOG(SSS))
      IF (IJPROJ.GT.12)BSLOPE=6.D0*(1.D0+0.065D0*LOG(SSS))
      GSH=1.D0/(2.D0*BSLOPE*0.038938D0)
      ROSH=-0.43D0
      IF (IJPROJ.LE.12)THEN
      IF(ECM.GT.3.0D0.AND.ECM.LE.50.D0) ROSH=-0.63D0+
     +                               0.175D0*LOG(ECM)
      IF(ECM.GT.50.) ROSH=0.1D0
      ENDIF
      IF (IJPROJ.GT.12) ROSH=0.01D0
      WRITE(6, 1030)SIGSH,ROSH,GSH,BSLOPE,ECM
 1020 FORMAT(3F10.5)
 1030 FORMAT(' PARAMETERS OF THE NN AMPLITUDE SIG,RO,G,BSLOPE,ECM ' /5
     +(1PE12.5))
      CALL TITLE(NA,NB,NCA,NCB)
      WRITE(6,*)' vor PREVIO(RA,RB,NSTB,BMAX,BSTEP,SIG,RO,G)',
     &RASH,RBSH,NSITEB,BMAX,BSTEP,SIGSH,ROSH,GSH
      CALL PREVIO(RASH,RBSH,NSITEB,BMAX,BSTEP,SIGSH,ROSH,GSH)
      WRITE(6,*)' SHMAKI: RASH, RBSH = ',RASH,RBSH
      CALL PROFB(BSTEP,NSTATB,NA,RASH,NB,RBSH,BSITE,NSITEB)
      WRITE(37,'(4I5,F15.2,E15.5,F15.2)')
     * NA,NCA,NB,NCB,SIGSH,PPN,SIGLAU
      RETURN
      END
      DOUBLE PRECISION FUNCTION SIGSDS(S)
C                  Approximate sigma_SD in p-p (mb)      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      SIGSDS=4.D0+1.2D0*LOG10(S)
      RETURN
      END
*-- Author :
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE SHMAKF(NA,NCA,NB,NCB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C***     SHMAKOV INITIALIZATION
C***     To be called for each material separately
C***     Input of Shmakov data from file:
C***           1) 'NUCLEUS'      IT, ITZ
C***           2) BMAX, BSTEP, RASH, RBSH
C***           3) BSITE (energy and projectile dependence)
C-----------------------
      CHARACTER*10 BNUC
*KEEP,DSHM.
      COMMON /DSHM/ RASH,RBSH,BMAX,BSTEP,SIGSH,ROSH,GSH,
     *              BSITE(0:1,200),NSTATB,NSITEB
      COMMON /DSHMS/ SIGSHS
*KEEP,RTAR.
      COMMON /RTAR/ RTARNU
*KEEP,NUCC.
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
*KEEP,DTUMAT.
      COMMON /DTUMAT/ BSITEN(200,24,50),BSITEM(200,24,50),
     +                RPROJJ(50),RTARGG(50),BSTEPP(50),BMAXX(50),
     +                NTAXX(50),NZTAXX(50),NPRXX(50),NZPRXX(50)
*KEEP,DAMP.
      COMPLEX*16 CA,CI
      COMMON /DAMP/   CA,CI,GA
*KEND.
      COMMON /KKMATU/KKMATO,KKMATA,IPOO,IPOA,IPZOO,IPZOA
     *,IBPROO,IBPROA,IREADO
      DIMENSION HELP(200)
      DATA MATNUM /0/
C--------------------------------------
      KKMATO=0
      IPOO=0
      IPZOO=0
      IBPROO=0
      IREADO=0
      REWIND 47
      MATNUM=MATNUM + 1
      IF(MATNUM.GT.50) THEN
      WRITE(6,'(2A,I3/A)')
     &  ' Too large number of materials requested for Glauber',
     &  ' initialization in SHMAKF / MATNUM=',MATNUM,
     &  ' execution stopped in SHMAKF'
      STOP
      ENDIF
      WRITE(6,'(A,I3,A)') ' Read Glauber data for material no.',MATNUM,
     &                    ' from unit 47'
C----------------------------------Read Glauber data----------------
      DO 72 I=1,100000
        READ(47,'(A10)',END=79) BNUC
        IF(BNUC.EQ.' NUCLEUS  ') THEN
          BACKSPACE 47
          READ(47,'(A10,4I10)') BNUC,NPRX,NZPRX,NTAX,NZTAX
          IF(NB.EQ.NTAX.AND.NCB.EQ.NZTAX.AND.NA.EQ.NPRX.
     &	  AND.NCA.EQ.NZPRX) THEN
          NTAXX(MATNUM)=NTAX
          NZTAXX(MATNUM)=NZTAX
	  NPRXX(MATNUM)=NPRX
	  NZPRXX(MATNUM)=NZPRX
        READ(47,'(4F10.5)') BMAXX(MATNUM),BSTEPP(MATNUM),
     &                        RPROJJ(MATNUM),RTARGG(MATNUM)
        DO 170 IE=1,24
            READ(47,'(5E16.8)') (BSITEN(IDA,IE,MATNUM),IDA=1,200)
 170      CONTINUE
          IF(NPRX.EQ.1)THEN
            DO 171 IE=1,24
              READ(47,'(5E16.8)') (BSITEM(IDA,IE,MATNUM),IDA=1,200)
 171        CONTINUE
          ENDIF
        GOTO 78
        ENDIF
      ENDIF
 72   CONTINUE
 79   CONTINUE
      WRITE(6,'(A)') ' GLAUBER DATA NOT FOUND'
      STOP
 
 78   CONTINUE
      DO 80 I=1,200
      HELP(I)=I*BSTEPP(MATNUM)
      WRITE (6,1040) HELP(I),(BSITEN(I,IE,MATNUM),IE=1,24)
 1040   FORMAT (F10.4,10(1PE12.4)/10(1PE12.4))
   80 CONTINUE
C-------------------------------------------------------------------
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE SHMAKO(NA,NB,B,INTT,INTA,INTB,JS,JT,PPN,KKMAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER(NAMX=248)
      PARAMETER (INTMX=2488,INTMD=252)
      DIMENSION JS(NAMX),JT(NAMX)
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEEP,NUCKOO.
      COMMON /NUCKOO/ PKOO(3,INTMX),TKOO(3,INTMX),PPOO(3,INTMX),
     +TPOO(3,INTMX)
*KEEP,DSHM.
      COMMON /DSHM/ RASH,RBSH,BMAX,BSTEP,SIGSH,ROSH,GSH,
     *              BSITE(0:1,200),NSTATB,NSITEB
      COMMON /DSHMS/ SIGSHS
*KEEP,NUCC.
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
*KEEP,INTMX.
*KEEP,SHMAKL.
C     INCLUDE (SHMAKL)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON/SHMAKL/JSSH(INTMX),JTSH(INTMX),INTER1(INTMX),INTER2(INTMX)
*KEEP,DTUMAT.
      COMMON /DTUMAT/ BSITEN(200,24,50),BSITEM(200,24,50),
     +                RPROJJ(50),RTARGG(50),BSTEPP(50),BMAXX(50),
     +                NTAXX(50),NZTAXX(50),NPRXX(50),NZPRXX(50)
*KEEP,RTAR.
      COMMON /RTAR/ RTARNU
*KEEP,DAMP.
      COMPLEX*16 CA,CI
      COMMON /DAMP/   CA,CI,GA
*KEEP,RPTSHM.
      COMMON /RPTSHM/ RPROJ,RTARG,BIMPAC
*KEND.
      COMMON /KKMATU/KKMATO,KKMATA,IPOO,IPOA,IPZOO,IPZOA
     *,IBPROO,IBPROA,IREADO
      COMMON /KGLAUB/JGLAUB
C------------------------------------------
*---use predefined values for BSITE array and redefine necessary parameters
*---if initialization by SHMAKI (KKMAT=0) all information directly available
      KKMATA=KKMAT
      IPOA=IP
      IPZOA=IPZ
      IF(IP.GE.2)IJPROJ=1 
      IF(IP.GE.2)IBPROJ=1 
      IBPROA=IBPROJ
      IF(KKMAT.GT.0) THEN
      RASH=RPROJJ(KKMAT)
      RBSH=RTARGG(KKMAT)
      RPROJ=RASH
      RTARG=RBSH
      RTARNU=RBSH
      BMAX=BMAXX(KKMAT)
      BSTEP=BSTEPP(KKMAT)
C     SIGSH=DSHPTO(IJPROJ,PPN)/10.
C     IF(NA.GT.2.D0)
C    *SIGSH=(DSHPTO(IJPROJ,PPN)-SIGSDS(PPN*2.D0))/10.D0
      IF(JGLAUB.EQ.2)THEN
        SIGSHS=(DSHPTO(IJPROJ,PPN)-SIGSDS(PPN*2.D0))/10.D0
      ELSEIF(JGLAUB.EQ.1)THEN
        SIGSHS=(DSHPTO(IJPROJ,PPN))/10.D0
      ENDIF
      IF(IJPROJ.EQ.5)THEN
        SIGSHS=(DSHPTO(1,PPN)-SIGSDS(PPN*2.D0))/10.D0
      ENDIF
      SIGSH=SIGSHS
      BSLOPE=10.D0
      SSS=2.*PPN
      ECM=SQRT(SSS)
      IF (IJPROJ.LE.12)BSLOPE=8.5D0*(1.D0+0.065D0*LOG(SSS))
      IF (IJPROJ.GT.12)BSLOPE=6.D0*(1.D0+0.065D0*LOG(SSS))
*     GSH=1.6D0
      GSH=1.D0/(2.D0*BSLOPE*0.038938D0)
      ROSH=-0.43D0
      IF (IBPROJ.LE.12)THEN
        IF(ECM.GT.3.0D0.AND.ECM.LE.50.D0) ROSH=-0.63D0+
     +                                 0.175D0*LOG(ECM)
        IF(ECM.GT.50.D0) ROSH=0.1
      ENDIF
      IF (IJPROJ.GT.12) ROSH=0.01
      RCA=GSH*SIGSHS/6.2831854D0
      FCA=-GSH*SIGSHS*ROSH/6.2831854D0
      CA=CMPLX(RCA,FCA)
      GA=GSH
      IREAD=0
      WU10=SQRT(10.0)
      DO 73 IPO=1,24
        PPO=WU10**(IPO+1) + 10.
        IF(PPN.LE.PPO) THEN
          IREAD=IREAD + IPO - 1
          GOTO 74
        ENDIF
 73     CONTINUE
 74     CONTINUE
        IF(IREAD.LE.0)IREAD=1
        IF(IREAD.GE.25)IREAD=24
C       old (before 10/98 IF, which might be correct for HEMAS interface
C       INTER_DPM calling DPMEVT calling KKINC (to be checked)
C       IF (KKMAT.NE.KKMATO.AND.IP.NE.IPOO.AND.
C    *  IJPROJ.NE.IBPROO.AND.IPZ.NE.IPZOO )THEN
        IF ((KKMAT.NE.KKMATO).OR.(IP.NE.IPOO).OR.
     *  (IBPROJ.NE.IBPROO).OR.(IPZ.NE.IPZOO).OR.(IREAD.NE.IREADO) )THEN
          IF(IBPROJ.NE.0) THEN
            DO 75 II=1,200
              BSITE(1,II)=BSITEN(II,IREAD,KKMAT)
 75         CONTINUE
	    IBPROO=IBPROJ
C           ADDED 10/98
            KKMATO=KKMAT
            IREADO=IREAD
          ELSE
            DO 76 II=1,200
              BSITE(1,II)=BSITEM(II,IREAD,KKMAT)
 76         CONTINUE
	    IBPROO=IBPROJ
C           ADDED 10/98
            KKMATO=KKMAT
            IREADO=IREAD
          ENDIF
        ENDIF
      ENDIF
      IF(IPEV.GE.6) THEN
      WRITE(6,'(A)')  ' SHMAKO - BEFORE DIAGR'
      IF(IPRI.GE.6) WRITE(6, 1030)SIGSH,ROSH,GSH,BSLOPE,ECM
 1030 FORMAT(' PARAMETERS OF THE NN AMPLITUDE SIG,RO,G,BSLOPE,ECM ' /5
     +(1PE12.5))
      ENDIF
      CALL DIAGR(NA,NB,B,JS,JT,INTT,INTA,INTB)
      IF(IPEV.GE.6) WRITE(6,1000)NA,NB,B,INTT,INTA,INTB,JS(1),JT(1),
     +              PKOO(1,1),TKOO(1,1)
 1000 FORMAT(' NA,NB,B,INTT,INTA,INTB,JS(1),JT(1),PKOO(1,1),TKOO(1,1)
     + IN SHMAKO '/2I5,F10.4,5I6,2F10.3)
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE TITLE(NA,NB,NCA,NCB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      WRITE(6,1000)NA,NCA,NB,NCB
 1000 FORMAT(//10X,39HGLAUBER S APPROACH IS USED TO SIMULATE ,
     +26HNUCLEUS-NUCLEUS COLLISIONS/ 24X,
     +40HTHE CALCULATION NAS BEEN CARRIED OUT FOR/
     +27H PROJECTED NUCLEI WITH A = ,I5,12H CHARGE A = ,I5/
     +24H TARGET NUCLEI WITH B = ,I5,12H CHARGE B = ,I5///)
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE CONUCL(X,N,R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C     -----------------------------------
      COMMON /KGLAUB/JGLAUB
      DIMENSION X(3,N),WD(4),RD(3)
      LOGICAL ISTART
      DATA SQR2/1.414216D0/
      DATA PDIF/0.545D0/,R2MIN/0.16D0/
      DATA WD/0.D0,0.178D0,0.465D0,1.D0/
      DATA RD/2.09D0,0.935D0,0.697D0/
C     DATA RC12/1.6976D0/,WC12/0.4444444D0/
      AAN=N
      IF (N.EQ.1)THEN
        GOTO 10
      ELSEIF (N.EQ.2)THEN
        GOTO 20
      ELSEIF (N.EQ.4) THEN
        GOTO 60
      ELSEIF (N.GE.12) THEN
        GOTO 110
      ELSE
        GO TO 110
      ENDIF
   10 CONTINUE
      X(1,1)=0
      X(2,1)=0
      X(3,1)=0
      RETURN
   20 CONTINUE
      EPS=RNDM(V)
      DO 30 I=1,3
        IF ((EPS.GE.WD(I)).AND.(EPS.LE.WD(I+1)))                 GOTO 40
 
   30 CONTINUE
   40 CONTINUE
      DO 50 J=1,3
        CALL RANNOR(X1,X2)
        X(J,1)=RD(I)*X1
        X(J,2)=-X(J,1)
   50 CONTINUE
      RETURN
   60 CONTINUE
      SIGMA=R/SQR2
   70 CONTINUE
      ISTART=.TRUE.
      CALL RANNOR(X3,X4)
      DO 100 I=1,N
        CALL RANNOR(X1,X2)
        X(1,I)=SIGMA*X1
        X(2,I)=SIGMA*X2
        IF (ISTART)                                              GOTO 80
 
        X(3,I)=SIGMA*X4
        CALL RANNOR(X3,X4)
                                                                 GOTO 90
   80   CONTINUE
        X(3,I)=SIGMA*X3
   90   CONTINUE
        ISTART=.NOT.ISTART
  100 CONTINUE
      RETURN
  110 CONTINUE
      RMAX=R+4.605*PDIF
      DO 140 I=1,N
  120   CONTINUE
        RAD=RMAX*(RNDM(V))**0.3333333
        CT=1.-2.*RNDM(V)
        FI=6.283185D0*RNDM(V)
        ST=SQRT(1.-CT*CT)
        X(1,I)=RAD*ST*COS(FI)
        X(2,I)=RAD*ST*SIN(FI)
        X(3,I)=RAD*CT
        RR=SQRT(X(1,I)**2+X(2,I)**2+X(3,I)**2)
        IF(JGLAUB.EQ.2)F=1./(1.+EXP((RR-R)/PDIF))
        IF(JGLAUB.EQ.1)THEN
	  IF(N.GE.11.OR.N.LE.18)THEN
	    RR0=R*R/(2.5D0-4.D0/AAN)
            F=(1.D0+(AAN-4.D0)*RR**2/(6.D0*RR0))*6.D0/(AAN-4.D0)
     *        *EXP(-RR**2/RR0+(AAN-10.D0)/(AAN-4.D0))
C           F=(1.D0+(AAN-4.D0)*RR**2/(6.D0*R0**2))
C    *        *EXP(-RR**2/R0**2)
C    PDIF=0.513
C           F=(1.D0-0.051D0*RR**2/R**2)/(1.+EXP((RR-R)/PDIF))
	  ELSEIF(N.GE.9.AND.N.LE.10)THEN
	    RR0=R*R/(2.5D0-4.D0/AAN)
            F=(1.D0+(AAN-4.D0)*RR**2/(6.D0*RR0))
     *        *EXP(-RR**2/RR0)
	  ELSE
            F=1./(1.+EXP((RR-R)/PDIF))
          ENDIF
        ENDIF
        IF (RNDM(V).GT.F)                                      GOTO 120
 
        IF (I.LT.2)                                            GOTO 140
 
        I1=I-1
        DO 130 I2=1,I1
          DIST2=(X(1,I)-X(1,I2))**2+(X(2,I)-X(2,I2))**2+(X(3,I)-X(3,I2))
     *    **2
          IF (DIST2.LE.R2MIN)                                   GOTO 120
 
  130   CONTINUE
  140 CONTINUE
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE MODB(BSITE,N,BSTEP,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
*KEEP,ZENTRA.
      COMMON /ZENTRA/ ICENTR
      COMMON /DSHM/ RASH,RBSH,BMAX,BSTAP,SIGSH,ROSH,GSH,
     *              BSITA(0:1,200),NSTATB,NSITEB
      COMMON /DSHMS/ SIGSHS
*KEND.
      DIMENSION BSITE(0:1,N)
      LOGICAL LEFT
C**** BSITE(1)=0. , BSITE(N)=1.
C
C                            Central collisions for lead-lead
C
C     WRITE(6,*)'RASH,RBSH,BMAX,BSTAP,SIGSH,ROSH,GSH,ICENTR,N',
C    *RASH,RBSH,BMAX,BSTAP,SIGSH,ROSH,GSH,ICENTR,N
      IF (ICENTR.EQ.1) THEN
        IF(RASH.EQ.RBSH)THEN
	IF(RASH.LE.15.D0)THEN
          BB=RNDM(V)*9.D0
          B=SQRT(BB)
          RETURN
        ELSEIF(RASH.GT.5.D0)THEN
          BB=RNDM(V)*30.D0
          B=SQRT(BB)
          RETURN
	ENDIF
        ELSEIF(RASH.LT.RBSH)THEN
          BB=RNDM(V)*(RBSH-RASH+3.D0)**2
          B=SQRT(BB)
          RETURN
        ELSEIF(RASH.GT.RBSH)THEN
          BB=RNDM(V)*(RASH-RBSH+3.D0)**2
          B=SQRT(BB)
          RETURN
        ENDIF
      ENDIF
C
      Y=RNDM(V)
      I0=1
      I2=N
   10 CONTINUE
      I1=(I0+I2)/2
      LEFT=((BSITE(1,I0)-Y)*(BSITE(1,I1)-Y)).LT.0.D0
      IF(LEFT)                                                  GO TO 20
      I0=I1
                                                                GO TO 30
   20 CONTINUE
      I2=I1
   30 CONTINUE
      IF(I2-I0-2)40,50,60
   40 CONTINUE
      I1=I2+1
      IF(I1.GT.N)I1=I0-1
                                                                GO TO 70
   50 CONTINUE
      I1=I0+1
                                                                GO TO 70
   60 CONTINUE
                                                                GO TO 10
   70 CONTINUE
      X0=(I0-1)*BSTEP
      X1=(I1-1)*BSTEP
      X2=(I2-1)*BSTEP
      Y0=BSITE(1,I0)
      Y1=BSITE(1,I1)
      Y2=BSITE(1,I2)
   80 CONTINUE
      B=X0*(Y-Y1)*(Y-Y2)/((Y0-Y1)*(Y0-Y2)+1.E-9)+ X1*(Y-Y0)*(Y-Y2)
     +/((Y1-Y0)*(Y1-Y2)+1.E-9)+ X2*(Y-Y0)*(Y-Y1)/((Y2-Y0)*(Y2-Y1)+1.E-9)
**sr 14.4.98
      B = B+0.5D0*BSTEP
      IF (B.LT.0.0D0) B = X1
      IF (B.GT.BMAX)  B = BMAX
**
 
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE DIAGR(NA,NB,B,JS,JT,INT,INTA,INTB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
*---
*   -sample new impct parameter B
*   -sample nucleon configuration for projectile and target
*   -sample number of collisions  INT, INTA, INTB
*---
*KEEP,INTMX.
      PARAMETER (INTMX=2488,INTMD=252)
*KEEP,SHMAKL.
C     INCLUDE (SHMAKL)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON/SHMAKL/JSSH(INTMX),JTSH(INTMX),INTER1(INTMX),INTER2(INTMX)
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEEP,DSHM.
      COMMON /DSHM/ RASH,RBSH,BMAX,BSTEP,SIGSH,ROSH,GSH,
     *              BSITE(0:1,200),NSTATB,NSITEB
      COMMON /DSHMS/ SIGSHS
*KEEP,NUCKOO.
      COMMON /NUCKOO/ PKOO(3,INTMX),TKOO(3,INTMX),PPOO(3,INTMX),
     +TPOO(3,INTMX)
*KEEP,DAMP.
C     COMPLEX*16 CA,CI
      DOUBLE COMPLEX CA,CI
      COMMON /DAMP/   CA,CI,GA
*KEND.
      COMMON /KKMATU/KKMATO,KKMATA,IPOO,IPOA,IPZOO,IPZOA
     *,IBPROO,IBPROA,IREADO
      COMMON /FLUCTU/IFLUCT
      COMMON /FLUARR/FLUSI(1000),FLUIX(1000),FLUIXX(1000)
      PARAMETER(NAMX=248)
      DIMENSION JS(NAMX),JT(NAMX)
C     COMPLEX*16 C
      DOUBLE COMPLEX C
      DATA ICNT/0/
      DATA INTCO/0/
C------------------------------
**sr 14.4.98
      CALL MODB(BSITE,NSITEB,BSTEP,B)
      INTCO=0
**
      DO 10 I=1,NA
   10 JS(I)=0
      DO 20 I=1,NB
   20 JT(I)=0
C--------
   30 INT=0
      INTA=0
      INTB=0
      INTCO=INTCO+1
      IF(INTCO.GE.500)THEN
        INTCO=0
        CALL MODB(BSITE,NSITEB,BSTEP,B)
      ENDIF
C--------
C     IF (KKMATO.NE.KKMATA.AND.IPOO.NE.IPOA.AND.
C    *  IBPROO.NE.IBPROA.AND.IPZOO.NE.IPZOA )THEN
        CALL CONUCL(PKOO,NA,RASH)
        CALL SORTIN(PKOO,NA)
        CALL CONUCL(TKOO,NB,RBSH)
        CALL SORT(TKOO,NB)
        KKMATO=KKMATA
	IPOO=IPOA
	IPZOO=IPZOA
	IBPROO=IBPROA
C     ELSEIF (KKMATO.EQ.KKMATA.AND.IPOO.EQ.IPOA.AND.
C    *                              IPZOO.EQ.IPZOA)THEN
C       IF(MOD(ICNT,5).EQ.0) THEN
C         CALL CONUCL(PKOO,NA,RASH)
C         CALL SORTIN(PKOO,NA)
C         CALL CONUCL(TKOO,NB,RBSH)
C         CALL SORT(TKOO,NB)
C       ENDIF
C       ICNT=ICNT+1
C     ENDIF
C
      IF(IPEV.GE.6) WRITE (6,1000)ICNT,PKOO(1,1),TKOO(1,1)
 1000 FORMAT (' 111 FORM IN DIAGR ICNT,PKOO(1,1),TKOO(1,1) ',I6,2F10.3)
C--------
      DO 40 I=1,NA
        X1=B-PKOO(1,I)
        X2=-PKOO(2,I)
        IF(IFLUCT.EQ.0)THEN
          AFLUK=1.
        ELSEIF(IFLUCT.EQ.1)THEN
          IFUK=(RNDM(V)+0.001)*1000.
          AFLUK=FLUIXX(IFUK)
        ENDIF
        DO 40 J=1,NB
          Q1=X1+TKOO(1,J )
          Q2=X2+TKOO(2,J )
          XY=GA*(Q1*Q1+Q2*Q2)
          IF(XY.GT.15.D0)                                    GO TO 40
          E=EXP(-XY)
          C=CI-CA*E*AFLUK
          AR=REAL(REAL(C))
          AI=IMAG(C)
          P=AR*AR+AI*AI
          IF(RNDM(V).LT.P)                                  GO TO 40
          INT=INT+1
          IF(INT.GT.INTMX)                                   GO TO 40
          JS(I)=JS(I)+1
          JT(J)=JT(J)+1
          INTER1(INT)=I
          INTER2(INT)=J
   40 CONTINUE
      DO 50 I=1,NA
        IF (JS(I).NE.0) INTA=INTA+1
   50 CONTINUE
      DO 60 J=1,NB
        IF (JT(J).NE.0) INTB=INTB+1
   60 CONTINUE
      IF(IPEV.GE.6) THEN
        WRITE(6,'(A)')
     +  ' DIAGR - AFTER 30 CONTINUE: ICNT, INT, B, NA,RA, NB,RB'
      WRITE(6,'(I10,I5,1PE11.3,2(I5,1PE11.3))') ICNT, INT, B, NA,RASH,
     +  NB,RBSH
      ENDIF
      IF(INT.EQ.0)                                              GO TO 30
      RETURN
      END
C----------------------------------------------------------------
       SUBROUTINE FLUINI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C      INITIALIZE THE ARRAY FOR CROSS SECTION FLUCTUATIONS
C
      COMMON /FLUARR/ FLUSI(1000),FLUIX(1000),FLUIXX(1000)
      DX=0.003D0
      B=0.893D0
      N=6
      A=0.1
      OM=1.1
      FLUSU=0
      FLUSUU=0
      DO 1 I=1,1000
        X=I*DX
        FLUIX(I)=X
        FLUS=((X-B)/(OM*B))**N
        IF(FLUS.LE.20.)THEN
          FLUSI(I)=(X/B)*EXP(-((X-B)/(OM*B))**N)/(X/B+A)
        ELSE
          FLUSI(I)=0.
        ENDIF
        FLUSU=FLUSU+FLUSI(I)
    1 CONTINUE
      DO 2 I=1,1000
        FLUSUU=FLUSUU+FLUSI(I)/FLUSU
        FLUSI(I)=FLUSUU
    2 CONTINUE
      WRITE(6,3)
    3 FORMAT(' FLUCTUATIONS')
      CALL PLOT(FLUIX,FLUSI,1000,1,1000,0.D0,0.06D0,0.D0,0.01D0)
      DO 5 I=1,1000
        AF=I*0.001D0
        DO 6 J=1,1000
          IF(AF.LE.FLUSI(J))THEN
            FLUIXX(I)=FLUIX(J)
            GO TO 7
          ENDIF
    6   CONTINUE
    7   CONTINUE
    5 CONTINUE
      FLUIXX(1)=FLUIX(1)
      FLUIXX(1000)=FLUIX(1000)
      RETURN
      END
*-- Author :
C
C--------------------------------------------------------------------
C
C     FILE TECALBAM
C
C
C--------------------------------------------------------------------
c     SUBROUTINE TECALB
C     SUBROUTINE TECALBAM
C
C------------------------------------------------------------------
C
C               DTNTCBI.FOR
C
C------------------------------------------------------------------
C
      SUBROUTINE CALBAM(NNCH,I1,I2,IFB11,IFB22,IFB33,IFB44,
     *                  AMCH,NOBAM,IHAD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C--------------------------------------------------------------------
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
C            = 10 q -- q -- q Jet Capella Kopeliovich 
C                or aq -- aq -- aq flavors IFB11,IFB22,IFB33
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
      DATA ISYMM/0/
*KEND.
C---------------------
C     IPCO=6
C-------------------------------------------------------------------------------C  
C                       Symmetrize JETSET and BAMJET at small chain masses
C                       for NOBAM=4 or 6
C
C-------------------------------------------------------------------------------
      IF(NOBAM.EQ.4.AND.ISYMM.EQ.1)THEN
        IFB4=IFB44
        IF (AMCH.LT.3.D0)THEN
          RR=RNDM(V)
          IF (RR.LT.0.33333D0)THEN
            IFB1=IFB11
            IFB2=IFB22
            IFB3=IFB33
          ELSEIF (RR.GT.0.666666D0)THEN
            IFB1=IFB22
            IFB2=IFB11
            IFB3=IFB33
          ELSE
            IFB1=IFB33
            IFB2=IFB22
            IFB3=IFB11
          ENDIF
        ELSEIF(AMCH.GT.7.D0)THEN
          IFB1=IFB11
          IFB2=IFB22
          IFB3=IFB33
          IFB4=IFB44
        ELSE
          SSSS=(7.D0-AMCH)/4.D0
          RRR=RNDM(VV)
          IF(RRR.LT.1.D0-SSSS)THEN
            IFB1=IFB11
            IFB2=IFB22
            IFB3=IFB33
          ELSE
            RR=RNDM(V)
            IF (RR.LT.0.33333D0)THEN
              IFB1=IFB11
              IFB2=IFB22
              IFB3=IFB33
            ELSEIF (RR.GT.0.666666D0)THEN
              IFB1=IFB22
              IFB2=IFB11
              IFB3=IFB33
            ELSE
              IFB1=IFB33
              IFB2=IFB22
              IFB3=IFB11
            ENDIF
          ENDIF
        ENDIF
      ELSEIF(NOBAM.EQ.6.AND.ISYMM.EQ.1)THEN
        IFB4=IFB44
        IF (AMCH.LT.3.D0)THEN
          RR=RNDM(V)
          IF (RR.LT.0.33333D0)THEN
            IFB1=IFB11
            IFB2=IFB22
            IFB3=IFB33
          ELSEIF (RR.GT.0.666666D0)THEN
            IFB3=IFB22
            IFB1=IFB11
            IFB2=IFB33
          ELSE
            IFB3=IFB11
            IFB2=IFB22
            IFB1=IFB33
          ENDIF
        ELSEIF(AMCH.GT.7.D0)THEN
          IFB1=IFB11
          IFB2=IFB22
          IFB3=IFB33
          IFB4=IFB44
        ELSE
          SSSS=(7.D0-AMCH)/4.D0
          RRR=RNDM(VV)
          IF(RRR.LT.1.D0-SSSS)THEN
            IFB1=IFB11
            IFB2=IFB22
            IFB3=IFB33
          ELSE
            RR=RNDM(V)
            IF (RR.LT.0.33333D0)THEN
              IFB1=IFB11
              IFB2=IFB22
              IFB3=IFB33
            ELSEIF (RR.GT.0.666666D0)THEN
              IFB3=IFB22
              IFB1=IFB11
              IFB2=IFB33
            ELSE
              IFB1=IFB33
              IFB2=IFB22
              IFB3=IFB11
            ENDIF
          ENDIF
        ENDIF
      ELSE
        IFB1=IFB11
        IFB2=IFB22
        IFB3=IFB33
        IFB4=IFB44
      ENDIF
C------------------------------------------------------------------------------
      IF(IPCO.GE.6)THEN
        WRITE (6,1000)NNCH,I1,I2,IFB1,IFB2,IFB3,IFB4,AMCH,NOBAM,IHAD
 1000 FORMAT(' CALBAM:NNCH,I1,I2,IFB1,IFB2,IFB3,IFB4,AMCH,NOBAM,IHAD' /7
     +I5,F10.3,2I5)
      ENDIF
      IF(NOBAM.EQ.10)THEN
        CALL DBAMJE(IHAD,IFB1,IFB2,IFB3,IFB4,AMCH,NOBAM)
      ENDIF
      ITURN=0
      IF (NNCH) 10,30,20
   10 CONTINUE
      IF(NNCH.EQ.-99) THEN
        IHAD=0
        RETURN
      ENDIF
      ICH1=I1
                                                                GO TO 50
   20 CONTINUE
      IF(NNCH.EQ.99) THEN
        IHAD=0
        RETURN
      ENDIF
      ICH1=I2
                                                                GO TO 50
   30 CONTINUE
C***  ITURN=0                       HJM 24/01/91
      IF (IFB1.LE.6)                                            GO TO 40
      ITURN=1
      IF(NOBAM.EQ.3) CALL DBAMJE(IHAD,IFB2,IFB1,IFB3,IFB4,AMCH,NOBAM)
      IF(NOBAM.EQ.4)THEN
        ITURN=0
        CALL DBAMJE(IHAD,IFB1,IFB2,IFB3,IFB4,AMCH,NOBAM)
      END IF
      IF(NOBAM.EQ.6) CALL DBAMJE(IHAD,IFB3,IFB1,IFB2,IFB4,AMCH,4)
      IF(NOBAM.EQ.5)THEN
        CALL DBAMJE(IHAD,IFB3,IFB4,IFB1,IFB2,AMCH,NOBAM)
      ENDIF
                                                                GO TO 60
   40 CONTINUE
 
      IF (NOBAM.EQ.3.OR.NOBAM.EQ.4.OR.NOBAM.EQ.5) THEN
        CALL DBAMJE(IHAD,IFB1,IFB2,IFB3,IFB4,AMCH,NOBAM)
      ELSEIF(NOBAM.EQ.6)THEN
        ITURN=1
        CALL DBAMJE(IHAD,IFB3,IFB1,IFB2,IFB4,AMCH,4)
      END IF
                                                                GO TO 60
   50 CONTINUE
      CALL DPARJE(IHAD,ICH1)
   60 CONTINUE
C     CALL DECAY(IHAD,2)
   70 CONTINUE
C*** TURN CHAINS AROUND IF NESSESARY
      IF (ITURN.EQ.0)                                          GO TO 100
C*** TURN JET AROUND
      DO 80 I=1,IHAD
        PZF(I)=-PZF(I)
   80 CONTINUE
   90 CONTINUE
  100 CONTINUE
      IF (IPCO.GE.6)THEN
        DO 1244 I=1,IHAD
          WRITE(6,1245)I,PZF(I),HEF(I),ANF(I)
 1245     FORMAT(I5,2F10.4,A8)
 1244   CONTINUE
      ENDIF
      RETURN
      END
*-- Author :
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
      SUBROUTINE DBAMJE(IHAD,KFA1,KFA2,KFA3,KFA4,AE0,IOPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C*****IHAD=NUMBER OF FINAL HADRONS AND HADRON RESONANCES
C*****AE0=INITIAL ENERGY IN GEV
C*****KFA=INITIAL QUARK FLAVOUR
C*****KFD1,KFD2=FLAVOUR CONTENTS OF A INITIAL DIQUARK
C*****IOPT=1,2,3,4 MEANS: SINGLE QUARK JET,SINGLE DIQUARK JET,
C*****COMPLETE QUARK ANTIQUARK TWOJET EVENT,COMPLETE QUARK-DIQUARK TWO
C*****JET EVENT
C*****IOPT=10 q -- q -- q Jet Capella Kopeliovich (only jetset defined)
C*****COMMON/FINPAR/ CONTAINS THE MOMENTA,ENERGIES AND QUANTUM NUMBERS
C*****OF THE CREATED HADRONS
C*****IV IS THE ACTUAL VERTEX,IV=1,4,5,6,9,10 ARE MESON VERTIZES
C*****IV=2,3,7,8 ARE BARYON VERTIZES
C*****LA=1 MEANS CUT-OFF
C*****LL=0,1 MEANS QUARK JET,ANTIQUARK JET(DIQUARK JET,ANTIDIQUARK JET)
C*****COMMON/REMAIN/ CONTAINS REST JET ENERGY,MOMENTA AND QUANTUMNUMBERS
C------------------
*KEEP,DFINPA.
      CHARACTER*8 ANF
      PARAMETER (NFIMAX=249)
      COMMON /DFINPA/ ANF(NFIMAX),PXF(NFIMAX),PYF(NFIMAX),PZF(NFIMAX),
     +HEF(NFIMAX),AMF(NFIMAX), ICHF(NFIMAX),IBARF(NFIMAX),NREF(NFIMAX)
      COMMON /DFINPZ/IORMO(NFIMAX),IDAUG1(NFIMAX),IDAUG2(NFIMAX),
     * ISTATH(NFIMAX)
*KEEP,DINPDA.
      COMMON /DINPDA/ IMPS(6,6),IMVE(6,6),IB08(6,21),IB10(6,21), IA08
     +(6,21),IA10(6,21), A1,B1,B2,B3,LT,LE,BET,AS,B8,AME,DIQ,ISU
*KEND.
      CHARACTER*8 ANAME
      COMMON /DPAR/ ANAME(210),AM(210),GA(210),TAU(210),ICH(210), IBAR
     +(210),K1(210),K2(210)
      COMMON/DREMAI/ RPXR,RPYR,RPZR,RER,KR1R,KR2R
      COMMON /BAMCO/  NVDD
      COMMON/IFRAGM/IFRAG
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*--------------------------- S.Roesler 5/27/93
      COMMON /DIFFRA/ ISINGD,IDIFTP,IOUDIF,IFLAGD
*
      DIMENSION RPX(100),RPY(100),RE(100)
      DIMENSION KFR1(100),KFR2(100),IV(100)
      PARAMETER (PIMASS=0.15D0)
C
      IF(IPCO.GE.6)LT=1
C-------------------------------------------------------------------
      IF (LT.EQ.1)WRITE(6, 1000)IHAD,KFA1,KFA2,KFA3,KFA4,AE0,IOPT
 1000 FORMAT (5I5,E12.4,I5,
     +      ' BAMJET,IHAD,KFA1,KFA2,KFA3,KFA4,AE0,IOPT')
C------------------------------------------------------------------
C
C                            JETSET-7.3 FRAGMENTATION j.r.6/93
C
      IF(IOPT.EQ.10)THEN
        IREJ=0
        CALL BAMLUN(IHAD,KFA1,KFA2,KFA3,KFA4,AE0,IOPT,IREJ)
        IF(IREJ.EQ.1)THEN
          RETURN
        ENDIF
        RETURN
      ENDIF
      IF(IFRAG.EQ.1.OR.IFRAG.EQ.2.OR.IFRAG.GE.10)THEN
        IREJ=0
        CALL BAMLUN(IHAD,KFA1,KFA2,KFA3,KFA4,AE0,IOPT,IREJ)
        IF(IREJ.EQ.1)THEN
          RETURN
        ENDIF
        RETURN
      ENDIF
C_________________________________________________________________
      AS=0.5
      B8=0.4
      A1=0.88
*-------------------------- S.Roesler 5/26/93
      B1=3.
      B2=3.
*
      B3=8.0
C     BET=9.5
      BET=12.0
      IF(NVDD.EQ.15) THEN
        A1=0.99
        B3=2.0
        BET=3.
      ENDIF
*-------------------------- S.Roesler 5/26/93
*                           diffractive chains
*
      IF(IFLAGD.EQ.1)THEN
         A1=0.88
         B1=6.
         B2=9.
         B3=4.0
      ENDIF
*
C
      ITRY = 0
   10 CONTINUE
C
C  avoid hang ups
      ITRY = ITRY+1
      IF(ITRY.GT.10000) THEN
        WRITE(6,'(/1X,A)') 'DBAMJE:ERROR: FRAGMENTATION IMPOSSIBLE'
        WRITE(6,'(1X,A,5I5,E12.3,I5)')
     &    'DBAMJE:IHAD,KFA1,KFA2,KFA3,KFA4,AE0,IOPT ',
     &    IHAD,KFA1,KFA2,KFA3,KFA4,AE0,IOPT
        STOP
      ENDIF
C
      DO 20 I=1,100
        KFR1(I)=0
        KFR2(I)=0
   20 CONTINUE
   30 CONTINUE
      IYY=0
      IHAD=0
      IT=0
      E0=AE0/2.
C  low mass asymmetric fragmentation (r.e. 12/93)
      IF(ITRY.GT.100) THEN
        IF(IOPT.EQ.3) THEN
          IF(KFA1.GT.(KFA2-6)) THEN
            E0=AE0-MAX(AE0*0.1,PIMASS)
          ELSE
            E0=MAX(AE0*0.1,PIMASS)
          ENDIF
        ENDIF
      ENDIF
C
      IF(IOPT.EQ.1.OR.IOPT.EQ.2) E0=AE0
      LL=0
      PGX=0.
      IF(KFA1.GT.6.AND.IOPT.EQ.1) LL=1
      IF(KFA1.LE.6.AND.IOPT.EQ.2) LL=1
      IF(KFA1.GT.6.AND.IOPT.EQ.4) LL=1
      PGY=0.
      PGZ=0.
      RPX0=0.
      RPY0=0.
      DO 130 I=1,100
        LA=0
        IT=IT+1
        J=IT-1
   40   CONTINUE
C*****CUT OFF TASK
C       CALL DABBRC(IT,LL,LA,LT,E0,PGX,PGY,PGZ,KFR1,KFR2,RE, KR1R,KR2R,
C    +  KR1L,KR2L,RPX,RPY,RPZ,RPXR,RPYR,RPZR,RPXL,RPYL,RPZL, RER,REL,IV,
C    +  B1,B2,KFA1,KFA2,KFA3,KFA4,IOPT,IYY)
        IF(LA.EQ.0)                                               GOTO60
        IT=IT-1
        IF(IOPT.EQ.3.AND.LL.EQ.0)                                GOTO 50
        IF(IOPT.EQ.4.AND.KFA1.LE.6.AND.LL.EQ.0)                  GOTO 50
        IF(IOPT.EQ.4.AND.KFA1.GT.6.AND.LL.EQ.1)                  GOTO 50
        IF(IOPT.EQ.5.AND.LL.EQ.0)                                GOTO 50
C  asymmetric fragmentation (r.e. 12/93)
        E0 = AE0/2.
                                                                GOTO 140
   50   CONTINUE
C  asymmetric fragmentation (r.e. 12/93)
        E0 = AE0-E0
C
        IYY=1
        LL=1
        IF(IOPT.EQ.4.AND.KFA1.GT.6) LL=0
        IAR=IT
                                                                 GOTO120
   60   CONTINUE
C*****CHOICE OF THE VERTEX
C       CALL DVERTE(IT,LT,LL,KFA1,E0,IV,AME,IOPT)
C*****CHOICE OF THE FLAVOUR
C       CALL DFLAVO(IT,LT,LL,E0,IV,RE,KFR1,KFR2,ISU,BET,KFA1,KFA2,KFA3,
C    +  KFA4,IOPT)
C*****CLASSIFICATION OF THEPARTICLES
C       CALL DHKLAS(IT,LT,LA,LL,KFR1,KFR2,KR1R,KR2R,KR1L,KR2L,IV,IMPS,
C    +  IMVE,IB08,IA08,IB10,IA10,AS,B8,KFA1,KFA2,KFA3,KFA4,IOPT)
        IF(IT.EQ.1)RX=E0
        IF(IT.GT.1)RX=RE(J)
        IF(AMF(IT).GT.RX)                                        GOTO 10
        IF(AMF(IT).LE.RX)                                        GOTO 70
        LA=1
                                                                  GOTO40
   70   CONTINUE
        IHAD=IHAD+1
        IF(LT.EQ.0)                                               GOTO80
        WRITE(6, 1070)IHAD
   80   CONTINUE
C*****CHOICE OF THE ENERGY
C       CALL DENERG(IT,IV,RE,HMA,HE,E0,A1)
C*****CHOICE OF THE MOMENTUM
C       CALL IMPULD(HE,HMA,HPS,HPX,HPY,HPZ,LT,LL,B3)
        IF(IT.GT.1)                                               GOTO90
        RPX(IT)=RPX0-HPX
        RPY(IT)=RPY0-HPY
                                                                 GOTO100
   90   RPX(IT)=RPX(J)-HPX
        RPY(IT)=RPY(J)-HPY
  100   CONTINUE
        IF (IOPT.EQ.1.AND.LL.EQ.1)HPZ=-HPZ
        IF(IOPT.EQ.2.AND.LL.EQ.1) HPZ=-HPZ
        IF(IOPT.EQ.4.AND.KFA1.GT.6) HPZ=-HPZ
        IF(IOPT.EQ.5) HPZ=-HPZ
        PGX=PGX+HPX
        PGY=PGY+HPY
        PGZ=PGZ+HPZ
        PXF(IT)=HPX
        PYF(IT)=HPY
        PZF(IT)=HPZ
        IF(LT.EQ.0)                                              GOTO110
        WRITE(6, 1010)PGX,PGY,PGZ
 1010 FORMAT(1H0,12HPGX,PGY,PGZ=,3F8.4)
  110   CONTINUE
  120   CONTINUE
  130 CONTINUE
  140 CONTINUE
      IF(IOPT.EQ.1.OR.IOPT.EQ.2)                                GOTO 150
C*****PUT THE RIGHT AND LEFT JET TOGETHER
C     CALL DVEREI(IT,LA,LT,RER,REL,RPXR,RPYR,RPZR,RPXL,RPYL,RPZL,
C    *KR1R,KR2R,KR1L,KR2L,IHAD,LL,KFR1,KFR2,IMPS,IMVE,IB08,IA08,
C    *IB10,IA10,B3,AS,B8,IAR,KFA1,KFA2,KFA3,KFA4,IOPT)
      IF(LA.EQ.3)                                                GOTO 10
      IF(LA.EQ.2)                                                GOTO 10
  150 CONTINUE
      IF(IOPT.EQ.3.OR.IOPT.EQ.4.OR.IOPT.EQ.5)                   GOTO 160
      IF(LL.EQ.0)                                               GOTO 160
      RPXR=RPXL
      RPYR=RPYL
      RPZR=RPZL
      RER=REL
      KR1R=KR1L
      KR2R=KR2L
  160 CONTINUE
      IF(LE.EQ.0)                                               GOTO 180
      WRITE(6, 1030)
      DO170 J=1,IT
        WRITE(6, 1020)NREF(J),ANF(J),AMF(J),ICHF(J),
     +  IBARF(J),PXF(J),PYF(J),
     +  PZF(J),HEF(J)
 1020 FORMAT(2X,I3,A6,F6.3,2I4,4F8.4)
 1030 FORMAT(2X,'NF,NAME,MASS,IQ,IB,PX,PY,PZ,E')
  170 CONTINUE
  180 CONTINUE
 1040 FORMAT(1H0,38HNUMBER OF EVENTS WITH PREST GT. EREST=,I4, /,
     +21HNUMBER OF ALL EVENTS=,I4)
 1050 FORMAT(1H0,'NUMBER OF EVENTS WITH ONLY ONE PARTICLE=',I4)
C*****TEST OF THE CONSERVATION LAWS
C      CALL TERHAL(IT,LE,KFA1,KFA2,IOPT)
      IF(LT.EQ.0)                                                GOTO190
      WRITE(6, 1060)IHAD
 1060 FORMAT(1H0,' MULTIPLIZITAET=',I3)
 1070 FORMAT(1H0,13HHADRONANZAHL=,I3)
  190 CONTINUE
C     IF (IHAD.EQ.2)THEN
C       DO 200 I=1,IHAD
C         PZF(I)=-PZF(I)
C 200   CONTINUE
C     ENDIF
C
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE INDEXD(KA,KB,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      KP=KA*KB
      KS=KA+KB
      IF(KP.EQ.1)IND=1
      IF(KP.EQ.2)IND=2
      IF(KP.EQ.3)IND=3
      IF(KP.EQ.4.AND.KS.EQ.5)IND=4
      IF(KP.EQ.5)IND=5
      IF(KP.EQ.6.AND.KS.EQ.7)IND=6
      IF(KP.EQ.4.AND.KS.EQ.4)IND=7
      IF(KP.EQ.6.AND.KS.EQ.5)IND=8
      IF(KP.EQ.8)IND=9
      IF(KP.EQ.10)IND=10
      IF(KP.EQ.12.AND.KS.EQ.8)IND=11
      IF(KP.EQ.9)IND=12
      IF(KP.EQ.12.AND.KS.EQ.7)IND=13
      IF(KP.EQ.15)IND=14
      IF(KP.EQ.18)IND=15
      IF(KP.EQ.16)IND=16
      IF(KP.EQ.20)IND=17
      IF(KP.EQ.24)IND=18
      IF(KP.EQ.25)IND=19
      IF(KP.EQ.30)IND=20
      IF(KP.EQ.36)IND=21
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      DOUBLE PRECISION FUNCTION DBETA(X1,X2,BET)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      AX=0.0
      BETX1=BET*X1
      IF(BETX1.LT.70.) AX=-1./BET**2*(BETX1+1.)*EXP(-BETX1)
      AY=1./BET**2*(BET*X2+1.)*EXP(-BET*X2)
      DBETA=AX+AY
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE DDRELA(X,Y,Z,COTE,SITE,COPS,SIPS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      X1=COPS*X-SIPS*COTE*Y+SIPS*SITE*Z
      X2=SIPS*X+COPS*COTE*Y-COPS*SITE*Z
      X3=SITE*Y+COTE*Z
      X=X1
      Y=X2
      Z=X3
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C      SUBROUTINE DPOLI(CS,SI)
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     SAVE
C      U=RNDM(V)
C     CS=RNDM(VV)
C     IF (U.LT.0.5) CS=-CS
C     SI=SQRT(1.-CS*CS+1.E-10)
C     RETURN
C     END
C     SUBROUTINE ALTRAF(GA,BGA,CX,CY,CZ,COD,COF,SIF,PC,EC,P,PX,PY,PZ,E)
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     SAVE
C     BGX=BGA*CX
C     BGY=BGA*CY
C     BGZ=BGA*CZ
C     COD2=COD*COD
C     IF (COD2.GT.0.999999) COD2=0.999999
C     SID=SQRT(1.-COD2)*PC
C     PCX=SID*COF
C     PCY=SID*SIF
C     PCZ=COD*PC
C     EP=PCX*BGX+PCY*BGY+PCZ*BGZ
C     PE=EP/(GA+1.)+EC
C     PX=PCX+BGX*PE
C     PY=PCY+BGY*PE
C     PZ=PCZ+BGZ*PE
C     P=SQRT(PX*PX+PY*PY+PZ*PZ)
C     PM=1./P
C     PX=PX*PM
C     PY=PY*PM
C     PZ=PZ*PM
C     E=GA*EC+EP
C     RETURN
C     END
C     SUBROUTINE ROTAT(PX,PY,PZ,PXN,PYN,PZN,COTE,SITE,COPS,SIPS)
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     SAVE
C     PXN=-PX*SIPS-PY*COTE*COPS+PZ*SITE*COPS
C     PYN=PX*COPS-PY*COTE*SIPS+PZ*SITE*SIPS
C     PZN=PY*SITE+PZ*COTE
C     RETURN
C     END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
*=== threpd ===========================================================*
*
      SUBROUTINE DTHREP(UMO,ECM1,ECM2,ECM3,PCM1,PCM2,PCM3,COD1,COF1,
     &SIF1,COD2,COF2,SIF2,COD3,COF3,SIF3,AM1,AM2,AM3)
 
*$ CREATE DBLPRC.ADD
*COPY DBLPRC
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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
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
*----------------------------------------------------------------------*
*  Threpd89: slight revision by A. Ferrari                             *
*  Last change on   11-oct-93   by    Alfredo Ferrari, INFN - Milan    *
*----------------------------------------------------------------------*
*
      DIMENSION F(5),XX(5)
C***THREE PARTICLE DECAY IN THE CM - SYSTEM
      COMMON /DGAMRE/ REDU,AMO,AMM(15 )
      COMMON/DDREI/UUMO,AAM1,AAM2,AAM3,S22,UMO2,
     *AM11,AM22,AM33,S2SUP,S2SAP(2)
C     COMMON/PRUNT/ISYS
      COMMON/PRITT/ISYS
      SAVE EPS
      DATA EPS/AZRZRZ/
*
      UMOO=UMO+UMO
C***S1, S2, S3 ARE THE INVARIANT MASSES OF THE PARTICLES 1, 2, 3
C***J. VON NEUMANN - RANDOM - SELECTION OF S2
C***CALCULATION OF THE MAXIMUM OF THE S2 - DISTRIBUTION
      UUMO=UMO
      AAM1=AM1
      AAM2=AM2
      AAM3=AM3
      GU=(AM2+AM3)**2
      GO=(UMO-AM1)**2
*     UFAK=1.0000000000001D0
*     IF (GU.GT.GO) UFAK=0.9999999999999D0
      IF (GU.GT.GO) THEN
         UFAK=ONEMNS
      ELSE
         UFAK=ONEPLS
      END IF
      OFAK=2.D0-UFAK
      GU=GU*UFAK
      GO=GO*OFAK
      DS2=(GO-GU)/99.D0
      AM11=AM1*AM1
      AM22=AM2*AM2
      AM33=AM3*AM3
      UMO2=UMO*UMO
      RHO2=0.D0
      S22=GU
      DO 124 I=1,100
         S21=S22
         S22=GU+(I-1.D0)*DS2
         RHO1=RHO2
         RHO2=DXLAMB(S22,UMO2,AM11)*DXLAMB(S22,AM22,AM33)/
     *                                             (S22+EPS)
         IF(RHO2.LT.RHO1) GO TO 125
  124 CONTINUE
  125 S2SUP=(S22-S21)*.5D0+S21
      SUPRHO=DXLAMB(S2SUP,UMO2,AM11)*DXLAMB(S2SUP,AM22,AM33)/
     *                                           (S2SUP+EPS)
      SUPRHO=SUPRHO*1.05D0
      XO=S21-DS2
      IF (GU.LT.GO.AND.XO.LT.GU) XO=GU
      IF (GU.GT.GO.AND.XO.GT.GU) XO=GU
      XX(1)=XO
      XX(3)=S22
      X1=(XO+S22)*0.5D0
      XX(2)=X1
      F(3)=RHO2
      F(1)=DXLAMB(XO,UMO2,AM11)*DXLAMB(XO,AM22,AM33)/(XO+EPS)
      F(2)=DXLAMB(X1,UMO2,AM11)*DXLAMB(X1,AM22,AM33)/(X1+EPS)
      DO 126 I=1,16
         X4=(XX(1)+XX(2))*0.5D0
         X5=(XX(2)+XX(3))*0.5D0
         F(4)=DXLAMB(X4,UMO2,AM11)*DXLAMB(X4,AM22,AM33)/
     *                                               (X4+EPS)
         F(5)=DXLAMB(X5,UMO2,AM11)*DXLAMB(X5,AM22,AM33)/
     *                                               (X5+EPS)
         XX(4)=X4
         XX(5)=X5
         DO 128 II=1,5
            IA=II
            DO 128 III=IA,5
               IF (F (II).GE.F (III)) GO TO 128
               FH=F(II)
               F(II)=F(III)
               F(III)=FH
               FH=XX(II)
               XX(II)=XX(III)
               XX(III)=FH
128      CONTINUE
         SUPRHO=F(1)
         S2SUP=XX(1)
         DO 129 II=1,3
            IA=II
            DO 129 III=IA,3
               IF (XX(II).GE.XX(III)) GO TO 129
               FH=F(II)
               F(II)=F(III)
               F(III)=FH
               FH=XX(II)
               XX(II)=XX(III)
               XX(III)=FH
129      CONTINUE
126   CONTINUE
      AM23=(AM2+AM3)**2
      ITH=0
      REDU=2.D0
    1 CONTINUE
      ITH=ITH+1
      IF (ITH.GT.200) REDU=-9.D0
      IF (ITH.GT.200) GO TO 400
      C=RNDM(C)
*     S2=AM23+C*((UMO-AM1)**2-AM23)
      S2=AM23+C*(UMO-AM1-AM2-AM3)*(UMO-AM1+AM2+AM3)
      Y=RNDM(Y)
      Y=Y*SUPRHO
      RHO=DXLAMB(S2,UMO2,AM11)*DXLAMB(S2,AM22,AM33)/S2
      IF(Y.GT.RHO) GO TO 1
C***RANDOM SELECTION OF S3 AND CALCULATION OF S1
      S1=RNDM(S1)
      S1=S1*RHO+AM11+AM22-(S2-UMO2+AM11)*(S2+AM22-AM33)/(2.D0*S2)-
     &RHO*.5D0
      S3=UMO2+AM11+AM22+AM33-S1-S2
      ECM1=(UMO2+AM11-S2)/UMOO
      ECM2=(UMO2+AM22-S3)/UMOO
      ECM3=(UMO2+AM33-S1)/UMOO
      PCM1=SQRT((ECM1+AM1)*(ECM1-AM1))
      PCM2=SQRT((ECM2+AM2)*(ECM2-AM2))
      PCM3=SQRT((ECM3+AM3)*(ECM3-AM3))
      CALL DSFECF(SFE,CFE)
C***TH IS THE ANGLE BETWEEN PARTICLES 1 AND 2
C***TH1, TH2 ARE THE ANGLES BETWEEN PARTICLES 1, 2 AND THE DIRECTION OF
      PCM12 = PCM1 * PCM2
      IF ( PCM12 .LT. ANGLSQ ) GO TO 200
      COSTH=(ECM1*ECM2+0.5D+00*(AM11+AM22-S1))/PCM12
      GO TO 300
 200  CONTINUE
         UW=RNDM(UW)
         COSTH=(UW-0.5D+00)*2.D+00
 300  CONTINUE
*     IF(ABS(COSTH).GT.0.9999999999999999D0)
*    &COSTH=SIGN(0.9999999999999999D0,COSTH)
      IF(ABS(COSTH).GT.ONEONE)
     &COSTH=SIGN(ONEONE,COSTH)
      IF (REDU.LT.1.D+00) RETURN
      COSTH2=(PCM3*PCM3+PCM2*PCM2-PCM1*PCM1)/(2.D+00*PCM2*PCM3)
*     IF(ABS(COSTH2).GT.0.9999999999999999D0)
*    &COSTH2=SIGN(0.9999999999999999D0,COSTH2)
      IF(ABS(COSTH2).GT.ONEONE)
     &COSTH2=SIGN(ONEONE,COSTH2)
      SINTH2=SQRT((ONEONE-COSTH2)*(ONEONE+COSTH2))
      SINTH =SQRT((ONEONE-COSTH)*(ONEONE+COSTH))
      SINTH1=COSTH2*SINTH-COSTH*SINTH2
      COSTH1=COSTH*COSTH2+SINTH2*SINTH
C***RANDOM SELECTION OF THE SPHERICAL COORDINATES OF THE DIRECTION OF PA
C***CFE, SFE ARE COS AND SIN OF THE ROTATION ANGLE OF THE SYSTEM 1, 2 AR
C***THE DIRECTION OF PARTICLE 3
C***CALCULATION OF THE SPHERICAL COORDINATES OF PARTICLES 1, 2
      CX11=-COSTH1
      CY11=SINTH1*CFE
      CZ11=SINTH1*SFE
      CX22=-COSTH2
      CY22=-SINTH2*CFE
      CZ22=-SINTH2*SFE
      CALL DSFECF(SIF3,COF3)
      COD3=TWOTWO*RNDM(COD3)-ONEONE
      SID3=SQRT((1.D+00-COD3)*(1.D+00+COD3))
    2 FORMAT(5F20.15)
      COD1=CX11*COD3+CZ11*SID3
      CHLP=(ONEONE-COD1)*(ONEONE+COD1)
      IF(CHLP.LT.1.D-14)WRITE(ISYS,2)COD1,COF3,SID3,
     &CX11,CZ11
      SID1=SQRT(CHLP)
      COF1=(CX11*SID3*COF3-CY11*SIF3-CZ11*COD3*COF3)/SID1
      SIF1=(CX11*SID3*SIF3+CY11*COF3-CZ11*COD3*SIF3)/SID1
      COD2=CX22*COD3+CZ22*SID3
      SID2=SQRT((ONEONE-COD2)*(ONEONE+COD2))
      COF2=(CX22*SID3*COF3-CY22*SIF3-CZ22*COD3*COF3)/SID2
      SIF2=(CX22*SID3*SIF3+CY22*COF3-CZ22*COD3*SIF3)/SID2
 400  CONTINUE
* === Energy conservation check: === *
      EOCHCK = UMO - ECM1 - ECM2 - ECM3
*     SID1 = SQRT ( ( ONEONE - COD1 ) * ( ONEONE + COD1 ) )
*     SID2 = SQRT ( ( ONEONE - COD2 ) * ( ONEONE + COD2 ) )
*     SID3 = SQRT ( ( ONEONE - COD3 ) * ( ONEONE + COD3 ) )
      PZCHCK = PCM1 * COD1 + PCM2 * COD2 + PCM3 * COD3
      PXCHCK = PCM1 * COF1 * SID1 + PCM2 * COF2 * SID2
     &       + PCM3 * COF3 * SID3
      PYCHCK = PCM1 * SIF1 * SID1 + PCM2 * SIF2 * SID2
     &       + PCM3 * SIF3 * SID3
      EOCMPR = 1.D-12 * UMO
      IF ( ABS (EOCHCK) + ABS (PXCHCK) + ABS (PYCHCK) + ABS (PZCHCK)
     &     .GT. EOCMPR ) THEN
         WRITE(LUNERR,*)
     &   ' *** Threpd: energy/momentum conservation failure! ***',
     &   EOCHCK,PXCHCK,PYCHCK,PZCHCK
         WRITE(LUNERR,*)' *** SID1,SID2,SID3',SID1,SID2,SID3
      END IF
      RETURN
      END
*=== xlamb ============================================================*
*
      DOUBLE PRECISION FUNCTION DXLAMB(X,Y,Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      DATA IDGB/0/
      COMMON /DGAMRE/ REDU,AMO,AMM(15 )
      COMMON/DDREI/TEST(12)
      YZ=Y-Z
      DXLAMB=X*X-2.D0*X*(Y+Z)+YZ*YZ
      XLAM =DXLAMB
      IF (IDGB.LE.0) GO TO 11
      IF(DXLAMB.GT.1.D-12) GOTO 11
      WRITE(6,12)
      WRITE(6,10) XLAM,X,Y,Z,TEST
      WRITE(6,13)
   12 FORMAT(/,10X,' DXLAMB PRINT')
   13 FORMAT(10X,60(1H*))
   10 FORMAT(4E20.8,'DXLAMB',/,12F10.5)
   11 CONTINUE
      IF(DXLAMB.LE.0.D0)DXLAMB=ABS(DXLAMB)
      DXLAMB=SQRT(DXLAMB)
      RETURN
      END
 
*-- Author :
C     DOUBLE PRECISION FUNCTION DXLAMB(X,Y,Z)
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     SAVE
C     YZ=Y-Z
C     DXLAMB=X*X-2.*X*(Y+Z)+YZ*YZ
C     IF(DXLAMB.LE.0.)DXLAMB=ABS(DXLAMB)
C     DXLAMB=SQRT(DXLAMB)
C     RETURN
C     END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE STRAFO(GAM,BGAM,CX,CY,CZ,COD,COF,SIF,P,ECM,
     1PL,CXL,CYL,CZL,EL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C     LORENTZ TRANSFORMATION INTO THE LAB - SYSTEM
      SID=SQRT((1.-COD)*(1.+COD)+1.E-22)
      SIF=SQRT((1.-COF)*(1.+COF)+1.E-22)
      PLX=P*SID*COF
      PLY=P*SID*SIF
      PCMZ=P*COD
      PLZ=GAM*PCMZ+BGAM*ECM
      PL=SQRT(PLX*PLX+PLY*PLY+PLZ*PLZ)
      EL=GAM*ECM+BGAM*PCMZ
C     ROTATION INTO THE ORIGINAL DIRECTION
      COZ=PLZ/PL
      IF(COZ.GE.1.)COZ=0.999999999999
      SIZ=SQRT((1.-COZ)*(1.+COZ))
      CALL DRTRAN(CX,CY,CZ,COZ,SIZ,SIF,COF,CXL,CYL,CZL)
      RETURN
      END
*-- Author :
      SUBROUTINE DRTRAN(XO,YO,ZO,CDE,SDE,SFE,CFE,X,Y,Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      IF (ABS(XO)-0.0001) 10,10,30
   10 IF (ABS(YO)-0.0001) 20,20,30
   20 CONTINUE
      X=SDE*CFE
      Y=SDE*SFE
      Z=CDE*ZO
      RETURN
   30 CONTINUE
      XI=SDE*CFE
      YI=SDE*SFE
      ZI=CDE
      A=SQRT(XO**2+YO**2)
      X=-YO*XI/A-ZO*XO*YI/A+XO*ZI
      Y=XO*XI/A-ZO*YO*YI/A+YO*ZI
      Z=A*YI+ZO*ZI
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE DCHANT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      CHARACTER*8 ZKNAME
C      COMMON/DDECAC/ ZKNAME(540),NZK(540,3),WT(540)
 
      PARAMETER (IDMAX9=602)
C      CHARACTER*8 ZKNAME
      COMMON/DDECAC/ ZKNAME(IDMAX9),WT(IDMAX9),NZK(IDMAX9,3)
 
 
      CHARACTER*8 ANAME
      COMMON/DPAR/ANAME(210),AM(210),GA(210),TAU(210),ICH(210),IBAR(210)
     *,K1(210),K2(210)
      DIMENSION HWT(602)
C  CHANGE OF WEIGHTS WT FROM ABSOLUT VALUES INTO THE SUM OF WT OF A DEC.
      DO 10 J=1,602
   10 HWT(J)=0.
      DO 30 I=1,210
        IK1=K1(I)
        IK2=K2(I)
        HV=0.
        DO 20 J=IK1,IK2
          HV=HV+WT(J)
          HWT(J)=HV
C         IF(HWT(J).GT.1.) WRITE(6,1000) HWT(J),J,I,IK1
   20   CONTINUE
C1000 FORMAT(2X,15H ERROR IN HWT =,1F10.5,8H J,I,K1=,3I5)
   30 CONTINUE
      DO 40 J=1,602
   40 WT(J)=HWT(J)
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE DTWOPD(UMO,ECM1,ECM2,PCM1,PCM2,COD1,COF1,SIF1, COD2,
     +COF2,SIF2,AM1,AM2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C*****TWO PARTICLE DECAY IN THE CM - SYSTEM
C
      IF(UMO.LT.(AM1+AM2)) THEN
        WRITE(6,'(/,A/A,3(1PE12.4))')
     +  ' INCONSISTENT CALL OF TWOPAD / EXECUTION STOPPED',
     +  ' UMO, AM1, AM2 :', UMO, AM1, AM2
        STOP
      ENDIF
C
      ECM1=((UMO-AM2)*(UMO+AM2) + AM1*AM1)/(2.*UMO)
      ECM2=UMO-ECM1
      PCM1=SQRT((ECM1-AM1)*(ECM1+AM1))
      PCM2=PCM1
      CALL DSFECF(SIF1,COF1)
      COD1=2.*RNDM(X)-1.
      COD2=-COD1
      COF2=-COF1
      SIF2=-SIF1
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE DSFECF(SFE,CFE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
   10 X=RNDM(V)
      Y=RNDM(V)
      XX=X*X
      YY=Y*Y
      XY=XX+YY
      IF(XY.GT.1)                                                 GOTO10
      CFE=(XX-YY)/XY
      SFE=2.*X*Y/XY
      IF(RNDM(V).LT.0.5D0)                                      GOTO20
      RETURN
   20 SFE=-SFE
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE DDATES
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      CHARACTER*8 ZKNAME,Z
C      COMMON /DDECAC/ ZKNAME(540),NZK(540,3),WT(540)
 
      PARAMETER (IDMAX9=602)
C      CHARACTER*8 ZKNAME
      COMMON/DDECAC/ ZKNAME(IDMAX9),WT(IDMAX9),NZK(IDMAX9,3)
 
 
      CHARACTER*8 ANAME
      COMMON /DPAR/ ANAME(210),AM(210),GA(210),TAU(210), ICH(210),IBAR
     +(210),K1(210),K2(210)
C----------------------
      DIMENSION ICHAR(210)
      EQUIVALENCE (ICH(1),ICHAR(1))
      DIMENSION Z(3)
C---------------------------
      WRITE(6, 1000)
 1000 FORMAT(1H1,'   ')
      WRITE(6, 1010)
 1010 FORMAT(///' TABLE OF USED PARTICLES AND RESONANCES (I)',//
     +' I = NUMBER OF PARTICLE OR RESONANCE',/
     +' IPDG = P D G NUMBER OF PARTICLE OR RESONANCE',/
     +' ANAME = NAME OF I'/, ' AM = MASS OF I  (GEV)',/
     +' GA = WIDTH OF I  (GEV)',/ ' TAU = LIFE TIME OF I  (SEC.)',/
     +' ICH = ELECTRIC CHARGE OF I, IBAR = BARYONIC CHARGE OF I',/' ', '
     +K1 = FIRST DECAY CHANNEL NUMBER, K2 = LAST DECAY CHANNEL NUMBER OF
     +I')
 
 
      WRITE(6, 1020)
 1020 FORMAT(///
     +'   I  ANAME   AM         GA         TAU       ICH IBAR K1  K2'/)
      JOO=210
      DO 10 I=1,JOO
        IPDG=MPDGHA(I)
        WRITE(6, 1030)I,IPDG,ANAME(I),AM(I),
     +  GA(I),TAU(I),ICH(I),IBAR(I), K1
     +  (I),K2(I)
 1030 FORMAT (1I4,I6,2X,1A8,3E11.4,4I4)
        IF(I.EQ.43) WRITE(6, 1000)
        IF(I.EQ.43) WRITE(6, 1020)
        IF(I.EQ.99) WRITE(6, 1000)
        IF(I.EQ.99) WRITE(6, 1020)
        IF(I.EQ.155) WRITE(6, 1000)
        IF(I.EQ.155) WRITE(6, 1020)
   10 CONTINUE
      WRITE(6, 1000)
      WRITE(6, 1040)
 1040 FORMAT(///' DECAY CHANNELS OF PARTICLES AND RESONANCES',//)
      WRITE(6, 1050)
 1050 FORMAT(' ANAME = PARTICLE AND RESONANCE NAME'/,
     +' DNAME = DECAY CHANNEL NAME'/, ' J = DECAY CHANNEL NUMBER'/,
     +' I = NUMBER OF DECAYING PARTICLE'/,
     +' WT = SUM OF DECAY CHANNEL WEIGHTS FROM K1(I) UP TO J'/,
     +' NZK = PROGRAM INTERNAL NUMBERS OF DECAY PRODUCTS')
 
      WRITE(6, 1060)
 1060 FORMAT(///'   I     J      ANAME       DNAME                DECAY
     +PRODUCTS            WT       NZK'/)
      DO 60 I=1,JOO
        IK1=K1(I)
        IK2=K2(I)
        IF (IK1.LE.0)                                           GO TO 60
        DO 50 IK=IK1,IK2
          I1=NZK(IK,1)
          I2=NZK(IK,2)
          I3=NZK(IK,3)
          IF (I1.LE.0) I1=29
          IF (I2.LE.0) I2=29
          IF (I3.LE.0) I3=29
          J1=I1
          J2=I2
          J3=I3
          Z(1)=ANAME(I1)
          Z(2)=ANAME(I2)
          Z(3)=ANAME(I3)
          WRITE(6, 1070)I,IK,ANAME(I),ZKNAME(IK),(Z(J),J=1,3),WT(IK),J1,J2,
     +    J3
 1070 FORMAT(2I5,' DECAY OF ',1A8,' (CHANNEL: ',1A6,' ) TO ',3(1A6,2X),
     +1F8.4,3I5)
          AMTEST=AM(I)-AM(J1)-AM(J2)-AM(J3)
          IBTEST=IBAR(I)-IBAR(J1)-IBAR(J2)-IBAR(J3)
          ICTEST=ICHAR(I)-ICHAR(J1)-ICHAR(J2)-ICHAR(J3)
          IF (AMTEST) 20,30,30
   20     MTEST=1
                                                                GO TO 40
   30     MTEST=0
   40     CONTINUE
          IF (MTEST+IBTEST**2+ICTEST**2.NE.0) WRITE(6, 1080)AMTEST,
     +    IBTEST,
     +    ICTEST
 1080 FORMAT (' ***** ERROR IN MASS, BAR.CH. OR E.CH. ',F10.5,2I6)
          IF(IK.EQ.27) WRITE(6, 1000)
          IF(IK.EQ.27) WRITE(6, 1060)
          IF(IK.EQ.62) WRITE(6, 1000)
          IF(IK.EQ.62) WRITE(6, 1060)
          IF(IK.EQ.101) WRITE(6, 1000)
          IF(IK.EQ.101) WRITE(6, 1060)
          IF(IK.EQ.144) WRITE(6, 1000)
          IF(IK.EQ.144) WRITE(6, 1060)
          IF(IK.EQ.183) WRITE(6, 1000)
          IF(IK.EQ.183) WRITE(6, 1060)
          IF(IK.EQ.222) WRITE(6, 1000)
          IF(IK.EQ.222) WRITE(6, 1060)
          IF(IK.EQ.261) WRITE(6, 1000)
          IF(IK.EQ.261) WRITE(6, 1060)
          IF(IK.EQ.300) WRITE(6, 1000)
          IF(IK.EQ.300) WRITE(6, 1060)
          IF(IK.EQ.362) WRITE(6, 1000)
          IF(IK.EQ.362) WRITE(6, 1060)
          IF(IK.EQ.401) WRITE(6, 1000)
          IF(IK.EQ.401) WRITE(6, 1060)
          IF(IK.EQ.440) WRITE(6, 1000)
          IF(IK.EQ.440) WRITE(6, 1060)
          IF(IK.EQ.479) WRITE(6, 1000)
          IF(IK.EQ.479) WRITE(6, 1060)
          IF(IK.EQ.518) WRITE(6, 1000)
          IF(IK.EQ.518) WRITE(6, 1060)
   50   CONTINUE
   60 CONTINUE
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE DDATAR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      COMMON/DINPDA/IMPS(6,6),IMVE(6,6),IB08(6,21),IB10(6,21),
     *IA08(6,21),IA10(6,21),A1,B1,B2,B3,LT,LB,BET,AS,B8,AME,DIQ,ISU
      DIMENSION IV(36),IP(36),IB(126),IBB(126),IA(126),IAA(126)
C     DEFINE THE FIELDS FOR PARTICLE CLASSIFICATION
C     IMPS=PSEUDO SCALAR MESONS (SPIN=0)
C     IMVE=VECTOR MESONS (SPIN=1)
C     IB08(IA08)=BARYONS (ANTIBARYONS) (SPIN=1/2)
C     IB10(IA10)=BARYONS (ANTIBARYONS) (SPIN=3/2)
      DATA IP/
     *23,14,16,116,0,0,13,23,25,117,0,0,15,24,31,120,0,0,119,118,121,
     *122,14*0/
      L=0
      DO 20 I=1,6
        DO 10 J=1,6
          L=L+1
          IMPS(I,J)=IP(L)
   10   CONTINUE
   20 CONTINUE
      DATA IV/
     *33,34,38,123,0,0,32,33,39,124,0,0,36,37,96,127,0,0,126,125,128,
     *129,14*0/
      L=0
      DO 40 I=1,6
        DO 30 J=1,6
          L=L+1
          IMVE(I,J)=IV(L)
   30   CONTINUE
   40 CONTINUE
      DATA IB/
     *0,1,21,140,0,0,8,22,137,0,0,97,138,0,0,146,5*0,
     *1,8,22,137,0,0,0,20,142,0,0,98,139,0,0,147,5*0,
     *21,22,97,138,0,0,20,98,139,0,0,0,145,0,0,148,5*0,
     *140,137,138,146,0,0,142,139,147,0,0,145,148,50*0/
      L=0
      DO 60 I=1,6
        DO 50 J=1,21
          L=L+1
          IB08(I,J)=IB(L)
   50   CONTINUE
   60 CONTINUE
      DATA IBB/
     *53,54,104,161,0,0,55,105,162,0,0,107,164,0,0,167,5*0,
     *54,55,105,162,0,0,56,106,163,0,0,108,165,0,0,168,5*0,
     *104,105,107,164,0,0,106,108,165,0,0,109,166,0,0,169,5*0,
     *161,162,164,167,0,0,163,165,168,0,0,166,169,0,0,170,47*0/
      L=0
      DO 80 I=1,6
        DO 70 J=1,21
          L=L+1
          IB10(I,J)=IBB(L)
   70   CONTINUE
   80 CONTINUE
      DATA IA/
     *0,2,99,152,0,0,9,100,149,0,0,102,150,0,0,158,5*0,
     *2,9,100,149,0,0,0,101,154,0,0,103,151,0,0,159,5*0,
     *99,100,102,150,0,0,101,103,151,0,0,0,157,0,0,160,5*0,
     *152,149,150,158,0,0,154,151,159,0,0,157,160,50*0/
      L=0
      DO 100 I=1,6
        DO 90 J=1,21
          L=L+1
          IA08(I,J)=IA(L)
   90   CONTINUE
  100 CONTINUE
      DATA IAA/
     *67,68,110,171,0,0,69,111,172,0,0,113,174,0,0,177,5*0,
     *68,69,111,172,0,0,70,112,173,0,0,114,175,0,0,178,5*0,
     *110,111,113,174,0,0,112,114,175,0,0,115,176,0,0,179,5*0,
     *171,172,174,177,0,0,173,175,178,0,0,176,179,0,0,180,47*0/
      L=0
      DO 120 I=1,6
        DO 110 J=1,21
          L=L+1
          IA10(I,J)=IAA(L)
  110   CONTINUE
  120 CONTINUE
C     DEFINE THE FREE PARAMETERS FOR THE MONTE-CARLO PROGRAMMES BAMJET
C     AND PARJET
      A1=0.88
      B3=8.0
      B1=8.0
      B2=8.0
      ISU=4
c     BET=8.0
      BET=9.5
      BET=12.
      AS=0.50
      AME=0.95
      B8=0.40
      DIQ=0.375
      LT=0
      LB=0
C     THE FOLLOWING ARE THE PARAMETERS USED IN ABR8611
      A1=0.95
      A1=0.50
      A1=0.88
      B3=8.
      B1=6.
      B1=3.
      B2=6.
      B2=3.
      AS=0.25
*     AME=0.95
      B8=0.33
C     WRITE (6,123)A1,B3,B1,B2,ISU,BET,AS,AME,LT,LB,B8,DIQ
C 123 FORMAT (' DATAR3 INITIALIZATION:PARAMETERS SET LIKE IN ABR8611'/
C    *' A1    = ',F10.4/
C    *' B3    = ',F10.4/
C    *' B1    = ',F10.4/
C    *' B2    = ',F10.4/
C    *' ISU   = ',I10/
C    *' BET   = ',F10.4/
C    *' AS    = ',F10.4/
C    *' AME   = ',F10.4/
C    *' LT    = ',I10/
C    *' LB    = ',I10/
C    *' B8    = ',F10.4/
C    *' DIQ   = ',F10.4)
      RETURN
      END
*-- Author :
C
C
C--------------------------------------------------------------------
C
C     FILE TECALBAM
C
C
C--------------------------------------------------------------------
      SUBROUTINE TECALB
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C     SUBROUTINE TECALBAM
C
C     TWO CHAIN FRAGMENTATION MODEL FOR PARTICLE PRODUCTION
C     TEST OF  CALBAM ROUTINE CALLING BAMJET
C                 JUNE   1987, J.RANFT
C********************************************************************
C
C      OPTIONS
C
C      JNI=7 BAMJET DISTRIBUTIONS ANALOG TO JNI=2
C
C      IP=1=P,  IP=2=AP,  IP=8=N,  IP=9=AN,  IP=13=PI+,  IP=14=PI-,
C      IP=15=K+,  IP=16=K-,
C
C**********************************************************************
C
*KEEP,DFINPA.
      CHARACTER*8 ANF
      PARAMETER (NFIMAX=249)
      COMMON /DFINPA/ ANF(NFIMAX),PXF(NFIMAX),PYF(NFIMAX),PZF(NFIMAX),
     +HEF(NFIMAX),AMF(NFIMAX), ICHF(NFIMAX),IBARF(NFIMAX),NREF(NFIMAX)
      COMMON /DFINPZ/IORMO(NFIMAX),IDAUG1(NFIMAX),IDAUG2(NFIMAX),
     * ISTATH(NFIMAX)
*KEEP,DINPDA.
      COMMON /DINPDA/ IMPS(6,6),IMVE(6,6),IB08(6,21),IB10(6,21), IA08
     +(6,21),IA10(6,21), A1,B1,B2,B3,LT,LE,BET,AS,B8,AME,DIQ,ISU
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEND.
      COMMON/JNI/JNI
      COMMON /DIFF/ IDIFF
      COMMON /DKPL/UPLO,IPQ
      COMMON /DINV/PNUC(3),INUCVT
      COMMON/CAPKOP/XX1,XX3
C----------------------------------------------------------------------
      WRITE (6,1000)
 1000 FORMAT (' ##############################################'/
     +' PROGRAM TECABAPT'/
     +' ######################################################')
C     CALL PRIBLO
C     CALL DATAR3
C     CALL CHANWT
      INIT=0
      LT=0
      IT=1
      IP=1
      IPRI=0
   10 CONTINUE
      READ(5,1010)JNI,NEVT,IP,NCASES,POO,AOO,ZNUC
 1010 FORMAT(4I10,3F10.2)
      IF (IT.EQ.0) IT=1
      IF (JNI.LE.0)                                            GO TO 120
      WRITE(6, 1010)JNI,NEVT,IP,NCASES,POO,AOO,ZNUC,IT
      IF (JNI.LT.0)STOP
C**********   JNI SELECTS OPTION   *************************************
      JNI=7
      GO TO (20,30,40,50,60,70,80,100),JNI
   20 CONTINUE
   30 CONTINUE
   40 CONTINUE
   50 CONTINUE
   60 CONTINUE
   70 CONTINUE
                                                               GO TO 110
   80 CONTINUE
C*** JNI=7 CALCUL. BAMJET EVENTS
C*** WITH ENERGY POO     (CMS ENERGY)
      IPRI=0
      LT=0
      INIT=0
      IPQ=1
      XX1=0.9
      XX3=1.
      UPLO=POO
      III=1
      CALL DISTCM(1,IPQ,POO,IPQ,IPQ)
C     CALL DISRES(1,IPQ,POO,IPQ,IPQ)
      DO 90 L=1, NEVT
        IF (IP.EQ.103)CALL CALBAM(0,1,1,1,7,1,1,POO,3,NHAD)
        IF (IP.EQ.109)CALL CALBAM(0,1,1,7,1,1,1,POO,3,NHAD)
        IF (IP.EQ.104)CALL CALBAM(0,1,1,1,2,2,1,POO,4,NHAD)
        IF (IP.EQ.1010)CALL CALBAM(0,1,1,1,2,3,1,POO,4,NHAD)
        IF (IP.EQ.105)CALL CALBAM(0,1,1,1,2,7,8,POO,5,NHAD)
        IF (IP.EQ.1011)CALL CALBAM(0,1,1,7,7,1,1,POO,5,NHAD)
        IF (IP.EQ.106)CALL CALBAM(0,1,1,1,1,1,1,POO,6,NHAD)
        IF (IP.EQ.1012)CALL CALBAM(0,1,1,7,7,7,1,POO,6,NHAD)
        IF (IP.EQ.1050)CALL CALBAM(0,1,1,1,1,2,1,POO,10,NHAD)
C
        CALL DDECAY(NHAD,2)
        CALL DISTCM(2,NHAD,POO,IPQ,NCASES)
   90 CONTINUE
C
      WRITE(6, 1020)POO,IP,NCASES
 1020 FORMAT (' BAMJET (POO,IP,NCASES) = ',1F10.2,2I10)
      CALL DISTCM(3,NEVT,POO,IPQ,NCASES)
C     CALL DISRES(3,III*NEVT,POO,IPQ,NCASES)
                                                               GO TO 110
  100 CONTINUE
  110 CONTINUE
                                                                GO TO 10
  120 CONTINUE
      RETURN
      END
*-- Author :
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE DISTCM(IOP,NHAD,POLAB,KPROJ,KTARG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C*** 1=P,  2=N, 3=PI+, 4=PI-, 5=PIO, 6=GAM+HYP, 7=K, 8=ANUC, 9=CHARGED
C*** 10=TOT, 11=TOTHAD
C
      CHARACTER*8 ANH
      PARAMETER (NFIMAX=249)
      COMMON /DFINPA/ ANH(NFIMAX),PX(NFIMAX),PY(NFIMAX),PZ(NFIMAX),
     +HE(NFIMAX),AM(NFIMAX), ICH(NFIMAX),IBAR(NFIMAX),NR(NFIMAX)
      COMMON /DFINPZ/IORMO(NFIMAX),IDAUG1(NFIMAX),IDAUG2(NFIMAX),
     * ISTATH(NFIMAX)
C---------------------
      CHARACTER*8 ANAME
      COMMON /DPAR/ ANAME(210),AAM(210),GA(210),TAU(210),IICH(210),
     +IIBAR(210),K1(210),K2(210)
C------------------
      COMMON /HISTO / XMULT(100,10),YMULT(100,10),XXFL(50,20), YXFL
     +(50,20),XYL(50,20),YYL(50,20), YYLPS(50,20),PTP(50,20),PTY(50,20),
     +FILL(6000)
      DIMENSION AVMULT(12,30),AVE(12,30),INDX(25),MU(12,30), AKNO
     +(100,2),XKNO(100,2),AKE(12,30),AASO(12,30)
      COMMON /DKPL/UPLO,KPL
      COMMON/JNI/JNI
C
      DATA IPRIOP/1/
      DATA INDX/1,8,10,10,10,10,7,2,7,10,10,7,3,4,5,6,
     *11,12,7,13,14,15,16,17,18/
C-----------------------------------------------------------------------
      GO TO (10,60,100),IOP
   10 CONTINUE
      KPL=1
      AVPT=0.
      NAVPT=0
      DXFL=0.04
      PO=POLAB
C     EEO=SQRT(PO**2+AAM(NHAD)**2)
C     UMO=SQRT(AAM(KTARG)**2+AAM(KPROJ)**2+2.*AAM(KTARG)*EEO)
      IF(JNI.EQ.7) UMO=POLAB
      UMO=POLAB
      EEO=UMO
      PO=EEO/2.D0
      WRITE(6, 1000)EEO,PO,NHAD
 1000 FORMAT (' EEO',F10.2,F10.2,I10)
      DY=0.2
      DPT=0.10
      DO 20 I=1,10
        AVMULT(KPL,I)=1.E-18
        AVE(KPL,I)=0.
        DO 20 J=1,100
          XMULT(J,I)=J-1
          YMULT(J,I)=1.D-18
          AKNO(J,1)=1.D-18
          AKNO(J,2)=1.D-18
   20 CONTINUE
      WRITE(6, 1000)EEO,PO,NHAD
      DO 30 I=1,20
        DO 30 J=1,50
          XXFL(J,I)=J*DXFL -1.
          YXFL(J,I)=1.D-18
          XYL(J,I)=-5.0+J*DY
          YYL(J,I)=1.D-18
          YYLPS(J,I)=1.E-18
          PTP(J,I)=J*DPT
          PTY(J,I)=1.D-18
   30 CONTINUE
      WRITE(6, 1000)EEO,PO,NHAD
   40 CONTINUE
      DO 50 I=1,30
        AVE(KPL,I)=1.D-18
        AVMULT(KPL,I)=1.D-18
        MU(KPL,I)=0
        AASO(KPL,I)=UPLO
   50 CONTINUE
      WRITE(6, 1000)EEO,PO,NHAD
      RETURN
C
   60 CONTINUE
C      WRITE(6, 1000)EEO,PO,NHAD
      AVMULT(KPL,30)=AVMULT(KPL,30)+NHAD
      NNHAD=NHAD+1
      IF (NNHAD.GT.100) NNHAD=100
      YMULT(NNHAD,10)=YMULT(NNHAD,10)+1.
      DO 70 I=1,30
        MU(KPL,I)=0
   70 CONTINUE
      EETOT=0.D0
      DO 71 I=1,NHAD
	IF(IBAR(I).NE.500)THEN
C         WRITE(6,*)I,ANH(I),HE(I),AM(I),IBAR(I),NR(I)
	  EETOT=EETOT+HE(I)
        ENDIF
   71 CONTINUE
      IF(EETOT.GT.POLAB+1.D-6)THEN
	WRITE(6,*)' eetot.gt.polab ',EETOT,POLAB
      ENDIF
C     WRITE(6,*)' eetot ',EETOT
      DO 80 I=1,NHAD
	IF(IBAR(I).NE.500)THEN
        NRE=NR(I)
        IF (NRE.GT.25) NRE=3
        IF (NRE.LT. 1) NRE=3
        NI=INDX(NRE)
        IF (NRE.EQ.28)NI=8
        AVE(KPL,NRE)=AVE(KPL,NRE)+HE(I)
        AVE(KPL,30)=AVE(KPL,30)+HE(I)
        IF (NI.NE.6) AVE(KPL,29)=AVE(KPL,29)+HE(I)
        AVMULT(KPL,NRE)=AVMULT(KPL,NRE)+1.
        IF (NI.NE.6) AVMULT(KPL,29)=AVMULT(KPL,29)+1.
        MU(KPL,NI)=MU(KPL,NI)+1
        IF (ICH(I).NE.0)MU(KPL,9)=MU(KPL,9)+1
C    TOTAL=30   TOTAL-GAMMA=29   ANTIHYP=28
C    CHARGED=27
        IF (ICH(I).NE.0)AVE(KPL,27)=AVE(KPL,27)+HE(I)
        IF (ICH(I).NE.0)AVMULT(KPL,27)=AVMULT(KPL,27)+1
C       XFL=PZ(I)/PO
        XFL=(PZ(I)/ABS(PZ(I)))*HE(I)/PO
        IXFL=XFL/DXFL+26.
        IF (IXFL.LT. 1) IXFL=1
        IF (IXFL.GT.50) IXFL=50
C       XXXFL=SQRT(XFL**2+(AM(I)+0.3)**2/PO**2)
        XXXFL=ABS(XFL)
        IF (ICH(I).NE.0)YXFL(IXFL,9)=YXFL(IXFL,9)+XXXFL
        YXFL(IXFL,NI)=YXFL(IXFL,NI)+XXXFL
        YXFL(IXFL,10)=YXFL(IXFL,10)+XXXFL
        PTT=PX(I)**2+PY(I)**2
        YL=0.5*LOG(ABS((HE(I)+PZ(I)+1.E-10)/(HE(I)-PZ(I)+1.E-10)))
        YLPS=LOG(ABS((PZ(I)+SQRT(PZ(I)**2+PTT))/SQRT(PTT)+1.E-18))
        IYLPS=(YLPS+5.0)/DY
        IF (IYLPS.LT.1)IYLPS=1
        IF (IYLPS.GT.50)IYLPS=50
        YYLPS(IYLPS,NI)=YYLPS(IYLPS,NI)+1.
        YYLPS(IYLPS,10)=YYLPS(IYLPS,10)+1.
        IF (ICH(I).NE.0)YYLPS(IYLPS,9)=YYLPS(IYLPS,9)+1.
        IYL=(YL+5.0)/DY
        IF (IYL.LT.1) IYL=1
        IF (IYL.GT.50) IYL=50
        IF (ICH(I).NE.0)YYL(IYL,9)=YYL(IYL,9)+1.
        YYL(IYL,NI)=YYL(IYL,NI)+1.
        YYL(IYL,10)=YYL(IYL,10)+1.
        PT=SQRT(PTT)+0.001
        AVPT=AVPT+PT
        NAVPT=NAVPT+1
        IPT=PT/DPT+1.
        IF (IPT.LT.1)IPT=1
        IF (IPT.GT.50) IPT=50
        IF (ICH(I).NE.0)PTY(IPT,9)=PTY(IPT,9)+1./PT
        PTY(IPT,NI)=PTY(IPT,NI)+1./PT
        PTY(IPT,10)=PTY(IPT,10)+1./PT
        ENDIF
   80 CONTINUE
      DO 90 I=1,9
        IM=MU(KPL,I)+1
        IF (IM.GT.100)IM=100
        YMULT(IM,I)=YMULT(IM,I)+1.
   90 CONTINUE
      RETURN
C------------------------------------------------
  100 CONTINUE
      WRITE(6, 1000)EEO,PO,NHAD
C1020 FORMAT (' AVMULT=',11F10.5/,' AVE=',11F10.5)
      DO 110 I=1,30
        AVMULT(KPL,I)=AVMULT(KPL,I)/NHAD
        AVE(KPL,I)=AVE(KPL,I)/NHAD
  110 CONTINUE
      AVPT=AVPT/NAVPT
      WRITE (6,1030)AVPT,NAVPT
 1030 FORMAT (' AVERAGE PT= ',F12.4,I10)
      WRITE(6, 1040)
 1040 FORMAT(' PARTICLE REF,CHAR,IBAR, MASS      AVERAGE',
     +' ENERGY, MULTIPLICITY, INELASTICITY')
      DO 120 I=1,30
        AKE(KPL,I)=AVE(KPL,I)/EEO
        WRITE(6, 1050)ANAME(I),I,IICH(I),IIBAR(I),
     +   AAM(I), AVE(KPL,I),AVMULT
     +  (KPL,I),AKE(KPL,I)
 1050 FORMAT (' ',A8,3I5,F10.3,3F18.6)
  120 CONTINUE
      DO 130 I=1,10
        DO 130 J=1,100
          YMULT(J,I)=YMULT(J,I)/NHAD
  130 CONTINUE
      DO 140 I=1,20
        DO 140 J=1,50
          YXFL(J,I)=YXFL(J,I)/(NHAD*DXFL)
          YY L(J,I)=YY L(J,I)/(NHAD*DY)
          YYLPS(J,I)=YYLPS(J,I)/(NHAD*DY)
          PTY(J,I)=PTY(J,I)/(NHAD*DPT)
  140 CONTINUE
  150 CONTINUE
      WRITE(6, 1060)
 1060 FORMAT('1 RAPIDITY DISTRIBUTION')
      DO 160 J=1,50
        WRITE(6, 1070)XYL(J,1),(YYL(J,I),I=1,10)
 1070 FORMAT (F10.2,10E11.3)
  160 CONTINUE
      DO 161 J=1,50
        WRITE(6, 1070)XYL(J,1),(YYL(J,I),I=11,20)
  161 CONTINUE
      WRITE(6, 1060)
      WRITE(6,*)' 1=*=p, 2=n, 3=pi+, 4=pi-, 5=K+, 6=K-,',
     &' 8=ap, 9=chargd., 10=Z=all, 11=+=Lamb, 12=A=aLam,',
     &' 13=O=sig-, 14=B=Sig+, 15=C=Sig0, 16=D=pi0'
      CALL PLOT(XYL,YYL,1000,20,50,-5.D0,DY,0.D0,0.1D0)
      WRITE(6,*)' 1=*=p, 2=n, 3=pi+, 4=pi-, 5=K+, 6=K-,',
     &' 8=ap, 9=chargd., 10=Z=all, 11=+=Lamb, 12=A=aLam,',
     &' 13=O=sig-, 14=B=Sig+, 15=C=Sig0, 16=D=pi0'
      CALL PLOT(XYL,YYLPS,1000,20,50,-5.D0,DY,0.D0,0.1D0)
C     IF (IPRIOP.EQ.1) GO TO 1423
      WRITE(6, 1080)
 1080 FORMAT ('1  LONG MOMENTUM (SCALED) DISTRIBUTION')
      DO 170 J=1,50
        WRITE(6, 1070)XXFL(J,1),(YXFL(J,I),I=1,10)
  170 CONTINUE
      DO 171 J=1,50
        WRITE(6, 1070)XXFL(J,1),(YXFL(J,I),I=11,20)
  171 CONTINUE
  180 CONTINUE
      WRITE(6, 1080)
      WRITE(6,*)' 1=*=p, 2=n, 3=pi+, 4=pi-, 5=K+, 6=K-,',
     &' 8=ap, 9=chargd., 10=Z=all, 11=+=Lamb, 12=A=aLam,',
     &' 13=O=sig-, 14=B=Sig+, 15=C=Sig0, 16=D=pi0'
      CALL PLOT(XXFL,YXFL,1000,20,50,-1.D0,DXFL,0.D0,0.05D0)
      WRITE(6, 1090)
 1090 FORMAT ('1 MULTIPLICITY DISTRIBUTIONS')
      SIMUL=0.
      SUMUL=0.
      DO 190 J=1,100
        SUMUL=SUMUL+YMULT(J,10)
        SIMUL=SIMUL+YMULT(J,9)
  190 CONTINUE
      WRITE(6, 1100)(XMULT(J,1),YMULT(J,9),YMULT(J,10),J=1,100)
 1100 FORMAT(F6.1,2E12.4,F6.1,2E12.4,F6.1,2E12.4,F6.1,2E12.4)
      WRITE(6, 1090)
      WRITE(6,*)' 1=*=p, 2=n, 3=pi+, 4=pi-, 5=K+, 6=K-,',
     &' 8=ap, 9=chargd., 10=Z=all, 11=+=Lamb, 12=A=aLam,',
     &' 13=O=sig-, 14=B=Sig+, 15=C=Sig0, 16=D=pi0'
      CALL PLOT(XMULT,YMULT,1000,10,100,0.D0,1.D0,0.D0,0.01D0)
      DO 200 I=1,100
        XKNO(I,1)=I/AVMULT(KPL,30)
        XKNO(I,2)=I/AVMULT(KPL,27)
        AKNO(I,1)=YMULT(I,10)*AVMULT(KPL,30)/SUMUL
        AKNO(I,2)=YMULT(I,9)*AVMULT(KPL,27)/SIMUL
        AKNO(I,1)=LOG10(AKNO(I,1)+1.D-9)
        AKNO(I,2)=LOG10(AKNO(I,2)+1.D-9)
  200 CONTINUE
      WRITE(6, 1110)
 1110 FORMAT ('1 KNO MULTIPLICITY DISTRIBUTIONS')
      WRITE(6,*)' 1=*=p, 2=n, 3=pi+, 4=pi-, 5=K+, 6=K-,',
     &' 8=ap, 9=chargd., 10=Z=all, 11=+=Lamb, 12=A=aLam,',
     &' 13=O=sig-, 14=B=Sig+, 15=C=Sig0, 16=D=pi0'
      CALL PLOT(XKNO,AKNO,200,2,100,0.D0,0.08D0,-4.D0,0.05D0)
      DO 210 I=1,10
        DO 210 J=1,100
          YMULT(J,I)=LOG10(YMULT(J,I))
  210 CONTINUE
      DO 220 I=1,20
        DO 220 J=1,50
          YXFL(J,I)=LOG10(ABS(YXFL(J,I)+1.D-8))
          YYL(J,I)=LOG10(YYL(J,I)+1.D-8)
          PTY(J,I)=LOG10(PTY(J,I)+1.D-8)
  220 CONTINUE
  230 CONTINUE
      WRITE(6, 1060)
      WRITE(6,*)' 1=*=p, 2=n, 3=pi+, 4=pi-, 5=K+, 6=K-,',
     &' 8=ap, 9=chargd., 10=Z=all, 11=+=Lamb, 12=A=aLam,',
     &' 13=O=sig-, 14=B=Sig+, 15=C=Sig0, 16=D=pi0'
      CALL PLOT(XYL,YYL,1000,20,50,-5.D0,DY,-3.5D0,0.05D0)
      DO 240 J=1,50
        WRITE(6, 1070)XXFL(J,1),(YXFL(J,I),I=1,10)
  240 CONTINUE
      DO 241 J=1,50
        WRITE(6, 1070)XXFL(J,1),(YXFL(J,I),I=11,20)
  241 CONTINUE
      WRITE(6, 1080)
      WRITE(6,*)' 1=*=p, 2=n, 3=pi+, 4=pi-, 5=K+, 6=K-,',
     &' 8=ap, 9=chargd., 10=Z=all, 11=+=Lamb, 12=A=aLam,',
     &' 13=O=sig-, 14=B=Sig+, 15=C=Sig0, 16=D=pi0'
      CALL PLOT(XXFL,YXFL,1000,20,50,-1.D0,DXFL,-4.5D0,0.05D0)
      WRITE(6,1120)
 1120 FORMAT ('1 PT DISTRIBUTION DN/PTDPT')
      CALL PLOT(PTP,PTY,1000,20,50,0.D0,DPT,-2.0D0,0.05D0)
      IF (IPRIOP.EQ.1)                                         GO TO 250
      WRITE(6,1090)
      WRITE(6,*)' 1=*=p, 2=n, 3=pi+, 4=pi-, 5=K+, 6=K-,',
     &' 8=ap, 9=chargd., 10=Z=all, 11=+=Lamb, 12=A=aLam,',
     &' 13=O=sig-, 14=B=Sig+, 15=C=Sig0, 16=D=pi0'
      CALL PLOT(XMULT,YMULT,1000,10,100,0.D0,1.D0, -3.5D0,0.05D0)
  250 CONTINUE
      IF (KPL.NE.12)                                           GO TO 270
      DO 260 I=1,12
        DO 260 J=1,30
          AASO(I,J)=LOG10(AASO(I,J)+1.D-18)
          AVMULT(I,J)=LOG10(AVMULT(I,J)+1.D-18)
          AKE(I,J)=LOG10(AKE(I,J)+1.D-18)
  260 CONTINUE
      WRITE(6,*)' 1=*=p, 2=n, 3=pi+, 4=pi-, 5=K+, 6=K-,',
     &' 8=ap, 9=chargd., 10=Z=all, 11=+=Lamb, 12=A=aLam,',
     &' 13=O=sig-, 14=B=Sig+, 15=C=Sig0, 16=D=pi0'
      CALL PLOT(AASO,AVMULT,360,30,12,0.D0,0.1D0,-3.D0,0.05D0)
      WRITE(6,*)' 1=*=p, 2=n, 3=pi+, 4=pi-, 5=K+, 6=K-,',
     &' 8=ap, 9=chargd., 10=Z=all, 11=+=Lamb, 12=A=aLam,',
     &' 13=O=sig-, 14=B=Sig+, 15=C=Sig0, 16=D=pi0'
      CALL PLOT(AASO,AKE,360,30,12,0.D0,0.1D0,-5.D0,0.05D0)
  270 CONTINUE
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE DDECAY(IHAD,ISTAB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C------------------
*KEEP,DFINPA.
      CHARACTER*8 ANF
      PARAMETER (NFIMAX=249)
      COMMON /DFINPA/ ANF(NFIMAX),PXF(NFIMAX),PYF(NFIMAX),PZF(NFIMAX),
     +HEF(NFIMAX),AMF(NFIMAX), ICHF(NFIMAX),IBARF(NFIMAX),NREF(NFIMAX)
      COMMON /DFINPZ/IORMO(NFIMAX),IDAUG1(NFIMAX),IDAUG2(NFIMAX),
     * ISTATH(NFIMAX)
*KEND.
      CHARACTER*8 ZKNAME
      CHARACTER*8 ANAME
C      COMMON/DDECAC/ ZKNAME(540),NZK(540,3),WT(540)
 
      PARAMETER (IDMAX9=602)
C      CHARACTER*8 ZKNAME
      COMMON/DDECAC/ ZKNAME(IDMAX9),WT(IDMAX9),NZK(IDMAX9,3)
 
 
      COMMON/DPAR/ANAME(210),AM(210),GA(210),TAU(210),ICH(210),IBAR(210)
     *,K1(210),K2(210)
      COMMON/DMETLS/ CXS(149),CYS(149),CZS(149),ELS(149),
     *PLS(149),IS,ITS(149)
      COMMON/DDRE/ TEST(12)
      COMMON/PRITT/ISYS
C------------------
      ISYS=6
      DO 10 I=1,IHAD
        ITS(I)=NREF(I)
        PLS(I)=SQRT(PXF(I)**2+PYF(I)**2+PZF(I)**2)
        IF(PLS(I).NE.0.)CXS(I)=PXF(I)/PLS(I)
        IF(PLS(I).NE.0.)CYS(I)=PYF(I)/PLS(I)
        IF(PLS(I).NE.0.)CZS(I)=PZF(I)/PLS(I)
        ELS(I)=HEF(I)
   10 CONTINUE
      IST=IHAD
      IR=0
   20 CONTINUE
C*****TEST STABLE OR UNSTABLE
C   ISTAB=1/2/3 MEANS  STRONG + WEAK DECAYS / ONLY STRONG DECAYS /
C   STRONG DECAYS + WEAK DECAYS FOR CHARMED PARTICLES AND TAU LEPTONS
      IF(ISTAB.EQ.1)                                             GOTO 30
      IF(ISTAB.EQ.2)                                             GOTO 50
      IF(ISTAB.EQ.3)                                             GOTO 40
   30 IF(ITS(IST).EQ.135.OR.ITS(IST).EQ.136)                     GOTO 60
      IF(ITS(IST).GE.1.AND.ITS(IST).LE.7)                        GOTO 60
                                                                 GOTO 70
   40 IF(ITS(IST).GE.1.AND.ITS(IST).LE.23)                       GOTO 60
      IF(ITS(IST).GE. 97.AND.ITS(IST).LE.103)                    GOTO 60
C*    IF(ITS(IST).EQ.109.AND.ITS(IST).EQ.115) GOTO 202
      IF(ITS(IST).EQ.109.OR.ITS(IST).EQ.115)                     GOTO 60
      IF(ITS(IST).GE.133.AND.ITS(IST).LE.136)                    GOTO 60
                                                                 GOTO 70
   50 IF(ITS(IST).GE. 1.AND.ITS(IST).LE. 30)                     GOTO 60
      IF(ITS(IST).GE. 97.AND.ITS(IST).LE.103)                    GOTO 60
      IF(ITS(IST).GE.115.AND.ITS(IST).LE.122)                    GOTO 60
      IF(ITS(IST).GE.131.AND.ITS(IST).LE.136)                    GOTO 60
      IF(ITS(IST).EQ.109)                                        GOTO 60
      IF(ITS(IST).GE.137.AND.ITS(IST).LE.160)                    GOTO 60
                                                                 GOTO 70
   60 IR=IR+1
      IF (IR.GT.NFIMAX)THEN
        WRITE (6,1000)IR,NFIMAX
 1000 FORMAT(' DECAY IR.GT.NFIMAX RETURN ',2I10)
        RETURN
      ENDIF
      NREF(IR)=ITS(IST)
      ITT=ITS(IST)
      AMF(IR)=AM(ITT)
      ANF(IR)=ANAME(ITT)
      ICHF(IR)=ICH(ITT)
      IBARF(IR)=IBAR(ITT)
      HEF(IR)=ELS(IST)
      PXF(IR)=CXS(IST)*PLS(IST)
      PYF(IR)=CYS(IST)*PLS(IST)
      PZF(IR)=CZS(IST)*PLS(IST)
      IST=IST-1
      IF(IST.GE.1)                                               GOTO 20
                                                                 GOTO140
   70 IT=ITS(IST)
      GAM=ELS(IST)/AM(IT)
      BGAM=PLS(IST)/AM(IT)
      ECO=AM(IT)
      KZ1=K1(IT)
   80 CONTINUE
      VV=RNDM(VW)-1.D-17
      IIK=KZ1-1
   90 IIK=IIK+1
      IF (VV.GT.WT(IIK))                                        GO TO 90
C  IIK IS THE DECAY CHANNEL
      IT1=NZK(IIK,1)
      IT2=NZK(IIK,2)
      IF (IT2-1.LT.0)                                          GO TO 120
      IT3=NZK(IIK,3)
C  IT1,IT2, IT3 ARE THE PRODUCED PARTICLES FROM  IT
      IF(IT3.EQ.0)                                             GO TO 100
      CALL DTHREP(ECO,ECM1,ECM2,ECM3,PCM1,PCM2,PCM3,COD1,COF1,SIF1,
     *COD2,COF2,SIF2,COD3,COF3,SIF3,AM(IT1),AM(IT2),AM(IT3))
                                                               GO TO 110
  100 CALL DTWOPD(ECO,ECM1,ECM2,PCM1,PCM2,COD1,COF1,SIF1,COD2,COF2,SIF2,
     +AM(IT1),AM(IT2))
  110 CONTINUE
  120 CONTINUE
      ITS(IST  )=IT1
      IF (IT2-1.LT.0)                                          GO TO 130
      ITS(IST+1)  =IT2
      ITS(IST+2)=IT3
      RX=CXS(IST)
      RY=CYS(IST)
      RZ=CZS(IST)
      CALL DTRAFO(GAM,BGAM,RX,RY,RZ,COD1,COF1,SIF1,PCM1,ECM1,
     *PLS(IST),CXS(IST),CYS(IST),CZS(IST),ELS(IST))
      IST=IST+1
      CALL DTRAFO(GAM,BGAM,RX,RY,RZ,COD2,COF2,SIF2,PCM2,ECM2,
     *PLS(IST),CXS(IST),CYS(IST),CZS(IST),ELS(IST))
      IF (IT3.LE.0)                                            GO TO 130
      IST=IST+1
      CALL DTRAFO(GAM,BGAM,RX,RY,RZ,COD3,COF3,SIF3,PCM3,ECM3,
     *PLS(IST),CXS(IST),CYS(IST),CZS(IST),ELS(IST))
  130 CONTINUE
                                                                GO TO 20
  140 CONTINUE
      IF(IR.GT.7998) WRITE(ISYS,1010)
 1010 FORMAT(2X,'  NUMBER OF STAB. FINAL PART. IS GREATER THAN 7998')
      IHAD=IR
      RETURN
      END
*-- Author :
C--------------------------------------------------------------------
C
C     FILE SHMAK
C
C--------------------------------------------------------------------
      SUBROUTINE SHMAK(ICASE,NN,NNA,NNB,NA,NB,UMO,BIMP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
*  scoring of unbiased Glauber events sampled in KKEVT
*  (reduction of interactions possible in case event rejection because
*   of limitations from kinematics)
*
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEND.
      PARAMETER(NAMX=248)
      DIMENSION FNUA(NAMX),FNUB(NAMX),FNUT(NAMX)
      DIMENSION ANN(NAMX)
      DIMENSION XB(200),BIMPP(200)
C-------------------------------
C---------------------------------------------------------------
C
C                  plot impact parameter distribution
C
C----------------------------------------------------------------
C     DO 7784 II=1,200
C       BIMPP(II)=0.D0
CXB(II)=0.1D0*II
C7784 CONTINUE
C     IP=207
C     IT=207
C     KKMAT=1
C     PPROJ=PPN
C     DO 7785 II=1,100000
C     CALL SHMAKO(IP,IT,BIMP,NN,NP,NT,JSSH,JTSH,PPROJ,KKMAT)	
C     IF(II.LE.1000)WRITE(6,*)' IP,IT,BIMP,NN,NP,NT ',
C    * IP,IT,BIMP,NN,NP,NT      
C     IB=BIMP/0.1D0+1.D0
C     IF(IB.GE.200)IB=200
C     BIMPP(IB)=BIMPP(IB)+1
C7785 CONTINUE 
C     WRITE(6,*)' Impact parameter distribution B,BIMPP'
C     DO 7786 II=1,200
C     WRITE(6,*)XB(II),BIMPP(II)
C7786 CONTINUE     
C---------------------------------------------------------------
C
C                  plot impact parameter distribution
C
C---------------------------------------------------------------
      GO TO (10,30,40),ICASE
   10 CONTINUE
      DO 7784 II=1,200
        BIMPP(II)=0.D0
        XB(II)=0.1D0*II
 7784 CONTINUE
      BNUT=0.
      BNUA=0.
      BNUB=0.
      BNVV=0.
      BNSV=0.
      BNVS=0.
      BNSS=0.
      DO 20 I=1,NAMX
	ANN(I)=I
        FNU A(I)=0.
        FNU B(I)=0.
        FNU T(I)=0.
   20 CONTINUE
      ANUSD=0.D0
      RETURN
   30 CONTINUE
      IB=BIMP/0.1D0+1.D0
      IF(IB.GE.200)IB=200
      BIMPP(IB)=BIMPP(IB)+1
C     Calculate fraction of diffractive events
      IF((NN.EQ.1).AND.(NNA.EQ.1).AND.(NNB.EQ.1))THEN
        CALL SIHNDI(UMO,1,1,SINGDIF,SIGDIH) 
	SIGABS=SIINEL(1,1,UMO)
        ANUSD=ANUSD + SINGDIF/SIGABS
      ENDIF
      INTT=NN
      IF (INTT.GT.NAMX)INTT=NAMX
      NUA=NNA
      IF (NUA.GT.NAMX) NUA=NAMX
      NUB=NNB
      IF (NUB.GT.NAMX) NUB=NAMX
      FNUA(NUA)=FNUA(NUA)+1.
      FNUT(INTT)=FNUT(INTT)+1.
      FNUB(NUB)=FNUB(NUB)+1.
      BNUT=BNUT+NN
      BNUA=BNUA+NNA
      BNUB=BNUB+NNB
      IF(NNB.GE.NNA) THEN
        NNVV=NNA
        NNSV=NNB-NNA
        NNVS=0
        NNSS=NN-NNB
      ELSE
        NNVV=NNB
        NNSV=0
        NNVS=NNA-NNB
        NNSS=NN-NNA
      ENDIF
      BNVV=BNVV + NNVV
      BNSV=BNSV + NNSV
      BNVS=BNVS + NNVS
      BNSS=BNSS + NNSS
      RETURN
   40 CONTINUE
      IF(NN.EQ.0)THEN
	WRITE(6,*)' shmak(3,NN,... ) NN= ',NN
	RETURN
      ENDIF
C     WRITE(6,*)' Impact parameter distribution B,BIMPP'
C     WRITE(6,*)(XB(II),BIMPP(II),II=1,200)
C     WRITE(6,*)' Impact parameter distribution B,BIMPP'
C       CALL PLOT(XB,BIMPP,50,1,50,0.,0.4D0,0.,10.D0)
      ANUSD=ANUSD/NN
      BNUT=BNUT/NN
      BNUA=BNUA/NN
      BNUB=BNUB/NN
      BNVV=BNVV/NN
      BNSV=BNSV/NN
      BNVS=BNVS/NN
      BNSS=BNSS/NN
      WRITE(6,'(1H1,50(1H*))')
      WRITE(6,'(/10X,A/)') ' OUTPUT FROM SHMAK  all events before',
     *' diffraction modification'
      WRITE(6,'(50(1H*))')
      WRITE(6,'(A,I10)') ' NUMBER OF TOTALLY SAMPLED GLAUBER EVENTS',NN
      WRITE(6, 1000) BNUT,BNUA,BNUB
      WRITE(6,*)' Fraction of diffractive evnts: ',ANUSD
 1000 FORMAT('  AVERAGE NO OF COLLISIONS BNUT,BNUA,BNUB=',3F10.3)
      WRITE(6,'(/A)') ' AVERAGE NUMBERS OF DIFFERENT COLLISION TYPES'
      WRITE(6,'(4(5X,A,F8.2/))') ' VAL-VAL:',BNVV, ' SEA-VAL:',BNSV,
     +' VAL-SEA:',BNVS, ' SEA-SEA:',BNSS
      IF(IPRI.GE.1) THEN
        DNNA=NA/50+1
        DNNB=NB/50+1
        DNNT=2.*DNNB
        WRITE(6,1010)
 1010   FORMAT (' FNUA')
        WRITE(6,1040) FNUA
        DO 323 I=1,NAMX
          FNU A(I)=LOG10(FNU A(I)+1.D-5)
 323    CONTINUE
        CALL PLOT(ANN,FNU A,NAMX,1,NAMX,0.D0,DNNA,0.D0,0.05D0)
        WRITE(6,1020 )
        WRITE(6,1040) FNUB
 1020   FORMAT (' FNUB')
        DO 324 I=1,NAMX
          FNU B(I)=LOG10(FNU B(I)+1.D-5)
 324    CONTINUE
        CALL PLOT(ANN,FNU B,NAMX,1,NAMX,0.D0,DNNB,0.D0,0.05D0)
        WRITE(6,1030 )
        WRITE(6,1040) FNUT
 1030    FORMAT (' FNUT')
        DO 325 I=1,NAMX
          FNU T(I)=LOG10(FNU T(I)+1.E-5)
 325     CONTINUE
        CALL PLOT(ANN,FNU T,NAMX,1,NAMX,0.D0,DNNT,0.D0,0.05D0)
 1040    FORMAT (10F12.2)
      ENDIF
      RETURN
      END
C--------------------------------------------------------------------
C
C     FILE SHMAK1
C
C--------------------------------------------------------------------
      SUBROUTINE SHMAK1(ICASE,NN,NNA,NNB,NA,NB,UMO,BIMP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
*  scoring of unbiased Glauber events sampled in KKEVT
*  (reduction of interactions possible in case event rejection because
*   of limitations from kinematics)
*
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEND.
      PARAMETER(NAMX=248)
      DIMENSION FNUA(NAMX),FNUB(NAMX),FNUT(NAMX)
      DIMENSION ANN(NAMX)
      DIMENSION XB(200),BIMPP(200)
C-------------------------------
      GO TO (10,30,40),ICASE
   10 CONTINUE
      DO 7784 II=1,200
        BIMPP(II)=0.D0
        XB(II)=0.1D0*II
 7784 CONTINUE
      BNUT=0.
      BNUA=0.
      BNUB=0.
      ANUSD=0.D0
      DO 20 I=1,NAMX
        ANN(I)=I
        FNU A(I)=0.
        FNU B(I)=0.
        FNU T(I)=0.
   20 CONTINUE
      RETURN
   30 CONTINUE
      IB=BIMP/0.1D0+1.D0
      IF(IB.GE.200)IB=200
      BIMPP(IB)=BIMPP(IB)+1
C     Calculate fraction of diffractive events
      IF((NN.EQ.1).AND.(NNA.EQ.1).AND.(NNB.EQ.1))THEN
        CALL SIHNDI(UMO,1,1,SINGDIF,SIGDIH) 
	SIGABS=SIINEL(1,1,UMO)
        ANUSD=ANUSD + SINGDIF/SIGABS
      ENDIF
      INTT=NN
      IF (INTT.GT.NAMX)INTT=NAMX
      NUA=NNA
      IF (NUA.GT.NAMX) NUA=NAMX
      NUB=NNB
      IF (NUB.GT.NAMX) NUB=NAMX
      FNUA(NUA)=FNUA(NUA)+1.
      FNUT(INTT)=FNUT(INTT)+1.
      FNUB(NUB)=FNUB(NUB)+1.
      BNUT=BNUT+NN
      BNUA=BNUA+NNA
      BNUB=BNUB+NNB
      RETURN
   40 CONTINUE
      IF(NN.EQ.0)THEN
	WRITE(6,*)' shmak1(3,NN,... ) NN= ',NN
	RETURN
      ENDIF
C     WRITE(6,*)' Impact parameter distribution B,BIMPP'
C     WRITE(6,*)(XB(II),BIMPP(II),II=1,200)
C     WRITE(6,*)' Impact parameter distribution B,BIMPP'
C       CALL PLOT(XB,BIMPP,50,1,50,0.D0,0.4D0,0.D0,10.D0)
      BNUT=BNUT/NN
      BNUA=BNUA/NN
      BNUB=BNUB/NN
      ANUSD=ANUSD/NN
      WRITE(6,'(1H1,50(1H*))')
      WRITE(6,'(/10X,A/)') ' OUTPUT FROM SHMAK1 after modification',
     *' of Glauber events for diffractive cross section'
      WRITE(6,'(50(1H*))')
      WRITE(6,'(A,I10)') ' NUMBER OF TOTALLY SAMPLED GLAUBER EVENTS',NN
      WRITE(6, 1000) BNUT,BNUA,BNUB
 1000 FORMAT('  AVERAGE NO OF COLLISIONS BNUT,BNUA,BNUB=',3F10.3)
      WRITE(6,*)' Fraction of diffractive evnts: ',ANUSD
C     WRITE(6,'(/A)') ' AVERAGE NUMBERS OF DIFFERENT COLLISION TYPES'
C     WRITE(6,'(4(5X,A,F8.2/))') ' VAL-VAL:',BNVV, ' SEA-VAL:',BNSV,
C    +' VAL-SEA:',BNVS, ' SEA-SEA:',BNSS
      IF(IPRI.GE.1) THEN
        DNNA=NA/50+1
        DNNB=NB/50+1
        DNNT=2.*DNNB
        WRITE(6,1010)
 1010   FORMAT (' FNUA')
        WRITE(6,1040) FNUA
        DO 323 I=1,NAMX
          FNU A(I)=LOG10(FNU A(I)+1.D-5)
323     CONTINUE
        CALL PLOT(ANN,FNU A,NAMX,1,NAMX,0.D0,DNNA,0.D0,0.05D0)
        WRITE(6,1020 )
        WRITE(6,1040) FNUB
 1020   FORMAT (' FNUB')
        DO 324 I=1,NAMX
          FNU B(I)=LOG10(FNU B(I)+1.D-5)
324     CONTINUE
        CALL PLOT(ANN,FNU B,NAMX,1,NAMX,0.D0,DNNB,0.D0,0.05D0)
        WRITE(6,1030 )
        WRITE(6,1040) FNUT
1030    FORMAT (' FNUT')
        DO 325 I=1,NAMX
          FNU T(I)=LOG10(FNU T(I)+1.E-5)
325     CONTINUE
        CALL PLOT(ANN,FNU T,NAMX,1,NAMX,0.D0,DNNT,0.D0,0.05D0)
1040    FORMAT (10F12.2)
      ENDIF
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE PREVIO(RA,RB,NSTB,BMAX,BSTEP,SIG,RO,G)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      COMPLEX*16 CA,CI
      COMMON/DAMP/CA,CI,GA
C     WRITE(6,*)' PREVIO(RA,RB,NSTB,BMAX,BSTEP,SIG,RO,G)',
C    &RA,RB,NSTB,BMAX,BSTEP,SIG,RO,G
C     WRITE(6,*)' /CA,CI,GA/ ',CA,CI,GA
C     WRITE(6,*)' PREVIO: RA, RB = ',RA,RB
      BMAX=4.*(RA+RB)
C     WRITE(6,*)' PREVIO(RA,RB,NSTB,BMAX,BSTEP,SIG,RO,G)',
C    &RA,RB,NSTB,BMAX,BSTEP,SIG,RO,G
C     WRITE(6,*)' /CA,CI,GA/ ',CA,CI,GA
C     WRITE(6,*)' PREVIO: RA, RB BMAX= ',RA,RB,BMAX
      BSTEP=BMAX/(NSTB-1)
C     WRITE(6,*)' PREVIO(RA,RB,NSTB,BMAX,BSTEP,SIG,RO,G)',
C    &RA,RB,NSTB,BMAX,BSTEP,SIG,RO,G
C     WRITE(6,*)' /CA,CI,GA/ ',CA,CI,GA
C     WRITE(6,*)' PREVIO: RA, RB ,BSTEP= ',RA,RB,BSTEP
      BSTEP=0.15D0
C     WRITE(6,*)' PREVIO(RA,RB,NSTB,BMAX,BSTEP,SIG,RO,G)',
C    &RA,RB,NSTB,BMAX,BSTEP,SIG,RO,G
C     WRITE(6,*)' /CA,CI,GA/ ',CA,CI,GA
C     WRITE(6,*)' PREVIO: RA, RB ,BSTEP= ',RA,RB,BSTEP
      GA=G
C     WRITE(6,*)' PREVIO(RA,RB,NSTB,BMAX,BSTEP,SIG,RO,G)',
C    &RA,RB,NSTB,BMAX,BSTEP,SIG,RO,G
C     WRITE(6,*)' /CA,CI,GA/ ',CA,CI,GA
C     WRITE(6,*)' PREVIO: RA, RB ,GA= ',RA,RB,GA
      RCA=GA*SIG/6.2831854
C     WRITE(6,*)' PREVIO(RA,RB,NSTB,BMAX,BSTEP,SIG,RO,G)',
C    &RA,RB,NSTB,BMAX,BSTEP,SIG,RO,G
C     WRITE(6,*)' /CA,CI,GA/ ',CA,CI,GA
C     WRITE(6,*)' PREVIO: RA, RB ,RCA= ',RA,RB,RCA
      FCA=-GA*SIG*RO/6.2831854
C     WRITE(6,*)' PREVIO(RA,RB,NSTB,BMAX,BSTEP,SIG,RO,G)',
C    &RA,RB,NSTB,BMAX,BSTEP,SIG,RO,G
C     WRITE(6,*)' /CA,CI,GA/ ',CA,CI,GA
C     WRITE(6,*)' PREVIO: RA, RB ,FCA= ',RA,RB,FCA
      CA=CMPLX(RCA,FCA)
C     WRITE(6,*)' PREVIO(RA,RB,NSTB,BMAX,BSTEP,SIG,RO,G)',
C    &RA,RB,NSTB,BMAX,BSTEP,SIG,RO,G
C     WRITE(6,*)' /CA,CI,GA/ ',CA,CI,GA
C     WRITE(6,*)' PREVIO: RA, RB,CA = ',RA,RB,CA
      CI=(1.D0,0.D0)
      WRITE(6,*)' PREVIO(RA,RB,NSTB,BMAX,BSTEP,SIG,RO,G)',
     &RA,RB,NSTB,BMAX,BSTEP,SIG,RO,G
      WRITE(6,*)' /CA,CI,GA/ ',CA,CI,GA
      WRITE(6,*)' PREVIO: RA, RB ,CI= ',RA,RB,CI
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE PROFB(BSTEP,NSTAT,NA,RA,NB,RB,BSITE,NSITEB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C  THE PROGRAM CALCULATES THE PROFIL-FUNCTION AND FILLS
C  THE ARRAY BSITE TO APPROXIMATE THE B-DISTRIBUTION.
C-------------------
      PARAMETER (INTMX=2488,INTMD=252)
*KEEP,DROPPT.
      LOGICAL INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA, IPADIS,
     +ISHMAL,LPAULI
      COMMON /DROPPT/ INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA,
     +IPADIS,ISHMAL,LPAULI
*KEEP,NUCKOO.
      COMMON /NUCKOO/ PKOO(3,INTMX),TKOO(3,INTMX),PPOO(3,INTMX),
     +TPOO(3,INTMX)
*KEEP,DAMP.
C     COMPLEX*16 CA,CI
      DOUBLE COMPLEX CA,CI
      COMMON /DAMP/   CA,CI,GA
*KEND.
      DIMENSION BSITE(0:1,NSITEB)
C--------
      DIMENSION HELPP(200)
      DIMENSION HELP(200)
      DIMENSION BS(200)
      COMMON /SIGLA/SIGLAU
C     COMPLEX*16 C
      DOUBLE COMPLEX C
      DATA IRW /0/
C--------
      WRITE(6,*)' PROFB: RA, RB = ',RA,RB
      WRITE(6, 1000)BSTEP,NSTAT,NA,RA,NB,RB,IRW,NSITEB
 1000 FORMAT (' PROFB',E15.5,2I10,F15.5,I10,E15.5,2I10)
      NS=NSTAT
      NSITE=NSITEB-1
      BST=BSTEP
      DO 10 I=1,NSITEB
        BS(I)=0.
   10 CONTINUE
      DO 40 I=1,NS
      CALL CONUCL(TKOO,NB,RB)
C     CALL SORTIN(TKOO,NB)
      CALL CONUCL(PKOO,NA,RA)
C     CALL SORT(PKOO,NA)
        DO 40 I3=1,NSITE
          B=I3*BST
          PI=1.
          DO 30 I1=1,NA
          X1=B-PKOO(1,I1)
	  IF(PI.LT.1.D-100)GO TO 31
          X2=-PKOO(2,I1)
            DO 32 I2=1,NB
            Q1=X1+TKOO(1,I2)
            Q2=X2+TKOO(2,I2)
              XY=GA*(Q1*Q1+Q2*Q2)
C
              IF(XY.GT.15.)                                     GO TO 20
              E=EXP(-XY)
              C=CI-CA*E
              AR=REAL(REAL(C))
              AI=IMAG(C)
              P=AR*AR+AI*AI
C	      WRITE(6,'(A,5E13.3,3I6)')' PROFB:Pi,P,AR,AI,Ei,I3,I2,I1',
C    *PI,P,AR,AI,Ei,I3,I2,I1
              PI=PI*P
   20         CONTINUE
   32     CONTINUE
   31     CONTINUE
   30     CONTINUE
          BS(I3+1)=BS(I3+1)+1.-PI
   40 CONTINUE
      BS(1)=BS(2)
      SUMB=0.
      DO 50 I=1,NSITEB
	HELPP(I)=BS(I)/NS
        BS(I)=BS(I)*(I-1)*BST/NS
        SUMB=SUMB+BS(I)
   50 CONTINUE
      BSITE(1,1)=0.
      DO 60 I=2,NSITEB
        BSITE(1,I)=BS(I)/SUMB+BSITE(1,I-1)
   60 CONTINUE
      DO 70 I=1,NSITEB
        HELP(I)=I*BST
   70 CONTINUE
      SUMB=SUMB*BST*6.2831854
      SIGLAU=SUMB*10.
      WRITE(6,1020) SUMB
 1020 FORMAT(/5X,7HSIGMA =,F7.3)
      IF(IRW.GE.1) RETURN
      IF(ISHMAL) THEN
      DO 80 I=1,200
        WRITE (6,1030) HELP(I),HELPP(I),BS(I),BSITE(1,I)
 1030 FORMAT (F10.4,3E15.5)
   80 CONTINUE
        CALL PLOT(HELP,BSITE,50,1,50,0.D0,0.5D0,0.D0,0.01D0)
        CALL PLOT(HELP,BS   ,50,1,50,0.D0,0.5D0,0.D0,0.07D0)
      ENDIF
      RETURN
      END
*-- Author :
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE DPARJE(IHAD,I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
*KEEP,DFINPA.
      CHARACTER*8 ANF
      PARAMETER (NFIMAX=249)
      COMMON /DFINPA/ ANF(NFIMAX),PXF(NFIMAX),PYF(NFIMAX),PZF(NFIMAX),
     +HEF(NFIMAX),AMF(NFIMAX), ICHF(NFIMAX),IBARF(NFIMAX),NREF(NFIMAX)
      COMMON /DFINPZ/IORMO(NFIMAX),IDAUG1(NFIMAX),IDAUG2(NFIMAX),
     * ISTATH(NFIMAX)
*KEEP,DINPDA.
      COMMON /DINPDA/ IMPS(6,6),IMVE(6,6),IB08(6,21),IB10(6,21), IA08
     +(6,21),IA10(6,21), A1,B1,B2,B3,LT,LE,BET,AS,B8,AME,DIQ,ISU
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEND.
      CHARACTER*8 ANAME
      COMMON /DPAR/ ANAME(210),AM(210),GA(210),TAU(210), ICH(210),IBAR
     +(210),K1(210),K2(210)
      IHAD=1
      NREF(1)=I
      PXF(1)=0.
      PYF(1)=0.
      PZF(1)=0.
      HEF(1)=AM(I)
      AMF(1)=AM(I)
      ICHF(1)=ICH(I)
      IBARF(1)=IBAR(I)
      ANF(1)=ANAME(I)
      IF (IPCO.GE.6)THEN
        WRITE(6,1000)IHAD,I,PXF(1),PYF(1),PZF(1),HEF(1),AMF(1)
 1000 FORMAT(' PARJET: IHAD,I,PXF(1),PYF(1),PZF(1),HEP(1),AMF(1)'/ 2I5,5
     +F10.3)
      ENDIF
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE SORT(A,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      DIMENSION A(3,N)
      M=N
   10 CONTINUE
      M=N-1
      IF(M.LE.0) RETURN
      L=0
      DO 20 I=1,M
        J=I+1
        IF (A(3,I).LE.A(3,J))                                   GO TO 20
        B=A(3,I)
        C=A(1,I)
        D=A(2,I)
        A(3,I)=A(3,J)
        A(2,I)=A(2,J)
        A(1,I)=A(1,J)
        A(3,J)=B
        A(1,J)=C
        A(2,J)=D
        L=1
   20 CONTINUE
      IF(L.EQ.1)                                                GO TO 10
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE SORTIN(A,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      DIMENSION A(3,N)
      M=N
   10 CONTINUE
      M=N-1
      IF(M.LE.0) RETURN
      L=0
      DO 20 I=1,M
        J=I+1
        IF (A(3,I).GE.A(3,J))                                   GO TO 20
        B=A(3,I)
        C=A(1,I)
        D=A(2,I)
        A(3,I)=A(3,J)
        A(2,I)=A(2,J)
        A(1,I)=A(1,J)
        A(3,J)=B
        A(1,J)=C
        A(2,J)=D
        L=1
   20 CONTINUE
      IF(L.EQ.1)                                                GO TO 10
      RETURN
      END
*
*=== blkdt6 ===========================================================*
*==                                                                    *
      BLOCK DATA BLKD46
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
*$ CREATE DBLPRC.ADD
*COPY DBLPRC
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
 
 
       CHARACTER*8 ANAME
       COMMON /DPAR/ ANAME(210),AM(210),GA(210),TAU(210),
     +               ICH(210),IBAR(210),K1(210),K2(210)
* / Part /
*     datas     datas    datas      datas     datas                    *
*     ---------------------------------------------                    *
*
*
*     Particle  masses Engel version JETSET compatible                 *
*                                                                      *
      DATA (AM(K),K=1,85) /
     &   .9383D+00, .9383D+00,  AMELCT  ,  AMELCT  , .0000D+00,
     &   .0000D+00, .0000D+00, .9396D+00, .9396D+00, AMMUON   ,
     &   AMMUON   , .4977D+00, .1396D+00, .1396D+00, .4936D+00,
     &   .4936D+00, .1116D+01, .1116D+01, .4977D+00, .1197D+01,
     &   .1189D+01, .1193D+01, .1350D+00, .4977D+00, .4977D+00,
     &   .0000D+00, .0000D+00, .0000D+00, .0000D+00, .0000D+00,
     &   .5488D+00, .7669D+00, .7700D+00, .7669D+00, .7820D+00,
     &   .8921D+00, .8962D+00, .8921D+00, .8962D+00, .1300D+01,
     &   .1300D+01, .1300D+01, .1300D+01, .1421D+01, .1421D+01,
     &   .1421D+01, .1421D+01, .1383D+01, .1384D+01, .1387D+01,
     &   .1820D+01, .2030D+01, .1231D+01, .1232D+01, .1233D+01,
     &   .1234D+01, .1675D+01, .1675D+01, .1675D+01, .1675D+01,
     &   .1500D+01, .1500D+01, .1515D+01, .1515D+01, .1775D+01,
     &   .1775D+01, .1231D+01, .1232D+01, .1233D+01, .1234D+01,
     &   .1675D+01, .1675D+01, .1675D+01, .1675D+01, .1515D+01,
     &   .1515D+01, .2500D+01, .4890D+00, .4890D+00, .4890D+00,
     &   .1300D+01, .1300D+01, .1300D+01, .1300D+01, .2200D+01  /
      DATA (AM(K),K=86,183) /
     &   .2200D+01, .2200D+01, .2200D+01, .1700D+01, .1700D+01,
     &   .1700D+01, .1700D+01, .1820D+01, .2030D+01, .9575D+00,
     &   .1019D+01, .1315D+01, .1321D+01, .1189D+01, .1193D+01,
     &   .1197D+01, .1315D+01, .1321D+01, .1383D+01, .1384D+01,
     &   .1387D+01, .1532D+01, .1535D+01, .1672D+01, .1383D+01,
     &   .1384D+01, .1387D+01, .1532D+01, .1535D+01, .1672D+01,
     &   .1865D+01, .1869D+01, .1869D+01, .1865D+01, .1969D+01,
     &   .1969D+01, .2980D+01, .2007D+01, .2010D+01, .2010D+01,
     &   .2007D+01, .2113D+01, .2113D+01, .3686D+01, .3097D+01,
     &   .1777D+01, .1777D+01, .0000D+00, .0000D+00, .0000D+00,
     &   .0000D+00, .2285D+01, .2460D+01, .2460D+01, .2452D+01,
     &   .2453D+01, .2454D+01, .2560D+01, .2560D+01, .2730D+01,
     &   .3610D+01, .3610D+01, .3790D+01, .2285D+01, .2460D+01,
     &   .2460D+01, .2452D+01, .2453D+01, .2454D+01, .2560D+01,
     &   .2560D+01, .2730D+01, .3610D+01, .3610D+01, .3790D+01,
     &   .2490D+01, .2490D+01, .2490D+01, .2610D+01, .2610D+01,
     &   .2770D+01, .3670D+01, .3670D+01, .3850D+01, .4890D+01,
     &   .2490D+01, .2490D+01, .2490D+01, .2610D+01, .2610D+01,
     &   .2770D+01, .3670D+01, .3670D+01, .3850D+01, .4890D+01,
     &   .1250D+01, .1250D+01, .1250D+01  /
      DATA ( AM ( I ), I = 184,210 ) /
     & 1.44000000000000D+00, 1.44000000000000D+00, 1.30000000000000D+00,
     & 1.30000000000000D+00, 1.30000000000000D+00, 1.40000000000000D+00,
     & 1.46000000000000D+00, 1.46000000000000D+00, 1.46000000000000D+00,
     & 1.46000000000000D+00, 1.60000000000000D+00, 1.60000000000000D+00,
     & 1.66000000000000D+00, 1.66000000000000D+00, 1.66000000000000D+00,
     & 1.66000000000000D+00, 1.66000000000000D+00, 1.66000000000000D+00,
     & 1.95000000000000D+00, 1.95000000000000D+00, 1.95000000000000D+00,
     & 1.95000000000000D+00, 2.25000000000000D+00, 2.25000000000000D+00,
     & 1.44000000000000D+00, 1.44000000000000D+00, 0.00000000000000D+00/
*                                                                      *
*     Particle  mean lives                                             *
*                                                                      *
      DATA (TAU(K),K=1,183) /
     &   .1000D+19, .1000D+19, .1000D+19, .1000D+19, .1000D+19,
     &   .1000D+19, .1000D+19, .9180D+03, .9180D+03, .2200D-05,
     &   .2200D-05, .5200D-07, .2600D-07, .2600D-07, .1200D-07,
     &   .1200D-07, .2600D-09, .2600D-09, .9000D-10, .1500D-09,
     &   .8000D-10, .5000D-14, .8000D-16, .0000D+00, .0000D+00,
     &   70*.0000D+00,
     &   .0000D+00, .3000D-09, .1700D-09, .8000D-10, .1000D-13,
     &   .1500D-09, .3000D-09, .1700D-09, .0000D+00, .0000D+00,
     &   .0000D+00, .0000D+00, .0000D+00, .1000D-09, .0000D+00,
     &   .0000D+00, .0000D+00, .0000D+00, .0000D+00, .1000D-09,
     &   .0000D+00, .0000D+00, .0000D+00, .0000D+00, .0000D+00,
     &   .0000D+00, .0000D+00, .0000D+00, .0000D+00, .0000D+00,
     &   .0000D+00, .0000D+00, .0000D+00, .0000D+00, .0000D+00,
     &   .9000D-11, .9000D-11, .9000D-11, .9000D-11, .1000D+19,
     &   .1000D+19, .0000D+00, .0000D+00, .0000D+00, .0000D+00,
     &   40*.0000D+00,
     &   .0000D+00, .0000D+00, .0000D+00  /
      DATA ( TAU ( I ), I = 184,210 ) /
     & 0.00000000000000D+00, 0.00000000000000D+00, 0.00000000000000D+00,
     & 0.00000000000000D+00, 0.00000000000000D+00, 0.00000000000000D+00,
     & 0.00000000000000D+00, 0.00000000000000D+00, 0.00000000000000D+00,
     & 0.00000000000000D+00, 0.00000000000000D+00, 0.00000000000000D+00,
     & 0.00000000000000D+00, 0.00000000000000D+00, 0.00000000000000D+00,
     & 0.00000000000000D+00, 0.00000000000000D+00, 0.00000000000000D+00,
     & 0.00000000000000D+00, 0.00000000000000D+00, 0.00000000000000D+00,
     & 0.00000000000000D+00, 0.00000000000000D+00, 0.00000000000000D+00,
     & 0.00000000000000D+00, 0.00000000000000D+00, 0.00000000000000D+00/
*                                                                      *
*     Resonance width Gamma in GeV                                     *
*                                                                      *
      DATA (GA(K),K=  1,85) /
     &    30*.0000D+00,
     &   .8500D-06, .1520D+00, .1520D+00, .1520D+00, .1000D-01,
     &   .7900D-01, .7900D-01, .7900D-01, .7900D-01, .4500D+00,
     &   .4500D+00, .4500D+00, .4500D+00, .1080D+00, .1080D+00,
     &   .1080D+00, .1080D+00, .5000D-01, .5000D-01, .5000D-01,
     &   .8500D-01, .1800D+00, .1150D+00, .1150D+00, .1150D+00,
     &   .1150D+00, .2000D+00, .2000D+00, .2000D+00, .2000D+00,
     &   .2000D+00, .2000D+00, .1000D+00, .1000D+00, .2000D+00,
     &   .2000D+00, .1150D+00, .1150D+00, .1150D+00, .1150D+00,
     &   .2000D+00, .2000D+00, .2000D+00, .2000D+00, .1000D+00,
     &   .1000D+00, .2000D+00, .1000D+00, .1000D+00, .1000D+00,
     &   .1000D+00, .1000D+00, .1000D+00, .1000D+00, .2000D+00  /
      DATA (GA(K),K= 86,183) /
     &   .2000D+00, .2000D+00, .2000D+00, .1500D+00, .1500D+00,
     &   .1500D+00, .1500D+00, .8500D-01, .1800D+00, .2000D-02,
     &   .4000D-02, .0000D+00, .0000D+00, .0000D+00, .0000D+00,
     &   .0000D+00, .0000D+00, .0000D+00, .3400D-01, .3400D-01,
     &   .3600D-01, .9000D-02, .9000D-02, .0000D+00, .3400D-01,
     &   .3400D-01, .3600D-01, .9000D-02, .9000D-02, .0000D+00,
     &   .0000D+00, .0000D+00, .0000D+00, .0000D+00, .0000D+00,
     &   .0000D+00, .0000D+00, .5000D-02, .2000D-02, .2000D-02,
     &   .5000D-02, .2000D-02, .2000D-02, .2000D-03, .7000D-03,
     &   50*.0000D+00,
     &   .3000D+00, .3000D+00, .3000D+00  /
      DATA ( GA ( I ), I = 184,210 ) /
     & 2.00000000000000D-01, 2.00000000000000D-01, 3.00000000000000D-01,
     & 3.00000000000000D-01, 3.00000000000000D-01, 2.70000000000000D-01,
     & 2.50000000000000D-01, 2.50000000000000D-01, 2.50000000000000D-01,
     & 2.50000000000000D-01, 1.50000000000000D-01, 1.50000000000000D-01,
     & 1.00000000000000D-01, 1.00000000000000D-01, 1.00000000000000D-01,
     & 1.00000000000000D-01, 1.00000000000000D-01, 1.00000000000000D-01,
     & 6.00000000000000D-02, 6.00000000000000D-02, 6.00000000000000D-02,
     & 6.00000000000000D-02, 5.50000000000000D-02, 5.50000000000000D-02,
     & 2.00000000000000D-01, 2.00000000000000D-01, 0.00000000000000D+00/
*                                                                      *
*     Particle  names                                                  *
*                                                                      *
*     S+1385+Sigma+(1385)    L02030+Lambda0(2030)                      *
*     Rho77=Rho(770) Om783=Omega(783) K*14=K*(1420) and so on          *
*     designation N*@@ means N*@1(@2)                                  *
*                                                                      *
*                                                                      *
      DATA (ANAME(K),K=1,85) /
     &  'P       ','AP      ','E-      ','E+      ','NUE     ',
     &  'ANUE    ','GAM     ','NEU     ','ANEU    ','MUE+    ',
     &  'MUE-    ','K0L     ','PI+     ','PI-     ','K+      ',
     &  'K-      ','LAM     ','ALAM    ','K0S     ','SIGM-   ',
     &  'SIGM+   ','SIGM0   ','PI0     ','K0      ','AK0     ',
     &  'BLANK   ','BLANK   ','BLANK   ','BLANK   ','BLANK   ',
     &  'ETA550  ','RHO+77  ','RHO077  ','RHO-77  ','OM0783  ',
     &  'K*+892  ','K*0892  ','K*-892  ','AK*089  ','KA+125  ',
     &  'KA0125  ','KA-125  ','AKA012  ','K*+142  ','K*0142  ',
     &  'K*-142  ','AK*014  ','S+1385  ','S01385  ','S-1385  ',
     &  'L01820  ','L02030  ','N*++12  ','N*+ 12  ','N*012   ',
     &  'N*-12   ','N*++16  ','N*+16   ','N*016   ','N*-16   ',
     &  'N*+14   ','N*014   ','N*+15   ','N*015   ','N*+18   ',
     &  'N*018   ','AN--12  ','AN*-12  ','AN*012  ','AN*+12  ',
     &  'AN--16  ','AN*-16  ','AN*016  ','AN*+16  ','AN*-15  ',
     &  'AN*015  ','DE*=24  ','RPI+49  ','RPI049  ','RPI-49  ',
     &  'PIN++   ','PIN+0   ','PIN+-   ','PIN-0   ','PPPI    ' /
      DATA (ANAME(K),K=86,183) /
     &  'PNPI    ','APPPI   ','APNPI   ','K+PPI   ','K-PPI   ',
     &  'K+NPI   ','K-NPI   ','S+1820  ','S-2030  ','ETA*    ',
     &  'PHI     ','TETA0   ','TETA-   ','ASIG-   ','ASIG0   ',
     &  'ASIG+   ','ATETA0  ','ATETA+  ','SIG*+   ','SIG*0   ',
     &  'SIG*-   ','TETA*0  ','TETA*   ','OMEGA-  ','ASIG*-  ',
     &  'ASIG*0  ','ASIG*+  ','ATET*0  ','ATET*+  ','OMEGA+  ',
     &  'D0      ','D+      ','D-      ','AD0     ','DS+     ',
     &  'DS-     ','ETAC    ','D*0     ','D*+     ','D*-     ',
     &  'AD*0    ','DS*+    ','DS*-    ','CHI1C   ','JPSI    ',
     &  'TAU+    ','TAU-    ','NUET    ','ANUET   ','NUEM    ',
     &  'ANUEM   ','LAMC+   ','XIC+    ','XIC0    ','SIGC++  ',
     &  'SIGC+   ','SIGC0   ','S+      ','S0      ','T0      ',
     &  'XU++    ','XD+     ','XS+     ','ALAMC-  ','AXIC-   ',
     &  'AXIC0   ','ASIGC-- ','ASIGC-  ','ASIGC0  ','AS-     ',
     &  'AS0     ','AT0     ','AXU--   ','AXD-    ','AXS     ',
     &  'C1*++   ','C1*+    ','C1*0    ','S*+     ','S*0     ',
     &  'T*0     ','XU*++   ','XD*+    ','XS*+    ','TETA++  ',
     &  'AC1*--  ','AC1*-   ','AC1*0   ','AS*-    ','AS*0    ',
     &  'AT*0    ','AXU*--  ','AXD*-   ','AXS*-   ','ATET--  ',
     &  'RO      ','R+      ','R-      '  /
      DATA (    ANAME ( I ), I = 184,210 ) /
     &'AN*-14  ','AN*014  ','PI+130  ','PI0130  ','PI-130  ','F01400  ',
     &'K*+146  ','K*-146  ','K*0146  ','AK0146  ','L01600  ','AL0160  ',
     &'S+1660  ','S01660  ','S-1660  ','AS-166  ','AS0166  ','AS+166  ',
     &'X01950  ','X-1950  ','AX0195  ','AX+195  ','OM-225  ','AOM+22  ',
     &'N*+14   ','N*014   ','BLANK   '/
*                                                                      *
*     Charge of particles and resonances                               *
*                                                                      *
      DATA ( ICH ( I ), I =   1,210 ) /
     &  1, -1, -1,  1,  0,  0,  0,  0,  0,  1, -1,  0,  1, -1,  1,
     & -1,  0,  0,  0, -1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &  0,  1,  0, -1,  0,  1,  0, -1,  0,  1,  0, -1,  0,  1,  0,
     & -1,  0,  1,  0, -1,  0,  0,  2,  1,  0, -1,  2,  1,  0, -1,
     &  1,  0,  1,  0,  1,  0, -2, -1,  0,  1, -2, -1,  0,  1, -1,
     &  0,  1,  1,  0, -1,  2,  1,  0, -1,  2,  1,  0, -1,  2,  0,
     &  1, -1,  1, -1,  0,  0,  0, -1, -1,  0,  1,  0,  1,  1,  0,
     & -1,  0, -1, -1, -1,  0,  1,  0,  1,  1,  0,  1, -1,  0,  1,
     & -1,  0,  0,  1, -1,  0,  1, -1,  0,  0,  1, -1,  0,  0,  0,
     &  0,  1,  1,  0,  2,  1,  0,  1,  0,  0,  2,  1,  1, -1, -1,
     &  0, -2, -1,  0, -1,  0,  0, -2, -1, -1,  2,  1,  0,  1,  0,
     &  0,  2,  1,  1,  2, -2, -1,  0, -1,  0,  0, -2, -1, -1, -2,
     &  0,  1, -1, -1,  0,  1,  0, -1,  0,  1, -1,  0,  0,  0,  0,
     &  1,  0, -1, -1,  0,  1,  0, -1,  0,  1, -1,  1,  1,  0,  0/
*                                                                      *
*     Particle  baryonic charges                                       *
*                                                                      *
      DATA ( IBAR ( I ), I =   1,210 ) /
     &  1, -1,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     &  0,  1, -1,  0,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,
     &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
     &  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
     & -1,  2,  0,  0,  0,  1,  1,  1,  1,  2,  2,  0,  0,  1,  1,
     &  1,  1,  1,  1,  0,  0,  1,  1, -1, -1, -1, -1, -1,  1,  1,
     &  1,  1,  1,  1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,
     &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1,
     & -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,
     &  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
     &  0,  0,  0, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1,
     &  1,  1,  1, -1, -1, -1,  1,  1, -1, -1,  1, -1,  1,  1,  0/
*                                                                      *
*     First number of decay channels used for resonances               *
*     and decaying particles                                           *
*                                                                      *
      DATA K1/   1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 16, 17,
     &  18, 24, 30, 34, 38, 40, 41, 43, 44, 136, 138, 330, 327, 328,
     &   2*330, 46, 51, 52, 54, 55, 58,
     &  60, 62, 64, 66, 68, 70, 72, 74, 82, 90, 98, 106, 109, 112, 114,
     & 123, 140, 141, 143, 145, 146, 150, 157, 164, 168, 174, 180, 187,
     & 194, 202, 210, 211, 213, 215, 216, 220, 227, 234, 238, 245, 252,
     & 254, 255, 256, 257, 259, 262, 265, 267, 269, 272, 276, 279, 282,
     & 286, 290, 293, 299, 331, 335, 339, 340, 341, 343, 344, 345, 346,
     & 347, 350, 353, 356, 358, 360, 363, 366, 369, 372, 374, 376, 379,
     & 383, 385, 387, 391, 394, 397, 400, 402, 405, 408, 410, 412, 414,
     & 417, 420, 425, 430, 431, 432, 433, 434, 448, 452, 457, 458, 459,
     & 460, 461, 462, 466, 468, 470, 472, 486, 490, 495, 496, 497, 498,
     & 499, 500, 504, 506, 508, 510, 511, 512, 513, 514, 515, 516, 517,
     & 518, 519, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 534,
     & 537, 539, 541, 547, 553, 558, 563, 568, 572, 573, 574, 575, 576,
     & 577, 578, 579, 580, 581, 582, 583, 584, 585, 586, 587, 588, 589,
     & 590, 596, 602 /
*                                                                      *
*     Last number of decay channels used for resonances                *
*     and decaying particles                                           *
*                                                                      *
      DATA K2/   1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 15, 16, 17,
     & 23, 29, 31, 35, 39, 40, 42, 43, 45, 137, 139, 330, 327, 328,
     & 2* 330, 50, 51, 53, 54, 57,
     & 59, 61, 63, 65, 67, 69, 71, 73, 81, 89, 97, 105, 108, 111, 113,
     & 122, 135, 140, 142, 144, 145, 149, 156, 163, 167, 173, 179, 186,
     & 193, 201, 209, 210, 212, 214, 215, 219, 226, 233, 237, 244, 251,
     & 253, 254, 255, 256, 258, 261, 264, 266, 268, 271, 275, 278, 281,
     & 285, 289, 292, 298, 307, 334, 338, 339, 340, 342, 343, 344, 345,
     & 346, 349, 352, 355, 357, 359, 362, 365, 368, 371, 373, 375, 378,
     & 382, 384, 386, 390, 393, 396, 399, 401, 404, 407, 409, 411, 413,
     & 416, 419, 424, 429, 430, 431, 432, 433, 447, 451, 456, 457, 458,
     & 459, 460, 461, 465, 467, 469, 471, 485, 489, 494, 495, 496, 497,
     & 498, 499, 503, 505, 507, 509, 510, 511, 512, 513, 514, 515, 516,
     & 517, 518, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 533,
     & 536, 538, 540, 546, 552, 557, 562, 567, 571, 572, 573, 574, 575,
     & 576, 577, 578, 579, 580, 581, 582, 583, 584, 585, 586, 587, 588,
     & 589, 595, 601, 602 /
*                                                                      *
*
       END
*
 
*=== blkdt7 ===========================================================*
*==                                                                    *
      BLOCK DATA BLKD47
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
* * Block data 7 (ex 2)                                                *
C      INCLUDE '(DECAYC)'
C***********************************************************************
      PARAMETER (IDMAX9=602)
      CHARACTER*8 ZKNAME
      COMMON/DDECAC/ ZKNAME(IDMAX9),WT(IDMAX9),NZK(IDMAX9,3)
 
 
*                                                                      *
*     Name of decay channel                                            *
*                                                                      *
*                                                                      *
*   Designation N*@ means N*@1(1236)                                   *
*   @1=# means ++,  @1 = = means --                                    *
*   Designation  P+/0/- means Pi+/Pi0/Pi- , respectively               *
*                                                                      *
      DATA (ZKNAME(K),K=  1, 85) /
     &  'P       ','AP      ','E-      ','E+      ','NUE     ',
     &  'ANUE    ','GAM     ','PE-NUE  ','APEANU  ','EANUNU  ',
     &  'E-NUAN  ','3PI0    ','PI+-0   ','PIMUNU  ','PIE-NU  ',
     &  'MU+NUE  ','MU-NUE  ','MU+NUE  ','PI+PI0  ','PI++-   ',
     &  'PI+00   ','M+P0NU  ','E+P0NU  ','MU-NU   ','PI-0    ',
     &  'PI+--   ','PI-00   ','M-P0NU  ','E-P0NU  ','PPI-    ',
     &  'NPI0    ','PD-NUE  ','PM-NUE  ','APPI+   ','ANPI0   ',
     &  'APE+NU  ','APM+NU  ','PI+PI-  ','PI0PI0  ','NPI-    ',
     &  'PPI0    ','NPI+    ','LAGA    ','GAGA    ','GAE+E-  ',
     &  'GAGA    ','GAGAP0  ','PI000   ','PI+-0   ','PI+-GA  ',
     &  'PI+0    ','PI+-    ','PI00    ','PI-0    ','PI+-0   ',
     &  'PI+-    ','PI0GA   ','K+PI0   ','K0PI+   ','KOPI0   ',
     &  'K+PI-   ','K-PI0   ','AK0PI-  ','AK0PI0  ','K-PI+   ',
     &  'K+PI0   ','K0PI+   ','K0PI0   ','K+PI-   ','K-PI0   ',
     &  'K0PI-   ','AK0PI0  ','K-PI+   ','K+PI0   ','K0PI+   ',
     &  'K+89P0  ','K08PI+  ','K+RO77  ','K0RO+7  ','K+OM07  ',
     &  'K+E055  ','K0PI0   ','K+PI+   ','K089P0  ','K+8PI-  '  /
      DATA (ZKNAME(K),K= 86,170) /
     &  'K0R077  ','K+R-77  ','K+R-77  ','K0OM07  ','K0E055  ',
     &  'K-PI0   ','K0PI-   ','K-89P0  ','AK08P-  ','K-R077  ',
     &  'AK0R-7  ','K-OM07  ','K-E055  ','AK0PI0  ','K-PI+   ',
     &  'AK08P0  ','K-8PI+  ','AK0R07  ','AK0OM7  ','AK0E05  ',
     &  'LA0PI+  ','SI0PI+  ','SI+PI0  ','LA0PI0  ','SI+PI-  ',
     &  'SI-PI+  ','LA0PI-  ','SI0PI-  ','NEUAK0  ','PK-     ',
     &  'SI+PI-  ','SI0PI0  ','SI-PI+  ','LA0ET0  ','S+1PI-  ',
     &  'S-1PI+  ','SO1PI0  ','NEUAK0  ','PK-     ','LA0PI0  ',
     &  'LA0OM0  ','LA0RO0  ','SI+RO-  ','SI-RO+  ','SI0RO0  ',
     &  'LA0ET0  ','SI0ET0  ','SI+PI-  ','SI-PI+  ','SI0PI0  ',
     &  'K0S     ','K0L     ','K0S     ','K0L     ','P PI+   ',
     &  'P PI0   ','N PI+   ','P PI-   ','N PI0   ','N PI-   ',
     &  'P PI+   ','N*#PI0  ','N*+PI+  ','PRHO+   ','P PI0   ',
     &  'N PI+   ','N*#PI-  ','N*+PI0  ','N*0PI+  ','PRHO0   ',
     &  'NRHO+   ','P PI-   ','N PI0   ','N*+PI-  ','N*0PI0  ',
     &  'N*-PI+  ','PRHO-   ','NRHO0   ','N PI-   ','N*0PI-  ',
     &  'N*-PI0  ','NRHO-   ','PETA0   ','N*#PI-  ','N*+PI0  '  /
      DATA (ZKNAME(K),K=171,255) /
     &  'N*0PI+  ','PRHO0   ','NRHO+   ','NETA0   ','N*+PI-  ',
     &  'N*0PI0  ','N*-PI+  ','PRHO-   ','NRHO0   ','P PI0   ',
     &  'N PI+   ','N*#PI-  ','N*+PI0  ','N*0PI+  ','PRHO0   ',
     &  'NRHO+   ','P PI-   ','N PI0   ','N*+PI-  ','N*0PI0  ',
     &  'N*-PI+  ','PRHO-   ','NRHO0   ','P PI0   ','N PI+   ',
     &  'PRHO0   ','NRHO+   ','LAMK+   ','S+ K0   ','S0 K+   ',
     &  'PETA0   ','P PI-   ','N PI0   ','PRHO-   ','NRHO0   ',
     &  'LAMK0   ','S0 K0   ','S- K+   ','NETA/   ','APPI-   ',
     &  'APPI0   ','ANPI-   ','APPI+   ','ANPI0   ','ANPI+   ',
     &  'APPI-   ','AN*=P0  ','AN*-P-  ','APRHO-  ','APPI0   ',
     &  'ANPI-   ','AN*=P+  ','AN*-P0  ','AN*0P-  ','APRHO0  ',
     &  'ANRHO-  ','APPI+   ','ANPI0   ','AN*-P+  ','AN*0P0  ',
     &  'AN*+P-  ','APRHO+  ','ANRHO0  ','ANPI+   ','AN*0P+  ',
     &  'AN*+P0  ','ANRHO+  ','APPI0   ','ANPI-   ','AN*=P+  ',
     &  'AN*-P0  ','AN*0P-  ','APRHO0  ','ANRHO-  ','APPI+,  ',
     &  'ANPI0   ','AN*-P+  ','AN*0P0  ','AN*+P-  ','APRHO+  ',
     &  'ANRHO0  ','PN*014  ','NN*=14  ','PI+0    ','PI+-    '  /
      DATA (ZKNAME(K),K=256,340) /
     &  'PI-0    ','P+0     ','N++     ','P+-     ','P00     ',
     &  'N+0     ','N+-     ','N00     ','P-0     ','N-0     ',
     &  'P--     ','PPPI0   ','PNPI+   ','PNPI0   ','PPPI-   ',
     &  'NNPI+   ','APPPI0  ','APNPI+  ','ANNPI0  ','ANPPI-  ',
     &  'APNPI0  ','APPPI-  ','ANNPI-  ','K+PPI0  ','K+NPI+  ',
     &  'K0PPI0  ','K-PPI0  ','K-NPI+  ','AKPPI-  ','AKNPI0  ',
     &  'K+NPI0  ','K+PPI-  ','K0PPI0  ','K0NPI+  ','K-NPI0  ',
     &  'K-PPI-  ','AKNPI-  ','PAK0    ','SI+PI0  ','SI0PI+  ',
     &  'SI+ETA  ','S+1PI0  ','S01PI+  ','NEUK-   ','LA0PI-  ',
     &  'SI-OM0  ','LA0RO-  ','SI0RO-  ','SI-RO0  ','SI-ET0  ',
     &  'SI0PI-  ','SI-0    ','BLANC   ','BLANC   ','BLANC   ',
     &  'BLANC   ','BLANC   ','BLANC   ','BLANC   ','BLANC   ',
     &  'BLANC   ','BLANC   ','BLANC   ','BLANC   ','BLANC   ',
     &  'BLANC   ','BLANC   ','BLANC   ','BLANC   ','BLANC   ',
     &  'BLANC   ','BLANC   ','BLANC   ','BLANC   ','BLANC   ',
     &  'EPI+-   ','EPI00   ','GAPI+-  ','GAGA*   ','K+-     ',
     &  'KLKS    ','PI+-0   ','EGA     ','LPI0    ','LPI     '  /
      DATA (ZKNAME(K),K=341,425) /
     &  'APPI0   ','ANPI-   ','ALAGA   ','ANPI    ','ALPI0   ',
     &  'ALPI+   ','LAPI+   ','SI+PI0  ','SI0PI+  ','LAPI0   ',
     &  'SI+PI-  ','SI-PI+  ','LAPI-   ','SI-PI0  ','SI0PI-  ',
     &  'TE0PI0  ','TE-PI+  ','TE0PI-  ','TE-PI0  ','TE0PI   ',
     &  'TE-PI   ','LAK-    ','ALPI-   ','AS-PI0  ','AS0PI-  ',
     &  'ALPI0   ','AS+PI-  ','AS-PI+  ','ALPI+   ','AS+PI0  ',
     &  'AS0PI+  ','AT0PI0  ','AT+PI-  ','AT0PI+  ','AT+PI0  ',
     &  'AT0PI   ','AT+PI   ','ALK+    ','K-PI+   ','K-PI+0  ',
     &  'K0PI+-  ','K0PI0   ','K-PI++  ','AK0PI+  ','K+PI--  ',
     &  'K0PI-   ','K+PI-   ','K+PI-0  ','AKPI-+  ','AK0PI0  ',
     &  'ETAPIF  ','K++-    ','K+AK0   ','ETAPI-  ','K--+    ',
     &  'K-K0    ','PI00    ','PI+-    ','GAGA    ','D0PI0   ',
     &  'D0GA    ','D0PI+   ','D+PI0   ','DFGA    ','AD0PI-  ',
     &  'D-PI0   ','D-GA    ','AD0PI0  ','AD0GA   ','F+GA    ',
     &  'F+GA    ','F-GA    ','F-GA    ','PSPI+-  ','PSPI00  ',
     &  'PSETA   ','E+E-    ','MUE+-   ','PI+-0   ','M+NN    ',
     &  'E+NN    ','RHO+NT  ','PI+ANT  ','K*+ANT  ','M-NN    '  /
      DATA (ZKNAME(K),K=426,510) /
     &  'E-NN    ','RHO-NT  ','PI-NT   ','K*-NT   ','NUET    ',
     &  'ANUET   ','NUEM    ','ANUEM   ','SI+ETA  ','SI+ET*  ',
     &  'PAK0    ','TET0K+  ','SI*+ET  ','N*+AK0  ','N*++K-  ',
     &  'LAMRO+  ','SI0RO+  ','SI+RO0  ','SI+OME  ','PAK*0   ',
     &  'N*+AK*  ','N*++K*  ','SI+AK0  ','TET0PI  ','SI+AK*  ',
     &  'TET0RO  ','SI0AK*  ','SI+K*-  ','TET0OM  ','TET-RO  ',
     &  'SI*0AK  ','C0+PI+  ','C0+PI0  ','C0+PI-  ','A+GAM   ',
     &  'A0GAM   ','TET0AK  ','TET0K*  ','OM-RO+  ','OM-PI+  ',
     &  'C1++AK  ','A+PI+   ','C0+AK0  ','A0PI+   ','A+AK0   ',
     &  'T0PI+   ','ASI-ET  ','ASI-E*  ','APK0    ','ATET0K  ',
     &  'ASI*-E  ','AN*-K0  ','AN*--K  ','ALAMRO  ','ASI0RO  ',
     &  'ASI-RO  ','ASI-OM  ','APK*0   ','AN*-K*  ','AN*--K  ',
     &  'ASI-K0  ','ATETPI  ','ASI-K*  ','ATETRO  ','ASI0K*  ',
     &  'ASI-K*  ','ATE0OM  ','ATE+RO  ','ASI*0K  ','AC-PI-  ',
     &  'AC-PI0  ','AC-PI+  ','AA-GAM  ','AA0GAM  ','ATET0K  ',
     &  'ATE0K*  ','AOM+RO  ','AOM+PI  ','AC1--K  ','AA-PI-  ',
     &  'AC0-K0  ','AA0PI-  ','AA-K0   ','AT0PI-  ','C1++GA  '  /
      DATA (ZKNAME(K),K=511,540) /
     &  'C1++GA  ','C10GAM  ','S+GAM   ','S0GAM   ','T0GAM   ',
     &  'XU++GA  ','XD+GAM  ','XS+GAM  ','A+AKPI  ','T02PI+  ',
     &  'C1++2K  ','AC1--G  ','AC1-GA  ','AC10GA  ','AS-GAM  ',
     &  'AS0GAM  ','AT0GAM  ','AXU--G  ','AXD-GA  ','AXS-GA  ',
     &  'AA-KPI  ','AT02PI  ','AC1--K  ','RH-PI+  ','RH+PI-  ',
     &  'RH3PI0  ','RH0PI+  ','RH+PI0  ','RH0PI-  ','RH-PI0  '  /
      DATA (ZKNAME(I),I=541,602)/
     & 'APETA ','AN=P+ ','AN-PO ','ANOPO ','APRHO0','ANRHO-','ANETA ',
     & 'AN-P+ ','AN0PO ','AN+P- ','APRHO+','ANRHO0','RH0PI+','RH+PI0',
     & '3PI+00','3PI-++','F0PI+ ','RH+PI-','RH0PI0','3PI000','3PI0+-',
     & 'F0PI0 ','RH0PI-','RH-PI0','3PI-00','3PI--+','F0PI- ','PI0PI0',
     & 'PI+PI-','K+K-  ','K0AK0 ','L01600','AL0160','K*+146','K*-146',
     & 'K*0146','AK0146','S+1660','S01660','S-1660','AS-166','AS0166',
     & 'AS+166','X01690','X-1690','AX0169','AX+169','OM-225','AOM+22',
     & 'N*PPI0','N*NPI+','N*P2P0','N*PP+-','N*D+P0','N*D0P+','N*NPI0',
     & 'N*PPI-','N*N2P0','N*NP+-','N*D+P-','N*D0P0','BLANK '/
*                                                                      *
*     Weight of decay channel                                          *
*                                                                      *
      DATA (WT(K),K=  1, 85) /
     &   .1000D+01, .1000D+01, .1000D+01, .1000D+01, .1000D+01,
     &   .1000D+01, .1000D+01, .1000D+01, .1000D+01, .1000D+01,
     &   .1000D+01, .2100D+00, .1200D+00, .2700D+00, .4000D+00,
     &   .1000D+01, .1000D+01, .6400D+00, .2100D+00, .6000D-01,
     &   .2000D-01, .3000D-01, .4000D-01, .6400D+00, .2100D+00,
     &   .6000D-01, .2000D-01, .3000D-01, .4000D-01, .6400D+00,
     &   .3600D+00, .0000D+00, .0000D+00, .6400D+00, .3600D+00,
     &   .0000D+00, .0000D+00, .6900D+00, .3100D+00, .1000D+01,
     &   .5200D+00, .4800D+00, .1000D+01, .9900D+00, .1000D-01,
     &   .3800D+00, .3000D-01, .3000D+00, .2400D+00, .5000D-01,
     &   .1000D+01, .1000D+01, .0000D+00, .1000D+01, .9000D+00,
     &   .1000D-01, .9000D-01, .3300D+00, .6700D+00, .3300D+00,
     &   .6700D+00, .3300D+00, .6700D+00, .3300D+00, .6700D+00,
     &   .3300D+00, .6700D+00, .3300D+00, .6700D+00, .3300D+00,
     &   .6700D+00, .3300D+00, .6700D+00, .1900D+00, .3800D+00,
     &   .9000D-01, .2000D+00, .3000D-01, .4000D-01, .5000D-01,
     &   .2000D-01, .1900D+00, .3800D+00, .9000D-01, .2000D+00  /
      DATA (WT(K),K= 86,170) /
     &   .3000D-01, .4000D-01, .5000D-01, .2000D-01, .1900D+00,
     &   .3800D+00, .9000D-01, .2000D+00, .3000D-01, .4000D-01,
     &   .5000D-01, .2000D-01, .1900D+00, .3800D+00, .9000D-01,
     &   .2000D+00, .3000D-01, .4000D-01, .5000D-01, .2000D-01,
     &   .8800D+00, .6000D-01, .6000D-01, .8800D+00, .6000D-01,
     &   .6000D-01, .8800D+00, .1200D+00, .1900D+00, .1900D+00,
     &   .1600D+00, .1600D+00, .1700D+00, .3000D-01, .3000D-01,
     &   .3000D-01, .4000D-01, .1000D+00, .1000D+00, .2000D+00,
     &   .1200D+00, .1000D+00, .4000D-01, .4000D-01, .5000D-01,
     &   .7500D-01, .7500D-01, .3000D-01, .3000D-01, .4000D-01,
     &   .5000D+00, .5000D+00, .5000D+00, .5000D+00, .1000D+01,
     &   .6700D+00, .3300D+00, .3300D+00, .6700D+00, .1000D+01,
     &   .2500D+00, .2700D+00, .1800D+00, .3000D+00, .1700D+00,
     &   .8000D-01, .1800D+00, .3000D-01, .2400D+00, .2000D+00,
     &   .1000D+00, .8000D-01, .1700D+00, .2400D+00, .3000D-01,
     &   .1800D+00, .1000D+00, .2000D+00, .2500D+00, .1800D+00,
     &   .2700D+00, .3000D+00, .5000D+00, .3000D+00, .1250D+00  /
C    &   .2700D+00, .3000D+00, .4000D+00, .2000D+00, .1250D+00  /
C    &   .7500D-01, .7500D-01, .1250D+00, .4000D+00, .7500D-01,
C    &   .1250D+00, .2000D+00, .1250D+00, .7500D-01, .1800D+00,
      DATA (WT(K),K=171,255) /
     &   .7500D-01, .0000D+00, .0000D+00, .5000D+00, .7500D-01,
     &   .1250D+00, .3000D+00, .0000D+00, .0000D+00, .1800D+00,
     &   .3700D+00, .1300D+00, .8000D-01, .4000D-01, .7000D-01,
     &   .1300D+00, .3700D+00, .1800D+00, .4000D-01, .8000D-01,
     &   .1300D+00, .1300D+00, .7000D-01, .7000D-01, .1300D+00,
     &   .2300D+00, .4700D+00, .5000D-01, .2000D-01, .1000D-01,
     &   .2000D-01, .1300D+00, .7000D-01, .4700D+00, .2300D+00,
     &   .5000D-01, .1000D-01, .2000D-01, .2000D-01, .1000D+01,
     &   .6700D+00, .3300D+00, .3300D+00, .6700D+00, .1000D+01,
     &   .2500D+00, .2700D+00, .1800D+00, .3000D+00, .1700D+00,
     &   .8000D-01, .1800D+00, .3000D-01, .2400D+00, .2000D+00,
     &   .1000D+00, .8000D-01, .1700D+00, .2400D+00, .3000D-01,
     &   .1800D+00, .1000D+00, .2000D+00, .2500D+00, .1800D+00,
     &   .2700D+00, .3000D+00, .1800D+00, .3700D+00, .1300D+00,
     &   .8000D-01, .4000D-01, .7000D-01, .1300D+00, .3700D+00,
     &   .1800D+00, .4000D-01, .8000D-01, .1300D+00, .1300D+00,
     &   .7000D-01, .5000D+00, .5000D+00, .1000D+01, .1000D+01  /
      DATA (WT(K),K=256,340) /
     &   .1000D+01, .8000D+00, .2000D+00, .6000D+00, .3000D+00,
     &   .1000D+00, .6000D+00, .3000D+00, .1000D+00, .8000D+00,
     &   .2000D+00, .3300D+00, .6700D+00, .6600D+00, .1700D+00,
     &   .1700D+00, .3200D+00, .1700D+00, .3200D+00, .1900D+00,
     &   .3300D+00, .3300D+00, .3400D+00, .3000D+00, .5000D-01,
     &   .6500D+00, .3800D+00, .1200D+00, .3800D+00, .1200D+00,
     &   .3800D+00, .1200D+00, .3800D+00, .1200D+00, .3000D+00,
     &   .5000D-01, .6500D+00, .3800D+00, .2500D+00, .2500D+00,
     &   .2000D-01, .5000D-01, .5000D-01, .2000D+00, .2000D+00,
     &   .1200D+00, .1000D+00, .7000D-01, .7000D-01, .1400D+00,
     &   .5000D-01, .5000D-01, .1000D+01, .1000D+01, .1000D+01,
     &   .1000D+01, .1000D+01, .1000D+01, .1000D+01, .1000D+01,
     &   .1000D+01, .1000D+01, .1000D+01, .1000D+01, .1000D+01,
     &   .1000D+01, .1000D+01, .1000D+01, .1000D+01, .1000D+01,
     &   .1000D+01, .1000D+01, .1000D+01, .1000D+01, .1000D+01,
     &   .4800D+00, .2400D+00, .2600D+00, .2000D-01, .4700D+00,
     &   .3500D+00, .1500D+00, .3000D-01, .1000D+01, .1000D+01  /
      DATA (WT(K),K=341,425) /
     &   .5200D+00, .4800D+00, .1000D+01, .1000D+01, .1000D+01,
     &   .1000D+01, .9000D+00, .5000D-01, .5000D-01, .9000D+00,
     &   .5000D-01, .5000D-01, .9000D+00, .5000D-01, .5000D-01,
     &   .3300D+00, .6700D+00, .6700D+00, .3300D+00, .2500D+00,
     &   .2500D+00, .5000D+00, .9000D+00, .5000D-01, .5000D-01,
     &   .9000D+00, .5000D-01, .5000D-01, .9000D+00, .5000D-01,
     &   .5000D-01, .3300D+00, .6700D+00, .6700D+00, .3300D+00,
     &   .2500D+00, .2500D+00, .5000D+00, .1000D+00, .5000D+00,
     &   .1600D+00, .2400D+00, .7000D+00, .3000D+00, .7000D+00,
     &   .3000D+00, .1000D+00, .5000D+00, .1600D+00, .2400D+00,
     &   .3000D+00, .4000D+00, .3000D+00, .3000D+00, .4000D+00,
     &   .3000D+00, .4900D+00, .4900D+00, .2000D-01, .5500D+00,
     &   .4500D+00, .6800D+00, .3000D+00, .2000D-01, .6800D+00,
     &   .3000D+00, .2000D-01, .5500D+00, .4500D+00, .9000D+00,
     &   .1000D+00, .9000D+00, .1000D+00, .6000D+00, .3000D+00,
     &   .1000D+00, .1000D+00, .1000D+00, .8000D+00, .2800D+00,
     &   .2800D+00, .3500D+00, .7000D-01, .2000D-01, .2800D+00  /
      DATA (WT(K),K=426,510) /
     &   .2800D+00, .3500D+00, .7000D-01, .2000D-01, .1000D+01,
     &   .1000D+01, .1000D+01, .1000D+01, .2000D-01, .3000D-01,
     &   .7000D-01, .2000D-01, .2000D-01, .4000D-01, .1300D+00,
     &   .7000D-01, .6000D-01, .6000D-01, .2000D+00, .1400D+00,
     &   .4000D-01, .1000D+00, .2500D+00, .3000D-01, .3000D+00,
     &   .4200D+00, .2200D+00, .3500D+00, .1900D+00, .1600D+00,
     &   .8000D-01, .1000D+01, .1000D+01, .1000D+01, .1000D+01,
     &   .1000D+01, .3700D+00, .2000D+00, .3600D+00, .7000D-01,
     &   .5000D+00, .5000D+00, .5000D+00, .5000D+00, .5000D+00,
     &   .5000D+00, .2000D-01, .3000D-01, .7000D-01, .2000D-01,
     &   .2000D-01, .4000D-01, .1300D+00, .7000D-01, .6000D-01,
     &   .6000D-01, .2000D+00, .1400D+00, .4000D-01, .1000D+00,
     &   .2500D+00, .3000D-01, .3000D+00, .4200D+00, .2200D+00,
     &   .3500D+00, .1900D+00, .1600D+00, .8000D-01, .1000D+01,
     &   .1000D+01, .1000D+01, .1000D+01, .1000D+01, .3700D+00,
     &   .2000D+00, .3600D+00, .7000D-01, .5000D+00, .5000D+00,
     &   .5000D+00, .5000D+00, .5000D+00, .5000D+00, .1000D+01  /
      DATA (WT(K),K=511,540) /
     &   .1000D+01, .1000D+01, .1000D+01, .1000D+01, .1000D+01,
     &   .1000D+01, .1000D+01, .1000D+01, .3000D+00, .3000D+00,
     &   .4000D+00, .1000D+01, .1000D+01, .1000D+01, .1000D+01,
     &   .1000D+01, .1000D+01, .1000D+01, .1000D+01, .1000D+01,
     &   .3000D+00, .3000D+00, .4000D+00, .3300D+00, .3300D+00,
     &   .3400D+00, .5000D+00, .5000D+00, .5000D+00, .5000D+00  /
C
      DATA (WT(I),I=541,602) / .0D+00, .3334D+00, .2083D+00, 2*.125D+00,
     & .2083D+00, .0D+00, .125D+00, .2083D+00, .3334D+00, .2083D+00,
     & .125D+00,  0.2D+00, 0.2D+00, 0.3D+00, 0.3D+00, 0.0D+00, 0.2D+00,
     & 0.2D+00, 0.3D+00, 0.3D+00, 0.0D+00, 0.2D+00, 0.2D+00, 0.3D+00,
     & 0.3D+00, 0.0D+00, 0.31D+00, 0.62D+00, 0.035D+00, 0.035D+00,
     & 18*1.D+00, 0.5D+00, 0.16D+00, 2*0.12D+00, 2*0.05D+00, 0.5D+00,
     & 0.16D+00, 2*0.12D+00, 2*0.05D+00, 1.D+00 /
*
*     Particle numbers in decay channel                                *
*                                                                      *
      DATA (NZK(K,1),K=  1,170) /
     &     1,   2,   3,   4,   5,   6,   7,   1,   2,   4,
     &     3,  23,  13,  13,  13,  10,  11,  10,  13,  13,
     &    13,  10,   4,  11,  14,  14,  14,  11,   3,   1,
     &     8,   1,   1,   2,   9,   2,   2,  13,  23,   8,
     &     1,   8,  17,   7,   7,   7,  23,  23,  13,  13,
     &    13,  13,  23,  14,  13,  13,  23,  15,  24,  24,
     &    15,  16,  25,  25,  16,  15,  24,  24,  15,  16,
     &    24,  25,  16,  15,  24,  36,  37,  15,  24,  15,
     &    15,  24,  15,  37,  36,  24,  15,  24,  24,  16,
     &    24,  38,  39,  16,  25,  16,  16,  25,  16,  39,
     &    38,  25,  16,  25,  25,  17,  22,  21,  17,  21,
     &    20,  17,  22,   8,   1,  21,  22,  20,  17,  48,
     &    50,  49,   8,   1,  17,  17,  17,  21,  20,  22,
     &    17,  22,  21,  20,  22,  19,  12,  19,  12,   1,
     &     1,   8,   1,   8,   8,   1,  53,  54,   1,   1,
     &     8,  53,  54,  55,   1,   8,   1,   8,  54,  55,
     &    56,   1,   8,   8,  55,  56,   8,   1,  53,  54  /
      DATA (NZK(K,1),K=171,340) /
     &    55,   1,   8,   8,  54,  55,  56,   1,   8,   1,
     &     8,  53,  54,  55,   1,   8,   1,   8,  54,  55,
     &    56,   1,   8,   1,   8,   1,   8,  17,  21,  22,
     &     1,   1,   8,   1,   8,  17,  22,  20,   8,   2,
     &     2,   9,   2,   9,   9,   2,  67,  68,   2,   2,
     &     9,  67,  68,  69,   2,   9,   2,   9,  68,  69,
     &    70,   2,   9,   9,  69,  70,   9,   2,   9,  67,
     &    68,  69,   2,   9,   2,   9,  68,  69,  70,   2,
     &     9,   1,   8,  13,  13,  14,   1,   8,   1,   1,
     &     8,   8,   8,   1,   8,   1,   1,   1,   1,   1,
     &     8,   2,   2,   9,   9,   2,   2,   9,  15,  15,
     &    24,  16,  16,  25,  25,  15,  15,  24,  24,  16,
     &    16,  25,   1,  21,  22,  21,  48,  49,   8,  17,
     &    20,  17,  22,  20,  20,  22,  20,   0,   0,   0,
     &     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     &     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     &    31,  31,  13,   7,  15,  12,  13,  31,  17,  17  /
      DATA (NZK(K,1),K=341,510) /
     &     2,   9,  18,   9,  18,  18,  17,  21,  22,  17,
     &    21,  20,  17,  20,  22,  97,  98,  97,  98,  97,
     &    98,  17,  18,  99, 100,  18, 101,  99,  18, 101,
     &   100, 102, 103, 102, 103, 102, 103,  18,  16,  16,
     &    24,  24,  16,  25,  15,  24,  15,  15,  25,  25,
     &    31,  15,  15,  31,  16,  16,  23,  13,   7, 116,
     &   116, 116, 117, 117, 119, 118, 118, 119, 119, 120,
     &   120, 121, 121, 130, 130, 130,   4,  10,  13,  10,
     &     4,  32,  13,  36,  11,   3,  34,  14,  38, 133,
     &   134, 135, 136,  21,  21,   1,  97, 104,  54,  53,
     &    17,  22,  21,  21,   1,  54,  53,  21,  97,  21,
     &    97,  22,  21,  97,  98, 105, 137, 137, 137, 138,
     &   139,  97,  97, 109, 109, 140, 138, 137, 139, 138,
     &   145,  99,  99,   2, 102, 110,  68,  67,  18, 100,
     &    99,  99,   2,  68,  67,  99, 102,  99, 102, 100,
     &    99, 102, 103, 111, 149, 149, 149, 150, 151, 113,
     &   113, 115, 115, 152, 150, 149, 151, 150, 157, 140  /
      DATA (NZK(K,1),K=511,540) /
     &   141, 142, 143, 144, 145, 146, 147, 148, 138, 145,
     &   140, 152, 153, 154, 155, 156, 157, 158, 159, 160,
     &   150, 157, 152,  34,  32,  33,  33,  32,  33,  34  /
      DATA (NZK(I,1),I=541,602) /  2, 67, 68, 69,  2,  9,  9, 68, 69,
     & 70,  2,  9, 33, 32, 13, 14, 189, 32, 34, 23, 23, 189, 33, 34, 14,
     & 14, 189, 23, 13, 15, 24,  36,  38,  37,  39, 194, 195, 196, 197,
     & 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 1, 8, 1, 1, 54,
     & 55, 8, 1, 8, 8, 54, 55, 210/
      DATA (NZK(K,2),K=  1,170) /
     &     0,   0,   0,   0,   0,   0,   0,   3,   4,   6,
     &     5,  23,  14,  11,   3,   5,   5,   5,  23,  13,
     &    23,  23,  23,   5,  23,  13,  23,  23,  23,  14,
     &    23,   3,  11,  13,  23,   4,  10,  14,  23,  14,
     &    23,  13,   7,   7,   4,   7,   7,  23,  14,  14,
     &    23,  14,  23,  23,  14,  14,   7,  23,  13,  23,
     &    14,  23,  14,  23,  13,  23,  13,  23,  14,  23,
     &    14,  23,  13,  23,  13,  23,  13,  33,  32,  35,
     &    31,  23,  14,  23,  14,  33,  34,  35,  31,  23,
     &    14,  23,  14,  33,  34,  35,  31,  23,  13,  23,
     &    13,  33,  32,  35,  31,  13,  13,  23,  23,  14,
     &    13,  14,  14,  25,  16,  14,  23,  13,  31,  14,
     &    13,  23,  25,  16,  23,  35,  33,  34,  32,  33,
     &    31,  31,  14,  13,  23,   0,   0,   0,   0,  13,
     &    23,  13,  14,  23,  14,  13,  23,  13,  78,  23,
     &    13,  14,  23,  13,  79,  78,  14,  23,  14,  23,
     &    13,  80,  79,  14,  14,  23,  80,  31,  14,  23  /
      DATA (NZK(K,2),K=171,340) /
     &    13,  79,  78,  31,  14,  23,  13,  80,  79,  23,
     &    13,  14,  23,  13,  79,  78,  14,  23,  14,  23,
     &    13,  80,  79,  23,  13,  33,  32,  15,  24,  15,
     &    31,  14,  23,  34,  33,  24,  24,  15,  31,  14,
     &    23,  14,  13,  23,  13,  14,  23,  14,  80,  23,
     &    14,  13,  23,  14,  79,  80,  13,  23,  13,  23,
     &    14,  78,  79,  13,  13,  23,  78,  23,  14,  13,
     &    23,  14,  79,  80,  13,  23,  13,  23,  14,  78,
     &    79,  62,  61,  23,  14,  23,  13,  13,  13,  23,
     &    13,  13,  23,  14,  14,  14,   1,   8,   8,   1,
     &     8,   1,   8,   8,   1,   8,   1,   8,   1,   8,
     &     1,   1,   8,   1,   8,   8,   1,   1,   8,   8,
     &     1,   8,  25,  23,  13,  31,  23,  13,  16,  14,
     &    35,  34,  34,  33,  31,  14,  23,   0,   0,   0,
     &     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     &     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     &    13,  23,  14,   7,  16,  19,  14,   7,  23,  14  /
      DATA (NZK(K,2),K=341,510) /
     &    23,  14,   7,  13,  23,  13,  13,  23,  13,  23,
     &    14,  13,  14,  23,  14,  23,  13,  14,  23,  14,
     &    23,  16,  14,  23,  14,  23,  14,  13,  13,  23,
     &    13,  23,  14,  13,  23,  13,  23,  15,  13,  13,
     &    13,  23,  13,  13,  14,  14,  14,  14,  14,  23,
     &    13,  16,  25,  14,  15,  24,  23,  14,   7,  23,
     &     7,  13,  23,   7,  14,  23,   7,  23,   7,   7,
     &     7,   7,   7,  13,  23,  31,   3,  11,  14, 135,
     &     5, 134, 134, 134, 136,   6, 133, 133, 133,   0,
     &     0,   0,   0,  31,  95,  25,  15,  31,  95,  16,
     &    32,  32,  33,  35,  39,  39,  38,  25,  13,  39,
     &    32,  39,  38,  35,  32,  39,  13,  23,  14,   7,
     &     7,  25,  37,  32,  13,  25,  13,  25,  13,  25,
     &    13,  31,  95,  24,  16,  31,  24,  15,  34,  34,
     &    33,  35,  37,  37,  36,  24,  14,  37,  34,  37,
     &    36,  35,  34,  37,  14,  23,  13,   7,   7,  24,
     &    39,  34,  14,  24,  14,  24,  14,  24,  14,   7  /
      DATA (NZK(K,2),K=511,540) /
     &     7,   7,   7,   7,   7,   7,   7,   7,  25,  13,
     &    25,   7,   7,   7,   7,   7,   7,   7,   7,   7,
     &    24,  14,  24,  13,  14,  23,  13,  23,  14,  23  /
      DATA (NZK(I,2),I=541,602) / 31, 13, 23, 14, 79, 80, 31, 13, 23,
     & 14, 78, 79, 13, 23, 23, 13, 13, 14, 13, 23, 13, 23, 14, 23, 23,
     & 14, 14, 23, 14, 16, 25,
     & 4*23, 14*0, 23, 13, 23, 13, 23, 13, 23, 14,
     & 23, 13, 14, 23,  0 /
      DATA (NZK(K,3),K=  1,170) /
     &     0,   0,   0,   0,   0,   0,   0,   5,   6,   5,
     &     6,  23,  23,   5,   5,   0,   0,   0,   0,  14,
     &    23,   5,   5,   0,   0,  14,  23,   5,   5,   0,
     &     0,   5,   5,   0,   0,   5,   5,   0,   0,   0,
     &     0,   0,   0,   0,   3,   0,   7,  23,  23,   7,
     &     0,   0,   0,   0,  23,   0,   0,   0,   0,   0,
     &     110*0   /
      DATA (NZK(K,3),K=171,340) /
     &     80*0,
     &     0,   0,   0,   0,   0,   0,  23,  13,  14,  23,
     &    23,  14,  23,  23,  23,  14,  23,  13,  23,  14,
     &    13,  23,  13,  23,  14,  23,  14,  14,  23,  13,
     &    13,  23,  13,  14,  23,  23,  14,  23,  13,  23,
     &    14,  14,   0,   0,   0,   0,   0,   0,   0,   0,
     &     30*0,
     &    14,  23,   7,   0,   0,   0,  23,   0,   0,   0  /
      DATA (NZK(K,3),K=341,510) /
     &     30*0,
     &     0,   0,   0,   0,   0,   0,   0,   0,   0,  23,
     &    14,   0,  13,   0,  14,   0,   0,  23,  13,   0,
     &     0,  15,   0,   0,  16,   0,   0,   0,   0,   0,
     &     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     &     0,   0,   0,  14,  23,   0,   0,   0,  23, 134,
     &   134,   0,   0,   0, 133, 133,   0,   0,   0,   0,
     &     80*0  /
      DATA (NZK(K,3),K=511,540) /
     &     0,   0,   0,   0,   0,   0,   0,   0,  13,  13,
     &    25,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     &    14,  14,  24,   0,   0,   0,   0,   0,   0,   0  /
      DATA (NZK(I,3),I=541,602) / 12*0, 2*0, 23, 13, 0, 2*0, 23, 14, 0,
     & 2*0, 23, 13, 0, 4*0, 18*0, 2*0, 23, 14, 2*0, 2*0, 23, 14, 2*0, 0/
*=                                               end.block.blkdt7      *
*                                                                      *
      END
 
*
*===xsglau=============================================================*
*
      SUBROUTINE XSGLAU(NA,NB,IJPROJ,NTARG)

************************************************************************
* Total, elastic, quasi-elastic, inelastic cross sections according to *
* Glauber's approach.                                                  *
*  NA / NB     mass numbers of proj./target nuclei                     *
*  IJPROJ      bamjet-index of projectile (=1 in case of proj.nucleus) *
*  ECMI kinematical variables   E_cm                      *
*  IE       indices of energy 
*  NTARG       index of target nucleus set o NTARG=1 here             *
* This version dated 17.3.98  is written by S. Roesler mod by J.R.    *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)

      COMPLEX*16 CZERO,CONE,CTWO
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,
     &           ONETHI=ONE/THREE,TINY25=1.0D-25)
      PARAMETER (TWOPI  = 6.283185307179586454D+00,
     &           PI     = TWOPI/TWO,
     &           GEV2MB = 0.38938D0,
     &           GEV2FM = 0.1972D0,
     &           ALPHEM = ONE/137.0D0,
* proton mass
     &           AMP    = 0.938D0,
     &           AMP2   = AMP**2,
* approx. nucleon radius
     &           RNUCLE = 1.12D0,
* number of bins in b-space
     &           KSITEB = 200    )

      CHARACTER*8  ANAME
      COMMON /DPAR/  ANAME(210),AAM(210),GAM(210),TAU(210),IICH(210),
     &               IIBAR(210),KA1(210),KA2(210)

      PARAMETER (NCOMPX=1,NEB=50) 
      COMMON /DSHMM/   RASH,RBSH(NCOMPX),BMAX(NCOMPX),BSTEP(NCOMPX),
     &                SIGSH,ROSH,GSH,BSITE(0:NEB,NCOMPX,KSITEB),
     &                NSTATB,NSITEB
      COMMON /GLABER/ ECMNN(NEB),ECMNOW,
     &                XSTOT(NEB),XSELA(NEB),
     &                XSQEP(NEB),XSQET(NEB),
     &                XSQE2(NEB),XSPRO(NEB),
     &                XETOT(NEB),XEELA(NEB),
     &                XEQEP(NEB),XEQET(NEB),
     &                XEQE2(NEB),XEPRO(NEB),
     &                BSLOPE,ELABB(NEB)

      COMMON /VDMPAR/ RL2,EPSPOL,INTRGE(2),IDPDF,MODEGA,ISHAD(3)
      COMMON /GLAPAR/ JSTATB

      COMPLEX*16 C,CA,CI
      COMMON /DAMP/   CA,CI,GA
      COMMON /XSECNU/ECMUU,ECMOO,NGRITT,NEVTT
      COMMON /KGLAUB/JGLAUB

      PARAMETER (MAXNCL = 210)
      COMPLEX*16 PP11,PP12,PP21,PP22,
     &           OMPP11,OMPP12,OMPP21,OMPP22
      DIMENSION COOP1(3,MAXNCL),COOT1(3,MAXNCL),
     &          COOP2(3,MAXNCL),COOT2(3,MAXNCL),
     &          BPROD(KSITEB),SIGSHH(NEB),
     &          SIGTO(NEB),SIGEL(NEB),SIGIN(NEB),SIGSD(NEB),SIGDIF(NEB)

      JGLAUB=1
      Write(6,*)' XSGLAU(NA,NB,IJPROJ,NTARG)',
     &NA,NB,IJPROJ,NTARG
      WRITE(6,*)'/XSECNU/ECMUU,ECMOO,NGRITT,NEVTT',
     &ECMUU,ECMOO,NGRITT,NEVTT 

      CZERO  = DCMPLX(ZERO,ZERO)
      CONE   = DCMPLX(ONE,ZERO)
      CTWO   = DCMPLX(TWO,ZERO)

* re-define kinematics
      EC000=ECMUU
      DELLOG=(LOG10(ECMOO)-LOG10(ECMUU))/(NGRITT-1)
      DELDEL=10.D0**DELLOG
      EC111=ECMUU/DELDEL
      DO 1123 IEEE=1,NGRITT
      IE=IEEE
      EC111=DELDEL*EC111
      S=EC111**2
      ECMNN(IE) = EC111
      WRITE(6,*)'IE,EC111,S',IE,EC111,S

* parameters determining statistics in evaluating Glauber-xsection
      JSTATB=NEVTT
      NSTATB = JSTATB
      NSITEB = KSITEB

* set up interaction geometry (common /DSHM/)
*  projectile/target radii
         RASH = RNUCLE*DBLE(NA)**ONETHI
      RBSH(NTARG) = RNUCLE*DBLE(NB)**ONETHI
      IF(JGLAUB.EQ.1)THEN
        IF(NA.EQ.9)RASH=2.52D0
        IF(NA.EQ.10)RASH=2.45D0
        IF(NA.EQ.11)RASH=2.37D0
        IF(NA.EQ.12)RASH=2.45D0
        IF(NA.EQ.13)RASH=2.44D0
        IF(NA.EQ.14)RASH=2.55D0
        IF(NA.EQ.15)RASH=2.58D0
        IF(NA.EQ.16)RASH=2.71D0
        IF(NA.EQ.17)RASH=2.66D0
        IF(NA.EQ.18)RASH=2.71D0
        IF(NB.EQ.9)RBSH(NTARG)=2.52D0
        IF(NB.EQ.10)RBSH(NTARG)=2.45D0
        IF(NB.EQ.11)RBSH(NTARG)=2.37D0
        IF(NB.EQ.12)RBSH(NTARG)=2.45D0
        IF(NB.EQ.13)RBSH(NTARG)=2.44D0
        IF(NB.EQ.14)RBSH(NTARG)=2.55D0
        IF(NB.EQ.15)RBSH(NTARG)=2.58D0
	IF(NB.EQ.16)RBSH(NTARG)=2.71D0
        IF(NB.EQ.17)RBSH(NTARG)=2.66D0
        IF(NB.EQ.18)RBSH(NTARG)=2.71D0
       ENDIF
*  maximum impact-parameter
      BMAX(NTARG) = 4.0D0*(RASH+RBSH(NTARG))
      BSTEP(NTARG)= BMAX(NTARG)/DBLE(NSITEB-1)

* slope, rho ( Re(f(0))/Im(f(0)) )
      IF (IJPROJ.LE.12) THEN
            BSLOPE = 8.5D0*(1.0D0+0.065D0*LOG(S))
         IF (ECMNN(IE).LE.3.0D0) THEN
            ROSH = -0.43D0
         ELSEIF ((ECMNN(IE).GT.3.0D0).AND.(ECMNN(IE).LE.50.D0)) THEN
            ROSH = -0.63D0+0.175D0*LOG(ECMNN(IE))
         ELSEIF (ECMNN(IE).GT.50.0D0) THEN
            ROSH = 0.1D0
         ENDIF
      ELSE
         BSLOPE = 6.0D0*(1.0D0+0.065D0*LOG(S))
         ROSH   = 0.01D0
      ENDIF

* projectile-nucleon xsection (in fm)
         ELAB  = (S-AAM(IJPROJ)**2-AMP2)/(TWO*AMP)
	 ELABB(IE)=ELAB/1000.
         PLAB  = SQRT( (ELAB-AAM(IJPROJ))*(ELAB+AAM(IJPROJ)) )
C        SIGSH = SHNTOT(IJPROJ,1,ZERO,PLAB)/10.0D0
         SIGSH = DSHPTO(IJPROJ,PLAB)/10.D0
	 SIGSHH(IE)=SIGSH*10.D0
      WRITE(6,*)' NSTATB,NSITEB,RASH,RBSH(NTARG),BMAX(NTARG),
     &BSLOPE,ROSH,SIGSH,ECM ELAB',
     & NSTATB,NSITEB,RASH,RBSH(NTARG),BMAX(NTARG),
     &BSLOPE,ROSH,SIGSH,EC111,ELAB 
* initializations
      DO 10 I=1,NSITEB
         BSITE( 0,NTARG,I) = ZERO
         BSITE(IE,NTARG,I) = ZERO
         BPROD(I) = ZERO
   10 CONTINUE
      STOT  = ZERO
      STOT2 = ZERO
      SELA  = ZERO
      SELA2 = ZERO
      SQEP  = ZERO
      SQEP2 = ZERO
      SQET  = ZERO
      SQET2 = ZERO
      SQE2  = ZERO
      SQE22 = ZERO
      SPRO  = ZERO
      SPRO2 = ZERO
      FACN   = ONE/DBLE(NSTATB)

      IPNT = 0
      RPNT = ZERO

C------------------------------------------------------

* cross sections averaged over NSTATB nucleon configurations
      DO 11 IS=1,NSTATB
         STOTN = ZERO
         SELAN = ZERO
         SQEPN = ZERO
         SQETN = ZERO
         SQE2N = ZERO
         SPRON = ZERO
         CALL CONUCLX(COOP1,NA,RASH,0)
         CALL CONUCLX(COOT1,NB,RBSH(NTARG),1)
         CALL CONUCLX(COOP2,NA,RASH,0)
         CALL CONUCLX(COOT2,NB,RBSH(NTARG),1)

*  integration over impact parameter B
         DO 12 IB=1,NSITEB-1
            STOTB = ZERO
            SELAB = ZERO
            SQEPB = ZERO
            SQETB = ZERO
            SQE2B = ZERO
            SPROB = ZERO
            SDIR  = ZERO
            B     = DBLE(IB)*BSTEP(NTARG)
            FACB  = 10.0D0*TWOPI*B*BSTEP(NTARG)

*   integration over M_V^2 for photon-proj.
C           DO 14 IM=1,JPOINT
               PP11 = CONE
               PP12 = CONE
               PP21 = CONE
               PP22 = CONE
               SHI  = ZERO
               FACM = ONE
               DCOH = 1.0D10

C------------------------------------------------------------

               GSH = 10.0D0/(TWO*BSLOPE*GEV2MB)
*    common /DAMP/
               GA  = GSH
               RCA = GA*SIGSH/TWOPI
               FCA = -ROSH*RCA
               CA  = DCMPLX(RCA,FCA)
               CI  = CONE

               DO 15 INA=1,NA
                  KK1  = 1
                  KK2  = 1
                  DO 16 INB=1,NB

                     X11 = B+COOT1(1,INB)-COOP1(1,INA)
                     Y11 =   COOT1(2,INB)-COOP1(2,INA)
                     XY11 = GA*(X11*X11+Y11*Y11)
                     X12 = B+COOT2(1,INB)-COOP1(1,INA)
                     Y12 =   COOT2(2,INB)-COOP1(2,INA)
                     XY12 = GA*(X12*X12+Y12*Y12)
                     X21 = B+COOT1(1,INB)-COOP2(1,INA)
                     Y21 =   COOT1(2,INB)-COOP2(2,INA)
                     XY21 = GA*(X21*X21+Y21*Y21)
                     X22 = B+COOT2(1,INB)-COOP2(1,INA)
                     Y22 =   COOT2(2,INB)-COOP2(2,INA)
                     XY22 = GA*(X22*X22+Y22*Y22)
                     IF (XY11.LE.15.0D0) THEN
                        C = CONE-CA*EXP(-XY11)
                        AR = DBLE(PP11)
                        AI = DIMAG(PP11)
                        IF (ABS(AR).LT.TINY25) AR = ZERO
                        IF (ABS(AI).LT.TINY25) AI = ZERO
                        PP11 = DCMPLX(AR,AI)
                        PP11 = PP11*C
                        AR  = DBLE(C)
                        AI  = DIMAG(C)
                        SHI = SHI+LOG(AR*AR+AI*AI)
                     ENDIF
                     IF (XY12.LE.15.0D0) THEN
                        C = CONE-CA*EXP(-XY12)
                        AR = DBLE(PP12)
                        AI = DIMAG(PP12)
                        IF (ABS(AR).LT.TINY25) AR = ZERO
                        IF (ABS(AI).LT.TINY25) AI = ZERO
                        PP12 = DCMPLX(AR,AI)
                        PP12 = PP12*C
                     ENDIF
                     IF (XY21.LE.15.0D0) THEN
                        C = CONE-CA*EXP(-XY21)
                        AR = DBLE(PP21)
                        AI = DIMAG(PP21)
                        IF (ABS(AR).LT.TINY25) AR = ZERO
                        IF (ABS(AI).LT.TINY25) AI = ZERO
                        PP21 = DCMPLX(AR,AI)
                        PP21 = PP21*C
                     ENDIF
                     IF (XY22.LE.15.0D0) THEN
                        C = CONE-CA*EXP(-XY22)
                        AR = DBLE(PP22)
                        AI = DIMAG(PP22)
                        IF (ABS(AR).LT.TINY25) AR = ZERO
                        IF (ABS(AI).LT.TINY25) AI = ZERO
                        PP22 = DCMPLX(AR,AI)
                        PP22 = PP22*C
                     ENDIF
   16             CONTINUE
   15          CONTINUE

               OMPP11 = CZERO
               OMPP21 = CZERO
               OMPP11 = OMPP11+(CONE-PP11)
               OMPP21 = OMPP21+(CONE-PP21)
               OMPP12 = CZERO
               OMPP22 = CZERO
               OMPP12 = OMPP12+(CONE-PP12)
               OMPP22 = OMPP22+(CONE-PP22)

               STOTM = DBLE(OMPP11+OMPP22)
               SELAM = DBLE(OMPP11*DCONJG(OMPP22))
               SPROM = ONE-EXP(SHI)
               SQEPM = DBLE(OMPP11*DCONJG(OMPP21))-SELAM
               SQETM = DBLE(OMPP11*DCONJG(OMPP12))-SELAM
               SQE2M = DBLE(OMPP11*DCONJG(OMPP11))-SELAM-SQEPM-SQETM

               STOTB = STOTB+FACM*STOTM
               SELAB = SELAB+FACM*SELAM
               IF (NB.GT.1) SQEPB = SQEPB+FACM*SQEPM
               IF (NA.GT.1) SQETB = SQETB+FACM*SQETM
               IF ((NA.GT.1).AND.(NB.GT.1)) SQE2B = SQE2B+FACM*SQE2M
               SPROB = SPROB+FACM*SPROM

C  14       CONTINUE

            STOTN = STOTN+FACB*STOTB
            SELAN = SELAN+FACB*SELAB
            SQEPN = SQEPN+FACB*SQEPB
            SQETN = SQETN+FACB*SQETB
            SQE2N = SQE2N+FACB*SQE2B
            SPRON = SPRON+FACB*SPROB
            BPROD(IB+1)= BPROD(IB+1)+FACN*FACB*SPROB

   12    CONTINUE

         STOT  = STOT +FACN*STOTN
         STOT2 = STOT2+FACN*STOTN**2
         SELA  = SELA +FACN*SELAN
         SELA2 = SELA2+FACN*SELAN**2
         SQEP  = SQEP +FACN*SQEPN
         SQEP2 = SQEP2+FACN*SQEPN**2
         SQET  = SQET +FACN*SQETN
         SQET2 = SQET2+FACN*SQETN**2
         SQE2  = SQE2 +FACN*SQE2N
         SQE22 = SQE22+FACN*SQE2N**2
         SPRO  = SPRO +FACN*SPRON
         SPRO2 = SPRO2+FACN*SPRON**2

   11 CONTINUE

* final cross sections
* 1) total
      XSTOT(IE) = STOT
* 2) elastic
      XSELA(IE) = SELA
* 3) quasi-el.: A+B-->A+X (excluding 2)
      XSQEP(IE) = SQEP
* 4) quasi-el.: A+B-->X+B (excluding 2)
      XSQET(IE) = SQET
* 5) quasi-el.: A+B-->X (excluding 2-4)
      XSQE2(IE) = SQE2
* 6) production (= STOT-SELA-SQEP-SQET-SQE2!)
      XSPRO(IE) = SPRO
      WRITE(6,*)' STOT,SELA ,SQEP,SQET,SQE2,SPRO ',
     & STOT,SELA ,SQEP,SQET,SQE2,SPRO
*  stat. errors
      XETOT(IE) = SQRT(ABS(STOT2-STOT**2)/DBLE(NSTATB-1))
      XEELA(IE) = SQRT(ABS(SELA2-SELA**2)/DBLE(NSTATB-1))
      XEQEP(IE) = SQRT(ABS(SQEP2-SQEP**2)/DBLE(NSTATB-1))
      XEQET(IE) = SQRT(ABS(SQET2-SQET**2)/DBLE(NSTATB-1))
      XEQE2(IE) = SQRT(ABS(SQE22-SQE2**2)/DBLE(NSTATB-1))
      XEPRO(IE) = SQRT(ABS(SPRO2-SPRO**2)/DBLE(NSTATB-1))
      WRITE(6,*)' XETOT,XEELA,XEQEP,XEQET,XEQE2,XEPRO ',
     & XETOT(IE),XEELA(IE),XEQEP(IE),
     &XEQET(IE),XEQE2(IE),XEPRO(IE)
1123  CONTINUE
      DO 19 I=2,NSITEB
         BSITE(IE,NTARG,I) = BPROD(I)/SPRO+BSITE(IE,NTARG,I-1)
         IF (IE.EQ.1)
     &      BSITE(0,NTARG,I) = BPROD(I)/SPRO+BSITE(0,NTARG,I-1)
   19 CONTINUE
      WRITE (6,*)' ECMNN,ELABB,SIGSHH,SIGTO,SIGEL,SIGIN,SIGSD'
C    &          SIGTO(NEB),SIGEL(NEB),SIGIN(NEB),SIGSD(NEB),SIGDIF(NEB)
      DO 129 I=1,NGRITT
        SIGTO(I)=DSHNTO(1,1,ECMNN(I))
        SIGEL(I)=DSHNEL(1,1,ECMNN(I))
        SIGIN(I)=SIINEL(1,1,ECMNN(I))
        SIGSD(I)=SIPPSD(ECMNN(I))
	CALL SIHNDI(ECMNN(I),1,1,SIGDIF(I),SIGDIH)
        WRITE (6,'(2F18.4,6F11.3)')ECMNN(I),ELABB(I),SIGSHH(I),
     &  SIGTO(I),SIGEL(I),SIGIN(I),SIGSD(I),SIGDIF(I) 
 129  CONTINUE      
      WRITE (6,*)' ECMNN,ELABB,XSQEP,XEQEP,XSQET,XEQET,XSQE2,XEQE2'
      DO 139 I=1,NGRITT
      WRITE (6,'(2F18.4,6F11.3)')ECMNN(I),ELABB(I),XSQEP(I),XEQEP(I),
     * XSQET(I),XEQET(I),XSQE2(I),XEQE2(I)
 139  CONTINUE      
      WRITE (6,*)' ECMNN,ELABB,XSTOT,XETOT,XSELA,XEELA,XSPRO,XEPRO'
      DO 119 I=1,NGRITT
      WRITE (6,'(2F18.4,6F11.3)')ECMNN(I),ELABB(I),XSTOT(I),XETOT(I),
     * XSELA(I),XEELA(I),XSPRO(I),XEPRO(I)
 119  CONTINUE      

      RETURN
      END


      SUBROUTINE CONUCLX(COOP1,NA,RASH,I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (MAXNCL = 210)
      DIMENSION COOP1(3,MAXNCL)
      CALL CONUCL(COOP1,NA,RASH)
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE DBKLAS(I,J,K,I8,I10)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C*** I,J,K QUARK FLAVOURS U,D,S,C=1,2,3,4
C*** AQ = -Q
C*** I8,I10 BARYON INDICES
*KEEP,DINPDA.
      COMMON /DINPDA/ IMPS(6,6),IMVE(6,6),IB08(6,21),IB10(6,21), IA08
     +(6,21),IA10(6,21), A1,B1,B2,B3,LT,LE,BET,AS,B8,AME,DIQ,ISU
*KEND.
      IF (I) 20,20,10
C*** BARYON
   10 CONTINUE
      CALL INDEXD(J,K,IND)
      I8=IB08(I,IND)
      I10=IB10(I,IND)
      IF (I8.LE.0) I8=I10
      RETURN
   20 CONTINUE
C*** ANTIBARYONS
      II=IABS(I)
      JJ=IABS(J)
      KK=IABS(K)
      CALL INDEXD(JJ,KK,IND)
      I8=IA08(II,IND)
      I10=IA10(II,IND)
      IF (I8.LE.0) I8=I10
      RETURN
      END
C-----------------------------------------------------------

      DOUBLE PRECISION FUNCTION SIPPSD(ECM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     Single Diffraction cross section in p-p collision
C     Tables calculated with DPMJET-II.4.2
      DIMENSION EC(30),SD(30)
      DATA EC /0.D0, 5.D0, 20.D0, 50.D0, 100.D0,
     *       200.D0, 500.D0, 1000.D0, 1500.D0, 2000.D0,
     *      3000.D0, 4000.D0, 6000.D0, 8000.D0, 10000.D0,
     *     15000.D0, 20000.D0, 30000.D0, 40000.D0, 60000.D0, 
     *     80000.D0, 100000.D0, 150000.D0, 200000.D0, 300000.D0,
     *    400000.D0, 600000.D0, 800000.D0, 1000000.D0, 2000000.D0/
      DATA SD /0.D0, 0.D0, 5.00D0, 6.14D0, 6.93D0,
     *       7.64D0, 8.43D0, 8.87D0, 9.07D0, 9.17D0,
     *       9.33D0, 9.40D0, 9.49D0, 9.56D0, 9.58D0,
     *       9.69D0, 9.72D0, 9.82D0, 9.85D0, 9.97D0,
     *      10.02D0, 10.03D0, 10.13D0, 10.16D0, 10.25D0,
     *      10.28D0, 10.39D0, 10.42D0, 10.43D0, 10.53D0/
      II=1
      DO 1 I=1,30
        IF((ECM.GE.EC(I)).AND.(ECM.LT.EC(I+1)))THEN
          II=I
          DEL=(ECM-EC(I))*(SD(I+1)-SD(I))/(EC(I+1)-EC(I))
          SIPPSD=SD(I)+DEL
          RETURN
        ENDIF
 1    CONTINUE  
      SIPPSD=0.D0
      RETURN
      END
      DOUBLE PRECISION FUNCTION SIINEL(KPROJ,KTARG,UMO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     Inelastic cross section
      SIINEL=DSHNTO(KPROJ,KTARG,UMO)-DSHNEL(KPROJ,KTARG,UMO)
      RETURN
      END
C---------------------------------------------------------------
C                 was dpmsicha.f
C---------------------------------------------------------------
*$ CREATE PHNSCH.FOR
*COPY PHNSCH
*
*=== phnsch ===========================================================*
*
      DOUBLE PRECISION FUNCTION PHNSCH ( KP, KTARG, PLAB )

C      INCLUDE '(DBLPRC)'
C      INCLUDE '(DIMPAR)'
C      INCLUDE '(IOUNIT)'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Probability for Hadron Nucleon Single CHain interactions:        *
*                                                                      *
*     Created on 30 december 1993  by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 04-jan-94     by    Alfredo Ferrari               *
*                                                                      *
*             modified by J.R.for use in DTUNUC  6.1.94                *
*                                                                      *
*     Input variables:                                                 *
*                      Kp = hadron projectile index (Part numbering    *
*                           scheme)                                    *
*                   Ktarg = target nucleon index (1=proton, 8=neutron) *
*                    Plab = projectile laboratory momentum (GeV/c)     *
*     Output variable:                                                 *
*                  Phnsch = probability per single chain (particle     *
*                           exchange) interactions                     *
*                                                                      *
*----------------------------------------------------------------------*
*
C      INCLUDE '(PAPROP)'
C      INCLUDE '(PART2)'
C      INCLUDE '(QQUARK)'

*$ CREATE DBLPRC.ADD
*COPY DBLPRC
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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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


*$ CREATE PAPROP.ADD
*COPY PAPROP
*
*=== paprop ===========================================================*
*
*----------------------------------------------------------------------*
*     include file: paprop copy                   created 26/11/86 by p*
*     changes: on  16 december 1992 by Alfredo Ferrari                 *
*     included in the following subroutines or functions: not updated  *
*                                                                      *
*     description of the common block(s) and variable(s)               *
*                                                                      *
*     /paprop/ contains particle properties                            *
*        btype  = literal name of the particle                         *
*        am     = particle mass in gev                                 *
*        ichrge = electric charge of the particle                      *
*        iscore = explanations for the scored distribution             *
*        genpar = names of the generalized particles                   *
*        ijdisc = list of the particle types to be discarded           *
*        thalf  = half life of the particle in sec                     *
*        biasdc = decay biasing factors                                *
*        biasin = inelastic interaction biasing factors                *
*        lhadro = flag for hadrons                                     *
*        jspinp = particle spin (in units of 1/2)                      *
*        lbsdcy = logical flag for biased decay: if .true. the biasing *
*                 factor is used as an upper limit to the decay length *
*        lprbsd = logical flag for biased decay: if .true. the biasing *
*                 factor is applied only to primaries                  *
*        lprbsi = logical flag for inelastic interaction biasing: if   *
*                 .true. the biasing factor is applied only to prima-  *
*                 ries                                                 *
*                                                                      *
*----------------------------------------------------------------------*
*
C     LOGICAL LHADRO, LBSDCY, LPRBSD, LPRBSI
C     CHARACTER*8 BTYPE,GENPAR
C     COMMON / PAPROP / AM  (NALLWP), AMDISC (NALLWP), THALF  (NALLWP),
C    &               BIASDC (NALLWP), BIASIN (NALLWP), ICHRGE (NALLWP),
C    &               ISCORE     (10), IJDISC (NALLWP), LHADRO (NALLWP),
C    &               JSPINP (NALLWP), LBSDCY (NALLWP), LPRBSD, LPRBSI
C     COMMON / CHPPRP / BTYPE  (NALLWP), GENPAR (30)

      DIMENSION ICHRGE(39),AM(39)

*$ CREATE PART2.ADD
*COPY PART2
*
*=== part2 ============================================================*
*
*----------------------------------------------------------------------*
*     Include file: part2 copy        Revised on 20-7-90 by A. Ferrari *
*     Note: see also part copy and part3 copy                          *
*     Changes: none                                                    *
*     Included in the following subroutines or functions: not updated  *
*                                                                      *
*     Description of the common block(s) and variable(s)               *
*                                                                      *
*         Kptoip = conversion from part to paprop numbering            *
*         Iptokp = conversion from paprop to part numbering            *
*                                                                      *
*----------------------------------------------------------------------*
*
C     CHARACTER*8  ANAME
C     COMMON / PART / AAM    (IDMAXP), GA     (IDMAXP), TAU    (IDMAXP),
C    &                AAMDSC (IDMAXP), ZMNABS (IDMAXP), ATNMNA (IDMAXP),
C    &                IICH   (IDMAXP), IIBAR  (IDMAXP), K1     (IDMAXP),
C    &                K2     (IDMAXP), KPTOIP (IDMAXP), IPTOKP (NALLWP)
C     COMMON / CHPART / ANAME (IDMAXP)

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
       COMMON /DPAR/ ANAME(210),AAM(210),GA(210),TAU(210),
     +               IICH(210),IIBAR(210),K1(210),K2(210)
       DIMENSION KPTOIP(210),IPTOKP(39)
C      DATA KPTOIP/1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
C    +             11,12,13,14,15,16,17,18,19,20,
C    +             21,22,23,24,25, 0, 0, 0, 0, 0,
C    +             60*0,
C    +              0, 0, 0, 0, 0, 0,34,36,31,32,
C    +             33,35,37, 0, 0, 0, 0, 0,38, 0,
C    +              0, 0, 0, 0,39, 0, 0, 0, 0, 0,
C    +             90*0/
C      DATA IPTOKP/1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
C    +             11,12,13,14,15,16,17,18,19,20,
C    +             21,22,23,24,25, 0, 0, 0, 0, 0,
C    +             99,100,101,97,102,98,103,109,115/
*                                                                      *
*     Conversion from part to paprop numbering                         *
*                                                                      *
      DATA KPTOIP / 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
     & 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 66*0,
     & 34, 36, 31, 32, 33, 35, 37, 5*0, 38, 5*0, 39, 19*0, 27, 28, 74*0/
*                                                                      *
*     Conversion from paprop to part numbering                         *
*                                                                      *
      DATA IPTOKP / 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
     & 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 99,
     & 100, 101, 97, 102, 98, 103, 109, 115 /
       

*$ CREATE QQUARK.ADD
*COPY QQUARK
*
*=== qquark ===========================================================*
*
*----------------------------------------------------------------------*
*                                                                      *
*     Created on    6 february 1991    by        Alfredo Ferrari       *
*                                                  INFN - Milan        *
*                                                                      *
*     Last change  on  6 february 1991  by       Alfredo Ferrari       *
*                                                                      *
*     Included in the following routines :                             *
*                                                                      *
*                     COREVT                                           *
*                     CORRIN                                           *
*                     HADEVV                                           *
*                     HADEVT                                           *
*                     NUCEVV                                           *
*                     NUCEVT                                           *
*                                                                      *
*     Quark content of particles:                                      *
*          index   quark   el. charge  bar. charge  isospin  isospin3  *
*              1 = u          2/3          1/3        1/2       1/2    *
*             -1 = ubar      -2/3         -1/3        1/2      -1/2    *
*              2 = d         -1/3          1/3        1/2      -1/2    *
*             -2 = dbar       1/3         -1/3        1/2       1/2    *
*              3 = s         -1/3          1/3         0         0     *
*             -3 = sbar       1/3         -1/3         0         0     *
*              4 = c          2/3          1/3         0         0     *
*             -4 = cbar      -2/3         -1/3         0         0     *
*              5 = b         -1/3          1/3         0         0     *
*             -5 = bbar       1/3         -1/3         0         0     *
*              6 = t          2/3          1/3         0         0     *
*             -6 = tbar      -2/3         -1/3         0         0     *
*                                                                      *
*         Mquark = particle quark composition (Paprop numbering)       *
*         Iqechr = electric charge ( in 1/3 unit )                     *
*         Iqbchr = baryonic charge ( in 1/3 unit )                     *
*         Iqichr = isospin ( in 1/2 unit ), z component                *
*         Iqschr = strangeness                                         *
*         Iqcchr = charm                                               *
*         Iquchr = beauty                                              *
*         Iqtchr = ......                                              *
*                                                                      *
*----------------------------------------------------------------------*
*
      COMMON / QQUARK / IQECHR (-6:6), IQBCHR (-6:6), IQICHR (-6:6),
     &                  IQSCHR (-6:6), IQCCHR (-6:6), IQUCHR (-6:6),
     &                  IQTCHR (-6:6), MQUARK (3,39)
C
      DIMENSION SIEAPP (11), SITAPP (16), PLAETB (16)
      DIMENSION SGTCOE (5,33), PLALIM (2,33), IHLP (NALLWP)
      DIMENSION SGTCO1(5,10),SGTCO2(5,8),SGTCO3(5,15) 
      SAVE PLAETB, SIEAPP, SITAPP, SGTCOE, PLALIM, IHLP
      SAVE IQFSC1, IQFSC2, IQBSC1, IQBSC2
      EQUIVALENCE (SGTCO1(1,1),SGTCOE(1,1))
      EQUIVALENCE (SGTCO2(1,1),SGTCOE(1,11))
      EQUIVALENCE (SGTCO3(1,1),SGTCOE(1,19))
*  1=baryon, 2=pion, 3=kaon, 4=antibaryon:
      DATA IHLP/1,4,5*0,1,4,2*0,3,2*2,2*3,1,4,3,3*1,2,
     &    2*3, 2, 4*0, 3*4, 1, 4, 1, 4, 1, 4 /
C     DATA ( ( SGTCOE (J,I), J=1,5 ), I=1,10 ) /
      DATA  SGTCO1  /
* 1st reaction: gamma p total
     &0.147 D+00, ZERZER  , ZERZER   , 0.0022D+00, -0.0170D+00,
* 2nd reaction: gamma d total
     &0.300 D+00, ZERZER  , ZERZER   , 0.0095D+00, -0.057 D+00,
* 3rd reaction: pi+ p total
     &16.4  D+00, 19.3D+00, -0.42D+00, 0.19  D+00, ZERZER     ,
* 4th reaction: pi- p total
     &33.0  D+00, 14.0D+00, -1.36D+00, 0.456 D+00, -4.03  D+00,
* 5th reaction: pi+/- d total
     &56.8  D+00, 42.2D+00, -1.45D+00, 0.65  D+00, -5.39  D+00,
* 6th reaction: K+ p total
     &18.1  D+00, ZERZER  , ZERZER   , 0.26  D+00, -1.0   D+00,
* 7th reaction: K+ n total
     &18.7  D+00, ZERZER  , ZERZER   , 0.21  D+00, -0.89  D+00,
* 8th reaction: K+ d total
     &34.2  D+00, 7.9 D+00, -2.1 D+00, 0.346 D+00, -0.99  D+00,
* 9th reaction: K- p total
     &32.1  D+00, ZERZER  , ZERZER   , 0.66  D+00, -5.6   D+00,
* 10th reaction: K- n total
     &25.2  D+00, ZERZER  , ZERZER   , 0.38  D+00, -2.9   D+00/
C     DATA ( ( SGTCOE (J,I), J=1,5 ), I=11,18 ) /
      DATA  SGTCO2  /
* 11th reaction: K- d total
     &57.6  D+00, ZERZER  , ZERZER   , 1.17  D+00, -9.5   D+00,
* 12th reaction: p p total
     &48.0  D+00, ZERZER  , ZERZER   , 0.522 D+00, -4.51  D+00,
* 13th reaction: p n total
     &47.30 D+00, ZERZER  , ZERZER   , 0.513 D+00, -4.27  D+00,
* 14th reaction: p d total
     &91.3  D+00, ZERZER  , ZERZER   , 1.05  D+00, -8.8   D+00,
* 15th reaction: pbar p total
     &38.4  D+00, 77.6D+00, -0.64D+00, 0.26  D+00, -1.2   D+00,
* 16th reaction: pbar n total
     &ZERZER    ,133.6D+00, -0.70D+00, -1.22 D+00, 13.7   D+00,
* 17th reaction: pbar d total
     &112.  D+00, 125.D+00, -1.08D+00, 1.14  D+00, -12.4  D+00,
* 18th reaction: Lamda p total
     &30.4  D+00, ZERZER  , ZERZER   , ZERZER    , 1.6    D+00/
C     DATA ( ( SGTCOE (J,I), J=1,5 ), I=19,33 ) /
      DATA SGTCO3  /
* 19th reaction: pi+ p elastic
     &ZERZER    , 11.4D+00, -0.4 D+00, 0.079 D+00, ZERZER     ,
* 20th reaction: pi- p elastic
     &1.76  D+00, 11.2D+00, -0.64D+00, 0.043 D+00, ZERZER     ,
* 21st reaction: K+ p elastic
     &5.0   D+00, 8.1 D+00, -1.8 D+00, 0.16  D+00, -1.3   D+00,
* 22nd reaction: K- p elastic
     &7.3   D+00, ZERZER  , ZERZER   , 0.29  D+00, -2.40  D+00,
* 23rd reaction: p p elastic
     &11.9  D+00, 26.9D+00, -1.21D+00, 0.169 D+00, -1.85  D+00,
* 24th reaction: p d elastic
     &16.1  D+00, ZERZER  , ZERZER   , 0.32  D+00, -3.4   D+00,
* 25th reaction: pbar p elastic
     &10.2  D+00, 52.7D+00, -1.16D+00, 0.125 D+00, -1.28  D+00,
* 26th reaction: pbar p elastic bis
     &10.6  D+00, 53.1D+00, -1.19D+00, 0.136 D+00, -1.41  D+00,
* 27th reaction: pbar n elastic
     &36.5  D+00, ZERZER  , ZERZER   , ZERZER    , -11.9  D+00,
* 28th reaction: Lamda p elastic
     &12.3  D+00, ZERZER  , ZERZER   , ZERZER    , -2.4   D+00,
* 29th reaction: K- p ela bis
     &7.24  D+00, 46.0D+00, -4.71D+00, 0.279 D+00, -2.35  D+00,
* 30th reaction: pi- p cx
     &ZERZER    ,0.912D+00, -1.22D+00, ZERZER    , ZERZER     ,
* 31st reaction: K- p cx
     &ZERZER    , 3.39D+00, -1.75D+00, ZERZER    , ZERZER     ,
* 32nd reaction: K+ n cx
     &ZERZER    , 7.18D+00, -2.01D+00, ZERZER    , ZERZER     ,
* 33rd reaction: pbar p cx
     &ZERZER    , 18.8D+00, -2.01D+00, ZERZER    , ZERZER     /
*
      DATA PLALIM /
*          gamma p tot     ,     gamma d tot     ,       pi+ p tot     ,
     & 3.0D+00, 183.D+00, 2.0D+00, 17.8D+00, 4.0D+00, 340.D+00,
*            pi- p tot     ,     pi+/- d tot     ,        K+ p tot     ,
     & 2.5D+00, 370.D+00, 2.5D+00, 370.D+00, 2.0D+00, 310.D+00,
*             K+ n tot     ,        K+ d tot     ,        K- p tot     ,
     & 2.0D+00, 310.D+00, 2.0D+00, 310.D+00, 3.0D+00, 310.D+00,
*             K- n tot     ,        K- d tot     ,         p p tot     ,
     & 1.8D+00, 310.D+00, 3.0D+00, 310.D+00, 3.0D+00, 2100.D+00,
*              p n tot     ,         p d tot     ,      pbar p tot     ,
     & 3.0D+00, 370.D+00, 3.0D+00, 370.D+00, 5.0D+00, 1.73D+06,
*           pbar n tot     ,      pbar d tot     ,     Lamda p tot     ,
     & 1.1D+00, 280.D+00, 2.0D+00, 280.D+00, 0.6D+00, 21.D+00,
*            pi+ p ela     ,       pi- p ela     ,        K+ p ela     ,
     & 2.0D+00, 200.D+00, 2.0D+00, 360.D+00, 2.0D+00, 175.D+00,
*             K- p ela     ,         p p ela     ,         p d ela     ,
     & 3.0D+00, 175.D+00, 3.0D+00, 2100.D+00, 2.0D+00, 384.D+00,
*           pbar p ela     ,      pbar p ela bis ,      pbar n ela     ,
     & 5.0D+00, 1.73D+06, 2.0D+00, 1.59D+05, 1.1D+00, 5.55D+00,
*          Lamda p ela     ,        K- p ela bis ,       pi- p cx      ,
     & 0.6D+00, 24.D+00, 2.0D+00, 175.D+00, 3.5D+00, 200.D+00,
*             K- p cx      ,        K+ n cx      ,      pbar p cx      /
     & 2.0D+00, 40.D+00, 2.0D+00, 12.8D+00, 3.0D+00, 350.D+00/
*  Momenta for which tabulated data exist for elastic/total pbar p
      DATA PLAETB / 0.1D+00, 0.2D+00,
     &      0.3D+00, 0.4D+00, 0.5D+00, 0.6D+00, 0.8D+00, 1.D+00,
     &      1.2D+00, 1.5D+00, 2. D+00, 2.5D+00, 3. D+00, 4.D+00,
     &      4.5D+00, 5. D+00 /
*  Tabulated data for pbar p elastic:
*  The two lowest energy points are educated guesses:
      DATA SIEAPP / 142.D+00, 95.1D+00,
     &      75.0D+00, 70.0D+00, 62.0D+00, 57.0D+00, 48.0D+00,
     &      44.5D+00, 43.5D+00, 38.0D+00, 33.0D+00 /
*  Tabulated data for pbar p total cross section:
      DATA SITAPP /1129.D+00, 424.D+00,
     &      239.D+00, 195.D+00, 172.D+00, 150.D+00, 124.D+00,
     &      117.D+00, 109.D+00, 100.D+00, 90.2D+00, 81.5D+00,
     &      78.0D+00, 72.0D+00, 67.0D+00, 64.8D+00 /
*
*  +-------------------------------------------------------------------*
         ICHRGE(KTARG)=IICH(KTARG)
         AM    (KTARG)=AAM (KTARG)
*  |  Check for pi0 (d-dbar)
      IF ( KP .NE. 26 ) THEN
         IP  = KPTOIP (KP)
         IF(IP.EQ.0)IP=1 
         ICHRGE(IP)=IICH(KP)
         AM    (IP)=AAM (KP)
*  |
*  +-------------------------------------------------------------------*
*  |
      ELSE
         IP = 23
         ICHRGE(IP)=0
      END IF
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |  No such interactions for baryon-baryon
      IF ( IIBAR (KP) .GT. 0 ) THEN
         PHNSCH = ZERZER
         RETURN
*  |
*  +-------------------------------------------------------------------*
*  |  No "annihilation" diagram possible for K+ p/n
      ELSE IF ( IP .EQ. 15 ) THEN
         PHNSCH = ZERZER
         RETURN
*  |
*  +-------------------------------------------------------------------*
*  |  No "annihilation" diagram possible for K0 p/n
      ELSE IF ( IP .EQ. 24 ) THEN
         PHNSCH = ZERZER
         RETURN
*  |
*  +-------------------------------------------------------------------*
*  |  No "annihilation" diagram possible for Omebar p/n
      ELSE IF ( IP .GE. 38 ) THEN
         PHNSCH = ZERZER
         RETURN
      END IF
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |  If the momentum is larger than 50 GeV/c, compute the single
*  |  chain probability at 50 GeV/c and extrapolate to the present
*  |  momentum according to 1/sqrt(s)
*  |  sigma = sigma_sch (50) * sqrt (s(50)/s) + sigma_dch
*  |  P_sch (50) = sigma_sch (50) / ( sigma_dch + sigma_sch (50) )
*  |  sigma_dch / sigma_sch (50) = 1 / P_sch (50) - 1
*  |  sigma_dch / sigma_sch = 1 / P_sch - 1 = ( 1 / P_sch (50) - 1 )
*  |                        x sqrt(s/s(50))
*  |  P_sch = 1 / [ ( 1 / P_sch (50) - 1 ) x sqrt(s/s(50)) + 1 ]
      IF ( PLAB .GT. 50.D+00 ) THEN
         PLA    = 50.D+00
         AMPSQ  = AM (IP)**2
         AMTSQ  = AM (KTARG)**2
         EPROJ  = SQRT ( PLAB**2 + AMPSQ )
         UMOSQ  = AMPSQ + AMTSQ + TWOTWO * AM (KTARG) * EPROJ
         EPROJ  = SQRT ( PLA**2 + AMPSQ )
         UMO50  = AMPSQ + AMTSQ + TWOTWO * AM (KTARG) * EPROJ
         UMORAT = SQRT ( UMOSQ / UMO50 )
*  |
*  +-------------------------------------------------------------------*
*  |  P < 3 GeV/c
      ELSE IF ( PLAB .LT. 3.D+00 ) THEN
         PLA    = 3.D+00
         AMPSQ  = AM (IP)**2
         AMTSQ  = AM (KTARG)**2
         EPROJ  = SQRT ( PLAB**2 + AMPSQ )
         UMOSQ  = AMPSQ + AMTSQ + TWOTWO * AM (KTARG) * EPROJ
         EPROJ  = SQRT ( PLA**2 + AMPSQ )
         UMO50  = AMPSQ + AMTSQ + TWOTWO * AM (KTARG) * EPROJ
         UMORAT = SQRT ( UMOSQ / UMO50 )
*  |
*  +-------------------------------------------------------------------*
*  |  P < 50 GeV/c
      ELSE
         PLA    = PLAB
         UMORAT = ONEONE
      END IF
*  |
*  +-------------------------------------------------------------------*
      ALGPLA = LOG (PLA)
*  +-------------------------------------------------------------------*
*  |  Pions:
      IF ( IHLP (IP) .EQ. 2 ) THEN
         ACOF = SGTCOE (1,3)
         BCOF = SGTCOE (2,3)
         ENNE = SGTCOE (3,3)
         CCOF = SGTCOE (4,3)
         DCOF = SGTCOE (5,3)
*  |  Compute the pi+ p total cross section:
         SPPPTT = ACOF + BCOF * PLA**ENNE + CCOF * ALGPLA**2
     &          + DCOF * ALGPLA
         ACOF = SGTCOE (1,19)
         BCOF = SGTCOE (2,19)
         ENNE = SGTCOE (3,19)
         CCOF = SGTCOE (4,19)
         DCOF = SGTCOE (5,19)
*  |  Compute the pi+ p elastic cross section:
         SPPPEL = ACOF + BCOF * PLA**ENNE + CCOF * ALGPLA**2
     &          + DCOF * ALGPLA
*  |  Compute the pi+ p inelastic cross section:
         SPPPIN = SPPPTT - SPPPEL
         ACOF = SGTCOE (1,4)
         BCOF = SGTCOE (2,4)
         ENNE = SGTCOE (3,4)
         CCOF = SGTCOE (4,4)
         DCOF = SGTCOE (5,4)
*  |  Compute the pi- p total cross section:
         SPMPTT = ACOF + BCOF * PLA**ENNE + CCOF * ALGPLA**2
     &          + DCOF * ALGPLA
         ACOF = SGTCOE (1,20)
         BCOF = SGTCOE (2,20)
         ENNE = SGTCOE (3,20)
         CCOF = SGTCOE (4,20)
         DCOF = SGTCOE (5,20)
*  |  Compute the pi- p elastic cross section:
         SPMPEL = ACOF + BCOF * PLA**ENNE + CCOF * ALGPLA**2
     &          + DCOF * ALGPLA
*  |  Compute the pi- p inelastic cross section:
         SPMPIN = SPMPTT - SPMPEL
         SIGDIA = SPMPIN - SPPPIN
*  |  +----------------------------------------------------------------*
*  |  |  Charged pions: besides isospin consideration it is supposed
*  |  |                 that (pi+ n)el is almost equal to (pi- p)el
*  |  |                 and  (pi+ p)el "    "     "    "  (pi- n)el
*  |  |                 and all are almost equal among each others
*  |  |                 (reasonable above 5 GeV/c)
         IF ( ICHRGE (IP) .NE. 0 ) THEN
            KHELP = KTARG / 8
            JREAC = 3 + IP - 13 + ICHRGE (IP) * KHELP
            ACOF = SGTCOE (1,JREAC)
            BCOF = SGTCOE (2,JREAC)
            ENNE = SGTCOE (3,JREAC)
            CCOF = SGTCOE (4,JREAC)
            DCOF = SGTCOE (5,JREAC)
*  |  |  Compute the total cross section:
            SHNCTT = ACOF + BCOF * PLA**ENNE + CCOF * ALGPLA**2
     &             + DCOF * ALGPLA
            JREAC = 19 + IP - 13 + ICHRGE (IP) * KHELP
            ACOF = SGTCOE (1,JREAC)
            BCOF = SGTCOE (2,JREAC)
            ENNE = SGTCOE (3,JREAC)
            CCOF = SGTCOE (4,JREAC)
            DCOF = SGTCOE (5,JREAC)
*  |  |  Compute the elastic cross section:
            SHNCEL = ACOF + BCOF * PLA**ENNE + CCOF * ALGPLA**2
     &             + DCOF * ALGPLA
*  |  |  Compute the inelastic cross section:
            SHNCIN = SHNCTT - SHNCEL
*  |  |  Number of diagrams:
            NDIAGR = 1 + IP - 13 + ICHRGE (IP) * KHELP
*  |  |  Now compute the chain end (anti)quark-(anti)diquark
            IQFSC1 = 1 + IP - 13
            IQFSC2 = 0
            IQBSC1 = 1 + KHELP
            IQBSC2 = 1 + IP - 13
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |  pi0: besides isospin consideration it is supposed that the
*  |  |       elastic cross section is not very different from
*  |  |       pi+ p and/or pi- p (reasonable above 5 GeV/c)
         ELSE
            KHELP  = KTARG / 8
            K2HLP  = ( KP - 23 ) / 3
*  |  |  Number of diagrams:
*  |  |  For u ubar (k2hlp=0):
*           NDIAGR = 2 - KHELP
*  |  |  For d dbar (k2hlp=1):
*           NDIAGR = 2 + KHELP - K2HLP
            NDIAGR = 2 + KHELP * ( 2 * K2HLP - 1 ) - K2HLP
            SHNCIN = HLFHLF * ( SPPPIN + SPMPIN )
*  |  |  Now compute the chain end (anti)quark-(anti)diquark
            IQFSC1 = 1 + K2HLP
            IQFSC2 = 0
            IQBSC1 = 1 + KHELP
            IQBSC2 = 2 - K2HLP
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  |                                                   end pi's
*  +-------------------------------------------------------------------*
*  |  Kaons:
      ELSE IF ( IHLP (IP) .EQ. 3 ) THEN
         ACOF = SGTCOE (1,6)
         BCOF = SGTCOE (2,6)
         ENNE = SGTCOE (3,6)
         CCOF = SGTCOE (4,6)
         DCOF = SGTCOE (5,6)
*  |  Compute the K+ p total cross section:
         SKPPTT = ACOF + BCOF * PLA**ENNE + CCOF * ALGPLA**2
     &          + DCOF * ALGPLA
         ACOF = SGTCOE (1,21)
         BCOF = SGTCOE (2,21)
         ENNE = SGTCOE (3,21)
         CCOF = SGTCOE (4,21)
         DCOF = SGTCOE (5,21)
*  |  Compute the K+ p elastic cross section:
         SKPPEL = ACOF + BCOF * PLA**ENNE + CCOF * ALGPLA**2
     &          + DCOF * ALGPLA
*  |  Compute the K+ p inelastic cross section:
         SKPPIN = SKPPTT - SKPPEL
         ACOF = SGTCOE (1,9)
         BCOF = SGTCOE (2,9)
         ENNE = SGTCOE (3,9)
         CCOF = SGTCOE (4,9)
         DCOF = SGTCOE (5,9)
*  |  Compute the K- p total cross section:
         SKMPTT = ACOF + BCOF * PLA**ENNE + CCOF * ALGPLA**2
     &          + DCOF * ALGPLA
         ACOF = SGTCOE (1,22)
         BCOF = SGTCOE (2,22)
         ENNE = SGTCOE (3,22)
         CCOF = SGTCOE (4,22)
         DCOF = SGTCOE (5,22)
*  |  Compute the K- p elastic cross section:
         SKMPEL = ACOF + BCOF * PLA**ENNE + CCOF * ALGPLA**2
     &          + DCOF * ALGPLA
*  |  Compute the K- p inelastic cross section:
         SKMPIN = SKMPTT - SKMPEL
         SIGDIA = HLFHLF * ( SKMPIN - SKPPIN )
*  |  +----------------------------------------------------------------*
*  |  |  Charged Kaons: actually only K-
         IF ( ICHRGE (IP) .NE. 0 ) THEN
            KHELP = KTARG / 8
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Proton target:
            IF ( KHELP .EQ. 0 ) THEN
               SHNCIN = SKMPIN
*  |  |  |  Number of diagrams:
               NDIAGR = 2
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Neutron target: besides isospin consideration it is supposed
*  |  |  |              that (K- n)el is almost equal to (K- p)el
*  |  |  |              (reasonable above 5 GeV/c)
            ELSE
               ACOF = SGTCOE (1,10)
               BCOF = SGTCOE (2,10)
               ENNE = SGTCOE (3,10)
               CCOF = SGTCOE (4,10)
               DCOF = SGTCOE (5,10)
*  |  |  |  Compute the total cross section:
               SHNCTT = ACOF + BCOF * PLA**ENNE + CCOF * ALGPLA**2
     &                + DCOF * ALGPLA
*  |  |  |  Compute the elastic cross section:
               SHNCEL = SKMPEL
*  |  |  |  Compute the inelastic cross section:
               SHNCIN = SHNCTT - SHNCEL
*  |  |  |  Number of diagrams:
               NDIAGR = 1
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  Now compute the chain end (anti)quark-(anti)diquark
            IQFSC1 = 3
            IQFSC2 = 0
            IQBSC1 = 1 + KHELP
            IQBSC2 = 2
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |  K0's: (actually only K0bar)
         ELSE
            KHELP  = KTARG / 8
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Proton target: (K0bar p)in supposed to be given by
*  |  |  |                 (K- p)in - Sig_diagr
            IF ( KHELP .EQ. 0 ) THEN
               SHNCIN = SKMPIN - SIGDIA
*  |  |  |  Number of diagrams:
               NDIAGR = 1
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Neutron target: (K0bar n)in supposed to be given by
*  |  |  |                 (K- n)in + Sig_diagr
*  |  |  |              besides isospin consideration it is supposed
*  |  |  |              that (K- n)el is almost equal to (K- p)el
*  |  |  |              (reasonable above 5 GeV/c)
            ELSE
               ACOF = SGTCOE (1,10)
               BCOF = SGTCOE (2,10)
               ENNE = SGTCOE (3,10)
               CCOF = SGTCOE (4,10)
               DCOF = SGTCOE (5,10)
*  |  |  |  Compute the total cross section:
               SHNCTT = ACOF + BCOF * PLA**ENNE + CCOF * ALGPLA**2
     &                + DCOF * ALGPLA
*  |  |  |  Compute the elastic cross section:
               SHNCEL = SKMPEL
*  |  |  |  Compute the inelastic cross section:
               SHNCIN = SHNCTT - SHNCEL + SIGDIA
*  |  |  |  Number of diagrams:
               NDIAGR = 2
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  Now compute the chain end (anti)quark-(anti)diquark
            IQFSC1 = 3
            IQFSC2 = 0
            IQBSC1 = 1
            IQBSC2 = 1 + KHELP
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  |                                                   end Kaon's
*  +-------------------------------------------------------------------*
*  |  Antinucleons:
      ELSE IF ( IHLP (IP) .EQ. 4 .AND. IP .LE. 9 ) THEN
*  |  For momenta between 3 and 5 GeV/c the use of tabulated data
*  |  should be implemented!
         ACOF = SGTCOE (1,15)
         BCOF = SGTCOE (2,15)
         ENNE = SGTCOE (3,15)
         CCOF = SGTCOE (4,15)
         DCOF = SGTCOE (5,15)
*  |  Compute the pbar p total cross section:
         SAPPTT = ACOF + BCOF * PLA**ENNE + CCOF * ALGPLA**2
     &          + DCOF * ALGPLA
         IF ( PLA .LT. FIVFIV ) THEN
            JREAC = 26
         ELSE
            JREAC = 25
         END IF
         ACOF = SGTCOE (1,JREAC)
         BCOF = SGTCOE (2,JREAC)
         ENNE = SGTCOE (3,JREAC)
         CCOF = SGTCOE (4,JREAC)
         DCOF = SGTCOE (5,JREAC)
*  |  Compute the pbar p elastic cross section:
         SAPPEL = ACOF + BCOF * PLA**ENNE + CCOF * ALGPLA**2
     &          + DCOF * ALGPLA
*  |  Compute the pbar p inelastic cross section:
         SAPPIN = SAPPTT - SAPPEL
         ACOF = SGTCOE (1,12)
         BCOF = SGTCOE (2,12)
         ENNE = SGTCOE (3,12)
         CCOF = SGTCOE (4,12)
         DCOF = SGTCOE (5,12)
*  |  Compute the p p total cross section:
         SPPTOT = ACOF + BCOF * PLA**ENNE + CCOF * ALGPLA**2
     &          + DCOF * ALGPLA
         ACOF = SGTCOE (1,23)
         BCOF = SGTCOE (2,23)
         ENNE = SGTCOE (3,23)
         CCOF = SGTCOE (4,23)
         DCOF = SGTCOE (5,23)
*  |  Compute the p p elastic cross section:
         SPPELA = ACOF + BCOF * PLA**ENNE + CCOF * ALGPLA**2
     &          + DCOF * ALGPLA
*  |  Compute the K- p inelastic cross section:
         SPPINE = SPPTOT - SPPELA
         SIGDIA = ( SAPPIN - SPPINE ) / FIVFIV
         KHELP  = KTARG / 8
*  |  +----------------------------------------------------------------*
*  |  |  Pbar:
         IF ( ICHRGE (IP) .NE. 0 ) THEN
            NDIAGR = 5 - KHELP
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Proton target:
            IF ( KHELP .EQ. 0 ) THEN
*  |  |  |  Number of diagrams:
               SHNCIN = SAPPIN
               PUUBAR = 0.8D+00
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Neutron target: it is supposed that (ap n)el is almost equal
*  |  |  |                  to (ap p)el (reasonable above 5 GeV/c)
            ELSE
               ACOF = SGTCOE (1,16)
               BCOF = SGTCOE (2,16)
               ENNE = SGTCOE (3,16)
               CCOF = SGTCOE (4,16)
               DCOF = SGTCOE (5,16)
*  |  |  |  Compute the total cross section:
               SHNCTT = ACOF + BCOF * PLA**ENNE + CCOF * ALGPLA**2
     &                + DCOF * ALGPLA
*  |  |  |  Compute the elastic cross section:
               SHNCEL = SAPPEL
*  |  |  |  Compute the inelastic cross section:
               SHNCIN = SHNCTT - SHNCEL
               PUUBAR = HLFHLF
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  Now compute the chain end (anti)quark-(anti)diquark
*  |  |  there are different possibilities, make a random choiche:
            IQFSC1 = -1
            RNCHEN = RNDM (RNCHEN)
            IF ( RNCHEN .LT. PUUBAR ) THEN
               IQFSC2 = -2
            ELSE
               IQFSC2 = -1
            END IF
            IQBSC1 = -IQFSC1 + KHELP
            IQBSC2 = -IQFSC2
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |  nbar:
         ELSE
            NDIAGR = 4 + KHELP
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Proton target: (nbar p)in supposed to be given by
*  |  |  |                 (pbar p)in - Sig_diagr
            IF ( KHELP .EQ. 0 ) THEN
               SHNCIN = SAPPIN - SIGDIA
               PDDBAR = HLFHLF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Neutron target: (nbar n)el is supposed to be equal to
*  |  |  |                  (pbar p)el (reasonable above 5 GeV/c)
            ELSE
*  |  |  |  Compute the total cross section:
               SHNCTT = SAPPTT
*  |  |  |  Compute the elastic cross section:
               SHNCEL = SAPPEL
*  |  |  |  Compute the inelastic cross section:
               SHNCIN = SHNCTT - SHNCEL
               PDDBAR = 0.8D+00
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  Now compute the chain end (anti)quark-(anti)diquark
*  |  |  there are different possibilities, make a random choiche:
            IQFSC1 = -2
            RNCHEN = RNDM (RNCHEN)
            IF ( RNCHEN .LT. PDDBAR ) THEN
               IQFSC2 = -1
            ELSE
               IQFSC2 = -2
            END IF
            IQBSC1 = -IQFSC1 + KHELP - 1
            IQBSC2 = -IQFSC2
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  |
*  +-------------------------------------------------------------------*
*  |  Others: not yet implemented
      ELSE
         SIGDIA = ZERZER
         SHNCIN = ONEONE
         NDIAGR = 0
         PHNSCH = ZERZER
         RETURN
      END IF
*  |                                                   end others
*  +-------------------------------------------------------------------*
      PHNSCH = NDIAGR * SIGDIA / SHNCIN
      IQECHC = IQECHR (IQFSC1) + IQECHR (IQFSC2) + IQECHR (IQBSC1)
     &       + IQECHR (IQBSC2)
      IQBCHC = IQBCHR (IQFSC1) + IQBCHR (IQFSC2) + IQBCHR (IQBSC1)
     &       + IQBCHR (IQBSC2)
      IQECHC = IQECHC / 3
      IQBCHC = IQBCHC / 3
      IQSCHC = IQSCHR (IQFSC1) + IQSCHR (IQFSC2) + IQSCHR (IQBSC1)
     &       + IQSCHR (IQBSC2)
      IQSPRO = IQSCHR (MQUARK(1,IP)) + IQSCHR (MQUARK(2,IP))
     &       + IQSCHR (MQUARK(3,IP))
*  +-------------------------------------------------------------------*
*  |  Consistency check:
      IF ( PHNSCH .LE. ZERZER .OR. PHNSCH .GT. ONEONE ) THEN
         WRITE (LUNOUT,*)' *** Phnsch,kp,ktarg,pla',
     &                         PHNSCH,KP,KTARG,PLA,' ****'
         WRITE (LUNERR,*)' *** Phnsch,kp,ktarg,pla',
     &                         PHNSCH,KP,KTARG,PLA,' ****'
         PHNSCH = MAX ( PHNSCH, ZERZER )
         PHNSCH = MIN ( PHNSCH, ONEONE )
      END IF
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |  Consistency check:
      IF ( IQSPRO .NE. IQSCHC .OR. ICHRGE (IP) + ICHRGE (KTARG)
     &     .NE. IQECHC .OR. IIBAR (KP) + IIBAR (KTARG) .NE. IQBCHC) THEN
         WRITE (LUNOUT,*)
     &' *** Phnsch,iqspro,iqschc,ichrge,iqechc,ibar,iqbchc,ktarg',
     &      IQSPRO,IQSCHC,ICHRGE(IP),IQECHC,IIBAR(KP),IQBCHC,KTARG
         WRITE (LUNERR,*)
     &' *** Phnsch,iqspro,iqschc,ichrge,iqechc,ibar,iqbchc,ktarg',
     &      IQSPRO,IQSCHC,ICHRGE(IP),IQECHC,IIBAR(KP),IQBCHC,KTARG
      END IF
*  |
*  +-------------------------------------------------------------------*
*  P_sch = 1 / [ ( 1 / P_sch (50) - 1 ) x sqrt(s/s(50)) + 1 ]
      IF ( UMORAT .GT. ONEPLS ) PHNSCH = ONEONE / ( ( ONEONE / PHNSCH
     &                                 - ONEONE ) * UMORAT + ONEONE )
      RETURN
*
      ENTRY SCHQUA ( JQFSC1, JQFSC2, JQBSC1, JQBSC2 )
      SCHQUA = ONEONE
      JQFSC1 = IQFSC1
      JQFSC2 = IQFSC2
      JQBSC1 = IQBSC1
      JQBSC2 = IQBSC2
*=== End of function Phnsch ===========================================*
      RETURN
      END
*
*=== qprop ============================================================*
*
      BLOCK DATA QPROP
*----------------------------------------------------------------------*
*                                                                      *
*     Created on    6 february 1991    by        Alfredo Ferrari       *
*                                                  INFN - Milan        *
*                                                                      *
*     Last change  on  6 february 1991  by       Alfredo Ferrari       *
*                                                                      *
*     Included in the following routines :                             *
*                                                                      *
*                     COREVT                                           *
*                     CORRIN                                           *
*                     HADEVV                                           *
*                     HADEVT                                           *
*                     NUCEVV                                           *
*                     NUCEVT                                           *
*                                                                      *
*     Quark content of particles:                                      *
*          index   quark   el. charge  bar. charge  isospin  isospin3  *
*              1 = u          2/3          1/3        1/2       1/2    *
*             -1 = ubar      -2/3         -1/3        1/2      -1/2    *
*              2 = d         -1/3          1/3        1/2      -1/2    *
*             -2 = dbar       1/3         -1/3        1/2       1/2    *
*              3 = s         -1/3          1/3         0         0     *
*             -3 = sbar       1/3         -1/3         0         0     *
*              4 = c          2/3          1/3         0         0     *
*             -4 = cbar      -2/3         -1/3         0         0     *
*              5 = b         -1/3          1/3         0         0     *
*             -5 = bbar       1/3         -1/3         0         0     *
*              6 = t          2/3          1/3         0         0     *
*             -6 = tbar      -2/3         -1/3         0         0     *
*                                                                      *
*         Mquark = particle quark composition (Paprop numbering)       *
*         Iqechr = electric charge ( in 1/3 unit )                     *
*         Iqbchr = baryonic charge ( in 1/3 unit )                     *
*         Iqichr = isospin ( in 1/2 unit ), z component                *
*         Iqschr = strangeness                                         *
*         Iqcchr = charm                                               *
*         Iquchr = beauty                                              *
*         Iqtchr = ......                                              *
*                                                                      *
*----------------------------------------------------------------------*
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON / QQUARK / IQECHR (-6:6), IQBCHR (-6:6), IQICHR (-6:6),
     &                  IQSCHR (-6:6), IQCCHR (-6:6), IQUCHR (-6:6),
     &                  IQTCHR (-6:6), MQUARK (3,39)

*
* / Qquark / 
      DATA IQECHR / -2, 1, -2, 1, 1, -2, 0, 2, -1, -1, 2, -1, 2 /
      DATA IQBCHR / 6*-1, 0, 6*1 /
      DATA IQICHR / 4*0, 1, -1, 0, 1, -1, 4*0 /
      DATA IQSCHR / 3*0, 1, 5*0, -1, 3*0 /
      DATA IQCCHR / 2*0, -1, 7*0, 1, 2*0 /
      DATA IQUCHR / 0, 1, 9*0, -1, 0 /
      DATA IQTCHR / -1, 11*0, 1 /
      DATA MQUARK /                1,1,2,    -1,-1,-2,
     *   0,0,0,       0,0,0,       0,0,0,       0,0,0,       0,0,0,
     *   1,2,2,    -1,-2,-2,       0,0,0,       0,0,0,       0,0,0,
     *  1,-2,0,      2,-1,0,      1,-3,0,      3,-1,0,
     *   1,2,3,    -1,-2,-3,       0,0,0,
     *   2,2,3,     1,1,3,     1,2,3,     1,-1,0,
     *   2,-3,0,    3,-2,0,    2,-2,0,    0,0,0,
     *   0,0,0,       0,0,0,       0,0,0,
     *  -1,-1,-3,    -1,-2,-3,    -2,-2,-3,
     *   1,3,3,      -1,-3,-3,     2,3,3,      -2,-3,-3,
     *   3,3,3,      -3,-3,-3 /

      END
C******************************************************************
      SUBROUTINE SELPTS( PTXSQ1,PTYSQ1,
     +PLQ1,EQ1,PTXSA2,
     +PTYSA2,PLAQ2,EAQ2, AMCH1,IREJ,IKVALA,PTTQ1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  SELECT PT VALUES FOR A SINGLE CHAIN SYSTEM
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
      QTXSA2=PTXSA2
      QTYSQ1=PTYSQ1
      QTYSA2=PTYSA2
      QLQ1=PLQ1
      QLAQ2=PLAQ2
      QEQ1=EQ1
      QEAQ2=EAQ2
C                                        ----------------
      IANFA=0
      ITAGPT=0
C                           changed from 3.  j.r.21.8.93
      B33=3.00
      IF (IKVALA.EQ.1)B33=6.0
      ICOUNT=0
      IREJ=0
   10 CONTINUE
      ICOUNT=ICOUNT+1
      IF (ICOUNT.EQ.10)THEN
        IREJ=1
C                            REJECT EVENT
        RETURN
      ENDIF
      IF (ICOUNT.GE.1)THEN
        HPS=HPS*0.9
      PTXSQ1=QTXSQ1+HPS*CFE
      PTYSQ1=QTYSQ1+HPS*SFE
      PTXSA2=QTXSA2-HPS*CFE
      PTYSA2=QTYSA2-HPS*SFE
        GO TO 111       
      ENDIF
      B33=2.*B33
C
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
      PTXSA2=QTXSA2-HPS*CFE
      PTYSA2=QTYSA2-HPS*SFE
  111 CONTINUE
C                                      -----------------
C
      IF (IPEV.GE.6)WRITE(6,1000)PTXSQ1,PTYSQ1,
     +PTXSA2,PTYSA2
 1000 FORMAT (' PT S  ',8F12.6)
C                           KINEMATICS OF THE TWO CHAINS Q1-AQ2,AQ1-Q2
      PTTQ1=PTXSQ1**2+PTYSQ1**2
      IF((EQ1**2.LE.PTTQ1))            GO TO 10
C                                      
       IANFA2=0
      ITAGP2=0
      B33=3.00
      IF (IKVALA.EQ.1)B33=6.0
      ICOUN2=0
      IREJ=0
   12 CONTINUE
      ICOUN2=ICOUN2+1
      IF (ICOUN2.EQ.12)THEN
        IREJ=1
C                            REJECT EVENT
        RETURN
      ENDIF
C                                      -----------------
C
      IF (IPEV.GE.6)WRITE(6,1000)PTXSQ1,PTYSQ1,
     +PTXSA2,PTYSA2
C                           KINEMATICS OF THE TWO CHAINS Q1-AQ2,AQ1-Q2
      PTTA2=PTXSA2**2+PTYSA2**2
      IF((EAQ2**2.LE.PTTA2))            GO TO 12
 
C
      IF(IP.GE.1)GO TO 1779
        PLQ1=SQRT(EQ1**2-PTTQ1)
        PLAQ2=-SQRT(EAQ2**2-PTTA2)
 1779 CONTINUE
C-----------
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
      RETURN
      END
