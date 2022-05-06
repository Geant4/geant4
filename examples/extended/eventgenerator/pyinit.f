 
C*********************************************************************
 
C...PYINIT
C...Initializes the generation procedure; finds maxima of the
C...differential cross-sections to be used for weighting.
      SUBROUTINE PYINIT(FRAME,BEAM,TARGET,WIN)
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYDAT4/CHAF(500,2)
      CHARACTER CHAF*16
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
      COMMON/PYPUED/IUED(0:99),RUED(0:99)
      SAVE /PYDAT1/,/PYDAT2/,/PYDAT3/,/PYDAT4/,/PYSUBS/,/PYPARS/,
     &/PYINT1/,/PYINT2/,/PYINT5/,/PYPUED/
C...Local arrays and character variables.
      DIMENSION ALAMIN(20),NFIN(20)
      CHARACTER*(*) FRAME,BEAM,TARGET
      CHARACTER CHFRAM*12,CHBEAM*12,CHTARG*12,CHLH(2)*6
 
C...Interface to PDFLIB.
      COMMON/W50511/NPTYPE,NGROUP,NSET,MODE,NFL,LO,TMAS
      COMMON/W50512/QCDL4,QCDL5
      SAVE /W50511/,/W50512/
      DOUBLE PRECISION VALUE(20),TMAS,QCDL4,QCDL5
      CHARACTER*20 PARM(20)
      DATA VALUE/20*0D0/,PARM/20*' '/
 
C...Data:Lambda and n_f values for parton distributions..
      DATA ALAMIN/0.177D0,0.239D0,0.247D0,0.2322D0,0.248D0,0.248D0,
     &0.192D0,0.326D0,2*0.2D0,0.2D0,0.2D0,0.29D0,0.2D0,0.4D0,5*0.2D0/,
     &NFIN/20*4/
      DATA CHLH/'lepton','hadron'/
C...Check that BLOCK DATA PYDATA has been loaded.
      CALL PYCKBD
 
C...Reset MINT and VINT arrays. Write headers.
      MSTI(53)=0
      DO 100 J=1,400
        MINT(J)=0
        VINT(J)=0D0
  100 CONTINUE
      IF(MSTU(12).NE.12345) CALL PYLIST(0)
      IF(MSTP(122).GE.1) WRITE(MSTU(11),5100)
 
 
C...Reset error counters.
      MSTU(23)=0
      MSTU(27)=0
      MSTU(30)=0
 
 
C...Reset processes that should not be on.
      MSUB(96)=0
      MSUB(97)=0
C...Select global FSR/ISR/UE parameter set = 'tune' 
C...See routine PYTUNE for details
      IF (MSTP(5).NE.0) THEN
        MSTP5=MSTP(5)
        CALL PYTUNE(MSTP5)
      ENDIF

C...Call user process initialization routine.
      IF(FRAME(1:1).EQ.'u'.OR.FRAME(1:1).EQ.'U') THEN
        MSEL=0
        CALL UPINIT
        MSEL=0
      ENDIF
 
C...Maximum 4 generations; set maximum number of allowed flavours.
      MSTP(1)=MIN(4,MSTP(1))
      MSTU(114)=MIN(MSTU(114),2*MSTP(1))
      MSTP(58)=MIN(MSTP(58),2*MSTP(1))
 
C...Sum up Cabibbo-Kobayashi-Maskawa factors for each quark/lepton.
      DO 120 I=-20,20
        VINT(180+I)=0D0
        IA=IABS(I)
        IF(IA.GE.1.AND.IA.LE.2*MSTP(1)) THEN
          DO 110 J=1,MSTP(1)
            IB=2*J-1+MOD(IA,2)
            IF(IB.GE.6.AND.MSTP(9).EQ.0) GOTO 110
            IPM=(5-ISIGN(1,I))/2
            IDC=J+MDCY(IA,2)+2
            IF(MDME(IDC,1).EQ.1.OR.MDME(IDC,1).EQ.IPM) VINT(180+I)=
     &      VINT(180+I)+VCKM((IA+1)/2,(IB+1)/2)
  110     CONTINUE
        ELSEIF(IA.GE.11.AND.IA.LE.10+2*MSTP(1)) THEN
          VINT(180+I)=1D0
        ENDIF
  120 CONTINUE
C...Initialize parton distributions: PDFLIB.
      IF(MSTP(52).EQ.2) THEN
        PARM(1)='NPTYPE'
        VALUE(1)=1
        PARM(2)='NGROUP'
        VALUE(2)=MSTP(51)/1000
        PARM(3)='NSET'
        VALUE(3)=MOD(MSTP(51),1000)
        PARM(4)='TMAS'
        VALUE(4)=PMAS(6,1)
        CALL PDFSET(PARM,VALUE)
        MINT(93)=1000000+MSTP(51)
      ENDIF
C...Choose Lambda value to use in alpha-strong.
      MSTU(111)=MSTP(2)
      IF(MSTP(3).GE.2) THEN
        ALAM=0.2D0
        NF=4
        IF(MSTP(52).EQ.1.AND.MSTP(51).GE.1.AND.MSTP(51).LE.20) THEN
          ALAM=ALAMIN(MSTP(51))
          NF=NFIN(MSTP(51))
        ELSEIF(MSTP(52).EQ.2.AND.NFL.EQ.5) THEN
          ALAM=QCDL5
          NF=5
        ELSEIF(MSTP(52).EQ.2) THEN
          ALAM=QCDL4
          NF=4
        ENDIF
        PARP(1)=ALAM
        PARP(61)=ALAM
        PARP(72)=ALAM
        PARU(112)=ALAM
        MSTU(112)=NF
        IF(MSTP(3).EQ.3) PARJ(81)=ALAM
      ENDIF
 
C...Initialize the UED masses and widths
      IF (IUED(1).EQ.1) CALL PYXDIN
C...Initialize the SUSY generation: couplings, masses,
C...decay modes, branching ratios, and so on.
      CALL PYMSIN
C...Initialize widths and partial widths for resonances.
      CALL PYINRE
C...Set Z0 mass and width for e+e- routines.
      PARJ(123)=PMAS(23,1)
      PARJ(124)=PMAS(23,2)
C...Identify beam and target particles and frame of process.
      CHFRAM=FRAME(1:len(FRAME))//' '
      CHBEAM=BEAM(1:len(BEAM))//' '
      CHTARG=TARGET(1:len(TARGET))//' '
      CALL PYINBM(CHFRAM,CHBEAM,CHTARG,WIN)
      IF(MINT(65).EQ.1) GOTO 170
C...For gamma-p or gamma-gamma allow many (3 or 6) alternatives.
C...For e-gamma allow 2 alternatives.
      MINT(121)=1
      IF(MSTP(14).EQ.10.AND.(MSEL.EQ.1.OR.MSEL.EQ.2)) THEN
        IF((MINT(11).EQ.22.OR.MINT(12).EQ.22).AND.
     &  (IABS(MINT(11)).GT.100.OR.IABS(MINT(12)).GT.100)) MINT(121)=3
        IF(MINT(11).EQ.22.AND.MINT(12).EQ.22) MINT(121)=6
        IF((MINT(11).EQ.22.OR.MINT(12).EQ.22).AND.
     &  (IABS(MINT(11)).EQ.11.OR.IABS(MINT(12)).EQ.11)) MINT(121)=2
      ELSEIF(MSTP(14).EQ.20.AND.(MSEL.EQ.1.OR.MSEL.EQ.2)) THEN
        IF((MINT(11).EQ.22.OR.MINT(12).EQ.22).AND.
     &  (IABS(MINT(11)).GT.100.OR.IABS(MINT(12)).GT.100)) MINT(121)=3
        IF(MINT(11).EQ.22.AND.MINT(12).EQ.22) MINT(121)=9
      ELSEIF(MSTP(14).EQ.25.AND.(MSEL.EQ.1.OR.MSEL.EQ.2)) THEN
        IF((MINT(11).EQ.22.OR.MINT(12).EQ.22).AND.
     &  (IABS(MINT(11)).GT.100.OR.IABS(MINT(12)).GT.100)) MINT(121)=2
        IF(MINT(11).EQ.22.AND.MINT(12).EQ.22) MINT(121)=4
      ELSEIF(MSTP(14).EQ.30.AND.(MSEL.EQ.1.OR.MSEL.EQ.2)) THEN
        IF((MINT(11).EQ.22.OR.MINT(12).EQ.22).AND.
     &  (IABS(MINT(11)).GT.100.OR.IABS(MINT(12)).GT.100)) MINT(121)=4
        IF(MINT(11).EQ.22.AND.MINT(12).EQ.22) MINT(121)=13
      ENDIF
      MINT(123)=MSTP(14)
      IF((MSTP(14).EQ.10.OR.MSTP(14).EQ.20.OR.MSTP(14).EQ.25.OR.
     &MSTP(14).EQ.30).AND.MSEL.NE.1.AND.MSEL.NE.2) MINT(123)=0
      IF(MSTP(14).GE.11.AND.MSTP(14).LE.19) THEN
        IF(MSTP(14).EQ.11) MINT(123)=0
        IF(MSTP(14).EQ.12.OR.MSTP(14).EQ.14) MINT(123)=5
        IF(MSTP(14).EQ.13.OR.MSTP(14).EQ.17) MINT(123)=6
        IF(MSTP(14).EQ.15) MINT(123)=2
        IF(MSTP(14).EQ.16.OR.MSTP(14).EQ.18) MINT(123)=7
        IF(MSTP(14).EQ.19) MINT(123)=3
      ELSEIF(MSTP(14).GE.21.AND.MSTP(14).LE.24) THEN
        IF(MSTP(14).EQ.21) MINT(123)=0
        IF(MSTP(14).EQ.22.OR.MSTP(14).EQ.23) MINT(123)=4
        IF(MSTP(14).EQ.24) MINT(123)=1
      ELSEIF(MSTP(14).GE.26.AND.MSTP(14).LE.29) THEN
        IF(MSTP(14).EQ.26.OR.MSTP(14).EQ.28) MINT(123)=8
        IF(MSTP(14).EQ.27.OR.MSTP(14).EQ.29) MINT(123)=9
      ENDIF
C...Set up kinematics of process.
      CALL PYINKI(0)
 
C...Set up kinematics for photons inside leptons.
      IF(MINT(141).NE.0.OR.MINT(142).NE.0) CALL PYGAGA(1,WTGAGA)
 
C...Precalculate flavour selection weights.
      CALL PYKFIN
C...Loop over gamma-p or gamma-gamma alternatives.
      CKIN3=CKIN(3)
      MSAV48=0
      DO 160 IGA=1,MINT(121)
        CKIN(3)=CKIN3
        MINT(122)=IGA
 
C...Select partonic subprocesses to be included in the simulation.
        CALL PYINPR
        MINT(101)=1
        MINT(102)=1
        MINT(103)=MINT(11)
        MINT(104)=MINT(12)
 
C...Count number of subprocesses on.
        MINT(48)=0
        DO 130 ISUB=1,500
          IF(MINT(50).EQ.0.AND.ISUB.GE.91.AND.ISUB.LE.96.AND.
     &    MSUB(ISUB).EQ.1.AND.MINT(121).GT.1) THEN
            MSUB(ISUB)=0
          ELSEIF(MINT(50).EQ.0.AND.ISUB.GE.91.AND.ISUB.LE.96.AND.
     &    MSUB(ISUB).EQ.1) THEN
            WRITE(MSTU(11),5200) ISUB,CHLH(MINT(41)),CHLH(MINT(42))
            CALL PYSTOP(1)
          ELSEIF(MSUB(ISUB).EQ.1.AND.ISET(ISUB).EQ.-1) THEN
            WRITE(MSTU(11),5300) ISUB
            CALL PYSTOP(1)
          ELSEIF(MSUB(ISUB).EQ.1.AND.ISET(ISUB).LE.-2) THEN
            WRITE(MSTU(11),5400) ISUB
            CALL PYSTOP(1)
          ELSEIF(MSUB(ISUB).EQ.1) THEN
            MINT(48)=MINT(48)+1
          ENDIF
  130   CONTINUE
 
C...Stop or raise warning flag if no subprocesses on.
        IF(MINT(121).EQ.1.AND.MINT(48).EQ.0) THEN
          IF(MSTP(127).NE.1) THEN
            WRITE(MSTU(11),5500)
            CALL PYSTOP(1)
          ELSE
            WRITE(MSTU(11),5700)
            MSTI(53)=1
          ENDIF
        ENDIF
        MINT(49)=MINT(48)-MSUB(91)-MSUB(92)-MSUB(93)-MSUB(94)
        MSAV48=MSAV48+MINT(48)
 
C...Reset variables for cross-section calculation.
        DO 150 I=0,500
          DO 140 J=1,3
            NGEN(I,J)=0
            XSEC(I,J)=0D0
  140     CONTINUE
  150   CONTINUE
 
C...Find parametrized total cross-sections.
        CALL PYXTOT
        VINT(318)=VINT(317)
 
C...Maxima of differential cross-sections.
        IF(MSTP(121).LE.1) CALL PYMAXI
 
C...Initialize possibility of pileup events.
        IF(MINT(121).GT.1) MSTP(131)=0
        IF(MSTP(131).NE.0) CALL PYPILE(1)
 
C...Initialize multiple interactions with variable impact parameter.
        IF(MINT(50).EQ.1) THEN
          PTMN=PARP(82)*(VINT(1)/PARP(89))**PARP(90)
          IF(MOD(MSTP(81),10).EQ.0.AND.(CKIN(3).GT.PTMN.OR.
     &    ((MSEL.NE.1.AND.MSEL.NE.2)))) MSTP(82)=MIN(1,MSTP(82))
          IF((MINT(49).NE.0.OR.MSTP(131).NE.0).AND.MSTP(82).GE.2) THEN
            MINT(35)=1
            CALL PYMULT(1)
            MINT(35)=3
            CALL PYMIGN(1)
          ENDIF
        ENDIF
 
C...Save results for gamma-p and gamma-gamma alternatives.
        IF(MINT(121).GT.1) CALL PYSAVE(1,IGA)
  160 CONTINUE
 
C...Initialization finished.
      IF(MSAV48.EQ.0) THEN
        IF(MSTP(127).NE.1) THEN
          WRITE(MSTU(11),5500)
          CALL PYSTOP(1)
        ELSE
          WRITE(MSTU(11),5700)
          MSTI(53)=1
        ENDIF
      ENDIF
  170 IF(MSTP(122).GE.1) WRITE(MSTU(11),5600)
 
C...Formats for initialization information.
 5100 FORMAT('1',18('*'),1X,'PYINIT: initialization of PYTHIA ',
     &'routines',1X,17('*'))
 5200 FORMAT(1X,'Error: process number ',I3,' not meaningful for ',A6,
     &'-',A6,' interactions.'/1X,'Execution stopped!')
 5300 FORMAT(1X,'Error: requested subprocess',I4,' not implemented.'/
     &1X,'Execution stopped!')
 5400 FORMAT(1X,'Error: requested subprocess',I4,' not existing.'/
     &1X,'Execution stopped!')
 5500 FORMAT(1X,'Error: no subprocess switched on.'/
     &1X,'Execution stopped.')
 5600 FORMAT(/1X,22('*'),1X,'PYINIT: initialization completed',1X,
     &22('*'))
 5700 FORMAT(1X,'Error: no subprocess switched on.'/
     &1X,'Execution will stop if you try to generate events.')
 
      RETURN
      END
