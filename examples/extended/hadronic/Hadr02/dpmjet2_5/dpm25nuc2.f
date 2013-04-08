C----------------------------------------------------------------------
C
C    FILE DPMNUC2.FOR
C
C----------------------------------------------------------------------
*
      SUBROUTINE PRIMPT(MPO,ECM)
C
C  SELECT PRIMORDIAL PT FOR HARD SCATTERED PARTONS
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER(MSH=250)
      COMMON /ABRHRD/XH1(MSH),XH2(MSH),IJHI1(MSH),IJHI2(MSH),
     *IJHF1(MSH),IJHF2(MSH),PHARD1(MSH,4),PHARD2(MSH,4)
      COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
      B33=4.0+3.0/LOG10(ECM+10.D0)
      DO 20 I=1,MPO
        ES=-2./(B33**2)*LOG(RNDM(V)*RNDM(U))
        HPS=SQRT(ES*ES+2.*ES*0.94)
        CALL DSFECF(SFE,CFE)
        PTXSQ1=HPS*CFE
        PTYSQ1=HPS*SFE
        PTXSA1=-PTXSQ1
        PTYSA1=-PTYSQ1
        IF (IOUXEV.GE.6)WRITE(6,115)PTXSQ1,PTYSQ1,PTXSA1,PTYSA1
  115   FORMAT (' PT S  ',8F12.6)
        PHARD1(I,1)=PHARD1(I,1)+PTXSQ1
        PHARD1(I,2)=PHARD1(I,2)+PTYSQ1
        PHARD2(I,1)=PHARD2(I,1)+PTXSA1
        PHARD2(I,2)=PHARD2(I,2)+PTYSA1
        DE1=SQRT(PHARD1(I,1)**2+PHARD1(I,2)**2+PHARD1(I,3)**2)
     *                                        -PHARD1(I,4)
        DE2=SQRT(PHARD2(I,1)**2+PHARD2(I,2)**2+PHARD2(I,3)**2)
     *                                        -PHARD2(I,4)
        PHARD1(I,4)=PHARD1(I,4)+DE1
        PHARD2(I,4)=PHARD2(I,4)+DE2
        DX1=2.*DE1/ECM
        DX2=2.*DE2/ECM
        XH1(I)=XH1(I)+DX1
        XH2(I)=XH2(I)+DX2
   20 CONTINUE
      RETURN
      END

C______________________________________________________________________
C
C****************************************************************8**
      SUBROUTINE SELPTH(PQUAR,PAQUAR,TQUAR,TAQUAR,ECM,
     *                 PTXSQ1,PTYSQ1,PLQ1,EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1,
     *                 PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,PTYSA2,PLAQ2,EAQ2,
     *                 AMCH1,AMCH2,IREJ,IKVALA,pttq1,ptta1,pttq2,ptta2)
C  SELECT PT VALUES FOR A TWO CHAIN SYSTEM
C                            SELECT SEA QUARK AND ANTIQUARK PT-VALUES
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      DIMENSION PQUAR(4),TQUAR(4),PAQUAR(4),TAQUAR(4)
      COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
        B33=1.55
      B33=4.0+3.0/LOG10(ECM+10.D0)
        IF (IKVALA.EQ.1)B33=8.0
        ICOUNT=0
        IREJ=0
    1   CONTINUE
        B33=2.*B33
        ICOUNT=ICOUNT+1
        IF (ICOUNT.EQ.4)THEN
          IREJ=1
C                            REJECT EVENT
          RETURN
        ENDIF
        ES=-2./(B33**2)*LOG(ABS(RNDM(V)*RNDM(U))+1.E-18)
        HPS=SQRT(ES*ES+2.*ES*0.94)
        CALL DSFECF(SFE,CFE)
        PTXSQ1=HPS*CFE+PQUAR(1)
        PTYSQ1=HPS*SFE+PQUAR(2)
        PTXSA1=-PTXSQ1+PAQUAR(1)
        PTYSA1=-PTYSQ1+PAQUAR(2)
        ES=-2./(B33**2)*LOG(ABS(RNDM(V)*RNDM(U))+1.E-18)
        HPS=SQRT(ES*ES+2.*ES*0.94)
        CALL DSFECF(SFE,CFE)
        PTXSQ2=HPS*CFE+TQUAR(1)
        PTYSQ2=HPS*SFE+TQUAR(2)
        PTXSA2=-PTXSQ2+TAQUAR(1)
        PTYSA2=-PTYSQ2+TAQUAR(2)
        IF (IOUXEV.GE.6)WRITE(6,115)PTXSQ1,PTYSQ1,PTXSA1,PTYSA1
     *                             ,PTXSQ2,PTYSQ2,PTXSA2,PTYSA2
  115   FORMAT (' PT S  ',8F12.6)
C                           KINEMATICS OF THE TWO CHAINS Q1-AQ2,AQ1-Q2
        PTTQ1=PTXSQ1**2+PTYSQ1**2
        PTTA1=PTXSA1**2+PTYSA1**2
        PTTQ2=PTXSQ2**2+PTYSQ2**2
        PTTA2=PTXSA2**2+PTYSA2**2
        EQ1=PQUAR(4)
        EAQ1=PAQUAR(4)
        EQ2=TQUAR(4)
        EAQ2=TAQUAR(4)
        IF((EQ1**2.LE.PTTQ1).OR.(EQ2**2.LE.PTTQ2)
     *           .OR.(EAQ1**2.LE.PTTA1).OR.(EAQ2**2.LE.PTTA2))THEN
          GO TO 1
        ENDIF
        PLQ1=SQRT(EQ1**2-PTTQ1+1.E-6)*PQUAR(3)/ABS(PQUAR(3))
        PLQ2=SQRT(EQ2**2-PTTQ2+1.E-6)*TQUAR(3)/ABS(TQUAR(3))
        PLAQ1=SQRT(EAQ1**2-PTTA1+1.E-6)*PAQUAR(3)/ABS(PAQUAR(3))
        PLAQ2=SQRT(EAQ2**2-PTTA2+1.E-6)*TAQUAR(3)/ABS(TAQUAR(3))
C                          CHAIN 1: Q1-AQ2     CHAIN2:  AQ1-Q2
        AMCH1=SQRT((EQ1+EAQ2)**2-(PTXSQ1+PTXSA2)**2
     *       -(PTYSQ1+PTYSA2)**2-(PLQ1+PLAQ2)**2)
        AMCH2=SQRT((EQ2+EAQ1)**2-(PTXSQ2+PTXSA1)**2
     *       -(PTYSQ2+PTYSA1)**2-(PLQ2+PLAQ1)**2)
        RETURN
        END
C
C****************************************************************8**
      SUBROUTINE XPTFL(NHARD,NSEA,IREG,XMAX1,XMAX2)
C  PARTON LEVEL COLLISION EVENTS (X,PT FLAVOR)
C  THE ROUTINE XPTFL CALLS ONE EVENT (in DTUJET)
C  in DTUNUC, DPMJET it calls the multiple soft and hard chains
C   in one elementary collision 
C
C   IJPROJ,IJTAR: PROJECTILE AND TARGET PARTICLE OF THE REACTION
C                 1=PROTON,  2=ANTIPROTON
C   IJPVAL,IJTVAL =0 VALENCE QUARKS OF PROJECTILE OR TARGET NOT INVOLVED
C                                                     IN HARD SCATTERING
C   IJPVAL,IJTVAL =1 VALENCE QUARKS OF PROJECTILE OR TARGET  INVOLVED
C                                                     IN HARD SCATTERING
C
C   PARTEV VERSIONS CONTROLLED BY NVERS
C              NVERS=1: ALL HARD PARTONS CONSIDERED TO BE GLUONS
C                       soft x values by rejection
C              NVERS=2: ALL HARD PARTONS CONSIDERED TO BE GLUONS
C                       soft x values by P.Aurenche P.Maire method
C
C THE RESULTS (HARD SCATTERING) ARE IN COMMON /ABRHRD/
C               XH1(I),XH2(I):     X-VALUES OF INITIAL PARTONS
C               IJHI1(I),IJHI2(I): FLAVOR OF INITIAL PARTON
C                                  0            GLUON
C                                  1,2          VALENCE U,D QUARKS
C                                  11,12,13,14  SEA UDSC-QUARKS
C                                  NEGATIVE     ANTI S OR V QUARKS
C               IJHF1(I),IJHF2(I): FLAVOR OF FINAL STATE PARTONS
C               PHARD1(I,J),PHARD2(I,J): FINAL PART. MOMENTUM AND ENERGY
C                                J=1   PX
C                                 =2   PY
C                                 =3   PZ
C                                 =4   ENERGY  (MASSLESS PARTONS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( ONE=1.D0,ONEH=.5D0, ZERO=0.D0)
      PARAMETER (UMMM=0.3D0)
      PARAMETER (SMMM=0.5D0)
      PARAMETER (CMMM=1.3D0)
      PARAMETER(MSH=250)
      COMMON /ABRHRD/XH1(MSH),XH2(MSH),IJHI1(MSH),IJHI2(MSH),
     *IJHF1(MSH),IJHF2(MSH),PHARD1(MSH,4),PHARD2(MSH,4)
      COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
c     addded as needed:
      COMMON /SINGDI/SILMSD,SIGDI
C repl:COMMON/DIFFRA/ISINGD,IDUBLD,SDFRAC
C repl:COMMON /COLLIS/ECM,S,IJPROJ,IJTAR,PTTHR,IOPHRD,IJPRLU,IJTALU,PTTHR2
C     COMMON /COLLIS/S,IJPROJ,IJTAR,PTTHR,IOPHRD,IJPRLU,IJTALU,PTTHR2
      COMMON/COLLIS/S,IJPROJ,IJTAR,PTTHR,PTTHR2,IOPHRD,IJPRLU,IJTALU
      CHARACTER*80 TITLE
      CHARACTER*8 PROJTY,TARGTY
C     COMMON /USER/TITLE,PROJTY,TARGTY,CMENER,ISTRUF
C    &            ,ISINGD,IDUBLD,SDFRAC,PTLAR
      COMMON /USER1/TITLE,PROJTY,TARGTY
      COMMON /USER2/CMENER,SDFRAC,PTLAR,ISTRUF,ISINGD,IDUBLD
*
      COMMON /COLLE/ NEVHAD,NVERS,IHADRZ,NFILE
      COMMON /POMTYP/IPIM,ICON,ISIG,LMAX,MMAX,NMAX,DEFEL,DIFNU
      common/IPIMM/IPIMO
C-------------------------------------------------------------

      PARAMETER (INTMX=2488)
      COMMON /ABRSOF/XSQ1(INTMX),XSAQ1(INTMX),XSQ2(INTMX),XSAQ2(INTMX),
     *        IJSQ1(INTMX),IJSAQ1(INTMX),IJSQ2(INTMX),IJSAQ2(INTMX),
     *        AMCCH1(INTMX),AMCCH2(INTMX),GAMCH1(INTMX),GAMCH2(INTMX),
     *        BGCH1(INTMX),BGCH2(INTMX),THECH1(INTMX),THECH2(INTMX),
     *        BGXCH1(INTMX),BGYCH1(INTMX),BGZCH1(INTMX),
     *        BGXCH2(INTMX),BGYCH2(INTMX),BGZCH2(INTMX),
     *        NCH1(INTMX),NCH2(INTMX),IJCH1(INTMX),IJCH2(INTMX),
     *  PSOFA1(INTMX,4),PSOFA2(INTMX,4),PSOFB1(INTMX,4),PSOFB2(INTMX,4)
     * ,JHKKPZ(INTMX),JHKKTZ(INTMX),JHKKSX(INTMX),JHKKS1(INTMX)
      COMMON /ABRJT/XJQ1(INTMX),XJAQ1(INTMX),XJQ2(INTMX),XJAQ2(INTMX),
     *        IJJQ1(INTMX),IJJAQ1(INTMX),IJJQ2(INTMX),IJJAQ2(INTMX),
     *        AMJCH1(INTMX),AMJCH2(INTMX),GAMJH1(INTMX),GAMJH2(INTMX),
     *        BGJH1(INTMX),BGJH2(INTMX),THEJH1(INTMX),THEJH2(INTMX),
     *        BGXJH1(INTMX),BGYJH1(INTMX),BGZJH1(INTMX),
     *        BGXJH2(INTMX),BGYJH2(INTMX),BGZJH2(INTMX),
     *  PJETA1(INTMX,4),PJETA2(INTMX,4),PJETB1(INTMX,4),PJETB2(INTMX,4)
     * ,JHKKPH(INTMX),JHKKTH(INTMX),JHKKEX(INTMX),JHKKE1(INTMX)
      COMMON /NUCJTN/NONUJ1,NONUJT,NONUS1,NONUST
C-------------------------------------------------------------
C     COMMON /SVSWAP/ISVSWP,ISVSWT,JSVSWP,JSVSWT,XPSVSW,XTSVSW
      COMMON /VALHVG/PHPVAL(4),PHTVAL(4),IJVGP,IJVGT,IVALHP,IVALHT
      COMMON /PTLARG/ XSMAX
      COMMON /GLUSPL/NUGLUU,NSGLUU
      COMMON /SEASU3/SEASQ
      COMMON/VVDIFF/NVALCH,NVALDI,NSOFVD,IDIFTP,AMCHDD,NVADUD
      COMMON/INTNEZ/NDZ,NZD
C
C  THE COMMON BLOCK /PART/ DEFINES THE PARTICLE PROPERTIES AS USED IN
C  BAMJET AND DECAY
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
C  THE COMMON BLOCK /INPDAT/ DEFINES QUANTITIES NEEDED IN THE BAMJET
C  CHAIN DECAY CODE
*KEEP,DINPDA.
      COMMON /DINPDA/ IMPS(6,6),IMVE(6,6),IB08(6,21),IB10(6,21), IA08
     +(6,21),IA10(6,21), A1,B1,B2,B3,LT,LE,BET,AS,B8,AME,DIQ,ISU
      COMMON /LMMAXI/ LMMAX
      PARAMETER (NSTRMX=50)
      COMMON/SKINE1/NGLUEF,NNN,GL(NSTRMX),GR(NSTRMX),VL,VR,WL,WR,
     *              PTGL(2,NSTRMX),PTVL(2),PTWL(2),
     *              PTGR(2,NSTRMX),PTVR(2),PTWR(2)
      COMMON /DROPJJ/DROPJT,DROPVA
C     COMMON /XYTES/XTEST(50),XYTEST(0:11,50)
      DATA NCMPO/0/
      DATA INICHA/0/
C     COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
      IF(XMAX1.LE.0.D0.OR.XMAX2.LE.0.D0)THEN
      WRITE(6,'(A,3I5,2F10.4)')' XPTFL(',NHARD,NSEA,IREG,XMAX1,XMAX2
	NHARD=0
	NSEA=0
	IREG=0
	RETURN
      ENDIF
      IJPVAL=0
      IJTVAL=0
      NONUJ1=NONUJT+1
      NONUS1=NONUST+1
      IREG=0
      IOUXEV=IPEV
      IOUTPA=IPPA
      IOPTPO=IPRI
      IOUCOl=IPCO
C     NDZ=0
C     NZD=0
C     to keep identical commons
      ECM=CMENER
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C                     Initialize Charm selection at hard chain ends
C
      IF(INICHA.EQ.0)THEN
        INICHA=1
        PCCC=0.333*(UMMM/(CMMM*LOG(CMMM/0.2)))**2
        WRITE(6,4567)PCCC
 4567   FORMAT(' Charm at hard chain ends XPTFL: PCCC ',1F10.5)
      ENDIF
C
C----------------------------------------------------------------------
*
      ZXC=5000.
C  CALL NSOFT-NHARD EVENT
      NHARD=0
      NSEA=0
      NVAL=0
      IF (IOUXEV.GE.6)WRITE (6,107)IOUXEV,NHARD,NSEA,NVAL,NVERS
  107 FORMAT ('  XPTFL  IOUXEV,NHARD,NSEA,NVAL,NVERS = ',6I10)
      NC1000=0
 1000 CONTINUE
      DROPVA=0.
      NC1000=NC1000+1
      IF (IOUXEV.GE.1.AND.MOD(NC1000,20).EQ.0)WRITE(6,1100)NC1000
 1100 FORMAT(' REJECTION IN XPTFL ',I10)
      NNPO=0
      IF (IPIM.NE.2)THEN
        CALL SAMPLM(LPO,MPO,NPO)
        NPOLO=0
        NPODD=0
      ELSEIF(IPIM.EQ.2)THEN
      IF (IOUXEV.GE.6)WRITE(6,'(A)')' XPTFL call SAMPLX'

        CALL SAMPLX(LPO,MPO,NPO,NPODD,NPOLO)
        NCMPO=NCMPO+MPO
C       WRITE(6,*) ' NCMPO,MPO', NCMPO,MPO
      ENDIF
C     EACH HARD SCATTERING ALSO GETS SOFT COLOR MIXUP
      LPO = LPO + MPO
      IF (IOUXEV.GE.6)WRITE(6,107)IOUXEV,NHARD,NSEA,NVAL,NVERS
      NHARD=MPO
      IF (IOUXEV.GE.1)WRITE(6,101)LPO,MPO,NPO ,NNPO
  101 FORMAT (' XPTFL SAMPLM-LPO,MPO,NPO,NNPO= ',4I10
     *          /' NEXT CALL SELHRD')
C  CALL HARD PARTON EVENT
      NAX12=0
 2628   CONTINUE
      HAX1=0.
      HAXX1=0.
      HAX2=0.
      HAXX2=0.
 7717 CONTINUE
      IF(MPO.GE.1)THEN
      IF (IOUXEV.GE.6)WRITE(6,'(A)')' XPTFL call SELHRD'
        CALL SELHRD(MPO,IJPVAL,IJTVAL,PTTHR2)
        NHARD = MPO
 7727   CONTINUE
        IF (IOUXEV.GE.6)WRITE (6,107)IOUXEV,NHARD,NSEA,NVAL,NVERS
        DO 10 I=1,MPO
          HAX1=HAX1+XH1(I)
          HAX2=HAX2+XH2(I)
          IF(IOUXEV.GE.1.AND.XH1(I).LT.0 )WRITE(6,7788) I,XH1(I)
          IF(IOUXEV.GE.1.AND.XH2(I).LT.0 )WRITE(6,7787) I,XH2(I)
 7788     FORMAT(' XPTFL: XH1(',I5,') =',E12.5)
 7787     FORMAT(' XPTFL: XH2(',I5,') =',E12.5)
   10   CONTINUE
      ENDIF
      SOXM=0.
      SOX1=XMAX1-HAX1
      SOX2=XMAX2-HAX2
      IF (SOX1.LT.SOXM .OR. SOX2.LT.SOXM)THEN
        IF (IOUTPA.GE.1)WRITE (6,2510)HAX1,HAX2,XMAX1,XMAX2,MPO
 2510   FORMAT(' REJECT HAX1,HAX2.GT.1 HAX1,HAX2,XMAX1,XMAX2MPO='
     *         ,4F10.3,I10)
        NAX12=NAX12+1
        IF(MOD(NAX12,2).EQ.0)THEN
          GO TO 1000
        ELSE
          GO TO 2628
        ENDIF
      ENDIF
      NC1002=0
      GO TO 1002
 1001 CONTINUE
C     NDZ=0
C     NZD=0
      IF (LPO.GT.1) THEN
        IF(MOD(NC1002,6).EQ.0) THEN
          LPO=LPO-1
        ELSEIF( NC1002.GT.50 ) THEN
          IF(IOUXEV.GE.3)WRITE(6,9874) MPO
 9874     FORMAT(' XPTFL: 1001 SOFT X REJECTION TO 1000, MPO=',I5)
          GO TO 1000
        ENDIF
        NC1002=NC1002 + 1
      ENDIF
 1002 CONTINUE
      SOXUS1=0
      SOXUS2=0
      IF (IOUXEV.GE.3)WRITE (6,105)SOXUS1,SOXUS2,SOX1,SOX2,HAX1,HAX2
  105     FORMAT('XPTFL SOXUS1,SOXUS2,SOX1,SOX2,HAX1,HAX2 ',6F10.6)
      NHARD=MPO
      LPASOF=0
      IF (LPO.GT.1)THEN
        IF (NVERS.EQ.1) THEN
          UNOGLU=5.
          UNOVAL=2.
      IF (IOUXEV.GE.6)WRITE(6,*)' XPTFL call XPTFL1,NSEA,NVAL',NSEA,NVAL
          CALL XPTFL1(NHARD,NSEA,NVAL,SOXUS1,SOXUS2,SOX1,SOX2,HAX1,HAX2,
     *          LPO,MPO,NPO,LPASOF,IJPVAL,IJTVAL,RJ1000,XMAX1,XMAX2)
          IF (RJ1000.EQ.1.D0) THEN
           IF (IOUXEV.GE.6) THEN
            WRITE(6,*)'REJECTION TO 1001 AFTER XPTFL1 RJ1000=',RJ1000
           ENDIF
C	   IREG=1
C	   RETURN
           GO TO 1001
          ENDIF
        ENDIF
      ENDIF
 2020 CONTINUE
C
      NSEA=LPASOF
      IF (IOUXEV.GE.3)THEN
       WRITE (6,1303)NSEA
 1303  FORMAT ('  XPTFL (after xptfl1/2): NSEA=',I10,
     *'ii,ijsq1,ijsaq1,ijsq2,ijsaq2,amcch1,amcch2,...')
       DO 305 II=1,NSEA
        WRITE(6,304)II,
     *               IJSQ1(II),IJSAQ1(II),IJSQ2(II),IJSAQ2(II),
     *               AMCCH1(II),AMCCH2(II),GAMCH1(II),GAMCH2(II),
     *               BGXCH1(II),BGYCH1(II),BGZCH1(II),
     *               BGXCH2(II),BGYCH2(II),BGZCH2(II),
     *               NCH1(II),NCH2(II),IJCH1(II),IJCH2(II),
     *               (PSOFA1(II,JU),PSOFA2(II,JU),PSOFB1(II,JU),
     *               PSOFB2(II,JU),JU=1,4)
  304   FORMAT(5I4,6E18.8/4E18.8,4I4,2E18.8/7E18.8/7E18.8)
  305  CONTINUE
      ENDIF
C
C  END LOOP OVER SOFT SEA-SEA CHAINS--------------------------------
C
C                                X-VALUES REMAINING FOR VALENCE CHAINS
C
      SOXVA2=SOX2-SOXUS2
      SOXVA1=SOX1-SOXUS1
      IF(SOXVA1.LT.0.0D0.OR.SOXVA2.LT.0.0D0) THEN
       IF(IOUXEV.GE.6) THEN
        WRITE(6,*) '  XPTFL: REJECTION TO 1001 DUE TO SOXVA1/2 < 0.1'
     *  ,SOXVA1,SOXVA2
       ENDIF
       GOTO 1001
      ENDIF
C
C  PARTEV VERSIONS
C
      IF ((NVERS.EQ.1.OR.NVERS.EQ.2).AND.MPO.GE.1) THEN
C   PARTEV VERSION 1 : ALL HARD PARTONS CONSIDERED TO BE GLUONS
C                      IN AP-P ALSO VALENCE GLUON SCATTERING TREATED
C                      CHAINS FRAGMENTING:
C                          -SOFT VALENCE CHAINS
C                          -SOFT SEA CHAINS
C                          -SPLIT EACHHARD GLUON INTO Q-AQ PAIR
C                          -GLUON-GLUON BECOMES TWO Q-AQ CHAINS
C
        I=0
        NONUJY=NONUJT+1
        DO 301 IXNUJT=1,MPO
          I=I+1
          NONUJT=NONUJT+1
C                                 FIRST SPLIT GLUON MOMENTUM
          IC302=0
  302     CONTINUE
          IC302=IC302+1
          IF (IOUXEV.GE.3.AND.MOD(IC302,12).EQ.0)WRITE(6,1302)IC302
 1302 FORMAT(' REJECTION IN XPTFL 302 HARD GLUON SPLIT ',I10)
C                                 REJECT TOTAL EVENT FOR IC302=12
          IF (IC302.EQ.12) GO TO 1001
          XXXG1=(RNDM(V))**0.50
          XXXG2=(RNDM(U))**0.50
          IF (NUGLUU.EQ.1) THEN
            XXXG1=0.999999999999D0
            XXXG2=0.000000000001D0 
          ENDIF
          XJQ1(NONUJT)=XH1(I)
          XJQ2(NONUJT)=XH2(I)
          DO 303 J=1,3
            PJETA1(NONUJT,J)=PHARD1(I,J)*XXXG1
            PJETB1(NONUJT,J)=PHARD2(I,J)*XXXG2
            PJETA2(NONUJT,J)=PHARD2(I,J)*(1.-XXXG2)
            PJETB2(NONUJT,J)=PHARD1(I,J)*(1.-XXXG1)
  303     CONTINUE
          PJETA1(NONUJT,4)=SQRT(PJETA1(NONUJT,1)**2+
     *                          PJETA1(NONUJT,2)**2
     *                         +PJETA1(NONUJT,3)**2)
          PJETB1(NONUJT,4)=SQRT(PJETB1(NONUJT,1)**2+
     *                          PJETB1(NONUJT,2)**2
     *                         +PJETB1(NONUJT,3)**2)
          PJETA2(NONUJT,4)=SQRT(PJETA2(NONUJT,1)**2+
     *                          PJETA2(NONUJT,2)**2
     *                         +PJETA2(NONUJT,3)**2)
          PJETB2(NONUJT,4)=SQRT(PJETB2(NONUJT,1)**2+
     *                          PJETB2(NONUJT,2)**2
     *                         +PJETB2(NONUJT,3)**2)
C                                  MASSES OF SUBCHAINS
          AMJCH1(NONUJT)=SQRT((PJETA1(NONUJT,4)+
     *                         PJETA2(NONUJT,4))**2
     *                       -(PJETA1(NONUJT,1)+
     *                         PJETA2(NONUJT,1))**2
     *                       -(PJETA1(NONUJT,2)+
     *                         PJETA2(NONUJT,2))**2
     *                       -(PJETA1(NONUJT,3)+
     *                         PJETA2(NONUJT,3))**2)
          AMJCH2(NONUJT)=SQRT((PJETB1(NONUJT,4)+
     *                         PJETB2(NONUJT,4))**2
     *                       -(PJETB1(NONUJT,1)+
     *                         PJETB2(NONUJT,1))**2
     *                       -(PJETB1(NONUJT,2)+
     *                         PJETB2(NONUJT,2))**2
     *                       -(PJETB1(NONUJT,3)+
     *                         PJETB2(NONUJT,3))**2)
C                                  FLAVORS OF QUARKS
          AI=I
          BI=I+I
          IPJQ1=1.D0+RNDM(QA1)*(2.D0+SEASQ)
          IF(RNDM(V3).LT.PCCC)IPJQ1=4
          IPJAQ1=-IPJQ1
          IPJQ2=1.D0+RNDM(QB1)*(2.D0+SEASQ)
          IF(RNDM(V4).LT.PCCC)IPJQ2=4
          IPJAQ2=-IPJQ2
          IF (IOUXEV.GE.6)WRITE (6,113)IPJQ1,IPJQ2
  113     FORMAT(' IPJQ1,IPJQ2 ',2I10)
C                               REJECT SPLITTING FOR SMALL CHAIN MASSE
          IFPS1=IMPS(IPJQ2,IPJQ1)
          IFV1=IMVE(IPJQ2,IPJQ1)
          AMPS1=AAM(IFPS1)
          AMV1=AAM(IFV1)
          AMFF1=AMV1+0.3
C
          IFPS2=IMPS(IPJQ1,IPJQ2)
          IFV2=IMVE(IPJQ1,IPJQ2)
          AMPS2=AAM(IFPS2)
          AMV2=AAM(IFV2)
          AMFF2=AMV2+0.3
C
          IF(NUGLUU.EQ.0.AND.
     *    ((AMJCH1(NONUJT).LE.AMFF1).OR.
     *     (AMJCH2(NONUJT).LE.AMFF2))) GO TO 302
C
          GAMJH1(NONUJT)=(PJETA1(NONUJT,4)+
     *                    PJETA2(NONUJT,4))/AMJCH1(NONUJT)
          BGXJH1(NONUJT)=(PJETA1(NONUJT,1)+
     *                    PJETA2(NONUJT,1))/AMJCH1(NONUJT)
          BGYJH1(NONUJT)=(PJETA1(NONUJT,2)+
     *                    PJETA2(NONUJT,2))/AMJCH1(NONUJT)
          BGZJH1(NONUJT)=(PJETA1(NONUJT,3)+
     *                    PJETA2(NONUJT,3))/AMJCH1(NONUJT)
          GAMJH2(NONUJT)=(PJETB1(NONUJT,4)+
     *                    PJETB2(NONUJT,4))/AMJCH2(NONUJT)
          BGXJH2(NONUJT)=(PJETB1(NONUJT,1)+
     *                    PJETB2(NONUJT,1))/AMJCH2(NONUJT)
          BGYJH2(NONUJT)=(PJETB1(NONUJT,2)+
     *                    PJETB2(NONUJT,2))/AMJCH2(NONUJT)
          BGZJH2(NONUJT)=(PJETB1(NONUJT,3)+
     *                    PJETB2(NONUJT,3))/AMJCH2(NONUJT)
          IJJQ1 (NONUJT)=IPJQ1
          IJJAQ1(NONUJT)=IPJAQ1
          IJJQ2 (NONUJT)=IPJQ2
C         CHange r.e.21.4.94 flavor conservation
C         IJJAQ2(NONUJT)=IPJAQ2
          IJJAQ2(NONUJT)=-IPJQ1
  301   CONTINUE
  403 FORMAT (I10)
      DO 405 II=NONUJY,NONUJT
      IF (IOUTPA.GE.3)
     *  WRITE(6,404)II,
     *               IJJQ1(II),IJJAQ1(II),IJJQ2(II),IJJAQ2(II),
     *               AMJCH1(II),AMJCH2(II),GAMJH1(II),GAMJH2(II),
     *               BGXJH1(II),BGYJH1(II),BGZJH1(II),
     *               BGXJH2(II),BGYJH2(II),BGZJH2(II),
     *               (PJETA1(II,JU),PJETA2(II,JU),PJETB1(II,JU),
     *               PJETB2(II,JU),JU=1,4)
  404   FORMAT(5I4,6E18.8/4E18.8,2E18.8/7E18.8/7E18.8)
  405 CONTINUE
      ENDIF
      IF (IOUXEV.GE.6)WRITE (6,107)IOUXEV,NHARD,NSEA,NVAL,NVERS
      RETURN
      END
C-----------------------------------------------------------------     
C-----------------------------------------------------------------     
C-----------------------------------------------------------------     
      SUBROUTINE XPTFL1(NHARD,NSEA,NVAL,SOXUS1,SOXUS2,SOX1,SOX2,HAX1,
     *      HAX2,LPO,MPO,NPO,LPASOF,IJPVAL,IJTVAL,RJ1000,XMAX1,XMAX2)
C  PARTON LEVEL COLLISION EVENTS (X,PT FLAVOR)
C  THE ROUTINE XPTFL CALLS ONE EVENT
C
C   IJPROJ,IJTAR: PROJECTILE AND TARGET PARTICLE OF THE REACTION
C                 1=PROTON,  2=ANTIPROTON
C   IJPVAL,IJTVAL =0 VALENCE QUARKS OF PROJECTILE OR TARGET NOT INVOLVED
C                                                     IN HARD SCATTERING
C   IJPVAL,IJTVAL =1 VALENCE QUARKS OF PROJECTILE OR TARGET  INVOLVED
C                                                     IN HARD SCATTERING
C
C   PARTEV VERSIONS CONTROLLED BY NVERS
C              NVERS=1: ALL HARD PARTONS CONSIDERED TO BE GLUONS
C                       soft x values by rejection
C              NVERS=2: ALL HARD PARTONS CONSIDERED TO BE GLUONS
C                       soft x values by P.Aurenche P.Maire method
C
C THE RESULTS (HARD SCATTERING) ARE IN COMMON /ABRHRD/
C               XH1(I),XH2(I):     X-VALUES OF INITIAL PARTONS
C               IJHI1(I),IJHI2(I): FLAVOR OF INITIAL PARTON
C                                  0            GLUON
C                                  1,2          VALENCE U,D QUARKS
C                                  11,12,13,14  SEA UDSC-QUARKS
C                                  NEGATIVE     ANTI S OR V QUARKS
C               IJHF1(I),IJHF2(I): FLAVOR OF FINAL STATE PARTONS
C               PHARD1(I,J),PHARD2(I,J): FINAL PART. MOMENTUM AND ENERGY
C                                J=1   PX
C                                 =2   PY
C                                 =3   PZ
C                                 =4   ENERGY  (MASSLESS PARTONS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER ( ONE=1.D0,ONEH=.5D0, ZERO=0.D0)
      PARAMETER (UMMM=0.3D0)
      PARAMETER (SMMM=0.5D0)
      PARAMETER (CMMM=1.3D0)
      PARAMETER(MSH=250)
      PARAMETER (INTMD=252)
*KEEP,NUCC.
      COMMON /NUCC/   IT,ITZ,IP,IPZ,MJPROJ,IBPROJ,IJTARG,IBTARG
      COMMON /ABRHRD/XH1(MSH),XH2(MSH),IJHI1(MSH),IJHI2(MSH),
     *IJHF1(MSH),IJHF2(MSH),PHARD1(MSH,4),PHARD2(MSH,4)
      COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
      COMMON/INTNEZ/NDZ,NZD
C repl:COMMON /COLLIS/ECM,S,IJPROJ,IJTAR,PTTHR,IOPHRD,IJPRLU,IJTALU,PTTHR2
C     COMMON /COLLIS/S,IJPROJ,IJTAR,PTTHR,IOPHRD,IJPRLU,IJTALU,PTTHR2
      COMMON/COLLIS/S,IJPROJ,IJTAR,PTTHR,PTTHR2,IOPHRD,IJPRLU,IJTALU
      CHARACTER*80 TITLE
      CHARACTER*8 PROJTY,TARGTY
C     COMMON /USER/TITLE,PROJTY,TARGTY,CMENER,ISTRUF
C    &            ,ISINGD,IDUBLD,SDFRAC,PTLAR
      COMMON /USER1/TITLE,PROJTY,TARGTY
      COMMON /USER2/CMENER,SDFRAC,PTLAR,ISTRUF,ISINGD,IDUBLD
*
      COMMON /COLLE/ NEVHAD,NVERS,IHADRZ,NFILE
      COMMON /IPIMM/IPIM
      COMMON /SEASU3/SEASQ
      COMMON /POMTYP/ IPOM2,IPOM1,IPOSOM(4),APOSOM(2)
      COMMON /DIQUAX/AMEDD,IDIQUA,IDIQUU
C-------------------------------------------------------------

      PARAMETER (INTMX=2488)
      COMMON /ABRSOF/XSQ1(INTMX),XSAQ1(INTMX),XSQ2(INTMX),XSAQ2(INTMX),
     *        IJSQ1(INTMX),IJSAQ1(INTMX),IJSQ2(INTMX),IJSAQ2(INTMX),
     *        AMCCH1(INTMX),AMCCH2(INTMX),GAMCH1(INTMX),GAMCH2(INTMX),
     *        BGCH1(INTMX),BGCH2(INTMX),THECH1(INTMX),THECH2(INTMX),
     *        BGXCH1(INTMX),BGYCH1(INTMX),BGZCH1(INTMX),
     *        BGXCH2(INTMX),BGYCH2(INTMX),BGZCH2(INTMX),
     *        NCH1(INTMX),NCH2(INTMX),IJCH1(INTMX),IJCH2(INTMX),
     *  PSOFA1(INTMX,4),PSOFA2(INTMX,4),PSOFB1(INTMX,4),PSOFB2(INTMX,4)
     * ,JHKKPZ(INTMX),JHKKTZ(INTMX),JHKKSX(INTMX),JHKKS1(INTMX)
      COMMON /ABRJT/XJQ1(INTMX),XJAQ1(INTMX),XJQ2(INTMX),XJAQ2(INTMX),
     *        IJJQ1(INTMX),IJJAQ1(INTMX),IJJQ2(INTMX),IJJAQ2(INTMX),
     *        AMJCH1(INTMX),AMJCH2(INTMX),GAMJH1(INTMX),GAMJH2(INTMX),
     *        BGJH1(INTMX),BGJH2(INTMX),THEJH1(INTMX),THEJH2(INTMX),
     *        BGXJH1(INTMX),BGYJH1(INTMX),BGZJH1(INTMX),
     *        BGXJH2(INTMX),BGYJH2(INTMX),BGZJH2(INTMX),
     *  PJETA1(INTMX,4),PJETA2(INTMX,4),PJETB1(INTMX,4),PJETB2(INTMX,4)
     * ,JHKKPH(INTMX),JHKKTH(INTMX),JHKKEX(INTMX),JHKKE1(INTMX)
*KEEP,ABRZD.
      COMMON /ABRZD/ AMCZD1(INTMD),AMCZD2(INTMD),
     +GACZD1(INTMD),GACZD2(INTMD),
     +BGXZD1(INTMD),BGYZD1(INTMD),BGZZD1(INTMD), 
     +BGXZD2(INTMD),BGYZD2(INTMD),
     +BGZZD2(INTMD), NCHZD1(INTMD),NCHZD2(INTMD),
     +IJCZD1(INTMD),IJCZD2(INTMD),
     +PQZDA1(INTMD,4),PQZDA2(INTMD,4), PQZDB1(INTMD,4),
     +PQZDB2(INTMD,4),
     +IPCQ(INTMD),ITCQ(INTMD),ITCQ2(INTMD),IPCAQ(INTMD),
     +ITCAQ(INTMD),ITCAQ2(INTMD)
     +,IZDSS(INTMD)
*KEEP,ABRDZ.
      COMMON /ABRDZ/ AMCDZ1(INTMD),AMCDZ2(INTMD),
     +GACDZ1(INTMD),GACDZ2(INTMD),
     +BGXDZ1(INTMD),BGYDZ1(INTMD),BGZDZ1(INTMD), 
     +BGXDZ2(INTMD),BGYDZ2(INTMD),
     +BGZDZ2(INTMD), NCHDZ1(INTMD),NCHDZ2(INTMD),
     +IJCDZ1(INTMD),IJCDZ2(INTMD),
     +PQDZA1(INTMD,4),PQDZA2(INTMD,4), PQDZB1(INTMD,4),
     +PQDZB2(INTMD,4),
     +IPZQ(INTMD),IPZQQ2(INTMD),ITZQ(INTMD),IPZAQ(INTMD),
     +IZAQQ2(INTMD),ITZAQ(INTMD)
     +,IDZSS(INTMD)
C-------------------
      COMMON /NUCJTN/NONUJ1,NONUJT,NONUS1,NONUST
C-------------------------------------------------------------
C     COMMON /SVSWAP/ISVSWP,ISVSWT,JSVSWP,JSVSWT,XPSVSW,XTSVSW
      COMMON /VALHVG/PHPVAL(4),PHTVAL(4),IJVGP,IJVGT,IVALHP,IVALHT
      COMMON /PTLARG/XSMAX
      COMMON /GLUSPL/NUGLUU,NSGLUU
C
C  THE COMMON BLOCK /PART/ DEFINES THE PARTICLE PROPERTIES AS USED IN
C  BAMJET AND DECAY
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

C  THE COMMON BLOCK /INPDAT/ DEFINES QUANTITIES NEEDED IN THE BAMJET
C  CHAIN DECAY CODE
*KEEP,DINPDA.
      COMMON /DINPDA/ IMPS(6,6),IMVE(6,6),IB08(6,21),IB10(6,21), IA08
     +(6,21),IA10(6,21), A1,B1,B2,B3,LT,LE,BET,AS,B8,AME,DIQ,ISU
      COMMON /LMMAXI/ LMMAX
      PARAMETER (NSTRMX=50)
      COMMON/SKINE1/NGLUEF,NNN,GL(NSTRMX),GR(NSTRMX),VL,VR,WL,WR,
     *              PTGL(2,NSTRMX),PTVL(2),PTWL(2),
     *              PTGR(2,NSTRMX),PTVR(2),PTWR(2)
      COMMON /PCHARM/PC
      COMMON /SEAQXX/ SEAQX,SEAQXN
C      DATA INIPRI/0/
       DATA INICHA/0/
       DATA JTSP /0/
C      to keep identical commons
       ECM=CMENER
*
      IF(IOUXEV.GE.4)WRITE(6,*)'XPTFL1:entry:NDZ,NZD,NNDZ,NNZD,NHARD,',
     *'NSEA,NVAL'
     *,NDZ,NZD,NNDZ,NNZD,NHARD,NSEA,NVAL
      NOSTIN=NONUST
      GO TO 1179
 1199 CONTINUE
      NDZ=NDZ-NNDZ
      NZD=NZD-NNZD
C                       change 27.10.96
C     NDZ=0
C     NZD=0
      IF (IOUXEV.GE.6)WRITE (6,*)'XPTFL1: 1199 ndz nzd nndz nnzd'
     *,NDZ,NZD,NNDZ,NNZD
      IF (IOUXEV.GE.6)WRITE (6,107)IOUXEV,NHARD,LPO,NZD,NDZ,LPASOF
 1179 CONTINUE
      NONUST=NOSTIN
C     NDZ=0
C     NZD=0
      NNDZ=0
      NNZD=0
      IPIM  =0
      RJ1000=0.D0
      NHARD=MPO
      UNOGLU=5.
      UNOVAL=2.
      SOX1=XMAX1-HAX1
      SOX2=XMAX2-HAX2
        LPASOF=0
C       LLPPOO=LPO-1+MPO
        LLPPOO=LPO-1
        IF (LLPPOO.LE.0) GO TO 2020
        ALPO=LPO
C       ALPO=LPO+MPO
        NSEA=0
          IF(IPOM1.EQ.48.AND.IPOM2.EQ.2.AND.ECM.LT.20.D0)THEN
            XPTHRO=1.5*LOG10(ECM/2000.)+5.
            XPTHRO=1.5*LOG10(ECM/200.)+3.5
          ELSEIF(IPOM1.EQ.48.AND.IPOM2.EQ.2.AND.ECM.GE.20.D0)THEN
            XPTHRO=5.0
          ENDIF
          IF(IPOM1.EQ.11.AND.IPOM2.EQ.5)XPTHRO=15.
          IF(IPOM1.EQ.5.AND.IPOM2.EQ.5)XPTHRO=20.
          IF (IPIM.EQ.2)XPTHRO=2.
	  XPTHRO=8.
	  IF(ISTRUF.EQ.15) XPTHRO=5.
	  IF(ISTRUF.EQ.22) XPTHRO=8.
          IF( JTSP.EQ.0 ) THEN
            WRITE(6,*)' XPTFL1: XPTHRO=',XPTHRO
            JTSP=1
          ENDIF
C                                 j.r.5.12.94
C         XPTHR=XPTHRO/ECM
          XPTHR=1.5*XPTHRO/(ECM**1.5*14.)
C------------------------------------------------------------
C                                  j.r.11.5.94---16.5.94
C	  IF(IP.EQ.1)XPTHR=1.5*XPTHRO/ECM**2
 	  IF(IP.EQ.1)XPTHR=1.5*XPTHRO/(ECM**1.5*14.)
C------------------------------------------------------------
          XPTHR2=1.8
          IF (XPTHR2.GT.XPTHRO)XPTHR2=XPTHRO
C                                     j.r,5.12.94
C         XSTHR2=XPTHR2/ECM
          XSTHR2=1.5*XPTHR2/(ECM**1.5*14.)
C------------------------------------------------------------
C                                  j.r.11.5.94---16.5.94
C	  IF(IP.EQ.1)XSTHR2=1.5*XPTHR2/ECM**2
 	  IF(IP.EQ.1)XSTHR2=1.5*XPTHR2/(ECM**1.5*14.)
C------------------------------------------------------------
C         IF(INIPRI.EQ.0)THEN
C           INIPRI=1
            ALOX1=LOG(SOX1/XPTHR)
            ALOX2=LOG(SOX2/XPTHR)
            ALOOX1=1.+ALOX1
            ALOOX2=1.+ALOX2
            ALOOO1=1./ALOOX1
            ALOOO2=1./ALOOX2
          IF( JTSP.EQ.1 ) THEN
            WRITE(6,9753)XPTHRO,XPTHR,XSTHR2
 9753       FORMAT(' XPTFL1: XPTHRO,XPTHR,XSTHR2= ',3E15.5)
            JTSP=2
          ENDIF
C         ENDIF
C       one pair of soft chains for each hard pomeron
C       IF(NPO.EQ.1)LLPPOO=LPO-2
C       IF(NPO.EQ.1)LLPPOO=LPO
C----------------------------------------------------------------------
C                     Initialize Charm selection at soft chain ends
C
      IF(INICHA.EQ.0)THEN
        GM=2.140
        X2=UMMM
	BETOO=7.5D0
      ENDIF
        RX=XPTHRO
        X1=RX
        BETCHA=BETOO+1.3-LOG10(ECM)
        PU=DBETA(X1,X2,BETCHA)
        X2=SMMM
        PS=DBETA(X1,X2,BETCHA)
        X2=CMMM
        PC=DBETA(X1,X2,BETCHA)
C       PU1=PU/(2*PU+PS+PC)
C       PS1=PS/(2*PU+PS+PC)
        PC1=PC/(2*PU+PS+PC)
        PC=PC1
        PU1=PU/(2*PU+PS+PC)
        PS1=PS/(2*PU+PS+PC)
      IF(INICHA.EQ.0)THEN
        INICHA=1
        WRITE(6,4567)PC,BETCHA,PU1,PS1
 4567   FORMAT(' Charm at chain ends XPTFL1: PC,BETCHA,PU,PS ',4F10.5)
      ENDIF
C----------------------------------------------------------------------
C
        DO 20 I=1,LLPPOO
C         JSVSWP=0
C         JSVSWT=0
          AI=I-1
C         XPTHRX=XPTHR-0.5*AI/ECM
          XPTHRX=XPTHR-0.5*AI/ECM**2
C------------------------------------------------------------
C                                  j.r.11.5.94---16.5.94
 	  IF(IP.EQ.1)XPTHRX=XPTHR-0.5*AI/ECM**2
C------------------------------------------------------------
C         IF (XPTHRX.LT.2.D0/ECM)XPTHRX=2./ECM
          IF (XPTHRX.LT.4.D0/ECM**2)XPTHRX=4./ECM**2
C------------------------------------------------------------
C                                  j.r.11.5.94---16.5.94
          IF(IP.EQ.1.AND.XPTHRX.LT.4.D0/ECM**2)XPTHRX=4./ECM**2
C------------------------------------------------------------
C
C  LOOP OVER (LPO-1) SOFT SEA-SEA CHAINS WITH GLUONS AT CHAIN ENDS
C  GLUONS WILL BE SPLIT INTO QUARK-ANTIQUARK PAIRS
C  AND TWO Q-AQ CHAINS FORMED PER COLLISION
C
C                            GLUON X-VALUES
C                                             CHANGE J.R.21.5.90
C---------------------------------------------------------------
          NCOGLU=0
 5577     CONTINUE
          NCOGLU=NCOGLU+1
	  IF(NCOGLU.GE.6)THEN
C            REJECT XGLU values too large	    
C                                 REJECT THE TOTAL EVENT
            IF (IOUXEV.GE.6)WRITE (6,*)' REJECT  EVENT XGLU-VALUES'
            LPO=LPO-1
            SOXUS1=0.
            SOXUS2=0.
            GO TO 1199
	  ENDIF
          IF (RNDM(V1).LT.ALOOO1)THEN
            XGLU1=RNDM(A2)*(XPTHRX-XSTHR2)+XSTHR2
          ELSE
   25       CONTINUE
C................................................................a
	    IF(SEAQX.LE.0.75D0)THEN
              XGLU1=SAMPEX(XPTHRX,SOX1)
	    ELSEIF(SEAQX.GT.0.75D0)THEN
              XGLU1=SAMPEY(XPTHRX,SOX1)
            ENDIF
C................................................................
          ENDIF
          IF (RNDM(V3).LT.ALOOO2)THEN
            XGLU2=RNDM(A4)*(XPTHRX-XSTHR2)+XSTHR2
          ELSE
   26       CONTINUE
C................................................................
	    IF(SEAQX.LE.0.75D0)THEN
              XGLU2=SAMPEX(XPTHRX,SOX1)
	    ELSEIF(SEAQX.GT.0.75D0)THEN
              XGLU2=SAMPEY(XPTHRX,SOX1)
            ENDIF
C................................................................
C                                   CHANGE 18.6.90 PREVENT EVENT LOOP
          ENDIF
C---------------------------------------------------------------
C                              SPLIT GLUON INTO TWO SEA QUARKS
C                              FLAVORS OF SEA QUARKS
          IF(IOUXEV.GE.6)WRITE (6,109) XGLU1,XGLU2
C                 Are these xglu values allowed
          IF(XGLU1+SOXUS1.GT.SOX1.OR.XGLU2+SOXUS2.GT.SOX2)GO TO 5577
  109     FORMAT (' XPTFL1  XGLU1,XGLU2 ',2F10.6)
          AI=I
          BI=I+I
          IPSQ1=1.D0+RNDM(QA1)*(2.D0+SEASQ)
          IF(RNDM(W1).LT.PC)IPSQ1=4
          IPSAQ1=-IPSQ1
          IPSQ2=1.D0+RNDM(QB1)*(2.D0+SEASQ)
          IF(RNDM(W2).LT.PC)IPSQ2=4
          IPSAQ2=-IPSQ2
          IF (IOUXEV.GE.6)WRITE (6,113)IPSQ1,IPSQ2
  113     FORMAT(' XPTFL1  IPSQ1,IPSQ2 ',2I10)
C                              X-FRAXTIONS OF SEA QUARKS
C------------------------------------------------------j.r.29.4.93
          IF(IPSQ1.LE.2)THEN
            XPSQ1=(0.2+(0.36*RNDM(A1))**0.50)*XGLU1
            XPSAQ1=XGLU1-XPSQ1
          ELSEIF(IPSQ1.EQ.3)THEN
            BSQ=0.7/ECM
            XSTHR=2./ECM
	    ICOXS1=0
 5588       CONTINUE	    
	    ICOXS1=ICOXS1+1
	    IF(ICOXS1.GT.8)THEN
C             REJECT XPSQ1 values too large	    
C                                 REJECT THE TOTAL EVENT
              IF (IOUXEV.GE.6)WRITE (6,*)' REJECT  EVENT XPSQ1-VALUES'
              LPO=LPO-1
              SOXUS1=0.
              SOXUS2=0.
	      IF(IOUXEV.GE.4) WRITE(6,*)' xptfl1 LPO,SOXUS1,SOXUS2 reject ',
     *	      LPO,SOXUS1,SOXUS2
              GO TO 1199
	    ENDIF
            XPSQ1=SAMPXB(XSTHR+BSQ,0.9D0,BSQ)
	    IF(XPSQ1.GE.XGLU1)GO TO 5588
	    XPSAQ1=XGLU1-XPSQ1
C           XPSAQ1=SAMPXB(XSTHR+BSQ,0.9D0,BSQ)
          ELSEIF(IPSQ1.EQ.4)THEN
            BCQ=2./ECM
            XSTHR=2./ECM
            XPSQ1=SAMPXB(XSTHR+BCQ,0.9D0,BCQ)
            XPSAQ1=SAMPXB(XSTHR+BCQ,0.9D0,BCQ)
          ENDIF
          IF(IPSQ2.LE.2)THEN
            XPSQ2=(0.2+(0.36*RNDM(B1))**0.50)*XGLU2
            XPSAQ2=XGLU2-XPSQ2
          ELSEIF(IPSQ2.EQ.3)THEN
            BSQ=0.7/ECM
            XSTHR=2./ECM
	    ICOXS2=0
 5599       CONTINUE	    
	    ICOXS2=ICOXS2+1
	    IF(ICOXS2.GT.8)THEN
C             REJECT XPSQ2 values too large	    
C                                 REJECT THE TOTAL EVENT
              IF (IOUXEV.GE.6)WRITE (6,*)' REJECT  EVENT XPSQ2-VALUES'
              LPO=LPO-1
              SOXUS1=0.
              SOXUS2=0.
              GO TO 1199
	    ENDIF
            XPSQ2=SAMPXB(XSTHR+BSQ,0.9D0,BSQ)
	    IF(XPSQ2.GE.XGLU2)GO TO 5599
	    XPSAQ2=XGLU2-XPSQ2
C           XPSAQ2=SAMPXB(XSTHR+BSQ,0.9D0,BSQ)
          ELSEIF(IPSQ2.EQ.4)THEN
            BCQ=2./ECM
            XSTHR=2./ECM
            XPSQ2=SAMPXB(XSTHR+BCQ,0.9D0,BCQ)
            XPSAQ2=SAMPXB(XSTHR+BCQ,0.9D0,BCQ)
          ENDIF
C------------------------------------------------------j.r.29.4.93
      IF (IOUXEV.GE.6)WRITE (6,107)IOUXEV,NHARD,LPO,NZD,NDZ,LPASOF
  107 FORMAT ('  XPTFL1  IOUXEV,NHARD,LPO,NZD,NDZ,LPASOF = ',6I10)
          IF(IOUXEV.GE.6)WRITE(6,114) XPSQ1,XPSAQ1,XPSQ2,XPSAQ2
  114     FORMAT('  XPSQ1,XPSAQ1,XPSQ2,XPSAQ2 ',4F12.6)
C  ------------------------------------------------------------------
C                        define eventually D-Z chains (sea-diquark--sea)
C
       IREJDZ=0
       IREJZD=0
       NDIQDZ=0
       NDIQZD=0
C      AME=0.9
       IF(RNDM(V).GT.2.D0*AMEDD-1.D0)THEN
         IF(IDIQUU.EQ.1)THEN
           IF(IOUXEV.GE.3)WRITE(6,*)' XPTFL1 call DIQDZZ ',
     *	   'LPO,AMEDD',LPO,AMEDD
           CALL DIQDZZ(ECM,XPSQ1,XPSAQ1,XPSQ2,XPSAQ2,IPSQ1,IPSAQ1,
     *               IPSQ2,IPSAQ2,IREJDZ)
           IF(IREJDZ.EQ.1)THEN
             IF (IOUXEV.GE.4)WRITE (6,'(2A,4I5)')'DIQDZZ1 ndz nzd nndz '
     *       ,'nnzd XPTFL1',NDZ,NZD,NNDZ,NNZD
	   ENDIF
           IF(IREJDZ.EQ.0) THEN
             NNDZ=NNDZ+1
             IF (IOUXEV.GE.3)WRITE (6,'(2A,4I5)')' DIQDZZ0 ndz nzd nndz'
     *       ,' nnzd XPTFL1',NDZ,NZD,NNDZ,NNZD
             NDIQDZ=1
C                                   TEST ARE THESE X VALUES ALLOWED
             SOXUS1=SOXUS1+XPSQ1+XPSAQ1
             SOXUS2=SOXUS2+XPSQ2+XPSAQ2
             IF(IOUXEV.GE.3)WRITE (6,*)' SOXUS1,SOXUS2,SOX1,SOX2 ',
     *	     'HAX1,HAX2 after call diqdzz ',
     *       SOXUS1,SOXUS2,SOX1,SOX2,
     *	     HAX1,HAX2
             IF ((SOXUS1.GT.SOX1).OR.(SOXUS2.GT.SOX2)) THEN
C                                   REJECT THE TOTAL EVENT
               IF (IOUXEV.GE.3)WRITE (6,106)
               RJ1000=1.D0
               NDZ=NDZ-NNDZ
               NZD=NZD-NNZD
C                          change 27.10.96
	       NNDZ=0
	       NNZD=0
C                          change 27.10.96
               NDIQDZ=0 
               LPO=LPO-1
               SOXUS1=0.
               SOXUS2=0.
               NONUST=NOSTIN
               IF (IOUXEV.GE.3)WRITE (6,*)' RETURN ndz nzd '
     *         ,'nndz,nnzd,LPO soxus.GT.sox  DIQDZZ0',
     *         NDZ,NZD,NNDZ,NNZD,LPO
               RETURN
             ENDIF
C            GO TO 20
           ENDIF
         ENDIF
       ENDIF
       IF(RNDM(V).GT.2.D0*AMEDD-1.D0.AND.NDIQDZ.EQ.0)THEN
         IF(IDIQUU.EQ.1)THEN
           IF(IOUXEV.GE.3)WRITE(6,*)' XPTFL1 call DIQZZD ',
     *	   'LPO,AMEDD',LPO,AMEDD
           CALL DIQZZD(ECM,XPSQ1,XPSAQ1,XPSQ2,XPSAQ2,IPSQ1,IPSAQ1,
     *               IPSQ2,IPSAQ2,IREJZD)
           IF(IREJZD.EQ.1)THEN
             IF (IOUXEV.GE.3)WRITE (6,'(2A,4I5)')' DIQZZD1 ndz nzd nndz'
     *       ,' nnzd XPTFL1',NDZ,NZD,NNDZ,NNZD
C	     NZD=NZD-1
	   ENDIF
           IF(IREJZD.EQ.0) THEN
             NNZD=NNZD+1
             IF (IOUXEV.GE.3)WRITE (6,'(2A,4I5)')' DIQZZD0 ndz nzd '
     *       ,'nndz,nnzd XPTFL1',NDZ,NZD,NNDZ,NNZD
             NDIQZD=1
C                                   TEST ARE THESE X VALUES ALLOWED
             SOXUS1=SOXUS1+XPSQ1+XPSAQ1
             SOXUS2=SOXUS2+XPSQ2+XPSAQ2
             IF(IOUXEV.GE.3)WRITE (6,*)' SOXUS1,SOXUS2,SOX1,SOX2 ,',
     *	     'HAX1,HAX2 after call diqzzd0',
     *       SOXUS1,SOXUS2,SOX1,SOX2,
     *	     HAX1,HAX2
             IF ((SOXUS1.GT.SOX1).OR.(SOXUS2.GT.SOX2)) THEN
C                                   REJECT THE TOTAL EVENT
               IF (IOUXEV.GE.3)WRITE (6,106)
               RJ1000=1.D0
               NZD=NZD-NNZD
               NDZ=NDZ-NNDZ
C                          change 27.10.96
C	       NDZ=0
C	       NZD=0
	       NNDZ=0
	       NNZD=0
C                          change 27.10.96
               NDIQZD=0
               LPO=LPO-1
               SOXUS1=0.
               SOXUS2=0.
               NONUST=NOSTIN
               IF (IOUXEV.GE.3)WRITE (6,*)' RETURN2 ndz nzd '
     *         ,'nndz,nnzd,LPO SOXUS.GT.SOX',
     *         'diqzzd0',NDZ,NZD,NNDZ,NNZD,LPO
               RETURN
             ENDIF
C            GO TO 20
           ENDIF
         ENDIF
       ENDIF
       AME=0.95
       PTXSQ1=0  
       PTYSQ1=0
       PTXSA1=0
       PTYSA1=0 
       PTXSQ2=0  
       PTYSQ2=0
       PTXSA2=0
       PTYSA2=0 
       PLQ1 = XPSQ1 *ECM/2.
       EQ1  = XPSQ1 *ECM/2. 
       PLAQ1= XPSAQ1*ECM/2.
       EAQ1 = XPSAQ1*ECM/2.
       PLQ2 =-XPSQ2 *ECM/2.
       EQ2  = XPSQ2 *ECM/2. 
       PLAQ2=-XPSAQ2*ECM/2.
       EAQ2 = XPSAQ2*ECM/2.
C  ------------------------------------------------------------------

C                            SELECT SEA QUARK AND ANTIQUARK PT-VALUES
C
          IKVALA=2
          NSELPT=1
       IF(IOUXEV.GE.6)WRITE(6,'(A)')' XPTFL1 call SELPT'
          CALL     SELPT(
     *                 PTXSQ1,PTYSQ1,PLQ1,EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1,
     *                 PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,PTYSA2,PLAQ2,EAQ2,
     *                 AMCH1,AMCH2,IREJ,IKVALA,pttq1,ptta1,PTTQ2,PTTA2,
     *                 NSELPT)
C
          IF (IREJ.EQ.1) THEN
           IF(IOUXEV.GE.6)WRITE(6,*)'  XPTFL1: --> 9922 IREJ=',IREJ
           IF(IOUXEV.GE.6)WRITE(6,'(A,6I5)')
     *' XPTFL1:NDZ,NNDZ,NDIQDZ,NZD,NNZD,NDIQZD ',
     *NDZ,NNDZ,NDIQDZ,NZD,NNZD,NDIQZD 
           IF(NDIQDZ.EQ.1)THEN
             NDZ=NDZ-1
             NDIQDZ=0
             NNDZ=NNDZ-1
           ENDIF
           IF(NDIQZD.EQ.1)THEN
             NZD=NZD-1
             NDIQZD=0
             NNZD=NNZD-1
           ENDIF
           GO TO 9922
          ENDIF
      IF (IOUXEV.GE.6)WRITE (6,*)'IOUXEV,NHARD,LPO,NZD,NDZ,LPASOF',
     *IOUXEV,NHARD,LPO,NZD,NDZ,LPASOF
C
C  REPLACE SMALL MASS CHAINS BY PSEUDOSCALAR OR VECTOR MESONS
C  FIRST FOR CHAIN 1
C
          IFPS1=IMPS(IPSQ2,IPSQ1)
          IFV1=IMVE(IPSQ2,IPSQ1)
          AMPS1=AAM(IFPS1)
          AMV1=AAM(IFV1)
          NNCH1=0
          AMFF1=AMV1+0.3
          IF(IOUXEV.GE.3)WRITE(6,102)AMCH1,AMPS1,AMV1,IFPS1,IFV1
  102     FORMAT(' AMCH1,AMPS1,AMV1,IFPS1,IFV1 ',3F12.4,2I10)
          IF(AMCH1.LT.AMFF1) THEN
           IF(IOUXEV.GE.6)WRITE(6,*)'  XPTFL1: --> 9922 AMCH1 < AMFF1'
           IF(IOUXEV.GE.6)WRITE(6,'(A,6I5)')
     *' XPTFL1:NDZ,NNDZ,NDIQDZ,NZD,NNZD,NDIQZD ',
     *NDZ,NNDZ,NDIQDZ,NZD,NNZD,NDIQZD 
           IF(NDIQDZ.EQ.1)THEN
             NDZ=NDZ-1
             NDIQDZ=0
             NNDZ=NNDZ-1
           ENDIF
           IF(NDIQZD.EQ.1)THEN
             NZD=NZD-1
             NDIQZD=0
             NNZD=NNZD-1
           ENDIF
           GO TO 9922
          ENDIF
          IF (AMCH1.LT.AMV1)THEN
C                                  PRODUCE PSEUDOSCALAR
            IJNCH1=IFPS1
            NNCH1=-1
C                                   CORRECT KINEMATICS
            XPSQ1=XPSQ1*AMPS1/AMCH1
            XPSAQ2=XPSAQ2*AMPS1/AMCH1
            AMCH1=AMPS1
C                                   GO TO REDO THE KINEMATICS
          ELSEIF(AMCH1.LT.AMFF1) THEN
C                                   PRODUCE VECTOR MESON
            IJNCH1=IFV1
            NNCH1=1
C                                   CORRECT KINEMATICS
            XPSQ1=XPSQ1*AMV1/AMCH1
            XPSAQ2=XPSAQ2*AMV1/AMCH1
            AMCH1=AMV1
C                                   GO TO REDO THE KINEMATICS
          ELSE
C                                   NO CORRECTIONS BUT DO  CHAIN 2
            GO TO 31
          ENDIF
C                                   CORRECT KINEMATICS FOR CHAIN 1

          EQ1=XPSQ1*ECM/2.
          EAQ2=XPSAQ2*ECM/2.
          IF(    (EQ1**2.LT.PTTQ1)
     *     .OR.(EAQ2**2.LT.PTTA2)) THEN
           IF(IOUXEV.GE.6)WRITE(6,*)'  XPTFL1: --> 9922 EQ^2 < PT'
     *     ,'EQ1 PTTQ1 EAQ2 PTTA2',EQ1,PTTQ1,EAQ2,PTTA2
           IF(IOUXEV.GE.6)WRITE(6,'(A,6I5)')
     *' XPTFL1:NDZ,NNDZ,NDIQDZ,NZD,NNZD,NDIQZD ',
     *NDZ,NNDZ,NDIQDZ,NZD,NNZD,NDIQZD 
           IF(NDIQDZ.EQ.1)THEN
             NDZ=NDZ-1
             NDIQDZ=0
             NNDZ=NNDZ-1
           ENDIF
           IF(NDIQZD.EQ.1)THEN
             NZD=NZD-1
             NDIQZD=0
             NNZD=NNZD-1
           ENDIF
           GO TO 9922
          ENDIF
          PLQ1=SQRT(EQ1**2-PTTQ1)
          PLAQ2=-SQRT(EAQ2**2-PTTA2)
   31     CONTINUE
      IF (IOUXEV.GE.6)WRITE (6,107)IOUXEV,NHARD,LPO,NZD,NDZ,LPASOF
C
C  REPLACE SMALL MASS CHAINS BY PSEUDOSCALAR OR VECTOR MESONS
C  SECOND FOR CHAIN 2
C
          IFPS2=IMPS(IPSQ1,IPSQ2)
          IFV2=IMVE(IPSQ1,IPSQ2)
          AMPS2=AAM(IFPS2)
          AMV2=AAM(IFV2)
          NNCH2=0
          AMFF2=AMV2+0.3
          IF(IOUXEV.GE.3)WRITE(6,103)AMCH2,AMPS2,AMV2,IFPS2,IFV2
  103     FORMAT(' AMCH2,AMPS2,AMV2,IFPS2,IFV2 ',3F12.4,2I10)
          IF(AMCH2.LT.AMFF2) THEN
           IF(IOUXEV.GE.6)WRITE(6,*)'  XPTFL1: --> 9922 AMCH2 < AMFF2'
           IF(IOUXEV.GE.6)WRITE(6,'(A,6I5)')
     *' XPTFL1:NDZ,NNDZ,NDIQDZ,NZD,NNZD,NDIQZD ',
     *NDZ,NNDZ,NDIQDZ,NZD,NNZD,NDIQZD 
           IF(NDIQDZ.EQ.1)THEN
             NDZ=NDZ-1
             NDIQDZ=0
             NNDZ=NNDZ-1
           ENDIF
           IF(NDIQZD.EQ.1)THEN
             NZD=NZD-1
             NDIQZD=0
             NNZD=NNZD-1
           ENDIF
           GO TO 9922
          ENDIF
          IF (AMCH2.LT.AMV2)THEN
C                                PRODUCE PSEUDO SCALAR
            IJNCH2=IFPS2
            NNCH2=-1
C                                   CORRECT KINEMATICS
            XPSQ2=XPSQ2*AMPS2/AMCH2
            XPSAQ1=XPSAQ1*AMPS2/AMCH2
            AMCH2=AMPS2
C                                   GO TO REDO THE KINEMATICS
          ELSEIF(AMCH2.LT.AMFF2) THEN
C                                   PRODUCE VECTOR MESON
            IJNCH2=IFV2
            NNCH2=1
C                                   CORRECT KINEMATICS
            XPSQ2=XPSQ2*AMV2/AMCH2
            XPSAQ1=XPSAQ1*AMV2/AMCH2
            AMCH2=AMV2
C                                   GO TO REDO THE KINEMATICS
          ELSE
C                                   NO CORRECTIONS
            GO TO 32
          ENDIF
C                                   CORRECT KINEMATICS FOR CHAIN 2

          EQ2=XPSQ2*ECM/2.
          EAQ1=XPSAQ1*ECM/2.
          IF(    (EQ2**2.LT.PTTQ2)
     *       .OR.(EAQ1**2.LT.PTTA1)) THEN
           IF(IOUXEV.GE.6)WRITE(6,*)'  XPTFL1: --> 9922 EQ^2 < PT'
     *     ,'EQ2 PTTQ2 EAQ1 PTTA1',EQ2,PTTQ2,EAQ1,PTTA1
           IF(IOUXEV.GE.6)WRITE(6,'(A,6I5)')
     *' XPTFL1:NDZ,NNDZ,NDIQDZ,NZD,NNZD,NDIQZD ',
     *NDZ,NNDZ,NDIQDZ,NZD,NNZD,NDIQZD 
           IF(NDIQDZ.EQ.1)THEN
             NDZ=NDZ-1
             NDIQDZ=0
             NNDZ=NNDZ-1
           ENDIF
           IF(NDIQZD.EQ.1)THEN
             NZD=NZD-1
             NDIQZD=0
             NNZD=NNZD-1
           ENDIF
           GO TO 9922
          ENDIF
          PLQ2=-SQRT(EQ2**2-PTTQ2)
          PLAQ1=SQRT(EAQ1**2-PTTA1)
   32     CONTINUE
C                                   TEST ARE THESE X VALUES ALLOWED
          IF(NDIQDZ.EQ.0.AND.NDIQZD.EQ.0)THEN
            SOXUS1=SOXUS1+XPSQ1+XPSAQ1
            SOXUS2=SOXUS2+XPSQ2+XPSAQ2
          ENDIF
          IF(IOUXEV.GE.3)WRITE (6,105)SOXUS1,SOXUS2,SOX1,SOX2,HAX1,HAX2
  105     FORMAT('XPTFL1 SOXUS1,SOXUS2,SOX1,SOX2,HAX1,HAX2 ',6F10.6)
          IF ((SOXUS1.GT.SOX1).OR.(SOXUS2.GT.SOX2)) THEN
C                                   REJECT THE TOTAL EVENT
            IF (IOUXEV.GE.6)WRITE (6,106)
  106       FORMAT(' REJECT THE EVENT  SEA X-VALUES')
C                                              j.r.10.5.94
            LPO=LPO-1
            SOXUS1=0.
            SOXUS2=0.
            GO TO 1199
          ENDIF
          GO TO 9923
 9922     CONTINUE
CC          LPO=LPO-1
            SOXUS1=0.
            SOXUS2=0.
            GO TO 1199
C                                   VALENCE-SEA SWAP-------------------
C         IF (JSVSWP.EQ.1)ISVSWP=0
C         IF (JSVSWT.EQ.1)ISVSWT=0
          GO TO 22
 9923     CONTINUE
C  NOW WE HAVE AN ACCEPTABLE SOFT GLUON-GLUON EVENT
C  AND PUT IT INTO THE HISTOGRAM
C
      IF (IOUXEV.GE.6)WRITE (6,107)IOUXEV,NHARD,LPO,NZD,NDZ,LPASOF
          LPASOF=LPASOF+1
          II=LPASOF
          NONUST=NONUST+1
          II=NONUST
          XSQ1(II)=XPSQ1
          XSAQ1(II)=XPSAQ1
          XSQ2(II)=XPSQ2
          XSAQ2(II)=XPSAQ2
          IJSQ1(II)=IPSQ1
          IJSAQ1(II)=IPSAQ1
          IJSQ2(II)=IPSQ2
          IJSAQ2(II)=IPSAQ2
          AMCCH1(II)=AMCH1
          AMCCH2(II)=AMCH2
          GAMCH1(II)=(EQ1+EAQ2)/AMCH1
          BGXCH1(II)=(PTXSQ1+PTXSA2)/AMCH1
          BGYCH1(II)=(PTYSQ1+PTYSA2)/AMCH1
          BGZCH1(II)=(PLQ1+PLAQ2)/AMCH1
          GAMCH2(II)=(EQ2+EAQ1)/AMCH2
          BGXCH2(II)=(PTXSQ2+PTXSA1)/AMCH2
          BGYCH2(II)=(PTYSQ2+PTYSA1)/AMCH2
          BGZCH2(II)=(PLAQ1+PLQ2)/AMCH2
          NCH1(II)=NNCH1
          NCH2(II)=NNCH2
          IF (IREJDZ.EQ.0.AND.NDIQDZ.EQ.1)THEN
            NCH1(II)=88
            NCH2(II)=88
          ENDIF
          IF (IREJZD.EQ.0.AND.NDIQZD.EQ.1)THEN
            NCH1(II)=88
            NCH2(II)=88
          ENDIF
          IF(NDIQDZ.EQ.1.AND.NDZ.GT.0)IDZSS(NDZ)=II
          IF(NDIQZD.EQ.1.AND.NZD.GT.0)IZDSS(NZD)=II
          IJCH1(II)=IJNCH1
          IJCH2(II)=IJNCH2
          PSOFA1(II,1)=PTXSQ1
          PSOFA1(II,2)=PTYSQ1
          PSOFA1(II,3)=PLQ1
          PSOFA1(II,4)=EQ1
          PSOFA2(II,1)=PTXSA2
          PSOFA2(II,2)=PTYSA2
          PSOFA2(II,3)=PLAQ2
          PSOFA2(II,4)=EAQ2
          PSOFB1(II,1)=PTXSQ2
          PSOFB1(II,2)=PTYSQ2
          PSOFB1(II,3)=PLQ2
          PSOFB1(II,4)=EQ2
          PSOFB2(II,1)=PTXSA1
          PSOFB2(II,2)=PTYSA1
          PSOFB2(II,3)=PLAQ1
          PSOFB2(II,4)=EAQ1
          IF (IOUXEV.GE.3)WRITE(6,104)II,
     *                 XSQ1(II),XSAQ1(II),XSQ2(II),XSAQ2(II),
     *                 IJSQ1(II),IJSAQ1(II),IJSQ2(II),IJSAQ2(II),
     *                 AMCCH1(II),AMCCH2(II),GAMCH1(II),GAMCH2(II),
     *                 BGCH1(II),BGCH2(II),THECH1(II),THECH2(II),
     *                 BGXCH1(II),BGYCH1(II),BGZCH1(II),
     *                 BGXCH2(II),BGYCH2(II),BGZCH2(II),
     *                 NCH1(II),NCH2(II),IJCH1(II),IJCH2(II),
     *               (PSOFA1(II,JU),PSOFA2(II,JU),PSOFB1(II,JU),
     *               PSOFB2(II,JU),JU=1,4)
  104     FORMAT(I10,4F12.7,4I5/10X,8F12.6/10X,6F12.6,4I5/8F15.5/8F15.5)
   22     CONTINUE
   20   CONTINUE
        IF (IOUXEV.GE.6) WRITE(6,*)'  LPASOF =',LPASOF
 2020   CONTINUE
      IF (IOUXEV.GE.4)WRITE (6,*)'END XPTFL1',
     * '  IOUXEV,NHARD,LPO,NZD,NDZ,LPASOF,IREJ',
     * IOUXEV,NHARD,LPO,NZD,NDZ,LPASOF,IREJ
      RETURN
      END
C*****************************************************************
      SUBROUTINE PTVAL(XP,XXP,XXT,XT,ECM,
     *                 PTXVQ1,PTYVQ1,PLQ1,EQ1,PTXVA1,PTYVA1,PLAQ1,EAQ1,
     *                 PTXVQ2,PTYVQ2,PLQ2,EQ2,PTXVA2,PTYVA2,PLAQ2,EAQ2,
     *                 AMCH1,AMCH2,IREJ,IKVALA)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      COMMON /COLLE/ NEVHAD,NVERS,IHADRZ,NFILE
      PARAMETER (NSTRMX=50)
      COMMON/SKINE1/NGLUEF,NNN,GL(NSTRMX),GR(NSTRMX),VL,VR,WL,WR,
     *              PTGL(2,NSTRMX),PTVL(2),PTWL(2),
     *              PTGR(2,NSTRMX),PTVR(2),PTWR(2)
      PTXVQ1 = PTVL(1)
      PTYVQ1 = PTVL(2)
      PTXVA1 = PTWL(1)
      PTYVA1 = PTWL(2)
      PTXVQ2 = PTVR(1)
      PTYVQ2 = PTVR(2)
      PTXVA2 = PTWR(1)
      PTYVA2 = PTWR(2)
      EQ1    = XP*ECM/2.
      EQ2    = XT*ECM/2.
      EAQ1   = XXP*ECM/2.
      EAQ2   = XXT*ECM/2.
C     PLQ1   = SQRT(EQ1**2-PTXVQ1**2-PTYVQ1**2)
C     PLQ2   = SQRT(EQ2**2-PTXVQ2**2-PTYVQ2**2)
C     PLAQ1  = SQRT(EAQ1**2-PTXVA1**2-PTYVA1**2)
C     PLAQ2  = SQRT(EAQ2**2-PTXVA2**2-PTYVA2**2)
      PLQ1   = EQ1
      PLQ2   = -EQ2
      PLAQ1  = EAQ1
      PLAQ2  = -EAQ2
      AMCH1=SQRT(XP*XXT*ECM*ECM-(PTXVQ1+PTXVA2)**2
     *       -(PTYVQ1+PTYVA2)**2)
C     AMCH1=SQRT((EQ1+EAQ2)**2-(PTXVQ1+PTXVA2)**2
C    *       -(PTYVQ1+PTYVA2)**2-(PLQ1+PLAQ2)**2)
      AMCH2=SQRT(XT*XXP*ECM*ECM-(PTXVQ2+PTXVA1)**2
     *       -(PTYVQ2+PTYVA1)**2)
C     AMCH2=SQRT((EQ2+EAQ1)**2-(PTXVQ2+PTXVA1)**2
C    *       -(PTYVQ2+PTYVA1)**2-(PLQ2+PLAQ1)**2)
      RETURN
      END
      SUBROUTINE KKEVT(NHKKH1,EPN,PPN,KKMAT,IREJ)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      COMMON/INTNEZ/NDZ,NZD
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
C      /
C INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C------------------
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
*KEEP,SHMAKL.
C     INCLUDE (SHMAKL)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON/SHMAKL/JSSH(INTMX),JTSH(INTMX),INTER1(INTMX),INTER2(INTMX)
*KEEP,NUCC.
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
*KEEP,RPTSHM.
      COMMON /RPTSHM/ RPROJ,RTARG,BIMPAC
*KEEP,NSHMAK.
      COMMON /NSHMAK/ NNSHMA,NPSHMA,NTSHMA,NSHMAC,NSHMA2
*KEEP,ZENTRA.
      COMMON /ZENTRA/ ICENTR
*KEEP,NUCIMP.
      COMMON /NUCIMP/ PRMOM(5,248),TAMOM(5,248), PRMFEP,PRMFEN,TAMFEP,
     +TAMFEN, PREFEP,PREFEN,TAEFEP,TAEFEN, PREPOT(210),TAEPOT(210),
     +PREBIN,TAEBIN,FERMOD,ETACOU
*KEEP,DROPPT.
      LOGICAL INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA, IPADIS,
     +ISHMAL,LPAULI
      COMMON /DROPPT/ INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA,
     +IPADIS,ISHMAL,LPAULI
*KEEP,NNCMS.
      COMMON /NNCMS/  GAMCM,BGCM,UMO,PCM,EPROJ,PPROJ
*KEEP,NUCPOS.
      COMMON /NUCPOS/INVVP(248),INVVT(248),INVSP(248),INVST(248), NUVV,
     +NUVS,NUSV,NUSS,INSVP(248),INSVT(248),INSSP(248),INSST(248), ISVEAP
     +(248),ISVEAT(248),ISSEAP(248),ISSEAT(248), IVSEAP(248),IVSEAT
     +(248), ISLOSP(248),ISLOST(248),INOOP(248),INOOT(248),NUOO
*KEEP,TAUFO.
      COMMON /TAUFO/  TAUFOR,KTAUGE,ITAUVE,INCMOD
      COMMON /EVAPPP/IEVAP
*KEEP,RTAR.
      COMMON /RTAR/ RTARNU
*KEEP,INNU.
      COMMON /INNU/INUDEC
*KEEP,HADTHR.
      COMMON /HADTHR/ EHADTH,INTHAD
*KEEP,DINPDA.
      COMMON /DINPDA/ IMPS(6,6),IMVE(6,6),IB08(6,21),IB10(6,21), IA08
     +(6,21),IA10(6,21), A1,B1,B2,B3,LT,LE,BET,AS,B8,AME,DIQ,ISU
*KEEP,FERMI.
      COMMON /FERMI/ PQUAR(4,248),PAQUAR(4,248), TQUAR(4,248),TAQUAR
     +(4,248)
*KEEP,KETMAS.
      COMMON /KETMAS/ AM8(2),AM10(2),IB88(2),IB1010(2),AMCH1N,AMCH2N
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
*KEEP,REJEC.
      COMMON /REJEC/ IRCO1,IRCO2,IRCO3,IRCO4,IRCO5, IRSS11,IRSS12,
     +IRSS13,IRSS14, IRSV11,IRSV12,IRSV13,IRSV14, IRVS11,IRVS12,IRVS13,
     +IRVS14, IRVV11,IRVV12,IRVV13,IRVV14
*KEEP,PROJK.
      COMMON /PROJK/ IPROJK
*KEEP,TANUIN.
      COMMON /TANUIN/ TASUMA,TASUBI,TABI,TAMASU,TAMA,TAIMMA
*KEND.
      COMMON /DIFFRA/ ISINGD,IDIFTP,IOUDIF,IFLAGD
*
      LOGICAL LSEADI
      COMMON /SEADIQ/ LSEADI
      COMMON /EVFLAG/NUMEV
      COMMON /DIQUAX/AMEDD,IDIQUA,IDIQUU
C
C-----------------------------------------------------------------------
C     PARAMETER (INTMX=2488)
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
*KEEP,ABRSS.
C     INCLUDE (ABRSS)
      COMMON /ABRSS/ AMCSS1(INTMX),AMCSS2(INTMX), GACSS1(INTMX),GACSS2
     +(INTMX), BGXSS1(INTMX),BGYSS1(INTMX),BGZSS1(INTMX), BGXSS2(INTMX),
     +BGYSS2(INTMX),BGZSS2(INTMX), NCHSS1(INTMX),NCHSS2(INTMX), IJCSS1
     +(INTMX),IJCSS2(INTMX), PQSSA1(INTMX,4),PQSSA2(INTMX,4), PQSSB1
     +(INTMX,4),PQSSB2(INTMX,4)
      COMMON /NUCJTN/NONUJ1,NONUJT,NONUS1,NONUST
      COMMON /XSVTHR/ XSTHR,XVTHR,XDTHR,XSSTHR
      COMMON /MINIJ/IMINIJ,NOMJE,NOMJER,NREJEV,NOMJT,NOMJTR
      COMMON /NCSHXX/NCOUXH,NCOUXT
      COMMON/INTNEU/NDZSU,NZDSU
      COMMON /VXSVD/VXSP(50),VXST(50),VXSAP(50),VXSAT(50),
     *              VXVP(50),VXVT(50),VXDP(50),VXDT(50),
     *      NXSP,NXST,NXSAP,NXSAT,NXVP,NXVT,NXDP,NXDT
      COMMON /NPARTT/NPART
      COMMON /ZSEA/ZSEAAV,ZSEASU,ANZSEA
      DATA NCOUSH/0/
      DATA NCOUST/0/
      DATA SUNDZ /0/
      DATA SUNZD /0/
      ANZSEA=0.D0
      ZSEASU=0.D0
      ZSEAAV=0.D0
C
      DO 5533 JJ=1,INTMX
        NCHSS1(JJ)=0
        NCHSS2(JJ)=0
 5533 CONTINUE
C                                    smoth transition between HADRIN and DPM
       EHADTW=EHADTH-RNDM(V)*2.D0
C*******************************************************************"
C
C     KINEMATICS
C
C********************************************************************
C
      IREJ = 0
*
      AAM(26)=AAM(23)
C
      KPROJ=1
      IF(IJPROJ.NE.0) KPROJ=IJPROJ
      KTARG=1
      ATNUC=IT
      ITN=IT-ITZ
      APNUC=IP
      IPN=IP-IPZ
      AMPROJ =AAM(KPROJ)
      AMTAR =AAM(KTARG)
*                             nucleon-nucleon cms
C     IBPROJ=1
      EPROJ=EPN
      PPROJ = SQRT((EPN-AMPROJ)*(EPN+AMPROJ))
      UMO = SQRT(AMPROJ**2 + AMTAR**2 + 2.*AMTAR*EPROJ)
      GAMCM = (EPROJ+AMTAR)/UMO
      BGCM=PPROJ/UMO
      ECM=UMO
      PCM=GAMCM*PPROJ - BGCM*EPROJ
C
      IF(IPEV.GE.1) PRINT 1000,IP,IPZ,IT,ITZ,IJPROJ,IBPROJ, EPROJ,PPROJ,
     +AMPROJ,AMTAR,UMO,GAMCM,BGCM
 1000 FORMAT(' ENTRY KKEVT'/ '    IP,IPZ,IT,ITZ,IJPROJ,IBPROJ',6I5/
     +'    EPROJ,PPROJ,AMPROJ,AMTAR,UMO,GAMCM,BGCM'/10E12.3)
 
C
C****                            CHANGE PARAMETERS FROM COMMON \INPDAT\
      AS=0.5
      B8=0.4
C                                CHAIN PT BIGGER THAN PARTICLE PT
      N9483=0
*  entry after rejection of an event because of kinematical reasons
*  several trials are made to realize a sampled Glauber event
   10 CONTINUE
      NDZ=0
      NZD=0
      N9483=N9483+1
      IF (MOD(N9483,125000).EQ.0) THEN
        WRITE(6,'(A,I5,A,I5,A)') ' KKEVT: Glauber event',NUMEV,
     +  ' rejected after', N9483, ' trials'
        WRITE(6, 1010) NN,NP,NT
        WRITE(6,1020) IRCO1,IRCO2,IRCO3,IRCO4,IRCO5, IRSS11,IRSS12,
     +  IRSS13,IRSS14, IRSV11,IRSV12,IRSV13,IRSV14, IRVS11,IRVS12,
     +  IRVS13,IRVS14, IRVV11,IRVV12,IRVV13,IRVV14
        N9483=1
        GO TO 20
      ELSEIF(N9483.GT.1) THEN
                                                                 GOTO 30
      ENDIF
 1010 FORMAT (5X,' N9483 LOOP - NN, NP, NT',5I10)
 1020 FORMAT (5X,' N9483 LOOP - REJECTION COUNTERS ',/,5I8/2(8I8/))
C
C***************************************************************
C
C     SAMPLE NUMBERS OF COLLISION A LA SHMAKOV---------------
C
C     TOTAL NUMBER OF INTERACTIONS = NN
C     NUMBER OF INTERACTING NUCLEONS
C                  FROM PROJECTILE = NP
C                  FROM TARGET = NT
C
      NCOUSH=NCOUSH+1
      NCOUXH=NCOUSH
      GO TO 2077
   20 CONTINUE
   22 CONTINUE
      NCOUST=NCOUST+1
      NCOUXT=NCOUST
 2077 CONTINUE
      CALL SHMAKO(IP,IT,BIMP,NN,NP,NT,JSSH,JTSH,PPROJ,KKMAT)
C     WRITE(6,*)' IP,IT,BIMP,NN,NP,NT ',IP,IT,BIMP,NN,NP,NT
      NPART=NP+NT
************    score characteristics of all sampled Glauber events
      CALL SHMAK(2,NN,NP,NT,IP,IT,ECM,BIMP)
      NSHMAC=NSHMAC+1
      IF ((ISINGD.GE.2).AND.((NT.NE.1).OR.(NN.NE.1))) GOTO 22
*
      IF (NN.GT.INTMX) THEN
        WRITE (6,1030)NN,NP,NT
 1030 FORMAT (' NN.GT.INTMX  SHMAKO SET TO INTMX ',3I10)
        NN=INTMX
      ENDIF
      NNSHMA=NN
      NPSHMA=NP
      NTSHMA=NT
C                    CENTRAL PRODUCTION FOR ICENTR.EQ.1 or 2
      IF (IP.LT.IT.AND.IT.LE.150)THEN
        IF(IP.LE.8)THEN
          IF ((ICENTR.EQ.1.OR.ICENTR.EQ.2).AND.NP.LT.IP-1)GO TO 20
        ELSEIF(IP.LE.16)THEN
          IF ((ICENTR.EQ.1.OR.ICENTR.EQ.2).AND.NP.LT.IP-2)GO TO 20
        ELSEIF(IP.LT.32)THEN
          IF ((ICENTR.EQ.1.OR.ICENTR.EQ.2).AND.NP.LT.IP-3)GO TO 20
        ELSEIF(IP.GE.32)THEN
C      Example S-Ag	
          IF ((ICENTR.EQ.1.OR.ICENTR.EQ.2).AND.NP.LT.IP-1)GO TO 20
        ENDIF
      ELSEIF (IP.LT.IT.AND.IT.GT.150)THEN
        IF(IP.LE.8)THEN
          IF ((ICENTR.EQ.1.OR.ICENTR.EQ.2).AND.NP.LT.IP-1)GO TO 20
        ELSEIF(IP.LE.16)THEN
          IF ((ICENTR.EQ.1.OR.ICENTR.EQ.2).AND.NP.LT.IP-2)GO TO 20
        ELSEIF(IP.LT.32)THEN
          IF ((ICENTR.EQ.1.OR.ICENTR.EQ.2).AND.NP.LT.IP-3)GO TO 20
        ELSEIF(IP.GE.32)THEN
C      Example S-Au	
          IF ((ICENTR.EQ.1.OR.ICENTR.EQ.2).AND.NP.LT.IP)GO TO 20
        ENDIF
      ELSEIF(IP.EQ.IT)THEN
        IF ((ICENTR.EQ.1.OR.ICENTR.EQ.2).AND.IP.EQ.32)THEN
C      Example S-S	
           IF(NP.LT.22.OR.NT.LT.22)                             GO TO 20
        ELSEIF ((ICENTR.EQ.1.OR.ICENTR.EQ.2).
     *AND.(UMO.GT.100.).AND.(NP.LT.IP-IP/10))THEN
C      Example Pb-Pb   central like at RHIC,LHC  UMO .GT.100	
                     GO TO 20
        ELSEIF ((ICENTR.EQ.1.OR.ICENTR.EQ.2).
     *AND.(UMO.LT.100.).AND.(NP.LT.IP-IP/4))THEN
C      Example Pb-Pb   central like at SPS umo .LT.100	
                     GO TO 20
        ELSEIF ((ICENTR.EQ.3).AND.NP.LT.IP-2*IP/3)THEN
C      Example Pb-Pb less central	
                     GO TO 20
        ENDIF
      ELSEIF(ABS(IP-IT).LT.3)THEN
        IF ((ICENTR.EQ.1.OR.ICENTR.EQ.2).AND.NP.LT.IP-IP/8)GO TO 20
      ENDIF
      BIMPAC=BIMP
C                    PERIPHERAL PRODUCTION FOR ICENTR.EQ.10
          IF (ICENTR.EQ.10.AND.NP.GT.6)                        GO TO 20
C------------------------------------------------------------
C     Drop diffractive collisions out of the Glauber
C     cascade in nucleus-nucleus collisions (For NN > 1 only)
      IF((ISINGD.LE.1).AND.(NN.GE.2).AND.(IP.GE.2).AND.(IT.GE.2).AND.
     *(IP.LE.200))THEN
        CALL DROPDI(NN,NP,NT,ECM)
        IF (IPEV.GE.3) THEN
        WRITE(6, 1040) IP,IPZ,IT,ITZ,EPROJ,PPROJ,NN,NP,NT
        WRITE (6,'(/A,2I5,1PE10.2,3I5)') ' KKEVT: IP,IT,BIMP,NN,NP,NT ',
     +  IP,IT,BIMP,NN,NP,NT
        WRITE (6,'(/2A)')
     +  ' KKEVT: JSSH(KKK),JTSH(KKK),INTER1(KKK),INTER2(KKK),',
     +  ' PKOO(3,KKK),TKOO(3,KKK)'
        ITUM=MAX(IT,IP,NN)
        DO 4011 KKK=1,ITUM
          WRITE (6,'(4I5,6(1PE11.3))') JSSH(KKK),JTSH(KKK),INTER1(KKK),
     +    INTER2(KKK), PKOO(1,KKK),PKOO(2,KKK),PKOO(3,KKK), TKOO(1,KKK),
     +    TKOO(2,KKK),TKOO(3,KKK)
 4011   CONTINUE
        ENDIF
      ENDIF
************    score characteristics of all sampled Glauber events
      CALL SHMAK1(2,NN,NP,NT,IP,IT,ECM,BIMP)
      NSHMA2=NSHMA2+1
C------------------------------------------------------------
C
*  entry for repeated trial to realize a sampled Glauber event
   30 CONTINUE
      IF (IPEV.GE.3) THEN
        WRITE(6, 1040) IP,IPZ,IT,ITZ,EPROJ,PPROJ,NN,NP,NT
 1040 FORMAT ('   752 FORM ',4I10,2F10.3,5I10)
        WRITE (6,'(/A,2I5,1PE10.2,3I5)') ' KKEVT: IP,IT,BIMP,NN,NP,NT ',
     +  IP,IT,BIMP,NN,NP,NT
        WRITE (6,'(/2A)')
     +  ' KKEVT: JSSH(KKK),JTSH(KKK),INTER1(KKK),INTER2(KKK),',
     +  ' PKOO(3,KKK),TKOO(3,KKK)'
        ITUM=MAX(IT,IP,NN)
        DO 40 KKK=1,ITUM
          WRITE (6,'(4I5,6(1PE11.3))') JSSH(KKK),JTSH(KKK),INTER1(KKK),
     +    INTER2(KKK), PKOO(1,KKK),PKOO(2,KKK),PKOO(3,KKK), TKOO(1,KKK),
     +    TKOO(2,KKK),TKOO(3,KKK)
 
   40   CONTINUE
      ENDIF
C
C-----------------------------------------------------------------------
C                  STORE PROJECTILE HADRON/NUCLEONS INTO /HKKEVT/
C                            - PROJECTILE CENTRE AT ORIGIN IN SPACE
C                            - TARGET SHIFTED IN X DIRECTION BY
C                                             IMPACT PARAMETER 'BIMP'
C
C                            - SAMPLING OF NUCLEON TYPES
C                            - CONSISTENCY CHECK
C                              FOR SAMPLED P/N NUMBERS
C                            - INTERACTING PROJECTILES ISTHKK=11
C                              NONINTERACTING ...      ISTHKK=13
C                            - FERMI MOMENTA IN CORRESP. REST SYSTEM
C-----------
      NHKK=0
C
      NCPP=0
      NCPN=0
C                      DEFINE FERMI MOMENTA/ENERGIES FOR PROJECTILE
C
      PXFE=0.0
      PYFE=0.0
      PZFE=0.0
      DO 50 KKK=1,IP
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IF (JSSH(KKK).GT.0) THEN
          ISTHKK(NHKK)=11
        ELSE
          ISTHKK(NHKK)=13
        ENDIF
C*
C                                       CHANGED 28.2.90.J.R.
C       IF(IJPROJ.EQ.0) THEN
        IF(IP.GE.2) THEN
C
          FRPNEU=FLOAT(IPN)/APNUC
          SAMTES=RNDM(v)
          IF(SAMTES.LT.FRPNEU.AND.NCPN.LT.IPN) THEN
            KPROJ=8
            NCPN=NCPN + 1
          ELSEIF(SAMTES.GE.FRPNEU.AND.NCPP.LT.IPZ) THEN
            KPROJ=1
            NCPP=NCPP + 1
          ELSEIF(NCPN.LT.IPN) THEN
            KPROJ=8
            NCPN=NCPN + 1
          ELSEIF(NCPP.LT.IPZ) THEN
            KPROJ=1
            NCPP=NCPP + 1
          ENDIF
C
          IF(KPROJ.EQ.1) THEN
            PFERM = PRMFEP
          ELSE
            PFERM = PRMFEN
          ENDIF
C         CALL FER4M(PFERM,FPX,FPY,FPZ,FE,KPROJ)
          CALL FER4MP(IP,PFERM,FPX,FPY,FPZ,FE,KPROJ)
          PXFE=PXFE + FPX
          PYFE=PYFE + FPY
          PZFE=PZFE + FPZ
          PHKK(1,NHKK)=FPX
          PHKK(2,NHKK)=FPY
          PHKK(3,NHKK)=FPZ
          PHKK(4,NHKK)=FE
          PHKK(5,NHKK)=AAM(KPROJ)
C
        ELSE
          KPROJ=IJPROJ
          PHKK(1,NHKK)=0.
          PHKK(2,NHKK)=0.
          PHKK(3,NHKK)=0.
          PHKK(4,NHKK)=AAM(KPROJ)
          PHKK(5,NHKK)=AAM(KPROJ)
        ENDIF
C
        KKPROJ(KKK)=KPROJ
        IDHKK(NHKK)=MPDGHA(KPROJ)
        JMOHKK(1,NHKK)=0
        JMOHKK(2,NHKK)=0
        JDAHKK(1,NHKK)=0
        JDAHKK(2,NHKK)=0
C
        PHKK(5,NHKK)=AAM(KPROJ)
        VHKK(1,NHKK)=PKOO(1,KKK)*1.E-12
        VHKK(2,NHKK)=PKOO(2,KKK)*1.E-12
        VHKK(3,NHKK)=PKOO(3,KKK)*1.E-12
        VHKK(4,NHKK)=0.
        WHKK(1,NHKK)=PKOO(1,KKK)*1.E-12
        WHKK(2,NHKK)=PKOO(2,KKK)*1.E-12
        WHKK(3,NHKK)=PKOO(3,KKK)*1.E-12
        WHKK(4,NHKK)=0.
        JHKKNP(KKK)=NHKK
C
        IF (IPHKK.GE.2) WRITE(6,1050) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     +  JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +  (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
 1050 FORMAT (I6,I4,5I6,9E10.2)
C
   50 CONTINUE
C                          balance Sampled Fermi momenta
      IF(IP.GE.2) THEN
        PXFE=PXFE/IP
        PYFE=PYFE/IP
        PZFE=PZFE/IP
        DO 60 KKK=1,IP
          IHKK=KKK
          PHKK(1,IHKK)=PHKK(1,IHKK) - PXFE
          PHKK(2,IHKK)=PHKK(2,IHKK) - PYFE
          PHKK(3,IHKK)=PHKK(3,IHKK) - PZFE
          PHKK(4,IHKK)=SQRT(PHKK(5,IHKK)** 2+ PHKK(1,IHKK)** 2+ PHKK
     +    (2,IHKK)** 2+ PHKK(3,IHKK)**2)
        IF (IPHKK.GE.2) WRITE(6,1050) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +  JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +  (PHKK(KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
   60   CONTINUE
      ENDIF
C
C-----------------------------------------------------------------------
C                  STORE TARGET HADRON/NUCLEONS INTO /HKKEVT/
C                            - PROJECTILE CENTRE AT ORIGIN IN SPACE
C                            - TARGET SHIFTED IN X DIRECTION BY
C                                             IMPACT PARAMETER 'BIMP'
C
C                            - SAMPLING OF NUCLEON TYPES
C                            - CONSISTENCY CHECK
C                              FOR SAMPLED P/N NUMBERS
C                            - INTERACTING TARGETS     ISTHKK=12
C                              NONINTERACTING ...      ISTHKK=14
C-----------
C---------------------
      NHADRI=0
      NCTP=0
      NCTN=0
C
      TXFE=0.0
      TYFE=0.0
      TZFE=0.0
      DO 70 KKK=1,IT
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IF (JTSH(KKK).GT.0) THEN
          ISTHKK(NHKK)=12
          NHADRI=NHADRI+1
          IF (NHADRI.EQ.1) IHTAWW=NHKK
          IF (EPN.LE.EHADTW) THEN
            IF (NHADRI.GT.1) ISTHKK(NHKK)=14
          ENDIF
        ELSE
          ISTHKK(NHKK)=14
        ENDIF
        IF(IT.GE.2)THEN
        FRTNEU=FLOAT(ITN)/ATNUC
        SAMTES=RNDM(V)
        IF(SAMTES.LT.FRTNEU.AND.NCTN.LT.ITN) THEN
          KTARG=8
          NCTN=NCTN + 1
        ELSEIF(SAMTES.GE.FRTNEU.AND.NCTP.LT.ITZ) THEN
          KTARG=1
          NCTP=NCTP + 1
        ELSEIF(NCTN.LT.ITN) THEN
          KTARG=8
          NCTN=NCTN + 1
        ELSEIF(NCTP.LT.ITZ) THEN
          KTARG=1
          NCTP=NCTP + 1
        ENDIF
C
        IF(KTARG.EQ.1) THEN
          PFERM = TAMFEP
        ELSE
          PFERM = TAMFEN
        ENDIF
C       CALL FER4M(PFERM,FPX,FPY,FPZ,FE,KTARG)
        CALL FER4MT(IT,PFERM,FPX,FPY,FPZ,FE,KTARG)
        TXFE=TXFE + FPX
        TYFE=TYFE + FPY
        TZFE=TZFE + FPZ
        PHKK(1,NHKK)=FPX
        PHKK(2,NHKK)=FPY
        PHKK(3,NHKK)=FPZ
        PHKK(4,NHKK)=FE
        PHKK(5,NHKK)=AAM(KTARG)
        ELSE
        PHKK(1,NHKK)=0. 
        PHKK(2,NHKK)=0. 
        PHKK(3,NHKK)=0. 
        PHKK(4,NHKK)=AAM(KTARG)
        PHKK(5,NHKK)=AAM(KTARG)
        ENDIF
C
        KKTARG(KKK)=KTARG
        IDHKK(NHKK)=MPDGHA(KTARG)
        JMOHKK(1,NHKK)=0
        JMOHKK(2,NHKK)=0
        JDAHKK(1,NHKK)=0
        JDAHKK(2,NHKK)=0
        VHKK(1,NHKK)=(TKOO(1,KKK)+BIMP)*1.E-12
        VHKK(2,NHKK)=TKOO(2,KKK)*1.E-12
        VHKK(3,NHKK)=TKOO(3,KKK)*1.E-12
        VHKK(4,NHKK)=0.
        WHKK(1,NHKK)=(TKOO(1,KKK)+BIMP)*1.E-12
        WHKK(2,NHKK)=TKOO(2,KKK)*1.E-12
        WHKK(3,NHKK)=TKOO(3,KKK)*1.E-12
        WHKK(4,NHKK)=0.
        JHKKNT(KKK)=NHKK
C
        IF (IPHKK.GE.2) WRITE(6,1050) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     +  JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +  (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
C
   70 CONTINUE
C                          balance Sampled Fermi momenta
      IF(IT.GE.2) THEN
        TASUMA=ITZ*AAM(1) + (IT-ITZ)*AAM(8)
        TASUBI=0.0
        TAMASU=0.0
        TXFE=TXFE/IT
        TYFE=TYFE/IT
        TZFE=TZFE/IT
        DO 80 KKK=1,IT
          IHKK=KKK + IP
          PHKK(1,IHKK)=PHKK(1,IHKK) - TXFE
          PHKK(2,IHKK)=PHKK(2,IHKK) - TYFE
          PHKK(3,IHKK)=PHKK(3,IHKK) - TZFE
          PHKK(4,IHKK)=SQRT(PHKK(5,IHKK)** 2+ PHKK(1,IHKK)** 2+ PHKK
     +    (2,IHKK)** 2+ PHKK(3,IHKK)**2)
          ITSEC=MCIHAD(IDHKK(IHKK))
          TASUBI=TASUBI + PHKK(4,IHKK) - PHKK(5,IHKK) - TAEPOT(ITSEC)
          TAMASU=TAMASU + PHKK(4,IHKK) - TAEPOT(ITSEC)
        IF (IPHKK.GE.2) WRITE(6,1050) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +  JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +  (PHKK(KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
   80   CONTINUE
C***         definition of initial state
        TABI=-EBIND(IT,ITZ)
        TAMA=(IT-ITZ)*AAM(8) + ITZ*AAM(1) + TABI
        TAIMMA=TAMA - TAMASU
      ENDIF
C
      IF(IPEV.GT.3) THEN
        WRITE(6,'(/A/5X,A/5X,4(1PE11.3))') ' KKEVT: FERMI MOMENTA',
     +  'PRMFEP,PRMFEN, TAMFEP,TAMFEN', PRMFEP,PRMFEN, TAMFEP,TAMFEN
 
      ENDIF
C-----------------------------------------------------------------------
      IFLAGD = 0
C-----------------------------------------------------------------
C
C                                   Test j.r. 2/94
C
C----------------------------------------------------------------
	IQQDD=0
      IF(IPEV.GT.3) THEN
	WRITE(6,'(A,4I5)')' KKEVT before SDIFF',NP,NT,NN,ISINGD
      ENDIF
      IF ((NP.EQ.1).AND.(NT.EQ.1).AND.(NN.EQ.1)
C         Diffraction als for A-B collisions j.r. 2.2.99
C    &.AND.((IP.EQ.1).OR.(IT.EQ.1))
     &.AND.(EPN.GT.EHADTW))
     &   CALL SDIFF(EPROJ,PPROJ,KPROJ,NHKKH1,IQQDD)
C----------------------------------------------------------------
      IF (IFLAGD.EQ.1) RETURN
*
C----------------------------------------------------------------------
C
C********************************    SAMPLE THE X FRACTIONS
C                                    OF INTERACTING NUCLEONS / HADRONS
C
      IF (EPN.LE.EHADTW) THEN
 7107   CONTINUE
        ITTA=MCIHAD(IDHKK(IHTAWW))
      IF(IPEV.GT.1) THEN
        WRITE(6,'(A,I5,2F10.3)')' HADRIN CALL, IREJFO=',IREJFO, EHADTW
     *	,EHADTH
      ENDIF   
        CALL HADHAD(EPN,PPN,NHKKH1,IHTAWW,ITTA,IREJFO)
        IF(IREJFO.EQ.1) GO TO 7107
C                                  Transform into cms
	DO 111 I=NHKKH1+1,NHKK
	  PZNN=PHKK(3,I)
	  ENN=PHKK(4,I)
	  PHKK(3,I)=GAMCM*PZNN-BGCM*ENN
	  PHKK(4,I)=GAMCM*ENN-BGCM*PZNN
  111   CONTINUE
	
                                                           GO TO 110
      ENDIF
C-----------------------------------------------------------------------
C
C*********************  SAMPLE THE FLAVORS OF INTERACTING
C                       PROJECTILE AND TARGET HADRONS/NUCLEONS
C                       first run sea quarks
      CALL FLKSAA(NN,ECM)
C
      CALL XKSAMP(NN,ECM)
      DO 90 N=1,NSS
        INLOSS(N)=.TRUE.
   90 CONTINUE
      IF(IPEV.GE.6)WRITE(6,*)' KKEVT ,NSS,NVS,NSV,NVV,NDS,NSD,NDV,NVD ',
     *' after XKSAMP ',
     *NSS,NVS,NSV,NVV,NDS,NSD,NDV,NVD
C
      IF (IPEV.GE.6) THEN
        ITUM=MAX0(IP,IT,NN)
        WRITE(6,'(A,I10)')' KKEVT ITUM loop limit',ITUM
        WRITE(6,'(A,2A)') ' KKEVT (AFTER XKSAMP):',
     +  ' JSSH, JSSHS, JTSH, JTSHS, INTER1, INTER2',
     +  ' PKOO(1),PKOO(2),PKOO(3), TKOO(1),TKOO(2),TKOO(3)'
        DO 100 KKK=1,ITUM
          WRITE (6,'(6I5,6(1PE11.3))') JSSH(KKK),JSSHS(KKK),JTSH(KKK),
     +    JTSHS(KKK), INTER1(KKK),INTER2(KKK), PKOO(1,KKK),PKOO(2,KKK),
     +    PKOO(3,KKK), TKOO(1,KKK),TKOO(2,KKK),TKOO(3,KKK)
 
 
  100   CONTINUE
      ENDIF
C-----------------------------------------------------------------------
C
C*********************  SAMPLE THE FLAVORS OF INTERACTING
C                       PROJECTILE AND TARGET HADRONS/NUCLEONS
C                       second run valence quarks
      IF(IPEV.GE.6)WRITE(6,*)' KKEVT ,NSS,NVS,NSV,NVV,NDS,NSD,NDV,NVD ',
     *NSS,NVS,NSV,NVV,NDS,NSD,NDV,NVD
      IF(IPEV.GE.6)WRITE(6,'(A)')' KKEVT before flksam'
      CALL FLKSAM
C     IPEV=8
      IF(IPEV.GE.6)WRITE(6,'(A)')' KKEVT after flksam'
      
 1012 FORMAT(' XKSAMP:',
     +' I,XPVQ(I),XPVD(I),IFROVP(I),ITOVP(I),ZUOVP(I),KKPROJ(I)')
 1022 FORMAT(I5,2E15.5,2I5,L5,I5)
 1032 FORMAT(' XKSAMP :  I,XPSQ(I),XPSAQ(I),IFROSP(I),ZUOSP(I)')
 1042 FORMAT(I5,2E15.5,I5,L5)
 1060 FORMAT(' XKSAMP :  I,XTSQ(I),XTSAQ(I),IFROST(I),ZUOST(I)')
 1052 FORMAT(' XKSAMP:',
     +' I,XTVQ(I),XTVD(I),IFROVT(I),ITOVT(I),ZUOVT(I),KKTARG(I)')
C     COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
C    +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
C    +JHKKNT
C    +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
C    +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
C    &                MHKKHH(INTMX),
C    +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
C
      DO 511 I=1,IXPV
        IIPV=1+XPVQ(I)/0.02D0
	VXVP(IIPV)=VXVP(IIPV)+1.D0
        IIPD=1+XPVD(I)/0.02D0
	VXDP(IIPD)=VXDP(IIPD)+1.D0
	NXVP=NXVP+1
	NXDP=NXDP+1
  511 CONTINUE      
      DO 521 I=1,IXPS
        IIPS=1+XPSQ(I)/0.02D0
	VXSP(IIPS)=VXSP(IIPS)+1.D0
        IIPA=1+XPSAQ(I)/0.02D0
	VXSAP(IIPA)=VXSAP(IIPA)+1.D0
	NXSP=NXSP+1
	NXSAP=NXSAP+1
  521 CONTINUE      
      DO 531 I=1,IXTV
        IITV=1+XTVQ(I)/0.02D0
	VXVT(IITV)=VXVT(IITV)+1.D0
        IITD=1+XTVD(I)/0.02D0
	VXDT(IITD)=VXDT(IITD)+1.D0
	NXVT=NXVT+1
	NXDT=NXDT+1
  531 CONTINUE      
      DO 541 I=1,IXTS
        IITS=1+XTSQ(I)/0.02D0
	VXST(IITS)=VXST(IITS)+1.D0
        IITA=1+XTSAQ(I)/0.02D0
	VXSAT(IITA)=VXSAT(IITA)+1.D0
	NXST=NXST+1
	NXSAT=NXSAT+1
  541 CONTINUE      
      IF(IPCO.GE.1)THEN
        WRITE(6,'(A)')
     +  ' XKSAMP :  FINAL X-VALUES AFTER POTENTIAL CORRECTION'
        WRITE(6,1012)
        DO 510 I=1,IXPV
          WRITE(6,1022) I,XPVQ(I),XPVD(I),IFROVP(I),ITOVP(I),ZUOVP(I)
          WRITE(6,*)' I(1-IXPV),IPVQ(I),IPPV1(I),IPPV2(I)JHKKPV(I) ',
     *    I,IPVQ(I),IPPV1(I),IPPV2(I),JHKKPV(I)
  510   CONTINUE
        WRITE(6,1032)
        DO 520 I=1,IXPS
          WRITE(6,1042) I,XPSQ(I),XPSAQ(I),IFROSP(I),ZUOSP(I)
          WRITE(6,*)' I(1-IXPS),IPSQ(I),IPSAQ(I ),JHKKPS(I) ',
     *    I,IPSQ(I),IPSAQ(I),JHKKPS(I)
  520   CONTINUE
        WRITE(6,1052)
        DO 530 I=1,IXTV
          WRITE(6,1022) I,XTVQ(I),XTVD(I),IFROVT(I),ITOVT(I),ZUOVT(I)
          WRITE(6,*)' I(1-IXTV),ITVQ(I),ITTV1(I),ITTV2(I),JHKKTV(I) ',
     *    I,ITVQ(I),ITTV1(I),ITTV2(I),JHKKTV(I)
  530   CONTINUE
        WRITE(6,1060)
        DO 540 I=1,IXTS
          WRITE(6,1042) I,XTSQ(I),XTSAQ(I),IFROST(I),ZUOST(I)
          WRITE(6,*)' I(1-IXTS),ITSQ(I),ITSAQ(I),JHKKTS(I) ',
     *    I,ITSQ(I),ITSAQ(I),JHKKTS(I)
  540   CONTINUE
      ENDIF
      IF(IPEV.GE.6)WRITE(6,'(A,6I5)')
     *' XKSAMP NSV,NDV,NVS,NVD',
     +    NSV,NDV,NVS,NVD
C     IPEV=-1
C
C-------------------------
C                 TRANSFORM MOMENTA OF INTERACTING NUCLEONS
C                 (INCLUDING FERMI MOMENTA FROM NUCLEUS REST FRAMES)
C                 INTO NUCLEON-NUCLEON CMS (DEFINED WITHOUT FERMI MOM.
      IF(IPEV.GE.6)WRITE(6,'(A)')' KKEVT before NUCMOM'
      DO 7745 IHKK=1,NHKK
        IF (IPHKK.GE.2) WRITE(6,1050) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +  JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +  (PHKK(KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
 7745 CONTINUE
      CALL NUCMOM
      IF(IPEV.GE.6)WRITE(6,'(A)')' KKEVT after NUCMOM'
      NONUST=0
      NONUJT=0
      NOMJE=0
      NOMJER=0
      IF(IPEV.GE.6)WRITE(6,*)' KKEVT ,NSS,NVS,NSV,NVV,NDS,NSD,NDV,NVD ',
     *NSS,NVS,NSV,NVV,NDS,NSD,NDV,NVD
C
C-----------------------------------------------------------------------
C
C          NOTE: THE FOLLOWING TREATMENT OF CHAIN SYSTEMS
C          ****  GENERALLY STARTS FROM THE NUCLEON-NUCLEON CMS
C                (DEFINED WITHOUT FERMI MOMENTA)
C-------------------------    TREATMENT OF SEA-SEA CHAIN SYSTEMS
      IF(IPEV.GE.6)WRITE(6,*)' KKEVT before KKEVSS, NSS',NSS
      IF(NSS.GT.0) CALL KKEVSS
      IF(IPEV.GE.6)WRITE(6,'(A)')' KKEVT after KKEVSS'
C
C-----------------  TREATMENT OF sea diquark - sea CHAIN SYSTEMS
      IF(IPEV.GE.6)WRITE(6,*)' KKEVT before KKEVDS, NDS',NDS
C     IF(NDS.GT.0 .AND. LSEADI) THEN
      IF(NDS.GT.0) THEN
      IF(IPEV.GE.6)WRITE(6,'(A)')' KKEVT before KKEVDS'
      IF(IDIQUA.EQ.1)  CALL KKEVDS(IREJDS)
      IF(IPEV.GE.6)WRITE(6,'(A)')' KKEVT after KKEVDS'
C       IPEV=1
        IF (IREJDS.EQ.1)                                        GO TO 10
      ENDIF
C
C
C-----------------  TREATMENT OF sea  - sea-diquark CHAIN SYSTEMS
      IF(IPEV.GE.6)WRITE(6,*)' KKEVT before KKEVSD NSD',NSD
C     IF(NSD.GT.0 .AND. LSEADI) THEN
      IF(NSD.GT.0) THEN
      IF(IPEV.GE.6)WRITE(6,'(A)')' KKEVT before KKEVSD'
      IF(IDIQUA.EQ.1)  CALL KKEVSD(IREJSD)
      IF(IPEV.GE.6)WRITE(6,'(A)')' KKEVT after KKEVSD'
C       IPEV=1
        IF (IREJSD.EQ.1)                                        GO TO 10
      ENDIF
C
C
C-------------------------    TREATMENT OF SEA-VALENCE CHAIN SYSTEMS
C     IPEV=6
      IF(IPEV.GE.6)WRITE(6,*)' KKEVT before KKEVSV, NSV',NSV
      IF(NSV.GT.0) THEN
      IF(IPEV.GE.6)WRITE(6,'(A)')' KKEVT before KKEVSV'
        CALL KKEVSV(IREJSV)
      IF(IPEV.GE.6)WRITE(6,'(A)')' KKEVT after KKEVSV'
C       IPEV=1
        IF (IREJSV.EQ.1)                                        GO TO 10
      ENDIF
C
C-----------------  TREATMENT OF sea diquark - VALENCE CHAIN SYSTEMS
      IF(IPEV.GE.6)WRITE(6,*)' KKEVT before KKEVDV, NDV',NDV
C     IF(NDV.GT.0 .AND. LSEADI) THEN
      IF(NDV.GT.0) THEN
      IF(IPEV.GE.6)WRITE(6,'(A)')' KKEVT before KKEVDV'
      IF(IDIQUA.EQ.1)  CALL KKEVDV(IREJDV)
      IF(IPEV.GE.6)WRITE(6,'(A)')' KKEVT after KKEVDV'
C       IPEV=1
        IF (IREJDV.EQ.1)                                        GO TO 10
      ENDIF
C
C-------------------------    TREATMENT OF VALENCE-SEA CHAIN SYSTEMS
C     IPEV=6
      IF(IPEV.GE.6)WRITE(6,*)' KKEVT before KKEVVS, NVS',NVS
      IF(NVS.GT.0) THEN
      IF(IPEV.GE.6)WRITE(6,'(A)')' KKEVT before KKEVVS'
        CALL KKEVVS(IREJVS)
      IF(IPEV.GE.6)WRITE(6,'(A)')' KKEVT after KKEVVS'
C       IPEV=1
        IF (IREJVS.EQ.1)                                        GO TO 10
      ENDIF
C
C-----------------  TREATMENT OF valence - sea diquark CHAIN SYSTEMS
      IF(IPEV.GE.6)WRITE(6,*)' KKEVT before KKEVVD,NVD',NVD
C     IF(NVD.GT.0 .AND. LSEADI) THEN
      IF(NVD.GT.0) THEN
      IF(IPEV.GE.6)WRITE(6,'(A)')' KKEVT before KKEVVD'
      IF(IDIQUA.EQ.1)  CALL KKEVVD(IREJVD)
      IF(IPEV.GE.6)WRITE(6,'(A)')' KKEVT after KKEVVD'
C       IPEV=1
        IF (IREJVD.EQ.1)                                        GO TO 10
      ENDIF
C
C-------------------    TREATMENT OF VALENCE-VALENCE CHAIN SYSTEMS
C
      IF(IPEV.GE.6)WRITE(6,*)' KKEVT before KKEVVV, NVV',NVV
      CALL KKEVVV(IREJVV,IBPROJ)
      IF(IPEV.GE.6)WRITE(6,'(A)')' KKEVT after KKEVVV'
C    &            IPVQ,IPPV1,IPPV2,ITVQ,ITTV1,ITTV2)
      IF (IREJVV.EQ.1)                                          GO TO 10
C     DO 5004 IHKK=1,NHKK
C         IF (IPHKK.GE.1) WRITE(6,5001)
C    *      IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),JMOHKK(2,IHKK),
C    &      JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK(KHKK,IHKK),KHKK=1,5),
C    &      (VHKK(KHKK,IHKK),KHKK=1,4)
C5004 CONTINUE
C

C
C-------------------    TREATMENT OF HARD CHAIN SYSTEMS
C
      IF(IPEV.GE.6)WRITE(6,*)' KKEVT before KKEVHH, NHH',NHH
      IF (IMINIJ.EQ.1) CALL KKEVHH
      IF(IPEV.GE.6)WRITE(6,'(A)')' KKEVT after KKEVHH'
      IF(IPEV.GE.6)WRITE(6,*)' KKEVT before KKEVZZ, NZZ',NZZ
      CALL KKEVZZ
      IF(IPEV.GE.6)WRITE(6,'(A)')' KKEVT after KKEVZZ'
      NOMJT=NOMJT+NOMJE
      NOMJTR=NOMJTR+NOMJER
      IEVI=IEVI+1
C     IF(IEVI.LE.20)WRITE (6,7233)NOMJE,NOMJER,NREJEV,NOMJT,NOMJTR
C7233 FORMAT (  '   NOMJE,NOMJER,NREJEV,NOMJT,NOMJTR ',5I10)
        DO 7787 III=1,NONUJT
          IF (IJJQ1(III).EQ.0.OR.IJJAQ1(III).EQ.0)THEN
            WRITE (6,7786)III,JHKKEX(III),IJJQ1(III),IJJAQ1(III),
     *                    AMJCH1(III)
 7786       FORMAT(' KKEVHH: III,JHKKEX,IJJQ1,IJJAQ1,AMCH1 ',4I10,F10.3)
            JHKKEX(III)=0
          ENDIF
 7787   CONTINUE
C-------------------------------------------------------------------
C
C                        COMBINE JETS
C
C------------------------------------------------------------------
C     IF(IPEV.GE.1)WRITE(6,*)' KKEVT before KKEVCC',LCOMBI
C     IF(LCOMBI) CALL KKEVCC(NN,IP,IT)
C     IF(IPEV.GE.6)WRITE(6,'(A)')' KKEVT after KKEVCC'
C----------------------------------------------------------------------
C          - TEST OF ENERGY MOMENTUM CONSERVATION ON NIVEAU OF CHAINS
C                                      AND ON NIVEAU OF CHAIN  ENDS
      IF(IPEV.GE.6)WRITE(6,'(A)')' KKEVT before EVTEST'
C     IF(IPEV.GE.1)CALL EVTEST(IREJ)
      CALL EVTEST(IREJ)
        IF (IPEV.GE.1) THEN
        WRITE(6,'(/A/)') ' KKEVT: FINAL LIST OF ENTRIES TO /HKKEVT/'
          DO 121 IHKK=1,NHKK
          WRITE(6,1050) IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),
     +    JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK
     +    (KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
 
  121     CONTINUE
        ENDIF
      IF (IREJ.EQ.1)THEN 
        IF(IPEV.GE.1)WRITE(6,'(A)')' EVTEST REJECTION would be '
 	IREJ=0
        IF (IREJ.EQ.1)GO TO 10
      ENDIF
      IF(IPEV.GE.1)WRITE(6,'(A)')' KKEVT after EVTEST'
C
C-----------------------------------------------------------------------
C                        CONSTRUCT HISTOGRAMS FROM PARTONS ON CHAIN ENDS
C
C     IF(IPADIS) CALL DISTPA(2)
C
C-----------------------------------------------------------------------
C               -  HADRONIZATION OF CHAINS  (HADRKK)
C                  AND BACK TRANSFORMATION FROM NN-CMS INTO LAB (HADRKK)
C
      IF(IPEV.GE.1)WRITE(6,'(A)')' KKEVT long before HADRKK'
      IF(IHADA.OR.IHADSS.OR.IHADSV.OR.IHADVS.OR.IHADVV) THEN
      IF(IPEV.GE.1)WRITE(6,'(A)')' KKEVT before HADRKK'
        CALL HADRKK(NHKKH1,PPN)
      IF(IPEV.GE.1)WRITE(6,'(A)')' KKEVT after HADRKK'
      ENDIF
C
  110 CONTINUE
C
C                           Correct HKKEVT COMMON
C                           ONLY FOR RUNS WITHOUT EVAPORATION
C     IF((IEVAP.EQ.0).AND.(KTAUGE.EQ.0))CALL CORRCO
C     not for runs with hadhad
      IF (EPN.GE.EHADTW) THEN
        CALL CORRCO
      ENDIF
C
C        Trigger ICENTR=8 NCH > 5 (NA35 p-S data)
      IIICH=0
      IF (ICENTR.EQ.8)THEN
        DO 128 IHKK=1,NHKK
	  IF(ISTHKK(IHKK).EQ.1)THEN
            NRHKK=MCIHAD(IDHKK(IHKK))
            ICHHKK=IICH(NRHKK)
            IF(ICHHKK.NE.0)THEN
C             PTT=PHKK(1,IHKK)**2+PHKK(2,IHKK)**2+0.000001
C             AMT=SQRT(PTT+PHKK(5,IHKK)**2)
C             YL=LOG((ABS(PHKK(3,IHKK) + PHKK(4,IHKK)))/AMT+1.E-18)
C             IF(YL.GT.0.6D0.AND.YL.LT.7.D0)THEN
	        IIICH=IIICH+1
C             ENDIF
	    ENDIF
	  ENDIF
  128   CONTINUE
         IF(IIICH.LE.14)THEN
          WRITE(6,*)' reject ',IIICH
          GO TO 22
	ENDIF
          WRITE(6,*)' no reject ',IIICH
      ENDIF
C
      IF (IPEV.GE.1) THEN
        WRITE(6,'(/A/)') ' KKEVT: FINAL LIST OF ENTRIES TO /HKKEVT/'
        DO 120 IHKK=1,NHKK
          WRITE(6,1050) IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),
     +    JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK
     +    (KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
 
  120   CONTINUE
      ENDIF
C
      SUNDZ=SUNDZ+NDZ
      SUNZD=SUNZD+NZD
      NZDSU=SUNZD
      NDZSU=SUNDZ
      IF(IPEV.GE.6)WRITE(6,*)' END KKEVT NZD,NZDSU,NDZ,NDZSU',
     * NZD,NZDSU,NDZ,NDZSU
      RETURN
      END
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE KKEVVV(IREJVV,NBPROJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C    &                  IPVQ,IPPV1,IPPV2,ITVQ,ITTV1,ITTV2)
C
C-------------------    TREATMENT OF VALENCE-VALENCE CHAIN SYSTEMS
C
C
*KEEP,NNCMS.
      COMMON /NNCMS/  GAMCM,BGCM,UMO,PCM,EPROJ,PPROJ
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
*KEEP,TRAFOP.
      COMMON /TRAFOP/ GAMP,BGAMP,BETP
*KEEP,NUCC.
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
*KEEP,NUCIMP.
      COMMON /NUCIMP/ PRMOM(5,248),TAMOM(5,248), PRMFEP,PRMFEN,TAMFEP,
     +TAMFEN, PREFEP,PREFEN,TAEFEP,TAEFEN, PREPOT(210),TAEPOT(210),
     +PREBIN,TAEBIN,FERMOD,ETACOU
*KEEP,ABRVV.
      COMMON /ABRVV/ AMCVV1(248),AMCVV2(248),GACVV1(248),GACVV2(248),
     +BGXVV1(248),BGYVV1(248),BGZVV1(248), BGXVV2(248),BGYVV2(248),
     +BGZVV2(248), NCHVV1(248),NCHVV2(248),IJCVV1(248),IJCVV2(248),
     +PQVVA1(248,4),PQVVA2(248,4), PQVVB1(248,4),PQVVB2(248,4)
*KEEP,DROPPT.
      LOGICAL INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA, IPADIS,
     +ISHMAL,LPAULI
      COMMON /DROPPT/ INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA,
     +IPADIS,ISHMAL,LPAULI
*KEEP,NUCPOS.
      COMMON /NUCPOS/INVVP(248),INVVT(248),INVSP(248),INVST(248), NUVV,
     +NUVS,NUSV,NUSS,INSVP(248),INSVT(248),INSSP(248),INSST(248), ISVEAP
     +(248),ISVEAT(248),ISSEAP(248),ISSEAT(248), IVSEAP(248),IVSEAT
     +(248), ISLOSP(248),ISLOST(248),INOOP(248),INOOT(248),NUOO
*KEEP,TAUFO.
      COMMON /TAUFO/  TAUFOR,KTAUGE,ITAUVE,INCMOD
*KEEP,RTAR.
      COMMON /RTAR/ RTARNU
*KEEP,INNU.
      COMMON /INNU/INUDEC
*KEEP,DINPDA.
      COMMON /DINPDA/ IMPS(6,6),IMVE(6,6),IB08(6,21),IB10(6,21), IA08
     +(6,21),IA10(6,21), A1,B1,B2,B3,LT,LE,BET,AS,B8,AME,DIQ,ISU
*KEEP,FERMI.
      COMMON /FERMI/ PQUAR(4,248),PAQUAR(4,248), TQUAR(4,248),TAQUAR
     +(4,248)
*KEEP,KETMAS.
      COMMON /KETMAS/ AM8(2),AM10(2),IB88(2),IB1010(2),AMCH1N,AMCH2N
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
*KEEP,REJEC.
      COMMON /REJEC/ IRCO1,IRCO2,IRCO3,IRCO4,IRCO5, IRSS11,IRSS12,
     +IRSS13,IRSS14, IRSV11,IRSV12,IRSV13,IRSV14, IRVS11,IRVS12,IRVS13,
     +IRVS14, IRVV11,IRVV12,IRVV13,IRVV14
*KEEP,PROJK.
      COMMON /PROJK/ IPROJK
      COMMON/RPTSHM/RPROJ,RTARG,BIMPAC
*KEND.
C
C-----------------------------------------------------------------------
C     PARAMETER (INTMX=3988)
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
      COMMON /NUCJTN/NONUJ1,NONUJT,NONUS1,NONUST
      COMMON /XSVTHR/ XSTHR,XVTHR,XDTHR,XSSTHR
      COMMON /MINIJ/IMINIJ,NOMJE,NOMJER,NREJEV,NOMJT,NOMJTR
C-----------------------------------------------------------------
C-------------------------------Single Chain Option-----------
      COMMON /SINCHA/ISICHA
      COMMON /ZSEA/ZSEAAV,ZSEASU,ANZSEA
      DATA ISCH/0/
      DATA NZSEA/0/
      ISCH=ISCH+1
      ZERO=0 
C------------------------------Single Chain Option--------------------
      IREJVV=0
C
      IMINIJ=1
      DO 20 N=1,NVV
C
C---------------------------drop recombined chain pairs
        IF(NCHVV1(N).EQ.99.AND.NCHVV2(N).EQ.99)GO TO 20
C
C
C
        IXVPR=INTVV1(N)
        INUCPR=IFROVP(IXVPR)
C
        IXVTA=INTVV2(N)
        INUCTA=IFROVT(IXVTA)
C
        XMAX1=XPVQ(IXVPR)+XPVD(IXVPR)-XVTHR-XDTHR
        XMAX2=XTVQ(IXVTA)+XTVD(IXVTA)-XVTHR-XDTHR
C
C
        IF(IPEV.GE.1)WRITE(6,'(A,2I5,3F9.3)')' KKEVVV,bef xptfl:n,nvv'
     *  ,N,NVV,XMAX1,XMAX2
        IF (IMINIJ.EQ.1)THEN
      	  CALL XPTFL(NHARD,NSEA,IREG,XMAX1,XMAX2)
C	  NZSEA=NZSEA+1
C	  ANZSEA=NZSEA
	  ANZSEA=ANZSEA+1.D0
	  ZSEASU=ZSEASU+NSEA
	  ZSEAAV=ZSEASU/ANZSEA
        ENDIF
        IF(IPEV.GE.1)WRITE(6,'(A,3I10)')' VV,xptfl:nhard,nsea,ireg '
     *  ,NHARD,NSEA,IREG
	IF(IREG.EQ.1)NHARD=0
	IF(IREG.EQ.1)NSEA=0
        NOMJE=NOMJE+NHARD
C
C                        SUBTRACT HARD SCATTERED X VALUES FROM DIQUARKS
C
        IF (NHARD.GE.1.AND.IMINIJ.EQ.1)THEN
        DO 71 IXX=NONUJ1,NONUJT
          JHKKPH(IXX)=IXVPR
          JHKKEX(IXX)=0
          JHKKE1(IXX)=0
          IF (XPVQ(IXVPR)-XJQ1(IXX).GT.XVTHR)THEN
            XPVQ(IXVPR)=XPVQ(IXVPR)-XJQ1(IXX)
            JHKKE1(IXX)=1
          ELSEIF (XPVD(IXVPR)-XJQ1(IXX).GT.XDTHR)THEN
            XPVD(IXVPR)=XPVD(IXVPR)-XJQ1(IXX)
            JHKKE1(IXX)=2
          ENDIF
   71   CONTINUE
        ENDIF
C
        IF(IPEV.GE.1)THEN
          PQP=GAMCM*PRMOM(3,INUCPR)+BGCM*PRMOM(4,INUCPR)
          PQE=GAMCM*PRMOM(4,INUCPR)+BGCM*PRMOM(3,INUCPR)
          PQPQ=GAMCM*PVQPZ+BGCM*PVQE
          PQEQ=GAMCM*PVQE+BGCM*PVQPZ
          PQPD=GAMCM*PVDQPZ+BGCM*PVDQE
          PQED=GAMCM*PVDQE+BGCM*PVDQPZ
        WRITE(6,1655)PRMOM(3,INUCPR),PRMOM(4,INUCPR),PQP,PQE,
     +  XPVQ(IXVPR),XPVD(IXVPR),IXVPR
        WRITE(6,1656)PVQPZ,PVQE,PQPQ,PQEQ
        WRITE(6,1657)PVDQPZ,PVDQE,PQPD,PQED
        ENDIF 
C
C
C                        SUBTRACT HARD SCATTERED X VALUES FROM DIQUARKS
C
        IF (NHARD.GE.1.AND.IMINIJ.EQ.1)THEN
        DO 771 IXX=NONUJ1,NONUJT
          JHKKTH(IXX)=IXVTA
          IF (JHKKE1(IXX).EQ.0)THEN
            JHKKEX(IXX)=0
            GO TO 771
          ENDIF
          IF (XTVQ(IXVTA)-XJQ2(IXX).GT. XVTHR) THEN
            XTVQ(IXVTA)=XTVQ(IXVTA)-XJQ2(IXX)
            JHKKEX(IXX)=1
            NOMJER=NOMJER+1
          ELSEIF(XTVD(IXVTA)-XJQ2(IXX).GT.XDTHR)THEN
            XTVD(IXVTA)=XTVD(IXVTA)-XJQ2(IXX)
            JHKKEX(IXX)=1
            NOMJER=NOMJER+1
          ELSE
            JHKKEX(IXX)=0
            IF (JHKKE1(IXX).EQ.1)THEN
              XPVQ(IXVPR)=XPVQ(IXVPR)+XJQ1(IXX)
            ELSEIF(JHKKE1(IXX).EQ.2)THEN
              XPVD(IXVPR)=XPVD(IXVPR)+XJQ1(IXX)
            ENDIF
          ENDIF
  771   CONTINUE
        ENDIF
C
C                        SUBTRACT SECONDARY SEA X VALUES FROM DIQUARKS
C
        IF (NSEA.GE.1)THEN
        IF(IPEV.GE.1)WRITE(6,'(A,3I10)')' VV,NSEA:NONUS1,NONUST '
     * ,NSEA,NONUS1,NONUST
        DO 271 IXX=NONUS1,NONUST
          JHKKPZ(IXX)=IXVPR
          JHKKSX(IXX)=0
          JHKKS1(IXX)=0
          IF (XPVQ(IXVPR)-XSQ1(IXX)-XSAQ1(IXX).GT.XVTHR)THEN
            XPVQ(IXVPR)=XPVQ(IXVPR)-XSQ1(IXX)-XSAQ1(IXX)
            JHKKS1(IXX)=1
          ELSEIF (XPVD(IXVPR)-XSQ1(IXX)-XSAQ1(IXX).GT.XDTHR)THEN
            XPVD(IXVPR)=XPVD(IXVPR)-XSQ1(IXX)-XSAQ1(IXX)
            JHKKS1(IXX)=2
          ENDIF
  271   CONTINUE
        ENDIF
C
        IXVTA=INTVV2(N)
        INUCTA=IFROVT(IXVTA)
C
C                        SUBTRACT SECONDARY SEA X VALUES FROM DIQUARKS
C
        IF (NSEA.GE.1)THEN
        DO 2771 IXX=NONUS1,NONUST
          JHKKTZ(IXX)=IXVTA
          IF (JHKKS1(IXX).EQ.0)THEN
            JHKKSX(IXX)=0
            GO TO 2775
          ENDIF
          IF (XTVQ(IXVTA)-XSQ2(IXX)-XSAQ2(IXX).GT. XVTHR) THEN
            XTVQ(IXVTA)=XTVQ(IXVTA)-XSQ2(IXX)-XSAQ2(IXX)
            JHKKSX(IXX)=1
C           NOMJER=NOMJER+1
          ELSEIF(XTVD(IXVTA)-XSQ2(IXX)-XSAQ2(IXX).GT.XDTHR)THEN
            XTVD(IXVTA)=XTVD(IXVTA)-XSQ2(IXX)-XSAQ2(IXX)
            JHKKSX(IXX)=1
C           NOMJER=NOMJER+1
          ELSE
            JHKKSX(IXX)=0
            IF (JHKKS1(IXX).EQ.1)THEN
              XPVQ(IXVPR)=XPVQ(IXVPR)+XSQ1(IXX)+XSAQ1(IXX)
            ELSEIF(JHKKS1(IXX).EQ.2)THEN
              XPVD(IXVPR)=XPVD(IXVPR)+XSQ1(IXX)+XSAQ1(IXX)
            ENDIF
          ENDIF
 2775   CONTINUE
        IF(IPEV.GE.1)WRITE(6,'(A,3I10)')' VV,ixx:jhkksx,jhkks1, '
     * ,IXX,JHKKSX(IXX),JHKKS1(IXX)
 2771   CONTINUE
        ENDIF
Cx
C
        XMAX1=XPVQ(IXVPR)+XPVD(IXVPR)-XVTHR-XDTHR
        XMAX2=XTVQ(IXVTA)+XTVD(IXVTA)-XVTHR-XDTHR
C
C
        IF(IPEV.GE.1)WRITE(6,'(A,2I5,3F9.3)')' KKEVVV,aft xptfl:n,nvv'
     *  ,N,NVV,XMAX1,XMAX2
C-----------------------------------------------------------------------
      IREJVV=0
C
C---------------------------drop recombined chain pairs
        IF(NCHVV1(N).EQ.99.AND.NCHVV2(N).EQ.99)GO TO 20
C
C***             4-MOMENTA OF PROJECTILE QUARK-DIQUARK PAIRS IN NN-CMS
        IXVPR=INTVV1(N)
        INUCPR=IFROVP(IXVPR)
C-------------------------------------Single Chain Option------
C                  NSICHA=0 : Normal 2 chain event
C                  NSICHA=1 : Single chain event
C     single chain  Meson  -Bayon  : Chain 1 ( q-qq  ) remains
C     single chain  Abaryon-Baryon : Chain 2 (aqaq-qq) remains
        NSICHA=0
        IF (ISICHA.EQ.1) THEN
          IF (NBPROJ.LE.0) THEN
C                                     Projectile particle
            IS1=INTVV1(N)
            IS2=INTVV2(N) 
            KHPROJ=KKPROJ(IS1)
C
C                                     Projectile Momentum(lab)
            PQP=GAMCM*PRMOM(3,INUCPR)+BGCM*PRMOM(4,INUCPR)
            PHPROJ=PQP
C                                     Original flavors
            IIQP=IPVQ(IS1)
            IIDP1=IPPV1(IS1)
            IIDP2=IPPV2(IS1)
            IIQT=ITVQ(IS2)
            IIDT1=ITTV1(IS2)
            IIDT2=ITTV2(IS2)
C                                    Target     particle
            IITSUM=IIQT+IIDT1+IIDT2
            IF(IITSUM.EQ.4)KHTARG=1
            IF(IITSUM.EQ.5)KHTARG=8
C           KHTARG=KKTARG(IS2)
C                                     Single chain probability  
            SICHAP=PHNSCH(KHPROJ,KHTARG,PHPROJ)
            IF (RNDM(V).LE.SICHAP)NSICHA=1
C           IF (NBPROJ.EQ.-1)NSICHA=0
C                                     Single chain quark flavors
            AAAAA=SCHQUA(JQFSC1,JQFSC2,JQBSC1,JQBSC2)
            IF(ISCH.LE.20)
     +      WRITE(6,'(A,3I5,2F10.3,10I5)')' KKEVVV Single chain ',
     +          NSICHA,KHPROJ,KHTARG,PHPROJ,SICHAP,
     +          IIQP,IIDP1,IIDP2,IIQT,IIDT1,IIDT2,
     +          JQFSC1,JQFSC2,JQBSC1,JQBSC2 
            IF(NBPROJ.EQ.0.AND.NSICHA.EQ.1)THEN
C                         Correct x-values and quark flavors
              NCHVV2(N)=99
              XPVQ(IXVPR)=XPVQ(IXVPR)+XPVD(IXVPR)
              XPVD(IXVPR)=0
              XTVD(IXVTA)=XTVD(IXVTA)+XTVQ(IXVTA)
              XTVQ(IXVTA)=0
C             IPVQ(IS1)=JQFSC1
              ITTV1(IS2)=JQBSC1
              ITTV2(IS2)=JQBSC2
            ELSEIF(NBPROJ.EQ.-1.AND.NSICHA.EQ.1)THEN
C                         Correct x-values and quark flavors
C                                                       ?
              XPVD(IXVPR)=XPVQ(IXVPR)+XPVD(IXVPR)
              XPVQ(IXVPR)=0
              XTVD(IXVTA)=XTVD(IXVTA)+XTVQ(IXVTA)
              XTVQ(IXVTA)=0
              NCHVV1(N)=99
              IPPV1(IS1)=JQFSC1
              IPPV2(IS1)=JQFSC2
              ITTV1(IS2)=JQBSC1
              ITTV2(IS2)=JQBSC2
            ENDIF         
          ENDIF 
        ENDIF
C-------------------------------------Single Chain Option------
C
        PVQPX=XPVQ(IXVPR)*PRMOM(1,INUCPR)
        PVQPY=XPVQ(IXVPR)*PRMOM(2,INUCPR)
        PVQPZ=XPVQ(IXVPR)*PRMOM(3,INUCPR)
        PVQE =XPVQ(IXVPR)*PRMOM(4,INUCPR)
C
        PVDQPX=XPVD(IXVPR)*PRMOM(1,INUCPR)
        PVDQPY=XPVD(IXVPR)*PRMOM(2,INUCPR)
        PVDQPZ=XPVD(IXVPR)*PRMOM(3,INUCPR)
        PVDQE =XPVD(IXVPR)*PRMOM(4,INUCPR)
        IF(IPEV.GE.1)THEN
          PQP=GAMCM*PRMOM(3,INUCPR)+BGCM*PRMOM(4,INUCPR)
          PQE=GAMCM*PRMOM(4,INUCPR)+BGCM*PRMOM(3,INUCPR)
          PQPQ=GAMCM*PVQPZ+BGCM*PVQE
          PQEQ=GAMCM*PVQE+BGCM*PVQPZ
          PQPD=GAMCM*PVDQPZ+BGCM*PVDQE
          PQED=GAMCM*PVDQE+BGCM*PVDQPZ
        WRITE(6,1655)PRMOM(3,INUCPR),PRMOM(4,INUCPR),PQP,PQE,
     +  XPVQ(IXVPR),XPVD(IXVPR),IXVPR
 1655     FORMAT(' vv PQP,PQE ',6F12.5,I5)
        WRITE(6,1656)PVQPZ,PVQE,PQPQ,PQEQ
 1656     FORMAT(' vv PQPQ,PQEQ ',4F12.5)
        WRITE(6,1657)PVDQPZ,PVDQE,PQPD,PQED
 1657     FORMAT(' vv PQPD,PQED ',4F12.5)
        ENDIF 
C
C***                 4-MOMENTA OF TARGET QUARK-DIQUARK PAIRS IN NN-CMS
        IXVTA=INTVV2(N)
        INUCTA=IFROVT(IXVTA)
C
        TVQPX=XTVQ(IXVTA)*TAMOM(1,INUCTA)
        TVQPY=XTVQ(IXVTA)*TAMOM(2,INUCTA)
        TVQPZ=XTVQ(IXVTA)*TAMOM(3,INUCTA)
        TVQE =XTVQ(IXVTA)*TAMOM(4,INUCTA)
 
C
        TVDQPX=XTVD(IXVTA)*TAMOM(1,INUCTA)
        TVDQPY=XTVD(IXVTA)*TAMOM(2,INUCTA)
        TVDQPZ=XTVD(IXVTA)*TAMOM(3,INUCTA)
        TVDQE =XTVD(IXVTA)*TAMOM(4,INUCTA)
        IF(IPEV.GE.1)THEN
          TQP=GAMCM*TAMOM(3,INUCTA)+BGCM*TAMOM(4,INUCTA)
          TQE=GAMCM*TAMOM(4,INUCTA)+BGCM*TAMOM(3,INUCTA)
          TQPQ=GAMCM*TVQPZ+BGCM*TVQE
          TQEQ=GAMCM*TVQE+BGCM*TVQPZ
          TQPD=GAMCM*TVDQPZ+BGCM*TVDQE
          TQED=GAMCM*TVDQE+BGCM*TVDQPZ
        WRITE(6,1455)TAMOM(3,INUCTA),TAMOM(4,INUCTA),TQP,TQE
 1455     FORMAT(' vv TQP,TQE ',4F12.5)
        WRITE(6,1456)TVQPZ,TVQE,TQPQ,TQEQ
 1456     FORMAT(' vv TQPQ,TQEQ ',4F12.5)
        WRITE(6,1457)TVDQPZ,TVDQE,TQPD,TQED
 1457     FORMAT(' vv TQPD,TQED ',4F12.5)
          WRITE(6,1355)XPVQ(IXVPR),XPVD(IXVPR),XTVQ(IXVTA),
     *          XTVD(IXVTA),PRMOM(4,INUCPR),TAMOM(4,INUCTA)
 1355     FORMAT(' VV  xpq.xpd,xtq,xtd,ep,et ',6F12.5)
        ENDIF 
C                                               j.r.6.5.93
C
C                     multiple scattering of valence quark chain ends
C
      IF(IT.GT.1)THEN
      ITNU=IP+INUCTA
      RTIX=(VHKK(1,ITNU))*1.E12-BIMPAC*0.1
      RTIY=VHKK(2,ITNU)*1.E12
      RTIZ=VHKK(3,ITNU)*1.E12
      CALL CROMSC(PVQPX,PVQPY,PVQPZ,PVQE,RTIX,RTIY,RTIZ,
     *            PVQNX,PVQNY,PVQNZ,PVQNE,1)
      PVQPX=PVQNX
      PVQPY=PVQNY
      PVQPZ=PVQNZ
      PVQE=PVQNE
      CALL CROMSC(PVDQPX,PVDQPY,PVDQPZ,PVDQE,RTIX,RTIY,RTIZ,
     *            PVDQNX,PVDQNY,PVDQNZ,PVDQNE,2)
      PVDQPX=PVDQNX
      PVDQPY=PVDQNY
      PVDQPZ=PVDQNZ
      PVDQE=PVDQNE
C                                                ---------
 
C                                               j.r.6.5.93
C
C                     multiple scattering of valence quark chain ends
C
      ITNU=IP+INUCTA
      RTIX=(VHKK(1,ITNU))*1.E12-BIMPAC*0.1
      RTIY=VHKK(2,ITNU)*1.E12
      RTIZ=VHKK(3,ITNU)*1.E12
      CALL CROMSC(TVQPX,TVQPY,TVQPZ,TVQE,RTIX,RTIY,RTIZ,
     *            TVQNX,TVQNY,TVQNZ,TVQNE,3)
      TVQPX=TVQNX
      TVQPY=TVQNY
      TVQPZ=TVQNZ
      TVQE=TVQNE
      CALL CROMSC(TVDQPX,TVDQPY,TVDQPZ,TVDQE,RTIX,RTIY,RTIZ,
     *            TVDQNX,TVDQNY,TVDQNZ,TVDQNE,4)
      TVDQPX=TVDQNX
      TVDQPY=TVDQNY
      TVDQPZ=TVDQNZ
      TVDQE=TVDQNE
      ENDIF
 
C                                                j.r.10.5.93
       IF(IP.GE.1)GO TO 1779
        PVQPZ2=PVQE**2-PVQPX**2-PVQPY**2
        IF(PVQPZ2.GE.0.)THEN
          PVQPZ=SQRT(PVQPZ2)
        ELSE
          PVQPX=0.
          PVQPY=0.
          PVQPZ=PVQE
        ENDIF
C
        PDQPZ2=PVDQE**2-PVDQPX**2-PVDQPY**2
        IF(PDQPZ2.GE.0.)THEN
          PVDQPZ=SQRT(PDQPZ2)
        ELSE
          PVDQPX=0.
          PVDQPY=0.
          PVDQPZ=PVDQE
        ENDIF
C
        TVQPZ2=TVQE**2-TVQPX**2-TVQPY**2
        IF(TVQPZ2.GE.0.)THEN
          TVQPZ=-SQRT(TVQPZ2)
        ELSE
          TVQPX=0.
          TVQPY=0.
          TVQPZ=TVQE
        ENDIF
C
        TDQPZ2=TVDQE**2-TVDQPX**2-TVDQPY**2
        IF(TDQPZ2.GE.0.)THEN
          TVDQPZ=-SQRT(TDQPZ2)
        ELSE
          TVDQPX=0.
          TVDQPY=0.
          TVDQPZ=PVDQE
        ENDIF
 1779  CONTINUE
C                                            ----------------
***  SAMPLE PARTON-PT VALUES / DETERMINE PARTON 4-MOMENTA AND CHAIN MAS
C***                            IN THE REST FRAME DEFINED ABOVE
C
C
C                                      changej.r.6.5.93
        PTXSQ1=0.
        PTXSA1=0.
        PTXSQ2=0.
        PTXSA2=0.
        PTYSQ1=0.
        PTYSA1=0.
        PTYSQ2=0.
        PTYSA2=0.
        PTXSQ1=PVQPX
        PTXSA1=PVDQPX
        PTXSQ2=TVQPX
        PTXSA2=TVDQPX
        PTYSQ1=PVQPY
        PTYSA1=PVDQPY
        PTYSQ2=TVQPY
        PTYSA2=TVDQPY
        PLQ1=PVQPZ
        PLAQ1=PVDQPZ
        PLQ2=TVQPZ
        PLAQ2=TVDQPZ
        EQ1=PVQE
        EAQ1=PVDQE
        EQ2=TVQE
        EAQ2=TVDQE
C                                       ---------------
        IF(IPEV.GE.1) THEN
          WRITE(6,1050)  PTXSQ1,PTYSQ1,
     +    PLQ1,EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,
     +    PTXSA2,PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1
          BPLQ1=GAMCM*PLQ1+BGCM*EQ1
          BEQ1=GAMCM*EQ1+BGCM*PLQ1
          BPLAQ1=GAMCM*PLAQ1+BGCM*EAQ1
          BEAQ1=GAMCM*EAQ1+BGCM*PLAQ1
          BPLQ2=GAMCM*PLQ2+BGCM*EQ2
          BEQ2=GAMCM*EQ2+BGCM*PLQ2
          BPLAQ2=GAMCM*PLAQ2+BGCM*EAQ2
          BEAQ2=GAMCM*EAQ2+BGCM*PLAQ2
          WRITE(6,1050)  PTXSQ1,PTYSQ1,
     +    BPLQ1,BEQ1,PTXSA1,PTYSA1,BPLAQ1,BEAQ1,
     +     PTXSQ2,PTYSQ2,BPLQ2,BEQ2,
     +    PTXSA2,PTYSA2,BPLAQ2,BEAQ2,
     +     AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1
        ENDIF
        IKVALA=1
        NSELPT=0
        IF(NSICHA.EQ.0)THEN
          IF(NBPROJ.GE.0)THEN
       IF(IOUXEV.GE.6)WRITE(6,'(A)')' KKEVVV call SELPT'
            CALL SELPT( PTXSQ1,PTYSQ1,PLQ1,EQ1,
     +       PTXSA1,PTYSA1,PLAQ1,EAQ1,
     +       PTXSQ2,PTYSQ2,PLQ2,EQ2,
     +       PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +       AMCH1,AMCH2,
     +       IREJ,IKVALA,PTTQ1,PTTA1,
     *       PTTQ2,PTTA2,
     +       NSELPT)
          ENDIF
          IF(NBPROJ.EQ.-1)THEN
       IF(IOUXEV.GE.6)WRITE(6,'(A)')' KKEVVV call SELPT'
            CALL SELPT( PTXSQ1,PTYSQ1,PLQ1,EQ1,
     +       PTXSA1,PTYSA1,PLAQ1,EAQ1,
     +       PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +       PTXSQ2,PTYSQ2,PLQ2,EQ2,
     +       AMCH1,AMCH2,
     +       IREJ,IKVALA,PTTQ1,PTTA1,
     *       PTTQ2,PTTA2,
     +       NSELPT)
          ENDIF
        ENDIF
C-------------------------------------Single Chain Option------
        IF(NSICHA.EQ.1.AND.NBPROJ.EQ.0)THEN
          CALL SELPTS( PTXSQ1,PTYSQ1,
     +        PLQ1,EQ1,PTXSA2,
     +        PTYSA2,PLAQ2,EAQ2, AMCH1,IREJ,IKVALA,PTTQ1)
        ENDIF
        IF(NSICHA.EQ.1.AND.NBPROJ.EQ.-1)THEN
          CALL SELPTS( PTXSA1,PTYSA1,
     +        PLAQ1,EAQ1,PTXSA2,
     +        PTYSA2,PLAQ2,EAQ2, AMCH2,IREJ,IKVALA,PTTA1)
        ENDIF
C-------------------------------------Single Chain Option------
        IF(IPEV.GE.1) THEN
          WRITE(6,1050)  PTXSQ1,PTYSQ1,
     +    PLQ1,EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,
     +    PTXSA2,PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1
        ENDIF
C
        IF (IREJ.EQ.1) THEN
          IRVV13=IRVV13 + 1
          IF(IPEV.GE.1) THEN
            WRITE(6,1100) IRVV13
            WRITE(6,1050)  PTXSQ1,
     +      PTYSQ1,PLQ1,EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,PTYSQ2,
     +      PLQ2,EQ2,PTXSA2,PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2,IREJ,IKVALA,
     +      PTTQ1,PTTA1
          ENDIF
                                                                GO TO 10
        ENDIF
C
C***  4-MOMENTA OF CHAINS IN THIS FRAME
C
        IF(NBPROJ.GE.0)THEN
          PTXCH1=PTXSQ1 + PTXSA2
          PTYCH1=PTYSQ1 + PTYSA2
          PTZCH1=PLQ1 + PLAQ2
          ECH1=EQ1 + EAQ2
          PTXCH2=PTXSQ2 + PTXSA1
          PTYCH2=PTYSQ2 + PTYSA1
          PTZCH2=PLQ2 + PLAQ1
          ECH2=EQ2 + EAQ1
        ENDIF
        IF(NBPROJ.EQ.-1)THEN
          PTXCH1=PTXSQ1 + PTXSQ2
          PTYCH1=PTYSQ1 + PTYSQ2
          PTZCH1=PLQ1 + PLQ2
          ECH1=EQ1 + EQ2
          PTXCH2=PTXSA2 + PTXSA1
          PTYCH2=PTYSA2 + PTYSA1
          PTZCH2=PLAQ2 + PLAQ1
          ECH2=EAQ2 + EAQ1
        ENDIF
        AMMM=SQRT((ECH1+ECH2)**2-(PTXCH1+PTXCH2)**2
     +            -(PTYCH1+PTYCH2)**2-(PTZCH1+PTZCH2)**2)
C
        IF (IPEV.GE.6)WRITE(6,1040) IREJ,
     +  AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1, AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2
 
C
C  REPLACE SMALL MASS CHAINS BY OCTETT OR DECUPLETT BARYONS
C  FIRST FOR CHAIN 1
C                    PROJ VAL-QUARK  - TAR DIQUARK   FOR MESONS/BARYONS
C                    PROJ VAL-AQUARK - TAR QUARK     FOR ANTIBARYONS
C
C     IF(NBPROJ.LE.100)GO TO 5559
      IF(NSICHA.EQ.0)THEN 
        IF(NBPROJ.GE.0) THEN
          CALL COBCMA(IPVQ(IXVPR),ITTV1(IXVTA),ITTV2(IXVTA), IJNCH1,
     +    NNCH1,IREJ,AMCH1,AMCH1N,1)
        ELSE
          CALL COMCMA(ITVQ(IXVTA),IPVQ(IXVPR), IJNCH1,NNCH1,IREJ,AMCH1,
     +    AMCH1N)
        ENDIF
C***                      MASS BELOW OCTETT BARYON MASS (MESON/BARYON)
C***                      MASS BELOW PSEUDOSCALAR MASS  (ANTIBARYON)
        IF(IREJ.EQ.1) THEN
          IRVV11=IRVV11 + 1
          IF(IPEV.GE.1) THEN
            WRITE(6,1110) IRVV11
            WRITE(6,1060) IPVQ(IXVPR),ITTV1(IXVTA),ITTV2(IXVTA), IJNCH1,
     +      NNCH1,IREJ, XPVQ(IXVPR),XPVD(IXVPR),XPVQCM,XPVDCM, XTVQ
     +      (IXVTA),XTVD(IXVTA),AMCH1,AMCH1N
 
          ENDIF
                                                                 GOTO 10
        ENDIF
C                                 CORRECT KINEMATICS FOR CHAIN 1
C***                MOMENTUM CORRECTION FOR CHANGED MASS OF CHAIN 1
        IF(NNCH1.NE.0) THEN
          IF(NBPROJ.GE.0) THEN
            CALL CORMOM(AMCH1,AMCH2,AMCH1N,AMCH2N,
     +      PTXSQ1,PTYSQ1,PLQ1,EQ1,
     +      PTXSA1,PTYSA1,PLAQ1,EAQ1,
     +      PTXSQ2,PTYSQ2,PLQ2,EQ2,
     +      PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +      PTXCH1,PTYCH1,PTZCH1,ECH1,
     +      PTXCH2,PTYCH2,PTZCH2,ECH2,IREJ)
          AMCH2=AMCH2N
          ELSE
            CALL CORMOM(AMCH1,AMCH2,AMCH1N,AMCH2N, 
     +      PTXSQ1,PTYSQ1,PLQ1,EQ1,
     +      PTXSA1,PTYSA1,PLAQ1,EAQ1,
     +      PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +      PTXSQ2,PTYSQ2,PLQ2,EQ2,
     +      PTXCH1,PTYCH1,PTZCH1,ECH1,
     +      PTXCH2,PTYCH2,PTZCH2,ECH2,IREJ)
            AMCH2=AMCH2N
          ENDIF
          IF(IREJ.EQ.1)THEN
	    IF(IPEV.GE.1)WRITE(6,'(A)')' vv cormom rej'
	  GO TO 10
	  ENDIF
C
          IF (IPEV.GE.1)WRITE(6,1040) IREJ,
     +    AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1, AMCH2,PTXCH2,PTYCH2,PTZCH2,
     +    ECH2
        ENDIF
C
C  SECOND FOR CHAIN 2:
C             PROJ VAL-DIQUARK  - TAR VAL-QUARK     FOR INC. BARYONS
C             PROJ VAL-AQUARK   - TAR VAL-QUARK     FOR INC. MESONS
C             PROJ VAL-ADIQUARK - TAR VAL-DIQUARK   FOR INC. ANTIBARYON
C
C	IF(NBPROJ.LE.100)GO TO 5557
        IF(NBPROJ.GT.0) THEN
          CALL COBCMA(ITVQ(IXVTA),IPPV1(IXVPR),IPPV2(IXVPR), IJNCH2,
     +    NNCH2,IREJ,AMCH2,AMCH2N,2)
        ELSEIF(NBPROJ.EQ.0) THEN
          CALL COMCMA(ITVQ(IXVTA),IPPV1(IXVPR), IJNCH2,NNCH2,IREJ,AMCH2,
     +    AMCH2N)
        ELSE
          CALL COMCM2(ITTV1(IXVTA),ITTV2(IXVTA), IPPV1(IXVPR),IPPV2
     +    (IXVPR), NNCH2,IREJ,AMCH2)
 
        ENDIF
C***                  MASS BELOW OCTETT BARYON/PSEUDOSCALAR MESON MASS
C                     OR INCONSISTENT QUARK FLAVORS FOR ANTIBARYONS
C
C5557   CONTINUE

        IF(IREJ.EQ.1) THEN
          IRVV12=IRVV12 + 1
          IF(IPEV.GE.1) THEN
            WRITE(6,1120) IRVV12
            WRITE(6,1080) IPPV1(IXVPR),IPPV2(IXVPR),ITTV1(IXVTA), ITTV2
     +      (IXVTA),IJNCH2,NNCH2,IREJ, XPVQ(IXVPR),XPVD(IXVPR),XPVQCM,
     +      XPVDCM, XTVQ(IXVTA),XTVD(IXVTA),XTVQCM,XTVDCM, AMCH2,AMCH2N
          ENDIF
                                                                 GOTO 10
        ENDIF
       ENDIF
C5559   CONTINUE
       IF(NSICHA.EQ.0)THEN
        IF(NNCH2.NE.0) THEN
          AMCH2=AMCH2N
         AMMM=SQRT((ECH1+ECH2)**2-(PTXCH1+PTXCH2)**2
     +            -(PTYCH1+PTYCH2)**2-(PTZCH1+PTZCH2)**2)
        EEE=ECH1+ECH2
        PXXX=PTXCH1+PTXCH2
        PYYY=PTYCH1+PTYCH2
        PZZZ=PTZCH1+PTZCH2
        GAMMM=EEE/(AMMM+1.E-4)
        BGGGX=PXXX/(AMMM+1.E-4)
        BGGGY=PYYY/(AMMM+1.E-4)
        BGGGZ=PZZZ/(AMMM+1.E-4)
C                        TRANSFORM BOTH CHAINS  INTO TWO CHAIN-CMS
C
C                                   4-MOMENTA OF CHAINS
        CALL DALTRA(GAMMM,-BGGGX,-BGGGY,-BGGGZ,
     +  PTXCH1,PTYCH1,PTZCH1,ECH1,
     +  PPPCH1, QTXCH1,QTYCH1,QTZCH1,QECH1)
C
        CALL DALTRA(GAMMM,-BGGGX,-BGGGY,-BGGGZ,
     +  PTXCH2,PTYCH2,PTZCH2,ECH2,
     +  PPPCH2, QTXCH2,QTYCH2,QTZCH2,QECH2)
C
C                                IF AMCH2 CHANGED IN COBCMA/COMCMA
C                                CORRESPONDING REPLACEMENT IN CORMOM
          NORIG=21
C	IF(NBPROJ.LE.100)GO TO 5558
          CALL CORVAL(AMMM,IREJ,AMCH1,AMCH2, QTXCH1,QTYCH1,QTZCH1,QECH1,
     +    QTXCH2,QTYCH2,QTZCH2,QECH2,NORIG)
C5558   CONTINUE
C                        TRANSFORM BOTH CHAINS  INTO TWO CHAIN-CMS
C
C                                   4-MOMENTA OF CHAINS
        CALL DALTRA(GAMMM,BGGGX,BGGGY,BGGGZ, QTXCH1,QTYCH1,QTZCH1,QECH1,
     +  PPPCH1, PTXCH1,PTYCH1,PTZCH1,ECH1)
C
        CALL DALTRA(GAMMM,BGGGX,BGGGY,BGGGZ, QTXCH2,QTYCH2,QTZCH2,QECH2,
     +  PPPCH2, PTXCH2,PTYCH2,PTZCH2,ECH2)
C
C
          IF(IPEV.GE.3) THEN
            WRITE(6,'(A/3(1PE15.4),3I5)')
     +      ' VV - CALL CORVAL: AMMM, AMCH1, AMCH2, NNCH1, NNCH2, IREJ',
     +      AMMM, AMCH1, AMCH2, NNCH1, NNCH2, IREJ
          ENDIF
          IF(IREJ.EQ.1) THEN
	    IF(IPEV.GE.1)WRITE(6,'(A)')' vv14 rej.'
C                           AMCH1N + AMCH2N > AMMM - 0.2
C                           REJECT EVENT
            IRVV14=IRVV14+1
                                                                 GOTO 10
          ENDIF
        ENDIF
       ENDIF
       QTXCH1=PTXCH1
       QTYCH1=PTYCH1
       QTZCH1=PTZCH1
       QECH1=ECH1
       QTXCH2=PTXCH2
       QTYCH2=PTYCH2
       QTZCH2=PTZCH2
       QECH2=ECH2
       PQVVA1(N,1)=PTXSQ1
       PQVVA1(N,2)=PTYSQ1
       PQVVA1(N,3)=PLQ1
       PQVVA1(N,4)=EQ1
       PQVVA2(N,1)=PTXSA2
       PQVVA2(N,2)=PTYSA2
       PQVVA2(N,3)=PLAQ2
       PQVVA2(N,4)=EAQ2
       PQVVB1(N,1)=PTXSQ2
       PQVVB1(N,2)=PTYSQ2
       PQVVB1(N,3)=PLQ2
       PQVVB1(N,4)=EQ2
       PQVVB2(N,1)=PTXSA1
       PQVVB2(N,2)=PTYSA1
       PQVVB2(N,3)=PLAQ1
       PQVVB2(N,4)=EAQ1
C-------------------
C                                PUT V-V CHAIN ENDS AND CHAINS
C                                                 INTO /HKKEVT/
C                                 MOMENTA IN NN-CMS
C                                 POSITION OF ORIGINAL NUCLEONS
C
C                                 FLAG FOR VV-CHAIN ENDS
C                                            PROJECTILE: ISTHKK=121
C                                            TARGET:     ISTHKK=122
C                                      FOR VV-CHAINS     ISTHKK=3
        IHKKPD=JHKKPV(IXVPR)
        IHKKPO=IHKKPD - 1
C***                                 hjm 27/08/90
        IF(NBPROJ.GE.0) THEN
          IHKKTD=JHKKTV(IXVTA)
          IHKKTO=IHKKTD - 1
        ELSE
          IHKKTO=JHKKTV(IXVTA)
          IHKKTD=IHKKTO - 1
        ENDIF
C
C
        IF (IPEV.GT.3) THEN
          WRITE(6,1000) IXVPR,INUCPR,IHKKPO,IHKKPD
 1000 FORMAT (' IXVPR,INUCPR,IHKKPO,IHKKPD ',5I5)
          WRITE(6,1010) IXVTA,INUCTA,IHKKTO,IHKKTD
 1010 FORMAT (' IXVTA,INUCTA,IHKKTO,IHKKTD ',5I5)
        ENDIF
C                                     CHAIN 1 PROJECTILE QUARK
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        ISTHKK(NHKK)=121
        IDHKK(NHKK)=IDHKK(IHKKPO)
        JMOHKK(1,NHKK)=IHKKPO
        JMOHKK(2,NHKK)=JMOHKK(1,IHKKPO)
        JDAHKK(1,NHKK)=NHKK+2
        JDAHKK(2,NHKK)=NHKK+2
        PHKK(1,NHKK)=PQVVA1(N,1)
        PHKK(2,NHKK)=PQVVA1(N,2)
        PHKK(3,NHKK)=PQVVA1(N,3)
        PHKK(4,NHKK)=PQVVA1(N,4)
        PHKK(5,NHKK)=0.
      CALL QINNUC(XXPP,YYPP)
        VHKK(1,NHKK)=VHKK(1,IHKKPO)+XXPP
        VHKK(2,NHKK)=VHKK(2,IHKKPO)+YYPP
        VHKK(3,NHKK)=VHKK(3,IHKKPO)
        VHKK(4,NHKK)=VHKK(4,IHKKPO)
C
        IF (IPHKK.GE.2) WRITE(6,1020) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     +  JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +  (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
 1020 FORMAT (I6,I4,5I6,9E10.2)
C                                     CHAIN 1 TARGET DIQUARK
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        ISTHKK(NHKK)=122
        IDHKK(NHKK)=IDHKK(IHKKTD)
        JMOHKK(1,NHKK)=IHKKTD
        JMOHKK(2,NHKK)=JMOHKK(1,IHKKTD)
        JDAHKK(1,NHKK)=NHKK+1
        JDAHKK(2,NHKK)=NHKK+1
        PHKK(1,NHKK)=PQVVA2(N,1)
        PHKK(2,NHKK)=PQVVA2(N,2)
        PHKK(3,NHKK)=PQVVA2(N,3)
        PHKK(4,NHKK)=PQVVA2(N,4)
        PHKK(5,NHKK)=0.
      CALL QINNUC(XXPP,YYPP)
        VHKK(1,NHKK)=VHKK(1,IHKKTD)+XXPP
        VHKK(2,NHKK)=VHKK(2,IHKKTD)+YYPP
        VHKK(3,NHKK)=VHKK(3,IHKKTD)
        VHKK(4,NHKK)=VHKK(4,IHKKTD)
C
        IF (IPHKK.GE.2) WRITE(6,1020) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     +  JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +  (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
C
C                                     CHAIN 1 BEFORE FRAGMENTATION
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        ISTHKK(NHKK)=3
        IDHKK(NHKK)=88888+NNCH1
	IF(NCHVV1(N).EQ.99)IDHKK(NHKK)=77777
        JMOHKK(1,NHKK)=NHKK-2
        JMOHKK(2,NHKK)=NHKK-1
        PHKK(1,NHKK)=QTXCH1
        PHKK(2,NHKK)=QTYCH1
        PHKK(3,NHKK)=QTZCH1
        PHKK(4,NHKK)=QECH1
        PHKK(5,NHKK)=AMCH1
C                             POSITION OF CREATED CHAIN IN LAB
C                             =POSITION OF TARGET NUCLEON
C                             TIME OF CHAIN CREATION IN LAB
C                             =TIME OF PASSAGE OF PROJECTILE
C                              NUCLEUS AT POSITION OF TAR. NUCLEUS
        VHKK(1,NHKK)= VHKK(1,NHKK-1)
        VHKK(2,NHKK)= VHKK(2,NHKK-1)
        VHKK(3,NHKK)= VHKK(3,NHKK-1)
        VHKK(4,NHKK)=VHKK(3,NHKK)/BETP-VHKK(3,NHKK-2)/BGAMP
        MHKKVV(N)=NHKK
        IF (IPROJK.EQ.1)THEN
          WHKK(1,NHKK)= VHKK(1,NHKK-2)
          WHKK(2,NHKK)= VHKK(2,NHKK-2)
          WHKK(3,NHKK)= VHKK(3,NHKK-2)
          WHKK(4,NHKK)=VHKK(4,NHKK)*GAMP-VHKK(3,NHKK)*BGAMP
          IF (IPHKK.GE.2) WRITE(6,1020) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     +    JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +    (PHKK(KHKK,NHKK),KHKK=1,5), (WHKK(KHKK,NHKK),KHKK=1,4)
 
        ENDIF
C
        IF (IPHKK.GE.1) WRITE(6,1020) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     +  JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +  (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
C
C
C                                     CHAIN 2 PROJECTILE DIQUARK
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        ISTHKK(NHKK)=121
        IDHKK(NHKK)=IDHKK(IHKKPD)
        JMOHKK(1,NHKK)=IHKKPD
        JMOHKK(2,NHKK)=JMOHKK(1,IHKKPD)
        JDAHKK(1,NHKK)=NHKK+2
        JDAHKK(2,NHKK)=NHKK+2
        PHKK(1,NHKK)=PQVVB1(N,1)
        PHKK(2,NHKK)=PQVVB1(N,2)
        PHKK(3,NHKK)=PQVVB1(N,3)
        PHKK(4,NHKK)=PQVVB1(N,4)
        PHKK(5,NHKK)=0.
      CALL QINNUC(XXPP,YYPP)
        VHKK(1,NHKK)=VHKK(1,IHKKPD)+XXPP
        VHKK(2,NHKK)=VHKK(2,IHKKPD)+YYPP
        VHKK(3,NHKK)=VHKK(3,IHKKPD)
        VHKK(4,NHKK)=VHKK(4,IHKKPD)
C
        IF (IPHKK.GE.2) WRITE(6,1020) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     +  JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +  (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
C                                     CHAIN 2 TARGET QUARK
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        ISTHKK(NHKK)=122
        IDHKK(NHKK)=IDHKK(IHKKTO)
        JMOHKK(1,NHKK)=IHKKTO
        JMOHKK(2,NHKK)=JMOHKK(1,IHKKTO)
        JDAHKK(1,NHKK)=NHKK+1
        JDAHKK(2,NHKK)=NHKK+1
        PHKK(1,NHKK)=PQVVB2(N,1)
        PHKK(2,NHKK)=PQVVB2(N,2)
        PHKK(3,NHKK)=PQVVB2(N,3)
        PHKK(4,NHKK)=PQVVB2(N,4)
        PHKK(5,NHKK)=0.
      CALL QINNUC(XXPP,YYPP)
        VHKK(1,NHKK)=VHKK(1,IHKKTO)+XXPP
        VHKK(2,NHKK)=VHKK(2,IHKKTO)+YYPP
        VHKK(3,NHKK)=VHKK(3,IHKKTO)
        VHKK(4,NHKK)=VHKK(4,IHKKTO)
C
        IF (IPHKK.GE.2) WRITE(6,1020) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     +  JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +  (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
C
C                                     CHAIN 2 BEFORE FRAGMENTATION
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        ISTHKK(NHKK)=3
        IDHKK(NHKK)=88888+NNCH2
	IF(NCHVV2(N).EQ.99)IDHKK(NHKK)=77777
        JMOHKK(1,NHKK)=NHKK-2
        JMOHKK(2,NHKK)=NHKK-1
        PHKK(1,NHKK)=QTXCH2
        PHKK(2,NHKK)=QTYCH2
        PHKK(3,NHKK)=QTZCH2
        PHKK(4,NHKK)=QECH2
        PHKK(5,NHKK)=AMCH2
C                             POSITION OF CREATED CHAIN IN LAB
C                             =POSITION OF TARGET NUCLEON
C                             TIME OF CHAIN CREATION IN LAB
C                             =TIME OF PASSAGE OF PROJECTILE
C                              NUCLEUS AT POSITION OF TAR. NUCLEUS
        VHKK(1,NHKK)= VHKK(1,NHKK-1)
        VHKK(2,NHKK)= VHKK(2,NHKK-1)
        VHKK(3,NHKK)= VHKK(3,NHKK-1)
        VHKK(4,NHKK)=VHKK(3,NHKK)/BETP-VHKK(3,NHKK-2)/BGAMP
        MHKKVV(N)=NHKK
        IF (IPROJK.EQ.1)THEN
          WHKK(1,NHKK)= VHKK(1,NHKK-2)
          WHKK(2,NHKK)= VHKK(2,NHKK-2)
          WHKK(3,NHKK)= VHKK(3,NHKK-2)
          WHKK(4,NHKK)=VHKK(4,NHKK)*GAMP-VHKK(3,NHKK)*BGAMP
          IF (IPHKK.GE.1) WRITE(6,1020) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     +    JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +    (PHKK(KHKK,NHKK),KHKK=1,5), (WHKK(KHKK,NHKK),KHKK=1,4)
 
        ENDIF
C
        IF (IPHKK.GE.1) WRITE(6,1020) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     +  JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +  (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
C
C
C  NOW WE HAVE AN ACCEPTABLE VALENCE-VALENCE EVENT
C  AND PUT IT INTO THE HISTOGRAM
C
        AMCVV1(N)=AMCH1
        AMCVV2(N)=AMCH2
        IF(AMCH1.GT.0.D0)THEN
        GACVV1(N)=QECH1/AMCH1
        BGXVV1(N)=QTXCH1/AMCH1
        BGYVV1(N)=QTYCH1/AMCH1
        BGZVV1(N)=QTZCH1/AMCH1
        ELSE
        GACVV1(N)=0
        BGXVV1(N)=0
        BGYVV1(N)=0
        BGZVV1(N)=0
        ENDIF
        IF(AMCH2.GT.0.D0)THEN
        GACVV2(N)=QECH2/AMCH2
        BGXVV2(N)=QTXCH2/AMCH2
        BGYVV2(N)=QTYCH2/AMCH2
        BGZVV2(N)=QTZCH2/AMCH2
        ELSE
        GACVV2(N)=0
        BGXVV2(N)=0
        BGYVV2(N)=0
        BGZVV2(N)=0
        ENDIF
       IF(NSICHA.EQ.0)THEN
        NCHVV1(N)=NNCH1
        NCHVV2(N)=NNCH2
       ENDIF
C-------------------------------------Single Chain Option------
       IF(NSICHA.EQ.1.AND.IBPROJ.EQ.0)THEN
         NCHVV1(N)=0
         NCHVV2(N)=99
       ENDIF
       IF(NSICHA.EQ.1.AND.IBPROJ.EQ.-1)THEN
         NCHVV1(N)=99
         NCHVV2(N)=0
       ENDIF
C-------------------------------------Single Chain Option------
C-----------------------------------------------------------------
        IJCVV1(N)=IJNCH1
        IJCVV2(N)=IJNCH2
C
        IF (IPEV.GE.6)WRITE(6,1030) N, XPVQ(IXVPR),XPVD(IXVPR),XTVQ
     +  (IXVTA),XTVD(IXVTA), IPVQ(IXVPR),IPPV1(IXVPR),IPPV2(IXVPR), ITVQ
     +  (IXVTA),ITTV1(IXVTA),ITTV2(IXVTA), AMCVV1(N),AMCVV2(N),GACVV1
     +  (N),GACVV2(N), BGXVV1(N),BGYVV1(N),BGZVV1(N), BGXVV2(N),BGYVV2
     +  (N),BGZVV2(N), NCHVV1(N),NCHVV2(N),IJCVV1(N),IJCVV2(N), (PQVVA1
     +  (N,JU),PQVVA2(N,JU),PQVVB1(N,JU), PQVVB2(N,JU),JU=1,4)
 
 
 
 
   20 CONTINUE
C--------------------------------------------------------
      RETURN
   10   CONTINUE
C                                     EVENT REJECTED
C                                     START A NEW ONE
        IREJVV=1
      RETURN
C
 1030 FORMAT(I10,4F12.7,6I5/10X,4F12.6/10X,6F12.6,4I5/8F15.5/8F15.5)
 1040 FORMAT (' VV IREJ ',I10/
     +' AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',5F12.5/
     +' AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ',5F12.5)
 1050 FORMAT(' VV',4(4E12.4/),2E12.4/2I5/4E12.4)
 1060 FORMAT(' VV',6I5/6E12.4/2E12.4)
 1070 FORMAT(' VV',5I5/2(4E12.4/),2E12.4)
 1080 FORMAT(' VV',7I5/2(4E12.4/),2E12.4)
 1090 FORMAT(' VV',4I5/6E12.4/2E12.4)
 1100 FORMAT(' KKEVT - IRVV13=',I5)
 1110 FORMAT(' KKEVT - IRVV11=',I5)
 1120 FORMAT(' KKEVT - IRVV12=',I5)
C
      END
*-- Author :
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE KKEVSS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C-------------------------    TREATMENT OF SEA-SEA CHAIN SYSTEMS
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
*KEEP,ABRSS.
C     INCLUDE (ABRSS)
      COMMON /ABRSS/ AMCSS1(INTMX),AMCSS2(INTMX), GACSS1(INTMX),GACSS2
     +(INTMX), BGXSS1(INTMX),BGYSS1(INTMX),BGZSS1(INTMX), BGXSS2(INTMX),
     +BGYSS2(INTMX),BGZSS2(INTMX), NCHSS1(INTMX),NCHSS2(INTMX), IJCSS1
     +(INTMX),IJCSS2(INTMX), PQSSA1(INTMX,4),PQSSA2(INTMX,4), PQSSB1
     +(INTMX,4),PQSSB2(INTMX,4)
*KEEP,TRAFOP.
      COMMON /TRAFOP/ GAMP,BGAMP,BETP
*KEEP,NUCIMP.
      COMMON /NUCIMP/ PRMOM(5,248),TAMOM(5,248), PRMFEP,PRMFEN,TAMFEP,
     +TAMFEN, PREFEP,PREFEN,TAEFEP,TAEFEN, PREPOT(210),TAEPOT(210),
     +PREBIN,TAEBIN,FERMOD,ETACOU
*KEEP,DROPPT.
      LOGICAL INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA, IPADIS,
     +ISHMAL,LPAULI
      COMMON /DROPPT/ INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA,
     +IPADIS,ISHMAL,LPAULI
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEEP,REJEC.
      COMMON /REJEC/ IRCO1,IRCO2,IRCO3,IRCO4,IRCO5, IRSS11,IRSS12,
     +IRSS13,IRSS14, IRSV11,IRSV12,IRSV13,IRSV14, IRVS11,IRVS12,IRVS13,
     +IRVS14, IRVV11,IRVV12,IRVV13,IRVV14
*KEEP,PROJK.
      COMMON /PROJK/ IPROJK
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
      COMMON/RPTSHM/RPROJ,RTARG,BIMPAC
*KEND.
C===================================================================

      COMMON /ABRJT/XJQ1(INTMX),XJAQ1(INTMX),XJQ2(INTMX),XJAQ2(INTMX),
     *        IJJQ1(INTMX),IJJAQ1(INTMX),IJJQ2(INTMX),IJJAQ2(INTMX),
     *        AMJCH1(INTMX),AMJCH2(INTMX),GAMJH1(INTMX),GAMJH2(INTMX),
     *        BGJH1(INTMX),BGJH2(INTMX),THEJH1(INTMX),THEJH2(INTMX),
     *        BGXJH1(INTMX),BGYJH1(INTMX),BGZJH1(INTMX),
     *        BGXJH2(INTMX),BGYJH2(INTMX),BGZJH2(INTMX),
     *  PJETA1(INTMX,4),PJETA2(INTMX,4),PJETB1(INTMX,4),PJETB2(INTMX,4)
     * ,JHKKPH(INTMX),JHKKTH(INTMX),JHKKEX(INTMX),JHKKE1(INTMX)
      COMMON /NUCJTN/NONUJ1,NONUJT,NONUS1,NONUST
      COMMON /XSVTHR/ XSTHR,XVTHR,XDTHR,XSSTHR
      COMMON /MINIJ/IMINIJ,NOMJE,NOMJER,NREJEV,NOMJT,NOMJTR
      COMMON /ABRSOF/XSQ1(INTMX),XSAQ1(INTMX),XSQ2(INTMX),XSAQ2(INTMX),
     *        IJSQ1(INTMX),IJSAQ1(INTMX),IJSQ2(INTMX),IJSAQ2(INTMX),
     *        AMCCH1(INTMX),AMCCH2(INTMX),GAMCH1(INTMX),GAMCH2(INTMX),
     *        BGCH1(INTMX),BGCH2(INTMX),THECH1(INTMX),THECH2(INTMX),
     *        BGXCH1(INTMX),BGYCH1(INTMX),BGZCH1(INTMX),
     *        BGXCH2(INTMX),BGYCH2(INTMX),BGZCH2(INTMX),
     *        NCH1(INTMX),NCH2(INTMX),IJCH1(INTMX),IJCH2(INTMX),
     *  PSOFA1(INTMX,4),PSOFA2(INTMX,4),PSOFB1(INTMX,4),PSOFB2(INTMX,4)
     * ,JHKKPZ(INTMX),JHKKTZ(INTMX),JHKKSX(INTMX),JHKKS1(INTMX)
      COMMON /ZSEA/ZSEAAV,ZSEASU,ANZSEA
C
      XSOTHR=XSTHR
      IF(NREJEV.GE.0)XSOTHR=0.
      IMINIJ=1
C
      IF(IPEV.GE.4)WRITE(6,*)' KKEVSS:NSS ',NSS 
C-----------------------------------------------------------------------
      DO 20 N=1,NSS
C---------------------------drop recombined chain pairs
        IF(IPEV.GE.4)WRITE(6,*)' KKEVSS:NCHSS1(N),NCHSS2(N)',
     *    NCHSS1(N),NCHSS2(N)
        IF(NCHSS1(N).EQ.99.AND.NCHSS2(N).EQ.99)GO TO 20
C
C
C***                 4-MOMENTA OF PROJECTILE SEA-QUARK PAIRS IN NN-CMS
        IXSPR=INTSS1(N)
        IXSPR1=INTSS1(N+1)
        IF (N.EQ.NSS)THEN
          IF (NSS.GE.2)THEN
            IXSPR1=INTSS1(N-1)
          ELSE
            IXSPR1=IXSPR
          ENDIF
        ENDIF
        INUCPR=IFROSP(IXSPR)
        JNUCPR=ITOVP(INUCPR)
C
C                        SUBTRACT HARD SCATTERED X VALUES FROM DIQUARKS
C
        IIFROP=IFROSP(IXSPR)
        IXVPR=ITOVP(IIFROP)
C
        IXSTA=INTSS2(N)
        IIFROT=IFROST(IXSTA)
        IXVTA=ITOVT(IIFROT)
C
        XMAX1=XPSQ(IXSPR)+XPSAQ(IXSPR)+XPVQ(IXVPR)+XPVD(IXVPR)
     *        -2.D0*XSOTHR-XVTHR-XDTHR            
        XMAX2=XTSQ(IXSTA)+XTSAQ(IXSTA)+XTVQ(IXVTA)+XTVD(IXVTA)
     *        -2.D0*XSOTHR-XVTHR-XDTHR            
C
        IF(IPEV.GE.1)WRITE(6,'(A,2I5,4F9.3/A,2I5,4F9.3/A,3F9.3)')
     *'IXSPR,IXVPR,XPSQ(IXSPR),XPSAQ(IXSPR),XPVQ(IXVPR),XPVD(IXVPR)'
     *,IXSPR,IXVPR,XPSQ(IXSPR),XPSAQ(IXSPR),XPVQ(IXVPR),XPVD(IXVPR),
     *'IXSTA,IXVTA,XTSQ(IXSTA),XTSAQ(IXSTA),XTVQ(IXVTA),XTVD(IXVTA)'
     *,IXSTA,IXVTA,XTSQ(IXSTA),XTSAQ(IXSTA),XTVQ(IXVTA),XTVD(IXVTA),
     *'XSOTHR,XVTHR,XDTHR'
     *,XSOTHR,XVTHR,XDTHR
        IF(IPEV.GE.1)WRITE(6,'(A,2I5,2F9.3)')' KKEVSS,bef xptfl:n,nss'
     *  ,N,NSS,XMAX1,XMAX2
        IF (IMINIJ.EQ.1)THEN
      	  CALL XPTFL(NHARD,NSEA,IREG,XMAX1,XMAX2)
C	  NZSEA=NZSEA+1
C	  ANZSEA=NZSEA
	  ANZSEA=ANZSEA+1.D0
	  ZSEASU=ZSEASU+NSEA
	  ZSEAAV=ZSEASU/ANZSEA
        ENDIF
        IF(IPEV.GE.1)WRITE(6,'(A,3I10)')' SS,xptfl:nhard,nsea,ireg '
     *  ,NHARD,NSEA,IREG
	IF(IREG.EQ.1)NHARD=0
	IF(IREG.EQ.1)NSEA=0
C
        NOMJE=NOMJE+NHARD
C
        IF (NHARD.GE.1.AND.IMINIJ.EQ.1)THEN
        DO 71 IXX=NONUJ1,NONUJT
          JHKKPH(IXX)=IXVPR
          JHKKEX(IXX)=0
          JHKKE1(IXX)=0
          IF (XPSQ(IXSPR)-XJQ1(IXX).GE.XSOTHR) THEN
            XPSQ(IXSPR)=XPSQ(IXSPR)-XJQ1(IXX)
            JHKKE1(IXX)=1
          ELSEIF (XPSAQ(IXSPR)-XJQ1(IXX).GE.XSOTHR) THEN
            XPSAQ(IXSPR)=XPSAQ(IXSPR)-XJQ1(IXX)
            JHKKE1(IXX)=2
          ELSEIF (XPSAQ(IXSPR)-XJQ1(IXX)/2..GE.XSOTHR.AND.
     *            XPSQ(IXSPR)-XJQ1(IXX)/2..GE.XSOTHR) THEN
            XPSQ(IXSPR)=XPSQ(IXSPR)-XJQ1(IXX)/2.
            XPSAQ(IXSPR)=XPSAQ(IXSPR)-XJQ1(IXX)/2.
            JHKKE1(IXX)=5
          ELSEIF (XPSQ(IXSPR1)-XJQ1(IXX).GE.XSOTHR) THEN
            XPSQ(IXSPR1)=XPSQ(IXSPR1)-XJQ1(IXX)
            JHKKE1(IXX)=6
          ELSEIF (XPSAQ(IXSPR1)-XJQ1(IXX).GE.XSOTHR) THEN
            XPSAQ(IXSPR1)=XPSAQ(IXSPR1)-XJQ1(IXX)
            JHKKE1(IXX)=7
          ELSEIF (XPSAQ(IXSPR1)-XJQ1(IXX)/2..GE.XSOTHR.AND.
     *            XPSQ(IXSPR1)-XJQ1(IXX)/2..GE.XSOTHR) THEN
            XPSQ(IXSPR1)=XPSQ(IXSPR1)-XJQ1(IXX)/2.
            XPSAQ(IXSPR1)=XPSAQ(IXSPR1)-XJQ1(IXX)/2.
            JHKKE1(IXX)=8
          ELSEIF (XPVQ(IXVPR)-XJQ1(IXX).GE.XVTHR) THEN
            XPVQ(IXVPR)=XPVQ(IXVPR)-XJQ1(IXX)
            JHKKE1(IXX)=3
          ELSEIF(XPVD(IXVPR)-XJQ1(IXX).GE.XDTHR) THEN
            XPVD(IXVPR)=XPVD(IXVPR)-XJQ1(IXX)
            JHKKE1(IXX)=4
          ENDIF
   71   CONTINUE
        ENDIF
C
        IXSTA=INTSS2(N)
        IXSTA1=INTSS2(N+1)
        IF (N.EQ.NSS)THEN
          IF (NSS.GE.2)THEN
            IXSTA1=INTSS2(N-1)
          ELSE
            IXSTA1=IXSTA
          ENDIF
        ENDIF
        INUCTA=IFROST(IXSTA)
        JNUCTA=ITOVT(INUCTA)
C
C
C                        SUBTRACT HARD SCATTERED X VALUES FROM DIQUARKS
C
        IIFROT=IFROST(IXSTA)
        IXVTA=ITOVT(IIFROT)
C
        IF (NHARD.GE.1.AND.IMINIJ.EQ.1) THEN
        DO 771 IXX=NONUJ1,NONUJT
          JHKKTH(IXX)=IXVTA
          IF(JHKKE1(IXX).EQ.0) THEN
            JHKKEX(IXX)=0
            GO TO 771
          ENDIF
          IF (XTSQ(IXSTA)-XJQ2(IXX).GE.XSOTHR) THEN
            XTSQ(IXSTA)=XTSQ(IXSTA)-XJQ2(IXX)
            JHKKEX(IXX)=1
            NOMJER=NOMJER+1
          ELSEIF (XTSAQ(IXSTA)-XJQ2(IXX).GE.XSOTHR) THEN
            XTSAQ(IXSTA)=XTSAQ(IXSTA)-XJQ2(IXX)
            JHKKEX(IXX)=1
            NOMJER=NOMJER+1
          ELSEIF (XTSAQ(IXSTA)-XJQ2(IXX)/2..GE.XSOTHR.AND.
     *            XTSQ(IXSTA)-XJQ2(IXX)/2..GE.XSOTHR) THEN
            XTSAQ(IXSTA)=XTSAQ(IXSTA)-XJQ2(IXX)/2.
            XTSQ(IXSTA)=XTSQ(IXSTA)-XJQ2(IXX)/2.
            JHKKEX(IXX)=1
            NOMJER=NOMJER+1
          ELSEIF (XTSQ(IXSTA1)-XJQ2(IXX).GE.XSOTHR) THEN
            XTSQ(IXSTA1)=XTSQ(IXSTA1)-XJQ2(IXX)
            JHKKEX(IXX)=1
            NOMJER=NOMJER+1
          ELSEIF (XTSAQ(IXSTA1)-XJQ2(IXX).GE.XSOTHR) THEN
            XTSAQ(IXSTA1)=XTSAQ(IXSTA1)-XJQ2(IXX)
            JHKKEX(IXX)=1
            NOMJER=NOMJER+1
          ELSEIF (XTSAQ(IXSTA1)-XJQ2(IXX)/2..GE.XSOTHR.AND.
     *            XTSQ(IXSTA1)-XJQ2(IXX)/2..GE.XSOTHR) THEN
            XTSAQ(IXSTA1)=XTSAQ(IXSTA1)-XJQ2(IXX)/2.
            XTSQ(IXSTA1)=XTSQ(IXSTA1)-XJQ2(IXX)/2.
            JHKKEX(IXX)=1
            NOMJER=NOMJER+1
          ELSEIF (XTVQ(IXVTA)-XJQ2(IXX).GE.XVTHR) THEN
            XTVQ(IXVTA)=XTVQ(IXVTA)-XJQ2(IXX)
            JHKKEX(IXX)=1
            NOMJER=NOMJER+1
          ELSEIF(XTVD(IXVTA)-XJQ2(IXX).GE.XDTHR) THEN
            XTVD(IXVTA)=XTVD(IXVTA)-XJQ2(IXX)
            JHKKEX(IXX)=1
            NOMJER=NOMJER+1
          ELSE
            JHKKEX(IXX)=0
            IF (JHKKE1(IXX).EQ.1) THEN
              XPSQ(IXSPR)=XPSQ(IXSPR)+XJQ1(IXX)
            ELSEIF (JHKKE1(IXX).EQ.2) THEN
              XPSAQ(IXSPR)=XPSAQ(IXSPR)+XJQ1(IXX)
            ELSEIF (JHKKE1(IXX).EQ.3) THEN
              XPVQ(IXVPR)=XPVQ(IXVPR)+XJQ1(IXX)
            ELSEIF (JHKKE1(IXX).EQ.4) THEN
              XPVD(IXVPR)=XPVD(IXVPR)+XJQ1(IXX)
            ELSEIF (JHKKE1(IXX).EQ.5) THEN
              XPSQ(IXSPR)=XPSQ(IXSPR)+XJQ1(IXX)/2.
              XPSAQ(IXSPR)=XPSAQ(IXSPR)+XJQ1(IXX)/2.
            ELSEIF (JHKKE1(IXX).EQ.6) THEN
              XPSQ(IXSPR1)=XPSQ(IXSPR1)+XJQ1(IXX)
            ELSEIF (JHKKE1(IXX).EQ.7) THEN
              XPSAQ(IXSPR1)=XPSAQ(IXSPR1)+XJQ1(IXX)
            ELSEIF (JHKKE1(IXX).EQ.8) THEN
              XPSQ(IXSPR1)=XPSQ(IXSPR1)+XJQ1(IXX)/2.
              XPSAQ(IXSPR1)=XPSAQ(IXSPR1)+XJQ1(IXX)/2.
            ENDIF
          ENDIF
  771   CONTINUE
        ENDIF
C
C                        SUBTRACT SECONDARY SEA X VALUES FROM DIQUARKS
C
        IF (NSEA.GE.1)THEN
        DO 271 IXX=NONUS1,NONUST
          JHKKPZ(IXX)=IXVPR
          JHKKSX(IXX)=0
          JHKKS1(IXX)=0
          IF (XPSQ(IXSPR)-XSQ1(IXX)-XSAQ1(IXX).GT.XVTHR)THEN
            XPSQ(IXSPR)=XPSQ(IXSPR)-XSQ1(IXX)-XSAQ1(IXX)
            JHKKS1(IXX)=3
          ELSEIF (XPSAQ(IXSPR)-XSQ1(IXX)-XSAQ1(IXX).GT.XVTHR)THEN
            XPSAQ(IXSPR)=XPSAQ(IXSPR)-XSQ1(IXX)-XSAQ1(IXX)
            JHKKS1(IXX)=4
          ELSEIF (XPVQ(IXVPR)-XSQ1(IXX)-XSAQ1(IXX).GT.XVTHR)THEN
            XPVQ(IXVPR)=XPVQ(IXVPR)-XSQ1(IXX)-XSAQ1(IXX)
            JHKKS1(IXX)=1
          ELSEIF (XPVD(IXVPR)-XSQ1(IXX)-XSAQ1(IXX).GT.XDTHR)THEN
            XPVD(IXVPR)=XPVD(IXVPR)-XSQ1(IXX)-XSAQ1(IXX)
            JHKKS1(IXX)=2
          ENDIF
  271   CONTINUE
        ENDIF
C
        INUCTA=IFROVT(IXVTA)
C
C                        SUBTRACT SECONDARY SEA X VALUES FROM DIQUARKS
C
        IF (NSEA.GE.1)THEN
        DO 2771 IXX=NONUS1,NONUST
          JHKKTZ(IXX)=IXVTA
          IF (JHKKS1(IXX).EQ.0)THEN
            JHKKSX(IXX)=0
            GO TO 2771
          ENDIF
          IF (XTSQ(IXSTA)-XSQ2(IXX)-XSAQ2(IXX).GT. XVTHR) THEN
            XTSQ(IXSTA)=XTSQ(IXSTA)-XSQ2(IXX)-XSAQ2(IXX)
            JHKKSX(IXX)=1
C           NOMJER=NOMJER+1
          ELSEIF (XTSAQ(IXSTA)-XSQ2(IXX)-XSAQ2(IXX).GT. XVTHR) THEN
            XTSAQ(IXSTA)=XTSAQ(IXSTA)-XSQ2(IXX)-XSAQ2(IXX)
            JHKKSX(IXX)=1
C           NOMJER=NOMJER+1
          ELSEIF (XTVQ(IXVTA)-XSQ2(IXX)-XSAQ2(IXX).GT. XVTHR) THEN
            XTVQ(IXVTA)=XTVQ(IXVTA)-XSQ2(IXX)-XSAQ2(IXX)
            JHKKSX(IXX)=1
C           NOMJER=NOMJER+1
          ELSEIF(XTVD(IXVTA)-XSQ2(IXX)-XSAQ2(IXX).GT.XDTHR)THEN
            XTVD(IXVTA)=XTVD(IXVTA)-XSQ2(IXX)-XSAQ2(IXX)
            JHKKSX(IXX)=1
C           NOMJER=NOMJER+1
          ELSE
            JHKKSX(IXX)=0
            IF (JHKKS1(IXX).EQ.1)THEN
              XPVQ(IXVPR)=XPVQ(IXVPR)+XSQ1(IXX)+XSAQ1(IXX)
            ELSEIF(JHKKS1(IXX).EQ.2)THEN
              XPVD(IXVPR)=XPVD(IXVPR)+XSQ1(IXX)+XSAQ1(IXX)
            ELSEIF(JHKKS1(IXX).EQ.3)THEN
              XPSQ(IXSPR)=XPSQ(IXSPR)+XSQ1(IXX)+XSAQ1(IXX)
            ELSEIF(JHKKS1(IXX).EQ.4)THEN
              XPSAQ(IXSPR)=XPSAQ(IXSPR)+XSQ1(IXX)+XSAQ1(IXX)
            ENDIF
          ENDIF
 2771   CONTINUE
        ENDIF
C
C
        XMAX1=XPSQ(IXSPR)+XPSAQ(IXSPR)+XPVQ(IXVPR)+XPVD(IXVPR)
     *        -2.D0*XSOTHR-XVTHR-XDTHR            
        XMAX2=XTSQ(IXSTA)+XTSAQ(IXSTA)+XTVQ(IXVTA)+XTVD(IXVTA)
     *        -2.D0*XSOTHR-XVTHR-XDTHR            
C
        IF(IPEV.GE.1)WRITE(6,'(A,2I5,4F9.3/A,2I5,4F9.3/A,3F9.3)')
     *'IXSPR,IXVPR,XPSQ(IXSPR),XPSAQ(IXSPR),XPVQ(IXVPR),XPVD(IXVPR)'
     *,IXSPR,IXVPR,XPSQ(IXSPR),XPSAQ(IXSPR),XPVQ(IXVPR),XPVD(IXVPR),
     *'IXSTA,IXVTA,XTSQ(IXSTA),XTSAQ(IXSTA),XTVQ(IXVTA),XTVD(IXVTA)'
     *,IXSTA,IXVTA,XTSQ(IXSTA),XTSAQ(IXSTA),XTVQ(IXVTA),XTVD(IXVTA),
     *'XSOTHR,XVTHR,XDTHR'
     *,XSOTHR,XVTHR,XDTHR
        IF(IPEV.GE.1)WRITE(6,'(A,2I5,2F9.3)')' KKEVSS,aft xptfl:n,nss'
     *  ,N,NSS,XMAX1,XMAX2
C===================================================================
C-----------------------------------------------------------------------
C     DO 20 N=1,NSS
C---------------------------drop recombined chain pairs
C       IF(NCHSS1(N).EQ.99.AND.NCHSS2(N).EQ.99)GO TO 20
C***                 4-MOMENTA OF PROJECTILE SEA-QUARK PAIRS IN NN-CMS
        IXSPR=INTSS1(N)
        INUCPR=IFROSP(IXSPR)
        JNUCPR=ITOVP(INUCPR)
C
        PSQPX=XPSQ(IXSPR)*PRMOM(1,INUCPR)
        PSQPY=XPSQ(IXSPR)*PRMOM(2,INUCPR)
        PSQPZ=XPSQ(IXSPR)*PRMOM(3,INUCPR)
        PSQE=XPSQ(IXSPR)*PRMOM(4,INUCPR)
        PSAQPX=XPSAQ(IXSPR)*PRMOM(1,INUCPR)
        PSAQPY=XPSAQ(IXSPR)*PRMOM(2,INUCPR)
        PSAQPZ=XPSAQ(IXSPR)*PRMOM(3,INUCPR)
        PSAQE=XPSAQ(IXSPR)*PRMOM(4,INUCPR)
C
C***                 4-MOMENTA OF TARGET SEA-QUARK PAIRS IN NN-CMS
        IXSTA=INTSS2(N)
        INUCTA=IFROST(IXSTA)
        JNUCTA=ITOVT(INUCTA)
C
        TSQPX=XTSQ(IXSTA)*TAMOM(1,INUCTA)
        TSQPY=XTSQ(IXSTA)*TAMOM(2,INUCTA)
        TSQPZ=XTSQ(IXSTA)*TAMOM(3,INUCTA)
        TSQE=XTSQ(IXSTA)*TAMOM(4,INUCTA)
        TSAQPX=XTSAQ(IXSTA)*TAMOM(1,INUCTA)
        TSAQPY=XTSAQ(IXSTA)*TAMOM(2,INUCTA)
        TSAQPZ=XTSAQ(IXSTA)*TAMOM(3,INUCTA)
        TSAQE=XTSAQ(IXSTA)*TAMOM(4,INUCTA)
C                                               j.r.6.5.93
C
C                     multiple scattering of sea quark chain ends
C
      IF(IT.GT.1)THEN
      ITNU=IP+INUCTA
      RTIX=(VHKK(1,ITNU))*1.E12-BIMPAC*0.1
      RTIY=VHKK(2,ITNU)*1.E12
      RTIZ=VHKK(3,ITNU)*1.E12
      CALL CROMSC(PSQPX,PSQPY,PSQPZ,PSQE,RTIX,RTIY,RTIZ,
     *            PSQNX,PSQNY,PSQNZ,PSQNE,5)
      PSQPX=PSQNX
      PSQPY=PSQNY
      PSQPZ=PSQNZ
      PSQE=PSQNE
      CALL CROMSC(PSAQPX,PSAQPY,PSAQPZ,PSAQE,RTIX,RTIY,RTIZ,
     *            PSAQNX,PSAQNY,PSAQNZ,PSAQNE,6)
      PSAQPX=PSAQNX
      PSAQPY=PSAQNY
      PSAQPZ=PSAQNZ
      PSAQE=PSAQNE
C                                                ---------
 
C                                               j.r.6.5.93
C
C                     multiple scattering of sea quark chain ends
C
      ITNU=IP+INUCTA
      RTIX=(VHKK(1,ITNU))*1.E12-BIMPAC*0.1
      RTIY=VHKK(2,ITNU)*1.E12
      RTIZ=VHKK(3,ITNU)*1.E12
      CALL CROMSC(TSQPX,TSQPY,TSQPZ,TSQE,RTIX,RTIY,RTIZ,
     *            TSQNX,TSQNY,TSQNZ,TSQNE,7)
      TSQPX=TSQNX
      TSQPY=TSQNY
      TSQPZ=TSQNZ
      TSQE=TSQNE
      CALL CROMSC(TSAQPX,TSAQPY,TSAQPZ,TSAQE,RTIX,RTIY,RTIZ,
     *            TSAQNX,TSAQNY,TSAQNZ,TSAQNE,8)
      TSAQPX=TSAQNX
      TSAQPY=TSAQNY
      TSAQPZ=TSAQNZ
      TSAQE=TSAQNE
      ENDIF 
C                                                j.r.10.5.93
       IF(IP.GE.1)GO TO 1779
        PSQPZ2=PSQE**2-PSQPX**2-PSQPY**2
        IF(PSQPZ2.GE.0.)THEN
          PSQPZ=SQRT(PSQPZ2)
        ELSE
          PSQPX=0.
          PSQPY=0.
          PSQPZ=PSQE
        ENDIF
C
        PAQPZ2=PSAQE**2-PSAQPX**2-PSAQPY**2
        IF(PAQPZ2.GE.0.)THEN
          PSAQPZ=SQRT(PAQPZ2)
        ELSE
          PSAQPX=0.
          PSAQPY=0.
          PSAQPZ=PSAQE
        ENDIF
C
        TSQPZ2=TSQE**2-TSQPX**2-TSQPY**2
        IF(TSQPZ2.GE.0.)THEN
          TSQPZ=-SQRT(TSQPZ2)
        ELSE
          TSQPX=0.
          TSQPY=0.
          TSQPZ=TSQE
        ENDIF
C
        TAQPZ2=TSAQE**2-TSAQPX**2-TSAQPY**2
        IF(TAQPZ2.GE.0.)THEN
          TSAQPZ=-SQRT(TAQPZ2)
        ELSE
          TSAQPX=0.
          TSAQPY=0.
          TSAQPZ=PSAQE
        ENDIF
 1779  CONTINUE
C                                                ---------
C
C                                      changej.r.6.5.93
        PTXSQ1=0.
        PTXSA1=0.
        PTXSQ2=0.
        PTXSA2=0.
        PTYSQ1=0.
        PTYSA1=0.
        PTYSQ2=0.
        PTYSA2=0.
        PTXSQ1=PSQPX
        PTXSA1=PSAQPX
        PTXSQ2=TSQPX
        PTXSA2=TSAQPX
        PTYSQ1=PSQPY
        PTYSA1=PSAQPY
        PTYSQ2=TSQPY
        PTYSA2=TSAQPY
        PLQ1=PSQPZ
        PLAQ1=PSAQPZ
        PLQ2=TSQPZ
        PLAQ2=TSAQPZ
        EQ1=PSQE
        EAQ1=PSAQE
        EQ2=TSQE
        EAQ2=TSAQE
C                                       ---------------
          IF(IPEV.GE.1) THEN
            WRITE(6,1060) IRSS13
            WRITE(6,1070)  PTXSQ1,
     +      PTYSQ1,PLQ1,EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,PTYSQ2,
     +      PLQ2,EQ2,PTXSA2,PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2,IREJ,IKVALA,
     +      PTTQ1,PTTA1
          ENDIF
        IKVALA=0
        NSELPT=0
	IF(IP.EQ.1)NSELPT=1
       IF(IOUXEV.GE.6)WRITE(6,'(A)')' KKEVSS call SELPT'
        IF(NSELPT.EQ.1)CALL SELPT( PTXSQ1,PTYSQ1,PLQ1,
     +  EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,
     +  PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1,
     *       PTTQ2,PTTA2,
     +  NSELPT)
        IF(NSELPT.EQ.0)CALL SELPT4( PTXSQ1,PTYSQ1,PLQ1,
     +  EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,
     +  PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1,
     +  NSELPT)
          IF(IPEV.GE.1) THEN
            WRITE(6,1060) IRSS13
            WRITE(6,1070)  PTXSQ1,
     +      PTYSQ1,PLQ1,EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,PTYSQ2,
     +      PLQ2,EQ2,PTXSA2,PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2,IREJ,IKVALA,
     +      PTTQ1,PTTA1
          ENDIF
        IF (IREJ.EQ.1) THEN
          IRSS13=IRSS13 + 1
          IF(IPEV.GE.1) THEN
            WRITE(6,1060) IRSS13
            WRITE(6,1070) PTXSQ1,
     +      PTYSQ1,PLQ1,EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,PTYSQ2,
     +      PLQ2,EQ2,PTXSA2,PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2,IREJ,IKVALA,
     +      PTTQ1,PTTA1
          ENDIF
                                                                GO TO 10
        ENDIF
C
C***  4-MOMENTA OF CHAINS IN THIS FRAME
C
        PTXCH1=PTXSQ1 + PTXSA2
        PTYCH1=PTYSQ1 + PTYSA2
        PTZCH1=PLQ1 + PLAQ2
        ECH1=EQ1 + EAQ2
        PTXCH2=PTXSQ2 + PTXSA1
        PTYCH2=PTYSQ2 + PTYSA1
        PTZCH2=PLQ2 + PLAQ1
        ECH2=EQ2 + EAQ1
        AMMM=SQRT((ECH1+ECH2)**2-(PTXCH1+PTXCH2)**2
     +            -(PTYCH1+PTYCH2)**2-(PTZCH1+PTZCH2)**2)
C
*
        IF (IPEV.GE.6) WRITE(6,1040) IREJ,
     +  AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1, AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2
 
C
C  REPLACE SMALL MASS CHAINS BY PSEUDOSCALAR OR VECTOR MESONS
C      FIRST FOR CHAIN 1  (PROSQ - TARASQ,  I.E. QUARK-AQUARK)
C
        CALL COMCMA(IPSQ(IXSPR),ITSAQ(IXSTA), IJNCH1,NNCH1,IREJ,AMCH1,
     +  AMCH1N)
C***                                       MASS BELOW PSEUDOSCALAR MASS
        IF(IREJ.EQ.1) THEN
          IRSS11=IRSS11 + 1
          IF(IPEV.GE.1) THEN
            WRITE(6,1080) IRSS11
            WRITE(6,1100) IPSQ(IXSPR),ITSAQ(IXSTA),IJNCH1,NNCH1,IREJ,
     +      XPSQ(IXSPR),XPSAQ(IXSPR),XPSQCM,XPSACM, XTSQ(IXSTA),XTSAQ
     +      (IXSTA),XTSQCM,XTSACM, AMCH1,AMCH1N
 
          ENDIF
                                                                 GOTO 10
        ENDIF
C                                 CORRECT KINEMATICS FOR CHAIN 1
C***                MOMENTUM CORRECTION FOR CHANGED MASS OF CHAIN 1
        IF(NNCH1.NE.0)THEN
             CALL CORMOM(AMCH1,AMCH2,AMCH1N,AMCH2N, 
     +  PTXSQ1,PTYSQ1,PLQ1,EQ1,
     +  PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,PTYSA2,
     +  PLAQ2,EAQ2, PTXCH1,PTYCH1,PTZCH1,ECH1, PTXCH2,PTYCH2,PTZCH2,
     +  ECH2,IREJ)
          AMCH2=AMCH2N
        ENDIF
        IF(IREJ.EQ.1)GO TO 10
C
        IF(IPEV.GE.6)WRITE(6,1050) IREJ,
     +  AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1, AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2
 
C
C  REPLACE SMALL MASS CHAINS BY PSEUDOSCALAR OR VECTOR MESONS
C            SECOND FOR CHAIN 2 (PROSAQ - TARSQ,  I.E. AQUARK-QUARK)
C
        CALL COMCMA(ITSQ(IXSTA),IPSAQ(IXSPR), IJNCH2,NNCH2,IREJ,AMCH2,
     +  AMCH2N)
c  rejection of both s-s chains if mass of chain 2 too low
        IF(IREJ.EQ.1) THEN
          IRSS12=IRSS12 + 1
          IF(IPEV.GE.1) THEN
            WRITE(6,1090) IRSS12
            WRITE(6,1100) IPSAQ(IXSPR),ITSQ(IXSTA),IJNCH2,NNCH2,IREJ,
     +      XPSQ(IXSPR),XPSAQ(IXSPR),XPSQCM,XPSACM, XTSQ(IXSTA),XTSAQ
     +      (IXSTA),XTSQCM,XTSACM, AMCH2,AMCH2N
 
          ENDIF
                                                                 GOTO 10
        ENDIF
C                                if AMCH2 changed in COBCMA/COMCMA
C                                CORVAL corrects chain kinematics
C                                according to 2-body kinem.
C                                with fixed masses
        IF(NNCH2.NE.0) THEN
          AMCH2=AMCH2N
         AMMM=SQRT((ECH1+ECH2)**2-(PTXCH1+PTXCH2)**2
     +            -(PTYCH1+PTYCH2)**2-(PTZCH1+PTZCH2)**2)
        EEE=ECH1+ECH2
        PXXX=PTXCH1+PTXCH2
        PYYY=PTYCH1+PTYCH2
        PZZZ=PTZCH1+PTZCH2
        GAMMM=EEE/(AMMM+1.E-4)
        BGGGX=PXXX/(AMMM+1.E-4)
        BGGGY=PYYY/(AMMM+1.E-4)
        BGGGZ=PZZZ/(AMMM+1.E-4)
C                        TRANSFORM BOTH CHAINS  INTO TWO CHAIN-CMS
C
C                                   4-MOMENTA OF CHAINS
        CALL DALTRA(GAMMM,-BGGGX,-BGGGY,-BGGGZ,
     +  PTXCH1,PTYCH1,PTZCH1,ECH1,
     +  PPPCH1, QTXCH1,QTYCH1,QTZCH1,QECH1)
C
        CALL DALTRA(GAMMM,-BGGGX,-BGGGY,-BGGGZ,
     +  PTXCH2,PTYCH2,PTZCH2,ECH2,
     +  PPPCH2, QTXCH2,QTYCH2,QTZCH2,QECH2)
C
          NORIG=22
          CALL CORVAL(AMMM,IREJ,AMCH1,AMCH2, QTXCH1,QTYCH1,QTZCH1,QECH1,
     +    QTXCH2,QTYCH2,QTZCH2,QECH2,NORIG)
C                        TRANSFORM BOTH CHAINS  INTO TWO CHAIN-CMS
C
C                                   4-MOMENTA OF CHAINS
 
        CALL DALTRA(GAMMM,BGGGX,BGGGY,BGGGZ, QTXCH1,QTYCH1,QTZCH1,QECH1,
     +  PPPCH1, PTXCH1,PTYCH1,PTZCH1,ECH1)
C
        CALL DALTRA(GAMMM,BGGGX,BGGGY,BGGGZ, QTXCH2,QTYCH2,QTZCH2,QECH2,
     +  PPPCH2, PTXCH2,PTYCH2,PTZCH2,ECH2)
C
C
          IF(IPEV.GE.6) THEN
            WRITE(6,'(A/3(1PE15.4),3I5)')
     +      ' SS - CALL CORVAL: AMMM, AMCH1, AMCH2, NNCH1, NNCH2, IREJ',
     +      AMMM, AMCH1, AMCH2, NNCH1, NNCH2, IREJ
            WRITE(6,1050) IREJ, AMCH1,
     +      PTXCH1,PTYCH1,PTZCH1,ECH1, AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2
 
          ENDIF
          IF(IREJ.EQ.1) THEN
*                           AMCH1N + AMCH2N > AMMM - 0.2
*                           reject event
            IRSS14=IRSS14+1
                                                                 GOTO 10
          ENDIF
        ENDIF
       QTXCH1=PTXCH1
       QTYCH1=PTYCH1
       QTZCH1=PTZCH1
       QECH1=ECH1
       QTXCH2=PTXCH2
       QTYCH2=PTYCH2
       QTZCH2=PTZCH2
       QECH2=ECH2
       PQSSA1(N,1)=PTXSQ1
       PQSSA1(N,2)=PTYSQ1
       PQSSA1(N,3)=PLQ1
       PQSSA1(N,4)=EQ1
       PQSSA2(N,1)=PTXSA2
       PQSSA2(N,2)=PTYSA2
       PQSSA2(N,3)=PLAQ2
       PQSSA2(N,4)=EAQ2
       PQSSB1(N,1)=PTXSQ2
       PQSSB1(N,2)=PTYSQ2
       PQSSB1(N,3)=PLQ2
       PQSSB1(N,4)=EQ2
       PQSSB2(N,1)=PTXSA1
       PQSSB2(N,2)=PTYSA1
       PQSSB2(N,3)=PLAQ1
       PQSSB2(N,4)=EAQ1
C-------------------
 
C
C                                      PUT S-S CHAIN ENDS INTO /HKKEVT/
C                                      MOMENTA IN NN-CMS
C                                      POSITION OF ORIGINAL NUCLEONS
C
        IHKKPD=JHKKPS(IXSPR)
        IHKKPO=IHKKPD -1
        IHKKTD=JHKKTS(IXSTA)
        IHKKTO=IHKKTD - 1
        IF (IPEV.GT.3)WRITE(6,1000)IXSPR,INUCPR,JNUCPR,IHKKPO,IHKKPD
 1000 FORMAT (' IXSPR,INUCPR,JNUCPR,IHKKPO,IHKKPD ',5I5)
        IF (IPEV.GT.3)WRITE(6,1010)IXSTA,INUCTA,JNUCTA,IHKKTO,IHKKTD
 1010 FORMAT (' IXSTA,INUCTA,JNUCTA,IHKKTO,IHKKTD ',5I5)
C                                     CHAIN 1 PROJECTILE SEA-QUARK
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=131
        IDHKK(IHKK)=IDHKK(IHKKPO)
        JMOHKK(1,IHKK)=IHKKPO
        JMOHKK(2,IHKK)=JMOHKK(2,IHKKPO)
        JDAHKK(1,IHKK)=IHKK+2
        JDAHKK(2,IHKK)=IHKK+2
        PHKK(1,IHKK)=PQSSA1(N,1)
        PHKK(2,IHKK)=PQSSA1(N,2)
        PHKK(3,IHKK)=PQSSA1(N,3)
        PHKK(4,IHKK)=PQSSA1(N,4)
        PHKK(5,IHKK)=0.
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
        VHKK(1,IHKK)=VHKK(1,IHKKPO)+XXPP
        VHKK(2,IHKK)=VHKK(2,IHKKPO)+YYPP
        VHKK(3,IHKK)=VHKK(3,IHKKPO)
        VHKK(4,IHKK)=VHKK(4,IHKKPO)
        IF (IPHKK.GE.2) WRITE(6,1020) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +  JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +  (PHKK(KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
 
 1020 FORMAT (I6,I4,5I6,9E10.2)
C                                     CHAIN 1 TARGET SEA-QUARK
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=132
        IDHKK(IHKK)=IDHKK(IHKKTD)
        JMOHKK(1,IHKK)=IHKKTD
        JMOHKK(2,IHKK)=JMOHKK(2,IHKKTD)
        JDAHKK(1,IHKK)=IHKK+1
        JDAHKK(2,IHKK)=IHKK+1
        PHKK(1,IHKK)=PQSSA2(N,1)
        PHKK(2,IHKK)=PQSSA2(N,2)
        PHKK(3,IHKK)=PQSSA2(N,3)
        PHKK(4,IHKK)=PQSSA2(N,4)
        PHKK(5,IHKK)=0.
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
        VHKK(1,IHKK)=VHKK(1,IHKKTD)+XXPP
        VHKK(2,IHKK)=VHKK(2,IHKKTD)+YYPP
        VHKK(3,IHKK)=VHKK(3,IHKKTD)
        VHKK(4,IHKK)=VHKK(4,IHKKTD)
        IF (IPHKK.GE.2) WRITE(6,1020) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +  JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +  (PHKK(KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
 
C
C                                     CHAIN 1 BEFORE FRAGMENTATION
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=6
        IDHKK(IHKK)=88888+NNCH1
        JMOHKK(1,IHKK)=IHKK-2
        JMOHKK(2,IHKK)=IHKK-1
        PHKK(1,IHKK)=QTXCH1
        PHKK(2,IHKK)=QTYCH1
        PHKK(3,IHKK)=QTZCH1
        PHKK(4,IHKK)=QECH1
        PHKK(5,IHKK)=AMCH1
C                           POSITION OF CREATED CHAIN IN LAB
C                             =POSITION OF TARGET NUCLEON
C                             TIME OF CHAIN CREATION IN LAB
C                             =TIME OF PASSAGE OF PROJECTILE
C                              NUCLEUS AT POSITION OF TAR. NUCLEUS
        VHKK(1,NHKK)=                VHKK(1,NHKK-1)
        VHKK(2,NHKK)=                VHKK(2,NHKK-1)
        VHKK(3,NHKK)=                VHKK(3,NHKK-1)
        VHKK(4,NHKK)=VHKK(3,NHKK)/BETP-VHKK(3,NHKK-2)/BGAMP
        MHKKSS(N)=IHKK
        IF (IPROJK.EQ.1)THEN
          WHKK(1,NHKK)=                VHKK(1,NHKK-2)
          WHKK(2,NHKK)=                VHKK(2,NHKK-2)
          WHKK(3,NHKK)=                VHKK(3,NHKK-2)
          WHKK(4,NHKK)=VHKK(4,NHKK)*GAMP-VHKK(3,NHKK)*BGAMP
          IF (IPHKK.GE.2) WRITE(6,1020) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +    JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +    (PHKK(KHKK,IHKK),KHKK=1,5), (WHKK(KHKK,IHKK),KHKK=1,4)
 
        ENDIF
        IF (IPHKK.GE.2) WRITE(6,1020) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +  JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +  (PHKK(KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
 
C
C
C                                     CHAIN 2 PROJECTILE SEA-QUARK
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=131
        IDHKK(IHKK)=IDHKK(IHKKPD)
        JMOHKK(1,IHKK)=IHKKPD
        JMOHKK(2,IHKK)=JMOHKK(2,IHKKPD)
        JDAHKK(1,IHKK)=IHKK+2
        JDAHKK(2,IHKK)=IHKK+2
        PHKK(1,IHKK)=PQSSB1(N,1)
        PHKK(2,IHKK)=PQSSB1(N,2)
        PHKK(3,IHKK)=PQSSB1(N,3)
        PHKK(4,IHKK)=PQSSB1(N,4)
        PHKK(5,IHKK)=0.
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
        VHKK(1,IHKK)=VHKK(1,IHKKPD)+XXPP
        VHKK(2,IHKK)=VHKK(2,IHKKPD)+YYPP
        VHKK(3,IHKK)=VHKK(3,IHKKPD)
        VHKK(4,IHKK)=VHKK(4,IHKKPD)
        IF (IPHKK.GE.2) WRITE(6,1020) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +  JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +  (PHKK(KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
 
C                                     CHAIN 2 TARGET SEA-QUARK
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=132
        IDHKK(IHKK)=IDHKK(IHKKTO)
        JMOHKK(1,IHKK)=IHKKTO
        JMOHKK(2,IHKK)=JMOHKK(2,IHKKTO)
        JDAHKK(1,IHKK)=IHKK+1
        JDAHKK(2,IHKK)=IHKK+1
        PHKK(1,IHKK)=PQSSB2(N,1)
        PHKK(2,IHKK)=PQSSB2(N,2)
        PHKK(3,IHKK)=PQSSB2(N,3)
        PHKK(4,IHKK)=PQSSB2(N,4)
        PHKK(5,IHKK)=0.
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
        VHKK(1,IHKK)=VHKK(1,IHKKTO)+XXPP
        VHKK(2,IHKK)=VHKK(2,IHKKTO)+YYPP
        VHKK(3,IHKK)=VHKK(3,IHKKTO)
        VHKK(4,IHKK)=VHKK(4,IHKKTO)
        IF (IPHKK.GE.2) WRITE(6,1020) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +  JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +  (PHKK(KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
 
C
C                                     CHAIN 2 BEFORE FRAGMENTATION
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=6
        IDHKK(IHKK)=88888+NNCH2
        JMOHKK(1,IHKK)=IHKK-2
        JMOHKK(2,IHKK)=IHKK-1
        PHKK(1,IHKK)=QTXCH2
        PHKK(2,IHKK)=QTYCH2
        PHKK(3,IHKK)=QTZCH2
        PHKK(4,IHKK)=QECH2
        PHKK(5,IHKK)=AMCH2
C                            POSITION OF CREATED CHAIN IN LAB
C                             =POSITION OF TARGET NUCLEON
C                             TIME OF CHAIN CREATION IN LAB
C                             =TIME OF PASSAGE OF PROJECTILE
C                              NUCLEUS AT POSITION OF TAR. NUCLEUS
        VHKK(1,NHKK)=                VHKK(1,NHKK-1)
        VHKK(2,NHKK)=                VHKK(2,NHKK-1)
        VHKK(3,NHKK)=                VHKK(3,NHKK-1)
        VHKK(4,NHKK)=VHKK(3,NHKK)/BETP-VHKK(3,NHKK-2)/BGAMP
        MHKKSS(N)=IHKK
        IF (IPROJK.EQ.1)THEN
          WHKK(1,NHKK)=                VHKK(1,NHKK-2)
          WHKK(2,NHKK)=                VHKK(2,NHKK-2)
          WHKK(3,NHKK)=                VHKK(3,NHKK-2)
          WHKK(4,NHKK)=VHKK(4,NHKK)*GAMP-VHKK(3,NHKK)*BGAMP
          IF (IPHKK.GE.2) WRITE(6,1020) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +    JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +    (PHKK(KHKK,IHKK),KHKK=1,5), (WHKK(KHKK,IHKK),KHKK=1,4)
 
        ENDIF
        IF (IPHKK.GE.2) WRITE(6,1020) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +  JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +  (PHKK(KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
 
C
C  NOW WE HAVE AN ACCEPTABLE SEA--SEA  EVENT
C  AND PUT IT INTO THE HISTOGRAM
C
        AMCSS1(N)=AMCH1
        AMCSS2(N)=AMCH2
        GACSS1(N)=QECH1/AMCH1
        BGXSS1(N)=QTXCH1/AMCH1
        BGYSS1(N)=QTYCH1/AMCH1
        BGZSS1(N)=QTZCH1/AMCH1
        GACSS2(N)=QECH2/AMCH2
        BGXSS2(N)=QTXCH2/AMCH2
        BGYSS2(N)=QTYCH2/AMCH2
        BGZSS2(N)=QTZCH2/AMCH2
        NCHSS1(N)=NNCH1
        NCHSS2(N)=NNCH2
        IJCSS1(N)=IJNCH1
        IJCSS2(N)=IJNCH2
        IF (IPEV.GE.6)WRITE(6,1030)N, XPSQ(IXSPR),XPSAQ(IXSPR),XTSQ
     +  (IXSTA),XTSAQ(IXSTA), IPSQ(IXSPR),IPSAQ(IXSPR),ITSQ(IXSTA),ITSAQ
     +  (IXSTA), ITSAQ(IXSTA), AMCSS1(N),AMCSS2(N),GACSS1(N),GACSS2(N),
     +  BGXSS1(N),BGYSS1(N),BGZSS1(N), BGXSS2(N),BGYSS2(N),BGZSS2(N),
     +  NCHSS1(N),NCHSS2(N),IJCSS1(N),IJCSS2(N), (PQSSA1(N,JU),PQSSA2
     +  (N,JU),PQSSB1(N,JU), PQSSB2(N,JU),JU=1,4)
 
 
 
 
                                                                GO TO 20
C***                     TREATMENT OF REJECTED SEA-SEA INTERACTIONS
   10   CONTINUE
        INLOSS(N)=.FALSE.
        XPVD(JNUCPR)=XPVD(JNUCPR) + XPSQ(IXSPR) + XPSAQ(IXSPR)
        XTVD(JNUCTA)=XTVD(JNUCTA) + XTSAQ(IXSTA) + XTSQ(IXSTA)
   20 CONTINUE
C
 1030 FORMAT(' SS - 104', I10,4F12.7,5I5/10X,4F12.6/10X,6F12.6,4I5/8F15.
     +5/8F15.5)
 1040 FORMAT (' SS: IREJ ',I10/
     +'     AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',5F12.5/
     +'     AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ',5F12.5)
 1050 FORMAT (' SS: IREJ  ',I10/
     +'     AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',5F12.5/
     +'     AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ',5F12.5)
 1060 FORMAT(' KKEVSS - IRSS13=',I5)
 1070 FORMAT( ' SS - 8002',4(4E12.4/),2E12.4/2I5/4E12.4)
 1080 FORMAT(' KKEVSS - IRSS11=',I5)
 1090 FORMAT(' KKEVSS - IRSS12=',I5)
 1100 FORMAT(' SS - 8006', 5I5/2(4E12.4/),2E12.4)
      RETURN
      END
*-- Author :
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE KKEVVS(IREJVS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C                             TREATMENT OF VALENCE-SEA CHAIN SYSTEMS
C
C---------------------------------------------------------------------
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
*KEEP,TRAFOP.
      COMMON /TRAFOP/ GAMP,BGAMP,BETP
*KEEP,NUCC.
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
*KEEP,NUCIMP.
      COMMON /NUCIMP/ PRMOM(5,248),TAMOM(5,248), PRMFEP,PRMFEN,TAMFEP,
     +TAMFEN, PREFEP,PREFEN,TAEFEP,TAEFEN, PREPOT(210),TAEPOT(210),
     +PREBIN,TAEBIN,FERMOD,ETACOU
*KEEP,ABRVS.
      COMMON /ABRVS/ AMCVS1(248),AMCVS2(248),GACVS1(248),GACVS2(248),
     +BGXVS1(248),BGYVS1(248),BGZVS1(248), BGXVS2(248),BGYVS2(248),
     +BGZVS2(248), NCHVS1(248),NCHVS2(248),IJCVS1(248),IJCVS2(248),
     +PQVSA1(248,4),PQVSA2(248,4), PQVSB1(248,4),PQVSB2(248,4)
*KEEP,DROPPT.
      LOGICAL INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA, IPADIS,
     +ISHMAL,LPAULI
      COMMON /DROPPT/ INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA,
     +IPADIS,ISHMAL,LPAULI
*KEEP,NUCPOS.
      COMMON /NUCPOS/INVVP(248),INVVT(248),INVSP(248),INVST(248), NUVV,
     +NUVS,NUSV,NUSS,INSVP(248),INSVT(248),INSSP(248),INSST(248), ISVEAP
     +(248),ISVEAT(248),ISSEAP(248),ISSEAT(248), IVSEAP(248),IVSEAT
     +(248), ISLOSP(248),ISLOST(248),INOOP(248),INOOT(248),NUOO
*KEEP,TAUFO.
      COMMON /TAUFO/  TAUFOR,KTAUGE,ITAUVE,INCMOD
*KEEP,RTAR.
      COMMON /RTAR/ RTARNU
*KEEP,INNU.
      COMMON /INNU/INUDEC
*KEEP,DINPDA.
      COMMON /DINPDA/ IMPS(6,6),IMVE(6,6),IB08(6,21),IB10(6,21), IA08
     +(6,21),IA10(6,21), A1,B1,B2,B3,LT,LE,BET,AS,B8,AME,DIQ,ISU
*KEEP,FERMI.
      COMMON /FERMI/ PQUAR(4,248),PAQUAR(4,248), TQUAR(4,248),TAQUAR
     +(4,248)
*KEEP,KETMAS.
      COMMON /KETMAS/ AM8(2),AM10(2),IB88(2),IB1010(2),AMCH1N,AMCH2N
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
*KEEP,REJEC.
      COMMON /REJEC/ IRCO1,IRCO2,IRCO3,IRCO4,IRCO5, IRSS11,IRSS12,
     +IRSS13,IRSS14, IRSV11,IRSV12,IRSV13,IRSV14, IRVS11,IRVS12,IRVS13,
     +IRVS14, IRVV11,IRVV12,IRVV13,IRVV14
*KEEP,PROJK.
      COMMON /PROJK/ IPROJK
      COMMON/RPTSHM/RPROJ,RTARG,BIMPAC
 
*KEND.
C===============================================================

      COMMON /ABRJT/XJQ1(INTMX),XJAQ1(INTMX),XJQ2(INTMX),XJAQ2(INTMX),
     *        IJJQ1(INTMX),IJJAQ1(INTMX),IJJQ2(INTMX),IJJAQ2(INTMX),
     *        AMJCH1(INTMX),AMJCH2(INTMX),GAMJH1(INTMX),GAMJH2(INTMX),
     *        BGJH1(INTMX),BGJH2(INTMX),THEJH1(INTMX),THEJH2(INTMX),
     *        BGXJH1(INTMX),BGYJH1(INTMX),BGZJH1(INTMX),
     *        BGXJH2(INTMX),BGYJH2(INTMX),BGZJH2(INTMX),
     *  PJETA1(INTMX,4),PJETA2(INTMX,4),PJETB1(INTMX,4),PJETB2(INTMX,4)
     * ,JHKKPH(INTMX),JHKKTH(INTMX),JHKKEX(INTMX),JHKKE1(INTMX)
      COMMON /NUCJTN/NONUJ1,NONUJT,NONUS1,NONUST
      COMMON /XSVTHR/ XSTHR,XVTHR,XDTHR,XSSTHR
      COMMON /MINIJ/IMINIJ,NOMJE,NOMJER,NREJEV,NOMJT,NOMJTR
      COMMON /ABRSOF/XSQ1(INTMX),XSAQ1(INTMX),XSQ2(INTMX),XSAQ2(INTMX),
     *        IJSQ1(INTMX),IJSAQ1(INTMX),IJSQ2(INTMX),IJSAQ2(INTMX),
     *        AMCCH1(INTMX),AMCCH2(INTMX),GAMCH1(INTMX),GAMCH2(INTMX),
     *        BGCH1(INTMX),BGCH2(INTMX),THECH1(INTMX),THECH2(INTMX),
     *        BGXCH1(INTMX),BGYCH1(INTMX),BGZCH1(INTMX),
     *        BGXCH2(INTMX),BGYCH2(INTMX),BGZCH2(INTMX),
     *        NCH1(INTMX),NCH2(INTMX),IJCH1(INTMX),IJCH2(INTMX),
     *  PSOFA1(INTMX,4),PSOFA2(INTMX,4),PSOFB1(INTMX,4),PSOFB2(INTMX,4)
     * ,JHKKPZ(INTMX),JHKKTZ(INTMX),JHKKSX(INTMX),JHKKS1(INTMX)
      COMMON /ZSEA/ZSEAAV,ZSEASU,ANZSEA
      IREJVS=0
      IMINIJ=1
C
      DO 10 N=1,NVS
C---------------------------drop recombined chain pairs
        IF(NCHVS1(N).EQ.99.AND.NCHVS2(N).EQ.99)GO TO 10
C
C-----------------------------------------------------------
        IXVPR=INTVS1(N)
        INUCPR=IFROVP(IXVPR)
        JNUCPR=ITOVP(INUCPR)
C---
        IXSTA=INTVS2(N)
        INUCTA=IFROST(IXSTA)
        JNUCTA=ITOVT(INUCTA)
        IIFROT=IFROST(IXSTA)
        IXVTA=ITOVT(IIFROT)
C
        XMAX1=XPVQ(IXVPR)+XPVD(IXVPR)-XVTHR-XDTHR
        XMAX2=XTSQ(IXSTA)+XTSAQ(IXSTA)+XTVQ(IXVTA)+XTVD(IXVTA)
     *        -2.D0*XSTHR-XVTHR-XDTHR
C
        IF(IPEV.GE.1)WRITE(6,'(A,2I5,2F9.3)')' KKEVVS,bef xptfl:n,nvs'
     *  ,N,NVS,XMAX1,XMAX2
        IF (IMINIJ.EQ.1)THEN
      	  CALL XPTFL(NHARD,NSEA,IREG,XMAX1,XMAX2)
C	  NZSEA=NZSEA+1
C	  ANZSEA=NZSEA
	  ANZSEA=ANZSEA+1.D0
	  ZSEASU=ZSEASU+NSEA
	  ZSEAAV=ZSEASU/ANZSEA
        ENDIF
        IF(IPEV.GE.1)WRITE(6,'(A,3I10)')' VS,xptfl:nhard,nsea,ireg '
     *  ,NHARD,NSEA,IREG
	IF(IREG.EQ.1)NHARD=0
	IF(IREG.EQ.1)NSEA=0
        NOMJE=NOMJE+NHARD
C
C
C                        SUBTRACT HARD SCATTERED X VALUES FROM DIQUARKS
C
        IF (NHARD.GE.1.AND.IMINIJ.EQ.1)THEN
        DO 71 IXX=NONUJ1,NONUJT
          JHKKPH(IXX)=IXVPR
          JHKKEX(IXX)=0
          JHKKE1(IXX)=0
          IF (XPVQ(IXVPR)-XJQ1(IXX).GE.XVTHR) THEN
            XPVQ(IXVPR)=XPVQ(IXVPR)-XJQ1(IXX)
            JHKKE1(IXX)=1
          ELSEIF(XPVD(IXVPR)-XJQ1(IXX).GE.XDTHR) THEN
            XPVD(IXVPR)=XPVD(IXVPR)-XJQ1(IXX)
            JHKKE1(IXX)=2
          ENDIF
   71   CONTINUE
        ENDIF
C---
        IXSTA=INTVS2(N)
        INUCTA=IFROST(IXSTA)
        JNUCTA=ITOVT(INUCTA)
C
C                        SUBTRACT HARD SCATTERED X VALUES FROM DIQUARKS
C
        IIFROT=IFROST(IXSTA)
        IXVTA=ITOVT(IIFROT)
C
        IF (NHARD.GE.1.AND.IMINIJ.EQ.1) THEN
        DO 771 IXX=NONUJ1,NONUJT
          JHKKTH(IXX)=IXVTA
          IF (JHKKE1(IXX).EQ.0)THEN
            JHKKEX(IXX)=0
            GOTO 771
          ENDIF
          IF (XTSQ(IXSTA)-XJQ2(IXX).GE.XSTHR) THEN
            XTSQ(IXSTA)=XTSQ(IXSTA)-XJQ2(IXX)
            JHKKEX(IXX)=1
            NOMJER=NOMJER+1
          ELSEIF (XTSAQ(IXSTA)-XJQ2(IXX).GE.XSTHR) THEN
            XTSAQ(IXSTA)=XTSAQ(IXSTA)-XJQ2(IXX)
            JHKKEX(IXX)=1
            NOMJER=NOMJER+1
          ELSEIF (XTVQ(IXVTA)-XJQ2(IXX).GE.XVTHR) THEN
            XTVQ(IXVTA)=XTVQ(IXVTA)-XJQ2(IXX)
            JHKKEX(IXX)=1
            NOMJER=NOMJER+1
          ELSEIF (XTVD(IXVTA)-XJQ2(IXX).GE.XDTHR)THEN
            XTVD(IXVTA)=XTVD(IXVTA)-XJQ2(IXX)
            JHKKEX(IXX)=1
            NOMJER=NOMJER+1
          ELSE
            JHKKEX(IXX)=0
            IF (JHKKE1(IXX).EQ.1)THEN
              XPVQ(IXVPR)=XPVQ(IXVPR)+XJQ1(IXX)
            ELSEIF (JHKKE1(IXX).EQ.2)THEN
              XPVD(IXVPR)=XPVD(IXVPR)+XJQ1(IXX)
            ENDIF
          ENDIF
  771   CONTINUE
        ENDIF
C
C
C                        SUBTRACT SECONDARY SEA X VALUES FROM DIQUARKS
C
        IF (NSEA.GE.1)THEN
        DO 271 IXX=NONUS1,NONUST
          JHKKPZ(IXX)=IXVPR
          JHKKSX(IXX)=0
          JHKKS1(IXX)=0
          IF (XPVQ(IXVPR)-XSQ1(IXX)-XSAQ1(IXX).GT.XVTHR)THEN
            XPVQ(IXVPR)=XPVQ(IXVPR)-XSQ1(IXX)-XSAQ1(IXX)
            JHKKS1(IXX)=1
          ELSEIF (XPVD(IXVPR)-XSQ1(IXX)-XSAQ1(IXX).GT.XDTHR)THEN
            XPVD(IXVPR)=XPVD(IXVPR)-XSQ1(IXX)-XSAQ1(IXX)
            JHKKS1(IXX)=2
          ENDIF
  271   CONTINUE
        ENDIF
C
        INUCTA=IFROVT(IXVTA)
C
C                        SUBTRACT SECONDARY SEA X VALUES FROM DIQUARKS
C
        IF (NSEA.GE.1)THEN
        DO 2771 IXX=NONUS1,NONUST
          JHKKTZ(IXX)=IXVTA
          IF (JHKKS1(IXX).EQ.0)THEN
            JHKKSX(IXX)=0
            GO TO 2771
          ENDIF
          IF (XTSQ(IXSTA)-XSQ2(IXX)-XSAQ2(IXX).GT. XVTHR) THEN
            XTSQ(IXSTA)=XTSQ(IXSTA)-XSQ2(IXX)-XSAQ2(IXX)
            JHKKSX(IXX)=1
C           NOMJER=NOMJER+1
          ELSEIF (XTSAQ(IXSTA)-XSQ2(IXX)-XSAQ2(IXX).GT. XVTHR) THEN
            XTSAQ(IXSTA)=XTSAQ(IXSTA)-XSQ2(IXX)-XSAQ2(IXX)
            JHKKSX(IXX)=1
C           NOMJER=NOMJER+1
          ELSEIF (XTVQ(IXVTA)-XSQ2(IXX)-XSAQ2(IXX).GT. XVTHR) THEN
            XTVQ(IXVTA)=XTVQ(IXVTA)-XSQ2(IXX)-XSAQ2(IXX)
            JHKKSX(IXX)=1
C           NOMJER=NOMJER+1
          ELSEIF(XTVD(IXVTA)-XSQ2(IXX)-XSAQ2(IXX).GT.XDTHR)THEN
            XTVD(IXVTA)=XTVD(IXVTA)-XSQ2(IXX)-XSAQ2(IXX)
            JHKKSX(IXX)=1
C           NOMJER=NOMJER+1
          ELSE
            JHKKSX(IXX)=0
            IF (JHKKS1(IXX).EQ.1)THEN
              XPVQ(IXVPR)=XPVQ(IXVPR)+XSQ1(IXX)+XSAQ1(IXX)
            ELSEIF(JHKKS1(IXX).EQ.2)THEN
              XPVD(IXVPR)=XPVD(IXVPR)+XSQ1(IXX)+XSAQ1(IXX)
            ENDIF
          ENDIF
 2771   CONTINUE
        ENDIF
C
C
        XMAX1=XPVQ(IXVPR)+XPVD(IXVPR)-XVTHR-XDTHR
        XMAX2=XTSQ(IXSTA)+XTSAQ(IXSTA)+XTVQ(IXVTA)+XTVD(IXVTA)
     *        -2.D0*XSTHR-XVTHR-XDTHR
C
        IF(IPEV.GE.1)WRITE(6,'(A,2I5,2F9.3)')' KKEVVS,aft xptfl:n,nvs'
     *  ,N,NVS,XMAX1,XMAX2
        GO TO 1302
 2302   CONTINUE
C
C               TRY TO INCREASE THE X-FRACTION OF TARGET SEA QUARK
C               ANTIQUARK PAIR BY XSTHR  XTSQ(IXSTA) XTSAQ(IXSTA)
C               DECREASING TARGET DIQUARK XTVD(IXVTA)
C
C       IF(XTSUT(IXVTA).EQ.0..AND.
C    *     XTVD(IXVTA)-2.*XSTHR.GE.XDTHR) THEN
C         XTSQ(IXSTA)=XTSQ(IXSTA)+XSTHR
C         XTSAQ(IXSTA)=XTSAQ(IXSTA)+XSTHR
C         XTVD(IXVTA)=XTVD(IXVTA)-2.*XSTHR
C         IREJ=0
C       ELSE
C         GO TO 302
C         GO TO 10
C       ENDIF
 1302   CONTINUE
C===============================================================
C***          4-MOMENTA OF PROJECTILE QUARKS AND DIQUARK-PAIRS IN NN-CMS
        IXVPR=INTVS1(N)
        INUCPR=IFROVP(IXVPR)
        JNUCPR=ITOVP(INUCPR)
C
        PVQPX=XPVQ(IXVPR)*PRMOM(1,INUCPR)
        PVQPY=XPVQ(IXVPR)*PRMOM(2,INUCPR)
        PVQPZ=XPVQ(IXVPR)*PRMOM(3,INUCPR)
        PVQE=XPVQ(IXVPR)*PRMOM(4,INUCPR)
        PVDQPX=XPVD(IXVPR)*PRMOM(1,INUCPR)
        PVDQPY=XPVD(IXVPR)*PRMOM(2,INUCPR)
        PVDQPZ=XPVD(IXVPR)*PRMOM(3,INUCPR)
        PVDQE=XPVD(IXVPR)*PRMOM(4,INUCPR)
C
        IF(IPEV.GE.7) THEN
          WRITE(6,1000) PVQPX,PVQPY,PVQPZ,PVQE, PVDQPX,PVDQPY,PVDQPZ,
     +    PVDQE
 1000 FORMAT(' VS:  PVQPX,PVQPY,PVQPZ,PVQE',
     +' PVDQPX,PVDQPY,PVDQPZ,PVDQE',/4E15.5/15X,4E15.5)
        ENDIF
C
C***                 4-MOMENTA OF TARGET SEA-QUARK PAIRS IN NN-CMS
        IXSTA=INTVS2(N)
        INUCTA=IFROST(IXSTA)
        JNUCTA=ITOVT(INUCTA)
C
        TSQPX=XTSQ(IXSTA)*TAMOM(1,INUCTA)
        TSQPY=XTSQ(IXSTA)*TAMOM(2,INUCTA)
        TSQPZ=XTSQ(IXSTA)*TAMOM(3,INUCTA)
        TSQE=XTSQ(IXSTA)*TAMOM(4,INUCTA)
        TSAQPX=XTSAQ(IXSTA)*TAMOM(1,INUCTA)
        TSAQPY=XTSAQ(IXSTA)*TAMOM(2,INUCTA)
        TSAQPZ=XTSAQ(IXSTA)*TAMOM(3,INUCTA)
        TSAQE=XTSAQ(IXSTA)*TAMOM(4,INUCTA)
C                                               j.r.6.5.93
C
C                     multiple scattering of valence quark chain ends
C
      IF(IT.GT.1)THEN
      ITNU=IP+INUCTA
      RTIX=(VHKK(1,ITNU))*1.E12-BIMPAC*0.1
      RTIY=VHKK(2,ITNU)*1.E12
      RTIZ=VHKK(3,ITNU)*1.E12
      CALL CROMSC(PVQPX,PVQPY,PVQPZ,PVQE,RTIX,RTIY,RTIZ,
     *            PVQNX,PVQNY,PVQNZ,PVQNE,9)
      PVQPX=PVQNX
      PVQPY=PVQNY
      PVQPZ=PVQNZ
      PVQE=PVQNE
      CALL CROMSC(PVDQPX,PVDQPY,PVDQPZ,PVDQE,RTIX,RTIY,RTIZ,
     *            PVDQNX,PVDQNY,PVDQNZ,PVDQNE,10)
      PVDQPX=PVDQNX
      PVDQPY=PVDQNY
      PVDQPZ=PVDQNZ
      PVDQE=PVDQNE
C                                                ---------
C
        IF(IPEV.GE.7) THEN
          WRITE(6,1010) N,NVS,IXVPR,INUCPR,INUCPR,IXSTA,INUCTA,JNUCTA
 1010 FORMAT(' VS: N,NVS,IXVPR,INUCPR,INUCPR,IXSTA,INUCTA,JNUCTA'/ 8I5)
 
          WRITE(6,1020) TSQPX,TSQPY,TSQPZ,TSQE, TSAQPX,TSAQPY,TSAQPZ,
     +    TSAQE
 1020 FORMAT(' VS:  TSQPX,TSQPY,TSQPZ,TSQE',
     +' TSAQPX,TSAQPY,TSAQPZ,TSAQE',/4E15.5/15X,4E15.5)
        ENDIF
C                                               j.r.6.5.93
C
C                     multiple scattering of sea quark chain ends
C
      ITNU=IP+INUCTA
      RTIX=(VHKK(1,ITNU))*1.E12-BIMPAC*0.1
      RTIY=VHKK(2,ITNU)*1.E12
      RTIZ=VHKK(3,ITNU)*1.E12
      CALL CROMSC(TSQPX,TSQPY,TSQPZ,TSQE,RTIX,RTIY,RTIZ,
     *            TSQNX,TSQNY,TSQNZ,TSQNE,11)
      TSQPX=TSQNX
      TSQPY=TSQNY
      TSQPZ=TSQNZ
      TSQE=TSQNE
      CALL CROMSC(TSAQPX,TSAQPY,TSAQPZ,TSAQE,RTIX,RTIY,RTIZ,
     *            TSAQNX,TSAQNY,TSAQNZ,TSAQNE,12)
      TSAQPX=TSAQNX
      TSAQPY=TSAQNY
      TSAQPZ=TSAQNZ
      TSAQE=TSAQNE
      ENDIF  
C                                                ---------
C                                                j.r.10.5.93
       IF(IP.GE.1)GO TO 1779
        PVQPZ2=PVQE**2-PVQPX**2-PVQPY**2
        IF(PVQPZ2.GE.0.)THEN
          PVQPZ=SQRT(PVQPZ2)
        ELSE
          PVQPX=0.
          PVQPY=0.
          PVQPZ=PVQE
        ENDIF
C
        PDQPZ2=PVDQE**2-PVDQPX**2-PVDQPY**2
        IF(PDQPZ2.GE.0.)THEN
          PVDQPZ=SQRT(PDQPZ2)
        ELSE
          PVDQPX=0.
          PVDQPY=0.
          PVDQPZ=PVDQE
        ENDIF
C
        TSQPZ2=TSQE**2-TSQPX**2-TSQPY**2
        IF(TSQPZ2.GE.0.)THEN
          TSQPZ=-SQRT(TSQPZ2)
        ELSE
          TSQPX=0.
          TSQPY=0.
          TSQPZ=TSQE
        ENDIF
C
        TAQPZ2=TSAQE**2-TSAQPX**2-TSAQPY**2
        IF(TAQPZ2.GE.0.)THEN
          TSAQPZ=-SQRT(TAQPZ2)
        ELSE
          TSAQPX=0.
          TSAQPY=0.
          TSAQPZ=TSAQE
        ENDIF
 1779  CONTINUE
C                                            ----------------
C
C***  SAMPLE PARTON-PT VALUES / DETERMINE PARTON 4-MOMENTA AND CHAIN MAS
C***                            IN THE REST FRAME DEFINED ABOVE
C
        IKVALA=0
C                                      changej.r.6.5.93
        PTXSQ1=0.
        PTXSA1=0.
        PTXSQ2=0.
        PTXSA2=0.
        PTYSQ1=0.
        PTYSA1=0.
        PTYSQ2=0.
        PTYSA2=0.
        PTXSQ1=PVQPX
        PTXSA1=PVDQPX
        PTXSQ2=TSQPX
        PTXSA2=TSAQPX
        PTYSQ1=PVQPY
        PTYSA1=PVDQPY
        PTYSQ2=TSQPY
        PTYSA2=TSAQPY
        PLQ1=PVQPZ
        PLAQ1=PVDQPZ
        PLQ2=TSQPZ
        PLAQ2=TSAQPZ
        EQ1=PVQE
        EAQ1=PVDQE
        EQ2=TSQE
        EAQ2=TSAQE
C                                       ---------------
          IF(IPEV.GE.1) THEN
            WRITE(6,1140) IRVS13
            WRITE(6,1090)  PTXSQ1,
     +      PTYSQ1,PLQ1,EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,PTYSQ2,
     +      PLQ2,EQ2,PTXSA2,PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2,IREJ,IKVALA,
     +      PTTQ1,PTTA1
          ENDIF
        IKVALA=0
        NSELPT=1
        NSELPT=0
	IF(IP.EQ.1)NSELPT=1
       IF(IOUXEV.GE.6)WRITE(6,'(A)')' KKEVVS call SELPT'
        IF(NSELPT.EQ.1)CALL SELPT( PTXSQ1,PTYSQ1,PLQ1,
     +  EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,
     +  PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1,
     *       PTTQ2,PTTA2,
     +  NSELPT)
        IF(NSELPT.EQ.0)CALL SELPT4( PTXSQ1,PTYSQ1,PLQ1,
     +  EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,
     +  PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1,
     +  NSELPT)
        IF (IPEV.GE.1) WRITE(6,1070) IREJ
        IF (IREJ.EQ.1) THEN
          IRVS13=IRVS13 + 1
          IF(IPEV.GE.1) THEN
            WRITE(6,1140) IRVS13
            WRITE(6,1090) PTXSQ1,
     +      PTYSQ1,PLQ1,EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,PTYSQ2,
     +      PLQ2,EQ2,PTXSA2,PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2,IREJ,IKVALA,
     +      PTTQ1,PTTA1
          ENDIF
                                                                GO TO 20
        ENDIF
C
C***  4-MOMENTA OF CHAINS IN THIS FRAME
C
        PTXCH1=PTXSQ1 + PTXSA2
        PTYCH1=PTYSQ1 + PTYSA2
        PTZCH1=PLQ1 + PLAQ2
        ECH1=EQ1 + EAQ2
        PTXCH2=PTXSQ2 + PTXSA1
        PTYCH2=PTYSQ2 + PTYSA1
        PTZCH2=PLQ2 + PLAQ1
        ECH2=EQ2 + EAQ1
        AMMM=SQRT((ECH1+ECH2)**2-(PTXCH1+PTXCH2)**2
     +            -(PTYCH1+PTYCH2)**2-(PTZCH1+PTZCH2)**2)
C
C
        IF (IPEV.GE.6)WRITE(6,1070) IREJ,
     +  AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1, AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2
 
C
C  REPLACE SMALL MASS CHAINS BY PSEUDOSCALAR OR VECTOR MESONS OR OCTETT
C                                              OR DECUPLETT BARYONS
C  FIRST FOR CHAIN 1  (PROJ VAL-QUARK  -  TARGET SEA-AQUARK)
C
        CALL COMCMA(IPVQ(IXVPR),ITSAQ(IXSTA), IJNCH1,NNCH1,IREJ,AMCH1,
     +  AMCH1N)
C***                                       MASS BELOW PSEUDOSCALAR MASS
        IF(IREJ.EQ.1) THEN
          IRVS11=IRVS11 + 1
          IF(IPEV.GE.1) THEN
            WRITE(6,1150) IRVS11
            WRITE(6,1110) IPVQ(IXVPR),ITSAQ(IXSTA),IJNCH1,NNCH1,IREJ,
     +      XPVQ(IXVPR),XPVD(IXVPR),XPVQCM,XPVDCM, XTSQ(IXSTA),XTSAQ
     +      (IXSTA),XTSQCM,XTSACM, AMCH1,AMCH1N
 
          ENDIF
                                                                 GOTO 20
        ENDIF
C                                 CORRECT KINEMATICS FOR CHAIN 1
C***                MOMENTUM CORRECTION FOR CHANGED MASS OF CHAIN 1
        IF(NNCH1.NE.0) THEN
             CALL CORMOM(AMCH1,AMCH2,AMCH1N,AMCH2N,
     +  PTXSQ1,PTYSQ1,PLQ1,EQ1,
     +  PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,PTYSA2,
     +  PLAQ2,EAQ2, PTXCH1,PTYCH1,PTZCH1,ECH1, PTXCH2,PTYCH2,PTZCH2,
     +  ECH2,IREJ)
          AMCH2=AMCH2N
        ENDIF
        IF(IREJ.EQ.1)THEN
	  IF(IPEV.EQ.1)WRITE(6,'(A)')' VS CORMOM REJECTION'
	  GO TO 20
	ENDIF
C
        IF (IPEV.GE.6)WRITE(6,1080) IREJ,
     +  AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1, AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2
 
C
C  SECOND FOR CHAIN 2  (PROJ VAL-DIQUARK  -  TAR SEA-QUARK)
C
        CALL COBCMA(ITSQ(IXSTA),IPPV1(IXVPR),IPPV2(IXVPR), IJNCH2,NNCH2,
     +  IREJ,AMCH2,AMCH2N,2)
C***                            MASS BELOW OCTETT BARYON MASS
C
C                           AT PRESENT NO CORRECTION FOR CHAIN 2
        IF(IREJ.EQ.1) THEN
          IRVS12=IRVS12 + 1
          IF(IPEV.GE.1) THEN
            WRITE(6,1160) IRVS12
            WRITE(6,1100) IPPV1(IXVPR),IPPV2(IXVPR),ITSQ(IXSTA), IJNCH2,
     +      NNCH2,IREJ, XPVQ(IXVPR),XPVD(IXVPR),XPVQCM,XPVDCM, XTSQ
     +      (IXSTA),XTSAQ(IXSTA),AMCH2,AMCH2N
 
          ENDIF
                                                                 GOTO 20
        ENDIF
        IF(NNCH2.NE.0) THEN
          AMCH2=AMCH2N
         AMMM=SQRT((ECH1+ECH2)**2-(PTXCH1+PTXCH2)**2
     +            -(PTYCH1+PTYCH2)**2-(PTZCH1+PTZCH2)**2)
        EEE=ECH1+ECH2
        PXXX=PTXCH1+PTXCH2
        PYYY=PTYCH1+PTYCH2
        PZZZ=PTZCH1+PTZCH2
        GAMMM=EEE/(AMMM+1.E-4)
        BGGGX=PXXX/(AMMM+1.E-4)
        BGGGY=PYYY/(AMMM+1.E-4)
        BGGGZ=PZZZ/(AMMM+1.E-4)
C                        TRANSFORM BOTH CHAINS  INTO TWO CHAIN-CMS
C
C                                   4-MOMENTA OF CHAINS
        CALL DALTRA(GAMMM,-BGGGX,-BGGGY,-BGGGZ,
     +  PTXCH1,PTYCH1,PTZCH1,ECH1,
     +  PPPCH1, QTXCH1,QTYCH1,QTZCH1,QECH1)
C
        CALL DALTRA(GAMMM,-BGGGX,-BGGGY,-BGGGZ,
     +  PTXCH2,PTYCH2,PTZCH2,ECH2,
     +  PPPCH2, QTXCH2,QTYCH2,QTZCH2,QECH2)
C
 
C                                if AMCH2 changed in COBCMA/COMCMA
C                                CORVAL corrects chain kinematics
C                                according to 2-body kinem.
C                                with fixed masses
           NORIG=23
          CALL CORVAL(AMMM,IREJ,AMCH1,AMCH2, QTXCH1,QTYCH1,QTZCH1,QECH1,
     +    QTXCH2,QTYCH2,QTZCH2,QECH2,NORIG)
C                        TRANSFORM BOTH CHAINS  INTO TWO CHAIN-CMS
C
C                                   4-MOMENTA OF CHAINS
C	  IREJ=1
          IF(IREJ.EQ.1) THEN
	  IF(IPEV.GE.1)WRITE(6,'(A)')' vs14 rej. '
C                           AMCH1N + AMCH2N > AMMM - 0.2
C                           REJECT EVENT
            IRVS14=IRVS14+1
                                                                 GOTO 20
          ENDIF
 
        CALL DALTRA(GAMMM,BGGGX,BGGGY,BGGGZ, QTXCH1,QTYCH1,QTZCH1,QECH1,
     +  PPPCH1, PTXCH1,PTYCH1,PTZCH1,ECH1)
C
        CALL DALTRA(GAMMM,BGGGX,BGGGY,BGGGZ, QTXCH2,QTYCH2,QTZCH2,QECH2,
     +  PPPCH2, PTXCH2,PTYCH2,PTZCH2,ECH2)
C
 
C
          IF(IPEV.GE.6) THEN
            WRITE(6,'(A/3(1PE15.4),3I5)')
     +      ' VS - CALL CORVAL: AMMM, AMCH1, AMCH2, NNCH1, NNCH2, IREJ',
     +      AMMM, AMCH1, AMCH2, NNCH1, NNCH2, IREJ
            WRITE(6,1080) IREJ, AMCH1,
     +      PTXCH1,PTYCH1,PTZCH1,ECH1, AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2
 
          ENDIF
C
          IF(IREJ.EQ.1) THEN
	  IF(IPEV.GE.1)WRITE(6,'(A)')' vs14 rej. '
C                           AMCH1N + AMCH2N > AMMM - 0.2
C                           REJECT EVENT
            IRVS14=IRVS14+1
                                                                 GOTO 20
          ENDIF
        ENDIF
        QTXCH1=PTXCH1
       QTYCH1=PTYCH1
       QTZCH1=PTZCH1
       QECH1=ECH1
       QTXCH2=PTXCH2
       QTYCH2=PTYCH2
       QTZCH2=PTZCH2
       QECH2=ECH2
       PQVSA1(N,1)=PTXSQ1
       PQVSA1(N,2)=PTYSQ1
       PQVSA1(N,3)=PLQ1
       PQVSA1(N,4)=EQ1
       PQVSA2(N,1)=PTXSA2
       PQVSA2(N,2)=PTYSA2
       PQVSA2(N,3)=PLAQ2
       PQVSA2(N,4)=EAQ2
       PQVSB1(N,1)=PTXSQ2
       PQVSB1(N,2)=PTYSQ2
       PQVSB1(N,3)=PLQ2
       PQVSB1(N,4)=EQ2
       PQVSB2(N,1)=PTXSA1
       PQVSB2(N,2)=PTYSA1
       PQVSB2(N,3)=PLAQ1
       PQVSB2(N,4)=EAQ1
C-------------------
 
C
C                                      PUT V-S CHAIN ENDS INTO /HKKEVT/
C                                      MOMENTA IN NN-CMS
C                                      POSITION OF ORIGINAL NUCLEONS
C
        IHKKPD=JHKKPV(IXVPR )
        IHKKPO=JHKKPV(IXVPR )-1
        IHKKTD=JHKKTS(IXSTA )
        IHKKTO=JHKKTS(IXSTA )-1
        IF (IPEV.GT.3)WRITE(6,1030)IXVPR,INUCPR,JNUCPR,IHKKPO,IHKKPD
 1030 FORMAT (' VS: IXVPR,INUCPR,JNUCPR,IHKKPO,IHKKPD ',5I5)
        IF (IPEV.GT.3)WRITE(6,1040)IXSTA,INUCTA,JNUCTA,IHKKTO,IHKKTD
 1040 FORMAT (' VS: IXSTA,INUCTA,JNUCTA,IHKKTO,IHKKTD ',5I5)
C                                     CHAIN 1 PROJECTILE QUARK
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=121
        IDHKK(IHKK)=IDHKK(IHKKPO)
        JMOHKK(1,IHKK)=IHKKPO
        JMOHKK(2,IHKK)=JMOHKK(1,IHKKPO)
        JDAHKK(1,IHKK)=IHKK+2
        JDAHKK(2,IHKK)=IHKK+2
        PHKK(1,IHKK)=PQVSA1(N,1)
        PHKK(2,IHKK)=PQVSA1(N,2)
        PHKK(3,IHKK)=PQVSA1(N,3)
        PHKK(4,IHKK)=PQVSA1(N,4)
        PHKK(5,IHKK)=0.
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
        VHKK(1,IHKK)=VHKK(1,IHKKPO)+XXPP
        VHKK(2,IHKK)=VHKK(2,IHKKPO)+YYPP
        VHKK(3,IHKK)=VHKK(3,IHKKPO)
        VHKK(4,IHKK)=VHKK(4,IHKKPO)
        IF (IPHKK.GE.2) WRITE(6,1050) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +  JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +  (PHKK(KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
 
 1050 FORMAT (I6,I4,5I6,9E10.2)
C                                     CHAIN 1 TARGET SEA-QUARK
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=132
        IDHKK(IHKK)=IDHKK(IHKKTD)
        JMOHKK(1,IHKK)=IHKKTD
        JMOHKK(2,IHKK)=JMOHKK(1,IHKKTD)
        JDAHKK(1,IHKK)=IHKK+1
        JDAHKK(2,IHKK)=IHKK+1
        PHKK(1,IHKK)=PQVSA2(N,1)
        PHKK(2,IHKK)=PQVSA2(N,2)
        PHKK(3,IHKK)=PQVSA2(N,3)
        PHKK(4,IHKK)=PQVSA2(N,4)
        PHKK(5,IHKK)=0.
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
        VHKK(1,IHKK)=VHKK(1,IHKKTD)+XXPP
        VHKK(2,IHKK)=VHKK(2,IHKKTD)+YYPP
        VHKK(3,IHKK)=VHKK(3,IHKKTD)
        VHKK(4,IHKK)=VHKK(4,IHKKTD)
        IF (IPHKK.GE.2) WRITE(6,1050) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +  JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +  (PHKK(KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
 
C
C                                     CHAIN 1 BEFORE FRAGMENTATION
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=5
        IDHKK(IHKK)=88888+NNCH1
        JMOHKK(1,IHKK)=IHKK-2
        JMOHKK(2,IHKK)=IHKK-1
        PHKK(1,IHKK)=QTXCH1
        PHKK(2,IHKK)=QTYCH1
        PHKK(3,IHKK)=QTZCH1
        PHKK(4,IHKK)=QECH1
        PHKK(5,IHKK)=AMCH1
C                             POSITION OF CREATED CHAIN IN LAB
C                             =POSITION OF TARGET NUCLEON
C                             TIME OF CHAIN CREATION IN LAB
C                             =TIME OF PASSAGE OF PROJECTILE
C                              NUCLEUS AT POSITION OF TAR. NUCLEUS
        VHKK(1,NHKK)= VHKK(1,NHKK-1)
        VHKK(2,NHKK)= VHKK(2,NHKK-1)
        VHKK(3,NHKK)= VHKK(3,NHKK-1)
        VHKK(4,NHKK)=VHKK(3,NHKK)/BETP-VHKK(3,NHKK-2)/BGAMP
        MHKKVS(N)=IHKK
        IF (IPROJK.EQ.1)THEN
          WHKK(1,NHKK)= VHKK(1,NHKK-2)
          WHKK(2,NHKK)= VHKK(2,NHKK-2)
          WHKK(3,NHKK)= VHKK(3,NHKK-2)
          WHKK(4,NHKK)=VHKK(4,NHKK)*GAMP-VHKK(3,NHKK)*BGAMP
          IF (IPHKK.GE.2) WRITE(6,1050) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     +    JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +    (PHKK(KHKK,NHKK),KHKK=1,5), (WHKK(KHKK,NHKK),KHKK=1,4)
 
        ENDIF
        IF (IPHKK.GE.2) WRITE(6,1050) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     +  JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +  (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
C
C
C                                     CHAIN 2 PROJECTILE DIQUARK
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=121
        IDHKK(IHKK)=IDHKK(IHKKPD)
        JMOHKK(1,IHKK)=IHKKPD
        JMOHKK(2,IHKK)=JMOHKK(1,IHKKPD)
        JDAHKK(1,IHKK)=IHKK+2
        JDAHKK(2,IHKK)=IHKK+2
        PHKK(1,IHKK)=PQVSB1(N,1)
        PHKK(2,IHKK)=PQVSB1(N,2)
        PHKK(3,IHKK)=PQVSB1(N,3)
        PHKK(4,IHKK)=PQVSB1(N,4)
        PHKK(5,IHKK)=0.
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
        VHKK(1,IHKK)=VHKK(1,IHKKPD)+XXPP
        VHKK(2,IHKK)=VHKK(2,IHKKPD)+YYPP
        VHKK(3,IHKK)=VHKK(3,IHKKPD)
        VHKK(4,IHKK)=VHKK(4,IHKKPD)
        IF (IPHKK.GE.2) WRITE(6,1050) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +  JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +  (PHKK(KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
 
C                                     CHAIN 2 TARGET SEA-QUARK
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=132
        IDHKK(IHKK)=IDHKK(IHKKTO)
        JMOHKK(1,IHKK)=IHKKTO
        JMOHKK(2,IHKK)=JMOHKK(1,IHKKTO)
        JDAHKK(1,IHKK)=IHKK+1
        JDAHKK(2,IHKK)=IHKK+1
        PHKK(1,IHKK)=PQVSB2(N,1)
        PHKK(2,IHKK)=PQVSB2(N,2)
        PHKK(3,IHKK)=PQVSB2(N,3)
        PHKK(4,IHKK)=PQVSB2(N,4)
        PHKK(5,IHKK)=0.
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
        VHKK(1,IHKK)=VHKK(1,IHKKTO)+XXPP
        VHKK(2,IHKK)=VHKK(2,IHKKTO)+YYPP
        VHKK(3,IHKK)=VHKK(3,IHKKTO)
        VHKK(4,IHKK)=VHKK(4,IHKKTO)
        IF (IPHKK.GE.2) WRITE(6,1050) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +  JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +  (PHKK(KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
 
C
C                                     CHAIN 2 BEFORE FRAGMENTATION
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=5
        IDHKK(IHKK)=88888+NNCH2
        JMOHKK(1,IHKK)=IHKK-2
        JMOHKK(2,IHKK)=IHKK-1
        PHKK(1,IHKK)=QTXCH2
        PHKK(2,IHKK)=QTYCH2
        PHKK(3,IHKK)=QTZCH2
        PHKK(4,IHKK)=QECH2
        PHKK(5,IHKK)=AMCH2
C                             POSITION OF CREATED CHAIN IN LAB
C                             =POSITION OF TARGET NUCLEON
C                             TIME OF CHAIN CREATION IN LAB
C                             =TIME OF PASSAGE OF PROJECTILE
C                              NUCLEUS AT POSITION OF TAR. NUCLEUS
        VHKK(1,NHKK)= VHKK(1,NHKK-1)
        VHKK(2,NHKK)= VHKK(2,NHKK-1)
        VHKK(3,NHKK)= VHKK(3,NHKK-1)
        VHKK(4,NHKK)=VHKK(3,NHKK)/BETP-VHKK(3,NHKK-2)/BGAMP
        MHKKVS(N)=IHKK
        IF (IPROJK.EQ.1)THEN
          WHKK(1,NHKK)= VHKK(1,NHKK-2)
          WHKK(2,NHKK)= VHKK(2,NHKK-2)
          WHKK(3,NHKK)= VHKK(3,NHKK-2)
          WHKK(4,NHKK)=VHKK(4,NHKK)*GAMP-VHKK(3,NHKK)*BGAMP
          IF (IPHKK.GE.2) WRITE(6,1050) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +    JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +    (PHKK(KHKK,IHKK),KHKK=1,5), (WHKK(KHKK,IHKK),KHKK=1,4)
 
        ENDIF
        IF (IPHKK.GE.2) WRITE(6,1050) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +  JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +  (PHKK(KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
 
C
C  NOW WE HAVE AN ACCEPTABLE VALENCE-SEA EVENT
C  AND PUT IT INTO THE HISTOGRAM
C
        AMCVS1(N)=AMCH1
        AMCVS2(N)=AMCH2
        GACVS1(N)=QECH1/AMCH1
        BGXVS1(N)=QTXCH1/AMCH1
        BGYVS1(N)=QTYCH1/AMCH1
        BGZVS1(N)=QTZCH1/AMCH1
        GACVS2(N)=QECH2/AMCH2
        BGXVS2(N)=QTXCH2/AMCH2
        BGYVS2(N)=QTYCH2/AMCH2
        BGZVS2(N)=QTZCH2/AMCH2
        NCHVS1(N)=NNCH1
        NCHVS2(N)=NNCH2
        IJCVS1(N)=IJNCH1
        IJCVS2(N)=IJNCH2
        IF (IPEV.GE.6)WRITE(6,1060) N, XPVQ(IXVPR),XPVD(IXVPR),XTSQ
     +  (IXSTA),XTSAQ(IXSTA), IPVQ(IXVPR),IPPV1(IXVPR),IPPV2(IXVPR),
     +  ITSQ(IXSTA),ITSAQ(IXSTA), AMCVS1(N),AMCVS2(N),GACVS1(N),GACVS2
     +  (N), BGXVS1(N),BGYVS1(N),BGZVS1(N), BGXVS2(N),BGYVS2(N),BGZVS2
     +  (N), NCHVS1(N),NCHVS2(N),IJCVS1(N),IJCVS2(N), (PQVSA1(N,JU),
     +  PQVSA2(N,JU),PQVSB1(N,JU), PQVSB2(N,JU),JU=1,4)
 
 
 
 
   10   CONTINUE
	RETURN
   20 CONTINUE
C                                     EVENT REJECTED
C                                     START A NEW ONE
        IREJVS=1
C---------------------------------------------------------------------
C
      RETURN
C
 1060 FORMAT(I10,4F12.7,5I5/10X,4F12.6/10X,6F12.6,4I5/8F15.5/8F15.5)
 1070 FORMAT (' VS IREJ ',I10/
     +' AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',5F12.5/
     +' AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ',5F12.5)
 1080 FORMAT (' VS IREJ  ',I10/
     +' AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',5F12.5/
     +' AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ',5F12.5)
C
 1090 FORMAT(' VS', 4(4E12.4/),2E12.4/2I5/4E12.4)
 1100 FORMAT(' VS',6I5/6E12.4/2E12.4)
 1110 FORMAT(' VS ',5I5/2(4E12.4/),2E12.4)
 1120 FORMAT(' VS',7I5/2(4E12.4/),2E12.4)
 1130 FORMAT(' VS',4I5/6E12.4/2E12.4)
 1140 FORMAT(' KKEVT - IRVS13=',I5)
 1150 FORMAT(' KKEVT - IRVS11=',I5)
 1160 FORMAT(' KKEVT - IRVS12=',I5)
C
      END
*-- Author :
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE KKEVSV(IREJSV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C-------------------------    TREATMENT OF SEA-VALENCE CHAIN SYSTEMS
C
*KEEP,NNCMS.
      COMMON /NNCMS/  GAMCM,BGCM,UMO,PCM,EPROJ,PPROJ
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
*KEEP,TRAFOP.
      COMMON /TRAFOP/ GAMP,BGAMP,BETP
*KEEP,NUCIMP.
      COMMON /NUCIMP/ PRMOM(5,248),TAMOM(5,248), PRMFEP,PRMFEN,TAMFEP,
     +TAMFEN, PREFEP,PREFEN,TAEFEP,TAEFEN, PREPOT(210),TAEPOT(210),
     +PREBIN,TAEBIN,FERMOD,ETACOU
*KEEP,ABRSV.
      COMMON /ABRSV/ AMCSV1(248),AMCSV2(248),GACSV1(248),GACSV2(248),
     +BGXSV1(248),BGYSV1(248),BGZSV1(248), BGXSV2(248),BGYSV2(248),
     +BGZSV2(248), NCHSV1(248),NCHSV2(248),IJCSV1(248),IJCSV2(248),
     +PQSVA1(248,4),PQSVA2(248,4), PQSVB1(248,4),PQSVB2(248,4)
*KEEP,FERMI.
      COMMON /FERMI/ PQUAR(4,248),PAQUAR(4,248), TQUAR(4,248),TAQUAR
     +(4,248)
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
*KEEP,REJEC.
      COMMON /REJEC/ IRCO1,IRCO2,IRCO3,IRCO4,IRCO5, IRSS11,IRSS12,
     +IRSS13,IRSS14, IRSV11,IRSV12,IRSV13,IRSV14, IRVS11,IRVS12,IRVS13,
     +IRVS14, IRVV11,IRVV12,IRVV13,IRVV14
*KEEP,PROJK.
      COMMON /PROJK/ IPROJK
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
      COMMON/RPTSHM/RPROJ,RTARG,BIMPAC
*KEND.
C=============================================================

      COMMON /ABRJT/XJQ1(INTMX),XJAQ1(INTMX),XJQ2(INTMX),XJAQ2(INTMX),
     *        IJJQ1(INTMX),IJJAQ1(INTMX),IJJQ2(INTMX),IJJAQ2(INTMX),
     *        AMJCH1(INTMX),AMJCH2(INTMX),GAMJH1(INTMX),GAMJH2(INTMX),
     *        BGJH1(INTMX),BGJH2(INTMX),THEJH1(INTMX),THEJH2(INTMX),
     *        BGXJH1(INTMX),BGYJH1(INTMX),BGZJH1(INTMX),
     *        BGXJH2(INTMX),BGYJH2(INTMX),BGZJH2(INTMX),
     *  PJETA1(INTMX,4),PJETA2(INTMX,4),PJETB1(INTMX,4),PJETB2(INTMX,4)
     * ,JHKKPH(INTMX),JHKKTH(INTMX),JHKKEX(INTMX),JHKKE1(INTMX)
      COMMON /NUCJTN/NONUJ1,NONUJT,NONUS1,NONUST
      COMMON /XSVTHR/ XSTHR,XVTHR,XDTHR,XSSTHR
      COMMON /ABRSOF/XSQ1(INTMX),XSAQ1(INTMX),XSQ2(INTMX),XSAQ2(INTMX),
     *        IJSQ1(INTMX),IJSAQ1(INTMX),IJSQ2(INTMX),IJSAQ2(INTMX),
     *        AMCCH1(INTMX),AMCCH2(INTMX),GAMCH1(INTMX),GAMCH2(INTMX),
     *        BGCH1(INTMX),BGCH2(INTMX),THECH1(INTMX),THECH2(INTMX),
     *        BGXCH1(INTMX),BGYCH1(INTMX),BGZCH1(INTMX),
     *        BGXCH2(INTMX),BGYCH2(INTMX),BGZCH2(INTMX),
     *        NCH1(INTMX),NCH2(INTMX),IJCH1(INTMX),IJCH2(INTMX),
     *  PSOFA1(INTMX,4),PSOFA2(INTMX,4),PSOFB1(INTMX,4),PSOFB2(INTMX,4)
     * ,JHKKPZ(INTMX),JHKKTZ(INTMX),JHKKSX(INTMX),JHKKS1(INTMX)
      COMMON /MINIJ/IMINIJ,NOMJE,NOMJER,NREJEV,NOMJT,NOMJTR
      COMMON /ZSEA/ZSEAAV,ZSEASU,ANZSEA
C
      THMOD=1.
      IMINIJ=1
      IF(IP.GT.1)THMOD=20.
C     DO 201 N=1,NSV
C-------------------
      IREJSV=0
      IF(IPEV.GE.1)THEN
      WRITE(6,6589) NVV,NSV,NVS,NSS,NDV,NVD,NDS,NSD
 6589 FORMAT(' KKEVSV: NVV,NSV,NVS,NSS,NDV,NVD,NDS,NSD',8I5)
      ENDIF
      DO 10 N=1,NSV
C---------------------------drop recombined chain pairs
        IF(NCHSV1(N).EQ.99.OR.NCHSV2(N).EQ.99)GO TO 10
C
        IXSPR=INTSV1(N)
        INUCPR=IFROSP(IXSPR)
        JNUCPR=ITOVP(INUCPR)
        IF(IPEV.GE.1)THEN
          PQP=GAMCM*PRMOM(3,INUCPR)+BGCM*PRMOM(4,INUCPR)
          PQE=GAMCM*PRMOM(4,INUCPR)+BGCM*PRMOM(3,INUCPR)
          PQPQ=GAMCM*PSQPZ+BGCM*PSQE
          PQEQ=GAMCM*PSQE+BGCM*PSQPZ
          PQPD=GAMCM*PSAQPZ+BGCM*PSAQE
          PQED=GAMCM*PSAQE+BGCM*PSAQPZ
        WRITE(6,1655)PRMOM(3,INUCPR),PRMOM(4,INUCPR),PQP,PQE,
     +  XPSQ(IXSPR),XPSAQ(IXSPR),IXSPR
        WRITE(6,1656)PVQPZ,PVQE,PQPQ,PQEQ
        WRITE(6,1657)PVDQPZ,PVDQE,PQPD,PQED
        ENDIF 
C
C
C                        SUBTRACT HARD SCATTERED X VALUES FROM DIQUARKS
C
        IIFROP=IFROSP(IXSPR)
        IXVPR=ITOVP(IIFROP)
C
        IXVTA=INTSV2(N)
        INUCTA=IFROVT(IXVTA)
        JNUCTA=ITOVT(INUCTA)
C
C
        XMAX1=XPSQ(IXSPR)+XPSAQ(IXSPR)+XPVQ(IXVPR)+XPVD(IXVPR)
     *        -2.D0*XSTHR-XVTHR-XDTHR
        XMAX2=XTVQ(IXVTA)+XTVD(IXVTA)-XVTHR-XDTHR
C
        IF(IPEV.GE.1)WRITE(6,'(A,2I5,2F9.3)')' KKEVSV,bef xptfl:n,nsv'
     *  ,N,NSV,XMAX1,XMAX2
        IF (IMINIJ.EQ.1)THEN
      	  CALL XPTFL(NHARD,NSEA,IREG,XMAX1,XMAX2)
C	  NZSEA=NZSEA+1
C	  ANZSEA=NZSEA
	  ANZSEA=ANZSEA+1.D0
	  ZSEASU=ZSEASU+NSEA
	  ZSEAAV=ZSEASU/ANZSEA
        ENDIF
        IF(IPEV.GE.1)WRITE(6,'(A,3I10)')' SV,xptfl:nhard,nsea,ireg '
     *  ,NHARD,NSEA,IREG
	IF(IREG.EQ.1)NHARD=0
	IF(IREG.EQ.1)NSEA=0
        NOMJE=NOMJE+NHARD
C
C
        IF (NHARD.GE.1.AND.IMINIJ.EQ.1)THEN
        DO 71 IXX=NONUJ1,NONUJT
          JHKKPH(IXX)=IXVPR
          JHKKEX(IXX)=0
          JHKKE1(IXX)=0
          IF (XPSQ(IXSPR)-XJQ1(IXX).GE.THMOD*XSTHR) THEN
            XPSQ(IXSPR)=XPSQ(IXSPR)-XJQ1(IXX)
            JHKKE1(IXX)=1
          ELSEIF (XPSAQ(IXSPR)-XJQ1(IXX).GE.THMOD*XSTHR) THEN
            XPSAQ(IXSPR)=XPSAQ(IXSPR)-XJQ1(IXX)
            JHKKE1(IXX)=2
          ELSEIF (XPVQ(IXVPR)-XJQ1(IXX).GE.THMOD*XVTHR) THEN
            XPVQ(IXVPR)=XPVQ(IXVPR)-XJQ1(IXX)
            JHKKE1(IXX)=3
          ELSEIF(XPVD(IXVPR)-XJQ1(IXX).GE.THMOD*XDTHR) THEN
            XPVD(IXVPR)=XPVD(IXVPR)-XJQ1(IXX)
            JHKKE1(IXX)=4
          ENDIF
   71   CONTINUE
        ENDIF
C
        IXVTA=INTSV2(N)
        INUCTA=IFROVT(IXVTA)
        JNUCTA=ITOVT(INUCTA)
C
C                        SUBTRACT HARD SCATTERED X VALUES FROM DIQUARKS
C
        IF (NHARD.GE.1.AND.IMINIJ.EQ.1) THEN
        DO 771 IXX=NONUJ1,NONUJT
          JHKKTH(IXX)=IXVTA
          IF (JHKKE1(IXX).EQ.0)THEN
            JHKKEX(IXX)=0
            GO TO 771
          ENDIF
          IF (XTVQ(IXVTA)-XJQ2(IXX).GE.THMOD*XVTHR) THEN
            XTVQ(IXVTA)=XTVQ(IXVTA)-XJQ2(IXX)
            JHKKEX(IXX)=1
            NOMJER=NOMJER+1
          ELSEIF (XTVD(IXVTA)-XJQ2(IXX).GE.THMOD*XDTHR)THEN
            XTVD(IXVTA)=XTVD(IXVTA)-XJQ2(IXX)
            JHKKEX(IXX)=1
            NOMJER=NOMJER+1
          ELSE
            JHKKEX(IXX)=0
            IF (JHKKE1(IXX).EQ.1)THEN
              XPSQ(IXSPR)=XPSQ(IXSPR)+XJQ1(IXX)
            ELSEIF (JHKKE1(IXX).EQ.2)THEN
              XPSAQ(IXSPR)=XPSAQ(IXSPR)+XJQ1(IXX)
            ELSEIF (JHKKE1(IXX).EQ.3)THEN
              XPVQ(IXVPR)=XPVQ(IXVPR)+XJQ1(IXX)
            ELSEIF (JHKKE1(IXX).EQ.4)THEN
              XPVD(IXVPR)=XPVD(IXVPR)+XJQ1(IXX)
            ENDIF
          ENDIF
  771   CONTINUE
        ENDIF
C
C
C                        SUBTRACT SECONDARY SEA X VALUES FROM DIQUARKS
C
	IF(IPEV.GE.1)WRITE(6,'(A,2I10)')' sv: NONUS1,NONUST ',
     *	NONUS1,NONUST
        IF (NSEA.GE.1)THEN
        DO 271 IXX=NONUS1,NONUST
          JHKKPZ(IXX)=IXVPR
          JHKKSX(IXX)=0
          JHKKS1(IXX)=0
          IF (XPSQ(IXSPR)-XSQ1(IXX)-XSAQ1(IXX).GT.THMOD*XSTHR)THEN
            XPSQ(IXSPR)=XPSQ(IXSPR)-XSQ1(IXX)-XSAQ1(IXX)
            JHKKS1(IXX)=3
          ELSEIF (XPSAQ(IXSPR)-XSQ1(IXX)-XSAQ1(IXX).GT.THMOD*XSTHR)THEN
            XPSAQ(IXSPR)=XPSAQ(IXSPR)-XSQ1(IXX)-XSAQ1(IXX)
            JHKKS1(IXX)=4
          ELSEIF (XPVQ(IXVPR)-XSQ1(IXX)-XSAQ1(IXX).GT.THMOD*XVTHR)THEN
            XPVQ(IXVPR)=XPVQ(IXVPR)-XSQ1(IXX)-XSAQ1(IXX)
            JHKKS1(IXX)=1
          ELSEIF (XPVD(IXVPR)-XSQ1(IXX)-XSAQ1(IXX).GT.THMOD*XDTHR)THEN
            XPVD(IXVPR)=XPVD(IXVPR)-XSQ1(IXX)-XSAQ1(IXX)
            JHKKS1(IXX)=2
          ENDIF
	IF(IPEV.GE.1)WRITE(6,'(A,2I10)')' sv:JHKKS1(IXX), SX  ',
     *	JHKKS1(IXX),JHKKSX(IXX)
	IF(IPEV.GE.1)WRITE(6,'(A,I10)')' sv:IXSPR  ',
     *	IXSPR
	IF(IPEV.GE.1)WRITE(6,'(A,2F10.2)')' sv:XPSQ(IXSPR),SAQ',
     *	XPSQ(IXSPR),XPSAQ(IXSPR)
  271   CONTINUE
        ENDIF
C
        INUCTA=IFROVT(IXVTA)
C
C                        SUBTRACT SECONDARY SEA X VALUES FROM DIQUARKS
C
        IF (NSEA.GE.1)THEN
        DO 2771 IXX=NONUS1,NONUST
          JHKKTZ(IXX)=IXVTA
          IF (JHKKS1(IXX).EQ.0)THEN
            JHKKSX(IXX)=0
            GO TO 2771
          ENDIF
          IF (XTVQ(IXVTA)-XSQ2(IXX)-XSAQ2(IXX).GT. THMOD*XVTHR) THEN
            XTVQ(IXVTA)=XTVQ(IXVTA)-XSQ2(IXX)-XSAQ2(IXX)
            JHKKSX(IXX)=1
C           NOMJER=NOMJER+1
          ELSEIF(XTVD(IXVTA)-XSQ2(IXX)-XSAQ2(IXX).GT.THMOD*XDTHR)THEN
            XTVD(IXVTA)=XTVD(IXVTA)-XSQ2(IXX)-XSAQ2(IXX)
            JHKKSX(IXX)=1
C           NOMJER=NOMJER+1
          ELSE
            JHKKSX(IXX)=0
            IF (JHKKS1(IXX).EQ.1)THEN
              XPVQ(IXVPR)=XPVQ(IXVPR)+XSQ1(IXX)+XSAQ1(IXX)
            ELSEIF(JHKKS1(IXX).EQ.2)THEN
              XPVD(IXVPR)=XPVD(IXVPR)+XSQ1(IXX)+XSAQ1(IXX)
            ELSEIF(JHKKS1(IXX).EQ.3)THEN
              XPSQ(IXSPR)=XPSQ(IXSPR)+XSQ1(IXX)+XSAQ1(IXX)
            ELSEIF(JHKKS1(IXX).EQ.4)THEN
              XPSAQ(IXSPR)=XPSAQ(IXSPR)+XSQ1(IXX)+XSAQ1(IXX)
            ENDIF
          ENDIF
	IF(IPEV.GE.1)WRITE(6,'(A,2I10)')' sv:JHKKS1(IXX), SX  ',
     *	JHKKS1(IXX),JHKKSX(IXX)
	IF(IPEV.GE.1)WRITE(6,'(A,2F10.2)')' sv:XPSQ(IXSPR),SAQ',
     *	XPSQ(IXSPR),XPSAQ(IXSPR)
 2771   CONTINUE
        ENDIF
C
C
        XMAX1=XPSQ(IXSPR)+XPSAQ(IXSPR)+XPVQ(IXVPR)+XPVD(IXVPR)
     *        -2.D0*XSTHR-XVTHR-XDTHR
        XMAX2=XTVQ(IXVTA)+XTVD(IXVTA)-XVTHR-XDTHR
C
        IF(IPEV.GE.1)WRITE(6,'(A,2I5,2F9.3)')' KKEVSV,aft xptfl:n,nsv'
     *  ,N,NSV,XMAX1,XMAX2
        GO TO 1202
 2202   CONTINUE
C
C                     TRY TO INCREASE X-FRACTION OF PROJECTILE SEA
C                     QUARK ANTIQUARK PAIR BY XSTHR
C                     XPSQ(IXSPR), XPSAQ(IXSPR)
C                     DECREASING VALENCE DIQUARK XPVD(IXVPR)
C
C       IF (XPSUT(IXVPR).EQ.0..AND.
C    *      XPVD(IXVPR)-2.*XSTHR.GE.XDTHR) THEN
C         XPSQ(IXSPR)=XPSQ(IXSPR)+XSTHR
C         XPSAQ(IXSPR)=XPSAQ(IXSPR)+XSTHR
C         XPVD(IXVPR)=XPVD(IXVPR)-2.*XSTHR
C         IREJ=0
C       ELSE
C         GO TO 202
C         GO TO 10
C       ENDIF
 1202   CONTINUE
C============================================================
C-------------------
C     IREJSV=0
C     IF(IPEV.GE.1)THEN
C     WRITE(6,6589) NVV,NSV,NVS,NSS,NDV,NVD,NDS,NSD
C6589 FORMAT(' KKEVSV: NVV,NSV,NVS,NSS,NDV,NVD,NDS,NSD',8I5)
C     ENDIF
C     DO 10 N=1,NSV
C---------------------------drop recombined chain pairs
C       IF(NCHSV1(N).EQ.99.OR.NCHSV2(N).EQ.99)GO TO 10
      IF(IPEV.GE.1)THEN
      WRITE(6,6588)NCHSV1(N),NCHSV2(N)
 6588 FORMAT(' NCHSV1(N),NCHSV2(N)',2I5)
      ENDIF
C
C***                 4-MOMENTA OF PROJECTILE SEA-QUARK PAIRS IN NN-CMS
        IXSPR=INTSV1(N)
        INUCPR=IFROSP(IXSPR)
        JNUCPR=ITOVP(INUCPR)
C
        PRAMOM=SQRT(PRMOM(1,INUCPR)**2
     +  +PRMOM(2,INUCPR)**2
     +  +PRMOM(3,INUCPR)**2)
        IF(PRAMOM.EQ.0.)THEN
          XXQQ=1.
        ELSE
          XXQQ=PRMOM(4,INUCPR)/PRAMOM
        ENDIF
          XXQQ=1.
        PSQPX=XPSQ(IXSPR)*PRMOM(1,INUCPR)*XXQQ
        PSQPY=XPSQ(IXSPR)*PRMOM(2,INUCPR)*XXQQ
        PSQPZ=XPSQ(IXSPR)*PRMOM(3,INUCPR)*XXQQ
        PSQE=XPSQ(IXSPR)*PRMOM(4,INUCPR)
        PSAQPX=XPSAQ(IXSPR)*PRMOM(1,INUCPR)*XXQQ
        PSAQPY=XPSAQ(IXSPR)*PRMOM(2,INUCPR)*XXQQ
        PSAQPZ=XPSAQ(IXSPR)*PRMOM(3,INUCPR)*XXQQ
        PSAQE=XPSAQ(IXSPR)*PRMOM(4,INUCPR)
        IF(IPEV.GE.1)THEN
          PQP=GAMCM*PRMOM(3,INUCPR)+BGCM*PRMOM(4,INUCPR)
          PQE=GAMCM*PRMOM(4,INUCPR)+BGCM*PRMOM(3,INUCPR)
          PQPQ=GAMCM*PSQPZ+BGCM*PSQE
          PQEQ=GAMCM*PSQE+BGCM*PSQPZ
          PQPD=GAMCM*PSAQPZ+BGCM*PSAQE
          PQED=GAMCM*PSAQE+BGCM*PSAQPZ
C  DO III=1,200
        WRITE(6,1655)PRMOM(3,INUCPR),PRMOM(4,INUCPR),PQP,PQE,
     +  XPSQ(IXSPR),XPSAQ(IXSPR),IXSPR
C  ENDDO
 1655     FORMAT(' sv PQP,PQE ',6E12.3,I5)
C	  DO III=1,200
        WRITE(6,1656)PVQPZ,PVQE,PQPQ,PQEQ
C	  ENDDO
 1656     FORMAT(' sv PQPQ,PQEQ ',4E15.5)
C	  DO III=1,200
        WRITE(6,1657)PVDQPZ,PVDQE,PQPD,PQED
C	  ENDDO
 1657     FORMAT(' sv PQPD,PQED ',4E15.5)
        ENDIF 
C
C***                 4-MOMENTA OF TARGET QUARK-DIQUARK PAIRS IN NN-CMS
        IXVTA=INTSV2(N)
        INUCTA=IFROVT(IXVTA)
        JNUCTA=ITOVT(INUCTA)
C
        TAAMOM=SQRT(TAMOM(1,INUCPR)**2
     +  +TAMOM(2,INUCPR)**2
     +  +TAMOM(3,INUCPR)**2)
        IF(TAAMOM.EQ.0.)THEN
          XXQQ=1.
        ELSE
          XXQQ=TAMOM(4,INUCTA)/TAAMOM
        ENDIF
          XXQQ=1.
        TVQPX=XTVQ(IXVTA)*TAMOM(1,INUCTA)*XXQQ
        TVQPY=XTVQ(IXVTA)*TAMOM(2,INUCTA)*XXQQ
        TVQPZ=XTVQ(IXVTA)*TAMOM(3,INUCTA)*XXQQ
        TVQE=XTVQ(IXVTA)*TAMOM(4,INUCTA)
        TVDQPX=XTVD(IXVTA)*TAMOM(1,INUCTA)*XXQQ
        TVDQPY=XTVD(IXVTA)*TAMOM(2,INUCTA)*XXQQ
        TVDQPZ=XTVD(IXVTA)*TAMOM(3,INUCTA)*XXQQ
        TVDQE=XTVD(IXVTA)*TAMOM(4,INUCTA)
        IF(PSAQE.LT.0..OR.PSQE.LE.0..OR.TVDQE.LT.0..OR.TVQE.LT.0.)
     + THEN
C	  DO III=1,200
          WRITE(6,7799)PSQPX,PSQPY,PSQPZ,PSQE,
     +    PSAQPX,PSAQPY,PSAQPZ, PSAQE,
     +    TVQPX,TVQPY,TVQPZ,TVQE,
     +    TVDQPX,TVDQPY,TVDQPZ,TVDQE
 7799     FORMAT('PSQPX,PSQPY,PSQPZ,PSQE,PSAQPX,PSAQPY,PSAQPZ
     +    PSAQE,TVQPX,TVQPY,TVQPZ,TVQE,TVDQPX,TVDQPY,TVDQPZ,TVDQE',
     +  /4(4E15.5/))
          WRITE (6,7798)IXSPR,INUCPR,IXVTA,INUCTA,
     +     XPSQ(IXSPR),XPSAQ(IXSPR),XTVQ(IXVTA),XTVD(IXVTA),
     +     PRMOM(4,INUCPR),TAMOM(4,INUCTA)
 7798     FORMAT('IXSPR,INUCPR,IXVTA,INUCTA,
     +     XPSQ(IXSPR),XPSAQ(IXSPR),XTVQ(IXVTA),XTVD(IXVTA),
     +     PRMOM(4,INUCPR),TAMOM(4,INUCTA)'/4I10/4E15.5/2E15.5)
C	  ENDDO
        ENDIF
        IF(IPEV.GE.1)THEN
          TQP=GAMCM*TAMOM(3,INUCTA)+BGCM*TAMOM(4,INUCTA)
          TQE=GAMCM*TAMOM(4,INUCTA)+BGCM*TAMOM(3,INUCTA)
          TQPQ=GAMCM*TVQPZ+BGCM*TVQE
          TQEQ=GAMCM*TVQE+BGCM*TVQPZ
          TQPD=GAMCM*TVDQPZ+BGCM*TVDQE
          TQED=GAMCM*TVDQE+BGCM*TVDQPZ
C	  DO III=1,200
        WRITE(6,1455)TAMOM(3,INUCTA),TAMOM(4,INUCTA),TQP,TQE
 1455     FORMAT(' sv TQP,TQE ',4F12.5)
        WRITE(6,1456)TVQPZ,TVQE,TQPQ,TQEQ
 1456     FORMAT(' sv TQPQ,TQEQ ',4F12.5)
        WRITE(6,1457)TVDQPZ,TVDQE,TQPD,TQED
 1457     FORMAT(' sv TQPD,TQED ',4E15.5)
          WRITE(6,7799)PSQPX,PSQPY,PSQPZ,PSQE,
     +    PSAQPX,PSAQPY,PSAQPZ, PSAQE,
     +    TVQPX,TVQPY,TVQPZ,TVQE,
     +    TVDQPX,TVDQPY,TVDQPZ,TVDQE
          WRITE (6,7798)IXSPR,INUCPR,IXVTA,INUCTA,
     +     XPSQ(IXSPR),XPSAQ(IXSPR),XTVQ(IXVTA),XTVD(IXVTA),
     +     PRMOM(4,INUCPR),TAMOM(4,INUCTA)
C	  ENDDO
        ENDIF
C                                               j.r.6.5.93
C
C                     multiple scattering of valence quark chain ends
C
      IF(IT.GT.1)THEN
      ITNU=IP+INUCTA
      RTIX=(VHKK(1,ITNU))*1.E12-BIMPAC*0.1
      RTIY=VHKK(2,ITNU)*1.E12
      RTIZ=VHKK(3,ITNU)*1.E12
      CALL CROMSC(TVQPX,TVQPY,TVQPZ,TVQE,RTIX,RTIY,RTIZ,
     *            TVQNX,TVQNY,TVQNZ,TVQNE,13)
      TVQPX=TVQNX
      TVQPY=TVQNY
      TVQPZ=TVQNZ
      TVQE=TVQNE
      CALL CROMSC(TVDQPX,TVDQPY,TVDQPZ,TVDQE,RTIX,RTIY,RTIZ,
     *            TVDQNX,TVDQNY,TVDQNZ,TVDQNE,14)
      TVDQPX=TVDQNX
      TVDQPY=TVDQNY
      TVDQPZ=TVDQNZ
      TVDQE=TVDQNE
C                                               j.r.6.5.93
C
C                     multiple scattering of sea quark chain ends
C
      ITNU=IP+INUCTA
      RTIX=(VHKK(1,ITNU))*1.E12-BIMPAC*0.1
      RTIY=VHKK(2,ITNU)*1.E12
      RTIZ=VHKK(3,ITNU)*1.E12
      CALL CROMSC(PSQPX,PSQPY,PSQPZ,PSQE,RTIX,RTIY,RTIZ,
     *            PSQNX,PSQNY,PSQNZ,PSQNE,15)
      PSQPX=PSQNX
      PSQPY=PSQNY
      PSQPZ=PSQNZ
      PSQE=PSQNE
      CALL CROMSC(PSAQPX,PSAQPY,PSAQPZ,PSAQE,RTIX,RTIY,RTIZ,
     *            PSAQNX,PSAQNY,PSAQNZ,PSAQNE,16)
      PSAQPX=PSAQNX
      PSAQPY=PSAQNY
      PSAQPZ=PSAQNZ
      PSAQE=PSAQNE
      ENDIF
C                                                ---------
C
C                                                j.r.10.5.93
       IF(IP.GE.1) GO TO 1779
        PSQPZ2=PSQE**2-PSQPX**2-PSQPY**2
        IF(PSQPZ2.GE.0.)THEN
          PSQPZ=SQRT(PSQPZ2)
        ELSE
          PSQPX=0.
          PSQPY=0.
          PSQPZ=PSQE
        ENDIF
C
        PSAQP2=PSAQE**2-PSAQPX**2-PSAQPY**2
        IF(PSAQP2.GE.0.)THEN
          PSAQPZ=SQRT(PSAQP2)
        ELSE
          PSAQPX=0.
          PSAQPY=0.
          PSAQPZ=PSAQE
        ENDIF
C
        TVQPZ2=TVQE**2-TVQPX**2-TVQPY**2
        IF(TVQPZ2.GE.0.)THEN
          TVQPZ=-SQRT(TVQPZ2)
        ELSE
          TVQPX=0.
          TVQPY=0.
          TVQPZ=TVQE
        ENDIF
C
        TDQPZ2=TVDQE**2-TVDQPX**2-TVDQPY**2
        IF(TDQPZ2.GE.0.)THEN
          TVDQPZ=-SQRT(TDQPZ2)
        ELSE
          TVDQPX=0.
          TVDQPY=0.
          TVDQPZ=TVDQE
        ENDIF
 1779  CONTINUE
C                                            ----------------
 
C
C***  SAMPLE PARTON-PT VALUES / DETERMINE PARTON 4-MOMENTA AND CHAIN MAS
C***                            IN THE REST FRAME DEFINED ABOVE
C
        IKVALA=0
C                                      changej.r.6.5.93
        PTXSQ1=0.
        PTXSA1=0.
        PTXSQ2=0.
        PTXSA2=0.
        PTYSQ1=0.
        PTYSA1=0.
        PTYSQ2=0.
        PTYSA2=0.
        PTXSQ1=PSQPX
        PTXSA1=PSAQPX
        PTXSQ2=TVQPX
        PTXSA2=TVDQPX
        PTYSQ1=PSQPY
        PTYSA1=PSAQPY
        PTYSQ2=TVQPY
        PTYSA2=TVDQPY
        PLQ1=PSQPZ
        PLAQ1=PSAQPZ
        PLQ2=TVQPZ
        PLAQ2=TVDQPZ
        EQ1=PSQE
        EAQ1=PSAQE
        EQ2=TVQE
        EAQ2=TVDQE
C                                       ---------------
C                                       ---------------
          IF(IPEV.GE.1) THEN
C	    DO III=1,200
            WRITE(6,'(A,I5)') ' HAEVSV - IRSV13=',IRSV13
            WRITE(6,'(A/4(4E12.4/),2E12.4/2I5/4E12.4)')
     +      ' SV:   ...',
     +       PTXSQ1,PTYSQ1,PLQ1,EQ1,PTXSA1,PTYSA1,
     +      PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,
     +      PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +      AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1
            BPLQ1=GAMCM*PLQ1+BGCM*EQ1
            BEQ1=GAMCM*EQ1+BGCM*PLQ1
            BPLAQ1=GAMCM*PLAQ1+BGCM*EAQ1
            BEAQ1=GAMCM*EAQ1+BGCM*PLAQ1
            BPLQ2=GAMCM*PLQ2+BGCM*EQ2
            BEQ2=GAMCM*EQ2+BGCM*PLQ2
            BPLAQ2=GAMCM*PLAQ2+BGCM*EAQ2
            BEAQ2=GAMCM*EAQ2+BGCM*PLAQ2
            WRITE(6,'(A/4(4E12.4/),2E12.4/2I5/4E12.4)')
     +      ' SV:   ...',
     +       PTXSQ1,PTYSQ1,BPLQ1,BEQ1,PTXSA1,PTYSA1,
     +      BPLAQ1,BEAQ1, PTXSQ2,PTYSQ2,BPLQ2,BEQ2,
     +      PTXSA2,PTYSA2,BPLAQ2,BEAQ2,
     +      AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1
C	     ENDDO
          ENDIF
        IKVALA=0
        NSELPT=1
        NSELPT=0
	IF(IP.EQ.1)NSELPT=1
C	    DO III=1,200
       IF(IOUXEV.GE.6)WRITE(6,'(A)')' KKEVSV call SELPT'
C	    ENDDO
        IF(NSELPT.EQ.1)CALL SELPT( PTXSQ1,PTYSQ1,PLQ1,
     +  EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,
     +  PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1,
     *       PTTQ2,PTTA2,
     +  NSELPT)
        IF(NSELPT.EQ.0)CALL SELPT4( PTXSQ1,PTYSQ1,PLQ1,
     +  EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,
     +  PTYSA2,PLAQ2,EAQ2, AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1,
     +  NSELPT)
          IF(IPEV.GE.1) THEN
C	    DO III=1,200
             WRITE(6,'(A/4(4E12.4/),2E12.4/2I5/4E12.4)')
     +      ' SV:  ...',
     +       PTXSQ1,PTYSQ1,PLQ1,EQ1,PTXSA1,PTYSA1,
     +      PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +      AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1
C	    ENDDO
          ENDIF
 
        IF (IPEV.GE.1) WRITE(6,'(A/5X,I10)')
     +  'SV   ,IREJ ',
     +  IREJ
        IF (IREJ.EQ.1) THEN
          IRSV13=IRSV13 + 1
           IF(IPEV.GE.1) THEN
            WRITE(6,'(A,I5)') ' HAEVSV - IRSV13=',IRSV13
            WRITE(6,'(A/4(4E12.4/),2E12.4/2I5/4E12.4)')
     +      ' SV:   ...',
     +       PTXSQ1,PTYSQ1,PLQ1,EQ1,PTXSA1,PTYSA1,
     +      PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +      AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1
           ENDIF
                                                                GO TO 20
        ENDIF
C
C***  4-MOMENTA OF CHAINS IN THIS FRAME
C
        PTXCH1=PTXSQ1 + PTXSA2
        PTYCH1=PTYSQ1 + PTYSA2
        PTZCH1=PLQ1 + PLAQ2
        ECH1=EQ1 + EAQ2
        PTXCH2=PTXSQ2 + PTXSA1
        PTYCH2=PTYSQ2 + PTYSA1
        PTZCH2=PLQ2 + PLAQ1
        ECH2=EQ2 + EAQ1
        AMMM=SQRT((ECH1+ECH2)**2-(PTXCH1+PTXCH2)**2
     +            -(PTYCH1+PTYCH2)**2-(PTZCH1+PTZCH2)**2)
C
C
        IF (IPEV.GE.6) WRITE(6,'(A,I10/A,5F12.5/A,5F12.5)')
     +  ' SV: IREJ ',
     +  IREJ, '     AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',
     +  AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1,
     +  '     AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ', AMCH2,PTXCH2,PTYCH2,
     +  PTZCH2,ECH2
 
C
C  REPLACE SMALL MASS CHAINS BY PSEUDOSCALAR OR VECTOR MESONS OR OCTETT
C                                              OR DECUPLETT BARYONS
C  FIRST FOR CHAIN 1  (PROJ SEA-QUARK - TAR DIQUARK)
C
        CALL COBCMA(IPSQ(IXSPR),ITTV1(IXVTA),ITTV2(IXVTA), IJNCH1,NNCH1,
     +  IREJ,AMCH1,AMCH1N,1)
          IF(IPEV.GE.2) THEN
            WRITE(6,'(A,I5)') ' HAEVSV - IRSV11=',IRSV11
            WRITE(6,'(A,6I5/6E12.4/2E12.4)') ' SV:', IPSQ(IXSPR),ITTV1
     +      (IXVTA),ITTV2(IXVTA),IJNCH1,NNCH1,IREJ, XPSQ(IXSPR),XPSAQ
     +      (IXSPR),XPSQCM,XPSACM, XTVQ(IXVTA),XTVD(IXVTA),AMCH1,AMCH1N
          ENDIF
 
C***                            MASS BELOW OCTETT BARYON MASS
        IF(IREJ.EQ.1) THEN
	  IF(IPEV.GE.1)WRITE(6,'(A)')' sv11 rej.'
          IF(IPEV.GE.1) THEN
            WRITE(6,'(A,I5)') ' HAEVSV - IRSV11=',IRSV11
            WRITE(6,'(A,6I5/6E12.4/2E12.4)') ' SV:', IPSQ(IXSPR),ITTV1
     +      (IXVTA),ITTV2(IXVTA),IJNCH1,NNCH1,IREJ, XPSQ(IXSPR),XPSAQ
     +      (IXSPR),XPSQCM,XPSACM, XTVQ(IXVTA),XTVD(IXVTA),AMCH1,AMCH1N
          ENDIF
          IRSV11=IRSV11 + 1
                                                                 GOTO 20
        ENDIF
C                                 CORRECT KINEMATICS FOR CHAIN 1
C***                MOMENTUM CORRECTION FOR CHANGED MASS OF CHAIN 1
        IF(NNCH1.NE.0) THEN
              CALL CORMOM(AMCH1,AMCH2,AMCH1N,AMCH2N, 
     +  PTXSQ1,PTYSQ1,PLQ1,EQ1,
     +  PTXSA1,PTYSA1,PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,PTYSA2,
     +  PLAQ2,EAQ2, PTXCH1,PTYCH1,PTZCH1,ECH1, PTXCH2,PTYCH2,PTZCH2,
     +  ECH2,IREJ)
          AMCH2=AMCH2N
        ENDIF
C
        IF (IPEV.GE.6) WRITE(6,'(A,I10/A,5F12.5/A,5F12.5)')
     +  ' SV(2): IREJ ',
     +  IREJ, '     AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',
     +  AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1,
     +  '     AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ', AMCH2,PTXCH2,PTYCH2,
     +  PTZCH2,ECH2
          IF(IPEV.GE.2) THEN
            WRITE(6,'(A,I5)') ' HAEVSV - IRSV11=',IRSV11
            WRITE(6,'(A,6I5/6E12.4/2E12.4)') ' SV:', IPSQ(IXSPR),ITTV1
     +      (IXVTA),ITTV2(IXVTA),IJNCH1,NNCH1,IREJ, XPSQ(IXSPR),XPSAQ
     +      (IXSPR),XPSQCM,XPSACM, XTVQ(IXVTA),XTVD(IXVTA),AMCH1,AMCH1N
          ENDIF
        IF(IREJ.EQ.1) THEN
	  IF(IPEV.GE.1)WRITE(6,'(A)')' sv cormom rej.'
	  GO TO 20
	ENDIF
 
C
C  REPLACE SMALL MASS CHAINS BY PSEUDOSCALAR OR VECTOR MESONS
C  SECOND FOR CHAIN 2  XPSAQ(N)---XTVQ(ITTA)   ANTIQUARK-QUARK
C
        CALL COMCMA(ITVQ(IXVTA),IPSAQ(IXSPR), IJNCH2,NNCH2,IREJ,AMCH2,
     +  AMCH2N)
C
C                           AT PRESENT NO CORRECTION FOR CHAIN 2
        IF(IREJ.EQ.1) THEN
          IRSV12=IRSV12 + 1
          IF(IPEV.GE.1) THEN
            WRITE(6,'(A,I5)') ' HAEVSV - IRSV12=',IRSV12
            WRITE(6,'(A/5I5/2(4E12.4/),2E12.4)')
     +      ' SV: ITVQ(IXVTA),IPSAQ(IXSPR),IJNCH2,NNCH2,IREJ...', ITVQ
     +      (IXVTA),IPSAQ(IXSPR),IJNCH2,NNCH2,IREJ, XPSQ(IXSPR),XPSAQ
     +      (IXSPR),XPSQCM,XPSACM, XTVQ(IXVTA),XTVD(IXVTA),XTVQCM,
     +      XTVDCM, AMCH2,AMCH2N
 
          ENDIF
                                                                 GOTO 20
        ENDIF
          AMCH2=AMCH2N
C
        IF (IPEV.GE.2) WRITE(6,'(A,I10/A,5F12.5/A,5F12.5)')
     +  ' SV: IREJ ',
     +  IREJ, '     AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',
     +  AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1,
     +  '     AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ', AMCH2,PTXCH2,PTYCH2,
     +  PTZCH2,ECH2
          IF(IPEV.GE.1) THEN
             WRITE(6,'(A/4(4E12.4/),2E12.4/2I5/4E12.4)')
     +      ' SV:  ...',
     +       PTXSQ1,PTYSQ1,PLQ1,EQ1,PTXSA1,PTYSA1,
     +      PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +      AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1
 
          ENDIF
 
        IF(NNCH2.NE.0) THEN
          AMCH2=AMCH2N
C                                IF AMCH1 CHANGED IN COBCMA/COMCMA
C                                CORRESPONDING REPLACEMENT IN CORMOM
         AMMM=SQRT((ECH1+ECH2)**2-(PTXCH1+PTXCH2)**2
     +            -(PTYCH1+PTYCH2)**2-(PTZCH1+PTZCH2)**2)
        EEE=ECH1+ECH2
        PXXX=PTXCH1+PTXCH2
        PYYY=PTYCH1+PTYCH2
        PZZZ=PTZCH1+PTZCH2
        GAMMM=EEE/(AMMM+1.E-4)
        BGGGX=PXXX/(AMMM+1.E-4)
        BGGGY=PYYY/(AMMM+1.E-4)
        BGGGZ=PZZZ/(AMMM+1.E-4)
C                        TRANSFORM BOTH CHAINS  INTO TWO CHAIN-CMS
C
C                                   4-MOMENTA OF CHAINS
        CALL DALTRA(GAMMM,-BGGGX,-BGGGY,-BGGGZ,
     +  PTXCH1,PTYCH1,PTZCH1,ECH1,
     +  PPPCH1, QTXCH1,QTYCH1,QTZCH1,QECH1)
C
        CALL DALTRA(GAMMM,-BGGGX,-BGGGY,-BGGGZ,
     +  PTXCH2,PTYCH2,PTZCH2,ECH2,
     +  PPPCH2, QTXCH2,QTYCH2,QTZCH2,QECH2)
C
          NORIG=24
          CALL CORVAL(AMMM,
     +    IREJ,AMCH1,AMCH2, QTXCH1,QTYCH1,QTZCH1,QECH1,
     +    QTXCH2,QTYCH2,QTZCH2,QECH2,NORIG)
C                        TRANSFORM BOTH CHAINS  INTO TWO CHAIN-CMS
C
C                                   4-MOMENTA OF CHAINS
 
        CALL DALTRA(GAMMM,BGGGX,BGGGY,BGGGZ, QTXCH1,QTYCH1,QTZCH1,QECH1,
     +  PPPCH1, PTXCH1,PTYCH1,PTZCH1,ECH1)
C
        CALL DALTRA(GAMMM,BGGGX,BGGGY,BGGGZ, QTXCH2,QTYCH2,QTZCH2,QECH2,
     +  PPPCH2, PTXCH2,PTYCH2,PTZCH2,ECH2)
C
 
C
          IF(IPEV.GE.6) THEN
            WRITE(6,'(A/3(1PE15.4),3I5)')
     +      ' SV - CALL CORVAL: AMMM, AMCH1, AMCH2, NNCH1, NNCH2, IREJ',
     +      AMMM, AMCH1, AMCH2, NNCH1, NNCH2, IREJ
          ENDIF
          IF(IREJ.EQ.1) THEN
	    IF(IPEV.GE.1)WRITE(6,'(A)')' sv14 rej.'
C                           AMCH1N + AMCH2N > AMMM - 0.2
C                           REJECT EVENT
            IRSV14=IRSV14+1
                                                                 GOTO 20
          ENDIF
C                         12.5.95
C         GO TO 20
        ENDIF
C
       QTXCH1=PTXCH1
       QTYCH1=PTYCH1
       QTZCH1=PTZCH1
       QECH1=ECH1
       QTXCH2=PTXCH2
       QTYCH2=PTYCH2
       QTZCH2=PTZCH2
       QECH2=ECH2
       PQSVA1(N,1)=PTXSQ1
       PQSVA1(N,2)=PTYSQ1
       PQSVA1(N,3)=PLQ1
       PQSVA1(N,4)=EQ1
       PQSVA2(N,1)=PTXSA2
       PQSVA2(N,2)=PTYSA2
       PQSVA2(N,3)=PLAQ2
       PQSVA2(N,4)=EAQ2
       PQSVB1(N,1)=PTXSQ2
       PQSVB1(N,2)=PTYSQ2
       PQSVB1(N,3)=PLQ2
       PQSVB1(N,4)=EQ2
       PQSVB2(N,1)=PTXSA1
       PQSVB2(N,2)=PTYSA1
       PQSVB2(N,3)=PLAQ1
       PQSVB2(N,4)=EAQ1
C-------------------
 
C
C                                      PUT S-V CHAIN ENDS INTO /HKKEVT/
C                                      MOMENTA IN NN-CMS
C                                      POSITION OF ORIGINAL NUCLEONS
C
C                                 FLAG FOR SV-CHAIN ENDS
C                                            PROJECTILE: ISTHKK=131
C                                            TARGET:     ISTHKK=122
C                                      FOR SV-CHAINS     ISTHKK=4
C
        IHKKPD=JHKKPS(IXSPR )
        IHKKPO=JHKKPS(IXSPR )-1
        IHKKTD=JHKKTV(IXVTA )
        IHKKTO=JHKKTV(IXVTA )-1
        IF (IPEV.GT.3)WRITE(6,1000)IXSPR,INUCPR,JNUCPR,IHKKPO,IHKKPD
 1000 FORMAT (' IXSPR,INUCPR,JNUCPR,IHKKPO,IHKKPD ',5I5)
        IF (IPEV.GT.3)WRITE(6,1010)IXVTA,INUCTA,JNUCTA,IHKKTO,IHKKTD
 1010 FORMAT (' IXVTA,INUCTA,JNUCTA,IHKKTO,IHKKTD ',5I5)
C                                     CHAIN 1 PROJECTILE SEA-QUARK
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=131
        IDHKK(IHKK)=IDHKK(IHKKPO)
        JMOHKK(1,IHKK)=IHKKPO
        JMOHKK(2,IHKK)=JMOHKK(1,IHKKPO)
        JDAHKK(1,IHKK)=IHKK+2
        JDAHKK(2,IHKK)=IHKK+2
        PHKK(1,IHKK)=PQSVA1(N,1)
        PHKK(2,IHKK)=PQSVA1(N,2)
        PHKK(3,IHKK)=PQSVA1(N,3)
        PHKK(4,IHKK)=PQSVA1(N,4)
        PHKK(5,IHKK)=0.
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
        VHKK(1,IHKK)=VHKK(1,IHKKPO)+XXPP
        VHKK(2,IHKK)=VHKK(2,IHKKPO)+YYPP
        VHKK(3,IHKK)=VHKK(3,IHKKPO)
        VHKK(4,IHKK)=VHKK(4,IHKKPO)
        IF (IPHKK.GE.2) WRITE(6,1020) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +  JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +  (PHKK(KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
 
 1020 FORMAT (I6,I4,5I6,9E10.2)
C                                     CHAIN 1 TARGET DIQUARK
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=122
        IDHKK(IHKK)=IDHKK(IHKKTD)
        JMOHKK(1,IHKK)=IHKKTD
        JMOHKK(2,IHKK)=JMOHKK(1,IHKKTD)
        JDAHKK(1,IHKK)=IHKK+1
        JDAHKK(2,IHKK)=IHKK+1
        PHKK(1,IHKK)=PQSVA2(N,1)
        PHKK(2,IHKK)=PQSVA2(N,2)
        PHKK(3,IHKK)=PQSVA2(N,3)
        PHKK(4,IHKK)=PQSVA2(N,4)
        PHKK(5,IHKK)=0.
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
        VHKK(1,IHKK)=VHKK(1,IHKKTD)+XXPP
        VHKK(2,IHKK)=VHKK(2,IHKKTD)+YYPP
        VHKK(3,IHKK)=VHKK(3,IHKKTD)
        VHKK(4,IHKK)=VHKK(4,IHKKTD)
        IF (IPHKK.GE.2) WRITE(6,1020) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +  JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +  (PHKK(KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
 
C
C                                     CHAIN 1 BEFORE FRAGMENTATION
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=4
        IDHKK(IHKK)=88888+NNCH1
        JMOHKK(1,IHKK)=IHKK-2
        JMOHKK(2,IHKK)=IHKK-1
        PHKK(1,IHKK)=QTXCH1
        PHKK(2,IHKK)=QTYCH1
        PHKK(3,IHKK)=QTZCH1
        PHKK(4,IHKK)=QECH1
        PHKK(5,IHKK)=AMCH1
C                             POSITION OF CREATED CHAIN IN LAB
C                             =POSITION OF TARGET NUCLEON
C                             TIME OF CHAIN CREATION IN LAB
C                             =TIME OF PASSAGE OF PROJECTILE
C                              NUCLEUS AT POSITION OF TAR. NUCLEUS
        IF (IPEV.GT.3)WRITE(6,'(A,3E12.3)')' BETP,GAMP,BGAMP',
     *  BETP,GAMP,BGAMP
        VHKK(1,NHKK)= VHKK(1,NHKK-1)
        VHKK(2,NHKK)= VHKK(2,NHKK-1)
        VHKK(3,NHKK)= VHKK(3,NHKK-1)
        VHKK(4,NHKK)=VHKK(3,NHKK)/BETP-VHKK(3,NHKK-2)/BGAMP
        MHKKSV(N)=IHKK
        IF (IPROJK.EQ.1)THEN
          WHKK(1,NHKK)= VHKK(1,NHKK-2)
          WHKK(2,NHKK)= VHKK(2,NHKK-2)
          WHKK(3,NHKK)= VHKK(3,NHKK-2)
          WHKK(4,NHKK)=VHKK(4,NHKK)*GAMP-VHKK(3,NHKK)*BGAMP
          IF (IPHKK.GE.2) WRITE(6,1020) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +    JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +    (PHKK(KHKK,IHKK),KHKK=1,5), (WHKK(KHKK,IHKK),KHKK=1,4)
 
        ENDIF
        IF (IPHKK.GE.2) WRITE(6,1020) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +  JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +  (PHKK(KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
 
C
C
C                                     CHAIN 2 PROJECTILE SEA-ANTIQUARK
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=131
        IDHKK(IHKK)=IDHKK(IHKKPD)
        JMOHKK(1,IHKK)=IHKKPD
        JMOHKK(2,IHKK)=JMOHKK(1,IHKKPD)
        JDAHKK(1,IHKK)=IHKK+2
        JDAHKK(2,IHKK)=IHKK+2
        PHKK(1,IHKK)=PQSVB1(N,1)
        PHKK(2,IHKK)=PQSVB1(N,2)
        PHKK(3,IHKK)=PQSVB1(N,3)
        PHKK(4,IHKK)=PQSVB1(N,4)
        PHKK(5,IHKK)=0.
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
        VHKK(1,IHKK)=VHKK(1,IHKKPD)+XXPP
        VHKK(2,IHKK)=VHKK(2,IHKKPD)+YYPP
        VHKK(3,IHKK)=VHKK(3,IHKKPD)
        VHKK(4,IHKK)=VHKK(4,IHKKPD)
        IF (IPHKK.GE.2) WRITE(6,1020) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +  JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +  (PHKK(KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
 
C                                     CHAIN 2 TARGET QUARK
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=122
        IDHKK(IHKK)=IDHKK(IHKKTO)
        JMOHKK(1,IHKK)=IHKKTO
        JMOHKK(2,IHKK)=JMOHKK(1,IHKKTO)
        JDAHKK(1,IHKK)=IHKK+1
        JDAHKK(2,IHKK)=IHKK+1
        PHKK(1,IHKK)=PQSVB2(N,1)
        PHKK(2,IHKK)=PQSVB2(N,2)
        PHKK(3,IHKK)=PQSVB2(N,3)
        PHKK(4,IHKK)=PQSVB2(N,4)
        PHKK(5,IHKK)=0.
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
        VHKK(1,IHKK)=VHKK(1,IHKKTO)+XXPP
        VHKK(2,IHKK)=VHKK(2,IHKKTO)+YYPP
        VHKK(3,IHKK)=VHKK(3,IHKKTO)
        VHKK(4,IHKK)=VHKK(4,IHKKTO)
        IF (IPHKK.GE.2) WRITE(6,1020) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +  JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +  (PHKK(KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
 
C
C                                     CHAIN 2 BEFORE FRAGMENTATION
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=4
        IDHKK(IHKK)=88888+NNCH2
        JMOHKK(1,IHKK)=IHKK-2
        JMOHKK(2,IHKK)=IHKK-1
        PHKK(1,IHKK)=QTXCH2
        PHKK(2,IHKK)=QTYCH2
        PHKK(3,IHKK)=QTZCH2
        PHKK(4,IHKK)=QECH2
        PHKK(5,IHKK)=AMCH2
C                             POSITION OF CREATED CHAIN IN LAB
C                             =POSITION OF TARGET NUCLEON
C                             TIME OF CHAIN CREATION IN LAB
C                             =TIME OF PASSAGE OF PROJECTILE
C                              NUCLEUS AT POSITION OF TAR. NUCLEUS
        VHKK(1,NHKK)= VHKK(1,NHKK-1)
        VHKK(2,NHKK)= VHKK(2,NHKK-1)
        VHKK(3,NHKK)= VHKK(3,NHKK-1)
        VHKK(4,NHKK)=VHKK(3,NHKK)/BETP-VHKK(3,NHKK-2)/BGAMP
        MHKKSV(N)=IHKK
        IF (IPROJK.EQ.1)THEN
          WHKK(1,NHKK)= VHKK(1,NHKK-2)
          WHKK(2,NHKK)= VHKK(2,NHKK-2)
          WHKK(3,NHKK)= VHKK(3,NHKK-2)
          WHKK(4,NHKK)=VHKK(4,NHKK)*GAMP-VHKK(3,NHKK)*BGAMP
          IF (IPHKK.GE.2) WRITE(6,1020) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +    JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +    (PHKK(KHKK,IHKK),KHKK=1,5), (WHKK(KHKK,IHKK),KHKK=1,4)
 
        ENDIF
        IF (IPHKK.GE.2) WRITE(6,1020) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +  JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +  (PHKK(KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
 
C
C
C  NOW WE HAVE AN ACCEPTABLE SEA--VALENCE  EVENT
C  AND PUT IT INTO THE HISTOGRAM
C
        AMCSV1(N)=AMCH1
        AMCSV2(N)=AMCH2
        GACSV1(N)=QECH1/AMCH1
        BGXSV1(N)=QTXCH1/AMCH1
        BGYSV1(N)=QTYCH1/AMCH1
        BGZSV1(N)=QTZCH1/AMCH1
        GACSV2(N)=QECH2/AMCH2
        BGXSV2(N)=QTXCH2/AMCH2
        BGYSV2(N)=QTYCH2/AMCH2
        BGZSV2(N)=QTZCH2/AMCH2
        NCHSV1(N)=NNCH1
        NCHSV2(N)=NNCH2
        IJCSV1(N)=IJNCH1
        IJCSV2(N)=IJNCH2
        IF (IPEV.GE.2) WRITE(6,'(A/I10,4F12.7,5I5/10X,4F12.6/10X,6F12.6,
     +4I5/8F15.5/                8F15.5)') ' SV / FINAL PRINT',N, XPSQ
     +  (IXSPR),XPSAQ(IXSPR),XTVQ(IXVTA),XTVD(IXVTA), IPSQ(IXSPR),IPSAQ
     +  (IXSPR), ITVQ(IXVTA),ITTV1(IXVTA),ITTV2(IXVTA), AMCSV1(N),AMCSV2
     +  (N),GACSV1(N),GACSV2(N), BGXSV1(N),BGYSV1(N),BGZSV1(N), BGXSV2
     +  (N),BGYSV2(N),BGZSV2(N), NCHSV1(N),NCHSV2(N),IJCSV1(N),IJCSV2
     +  (N), (PQSVA1(N,JU),PQSVA2(N,JU),PQSVB1(N,JU), PQSVB2(N,JU),JU=1,
     +  4)
 
 
 
 
   10 CONTINUE
      RETURN
C
   20 CONTINUE
C                                     EVENT REJECTED
C                                     START A NEW ONE
      IREJSV=1
      RETURN
      END
C-------------------------------------------------------------------

C-------------------------------------------------------------------
      SUBROUTINE CROMSC(PX,PY,PZ,E,RX,RY,RZ,PXN,PYN,PZN,EN,IORIG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C                                                     j.r.6.5.93
C      parton with momentum components px,py,pz
C      at position rx,ry,rz in target nucleus
C      gets multiple scattering during travel through target
C
      COMMON /NNCMS/GAMCM,BGCM,UMO,PCM,EPROJ,PPROJ
      COMMON/RPTSHM/RPROJ,RTARG,BIMPAC
      COMMON/CRONIN/CRONCO,MKCRON
      COMMON/DPRIN/IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
C      IPRI=1
      IF(E.LE.0.D0)IPRI=1
      IF(MKCRON.EQ.0) THEN
        PXN=PX
        PYN=PY
        PZN=PZ
        EN=E
        RETURN
      ENDIF
      IF(IPRI.GE.1)
     * WRITE(6,'(A,7E12.3,I5)')
     *        ' CROMSC:PX,PY,PZ,E,RX,RY,RZ,IORIG',
     *          PX,PY,PZ,E,RX,RY,RZ,IORIG
C
C  first transform parton momenta px,py,pz,e back into target systen
C
      IF(IPRI.GE.1)
     * WRITE(6,'(A,7E12.3)')' CROMSC:GAMCM,BGCM',
     *                              GAMCM,BGCM 
      CALL SLTRAF(GAMCM,-BGCM,E,PZ,EL,PZL)
      IF(IPRI.GE.1)
     * WRITE(6,'(A,4E12.3)')' CROMSC:E,PZ,EL,PZL',
     *                               E,PZ,EL,PZL
C
C                       direction cosines of parton
C
      PP=PX**2+PY**2+PZL**2
      P=SQRT(PP)
      IF(P.LE.2.0)THEN
        PXN=PX
        PYN=PY
        PZN=PZ
        EN=E
        RETURN
      ENDIF
C
      CX=PX/P
      CY=PY/P
      CZ=PZL/P
      IF(IPRI.GE.1)
     * WRITE(6,'(A,4E12.3)')' CROMSC:P,CX,CY,CZ',
     *                               P,CX,CY,CZ
C
C     is position of parton within standard target nucleus (r=rtarg)
C
      RTESQ= RX**2+RY**2+RZ**2-RTARG**2
      IF(IPRI.GE.1)
     * WRITE(6,'(A,2E12.3)')' CROMSC:RTARG,RTESQ',
     *                               RTARG,RTESQ
      IF(RTESQ.GE.-0.001)THEN
        PXN=PX
        PYN=PY
        PZN=PZ
        EN=E
        RETURN
      ENDIF
C
C    calculate distance from point rx,ry,rz to surface of rtarg sphere
C                                             (origin:0,0,0)
C
      B=RTESQ
      A=CX*RX+CY*RY+CZ*RZ
C                      distance to surface ts
      TS=-A+SQRT(A**2-B)
      IF(IPRI.GE.1)
     * WRITE(6,'(A,3E12.3)')' CROMSC:A,B,TS',
     *                               A,B,TS
 
C
C     calculate multiple scattering angle
C
      THETO=CRONCO*SQRT(TS)/P
C     IF(IPRI.GE.0.AND.THETO.GT.0.10D0)
C    * WRITE(6,'(A,4E12.3)')' CROMSC:A,B,TS,THETO,truncate',
C    *                               A,B,TS,THETO
C     IF(THETO.GT.0.10D0)THETO=0.1
C       PXN=PX
C       PYN=PY
C       PZN=PZ
C       EN=E
C       RETURN
C     ENDIF
 1212 CONTINUE
C
C     Gaussian sampling of space angle
C
      CALL RANNOR(R1,R2)
      THETA=ABS(R1*THETO)
C     IF(THETA.GE.0.3D0)THEN
      IF(THETA.GE.0.9D0)THEN
      IF(IPRI.GE.1)
     * WRITE(6,'(A,4E12.3)')' CROMSC:A,B,TS,THETA,reject',
     *                               A,B,TS,THETA
        PXN=PX
        PYN=PY
        PZN=PZ
        EN=E
        RETURN
      ENDIF
      CALL DSFECF(SFE,CFE)
      CT=COS(THETA)
      ST=SIN(THETA)
      IF(IPRI.GE.1)
     * WRITE(6,'(A,2E12.3)')' CROMSC:THETO,THETA',
     *                               THETO,THETA
C
C       new direction cosines
C
      CALL DTRANS(CX,CY,CZ,CT,ST,CFE,SFE,CXN,CYN,CZN)
      IF(IPRI.GE.1)
     * WRITE(6,'(A,3E12.3)')' CROMSC:CXN,CYN,CZN',
     *                               CXN,CYN,CZN
C
C        new momenta in target system
C
      PXLN=CXN*P
      PYLN=CYN*P
      PZLN=CZN*P
      IF(IPRI.GE.1)
     * WRITE(6,'(A,3E12.3)')' CROMSC:PXLN,PYLN,PZLN',
     *                               PXLN,PYLN,PZLN
C
C        transformation back into cms
C
      PXN=PXLN
      PYN=PYLN
      IF(IPRI.GE.1)
     * WRITE(6,'(A,7E12.3)')' CROMSC:GAMCM,BGCM',
     *                              GAMCM,BGCM 
      CALL SLTRAF(GAMCM,BGCM,EL,PZLN,EN,PZN)
      IF(IPRI.GE.1)
     * WRITE(6,'(A,4E12.3)')' CROMSC:PXN,PYN,PZN,EN',
     *                               PXN,PYN,PZN,EN
      IF(ABS(E-EN).GT.0.2)THEN
        THETO=THETO/2.
        GO TO 1212
      ENDIF
C      IPRI=0
      IF (E.LE.0.)IPRI=0
      RETURN
      END

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE KKEVHH
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C-------------------------    TREATMENT OF Hard scattered CHAIN SYSTEMS
C
*KEEP,NUCC.
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
*KEEP,RPTSHM.
      COMMON /RPTSHM/ RPROJ,RTARG,BIMPAC
      PARAMETER (INTMX=2488,INTMD=252)
      COMMON /TRAFOP/ GAMP,BGAMP,BETP
      COMMON /DXQX/   XPVQ(248),XPVD(248),XTVQ(248),XTVD(248),
     *                XPSQ(INTMX),XPSAQ(INTMX),XTSQ(INTMX),XTSAQ(INTMX)
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
C--------------------
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX),
     *                IFROVT(248),ITOVT(248),IFROST(INTMX),
     *            JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),JHKKNT(248),
     *                JHKKPV(INTMX),JHKKPS(INTMX),
     *                JHKKTV(INTMX),JHKKTS(INTMX),
     *                MHKKVV(INTMX),MHKKSS(INTMX),
     &                MHKKVS(INTMX),MHKKSV(INTMX),
     +                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
C-------------------
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/ ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),ZUOST(INTMX),
     *               INTLO(INTMX),INLOSS(INTMX)
C-------------------
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
C--------------------
      COMMON /NUCIMP/ PRMOM(5,248),TAMOM(5,248),
     &                PRMFEP,PRMFEN,TAMFEP,TAMFEN,
     &                PREFEP,PREFEN,TAEFEP,TAEFEN,
     &        PREPOT(210),TAEPOT(210),PREBIN,TAEBIN,FERMOD,ETACOU
C--------------------
C--------------------
      LOGICAL INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA,
     *        IPADIS,ISHMAL,LPAULI
      COMMON /DROPPT/ INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA,
     *                IPADIS,ISHMAL,LPAULI
C-------------------
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
      COMMON /REJEC/ IRCO1,IRCO2,IRCO3,IRCO4,IRCO5,IRSS11,IRSS12,IRSS13,
     *               IRSS14,
     *               IRSV11,IRSV12,IRSV13,IRSV14,
     *               IRVS11,IRVS12,IRVS13,IRVS14,
     *               IRVV11,IRVV12,IRVV13,IRVV14
      COMMON /ABRHH/ AMCHH1(INTMX),AMCHH2(INTMX),
     *               GACHH1(INTMX),GACHH2(INTMX),
     *               BGXHH1(INTMX),BGYHH1(INTMX),BGZHH1(INTMX),
     *               BGXHH2(INTMX),BGYHH2(INTMX),BGZHH2(INTMX),
     *               NCHHH1(INTMX),NCHHH2(INTMX),
     *               IJCHH1(INTMX),IJCHH2(INTMX),
     *               PQHHA1(INTMX,4),PQHHA2(INTMX,4),
     *               PQHHB1(INTMX,4),PQHHB2(INTMX,4)
C-------------------
      PARAMETER (NMXHKK=  89998)
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),
     & VHKK(4,NMXHKK),WHKK(4,NMXHKK)
C                       WHKK(4,NMXHKK) GIVES POSITIONS AND TIMES IN
C                       PROJECTILE FRAME, THE CHAINS ARE CREATED ON
C                       THE POSITIONS OF THE PROJECTILE NUCLEONS
C                       IN THE PROJECTILE FRAME (TARGET NUCLEONS IN
C                       TARGET FRAME) BOTH POSITIONS ARE THREFORE NOT
C                       COMPLETELY CONSISTENT. THE TIMES IN THE
C                       PROJECTILE FRAME HOWEVER ARE OBTAINED BY
C                       LORENTZ TRANSFORMING FROM THE LAB SYSTEM.
      COMMON /PROJK/ IPROJK
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
C
C-----------------------------------------------------------------------
      COMMON /ABRJT/XJQ1(INTMX),XJAQ1(INTMX),XJQ2(INTMX),XJAQ2(INTMX),
     *        IJJQ1(INTMX),IJJAQ1(INTMX),IJJQ2(INTMX),IJJAQ2(INTMX),
     *        AMJCH1(INTMX),AMJCH2(INTMX),GAMJH1(INTMX),GAMJH2(INTMX),
     *        BGJH1(INTMX),BGJH2(INTMX),THEJH1(INTMX),THEJH2(INTMX),
     *        BGXJH1(INTMX),BGYJH1(INTMX),BGZJH1(INTMX),
     *        BGXJH2(INTMX),BGYJH2(INTMX),BGZJH2(INTMX),
     *  PJETA1(INTMX,4),PJETA2(INTMX,4),PJETB1(INTMX,4),PJETB2(INTMX,4)
     * ,JHKKPH(INTMX),JHKKTH(INTMX),JHKKEX(INTMX),JHKKE1(INTMX)
      COMMON /NUCJTN/NONUJ1,NONUJT,NONUS1,NONUST
      COMMON /XSVTHR/ XSTHR,XVTHR,XDTHR,XSSTHR
C----------------------------------------------------------------------
      DIMENSION IHKKQ(-6:6)
      DATA IHKKQ/-6,-5,-4,-3,-1,-2,0,2,1,3,4,5,6/
C
C----------------------------------------------------------------------
      DO 101 N=1,NONUJT
        IF (JHKKEX(N).EQ.1)THEN
C
          IXVPR=JHKKPH(N)
          IXVTA=JHKKTH(N)
          IHKKPO=JHKKPV(IXVPR)
          IHKKTO=JHKKTV(IXVTA)
C
C                               Cronin multiple scattering
      IF(IT.GT.1)THEN
      ITNU=IHKKTO
      RTIX=(VHKK(1,ITNU))*1.E12-BIMPAC*0.1
      RTIY=VHKK(2,ITNU)*1.E12
      RTIZ=VHKK(3,ITNU)*1.E12
      RTIR2=(RTIX**2+RTIY**2+RTIZ**2)
      IF(RTIR2.GT.RTARG**2)THEN
        IF(IPEV.GE.2)
     *  WRITE(6,774)RTARG,RTIX,RTIY,RTIZ,BIMPAC,IHKKTO,IXVTA
 774    FORMAT(' KKEVHH: RTARG,RTIX,RTIY,RTIZ,BIMPAC,IHKKTO,IXVTA'
     *  ,5E12.4,2I10)
        GO TO 779
      ENDIF
      PVQPX=PJETA1(N,1)
      PVQPY=PJETA1(N,2)
      PVQPZ=PJETA1(N,3)
      PVQE =PJETA1(N,4)
      CALL CROMSC(PVQPX,PVQPY,PVQPZ,PVQE,RTIX,RTIY,RTIZ,
     *            PVQNX,PVQNY,PVQNZ,PVQNE,20)
      PVDQPX=PJETA2(N,1)
      PVDQPY=PJETA2(N,2)
      PVDQPZ=PJETA2(N,3)
      PVDQE =PJETA2(N,4)
      CALL CROMSC(PVDQPX,PVDQPY,PVDQPZ,PVDQE,RTIX,RTIY,RTIZ,
     *            PVDQNX,PVDQNY,PVDQNZ,PVDQNE,21)
      AMTES2=((PVQNE+PVDQNE)**2-(PVQNX+PVDQNX)**2
     *       -(PVQNY+PVDQNY)**2-(PVQNZ+PVDQNZ)**2)
      IF(AMTES2.GE.AMJCH1(N)**2.OR.AMTES2.GE.25.D0)THEN 
      PJETA1(N,1)=PVQNX
      PJETA1(N,2)=PVQNY
      PJETA1(N,3)=PVQNZ
      PJETA1(N,4)=PVQNE
      PJETA2(N,1)=PVDQNX
      PJETA2(N,2)=PVDQNY
      PJETA2(N,3)=PVDQNZ
      PJETA2(N,4)=PVDQNE
      ENDIF
C                                  MASSES OF SUBCHAINS
          XMJCH1=SQRT((PJETA1(N,4)+
     *                         PJETA2(N,4))**2
     *                       -(PJETA1(N,1)+
     *                         PJETA2(N,1))**2
     *                       -(PJETA1(N,2)+
     *                         PJETA2(N,2))**2
     *                       -(PJETA1(N,3)+
     *                         PJETA2(N,3))**2)
         IF(XMJCH1.GE.AMJCH1(N))THEN
         AMJCH1(N)=XMJCH1
C
          GAMJH1(N)=(PJETA1(N,4)+
     *                    PJETA2(N,4))/AMJCH1(N)
          BGXJH1(N)=(PJETA1(N,1)+
     *                    PJETA2(N,1))/AMJCH1(N)
          BGYJH1(N)=(PJETA1(N,2)+
     *                    PJETA2(N,2))/AMJCH1(N)
          BGZJH1(N)=(PJETA1(N,3)+
     *                    PJETA2(N,3))/AMJCH1(N)
      ENDIF
      ENDIF
 779  CONTINUE
C                                                ---------
C
C
C                                      PUT h-h CHAIN ENDS INTO /HKKEVT/
C                                      MOMENTA IN NN-CMS
C                                      POSITION OF ORIGINAL NUCLEONS
c                                      flags for h-h chain ends
c                                        projectile: isthkk=151
c                                        target:     isthkk=152
c                                        h-h chains: isthkk=7
C
          IXVPR=JHKKPH(N)
          IXVTA=JHKKTH(N)
          IHKKPO=JHKKPV(IXVPR)
          IHKKTO=JHKKTV(IXVTA)
          IF (IPEV.GT.3)WRITE(6,5002)IXVPR,IHKKPO
 5002     FORMAT (' IXVPR,IHKKPO ',5I5)
          IF (IPEV.GT.3)WRITE(6,5003)IXVTA,IHKKTO
 5003     FORMAT (' IXVTA,IHKKTO ',5I5)
C                                     CHAIN 1 PROJECTILE SEA-QUARK
          NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
          IHKK=NHKK
          ISTHKK(IHKK)=151
          IDHKK(IHKK)=IHKKQ(IJJQ1(N))
          JMOHKK(1,IHKK)=IHKKPO
          JMOHKK(2,IHKK)=IHKKPO
          JDAHKK(1,IHKK)=IHKK+2
          JDAHKK(2,IHKK)=IHKK+2
          PHKK(1,IHKK)=PJETA1(N,1)
          PHKK(2,IHKK)=PJETA1(N,2)
          PHKK(3,IHKK)=PJETA1(N,3)
          PHKK(4,IHKK)=PJETA1(N,4)
          PHKK(5,IHKK)=0.
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
          VHKK(1,IHKK)=VHKK(1,IHKKPO)+XXPP
          VHKK(2,IHKK)=VHKK(2,IHKKPO)+YYPP
          VHKK(3,IHKK)=VHKK(3,IHKKPO)
          VHKK(4,IHKK)=VHKK(4,IHKKPO)
          IF (IPHKK.GE.2) WRITE(6,5001)
     *      IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),JMOHKK(2,IHKK),
     &      JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK(KHKK,IHKK),KHKK=1,5),
     &      (VHKK(KHKK,IHKK),KHKK=1,4)
 5001     FORMAT (I6,I4,5I6,9E10.2)
C                                     CHAIN 1 TARGET SEA-QUARK
          NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
          IHKK=NHKK
          ISTHKK(IHKK)=152
          IDHKK(IHKK)=IHKKQ(IJJAQ2(N))
          JMOHKK(1,IHKK)=IHKKTO
          JMOHKK(2,IHKK)=IHKKTO
          JDAHKK(1,IHKK)=IHKK+1
          JDAHKK(2,IHKK)=IHKK+1
          PHKK(1,IHKK)=PJETA2(N,1)
          PHKK(2,IHKK)=PJETA2(N,2)
          PHKK(3,IHKK)=PJETA2(N,3)
          PHKK(4,IHKK)=PJETA2(N,4)
          PHKK(5,IHKK)=0.
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
          VHKK(1,IHKK)=VHKK(1,IHKKTO)+XXPP
          VHKK(2,IHKK)=VHKK(2,IHKKTO)+YYPP
          VHKK(3,IHKK)=VHKK(3,IHKKTO)
          VHKK(4,IHKK)=VHKK(4,IHKKTO)
          IF (IPHKK.GE.2) WRITE(6,5001)
     *      IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),JMOHKK(2,IHKK),
     &      JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK(KHKK,IHKK),KHKK=1,5),
     &      (VHKK(KHKK,IHKK),KHKK=1,4)
C
C                                     CHAIN 1 BEFORE FRAGMENTATION
          NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
          IHKK=NHKK
          ISTHKK(IHKK)=7
          IDHKK(IHKK)=88888
          JMOHKK(1,IHKK)=IHKK-2
          JMOHKK(2,IHKK)=IHKK-1
          PHKK(1,IHKK)=PHKK(1,IHKK-2)+PHKK(1,IHKK-1)
          PHKK(2,IHKK)=PHKK(2,IHKK-2)+PHKK(2,IHKK-1)
          PHKK(3,IHKK)=PHKK(3,IHKK-2)+PHKK(3,IHKK-1)
          PHKK(4,IHKK)=PHKK(4,IHKK-2)+PHKK(4,IHKK-1)
          PHKK(5,IHKK)=AMJCH1(N)
C                             POSITION OF CREATED CHAIN IN LAB
C                             =POSITION OF TARGET NUCLEON
C                             TIME OF CHAIN CREATION IN LAB
C                             =TIME OF PASSAGE OF PROJECTILE
C                              NUCLEUS AT POSITION OF TAR. NUCLEUS
          VHKK(1,NHKK)=                VHKK(1,NHKK-1)
          VHKK(2,NHKK)=                VHKK(2,NHKK-1)
          VHKK(3,NHKK)=                VHKK(3,NHKK-1)
          VHKK(4,NHKK)=0.
          MHKKHH(N)=IHKK
          IF (IPROJK.EQ.1)THEN
            WHKK(1,NHKK)=                VHKK(1,NHKK-2)
            WHKK(2,NHKK)=                VHKK(2,NHKK-2)
            WHKK(3,NHKK)=                VHKK(3,NHKK-2)
            WHKK(4,NHKK)=VHKK(4,NHKK)*GAMP-VHKK(3,NHKK)*BGAMP
            IF (IPHKK.GE.2) WRITE(6,5001)
     *      IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),JMOHKK(2,IHKK),
     &      JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK(KHKK,IHKK),KHKK=1,5),
     &      (WHKK(KHKK,IHKK),KHKK=1,4)
          ENDIF
          IF (IPHKK.GE.1)THEN
	  WRITE(6,'(A)')' KKEVHH:'
	  WRITE(6,5001)
     *      IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),JMOHKK(2,IHKK),
     &      JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK(KHKK,IHKK),KHKK=1,5),
     &      (VHKK(KHKK,IHKK),KHKK=1,4)
	  ENDIF
C
C
C                                     CHAIN 2 PROJECTILE SEA-QUARK
          IIJJKK=0
          IF(IIJJKK.EQ.0)GO TO 33446	  
          NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
          IHKK=NHKK
          ISTHKK(IHKK)=151
          IDHKK(IHKK)=IHKKQ(IJJAQ1(N))
          JMOHKK(1,IHKK)=IHKKPO
          JMOHKK(2,IHKK)=IHKKPO
          JDAHKK(1,IHKK)=IHKK+2
          JDAHKK(2,IHKK)=IHKK+2
          PHKK(1,IHKK)=PJETB1(N,1)
          PHKK(2,IHKK)=PJETB1(N,2)
          PHKK(3,IHKK)=PJETB1(N,3)
          PHKK(4,IHKK)=PJETB1(N,4)
          PHKK(5,IHKK)=0.
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
          VHKK(1,IHKK)=VHKK(1,IHKKPO)+XXPP
          VHKK(2,IHKK)=VHKK(2,IHKKPO)+YYPP
          VHKK(3,IHKK)=VHKK(3,IHKKPO)
          VHKK(4,IHKK)=VHKK(4,IHKKPO)
          IF (IPHKK.GE.2) WRITE(6,5001)
     *      IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),JMOHKK(2,IHKK),
     &      JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK(KHKK,IHKK),KHKK=1,5),
     &      (VHKK(KHKK,IHKK),KHKK=1,4)
C                                     CHAIN 2 TARGET SEA-QUARK
          NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
          IHKK=NHKK
          ISTHKK(IHKK)=152
          IDHKK(IHKK)=IHKKQ(IJJQ2(N))
          JMOHKK(1,IHKK)=IHKKTO
          JMOHKK(2,IHKK)=IHKKTO
          JDAHKK(1,IHKK)=IHKK+1
          JDAHKK(2,IHKK)=IHKK+1
          PHKK(1,IHKK)=PJETB2(N,1)
          PHKK(2,IHKK)=PJETB2(N,2)
          PHKK(3,IHKK)=PJETB2(N,3)
          PHKK(4,IHKK)=PJETB2(N,4)
          PHKK(5,IHKK)=0.
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
          VHKK(1,IHKK)=VHKK(1,IHKKTO)+XXPP
          VHKK(2,IHKK)=VHKK(2,IHKKTO)+YYPP
          VHKK(3,IHKK)=VHKK(3,IHKKTO)
          VHKK(4,IHKK)=VHKK(4,IHKKTO)
          IF (IPHKK.GE.2) WRITE(6,5001)
     *      IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),JMOHKK(2,IHKK),
     &      JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK(KHKK,IHKK),KHKK=1,5),
     &      (VHKK(KHKK,IHKK),KHKK=1,4)
C
C                                     CHAIN 2 BEFORE FRAGMENTATION
          NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
          IHKK=NHKK
          ISTHKK(IHKK)=7
          IDHKK(IHKK)=88888
          JMOHKK(1,IHKK)=IHKK-2
          JMOHKK(2,IHKK)=IHKK-1
          PHKK(1,IHKK)=PHKK(1,IHKK-2)+PHKK(1,IHKK-1)
          PHKK(2,IHKK)=PHKK(2,IHKK-2)+PHKK(2,IHKK-1)
          PHKK(3,IHKK)=PHKK(3,IHKK-2)+PHKK(3,IHKK-1)
          PHKK(4,IHKK)=PHKK(4,IHKK-2)+PHKK(4,IHKK-1)
          PHKK(5,IHKK)=AMJCH2(N)
C                             POSITION OF CREATED CHAIN IN LAB
C                             =POSITION OF TARGET NUCLEON
C                             TIME OF CHAIN CREATION IN LAB
C                             =TIME OF PASSAGE OF PROJECTILE
C                              NUCLEUS AT POSITION OF TAR. NUCLEUS
          VHKK(1,NHKK)=                VHKK(1,NHKK-1)
          VHKK(2,NHKK)=                VHKK(2,NHKK-1)
          VHKK(3,NHKK)=                VHKK(3,NHKK-1)
          VHKK(4,NHKK)=VHKK(3,NHKK)/BETP-VHKK(3,NHKK-2)/BGAMP
          MHKKHH(N)=IHKK
          IF (IPROJK.EQ.1)THEN
            WHKK(1,NHKK)=                VHKK(1,NHKK-2)
            WHKK(2,NHKK)=                VHKK(2,NHKK-2)
            WHKK(3,NHKK)=                VHKK(3,NHKK-2)
            WHKK(4,NHKK)=VHKK(4,NHKK)*GAMP-VHKK(3,NHKK)*BGAMP
            IF (IPHKK.GE.2) THEN
	  WRITE(6,'(A)')' KKEVHH:'
	    WRITE(6,5001)
     *      IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),JMOHKK(2,IHKK),
     &      JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK(KHKK,IHKK),KHKK=1,5),
     &      (WHKK(KHKK,IHKK),KHKK=1,4)
	    ENDIF
          ENDIF
          IF (IPHKK.GE.2) THEN
	  WRITE(6,'(A)')' KKEVHH:'
	  WRITE(6,5001)
     *      IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),JMOHKK(2,IHKK),
     &      JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK(KHKK,IHKK),KHKK=1,5),
     &      (VHKK(KHKK,IHKK),KHKK=1,4)
	  ENDIF
33446     CONTINUE	  
C
C  NOW WE HAVE AN ACCEPTABLE HARD  EVENT
C  AND PUT IT INTO THE HISTOGRAM
C
        AMCHH1(N)=AMJCH1(N)
        AMCHH2(N)=AMJCH2(N)
        GACHH1(N)=GAMJH1(N)
        BGXHH1(N)=BGXJH1(N)
        BGYHH1(N)=BGYJH1(N)
        BGZHH1(N)=BGZJH1(N)
        GACHH2(N)=GAMJH2(N)
        BGXHH2(N)=BGXJH2(N)
        BGYHH2(N)=BGYJH2(N)
        BGZHH2(N)=BGZJH2(N)
	NNCH1=0
	NNCH2=0
	IJNCH1=0
	IJNCH2=0
        NCHHH1(N)=NNCH1
        NCHHH2(N)=NNCH2
        IJCHH1(N)=IJNCH1
        IJCHH2(N)=IJNCH2
        DO 1234 III=1,4
          PQHHA1(N,III)=PJETA1(N,III)
          PQHHA2(N,III)=PJETA2(N,III)
          PQHHB1(N,III)=PJETB1(N,III)
          PQHHB2(N,III)=PJETB2(N,III)
 1234   CONTINUE
        IF (IPEV.GE.6)WRITE(6,104)N,
     *               AMCHH1(N),AMCHH2(N),GACHH1(N),GACHH2(N),
     *               BGXHH1(N),BGYHH1(N),BGZHH1(N),
     *               BGXHH2(N),BGYHH2(N),BGZHH2(N),
     *               NCHHH1(N),NCHHH2(N),IJCHH1(N),IJCHH2(N)
      ENDIF
  101 CONTINUE
C
  104 FORMAT(' HH - 104',
     *       I10,4F12.7    /10X,6F12.6,4I5)
  211 FORMAT (' HH: AMMM,GAMMM,BGGGX,BGGGY,BGGGZ,IREJ ',5F12.5,I10/
     *        '     AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',5F12.5/
     *        '     AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ',5F12.5)
  212 FORMAT (' HH: AMMM,GAMMM,BGGGX,BGGGY,BGGGZ,IREJ || ',5F12.5,I10/
     *        '     AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',5F12.5/
     *        '     AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ',5F12.5)
 8001 FORMAT(' KKEVHH - IRHH13=',I5)
 8002 FORMAT( ' HH - 8002',5E12.4/4(4E12.4/),2E12.4/2I5/4E12.4)
 8003 FORMAT(' KKEVHH - IRHH11=',I5)
 8005 FORMAT(' KKEVHH - IRHH12=',I5)
 8006 FORMAT(' HH - 8006', 5I5/2(4E12.4/),2E12.4)
      RETURN
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE KKEVZZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C---------------    TREATMENT OF SUPPLEMENTARY SEA CHAIN SYSTEMS
C
*KEEP,NUCC.
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
*KEEP,RPTSHM.
      COMMON /RPTSHM/ RPROJ,RTARG,BIMPAC
      PARAMETER (INTMX=2488,INTMD=252)
      COMMON /TRAFOP/ GAMP,BGAMP,BETP
      COMMON /DXQX/   XPVQ(248),XPVD(248),XTVQ(248),XTVD(248),
     *                XPSQ(INTMX),XPSAQ(INTMX),XTSQ(INTMX),XTSAQ(INTMX)
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
C--------------------
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
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/ ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),ZUOST(INTMX),
     *               INTLO(INTMX),INLOSS(INTMX)
C-------------------
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
C--------------------
      COMMON /NUCIMP/ PRMOM(5,248),TAMOM(5,248),
     &                PRMFEP,PRMFEN,TAMFEP,TAMFEN,
     &                PREFEP,PREFEN,TAEFEP,TAEFEN,
     &        PREPOT(210),TAEPOT(210),PREBIN,TAEBIN,FERMOD,ETACOU
C--------------------
C--------------------
      LOGICAL INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA,
     *        IPADIS,ISHMAL,LPAULI
      COMMON /DROPPT/ INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA,
     *                IPADIS,ISHMAL,LPAULI
C-------------------
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
      COMMON /REJEC/ IRCO1,IRCO2,IRCO3,IRCO4,IRCO5,IRSS11,IRSS12,IRSS13,
     *               IRSS14,
     *               IRSV11,IRSV12,IRSV13,IRSV14,
     *               IRVS11,IRVS12,IRVS13,IRVS14,
     *               IRVV11,IRVV12,IRVV13,IRVV14
      COMMON /ABRZZ/ AMCZZ1(INTMX),AMCZZ2(INTMX),
     *               GACZZ1(INTMX),GACZZ2(INTMX),
     *               BGXZZ1(INTMX),BGYZZ1(INTMX),BGZZZ1(INTMX),
     *               BGXZZ2(INTMX),BGYZZ2(INTMX),BGZZZ2(INTMX),
     *               NCHZZ1(INTMX),NCHZZ2(INTMX),
     *               IJCZZ1(INTMX),IJCZZ2(INTMX),
     *               PQZZA1(INTMX,4),PQZZA2(INTMX,4),
     *               PQZZB1(INTMX,4),PQZZB2(INTMX,4)
C-------------------
      PARAMETER (NMXHKK=  89998)
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),
     & VHKK(4,NMXHKK),WHKK(4,NMXHKK)
C                       WHKK(4,NMXHKK) GIVES POSITIONS AND TIMES IN
C                       PROJECTILE FRAME, THE CHAINS ARE CREATED ON
C                       THE POSITIONS OF THE PROJECTILE NUCLEONS
C                       IN THE PROJECTILE FRAME (TARGET NUCLEONS IN
C                       TARGET FRAME) BOTH POSITIONS ARE THREFORE NOT
C                       COMPLETELY CONSISTENT. THE TIMES IN THE
C                       PROJECTILE FRAME HOWEVER ARE OBTAINED BY
C                       LORENTZ TRANSFORMING FROM THE LAB SYSTEM.
      COMMON /PROJK/ IPROJK
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
C
C-----------------------------------------------------------------------
      COMMON /ABRSOF/XSQ1(INTMX),XSAQ1(INTMX),XSQ2(INTMX),XSAQ2(INTMX),
     *        IJSQ1(INTMX),IJSAQ1(INTMX),IJSQ2(INTMX),IJSAQ2(INTMX),
     *        AMCCH1(INTMX),AMCCH2(INTMX),GAMCH1(INTMX),GAMCH2(INTMX),
     *        BGCH1(INTMX),BGCH2(INTMX),THECH1(INTMX),THECH2(INTMX),
     *        BGXCH1(INTMX),BGYCH1(INTMX),BGZCH1(INTMX),
     *        BGXCH2(INTMX),BGYCH2(INTMX),BGZCH2(INTMX),
     *        NCH1(INTMX),NCH2(INTMX),IJCH1(INTMX),IJCH2(INTMX),
     *  PSOFA1(INTMX,4),PSOFA2(INTMX,4),PSOFB1(INTMX,4),PSOFB2(INTMX,4)
     * ,JHKKPZ(INTMX),JHKKTZ(INTMX),JHKKSX(INTMX),JHKKS1(INTMX)
*KEEP,ABRZD.
      COMMON /ABRZD/ AMCZD1(INTMD),AMCZD2(INTMD),
     +GACZD1(INTMD),GACZD2(INTMD),
     +BGXZD1(INTMD),BGYZD1(INTMD),BGZZD1(INTMD), 
     +BGXZD2(INTMD),BGYZD2(INTMD),
     +BGZZD2(INTMD), NCHZD1(INTMD),NCHZD2(INTMD),
     +IJCZD1(INTMD),IJCZD2(INTMD),
     +PQZDA1(INTMD,4),PQZDA2(INTMD,4), PQZDB1(INTMD,4),
     +PQZDB2(INTMD,4),
     +IPCQ(INTMD),ITCQ(INTMD),ITCQ2(INTMD),IPCAQ(INTMD),
     +ITCAQ(INTMD),ITCAQ2(INTMD)
     +,IZDSS(INTMD)
*KEEP,ABRDZ.
      COMMON /ABRDZ/ AMCDZ1(INTMD),AMCDZ2(INTMD),
     +GACDZ1(INTMD),GACDZ2(INTMD),
     +BGXDZ1(INTMD),BGYDZ1(INTMD),BGZDZ1(INTMD), 
     +BGXDZ2(INTMD),BGYDZ2(INTMD),
     +BGZDZ2(INTMD), NCHDZ1(INTMD),NCHDZ2(INTMD),
     +IJCDZ1(INTMD),IJCDZ2(INTMD),
     +PQDZA1(INTMD,4),PQDZA2(INTMD,4), PQDZB1(INTMD,4),
     +PQDZB2(INTMD,4),
     +IPZQ(INTMD),IPZQQ2(INTMD),ITZQ(INTMD),IPZAQ(INTMD),
     +IZAQQ2(INTMD),ITZAQ(INTMD)
     +,IDZSS(INTMD)
C-------------------
      COMMON /NUCJTN/NONUJ1,NONUJT,NONUS1,NONUST
      COMMON /XSVTHR/ XSTHR,XVTHR,XDTHR,XSSTHR
      COMMON/INTNEZ/NDZ,NZD
C----------------------------------------------------------------------
      DIMENSION IHKKQ(-6:6)
      DATA IHKKQ/-6,-5,-4,-3,-1,-2,0,2,1,3,4,5,6/
C
C----------------------------------------------------------------------
      DO 1001 N=1,NONUST
         DO 1002 I=1,NDZ 
           IF(IDZSS(I).EQ.N.AND.NCH1(N).EQ.99)THEN
             NCHDZ1(I)=99
             NCHDZ2(I)=99
           ENDIF
           IF(IDZSS(I).EQ.N.AND.JHKKSX(N).NE.1)THEN
             NCHDZ1(I)=99
             NCHDZ2(I)=99
           ENDIF
           IF(IPEV.EQ.2)THEN
             WRITE(6,'(A,6I10)')' kkevzz:n,i,ndz,nchdz1,jhkksx,idzss'
     *       ,N,I,NDZ,NCHDZ1(I),JHKKSX(N),IDZSS(I)
           ENDIF
 1002    CONTINUE
         DO 1003 I=1,NZD 
           IF(IZDSS(I).EQ.N.AND.NCH1(N).EQ.99)THEN
             NCHZD1(I)=99
             NCHZD2(I)=99
           ENDIF
           IF(IZDSS(I).EQ.N.AND.JHKKSX(N).NE.1)THEN
             NCHZD1(I)=99
             NCHZD2(I)=99
           ENDIF
           IF(IPEV.EQ.2)THEN
             WRITE(6,'(A,6I10)')' kkevzz:n,i,nzd,nchzd1,jhkksx,izdss'
     *       ,N,I,NZD,NCHZD1(I),JHKKSX(N),IZDSS(I)
           ENDIF
 1003    CONTINUE
 1001 CONTINUE
      DO 101 N=1,NONUST
        IF(NCH1(N).EQ.88)GO TO 101
        IF(NCH2(N).EQ.88)GO TO 101
        IF (JHKKSX(N).EQ.1)THEN
          IXVPR=JHKKPZ(N)
          IXVTA=JHKKTZ(N)
          IHKKPO=JHKKPV(IXVPR)
          IHKKTO=JHKKTV(IXVTA)
C
C                               Cronin multiple scattering
C     IF(IT.GT.1.AND.N.GT.100)THEN
      IF(IT.GT.1)THEN
      ITNU=IHKKTO
      RTIX=(VHKK(1,ITNU))*1.E12-BIMPAC*0.1
      RTIY=VHKK(2,ITNU)*1.E12
      RTIZ=VHKK(3,ITNU)*1.E12
      RTIR2=(RTIX**2+RTIY**2+RTIZ**2)
      IF(RTIR2.GT.RTARG**2)THEN
        IF(IPEV.GE.2)
     *  WRITE(6,774)RTARG,RTIX,RTIY,RTIZ,BIMPAC,IHKKTO,IXVTA
 774    FORMAT(' KKEVZZ: RTARG,RTIX,RTIY,RTIZ,BIMPAC,IHKKTO,IXVTA'
     *  ,5E12.4,2I10)
        GO TO 779
      ENDIF
      IF(NCH1(N).EQ.0)THEN
      PVQPX=PSOFA1(N,1)
      PVQPY=PSOFA1(N,2)
      PVQPZ=PSOFA1(N,3)
      PVQE =PSOFA1(N,4)
      IF(PVQE.LE.0.D0)THEN
        PVQEN=SQRT(PVQPX**2+PVQPY**2+PVQPZ**2)
        WRITE(6,776)PVQE,PVQEN,N,NONUST
 776    FORMAT(' KKEVZZ: PVQE,PVQEN,N,NONUST ',2E12.4,2I5)
        PVQE=PVQEN
      ENDIF
      CALL CROMSC(PVQPX,PVQPY,PVQPZ,PVQE,RTIX,RTIY,RTIZ,
     *            PVQNX,PVQNY,PVQNZ,PVQNE,30)
      PVDQPX=PSOFA2(N,1)
      PVDQPY=PSOFA2(N,2)
      PVDQPZ=PSOFA2(N,3)
      PVDQE =PSOFA2(N,4)
      IF(PVDQE.LE.0.D0)THEN
        PVDQEN=SQRT(PVDQPX**2+PVDQPY**2+PVDQPZ**2)
        WRITE(6,778)PVDQE,PVDQEN,N,NONUST
 778    FORMAT(' KKEVZZ: PVDQE,PVDQEN,N,NONUST ',2E12.4,2I5)
        PVDQE=PVDQEN 
      ENDIF
      CALL CROMSC(PVDQPX,PVDQPY,PVDQPZ,PVDQE,RTIX,RTIY,RTIZ,
     *            PVDQNX,PVDQNY,PVDQNZ,PVDQNE,31)
      AMTES2=((PVQNE+PVDQNE)**2-(PVQNX+PVDQNX)**2
     *       -(PVQNY+PVDQNY)**2-(PVQNZ+PVDQNZ)**2)
      IF(AMTES2.GE.AMCCH1(N)**2.OR.AMTES2.GE.25.D0)THEN 
      PSOFA1(N,1)=PVQNX
      PSOFA1(N,2)=PVQNY
      PSOFA1(N,3)=PVQNZ
      PSOFA1(N,4)=PVQNE
      PSOFA2(N,1)=PVDQNX
      PSOFA2(N,2)=PVDQNY
      PSOFA2(N,3)=PVDQNZ
      PSOFA2(N,4)=PVDQNE
      ENDIF
C                                  MASSES OF SUBCHAINS
          XMCCH1=SQRT((PSOFA1(N,4)+
     *                         PSOFA2(N,4))**2
     *                       -(PSOFA1(N,1)+
     *                         PSOFA2(N,1))**2
     *                       -(PSOFA1(N,2)+
     *                         PSOFA2(N,2))**2
     *                       -(PSOFA1(N,3)+
     *                         PSOFA2(N,3))**2)
         IF(XMCCH1.GE.AMCCH1(N))THEN
          AMCCH1(N)=XMCCH1
C
          GAMCH1(N)=(PSOFA1(N,4)+
     *                    PSOFA2(N,4))/AMCCH1(N)
          BGXCH1(N)=(PSOFA1(N,1)+
     *                    PSOFA2(N,1))/AMCCH1(N)
          BGYCH1(N)=(PSOFA1(N,2)+
     *                    PSOFA2(N,2))/AMCCH1(N)
          BGZCH1(N)=(PSOFA1(N,3)+
     *                    PSOFA2(N,3))/AMCCH1(N)
      ENDIF
      ENDIF
      IF(NCH2(N).EQ.0)THEN
      PVQTX=PSOFB1(N,1)
      PVQTY=PSOFB1(N,2)
      PVQTZ=PSOFB1(N,3)
      PVQTE=PSOFB1(N,4)
      IF(PVQTE.LE.0.D0)THEN
        PVQTEN=SQRT(PVQTX**2+PVQTY**2+PVQTZ**2)
        WRITE(6,786)PVQTE,PVQTEN,N,NONUST
 786    FORMAT(' KKEVZZ: PVQTE,PVQTEN,N,NONUST ',2E12.4,2I5)
        PVQTE=PVQTEN
      ENDIF
      CALL CROMSC(PVQTX,PVQTY,PVQTZ,PVQTE,RTIX,RTIY,RTIZ,
     *            PVQNTX,PVQNTY,PVQNTZ,PVQNTE,32)
      PVDQTX=PSOFB2(N,1)
      PVDQTY=PSOFB2(N,2)
      PVDQTZ=PSOFB2(N,3)
      PVDQTE=PSOFB2(N,4)
      IF(PVDQTE.LE.0.D0)THEN
        PVDTEN=SQRT(PVDQTX**2+PVDQTY**2+PVDQTZ**2)
        WRITE(6,796)PVDQTE,PVDTEN,N,NONUST
 796    FORMAT(' KKEVZZ: PVQTE,PVQTEN,N,NONUST ',2E12.4,2I5)
        PVDQTE=PVDTEN
      ENDIF
      CALL CROMSC(PVDQTX,PVDQTY,PVDQTZ,PVDQTE,RTIX,RTIY,RTIZ,
     *            PVTQNX,PVTQNY,PVTQNZ,PVTQNE,33)
      AMTES2=((PVQNTE+PVTQNE)**2-(PVQNTX+PVTQNX)**2
     *       -(PVQNTY+PVTQNY)**2-(PVQNTZ+PVTQNZ)**2)
      IF(AMTES2.GE.AMCCH1(N)**2.OR.AMTES2.GE.25.D0)THEN 
      PSOFB1(N,1)=PVQNTX
      PSOFB1(N,2)=PVQNTY
      PSOFB1(N,3)=PVQNTZ
      PSOFB1(N,4)=PVQNTE
      PSOFB2(N,1)=PVTQNX
      PSOFB2(N,2)=PVTQNY
      PSOFB2(N,3)=PVTQNZ
      PSOFB2(N,4)=PVTQNE
      ENDIF
C                                  MASSES OF SUBCHAINS
          XMCCH2=SQRT((PSOFB1(N,4)+
     *                         PSOFB2(N,4))**2
     *                       -(PSOFB1(N,1)+
     *                         PSOFB2(N,1))**2
     *                       -(PSOFB1(N,2)+
     *                         PSOFB2(N,2))**2
     *                       -(PSOFB1(N,3)+
     *                         PSOFB2(N,3))**2)
         IF(XMCCH2.GE.AMCCH2(N))THEN
          AMCCH2(N)=XMCCH2
C
          GAMCH2(N)=(PSOFB1(N,4)+
     *                    PSOFB2(N,4))/AMCCH2(N)
          BGXCH2(N)=(PSOFB1(N,1)+
     *                    PSOFB2(N,1))/AMCCH2(N)
          BGYCH2(N)=(PSOFB1(N,2)+
     *                    PSOFB2(N,2))/AMCCH2(N)
          BGZCH2(N)=(PSOFB1(N,3)+
     *                    PSOFB2(N,3))/AMCCH2(N)
      ENDIF
      ENDIF
      ENDIF
 779  CONTINUE
C
C
C                                      PUT Z-Z CHAIN ENDS INTO /HKKEVT/
C                                      MOMENTA IN NN-CMS
C                                      POSITION OF ORIGINAL NUCLEONS
c                                      flags for Z-Z chain ends
c                                        projectile: isthkk=251
c                                        target:     isthkk=252
c                                        h-h chains: isthkk=9
C
          IXVPR=JHKKPZ(N)
          IXVTA=JHKKTZ(N)
          IHKKPO=JHKKPV(IXVPR)
          IHKKTO=JHKKTV(IXVTA)
          IF (IPEV.GT.3)WRITE(6,5002)IXVPR,IHKKPO
 5002     FORMAT (' IXVPR,IHKKPO ',5I5)
          IF (IPEV.GT.3)WRITE(6,5003)IXVTA,IHKKTO
 5003     FORMAT (' IXVTA,IHKKTO ',5I5)
C                                     CHAIN 1 PROJECTILE SEA-QUARK
          NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
          IHKK=NHKK
          ISTHKK(IHKK)=251
          IDHKK(IHKK)=IHKKQ(IJSQ1(N))
          JMOHKK(1,IHKK)=IHKKPO
          JMOHKK(2,IHKK)=IHKKPO
          JDAHKK(1,IHKK)=IHKK+2
          JDAHKK(2,IHKK)=IHKK+2
          PHKK(1,IHKK)=PSOFA1(N,1)
          PHKK(2,IHKK)=PSOFA1(N,2)
          PHKK(3,IHKK)=PSOFA1(N,3)
          PHKK(4,IHKK)=PSOFA1(N,4)
          PHKK(5,IHKK)=0.
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
          VHKK(1,IHKK)=VHKK(1,IHKKPO)+XXPP
          VHKK(2,IHKK)=VHKK(2,IHKKPO)+YYPP
          VHKK(3,IHKK)=VHKK(3,IHKKPO)
          VHKK(4,IHKK)=VHKK(4,IHKKPO)
          IF (IPHKK.GE.2) WRITE(6,5001)
     *      IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),JMOHKK(2,IHKK),
     &      JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK(KHKK,IHKK),KHKK=1,5),
     &      (VHKK(KHKK,IHKK),KHKK=1,4)
 5001     FORMAT (I6,I4,5I6,9E10.2)
C                                     CHAIN 1 TARGET SEA-QUARK
          NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
          IHKK=NHKK
          ISTHKK(IHKK)=252
          IDHKK(IHKK)=IHKKQ(IJSAQ2(N))
          JMOHKK(1,IHKK)=IHKKTO
          JMOHKK(2,IHKK)=IHKKTO
          JDAHKK(1,IHKK)=IHKK+1
          JDAHKK(2,IHKK)=IHKK+1
          PHKK(1,IHKK)=PSOFA2(N,1)
          PHKK(2,IHKK)=PSOFA2(N,2)
          PHKK(3,IHKK)=PSOFA2(N,3)
          PHKK(4,IHKK)=PSOFA2(N,4)
          PHKK(5,IHKK)=0.
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
          VHKK(1,IHKK)=VHKK(1,IHKKTO)+XXPP
          VHKK(2,IHKK)=VHKK(2,IHKKTO)+YYPP
          VHKK(3,IHKK)=VHKK(3,IHKKTO)
          VHKK(4,IHKK)=VHKK(4,IHKKTO)
          IF (IPHKK.GE.2) WRITE(6,5001)
     *      IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),JMOHKK(2,IHKK),
     &      JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK(KHKK,IHKK),KHKK=1,5),
     &      (VHKK(KHKK,IHKK),KHKK=1,4)
C
C                                     CHAIN 1 BEFORE FRAGMENTATION
          NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
          IHKK=NHKK
          ISTHKK(IHKK)=9
          IDHKK(IHKK)=88888+NCH1(N)
          JMOHKK(1,IHKK)=IHKK-2
          JMOHKK(2,IHKK)=IHKK-1
          PHKK(1,IHKK)=PHKK(1,IHKK-2)+PHKK(1,IHKK-1)
          PHKK(2,IHKK)=PHKK(2,IHKK-2)+PHKK(2,IHKK-1)
          PHKK(3,IHKK)=PHKK(3,IHKK-2)+PHKK(3,IHKK-1)
          PHKK(4,IHKK)=PHKK(4,IHKK-2)+PHKK(4,IHKK-1)
          PHKK(5,IHKK)=AMCCH1(N)
C                             POSITION OF CREATED CHAIN IN LAB
C                             =POSITION OF TARGET NUCLEON
C                             TIME OF CHAIN CREATION IN LAB
C                             =TIME OF PASSAGE OF PROJECTILE
C                              NUCLEUS AT POSITION OF TAR. NUCLEUS
          VHKK(1,NHKK)=                VHKK(1,NHKK-1)
          VHKK(2,NHKK)=                VHKK(2,NHKK-1)
          VHKK(3,NHKK)=                VHKK(3,NHKK-1)
          VHKK(4,NHKK)=0.
          MHKKHH(N)=IHKK
          IF (IPROJK.EQ.1)THEN
            WHKK(1,NHKK)=                VHKK(1,NHKK-2)
            WHKK(2,NHKK)=                VHKK(2,NHKK-2)
            WHKK(3,NHKK)=                VHKK(3,NHKK-2)
            WHKK(4,NHKK)=VHKK(4,NHKK)*GAMP-VHKK(3,NHKK)*BGAMP
            IF (IPHKK.GE.2) WRITE(6,5001)
     *      IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),JMOHKK(2,IHKK),
     &      JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK(KHKK,IHKK),KHKK=1,5),
     &      (WHKK(KHKK,IHKK),KHKK=1,4)
          ENDIF
          IF (IPHKK.GE.2) WRITE(6,5001)
     *      IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),JMOHKK(2,IHKK),
     &      JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK(KHKK,IHKK),KHKK=1,5),
     &      (VHKK(KHKK,IHKK),KHKK=1,4)
C
C
C                                     CHAIN 2 PROJECTILE SEA-QUARK
          NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
          IHKK=NHKK
          ISTHKK(IHKK)=251
          IDHKK(IHKK)=IHKKQ(IJSAQ1(N))
          JMOHKK(1,IHKK)=IHKKPO
          JMOHKK(2,IHKK)=IHKKPO
          JDAHKK(1,IHKK)=IHKK+2
          JDAHKK(2,IHKK)=IHKK+2
          PHKK(1,IHKK)=PSOFB1(N,1)
          PHKK(2,IHKK)=PSOFB1(N,2)
          PHKK(3,IHKK)=PSOFB1(N,3)
          PHKK(4,IHKK)=PSOFB1(N,4)
          PHKK(5,IHKK)=0.
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
          VHKK(1,IHKK)=VHKK(1,IHKKPO)+XXPP
          VHKK(2,IHKK)=VHKK(2,IHKKPO)+YYPP
          VHKK(3,IHKK)=VHKK(3,IHKKPO)
          VHKK(4,IHKK)=VHKK(4,IHKKPO)
          IF (IPHKK.GE.2) WRITE(6,5001)
     *      IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),JMOHKK(2,IHKK),
     &      JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK(KHKK,IHKK),KHKK=1,5),
     &      (VHKK(KHKK,IHKK),KHKK=1,4)
C                                     CHAIN 2 TARGET SEA-QUARK
          NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
          IHKK=NHKK
          ISTHKK(IHKK)=252
          IDHKK(IHKK)=IHKKQ(IJSQ2(N))
          JMOHKK(1,IHKK)=IHKKTO
          JMOHKK(2,IHKK)=IHKKTO
          JDAHKK(1,IHKK)=IHKK+1
          JDAHKK(2,IHKK)=IHKK+1
          PHKK(1,IHKK)=PSOFB2(N,1)
          PHKK(2,IHKK)=PSOFB2(N,2)
          PHKK(3,IHKK)=PSOFB2(N,3)
          PHKK(4,IHKK)=PSOFB2(N,4)
          PHKK(5,IHKK)=0.
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
          VHKK(1,IHKK)=VHKK(1,IHKKTO)+XXPP
          VHKK(2,IHKK)=VHKK(2,IHKKTO)+YYPP
          VHKK(3,IHKK)=VHKK(3,IHKKTO)
          VHKK(4,IHKK)=VHKK(4,IHKKTO)
          IF (IPHKK.GE.2) WRITE(6,5001)
     *      IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),JMOHKK(2,IHKK),
     &      JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK(KHKK,IHKK),KHKK=1,5),
     &      (VHKK(KHKK,IHKK),KHKK=1,4)
C
C                                     CHAIN 2 BEFORE FRAGMENTATION
          NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
          IHKK=NHKK
          ISTHKK(IHKK)=9
          IDHKK(IHKK)=88888+NCH2(N)
          JMOHKK(1,IHKK)=IHKK-2
          JMOHKK(2,IHKK)=IHKK-1
          PHKK(1,IHKK)=PHKK(1,IHKK-2)+PHKK(1,IHKK-1)
          PHKK(2,IHKK)=PHKK(2,IHKK-2)+PHKK(2,IHKK-1)
          PHKK(3,IHKK)=PHKK(3,IHKK-2)+PHKK(3,IHKK-1)
          PHKK(4,IHKK)=PHKK(4,IHKK-2)+PHKK(4,IHKK-1)
          PHKK(5,IHKK)=AMCCH2(N)
C                             POSITION OF CREATED CHAIN IN LAB
C                             =POSITION OF TARGET NUCLEON
C                             TIME OF CHAIN CREATION IN LAB
C                             =TIME OF PASSAGE OF PROJECTILE
C                              NUCLEUS AT POSITION OF TAR. NUCLEUS
          VHKK(1,NHKK)=                VHKK(1,NHKK-1)
          VHKK(2,NHKK)=                VHKK(2,NHKK-1)
          VHKK(3,NHKK)=                VHKK(3,NHKK-1)
          VHKK(4,NHKK)=VHKK(3,NHKK)/BETP-VHKK(3,NHKK-2)/BGAMP
          MHKKHH(N)=IHKK
          IF (IPROJK.EQ.1)THEN
            WHKK(1,NHKK)=                VHKK(1,NHKK-2)
            WHKK(2,NHKK)=                VHKK(2,NHKK-2)
            WHKK(3,NHKK)=                VHKK(3,NHKK-2)
            WHKK(4,NHKK)=VHKK(4,NHKK)*GAMP-VHKK(3,NHKK)*BGAMP
            IF (IPHKK.GE.2) WRITE(6,5001)
     *      IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),JMOHKK(2,IHKK),
     &      JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK(KHKK,IHKK),KHKK=1,5),
     &      (WHKK(KHKK,IHKK),KHKK=1,4)
          ENDIF
          IF (IPHKK.GE.2) WRITE(6,5001)
     *      IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),JMOHKK(2,IHKK),
     &      JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK(KHKK,IHKK),KHKK=1,5),
     &      (VHKK(KHKK,IHKK),KHKK=1,4)
C
C  NOW WE HAVE AN ACCEPTABLE HARD  EVENT
C  AND PUT IT INTO THE HISTOGRAM
C
        AMCZZ1(N)=AMCCH1(N)
        AMCZZ2(N)=AMCCH2(N)
        GACZZ1(N)=GAMCH1(N)
        BGXZZ1(N)=BGXCH1(N)
        BGYZZ1(N)=BGYCH1(N)
        BGZZZ1(N)=BGZCH1(N)
        GACZZ2(N)=GAMCH2(N)
        BGXZZ2(N)=BGXCH2(N)
        BGYZZ2(N)=BGYCH2(N)
        BGZZZ2(N)=BGZCH2(N)
        NCHZZ1(N)=NCH1(N)
        NCHZZ2(N)=NCH2(N)
        IJCZZ1(N)=IJCH1(N)
        IJCZZ2(N)=IJCH2(N)
        DO 1234 III=1,4
          PQZZA1(N,III)=PSOFA1(N,III)
          PQZZA2(N,III)=PSOFA2(N,III)
          PQZZB1(N,III)=PSOFB1(N,III)
          PQZZB2(N,III)=PSOFB2(N,III)
 1234   CONTINUE
        IF (IPEV.GE.6)WRITE(6,104)N,
     *               AMCZZ1(N),AMCZZ2(N),GACZZ1(N),GACZZ2(N),
     *               BGXZZ1(N),BGYZZ1(N),BGZZZ1(N),
     *               BGXZZ2(N),BGYZZ2(N),BGZZZ2(N),
     *               NCHZZ1(N),NCHZZ2(N),IJCZZ1(N),IJCZZ2(N)
      ENDIF
  101 CONTINUE
C
  104 FORMAT(' ZZ - 104',
     *       I10,4F12.7    /10X,6F12.6,           4I5              )
  211 FORMAT (' ZZ: AMMM,GAMMM,BGGGX,BGGGY,BGGGZ,IREJ ',5F12.5,I10/
     *        '     AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',5F12.5/
     *        '     AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ',5F12.5)
  212 FORMAT (' ZZ: AMMM,GAMMM,BGGGX,BGGGY,BGGGZ,IREJ || ',5F12.5,I10/
     *        '     AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',5F12.5/
     *        '     AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ',5F12.5)
 8001 FORMAT(' KKEVZZ - IRZZ13=',I5)
 8002 FORMAT( ' ZZ - 8002',5E12.4/4(4E12.4/),2E12.4/2I5/4E12.4)
 8003 FORMAT(' KKEVZZ - IRZZ11=',I5)
 8005 FORMAT(' KKEVZZ - IRZZ12=',I5)
 8006 FORMAT(' ZZ - 8006', 5I5/2(4E12.4/),2E12.4)
      RETURN
      END
C------------------------------------------------------------------------
      SUBROUTINE CORRCO
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
*KEEP,HKKEVT.
C     INCLUDE (HKKEVT)
      PARAMETER (NMXHKK= 89998)
C     PARAMETER (NMXHKK=25000)
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

      COMMON /EXTEVT/ IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),
     +                IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10)

      DO 1 I=1,NHKK
        IF (IDHKK(I).EQ.88888) THEN
          M1=I-2
          M2=I-1
          IF (JMOHKK(2,M1).EQ.0) THEN
            JM1=JMOHKK(1,M1)
            JMOHKK(2,M1)=JMOHKK(1,JM1)
          ENDIF
          IF (JMOHKK(2,M2).EQ.0) THEN
            JM2=JMOHKK(1,M2)
            JMOHKK(2,M2)=JMOHKK(1,JM2)
          ENDIF
        ENDIF
   1  CONTINUE

      DO 2 I=1,NHKK
        IF (IDHKK(I).EQ.88888) THEN
          M1=I-2
          M2=I-1
          M2M1=JMOHKK(2,M1)
          M2M2=JMOHKK(2,M2)
          IF (JDAHKK(1,M2M1).EQ.0) THEN
            JDAHKK(1,M2M1)=M1
          ELSE
            IF (JDAHKK(2,M2M1).EQ.0) THEN
            JDAHKK(2,M2M1)=M1
            ENDIF
          ENDIF
          IF (JDAHKK(1,M2M2).EQ.0) THEN
            JDAHKK(1,M2M2)=M2
          ELSE
            IF (JDAHKK(2,M2M2).EQ.0) THEN
            JDAHKK(2,M2M2)=M2
            ENDIF
          ENDIF
        ENDIF
        MO1=JMOHKK(1,M1)
        MO2=JMOHKK(1,M2)
        IF(JDAHKK(1,MO1).EQ.0)JDAHKK(1,MO1)=M1
	IF(JDAHKK(2,MO1).EQ.0)JDAHKK(2,MO1)=M1
        IF(JDAHKK(1,MO2).EQ.0)JDAHKK(1,MO2)=M2
        IF(JDAHKK(2,MO2).EQ.0)JDAHKK(2,MO2)=M2
   2  CONTINUE

      DO 3 I=1,NHKK
        IF (ISTHKK(I).EQ.11) THEN
          IF ((JDAHKK(1,I).EQ.0).AND.(JDAHKK(2,I).EQ.0)) THEN
            ISTHKK(I)=13
          ENDIF
        ENDIF
        IF (ISTHKK(I).EQ.12) THEN
          IF ((JDAHKK(1,I).EQ.0).AND.(JDAHKK(2,I).EQ.0)) THEN
          ISTHKK(I)=14
          ENDIF
        ENDIF
   3  CONTINUE

      RETURN
      END

C--------------------------------------------------------------      
      SUBROUTINE KKEVNU(NHKKH1,EPN,PPN,KKMAT,IREJ,ECM)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      COMMON/INTNEZ/NDZ,NZD
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
C      /
C INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C------------------
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
*KEEP,SHMAKL.
C     INCLUDE (SHMAKL)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON/SHMAKL/JSSH(INTMX),JTSH(INTMX),INTER1(INTMX),INTER2(INTMX)
*KEEP,NUCC.
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
*KEEP,RPTSHM.
      COMMON /RPTSHM/ RPROJ,RTARG,BIMPAC
*KEEP,NSHMAK.
      COMMON /NSHMAK/ NNSHMA,NPSHMA,NTSHMA,NSHMAC,NSHMA2
*KEEP,ZENTRA.
      COMMON /ZENTRA/ ICENTR
*KEEP,NUCIMP.
      COMMON /NUCIMP/ PRMOM(5,248),TAMOM(5,248), PRMFEP,PRMFEN,TAMFEP,
     +TAMFEN, PREFEP,PREFEN,TAEFEP,TAEFEN, PREPOT(210),TAEPOT(210),
     +PREBIN,TAEBIN,FERMOD,ETACOU
*KEEP,DROPPT.
      LOGICAL INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA, IPADIS,
     +ISHMAL,LPAULI
      COMMON /DROPPT/ INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA,
     +IPADIS,ISHMAL,LPAULI
*KEEP,NNCMS.
      COMMON /NNCMS/  GAMCM,BGCM,UMOJ,PCMJ,EPROJJ,PPROJJ
      COMMON /NUCCMS/ GACMS,BGCMS,GALAB,BGLAB,BLAB,UMO,PCM,EPROJ,PPROJ
*KEEP,NUCPOS.
      COMMON /NUCPOS/INVVP(248),INVVT(248),INVSP(248),INVST(248), NUVV,
     +NUVS,NUSV,NUSS,INSVP(248),INSVT(248),INSSP(248),INSST(248), ISVEAP
     +(248),ISVEAT(248),ISSEAP(248),ISSEAT(248), IVSEAP(248),IVSEAT
     +(248), ISLOSP(248),ISLOST(248),INOOP(248),INOOT(248),NUOO
*KEEP,TAUFO.
      COMMON /TAUFO/  TAUFOR,KTAUGE,ITAUVE,INCMOD
      COMMON /EVAPPP/IEVAP
      COMMON /NEUTYY/NEUTYP,NEUDEC
*KEEP,RTAR.
      COMMON /RTAR/ RTARNU
*KEEP,INNU.
      COMMON /INNU/INUDEC
*KEEP,HADTHR.
      COMMON /HADTHR/ EHADTH,INTHAD
*KEEP,DINPDA.
      COMMON /DINPDA/ IMPS(6,6),IMVE(6,6),IB08(6,21),IB10(6,21), IA08
     +(6,21),IA10(6,21), A1,B1,B2,B3,LT,LE,BET,AS,B8,AME,DIQ,ISU
*KEEP,FERMI.
      COMMON /FERMI/ PQUAR(4,248),PAQUAR(4,248), TQUAR(4,248),TAQUAR
     +(4,248)
*KEEP,KETMAS.
      COMMON /KETMAS/ AM8(2),AM10(2),IB88(2),IB1010(2),AMCH1N,AMCH2N
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
*KEEP,REJEC.
      COMMON /REJEC/ IRCO1,IRCO2,IRCO3,IRCO4,IRCO5, IRSS11,IRSS12,
     +IRSS13,IRSS14, IRSV11,IRSV12,IRSV13,IRSV14, IRVS11,IRVS12,IRVS13,
     +IRVS14, IRVV11,IRVV12,IRVV13,IRVV14
*KEEP,PROJK.
      COMMON /PROJK/ IPROJK
*KEEP,TANUIN.
      COMMON /TANUIN/ TASUMA,TASUBI,TABI,TAMASU,TAMA,TAIMMA
*KEND.
      COMMON /DIFFRA/ ISINGD,IDIFTP,IOUDIF,IFLAGD
*
      LOGICAL LSEADI
      COMMON /SEADIQ/ LSEADI
      COMMON /EVFLAG/NUMEV
      COMMON /DIQUAX/AMEDD,IDIQUA,IDIQUU
C
C-----------------------------------------------------------------------
C     PARAMETER (INTMX=2488)
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
      COMMON /NUCJTN/NONUJ1,NONUJT,NONUS1,NONUST
      COMMON /XSVTHR/ XSTHR,XVTHR,XDTHR,XSSTHR
      COMMON /MINIJ/IMINIJ,NOMJE,NOMJER,NREJEV,NOMJT,NOMJTR
      COMMON/PYJETS/NLU,NPAD,KLU(4000,5),PLU(4000,5),VLU(4000,5)
      COMMON/POL/POLARX(4),PMODUL
      COMMON /NEUREJ/ NONEUR
C     DIMENSION KKQ(1000,2),PPP(1000,5)
      DATA INIQEL /0/
C*******************************************************************"
C
C     KINEMATICS
C
C********************************************************************
C
      IREJ = 0
*
      AAM(26)=AAM(23)
C
      KPROJ=1
      IF(IJPROJ.NE.0) KPROJ=IJPROJ
      KPROJ=5
      IJPROJ=5
      KTARG=1
      ATNUC=IT
      ITN=IT-ITZ
      APNUC=IP
      IPN=IP-IPZ
      AMPROJ =AAM(KPROJ)
      AMTAR =AAM(KTARG)
*                             nucleon-nucleon cms
C     IBPROJ=1
      EPROJJ=EPN
      PPROJJ= SQRT((EPN-AMPROJ)*(EPN+AMPROJ))
C     PPROJJ= EPN
      UMOJ= SQRT(AMPROJ**2 + AMTAR**2 + 2.*AMTAR*EPROJJ)
C     UMOJ= SQRT( AMTAR**2 + 2.*AMTAR*EPROJJ)
      GAMCM = (EPROJJ+AMTAR)/UMOJ
      BGCM=PPROJJ/UMOJ
      ECM=UMOJ
      PCMJ=GAMCM*PPROJJ - BGCM*EPROJJ
C
      IF(IPEV.GE.1)WRITE(6, 1000)IP,IPZ,IT,ITZ,IJPROJ,IBPROJ,
     + EPROJJ,PPROJJ,
     +AMPROJ,AMTAR,UMO,GAMCM,BGCM
 1000 FORMAT(' ENTRY KKEVNU'/ '    IP,IPZ,IT,ITZ,IJPROJ,IBPROJ',6I5/
     +'    EPROJJ,PPROJJ,AMPROJ,AMTAR,UMO,GAMCM,BGCM'/10E12.3)
 
C
C****                            CHANGE PARAMETERS FROM COMMON \INPDAT\
      AS=0.5
      B8=0.4
C                                CHAIN PT BIGGER THAN PARTICLE PT
      N9483=0
*  entry after rejection of an event because of kinematical reasons
*  several trials are made to realize a sampled Glauber event
   10 CONTINUE
      NDZ=0
      NZD=0
      N9483=N9483+1
      IF (MOD(N9483,200).EQ.0) THEN
        WRITE(6,'(A,I5,A,I5,A)') ' KKEVT: Glauber event',NUMEV,
     +  ' rejected after', N9483, ' trials'
        WRITE(6, 1010) NN,NP,NT
        WRITE(6,1020) IRCO1,IRCO2,IRCO3,IRCO4,IRCO5, IRSS11,IRSS12,
     +  IRSS13,IRSS14, IRSV11,IRSV12,IRSV13,IRSV14, IRVS11,IRVS12,
     +  IRVS13,IRVS14, IRVV11,IRVV12,IRVV13,IRVV14
        N9483=1
        GO TO 20
      ELSEIF(N9483.GT.1) THEN
                                                                 GOTO 30
      ENDIF
 1010 FORMAT (5X,' N9483 LOOP - NN, NP, NT',5I10)
 1020 FORMAT (5X,' N9483 LOOP - REJECTION COUNTERS ',/,5I8/2(8I8/))
C
C***************************************************************
C
C     SAMPLE NUMBERS OF COLLISION A LA SHMAKOV---------------
C
C
C                 This is done here for a h-A event to obtain
C                the positions of the nucleons of A
C
C
C     TOTAL NUMBER OF INTERACTIONS = NN
C     NUMBER OF INTERACTING NUCLEONS
C                  FROM PROJECTILE = NP
C                  FROM TARGET = NT
C
   20 CONTINUE
   22 CONTINUE
      CALL SHMAKO(IP,IT,BIMP,NN,NP,NT,JSSH,JTSH,PPROJ,KKMAT)
      BIMPAC=BIMP
      NSHMAC=NSHMAC+1
      NNSHMA=NN
      NPSHMA=NP
      NTSHMA=NT
*  entry for repeated trial to realize a sampled Glauber event
   30 CONTINUE
      IF (IPEV.GE.2) THEN
        WRITE(6, 1040) IP,IPZ,IT,ITZ,EPROJ,PPROJ,NN,NP,NT
 1040 FORMAT ('   752 FORM ',4I10,2F10.3,5I10)
        WRITE (6,'(/A,2I5,1PE10.2,3I5)') ' KKEVT: IP,IT,BIMP,NN,NP,NT ',
     +  IP,IT,BIMP,NN,NP,NT
        WRITE (6,'(/2A)')
     +  ' KKEVT: JSSH(KKK),JTSH(KKK),INTER1(KKK),INTER2(KKK),',
     +  ' PKOO(3,KKK),TKOO(3,KKK)'
        ITUM=MAX(IT,IP)
        DO 40 KKK=1,ITUM
          WRITE (6,'(4I5,6(1PE11.3))') JSSH(KKK),JTSH(KKK),INTER1(KKK),
     +    INTER2(KKK), PKOO(1,KKK),PKOO(2,KKK),PKOO(3,KKK), TKOO(1,KKK),
     +    TKOO(2,KKK),TKOO(3,KKK)
 
   40   CONTINUE
      ENDIF
C
C-----------------------------------------------------------------------
C                  STORE PROJECTILE HADRON/NUCLEONS INTO /HKKEVT/
C                            - PROJECTILE CENTRE AT ORIGIN IN SPACE
C                            - TARGET SHIFTED IN X DIRECTION BY
C                                             IMPACT PARAMETER 'BIMP'
C
C                            - SAMPLING OF NUCLEON TYPES
C                            - CONSISTENCY CHECK
C                              FOR SAMPLED P/N NUMBERS
C                            - INTERACTING PROJECTILES ISTHKK=11
C                              NONINTERACTING ...      ISTHKK=13
C                            - FERMI MOMENTA IN CORRESP. REST SYSTEM
C-----------
      NHKK=0
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
C
      NCPP=0
      NCPN=0
C                      DEFINE FERMI MOMENTA/ENERGIES FOR PROJECTILE
C
      PXFE=0.0
      PYFE=0.0
      PZFE=0.0
      DO 50 KKK=1,IP
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
C       IF (JSSH(KKK).GT.0) THEN
          ISTHKK(NHKK)=11
C       ELSE
C         ISTHKK(NHKK)=13
C       ENDIF
C*
          KPROJ=IJPROJ
          PHKK(1,NHKK)=0.
          PHKK(2,NHKK)=0.
          PHKK(3,NHKK)=0.
          PHKK(4,NHKK)=AAM(KPROJ)
          PHKK(5,NHKK)=AAM(KPROJ)
C
        KKPROJ(KKK)=KPROJ
        IDHKK(NHKK)=MPDGHA(KPROJ)
        JMOHKK(1,NHKK)=0
        JMOHKK(2,NHKK)=0
        JDAHKK(1,NHKK)=0
        JDAHKK(2,NHKK)=0
C
        PHKK(5,NHKK)=AAM(KPROJ)
        VHKK(1,NHKK)=PKOO(1,KKK)*1.E-12
        VHKK(2,NHKK)=PKOO(2,KKK)*1.E-12
        VHKK(3,NHKK)=PKOO(3,KKK)*1.E-12
        VHKK(4,NHKK)=0.
        WHKK(1,NHKK)=PKOO(1,KKK)*1.E-12
        WHKK(2,NHKK)=PKOO(2,KKK)*1.E-12
        WHKK(3,NHKK)=PKOO(3,KKK)*1.E-12
        WHKK(4,NHKK)=0.
        JHKKNP(KKK)=NHKK
C
        IF (IPHKK.GE.2) WRITE(6,1050) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     +  JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +  (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
 1050 FORMAT (I6,I4,5I6,9E10.2)
C
   50 CONTINUE
C
C-----------------------------------------------------------------------
C                  STORE TARGET HADRON/NUCLEONS INTO /HKKEVT/
C                            - PROJECTILE CENTRE AT ORIGIN IN SPACE
C                            - TARGET SHIFTED IN X DIRECTION BY
C                                             IMPACT PARAMETER 'BIMP'
C
C                            - SAMPLING OF NUCLEON TYPES
C                            - CONSISTENCY CHECK
C                              FOR SAMPLED P/N NUMBERS
C                            - INTERACTING TARGETS     ISTHKK=12
C                              NONINTERACTING ...      ISTHKK=14
C-----------
C---------------------
      NHADRI=0
      NCTP=0
      NCTN=0
C
      TXFE=0.0
      TYFE=0.0
      TZFE=0.0
      DO 70 KKK=1,IT
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
C       IF (JTSH(KKK).GT.0) THEN
C         ISTHKK(NHKK)=12
C         NHADRI=NHADRI+1
C         IF (NHADRI.EQ.1) IHTAWW=NHKK
C         IF (EPN.LE.EHADTW) THEN
C           IF (NHADRI.GT.1) ISTHKK(NHKK)=14
C         ENDIF
C       ELSE
          ISTHKK(NHKK)=14
C       ENDIF
        IF(IT.GE.2)THEN
        FRTNEU=FLOAT(ITN)/ATNUC
        SAMTES=RNDM(v)
        IF(SAMTES.LT.FRTNEU.AND.NCTN.LT.ITN) THEN
          KTARG=8
          NCTN=NCTN + 1
        ELSEIF(SAMTES.GE.FRTNEU.AND.NCTP.LT.ITZ) THEN
          KTARG=1
          NCTP=NCTP + 1
        ELSEIF(NCTN.LT.ITN) THEN
          KTARG=8
          NCTN=NCTN + 1
        ELSEIF(NCTP.LT.ITZ) THEN
          KTARG=1
          NCTP=NCTP + 1
        ENDIF
C
        IF(KTARG.EQ.1) THEN
          PFERM = TAMFEP
        ELSE
          PFERM = TAMFEN
        ENDIF
C       CALL FER4M(PFERM,FPX,FPY,FPZ,FE,KTARG)
        CALL FER4MT(IT,PFERM,FPX,FPY,FPZ,FE,KTARG)
CWRITE(6,*)' Fermi PFERM;FPX,FPY,FPZ,FE '
C    & 	,PFERM,FPX,FPY,FPZ,FE
        TXFE=TXFE + FPX
        TYFE=TYFE + FPY
        TZFE=TZFE + FPZ
        PHKK(1,NHKK)=FPX
        PHKK(2,NHKK)=FPY
        PHKK(3,NHKK)=FPZ
        PHKK(4,NHKK)=FE
        PHKK(5,NHKK)=AAM(KTARG)
        ELSE
        PHKK(1,NHKK)=0. 
        PHKK(2,NHKK)=0. 
        PHKK(3,NHKK)=0. 
        PHKK(4,NHKK)=AAM(KTARG)
        PHKK(5,NHKK)=AAM(KTARG)
        ENDIF
C
        KKTARG(KKK)=KTARG
        IDHKK(NHKK)=MPDGHA(KTARG)
        JMOHKK(1,NHKK)=0
        JMOHKK(2,NHKK)=0
        JDAHKK(1,NHKK)=0
        JDAHKK(2,NHKK)=0
        VHKK(1,NHKK)=(TKOO(1,KKK)+BIMP)*1.E-12
        VHKK(2,NHKK)=TKOO(2,KKK)*1.E-12
        VHKK(3,NHKK)=TKOO(3,KKK)*1.E-12
        VHKK(4,NHKK)=0.
        WHKK(1,NHKK)=(TKOO(1,KKK)+BIMP)*1.E-12
        WHKK(2,NHKK)=TKOO(2,KKK)*1.E-12
        WHKK(3,NHKK)=TKOO(3,KKK)*1.E-12
        WHKK(4,NHKK)=0.
        JHKKNT(KKK)=NHKK
C
        IF (IPHKK.GE.2) WRITE(6,1050) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     +  JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +  (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
C
   70 CONTINUE
C                          balance Sampled Fermi momenta
      IF(IT.GE.2) THEN
        TASUMA=ITZ*AAM(1) + (IT-ITZ)*AAM(8)
        TASUBI=0.0
        TAMASU=0.0
        TXFE=TXFE/IT
        TYFE=TYFE/IT
        TZFE=TZFE/IT
        DO 80 KKK=1,IT
          IHKK=KKK + IP
          PHKK(1,IHKK)=PHKK(1,IHKK) - TXFE
          PHKK(2,IHKK)=PHKK(2,IHKK) - TYFE
          PHKK(3,IHKK)=PHKK(3,IHKK) - TZFE
          PHKK(4,IHKK)=SQRT(PHKK(5,IHKK)** 2+ PHKK(1,IHKK)** 2+ PHKK
     +    (2,IHKK)** 2+ PHKK(3,IHKK)**2)
          ITSEC=MCIHAD(IDHKK(IHKK))
          TASUBI=TASUBI + PHKK(4,IHKK) - PHKK(5,IHKK) - TAEPOT(ITSEC)
          TAMASU=TAMASU + PHKK(4,IHKK) - TAEPOT(ITSEC)
        IF (IPHKK.GE.2) WRITE(6,1050) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +  JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +  (PHKK(KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
   80   CONTINUE
C***         definition of initial state
        TABI=-EBIND(IT,ITZ)
        TAMA=(IT-ITZ)*AAM(8) + ITZ*AAM(1) + TABI
        TAIMMA=TAMA - TAMASU
      ENDIF
C
      IF(IPEV.GT.2) THEN
        WRITE(6,'(/A/5X,A/5X,4(1PE11.3))') ' KKEVT: FERMI MOMENTA',
     +  'PRMFEP,PRMFEN, TAMFEP,TAMFEN', PRMFEP,PRMFEN, TAMFEP,TAMFEN
 
      ENDIF
C-----------------------------------------------------------------------
      IFLAGD = 0
C-----------------------------------------------------------------------
C
      IF (IPEV.GE.6) THEN
        ITUM=MAX0(IP,IT,NN)
        WRITE(6,'(A,I10)')' KKEVT ITUM loop limit',ITUM
        WRITE(6,'(A,2A)') ' KKEVT (AFTER XKSAMP):',
     +  ' JSSH, JSSHS, JTSH, JTSHS, INTER1, INTER2',
     +  ' PKOO(1),PKOO(2),PKOO(3), TKOO(1),TKOO(2),TKOO(3)'
        DO 100 KKK=1,ITUM
          WRITE (6,'(6I5,6(1PE11.3))') JSSH(KKK),JSSHS(KKK),JTSH(KKK),
     +    JTSHS(KKK), INTER1(KKK),INTER2(KKK), PKOO(1,KKK),PKOO(2,KKK),
     +    PKOO(3,KKK), TKOO(1,KKK),TKOO(2,KKK),TKOO(3,KKK)
 
 
  100   CONTINUE
      ENDIF
C-----------------------------------------------------------------------
C                 TRANSFORM MOMENTA OF INTERACTING NUCLEONS
C                 (INCLUDING FERMI MOMENTA FROM NUCLEUS REST FRAMES)
C                 INTO NUCLEON-NUCLEON CMS (DEFINED WITHOUT FERMI MOM.
      IF(IPEV.GE.2)WRITE(6,'(A)')' KKEVT before NUCMOM'
      DO 7745 IHKK=1,NHKK
        IF (IPHKK.GE.2) WRITE(6,1050) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +  JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +  (PHKK(KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
 7745 CONTINUE
      CALL NUCMOM
      IF(IPEV.GE.2)THEN
C	DO IKI=1,200
	WRITE(6,'(A)')' KKEVNU after NUCMOM'
C	ENDDO
      ENDIF
      NONUST=0
      NONUJT=0
      NOMJE=0
      NOMJER=0
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
      INIQEL=INIQEL+1
      IF(INIQEL.EQ.1)CALL MASS_INI
      IF(IPEV.GE.2)THEN
C	DO IKI=1,200
	WRITE(6,'(A)')' KKEVNU after MASS_INI'
C	ENDDO
      ENDIF
      NHKKH1=NHKK
C-----------------------------------------------------------
C
C-----------------------------------------------------------
C
C                              Select target nucleon
C
C            NEUDEC < 9 qel_pol events
C
C.  INPUT  : LTYP = neutrino type (1,...,6)nue,anue,numu,anumu
C                                          nutau,anutau
C     for LTYP=1,3,5 Target is neutron
C     for LTYP=2,4,6 Target is proton
C
C-----------------------------------------------------------
C
C            NEUDEC > 10 gen_delta events
C
C.  INPUT  : LTYP = neutrino type (1,...,6)nue,anue,numu,anumu
C                                          nutau,anutau
C           The target nucleus is choosen randomly p or n
C
C-----------------------------------------------------------
	 IF(NEUDEC.LE.9)THEN
	   LTYP=NEUTYP
	   IF(LTYP.EQ.1.OR.LTYP.EQ.3.OR.LTYP.EQ.5)NUCTYP=2112
	   IF(LTYP.EQ.2.OR.LTYP.EQ.4.OR.LTYP.EQ.6)NUCTYP=2212
	 ELSEIF(NEUDEC.GE.10)THEN
	   LTYP=NEUTYP
	   NUCTYP=2112
	   RTYP=RNDM(V)*IT+1.
	   AITZ=ITZ
	   IF(RTYP.LE.AITZ)NUCTYP=2212
	 ENDIF
C            Neutrino energy is EPN in lab
  202      CONTINUE
	   IKTA=IT*RNDM(V)+2.
      IF(IPEV.GE.2)THEN
	WRITE(6,*)' NEUTYP,NUCTYP,IKTA,IDHKK(IKTA)',
     *              NEUTYP,NUCTYP,IKTA,IDHKK(IKTA)
      ENDIF
	   IF(IDHKK(IKTA).NE.NUCTYP) GO TO 202
	   ISTHKK(IKTA)=12
C          ENDIF
C	   IF(KKQ(III,1).EQ.1)THEN
	   IF(NUCTYP.EQ.2112)NUCTOP=2
	   IF(NUCTYP.EQ.2212)NUCTOP=1
           PLU21=PHKK(1,IKTA)
	   PLU22=PHKK(2,IKTA)
	   PLU23=PHKK(3,IKTA)
	   PLU24=PHKK(4,IKTA)
	   PLU25=PHKK(5,IKTA)
C          Call one qeldmo event
C	   CALL GEN_QEL(EPN,LTYP,PLU21,PLU22,PLU23,PLU24,PLU25)
C          Call one qel-pol event
	   IF(NEUDEC.LT.9)THEN
 	   CALL QEL_POL(EPN,LTYP,PLU21,PLU22,PLU23,PLU24,PLU25)
	   ELSEIF(NEUDEC.EQ.10)THEN
	      JINT=1
      IF(IPEV.GE.2)THEN
        WRITE(6,*)' CALL GEN_DELTA',EPN,LTYP,NUCTOP,JINT,
     &	      PLU21,PLU22,PLU23,PLU24,PLU25
      ENDIF
              CALL GEN_DELTA(EPN,LTYP,NUCTOP,JINT,
     &	      PLU21,PLU22,PLU23,PLU24,PLU25)
	   ELSEIF(NEUDEC.EQ.11)THEN
	      JINT=2
              CALL GEN_DELTA(EPN,LTYP,NUCTOP,JINT,
     &	      PLU21,PLU22,PLU23,PLU24,PLU25)
	   ELSEIF(NEUDEC.EQ.20)THEN
	     CALL FILENU(EPNN,LTYP,NUTYP,PLU21,PLU22,PLU23,
     &                      NHAD,IFLAG,LEND )
	     CALL ROTATE
	     EPN=EPNN
C
C            first initialize everything for energy EPN

             CALL LTINI(5,EPN,PPPN,EEECM)
C
C
      KPROJ=1
      IF(IJPROJ.NE.0) KPROJ=IJPROJ
      KPROJ=5
      IJPROJ=5
      KTARG=1
      ATNUC=IT
      ITN=IT-ITZ
      APNUC=IP
      IPN=IP-IPZ
      AMPROJ =AAM(KPROJ)
      AMTAR =AAM(KTARG)
*                             nucleon-nucleon cms
C     IBPROJ=1
      EPROJJ=EPN
      PPROJJ= SQRT((EPN-AMPROJ)*(EPN+AMPROJ))
      PPN=PPROJJ
C     PPROJJ= EPN
      UMOJ= SQRT(AMPROJ**2 + AMTAR**2 + 2.*AMTAR*EPROJJ)
C     UMOJ= SQRT( AMTAR**2 + 2.*AMTAR*EPROJJ)
      GAMCM = (EPROJJ+AMTAR)/UMOJ
      BGCM=PPROJJ/UMOJ
      GACMS=GAMCM
      BGCMS=BGCM
      UMO=UMOJ
      EPROJ=EPROJJ
      PPROJ=PPROJJ
C     COMMON /NUCCMS/ GACMS,BGCMS,GALAB,BGLAB,BLAB,UMO,PCM,EPROJ,PPROJ
      ECM=UMOJ
      PCMJ=GAMCM*PPROJJ - BGCM*EPROJJ
      PCM=PCMJ
C
      IF(IPEV.GE.1)WRITE(6,*)' EPN,PPROJJ,UMOJ,GAMCM,BGCM,PCMJ,ECM',
     &EPN,PPROJJ,UMOJ,GAMCM,BGCM,PCMJ,ECM

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C            pick out interacting nucleon and give the
C            Fermi momentum PLU21,PLU22,PLU23 to it
	     NUCTYP=NUTYP
	     NEUTYP=LTYP
  702        CONTINUE
	     IKTA=IT*RNDM(V)+2.
             IF(IPEV.GE.1)THEN
	       WRITE(6,*)' NEUTYP,NUCTYP,IKTA,IDHKK(IKTA)',
     *              NEUTYP,NUCTYP,IKTA,IDHKK(IKTA)
             ENDIF
	     IF(IDHKK(IKTA).NE.NUTYP) GO TO 702
	     ISTHKK(IKTA)=12
             PHKK(1,IKTA)=PLU21
	     PHKK(2,IKTA)=PLU22
	     PHKK(3,IKTA)=PLU23
               PHKK(4,IKTA)=SQRT(PHKK(5,IKTA)**2+ 
     +         PHKK(1,IKTA)**2+ PHKK(2,IKTA)**2+ PHKK(3,IKTA)**2)
C                Balance Fermi momenta
	     TXFE=0.0
             TYFE=0.0
             TZFE=0.0
             DO 704 KKK=1,IT
	       III=KKK+1
	       TXFE=TXFE+PHKK(1,III)
	       TYFE=TYFE+PHKK(2,III)
	       TZFE=TZFE+PHKK(3,III)
  704        CONTINUE
             TXFE=TXFE/(IT-1)
             TYFE=TYFE/(IT-1)
             TZFE=TZFE/(IT-1)
             DO 705 KKK=1,IT
               IHKK=KKK + 1
	       IF(IHKK.NE.IKTA)THEN
               PHKK(1,IHKK)=PHKK(1,IHKK) - TXFE
               PHKK(2,IHKK)=PHKK(2,IHKK) - TYFE
               PHKK(3,IHKK)=PHKK(3,IHKK) - TZFE
               PHKK(4,IHKK)=SQRT(PHKK(5,IHKK)**2+ 
     +         PHKK(1,IHKK)**2+ PHKK(2,IHKK)**2+ PHKK(3,IHKK)**2)
	       ENDIF
  705        CONTINUE
	   ENDIF
           IF(INIQEL.LE.20)THEN
C	     DO IKI=1,40
	     CALL PYLIST(1)
C	     ENDDO
           ENDIF
C
C                        Write events to file qeld.evt
C                       this is now done in dpmnuc6.f
C
C                              ADD particle to HKKEVT COMMON	 
      IIIMAX = 5
      IF(NEUDEC.GE.10)IIIMAX=7
      IF(NEUDEC.EQ.20)IIIMAX=NHAD
      IF(KLU(1,2).EQ.16.OR.KLU(1,2).EQ.-16)THEN
        IF(NEUDEC.EQ.1)THEN
	  IIIMAX = 8
        ENDIF
      ENDIF
      DO 200 III=4,IIIMAX
	 IF(KLU(III,1).EQ.1)THEN 
	   NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
	   ISTHKK(NHKK)=KLU(III,1)
	   IDHKK(NHKK)=KLU(III,2)
	   IF (ISTHKK(NHKK).EQ.15)ISTHKK(NHKK)=2
	   IF (ISTHKK(NHKK).EQ.11)ISTHKK(NHKK)=2
	   JMOHKK(1,NHKK)=IKTA
	   JMOHKK(2,NHKK)=0
	   JDAHKK(1,NHKK)=0
	   JDAHKK(2,NHKK)=0
	   PHKK(1,NHKK)=PLU(III,1)
	   PHKK(2,NHKK)=PLU(III,2)
	   PHKK(3,NHKK)=PLU(III,3)
	   PHKK(4,NHKK)=PLU(III,4)
C	   WRITE(6,*)'PHKK',NHKK,(PHKK(IKLO,NHKK),IKLO=1,4)
           NRHKK=MCIHAD(IDHKK(NHKK))
C                     drop Pauli blocking test (is in dpmqelpo)
C	   IF(NRHKK.EQ.1.OR.NRHKK.EQ.8)THEN
C	    IF(NRHKK.EQ.1)THEN
C     IF(PHKK(4,NHKK).LE.TAEFEP+AAM(NRHKK))THEN
C	       WRITE(6,*)' Pauli Blocking of p',PHKK(4,NHKK),TAEFEP
C	     ENDIF
C	    ENDIF
C	    IF(NRHKK.EQ.8)THEN
C     IF(PHKK(4,NHKK).LE.TAEFEN+AAM(NRHKK))THEN
C       WRITE(6,*)' Pauli Blocking of n',PHKK(4,NHKK),TAEFEN
C     ENDIF
C    ENDIF
	  IF(NRHKK.EQ.1.OR.NRHKK.EQ.8)THEN
	    IF(PHKK(4,NHKK)-AAM(NRHKK).LE.TAEPOT(NRHKK))THEN
	     ISTHKK(NHKK)=16
	    ENDIF
          ENDIF
C   ENDIF
C	   PHKK(5,NHKK)=AAM(NRHKK)
	   PHKK(5,NHKK)=PLU(III,5)
	   VHKK(1,NHKK)=VHKK(1,IKTA)
	   VHKK(2,NHKK)=VHKK(2,IKTA)
	   VHKK(3,NHKK)=VHKK(3,IKTA)
	   VHKK(4,NHKK)=VHKK(4,IKTA)
           IF(III.EQ.4)THEN
             WHKK(1,NHKK)=POLARX(1)
             WHKK(2,NHKK)=POLARX(2)
             WHKK(3,NHKK)=POLARX(3)
             WHKK(4,NHKK)=POLARX(4)
           ELSE
	   WHKK(1,NHKK)=WHKK(1,IKTA)
	   WHKK(2,NHKK)=WHKK(2,IKTA)
	   WHKK(3,NHKK)=WHKK(3,IKTA)
	   WHKK(4,NHKK)=WHKK(4,IKTA)
           ENDIF
C	 ENDIF
	ENDIF
  200 CONTINUE
  201 CONTINUE
         CALL BACKROT
           IF(INIQEL.LE.20)THEN
	     CALL PYLIST(1)
           ENDIF
C
C                                  Transform into cms
 	DO 111 I=NHKKH1+1,NHKK
 	  PZNN=PHKK(3,I)
 	  ENN=PHKK(4,I)
 	  PHKK(3,I)=GACMS*PZNN-BGCMS*ENN
 	  PHKK(4,I)=GACMS*ENN-BGCMS*PZNN
C	   WRITE(6,*)'PHKK',I,(PHKK(IKLO,I),IKLO=1,4)
  111   CONTINUE
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      IF (IPEV.GE.1) THEN
C	DO IKI=1,200
        WRITE(6,'(/A/)') ' KKEVT: FINAL LIST OF ENTRIES TO /HKKEVT/'
C	ENDDO
          DO 121 IHKK=1,NHKK
          WRITE(6,1050) IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),
     +    JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK
     +    (KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
 
  121     CONTINUE
        ENDIF
C
C
  110 CONTINUE
C
C
C
      RETURN
      END
      SUBROUTINE KKEVDI(NHKKH1,EPN,PPN,KKMAT,IREJ)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      COMMON/INTNEZ/NDZ,NZD
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
C      /
C INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C------------------
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
*KEEP,SHMAKL.
C     INCLUDE (SHMAKL)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON/SHMAKL/JSSH(INTMX),JTSH(INTMX),INTER1(INTMX),INTER2(INTMX)
*KEEP,NUCC.
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
*KEEP,RPTSHM.
      COMMON /RPTSHM/ RPROJ,RTARG,BIMPAC
*KEEP,NSHMAK.
      COMMON /NSHMAK/ NNSHMA,NPSHMA,NTSHMA,NSHMAC,NSHMA2
*KEEP,ZENTRA.
      COMMON /ZENTRA/ ICENTR
*KEEP,NUCIMP.
      COMMON /NUCIMP/ PRMOM(5,248),TAMOM(5,248), PRMFEP,PRMFEN,TAMFEP,
     +TAMFEN, PREFEP,PREFEN,TAEFEP,TAEFEN, PREPOT(210),TAEPOT(210),
     +PREBIN,TAEBIN,FERMOD,ETACOU
*KEEP,DROPPT.
      LOGICAL INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA, IPADIS,
     +ISHMAL,LPAULI
      COMMON /DROPPT/ INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA,
     +IPADIS,ISHMAL,LPAULI
*KEEP,NNCMS.
      COMMON /NNCMS/  GAMCM,BGCM,UMO,PCM,EPROJ,PPROJ
*KEEP,NUCPOS.
      COMMON /NUCPOS/INVVP(248),INVVT(248),INVSP(248),INVST(248), NUVV,
     +NUVS,NUSV,NUSS,INSVP(248),INSVT(248),INSSP(248),INSST(248), ISVEAP
     +(248),ISVEAT(248),ISSEAP(248),ISSEAT(248), IVSEAP(248),IVSEAT
     +(248), ISLOSP(248),ISLOST(248),INOOP(248),INOOT(248),NUOO
*KEEP,TAUFO.
      COMMON /TAUFO/  TAUFOR,KTAUGE,ITAUVE,INCMOD
      COMMON /EVAPPP/IEVAP
*KEEP,RTAR.
      COMMON /RTAR/ RTARNU
*KEEP,INNU.
      COMMON /INNU/INUDEC
*KEEP,HADTHR.
      COMMON /HADTHR/ EHADTH,INTHAD
*KEEP,DINPDA.
      COMMON /DINPDA/ IMPS(6,6),IMVE(6,6),IB08(6,21),IB10(6,21), IA08
     +(6,21),IA10(6,21), A1,B1,B2,B3,LT,LE,BET,AS,B8,AME,DIQ,ISU
*KEEP,FERMI.
      COMMON /FERMI/ PQUAR(4,248),PAQUAR(4,248), TQUAR(4,248),TAQUAR
     +(4,248)
*KEEP,KETMAS.
      COMMON /KETMAS/ AM8(2),AM10(2),IB88(2),IB1010(2),AMCH1N,AMCH2N
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
*KEEP,REJEC.
      COMMON /REJEC/ IRCO1,IRCO2,IRCO3,IRCO4,IRCO5, IRSS11,IRSS12,
     +IRSS13,IRSS14, IRSV11,IRSV12,IRSV13,IRSV14, IRVS11,IRVS12,IRVS13,
     +IRVS14, IRVV11,IRVV12,IRVV13,IRVV14
*KEEP,PROJK.
      COMMON /PROJK/ IPROJK
*KEEP,TANUIN.
      COMMON /TANUIN/ TASUMA,TASUBI,TABI,TAMASU,TAMA,TAIMMA
*KEND.
      COMMON /DIFFRA/ ISINGD,IDIFTP,IOUDIF,IFLAGD
*
      LOGICAL LSEADI
       CHARACTER*108  A109
      COMMON /SEADIQ/ LSEADI
      COMMON /EVFLAG/NUMEV
      COMMON /DIQUAX/AMEDD,IDIQUA,IDIQUU
C
C-----------------------------------------------------------------------
C     PARAMETER (INTMX=2488)
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
      COMMON /NUCJTN/NONUJ1,NONUJT,NONUS1,NONUST
      COMMON /XSVTHR/ XSTHR,XVTHR,XDTHR,XSSTHR
      COMMON /MINIJ/IMINIJ,NOMJE,NOMJER,NREJEV,NOMJT,NOMJTR
      COMMON /FELIRE/AMRECD,KJPRO
      DIMENSION PPPP(4),RMAX(5),NOMAX(5)
C*******************************************************************"
C
C     KINEMATICS
C
C********************************************************************
C
      IREJ = 0
*
      AAM(26)=AAM(23)
C
      KPROJ=1
      IF(IJPROJ.NE.0) KPROJ=IJPROJ
      KTARG=1
      ATNUC=IT
      ITN=IT-ITZ
      APNUC=IP
      IPN=IP-IPZ
      AMPROJ =AAM(KPROJ)
      AMTAR =AAM(KTARG)
*                             nucleon-nucleon cms
C     IBPROJ=1
      EPROJ=EPN
      PPROJ = SQRT((EPN-AMPROJ)*(EPN+AMPROJ))
      UMO = SQRT(AMPROJ**2 + AMTAR**2 + 2.*AMTAR*EPROJ)
      GAMCM = (EPROJ+AMTAR)/UMO
      BGCM=PPROJ/UMO
      ECM=UMO
      PCM=GAMCM*PPROJ - BGCM*EPROJ
C
      IF(IPEV.GE.1) PRINT 1000,IP,IPZ,IT,ITZ,IJPROJ,IBPROJ, EPROJ,PPROJ,
     +AMPROJ,AMTAR,UMO,GAMCM,BGCM
 1000 FORMAT(' ENTRY KKEVDI'/ '    IP,IPZ,IT,ITZ,IJPROJ,IBPROJ',6I5/
     +'    EPROJ,PPROJ,AMPROJ,AMTAR,UMO,GAMCM,BGCM'/10E12.3)
 
C
C****                            CHANGE PARAMETERS FROM COMMON \INPDAT\
      AS=0.5
      B8=0.4
C                                CHAIN PT BIGGER THAN PARTICLE PT
      N9483=0
*  entry after rejection of an event because of kinematical reasons
*  several trials are made to realize a sampled Glauber event
   10 CONTINUE
      NDZ=0
      NZD=0
      N9483=N9483+1
      IF (MOD(N9483,200).EQ.0) THEN
        WRITE(6,'(A,I5,A,I5,A)') ' KKEVT: Glauber event',NUMEV,
     +  ' rejected after', N9483, ' trials'
        WRITE(6, 1010) NN,NP,NT
        WRITE(6,1020) IRCO1,IRCO2,IRCO3,IRCO4,IRCO5, IRSS11,IRSS12,
     +  IRSS13,IRSS14, IRSV11,IRSV12,IRSV13,IRSV14, IRVS11,IRVS12,
     +  IRVS13,IRVS14, IRVV11,IRVV12,IRVV13,IRVV14
        N9483=1
        GO TO 20
      ELSEIF(N9483.GT.1) THEN
                                                                 GOTO 30
      ENDIF
 1010 FORMAT (5X,' N9483 LOOP - NN, NP, NT',5I10)
 1020 FORMAT (5X,' N9483 LOOP - REJECTION COUNTERS ',/,5I8/2(8I8/))
C
C***************************************************************
C
C     SAMPLE NUMBERS OF COLLISION A LA SHMAKOV---------------
C
C
C                 This is done here for a h-A event to obtain
C                the positions of the nucleons of A
C
C
C     TOTAL NUMBER OF INTERACTIONS = NN
C     NUMBER OF INTERACTING NUCLEONS
C                  FROM PROJECTILE = NP
C                  FROM TARGET = NT
C
   20 CONTINUE
   22 CONTINUE
      CALL SHMAKO(IP,IT,BIMP,NN,NP,NT,JSSH,JTSH,PPROJ,KKMAT)
      BIMPAC=BIMP
      NSHMAC=NSHMAC+1
      NNSHMA=NN
      NPSHMA=NP
      NTSHMA=NT
*  entry for repeated trial to realize a sampled Glauber event
   30 CONTINUE
      IF (IPEV.GE.3) THEN
        WRITE(6, 1040) IP,IPZ,IT,ITZ,EPROJ,PPROJ,NN,NP,NT
 1040 FORMAT ('   752 FORM ',4I10,2F10.3,5I10)
        WRITE (6,'(/A,2I5,1PE10.2,3I5)') ' KKEVT: IP,IT,BIMP,NN,NP,NT ',
     +  IP,IT,BIMP,NN,NP,NT
        WRITE (6,'(/2A)')
     +  ' KKEVT: JSSH(KKK),JTSH(KKK),INTER1(KKK),INTER2(KKK),',
     +  ' PKOO(3,KKK),TKOO(3,KKK)'
        ITUM=MAX(IT,IP)
        DO 40 KKK=1,ITUM
          WRITE (6,'(4I5,6(1PE11.3))') JSSH(KKK),JTSH(KKK),INTER1(KKK),
     +    INTER2(KKK), PKOO(1,KKK),PKOO(2,KKK),PKOO(3,KKK), TKOO(1,KKK),
     +    TKOO(2,KKK),TKOO(3,KKK)
 
   40   CONTINUE
      ENDIF
C
C-----------------------------------------------------------------------
C                  STORE PROJECTILE HADRON/NUCLEONS INTO /HKKEVT/
C                            - PROJECTILE CENTRE AT ORIGIN IN SPACE
C                            - TARGET SHIFTED IN X DIRECTION BY
C                                             IMPACT PARAMETER 'BIMP'
C
C                            - SAMPLING OF NUCLEON TYPES
C                            - CONSISTENCY CHECK
C                              FOR SAMPLED P/N NUMBERS
C                            - INTERACTING PROJECTILES ISTHKK=11
C                              NONINTERACTING ...      ISTHKK=13
C                            - FERMI MOMENTA IN CORRESP. REST SYSTEM
C-----------
      NHKK=0
C
      NCPP=0
      NCPN=0
C                      DEFINE FERMI MOMENTA/ENERGIES FOR PROJECTILE
C
      PXFE=0.0
      PYFE=0.0
      PZFE=0.0
      DO 50 KKK=1,IP
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
C       IF (JSSH(KKK).GT.0) THEN
          ISTHKK(NHKK)=11
C       ELSE
C         ISTHKK(NHKK)=13
C       ENDIF
C*
          KPROJ=IJPROJ
          PHKK(1,NHKK)=0.
          PHKK(2,NHKK)=0.
          PHKK(3,NHKK)=0.
          PHKK(4,NHKK)=AAM(KPROJ)
          PHKK(5,NHKK)=AAM(KPROJ)
C
        KKPROJ(KKK)=KPROJ
        IDHKK(NHKK)=MPDGHA(KPROJ)
        JMOHKK(1,NHKK)=0
        JMOHKK(2,NHKK)=0
        JDAHKK(1,NHKK)=0
        JDAHKK(2,NHKK)=0
C
        PHKK(5,NHKK)=AAM(KPROJ)
        VHKK(1,NHKK)=PKOO(1,KKK)*1.E-12
        VHKK(2,NHKK)=PKOO(2,KKK)*1.E-12
        VHKK(3,NHKK)=PKOO(3,KKK)*1.E-12
        VHKK(4,NHKK)=0.
        WHKK(1,NHKK)=PKOO(1,KKK)*1.E-12
        WHKK(2,NHKK)=PKOO(2,KKK)*1.E-12
        WHKK(3,NHKK)=PKOO(3,KKK)*1.E-12
        WHKK(4,NHKK)=0.
        JHKKNP(KKK)=NHKK
C
        IF (IPHKK.GE.2) WRITE(6,1050) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     +  JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +  (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
 1050 FORMAT (I6,I4,5I6,9E10.2)
C
   50 CONTINUE
C
C-----------------------------------------------------------------------
C                  STORE TARGET HADRON/NUCLEONS INTO /HKKEVT/
C                            - PROJECTILE CENTRE AT ORIGIN IN SPACE
C                            - TARGET SHIFTED IN X DIRECTION BY
C                                             IMPACT PARAMETER 'BIMP'
C
C                            - SAMPLING OF NUCLEON TYPES
C                            - CONSISTENCY CHECK
C                              FOR SAMPLED P/N NUMBERS
C                            - INTERACTING TARGETS     ISTHKK=12
C                              NONINTERACTING ...      ISTHKK=14
C-----------
C---------------------
      NHADRI=0
      NCTP=0
      NCTN=0
C
      TXFE=0.0
      TYFE=0.0
      TZFE=0.0
      DO 70 KKK=1,IT
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
C       IF (JTSH(KKK).GT.0) THEN
C         ISTHKK(NHKK)=12
C         NHADRI=NHADRI+1
C         IF (NHADRI.EQ.1) IHTAWW=NHKK
C         IF (EPN.LE.EHADTW) THEN
C           IF (NHADRI.GT.1) ISTHKK(NHKK)=14
C         ENDIF
C       ELSE
          ISTHKK(NHKK)=14
C       ENDIF
        IF(IT.GE.2)THEN
        FRTNEU=FLOAT(ITN)/ATNUC
        SAMTES=RNDM(v)
        IF(SAMTES.LT.FRTNEU.AND.NCTN.LT.ITN) THEN
          KTARG=8
          NCTN=NCTN + 1
        ELSEIF(SAMTES.GE.FRTNEU.AND.NCTP.LT.ITZ) THEN
          KTARG=1
          NCTP=NCTP + 1
        ELSEIF(NCTN.LT.ITN) THEN
          KTARG=8
          NCTN=NCTN + 1
        ELSEIF(NCTP.LT.ITZ) THEN
          KTARG=1
          NCTP=NCTP + 1
        ENDIF
C
        IF(KTARG.EQ.1) THEN
          PFERM = TAMFEP
        ELSE
          PFERM = TAMFEN
        ENDIF
C       CALL FER4M(PFERM,FPX,FPY,FPZ,FE,KTARG)
        CALL FER4MT(IT,PFERM,FPX,FPY,FPZ,FE,KTARG)
        TXFE=TXFE + FPX
        TYFE=TYFE + FPY
        TZFE=TZFE + FPZ
        PHKK(1,NHKK)=FPX
        PHKK(2,NHKK)=FPY
        PHKK(3,NHKK)=FPZ
        PHKK(4,NHKK)=FE
        PHKK(5,NHKK)=AAM(KTARG)
        ELSE
        PHKK(1,NHKK)=0. 
        PHKK(2,NHKK)=0. 
        PHKK(3,NHKK)=0. 
        PHKK(4,NHKK)=AAM(KTARG)
        PHKK(5,NHKK)=AAM(KTARG)
        ENDIF
C
        KKTARG(KKK)=KTARG
        IDHKK(NHKK)=MPDGHA(KTARG)
        JMOHKK(1,NHKK)=0
        JMOHKK(2,NHKK)=0
        JDAHKK(1,NHKK)=0
        JDAHKK(2,NHKK)=0
        VHKK(1,NHKK)=(TKOO(1,KKK))*1.E-12
        VHKK(2,NHKK)=TKOO(2,KKK)*1.E-12
        VHKK(3,NHKK)=TKOO(3,KKK)*1.E-12
        VHKK(4,NHKK)=0.
        WHKK(1,NHKK)=(TKOO(1,KKK)+BIMP)*1.E-12
        WHKK(2,NHKK)=TKOO(2,KKK)*1.E-12
        WHKK(3,NHKK)=TKOO(3,KKK)*1.E-12
        WHKK(4,NHKK)=0.
        JHKKNT(KKK)=NHKK
C
        IF (IPHKK.GE.2) WRITE(6,1050) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     +  JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +  (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
C
   70 CONTINUE
C                          balance Sampled Fermi momenta
      IF(IT.GE.2) THEN
        TASUMA=ITZ*AAM(1) + (IT-ITZ)*AAM(8)
        TASUBI=0.0
        TAMASU=0.0
        TXFE=TXFE/IT
        TYFE=TYFE/IT
        TZFE=TZFE/IT
        DO 80 KKK=1,IT
          IHKK=KKK + IP
          PHKK(1,IHKK)=PHKK(1,IHKK) - TXFE
          PHKK(2,IHKK)=PHKK(2,IHKK) - TYFE
          PHKK(3,IHKK)=PHKK(3,IHKK) - TZFE
          PHKK(4,IHKK)=SQRT(PHKK(5,IHKK)** 2+ PHKK(1,IHKK)** 2+ PHKK
     +    (2,IHKK)** 2+ PHKK(3,IHKK)**2)
          ITSEC=MCIHAD(IDHKK(IHKK))
          TASUBI=TASUBI + PHKK(4,IHKK) - PHKK(5,IHKK) - TAEPOT(ITSEC)
          TAMASU=TAMASU + PHKK(4,IHKK) - TAEPOT(ITSEC)
        IF (IPHKK.GE.2) WRITE(6,1050) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +  JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +  (PHKK(KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
   80   CONTINUE
C***         definition of initial state
        TABI=-EBIND(IT,ITZ)
        TAMA=(IT-ITZ)*AAM(8) + ITZ*AAM(1) + TABI
        TAIMMA=TAMA - TAMASU
      ENDIF
C
      IF(IPEV.GT.3) THEN
        WRITE(6,'(/A/5X,A/5X,4(1PE11.3))') ' KKEVT: FERMI MOMENTA',
     +  'PRMFEP,PRMFEN, TAMFEP,TAMFEN', PRMFEP,PRMFEN, TAMFEP,TAMFEN
 
      ENDIF
C-----------------------------------------------------------------------
      IFLAGD = 0
C-----------------------------------------------------------------------
C
      IF (IPEV.GE.6) THEN
        ITUM=MAX0(IP,IT,NN)
        WRITE(6,'(A,I10)')' KKEVT ITUM loop limit',ITUM
        WRITE(6,'(A,2A)') ' KKEVT (AFTER XKSAMP):',
     +  ' JSSH, JSSHS, JTSH, JTSHS, INTER1, INTER2',
     +  ' PKOO(1),PKOO(2),PKOO(3), TKOO(1),TKOO(2),TKOO(3)'
        DO 100 KKK=1,ITUM
          WRITE (6,'(6I5,6(1PE11.3))') JSSH(KKK),JSSHS(KKK),JTSH(KKK),
     +    JTSHS(KKK), INTER1(KKK),INTER2(KKK), PKOO(1,KKK),PKOO(2,KKK),
     +    PKOO(3,KKK), TKOO(1,KKK),TKOO(2,KKK),TKOO(3,KKK)
 
 
  100   CONTINUE
      ENDIF
C-----------------------------------------------------------------------
C                 TRANSFORM MOMENTA OF INTERACTING NUCLEONS
C                 (INCLUDING FERMI MOMENTA FROM NUCLEUS REST FRAMES)
C                 INTO NUCLEON-NUCLEON CMS (DEFINED WITHOUT FERMI MOM.
      IF(IPEV.GE.6)WRITE(6,'(A)')' KKEVT before NUCMOM'
      DO 7745 IHKK=1,NHKK
        IF (IPHKK.GE.2) WRITE(6,1050) IHKK,ISTHKK(IHKK),IDHKK(IHKK),
     +  JMOHKK(1,IHKK),JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),
     +  (PHKK(KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
 7745 CONTINUE
      CALL NUCMOM
      IF(IPEV.GE.6)WRITE(6,'(A)')' KKEVT after NUCMOM'
      NONUST=0
      NONUJT=0
      NOMJE=0
      NOMJER=0
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
      NHKKH1=NHKK
C-----------------------------------------------------------------------
C
C                 Read momentum shifts for nucleons from file
C                             diffnuc.evt
C
C                  KFORM=1 J.R.file diffnuc.evt
C                  KFORM=2 R.E.file diffnuc.evt
C----------------------------------------------------------------------
      KFORM=2
  214 CONTINUE
      IF(KFORM.EQ.1)THEN
        READ(29,'(I5,4E15.6)')NDIFFN,PPPP(1),PPPP(2),PPPP(3),PPPP(4)
      ELSEIF(KFORM.EQ.2)THEN
	READ(29,'(1X,I5,E12.4)')KJPRO,AMRECD
	WRITE(6,'(1X,I5,E12.4)')KJPRO,AMRECD
	READ(29,'(1X,I5)')IMIST
C	WRITE(6,'(1X,I5)')IMIST
	READ(29,'(1X,I5)')IMIST
C	WRITE(6,'(1X,I5)')IMIST
	READ(29,'(1X,I5,4E18.10)')IMIST,PPPP(1),PPPP(2),PPPP(3),PPPP(4)
	WRITE(6,'(1X,I5,4E18.10)')IMIST,PPPP(1),PPPP(2),PPPP(3),PPPP(4)
C                    The following READ were fine for CD
C	READ(29,'(1X,I5)')IMIST
C	WRITE(6,'(1X,I5)')IMIST
C	READ(29,'(1X,I5)')IMIST
C	WRITE(6,'(1X,I5)')IMIST
C	READ(29,'(1X,I5)')IMIST
C	WRITE(6,'(1X,I5)')IMIST
      ENDIF
C-----------------------------------------------------------------------
C
C               determine one of the target nucleons as involved
C                         in the diffractive scattering
C
C-----------------------------------------------------------------------
      RMAX(1)=0.
      RMAX(2)=0.
      RMAX(3)=0.
      RMAX(4)=0.
      RMAX(5)=0.
      NOMAX(1)=0
      NOMAX(2)=0
      NOMAX(3)=0
      NOMAX(4)=0
      NOMAX(5)=0
      DO 211 I=2,NHKK
	RRRN=SQRT(VHKK(1,I)**2+VHKK(2,I)**2)
	IF(RMAX(1).LT.RRRN)THEN
	  RMAX(1)=RRRN
	  NOMAX(1)=I
        ENDIF
  211 CONTINUE
      DO 212 I=2,NHKK
	IF(I.EQ.NOMAX(1))GO TO 212
	RRRN=SQRT(VHKK(1,I)**2+VHKK(2,I)**2)
	IF(RMAX(2).LT.RRRN)THEN
	  RMAX(2)=RRRN
	  NOMAX(2)=I
        ENDIF
  212 CONTINUE
      DO 213 I=2,NHKK
	IF(I.EQ.NOMAX(1))GO TO 213
	IF(I.EQ.NOMAX(2))GO TO 213
	RRRN=SQRT(VHKK(1,I)**2+VHKK(2,I)**2)
	IF(RMAX(3).LT.RRRN)THEN
	  RMAX(3)=RRRN
	  NOMAX(3)=I
        ENDIF
  213 CONTINUE
C-----------------------------------------------------------------------
C
C        have interaction with nucleon nomax(3)
C
C     READ(29,'(I5,4E15.6)')NDIFFN,PPPP(1),PPPP(2),PPPP(3),PPPP(4)
C-----------------------------------------------------------------------
C
      NWEPAU=0
 215  CONTINUE
      IF(NWEPAU.EQ.0)IINT=NOMAX(3)
      IF(NWEPAU.EQ.1)IINT=NOMAX(2)
      IF(NWEPAU.EQ.2)IINT=NOMAX(1)
      NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
      ISTHKK(NHKK)=1
      IDHKK(NHKK)=IDHKK(IINT)
      JMOHKK(1,NHKK)=IINT
      JMOHKK(2,NHKK)=0
      JDAHKK(1,NHKK)=0
      JDAHKK(2,NHKK)=0
      NRHKK=MCIHAD(IDHKK(NHKK))
      PHKK(1,NHKK)=PHKK(1,IINT)+PPPP(1)
      PHKK(2,NHKK)=PHKK(2,IINT)+PPPP(2)
      PHKK(3,NHKK)=PHKK(3,IINT)+PPPP(3)
      PHKK(4,NHKK)=SQRT(PHKK(1,NHKK)**2+PHKK(2,NHKK)**2+
     * PHKK(3,NHKK)**2+AAM(NRHKK)**2)
      PHKK(5,NHKK)=AAM(NRHKK)
      IF(NRHKK.EQ.-1.OR.NRHKK.EQ.-8)THEN
       IF(NRHKK.EQ.1)THEN
        IF(PHKK(4,NHKK).LE.TAEFEP+AAM(NRHKK))THEN
         WRITE(6,*)' Pauli Blocking of p',NWEPAU,PHKK(4,NHKK),TAEFEP
	 NWEPAU=NWEPAU+1
	 NHKK=NHKK-1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
C 	 IF(NWEPAU.LE.2)GO TO 215
         KFORM=2
         IF(KFORM.EQ.1)THEN
           AABBCC=0.
         ELSEIF(KFORM.EQ.2.AND.IREJ.EQ.0)THEN
C          The next 3 lines only for 6 (J/psi)
           READ(29,'(1X,I5)')KREPA
           READ(29,'(1X,I5)')KREPA
           READ(29,'(1X,I5)')KREPA
C
           READ(29,'(1X,I5)')KREPA
           DO 1975 KRE=1,KREPA
             READ(29,'(1X,A)')A109
 1975      CONTINUE
         ENDIF
         GO TO 214
        ENDIF
       ENDIF
       IF(NRHKK.EQ.8)THEN
        IF(PHKK(4,NHKK).LE.TAEFEN+AAM(NRHKK))THEN
         WRITE(6,*)' Pauli Blocking of n',NWEPAU,PHKK(4,NHKK),TAEFEN
	 NWEPAU=NWEPAU+1
	 NHKK=NHKK-1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' :NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
C	 IF(NWEPAU.LE.2)GO TO 215
       KFORM=2
       IF(KFORM.EQ.1)THEN
         AABBCC=0.
       ELSEIF(KFORM.EQ.2.AND.IREJ.EQ.0)THEN
C             The next 3 lines only for 6 (J/psi)
         READ(29,'(1X,I5)')KREPA
         READ(29,'(1X,I5)')KREPA
         READ(29,'(1X,I5)')KREPA
C
         READ(29,'(1X,I5)')KREPA
         DO 1976 KRE=1,KREPA
           READ(29,'(1X,A)')A109
 1976    CONTINUE
       ENDIF
	 GO TO 214
        ENDIF
       ENDIF
       IF(PHKK(4,NHKK)-AAM(NRHKK).LE.TAEPOT(NRHKK))THEN
        ISTHKK(NHKK)=16
       ENDIF
      ENDIF
      ISTHKK(IINT)=12
      IKTA=IINT
      VHKK(1,NHKK)=VHKK(1,IKTA)
      VHKK(2,NHKK)=VHKK(2,IKTA)
      VHKK(3,NHKK)=VHKK(3,IKTA)
      VHKK(4,NHKK)=VHKK(4,IKTA)
      WHKK(1,NHKK)=WHKK(1,IKTA)
      WHKK(2,NHKK)=WHKK(2,IKTA)
      WHKK(3,NHKK)=WHKK(3,IKTA)
      WHKK(4,NHKK)=WHKK(4,IKTA)
C
C-----------------------------------------------------------------------
        IF (IPEV.GE.1) THEN
        WRITE(6,'(/A/)') ' KKEVT: FINAL LIST OF ENTRIES TO /HKKEVT/'
          DO 121 IHKK=1,NHKK
          WRITE(6,1050) IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),
     +    JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK
     +    (KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
 
  121     CONTINUE
        ENDIF
C
C
  110 CONTINUE
C
C
C
      RETURN
      END
      SUBROUTINE LUINOL
C
C
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
C                           Prevent particles dacaying
C                 KOS
      KC=LUCOMP(310)
      MDCY(KC,1)=0
C                 PIO
      KC=LUCOMP(111)
      MDCY(KC,1)=0
C                 LAMBDA
      KC=LUCOMP(3122)
      MDCY(KC,1)=0
C                 ALAMBDA
      KC=LUCOMP(-3122)
      MDCY(KC,1)=0
C                 SIG+
      KC=LUCOMP(3222)
      MDCY(KC,1)=0
C                 ASIG+
      KC=LUCOMP(-3222)
      MDCY(KC,1)=0
C                 SIG-
      KC=LUCOMP(3112)
      MDCY(KC,1)=0
C                 ASIG-
      KC=LUCOMP(-3112)
      MDCY(KC,1)=0
C                 SIG0
C     KC=LUCOMP(3212)
C     MDCY(KC,1)=0
C                 ASIG0
C     KC=LUCOMP(-3212)
C     MDCY(KC,1)=0
C                 TET0
      KC=LUCOMP(3322)
      MDCY(KC,1)=0
C                 ATET0
      KC=LUCOMP(-3322)
      MDCY(KC,1)=0
C                 TET-
      KC=LUCOMP(3312)
      MDCY(KC,1)=0
C                 ATET-
      KC=LUCOMP(-3312)
      MDCY(KC,1)=0
C                 OMEGA-
      KC=LUCOMP(3334)
      MDCY(KC,1)=0
C                 AOMEGA-
      KC=LUCOMP(-3334)
      MDCY(KC,1)=0
C
      RETURN
      END
      SUBROUTINE TESTFILENU
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C      COMMON/BATLUND/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PHIROT/phr1,phr2,phr3
      INTEGER NFLAG(7)
      DO NEV = 1,100000
         CALL FILENU(EPN,LTYP,NUTYP,PLU21,PLU22,PLU23,NHAD,IFLAG,LEND)
         IF(LEND.EQ.1) GO TO 100
         NFLAG(IFLAG) = NFLAG(IFLAG) + 1
         write(6,150) (kw,K(kw,1),k(kw,2),(p(kw,j),j=1,5),kw=1,N)
         write(6,*)
C
C        Here rotates the event with the neutrino along +z
C
         CALL ROTATE
         write(6,150) (kw,K(kw,1),k(kw,2),(p(kw,j),j=1,5),kw=1,N)
         write(6,*)
c
c     Here returns the control to DPMJET
c
C
C        Here rotates the event back to lab frame
C
         CALL BACKROT
         write(6,150) (kw,K(kw,1),k(kw,2),(p(kw,j),j=1,5),kw=1,N)
         write(6,*)
      END DO
 100  CONTINUE
      WRITE(6,*) (NFLAG(J),J=1,7)
      STOP
 150  FORMAT(I5,2I5,5G10.3)
      END

      SUBROUTINE FILENU(EPN,LTYP,NUTYP,PLU21,PLU22,PLU23,NHO,
     $     IFLAG,LEND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C     COMMON/BATLUND/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON /CLOUT/ LUN
      LOGICAL FIRST
      DATA FIRST/.TRUE./
      DATA LUN /25/
      DATA INIT/0/
      SAVE FIRST
      OPEN (LUN,FILE='nuatm_new.dat',STATUS='OLD')
      IF(FIRST) THEN
         CALL READ_INI 
         FIRST = .FALSE.
      ENDIF
      INIT=INIT+1
      LEND = 0
      IFLAG = 0
      NHAD = 0
      READ (LUN, 10, ERR=1) NEV,  N, (V(1,J),J=1,3)      
      NHO=N
      DO L=1,N
	 READ (LUN, 15)LL, (K(L,J),J=1,5),(P(L,J),J=1,5)
         IF(L.GT.4.AND.K(L,1).EQ.1) NHAD = NHAD + 1
      ENDDO
      NONO=6
      IF(INIT.LE.20) THEN
         DO L=1,N
            WRITE(6, 15) L, (K(L,J),J=1,5),(P(L,J),J=1,5)
         ENDDO
      ENDIF
      EPN = P(1,4)
      LTYP = K(1,2)
      NUTYP = K(2,2)
      PLU21 = P(2,1)
      PLU21 = P(2,2)
      PLU21 = P(2,3)
      IF(N.EQ.4.OR.N.EQ.5) THEN
         IF(K(4,2).NE.K(1,2)) THEN
            IFLAG = 1                      ! quasi-elastic CC
         ELSE IF(K(4,2).EQ.K(1,2)) THEN
            IFLAG = 2                      ! quasi-elastic NC
         ENDIF
      ELSE IF(N.EQ.7) THEN
         IF(K(4,2).NE.K(1,2)) THEN
            IFLAG = 3                      ! delta resonance CC
         ELSE IF(K(4,2).EQ.K(1,2)) THEN
            IFLAG = 4                      ! delta resonance NC
         ENDIF
      ELSE IF(N.GT.7) THEN
         IF(K(4,2).NE.K(1,2)) THEN
            IFLAG = 5                      ! DIS CC
         ELSE IF(K(4,2).EQ.K(1,2)) THEN
            IFLAG = 6                      ! DIS NC
         ENDIF
      ELSE
         WRITE(6,*) N,NEV,K(1,2),K(4,2)
         IFLAG  = 7                        ! impossible
      ENDIF
      RETURN
1     LEND  = 1
10    FORMAT(1X,I7, 3X, I3, 2X, 3G12.4)
15    FORMAT(I5,5I7,5G12.4)
      END

      SUBROUTINE  READ_INI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      COMMON /CLOUT/ LUN
      REAL*4  RRAT(6),EMIN,VERS,AK
      CHARACTER*50  LINE
      DO J=1,10000
         READ(LUN, 10)   LINE
         IF (LINE(1:1) .EQ. '!')   GOTO 100
      ENDDO
 100  CONTINUE
      READ (LUN,  110) VERS, JCODE, JFLUX, JRAT, AK
      READ (LUN, 120) EMIN, (RRAT(J),J=1,6)
      RETURN
 10   FORMAT (A50)
 110  FORMAT(1X,F5.2,3X, I2, 3X, 2I2, 3X, F10.2)
 120  FORMAT(1X,F12.4, 3X, 6G12.4)
      END


      SUBROUTINE TESTROT1S(PI,PO,PHI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      DIMENSION ROT(3,3),PI(3),PO(3)
      ROT(1,1)=1.D0
      ROT(1,2)=0.D0
      ROT(1,3)=0.D0
      ROT(2,1)=0.D0
      ROT(2,2)=cos(phi)
      ROT(2,3)=-sin(phi)
      ROT(3,1)=0.D0
      ROT(3,2)=SIN(phi)
      ROT(3,3)=COS(phi)
      DO 140 J=1,3
        PO(J)=ROT(J,1)*PI(1)+ROT(J,2)*PI(2)+ROT(J,3)*PI(3)
  140 CONTINUE
      RETURN
      END


      SUBROUTINE TESTROT2S(PI,PO,PHI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      DIMENSION ROT(3,3),PI(3),PO(3)
      ROT(1,1)=0.D0
      ROT(1,2)=1.D0
      ROT(1,3)=0.D0
      ROT(2,1)=cos(phi)
      ROT(2,2)=0.D0
      ROT(2,3)=-sin(phi)
      ROT(3,1)=sin(phi)
      ROT(3,2)=0.D0
      ROT(3,3)=COS(phi)
      DO 140 J=1,3
        PO(J)=ROT(J,1)*PI(1)+ROT(J,2)*PI(2)+ROT(J,3)*PI(3)
  140 CONTINUE
      RETURN
      END

      SUBROUTINE TESTROT3S(PI,PO,PHI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      DIMENSION ROT(3,3),PI(3),PO(3)
      ROT(1,1)=0.D0
      ROT(2,1)=1.D0
      ROT(3,1)=0.D0
      ROT(1,2)=cos(phi)
      ROT(2,2)=0.D0
      ROT(3,2)=-sin(phi)
      ROT(1,3)=sin(phi)
      ROT(2,3)=0.D0
      ROT(3,3)=COS(phi)
      DO 140 J=1,3
        PO(J)=ROT(J,1)*PI(1)+ROT(J,2)*PI(2)+ROT(J,3)*PI(3)
  140 CONTINUE
      RETURN
      END

      SUBROUTINE TESTROT4S(PI,PO,PHI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      DIMENSION ROT(3,3),PI(3),PO(3)
      ROT(1,1)=1.D0
      ROT(2,1)=0.D0
      ROT(3,1)=0.D0
      ROT(1,2)=0.D0
      ROT(2,2)=cos(phi)
      ROT(3,2)=-sin(phi)
      ROT(1,3)=0.D0
      ROT(2,3)=SIN(phi)
      ROT(3,3)=COS(phi)
      DO 140 J=1,3
        PO(J)=ROT(J,1)*PI(1)+ROT(J,2)*PI(2)+ROT(J,3)*PI(3)
  140 CONTINUE
      RETURN
      END

      SUBROUTINE ROTATE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PHIROT/phr1,phr2,phr3
      DIMENSION PI(3),PO(3)
C
C     Rotate events so that neutrino goeas along +z
C
      phr1=atan(p(1,2)/p(1,3))
      DO kw=1,N
         pi(1)=p(kw,1)
         pi(2)=p(kw,2)
         pi(3)=p(kw,3)
         CALL TESTROT1S(PI,Po,PHR1)
         DO ll=1,3
            IF(abs(po(ll)).LT.1.D-07) po(ll)=0.
         END DO
         p(kw,1)=po(1)
         p(kw,2)=po(2)
         p(kw,3)=po(3)
      END DO
      phr2=atan(p(1,1)/p(1,3))
      DO kw=1,N
         pi(1)=p(kw,1)
         pi(2)=p(kw,2)
         pi(3)=p(kw,3)
         CALL TESTROT2S(Pi,Po,PHR2)
         DO ll=1,3
            IF(abs(po(ll)).LT.1.D-07) po(ll)=0.
         END DO
         p(kw,1)=po(1)
         p(kw,2)=po(2)
         p(kw,3)=po(3)
      END DO
      phr3 = 0
      IF(p(1,3).lt.0) THEN
         phr3 = -1.
         DO kw=1,N
            P(kw,3) = -P(kw,3)
         END DO
      ENDIF
      RETURN
      END

      SUBROUTINE BACKROT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PHIROT/phr1,phr2,phr3
      DIMENSION PI(3),PO(3)

c
c     Rotates back event to lab frame
c
      IF(phr3.EQ.-1.) THEN
         DO kw=1,N
            P(kw,3) = -P(kw,3)
         END DO
      END IF
      DO kw=1,N
         pi(1)=p(kw,1)
         pi(2)=p(kw,2)
         pi(3)=p(kw,3)
         CALL TESTROT3S(Pi,Po,PHR2)
         DO ll=1,3
            IF(abs(po(ll)).LT.1.D-07) po(ll)=0.
         END DO
         p(kw,1)=po(1)
         p(kw,2)=po(2)
         p(kw,3)=po(3)
      END DO
      DO kw=1,N
         pi(1)=p(kw,1)
         pi(2)=p(kw,2)
         pi(3)=p(kw,3)
         CALL TESTROT4S(Pi,Po,PHR1)
         DO ll=1,3
            IF(abs(po(ll)).LT.1.D-07) po(ll)=0.
         END DO
         p(kw,1)=po(1)
         p(kw,2)=po(2)
         p(kw,3)=po(3)
      END DO
      RETURN
      END

      SUBROUTINE BACKDPM
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
      COMMON/PHIROT/phr1,phr2,phr3
      DIMENSION PI(3),PO(3)

c
c     Rotates back event to lab frame
c
      IF(phr3.EQ.-1.) THEN
         DO kw=1,NHKK
         IF((ISTHKK(KW).EQ.-1).OR.
     *    (ISTHKK(KW).EQ.1).OR.
     *    (ISTHKK(KW).EQ.1001))THEN
            PHKK(3,kw) = -PHKK(3,kw)
         ENDIF
         END DO
      END IF
         DO kw=1,NHKK
         IF((ISTHKK(KW).EQ.-1).OR.
     *    (ISTHKK(KW).EQ.1).OR.
     *    (ISTHKK(KW).EQ.1001))THEN
         pi(1)=pHKK(1,kw)
         pi(2)=pHKK(2,kw)
         pi(3)=pHKK(3,kw)
         CALL TESTROT3S(Pi,Po,PHR2)
         DO ll=1,3
            IF(abs(po(ll)).LT.1.D-07) po(ll)=0.
         END DO
         pHKK(1,kw)=po(1)
         pHKK(2,kw)=po(2)
         pHKK(3,kw)=po(3)
         ENDIF
      END DO
         DO kw=1,NHKK
         IF((ISTHKK(KW).EQ.-1).OR.
     *    (ISTHKK(KW).EQ.1).OR.
     *    (ISTHKK(KW).EQ.1001))THEN
         pi(1)=pHKK(1,kw)
         pi(2)=pHKK(2,kw)
         pi(3)=pHKK(3,kw)
         CALL TESTROT4S(Pi,Po,PHR1)
         DO ll=1,3
            IF(abs(po(ll)).LT.1.D-07) po(ll)=0.
         END DO
         pHKK(1,kw)=po(1)
         pHKK(2,kw)=po(2)
         pHKK(3,kw)=po(3)
         ENDIF
      END DO
      RETURN
      END

      SUBROUTINE DROPDI(NN,NP,NT,ECM)
C     Drop diffractive collisions out of the Glauber
C     cascade in nuclear collisions (For NN > 1 only)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (INTMX=2488)
      COMMON/SHMAKL/JSSH(INTMX),JTSH(INTMX),INTER1(INTMX),INTER2(INTMX)
C     Fraction of diffractive collisions at given Energy
      FRACDIF=SIPPSD(ECM)/SIINEL(1,1,ECM)
      ANN=NN
      DANN=FRACDIF*ANN
      IDANN=DANN
      AIDANN=IDANN
      FDANN=DANN-AIDANN
      IF(RNDM(V).LT.FDANN)IDANN=IDANN+1
C     Total number of collisions NN is reduced by IDANN
C     WRITE(6,*)'NN,FRACDIF,FDANN,IDANN ',NN,FRACDIF,FDANN,IDANN
      NNNEW=NN-IDANN
      NPNEW=NP
      NTNEW=NT
      IF(IDANN.GT.0)THEN
        DO 1 I= NN-IDANN+1,NN
	  NI1=INTER1(I)
	  NI2=INTER2(I)
	  JSSH(NI1)=JSSH(NI1)-1
	  JTSH(NI2)=JTSH(NI2)-1
	  IF(JSSH(NI1).EQ.0)NPNEW=NPNEW-1
	  IF(JTSH(NI2).EQ.0)NTNEW=NTNEW-1
	  INTER1(I)=0
	  INTER2(I)=0
 1      CONTINUE	
      ENDIF	
C     WRITE (6,*)'NNNEW,NPNEW,NTNEW ',NNNEW,NPNEW,NTNEW
      NN=NNNEW
      NP=NPNEW
      NT=NTNEW
      RETURN
      END
