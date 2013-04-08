C****************************************************************
      SUBROUTINE KKEVLE(NHKKH1,EPN,PPN,KKMAT,IREJ)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      COMMON /NSHMAK/ NNSHMA,NPSHMA,NTSHMA,NSHMAC
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
      COMMON /DIQUAX/IDIQUA
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
      INTEGER NLU,KLU,LST,MDCY,MDME,KFDP
      REAL PLU,VLU,CUT,PARL,X,Y,W2,Q2,U,BRAT,ELAB
      COMMON/LUJETS/NLU,KLU(4000,5),PLU(4000,5),VLU(4000,5)
      COMMON /LEPTOU/CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      COMMON /NEUREJ/ NONEUR
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
      UMOJ= SQRT(AMPROJ**2 + AMTAR**2 + 2.*AMTAR*EPROJJ)
      GAMCM = (EPROJJ+AMTAR)/UMOJ
      BGCM=PPROJJ/UMOJ
      ECM=UMOJ
      PCMJ=GAMCM*PPROJJ - BGCM*EPROJJ
C
      IF(IPEV.GE.1) PRINT 1000,IP,IPZ,IT,ITZ,IJPROJ,IBPROJ, EPROJ,PPROJ,
     +AMPROJ,AMTAR,UMO,GAMCM,BGCM
 1000 FORMAT(' ENTRY KKEVNU'/ '    IP,IPZ,IT,ITZ,IJPROJ,IBPROJ',6I5/
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
        ELSEIF(SAMTES.GE.FRPNEU.AND.NCTP.LT.ITZ) THEN
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
        CALL FER4M(PFERM,FPX,FPY,FPZ,FE,KTARG)
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
	DO IKI=1,200
	WRITE(6,'(A)')' KKEVNU after NUCMOM'
	ENDDO
      ENDIF
      NONUST=0
      NONUJT=0
      NOMJE=0
      NOMJER=0
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
      NHKKH1=NHKK
C-----------------------------------------------------------
C
C
C                              Select target nucleon
C
C-----------------------------------------------------------
	   LTYP=NEUTYP
C	   IF(LTYP.EQ.1.OR.LTYP.EQ.3.OR.LTYP.EQ.5)NUCTYP=2112
C	   IF(LTYP.EQ.2.OR.LTYP.EQ.4.OR.LTYP.EQ.6)NUCTYP=2212
C            Neutrino energy is EPN in lab
  202      CONTINUE
	   IKTA=IT*RNDM(V)+2.
      IF(IPEV.GE.2)THEN
	WRITE(6,*)' NEUTYP,NUCTYP,IKTA,IDHKK(IKTA)',
     *              NEUTYP,NUCTYP,IKTA,IDHKK(IKTA)
      ENDIF
      IF(IDHKK(IKTA).EQ.2112) ITAR=2 
      IF(IDHKK(IKTA).EQ.2212) ITAR=1
      INU=NEUTYP
      ISTHKK(IKTA)=12
      ELAB=EPN
      CUT(10)=ELAB+1.
      CUT(11)=0.
      CUT(12)=ELAB+1.
      CUT(13)=0.
      CUT(14)=3.1314
      LST(22)=ITAR
      LST(17)=1
C     LST(19)=-1
      LST(8)=0
* GRID SUITABLE FOR FIXED TARGET < 300 GEV
      LST(19)=1
      PLU(1,1)=0.
      PLU(1,2)=0.
      PLU(1,3)=ELAB
      PLU(1,4)=ELAB
      PLU(1,5)=0.
      PLU(2,1)=PHKK(1,IKTA)
      PLU(2,2)=PHKK(2,IKTA)
      PLU(2,3)=PHKK(3,IKTA)
      PLU(2,4)=PHKK(4,IKTA)
      PLU(2,5)=PHKK(5,IKTA)
C          Call one lepto event
      IF(INIQEL.EQ.0)CALL LINIT(0,INU,ELAB,0.,2)
      INIQEL=INIQEL+1
      PLU(1,1)=0.
      PLU(1,2)=0.
      PLU(1,3)=ELAB
      PLU(1,4)=ELAB
      PLU(1,5)=0.
      PLU(2,1)=PHKK(1,IKTA)
      PLU(2,2)=PHKK(2,IKTA)
      PLU(2,3)=PHKK(3,IKTA)
      PLU(2,4)=PHKK(4,IKTA)
      PLU(2,5)=PHKK(5,IKTA)
      LST(22)=ITAR
      CALL LEPTO
C                     event in lab frame
      CALL LFRAME(3,1)
      IF(LST(21).NE.0)THEN
        CALL LULIST(1)
        WRITE(6,*)' event rejected '
        GO TO 10
      ENDIF
      IF(INIQEL.LE.100)WRITE(6,*)' Event ',INIQEL
      IF(INIQEL.LE.100)CALL LULIST(1)
C                        Write events to file lepto.evt
      IIII=0
      DO 205 III=1,NLU
        IF(KLU(III,1).EQ.1.OR.III.LE.2) THEN
          IIII=IIII+1
          WRITE(29,'(3I6,5F10.3)')IIII,KLU(III,1),KLU(III,2),
     *   (PLU(III,KK),KK=1,5)
        ENDIF
  205 CONTINUE
      IIII=-1
      WRITE(29,'(I6)')IIII
C                              ADD particle to HKKEVT COMMON	 
      DO 206 III=4,NLU
	IF(KLU(III,1).EQ.1)THEN
	   NHKK=NHKK+1
	   ISTHKK(NHKK)=1
	   IDHKK(NHKK)=KLU(III,2)
	   JMOHKK(1,NHKK)=IKTA
	   JMOHKK(2,NHKK)=0
	   JDAHKK(1,NHKK)=0
	   JDAHKK(2,NHKK)=0
	   PHKK(1,NHKK)=PLU(III,1)
	   PHKK(2,NHKK)=PLU(III,2)
	   PHKK(3,NHKK)=PLU(III,3)
	   PHKK(4,NHKK)=PLU(III,4)
           NRHKK=MCIHAD(IDHKK(NHKK))
	   IF(NRHKK.EQ.1.OR.NRHKK.EQ.8)THEN
	    IF(NRHKK.EQ.1)THEN
	     IF(PHKK(4,NHKK).LE.TAEFEP+AAM(NRHKK))THEN
	       WRITE(6,*)' Pauli Blocking of p',PHKK(4,NHKK),TAEFEP
	     ENDIF
	    ENDIF
	    IF(NRHKK.EQ.8)THEN
	     IF(PHKK(4,NHKK).LE.TAEFEN+AAM(NRHKK))THEN
	       WRITE(6,*)' Pauli Blocking of n',PHKK(4,NHKK),TAEFEN
	     ENDIF
	    ENDIF
	    IF(PHKK(4,NHKK)-AAM(NRHKK).LE.TAEPOT(NRHKK))THEN
	     ISTHKK(NHKK)=16
	    ENDIF
	   ENDIF
	   PHKK(5,NHKK)=AAM(NRHKK)
	   VHKK(1,NHKK)=VHKK(1,IKTA)
	   VHKK(2,NHKK)=VHKK(2,IKTA)
	   VHKK(3,NHKK)=VHKK(3,IKTA)
	   VHKK(4,NHKK)=VHKK(4,IKTA)
	   WHKK(1,NHKK)=WHKK(1,IKTA)
	   WHKK(2,NHKK)=WHKK(2,IKTA)
	   WHKK(3,NHKK)=WHKK(3,IKTA)
	   WHKK(4,NHKK)=WHKK(4,IKTA)
 	 ENDIF
  206 CONTINUE
  201 CONTINUE
C
C                                  Transform into cms
 	DO 111 I=NHKKH1+1,NHKK
 	  PZNN=PHKK(3,I)
 	  ENN=PHKK(4,I)
 	  PHKK(3,I)=GACMS*PZNN-BGCMS*ENN
 	  PHKK(4,I)=GACMS*ENN-BGCMS*PZNN
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
