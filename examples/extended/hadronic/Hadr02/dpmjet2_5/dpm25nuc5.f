      SUBROUTINE DIQSV(ECM,ITV,J,IREJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
*  define d-v chains (sea diquark - valence chains)
*   sq-q and saqsaq-qq chains instead of sq-qq and saq-q chains
      COMMON /ZSEA/ZSEAAV,ZSEASU,ANZSEA
       COMMON/POPCCK/PDBCK,PDBSE,PDBSEU,
     *  IJPOCK,IREJCK,ICK4,IHAD4,ICK6,IHAD6
     *,IREJSE,ISE4,ISE6,IREJS3,ISE43,ISE63,IREJS0
     *,IHADA4,IHADA6,IREJSA,ISEA4,ISEA6,IREJA3,
     *ISEA43,ISEA63,IREJAO
*KEEP,INTMX.
      PARAMETER (INTMX=2488,INTMD=252)
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
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
*KEEP,ABRDV.
      COMMON /ABRDV/ AMCDV1(248),AMCDV2(248),GACDV1(248),GACDV2(248),
     +BGXDV1(248),BGYDV1(248),BGZDV1(248), BGXDV2(248),BGYDV2(248),
     +BGZDV2(248), NCHDV1(248),NCHDV2(248),IJCDV1(248),IJCDV2(248),
     +PQDVA1(248,4),PQDVA2(248,4), PQDVB1(248,4),PQDVB2(248,4)
*KEND.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
      COMMON /SEASU3/SEASQ
      COMMON /XSEADI/ XSEACU,UNON,UNOM,UNOSEA, CVQ,CDQ,CSEA,SSMIMA,
     +SSMIMQ,VVMTHR
      COMMON /DIQREJ/IDIQRE(7),IDVRE(3),IVDRE(3),IDSRE(3),ISDRE(3),
     *IDZRE(3),IZDRE(3),IDIQRZ(7) 
      COMMON /DIQSUM/NDVUU,NDVUS,NDVSS,NVDUU,NVDUS,NVDSS,
     *               NDSUU,NDSUS,NDSSS,NSDUU,NSDUS,NSDSS,
     *               NDZUU,NDZUS,NDZSS,NZDUU,NZDUS,NZDSS 
     *              ,NADVUU,NADVUS,NADVSS,NAVDUU,NAVDUS,NAVDSS,
     *               NADSUU,NADSUS,NADSSS,NASDUU,NASDUS,NASDSS,
     *               NADZUU,NADZUS,NADZSS,NAZDUU,NAZDUS,NAZDSS 
C------------------------------------------------------------------
C     COMMON /PCHARM/PCCCC
      PARAMETER (UMMM=0.3D0)
      PARAMETER (SMMM=0.5D0)
      PARAMETER (CMMM=1.3D0)
      DATA PC/0.0001D0/
*KEND.
C----------
C
      DATA INICHA/0/
C----------------------------------------------------------------------
C                     Initialize Charm selection at soft chain ends
C
      IF(INICHA.EQ.0)THEN
        RX=8.D0
        X1=RX
        GM=2.140D0
        X2=UMMM
	BETOO=7.5D0
      ENDIF
      RX=8.D0
      X1=RX
      BETCHA=BETOO+1.3D0-LOG10(ECM)
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
      PC=PC1/7.0D0
      PU1=PU/(2*PU+PS+PC)
      PS1=PS/(2*PU+PS+PC)
      IF(INICHA.EQ.0)THEN
        INICHA=1
        WRITE(6,4567)PC,BETCHA,PU1,PS1,SEASQ
 4567   FORMAT(' Charm chain ends DIQSV: PC,BETCHA,PU,PS,SEASQ ',5F10.5)
      ENDIF
C----------------------------------------------------------------------
      IREJ=0
*  kinematics: is the mass of the adiquark-diquark chain big enough
*              to allow for fragmentation
      IF(IPHKK.GE.6)WRITE (6,'( A)') ' diqsv'
      IPSQ2(J)=1.D0+RNDM(v1)*(2.D0+2.D0*SEASQ)
        RR=RNDM(V)
	IF(RR.LT.PC)IPSQ2(J)=4
C------------------------------------------------------------------
      IPSAQ2(J)=-IPSQ2(J)
C---------------------------------------------------j.r.29.4.94
C                                   x**1.5 distr for sea diquarks
C                         number of projectile nucleon
      INUCPR=IFROSP(J)
C                          number of projectile diquark
      IITOP=ITOVP(INUCPR)
C                          diquark x
      XPDIQU=XPVD(IITOP)
C                          minimal value of diquark x
      XDTHR=CDQ/ECM
C
      XDFREE=XPDIQU-XDTHR
      XALL=XDFREE+XPSQ(J)+XPSAQ(J)-2.*XDTHR
      XDALT=XPVD(IITOP)
      XSALT=XPSQ(J)
      XAALT=XPSAQ(J) 
      IF(XALL.GE.0.)THEN
        RR1=RNDM(V1)
        RR2=RNDM(V2)
        RR3=RNDM(V3)
        SR123=RR1+RR2+RR3
        DX1=RR1*XALL/SR123
        DX2=RR2*XALL/SR123
        DX3=RR3*XALL/SR123
        XPVD(IITOP)=XDTHR+DX1
        XPSQ(J)=XDTHR+DX2
        XPSAQ(J)=XDTHR+DX3
      ENDIF
C--------------------------------------------------------------
      AMDVQ1=XPSQ(J)*XTVQ(ITV)*ECM**2
      AMDVQ2=XPSAQ(J)*XTVD(ITV)*ECM**2
      IDIQRE(1)=IDIQRE(1)+1
      IF(IPSQ(J).GE.3.AND.IPSQ2(J).GE.3)THEN
        IDIQRE(2)=IDIQRE(2)+1
C       IF(AMDVQ2.LE.9.0.OR.AMDVQ1.LE.2.30) THEN
        IF(AMDVQ2.LE.17.0D0.OR.AMDVQ1.LE.6.60D0) THEN
          IREJ=1
           IDIQRE(3)=IDIQRE(3)+1
           IDIQRE(2)=IDIQRE(2)-1
           IDIQRE(1)=IDIQRE(1)-1
          XPVD(IITOP)=XDALT
          XPSQ(J)=XSALT
          XPSAQ(J)=XAALT
          RETURN
        ENDIF
      ELSEIF(IPSQ(J).GE.3.OR.IPSQ2(J).GE.3)THEN
      IDIQRE(4)=IDIQRE(4)+1
C       IF(AMDVQ2.LE.7.3.OR.AMDVQ1.LE.1.90) THEN
        IF(AMDVQ2.LE.13.6D0.OR.AMDVQ1.LE.5.80D0) THEN
          IREJ=1
           IDIQRE(5)=IDIQRE(5)+1
           IDIQRE(4)=IDIQRE(4)-1
           IDIQRE(1)=IDIQRE(1)-1
          XPVD(IITOP)=XDALT
          XPSQ(J)=XSALT
          XPSAQ(J)=XAALT
          RETURN
        ENDIF
      ELSE
      IDIQRE(6)=IDIQRE(6)+1
C       IF(AMDVQ2.LE.6.70.OR.AMDVQ1.LE.1.50) THEN
        IF(AMDVQ2.LE.12.40D0.OR.AMDVQ1.LE.3.9D0) THEN
          IREJ=1
           IDIQRE(7)=IDIQRE(7)+1
           IDIQRE(6)=IDIQRE(6)-1
           IDIQRE(1)=IDIQRE(1)-1
          XPVD(IITOP)=XDALT
          XPSQ(J)=XSALT
          XPSAQ(J)=XAALT
          RETURN
        ENDIF
      ENDIF
      NDV=NDV+1
      NCHDV1(NDV)=0
      NCHDV2(NDV)=0
      INTDV1(NDV)=J
      INTDV2(NDV)=ITV
      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C      DEBUG SUBCHK
C      END DEBUG
      SUBROUTINE KKEVDV(IREJDV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C------------------ treatment of sea diquark - valence CHAIN SYSTEMS
      COMMON /ZSEA/ZSEAAV,ZSEASU,ANZSEA
       COMMON/POPCCK/PDBCK,PDBSE,PDBSEU,
     *  IJPOCK,IREJCK,ICK4,IHAD4,ICK6,IHAD6
     *,IREJSE,ISE4,ISE6,IREJS3,ISE43,ISE63,IREJS0
     *,IHADA4,IHADA6,IREJSA,ISEA4,ISEA6,IREJA3,
     *ISEA43,ISEA63,IREJAO
C
*KEEP,INTMX.

      PARAMETER (INTMX=2488,INTMD=252)
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
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
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
*KEEP,ABRDV.
      COMMON /ABRDV/ AMCDV1(248),AMCDV2(248),GACDV1(248),GACDV2(248),
     +BGXDV1(248),BGYDV1(248),BGZDV1(248), BGXDV2(248),BGYDV2(248),
     +BGZDV2(248), NCHDV1(248),NCHDV2(248),IJCDV1(248),IJCDV2(248),
     +PQDVA1(248,4),PQDVA2(248,4), PQDVB1(248,4),PQDVB2(248,4)
*KEEP,DXQX.
C     INCLUDE (XQXQ)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON /DXQX/ XPVQ(248),XPVD(248),XTVQ(248),XTVD(248), XPSQ
     +(INTMX),XPSAQ(INTMX),XTSQ(INTMX),XTSAQ(INTMX)
     *                ,XPSU(248),XTSU(248)
     *                ,XPSUT(248),XTSUT(248)
*KEEP,LOZUO.
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/ ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),ZUOST(INTMX),
     +INTLO(INTMX),INLOSS(INTMX)
      COMMON /DIQREJ/IDIQRE(7),IDVRE(3),IVDRE(3),IDSRE(3),ISDRE(3),
     *IDZRE(3),IZDRE(3),IDIQRZ(7) 
C  /LOZUO/
C       INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C------------------
*KEEP,TRAFOP.
      COMMON /TRAFOP/ GAMP,BGAMP,BETP
*KEEP,NUCIMP.
      COMMON /NUCIMP/ PRMOM(5,248),TAMOM(5,248), PRMFEP,PRMFEN,TAMFEP,
     +TAMFEN, PREFEP,PREFEN,TAEFEP,TAEFEN, PREPOT(210),TAEPOT(210),
     +PREBIN,TAEBIN,FERMOD,ETACOU
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
      COMMON /DIQSUM/NDVUU,NDVUS,NDVSS,NVDUU,NVDUS,NVDSS,
     *               NDSUU,NDSUS,NDSSS,NSDUU,NSDUS,NSDSS,
     *               NDZUU,NDZUS,NDZSS,NZDUU,NZDUS,NZDSS 
     *              ,NADVUU,NADVUS,NADVSS,NAVDUU,NAVDUS,NAVDSS,
     *               NADSUU,NADSUS,NADSSS,NASDUU,NASDUS,NASDSS,
     *               NADZUU,NADZUS,NADZSS,NAZDUU,NAZDUS,NAZDSS 
*KEND.
C-------------------
      IF(IPHKK.GE.6)WRITE (6,'( A)') ' kkevdv'
      IREJDV=0
      DO 10 N=1,NDV
C---------------------------drop recombined chain pairs
        IF(NCHDV1(N).EQ.99.AND.NCHDV2(N).EQ.99)GO TO 10
C
C***                 4-MOMENTA OF PROJECTILE SEA-QUARK PAIRS IN NN-CMS
        IXSPR=INTDV1(N)
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
C***                 4-MOMENTA OF TARGET QUARK-DIQUARK PAIRS IN NN-CMS
        IXVTA=INTDV2(N)
        INUCTA=IFROVT(IXVTA)
        JNUCTA=ITOVT(INUCTA)
C
        TVQPX=XTVQ(IXVTA)*TAMOM(1,INUCTA)
        TVQPY=XTVQ(IXVTA)*TAMOM(2,INUCTA)
        TVQPZ=XTVQ(IXVTA)*TAMOM(3,INUCTA)
        TVQE=XTVQ(IXVTA)*TAMOM(4,INUCTA)
        TVDQPX=XTVD(IXVTA)*TAMOM(1,INUCTA)
        TVDQPY=XTVD(IXVTA)*TAMOM(2,INUCTA)
        TVDQPZ=XTVD(IXVTA)*TAMOM(3,INUCTA)
        TVDQE=XTVD(IXVTA)*TAMOM(4,INUCTA)
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
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
     *            PSQNX,PSQNY,PSQNZ,PSQNE,51)
      PSQPX=PSQNX
      PSQPY=PSQNY
      PSQPZ=PSQNZ
      PSQE=PSQNE      
      CALL CROMSC(PSAQPX,PSAQPY,PSAQPZ,PSAQE,RTIX,RTIY,RTIZ,
     *            PSAQNX,PSAQNY,PSAQNZ,PSAQNE,52)
      PSAQPX=PSAQNX
      PSAQPY=PSAQNY
      PSAQPZ=PSAQNZ
      PSAQE=PSAQNE      
C                                                ---------
C                                               j.r.6.5.93
C
C                     multiple scattering of VALENCE quark chain ends
C
      ITNU=IP+INUCTA
      RTIX=(VHKK(1,ITNU))*1.E12-BIMPAC*0.1
      RTIY=VHKK(2,ITNU)*1.E12
      RTIZ=VHKK(3,ITNU)*1.E12
      CALL CROMSC(TVQPX,TVQPY,TVQPZ,TVQE,RTIX,RTIY,RTIZ,
     *            TVQNX,TVQNY,TVQNZ,TVQNE,53)
      TVQPX=TVQNX
      TVQPY=TVQNY
      TVQPZ=TVQNZ
      TVQE=TVQNE      
      CALL CROMSC(TVDQPX,TVDQPY,TVDQPZ,TVDQE,RTIX,RTIY,RTIZ,
     *            TVDQNX,TVDQNY,TVDQNZ,TVDQNE,54)
      TVDQPX=TVDQNX
      TVDQPY=TVDQNY
      TVDQPZ=TVDQNZ
      TVDQE=TVDQNE     
      ENDIF 
C                                                ---------

C
C                                                j.r.10.5.93
       IF(IP.GE.0)GO TO 1779
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

C                                                ---------
C                                      changej.r.6.5.93
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
C
C***  SAMPLE PARTON-PT VALUES / DETERMINE PARTON 4-MOMENTA AND CHAIN MAS
C***                            IN THE REST FRAME DEFINED ABOVE
C
C                                                change j.r.6.5.93
C                                               _________________
        IKVALA=0
          IF(IPEV.GE.2) THEN
            WRITE(6,'(A,I5)') ' KKEVDV - IRDV13=',IRDV13
            WRITE(6,'(A/4(4E12.4/),2E12.4/2I5/4E12.4)')
     +      ' DV:  ...', PTXSQ1,PTYSQ1,PLQ1,EQ1,PTXSA1,PTYSA1,
     +      PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +      AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1
         ENDIF
        IKVALA=0
        NSELPT=1
        NSELPT=0
	IF(IP.EQ.1)NSELPT=1
        IF(NSELPT.EQ.1)CALL SELPT( PTXSQ1,PTYSQ1,PLQ1,
     +             EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1,
     +             PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +             PTXSQ2,PTYSQ2,PLQ2,EQ2,
     +             AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1,
     * PTTQ2,PTTA2,
     * NSELPT)
        IF(NSELPT.EQ.0)CALL SELPT4( PTXSQ1,PTYSQ1,PLQ1,
     +             EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1,
     +             PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +             PTXSQ2,PTYSQ2,PLQ2,EQ2,
     +             AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1,NSELPT)
          IF(IPEV.GE.2) THEN
            WRITE(6,'(A,I5)') ' KKEVDV - IRDV13=',IRDV13
            WRITE(6,'(A/4(4E12.4/),2E12.4/2I5/4E12.4)')
     +      ' DV:  ...', PTXSQ1,PTYSQ1,PLQ1,EQ1,PTXSA1,PTYSA1,
     +      PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +      AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1
          ENDIF
        IF (IPEV.GE.7) WRITE(6,'(A,I10)')
     +  'DV  IREJ ', IREJ
        IF (IREJ.EQ.1) THEN
          IRDV13=IRDV13 + 1
          IF(IPEV.GE.1) THEN
            WRITE(6,'(A,I5)') ' KKEVDV - IRDV13=',IRDV13
            WRITE(6,'(A/4(4E12.4/),2E12.4/2I5/4E12.4)')
     +      ' DV:  ...', PTXSQ1,PTYSQ1,PLQ1,EQ1,PTXSA1,PTYSA1,
     +      PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +      AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1
          ENDIF
                                                              GO TO 20
        ENDIF
C
C***  4-MOMENTA OF CHAINS IN THIS FRAME
C
        PTXCH1=PTXSQ1 + PTXSQ2
        PTYCH1=PTYSQ1 + PTYSQ2
        PTZCH1=PLQ1 + PLQ2
        ECH1=EQ1 + EQ2
        PTXCH2=PTXSA2 + PTXSA1
        PTYCH2=PTYSA2 + PTYSA1
        PTZCH2=PLAQ2 + PLAQ1
        ECH2=EAQ2 + EAQ1
        AMMM=SQRT((ECH1+ECH2)**2-(PTXCH1+PTXCH2)**2
     +            -(PTYCH1+PTYCH2)**2-(PTZCH1+PTZCH2)**2) 
C
C
        IF (IPEV.GE.6) WRITE(6,'(A,I10/A,5F12.5/A,5F12.5)')
     +  ' DV: IREJ ',IREJ, '     AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',
     +  AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1,
     +  '     AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ', AMCH2,PTXCH2,PTYCH2,
     +  PTZCH2,ECH2
 
C
C  REPLACE SMALL MASS CHAINS BY PSEUDOSCALAR OR VECTOR MESONS OR OCTETT
C                                              OR DECUPLETT BARYONS
C  FIRST FOR CHAIN 1  (PROJ SEA-diquark - TAR QUARK)
C
        CALL COBCMA(IPSQ(IXSPR),IPSQ2(IXSPR),ITVQ(IXVTA), IJNCH1,NNCH1,
     +  IREJ,AMCH1,AMCH1N,1)
C***                            MASS BELOW OCTETT BARYON MASS
        IF(IREJ.EQ.1) THEN
          IRDV11=IRDV11 + 1
          IF(IPEV.GE.1) THEN
            WRITE(6,'(A,I5)') ' KKEVDV - IRDV11=',IRDV11
            WRITE(6,'(A,6I5/6E12.4/2E12.4)') ' DV:', IPSQ(IXSPR),ITTV1
     +      (IXVTA),ITTV2(IXVTA),IJNCH1,NNCH1,IREJ, XPSQ(IXSPR),XPSAQ
     +      (IXSPR),XPSQCM,XPSACM, XTVQ(IXVTA),XTVD(IXVTA),AMCH1,AMCH1N
 
          ENDIF
                                                                 GOTO 20
        ENDIF
C                                 CORRECT KINEMATICS FOR CHAIN 1
C***                MOMENTUM CORRECTION FOR CHANGED MASS OF CHAIN 1
        IF(NNCH1.NE.0)THEN
           CALL CORMOM(AMCH1,AMCH2,AMCH1N,AMCH2N,
     +         PTXSQ1,PTYSQ1,PLQ1,EQ1,
     +         PTXSA1,PTYSA1,PLAQ1,EAQ1,
     +         PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +         PTXSQ2,PTYSQ2,PLQ2,EQ2,
     +         PTXCH1,PTYCH1,PTZCH1,ECH1, PTXCH2,PTYCH2,PTZCH2,ECH2,
     +         IREJ) 
        AMCH2=AMCH2N 
        ENDIF
C
        IF (IPEV.GE.2) WRITE(6,'(A,I10/A,5F12.5/A,5F12.5)')
     +  ' DV(2): IREJ ',IREJ, '     AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',
     +  AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1,
     +  '     AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ', AMCH2,PTXCH2,PTYCH2,
     +  PTZCH2,ECH2
        IF(IREJ.EQ.1)THEN
	  IF(IPEV.GE.1)WRITE(6,'(A)')' dv cormom rej.'
	  GO TO 20
	ENDIF
        IF (AMCH2 .LT.3.)THEN
	  IF(IPEV.GE.1)WRITE(6,'(A,F10.2)')' dv amch2',AMCH2
	  GO TO 20
        ENDIF
C
C
C  no test for chain 2 / mass constraint from DIQVS
        IJNCH2=0
        NNCH2=0
       QTXCH1=PTXCH1
       QTYCH1=PTYCH1
       QTZCH1=PTZCH1
       QECH1=ECH1
       QTXCH2=PTXCH2
       QTYCH2=PTYCH2
       QTZCH2=PTZCH2
       QECH2=ECH2
       PQDVA1(N,1)=PTXSQ1
       PQDVA1(N,2)=PTYSQ1
       PQDVA1(N,3)=PLQ1
       PQDVA1(N,4)=EQ1
       PQDVA2(N,1)=PTXSQ2
       PQDVA2(N,2)=PTYSQ2
       PQDVA2(N,3)=PLQ2
       PQDVA2(N,4)=EQ2
       PQDVB1(N,1)=PTXSA2
       PQDVB1(N,2)=PTYSA2
       PQDVB1(N,3)=PLAQ2
       PQDVB1(N,4)=EAQ2
       PQDVB2(N,1)=PTXSA1
       PQDVB2(N,2)=PTYSA1
       PQDVB2(N,3)=PLAQ1
       PQDVB2(N,4)=EAQ1
C-------------------

C
C                                      PUT D-V CHAIN ENDS INTO /HKKEVT/
C                                      MOMENTA IN NN-CMS
C                                      POSITION OF ORIGINAL NUCLEONS
C
****  keep for the moment the old s-v notations
C                                 FLAG FOR DV-CHAIN ENDS
C                                            PROJECTILE: ISTHKK=131
C                                            TARGET:     ISTHKK=122
C                                      FOR DV-CHAINS     ISTHKK=4
C
        IHKKPD=JHKKPS(IXSPR )
        IHKKPO=JHKKPS(IXSPR )-1
        IHKKTD=JHKKTV(IXVTA )
        IHKKTO=JHKKTV(IXVTA )-1
        IF (IPEV.GT.3)WRITE(6,1000)IXSPR,INUCPR,JNUCPR,IHKKPO,IHKKPD
 1000 FORMAT (' IXSPR,INUCPR,JNUCPR,IHKKPO,IHKKPD ',5I5)
        IF (IPEV.GT.3)WRITE(6,1010)IXVTA,INUCTA,JNUCTA,IHKKTO,IHKKTD
 1010 FORMAT (' IXVTA,INUCTA,JNUCTA,IHKKTO,IHKKTD ',5I5)
C                                     CHAIN 1 PROJECTILE SEA-diquark
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=131
        IDHKK(IHKK)=IDHKK(IHKKPO)
        JMOHKK(1,IHKK)=IHKKPO
        JMOHKK(2,IHKK)=JMOHKK(1,IHKKPO)
        JDAHKK(1,IHKK)=IHKK+2
        JDAHKK(2,IHKK)=IHKK+2
        PHKK(1,IHKK)=PQDVA1(N,1)
        PHKK(2,IHKK)=PQDVA1(N,2)
        PHKK(3,IHKK)=PQDVA1(N,3)
        PHKK(4,IHKK)=PQDVA1(N,4)
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
C                                     CHAIN 1 TARGET QUARK
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=122
        IDHKK(IHKK)=IDHKK(IHKKTD)
        JMOHKK(1,IHKK)=IHKKTD
        JMOHKK(2,IHKK)=JMOHKK(1,IHKKTD)
        JDAHKK(1,IHKK)=IHKK+1
        JDAHKK(2,IHKK)=IHKK+1
        PHKK(1,IHKK)=PQDVA2(N,1)
        PHKK(2,IHKK)=PQDVA2(N,2)
        PHKK(3,IHKK)=PQDVA2(N,3)
        PHKK(4,IHKK)=PQDVA2(N,4)
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
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
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
        VHKK(1,NHKK)= VHKK(1,NHKK-1)
        VHKK(2,NHKK)= VHKK(2,NHKK-1)
        VHKK(3,NHKK)= VHKK(3,NHKK-1)
        VHKK(4,NHKK)=VHKK(3,NHKK)/BETP-VHKK(3,NHKK-2)/BGAMP
        MHKKDV(N)=IHKK
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
C                                   CHAIN 2 projectile sea antidiquark
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=131
        IDHKK(IHKK)=IDHKK(IHKKPD)
        JMOHKK(1,IHKK)=IHKKPD
        JMOHKK(2,IHKK)=JMOHKK(1,IHKKPD)
        JDAHKK(1,IHKK)=IHKK+2
        JDAHKK(2,IHKK)=IHKK+2
        PHKK(1,IHKK)=PQDVB1(N,1)
        PHKK(2,IHKK)=PQDVB1(N,2)
        PHKK(3,IHKK)=PQDVB1(N,3)
        PHKK(4,IHKK)=PQDVB1(N,4)
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
 
C                                     CHAIN 2 TARGET diquark
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=122
        IDHKK(IHKK)=IDHKK(IHKKTO)
        JMOHKK(1,IHKK)=IHKKTO
        JMOHKK(2,IHKK)=JMOHKK(1,IHKKTO)
        JDAHKK(1,IHKK)=IHKK+1
        JDAHKK(2,IHKK)=IHKK+1
        PHKK(1,IHKK)=PQDVB2(N,1)
        PHKK(2,IHKK)=PQDVB2(N,2)
        PHKK(3,IHKK)=PQDVB2(N,3)
        PHKK(4,IHKK)=PQDVB2(N,4)
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
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
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
        MHKKDV(N)=IHKK
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
*     sea diquark pair!
C  AND PUT IT INTO THE HISTOGRAM
C
        AMCDV1(N)=AMCH1
        AMCDV2(N)=AMCH2
        GACDV1(N)=QECH1/AMCH1
        BGXDV1(N)=QTXCH1/AMCH1
        BGYDV1(N)=QTYCH1/AMCH1
        BGZDV1(N)=QTZCH1/AMCH1
        GACDV2(N)=QECH2/AMCH2
        BGXDV2(N)=QTXCH2/AMCH2
        BGYDV2(N)=QTYCH2/AMCH2
        BGZDV2(N)=QTZCH2/AMCH2
        NCHDV1(N)=NNCH1
        NCHDV2(N)=NNCH2
        IJCDV1(N)=IJNCH1
        IJCDV2(N)=IJNCH2
        IF (IPEV.GE.6) WRITE(6,'(A/I10,4F12.7,5I5/10X,4F12.6/10X,6F12.6,
     +4I5/8F15.5/                8F15.5)') ' DV / FINAL PRINT',N
C    +, XPSQ
C    +  (IXSPR),XPSAQ(IXSPR),XTVQ(IXVTA),XTVD(IXVTA), IPSQ(IXSPR),IPSAQ
C    +  (IXSPR), ITVQ(IXVTA),ITTV1(IXVTA),ITTV2(IXVTA), AMCDV1(N),AMCDV2
C    +  (N),GACDV1(N),GACDV2(N), BGXDV1(N),BGYDV1(N),BGZDV1(N), BGXDV2
C    +  (N),BGYDV2(N),BGZDV2(N), NCHDV1(N),NCHDV2(N),IJCDV1(N),IJCDV2
C    +  (N), (PQDVA1(N,JU),PQDVA2(N,JU),PQDVB1(N,JU), PQDVB2(N,JU),JU=1,
C    +  4)
   10 CONTINUE
      RETURN
C
   20 CONTINUE
C                                     EVENT REJECTED
C                                     START A NEW ONE
      IREJDV=1
      ISSQQ=IPSQ(IXSPR)
      JSSQQ=IPSQ2(IXSPR)
      IF(ISSQQ.EQ.3.AND.JSSQQ.EQ.3)THEN
        IDVRE(3)=IDVRE(3)+1
      ELSEIF(ISSQQ.EQ.3.OR.JSSQQ.EQ.3)THEN
        IDVRE(2)=IDVRE(2)+1
      ELSE
        IDVRE(1)=IDVRE(1)+1
      ENDIF
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     DEBUG SUBCHK
C     END DEBUG
      SUBROUTINE HADRDV
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C-------------------------
C
C                       hadronize sea diquark - valence CHAINS
C
C                       ADD GENERATED HADRONS TO /ALLPAR/
C                          STARTING AT (NAUX + 1)
C                       AND TO /HKKEVT/ STARTING AT (NHKK + 1)
C
C---------------------------------------------------------
      COMMON /ZSEA/ZSEAAV,ZSEASU,ANZSEA
       COMMON/POPCCK/PDBCK,PDBSE,PDBSEU,
     *  IJPOCK,IREJCK,ICK4,IHAD4,ICK6,IHAD6
     *,IREJSE,ISE4,ISE6,IREJS3,ISE43,ISE63,IREJS0
     *,IHADA4,IHADA6,IREJSA,ISEA4,ISEA6,IREJA3,
     *ISEA43,ISEA63,IREJAO
*KEEP,INTMX.
      PARAMETER (INTMX=2488,INTMD=252)
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
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
*KEEP,ABRDV.
      COMMON /ABRDV/ AMCDV1(248),AMCDV2(248),GACDV1(248),GACDV2(248),
     +BGXDV1(248),BGYDV1(248),BGZDV1(248), BGXDV2(248),BGYDV2(248),
     +BGZDV2(248), NCHDV1(248),NCHDV2(248),IJCDV1(248),IJCDV2(248),
     +PQDVA1(248,4),PQDVA2(248,4), PQDVB1(248,4),PQDVB2(248,4)
*KEEP,LOZUO.
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/ ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),ZUOST(INTMX),
     +INTLO(INTMX),INLOSS(INTMX)
C  /LOZUO/
C       INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C------------------
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
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
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
      COMMON /DIQSUM/NDVUU,NDVUS,NDVSS,NVDUU,NVDUS,NVDSS,
     *               NDSUU,NDSUS,NDSSS,NSDUU,NSDUS,NSDSS,
     *               NDZUU,NDZUS,NDZSS,NZDUU,NZDUS,NZDSS 
     *              ,NADVUU,NADVUS,NADVSS,NAVDUU,NAVDUS,NAVDSS,
     *               NADSUU,NADSUS,NADSSS,NASDUU,NASDUS,NASDSS,
     *               NADZUU,NADZUS,NADZSS,NAZDUU,NAZDUS,NAZDSS 
C---------------------
      DIMENSION POJ(4),PAT(4)
      DATA NCALDV /0/
      IF(IPHKK.GE.6)WRITE (6,'( A)') ' hadrdv'
C-----------------------------------------------------------------------
      NCALDV=NCALDV+1
      DO 50 I=1,NDV
C-----------------------drop recombined chain pairs
        IF(NCHDV1(I).EQ.99.AND.NCHDV2(I).EQ.99) GO TO 50
        IS1=INTDV1(I)
        IS2=INTDV2(I)
	IF(IPCO.GE.3)WRITE(6,*)' hadrdv I IS1,IS2 ',I,IS1,IS2
C
        IF (IPCO.GE.6) WRITE (6,1000) IPSQ(IS1),IPSAQ(IS1),ITVQ(IS2),
     +  ITTV1(IS2),ITTV2(IS2), AMCDV1(I),AMCDV2(I),GACDV1(I),GACDV2(I),
     +  BGXDV1(I),BGYDV1(I),BGZDV1(I), BGXDV2(I),BGYDV2(I),BGZDV2(I),
     +  NCHDV1(I),NCHDV2(I),IJCDV1(I),IJCDV2(I), PQDVA1(I,4),PQDVA2
     +  (I,4),PQDVB1(I,4),PQDVB2(I,4)
 1000 FORMAT(10X,5I5,10F9.2/10X,4I5,4F12.4)
C
C++++++++++++++++++++++++++++++    CHAIN 1:  diquark-quark   +++++++++++
        IFB1=IPSQ(IS1)
        IFB2=IPSQ2(IS1)
        IFB3=ITVQ(IS2)
        DO 10 J=1,4
          POJ(J)=PQDVA1(I,J)
          PAT(J)=PQDVA2(I,J)
   10   CONTINUE
        IF((NCHDV1(I).NE.0.OR.NCHDV2(I).NE.0).AND.IP.NE.1)
     &  CALL SAPTRE(AMCDV1(I),GACDV1(I),BGXDV1(I),BGYDV1(I),BGZDV1(I),
     &              AMCDV2(I),GACDV2(I),BGXDV2(I),BGYDV2(I),BGZDV2(I))
C----------------------------------------------------------------
C       WRITE (6,1244) POJ,PAT
C1244   FORMAT ('  D-V QUARK-DIQUARK POJ,PAT ',8E12.3)
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
C               I= number of valence chain
C               Projectile Nr ippp= IFROVP(INTVS1(I))
C          No of Glauber sea q at Projectile JIPP=JSSHS(IPP)
       IF(IPCO.GE.3)WRITE(6,*)' INTDV1(I) ',INTDV1(I)
C      IPPP = IFROVP(INTVS1(I))
C      IPPP = IFROVP(INTDV1(I))
       IF(INTDV1(I).GE.1)THEN
       IPPPP = IFROSP(INTDV1(I))
       ELSE
       IPPPP=0
       ENDIF
       IF(IPCO.GE.3)WRITE(6,*)' IPPP,IPPPP ',IPPP,IPPPP
       IF(IPPPP.GE.1)THEN
       JIPP=JSSHS(IPPPP)
       ELSE
       JIPP=0
       ENDIF
C      JIPP=1
       IF(IPCO.GE.3)WRITE(6,*)' JIPP ',JIPP
C       IF(NCHVS2(I).EQ.0)THEN
       IF(IPCO.GE.3)WRITE(6,'(A,3I5)')'HADRVS: I,IPPP,JIPP ',
     *                     I,IPPP,JIPP
C       ENDIF
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
         IF(IFB1.LE.2.AND.IFB2.LE.2)THEN
	   NDVUU=NDVUU+1
	 ELSEIF((IFB1.EQ.3.AND.IFB2.LE.2).OR.
     *  	 (IFB2.EQ.3.AND.IFB1.LE.2))THEN
	   NDVUS=NDVUS+1
	 ELSEIF(IFB1.EQ.3.AND.IFB2.EQ.3)THEN
	   NDVSS=NDVSS+1
	 ENDIF  
         IF (NCHDV1(I).NE.0)
     *   CALL HADJET(NHAD,AMCDV1(I),POJ,PAT,GACDV1(I),
     *  BGXDV1(I), BGYDV1
     +  (I),BGZDV1(I),IFB1,IFB2,IFB3,IFB4, IJCDV1(I),
     *  IJCDV1(I),6,NCHDV1
     +  (I),11)
C--------------------------------------------------------------------
        ZSEAWU=RNDM(BB)*2.D0*ZSEAAV
        RSEACK=FLOAT(JIPP)*PDBSE +ZSEAWU*PDBSEU
        IF(IPCO.GE.1)WRITE(6,*)'HADJSE JIPP,RSEACK,PDBSE 1 dpmnuc5',
     +  JIPP,RSEACK,PDBSE
C--------------------------------------------------------------------
        IF(NCHDV1(I).EQ.0)THEN
C---------------------------------------------------------------------
          IREJSS=5
          IF(RNDM(V).LE.RSEACK)THEN
            IREJSS=2
            IF(AMCDV1(I).GT.2.3D0)THEN
              IREJSS=0
              CALL HADJSE(NHAD,AMCDV1(I),POJ,PAT,GACDV1(I),BGXDV1(I),
     *        BGYDV1
     +        (I),BGZDV1(I),IFB1,IFB2,IFB3,IFB4, IJCDV1(I),IJCDV1(I),6,
     *        NCHDV1
     +        (I),6,IREJSS,IISSQQ)
              IF(IPCO.GE.1)WRITE(6,'(2A,I5,F10.3,I5)')'HADJSE JIPP,',
     *        'RSEACK,IREJSS 1 dpmnuc5  ',
     +        JIPP,RSEACK,IREJSS
            ENDIF
            IF(IREJSS.GE.1)THEN
              IF(IREJSS.EQ.1)IREJSE=IREJSE+1
              IF(IREJSS.EQ.3)IREJS3=IREJS3+1
              IF(IREJSS.EQ.2)IREJS0=IREJS0+1
              CALL HADJET(NHAD,AMCDV1(I),POJ,PAT,GACDV1(I),
     *        BGXDV1(I), BGYDV1
     +        (I),BGZDV1(I),IFB1,IFB2,IFB3,IFB4, IJCDV1(I),
     *        IJCDV1(I),6,NCHDV1
     +        (I),11)
              IHAD6=IHAD6+1
            ENDIF
            IF(IREJSS.EQ.0)THEN
              IF(IISSQQ.EQ.3)THEN
                ISE63=ISE63+1
              ELSE
                ISE6=ISE6+1
              ENDIF
            ENDIF 
          ELSE
            CALL HADJET(NHAD,AMCDV1(I),POJ,PAT,GACDV1(I),
     *      BGXDV1(I), BGYDV1
     +      (I),BGZDV1(I),IFB1,IFB2,IFB3,IFB4, IJCDV1(I),
     *      IJCDV1(I),6,NCHDV1
     +      (I),11)
            IHAD6=IHAD6+1
          ENDIF
        ENDIF
C--------------------------------------------------------------------
          ACOUDV=ACOUDV+1
          NHKKAU=NHKK+1
        DO 20 J=1,NHAD
          IF (NHKK.EQ.NMXHKK) THEN
            WRITE (6,'(A,2I5/A)') ' HADRDV: NHKK.EQ.NMXHKK ',NHKK,NMXHKK
            RETURN
          ENDIF
C         NHKK=NHKK+1
          IF (NHKK.EQ.NMXHKK)THEN
            WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
            RETURN
          ENDIF
C
          EHECC=SQRT(PXF(J)**2+PYF(J)**2+PZF(J)**2+AMF(J)**2)
          IF (ABS(EHECC-HEF(J)).GT.0.001) THEN
C           WRITE(6,'(2A/3I5,3E15.6)')
C    &      ' HADRDV / CHAIN 1 : CORRECT INCONSISTENT ENERGY ',
C    *      '  NCALDV, NHKK,NREF(J), HEF(J),EHECC, AMF(J)',
C    *      NCALDV, NHKK,NREF(J), HEF(J),EHECC, AMF(J)
            HEF(J)=EHECC
          ENDIF
          ANNDV=ANNDV+1
          EEDV=EEDV+HEF(J)
          PTDV=PTDV+SQRT(PXF(J)**2+PYF(J)**2)
C                               PUT NN-CMS HADRONS INTO /HKKEVT/
          ISTIST=1
          IF(IBARF(J).EQ.500)ISTIST=2
          CALL HKKFIL(ISTIST,MPDGHA(NREF(J)),MHKKDV(I)-3,0,
     *                PXF(J),PYF(J),PZF(J),HEF(J),NHKKAU,IORMO(J),13)
          IF(IDHKK(NHKK).EQ.99999) WRITE (6,1030)NHKK,NREF(J), IDHKK
     +    (NHKK)
          IF (IPHKK.GE.2) WRITE(6,1010) NHKK, ISTHKK(NHKK),IDHKK(NHKK),
     +    JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +    (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
   20   CONTINUE
C       IF(NHAD.GT.0) THEN
C         JDAHKK(1,IMOHKK)=NHKKAU
C         JDAHKK(2,IMOHKK)=NHKK
C       ENDIF
C+++++++++++++++++++++++++++++   CHAIN 2: adiquark - diquark  +++++++++
        IFB1=IPSAQ(IS1)
        IFB2=IPSAQ2(IS1)
        IFB3=ITTV1(IS2)
        IFB4=ITTV2(IS2)
        IFB1=IABS(IFB1)+6
        IFB2=IABS(IFB2)+6
        DO 30 J=1,4
          POJ(J)=PQDVB2(I,J)
          PAT(J)=PQDVB1(I,J)
   30   CONTINUE
C
         IF(IFB1.LE.8.AND.IFB2.LE.8)THEN
	   NADVUU=NADVUU+1
	 ELSEIF((IFB1.EQ.9.AND.IFB2.LE.8).OR.
     *  	 (IFB2.EQ.9.AND.IFB1.LE.8))THEN
	   NADVUS=NADVUS+1
	 ELSEIF(IFB1.EQ.9.AND.IFB2.EQ.9)THEN
	   NADVSS=NADVSS+1
	 ENDIF  
        CALL HADJET(NHAD,AMCDV2(I),POJ,PAT,GACDV2(I),
     *  BGXDV2(I), BGYDV2
     +  (I),BGZDV2(I),IFB1,IFB2,IFB3,IFB4, IJCDV2(I),
     *  IJCDV2(I),5,NCHDV2
     +  (I),12)
C                                   ADD HADRONS/RESONANCES INTO
C                                   COMMON /ALLPAR/ STARTING AT NAUX
        NHKKAU=NHKK+1
        DO 40 J=1,NHAD
          IF (NHKK.EQ.NMXHKK) THEN
            WRITE (6,'(A,2I5/A)') ' HADRDV: NHKK.EQ.NMXHKK ', NHKK,
     +      NMXHKK
            RETURN
          ENDIF
C         NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
C
          EHECC=SQRT(PXF(J)**2+PYF(J)**2+PZF(J)**2+AMF(J)**2)
        IF (ABS(EHECC-HEF(J)).GT.0.001) THEN
C             WRITE(6,'(2A/3I5,3E15.6)')
C    &            ' HADRDV / CHAIN 2 : CORRECT INCONSISTENT ENERGY ',
C    *            '  NCALDV, NHKK,NREF(J), HEF(J),EHECC, AMF(J)',
C    *            NCALDV, NHKK,NREF(J), HEF(J),EHECC, AMF(J)
          HEF(J)=EHECC
          ENDIF
          ANNDV=ANNDV+1
          EEDV=EEDV+HEF(J)
          PTDV=PTDV+SQRT(PXF(J)**2+PYF(J)**2)
C                               PUT NN-CMS HADRONS INTO /HKKEVT/
          ISTIST=1
          IF(IBARF(J).EQ.500)ISTIST=2
          CALL HKKFIL(ISTIST,MPDGHA(NREF(J)),MHKKDV(I),0,
     *                PXF(J),PYF(J),PZF(J),HEF(J),NHKKAU,IORMO(J),14)
          IF(IDHKK(NHKK).EQ.99999) WRITE (6,1030)NHKK,NREF(J), IDHKK
     +    (NHKK)
          IF (IPHKK.GE.2) WRITE(6,1010) NHKK, ISTHKK(NHKK),IDHKK(NHKK),
     +    JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +    (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
   40   CONTINUE
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
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE DIQVS(ECM,IPV,J,IREJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
*  define v-d chains (valence - sea diquark chains)
*  q-sqsq and qq-saqsaq chains instead of qq-sq and q-saq chains
      COMMON /ZSEA/ZSEAAV,ZSEASU,ANZSEA
       COMMON/POPCCK/PDBCK,PDBSE,PDBSEU,
     *  IJPOCK,IREJCK,ICK4,IHAD4,ICK6,IHAD6
     *,IREJSE,ISE4,ISE6,IREJS3,ISE43,ISE63,IREJS0
     *,IHADA4,IHADA6,IREJSA,ISEA4,ISEA6,IREJA3,
     *ISEA43,ISEA63,IREJAO
*KEEP,INTMX.
      PARAMETER (INTMX=2488,INTMD=252)
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
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
*KEEP,ABRVD.
      COMMON /ABRVD/ AMCVD1(248),AMCVD2(248),GACVD1(248),GACVD2(248),
     +BGXVD1(248),BGYVD1(248),BGZVD1(248), BGXVD2(248),BGYVD2(248),
     +BGZVD2(248), NCHVD1(248),NCHVD2(248),IJCVD1(248),IJCVD2(248),
     +PQVDA1(248,4),PQVDA2(248,4), PQVDB1(248,4),PQVDB2(248,4)
*KEND.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
      COMMON/SEASU3/SEASQ
      COMMON /XSEADI/ XSEACU,UNON,UNOM,UNOSEA, CVQ,CDQ,CSEA,SSMIMA,
     +SSMIMQ,VVMTHR
      COMMON /DIQREJ/IDIQRE(7),IDVRE(3),IVDRE(3),IDSRE(3),ISDRE(3),
     *IDZRE(3),IZDRE(3),IDIQRZ(7) 
      COMMON /DIQSUM/NDVUU,NDVUS,NDVSS,NVDUU,NVDUS,NVDSS,
     *               NDSUU,NDSUS,NDSSS,NSDUU,NSDUS,NSDSS,
     *               NDZUU,NDZUS,NDZSS,NZDUU,NZDUS,NZDSS 
     *              ,NADVUU,NADVUS,NADVSS,NAVDUU,NAVDUS,NAVDSS,
     *               NADSUU,NADSUS,NADSSS,NASDUU,NASDUS,NASDSS,
     *               NADZUU,NADZUS,NADZSS,NAZDUU,NAZDUS,NAZDSS 
C-----------------------------------------------------------------------
C     COMMON /PCHARM/PCCCC
      PARAMETER (UMMM=0.3D0)
      PARAMETER (SMMM=0.5D0)
      PARAMETER (CMMM=1.3D0)
      DATA PC/0.0001D0/
*KEND.
C----------
C
      DATA INICHA/0/
C----------------------------------------------------------------------
C                     Initialize Charm selection at soft chain ends
C
      IF(INICHA.EQ.0)THEN
        RX=8.D0
        X1=RX
        GM=2.140D0
        X2=UMMM
	BETOO=7.5D0
      ENDIF
      RX=8.D0
      X1=RX
      BETCHA=BETOO+1.3D0-LOG10(ECM)
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
      PC=PC1/7.0D0
      PU1=PU/(2*PU+PS+PC)
      PS1=PS/(2*PU+PS+PC)
      IF(INICHA.EQ.0)THEN
        INICHA=1
        WRITE(6,4567)PC,BETCHA,PU1,PS1,SEASQ
 4567   FORMAT(' Charm chain ends DIQVS: PC,BETCHA,PU,PS,SEASQ ',5F10.5)
      ENDIF
C----------------------------------------------------------------------
      IF(IPHKK.GE.6)WRITE (6,'( A,2I10)') ' diqvs IPV,J ',IPV,J
      IREJ=0
*  kinematics: is the mass of both chains big enough
*              to allow for fragmentation
      ITSQ2(J)=1.D0+RNDM(v1)*(2.D0+2.D0*SEASQ)
        RR=RNDM(V)
	IF(RR.LT.PC)ITSQ2(J)=4
C-----------------------------------------------------------------------
      ITSAQ2(J)=-ITSQ2(J)
C---------------------------------------------------j.r.29.4.94
C                                   x**1.5 distr for sea diquarks
C                         number of target nucleon
      INUCTA=IFROST(J)
C                          number of target diquark
      IITOT=ITOVT(INUCTA)
C                          diquark x
      XTDIQU=XTVD(IITOT)
C                          minimal value of diquark x
      XDTHR=CDQ/ECM
C
      XDFREE=XTDIQU-XDTHR
      XALL=XDFREE+XTSQ(J)+XTSAQ(J)-2.*XDTHR
      XDALT=XTVD(IITOT)
      XSALT=XTSQ(J)
      XAALT=XTSAQ(J)
      IF(XALL.GE.0.)THEN
        RR1=RNDM(V1)
        RR2=RNDM(V2)
        RR3=RNDM(V3)
        SR123=RR1+RR2+RR3
        DX1=RR1*XALL/SR123
        DX2=RR2*XALL/SR123
        DX3=RR3*XALL/SR123
        XTVD(IITOT)=XDTHR+DX1
        XTSQ(J)=XDTHR+DX2
        XTSAQ(J)=XDTHR+DX3
      ENDIF
C--------------------------------------------------------------

      AMVDQ1=XTSQ(J)*XPVQ(IPV)*ECM**2
      AMVDQ2=XTSAQ(J)*XPVD(IPV)*ECM**2
      IDIQRE(1)=IDIQRE(1)+1
      IF(ITSQ(J).GE.3.AND.ITSQ2(J).GE.3)THEN
      IDIQRE(2)=IDIQRE(2)+1
C       IF(AMVDQ2.LE.9.0.OR.AMVDQ1.LE.2.30) THEN
        IF(AMVDQ2.LE.17.0D0.OR.AMVDQ1.LE.6.60D0) THEN
          IREJ=1
           IDIQRE(3)=IDIQRE(3)+1
           IDIQRE(2)=IDIQRE(2)-1
           IDIQRE(1)=IDIQRE(1)-1
          XTVD(IITOT)=XDALT
          XTSQ(J)=XSALT
          XTSAQ(J)=XAALT
          RETURN
        ENDIF
      ELSEIF(ITSQ(J).GE.3.OR.ITSQ2(J).GE.3)THEN
      IDIQRE(4)=IDIQRE(4)+1
C       IF(AMVDQ2.LE.7.3D0.OR.AMVDQ1.LE.1.90D0) THEN
        IF(AMVDQ2.LE.13.6.OR.AMVDQ1.LE.5.80) THEN
          IREJ=1
           IDIQRE(5)=IDIQRE(5)+1
           IDIQRE(4)=IDIQRE(4)-1
           IDIQRE(1)=IDIQRE(1)-1
          XTVD(IITOT)=XDALT
          XTSQ(J)=XSALT
          XTSAQ(J)=XAALT
          RETURN
        ENDIF
      ELSE
      IDIQRE(6)=IDIQRE(6)+1
C       IF(AMVDQ2.LE.6.70.OR.AMVDQ1.LE.1.50) THEN
        IF(AMVDQ2.LE.12.40D0.OR.AMVDQ1.LE.3.9D0) THEN
          IREJ=1
           IDIQRE(7)=IDIQRE(7)+1
           IDIQRE(6)=IDIQRE(6)-1
           IDIQRE(1)=IDIQRE(1)-1
          XTVD(IITOT)=XDALT
          XTSQ(J)=XSALT
          XTSAQ(J)=XAALT
          RETURN
        ENDIF
      ENDIF
      NVD=NVD+1
c     WRITE(6,'(A/5X,3F10.3,3I5/5X,3F10.3)')
c    +    ' DIQVS: AMVDQ1, XTSQ, XPVQ, IPV,J, NVD/ AMVDQ2, XTSAQ, XPVD',
c    +  AMVDQ1,XTSQ(J),XPVQ(IPV),IPV,J, NVD, AMVDQ2,XTSAQ(J),XPVD(IPV)
      NCHVD1(NVD)=0
      NCHVD2(NVD)=0
      INTVD1(NVD)=IPV
      INTVD2(NVD)=J
      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C      DEBUG SUBCHK
C      END DEBUG
      SUBROUTINE KKEVVD(IREJVD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C------------------ treatment of valence - sea diquark CHAIN SYSTEMS
C
      COMMON /ZSEA/ZSEAAV,ZSEASU,ANZSEA
       COMMON/POPCCK/PDBCK,PDBSE,PDBSEU,
     *  IJPOCK,IREJCK,ICK4,IHAD4,ICK6,IHAD6
     *,IREJSE,ISE4,ISE6,IREJS3,ISE43,ISE63,IREJS0
     *,IHADA4,IHADA6,IREJSA,ISEA4,ISEA6,IREJA3,
     *ISEA43,ISEA63,IREJAO
*KEEP,INTMX.
      PARAMETER (INTMX=2488,INTMD=252)
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
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
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
*KEEP,ABRVD.
      COMMON /ABRVD/ AMCVD1(248),AMCVD2(248),GACVD1(248),GACVD2(248),
     +BGXVD1(248),BGYVD1(248),BGZVD1(248), BGXVD2(248),BGYVD2(248),
     +BGZVD2(248), NCHVD1(248),NCHVD2(248),IJCVD1(248),IJCVD2(248),
     +PQVDA1(248,4),PQVDA2(248,4), PQVDB1(248,4),PQVDB2(248,4)
*KEEP,DXQX.
C     INCLUDE (XQXQ)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON /DXQX/ XPVQ(248),XPVD(248),XTVQ(248),XTVD(248), XPSQ
     +(INTMX),XPSAQ(INTMX),XTSQ(INTMX),XTSAQ(INTMX)
     *                ,XPSU(248),XTSU(248)
     *                ,XPSUT(248),XTSUT(248)
      COMMON /DIQREJ/IDIQRE(7),IDVRE(3),IVDRE(3),IDSRE(3),ISDRE(3),
     *IDZRE(3),IZDRE(3),IDIQRZ(7) 
*KEEP,LOZUO.
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/ ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),ZUOST(INTMX),
     +INTLO(INTMX),INLOSS(INTMX)
C  /LOZUO/
C       INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C------------------
*KEEP,TRAFOP.
      COMMON /TRAFOP/ GAMP,BGAMP,BETP
*KEEP,NUCIMP.
      COMMON /NUCIMP/ PRMOM(5,248),TAMOM(5,248), PRMFEP,PRMFEN,TAMFEP,
     +TAMFEN, PREFEP,PREFEN,TAEFEP,TAEFEN, PREPOT(210),TAEPOT(210),
     +PREBIN,TAEBIN,FERMOD,ETACOU
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
      COMMON /DIQSUM/NDVUU,NDVUS,NDVSS,NVDUU,NVDUS,NVDSS,
     *               NDSUU,NDSUS,NDSSS,NSDUU,NSDUS,NSDSS,
     *               NDZUU,NDZUS,NDZSS,NZDUU,NZDUS,NZDSS 
     *              ,NADVUU,NADVUS,NADVSS,NAVDUU,NAVDUS,NAVDSS,
     *               NADSUU,NADSUS,NADSSS,NASDUU,NASDUS,NASDSS,
     *               NADZUU,NADZUS,NADZSS,NAZDUU,NAZDUS,NAZDSS 
*KEND.
C-------------------
      IF(IPHKK.GE.6)WRITE (6,'( A)') ' kkevvd'
      IREJVD=0
      DO 10 N=1,NVD
C---------------------------drop recombined chain pairs
        IF(NCHVD1(N).EQ.99.AND.NCHVD2(N).EQ.99)GO TO 10
C
C***            4-MOMENTA OF projectile QUARK-DIQUARK PAIRS IN NN-CMS
        IXVPR=INTVD1(N)
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
C***                 4-MOMENTA OF TARGET QUARK-DIQUARK PAIRS IN NN-CMS
        IXSTA=INTVD2(N)
        INUCTA=IFROST(IXSTA)
        JNUCTA=ITOVT(INUCTA)
*
        TSQPX=XTSQ(IXSTA)*TAMOM(1,INUCTA)
        TSQPY=XTSQ(IXSTA)*TAMOM(2,INUCTA)
        TSQPZ=XTSQ(IXSTA)*TAMOM(3,INUCTA)
        TSQE=XTSQ(IXSTA)*TAMOM(4,INUCTA)
        TSAQPX=XTSAQ(IXSTA)*TAMOM(1,INUCTA)
        TSAQPY=XTSAQ(IXSTA)*TAMOM(2,INUCTA)
        TSAQPZ=XTSAQ(IXSTA)*TAMOM(3,INUCTA)
        TSAQE=XTSAQ(IXSTA)*TAMOM(4,INUCTA)
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C                                               j.r.6.5.93
C
C                     multiple scattering of VALENCE quark chain ends
C
      IF(IT.GT.1)THEN
      ITNU=IP+INUCTA
      RTIX=(VHKK(1,ITNU))*1.E12-BIMPAC*0.1
      RTIY=VHKK(2,ITNU)*1.E12
      RTIZ=VHKK(3,ITNU)*1.E12
      CALL CROMSC(PVQPX,PVQPY,PVQPZ,PVQE,RTIX,RTIY,RTIZ,
     *            PVQNX,PVQNY,PVQNZ,PVQNE,55)
      PVQPX=PVQNX
      PVQPY=PVQNY
      PVQPZ=PVQNZ
      PVQE=PVQNE      
      CALL CROMSC(PVDQPX,PVDQPY,PVDQPZ,PVDQE,RTIX,RTIY,RTIZ,
     *            PVDQNX,PVDQNY,PVDQNZ,PVDQNE,56)
      PVDQPX=PVDQNX
      PVDQPY=PVDQNY
      PVDQPZ=PVDQNZ
      PVDQE=PVDQNE      
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
     *            TSQNX,TSQNY,TSQNZ,TSQNE,57)
      TSQPX=TSQNX
      TSQPY=TSQNY
      TSQPZ=TSQNZ
      TSQE=TSQNE      
      CALL CROMSC(TSAQPX,TSAQPY,TSAQPZ,TSAQE,RTIX,RTIY,RTIZ,
     *            TSAQNX,TSAQNY,TSAQNZ,TSAQNE,58)
      TSAQPX=TSAQNX
      TSAQPY=TSAQNY
      TSAQPZ=TSAQNZ
      TSAQE=TSAQNE  
      ENDIF    
C                                                ---------
C                                                ---------
C                                                j.r.10.5.93
       IF(IP.GE.0)GO TO 1779
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
C
C***  SAMPLE PARTON-PT VALUES / DETERMINE PARTON 4-MOMENTA AND CHAIN MAS
C***                            IN THE REST FRAME DEFINED ABOVE
C
        IKVALA=0
         IF(IPEV.GE.2) THEN
            WRITE(6,'(A,I5)') ' KKEVVD - IRVD13=',IRVD13
            WRITE(6,'(A/4(4E12.4/),2E12.4/2I5/4E12.4)')
     +      ' VD:   ...', PTXSQ1,PTYSQ1,PLQ1,EQ1,PTXSA1,PTYSA1,
     +      PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +      AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1
          ENDIF
        IKVALA=0
        NSELPT=1
        NSELPT=0
	IF(IP.EQ.1)NSELPT=1
        IF(NSELPT.EQ.1)CALL SELPT(PTXSQ1,PTYSQ1,PLQ1,
     +             EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1,
     +             PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +             PTXSQ2,PTYSQ2,PLQ2,EQ2,
     +             AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1,
     * PTTQ2,PTTA2,
     * NSELPT)
        IF(NSELPT.EQ.0)CALL SELPT4(PTXSQ1,PTYSQ1,PLQ1,
     +             EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1,
     +             PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +             PTXSQ2,PTYSQ2,PLQ2,EQ2,
     +             AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1,NSELPT)
         IF(IPEV.GE.2) THEN
            WRITE(6,'(A,I5)') ' KKEVVD - IRVD13=',IRVD13
            WRITE(6,'(A/4(4E12.4/),2E12.4/2I5/4E12.4)')
     +      ' VD:   ...', PTXSQ1,PTYSQ1,PLQ1,EQ1,PTXSA1,PTYSA1,
     +      PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +      AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1
             WRITE(6,'(A,I5)') ' KKEVVD - IRVD13=',IRVD13
            WRITE(6,'(A/4(4E12.4/),2E12.4/2I5/4E12.4)')
     +      ' VD:  amch1,amch2 ',
     +          AMCH1,AMCH2
          ENDIF

        IF (IPEV.GE.7) WRITE(6,'(A/5X,I10)')
     +  'VD   IREJ ', IREJ
        IF (IREJ.EQ.1) THEN
          IRVD13=IRVD13 + 1
          IF(IPEV.GE.1) THEN
            WRITE(6,'(A,I5)') ' KKEVVD - IRVD13=',IRVD13
            WRITE(6,'(A/4(4E12.4/),2E12.4/2I5/4E12.4)')
     +      ' VD:   ...', PTXSQ1,PTYSQ1,PLQ1,EQ1,PTXSA1,PTYSA1,
     +      PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +      AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1
 
          ENDIF
                                                                GO TO 20
        ENDIF
C
C***  4-MOMENTA OF CHAINS IN THIS FRAME
C
        PTXCH1=PTXSQ1 + PTXSQ2
        PTYCH1=PTYSQ1 + PTYSQ2
        PTZCH1=PLQ1 + PLQ2
        ECH1=EQ1 + EQ2
        PTXCH2=PTXSA2 + PTXSA1
        PTYCH2=PTYSA2 + PTYSA1
        PTZCH2=PLAQ2 + PLAQ1
        ECH2=EAQ2 + EAQ1
        AMMM=SQRT((ECH1+ECH2)**2-(PTXCH1+PTXCH2)**2
     +            -(PTYCH1+PTYCH2)**2-(PTZCH1+PTZCH2)**2) 
C
        IF (IPEV.GE.6) WRITE(6,'(A,I10/A,5F12.5/A,5F12.5)')
     +  ' VD: IREJ ', IREJ, '     AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',
     +  AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1,
     +  '     AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ', AMCH2,PTXCH2,PTYCH2,
     +  PTZCH2,ECH2
 
C
C  REPLACE SMALL MASS CHAINS BY PSEUDOSCALAR OR VECTOR MESONS OR OCTETT
C                                              OR DECUPLETT BARYONS
C  FIRST FOR CHAIN 1  (PROJ quark - tar sea-diquark)
C
        CALL COBCMA(ITSQ(IXSTA),ITSQ2(IXSTA),IPVQ(IXVPR), IJNCH1,NNCH1,
     +  IREJ,AMCH1,AMCH1N,1)
C***                            MASS BELOW OCTETT BARYON MASS
        IF(IREJ.EQ.1) THEN
          IRVD11=IRVD11 + 1
          IF(IPEV.GE.1) THEN
            WRITE(6,'(A,I5)') ' KKEVVD - IRVD11=',IRVD11
            WRITE(6,'(A,6I5/6E12.4/2E12.4)') ' VD:', IPVQ(IXVPR),ITSQ
     +      (IXSTA),ITSQ2(IXSTA),IJNCH1,NNCH1,IREJ, XPVQ(IXVPR),XPVD
     +      (IXVPR),XPVQCM,XPVDCM,
     + XTSQ(IXSTA),XTSAQ(IXSTA),AMCH1,AMCH1N
 
          ENDIF
                                                              GOTO 20
        ENDIF
C                                 CORRECT KINEMATICS FOR CHAIN 1
C***                MOMENTUM CORRECTION FOR CHANGED MASS OF CHAIN 1
        IF(NNCH1.NE.0)THEN
           CALL CORMOM(AMCH1,AMCH2,AMCH1N,AMCH2N,
     +         PTXSQ1,PTYSQ1,PLQ1,EQ1,
     +         PTXSA1,PTYSA1,PLAQ1,EAQ1,
     +         PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +         PTXSQ2,PTYSQ2,PLQ2,EQ2,
     +         PTXCH1,PTYCH1,PTZCH1,ECH1, PTXCH2,PTYCH2,PTZCH2,ECH2,
     +         IREJ)
        AMCH2=AMCH2N
        ENDIF
        IF(IREJ.EQ.1)THEN
	  IF(IPEV.GE.1)WRITE(6,'(A)')' vd CORMOM rej.'
	  GO TO 20
        ENDIF
C
        IF (IPEV.GE.6) WRITE(6,'(A,5F12.5,I10/A,5F12.5/A,5F12.5)')
     +  ' VD(2): AMMM,GAMMM,BGGGX,BGGGY,BGGGZ,IREJ ', AMMM,GAMMM,BGGGX,
     +  BGGGY,BGGGZ,IREJ, '     AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',
     +  AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1,
     +  '     AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ', AMCH2,PTXCH2,PTYCH2,
     +  PTZCH2,ECH2
C
C  no test for chain 2 / mass constraint from DIQVS
        IJNCH2=0
        NNCH2=0
        QTXCH1=PTXCH1
       QTYCH1=PTYCH1
       QTZCH1=PTZCH1
       QECH1=ECH1
       QTXCH2=PTXCH2
       QTYCH2=PTYCH2
       QTZCH2=PTZCH2
       QECH2=ECH2
       PQVDA1(N,1)=PTXSQ1
       PQVDA1(N,2)=PTYSQ1
       PQVDA1(N,3)=PLQ1
       PQVDA1(N,4)=EQ1
       PQVDA2(N,1)=PTXSQ2
       PQVDA2(N,2)=PTYSQ2
       PQVDA2(N,3)=PLQ2
       PQVDA2(N,4)=EQ2
       PQVDB1(N,1)=PTXSA2
       PQVDB1(N,2)=PTYSA2
       PQVDB1(N,3)=PLAQ2
       PQVDB1(N,4)=EAQ2
       PQVDB2(N,1)=PTXSA1
       PQVDB2(N,2)=PTYSA1
       PQVDB2(N,3)=PLAQ1
       PQVDB2(N,4)=EAQ1
C
C
C
C                                      PUT D-V CHAIN ENDS INTO /HKKEVT/
C                                      MOMENTA IN NN-CMS
C                                      POSITION OF ORIGINAL NUCLEONS
C
****  keep for the moment the old v-s notations
C                                 FLAG FOR VD-CHAIN ENDS
C                                            PROJECTILE: ISTHKK=121
C                                            TARGET:     ISTHKK=132
C                                      FOR VD-CHAINS     ISTHKK=5
C
        IHKKPD=JHKKPV(IXVPR )
        IHKKPO=JHKKPV(IXVPR )-1
        IHKKTD=JHKKTS(IXSTA )
        IHKKTO=JHKKTS(IXSTA )-1
        IF (IPEV.GT.3)WRITE(6,1000)IXVPR,INUCPR,JNUCPR,IHKKPO,IHKKPD
 1000 FORMAT (' IXVPR,INUCPR,JNUCPR,IHKKPO,IHKKPD ',5I5)
        IF (IPEV.GT.3)WRITE(6,1010)IXSTA,INUCTA,JNUCTA,IHKKTO,IHKKTD
 1010 FORMAT (' IXSTA,INUCTA,JNUCTA,IHKKTO,IHKKTD ',5I5)
C                                     CHAIN 1 PROJECTILE SEA-diquark
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=121
        IDHKK(IHKK)=IDHKK(IHKKPO)
        JMOHKK(1,IHKK)=IHKKPO
        JMOHKK(2,IHKK)=JMOHKK(1,IHKKPO)
        JDAHKK(1,IHKK)=IHKK+2
        JDAHKK(2,IHKK)=IHKK+2
        PHKK(1,IHKK)=PQVDA1(N,1)
        PHKK(2,IHKK)=PQVDA1(N,2)
        PHKK(3,IHKK)=PQVDA1(N,3)
        PHKK(4,IHKK)=PQVDA1(N,4)
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
C                                     CHAIN 1 TARGET QUARK
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=132
        IDHKK(IHKK)=IDHKK(IHKKTD)
        JMOHKK(1,IHKK)=IHKKTD
        JMOHKK(2,IHKK)=JMOHKK(1,IHKKTD)
        JDAHKK(1,IHKK)=IHKK+1
        JDAHKK(2,IHKK)=IHKK+1
        PHKK(1,IHKK)=PQVDA2(N,1)
        PHKK(2,IHKK)=PQVDA2(N,2)
        PHKK(3,IHKK)=PQVDA2(N,3)
        PHKK(4,IHKK)=PQVDA2(N,4)
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
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
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
        MHKKVD(N)=IHKK
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
C                                   CHAIN 2 projectile sea antidiquark
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=121
        IDHKK(IHKK)=IDHKK(IHKKPD)
        JMOHKK(1,IHKK)=IHKKPD
        JMOHKK(2,IHKK)=JMOHKK(1,IHKKPD)
        JDAHKK(1,IHKK)=IHKK+2
        JDAHKK(2,IHKK)=IHKK+2
        PHKK(1,IHKK)=PQVDB1(N,1)
        PHKK(2,IHKK)=PQVDB1(N,2)
        PHKK(3,IHKK)=PQVDB1(N,3)
        PHKK(4,IHKK)=PQVDB1(N,4)
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
 
C                                     CHAIN 2 TARGET diquark
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=132
        IDHKK(IHKK)=IDHKK(IHKKTO)
        JMOHKK(1,IHKK)=IHKKTO
        JMOHKK(2,IHKK)=JMOHKK(1,IHKKTO)
        JDAHKK(1,IHKK)=IHKK+1
        JDAHKK(2,IHKK)=IHKK+1
        PHKK(1,IHKK)=PQVDB2(N,1)
        PHKK(2,IHKK)=PQVDB2(N,2)
        PHKK(3,IHKK)=PQVDB2(N,3)
        PHKK(4,IHKK)=PQVDB2(N,4)
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
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
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
        MHKKVD(N)=IHKK
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
*     sea diquark pair!
C  AND PUT IT INTO THE HISTOGRAM
C
        AMCVD1(N)=AMCH1
        AMCVD2(N)=AMCH2
        GACVD1(N)=QECH1/AMCH1
        BGXVD1(N)=QTXCH1/AMCH1
        BGYVD1(N)=QTYCH1/AMCH1
        BGZVD1(N)=QTZCH1/AMCH1
        GACVD2(N)=QECH2/AMCH2
        BGXVD2(N)=QTXCH2/AMCH2
        BGYVD2(N)=QTYCH2/AMCH2
        BGZVD2(N)=QTZCH2/AMCH2
        NCHVD1(N)=NNCH1
        NCHVD2(N)=NNCH2
        IJCVD1(N)=IJNCH1
        IJCVD2(N)=IJNCH2
        IF (IPEV.GE.2) WRITE(6,'(A/I10,4F12.7,5I5/10X,4F12.6/10X,6F12.6,
     +4I5/8F15.5/8F15.5/2I5)') ' VD / FINAL PRINT',N
C    +, XPVQ
C    +  (IXVPR),XPVD(IXVPR),XTSQ(IXSTA),XTSAQ(IXSTA), IPVQ(IXVPR),IPPV1
C    +  (IXVPR),IPPV2(IXVPR),ITSQ(IXSTA),ITSAQ(IXSTA), AMCVD1(N),AMCVD2
C    +  (N),GACVD1(N),GACVD2(N), BGXVD1(N),BGYVD1(N),BGZVD1(N), BGXVD2
C    +  (N),BGYVD2(N),BGZVD2(N), NCHVD1(N),NCHVD2(N),IJCVD1(N),IJCVD2
C    +  (N), (PQVDA1(N,JU),PQVDA2(N,JU),PQVDB1(N,JU), PQVDB2(N,JU),JU=1,
C    +  4),
C    +  IXVPR,IXSTA
   10 CONTINUE
      RETURN
C
   20 CONTINUE
C                                     EVENT REJECTED
C                                     START A NEW ONE
      IREJVD=1
      ISSQQ=ITSQ(IXSTA)
      JSSQQ=ITSQ2(IXSTA)
      IF(ISSQQ.EQ.3.AND.JSSQQ.EQ.3)THEN
        IVDRE(3)=IVDRE(3)+1
      ELSEIF(ISSQQ.EQ.3.OR.JSSQQ.EQ.3)THEN
        IVDRE(2)=IVDRE(2)+1
      ELSE
        IVDRE(1)=IVDRE(1)+1
      ENDIF
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     DEBUG SUBCHK
C     END DEBUG
      SUBROUTINE HADRVD
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C-------------------------
C
C                       hadronize sea diquark - valence CHAINS
C
C                       ADD GENERATED HADRONS TO /ALLPAR/
C                          STARTING AT (NAUX + 1)
C                       AND TO /HKKEVT/ STARTING AT (NHKK + 1)
C
C---------------------------------------------------------
      COMMON /ZSEA/ZSEAAV,ZSEASU,ANZSEA
       COMMON/POPCCK/PDBCK,PDBSE,PDBSEU,
     *  IJPOCK,IREJCK,ICK4,IHAD4,ICK6,IHAD6
     *,IREJSE,ISE4,ISE6,IREJS3,ISE43,ISE63,IREJS0
     *,IHADA4,IHADA6,IREJSA,ISEA4,ISEA6,IREJA3,
     *ISEA43,ISEA63,IREJAO
*KEEP,INTMX.
      PARAMETER (INTMX=2488,INTMD=252)
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
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
*KEEP,ABRVD.
      COMMON /ABRVD/ AMCVD1(248),AMCVD2(248),GACVD1(248),GACVD2(248),
     +BGXVD1(248),BGYVD1(248),BGZVD1(248), BGXVD2(248),BGYVD2(248),
     +BGZVD2(248), NCHVD1(248),NCHVD2(248),IJCVD1(248),IJCVD2(248),
     +PQVDA1(248,4),PQVDA2(248,4), PQVDB1(248,4),PQVDB2(248,4)
*KEEP,LOZUO.
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/ ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),ZUOST(INTMX),
     +INTLO(INTMX),INLOSS(INTMX)
C  /LOZUO/
C       INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C------------------
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
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
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
      COMMON /DIQSUM/NDVUU,NDVUS,NDVSS,NVDUU,NVDUS,NVDSS,
     *               NDSUU,NDSUS,NDSSS,NSDUU,NSDUS,NSDSS,
     *               NDZUU,NDZUS,NDZSS,NZDUU,NZDUS,NZDSS 
     *              ,NADVUU,NADVUS,NADVSS,NAVDUU,NAVDUS,NAVDSS,
     *               NADSUU,NADSUS,NADSSS,NASDUU,NASDUS,NASDSS,
     *               NADZUU,NADZUS,NADZSS,NAZDUU,NAZDUS,NAZDSS 
C---------------------
      DIMENSION POJ(4),PAT(4)
      DATA NCALVD /0/
C-----------------------------------------------------------------------
      IF(IPHKK.GE.6)WRITE (6,'( A)') ' hadrvd'
      NCALVD=NCALVD+1
      DO 50 I=1,NVD
C-----------------------drop recombined chain pairs
        IF(NCHVD1(I).EQ.99.AND.NCHVD2(I).EQ.99) GO TO 50
        IS1=INTVD1(I)
        IS2=INTVD2(I)
C
        IF (IPCO.GE.6) WRITE (6,1000) IPSQ(IS1),IPSAQ(IS1),ITVQ(IS2),
     +  ITTV1(IS2),ITTV2(IS2), AMCVD1(I),AMCVD2(I),GACVD1(I),GACVD2(I),
     +  BGXVD1(I),BGYVD1(I),BGZVD1(I), BGXVD2(I),BGYVD2(I),BGZVD2(I),
     +  NCHVD1(I),NCHVD2(I),IJCVD1(I),IJCVD2(I), PQVDA1(I,4),PQVDA2
     +  (I,4),PQVDB1(I,4),PQVDB2(I,4)
 1000 FORMAT(10X,5I5,10F9.2/10X,4I5,4F12.4)
C
C++++++++++++++++++++++++++++++    CHAIN 1:  quark-diquark   +++++++++++
        IFB1=IPVQ(IS1)
        IFB2=ITSQ(IS2)
        IFB3=ITSQ2(IS2)
        DO 10 J=1,4
          POJ(J)=PQVDA1(I,J)
          PAT(J)=PQVDA2(I,J)
   10   CONTINUE
        IF((NCHVD1(I).NE.0.OR.NCHVD2(I).NE.0).AND.IP.NE.1)
     &  CALL SAPTRE(AMCVD1(I),GACVD1(I),BGXVD1(I),
     *  BGYVD1(I),BGZVD1(I),
     &              AMCVD2(I),GACVD2(I),BGXVD2(I),
     *  BGYVD2(I),BGZVD2(I))
C----------------------------------------------------------------
C----------------------------------------------------------------
C       WRITE (6,1244) POJ,PAT
C1244   FORMAT ('  V-D QUARK-DIQUARK POJ,PAT ',8E12.3)
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
C               I= number of valence chain
C               Target     Nr itt = IFROVT(INTSV2(I))
C          No of Glauber sea q at Target     JITT=JTSHS(ITT)
C      ITTT = IFROVT(INTSV2(I))
C      IF(INTVD2(I).GE.1)THEN
C      ITTT = IFROVT(INTVD2(I))
C      ELSE
       ITTT=0
C      ENDIF
C      IF(ITTT.GE.1)THEN
C     JITT=JTSHS(ITTT)
C      ELSE
       JITT=0
C      ENDIF
C       IF(NCHSV1(I).EQ.0)THEN
C      WRITE(6,'(A,3I5)')'HADRSV: I,ITTT,JITT ',
C    *                     I,ITTT,JITT
C       ENDIF
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
         IF(IFB2.LE.2.AND.IFB3.LE.2)THEN
	   NVDUU=NVDUU+1
	 ELSEIF((IFB2.EQ.3.AND.IFB3.LE.2).OR.
     *  	 (IFB3.EQ.3.AND.IFB2.LE.2))THEN
	   NVDUS=NVDUS+1
	 ELSEIF(IFB2.EQ.3.AND.IFB3.EQ.3)THEN
	   NVDSS=NVDSS+1
	 ENDIF  
        IF((NCHVD1(I).NE.0))THEN
          CALL HADJET(NHAD,AMCVD1(I),POJ,PAT,GACVD1(I),
     *    BGXVD1(I), BGYVD1
     +    (I),BGZVD1(I),IFB1,IFB2,IFB3,IFB4, IJCVD1(I),
     *    IJCVD1(I),4,NCHVD1
     +    (I),13)
        ENDIF
C-----------------------------------------------------------------
        AACK=FLOAT(ICK4)/FLOAT(ICK4+IHAD4+1)
        IF((NCHVD1(I).EQ.0))THEN
          ZSEAWU=RNDM(BB)*2.D0*ZSEAAV
        RSEACK=FLOAT(JITT)*PDBSE +ZSEAWU*PDBSEU
          IF(IPCO.GE.1)WRITE(6,'(2A,I5,2F10.3)')'HADJSE JITT,',
     *    'RSEACK,PDBSE 2 dpmnuc5 ',
     +    JITT,RSEACK,PDBSE
          IREJSS=5
          IF(RNDM(V).LE.RSEACK)THEN
            IREJSS=2
            IF(AMCVD1(I).GT.2.3D0)THEN
              IREJSS=0
              CALL HADJSE(NHAD,AMCVD1(I),POJ,PAT,GACVD1(I),
     *        BGXVD1(I),
     *        BGYVD1
     +        (I),BGZVD1(I),IFB1,IFB2,IFB3,IFB4, IJCVD1(I),
     *        IJCVD1(I),4,
     *        NCHVD1
     +        (I),3,IREJSS,IISSQQ)
              IF(IPCO.GE.1)WRITE(6,'(2A,I5,F10.3,I5)')'HADJSE JITT,',
     *        'RSEACK,IREJSS 2 dpmnuc5 ',
     +        JITT,RSEACK,IREJSS
            ENDIF
            IF(IREJSS.GE.1)THEN
              IF(IREJSS.EQ.1)IREJSE=IREJSE+1
              IF(IREJSS.EQ.3)IREJS3=IREJS3+1
              IF(IREJSS.EQ.2)IREJS0=IREJS0+1
              CALL HADJET(NHAD,AMCVD1(I),POJ,PAT,GACVD1(I),
     *        BGXVD1(I), BGYVD1
     +        (I),BGZVD1(I),IFB1,IFB2,IFB3,IFB4,
     *        IJCVD1(I),IJCVD1(I),4,NCHVD1
     +        (I),13)
              IHAD4=IHAD4+1
            ENDIF
            IF(IREJSS.EQ.0)THEN
              IF(IISSQQ.EQ.3)THEN
                ISE43=ISE43+1
              ELSE
                ISE4=ISE4+1
              ENDIF
            ENDIF
          ELSE
            CALL HADJET(NHAD,AMCVD1(I),POJ,PAT,GACVD1(I),
     *      BGXVD1(I), BGYVD1
     +      (I),BGZVD1(I),IFB1,IFB2,IFB3,IFB4,
     *      IJCVD1(I),IJCVD1(I),4,NCHVD1
     +      (I),13)
            IHAD4=IHAD4+1
          ENDIF
        ENDIF

C-----------------------------------------------------------------
        ACOUVD=ACOUVD+1
        NHKKAU=NHKK+1
        DO 20 J=1,NHAD
            IF (NHKK.EQ.NMXHKK) THEN
            WRITE (6,'(A,2I5/A)') ' HADRVD: NHKK.EQ.NMXHKK ',NHKK,NMXHKK
            RETURN
          ENDIF
C         NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
C
          EHECC=SQRT(PXF(J)**2+PYF(J)**2+PZF(J)**2+AMF(J)**2)
          IF (ABS(EHECC-HEF(J)).GT.0.001) THEN
C           WRITE(6,'(2A/3I5,3E15.6)')
C    &            ' HADRVD / CHAIN 1 : CORRECT INCONSISTENT ENERGY ',
C    *            '  NCALVD, NHKK,NREF(J), HEF(J),EHECC, AMF(J)',
C    *            NCALVD, NHKK,NREF(J), HEF(J),EHECC, AMF(J)
          HEF(J)=EHECC
          ENDIF
          ANNVD=ANNVD+1
          EEVD=EEVD+HEF(J)
          PTVD=PTVD+SQRT(PXF(J)**2+PYF(J)**2)
C                               PUT NN-CMS HADRONS INTO /HKKEVT/
          ISTIST=1
          IF(IBARF(J).EQ.500)ISTIST=2
          CALL HKKFIL(ISTIST,MPDGHA(NREF(J)),MHKKVD(I)-3,0,
     *                PXF(J),PYF(J),PZF(J),HEF(J),NHKKAU,IORMO(J),15)
          IF(IDHKK(NHKK).EQ.99999) WRITE (6,1030)NHKK,NREF(J), IDHKK
     +    (NHKK)
          IF (IPHKK.GE.2) WRITE(6,1010) NHKK, ISTHKK(NHKK),IDHKK(NHKK),
     +    JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +    (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
   20   CONTINUE
C       IF(NHAD.GT.0) THEN
C         JDAHKK(1,IMOHKK)=NHKKAU
C         JDAHKK(2,IMOHKK)=NHKK
C       ENDIF
C+++++++++++++++++++++++++++++   CHAIN 2: diquark - adiquark  +++++++++
        IFB1=IPPV1(IS1)
        IFB2=IPPV2(IS1)
        IFB3=ITSAQ(IS2)
        IFB4=ITSAQ2(IS2)
        IFB3=IABS(IFB3)+6
        IFB4=IABS(IFB4)+6
        DO 30 J=1,4
          POJ(J)=PQVDB2(I,J)
          PAT(J)=PQVDB1(I,J)
   30   CONTINUE
C
        IF(AMCVD2(I).LT.2.3)THEN
          WRITE(6,'(A,F10.2,I5)')' HADRVD AMCVD2(I), I ',
     *    AMCVD2(I),I
          RETURN
        ENDIF
         IF(IFB3.LE.8.AND.IFB4.LE.8)THEN
	   NAVDUU=NAVDUU+1
	 ELSEIF((IFB3.EQ.9.AND.IFB4.LE.8).OR.
     *  	 (IFB4.EQ.9.AND.IFB3.LE.8))THEN
	   NAVDUS=NAVDUS+1
	 ELSEIF(IFB3.EQ.9.AND.IFB4.EQ.9)THEN
	   NAVDSS=NAVDSS+1
	 ENDIF  
        CALL HADJET(NHAD,AMCVD2(I),POJ,PAT,GACVD2(I),
     *  BGXVD2(I), BGYVD2
     +  (I),BGZVD2(I),IFB1,IFB2,IFB3,IFB4, IJCVD2(I),
     *  IJCVD2(I),5,NCHVD2
     +  (I),14)
C                                   ADD HADRONS/RESONANCES INTO
C                                   COMMON /ALLPAR/ STARTING AT NAUX
        NHKKAU=NHKK+1
        DO 40 J=1,NHAD
          IF (NHKK.EQ.NMXHKK) THEN
            WRITE (6,'(A,2I5/A)') ' HADRVD: NHKK.EQ.NMXHKK ',
     *      NHKK,
     +      NMXHKK
            RETURN
          ENDIF
C         NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
C
          EHECC=SQRT(PXF(J)**2+PYF(J)**2+PZF(J)**2+AMF(J)**2)
        IF (ABS(EHECC-HEF(J)).GT.0.001) THEN
C             WRITE(6,'(2A/3I5,3E15.6)')
C    &            ' HADRVD / CHAIN 2 : CORRECT INCONSISTENT ENERGY ',
C    *            '  NCALVD, NHKK,NREF(J), HEF(J),EHECC, AMF(J)',
C    *            NCALVD, NHKK,NREF(J), HEF(J),EHECC, AMF(J)
          HEF(J)=EHECC
          ENDIF
          ANNVD=ANNVD+1
          EEVD=EEVD+HEF(J)
          PTVD=PTVD+SQRT(PXF(J)**2+PYF(J)**2)
C                               PUT NN-CMS HADRONS INTO /HKKEVT/
          ISTIST=1
          IF(IBARF(J).EQ.500)ISTIST=2
          CALL HKKFIL(ISTIST,MPDGHA(NREF(J)),MHKKVD(I),0,
     *                PXF(J),PYF(J),PZF(J),HEF(J),NHKKAU,IORMO(J),16)
          IF(IDHKK(NHKK).EQ.99999) WRITE (6,1030)NHKK,NREF(J), IDHKK
     +    (NHKK)
          IF (IPHKK.GE.2) WRITE(6,1010) NHKK, ISTHKK(NHKK),IDHKK(NHKK),
     +    JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +    (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
   40   CONTINUE
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
C
C ---------------------------------------------------------------
C
      SUBROUTINE DIQDSS(ECM,ITS,IPS,IREJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
*  define d-s chains (sea diquark - sea chains)
*  sqsq-sq and saqsaq-saq chains instead of q-aq and aq-q chains
      COMMON /ZSEA/ZSEAAV,ZSEASU,ANZSEA
       COMMON/POPCCK/PDBCK,PDBSE,PDBSEU,
     *  IJPOCK,IREJCK,ICK4,IHAD4,ICK6,IHAD6
     *,IREJSE,ISE4,ISE6,IREJS3,ISE43,ISE63,IREJS0
     *,IHADA4,IHADA6,IREJSA,ISEA4,ISEA6,IREJA3,
     *ISEA43,ISEA63,IREJAO
*KEEP,INTMX.
      PARAMETER (INTMX=2488,INTMD=252)
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
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
*KEEP,ABRDS.
      COMMON /ABRDS/ AMCDS1(248),AMCDS2(248),GACDS1(248),GACDS2(248),
     +BGXDS1(248),BGYDS1(248),BGZDS1(248), BGXDS2(248),BGYDS2(248),
     +BGZDS2(248), NCHDS1(248),NCHDS2(248),IJCDS1(248),IJCDS2(248),
     +PQDSA1(248,4),PQDSA2(248,4), PQDSB1(248,4),PQDSB2(248,4)
*KEND.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
      COMMON/SEASU3/SEASQ
      COMMON /XSEADI/ XSEACU,UNON,UNOM,UNOSEA, CVQ,CDQ,CSEA,SSMIMA,
     +SSMIMQ,VVMTHR
      COMMON /DIQREJ/IDIQRE(7),IDVRE(3),IVDRE(3),IDSRE(3),ISDRE(3),
     *IDZRE(3),IZDRE(3),IDIQRZ(7) 
      COMMON /DIQSUM/NDVUU,NDVUS,NDVSS,NVDUU,NVDUS,NVDSS,
     *               NDSUU,NDSUS,NDSSS,NSDUU,NSDUS,NSDSS,
     *               NDZUU,NDZUS,NDZSS,NZDUU,NZDUS,NZDSS 
     *              ,NADVUU,NADVUS,NADVSS,NAVDUU,NAVDUS,NAVDSS,
     *               NADSUU,NADSUS,NADSSS,NASDUU,NASDUS,NASDSS,
     *               NADZUU,NADZUS,NADZSS,NAZDUU,NAZDUS,NAZDSS 
C-----------------------------------------------------------------
C     COMMON /PCHARM/PCCCC
      PARAMETER (UMMM=0.3D0)
      PARAMETER (SMMM=0.5D0)
      PARAMETER (CMMM=1.3D0)
      DATA PC/0.0001D0/
*KEND.
C----------
C
      DATA INICHA/0/
C----------------------------------------------------------------------
C                     Initialize Charm selection at soft chain ends
C
      IF(INICHA.EQ.0)THEN
        RX=8.D0
        X1=RX
        GM=2.140D0
        X2=UMMM
	BETOO=7.5D0
      ENDIF
      RX=8.D0
      X1=RX
      BETCHA=BETOO+1.3D0-LOG10(ECM)
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
      PC=PC1/7.0D0
      PU1=PU/(2*PU+PS+PC)
      PS1=PS/(2*PU+PS+PC)
      IF(INICHA.EQ.0)THEN
        INICHA=1
        WRITE(6,4567)PC,BETCHA,PU1,PS1,SEASQ
 4567   FORMAT(' Charm chain ends DIQDSS: PC,BETCHA,PU,PS,SEASQ',5F10.5)
      ENDIF
C----------------------------------------------------------------------
      IF(IPHKK.GE.6)WRITE (6,'( A)') ' diqdss'
      IREJ=0
*  kinematics: is the mass of the adiquark-diquark chain big enough
*              to allow for fragmentation
      IPSQ2(IPS)=1.D0+RNDM(v1)*(2.D0+2.D0*SEASQ)
      RR=RNDM(V)
      IF(RR.LT.PC)IPSQ2(IPS)=4
C-----------------------------------------------------------------
      IPSAQ2(IPS)=-IPSQ2(IPS)
C---------------------------------------------------j.r.29.4.94
C                                   x**1.5 distr for sea diquarks
C                         number of projectile nucleon
      INUCPR=IFROSP(IPS)
C                          number of projectile diquark
      IITOP=ITOVP(INUCPR)
C                          diquark x
      XPDIQU=XPVD(IITOP)
C                          minimal value of diquark x
      XDTHR=CDQ/ECM
C
      XDFREE=XPDIQU-XDTHR
      XALL=XDFREE+XPSQ(IPS)+XPSAQ(IPS)-2.*XDTHR
      XDALT=XPVD(IITOP)
      XSALT=XPSQ(IPS)
      XAALT=XPSAQ(IPS) 
      IF(XALL.GE.0.)THEN
        RR1=RNDM(V1)
        RR2=RNDM(V2)
        RR3=RNDM(V3)
        SR123=RR1+RR2+RR3
        DX1=RR1*XALL/SR123
        DX2=RR2*XALL/SR123
        DX3=RR3*XALL/SR123
        XPVD(IITOP)=XDTHR+DX1
        XPSQ(IPS)=XDTHR+DX2
        XPSAQ(IPS)=XDTHR+DX3
      ENDIF
C--------------------------------------------------------------
      AMDSQ1=XPSQ(IPS)*XTSQ(ITS)*ECM**2
      AMDSQ2=XPSAQ(IPS)*XTSAQ(ITS)*ECM**2
      IDIQRE(1)=IDIQRE(1)+1
      IF(IPSQ(IPS).GE.3.AND.IPSQ2(IPS).GE.3)THEN
        IDIQRE(2)=IDIQRE(2)+1
C       IF(AMDSQ2.LE.2.3.OR.AMDSQ1.LE.2.30) THEN
        IF(AMDSQ2.LE.6.6D0.OR.AMDSQ1.LE.6.60D0) THEN
          IREJ=1
          IDIQRE(3)=IDIQRE(3)+1
          IDIQRE(2)=IDIQRE(2)-1
          IDIQRE(1)=IDIQRE(1)-1
          XPVD(IITOP)=XDALT
          XPSQ(IPS)=XSALT
          XPSAQ(IPS)=XAALT
          RETURN
        ENDIF
      ELSEIF(IPSQ(IPS).GE.3.OR.IPSQ2(IPS).GE.3)THEN
        IDIQRE(4)=IDIQRE(4)+1
C       IF(AMDSQ2.LE.1.9.OR.AMDSQ1.LE.1.90) THEN
        IF(AMDSQ2.LE.5.8D0.OR.AMDSQ1.LE.5.80D0) THEN
          IREJ=1
          IDIQRE(5)=IDIQRE(5)+1
          IDIQRE(4)=IDIQRE(4)-1
          IDIQRE(1)=IDIQRE(1)-1
          XPVD(IITOP)=XDALT
          XPSQ(IPS)=XSALT
          XPSAQ(IPS)=XAALT
          RETURN
        ENDIF
      ELSE
        IDIQRE(6)=IDIQRE(6)+1
C       IF(AMDSQ2.LE.1.50.OR.AMDSQ1.LE.1.50) THEN
        IF(AMDSQ2.LE.3.9D0.OR.AMDSQ1.LE.3.9D0) THEN
          IREJ=1
          IDIQRE(7)=IDIQRE(7)+1
          IDIQRE(6)=IDIQRE(6)-1
          IDIQRE(1)=IDIQRE(1)-1
          XPVD(IITOP)=XDALT
          XPSQ(IPS)=XSALT
          XPSAQ(IPS)=XAALT
          RETURN
        ENDIF
      ENDIF
      NDS=NDS+1
      NCHDS1(NDS)=0
      NCHDS2(NDS)=0
      INTDS1(NDS)=IPS
      INTDS2(NDS)=ITS
C     WRITE(6,'(A/5X,3F10.3,3I5/5X,3F10.3)')
C    +' DIQDSS: AMDSQ1, XTSQ, XPSQ, ITS,IPS, NDS/ AMDSQ2, XTSAQ, XPSAQ',
C    +  AMDSQ1,XTSQ(ITS),XPSQ(IPS),ITS,IPS, NDS, AMDSQ2,
C    +XTSAQ(ITS),XPSAQ(IPS)
      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C      DEBUG SUBCHK
C      END DEBUG
      SUBROUTINE KKEVDS(IREJDS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C------------------ treatment of sea diquark - sea CHAIN SYSTEMS
C
      COMMON /ZSEA/ZSEAAV,ZSEASU,ANZSEA
       COMMON/POPCCK/PDBCK,PDBSE,PDBSEU,
     *  IJPOCK,IREJCK,ICK4,IHAD4,ICK6,IHAD6
     *,IREJSE,ISE4,ISE6,IREJS3,ISE43,ISE63,IREJS0
     *,IHADA4,IHADA6,IREJSA,ISEA4,ISEA6,IREJA3,
     *ISEA43,ISEA63,IREJAO
*KEEP,INTMX.
      PARAMETER (INTMX=2488,INTMD=252)
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
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
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
*KEEP,ABRDS.
      COMMON /ABRDS/ AMCDS1(248),AMCDS2(248),GACDS1(248),GACDS2(248),
     +BGXDS1(248),BGYDS1(248),BGZDS1(248), BGXDS2(248),BGYDS2(248),
     +BGZDS2(248), NCHDS1(248),NCHDS2(248),IJCDS1(248),IJCDS2(248),
     +PQDSA1(248,4),PQDSA2(248,4), PQDSB1(248,4),PQDSB2(248,4)
      COMMON /DIQREJ/IDIQRE(7),IDVRE(3),IVDRE(3),IDSRE(3),ISDRE(3),
     *IDZRE(3),IZDRE(3),IDIQRZ(7) 
*KEEP,DXQX.
C     INCLUDE (XQXQ)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON /DXQX/ XPVQ(248),XPVD(248),XTVQ(248),XTVD(248), XPSQ
     +(INTMX),XPSAQ(INTMX),XTSQ(INTMX),XTSAQ(INTMX)
     *                ,XPSU(248),XTSU(248)
     *                ,XPSUT(248),XTSUT(248)
*KEEP,LOZUO.
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/ ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),ZUOST(INTMX),
     +INTLO(INTMX),INLOSS(INTMX)
C  /LOZUO/
C       INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C------------------
*KEEP,TRAFOP.
      COMMON /TRAFOP/ GAMP,BGAMP,BETP
*KEEP,NUCIMP.
      COMMON /NUCIMP/ PRMOM(5,248),TAMOM(5,248), PRMFEP,PRMFEN,TAMFEP,
     +TAMFEN, PREFEP,PREFEN,TAEFEP,TAEFEN, PREPOT(210),TAEPOT(210),
     +PREBIN,TAEBIN,FERMOD,ETACOU
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
      COMMON /DIQSUM/NDVUU,NDVUS,NDVSS,NVDUU,NVDUS,NVDSS,
     *               NDSUU,NDSUS,NDSSS,NSDUU,NSDUS,NSDSS,
     *               NDZUU,NDZUS,NDZSS,NZDUU,NZDUS,NZDSS 
     *              ,NADVUU,NADVUS,NADVSS,NAVDUU,NAVDUS,NAVDSS,
     *               NADSUU,NADSUS,NADSSS,NASDUU,NASDUS,NASDSS,
     *               NADZUU,NADZUS,NADZSS,NAZDUU,NAZDUS,NAZDSS 
*KEND.
C-------------------
      IF(IPHKK.GE.6)WRITE (6,'( A)') ' kkevds'
      IREJDS=0
      DO 10 N=1,NDS
C---------------------------drop recombined chain pairs
        IF(NCHDS1(N).EQ.99.AND.NCHDS2(N).EQ.99)GO TO 10
C
C***                 4-MOMENTA OF PROJECTILE SEA-QUARK PAIRS IN NN-CMS
        IF(IPHKK.GE.7)WRITE(6,'(A,2I10)')' KKEVDS N,NDS',N,NDS
        IXSPR=INTDS1(N)
        IF(IPHKK.GE.7)WRITE(6,'(A,2I10)')' KKEVDS N,IXSPR',N,IXSPR
        INUCPR=IFROSP(IXSPR)
        JNUCPR=ITOVP(INUCPR)
        IF(IPHKK.GE.7)WRITE(6,'(A,2I10)')' KKEVDS INUCPR,JNUCPR',
     + INUCPR,JNUCPR
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
C***                 4-MOMENTA OF TARGET QUARK-AQUARK PAIRS IN NN-CMS
        IXSTA=INTDS2(N)
        IF(IPHKK.GE.7)WRITE(6,'(A,2I10)')' KKEVDS N,IXSTA',N,IXSTA
        INUCTA=IFROST(IXSTA)
        JNUCTA=ITOVT(INUCTA)
        IF(IPHKK.GE.7)WRITE(6,'(A,2I10)')' KKEVDS INUCTA,JNUCTA',
     + INUCTA,JNUCTA
C
        TSQPX=XTSQ(IXSTA)*TAMOM(1,INUCTA)
        TSQPY=XTSQ(IXSTA)*TAMOM(2,INUCTA)
        TSQPZ=XTSQ(IXSTA)*TAMOM(3,INUCTA)
        TSQE=XTSQ(IXSTA)*TAMOM(4,INUCTA)
        TSDQPX=XTSAQ(IXSTA)*TAMOM(1,INUCTA)
        TSDQPY=XTSAQ(IXSTA)*TAMOM(2,INUCTA)
        TSDQPZ=XTSAQ(IXSTA)*TAMOM(3,INUCTA)
        TSDQE=XTSAQ(IXSTA)*TAMOM(4,INUCTA)
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
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
     *            PSQNX,PSQNY,PSQNZ,PSQNE,59)
      PSQPX=PSQNX
      PSQPY=PSQNY
      PSQPZ=PSQNZ
      PSQE=PSQNE      
      CALL CROMSC(PSAQPX,PSAQPY,PSAQPZ,PSAQE,RTIX,RTIY,RTIZ,
     *            PSAQNX,PSAQNY,PSAQNZ,PSAQNE,60)
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
     *            TSQNX,TSQNY,TSQNZ,TSQNE,61)
      TSQPX=TSQNX
      TSQPY=TSQNY
      TSQPZ=TSQNZ
      TSQE=TSQNE      
      CALL CROMSC(TSDQPX,TSDQPY,TSDQPZ,TSDQE,RTIX,RTIY,RTIZ,
     *            TSDQNX,TSDQNY,TSDQNZ,TSDQNE,62)
      TSDQPX=TSDQNX
      TSDQPY=TSDQNY
      TSDQPZ=TSDQNZ
      TSDQE=TSDQNE   
      ENDIF    
C                                                ---------
C                                                j.r.10.5.93
       IF(IP.GE.0)GO TO 1779
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
        TDQPZ2=TSDQE**2-TSDQPX**2-TSDQPY**2
        IF(TDQPZ2.GE.0.)THEN
          TSDQPZ=-SQRT(TDQPZ2)
        ELSE
          TSDQPX=0.
          TSDQPY=0.
          TSDQPZ=TSDQE
        ENDIF
 1779  CONTINUE
C                                                ---------
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
        PTXSA2=TSDQPX
        PTYSQ1=PSQPY
        PTYSA1=PSAQPY
        PTYSQ2=TSQPY
        PTYSA2=TSDQPY
        PLQ1=PSQPZ
        PLAQ1=PSAQPZ
        PLQ2=TSQPZ
        PLAQ2=TSDQPZ
        EQ1=PSQE
        EAQ1=PSAQE
        EQ2=TSQE
        EAQ2=TSDQE
C                                       ---------------
C
        IKVALA=0 
         IF(IPEV.GE.2) THEN
            WRITE(6,'(A,I5)') ' KKEVDS - IRDS13=',IRDS13
            WRITE(6,'(A/4(4E12.4/),2E12.4/2I5/4E12.4)')
     +      ' DS:  ...',  PTXSQ1,PTYSQ1,PLQ1,EQ1,PTXSA1,PTYSA1,
     +      PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +      AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1
          ENDIF
        IKVALA=0
        NSELPT=1
        CALL SELPT( PTXSQ1,PTYSQ1,PLQ1,
     +             EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1,
     +             PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +             PTXSQ2,PTYSQ2,PLQ2,EQ2,
     +             AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1,
     * PTTQ2,PTTA2,
     * NSELPT) 
         IF(IPEV.GE.2) THEN
            WRITE(6,'(A,I5)') ' KKEVDS - IRDS13=',IRDS13
            WRITE(6,'(A/4(4E12.4/),2E12.4/2I5/4E12.4)')
     +      ' DS:  ...',  PTXSQ1,PTYSQ1,PLQ1,EQ1,PTXSA1,PTYSA1,
     +      PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +      AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1
          ENDIF

        IF (IPEV.GE.7) WRITE(6,'(A/5X,I10)')
     +  'DS   IREJ ', IREJ
        IF (IREJ.EQ.1) THEN
          IRDS13=IRDS13 + 1
          IF(IPEV.GE.2) THEN
            WRITE(6,'(A,I5)') ' KKEVDS - IRDS13=',IRDS13
            WRITE(6,'(A/4(4E12.4/),2E12.4/2I5/4E12.4)')
     +      ' DS:  ...',  PTXSQ1,PTYSQ1,PLQ1,EQ1,PTXSA1,PTYSA1,
     +      PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +      AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1
          ENDIF
                                                                GO TO 11
        ENDIF
C
C***  4-MOMENTA OF CHAINS IN THIS FRAME
C
        PTXCH1=PTXSQ1 + PTXSQ2
        PTYCH1=PTYSQ1 + PTYSQ2
        PTZCH1=PLQ1 + PLQ2
        ECH1=EQ1 + EQ2
        PTXCH2=PTXSA2 + PTXSA1
        PTYCH2=PTYSA2 + PTYSA1
        PTZCH2=PLAQ2 + PLAQ1
        ECH2=EAQ2 + EAQ1
        AMMM=SQRT((ECH1+ECH2)**2-(PTXCH1+PTXCH2)**2
     +            -(PTYCH1+PTYCH2)**2-(PTZCH1+PTZCH2)**2) 
C
C
        IF (IPEV.GE.6) WRITE(6,'(A,I10/A,5F12.5/A,5F12.5)')
     +  ' DS: IREJ ',IREJ, '     AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',
     +  AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1,
     +  '     AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ', AMCH2,PTXCH2,PTYCH2,
     +  PTZCH2,ECH2
C
C  REPLACE SMALL MASS CHAINS BY PSEUDOSCALAR OR VECTOR MESONS OR OCTETT
C                                              OR DECUPLETT BARYONS
C  FIRST FOR CHAIN 1  (PROJ SEA-diquark - TAR QUARK)
C
        CALL COBCMA(IPSQ(IXSPR),IPSQ2(IXSPR),ITSQ(IXSTA), IJNCH1,NNCH1,
     +  IREJ,AMCH1,AMCH1N,1)
C***                            MASS BELOW OCTETT BARYON MASS
        IF(IREJ.EQ.1) THEN
          IRDS11=IRDS11 + 1
          IF(IPEV.GE.2) THEN
            WRITE(6,'(A,I5)') ' KKEVDS - IRDS11=',IRDS11
            WRITE(6,'(A,6I5/6E12.4/2E12.4)') ' DS:', IPSQ(IXSPR),ITTV1
     +      (IXSTA),ITTV2(IXSTA),IJNCH1,NNCH1,IREJ, XPSQ(IXSPR),XPSAQ
     +      (IXSPR),XPSQCM,XPSACM, XTVQ(IXSTA),XTVD(IXSTA),AMCH1,AMCH1N
          ENDIF
                                                                 GOTO 11
        ENDIF
C                                 CORRECT KINEMATICS FOR CHAIN 1
C***                MOMENTUM CORRECTION FOR CHANGED MASS OF CHAIN 1
        IF(NNCH1.NE.0)THEN
           CALL CORMOM(AMCH1,AMCH2,AMCH1N,AMCH2N,
     +         PTXSQ1,PTYSQ1,PLQ1,EQ1,
     +         PTXSA1,PTYSA1,PLAQ1,EAQ1,
     +         PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +         PTXSQ2,PTYSQ2,PLQ2,EQ2,
     +         PTXCH1,PTYCH1,PTZCH1,ECH1, PTXCH2,PTYCH2,PTZCH2,ECH2,
     +         IREJ)
        AMCH2=AMCH2N
        ENDIF
        IF(IREJ.EQ.1)GO TO 11
C
        IF (IPEV.GE.6) WRITE(6,'(A,I10/A,5F12.5/A,5F12.5)')
     +  ' DS(2): IREJ ',IREJ, '     AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',
     +  AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1,
     +  '     AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ', AMCH2,PTXCH2,PTYCH2,
     +  PTZCH2,ECH2
C
C  REPLACE SMALL MASS CHAINS BY octet or decuplet baryons
C       SECOND FOR CHAIN 2 (proj sadiquark - tar saquark)
C
        CALL COBCMA(IPSAQ(IXSPR),IPSAQ2(IXSPR),ITSAQ(IXSTA),
     +              IJNCH2,NNCH2,IREJ,AMCH2,AMCH2N,2)
c  rejection of both s-s chains if mass of chain 2 too low
        IF(IREJ.EQ.1) THEN
          IRDS12=IRDS12 + 1
          IF(IPEV.GE.2) THEN
            WRITE(6,1090) IRDS12
            WRITE(6,1100) IPSAQ(IXSPR),IPSAQ2(IXSPR),ITSAQ(IXSTA),
     +                    IJNCH2,NNCH2,IREJ,
     +      XPSQ(IXSPR),XPSAQ(IXSPR),XPSQCM,XPSACM, XTSQ(IXSTA),XTSAQ
     +      (IXSTA),XTSQCM,XTSACM, AMCH2,AMCH2N
 1090       FORMAT(' KKEVDS - IRDS12=',I5)
 1100       FORMAT(' DS - 1100', 6I5/2(4E12.4/),2E12.4)
          ENDIF
                                                                 GOTO 11
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
          NORIG=51
          CALL CORVAL(AMMM,IREJ,AMCH1,AMCH2, QTXCH1,QTYCH1,QTZCH1,QECH1,
     +    QTXCH2,QTYCH2,QTZCH2,QECH2,NORIG)
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
     +      ' DS - CALL CORVAL: AMMM, AMCH1, AMCH2, NNCH1, NNCH2, IREJ',
     +      AMMM, AMCH1, AMCH2, NNCH1, NNCH2, IREJ
            WRITE(6,1050) IREJ, AMCH1,
     +      PTXCH1,PTYCH1,PTZCH1,ECH1, AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2
 1050 FORMAT (' DS: IREJ || ',I10/
     +        '     AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',5F12.5/
     +        '     AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ',5F12.5)
          ENDIF
          IF(IREJ.EQ.1) THEN
*                           AMCH1N + AMCH2N > AMMM - 0.2
*                           reject event
            IRDS14=IRDS14+1
                                                                 GOTO 11
          ENDIF
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
       PQDSA1(N,1)=PTXSQ1
       PQDSA1(N,2)=PTYSQ1
       PQDSA1(N,3)=PLQ1
       PQDSA1(N,4)=EQ1
       PQDSA2(N,1)=PTXSQ2
       PQDSA2(N,2)=PTYSQ2
       PQDSA2(N,3)=PLQ2
       PQDSA2(N,4)=EQ2
       PQDSB1(N,1)=PTXSA2
       PQDSB1(N,2)=PTYSA2
       PQDSB1(N,3)=PLAQ2
       PQDSB1(N,4)=EAQ2
       PQDSB2(N,1)=PTXSA1
       PQDSB2(N,2)=PTYSA1
       PQDSB2(N,3)=PLAQ1
       PQDSB2(N,4)=EAQ1
C-------------------
C
C                                      PUT D-S CHAIN ENDS INTO /HKKEVT/
C                                      MOMENTA IN NN-CMS
C                                      POSITION OF ORIGINAL NUCLEONS
C
****  keep for the moment the old s-v notations
C                                 FLAG FOR DS-CHAIN ENDS
C                                            PROJECTILE: ISTHKK=131
C                                            TARGET:     ISTHKK=122
C                                      FOR DS-CHAINS     ISTHKK=4
C
        IHKKPD=JHKKPS(IXSPR )
        IHKKPO=JHKKPS(IXSPR )-1
        IHKKTD=JHKKTS(IXSTA )
        IHKKTO=JHKKTS(IXSTA )-1
        IF (IPEV.GT.3)WRITE(6,1000)IXSPR,INUCPR,JNUCPR,IHKKPO,IHKKPD
 1000 FORMAT (' IXSPR,INUCPR,JNUCPR,IHKKPO,IHKKPD ',5I5)
        IF (IPEV.GT.3)WRITE(6,1010)IXSTA,INUCTA,JNUCTA,IHKKTO,IHKKTD
 1010 FORMAT (' IXSTA,INUCTA,JNUCTA,IHKKTO,IHKKTD ',5I5)
C                                     CHAIN 1 PROJECTILE SEA-diquark
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=131
        IDHKK(IHKK)=IDHKK(IHKKPO)
        JMOHKK(1,IHKK)=IHKKPO
        JMOHKK(2,IHKK)=JMOHKK(1,IHKKPO)
        JDAHKK(1,IHKK)=IHKK+2
        JDAHKK(2,IHKK)=IHKK+2
        PHKK(1,IHKK)=PQDSA1(N,1)
        PHKK(2,IHKK)=PQDSA1(N,2)
        PHKK(3,IHKK)=PQDSA1(N,3)
        PHKK(4,IHKK)=PQDSA1(N,4)
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
C                                     CHAIN 1 TARGET QUARK
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=122
        IDHKK(IHKK)=IDHKK(IHKKTD)
        JMOHKK(1,IHKK)=IHKKTD
        JMOHKK(2,IHKK)=JMOHKK(1,IHKKTD)
        JDAHKK(1,IHKK)=IHKK+1
        JDAHKK(2,IHKK)=IHKK+1
        PHKK(1,IHKK)=PQDSA2(N,1)
        PHKK(2,IHKK)=PQDSA2(N,2)
        PHKK(3,IHKK)=PQDSA2(N,3)
        PHKK(4,IHKK)=PQDSA2(N,4)
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
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
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
        VHKK(1,NHKK)= VHKK(1,NHKK-1)
        VHKK(2,NHKK)= VHKK(2,NHKK-1)
        VHKK(3,NHKK)= VHKK(3,NHKK-1)
        VHKK(4,NHKK)=VHKK(3,NHKK)/BETP-VHKK(3,NHKK-2)/BGAMP
        MHKKDS(N)=IHKK
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
C                                   CHAIN 2 projectile sea antidiquark
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=131
        IDHKK(IHKK)=IDHKK(IHKKPD)
        JMOHKK(1,IHKK)=IHKKPD
        JMOHKK(2,IHKK)=JMOHKK(1,IHKKPD)
        JDAHKK(1,IHKK)=IHKK+2
        JDAHKK(2,IHKK)=IHKK+2
        PHKK(1,IHKK)=PQDSB1(N,1)
        PHKK(2,IHKK)=PQDSB1(N,2)
        PHKK(3,IHKK)=PQDSB1(N,3)
        PHKK(4,IHKK)=PQDSB1(N,4)
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
 
C                                     CHAIN 2 TARGET diquark
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=122
        IDHKK(IHKK)=IDHKK(IHKKTO)
        JMOHKK(1,IHKK)=IHKKTO
        JMOHKK(2,IHKK)=JMOHKK(1,IHKKTO)
        JDAHKK(1,IHKK)=IHKK+1
        JDAHKK(2,IHKK)=IHKK+1
        PHKK(1,IHKK)=PQDSB2(N,1)
        PHKK(2,IHKK)=PQDSB2(N,2)
        PHKK(3,IHKK)=PQDSB2(N,3)
        PHKK(4,IHKK)=PQDSB2(N,4)
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
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
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
        MHKKDS(N)=IHKK
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
C  NOW WE HAVE AN ACCEPTABLE SEA--VALENCE  EVENT
*     sea diquark pair!
C  AND PUT IT INTO THE HISTOGRAM
C
        AMCDS1(N)=AMCH1
        AMCDS2(N)=AMCH2
        GACDS1(N)=QECH1/AMCH1
        BGXDS1(N)=QTXCH1/AMCH1
        BGYDS1(N)=QTYCH1/AMCH1
        BGZDS1(N)=QTZCH1/AMCH1
        GACDS2(N)=QECH2/AMCH2
        BGXDS2(N)=QTXCH2/AMCH2
        BGYDS2(N)=QTYCH2/AMCH2
        BGZDS2(N)=QTZCH2/AMCH2
        NCHDS1(N)=NNCH1
        NCHDS2(N)=NNCH2
        IJCDS1(N)=IJNCH1
        IJCDS2(N)=IJNCH2
        IF (IPEV.GE.6) WRITE(6,'(A/I10,4F12.7,5I5/10X,4F12.6/10X,6F12.6,
     +4I5/8F15.5/                8F15.5)') ' DS / FINAL PRINT',N
C    +, XPSQ
C    +  (IXSPR),XPSAQ(IXSPR),XTSQ(IXSTA),XTSAQ(IXSTA), IPSQ(IXSPR),IPSAQ
C    +  (IXSPR), ITSQ(IXSTA),ITTV1(IXSTA),ITTV2(IXSTA), AMCDS1(N),AMCDS2
C    +  (N),GACDS1(N),GACDS2(N), BGXDS1(N),BGYDS1(N),BGZDS1(N), BGXDS2
C    +  (N),BGYDS2(N),BGZDS2(N), NCHDS1(N),NCHDS2(N),IJCDS1(N),IJCDS2
C    +  (N), (PQDSA1(N,JU),PQDSA2(N,JU),PQDSB1(N,JU), PQDSB2(N,JU),JU=1,
C    +  4)
                                                                GO TO 20
C***                     TREATMENT OF REJECTED SEA-SEA INTERACTIONS
   11   CONTINUE
        NCHDS1(N)=99
        NCHDS2(N)=99
        XPVD(JNUCPR)=XPVD(JNUCPR) + XPSQ(IXSPR) + XPSAQ(IXSPR)
        XTVD(JNUCTA)=XTVD(JNUCTA) + XTSAQ(IXSTA) + XTSQ(IXSTA)
      ISSQQ=ABS(IPSAQ(IXSPR))
      JSSQQ=ABS(IPSAQ2(IXSPR))
      IF(ISSQQ.EQ.3.AND.JSSQQ.EQ.3)THEN
        IDSRE(3)=IDSRE(3)+1
      ELSEIF(ISSQQ.EQ.3.OR.JSSQQ.EQ.3)THEN
        IDSRE(2)=IDSRE(2)+1
      ELSE
        IDSRE(1)=IDSRE(1)+1
      ENDIF
   20   CONTINUE
   10 CONTINUE
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     DEBUG SUBCHK
C     END DEBUG
      SUBROUTINE HADRDS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C-------------------------
C
C                       hadronize sea diquark - valence CHAINS
C
C                       ADD GENERATED HADRONS TO /ALLPAR/
C                          STARTING AT (NAUX + 1)
C                       AND TO /HKKEVT/ STARTING AT (NHKK + 1)
C
C---------------------------------------------------------
      COMMON /ZSEA/ZSEAAV,ZSEASU,ANZSEA
       COMMON/POPCCK/PDBCK,PDBSE,PDBSEU,
     *  IJPOCK,IREJCK,ICK4,IHAD4,ICK6,IHAD6
     *,IREJSE,ISE4,ISE6,IREJS3,ISE43,ISE63,IREJS0
     *,IHADA4,IHADA6,IREJSA,ISEA4,ISEA6,IREJA3,
     *ISEA43,ISEA63,IREJAO
*KEEP,INTMX.
      PARAMETER (INTMX=2488,INTMD=252)
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
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
*KEEP,ABRDS.
      COMMON /ABRDS/ AMCDS1(248),AMCDS2(248),GACDS1(248),GACDS2(248),
     +BGXDS1(248),BGYDS1(248),BGZDS1(248), BGXDS2(248),BGYDS2(248),
     +BGZDS2(248), NCHDS1(248),NCHDS2(248),IJCDS1(248),IJCDS2(248),
     +PQDSA1(248,4),PQDSA2(248,4), PQDSB1(248,4),PQDSB2(248,4)
*KEEP,LOZUO.
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/ ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),ZUOST(INTMX),
     +INTLO(INTMX),INLOSS(INTMX)
C  /LOZUO/
C       INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C------------------
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
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
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
      COMMON /DIQSUM/NDVUU,NDVUS,NDVSS,NVDUU,NVDUS,NVDSS,
     *               NDSUU,NDSUS,NDSSS,NSDUU,NSDUS,NSDSS,
     *               NDZUU,NDZUS,NDZSS,NZDUU,NZDUS,NZDSS 
     *              ,NADVUU,NADVUS,NADVSS,NAVDUU,NAVDUS,NAVDSS,
     *               NADSUU,NADSUS,NADSSS,NASDUU,NASDUS,NASDSS,
     *               NADZUU,NADZUS,NADZSS,NAZDUU,NAZDUS,NAZDSS 
C---------------------
      DIMENSION POJ(4),PAT(4)
      DATA NCALDS /0/
C-----------------------------------------------------------------------
      IF(IPHKK.GE.6)WRITE (6,'( A)') ' hadrds'
      NCALDS=NCALDS+1
      DO 50 I=1,NDS
C-----------------------drop recombined chain pairs
        IF(NCHDS1(I).EQ.99.AND.NCHDS2(I).EQ.99) GO TO 50
        IS1=INTDS1(I)
        IS2=INTDS2(I)
C
        IF (IPCO.GE.6) WRITE (6,1000) IPSQ(IS1),IPSAQ(IS1),ITSQ(IS2),
     +  ITSAQ(IS2),ITTV2(IS2), AMCDS1(I),AMCDS2(I),GACDS1(I),GACDS2(I),
     +  BGXDS1(I),BGYDS1(I),BGZDS1(I), BGXDS2(I),BGYDS2(I),BGZDS2(I),
     +  NCHDS1(I),NCHDS2(I),IJCDS1(I),IJCDS2(I), PQDSA1(I,4),PQDSA2
     +  (I,4),PQDSB1(I,4),PQDSB2(I,4)
 1000 FORMAT(10X,5I5,10F9.2/10X,4I5,4F12.4)
C
C++++++++++++++++++++++++++++++    CHAIN 1:  diquark-quark   +++++++++++
        IFB1=IPSQ(IS1)
        IFB2=IPSQ2(IS1)
        IFB3=ITSQ(IS2)
        DO 10 J=1,4
          POJ(J)=PQDSA1(I,J)
          PAT(J)=PQDSA2(I,J)
   10   CONTINUE
        IF((NCHDS1(I).NE.0.OR.NCHDS2(I).NE.0).AND.IP.NE.1)
     &  CALL SAPTRE(AMCDS1(I),GACDS1(I),BGXDS1(I),BGYDS1(I),BGZDS1(I),
     &              AMCDS2(I),GACDS2(I),BGXDS2(I),BGYDS2(I),BGZDS2(I))
C----------------------------------------------------------------
C----------------------------------------------------------------
        IF(IPCO.GE.3)WRITE (6,1244) POJ,PAT
 1244   FORMAT ('  D-S QUARK-DIQUARK POJ,PAT ',8E12.3)
*       IF(AMCDS1(I).LT.1.6)THEN
*         IF(NCHDS1(I).EQ.0)THEN
*           WRITE(6,'(A,F10.2,5I5)')' HADRDS AMCDS1(I),NCHDS1(I),I ',
*    +                AMCDS1(I),NCHDS1(I),IJCDS1(I),I,IS1,IS2
*           RETURN
*         ENDIF
*       ENDIF
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
C               I= number of valence chain
C               Projectile Nr ipp = IFROVP(INTVS1(I))
C          No of Glauber sea q at Projectile JIPP=JSSHS(IPP)
C      IPPP = IFROVP(INTVS1(I))
C      WRITE(6,*)' INTVS1(I) INTDS1(I)',INTVS1(I),INTDS1(I)
C      IF(INTDS1(I).GE.1.AND.INTDS1(I).LE.INTMD)THEN
C      IPPP = IFROVP(INTDS1(I))
C      ELSE
C      WRITE(6,*)' HADRDS: INTDS1(I) ',INTDS1(I)
       IPPP=0
C      ENDIF
C      WRITE(6,*)' IPPP ',IPPP
C      IF(IPPP.GT.0)THEN 
C        JIPP=JSSHS(IPPP)
C      ELSEIF(IPPP.EQ.0)THEN
	 JIPP=1
C      ENDIF
C      WRITE(6,*)' JIPP ',JIPP
C       IF(NCHVS2(I).EQ.0)THEN
       IF(IPCO.GE.3)WRITE(6,'(A,3I5)')'HADRDS: I,IPPP,JIPP ',
     *                     I,IPPP,JIPP
C       ENDIF
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
         IF(IFB1.LE.2.AND.IFB2.LE.2)THEN
	   NDSUU=NDSUU+1
	 ELSEIF((IFB1.EQ.3.AND.IFB2.LE.2).OR.
     *  	 (IFB2.EQ.3.AND.IFB1.LE.2))THEN
	   NDSUS=NDSUS+1
	 ELSEIF(IFB1.EQ.3.AND.IFB2.EQ.3)THEN
	   NDSSS=NDSSS+1
	 ENDIF  
       IF((NCHDS1(I).NE.0))
     *  CALL HADJET(NHAD,AMCDS1(I),POJ,PAT,GACDS1(I),
     *  BGXDS1(I), BGYDS1
     +  (I),BGZDS1(I),IFB1,IFB2,IFB3,IFB4, IJCDS1(I),
     *  IJCDS1(I),6,NCHDS1
     +  (I),15)
C---------------------------------------------------------------
        AACK=FLOAT(ICK6)/FLOAT(ICK6+IHAD6+1)
        IF((NCHDS1(I).EQ.0))THEN
          ZSEAWU=RNDM(BB)*2.D0*ZSEAAV
        RSEACK=FLOAT(JITT)*PDBSE +ZSEAWU*PDBSEU
          IF(IPCO.GE.1)WRITE(6,'(2A,I5,2F10.3)')'HADJSE JIPP,',
     *    'RSEACK,PDBSE 3 dpmnuc5',
     +    JIPP,RSEACK,PDBSE
          IREJSS=5
          IF(RNDM(V).LE.RSEACK)THEN
            IREJSS=2
            IF(AMCDS1(I).GT.2.3D0)THEN
              IREJSS=0
              CALL HADJSE(NHAD,AMCDS1(I),POJ,PAT,GACDS1(I),BGXDS1(I),
     *        BGYDS1
     +        (I),BGZDS1(I),IFB1,IFB2,IFB3,IFB4, IJCDS1(I),IJCDS1(I),6,
     *        NCHDS1
     +        (I),6,IREJSS,IISSQQ)
              IF(IPCO.GE.1)WRITE(6,'(2A,I5,F10.3,I5)')'HADJSE JIPP,',
     *        'RSEACK,IREJSS 3 dpmnuc5 ',
     +        JIPP,RSEACK,IREJSS
            ENDIF
            IF(IREJSS.GE.1)THEN
              IF(IREJSS.EQ.1)IREJSE=IREJSE+1
              IF(IREJSS.EQ.3)IREJS3=IREJS3+1
              IF(IREJSS.EQ.2)IREJS0=IREJS0+1
              CALL HADJET(NHAD,AMCDS1(I),POJ,PAT,GACDS1(I),
     *        BGXDS1(I), BGYDS1
     +        (I),BGZDS1(I),IFB1,IFB2,IFB3,IFB4, IJCDS1(I),
     *        IJCDS1(I),6,NCHDS1
     +        (I),15)
              IHAD6=IHAD6+1
            ENDIF
            IF(IREJSS.EQ.0)THEN
              IF(IISSQQ.EQ.3)THEN
                ISE63=ISE63+1
              ELSE
                ISE6=ISE6+1
              ENDIF
            ENDIF
          ELSE
            CALL HADJET(NHAD,AMCDS1(I),POJ,PAT,GACDS1(I),
     *      BGXDS1(I), BGYDS1
     +      (I),BGZDS1(I),IFB1,IFB2,IFB3,IFB4, IJCDS1(I),
     *      IJCDS1(I),6,NCHDS1
     +      (I),15)
            IHAD6=IHAD6+1
          ENDIF
        ENDIF

C---------------------------------------------------------------
        ACOUDS=ACOUDS+1
        NHKKAU=NHKK+1
        DO 20 J=1,NHAD
          IF (NHKK.EQ.NMXHKK) THEN
            WRITE (6,'(A,2I5/A)') ' HADRDS: NHKK.EQ.NMXHKK ',NHKK,NMXHKK
            RETURN
          ENDIF
C         NHKK=NHKK+1
          IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
C 
          EHECC=SQRT(PXF(J)**2+PYF(J)**2+PZF(J)**2+AMF(J)**2)
        IF (ABS(EHECC-HEF(J)).GT.0.001) THEN
C           WRITE(6,'(2A/3I5,3E15.6)')
C    &            ' HADRDS / CHAIN 1 : CORRECT INCONSISTENT ENERGY ',
C    *            '  NCALDS, NHKK,NREF(J), HEF(J),EHECC, AMF(J)',
C    *            NCALDS, NHKK,NREF(J), HEF(J),EHECC, AMF(J)
          HEF(J)=EHECC
          ENDIF
          ANNDS=ANNDS+1
          EEDS=EEDS+HEF(J)
          PTDS=PTDS+SQRT(PXF(J)**2+PYF(J)**2)
C                               PUT NN-CMS HADRONS INTO /HKKEVT/
          ISTIST=1
          IF(IBARF(J).EQ.500)ISTIST=2
	  IF(IPCO.GE.3)WRITE(6,*)' HADRDS before HKKFIL'
          CALL HKKFIL(ISTIST,MPDGHA(NREF(J)),MHKKDS(I)-3,0,
     *                PXF(J),PYF(J),PZF(J),HEF(J),NHKKAU,IORMO(J),17)
          IF(IDHKK(NHKK).EQ.99999) WRITE (6,1030)NHKK,NREF(J), IDHKK
     +    (NHKK)
          IF (IPHKK.GE.2) WRITE(6,1010) NHKK, ISTHKK(NHKK),IDHKK(NHKK),
     +    JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +    (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
   20   CONTINUE
C       IF(NHAD.GT.0) THEN
C         JDAHKK(1,IMOHKK)=NHKKAU
C         JDAHKK(2,IMOHKK)=NHKK
C       ENDIF
C+++++++++++++++++++++++++++++   CHAIN 2: adiquark - aquark  +++++++++
        IFB1=IPSAQ(IS1)
        IFB2=IPSAQ2(IS1)
        IFB3=ITSAQ(IS2)
        IFB1=IABS(IFB1)+6
        IFB2=IABS(IFB2)+6
        IFB3=IABS(IFB3)+6
        DO 30 J=1,4
          POJ(J)=PQDSB2(I,J)
          PAT(J)=PQDSB1(I,J)
   30   CONTINUE
*       IF(AMCDS2(I).LT.1.6)THEN
*         IF(NCHDS2(I).EQ.0)THEN
C           WRITE(6,'(A,F10.2,5I5)')' HADRDS AMCDS2(I),NCHDS2(I),I ',
C    +                AMCDS2(I),NCHDS2(I),IJCDS2(I),I,IS1,IS2
*           RETURN
*         ENDIF
*       ENDIF
C
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
C               I= number of valence chain
C               Projectile Nr ipp = IFROVP(INTVS1(I))
C          No of Glauber sea q at Projectile JIPP=JSSHS(IPP)
C      WRITE(6,*)'HADRDS INTDS2(I) ',INTDS2(I)
C      IPPP = IFROVP(INTVS1(I))
C      IF(INTDS2(I).GE.1)THEN
C      ITTT = IFROVP(INTDS2(I))
C      ELSE
       ITTT=0
C      ENDIF
C      WRITE(6,*)' HADRDS 2 IPPP', IPPP
C      IF(ITTT.GT.0)THEN
C      JITT=JSSHS(ITTT)
C      ELSE
       JITT=0
C      ENDIF
C       IF(NCHVS2(I).EQ.0)THEN
C      WRITE(6,'(A,3I5)')'HADRDS: I,IPPP,JIPP ',
C    *                     I,IPPP,JIPP
C       ENDIF
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
         IF(IFB1.LE.8.AND.IFB2.LE.8)THEN
	   NADSUU=NADSUU+1
	 ELSEIF((IFB1.EQ.9.AND.IFB2.LE.8).OR.
     *  	 (IFB2.EQ.9.AND.IFB1.LE.8))THEN
	   NADSUS=NADSUS+1
	 ELSEIF(IFB1.EQ.9.AND.IFB2.EQ.9)THEN
	   NADSSS=NADSSS+1
	 ENDIF  
C	 WRITE(6,*)'NCHDS2(I)',NCHDS2(I)
C	WRITE(6,*)' before HADJET:AMCDS2(I),GACDS2(I),BGXDS2(I),',
C    *  AMCDS2(I),GACDS2(I),BGXDS2(I)
        IF((NCHDS2(I).NE.0))
     *  CALL HADJET(NHAD,AMCDS2(I),POJ,PAT,GACDS2(I),
     *  BGXDS2(I), BGYDS2
     +  (I),BGZDS2(I),IFB1,IFB2,IFB3,IFB4, IJCDS2(I),
     *  IJCDS2(I),6,NCHDS2
     +  (I),16)
C      WRITE(6,*)' after HADJET '
C-----------------------------------------------------------------
        AACK=FLOAT(ICK6)/FLOAT(ICK6+IHAD6+1)
        IF((NCHDS2(I).EQ.0))THEN
          ZSEAWU=RNDM(BB)*2.D0*ZSEAAV
        RSEACK=FLOAT(JITT)*PDBSE +ZSEAWU*PDBSEU
          IF(IPCO.GE.1)WRITE(6,'(2A,I5,2F10.3)')'HADJSE JIPP,',
     *    'RSEACK,PDBSE ',
     +    JIPP,RSEACK,PDBSE
          IREJSS=5
          IF(RNDM(V).LE.RSEACK)THEN
            IREJSS=2
            IF(AMCDS2(I).GT.2.3D0)THEN
              IREJSS=0
C	      WRITE(6,*)' before HADJASE '
              CALL HADJASE(NHAD,AMCDS2(I),POJ,PAT,GACDS2(I),BGXDS2(I),
     *        BGYDS2
     +        (I),BGZDS2(I),IFB1,IFB2,IFB3,IFB4, IJCDS2(I),IJCDS2(I),6,
     *        NCHDS2
     +        (I),6,IREJSS,IISSQQ)
              IF(IPCO.GE.1)WRITE(6,'(2A,I5,F10.3,I5)')'HADJSE JIPP,',
     *        'RSEACK,IREJSS ',
     +        JIPP,RSEACK,IREJSS
            ENDIF
            IF(IREJSS.GE.1)THEN
              IF(IREJSS.EQ.1)IREJSA=IREJSA+1
              IF(IREJSS.EQ.3)IREJA3=IREJA3+1
              IF(IREJSS.EQ.2)IREJA0=IREJA0+1
C	      WRITE(6,*)' before HADJET2 '
        CALL HADJET(NHAD,AMCDS2(I),POJ,PAT,GACDS2(I),
     *  BGXDS2(I), BGYDS2
     +  (I),BGZDS2(I),IFB1,IFB2,IFB3,IFB4, IJCDS2(I),
     *  IJCDS2(I),6,NCHDS2
     +  (I),16)
              IHADA6=IHADA6+1
            ENDIF
            IF(IREJSS.EQ.0)THEN
              IF(IISSQQ.EQ.3)THEN
                ISEA63=ISEA63+1
              ELSE
                ISEA6=ISEA6+1
              ENDIF
            ENDIF
          ELSE
C	      WRITE(6,*)' before HADJET3 '
        CALL HADJET(NHAD,AMCDS2(I),POJ,PAT,GACDS2(I),
     *  BGXDS2(I), BGYDS2
     +  (I),BGZDS2(I),IFB1,IFB2,IFB3,IFB4, IJCDS2(I),
     *  IJCDS2(I),6,NCHDS2
     +  (I),16)
            IHADA6=IHADA6+1
          ENDIF
        ENDIF
C--------------------------------------------------------------------
C                                   ADD HADRONS/RESONANCES INTO
C                                   COMMON /ALLPAR/ STARTING AT NAUX
        NHKKAU=NHKK+1
        DO 40 J=1,NHAD
          IF (NHKK.EQ.NMXHKK) THEN
            WRITE (6,'(A,2I5/A)') ' HADRDS: NHKK.EQ.NMXHKK ', NHKK,
     +      NMXHKK
            RETURN
          ENDIF
C         NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
C
          EHECC=SQRT(PXF(J)**2+PYF(J)**2+PZF(J)**2+AMF(J)**2)
        IF (ABS(EHECC-HEF(J)).GT.0.001) THEN
C             WRITE(6,'(2A/3I5,3E15.6)')
C    &            ' HADRDS / CHAIN 2 : CORRECT INCONSISTENT ENERGY ',
C    *            '  NCALDS, NHKK,NREF(J), HEF(J),EHECC, AMF(J)',
C    *            NCALDS, NHKK,NREF(J), HEF(J),EHECC, AMF(J)
          HEF(J)=EHECC
          ENDIF
          ANNDS=ANNDS+1
          EEDS=EEDS+HEF(J)
          PTDS=PTDS+SQRT(PXF(J)**2+PYF(J)**2)
C                               PUT NN-CMS HADRONS INTO /HKKEVT/
          ISTIST=1
          IF(IBARF(J).EQ.500)ISTIST=2
	  IF(IPCO.GE.3)WRITE(6,*)' HADRDS before  2 HKKFIL'
          CALL HKKFIL(ISTIST,MPDGHA(NREF(J)),MHKKDS(I),0,
     *                PXF(J),PYF(J),PZF(J),HEF(J),NHKKAU,IORMO(J),18)
          IF(IDHKK(NHKK).EQ.99999) WRITE (6,1030)NHKK,NREF(J), IDHKK
     +    (NHKK)
          IF (IPHKK.GE.2) WRITE(6,1010) NHKK, ISTHKK(NHKK),IDHKK(NHKK),
     +    JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +    (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
   40   CONTINUE
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
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE DIQSSD(ECM,ITS,IPS,IREJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
*  define s-d chains (sea - sea diquark chains)
*  sq-sqsq and saq-saqsaq chains instead of q-aq and aq-q chains
      COMMON /ZSEA/ZSEAAV,ZSEASU,ANZSEA
       COMMON/POPCCK/PDBCK,PDBSE,PDBSEU,
     *  IJPOCK,IREJCK,ICK4,IHAD4,ICK6,IHAD6
     *,IREJSE,ISE4,ISE6,IREJS3,ISE43,ISE63,IREJS0
     *,IHADA4,IHADA6,IREJSA,ISEA4,ISEA6,IREJA3,
     *ISEA43,ISEA63,IREJAO
*KEEP,INTMX.
      PARAMETER (INTMX=2488,INTMD=252)
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
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
*KEEP,ABRSD.
      COMMON /ABRSD/ AMCSD1(248),AMCSD2(248),GACSD1(248),GACSD2(248),
     +BGXSD1(248),BGYSD1(248),BGZSD1(248), BGXSD2(248),BGYSD2(248),
     +BGZSD2(248), NCHSD1(248),NCHSD2(248),IJCSD1(248),IJCSD2(248),
     +PQSDA1(248,4),PQSDA2(248,4), PQSDB1(248,4),PQSDB2(248,4)
*KEND.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
       COMMON/SEASU3/SEASQ
      COMMON /XSEADI/ XSEACU,UNON,UNOM,UNOSEA, CVQ,CDQ,CSEA,SSMIMA,
     +SSMIMQ,VVMTHR
      COMMON /DIQREJ/IDIQRE(7),IDVRE(3),IVDRE(3),IDSRE(3),ISDRE(3),
     *IDZRE(3),IZDRE(3),IDIQRZ(7) 
      COMMON /DIQSUM/NDVUU,NDVUS,NDVSS,NVDUU,NVDUS,NVDSS,
     *               NDSUU,NDSUS,NDSSS,NSDUU,NSDUS,NSDSS,
     *               NDZUU,NDZUS,NDZSS,NZDUU,NZDUS,NZDSS 
     *              ,NADVUU,NADVUS,NADVSS,NAVDUU,NAVDUS,NAVDSS,
     *               NADSUU,NADSUS,NADSSS,NASDUU,NASDUS,NASDSS,
     *               NADZUU,NADZUS,NADZSS,NAZDUU,NAZDUS,NAZDSS 
C----------------------------------------------------------------------
C     COMMON /PCHARM/PCCCC
      PARAMETER (UMMM=0.3D0)
      PARAMETER (SMMM=0.5D0)
      PARAMETER (CMMM=1.3D0)
      DATA PC/0.0001D0/
*KEND.
C----------
C
      DATA INICHA/0/
C----------------------------------------------------------------------
C                     Initialize Charm selection at soft chain ends
C
      IF(INICHA.EQ.0)THEN
        RX=8.D0
        X1=RX
        GM=2.140D0
        X2=UMMM
	BETOO=7.5D0
      ENDIF
      RX=8.D0
      X1=RX
      BETCHA=BETOO+1.3D0-LOG10(ECM)
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
      PC=PC1/7.0D0
      PU1=PU/(2*PU+PS+PC)
      PS1=PS/(2*PU+PS+PC)
      IF(INICHA.EQ.0)THEN
        INICHA=1
        WRITE(6,4567)PC,BETCHA,PU1,PS1,SEASQ
 4567   FORMAT(' Charm chain ends DIQSSD: PC,BETCHA,PU,PS,SEASQ',4F10.5)
      ENDIF
C----------------------------------------------------------------------
        RR=RNDM(V)
        IS=1.D0+RNDM(V1)*(2.D0+2.D0*SEASQ)
	IF(RR.LT.PC)IS=4
C----------------------------------------------------------------------
      IF(IPHKK.GE.6)WRITE (6,'( A)') ' diqssd'
      IREJ=0
*  kinematics: is the mass of both chains big enough
*              to allow for fragmentation
      ITSQ2(ITS)=1.D0+RNDM(v1)*(2.D0+2.D0*SEASQ)
        RR=RNDM(V)
	IF(RR.LT.PC)ITSQ2(ITS)=4
C----------------------------------------------------------------------
      ITSAQ2(ITS)=-ITSQ2(ITS)
C---------------------------------------------------j.r.29.4.94
C                                   x**1.5 distr for sea diquarks
C                         number of target nucleon
      INUCTA=IFROST(ITS)
C                          number of target diquark
      IITOT=ITOVT(INUCTA)
C                          diquark x
      XTDIQU=XTVD(IITOT)
C                          minimal value of diquark x
      XDTHR=CDQ/ECM
C
      XDFREE=XTDIQU-XDTHR
      XALL=XDFREE+XTSQ(ITS)+XTSAQ(ITS)-2.*XDTHR
      XDALT=XTVD(IITOT)
      XSALT=XTSQ(ITS)
      XAALT=XTSAQ(ITS)
      IF(XALL.GE.0.)THEN
        RR1=RNDM(V1)
        RR2=RNDM(V2)
        RR3=RNDM(V3)
        SR123=RR1+RR2+RR3
        DX1=RR1*XALL/SR123
        DX2=RR2*XALL/SR123
        DX3=RR3*XALL/SR123
        XTVD(IITOT)=XDTHR+DX1
        XTSQ(ITS)=XDTHR+DX2
        XTSAQ(ITS)=XDTHR+DX3
      ENDIF
C--------------------------------------------------------------
      AMSDQ1=XTSQ(ITS)*XPSQ(IPS)*ECM**2
      AMSDQ2=XTSAQ(ITS)*XPSAQ(IPS)*ECM**2
      IDIQRE(1)=IDIQRE(1)+1
      IF(ITSQ(ITS).GE.3.AND.ITSQ2(ITS).GE.3)THEN
      IDIQRE(2)=IDIQRE(2)+1
C       IF(AMSDQ2.LE.2.30.OR.AMSDQ1.LE.2.30) THEN
        IF(AMSDQ2.LE.6.60D0.OR.AMSDQ1.LE.6.60D0) THEN
          IREJ=1
           IDIQRE(3)=IDIQRE(3)+1
           IDIQRE(2)=IDIQRE(2)-1
           IDIQRE(1)=IDIQRE(1)-1
          XTVD(IITOT)=XDALT
          XTSQ(ITS)=XSALT
          XTSAQ(ITS)=XAALT
          RETURN
        ENDIF
      ELSEIF(ITSQ(ITS).GE.3.OR.ITSQ2(ITS).GE.3)THEN
      IDIQRE(4)=IDIQRE(4)+1
C        IF(AMSDQ2.LE.1.9.OR.AMSDQ1.LE.1.90) THEN
        IF(AMSDQ2.LE.5.8D0.OR.AMSDQ1.LE.5.80D0) THEN
          IREJ=1
           IDIQRE(5)=IDIQRE(5)+1
           IDIQRE(4)=IDIQRE(4)-1
           IDIQRE(1)=IDIQRE(1)-1
          XTVD(IITOT)=XDALT
          XTSQ(ITS)=XSALT
          XTSAQ(ITS)=XAALT
          RETURN
        ENDIF
      ELSE
      IDIQRE(6)=IDIQRE(6)+1
C        IF(AMSDQ2.LE.1.50.OR.AMSDQ1.LE.1.50) THEN
        IF(AMSDQ2.LE.3.9D0.OR.AMSDQ1.LE.3.9D0) THEN
          IREJ=1
           IDIQRE(7)=IDIQRE(7)+1
           IDIQRE(6)=IDIQRE(6)-1
           IDIQRE(1)=IDIQRE(1)-1
          XTVD(IITOT)=XDALT
          XTSQ(ITS)=XSALT
          XTSAQ(ITS)=XAALT
          RETURN
        ENDIF
      ENDIF
      NSD=NSD+1
c     WRITE(6,'(A/5X,3F10.3,3I5/5X,3F10.3)')
c    +' DIQVS: AMSDQ1, XTSQ, XPVQ, IPS,ITS, NSD/ AMSDQ2, XTSAQ, XPVD',
c    +AMSDQ1,XTSQ(ITS),XPVQ(IPS),IPS,ITS,NSD,AMSDQ2,XTSAQ(ITS),XPVD(IPS)
      NCHSD1(NSD)=0
      NCHSD2(NSD)=0
      INTSD1(NSD)=IPS
      INTSD2(NSD)=ITS
      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C      DEBUG SUBCHK
C      END DEBUG
      SUBROUTINE KKEVSD(IREJSD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C------------------ treatment of sea - sea diquark CHAIN SYSTEMS
C
      COMMON /ZSEA/ZSEAAV,ZSEASU,ANZSEA
       COMMON/POPCCK/PDBCK,PDBSE,PDBSEU,
     *  IJPOCK,IREJCK,ICK4,IHAD4,ICK6,IHAD6
     *,IREJSE,ISE4,ISE6,IREJS3,ISE43,ISE63,IREJS0
     *,IHADA4,IHADA6,IREJSA,ISEA4,ISEA6,IREJA3,
     *ISEA43,ISEA63,IREJAO
*KEEP,INTMX.
      PARAMETER (INTMX=2488,INTMD=252)
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
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
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
*KEEP,ABRSD.
      COMMON /ABRSD/ AMCSD1(248),AMCSD2(248),GACSD1(248),GACSD2(248),
     +BGXSD1(248),BGYSD1(248),BGZSD1(248), BGXSD2(248),BGYSD2(248),
     +BGZSD2(248), NCHSD1(248),NCHSD2(248),IJCSD1(248),IJCSD2(248),
     +PQSDA1(248,4),PQSDA2(248,4), PQSDB1(248,4),PQSDB2(248,4)
*KEEP,DXQX.
C     INCLUDE (XQXQ)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON /DXQX/ XPVQ(248),XPVD(248),XTVQ(248),XTVD(248), XPSQ
     +(INTMX),XPSAQ(INTMX),XTSQ(INTMX),XTSAQ(INTMX)
     *                ,XPSU(248),XTSU(248)
     *                ,XPSUT(248),XTSUT(248)
      COMMON /DIQREJ/IDIQRE(7),IDVRE(3),IVDRE(3),IDSRE(3),ISDRE(3),
     *IDZRE(3),IZDRE(3),IDIQRZ(7) 
*KEEP,LOZUO.
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/ ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),ZUOST(INTMX),
     +INTLO(INTMX),INLOSS(INTMX)
C  /LOZUO/
C       INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C------------------
*KEEP,TRAFOP.
      COMMON /TRAFOP/ GAMP,BGAMP,BETP
*KEEP,NUCIMP.
      COMMON /NUCIMP/ PRMOM(5,248),TAMOM(5,248), PRMFEP,PRMFEN,TAMFEP,
     +TAMFEN, PREFEP,PREFEN,TAEFEP,TAEFEN, PREPOT(210),TAEPOT(210),
     +PREBIN,TAEBIN,FERMOD,ETACOU
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
      COMMON /DIQSUM/NDVUU,NDVUS,NDVSS,NVDUU,NVDUS,NVDSS,
     *               NDSUU,NDSUS,NDSSS,NSDUU,NSDUS,NSDSS,
     *               NDZUU,NDZUS,NDZSS,NZDUU,NZDUS,NZDSS 
     *              ,NADVUU,NADVUS,NADVSS,NAVDUU,NAVDUS,NAVDSS,
     *               NADSUU,NADSUS,NADSSS,NASDUU,NASDUS,NASDSS,
     *               NADZUU,NADZUS,NADZSS,NAZDUU,NAZDUS,NAZDSS 
*KEND.
C-------------------
      IF(IPHKK.GE.6)WRITE (6,'( A)') ' kkevsd'
      IREJSD=0
      DO 10 N=1,NSD
C---------------------------drop recombined chain pairs
        IF(NCHSD1(N).EQ.99.AND.NCHSD2(N).EQ.99)GO TO 10
C
C***            4-MOMENTA OF projectile QUARK-DIQUARK PAIRS IN NN-CMS
        IXSPR=INTSD1(N)
        INUCPR=IFROSP(IXSPR)
        JNUCPR=ITOVP(INUCPR)
C
        PSQPX=XPSQ(IXSPR)*PRMOM(1,INUCPR)
        PSQPY=XPSQ(IXSPR)*PRMOM(2,INUCPR)
        PSQPZ=XPSQ(IXSPR)*PRMOM(3,INUCPR)
        PSQE=XPSQ(IXSPR)*PRMOM(4,INUCPR)
        PSDQPX=XPSAQ(IXSPR)*PRMOM(1,INUCPR)
        PSDQPY=XPSAQ(IXSPR)*PRMOM(2,INUCPR)
        PSDQPZ=XPSAQ(IXSPR)*PRMOM(3,INUCPR)
        PSDQE=XPSAQ(IXSPR)*PRMOM(4,INUCPR)
C
C***                 4-MOMENTA OF TARGET QUARK-DIQUARK PAIRS IN NN-CMS
        IXSTA=INTSD2(N)
        INUCTA=IFROST(IXSTA)
        JNUCTA=ITOVT(INUCTA)
*
        TSQPX=XTSQ(IXSTA)*TAMOM(1,INUCTA)
        TSQPY=XTSQ(IXSTA)*TAMOM(2,INUCTA)
        TSQPZ=XTSQ(IXSTA)*TAMOM(3,INUCTA)
        TSQE=XTSQ(IXSTA)*TAMOM(4,INUCTA)
        TSAQPX=XTSAQ(IXSTA)*TAMOM(1,INUCTA)
        TSAQPY=XTSAQ(IXSTA)*TAMOM(2,INUCTA)
        TSAQPZ=XTSAQ(IXSTA)*TAMOM(3,INUCTA)
        TSAQE=XTSAQ(IXSTA)*TAMOM(4,INUCTA)
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
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
     *            PSQNX,PSQNY,PSQNZ,PSQNE,63)
      PSQPX=PSQNX
      PSQPY=PSQNY
      PSQPZ=PSQNZ
      PSQE=PSQNE      
      CALL CROMSC(PSDQPX,PSDQPY,PSDQPZ,PSDQE,RTIX,RTIY,RTIZ,
     *            PSDQNX,PSDQNY,PSDQNZ,PSDQNE,64)
      PSDQPX=PSDQNX
      PSDQPY=PSDQNY
      PSDQPZ=PSDQNZ
      PSDQE=PSDQNE      
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
     *            TSQNX,TSQNY,TSQNZ,TSQNE,65)
      TSQPX=TSQNX
      TSQPY=TSQNY
      TSQPZ=TSQNZ
      TSQE=TSQNE      
      CALL CROMSC(TSAQPX,TSAQPY,TSAQPZ,TSAQE,RTIX,RTIY,RTIZ,
     *            TSAQNX,TSAQNY,TSAQNZ,TSAQNE,66)
      TSAQPX=TSAQNX
      TSAQPY=TSAQNY
      TSAQPZ=TSAQNZ
      TSAQE=TSAQNE  
      ENDIF    
C                                                ---------
C                                                j.r.10.5.93
       IF(IP.GE.0)GO TO 1779
        PSQPZ2=PSQE**2-PSQPX**2-PSQPY**2
        IF(PSQPZ2.GE.0.)THEN
          PSQPZ=SQRT(PSQPZ2)
        ELSE
          PSQPX=0.
          PSQPY=0.
          PSQPZ=PSQE
        ENDIF
C
        PDQPZ2=PSDQE**2-PSDQPX**2-PSDQPY**2
        IF(PDQPZ2.GE.0.)THEN
          PSDQPZ=SQRT(PDQPZ2)
        ELSE
          PSDQPX=0.
          PSDQPY=0.
          PSDQPZ=PSDQE
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
C                                                ---------
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
        PTXSA1=PSDQPX
        PTXSQ2=TSQPX
        PTXSA2=TSAQPX
        PTYSQ1=PSQPY
        PTYSA1=PSDQPY
        PTYSQ2=TSQPY
        PTYSA2=TSAQPY
        PLQ1=PSQPZ
        PLAQ1=PSDQPZ
        PLQ2=TSQPZ
        PLAQ2=TSAQPZ
        EQ1=PSQE
        EAQ1=PSDQE
        EQ2=TSQE
        EAQ2=TSAQE
C                                       ---------------
C
C                                               _________________
          IF(IPEV.GE.2) THEN
            WRITE(6,'(A,I5)') ' KKEVSD - IRSD13=',IRSD13
            WRITE(6,'(A/4(4E12.4/),2E12.4/2I5/4E12.4)')
     +      ' VD:  ...',  PTXSQ1,PTYSQ1,PLQ1,EQ1,PTXSA1,PTYSA1,
     +      PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +      AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1
          ENDIF
        IKVALA=0
        NSELPT=1
        CALL SELPT( PTXSQ1,PTYSQ1,PLQ1,
     +             EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1,
     +             PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +             PTXSQ2,PTYSQ2,PLQ2,EQ2,
     +             AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1,
     * PTTQ2,PTTA2,
     * NSELPT)
          IF(IPEV.GE.2) THEN
            WRITE(6,'(A,I5)') ' KKEVSD - IRSD13=',IRSD13
            WRITE(6,'(A/4(4E12.4/),2E12.4/2I5/4E12.4)')
     +      ' VD:  ...',  PTXSQ1,PTYSQ1,PLQ1,EQ1,PTXSA1,PTYSA1,
     +      PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +      AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1
          ENDIF
        IF (IPEV.GE.7) WRITE(6,'(A/5X,I10)')
     +  'SD   IREJ ', IREJ
        IF (IREJ.EQ.1) THEN
          IRSD13=IRSD13 + 1
          IF(IPEV.GE.2) THEN
            WRITE(6,'(A,I5)') ' KKEVSD - IRSD13=',IRSD13
            WRITE(6,'(A/4(4E12.4/),2E12.4/2I5/4E12.4)')
     +      ' VD:  ...',  PTXSQ1,PTYSQ1,PLQ1,EQ1,PTXSA1,PTYSA1,
     +      PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +      AMCH1,AMCH2,IREJ,IKVALA,PTTQ1,PTTA1
 
          ENDIF
                                                                GO TO 11
        ENDIF
C
C***  4-MOMENTA OF CHAINS IN THIS FRAME
C
        PTXCH1=PTXSQ1 + PTXSQ2
        PTYCH1=PTYSQ1 + PTYSQ2
        PTZCH1=PLQ1 + PLQ2
        ECH1=EQ1 + EQ2
        PTXCH2=PTXSA2 + PTXSA1
        PTYCH2=PTYSA2 + PTYSA1
        PTZCH2=PLAQ2 + PLAQ1
        ECH2=EAQ2 + EAQ1
        AMMM=SQRT((ECH1+ECH2)**2-(PTXCH1+PTXCH2)**2
     +            -(PTYCH1+PTYCH2)**2-(PTZCH1+PTZCH2)**2) 
C
C
        IF (IPEV.GE.6) WRITE(6,'(A,I10/A,5F12.5/A,5F12.5)')
     +  ' SD: IREJ ',IREJ, '     AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',
     +  AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1,
     +  '     AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ', AMCH2,PTXCH2,PTYCH2,
     +  PTZCH2,ECH2
 
C
C  REPLACE SMALL MASS CHAINS BY PSEUDOSCALAR OR VECTOR MESONS OR OCTETT
C                                              OR DECUPLETT BARYONS
C  FIRST FOR CHAIN 1  (PROJ quark - tar sea-diquark)
C
        CALL COBCMA(ITSQ(IXSTA),ITSQ2(IXSTA),IPSQ(IXSPR), IJNCH1,NNCH1,
     +  IREJ,AMCH1,AMCH1N,1)
C***                            MASS BELOW OCTETT BARYON MASS
        IF(IREJ.EQ.1) THEN
          IRSD11=IRSD11 + 1
          IF(IPEV.GE.2) THEN
            WRITE(6,'(A,I5)') ' KKEVSD - IRSD11=',IRSD11
            WRITE(6,'(A,6I5/6E12.4/2E12.4)') ' SD:', IPVQ(IXSPR),ITSQ
     +      (IXSTA),ITSQ2(IXSTA),IJNCH1,NNCH1,IREJ, XPVQ(IXSPR),XPVD
     +      (IXSPR),XPSQCM,XPSDCM, XTSQ(IXSTA),XTSAQ(IXSTA),AMCH1,AMCH1N
          ENDIF
                                                                 GOTO 11
        ENDIF
C                                 CORRECT KINEMATICS FOR CHAIN 1
C***                MOMENTUM CORRECTION FOR CHANGED MASS OF CHAIN 1
        IF(NNCH1.NE.0)THEN
           CALL CORMOM(AMCH1,AMCH2,AMCH1N,AMCH2N,
     +         PTXSQ1,PTYSQ1,PLQ1,EQ1,
     +         PTXSA1,PTYSA1,PLAQ1,EAQ1,
     +         PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +         PTXSQ2,PTYSQ2,PLQ2,EQ2,
     +         PTXCH1,PTYCH1,PTZCH1,ECH1, PTXCH2,PTYCH2,PTZCH2,ECH2,
     +         IREJ)
        AMCH2=AMCH2N
        ENDIF
        IF (IREJ.EQ.1)GO TO 11
C
        IF (IPEV.GE.6) WRITE(6,'(A,I10/A,5F12.5/A,5F12.5)')
     +  ' SD(2): IREJ ',IREJ, '     AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',
     +  AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1,
     +  '     AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ', AMCH2,PTXCH2,PTYCH2,
     +  PTZCH2,ECH2
C
C  REPLACE SMALL MASS CHAINS BY octet or decuplet baryons
C       SECOND FOR CHAIN 2 (proj saquark - tar sadiquark)
C
        CALL COBCMA(IPSAQ(IXSPR),ITSAQ(IXSTA),ITSAQ2(IXSTA),
     +              IJNCH2,NNCH2,IREJ,AMCH2,AMCH2N,2)
c  rejection of both s-s chains if mass of chain 2 too low
        IF(IREJ.EQ.1) THEN
          IRSD12=IRSD12 + 1
          IF(IPEV.GE.2) THEN
            WRITE(6,1090) IRSD12
            WRITE(6,1100) IPSAQ(IXSPR),ITSAQ(IXSTA),ITSAQ2(IXSTA),
     +                    IJNCH2,NNCH2,IREJ,
     +      XPSQ(IXSPR),XPSAQ(IXSPR),XPSQCM,XTSACM, XTSQ(IXSTA),XTSAQ
     +      (IXSTA),XTSQCM,XTSACM, AMCH2,AMCH2N
 1090       FORMAT(' KKEVSD - IRSD12=',I5)
 1100       FORMAT(' SD - 1100', 6I5/2(4E12.4/),2E12.4)
          ENDIF
                                                                 GOTO 11
        ENDIF
C                                if AMCH2 changed in COBCMA/COMCMA
C                                CORVAL corrects chain kinematics
C                                according to 2-body kinem.
C                                with fixed masses
        IF(NNCH2.NE.0) THEN
          AMCH2=AMCH2N
C                        TRANSFORM BOTH CHAINS  INTO TWO CHAIN-CMS
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
C-------------------
C                                   4-MOMENTA OF CHAINS
        CALL DALTRA(GAMMM,-BGGGX,-BGGGY,-BGGGZ, 
     +  PTXCH1,PTYCH1,PTZCH1,ECH1,
     +  PPPCH1, QTXCH1,QTYCH1,QTZCH1,QECH1)
C
        CALL DALTRA(GAMMM,-BGGGX,-BGGGY,-BGGGZ,
     +  PTXCH2,PTYCH2,PTZCH2,ECH2,
     +  PPPCH2, QTXCH2,QTYCH2,QTZCH2,QECH2)
C
          NORIG=52
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
     +      ' SD - CALL CORVAL: AMMM, AMCH1, AMCH2, NNCH1, NNCH2, IREJ',
     +      AMMM, AMCH1, AMCH2, NNCH1, NNCH2, IREJ
            WRITE(6,1050) IREJ, AMCH1,
     +      PTXCH1,PTYCH1,PTZCH1,ECH1, AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2
 1050 FORMAT (' SD: IREJ || ',I10/
     +        '     AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',5F12.5/
     +        '     AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ',5F12.5)
          ENDIF
          IF(IREJ.EQ.1) THEN
*                           AMCH1N + AMCH2N > AMMM - 0.2
*                           reject event
            IRSD14=IRSD14+1
                                                                 GOTO 11
          ENDIF
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
       PQSDA1(N,1)=PTXSQ1
       PQSDA1(N,2)=PTYSQ1
       PQSDA1(N,3)=PLQ1
       PQSDA1(N,4)=EQ1
       PQSDA2(N,1)=PTXSA2
       PQSDA2(N,2)=PTYSA2
       PQSDA2(N,3)=PLAQ2
       PQSDA2(N,4)=EAQ2
       PQSDB1(N,1)=PTXSQ2
       PQSDB1(N,2)=PTYSQ2
       PQSDB1(N,3)=PLQ2
       PQSDB1(N,4)=EQ2
       PQSDB2(N,1)=PTXSA1
       PQSDB2(N,2)=PTYSA1
       PQSDB2(N,3)=PLAQ1
       PQSDB2(N,4)=EAQ1
C-------------------
C
C                                      PUT D-S CHAIN ENDS INTO /HKKEVT/
C                                      MOMENTA IN NN-CMS
C                                      POSITION OF ORIGINAL NUCLEONS
C
****  keep for the moment the old v-s notations
C                                 FLAG FOR SD-CHAIN ENDS
C                                            PROJECTILE: ISTHKK=121
C                                            TARGET:     ISTHKK=132
C                                      FOR SD-CHAINS     ISTHKK=5
C
        IHKKPD=JHKKPS(IXSPR )
        IHKKPO=JHKKPS(IXSPR )-1
        IHKKTD=JHKKTS(IXSTA )
        IHKKTO=JHKKTS(IXSTA )-1
        IF (IPEV.GT.3)WRITE(6,1000)IXSPR,INUCPR,JNUCPR,IHKKPO,IHKKPD
 1000 FORMAT (' IXSPR,INUCPR,JNUCPR,IHKKPO,IHKKPD ',5I5)
        IF (IPEV.GT.3)WRITE(6,1010)IXSTA,INUCTA,JNUCTA,IHKKTO,IHKKTD
 1010 FORMAT (' IXSTA,INUCTA,JNUCTA,IHKKTO,IHKKTD ',5I5)
C                                     CHAIN 1 PROJECTILE SEA-diquark
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=121
        IDHKK(IHKK)=IDHKK(IHKKPO)
        JMOHKK(1,IHKK)=IHKKPO
        JMOHKK(2,IHKK)=JMOHKK(1,IHKKPO)
        JDAHKK(1,IHKK)=IHKK+2
        JDAHKK(2,IHKK)=IHKK+2
        PHKK(1,IHKK)=PQSDA1(N,1)
        PHKK(2,IHKK)=PQSDA1(N,2)
        PHKK(3,IHKK)=PQSDA1(N,3)
        PHKK(4,IHKK)=PQSDA1(N,4)
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
C                                     CHAIN 1 TARGET QUARK
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=132
        IDHKK(IHKK)=IDHKK(IHKKTD)
        JMOHKK(1,IHKK)=IHKKTD
        JMOHKK(2,IHKK)=JMOHKK(1,IHKKTD)
        JDAHKK(1,IHKK)=IHKK+1
        JDAHKK(2,IHKK)=IHKK+1
        PHKK(1,IHKK)=PQSDA2(N,1)
        PHKK(2,IHKK)=PQSDA2(N,2)
        PHKK(3,IHKK)=PQSDA2(N,3)
        PHKK(4,IHKK)=PQSDA2(N,4)
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
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
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
        MHKKSD(N)=IHKK
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
C                                   CHAIN 2 projectile sea antidiquark
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=121
        IDHKK(IHKK)=IDHKK(IHKKPD)
        JMOHKK(1,IHKK)=IHKKPD
        JMOHKK(2,IHKK)=JMOHKK(1,IHKKPD)
        JDAHKK(1,IHKK)=IHKK+2
        JDAHKK(2,IHKK)=IHKK+2
        PHKK(1,IHKK)=PQSDB1(N,1)
        PHKK(2,IHKK)=PQSDB1(N,2)
        PHKK(3,IHKK)=PQSDB1(N,3)
        PHKK(4,IHKK)=PQSDB1(N,4)
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
 
C                                     CHAIN 2 TARGET diquark
        NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
        IHKK=NHKK
        ISTHKK(IHKK)=132
        IDHKK(IHKK)=IDHKK(IHKKTO)
        JMOHKK(1,IHKK)=IHKKTO
        JMOHKK(2,IHKK)=JMOHKK(1,IHKKTO)
        JDAHKK(1,IHKK)=IHKK+1
        JDAHKK(2,IHKK)=IHKK+1
        PHKK(1,IHKK)=PQSDB2(N,1)
        PHKK(2,IHKK)=PQSDB2(N,2)
        PHKK(3,IHKK)=PQSDB2(N,3)
        PHKK(4,IHKK)=PQSDB2(N,4)
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
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
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
        MHKKSD(N)=IHKK
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
*     sea diquark pair!
C  AND PUT IT INTO THE HISTOGRAM
C
        AMCSD1(N)=AMCH1
        AMCSD2(N)=AMCH2
        GACSD1(N)=QECH1/AMCH1
        BGXSD1(N)=QTXCH1/AMCH1
        BGYSD1(N)=QTYCH1/AMCH1
        BGZSD1(N)=QTZCH1/AMCH1
        GACSD2(N)=QECH2/AMCH2
        BGXSD2(N)=QTXCH2/AMCH2
        BGYSD2(N)=QTYCH2/AMCH2
        BGZSD2(N)=QTZCH2/AMCH2
        NCHSD1(N)=NNCH1
        NCHSD2(N)=NNCH2
        IJCSD1(N)=IJNCH1
        IJCSD2(N)=IJNCH2
        IF (IPEV.GE.2) WRITE(6,'(A/I10,4F12.7,5I5/10X,4F12.6/10X,6F12.6,
     +4I5/8F15.5/8F15.5/2I5)') ' SD / FINAL PRINT',N
C    +, XPSQ
C    + (IXSPR),XPSAQ(IXSPR),XTSQ(IXSTA),XTSAQ(IXSTA), IPSQ(IXSPR),IPPV1
C    +  (IXSPR),IPPV2(IXSPR),ITSQ(IXSTA),ITSAQ(IXSTA), AMCSD1(N),AMCSD2
C    +  (N),GACSD1(N),GACSD2(N), BGXSD1(N),BGYSD1(N),BGZSD1(N), BGXSD2
C    +  (N),BGYSD2(N),BGZSD2(N), NCHSD1(N),NCHSD2(N),IJCSD1(N),IJCSD2
C    +  (N), (PQSDA1(N,JU),PQSDA2(N,JU),PQSDB1(N,JU), PQSDB2(N,JU),JU=1,
C    +  4),
C    +  IXSPR,IXSTA
                                                                GO TO 20
C***                     TREATMENT OF REJECTED SEA-SEA INTERACTIONS
   11   CONTINUE
        NCHSD1(N)=99
        NCHSD2(N)=99
        XPVD(JNUCPR)=XPVD(JNUCPR) + XPSQ(IXSPR) + XPSAQ(IXSPR)
        XTVD(JNUCTA)=XTVD(JNUCTA) + XTSAQ(IXSTA) + XTSQ(IXSTA)
      ISSQQ=ABS(ITSAQ(IXSTA))
      JSSQQ=ABS(ITSAQ2(IXSTA))
      IF(ISSQQ.EQ.3.AND.JSSQQ.EQ.3)THEN
        ISDRE(3)=ISDRE(3)+1
      ELSEIF(ISSQQ.EQ.3.OR.JSSQQ.EQ.3)THEN
        ISDRE(2)=ISDRE(2)+1
      ELSE
        ISDRE(1)=ISDRE(1)+1
      ENDIF
   20   CONTINUE
   10 CONTINUE
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     DEBUG SUBCHK
C     END DEBUG
      SUBROUTINE HADRSD
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C-------------------------
C
C                       hadronize sea diquark - valence CHAINS
C
C                       ADD GENERATED HADRONS TO /ALLPAR/
C                          STARTING AT (NAUX + 1)
C                       AND TO /HKKEVT/ STARTING AT (NHKK + 1)
C
C---------------------------------------------------------
      COMMON /ZSEA/ZSEAAV,ZSEASU,ANZSEA
       COMMON/POPCCK/PDBCK,PDBSE,PDBSEU,
     *  IJPOCK,IREJCK,ICK4,IHAD4,ICK6,IHAD6
     *,IREJSE,ISE4,ISE6,IREJS3,ISE43,ISE63,IREJS0
     *,IHADA4,IHADA6,IREJSA,ISEA4,ISEA6,IREJA3,
     *ISEA43,ISEA63,IREJAO
*KEEP,INTMX.
      PARAMETER (INTMX=2488,INTMD=252)
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,DIQI.
      COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
     +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
     +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
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
*KEEP,ABRSD.
      COMMON /ABRSD/ AMCSD1(248),AMCSD2(248),GACSD1(248),GACSD2(248),
     +BGXSD1(248),BGYSD1(248),BGZSD1(248), BGXSD2(248),BGYSD2(248),
     +BGZSD2(248), NCHSD1(248),NCHSD2(248),IJCSD1(248),IJCSD2(248),
     +PQSDA1(248,4),PQSDA2(248,4), PQSDB1(248,4),PQSDB2(248,4)
*KEEP,LOZUO.
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/ ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),ZUOST(INTMX),
     +INTLO(INTMX),INLOSS(INTMX)
C  /LOZUO/
C       INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C------------------
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
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
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
      COMMON /DIQSUM/NDVUU,NDVUS,NDVSS,NVDUU,NVDUS,NVDSS,
     *               NDSUU,NDSUS,NDSSS,NSDUU,NSDUS,NSDSS,
     *               NDZUU,NDZUS,NDZSS,NZDUU,NZDUS,NZDSS 
     *              ,NADVUU,NADVUS,NADVSS,NAVDUU,NAVDUS,NAVDSS,
     *               NADSUU,NADSUS,NADSSS,NASDUU,NASDUS,NASDSS,
     *               NADZUU,NADZUS,NADZSS,NAZDUU,NAZDUS,NAZDSS 
C---------------------
      DIMENSION POJ(4),PAT(4)
      DATA NCALSD /0/
C-----------------------------------------------------------------------
      IF(IPHKK.GE.6)WRITE (6,'( A)') ' hadrsd'
      NCALSD=NCALSD+1
      DO 50 I=1,NSD
C-----------------------drop recombined chain pairs
        IF(NCHSD1(I).EQ.99.AND.NCHSD2(I).EQ.99) GO TO 50
        IS1=INTSD1(I)
        IS2=INTSD2(I)
C
        IF (IPCO.GE.6) WRITE (6,1000) IPSQ(IS1),IPSAQ(IS1),ITVQ(IS2),
     +  ITTV1(IS2),ITTV2(IS2), AMCSD1(I),AMCSD2(I),GACSD1(I),GACSD2(I),
     +  BGXSD1(I),BGYSD1(I),BGZSD1(I), BGXSD2(I),BGYSD2(I),BGZSD2(I),
     +  NCHSD1(I),NCHSD2(I),IJCSD1(I),IJCSD2(I), PQSDA1(I,4),PQSDA2
     +  (I,4),PQSDB1(I,4),PQSDB2(I,4)
 1000 FORMAT(10X,5I5,10F9.2/10X,4I5,4F12.4)
C
C++++++++++++++++++++++++++++++    CHAIN 1:  quark-diquark   +++++++++++
        IFB1=IPSQ(IS1)
        IFB2=ITSQ(IS2)
        IFB3=ITSQ2(IS2)
        DO 10 J=1,4
          POJ(J)=PQSDA1(I,J)
          PAT(J)=PQSDA2(I,J)
   10   CONTINUE
        IF((NCHSD1(I).NE.0.OR.NCHSD2(I).NE.0).AND.IP.NE.1)
     &  CALL SAPTRE(AMCSD1(I),GACSD1(I),BGXSD1(I),BGYSD1(I),BGZSD1(I),
     &              AMCSD2(I),GACSD2(I),BGXSD2(I),BGYSD2(I),BGZSD2(I))
C----------------------------------------------------------------
C----------------------------------------------------------------
C       WRITE (6,1244) POJ,PAT
C1244   FORMAT ('  V-D QUARK-DIQUARK POJ,PAT ',8E12.3)
*       IF(AMCSD1(I).LT.1.6)THEN
*         IF(NCHSD1(I).EQ.0)THEN
*           WRITE(6,'(A,F10.2,5I5)')' HADRSD AMCDS1(I),NCHSD1(I),I ',
*    +                AMCSD1(I),NCHSD1(I),IJCSD1(I),I,IS1,IS2
*           RETURN
*         ENDIF
*       ENDIF
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
C               I= number of valence chain
C               Target     Nr itt = IFROVT(INTSV2(I))
C          No of Glauber sea q at Target     JITT=JTSHS(ITT)
C      ITTT = IFROVT(INTSV2(I))
C      IF(INTSD2(I).GE.1)THEN
C      ITTT = IFROVT(INTSD2(I))
C      ELSE
       ITTT=0
C      ENDIF
C      IF(ITTT.GE.1)THEN
C      JITT=JTSHS(ITTT)
C      ELSE
       JITT=0
C      ENDIF
C       IF(NCHSV1(I).EQ.0)THEN
C      WRITE(6,'(A,3I5)')'HADRSV: I,ITTT,JITT ',
C    *                     I,ITTT,JITT
C       ENDIF
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
         IF(IFB2.LE.2.AND.IFB3.LE.2)THEN
	   NSDUU=NSDUU+1
	 ELSEIF((IFB2.EQ.3.AND.IFB3.LE.2).OR.
     *  	 (IFB3.EQ.3.AND.IFB2.LE.2))THEN
	   NSDUS=NSDUS+1
	 ELSEIF(IFB2.EQ.3.AND.IFB3.EQ.3)THEN
	   NSDSS=NSDSS+1
	 ENDIF  
        IF((NCHSD1(I).NE.0))
     *  CALL HADJET(NHAD,AMCSD1(I),POJ,PAT,GACSD1(I),
     *  BGXSD1(I), BGYSD1
     +  (I),BGZSD1(I),IFB1,IFB2,IFB3,IFB4, IJCSD1(I),
     *  IJCSD1(I),4,NCHSD1
     +  (I),17)
C---------------------------------------------------------------
        AACK=FLOAT(ICK4)/FLOAT(ICK4+IHAD4+1)
        IF((NCHSD1(I).EQ.0))THEN
          ZSEAWU=RNDM(BB)*2.D0*ZSEAAV
        RSEACK=FLOAT(JITT)*PDBSE +ZSEAWU*PDBSEU
          IF(IPCO.GE.1)WRITE(6,'(2A,I5,2F10.3)')'HADJSE JITT,',
     *    'RSEACK,PDBSE 4 dpmnuc5 ',
     +    JITT,RSEACK,PDBSE
          IREJSS=5
          IF(RNDM(V).LE.RSEACK)THEN
            IREJSS=2
            IF(AMCSD1(I).GT.2.3D0)THEN
              IREJSS=0
              CALL HADJSE(NHAD,AMCSD1(I),POJ,PAT,GACSD1(I),BGXSD1(I),
     *        BGYSD1
     +        (I),BGZSD1(I),IFB1,IFB2,IFB3,IFB4, IJCSD1(I),IJCSD1(I),4,
     *        NCHSD1
     +        (I),3,IREJSS,IISSQQ)
              IF(IPCO.GE.1)WRITE(6,'(2A,I5,F10.3,I5)')'HADJSE JITT,',
     *        'RSEACK,IREJSS 4 dpmnuc5 NHAD',
     +        JITT,RSEACK,IREJSS,NHAD
            ENDIF
            IF(IREJSS.GE.1)THEN
              IF(IREJSS.EQ.1)IREJSE=IREJSE+1
              IF(IREJSS.EQ.3)IREJS3=IREJS3+1
              IF(IREJSS.EQ.2)IREJS0=IREJS0+1
              CALL HADJET(NHAD,AMCSD1(I),POJ,PAT,GACSD1(I),
     *        BGXSD1(I), BGYSD1
     +        (I),BGZSD1(I),IFB1,IFB2,IFB3,IFB4, IJCSD1(I),
     *        IJCSD1(I),4,NCHSD1
     +        (I),17)
              IHAD4=IHAD4+1
            ENDIF
            IF(IREJSS.EQ.0)THEN
              IF(IISSQQ.EQ.3)THEN
                ISE43=ISE43+1
              ELSE
                ISE4=ISE4+1
              ENDIF
            ENDIF
          ELSE
            CALL HADJET(NHAD,AMCSD1(I),POJ,PAT,GACSD1(I),
     *      BGXSD1(I), BGYSD1
     +      (I),BGZSD1(I),IFB1,IFB2,IFB3,IFB4, IJCSD1(I),
     *      IJCSD1(I),4,NCHSD1
     +      (I),17)
            IHAD4=IHAD4+1
          ENDIF
        ENDIF
C---------------------------------------------------------------
        ACOUSD=ACOUSD+1
        NHKKAU=NHKK+1
        DO 20 J=1,NHAD
          IF (NHKK.EQ.NMXHKK) THEN
            WRITE (6,'(A,2I5/A)') ' HADRSD: NHKK.EQ.NMXHKK ',NHKK,NMXHKK
            RETURN
          ENDIF
C         NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
C
          EHECC=SQRT(PXF(J)**2+PYF(J)**2+PZF(J)**2+AMF(J)**2)
          IF (ABS(EHECC-HEF(J)).GT.0.001) THEN
C           WRITE(6,'(2A/3I5,3E15.6)')
C    &            ' HADRSD / CHAIN 1 : CORRECT INCONSISTENT ENERGY ',
C    *            '  NCALSD, NHKK,NREF(J), HEF(J),EHECC, AMF(J)',
C    *            NCALSD, NHKK,NREF(J), HEF(J),EHECC, AMF(J)
          HEF(J)=EHECC
          ENDIF
          ANNSD=ANNSD+1
          EESD=EESD+HEF(J)
          PTSD=PTSD+SQRT(PXF(J)**2+PYF(J)**2)
C                               PUT NN-CMS HADRONS INTO /HKKEVT/
          ISTIST=1
          IF(IBARF(J).EQ.500)ISTIST=2
	  IF(IPCO.GE.3)WRITE(6,*)' HADRSD before HKKFIL J,NHAD',J,NHAD
          CALL HKKFIL(ISTIST,MPDGHA(NREF(J)),MHKKSD(I)-3,0,
     *                PXF(J),PYF(J),PZF(J),HEF(J),NHKKAU,IORMO(J),19)
          IF(IDHKK(NHKK).EQ.99999) WRITE (6,1030)NHKK,NREF(J), IDHKK
     +    (NHKK)
          IF (IPHKK.GE.2) WRITE(6,1010) NHKK, ISTHKK(NHKK),IDHKK(NHKK),
     +    JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +    (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
   20   CONTINUE
C	WRITE(6,*)' after 20 CONTINUE'
C       IF(NHAD.GT.0) THEN
C         JDAHKK(1,IMOHKK)=NHKKAU
C         JDAHKK(2,IMOHKK)=NHKK
C       ENDIF
C+++++++++++++++++++++++++++++   CHAIN 2: aquark - adiquark  +++++++++
        IFB1=IPSAQ(IS1)
        IFB2=ITSAQ(IS2)
        IFB3=ITSAQ2(IS2)
        IFB1=IABS(IFB1)+6
        IFB2=IABS(IFB2)+6
        IFB3=IABS(IFB3)+6
        DO 30 J=1,4
          POJ(J)=PQSDB2(I,J)
          PAT(J)=PQSDB1(I,J)
   30   CONTINUE
C
*       IF(AMCSD2(I).LT.1.6)THEN
*         IF(NCHSD2(I).EQ.0)THEN
*           WRITE(6,'(A,F10.2,5I5)')' HADRSD AMCSD2(I),NCHSD2(I),I ',
*    +                AMCSD2(I),NCHSD2(I),IJCSD2(I),I,IS1,IS2
*           RETURN
*         ENDIF
*       ENDIF
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
C               I= number of valence chain
C               Target     Nr itt = IFROVT(INTSV2(I))
C          No of Glauber sea q at Target     JITT=JTSHS(ITT)
C      ITTT = IFROVT(INTSV2(I))
C      WRITE(6,*)' INTSD2(I),I',INTSD2(I),I
C      IF(INTSD2(I).GE.1)THEN
C      ITTT = IFROVT(INTSD2(I))
C      ELSE
       ITTT=0
C      ENDIF
C      WRITE(6,*)' ITTT',ITTT
C      IF(ITTT.GE.1)THEN
C        JITT=JTSHS(ITTT)
C      ELSEIF(ITTT.EQ.0)THEN
	 JITT=0
C      ENDIF
C      WRITE(6,*)' JITT',JITT
C       IF(NCHSV1(I).EQ.0)THEN
C      WRITE(6,'(A,3I5)')'HADRSD: I,ITTT,JITT ',
C    *                     I,ITTT,JITT
C       ENDIF
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
         IF(IFB2.LE.8.AND.IFB3.LE.8)THEN
	   NASDUU=NASDUU+1
	 ELSEIF((IFB2.EQ.9.AND.IFB3.LE.8).OR.
     *  	 (IFB3.EQ.9.AND.IFB2.LE.8))THEN
	   NASDUS=NASDUS+1
	 ELSEIF(IFB2.EQ.9.AND.IFB3.EQ.9)THEN
	   NASDSS=NASDSS+1
	 ENDIF  
        IF((NCHSD2(I).NE.0))
     *  CALL HADJET(NHAD,AMCSD2(I),POJ,PAT,GACSD2(I),
     *  BGXSD2(I), BGYSD2
     +  (I),BGZSD2(I),IFB1,IFB2,IFB3,IFB4, IJCSD2(I),
     *  IJCSD2(I),4,NCHSD2
     +  (I),18)
C----------------------------------------------------------------
        AACK=FLOAT(ICK4)/FLOAT(ICK4+IHAD4+1)
        IF((NCHSD2(I).EQ.0))THEN
          ZSEAWU=RNDM(BB)*2.D0*ZSEAAV
        RSEACK=FLOAT(JITT)*PDBSE +ZSEAWU*PDBSEU
          IF(IPCO.GE.1)WRITE(6,'(2A,I5,2F10.3)')'HADJASE JITT,',
     *    'RSEACK,PDBSE ',
     +    JITT,RSEACK,PDBSE
          IREJSS=5
          IF(RNDM(V).LE.RSEACK)THEN
            IREJSS=2
            IF(AMCSD2(I).GT.2.3D0)THEN
              IREJSS=0
              CALL HADJASE(NHAD,AMCSD2(I),POJ,PAT,GACSD2(I),BGXSD2(I),
     *        BGYSD2
     +        (I),BGZSD2(I),IFB1,IFB2,IFB3,IFB4, IJCSD2(I),IJCSD2(I),4,
     *        NCHSD2
     +        (I),3,IREJSS,IISSQQ)
              IF(IPCO.GE.1)WRITE(6,'(2A,I5,F10.3,I5)')'HADJSE JITT,',
     *        'RSEACK,IREJSS ',
     +        JITT,RSEACK,IREJSS
            ENDIF
            IF(IREJSS.GE.1)THEN
              IF(IREJSS.EQ.1)IREJSA=IREJSA+1
              IF(IREJSS.EQ.3)IREJA3=IREJA3+1
              IF(IREJSS.EQ.2)IREJA0=IREJA0+1
        CALL HADJET(NHAD,AMCSD2(I),POJ,PAT,GACSD2(I),
     *  BGXSD2(I), BGYSD2
     +  (I),BGZSD2(I),IFB1,IFB2,IFB3,IFB4, IJCSD2(I),
     *  IJCSD2(I),4,NCHSD2
     +  (I),18)
              IHADA4=IHADA4+1
            ENDIF
            IF(IREJSS.EQ.0)THEN
              IF(IISSQQ.EQ.3)THEN
                ISEA43=ISEA43+1
              ELSE
                ISEA4=ISEA4+1
              ENDIF
            ENDIF
          ELSE
        CALL HADJET(NHAD,AMCSD2(I),POJ,PAT,GACSD2(I),
     *  BGXSD2(I), BGYSD2
     +  (I),BGZSD2(I),IFB1,IFB2,IFB3,IFB4, IJCSD2(I),
     *  IJCSD2(I),4,NCHSD2
     +  (I),18)
            IHADA4=IHADA4+1
          ENDIF
        ENDIF
C----------------------------------------------------------------
C                                   ADD HADRONS/RESONANCES INTO
C                                   COMMON /ALLPAR/ STARTING AT NAUX
        NHKKAU=NHKK+1
        DO 40 J=1,NHAD
          IF (NHKK.EQ.NMXHKK) THEN
            WRITE (6,'(A,2I5/A)') ' HADRSD: NHKK.EQ.NMXHKK ', NHKK,
     +      NMXHKK
            RETURN
          ENDIF
C         NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
C
          EHECC=SQRT(PXF(J)**2+PYF(J)**2+PZF(J)**2+AMF(J)**2)
        IF (ABS(EHECC-HEF(J)).GT.0.001) THEN
C             WRITE(6,'(2A/3I5,3E15.6)')
C    &            ' HADRSD / CHAIN 2 : CORRECT INCONSISTENT ENERGY ',
C    *            '  NCALSD, NHKK,NREF(J), HEF(J),EHECC, AMF(J)',
C    *            NCALSD, NHKK,NREF(J), HEF(J),EHECC, AMF(J)
          HEF(J)=EHECC
          ENDIF
          ANNSD=ANNSD+1
          EESD=EESD+HEF(J)
          PTSD=PTSD+SQRT(PXF(J)**2+PYF(J)**2)
C                               PUT NN-CMS HADRONS INTO /HKKEVT/
          ISTIST=1
          IF(IBARF(J).EQ.500)ISTIST=2
          CALL HKKFIL(ISTIST,MPDGHA(NREF(J)),MHKKSD(I),0,
     *                PXF(J),PYF(J),PZF(J),HEF(J),NHKKAU,IORMO(J),20)
          IF(IDHKK(NHKK).EQ.99999) WRITE (6,1030)NHKK,NREF(J), IDHKK
     +    (NHKK)
          IF (IPHKK.GE.2) WRITE(6,1010) NHKK, ISTHKK(NHKK),IDHKK(NHKK),
     +    JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +    (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
   40   CONTINUE
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
C----------------------------------------------------------------
C       formerly dpmdiqqq.f
C----------------------------------------------------------------
      SUBROUTINE DIQDZZ(ECM,XPSQ1,XPSAQ1,XPSQ2,XPSAQ2,
     *                  IPSQ1,IPSAQ1,IPSQ2,IPSAQ2,IREJDS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
*  define d-s chains (sea diquark - sea chains)
*  sqsq-sq and saqsaq-saq chains instead of q-aq and aq-q chains
      COMMON /ZSEA/ZSEAAV,ZSEASU,ANZSEA
       COMMON/POPCCK/PDBCK,PDBSE,PDBSEU,
     *  IJPOCK,IREJCK,ICK4,IHAD4,ICK6,IHAD6
     *,IREJSE,ISE4,ISE6,IREJS3,ISE43,ISE63,IREJS0
     *,IHADA4,IHADA6,IREJSA,ISEA4,ISEA6,IREJA3,
     *ISEA43,ISEA63,IREJAO
*KEEP,INTNEW.
      PARAMETER (INTMD=252)
      COMMON /INTNEZ/NDZ,NZD
      COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
C-------------------
*KEEP,ABRDS.
      COMMON /ABRDZ/ AMCDS1(INTMD),AMCDS2(INTMD),
     +GACDS1(INTMD),GACDS2(INTMD),
     +BGXDS1(INTMD),BGYDS1(INTMD),BGZDS1(INTMD), 
     +BGXDS2(INTMD),BGYDS2(INTMD),
     +BGZDS2(INTMD), NCHDS1(INTMD),NCHDS2(INTMD),
     +IJCDS1(INTMD),IJCDS2(INTMD),
     +PQDSA1(INTMD,4),PQDSA2(INTMD,4), 
     +PQDSB1(INTMD,4),PQDSB2(INTMD,4),
     +IPSQ(INTMD),IPSQQ2(INTMD),ITSQ(INTMD),
     +IPSAQ(INTMD),ISAQQ2(INTMD),ITSAQ(INTMD)
     +,IDZSS(INTMD)
C-------------------
*KEND.
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
      COMMON /DPAR/ ANAME(210),AM(210),GA(210),TAU(210),ICH(210),
     +IBAR(210),K1(210),K2(210)
C     COMMON /PCHARM/PCCCC
      COMMON/SEASU3/SEASQ
       COMMON /DIQREJ/IDIQRE(7),IDVRE(3),IVDRE(3),IDSRE(3),ISDRE(3),
     *IDZRE(3),IZDRE(3),IDIQRZ(7)
      COMMON /DIQSUM/NDVUU,NDVUS,NDVSS,NVDUU,NVDUS,NVDSS,
     *               NDSUU,NDSUS,NDSSS,NSDUU,NSDUS,NSDSS,
     *               NDZUU,NDZUS,NDZSS,NZDUU,NZDUS,NZDSS 
     *              ,NADVUU,NADVUS,NADVSS,NAVDUU,NAVDUS,NAVDSS,
     *               NADSUU,NADSUS,NADSSS,NASDUU,NASDUS,NASDSS,
     *               NADZUU,NADZUS,NADZSS,NAZDUU,NAZDUS,NAZDSS 
      PARAMETER (UMMM=0.3D0)
      PARAMETER (SMMM=0.5D0)
      PARAMETER (CMMM=1.3D0)
      DATA PC/0.0001D0/
*KEND.
C----------
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
      PC=PC1    
      PU1=PU/(2*PU+PS+PC)
      PS1=PS/(2*PU+PS+PC)
      IF(INICHA.EQ.0)THEN
        INICHA=1
        WRITE(6,4567)PC,BETCHA,PU1,PS1,SEASQ
 4567   FORMAT(' Charm chain ends DIQDZZ: PC,BETCHA,PU,PS,SEASQ',5F10.5)
      ENDIF
C----------------------------------------------------------------------
      IF(IPHKK.GE.3)WRITE (6,'( A)') ' diqdss'
      IREJDS=0
*  kinematics: is the mass of the adiquark-diquark chain big enough
*              to allow for fragmentation
      IPSQQ1=1.D0+RNDM(V1)*(2.D0+2.D0*SEASQ)
        RR=RNDM(V)
	IF(RR.LT.PC)IPSQQ1=4
C------------------
      ISAQQ1=-IPSQQ1
      AMDSQ1=XPSQ1*XPSQ2*ECM**2
      AMDSQ2=XPSAQ1*XPSAQ2*ECM**2
       IDIQRZ(1)=IDIQRZ(1)+1  
      IF(IPSQ1.EQ.3.AND.IPSQQ1.EQ.3)THEN
       IDIQRZ(2)=IDIQRZ(2)+1  
C       IF(AMDSQ2.LE.2.3.OR.AMDSQ1.LE.2.30) THEN
        IF(AMDSQ2.LE.6.6D0.OR.AMDSQ1.LE.6.60D0) THEN
       IDIQRZ(3)=IDIQRZ(3)+1  
       IDIQRZ(2)=IDIQRZ(2)-1  
       IDIQRZ(1)=IDIQRZ(1)-1  
          IREJDS=1
          RETURN
        ENDIF
      ELSEIF(IPSQ1.EQ.3.OR.IPSQQ1.EQ.3)THEN
       IDIQRZ(4)=IDIQRZ(4)+1  
C       IF(AMDSQ2.LE.1.9.OR.AMDSQ1.LE.1.90) THEN
        IF(AMDSQ2.LE.5.8D0.OR.AMDSQ1.LE.5.80D0) THEN
       IDIQRZ(5)=IDIQRZ(5)+1  
       IDIQRZ(4)=IDIQRZ(4)-1  
       IDIQRZ(1)=IDIQRZ(1)-1  
          IREJDS=1
          RETURN
        ENDIF
      ELSEIF(((IPSQ1.EQ.4).OR.(IPSQQ1.EQ.4)).AND.
     *       ((IPSQ1.EQ.3).OR.(IPSQQ1.EQ.3)))THEN
C       IF(AMDSQ2.LE.1.9.OR.AMDSQ1.LE.1.90) THEN
        IF(AMDSQ2.LE.30.8D0.OR.AMDSQ1.LE.30.80D0) THEN
          IREJDS=1
          RETURN
        ENDIF
      ELSEIF(IPSQ1.EQ.4.OR.IPSQQ1.EQ.4)THEN
C       IF(AMDSQ2.LE.1.9.OR.AMDSQ1.LE.1.90) THEN
        IF(AMDSQ2.LE.25.8.OR.AMDSQ1.LE.25.80) THEN
          IREJDS=1
          RETURN
        ENDIF
      ELSE
       IDIQRZ(6)=IDIQRZ(6)+1  
C        IF(AMDSQ2.LE.1.50.OR.AMDSQ1.LE.1.50) THEN
        IF(AMDSQ2.LE.3.9.OR.AMDSQ1.LE.3.9) THEN
       IDIQRZ(7)=IDIQRZ(7)+1  
       IDIQRZ(6)=IDIQRZ(6)-1  
       IDIQRZ(1)=IDIQRZ(1)-1  
          IREJDS=1
          RETURN
        ENDIF
      ENDIF
      NDZ=NDZ+1
      IF(NDZ.GE.INTMD)THEN
	IREJDS=1
	NDZ=NDZ-1
	RETURN
      ENDIF
C     WRITE(6,*)' DIQDZZ:IDIQRZ(1-7),NDZ ',(IDIQRZ(II),II=1,7),NDZ
      NCHDS1(NDZ)=0
      NCHDS2(NDZ)=0
C     WRITE(6,*)' DIQDZZ:NDZ,NCHDS1(NDZ),NCHDS2(NDZ) ',
C    * NDZ,NCHDS1(NDZ),NCHDS2(NDZ)
      IDZSS(NDZ)=0
C-------------------
C                                        KKEVDZ part
C-------------------
      IF(IPHKK.GE.3)WRITE (6,'( A,I10)') ' kkevdz',NDZ
      N=NDZ
C      DO 10 N=1,NDZ
C
C***                 4-MOMENTA OF PROJECTILE SEA-QUARK PAIRS IN NN-CMS
        IF(IPHKK.GE.7)WRITE(6,'(A,2I10)')' KKEVDZ N,NDZ',N,NDZ
C
        PRMOMZ=SQRT(ECM**2/4.-AM(1)**2)
        PSQPX=0.
        PSQPY=0.
        PSQPZ=XPSQ1*PRMOMZ
        PSQE=XPSQ1*ECM/2.
        PSAQPX=0.
        PSAQPY=0.
        PSAQPZ=XPSAQ1*PRMOMZ
        PSAQE=XPSAQ1*ECM/2.
C
C***                 4-MOMENTA OF TARGET QUARK-AQUARK PAIRS IN NN-CMS
C
        TSQPX=0.
        TSQPY=0.
        TSQPZ=-XPSQ2*PRMOMZ
        TSQE=XPSQ2*ECM/2.
        TSDQPX=0.
        TSDQPY=0.
        TSDQPZ=-XPSAQ2*PRMOMZ
        TSDQE=XPSAQ2*ECM/2.
C
C***           LORENTZ PARAMETER FOR CMS OF BOTH (sqsq-q) and
C              (saqsaq-diq)-SYSTEM
C              FROM PROJECTILE AND TARGET, RESP.
C
        PXXX=PSQPX + PSAQPX + TSQPX + TSDQPX
        PYYY=PSQPY + PSAQPY + TSQPY + TSDQPY
        PZZZ=PSQPZ + PSAQPZ + TSQPZ + TSDQPZ
        EEE =PSQE + PSAQE + TSQE + TSDQE
        PPTOTO=SQRT(PXXX**2+PYYY**2+PZZZ**2)
        AMMM=SQRT(ABS((EEE+PPTOTO)*(EEE-PPTOTO)))
        GAMMM=EEE/(AMMM+1.E-4)
        BGGGX=PXXX/(AMMM+1.E-4)
        BGGGY=PYYY/(AMMM+1.E-4)
        BGGGZ=PZZZ/(AMMM+1.E-4)
C
C***  SAMPLE PARTON-PT VALUES / DETERMINE PARTON 4-MOMENTA AND CHAIN MAS
C***                            IN THE REST FRAME DEFINED ABOVE
C
C       XPSQCM=XPSQ1/(XPSQ1+XPSAQ1)
C       XPSACM=1.0 - XPSQCM
C       XTSQCM=XPSQ2/(XPSQ2+XPSAQ2)
C       XTSACM=1.0 - XTSQCM
        XPSQCM=XPSQ1
        XPSACM=XPSAQ1
        XTSQCM=XPSQ2
        XTSACM=XPSAQ2
C
C***  SAMPLE PARTON-PT VALUES / DETERMINE PARTON 4-MOMENTA AND CHAIN MAS
C***                            IN THE REST FRAME DEFINED ABOVE
C
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
        IKVALA=0
        NSELPT=1
	IF(IPEV.GE.6)WRITE(6,'(A)')' KKEVDZ call selpt'
        CALL SELPT(
     +             PTXSQ1,PTYSQ1,PLQ1,EQ1,
     +             PTXSA1,PTYSA1,PLAQ1,EAQ1,
     +             PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +             PTXSQ2,PTYSQ2,PLQ2,EQ2,
     +             AMCH1,AMCH2,IREJDS,IKVALA,PTTQ1,PTTA1,
     +             PTTQ2,PTTA2,NSELPT)
C
        IF (IPEV.GE.7) WRITE(6,'(A/5X,5F12.5,I10)')
     +  'DS   AMMM,GAMMM,BGGGX,BGGGY,BGGGZ,IREJ ', AMMM,GAMMM,BGGGX,
     +  BGGGY,BGGGZ,IREJ
        IF (IREJDS.EQ.1) THEN
C         NDZ=NDZ-1
          IRDS13=IRDS13 + 1
          IF(IPEV.GE.2) THEN
            WRITE(6,'(A,I5)') ' KKEVDZ - IRDS13=',IRDS13
            WRITE(6,'(A/5E12.4/4(4E12.4/),2E12.4/2I5/4E12.4)')
     +      ' DS:  XPSQCM,XPSACM,XTSQCM,XTSACM,AMMM ...', XPSQCM,XPSACM,
     +      XTSQCM,XTSACM,AMMM, PTXSQ1,PTYSQ1,PLQ1,EQ1,PTXSA1,PTYSA1,
     +      PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +      AMCH1,AMCH2,IREJDS,IKVALA,PTTQ1,PTTA1
          ENDIF
                                                                GO TO 11
        ENDIF
C
C***  4-MOMENTA OF CHAINS IN THIS FRAME
C
        PTXCH1=PTXSQ1 + PTXSQ2
        PTYCH1=PTYSQ1 + PTYSQ2
        PTZCH1=PLQ1 + PLQ2
        ECH1=EQ1 + EQ2
        PTXCH2=PTXSA2 + PTXSA1
        PTYCH2=PTYSA2 + PTYSA1
        PTZCH2=PLAQ2 + PLAQ1
        ECH2=EAQ2 + EAQ1
C       WRITE(6,667)ECH1,ECH2,PTZCH1,PTZCH2
C 667 FORMAT(' DS ECH1,ECH2,PTZCH1,PTZCH2: ',4F10.3)
C
        IF (IPEV.GE.6) WRITE(6,'(A,5F12.5,I10/A,5F12.5/A,5F12.5)')
     +  ' DS: AMMM,GAMMM,BGGGX,BGGGY,BGGGZ,IREJDS ', AMMM,GAMMM,BGGGX,
     +  BGGGY,BGGGZ,IREJDS, '     AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',
     +  AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1,
     +  '     AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ', AMCH2,PTXCH2,PTYCH2,
     +  PTZCH2,ECH2

C
C  REPLACE SMALL MASS CHAINS BY PSEUDOSCALAR OR VECTOR MESONS OR OCTETT
C                                              OR DECUPLETT BARYONS
C  FIRST FOR CHAIN 1  (PROJ SEA-diquark - TAR QUARK)
C
        CALL ZOBCMA(IPSQ1,IPSQQ1,IPSQ2, IJNCH1,NNCH1,
     +  IREJDS,AMCH1,AMCH1N,1)
C***                            MASS BELOW OCTETT BARYON MASS
        IF(IREJDS.EQ.1) THEN
C         NDZ=NDZ-1
          IRDS11=IRDS11 + 1
                                                                 GOTO 11
        ENDIF
C                                 CORRECT KINEMATICS FOR CHAIN 1
C***                MOMENTUM CORRECTION FOR CHANGED MASS OF CHAIN 1
        IF(NNCH1.NE.0)
     +     CALL ZORMOM(AMMM,AMCH1,AMCH1N,AMCH2,
     +         XPSQ1,XPSAQ1,XPSAQ2,XPSQ2,
     +         PTXSQ1,PTYSQ1,PLQ1,EQ1,
     +         PTXSA1,PTYSA1,PLAQ1,EAQ1,
     +         PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +         PTXSQ2,PTYSQ2,PLQ2,EQ2,
     +         PTXCH1,PTYCH1,PTZCH1,ECH1, PTXCH2,PTYCH2,PTZCH2,ECH2,
     +                                                      IREJDS)
       IF(IREJDS.EQ.1)THEN
          GO TO 11
       ENDIF
C
        IF (IPEV.GE.6) WRITE(6,'(A,5F12.5,I10/A,5F12.5/A,5F12.5)')
     +  ' DS(2): AMMM,GAMMM,BGGGX,BGGGY,BGGGZ,IREJDS',AMMM,GAMMM,BGGGX,
     +  BGGGY,BGGGZ,IREJDS, '     AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',
     +  AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1,
     +  '     AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ', AMCH2,PTXCH2,PTYCH2,
     +  PTZCH2,ECH2
C
C  REPLACE SMALL MASS CHAINS BY octet or decuplet baryons
C       SECOND FOR CHAIN 2 (proj sadiquark - tar saquark)
C
        CALL ZOBCMA(IPSAQ1,ISAQQ1,IPSAQ2,
     +              IJNCH2,NNCH2,IREJDS,AMCH2,AMCH2N,2)
c  rejection of both s-s chains if mass of chain 2 too low
        IF(IREJDS.EQ.1) THEN
          IRDS12=IRDS12 + 1
C         NDZ=NDZ-1
          IF(IPEV.GE.2) THEN
            WRITE(6,1090) IRDS12
 1090       FORMAT(' KKEVDZ - IRDS12=',I5)
 1100       FORMAT(' DS - 1100', 6I5/2(4E12.4/),2E12.4)
          ENDIF
                                                                 GOTO 11
        ENDIF
C                                if AMCH2 changed in COBCMA/COMCMA
C                                ZORVAL corrects chain kinematics
C                                according to 2-body kinem.
C                                with fixed masses
        IF(NNCH2.NE.0) THEN
          AMCH2=AMCH2N
	  IORI=2
        CALL ZORVAL(AMMM,IREJDS,AMCH1,AMCH2, PTXCH1,PTYCH1,PTZCH1,ECH1,
     +    PTXCH2,PTYCH2,PTZCH2,ECH2,IORI)
          IF(IREJDS.EQ.1) THEN
*                           AMCH1N + AMCH2N > AMMM - 0.2
*                           reject event
C           NDZ=NDZ-1
            IRDS14=IRDS14+1
                                                                 GOTO 11
          ENDIF
C
          IF(IPEV.GE.6) THEN
            WRITE(6,'(A/3(1PE15.4),3I5)')
     +      ' DS - CALL ZORVAL: AMMM,AMCH1,AMCH2,NNCH1,NNCH2,IREJDS',
     +      AMMM, AMCH1, AMCH2, NNCH1, NNCH2, IREJDS
            WRITE(6,1050) AMMM,GAMMM,BGGGX,BGGGY,BGGGZ,IREJDS, AMCH1,
     +      PTXCH1,PTYCH1,PTZCH1,ECH1, AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2
 1050 FORMAT (' DS: AMMM,GAMMM,BGGGX,BGGGY,BGGGZ,IREJDS ',5F12.5,I10/
     +        '     AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',5F12.5/
     +        '     AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ',5F12.5)
          ENDIF
          IF(IREJDS.EQ.1) THEN
*                           AMCH1N + AMCH2N > AMMM - 0.2
*                           reject event
C           NDZ=NDZ-1
            IRDS14=IRDS14+1
                                                                 GOTO 11
          ENDIF
        ENDIF
C
C                           TRANSFORM BOTH CHAINS BACK INTO NN CMS
C
C                                   4-MOMENTA OF CHAINS
        QTXCH1=PTXCH1
        QTYCH1=PTYCH1
        QTZCH1=PTZCH1
        QECH1=ECH1
        QTXCH2=PTXCH2
        QTYCH2=PTYCH2
        QTZCH2=PTZCH2
        QECH2=ECH2
C       WRITE(6,887)QECH1,QECH2,QTZCH1,QTZCH2,AMCH1,AMCH2
  887 FORMAT( ' DS: QECH1,QECH2,QTZCH1,QTZCH2,AMCH1,AMCH2  ',6F10.2)

C                                   PARTONS AT ENDS OF CHAIN 1
C       CALL DALTRA(GAMMM,BGGGX,BGGGY,BGGGZ, PTXSQ1,PTYSQ1,PLQ1,EQ1,
C    +  PPPQ1, PQDSA1(N,1),PQDSA1(N,2),PQDSA1(N,3),PQDSA1(N,4) )
        PQDSA1(N,1)=PTXSQ1
        PQDSA1(N,2)=PTYSQ1
        PQDSA1(N,3)=PLQ1
        PQDSA1(N,4)=EQ1
        PQDSA2(N,1)=PTXSQ2
        PQDSA2(N,2)=PTYSQ2
        PQDSA2(N,3)=PLQ2
        PQDSA2(N,4)=EQ2

C                                   PARTONS AT ENDS OF CHAIN 2
        PQDSB2(N,1)=PTXSA2
        PQDSB2(N,2)=PTYSA2
        PQDSB2(N,3)=PLAQ2
        PQDSB2(N,4)=EAQ2

        PQDSB1(N,1)=PTXSA1
        PQDSB1(N,2)=PTYSA1
        PQDSB1(N,3)=PLAQ1
        PQDSB1(N,4)=EAQ1

C
C
C
C  NOW WE HAVE AN ACCEPTABLE SEA--VALENCE  EVENT
*     sea diquark pair!
C  AND PUT IT INTO THE HISTOGRAM
C
        IPSQ(N)=IPSQ1
        IPSQQ2(N)=IPSQQ1
        ITSQ(N)=IPSQ2
        IPSAQ(N)=IPSAQ1
        ISAQQ2(N)=ISAQQ1
        ITSAQ(N)=IPSAQ2
        AMCDS1(N)=AMCH1
        AMCDS2(N)=AMCH2
        GACDS1(N)=QECH1/AMCH1
        BGXDS1(N)=QTXCH1/AMCH1
        BGYDS1(N)=QTYCH1/AMCH1
        BGZDS1(N)=QTZCH1/AMCH1
        GACDS2(N)=QECH2/AMCH2
        BGXDS2(N)=QTXCH2/AMCH2
        BGYDS2(N)=QTYCH2/AMCH2
        BGZDS2(N)=QTZCH2/AMCH2
        NCHDS1(N)=NNCH1
        NCHDS2(N)=NNCH2
        IJCDS1(N)=IJNCH1
        IJCDS2(N)=IJNCH2
        IF (IPEV.GE.3) WRITE(6,'(A/I10,4F12.7,5I5/10X,4F12.6/10X,6F12.6,
     +4I5/8F15.5/                8F15.5)') ' DS / FINAL PRINT',N
                                                                GO TO 20
C***                     TREATMENT OF REJECTED SEA-SEA INTERACTIONS
   11   CONTINUE
        NCHDS1(N)=99
        NCHDS2(N)=99
C                               28.10.96
	NDZ=NDZ-1
	IF(NDZ.LT.0)THEN
	  NDZ=NDZ+1
	ENDIF
      ISSQQ=IPSQ1
      JSSQQ=IPSQQ1
      IF(ISSQQ.EQ.3.AND.JSSQQ.EQ.3)THEN
        IDZRE(3)=IDZRE(3)+1
      ELSEIF(ISSQQ.EQ.3.OR.JSSQQ.EQ.3)THEN
        IDZRE(2)=IDZRE(2)+1
      ELSE
        IDZRE(1)=IDZRE(1)+1
      ENDIF
   20   CONTINUE
   10 CONTINUE
C     WRITE(6,*)' DIQDZZ: IDZRE(1-3),NDZ ',(IDZRE(II),II=1,3),NDZ
C     WRITE(6,*)' DIQDZZ:NDZ,NCHDS1(NDZ),NCHDS2(NDZ) ',
C    * NDZ,NCHDS1(NDZ),NCHDS2(NDZ)
      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HADRDZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C-------------------------
C
C                       hadronize sea diquark - valence CHAINS
C
C                       ADD GENERATED HADRONS TO /ALLPAR/
C                          STARTING AT (NAUX + 1)
C                       AND TO /HKKEVT/ STARTING AT (NHKK + 1)
C
C---------------------------------------------------------
      COMMON /ZSEA/ZSEAAV,ZSEASU,ANZSEA
       COMMON/POPCCK/PDBCK,PDBSE,PDBSEU,
     *  IJPOCK,IREJCK,ICK4,IHAD4,ICK6,IHAD6
     *,IREJSE,ISE4,ISE6,IREJS3,ISE43,ISE63,IREJS0
     *,IHADA4,IHADA6,IREJSA,ISEA4,ISEA6,IREJA3,
     *ISEA43,ISEA63,IREJAO
*KEEP,INTNEW.
      PARAMETER (INTMD=252) 
      COMMON /INTNEZ/ NDZ,NZD
C-------------------
*KEEP,ABRDS.
      COMMON /ABRDZ/ AMCDS1(INTMD),AMCDS2(INTMD),
     +GACDS1(INTMD),GACDS2(INTMD),
     +BGXDS1(INTMD),BGYDS1(INTMD),BGZDS1(INTMD), 
     +BGXDS2(INTMD),BGYDS2(INTMD),
     +BGZDS2(INTMD), NCHDS1(INTMD),NCHDS2(INTMD),
     +IJCDS1(INTMD),IJCDS2(INTMD),
     +PQDSA1(INTMD,4),PQDSA2(INTMD,4), 
     +PQDSB1(INTMD,4),PQDSB2(INTMD,4),
     +IPSQ(INTMD),IPSQQ2(INTMD),ITSQ(INTMD),
     +IPSAQ(INTMD),ISAQQ2(INTMD),ITSAQ(INTMD)
     +,IDZSS(INTMD)
*KEEP,INTMX.
      PARAMETER (INTMX=2488)
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,DIQI.
C     COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
C    +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
C    +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
C    +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
*KEEP,DXQX.
C     INCLUDE (XQXQ)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON /DXQX/ XPVQ(248),XPVD(248),XTVQ(248),XTVD(248), XPSQ
     +(INTMX),XPSAQ(INTMX),XTSQ(INTMX),XTSAQ(INTMX)
     *                ,XPSU(248),XTSU(248)
     *                ,XPSUT(248),XTSUT(248)
 
C-------------------
*KEEP,LOZUO.
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/ ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),ZUOST(INTMX),
     +INTLO(INTMX),INLOSS(INTMX)
C  /LOZUO/
C       INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C------------------
*KEEP,HKKEVT.
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
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
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
      COMMON /DIQSUM/NDVUU,NDVUS,NDVSS,NVDUU,NVDUS,NVDSS,
     *               NDSUU,NDSUS,NDSSS,NSDUU,NSDUS,NSDSS,
     *               NDZUU,NDZUS,NDZSS,NZDUU,NZDUS,NZDSS 
     *              ,NADVUU,NADVUS,NADVSS,NAVDUU,NAVDUS,NAVDSS,
     *               NADSUU,NADSUS,NADSSS,NASDUU,NASDUS,NASDSS,
     *               NADZUU,NADZUS,NADZSS,NAZDUU,NAZDUS,NAZDSS 
*KEEP,INTNEW.
      COMMON /INTNEW/ NVV,NSV,NVS,NSS,NDV,NVD,NDS,NSD,
     +IXPV,IXPS,IXTV,IXTS, INTVV1(248),
     +INTVV2(248),INTSV1(248),INTSV2(248), INTVS1(248),INTVS2(248),
     +INTSS1(INTMX),INTSS2(INTMX),
     +INTDV1(248),INTDV2(248),INTVD1(248),INTVD2(248),
     +INTDS1(INTMD),INTDS2(INTMD),INTSD1(INTMD),INTSD2(INTMD)
C---------------------
      DIMENSION POJ(4),PAT(4)
      DATA NCALDS /0/
C     IPHKK=3
C-----------------------------------------------------------------------
      IF(IPHKK.GE.3)WRITE (6,'( A,4I10)') ' hadrdz',NDZ,NZD,
     *         NCHDS1(1),NCHDS2(1)
      NCALDS=NCALDS+1
      DO 50 I=1,NDZ
C-----------------------drop recombined chain pairs
        IF(NCHDS1(I).EQ.99.AND.NCHDS2(I).EQ.99) GO TO 50
        IS1=I
        IS2=I
C
        IF (IPCO.GE.6) WRITE (6,*)'IPSQ(IS1),IPSAQ(IS1),ITSQ(IS2)',
     +  'ITSAQ(IS2), AMCDS1(I),AMCDS2(I),GACDS1(I),GACDS2(I)',
     +  'BGXDS1(I),BGYDS1(I),BGZDS1(I), BGXDS2(I),BGYDS2(I),BGZDS2(I)',
     +  'NCHDS1(I),NCHDS2(I),IJCDS1(I),IJCDS2(I), PQDSA1(I,4),PQDSA2',
     +  '(I,4),PQDSB1(I,4),PQDSB2(I,4)',
     *  IPSQ(IS1),IPSAQ(IS1),ITSQ(IS2),
     +  ITSAQ(IS2), AMCDS1(I),AMCDS2(I),GACDS1(I),GACDS2(I),
     +  BGXDS1(I),BGYDS1(I),BGZDS1(I), BGXDS2(I),BGYDS2(I),BGZDS2(I),
     +  NCHDS1(I),NCHDS2(I),IJCDS1(I),IJCDS2(I), PQDSA1(I,4),PQDSA2
     +  (I,4),PQDSB1(I,4),PQDSB2(I,4)
 1000 FORMAT(10X,4I5,10F9.2/10X,4I5,4F12.4)
C
C
C++++++++++++++++++++++++++++++    CHAIN 1:  diquark-quark   +++++++++++
        IFB1=IPSQ(IS1)
        IFB2=IPSQQ2(IS1)
        IFB3=ITSQ(IS2)
        DO 10 J=1,4
          POJ(J)=PQDSA1(I,J)
          PAT(J)=PQDSA2(I,J)
   10   CONTINUE
        IF((NCHDS1(I).NE.0.OR.NCHDS2(I).NE.0).AND.IP.NE.1)
     &  CALL SAPTRE(AMCDS1(I),GACDS1(I),BGXDS1(I),BGYDS1(I),BGZDS1(I),
     &              AMCDS2(I),GACDS2(I),BGXDS2(I),BGYDS2(I),BGZDS2(I))
C----------------------------------------------------------------
C----------------------------------------------------------------
        IF (IPCO.GE.6)WRITE (6,1244) POJ,PAT
 1244   FORMAT ('  D-V QUARK-DIQUARK POJ,PAT ',8E12.3)
*       IF(AMCDS1(I).LT.1.6)THEN
*         IF(NCHDS1(I).EQ.0)THEN
*           WRITE(6,'(A,F10.2,5I5)')' HADRDZ AMCDS1(I),NCHDS1(I),I ',
*    +                AMCDS1(I),NCHDS1(I),IJCDS1(I),I,IS1,IS2
*           RETURN
*         ENDIF
*       ENDIF
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
C               I= number of valence chain
C               Projectile Nr ipp = IFROVP(INTVS1(I))
C          No of Glauber sea q at Projectile JIPP=JSSHS(IPP)
       IF(IPCO.GE.4)WRITE(6,*)' IPPP,INTVS1(I)',IPPP,INTVS1(I)
       IF(INTVS1(I).GT.0)THEN
         IPPP = IFROVP(INTVS1(I))
         IF(IPCO.GE.4)WRITE(6,*)' IPPP,INTVS1(I)',IPPP,INTVS1(I)
         JIPP=JSSHS(IPPP)
       ELSEIF(INTVS1(I).EQ.0)THEN
         JIPP=0
       ENDIF
       IF(IPCO.GE.4)WRITE(6,*)' JIPP,INTVS1(I)',JIPP,INTVS1(I)
C       IF(NCHVS2(I).EQ.0)THEN
C      WRITE(6,'(A,3I5)')'HADRVS: I,IPPP,JIPP ',
C    *                     I,IPPP,JIPP
C       ENDIF
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
         IF(IFB1.LE.2.AND.IFB2.LE.2)THEN
	   NDZUU=NDZUU+1
	 ELSEIF((IFB1.EQ.3.AND.IFB2.LE.2).OR.
     *  	 (IFB2.EQ.3.AND.IFB1.LE.2))THEN
	   NDZUS=NDZUS+1
	 ELSEIF(IFB1.EQ.3.AND.IFB2.EQ.3)THEN
	   NDZSS=NDZSS+1
	 ENDIF  
        IF((NCHDS1(I).NE.0))
     *  CALL HADJET(NHAD,AMCDS1(I),POJ,PAT,GACDS1(I),BGXDS1(I), BGYDS1
     +  (I),BGZDS1(I),IFB1,IFB2,IFB3,IFB4, IJCDS1(I),IJCDS1(I),6,NCHDS1
     +  (I),15)
C---------------------------------------------------------------------
        AACK=FLOAT(ICK6)/FLOAT(ICK6+IHAD6+1)
        IF((NCHDS1(I).EQ.0))THEN
          ZSEAWU=RNDM(BB)*2.D0*ZSEAAV
          RSEACK=FLOAT(JITT)*PDBSE +ZSEAWU*PDBSEU
          IF(IPCO.GE.1)WRITE(6,'(2A,I5,2F10.3)')'HADJSE JIPP,',
     *    'RSEACK,PDBSE 1 dpmdiqqq',
     +    JIPP,RSEACK,PDBSE
          IREJSS=5
          IF(RNDM(V).LE.RSEACK)THEN
            IREJSS=2
            IF(AMCDS1(I).GT.2.3D0)THEN
              IREJSS=0
              CALL HADJSE(NHAD,AMCDS1(I),POJ,PAT,GACDS1(I),BGXDS1(I),
     *        BGYDS1
     +        (I),BGZDS1(I),IFB1,IFB2,IFB3,IFB4, IJCDS1(I),IJCDS1(I),6,
     *        NCHDS1
     +        (I),6,IREJSS,IISSQQ)
              IF(IPCO.GE.1)WRITE(6,'(2A,I5,F10.3,I5)')'HADJSE JIPP,',
     *        'RSEACK,IREJSS 1 dpmdiqqq ',
     +        JIPP,RSEACK,IREJSS
            ENDIF
            IF(IREJSS.GE.1)THEN
              IF(IREJSS.EQ.1)IREJSE=IREJSE+1
              IF(IREJSS.EQ.3)IREJS3=IREJS3+1
              IF(IREJSS.EQ.2)IREJS0=IREJS0+1
        CALL HADJET(NHAD,AMCDS1(I),POJ,PAT,GACDS1(I),BGXDS1(I), BGYDS1
     +  (I),BGZDS1(I),IFB1,IFB2,IFB3,IFB4, IJCDS1(I),IJCDS1(I),6,NCHDS1
     +  (I),15)
              IHAD6=IHAD6+1
            ENDIF
            IF(IREJSS.EQ.0)THEN
              IF(IISSQQ.EQ.3)THEN
                ISE63=ISE63+1
              ELSE
                ISE6=ISE6+1
              ENDIF
            ENDIF
          ELSE
        CALL HADJET(NHAD,AMCDS1(I),POJ,PAT,GACDS1(I),BGXDS1(I), BGYDS1
     +  (I),BGZDS1(I),IFB1,IFB2,IFB3,IFB4, IJCDS1(I),IJCDS1(I),6,NCHDS1
     +  (I),15)
            IHAD6=IHAD6+1
          ENDIF
        ENDIF
C---------------------------------------------------------------------
        ACOUDZ=ACOUDZ+1
        NHKKAU=NHKK+1
        DO 20 J=1,NHAD
          IF (NHKK.EQ.NMXHKK) THEN
            WRITE (6,'(A,2I5/A)') ' HADRDZ: NHKK.EQ.NMXHKK ',NHKK,NMXHKK
            RETURN
          ENDIF
C         NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
C
          EHECC=SQRT(PXF(J)**2+PYF(J)**2+PZF(J)**2+AMF(J)**2)
        IF (ABS(EHECC-HEF(J)).GT.0.001) THEN
C           WRITE(6,'(2A/3I5,3E15.6)')
C    &            ' HADRDZ / CHAIN 1 : CORRECT INCONSISTENT ENERGY ',
C    *            '  NCALDS, NHKK,NREF(J), HEF(J),EHECC, AMF(J)',
C    *            NCALDS, NHKK,NREF(J), HEF(J),EHECC, AMF(J)
          HEF(J)=EHECC
          ENDIF
          ANNDZ=ANNDZ+1
          EEDZ=EEDZ+HEF(J)
          PTDZ=PTDZ+SQRT(PXF(J)**2+PYF(J)**2)
C                               PUT NN-CMS HADRONS INTO /HKKEVT/
          ISTIST=1
          IF(IBARF(J).EQ.500)ISTIST=2
          CALL HKKFIL(ISTIST,MPDGHA(NREF(J)),1,0,
     *                PXF(J),PYF(J),PZF(J),HEF(J),NHKKAU,IORMO(J),21)
          IF(IDHKK(NHKK).EQ.99999) WRITE (6,1030)NHKK,NREF(J), IDHKK
     +    (NHKK)
C         JMOHKK(1,NHKK)=MHKKSS(I)-3
          IF (IPHKK.GE.2) WRITE(6,1010) NHKK, ISTHKK(NHKK),IDHKK(NHKK),
     +    JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +    (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
   20   CONTINUE
C       IF(NHAD.GT.0) THEN
C         JDAHKK(1,IMOHKK)=NHKKAU
C         JDAHKK(2,IMOHKK)=NHKK
C       ENDIF
C+++++++++++++++++++++++++++++   CHAIN 2: adiquark - aquark  +++++++++
        IFB1=IPSAQ(IS1)
        IFB2=ISAQQ2(IS1)
        IFB3=ITSAQ(IS2)
        IFB1=IABS(IFB1)+6
        IFB2=IABS(IFB2)+6
        IFB3=IABS(IFB3)+6
        DO 30 J=1,4
          POJ(J)=PQDSB2(I,J)
          PAT(J)=PQDSB1(I,J)
   30   CONTINUE
*       IF(AMCDS2(I).LT.1.6)THEN
*         IF(NCHDS2(I).EQ.0)THEN
*           WRITE(6,'(A,F10.2,5I5)')' HADRDZ AMCDS2(I),NCHDS2(I),I ',
*    +                AMCDS2(I),NCHDS2(I),IJCDS2(I),I,IS1,IS2
*           RETURN
*         ENDIF
*       ENDIF
C
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
C               I= number of valence chain
C               Projectile Nr ipp = IFROVP(INTVS1(I))
C          No of Glauber sea q at Projectile JIPP=JSSHS(IPP)
       IF(INTVS1(I).GT.0)THEN
         IPPP = IFROVP(INTVS1(I))
         JITT=JSSHS(IPPP)
       ELSEIF(INTVS1(I).EQ.0)THEN
         JITT=0
       ENDIF
C       IF(NCHVS2(I).EQ.0)THEN
C      WRITE(6,'(A,3I5)')'HADRVS: I,IPPP,JIPP ',
C    *                     I,IPPP,JIPP
C       ENDIF
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
         IF(IFB1.LE.8.AND.IFB2.LE.8)THEN
	   NADZUU=NADZUU+1
	 ELSEIF((IFB1.EQ.9.AND.IFB2.LE.8).OR.
     *  	 (IFB2.EQ.9.AND.IFB1.LE.8))THEN
	   NADZUS=NADZUS+1
	 ELSEIF(IFB1.EQ.9.AND.IFB2.EQ.9)THEN
	   NADZSS=NADZSS+1
	 ENDIF  
        IF((NCHDS2(I).NE.0))
     *  CALL HADJET(NHAD,AMCDS2(I),POJ,PAT,GACDS2(I),BGXDS2(I), BGYDS2
     +  (I),BGZDS2(I),IFB1,IFB2,IFB3,IFB4, IJCDS2(I),IJCDS2(I),6,NCHDS2
     +  (I),16)
C-----------------------------------------------------------------------
        IF((NCHDS2(I).EQ.0))THEN
          ZSEAWU=RNDM(BB)*2.D0*ZSEAAV
          RSEACK=FLOAT(JITT)*PDBSE +ZSEAWU*PDBSEU
          IF(IPCO.GE.1)WRITE(6,'(2A,I5,2F10.3)')'HADJSE JIPP,',
     *    'RSEACK,PDBSE ',
     +    JIPP,RSEACK,PDBSE
          IREJSS=5
          IF(RNDM(V).LE.RSEACK)THEN
            IREJSS=2
            IF(AMCDS2(I).GT.2.3D0)THEN
              IREJSS=0
              CALL HADJASE(NHAD,AMCDS2(I),POJ,PAT,GACDS2(I),BGXDS2(I),
     *        BGYDS2
     +        (I),BGZDS2(I),IFB1,IFB2,IFB3,IFB4, IJCDS2(I),IJCDS2(I),6,
     *        NCHDS2
     +        (I),6,IREJSS,IISSQQ)
              IF(IPCO.GE.1)WRITE(6,'(2A,I5,F10.3,I5)')'HADJSE JIPP,',
     *        'RSEACK,IREJSS ',
     +        JIPP,RSEACK,IREJSS
            ENDIF
            IF(IREJSS.GE.1)THEN
              IF(IREJSS.EQ.1)IREJSA=IREJSA+1
              IF(IREJSS.EQ.3)IREJA3=IREJA3+1
              IF(IREJSS.EQ.2)IREJA0=IREJA0+1
        CALL HADJET(NHAD,AMCDS2(I),POJ,PAT,GACDS2(I),BGXDS2(I), BGYDS2
     +  (I),BGZDS2(I),IFB1,IFB2,IFB3,IFB4, IJCDS2(I),IJCDS2(I),6,NCHDS2
     +  (I),16)
              IHADA6=IHADA6+1
            ENDIF
            IF(IREJSS.EQ.0)THEN
              IF(IISSQQ.EQ.3)THEN
                ISEA63=ISEA63+1
              ELSE
                ISEA6=ISEA6+1
              ENDIF
            ENDIF
          ELSE
        CALL HADJET(NHAD,AMCDS2(I),POJ,PAT,GACDS2(I),BGXDS2(I), BGYDS2
     +  (I),BGZDS2(I),IFB1,IFB2,IFB3,IFB4, IJCDS2(I),IJCDS2(I),6,NCHDS2
     +  (I),16)
            IHADA6=IHADA6+1
          ENDIF
        ENDIF
C-----------------------------------------------------------------------
C                                   ADD HADRONS/RESONANCES INTO
C                                   COMMON /ALLPAR/ STARTING AT NAUX
        NHKKAU=NHKK+1
        DO 40 J=1,NHAD
          IF (NHKK.EQ.NMXHKK) THEN
            WRITE (6,'(A,2I5/A)') ' HADRDZ: NHKK.EQ.NMXHKK ', NHKK,
     +      NMXHKK
            RETURN
          ENDIF
C         NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
C
          EHECC=SQRT(PXF(J)**2+PYF(J)**2+PZF(J)**2+AMF(J)**2)
        IF (ABS(EHECC-HEF(J)).GT.0.001) THEN
C             WRITE(6,'(2A/3I5,3E15.6)')
C    &            ' HADRDZ / CHAIN 2 : CORRECT INCONSISTENT ENERGY ',
C    *            '  NCALDS, NHKK,NREF(J), HEF(J),EHECC, AMF(J)',
C    *            NCALDS, NHKK,NREF(J), HEF(J),EHECC, AMF(J)
          HEF(J)=EHECC
          ENDIF
C
          ANNDZ=ANNDZ+1
          EEDZ=EEDZ+HEF(J)
          PTDZ=PTDZ+SQRT(PXF(J)**2+PYF(J)**2)
C                               PUT NN-CMS HADRONS INTO /HKKEVT/
          ISTIST=1
          IF(IBARF(J).EQ.500)ISTIST=2
          CALL HKKFIL(ISTIST,MPDGHA(NREF(J)),1,0,
     *                PXF(J),PYF(J),PZF(J),HEF(J),NHKKAU,IORMO(J),22)
          IF(IDHKK(NHKK).EQ.99999) WRITE (6,1030)NHKK,NREF(J), IDHKK
     +    (NHKK)
C         JMOHKK(1,NHKK)=MHKKSS(I)
          IF (IPHKK.GE.2) WRITE(6,1010) NHKK, ISTHKK(NHKK),IDHKK(NHKK),
     +    JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +    (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
   40   CONTINUE
C       IF(NHAD.GT.0) THEN
C         JDAHKK(1,IMOHKK)=NHKKAU
C         JDAHKK(2,IMOHKK)=NHKK
C       ENDIF
   50 CONTINUE
C----------------------------------------------------------------
C
C     IPHKK=0
      RETURN
 1010 FORMAT (I6,I4,5I6,9E10.2)
 1020 FORMAT (' HADRKK J.GT.NAUMAX SKIP NEXT PARTICLES ',3I10)
 1030 FORMAT (' NHKK,IDHKK(NHKK)  ',3I10)
      END
C
C*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE DIQZZD(ECM,XPSQ1,XPSAQ1,XPSQ2,XPSAQ2,
     *                  IPSQ1,IPSAQ1,IPSQ2,IPSAQ2,IREJSD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
*  define s-d chains (sea - sea diquark chains)
*  sq-sqsq and saq-saqsaq chains instead of q-aq and aq-q chains
      COMMON /ZSEA/ZSEAAV,ZSEASU,ANZSEA
       COMMON/POPCCK/PDBCK,PDBSE,PDBSEU,
     *  IJPOCK,IREJCK,ICK4,IHAD4,ICK6,IHAD6
     *,IREJSE,ISE4,ISE6,IREJS3,ISE43,ISE63,IREJS0
     *,IHADA4,IHADA6,IREJSA,ISEA4,ISEA6,IREJA3,
     *ISEA43,ISEA63,IREJAO
*KEEP,INTNEW.
      PARAMETER (INTMD=252) 
      COMMON /INTNEZ/ NDZ,NZD
      COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
C-------------------
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
      COMMON /DPAR/ ANAME(210),AM(210),GA(210),TAU(210),ICH(210),
     +IBAR(210),K1(210),K2(210)
C------------------
*KEEP,ABRSD.
      COMMON /ABRZD/ AMCSD1(INTMD),AMCSD2(INTMD),
     +GACSD1(INTMD),GACSD2(INTMD),
     +BGXSD1(INTMD),BGYSD1(INTMD),BGZSD1(INTMD), 
     +BGXSD2(INTMD),BGYSD2(INTMD),
     +BGZSD2(INTMD), NCHSD1(INTMD),NCHSD2(INTMD),
     +IJCSD1(INTMD),IJCSD2(INTMD),
     +PQSDA1(INTMD,4),PQSDA2(INTMD,4), 
     +PQSDB1(INTMD,4),PQSDB2(INTMD,4),
     +IPSQ(INTMD),ITSQ(INTMD),ITSQ2(INTMD),
     +IPSAQ(INTMD),ITSAQ(INTMD),ITSAQ2(INTMD)
     +,IZDSS(INTMD)
*KEND.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
      COMMON/SEASU3/SEASQ
       COMMON /DIQREJ/IDIQRE(7),IDVRE(3),IVDRE(3),IDSRE(3),ISDRE(3),
     *IDZRE(3),IZDRE(3),IDIQRZ(7)
      COMMON /DIQSUM/NDVUU,NDVUS,NDVSS,NVDUU,NVDUS,NVDSS,
     *               NDSUU,NDSUS,NDSSS,NSDUU,NSDUS,NSDSS,
     *               NDZUU,NDZUS,NDZSS,NZDUU,NZDUS,NZDSS 
     *              ,NADVUU,NADVUS,NADVSS,NAVDUU,NAVDUS,NAVDSS,
     *               NADSUU,NADSUS,NADSSS,NASDUU,NASDUS,NASDSS,
     *               NADZUU,NADZUS,NADZSS,NAZDUU,NAZDUS,NAZDSS 
C     COMMON /PCHARM/PCCCC
      PARAMETER (UMMM=0.3D0)
      PARAMETER (SMMM=0.5D0)
      PARAMETER (CMMM=1.3D0)
      DATA PC/0.0001D0/
*KEND.
C----------
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
      PC=PC1
      PU1=PU/(2*PU+PS+PC)
      PS1=PS/(2*PU+PS+PC)
      IF(INICHA.EQ.0)THEN
        INICHA=1
        WRITE(6,4567)PC,BETCHA,PU1,PS1,SEASQ
 4567   FORMAT(' Charm chain ends DIQZZD: PC,BETCHA,PU,PS,SEASQ',5F10.5)
      ENDIF
C----------------------------------------------------------------------
      IF(IPHKK.GE.3)WRITE (6,'( A)') ' diqssd'
      IREJSD=0
*  kinematics: is the mass of both chains big enough
*              to allow for fragmentation
C     IPSQQ2=1
C     RND=RNDM(V)
C     IF(RND.GT.0.62) IPSQQ2=2
C     IF(RND.LT.0.24) IPSQQ2=3
      IPSQQ2=1.D0+RNDM(V1)*(2.D0+2.D0*SEASQ)
        RR=RNDM(V)
	IF(RR.LT.PC)IPSQQ2=4
      ISAQQ2=-IPSQQ2
      AMSDQ1=XPSQ2*XPSQ1*ECM**2
      AMSDQ2=XPSAQ2*XPSAQ1*ECM**2
       IDIQRZ(1)=IDIQRZ(1)+1  
C
      IF(IPSQ2.EQ.3.AND.IPSQQ2.EQ.3)THEN
       IDIQRZ(2)=IDIQRZ(2)+1  
C       IF(AMSDQ2.LE.2.30.OR.AMSDQ1.LE.2.30) THEN
        IF(AMSDQ2.LE.6.6D0.OR.AMSDQ1.LE.6.6D0) THEN
       IDIQRZ(3)=IDIQRZ(3)+1  
       IDIQRZ(2)=IDIQRZ(2)-1  
       IDIQRZ(1)=IDIQRZ(1)-1  
          IREJSD=1
          RETURN
        ENDIF
      ELSEIF(IPSQ2.EQ.3.OR.IPSQQ2.EQ.3)THEN
       IDIQRZ(4)=IDIQRZ(4)+1  
C        IF(AMSDQ2.LE.1.9.OR.AMSDQ1.LE.1.90) THEN
        IF(AMSDQ2.LE.5.8D0.OR.AMSDQ1.LE.5.80D0) THEN
       IDIQRZ(5)=IDIQRZ(5)+1  
       IDIQRZ(4)=IDIQRZ(4)-1  
       IDIQRZ(1)=IDIQRZ(1)-1  
          IREJSD=1
          RETURN
        ENDIF
      ELSEIF(((IPSQ2.EQ.4).OR.(IPSQQ2.EQ.4)).AND.
     *       ((IPSQ2.EQ.3).OR.(IPSQQ2.EQ.3)))THEN
C       IF(AMDSQ2.LE.1.9.OR.AMDSQ1.LE.1.90) THEN
        IF(AMSDQ2.LE.30.8D0.OR.AMSDQ1.LE.30.80D0) THEN
          IREJSD=1
          RETURN
        ENDIF
      ELSEIF(IPSQ2.EQ.4.OR.IPSQQ2.EQ.4)THEN
C        IF(AMSDQ2.LE.1.9.OR.AMSDQ1.LE.1.90) THEN
        IF(AMSDQ2.LE.25.8D0.OR.AMSDQ1.LE.25.80D0) THEN
          IREJSD=1
          RETURN
        ENDIF
      ELSE
       IDIQRZ(6)=IDIQRZ(6)+1  
C        IF(AMSDQ2.LE.1.50.OR.AMSDQ1.LE.1.50) THEN
        IF(AMSDQ2.LE.3.9D0.OR.AMSDQ1.LE.3.9D0) THEN
       IDIQRZ(7)=IDIQRZ(7)+1  
       IDIQRZ(6)=IDIQRZ(6)-1  
       IDIQRZ(1)=IDIQRZ(1)-1  
          IREJSD=1
          RETURN
        ENDIF
      ENDIF
      NZD=NZD+1
      IF(NZD.GE.INTMD)THEN
	IREJSD=1
	NZD=NZD-1
	RETURN
      ENDIF
C     WRITE(6,*)' DIQZZD:IDIQRZ(1-7),NZD ',(IDIQRZ(II),II=1,7),NZD
      NCHSD1(NZD)=0
      NCHSD2(NZD)=0
      IZDSS(NZD)=0
C-------------------
C                                 kkevzd part
C-------------------
      IF(IPHKK.GE.3)WRITE (6,'( A,3I10)') ' kkevzd',NZD,
     *            NCHSD1(1),NCHSD2(1)
      N=NZD
C      DO 10 N=1,NZD
C---------------------------drop recombined chain pairs
        IF(NCHSD1(N).EQ.99.AND.NCHSD2(N).EQ.99)GO TO 10
C
C***            4-MOMENTA OF projectile QUARK-DIQUARK PAIRS IN NN-CMS
C
        PRMOMZ=SQRT(ECM**2/4.-AM(1)**2)
        PSQPX=0.
        PSQPY=0.
        PSQPZ=XPSQ1*PRMOMZ
        PSQE=XPSQ1*ECM/2.
        PSDQPX=0.
        PSDQPY=0.
        PSDQPZ=XPSAQ1*PRMOMZ
        PSDQE=XPSAQ1*ECM/2.
C
C***                 4-MOMENTA OF TARGET QUARK-DIQUARK PAIRS IN NN-CMS
*
        TSQPX=0.
        TSQPY=0.
        TSQPZ=-XPSQ2*PRMOMZ
        TSQE=XPSQ2*ECM/2.
        TSAQPX=0.
        TSAQPY=0.
        TSAQPZ=-XPSAQ2*PRMOMZ
        TSAQE=XPSAQ2*ECM/2.
C
C***           LORENTZ PARAMETER FOR CMS OF BOTH (q-sqsq) and
C              (diq-saqsaq)-SYSTEM
C              FROM PROJECTILE AND TARGET, RESP.
C
        PXXX=TSQPX + TSAQPX + PSQPX + PSDQPX
        PYYY=TSQPY + TSAQPY + PSQPY + PSDQPY
        PZZZ=TSQPZ + TSAQPZ + PSQPZ + PSDQPZ
        EEE =TSQE + TSAQE + PSQE + PSDQE
        PPTOTO=SQRT(PXXX**2+PYYY**2+PZZZ**2)
        AMMM=SQRT(ABS((EEE+PPTOTO)*(EEE-PPTOTO)))
        GAMMM=EEE/(AMMM+1.E-4)
        BGGGX=PXXX/(AMMM+1.E-4)
        BGGGY=PYYY/(AMMM+1.E-4)
        BGGGZ=PZZZ/(AMMM+1.E-4)
C
C***  SAMPLE PARTON-PT VALUES / DETERMINE PARTON 4-MOMENTA AND CHAIN MAS
C***                            IN THE REST FRAME DEFINED ABOVE
C
        XTSQCM=XPSQ2
        XTSACM=XPSAQ2
        XPSQCM=XPSQ1
        XPSACM=XPSAQ1
C
C***  SAMPLE PARTON-PT VALUES / DETERMINE PARTON 4-MOMENTA AND CHAIN MAS
C***                            IN THE REST FRAME DEFINED ABOVE
C
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
        NSELPT=1
        IKVALA=0
	IF(IPEV.GE.6)WRITE(6,'(A)')' KKEVZD call selpt'
        CALL SELPT( PTXSQ1,PTYSQ1,PLQ1,
     +             EQ1,PTXSA1,PTYSA1,PLAQ1,EAQ1,
     +             PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +             PTXSQ2,PTYSQ2,PLQ2,EQ2,
     +             AMCH1,AMCH2,IREJSD,IKVALA,PTTQ1,PTTA1,
     +             PTTQ2,PTTA2,NSELPT)
c           WRITE(6,'(A,I5)') ' KKEVZD - IRSD13=',IRSD13
c           WRITE(6,'(A/5E12.4/4(4E12.4/),2E12.4/2I5/4E12.4)')
c    +      ' SD:  XPSQCM,XPSDCM,XTSQCM,XTSACM,AMMM,amch1,amch2 ',
c    +          XPVQCM,XPVDCM,
c    +      XTSQCM,XTSACM,AMMM, AMCH1,AMCH2
        IF (IPEV.GE.6) WRITE(6,'(A/5X,5F12.5,I10)')
     +  'SD   AMMM,GAMMM,BGGGX,BGGGY,BGGGZ,IREJSD ', AMMM,GAMMM,BGGGX,
     +  BGGGY,BGGGZ,IREJSD
        IF (IREJSD.EQ.1) THEN
          IRSD13=IRSD13 + 1
          IF(IPEV.GE.2) THEN
            WRITE(6,'(A,I5)') ' KKEVZD - IRSD13=',IRSD13
            WRITE(6,'(A/5E12.4/4(4E12.4/),2E12.4/2I5/4E12.4)')
     +      ' VD:  XPVQCM,XPVDCM,XTSQCM,XTSACM,AMMM ...', XPSQCM,XPSDCM,
     +      XTSQCM,XTSACM,AMMM, PTXSQ1,PTYSQ1,PLQ1,EQ1,PTXSA1,PTYSA1,
     +      PLAQ1,EAQ1, PTXSQ2,PTYSQ2,PLQ2,EQ2,PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +      AMCH1,AMCH2,IREJSD,IKVALA,PTTQ1,PTTA1

          ENDIF
                                                                GO TO 11
        ENDIF
C
C***  4-MOMENTA OF CHAINS IN THIS FRAME
C
        PTXCH1=PTXSQ1 + PTXSQ2
        PTYCH1=PTYSQ1 + PTYSQ2
        PTZCH1=PLQ1 + PLQ2
        ECH1=EQ1 + EQ2
        PTXCH2=PTXSA2 + PTXSA1
        PTYCH2=PTYSA2 + PTYSA1
        PTZCH2=PLAQ2 + PLAQ1
        ECH2=EAQ2 + EAQ1
C       WRITE(6,667)ECH1,ECH2,PTZCH1,PTZCH2
C 667 FORMAT(' SD ECH1,ECH2,PTZCH1,PTZCH2: ',4F10.3)
C
        IF (IPEV.GE.6) WRITE(6,'(A,5F12.5,I10/A,5F12.5/A,5F12.5)')
     +  ' SD: AMMM,GAMMM,BGGGX,BGGGY,BGGGZ,IREJSD ', AMMM,GAMMM,BGGGX,
     +  BGGGY,BGGGZ,IREJSD, '     AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',
     +  AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1,
     +  '     AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ', AMCH2,PTXCH2,PTYCH2,
     +  PTZCH2,ECH2

C
C  REPLACE SMALL MASS CHAINS BY PSEUDOSCALAR OR VECTOR MESONS OR OCTETT
C                                              OR DECUPLETT BARYONS
C  FIRST FOR CHAIN 1  (PROJ quark - tar sea-diquark)
C
        CALL ZOBCMA(IPSQ2,IPSQQ2,IPSQ1, IJNCH1,NNCH1,
     +  IREJSD,AMCH1,AMCH1N,1)
C***                            MASS BELOW OCTETT BARYON MASS
        IF(IREJSD.EQ.1) THEN
          IRSD11=IRSD11 + 1
                                                                 GOTO 11
        ENDIF
C                                 CORRECT KINEMATICS FOR CHAIN 1
C***                MOMENTUM CORRECTION FOR CHANGED MASS OF CHAIN 1
        IF(NNCH1.NE.0)
     +     CALL ZORMOM(AMMM,AMCH1,AMCH1N,AMCH2,
     +         XPSQ1,XPSAQ1,XPSAQ2,XPSQ2,
     +         PTXSQ1,PTYSQ1,PLQ1,EQ1,
     +         PTXSA1,PTYSA1,PLAQ1,EAQ1,
     +         PTXSA2,PTYSA2,PLAQ2,EAQ2,
     +         PTXSQ2,PTYSQ2,PLQ2,EQ2,
     +         PTXCH1,PTYCH1,PTZCH1,ECH1, PTXCH2,PTYCH2,PTZCH2,ECH2,
     +                                                      IREJSD)
       IF(IREJSD.EQ.1)THEN
          GO TO 11
       ENDIF
C
        IF (IPEV.GE.6) WRITE(6,'(A,5F12.5,I10/A,5F12.5/A,5F12.5)')
     +  ' SD(2): AMMM,GAMMM,BGGGX,BGGGY,BGGGZ,IREJSD',AMMM,GAMMM,BGGGX,
     +  BGGGY,BGGGZ,IREJSD, '     AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',
     +  AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1,
     +  '     AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ', AMCH2,PTXCH2,PTYCH2,
     +  PTZCH2,ECH2
C
C  REPLACE SMALL MASS CHAINS BY octet or decuplet baryons
C       SECOND FOR CHAIN 2 (proj saquark - tar sadiquark)
C
        CALL ZOBCMA(IPSAQ1,IPSAQ2,ISAQQ2,
     +              IJNCH2,NNCH2,IREJSD,AMCH2,AMCH2N,2)
c  rejection of both s-s chains if mass of chain 2 too low
        IF(IREJSD.EQ.1) THEN
          IRSD12=IRSD12 + 1
          IF(IPEV.GE.2) THEN
            WRITE(6,1090) IRSD12
            WRITE(6,1100) IPSAQ1,IPSAQ2,ISAQQ2,
     +                    IJNCH2,NNCH2,IREJSD,
     +      XPSQ1,XPSAQ1,XPSQCM,XTSACM, XPSQ2,XPSAQ2
     +      ,XTSQCM,XTSACM, AMCH2,AMCH2N
 1090       FORMAT(' KKEVZD - IRSD12=',I5)
 1100       FORMAT(' SD - 1100', 6I5/2(4E12.4/),2E12.4)
          ENDIF
                                                                 GOTO 11
        ENDIF
C                                if AMCH2 changed in COBCMA/COMCMA
C                                ZORVAL corrects chain kinematics
C                                according to 2-body kinem.
C                                with fixed masses
        IF(NNCH2.NE.0) THEN
          AMCH2=AMCH2N
	  IORI=1
        CALL ZORVAL(AMMM,IREJSD,AMCH1,AMCH2, PTXCH1,PTYCH1,PTZCH1,ECH1,
     +    PTXCH2,PTYCH2,PTZCH2,ECH2,IORI)
C
          IF(IPEV.GE.6) THEN
            WRITE(6,'(A/3(1PE15.4),3I5)')
     +      ' SD - CALL ZORVAL: AMMM,AMCH1,AMCH2,NNCH1,NNCH2,IREJSD',
     +      AMMM, AMCH1, AMCH2, NNCH1, NNCH2, IREJSD
            WRITE(6,1050) AMMM,GAMMM,BGGGX,BGGGY,BGGGZ,IREJSD, AMCH1,
     +      PTXCH1,PTYCH1,PTZCH1,ECH1, AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2
 1050 FORMAT (' SD: AMMM,GAMMM,BGGGX,BGGGY,BGGGZ,IREJSD ',5F12.5,I10/
     +        '     AMCH1,PTXCH1,PTYCH1,PTZCH1,ECH1 ',5F12.5/
     +        '     AMCH2,PTXCH2,PTYCH2,PTZCH2,ECH2 ',5F12.5)
          ENDIF
          IF(IREJSD.EQ.1) THEN
*                           AMCH1N + AMCH2N > AMMM - 0.2
*                           reject event
            IRSD14=IRSD14+1
                                                                 GOTO 11
          ENDIF
        ENDIF
C
C                           TRANSFORM BOTH CHAINS BACK INTO NN CMS
C
C                                   4-MOMENTA OF CHAINS
        QTXCH1=PTXCH1
        QTYCH1=PTYCH1
        QTZCH1=PTZCH1
        QECH1=ECH1

        QTXCH2=PTXCH2
        QTYCH2=PTYCH2
        QTZCH2=PTZCH2
        QECH2=ECH2
C       WRITE(6,887)QECH1,QECH2,QTZCH1,QTZCH2,AMCH1,AMCH2
C 887 FORMAT( ' SD: QECH1,QECH2,QTZCH1,QTZCH2,AMCH1,AMCH2  ',6F10.2)

C                                   PARTONS AT ENDS OF CHAIN 1
        PQSDA1(N,1)=PTXSQ1
        PQSDA1(N,2)=PTYSQ1
        PQSDA1(N,3)=PLQ1
        PQSDA1(N,4)=EQ1

        PQSDA2(N,1)=PTXSQ2
        PQSDA2(N,2)=PTYSQ2
        PQSDA2(N,3)=PLQ2
        PQSDA2(N,4)=EQ2

C                                   PARTONS AT ENDS OF CHAIN 2
        PQSDB2(N,1)=PTXSA2
        PQSDB2(N,2)=PTYSA2
        PQSDB2(N,3)=PLAQ2
        PQSDB2(N,4)=EAQ2

        PQSDB1(N,1)=PTXSA1
        PQSDB1(N,2)=PTYSA1
        PQSDB1(N,3)=PLAQ1
        PQSDB1(N,4)=EAQ1

C

C
C
C  NOW WE HAVE AN ACCEPTABLE SEA--VALENCE  EVENT
*     sea diquark pair!
C  AND PUT IT INTO THE HISTOGRAM
C
        IPSQ(N)=IPSQ1
        ITSQ(N)=IPSQ2
        ITSQ2(N)=IPSQQ2
        IPSAQ(N)=IPSAQ1
        ITSAQ(N)=IPSAQ2
        ITSAQ2(N)=ISAQQ2
        AMCSD1(N)=AMCH1
        AMCSD2(N)=AMCH2
        GACSD1(N)=QECH1/AMCH1
        BGXSD1(N)=QTXCH1/AMCH1
        BGYSD1(N)=QTYCH1/AMCH1
        BGZSD1(N)=QTZCH1/AMCH1
        GACSD2(N)=QECH2/AMCH2
        BGXSD2(N)=QTXCH2/AMCH2
        BGYSD2(N)=QTYCH2/AMCH2
        BGZSD2(N)=QTZCH2/AMCH2
        NCHSD1(N)=NNCH1
        NCHSD2(N)=NNCH2
        IJCSD1(N)=IJNCH1
        IJCSD2(N)=IJNCH2
        IF (IPEV.GE.2) WRITE(6,'(A/I10,4F12.7,5I5/10X,4F12.6/10X,6F12.6,
     +4I5/8F15.5/8F15.5/2I5)') ' SD / FINAL PRINT',N
                                                                GO TO 20
C***                     TREATMENT OF REJECTED SEA-SEA INTERACTIONS
   11   CONTINUE
        NCHSD1(N)=99
        NCHSD2(N)=99
C                                     28.10.96
	NZD=NZD-1
	IF(NZD.LT.0)THEN
	  NZD=NZD+1
	ENDIF
      ISSQQ=IPSQ2
      JSSQQ=IPSQQ2
      IF(ISSQQ.EQ.3.AND.JSSQQ.EQ.3)THEN
        IZDRE(3)=IZDRE(3)+1
      ELSEIF(ISSQQ.EQ.3.OR.JSSQQ.EQ.3)THEN
        IZDRE(2)=IZDRE(2)+1
      ELSE
        IZDRE(1)=IZDRE(1)+1
      ENDIF
   20   CONTINUE
   10 CONTINUE
C     WRITE(6,*)' DIQZZD: IZDRE(1-3),NZD ',(IZDRE(II),II=1,3),NZD
      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C ---------------------------------------------------------------
C ---------------------------------------------------------------
C ---------------------------------------------------------------

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HADRZD
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C-------------------------
C
C                       hadronize sea diquark - valence CHAINS
C
C                       ADD GENERATED HADRONS TO /ALLPAR/
C                          STARTING AT (NAUX + 1)
C                       AND TO /HKKEVT/ STARTING AT (NHKK + 1)
C
C---------------------------------------------------------
      COMMON /ZSEA/ZSEAAV,ZSEASU,ANZSEA
       COMMON/POPCCK/PDBCK,PDBSE,PDBSEU,
     *  IJPOCK,IREJCK,ICK4,IHAD4,ICK6,IHAD6
     *,IREJSE,ISE4,ISE6,IREJS3,ISE43,ISE63,IREJS0
     *,IHADA4,IHADA6,IREJSA,ISEA4,ISEA6,IREJA3,
     *ISEA43,ISEA63,IREJAO
      PARAMETER (INTMD=252)
      COMMON /INTNEZ/ NDZ,NZD
C-------------------
*KEEP,ABRSD.
      COMMON /ABRZD/ AMCSD1(INTMD),AMCSD2(INTMD),
     +GACSD1(INTMD),GACSD2(INTMD),
     +BGXSD1(INTMD),BGYSD1(INTMD),BGZSD1(INTMD), 
     +BGXSD2(INTMD),BGYSD2(INTMD),
     +BGZSD2(INTMD), NCHSD1(INTMD),NCHSD2(INTMD),
     +IJCSD1(INTMD),IJCSD2(INTMD),
     +PQSDA1(INTMD,4),PQSDA2(INTMD,4), 
     +PQSDB1(INTMD,4),PQSDB2(INTMD,4),
     +IPSQ(INTMD),ITSQ(INTMD),ITSQ2(INTMD),
     +IPSAQ(INTMD),ITSAQ(INTMD),ITSAQ2(INTMD)
     +,IZDSS(INTMD)
*KEEP,INTMX.
      PARAMETER (INTMX=2488)
*KEEP,IFROTO.
      COMMON /IFROTO/ IFROVP(248),ITOVP(248),IFROSP(INTMX), IFROVT(248),
     +ITOVT(248),IFROST(INTMX),JSSHS(INTMX),JTSHS(INTMX),JHKKNP(248),
     +JHKKNT
     +(248), JHKKPV(INTMX),JHKKPS(INTMX), JHKKTV(INTMX),JHKKTS(INTMX),
     +MHKKVV(INTMX),MHKKSS(INTMX), MHKKVS(INTMX),MHKKSV(INTMX),
     &                MHKKHH(INTMX),
     +MHKKDV(248),MHKKVD(248), MHKKDS(INTMD),MHKKSD(INTMD)
*KEEP,DIQI.
C     COMMON /DIQI/ IPVQ(248),IPPV1(248),IPPV2(248), ITVQ(248),ITTV1
C    +(248),ITTV2(248), IPSQ(INTMX),IPSQ2(INTMX),
C    +IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
C    +ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
*KEEP,DXQX.
C     INCLUDE (XQXQ)
* NOTE: INTMX set via INCLUDE(INTMX)
      COMMON /DXQX/ XPVQ(248),XPVD(248),XTVQ(248),XTVD(248), XPSQ
     +(INTMX),XPSAQ(INTMX),XTSQ(INTMX),XTSAQ(INTMX)
     *                ,XPSU(248),XTSU(248)
     *                ,XPSUT(248),XTSUT(248)
*KEEP,INTNEW.
C     COMMON /INTNEW/ NVV,NSV,NVS,NSS,NDV,NVD,NDS,NSD,
C    +IXPV,IXPS,IXTV,IXTS, INTVV1(248),
C    +INTVV2(248),INTSV1(248),INTSV2(248), INTVS1(248),INTVS2(248),
C    +INTSS1(INTMX),INTSS2(INTMX),
C    +INTDV1(248),INTDV2(248),INTVD1(248),INTVD2(248),
C    +INTDS1(INTMD),INTDS2(INTMD),INTSD1(INTMD),INTSD2(INTMD)
 
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
*KEEP,LOZUO.
      LOGICAL ZUOVP,ZUOSP,ZUOVT,ZUOST,INTLO,INLOSS
      COMMON /LOZUO/ ZUOVP(248),ZUOSP(INTMX),ZUOVT(248),ZUOST(INTMX),
     +INTLO(INTMX),INLOSS(INTMX)
C  /LOZUO/
C       INLOSS :  .FALSE. IF CORRESPONDING SEA-SEA INTERACTION
C                         REJECTED IN KKEVT
C------------------
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
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
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
      COMMON /DIQSUM/NDVUU,NDVUS,NDVSS,NVDUU,NVDUS,NVDSS,
     *               NDSUU,NDSUS,NDSSS,NSDUU,NSDUS,NSDSS,
     *               NDZUU,NDZUS,NDZSS,NZDUU,NZDUS,NZDSS 
     *              ,NADVUU,NADVUS,NADVSS,NAVDUU,NAVDUS,NAVDSS,
     *               NADSUU,NADSUS,NADSSS,NASDUU,NASDUS,NASDSS,
     *               NADZUU,NADZUS,NADZSS,NAZDUU,NAZDUS,NAZDSS 
C---------------------
*KEEP,INTNEW.
      COMMON /INTNEW/ NVV,NSV,NVS,NSS,NDV,NVD,NDS,NSD,
     +IXPV,IXPS,IXTV,IXTS, INTVV1(248),
     +INTVV2(248),INTSV1(248),INTSV2(248), INTVS1(248),INTVS2(248),
     +INTSS1(INTMX),INTSS2(INTMX),
     +INTDV1(248),INTDV2(248),INTVD1(248),INTVD2(248),
     +INTDS1(INTMD),INTDS2(INTMD),INTSD1(INTMD),INTSD2(INTMD)

      DIMENSION POJ(4),PAT(4)
      DATA NCALSD /0/
C     IPHKK=3
C-----------------------------------------------------------------------
      IF(IPHKK.GE.3)WRITE (6,'( A,4I10)') ' hadrzd',NDZ,NZD,
     *           NCHSD1(1),NCHSD2(1)
      NCALSD=NCALSD+1
      DO 50 I=1,NZD
C-----------------------drop recombined chain pairs
        IF(NCHSD1(I).EQ.99.AND.NCHSD2(I).EQ.99) GO TO 50
        IS1=I
        IS2=I
C
C       IF (IPCO.GE.6) WRITE (6,1000) IPSQ(IS1),IPSAQ(IS1),ITVQ(IS2),
C    +  ITTV1(IS2),ITTV2(IS2), AMCSD1(I),AMCSD2(I),GACSD1(I),GACSD2(I),
C    +  BGXSD1(I),BGYSD1(I),BGZSD1(I), BGXSD2(I),BGYSD2(I),BGZSD2(I),
C    +  NCHSD1(I),NCHSD2(I),IJCSD1(I),IJCSD2(I), PQSDA1(I,4),PQSDA2
C    +  (I,4),PQSDB1(I,4),PQSDB2(I,4)
 1000 FORMAT(10X,5I5,10F9.2/10X,4I5,4F12.4)
C
C++++++++++++++++++++++++++++++    CHAIN 1:  quark-diquark   +++++++++++
        IFB1=IPSQ(IS1)
        IFB2=ITSQ(IS2)
        IFB3=ITSQ2(IS2)
        DO 10 J=1,4
          POJ(J)=PQSDA1(I,J)
          PAT(J)=PQSDA2(I,J)
   10   CONTINUE
        IF((NCHSD1(I).NE.0.OR.NCHSD2(I).NE.0).AND.IP.NE.1)
     &  CALL SAPTRE(AMCSD1(I),GACSD1(I),BGXSD1(I),BGYSD1(I),BGZSD1(I),
     &              AMCSD2(I),GACSD2(I),BGXSD2(I),BGYSD2(I),BGZSD2(I))
C----------------------------------------------------------------
C----------------------------------------------------------------
C       WRITE (6,1244) POJ,PAT
C1244   FORMAT ('  V-D QUARK-DIQUARK POJ,PAT ',8E12.3)
*       IF(AMCSD1(I).LT.1.6)THEN
*         IF(NCHSD1(I).EQ.0)THEN
*           WRITE(6,'(A,F10.2,5I5)')' HADRZD AMCDS1(I),NCHSD1(I),I ',
*    +                AMCSD1(I),NCHSD1(I),IJCSD1(I),I,IS1,IS2
*           RETURN
*         ENDIF
*       ENDIF
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
C               I= number of valence chain
C               Target     Nr itt = IFROVT(INTSV2(I))
C          No of Glauber sea q at Target     JITT=JTSHS(ITT)
       IF(INTSV2(I).GT.0)THEN
         ITTT = IFROVT(INTSV2(I))
         JITT=JTSHS(ITTT)
       ELSEIF(INTSV2(I).EQ.0)THEN
         JITT=0
       ENDIF
C       IF(NCHSV1(I).EQ.0)THEN
C      WRITE(6,'(A,3I5)')'HADRSV: I,ITTT,JITT ',
C    *                     I,ITTT,JITT
C       ENDIF
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
         IF(IFB2.LE.2.AND.IFB3.LE.2)THEN
	   NZDUU=NZDUU+1
	 ELSEIF((IFB2.EQ.3.AND.IFB3.LE.2).OR.
     *  	 (IFB3.EQ.3.AND.IFB2.LE.2))THEN
	   NZDUS=NZDUS+1
	 ELSEIF(IFB2.EQ.3.AND.IFB3.EQ.3)THEN
	   NZDSS=NZDSS+1
	 ENDIF  
        IF((NCHSD1(I).NE.0))
     *  CALL HADJET(NHAD,AMCSD1(I),POJ,PAT,GACSD1(I),BGXSD1(I), BGYSD1
     +  (I),BGZSD1(I),IFB1,IFB2,IFB3,IFB4, IJCSD1(I),IJCSD1(I),4,NCHSD1
     +  (I),17)
C-------------------------------------------------------------------
        AACK=FLOAT(ICK4)/FLOAT(ICK4+IHAD4+1)
        IF((NCHSD1(I).EQ.0))THEN
          ZSEAWU=RNDM(BB)*2.D0*ZSEAAV
          RSEACK=FLOAT(JITT)*PDBSE +ZSEAWU*PDBSEU
          IF(IPCO.GE.1)WRITE(6,'(2A,I5,2F10.3)')'HADJSE JITT,',
     *    'RSEACK,PDBSE 2 dpmdiqqq ',
     +    JITT,RSEACK,PDBSE
          IREJSS=5
          IF(RNDM(V).LE.RSEACK)THEN
            IREJSS=2
            IF(AMCSD1(I).GT.2.3D0)THEN
              IREJSS=0
              CALL HADJSE(NHAD,AMCSD1(I),POJ,PAT,GACSD1(I),BGXSD1(I),
     *        BGYSD1
     +        (I),BGZSD1(I),IFB1,IFB2,IFB3,IFB4, IJCSD1(I),IJCSD1(I),4,
     *        NCHSD1
     +        (I),3,IREJSS,IISSQQ)
              IF(IPCO.GE.1)WRITE(6,'(2A,I5,F10.3,I5)')'HADJSE JITT,',
     *        'RSEACK,IREJSS 2 dpmdiqqq ',
     +        JITT,RSEACK,IREJSS
            ENDIF
            IF(IREJSS.GE.1)THEN
              IF(IREJSS.EQ.1)IREJSE=IREJSE+1
              IF(IREJSS.EQ.3)IREJS3=IREJS3+1
              IF(IREJSS.EQ.2)IREJS0=IREJS0+1
        CALL HADJET(NHAD,AMCSD1(I),POJ,PAT,GACSD1(I),BGXSD1(I), BGYSD1
     +  (I),BGZSD1(I),IFB1,IFB2,IFB3,IFB4, IJCSD1(I),IJCSD1(I),4,NCHSD1
     +  (I),17)
              IHAD4=IHAD4+1
            ENDIF
            IF(IREJSS.EQ.0)THEN
              IF(IISSQQ.EQ.3)THEN
                ISE43=ISE43+1
              ELSE
                ISE4=ISE4+1
              ENDIF
            ENDIF
          ELSE
        CALL HADJET(NHAD,AMCSD1(I),POJ,PAT,GACSD1(I),BGXSD1(I), BGYSD1
     +  (I),BGZSD1(I),IFB1,IFB2,IFB3,IFB4, IJCSD1(I),IJCSD1(I),4,NCHSD1
     +  (I),17)
            IHAD4=IHAD4+1
          ENDIF
        ENDIF
C-------------------------------------------------------------------
        ACOUZD=ACOUZD+1
        NHKKAU=NHKK+1
        DO 20 J=1,NHAD
          IF (NHKK.EQ.NMXHKK) THEN
            WRITE (6,'(A,2I5/A)') ' HADRZD: NHKK.EQ.NMXHKK ',NHKK,NMXHKK
            RETURN
          ENDIF
C         NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
C
          EHECC=SQRT(PXF(J)**2+PYF(J)**2+PZF(J)**2+AMF(J)**2)
          IF (ABS(EHECC-HEF(J)).GT.0.001) THEN
C           WRITE(6,'(2A/3I5,3E15.6)')
C    &            ' HADRZD / CHAIN 1 : CORRECT INCONSISTENT ENERGY ',
C    *            '  NCALSD, NHKK,NREF(J), HEF(J),EHECC, AMF(J)',
C    *            NCALSD, NHKK,NREF(J), HEF(J),EHECC, AMF(J)
          HEF(J)=EHECC
          ENDIF
          ANNZD=ANNZD+1
          EEZD=EEZD+HEF(J)
          PTZD=PTZD+SQRT(PXF(J)**2+PYF(J)**2)
C                               PUT NN-CMS HADRONS INTO /HKKEVT/
          ISTIST=1
          IF(IBARF(J).EQ.500)ISTIST=2
          CALL HKKFIL(ISTIST,MPDGHA(NREF(J)),1,0,
     *                PXF(J),PYF(J),PZF(J),HEF(J),NHKKAU,IORMO(J),23)
          IF(IDHKK(NHKK).EQ.99999) WRITE (6,1030)NHKK,NREF(J), IDHKK
     +    (NHKK)
C         JMOHKK(1,NHKK)=MHKKSS(I)-3
          IF (IPHKK.GE.2) WRITE(6,1010) NHKK, ISTHKK(NHKK),IDHKK(NHKK),
     +    JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +    (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
   20   CONTINUE
C       IF(NHAD.GT.0) THEN
C         JDAHKK(1,IMOHKK)=NHKKAU
C         JDAHKK(2,IMOHKK)=NHKK
C       ENDIF
C+++++++++++++++++++++++++++++   CHAIN 2: aquark - adiquark  +++++++++
        IFB1=IPSAQ(IS1)
        IFB2=ITSAQ(IS2)
        IFB3=ITSAQ2(IS2)
        IFB1=IABS(IFB1)+6
        IFB2=IABS(IFB2)+6
        IFB3=IABS(IFB3)+6
        DO 30 J=1,4
          POJ(J)=PQSDB2(I,J)
          PAT(J)=PQSDB1(I,J)
   30   CONTINUE
C
*       IF(AMCSD2(I).LT.1.6)THEN
*         IF(NCHSD2(I).EQ.0)THEN
*           WRITE(6,'(A,F10.2,5I5)')' HADRZD AMCSD2(I),NCHSD2(I),I ',
*    +                AMCSD2(I),NCHSD2(I),IJCSD2(I),I,IS1,IS2
*           RETURN
*         ENDIF
*       ENDIF
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
C               I= number of valence chain
C               Target     Nr itt = IFROVT(INTSV2(I))
C          No of Glauber sea q at Target     JITT=JTSHS(ITT)
       IF(INTSV2(I).GT.0)THEN
         ITTT = IFROVT(INTSV2(I))
         JITT=JTSHS(ITTT)
       ELSEIF(INTSV2(I).EQ.0)THEN
         JITT=0
       ENDIF
C       IF(NCHSV1(I).EQ.0)THEN
C      WRITE(6,'(A,3I5)')'HADRSV: I,ITTT,JITT ',
C    *                     I,ITTT,JITT
C       ENDIF
C------------------------------------------------------------------
C                                   check bookkeeping
C-----------------------------------------------------------------
         IF(IFB2.LE.8.AND.IFB3.LE.8)THEN
	   NAZDUU=NAZDUU+1
	 ELSEIF((IFB2.EQ.9.AND.IFB3.LE.8).OR.
     *  	 (IFB3.EQ.9.AND.IFB2.LE.8))THEN
	   NAZDUS=NAZDUS+1
	 ELSEIF(IFB2.EQ.9.AND.IFB3.EQ.9)THEN
	   NAZDSS=NAZDSS+1
	 ENDIF  
        IF((NCHSD2(I).NE.0))
     *  CALL HADJET(NHAD,AMCSD2(I),POJ,PAT,GACSD2(I),BGXSD2(I), BGYSD2
     +  (I),BGZSD2(I),IFB1,IFB2,IFB3,IFB4, IJCSD2(I),IJCSD2(I),4,NCHSD2
     +  (I),18)
C----------------------------------------------------------------------
        IF((NCHSD1(I).EQ.0))THEN
          ZSEAWU=RNDM(BB)*2.D0*ZSEAAV
          RSEACK=FLOAT(JITT)*PDBSE +ZSEAWU*PDBSEU
          IF(IPCO.GE.1)WRITE(6,'(2A,I5,2F10.3)')'HADJSE JITT,',
     *    'RSEACK,PDBSE ',
     +    JITT,RSEACK,PDBSE
          IREJSS=5
          IF(RNDM(V).LE.RSEACK)THEN
            IREJSS=2
            IF(AMCSD2(I).GT.2.3D0)THEN
              IREJSS=0
              CALL HADJASE(NHAD,AMCSD2(I),POJ,PAT,GACSD2(I),BGXSD2(I),
     *        BGYSD2
     +        (I),BGZSD2(I),IFB1,IFB2,IFB3,IFB4, IJCSD2(I),IJCSD2(I),4,
     *        NCHSD2
     +        (I),3,IREJSS,IISSQQ)
              IF(IPCO.GE.1)WRITE(6,'(2A,I5,F10.3,I5)')'HADJSE JITT,',
     *        'RSEACK,IREJSS ',
     +        JITT,RSEACK,IREJSS
            ENDIF
            IF(IREJSS.GE.1)THEN
              IF(IREJSS.EQ.1)IREJSA=IREJSA+1
              IF(IREJSS.EQ.3)IREJA3=IREJA3+1
              IF(IREJSS.EQ.2)IREJA0=IREJA0+1
        CALL HADJET(NHAD,AMCSD2(I),POJ,PAT,GACSD2(I),BGXSD2(I), BGYSD2
     +  (I),BGZSD2(I),IFB1,IFB2,IFB3,IFB4, IJCSD2(I),IJCSD2(I),4,NCHSD2
     +  (I),18)
              IHADA4=IHADA4+1
            ENDIF
            IF(IREJSS.EQ.0)THEN
              IF(IISSQQ.EQ.3)THEN
                ISEA43=ISEA43+1
              ELSE
                ISEA4=ISEA4+1
              ENDIF
            ENDIF
          ELSE
        CALL HADJET(NHAD,AMCSD2(I),POJ,PAT,GACSD2(I),BGXSD2(I), BGYSD2
     +  (I),BGZSD2(I),IFB1,IFB2,IFB3,IFB4, IJCSD2(I),IJCSD2(I),4,NCHSD2
     +  (I),18)
            IHADA4=IHADA4+1
          ENDIF
        ENDIF
C----------------------------------------------------------------------
C                                   ADD HADRONS/RESONANCES INTO
C                                   COMMON /ALLPAR/ STARTING AT NAUX
        NHKKAU=NHKK+1
        DO 40 J=1,NHAD
          IF (NHKK.EQ.NMXHKK) THEN
            WRITE (6,'(A,2I5/A)') ' HADRZD: NHKK.EQ.NMXHKK ', NHKK,
     +      NMXHKK
            RETURN
          ENDIF
C         NHKK=NHKK+1
         IF (NHKK.EQ.NMXHKK)THEN
           WRITE (6,'(A,2I5)')' XKSAMP:NHKK.EQ.NMXHKK ',NHKK,NMXHKK
           RETURN
         ENDIF
C
          EHECC=SQRT(PXF(J)**2+PYF(J)**2+PZF(J)**2+AMF(J)**2)
        IF (ABS(EHECC-HEF(J)).GT.0.001) THEN
C             WRITE(6,'(2A/3I5,3E15.6)')
C    &            ' HADRZD / CHAIN 2 : CORRECT INCONSISTENT ENERGY ',
C    *            '  NCALSD, NHKK,NREF(J), HEF(J),EHECC, AMF(J)',
C    *            NCALSD, NHKK,NREF(J), HEF(J),EHECC, AMF(J)
          HEF(J)=EHECC
          ENDIF
          ANNZD=ANNZD+1
          EEZD=EEZD+HEF(J)
          PTZD=PTZD+SQRT(PXF(J)**2+PYF(J)**2)
C                               PUT NN-CMS HADRONS INTO /HKKEVT/
          ISTIST=1
          IF(IBARF(J).EQ.500)ISTIST=2
          CALL HKKFIL(ISTIST,MPDGHA(NREF(J)),1,0,
     *                PXF(J),PYF(J),PZF(J),HEF(J),NHKKAU,IORMO(J),24)
          IF(IDHKK(NHKK).EQ.99999) WRITE (6,1030)NHKK,NREF(J), IDHKK
     +    (NHKK)
C         JMOHKK(1,NHKK)=MHKKSS(I)
          IF (IPHKK.GE.2) WRITE(6,1010) NHKK, ISTHKK(NHKK),IDHKK(NHKK),
     +    JMOHKK(1,NHKK),JMOHKK(2,NHKK), JDAHKK(1,NHKK),JDAHKK(2,NHKK),
     +    (PHKK(KHKK,NHKK),KHKK=1,5), (VHKK(KHKK,NHKK),KHKK=1,4)
 
   40   CONTINUE
C       IF(NHAD.GT.0) THEN
C         JDAHKK(1,IMOHKK)=NHKKAU
C         JDAHKK(2,IMOHKK)=NHKK
C       ENDIF
   50 CONTINUE
C----------------------------------------------------------------
C
C     IPHKK=0
      RETURN
 1010 FORMAT (I6,I4,5I6,9E10.2)
 1020 FORMAT (' HADRKK J.GT.NAUMAX SKIP NEXT PARTICLES ',3I10)
 1030 FORMAT (' NHKK,IDHKK(NHKK)  ',3I10)
      END
C ---------------------------------------------------------------
C ---------------------------------------------------------------
C ---------------------------------------------------------------
C
      SUBROUTINE ZOBCMA(IF1,IF2,IF3,IJNCH,NNCH,IREJ,AMCH,AMCHN,IKET)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C  REPLACE SMALL MASS BARYON CHAINS (AMCH)
C  BY OCTETT OR DECUPLETT BARYONS
C
C  MASS CORRECTED FOR NNCH.NE.0
C
C  IREJ=1: CHAIN GENERATION NOT ALLOWED BECAUSE OF TOO SMALL MASS
C          START FROM THE BEGINNING IN HAEVT
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
      COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
C------------------
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEEP,KETMAS.
      COMMON /KETMAS/ AM8(2),AM10(2),IB88(2),IB1010(2),AMCH1N,AMCH2N
*KEND.
C----------------
      CALL DBKLAS(IF1,IF2,IF3,IB8,IBB10)
C
      IF (IPEV.GE.6)WRITE(6,1000)IF1,IF2,IF3,IB8,IBB10
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
      IF(IPEV.GE.6) THEN
        WRITE(6,1010) AMCH,AMCHN,AM81,AM101
        WRITE(6,1020) IF1,IF2,IF3,IB8,IBB10,IJNCH,NNCH,IREJ
 1010 FORMAT(' COBCMA: AMCH,AMCHN,AM81,AM101', 4F13.4)
 1020 FORMAT(' COBCMA: IF1,IF2,IF3,IB8,IBB10,IJNCH,NNCH,IREJ',8I4)
      ENDIF
      RETURN
      END
*CMZ :  1.00/00 27/11/91  17.09.37  by  H.-J. M¶hring
*-- Author :
C
C++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE ZOMCMA(IFQ,IFAQ,IJNCH,NNCH,IREJ,AMCH,AMCHN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C  REPLACE SMALL MASS MESON CHAINS BY PSEUDOSCALAR OR VECTOR MESONS
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
C------------------
*KEEP,INPDAT.
      COMMON /INPDAT/ IMPS(6,6),IMVE(6,6),IB08(6,21),IB10(6,21), IA08
     +(6,21),IA10(6,21), A1,B1,B2,B3,LT,LE,BET,AS,B8,AME,DIQ,ISU
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEND.
      COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
C------------------------
      IIFAQ=IABS(IFAQ)
      IFPS=IMPS(IIFAQ,IFQ)
      IFV=IMVE(IIFAQ,IFQ)
      IF (IPEV.GE.6)WRITE (6,1000)IIFAQ,IFQ,IFPS,IFV
 1000 FORMAT (' COMCMA',5X,' IIPPAQ,ITQ,IFPS,IFV ',4I5)
      AMPS=AAM(IFPS)
      AMV=AAM(IFV)
      NNCH=0
      IJNCH=0
      IREJ=0
      AMFF=AMV+0.3
      IF(IPEV.GE.6) WRITE(6,1010) AMCH,AMPS,AMV,IFPS,IFV
 1010 FORMAT(' AMCH,AMPS,AMV,IFPS,IFV ',3F12.4,2I10)
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
C
      RETURN
      END
*CMZ :  1.00/00 27/11/91  17.09.38  by  H.-J. M¶hring
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C   17/10/89 910191458  MEMBER NAME  MCOMCM2  (KK89.S)      F77
      SUBROUTINE ZOMCM2(IQ1,IQ2,IAQ1,IAQ2,NNCH,IREJ,AMCH)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

C--------------------------------------------------------------------
C  (QQ)-(AQ AQ) CHAIN:
C                      CHECK QUANTUM NUMBERS AND
C                      CORRECT MASS IF NECESSARY
C             REJECT IF THERE IS NO CORRESPONDING PARTICLE
C                    OR TOO LOW MASS
C--------------------------------------------------------------------
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
C------------------
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEEP,INPDAT.
      COMMON /INPDAT/ IMPS(6,6),IMVE(6,6),IB08(6,21),IB10(6,21), IA08
     +(6,21),IA10(6,21), A1,B1,B2,B3,LT,LE,BET,AS,B8,AME,DIQ,ISU
*KEND.
      COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
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
*CMZ :  1.00/00 27/11/91  17.09.38  by  H.-J. M¶hring
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE ZORMOM(AMMM,AMCH1,AMCH1N,AMCH2N, XP,XPP,XTVQ,XTVD,
     +PQ1X,PQ1Y,PQ1Z,PQ1E,PA1X,PA1Y,PA1Z,PA1E, PQ2X,PQ2Y,PQ2Z,PQ2E,PA2X,
     +PA2Y,PA2Z,PA2E, PXCH1,PYCH1,PZCH1,ECH1, PXCH2,PYCH2,PZCH2,ECH2,
     +                                                         IREJ)
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
      COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
C------------------------------------
       IREJ=0
      FAK=AMCH1N/AMCH1
      AMCH1=AMCH1N
      XPSQOL=XP
      XP=XP*FAK
      XPP=XPP + XPSQOL - XP
      XTVDOL=XTVD
      XTVD=XTVD*FAK
      XTVQ=XTVQ + XTVDOL - XTVD
      XPPCM=XPP/(XP+XPP)
      XTVQCM=XTVQ/(XTVQ+XTVD)
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
C---                        NEW MOMENTUM OF CHAIN 1
C                           TO BE COMPENSATED BY CHAIN 2
      PXCH1=PQ1X+PA2X
      PYCH1=PQ1Y+PA2Y
      PZCH1=PQ1Z+PA2Z
      ECH1 =PQ1E+PA2E
      IF(ECH1.LE.AMCH1)THEN
        IREJ=1
C       WRITE(6,'(A)') ' ZORMOM: INCONSISTENT KINEMATICS'
        RETURN
      ENDIF
      PCH1 =SQRT(ABS((ECH1-AMCH1)*(ECH1+AMCH1)))+0.000001
      CXCH1=PXCH1/PCH1
      CYCH1=PYCH1/PCH1
      CZCH1=PZCH1/PCH1
C---                        NEW 4-MOMENTUM OF CHAIN 2
      PXCH2=-PXCH1
      PYCH2=-PYCH1
      PZCH2=-PZCH1
      ECH2 =AMMM - ECH1
      IF(ECH2.LE.PCH1)THEN
        IREJ=1
C       WRITE(6,'(A)') ' ZORMOM: INCONSISTENT KINEMATICS'
        RETURN
      ENDIF
      AMCH2N=SQRT(ABS((ECH2-PCH1)*(ECH2+PCH1)))
C---                        ENERGIES OF PARTONS FROM CHAIN 2
      PA1E=XPPCM*AMMM/2.
      PQ2E=XTVQCM*AMMM/2.
      IF(PCH1.GT.(PA1E+PQ2E)) THEN
        IREJ=1
C       WRITE(6,'(A)') ' ZORMOM: INCONSISTENT KINEMATICS'
        RETURN
      ENDIF
C---                         MOMENTUM COMPONENTS OF PARTONS FROM CHAIN 2
C                            WITH RESPECT TO THE MOMENTUM OF CHAIN 1 (Z)
      CT1=-(PCH1**2 + (PA1E-PQ2E)*(PA1E+PQ2E))/(2.0*PCH1*PA1E)
      if(abs(ct1).gt.1.0) then
C       write(6,'(5x,A/5x,4(1PE15.7))')
C    &        ' ZORMOM: PCH1,PA1E,PQ2E, CT1', PCH1,PA1E,PQ2E, CT1
      CT1=SIGN(0.999999999,CT1)
C       WRITE(6,'(A)') ' ZORMOM: INCONSISTENT KINEMATICS'
	IREJ=1
	RETURN
      endif
      ST1=SQRT(ABS((1.0+CT1)*(1.0-CT1)))
      CALL DSFECF(SFE,CFE)
      CALL DRTRAN(CXCH1,CYCH1,CZCH1,CT1,ST1,SFE,CFE,CXA1,CYA1,CZA1)
      PA1X=CXA1*PA1E
      PA1Y=CYA1*PA1E
      PA1Z=CZA1*PA1E
      PQ2X=PXCH2 - PA1X
      PQ2Y=PYCH2 - PA1Y
      PQ2Z=PZCH2 - PA1Z
C---
      IF(IPRI.GT.1) THEN
        PXSUM=PQ1X+PA1X+PQ2X+PA2X
        PYSUM=PQ1Y+PA1Y+PQ2Y+PA2Y
        PZSUM=PQ1Z+PA1Z+PQ2Z+PA2Z
        PESUM=PQ1E+PA1E+PQ2E+PA2E
        WRITE(6,'(A)') ' ZORMOM: KINEMATIC TEST FOR PARTONS'
        WRITE(6,'(A,1PE12.5)') ' AMMM',AMMM
        WRITE(6,'(A,4(1PE12.5))') ' PXSUM,PYSUM,PZSUM,PESUM', PXSUM,
     +  PYSUM,PZSUM,PESUM
      ENDIF
      RETURN
      END
*CMZ :  1.00/00 27/11/91  17.09.38  by  H.-J. M¶hring
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     DEBUG SUBCHK
C     END DEBUG
      SUBROUTINE ZORVAL(AMMM,IREJ,AMCH1,AMCH2, QTX1,QTY1,QZ1,QE1,QTX2,
     +QTY2,QZ2,QE2,IORI)
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
      COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
      IREJ=0
      IF(AMMM.LE.AMCH1+AMCH2+0.4) THEN
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
C                      ROTATE NEW CHAIN MOMENTA
C                      INTO DIRECTION OF CHAINS BEFORE CORRECTION
      GAM=(QE1+QE2)/AMMM
      BGX=(QTX1+QTX2)/AMMM
      BGY=(QTY1+QTY2)/AMMM
      BGZ=(QZ1+QZ2)/AMMM
C
      IF(ABS(GAM-1.).GT.1E-3) THEN
C       WRITE(6,'(A/5(1PE15.5)/15X,4(1PE15.4),I5)')
C    +  ' ZORVAL: INCONSISTENT KINEMATICS OF CHAINS 
C    +   AMMM, QE1, QTX1, QTY1, QZ1, QE2, QTX2, QTY2, QZ2', AMMM, QE1,
C    +  QTX1, QTY1, QZ1, QE2, QTX2, QTY2, QZ2,IORI
         IREJ=1
          RETURN
      ENDIF
C
      CALL DALTRA(GAM,-BGX,-BGY,-BGZ,PXK1,PYK1,PZK1,EK1,PPPCH1, QTX1,
     +QTY1,QZ1,QE1)
      CALL DALTRA(GAM,-BGX,-BGY,-BGZ,PXK2,PYK2,PZK2,EK2,PPPCH2, QTX2,
     +QTY2,QZ2,QE2)
      IF(IPRI.GT.1) THEN
        WRITE(6,'(2A)') ' ZORVAL - CORRECTION OF CHAIN MOMENTA',
     +  ' IF MASS OF CHAIN 2 HAD TO BE CHANGED'
      ENDIF
      RETURN
      END
