C***********************************************************************
      PROGRAM DPMJET
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
C     DO 7777 KK=1,NHKK
C     WRITE(6,'(2I4,I6,4I4,5F10.2,2I3,I2,I4)')KK,ISTHKK(KK),IDHKK(KK),
C    *JMOHKK(1,KK),JMOHKK(2,KK),JDAHKK(1,KK),JDAHKK(2,KK), 
C    *(PHKK(LL,KK),LL=1,5),IDRES(KK),IDXRES(KK),NOBAM(KK),IDBAM(KK)
C7777 CONTINUE      
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
C     COMMON /NUCCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
C     COMMON /NUCC/   JT,JTZ,JP,JPZ,JJPROJ,JBPROJ,JJTARG,JBTARG
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
      COMMON /NUCCC/   JT,JTZ,JP,JPZ,JJPROJ,JBPROJ,JJTARG,JBTARG
*KEEP,CMHICO.
      COMMON /CMHICO/ CMHIS
*KEEP,RESONA.
      COMMON /RESONA/ IRESO
*KEEP,TRAFOP.
      COMMON /TRAFOP/ GAMP,BGAMP,BETP
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEEP,REJEC.
      COMMON /REJEC/ IRCO1,IRCO2,IRCO3,IRCO4,IRCO5, IRSS11,IRSS12,
     +IRSS13,IRSS14, IRSV11,IRSV12,IRSV13,IRSV14, IRVS11,IRVS12,IRVS13,
     +IRVS14, IRVV11,IRVV12,IRVV13,IRVV14
*KEEP,DROPPT.
      LOGICAL INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA, IPADIS,
     +ISHMAL,LPAULI
      COMMON /DROPPT/ INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA,
     +IPADIS,ISHMAL,LPAULI
*KEEP,DNUN.
      COMMON /DNUN/   NN,NP,NT
*KEEP,NSHMAK.
      COMMON /NSHMAK/ NNSHMA,NPSHMA,NTSHMA,NSHMAC,NSHMA2
*KEEP,DSHM.
      COMMON /DSHM/ RASH,RBSH,BMAX,BSTEP,SIGSH,ROSH,GSH,
     *              BSITE(0:1,200),NSTATB,NSITEB
*KEND.
      COMMON /SEAQXX/ SEAQX,SEAQXN 
      COMMON /DIFFRA/ ISINGD,IDIFTP,IOUDIF,IFLAGD
      COMMON /EVFLAG/ NUMEV
      COMMON /NEUTYY/ NEUTYP,NEUDEC
      COMMON /FLUCTU/IFLUCT
       COMMON /DIQREJ/IDIQRE(7),IDVRE(3),IVDRE(3),IDSRE(3),ISDRE(3),
     *IDZRE(3),IZDRE(3),IDIQRZ(7)
       COMMON /INTNEU/NDZSU,NZDSU
       COMMON /HBOO/IHBOOK
       COMMON /FINAL/IFINAL
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
       COMMON/POPCCK/PDBCK,PDBSE,PDBSEU,
     *  IJPOCK,IREJCK,ICK4,IHAD4,ICK6,IHAD6
     *,IREJSE,ISE4,ISE6,IREJS3,ISE43,ISE63,IREJS0
     *,IHADA4,IHADA6,IREJSA,ISEA4,ISEA6,IREJA3,
     *ISEA43,ISEA63,IREJAO
       COMMON /INXDPM/INTDPM
       COMMON /NSTARI/NSTART
       COMMON /NCSHXX/NCOUXH,NCOUXT
      DIMENSION PP(4)
      COMMON /NUCROS/DSIGSU,DSIGMC,NDSIG
      COMMON /NNCMS/  GAMCM,BGCM,UMO,PCM,EPROJ,PPROJ
      COMMON /DIQSUM/NDVUU,NDVUS,NDVSS,NVDUU,NVDUS,NVDSS,
     *               NDSUU,NDSUS,NDSSS,NSDUU,NSDUS,NSDSS,
     *               NDZUU,NDZUS,NDZSS,NZDUU,NZDUS,NZDSS 
     *              ,NADVUU,NADVUS,NADVSS,NAVDUU,NAVDUS,NAVDSS,
     *               NADSUU,NADSUS,NADSSS,NASDUU,NASDUS,NASDSS,
     *               NADZUU,NADZUS,NADZSS,NAZDUU,NAZDUS,NAZDSS 
      COMMON /HDJASE/NHSE1,NHSE2,NHSE3,NHASE1,NHASE2,NHASE3
      COMMON /CASADI/CASAXX,ICASAD
      COMMON /VXSVD/VXSP(50),VXST(50),VXSAP(50),VXSAT(50),
     *              VXVP(50),VXVT(50),VXDP(50),VXDT(50),
     *	    NXSP,NXST,NXSAP,NXSAT,NXVP,NXVT,NXDP,NXDT
      DIMENSION VXSSS(50,6),VXVVV(50,6),XXXX(50,6)
      DIMENSION XB(200),BIMPP(200)
      PARAMETER (INTMX=2488,INTMD=252)
      COMMON/SHMAKL/JSSH(INTMX),JTSH(INTMX),INTER1(INTMX),INTER2(INTMX)
      COMMON /INFORE/IFREJ         
      CHARACTER*80 TITLED
      CHARACTER*8 PROJTY,TARGTY
      COMMON /USER1/TITLED,PROJTY,TARGTY
      COMMON /USER2/CMENER,SDFRAC,PTLAR,ISTRUF,ISINGX,IDUBLD
      COMMON /STRUFU/ISTRUM,ISTRUT
      COMMON /PTSAMP/ ISAMPT
      COMMON/COLLIS/S,IJPROX,IJTAR,PTTHR,PTTHR2,IOPHRD,IJPRLU,IJTALU
      COMMON /DROPJJ/DROPJT,DROPVA
      COMMON /POMTYP/IPIM,ICON,ISIG,LMAX,MMAX,NMAX,DIFEL,DIFNU
      COMMON/PSHOW/IPSHOW
      COMMON /ZENTRA/ ICENTR
      COMMON /EVAPPP/IEVAP
      COMMON /SEASU3/SEASQ
      COMMON /RECOM/IRECOM
      COMMON /TAUFO/  TAUFOR,KTAUGE,ITAUVE,INCMOD
      COMMON/POPCOR/PDB,AJSDEF
      COMMON /DIQUAX/AMEDD,IDIQUA,IDIQUU
      COMMON /COLLE/NEVHAD,NVERS,IHADRZ,NFILE
      COMMON /NUCIMP/ PRMOM(5,248),TAMOM(5,248), PRMFEP,PRMFEN,TAMFEP,
     +TAMFEN, PREFEP,PREFEN,TAEFEP,TAEFEN, PREPOT(210),TAEPOT(210),
     +PREBIN,TAEBIN,FERMOD,ETACOU
       COMMON /CRONIN/CRONCO,MKCRON
      COMMON /XSEADI/ XSEACU,UNON,UNOM,UNOSEA, CVQ,CDQ,CSEA,SSMIMA,
     +SSMIMQ,VVMTHR
      COMMON /SECINT/ISECIN
      COMMON /NDON/NDONE
C---------------------
*
      DATA NCOUNT/0/
C                               from DTUJET93
C                               from DTUJET93
C     ON DOUBLE PRECISION UNDERFLOW IGNORE
C     ON DOUBLE PRECISION OVERFLOW IGNORE
C     ON DOUBLE PRECISION INEXACT IGNORE
C     ON DOUBLE PRECISION ILLEGAL IGNORE
C     ON DOUBLE PRECISION DIV 0 IGNORE
C     ON REAL UNDERFLOW IGNORE
C     ON REAL OVERFLOW IGNORE
C     ON INTEGER OVERFLOW IGNORE
C     ON REAL INEXACT IGNORE
C     ON REAL ILLEGAL IGNORE
C     ON REAL DIV 0 IGNORE
C                               from DTUJET93
*---extended error handling on RISC 6000
C     include 'fexcp.h'
C     call signal(SIGTRAP,xl_trce)
C***********************************************************************a
C     OPEN(47,FILE='/u1/ranft/dtunuc44/GLAUBTAR.DAT',
C    *STATUS='UNKNOWN')
C     OPEN(47,FILE='/nfs/hptrack/user/ran/dtunuc44/GLAUBTAR.DAT',
C     OPEN(47,FILE='/user/ran/dtunuc44/GLAUBTAR.DAT',
C    *STATUS='UNKNOWN')
C     OPEN(47,FILE='/lapphp11_2/users/ranft/dtunuc44/GLAUBTAR.DAT',
C    *STATUS='UNKNOWN')
      OPEN(47,FILE='GLAUBTAR.DAT',
     *STATUS='UNKNOWN')
      OPEN(37,FILE='GLAUBCROSSPB.DAT',
     *STATUS='UNKNOWN')
C     OPEN( 2,FILE='HIBLD.DAT',STATUS='OLD')
      AAM(5)=0.001D0
      AAM(6)=0.001D0
      AAM(133)=0.001D0
      AAM(134)=0.001D0
      AAM(135)=0.001D0
      AAM(136)=0.001D0
C     Initialize x-distribtion survey
      NXSP=0
      NXST=0
      NXSAP=0
      NXSAT=0
      NXVP=0
      NXVT=0
      NXDP=0
      NXDT=0
      AXSP=0.
      AXST=0.
      AXSAP=0.
      AXSAT=0.
      AXVP=0.
      AXVT=0.
      AXDP=0.
      AXDT=0.
      DO 5271 II=1,50
        VXSP(II)=1.D-8
        VXST(II)=1.D-8
        VXSAP(II)=1.D-8
        VXSAT(II)=1.D-8
        VXVP(II)=1.D-8
        VXVT(II)=1.D-8
        VXDP(II)=1.D-8
        VXST(II)=1.D-8
 5271 CONTINUE      
C
* random number initialization for LEPTO
      CALL RLUXGO(LUX_LEVEL,ISEED,0,0)
      IDIQRE(1)=0
      IDIQRE(2)=0
      IDIQRE(3)=0
      IDIQRE(4)=0
      IDIQRE(5)=0
      IDIQRE(6)=0
      IDIQRE(7)=0
      IDIQRZ(1)=0
      IDIQRZ(2)=0
      IDIQRZ(3)=0
      IDIQRZ(4)=0
      IDIQRZ(5)=0
      IDIQRZ(6)=0
      IDIQRZ(7)=0
      IDVRE(1)=0
      IDVRE(2)=0
      IDVRE(3)=0
      IVDRE(1)=0
      IVDRE(2)=0
      IVDRE(3)=0
      IDSRE(1)=0
      IDSRE(2)=0
      IDSRE(3)=0
      ISDRE(1)=0
      ISDRE(2)=0
      ISDRE(3)=0
      IDZRE(1)=0
      IDZRE(2)=0
      IDZRE(3)=0
      IZDRE(1)=0
      IZDRE(2)=0
      IZDRE(3)=0
      NDVUU=0
      NDVUS=0
      NDVSS=0
      NVDUU=0
      NVDUS=0
      NVDSS=0
      NDSUU=0
      NDSUS=0
      NDSSS=0
      NSDUU=0
      NSDUS=0
      NSDSS=0
      NDZUU=0
      NDZUS=0
      NDZSS=0
      NZDUU=0
      NZDUS=0
      NZDSS=0
      NADVUU=0
      NADVUS=0
      NADVSS=0
      NAVDUU=0
      NAVDUS=0
      NAVDSS=0
      NADSUU=0
      NADSUS=0
      NADSSS=0
      NASDUU=0
      NASDUS=0
      NASDSS=0
      NADZUU=0
      NADZUS=0
      NADZSS=0
      NAZDUU=0
      NAZDUS=0
      NAZDSS=0
      NHSE1=0
      NHSE2=0
      NHSE3=0
      NHASE1=0
      NHASE2=0
      NHASE3=0
 1000 CONTINUE
      NCOUNT=NCOUNT+1
*      initialisation routine:
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                                     Parton pt distribution
      CALL PARPT(1,PT1,PT2,IPT,NEVT)
C-------------------------------------------------------------
C     CALL DMINIT(NCASES,MULTE,EPN,PPN,NCOUNT,IGLAUB)
      CALL DMINIT(NCASES,MULTE,EPNN,PPNN,NCOUNT,IGLAUB)
      EPN=EPNN
      PPN=PPNN
C----------------------------------------------------------------
C----------------------------------------------------------------
C---now starts the real work
      NHKKH1=1
      ISHMAL=.TRUE.
C
        TTOT=0.
        TMAX=0.
C     CALL TIMEL (TCPU)
C     TLIM=MAX(TCPU/100.,15.)
C     CALL TIMED(TDIFF)
      IIT=IT
      IITZ=ITZ
      IIP=IP
      IIPZ=IPZ
      IIPROJ=IJPROJ
      IITARG=IJTARG
      IF( IGLAUB.EQ.1) THEN
        KKMAT=0
      ELSE
        KKMAT=1
      ENDIF
C===============================================================      
C     Printout of important Parameters (defaults and input cards)       
C===============================================================      
      WRITE(6,*)' Printout of important Parameters before DPMJET run.'
     *,' Please note for DPMJET input all numbers are floating point!' 
      WRITE(6,*)'PROJPAR  ',IP,IPZ
      WRITE(6,*)'TARPAR   ',IT,ITZ
      WRITE(6,*)'MOMENTUM ',PPN
      WRITE(6,*)'ENERGY   ',EPN
      WRITE(6,*)'CMENERGY ',UMO
      WRITE(6,*)'NOFINALE ',IFINAL
      WRITE(6,*)'EVAPORAT ',IEVAP
      WRITE(6,*)'OUTLEVEL ',IPRI,IPEV,IPPA,IPCO,INIT,IPHKK
      AUAUAU=RD2OUT(ISEED1,ISEED2)
      WRITE(6,*)'RANDOMIZ ',ISEED1,ISEED2, ' Initial RNDM (RM48) seeds'
      WRITE(6,*)'STRUCFUN ',ISTRUF+100*ISTRUT
      WRITE(6,*)'SAMPT    ',ISAMPT
      WRITE(6,*)'SELHARD  ',0,IOPHRD, 0,DROPJT,PTTHR,PTTHR2 
      WRITE(6,*)'SIGMAPOM ',0,ISIG,IPIM+10*ICON,IMAX,MMAX,NMAX
      WRITE(6,*)'PSHOWER  ',IPSHOW
      WRITE(6,*)'CENTRAL  ',ICENTR
      WRITE(6,*)'CMHISTO  ',CMHIS
      WRITE(6,*)'SEASU3   ',SEASQ
      WRITE(6,*)'RECOMBIN ',IRECOM
      WRITE(6,*)'SINGDIFF ',ISINGD
      WRITE(6,*)'TAUFOR   ',TAUFOR,KTAUGE,ITAUVE
      WRITE(6,*)'POPCORN  ',PDB
      WRITE(6,*)'POPCORCK ',IJPOCK,PDBCK
      WRITE(6,*)'POPCORSE ',PDBSE,PDBSEU
      WRITE(6,*)'CASADIQU ',ICASAD,CASAXX
      WRITE(6,*)'DIQUARKS ',IDIQUA,IDIQUU,AMEDD
      WRITE(6,*)'HADRONIZ ',IHADRZ
      WRITE(6,*)'INTPT    ',INTPT
      WRITE(6,*)'PAULI    ',LPAULI
      WRITE(6,*)'FERMI    ',FERMP,FERMOD
      WRITE(6,*)'CRONINPT ',MKCRON,CRONCO
      WRITE(6,*)'SEADISTR ',XSEACU+0.95D0,UNON,UNOM,UNOSEA
      WRITE(6,*)'SEAQUARK ',SEAQX,SEAQXN
      WRITE(6,*)'SECINTER ',ISECIN
      WRITE(6,*)'XCUTS    ',CVQ,CDQ,CSEA,SSMIMA
      WRITE(6,*)'START    ',NCASES
      WRITE(6,*)' Printout of important Parameters before DPMJET run.'
     *,' Please note for DPMJET input all numbers are floating point!' 
C===============================================================      
C     Printout of important Parameters (defaults and input cards)       
C===============================================================      
      NCASET=NCASES/10
      DO 181 IIII=1,10
        NDONE=(IIII-1)*NCASET
        WRITE(6,1111)NDONE
        WRITE(6,'(A,4I5)')' KKINC call ',IIT,IITZ,IIP,IIPZ
        CALL TIMDAT
 1111   FORMAT(' NDONE= ',I10)
        DO 180 I=1,NCASET
	  NDONE=NDONE+1
C       WRITE(6,1111)NDONE
	  IF(MULTE.EQ.0)THEN
            EPN=EPNN
            PPN=PPNN
	  ELSEIF(MULTE.EQ.1)THEN
	    EPN=0.1D0*EPNN+RNDM(V)*(1.9D0*EPNN)
            NNPP=1
            IF(IJPROJ.NE.0) NNPP=IJPROJ
            PPN=SQRT((EPN-AAM(NNPP))*(EPN+AAM(NNPP)))

*                             nucleon-nucleon cms
C           IBPROJ=1
            EPROJ=EPN
            AMPROJ=AAM(NNPP)
            AMTAR=AAM(1)
            PPROJ = SQRT((EPN-AMPROJ)*(EPN+AMPROJ))
            UMO = SQRT(AMPROJ**2 + AMTAR**2 + 2.*AMTAR*EPROJ)
            CMENER=UMO
            IF(ISTRUT.EQ.1)THEN
              PTTHR=2.1+0.15*(LOG10(CMENER/50.))**3
              PTTHR2=PTTHR
            ELSEIF(ISTRUT.EQ.2)THEN
              PTTHR=2.5+0.12*(LOG10(CMENER/50.))**3
              PTTHR2=PTTHR
            ENDIF
            GAMCM = (EPROJ+AMTAR)/UMO
            BGCM=PPROJ/UMO
            ECM=UMO
            PCM=GAMCM*PPROJ - BGCM*EPROJ
C
C           PRINT 1033, EPROJ,PPROJ,
C    +      AMPROJ,AMTAR,UMO,GAMCM,BGCM,PCM
 1033       FORMAT(' CMS: ' , 
     +'     EPROJ,PPROJ,AMPROJ,AMTAR,UMO,GAMCM,BGCM,PCM'/8E22.13)
	  ENDIF
  765     CONTINUE
          NUMEV = I+(IIII-1)*NCASET
          IF ((I.EQ.486).OR.(I.EQ.803).OR.(I.EQ.1368).OR.
     &    (I.EQ.1465).OR.(I.EQ.1693).OR.(I.EQ.1808)) THEN
C           IPEV  = 7
C           IPCO  = 7
C           IPHKK = 7
          ENDIF
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
C         COMMON /NCOUCH/ ACOUVV,ACOUSS,ACOUSV,ACOUVS,
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
C         IIP=IP+1-IIII
C         WRITE(6,'(A,4I5)')' KKINC call ',IIT,IITZ,IIP,IIPZ
C         CALL KKINC(EPN,IIT,IITZ,IIP,IIPZ,IIPROJ,KKMAT,
C    *    IITARG,NHKKH1,IREJ)
      IF(INTDPM.EQ.0)THEN
        CALL KKINC(EPN,IIT,IITZ,IIP,IIPZ,IIPROJ,KKMAT,
     *   IITARG,NHKKH1,IREJ)
      ELSEIF(INTDPM.EQ.1)THEN
C       ELABLO=1.D0+7.D0*RNDM(V)
CELABT=10.D0**ELABLO
        ELABT=EPN/1000.D0
CIIIPRO=1
        IIIPRO=IIPROJ
	IIIP=IIP
	IIIPZ=IIPZ
	IIIT=IIT
	IIITZ=IITZ
CIIIP=1
CIIIPZ=1
CIF(RNDM(VV).LT.0.3D0 )THEN
C  IIIPRO=13
CELSEIF(RNDM(VVV).GT.0.7D0)THEN
C  IIIPRO=2
CENDIF
C       CALL DPMEVT(ELABT,IIIPRO,IIIP,IIIPZ,NHKKH1)
        CALL DPMEVT(ELABT,IIIPRO,IIIP,IIIPZ,IIIT,IIITZ,KKMAT,NHKKH1)
      ENDIF
C         CALL TIMDAT
          IF(IREJ.EQ.1)GO TO 765
C         IPEV  = 0
C         IPCO  = 0
C         IPHKK = 0
C                                   INITIALIZE COUNTERS
          BNNVV=BNNVV+ANNVV
          BNNSS=BNNSS+ANNSS
          BNNSV=BNNSV+ANNSV
          BNNVS=BNNVS+ANNVS
          BNNCC=BNNCC+ANNCC
          BNNDV=BNNDV+ANNDV
          BNNVD=BNNVD+ANNVD
          BNNDS=BNNDS+ANNDS
          BNNSD=BNNSD+ANNSD
          BNNHH=BNNHH+ANNHH
          BNNZZ=BNNZZ+ANNZZ
          BNNDI=BNNDI+ANNDI
          BNNZD=BNNZD+ANNZD
          BNNDZ=BNNDZ+ANNDZ
          BPTVV=BPTVV+PTVV
          BPTSS=BPTSS+PTSS
          BPTSV=BPTSV+PTSV
          BPTVS=BPTVS+PTVS
          BPTCC=BPTCC+PTCC
          BPTDV=BPTDV+PTDV
          BPTVD=BPTVD+PTVD
          BPTDS=BPTDS+PTDS
          BPTSD=BPTSD+PTSD
          BPTHH=BPTHH+PTHH
          BPTZZ=BPTZZ+PTZZ
          BPTDI=BPTDI+PTDI
          BPTZD=BPTDZ+PTDZ
          BPTDZ=BPTDZ+PTDZ
          BEEVV=BEEVV+EEVV
          BEESS=BEESS+EESS
          BEESV=BEESV+EESV
          BEEVS=BEEVS+EEVS
          BEECC=BEECC+EECC
          BEEDV=BEEDV+EEDV
          BEEVD=BEEVD+EEVD
          BEEDS=BEEDS+EEDS
          BEESD=BEESD+EESD
          BEEHH=BEEHH+EEHH
          BEEZZ=BEEZZ+EEZZ
          BEEDI=BEEDI+EEDI
          BEEZD=BEEZD+EEZD
          BEEDZ=BEEDZ+EEDZ
C         COMMON /NCOUCH/ ACOUVV,ACOUSS,ACOUSV,ACOUVS,
C    *                 ACOUZZ,ACOUHH,ACOUDS,ACOUSD,
C    *                 ACOUDZ,ACOUZD,ACOUDI
          BCOUVV=BCOUVV+ACOUVV
          BCOUSS=BCOUSS+ACOUSS
          BCOUSV=BCOUSV+ACOUSV
          BCOUVS=BCOUVS+ACOUVS
          BCOUZZ=BCOUZZ+ACOUZZ
          BCOUHH=BCOUHH+ACOUHH
          BCOUDS=BCOUDS+ACOUDS
          BCOUSD=BCOUSD+ACOUSD
          BCOUDZ=BCOUDZ+ACOUDZ
          BCOUZD=BCOUZD+ACOUZD
          BCOUDI=BCOUDI+ACOUDI
          BCOUDV=BCOUDV+ACOUDV
          BCOUVD=BCOUVD+ACOUVD
          BCOUCC=BCOUCC+ACOUCC
*
C         HOW LONG DID IT TAKE TO PROCESS THIS ONE?
C
C         CALL TIMED(TDIFF)
C         IF(TDIFF.GT.TMAX)TMAX=TDIFF
C         TTOT=TTOT+TDIFF
C         TMEAN=TTOT/FLOAT(I)
C
C         CONDITIONS FOR LOOP TERMINATION
C
C         CALL TIMEL(TLEFT)
C         IF ( TLEFT .LE. 3.*TMAX+TLIM )                       GO TO 190
  180   CONTINUE
C
C         WRITE(6,'(A,4I5)')' KKINC call ',IIT,IITZ,IIP,IIPZ
C     DO 7777 KK=1,NHKK
C     WRITE(6,'(2I4,I6,4I4,5F10.2,2I3,I2,I4)')KK,ISTHKK(KK),IDHKK(KK),
C    *JMOHKK(1,KK),JMOHKK(2,KK),JDAHKK(1,KK),JDAHKK(2,KK), 
C    *(PHKK(LL,KK),LL=1,5),IDRES(KK),IDXRES(KK),NOBAM(KK),IDBAM(KK)
C7777 CONTINUE      
C
  181 CONTINUE
C===============================================================      
C     Printout of important Parameters (defaults and input cards)       
C===============================================================      
      WRITE(6,*)' Printout of important Parameters after DPMJET run.'
     *,' Please note for DPMJET input all numbers are floating point!' 
      WRITE(6,*)'PROJPAR  ',IP,IPZ
      WRITE(6,*)'TARPAR   ',IT,ITZ
      WRITE(6,*)'MOMENTUM ',PPN
      WRITE(6,*)'ENERGY   ',EPN
      WRITE(6,*)'CMENERGY ',UMO
      WRITE(6,*)'NOFINALE ',IFINAL
      WRITE(6,*)'EVAPORAT ',IEVAP
      WRITE(6,*)'OUTLEVEL ',IPRI,IPEV,IPPA,IPCO,INIT,IPHKK
      AUAUAU=RD2OUT(ISEED1,ISEED2)
      WRITE(6,*)'RANDOMIZ ',ISEED1,ISEED2, ' Final RNDM (RM48) seeds'
      WRITE(6,*)'STRUCFUN ',ISTRUF+100*ISTRUT
      WRITE(6,*)'SAMPT    ',ISAMPT
      WRITE(6,*)'SELHARD  ',0,IOPHRD, 0,DROPJT,PTTHR,PTTHR2 
      WRITE(6,*)'SIGMAPOM ',0,ISIG,IPIM+10*ICON,IMAX,MMAX,NMAX
      WRITE(6,*)'PSHOWER  ',IPSHOW
      WRITE(6,*)'CENTRAL  ',ICENTR
      WRITE(6,*)'CMHISTO  ',CMHIS
      WRITE(6,*)'SEASU3   ',SEASQ
      WRITE(6,*)'RECOMBIN ',IRECOM
      WRITE(6,*)'SINGDIFF ',ISINGD
      WRITE(6,*)'TAUFOR   ',TAUFOR,KTAUGE,ITAUVE
      WRITE(6,*)'POPCORN  ',PDB
      WRITE(6,*)'POPCORCK ',IJPOCK,PDBCK
      WRITE(6,*)'POPCORSE ',PDBSE,PDBSEU
      WRITE(6,*)'CASADIQU ',ICASAD,CASAXX
      WRITE(6,*)'DIQUARKS ',IDIQUA,IDIQUU,AMEDD
      WRITE(6,*)'HADRONIZ ',IHADRZ
      WRITE(6,*)'INTPT    ',INTPT
      WRITE(6,*)'PAULI    ',LPAULI
      WRITE(6,*)'FERMI    ',FERMP,FERMOD
      WRITE(6,*)'CRONINPT ',MKCRON,CRONCO
      WRITE(6,*)'SEADISTR ',XSEACU+0.95D0,UNON,UNOM,UNOSEA
      WRITE(6,*)'SEAQUARK ',SEAQX,SEAQXN
      WRITE(6,*)'SECINTER ',ISECIN
      WRITE(6,*)'XCUTS    ',CVQ,CDQ,CSEA,SSMIMA
      WRITE(6,*)' Printout of important Parameters after DPMJET run.'
     *,' Please note for DPMJET input all numbers are floating point!' 
C===============================================================      
C     Printout of important Parameters (defaults and input cards)       
C===============================================================      

C                        OUTPUT of RNDM seeds
      AUAUAU=RD2OUT(ISEED1,ISEED2)
      WRITE (6,*)' Final RNDM seeds (RM48) ',ISEED1,ISEED2
      WRITE (6,*)' Final RNDM seeds (RM48) ',ISEED1,ISEED2
      WRITE (6,*)' Final RNDM seeds (RM48) ',ISEED1,ISEED2
      WRITE (6,*)' Final RNDM seeds (RM48) ',ISEED1,ISEED2
                                                               GO TO 200
  190 CONTINUE
C     WRITE(6,*)' STOPPED FOR CPUTIME LIMIT: ',I,' EVENTS ',
C    +'INSTEAD OF ',NCASES,' PRODUCED'
C     NCASES = I
  200 CONTINUE
C     WRITE (6,1090)TTOT,TMAX,TCPU,TLIM,TDIF,TMEAN,TLEFT
C1090 FORMAT (' TTOT,TMAX,TCPU,TLIM,TDIF,TMEAN.TLEFT '/7F10.2)
C
      IF(IPEV.GE.-1) THEN
      IF(IFREJ.EQ.1)THEN
        WRITE(6,1100) IRVV11,IRVV12,IRVV13,IRVV14, IRSV11,IRSV12,IRSV13,
     +  IRSV14, IRVS11,IRVS12,IRVS13,IRVS14, IRSS11,IRSS12,IRSS13,IRSS14
 1100 FORMAT (' REJECTION COUNTERS FROM KKEVT',/, 5X,' V-V CHAINS',4I6/
     +5X,' S-V CHAINS',4I6/ 5X,' V-S CHAINS',4I6/ 5X,' S-S CHAINS',4I6)
	WRITE(6,'(A,4I10)')' POPCCK/SE/S3/S0 rejections ',
     *	IREJCK,IREJSE,IREJS3,IREJS0
	WRITE(6,'(A,4I10)')' POPCCK/ASE/AS3/AS0 rejections ',
     *	IREJSA,IREJA3,IREJA0
	WRITE(6,'(2A,8I6)')' POPCCK ICK4,ICK6,IHAD4,IHAD6,ISE4,ISE6 ',
     *  'ISE43,ISE63 ', ICK4,ICK6,IHAD4,IHAD6,ISE4,ISE6,ISE43,ISE63
	WRITE(6,'(2A,8I6)')' POPCSAQ IHADA4,IHADA6,ISEA4,ISEA6 ',
     *  'ISEA43,ISEA63 ', IHADA4,IHADA6,ISEA4,ISEA6,ISEA43,ISEA63
      WRITE(6,*)   ' NDVUU,NDVUS,NDVSS,NVDUU,NVDUS,NVDSS',
     *             ' NDSUU,NDSUS,NDSSS,NSDUU,NSDUS,NSDSS',
     *             ' NDZUU,NDZUS,NDZSS,NZDUU,NZDUS,NZDSS' ,
     *               NDVUU,NDVUS,NDVSS,NVDUU,NVDUS,NVDSS,
     *               NDSUU,NDSUS,NDSSS,NSDUU,NSDUS,NSDSS,
     *               NDZUU,NDZUS,NDZSS,NZDUU,NZDUS,NZDSS  
      WRITE(6,*)   ' NADVUU,NADVUS,NADVSS,NAVDUU,NAVDUS,NAVDSS',
     *             ' NADSUU,NADSUS,NADSSS,NASDUU,NASDUS,NASDSS',
     *             ' NADZUU,NADZUS,NADZSS,NAZDUU,NAZDUS,NAZDSS' ,
     *               NADVUU,NADVUS,NADVSS,NAVDUU,NAVDUS,NAVDSS,
     *               NADSUU,NADSUS,NADSSS,NASDUU,NASDUS,NASDSS,
     *               NADZUU,NADZUS,NADZSS,NAZDUU,NAZDUS,NAZDSS  
      WRITE(6,*)' NHSE1,NHSE2,NHSE3,NHASE1,NHASE2,NHASE3 ',
     * NHSE1,NHSE2,NHSE3,NHASE1,NHASE2,NHASE3 
      ENDIF
      ENDIF
C
      CALL TIMDAT
C****************** PRINTOUT********************************
      IF(IFREJ.EQ.1)THEN
      WRITE(6,'(A/7I8)')' Diquark rejection IDIQRE(1-7),N,ss,su,ud',
     &                                   (IDIQRE(JJ),JJ=1,7)
      WRITE(6,'(A/7I8)')' Diquark rejection IDIQRZ(1-7),N,ss,su,ud',
     &                                   (IDIQRZ(JJ),JJ=1,7)
      WRITE(6,*)' Diquark rej. IDVRE(1-3),ud,us,ss ',(IDVRE(JJ),JJ=1,3)
      WRITE(6,*)' Diquark rej. IVDRE(1-3),ud,us,ss ',(IVDRE(JJ),JJ=1,3)
      WRITE(6,*)' Diquark rej. IDSRE(1-3),ud,us,ss ',(IDSRE(JJ),JJ=1,3)
      WRITE(6,*)' Diquark rej. ISDRE(1-3),ud,us,ss ',(ISDRE(JJ),JJ=1,3)
      WRITE(6,*)' Diquark rej. IDZRE(1-3),ud,us,ss ',(IDZRE(JJ),JJ=1,3)
      WRITE(6,*)' Diquark rej. IZDRE(1-3),ud,us,ss ',(IZDRE(JJ),JJ=1,3)
      WRITE(6,*)' NDZSU,NZDSU ',NDZSU,NZDSU
      ENDIF
      IF ((CMHIS.EQ.1.D0).AND.(IOUDIF.EQ.1))
     &    CALL DIADIF(3,NHKKH1)
C     Output of x-distribution survey
      IF(IFREJ.EQ.1)THEN
      WRITE(6,*)' Output of x-distribution survey',
     * 'VXSP(II),VXST(II),VXSAP(II),VXSAT(II),',
     *'VXVP(II),VXVT(II),VXDP(II),VXDT(II)' ,
     *NXSP,NXST,NXSAP,NXSAT,NXVP,NXVT,NXDP,NXDT
      DO 6671 II=1,50
        IF(NXSP.GE.1)VXSP(II)=50.D0*VXSP(II)/NXSP
        IF(NXST.GE.1)VXST(II)=50.D0*VXST(II)/NXST
        IF(NXSAP.GE.1)VXSAP(II)=50.D0*VXSAP(II)/NXSAP
        IF(NXSAT.GE.1)VXSAT(II)=50.D0*VXSAT(II)/NXSAT
        IF(NXVP.GE.1)VXVP(II)=50.D0*VXVP(II)/NXVP
        IF(NXVT.GE.1)VXVT(II)=50.D0*VXVT(II)/NXVT
        IF(NXDP.GE.1)VXDP(II)=50.D0*VXDP(II)/NXDP
        IF(NXDT.GE.1)VXDT(II)=50.D0*VXDT(II)/NXDT
	XXXXX=II*0.02D0-0.01D0
	XXXX(II,1)=XXXXX
	XXXX(II,2)=XXXXX
	XXXX(II,3)=XXXXX
	XXXX(II,4)=XXXXX
	XXXX(II,5)=XXXXX
	XXXX(II,6)=XXXXX
	FXVVV=(1.-XXXXX)**3/SQRT(XXXXX)
	FXDDD=2.D0*XXXXX**3.0D0/SQRT(1.D0-XXXXX)
	VXSSS(II,1)=LOG10(VXSP(II))
	VXSSS(II,2)=LOG10(VXST(II))
	VXSSS(II,3)=LOG10(VXSAP(II))
	VXSSS(II,4)=LOG10(VXSAT(II))
	VXVVV(II,1)=LOG10(VXVP(II))
	VXVVV(II,2)=LOG10(VXVT(II))
	VXVVV(II,3)=LOG10(VXDP(II))
	VXVVV(II,4)=LOG10(VXDT(II))
	VXVVV(II,5)=LOG10(FXVVV)
	VXVVV(II,6)=LOG10(FXDDD)
	AXSP=AXSP+0.02D0*VXSP(II)*XXXXX
	AXST=AXST+0.02D0*VXST(II)*XXXXX
	AXSAP=AXSAP+0.02D0*VXSAP(II)*XXXXX
	AXSAT=AXSAT+0.02D0*VXSAT(II)*XXXXX
	AXVP=AXVP+0.02D0*VXVP(II)*XXXXX
	AXVT=AXVT+0.02D0*VXVT(II)*XXXXX
	AXDP=AXDP+0.02D0*VXDP(II)*XXXXX
	AXDT=AXDT+0.02D0*VXDT(II)*XXXXX
	WRITE(6,*)VXSP(II),VXST(II),VXSAP(II),VXSAT(II),
     *  	  VXVP(II),VXVT(II),VXDP(II),VXDT(II)
 6671 CONTINUE      
      WRITE(6,*)
     *AXSP,AXST,AXSAP,AXSAT,AXVP,AXVT,AXDP,AXDT
      CALL PLOT(XXXX,VXSSS,200,4,50,0.D0,0.02D0,-3.D0,0.05D0)
      CALL PLOT(XXXX,VXVVV,300,6,50,0.D0,0.02D0,-3.D0,0.05D0)
      ENDIF
*
      IF(IPADIS) CALL DISTPA(3)
      CALL PARPT(3,PT1,PT2,IPT,NCASES)
      FRACXS=0.D0
      IF(NSTART.EQ.1)THEN
      IF(NCOUXH.GE.0)THEN
        FRACXS=FLOAT(NCOUXH)/(FLOAT(NCOUXH)+FLOAT(NCOUXT))
      ENDIF
      WRITE(6,*)' Fraction of x-sect: ',FRACXS,NCOUXH,NCOUXT
      ENDIF
      IF(NSTART.EQ.2)THEN
C                         print neutrino-nucleon cross section
	DSIGMC=0.
	IF(NDSIG.GE.1) DSIGMC=DSIGSU/NDSIG
	WRITE(6,*)' Neutrino-nucleon cross section DSIGMC,NDSIG ',
     &  DSIGMC,' *10**(-38) cm**2 ',NDSIG,' evts'
      ENDIF
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
C     IF(IFREJ.EQ.1)THEN
      IF(ISHMAL) CALL SHMAK(3,NSHMAC,NP,NT,IP,IT,UMO,BIMP)
      IF(ISHMAL) CALL SHMAK1(3,NSHMA2,NP,NT,IP,IT,UMO,BIMP)
C     ENDIF
      IF (IRESO.EQ.1)  CALL DISTRP(3,NCASES,PPN)
      IF (CMHIS.EQ.0.D0) CALL DISTR(3,NCASES,PPN,IDUMMY)
      IF (CMHIS.EQ.1.D0) CALL DISTRC(3,NCASES,PPN,IDUMMY)
      IF (CMHIS.EQ.2.D0) CALL DISTCO(3,NCASES,PPN,IDUMMY)
      IF (IRESO.EQ.1) CALL DISRES(3,NCASES,PPN)
C                                         HBOOK HISTOGRAMS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       IF(IHBOOK.EQ.1.AND.CMHIS.EQ.0.D0)THEN
         CALL PLOMB(5,PP,CHAR,XFXFXF,ITIF,IJPROJ)
       ENDIF
       IF(IHBOOK.EQ.1.AND.CMHIS.EQ.1.D0)THEN
         CALL PLOMBC(5,PP,CHAR,XFXFXF,ITIF,IJPROJ)
       ENDIF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL TIMDAT
C
C----------------------------------------------------------------
      GO TO 1000
      END
C
C*****************************************************************
C
      SUBROUTINE DMINIT(NCASES,MULTE,EPN,PPN,NCOUNT,IGLAUB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C     UNITARIZED HARD AND SOFT MULTICHAIN FRAGMENTATION MODEL
C
C                    CODE DPMJET-II.5
C
C********************************************************************
C
C
C     AUTHORS:
C               J. RANFT, (Johannes.Ranft@cern.ch) 
C
C
C
C     THE CODE
C                      D P M J E T
C
C IS  A CODE TO CALCULATE PARTICLE PRODUCTION
C
C   IN HADRON-NUCLEUS AND NUCLEUS-NUCLEUS COLLISIONS
C
C   USING THE DUAL PARTON MODEL AND THE FORMATION ZONE
C
C                                   INTRANUCLEAR CASCADE
C
C                          Minijets  and soft multijets (DTUJET like)
C
C
C     D P M J E T  PYTHIA-6.1 
C
C     REVISION HISTORY: OCTOBER 1989: FORMATION ZONE CASCADE
C                                     IMPLEMENTED IN TARGET NUCLEUS
C                                     ONLY. THIS IS SUFFICIENT FOR
C                                     HADRON-NUCLEUS COLLISIONS
C
C                       JANUARY 1990: FORMATION ZONE CASCADE IN PROJEC.
C                                     PAULI PRINCIPLE
C                                     ELASTIC COLLISIONS IN FORMATION
C                                     ZONE CASCADE
C
C                                     DUAL PARTON MODEL FOR H-A AND A-A
C
C
C                       JULY    1991: PROPERTIES OF RESIDUAL NUCLEI
C                                     WORKED OUT
C
C                        1992         HP-UX Version with minijets
C
C                                     and multiple soft jets
C
C                        1995/6  introduction of Evaporation
C                                and residual nuclei
C                          (with S.Roesler, A.Ferrari, P.Sala)
C
C                        1997    extension to sqet(s)=2000 TeV
C                   Using GRV94LO and CTEQ96 structure functions
C
C
C             Up to Version 2.3 string fragmentation using
C                         the jetset-7.3 code
C
C              Version 2.4  uses the double precision
C                           PYTHIA-6.1 code for jet
C                           fragmentation
C
C********************************************************************
C
C
C                    D P M J E T
C
C
C
C                         I-------I
C                         I BEGIN I
C                         I-------I
C                             I
C                             I
C               I---------------------------I   I----------------I
C               I                           I   I COMMON BLOCKS  I
C               I      INITIALISATION       I   I----------------I
C               I                           I   I   BLOCK  DATA  I
C               I---------------------------I   I----------------I
C                             I
C                             I-------------------<------------I
C                             I                                I
C               I---------------------------I                  I
C               I                           I                  I
C               I    READ A CONTROL CARD    I                  I
C               I                           I                  I
C       ----------------------------------------------         I
C      /        /         I                \          \        I
C     /        /          I                 \          \       I
C I------I  I------I  I------I           I------I  I------I    I
C I 100  I  I 200  I  I 300  I   . . .   I 3900 I  I 4000 I    I
C ITITLE I  IPROJPAR  ITARPARI           I      I  I STOP I    I
C I------I  I------I  I------I           I------I  I------I    I
C    I         I         I                   I        I        I
C    I         I         I                   I     I------I    I
C    I         I         I                   I     I STOP I    I
C    I         I         I                   I     I------I    I
C    I         I         I                   I                 I
C    I--->---------->------------->------------------>---------I
C
C
C
C
C
C***********************************************************************
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
*KEEP,PANAME.
C------------------
C
C     /PANAME/ CONTAINS PARTICLE NAMES
C        BTYPE  = LITERAL NAME OF THE PARTICLE
C
      CHARACTER*8 BTYPE
      COMMON /PANAME/ BTYPE(30)
*KEEP,DINPDA.
      COMMON /DINPDA/ IMPS(6,6),IMVE(6,6),IB08(6,21),IB10(6,21), IA08
     +(6,21),IA10(6,21), A1,B1,B2,B3,LT,LE,BET,AS,B8,AME,DIQ,ISU
*KEEP,FACTMO.
      COMMON /FACTMO/ IFACTO
*KEEP,TAUFO.
      COMMON /TAUFO/  TAUFOR,KTAUGE,ITAUVE,INCMOD
*KEEP,RPTSHM.
      COMMON /RPTSHM/ RPROJ,RTARG,BIMPAC
*KEEP,TRAFOP.
      COMMON /TRAFOP/ GAMP,BGAMP,BETP
*KEEP,REJEC.
      COMMON /REJEC/ IRCO1,IRCO2,IRCO3,IRCO4,IRCO5, IRSS11,IRSS12,
     +IRSS13,IRSS14, IRSV11,IRSV12,IRSV13,IRSV14, IRVS11,IRVS12,IRVS13,
     +IRVS14, IRVV11,IRVV12,IRVV13,IRVV14
*KEEP,DPRIN.
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
*KEEP,DNUN.
      COMMON /DNUN/   NN,NP,NT
*KEEP,NSHMAK.
      COMMON /NSHMAK/ NNSHMA,NPSHMA,NTSHMA,NSHMAC,NSHMA2
*KEEP,DSHM.
      COMMON /DSHM/ RASH,RBSH,BMAX,BSTEP,SIGSH,ROSH,GSH,
     *              BSITE(0:1,200),NSTATB,NSITEB
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
*KEEP,DROPPT.
      LOGICAL INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA, IPADIS,
     +ISHMAL,LPAULI
      COMMON /DROPPT/ INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA,
     +IPADIS,ISHMAL,LPAULI
*KEEP,XSEADI.
      COMMON /XSEADI/ XSEACU,UNON,UNOM,UNOSEA, CVQ,CDQ,CSEA,SSMIMA,
     +SSMIMQ,VVMTHR
*KEEP,NUCIMP.
      COMMON /NUCIMP/ PRMOM(5,248),TAMOM(5,248), PRMFEP,PRMFEN,TAMFEP,
     +TAMFEN, PREFEP,PREFEN,TAEFEP,TAEFEN, PREPOT(210),TAEPOT(210),
     +PREBIN,TAEBIN,FERMOD,ETACOU
*KEEP,COULO.
      COMMON/COULO/ICOUL
*KEEP,EDENS.
      COMMON/EDENS/IEDEN
*KEEP,PROJK.
      COMMON /PROJK/ IPROJK
*KEEP,INTMX.
      PARAMETER (LUNBER=14)
      PARAMETER (INTMX=2488,INTMD=252)
*KEND.
      COMMON /SEAQXX/ SEAQX,SEAQXN 
       COMMON /CRONIN/CRONCO,MKCRON
       LOGICAL LSEADI
       COMMON /SEADIQ/LSEADI
       COMMON /FINAL/IFINAL
       COMMON /RECOM/IRECOM
       COMMON /HBOO/IHBOOK
      COMMON /NEUTYY/ NEUTYP,NEUDEC
       COMMON /NSTARI/NSTART
       COMMON/POPCOR/PDB,AJSDEF
       COMMON/POPCCK/PDBCK,PDBSE,PDBSEU,
     *  IJPOCK,IREJCK,ICK4,IHAD4,ICK6,IHAD6
     *,IREJSE,ISE4,ISE6,IREJS3,ISE43,ISE63,IREJS0
     *,IHADA4,IHADA6,IREJSA,ISEA4,ISEA6,IREJA3,
     *ISEA43,ISEA63,IREJAO
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
      COMMON /DIQSUM/NDVUU,NDVUS,NDVSS,NVDUU,NVDUS,NVDSS,
     *               NDSUU,NDSUS,NDSSS,NSDUU,NSDUS,NSDSS,
     *               NDZUU,NDZUS,NDZSS,NZDUU,NZDUS,NZDSS 
     *              ,NADVUU,NADVUS,NADVSS,NAVDUU,NAVDUS,NAVDSS,
     *               NADSUU,NADSUS,NADSSS,NASDUU,NASDUS,NASDSS,
     *               NADZUU,NADZUS,NADZSS,NAZDUU,NAZDUS,NAZDSS 
      COMMON /HDJASE/NHSE1,NHSE2,NHSE3,NHASE1,NHASE2,NHASE3
C---------------------
      COMMON /DIFFRA/ ISINGD,IDIFTP,IOUDIF,IFLAGD
      COMMON /SEASU3/SEASQ
      COMMON /IFRAGM/IFRAG
      COMMON /FLUCTU/IFLUCT
      COMMON /DIQUAX/AMEDD,IDIQUA,IDIQUU
      COMMON /INXDPM/INTDPM
      COMMON /VXSVD/VXSP(50),VXST(50),VXSAP(50),VXSAT(50),
     *              VXVP(50),VXVT(50),VXDP(50),VXDT(50),
     *	    NXSP,NXST,NXSAP,NXSAT,NXVP,NXVT,NXDP,NXDT
*
*KEEP,NNCMS.
      COMMON /NNCMS/  GAMCM,BGCM,UMO,PCM,EPROJ,PPROJ
C---------------------
C                               from DTUJET93
      COMMON /XSECPT/ PTCUT,SIGS,DSIGH
      COMMON /KGLAUB/JGLAUB
      COMMON /XSECNU/ECMUU,ECMOO,NGRITT,NEVTT
      DIMENSION PPNPN(4)
*
*   **********************************************************************
*   *  DESCRIPTION OF THE COMMON BLOCK(S), VARIABLE(S) AND DECLARATIONS  *
*   **********************************************************************
*
*
*     the following two COMMON blocks are prepared
*
*            *** for normal user ***
*
*     who wants to use his own histogramming, detector progams ect.
*     without going into details
*
*  *********************************************************************
*  /USER/ contains the parameters, expected to be modified by normal user
*      TITLE is a litteral string TITLE printet in the OUTPUT
*      PROJTY resp. TARGTY specify the type of particle scattering
*            The projectile moves in positive z-direction.
*          (Particle type specifications numbers for scatterers are stored
*           in COMMON /BOOKLT/ in BLOCKDATA on the end of this file.
*           Our comlete particle and resonance numbering is given in
*           the file DTUTCB in BLOCK DATA partic DATA ANAME
*           Also a list of our particle numbering
*           is obtained running the code word PARTICLE)
*      CMENERGY the center of mass energy in GeV
*      ISTRUF specifies the structure function as
C              ISTRUF=21:  GLUECK,REYA,VOGT GRV94LO with K=1.
C              ISTRUF=22:  GLUECK,REYA,VOGT GRV98LO with K=2.
C              ISTRUF=23:  CTEQ Collab. CTEQ96 with K=2.
*      ISINGD (ISINGX)specifies what is done with diffractive events
*              ISINGD=0: Single diffraction surpressed
*              ISINGD=1: Single diffraction included to fraction SDFRAC
*              ISINGD=2: Only single diffraction with target excited
*              ISINGD=3: Only single diffraction with projectile excited
*      IDUBLD specifies what is done with double diffractive events
*              ISINGD=0: Double diffraction included
*              ISINGD=1: Only double diffraction
*      SDFRAC see ISINGD
*      PTLAR  cutoff parameter requiring minijet of given size??

*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CHARACTER*80 TITLED
      CHARACTER*8 PROJTY,TARGTY
C     COMMON /USER/TITLED,PROJTY,TARGTY,CMENER,ISTRUF
C    &            ,ISINGX,IDUBLD,SDFRAC,PTLAR
      COMMON /USER1/TITLED,PROJTY,TARGTY
      COMMON /USER2/CMENER,SDFRAC,PTLAR,ISTRUF,ISINGX,IDUBLD
*
*  *********************************************************************
*
*
*     the following alphabetically ordered COMMON blocks are discribed
*
*                  *** for insider use ***
*
*     (parameters of interest to outsiders are explained below
*      in the description of "input card" input)
*
*
*  *********************************************************************
*     /COLLE/           contains the input specifying the MC. run
*        NEVHAD = is the number of events
C                 in older version NCASES
*        NVERS  = 1  all hard partons considered to be gluons
*                    rejection method for soft x-selection
*        NVERS  = 2  all hard partons considered to be gluons
*                    Aurenche-Maire method for soft x and pt selection
*        IHADRZ =
*        NFILE  =
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      COMMON /COLLE/NEVHAD,NVERS,IHADRZ,NFILE
*
*  *********************************************************************
*     /COLLIS/       contains the input specifying the considered event
C        ECM dropped as now in /USER/CMENER
*        S = is the Mandelstam s variable (=ECM**2)
*        IJPROJ,IJTARG = specifies the projectile rsp. target Q.N.
*        PTTHR  = the minimum pt still hard
*        PTTHR2 = the  pt of the first sampled hard scattering
*        IOPHRD = the option chosen for the hard scatterring
*        IJPRLU,IJTALU =
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      COMMON/COLLIS/S,IJPROX,IJTAR,PTTHR,PTTHR2,IOPHRD,IJPRLU,IJTALU
*
*  *********************************************************************
*     /BOOKLT/ contains the final  particle names and PPDB-numbers
*        BTYPEX = literal name of the particle
*        NBOOK  = the number of the particle
*                 proposed in the particle data booklet (90)
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CHARACTER*8 BTYPEX
      COMMON/BOOKLT/BTYPEX(30),NBOOK(30)
*  *********************************************************************
*     /POLMN/          stores arrays describing probabilities of parton
* split in /..0 et al.  configurations determined in
*        AILMN = "I(L,M,N)" amplidude matrix with L soft, M hard,
*                       N trippel pomeron exchanges used for ISIG=2,4..
*        PLMN   =       probability of L soft and M hard cut Pomeron
*                       and N cut tripple Pomeron (or similar)
*        PLMNCU =       cummulative PLMN
*        PDIFR =        probability that quasielastic is diffractive
*        PSOFT =        sum of PLMN(l,m=0,n)
*        PHARD =        1-AAAH
*        ALFAH =        SIGHIN/SIGIN
*        BETAH =        1-ALFAH
*        SIGTOT,SIGQEL,SIGEL,SIGINE,SIGHIN,SIGD,SIGDD =
*                       total, quasielastic, elastic, inelastic,
*                       hard inelsatic, single & double diffractiv X-sec
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       COMMON /POLMN0/PDIFR,PHARD,PSOFT,ALFAH,BETAH,
     *              SIGTOT,SIGQEL,SIGEL,SIGINE,SIGHIN,SIGD,SIGDD
*
*  *********************************************************************
*     /POMTYP/ contains parameters determining X-sections
*        IPIM,ICON,ISIG,DIFEL,DIFNU    resp.
*        IPIM,ICON,ISIG,LMAX,MMAX,NMAX as described at "CODEWD=SIGMAPOM"
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      COMMON /POMTYP/IPIM,ICON,ISIG,LMAX,MMAX,NMAX,DIFEL,DIFNU
*
*  *********************************************************************
*     various smaller commons
*     in alphabetical order
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C-dtp COMMON /DIFFRA/ISINGD,IDUBLD,SDFRAC dropped as in /USER/
      COMMON /DROPJJ/DROPJT,DROPVA
      COMMON /GLUSPL/NUGLUU,NSGLUU
C-dtu90      COMMON /PTLARG/PTLAR dropped as in /USER/
C-dtp      COMMON /PTLARG/PTLAR, XSMAX contains addition kept:
      COMMON /PTLARG/XSMAX
C-dtp addition next line:
      COMMON /PTSAMP/ ISAMPT
      COMMON /STARS/ISTAR2,ISTAR3
C-dtu90    COMMON /STRUFU/ISTRUF dropped as in /USER/
C-dtp      COMMON /STRUFU/ISTRUF,ISTRUM contains addition kept:
      COMMON /STRUFU/ISTRUM,ISTRUT
      COMMON /CUTOFN/NCUTOX
      COMMON/PSHOW/IPSHOW
C     COMMON/HARLUN/IHARLU,QLUN
      COMMON /HARLUN/ QLUN,IHARLU
      COMMON /POMTAB/IPOMTA   
      COMMON /SINCHA/ISICHAa
* evaporation module
      COMMON /EVAPPP/IEVAP
*$ CREATE PAREVT.ADD
      PARAMETER ( FRDIFF = 0.2D+00 )
      PARAMETER ( ETHSEA = 1.0D+00 )
 
      LOGICAL LDIFFR, LINCTV, LEVPRT, LHEAVY, LDEEXG, LGDHPR, LPREEX,
     &        LHLFIX, LPRFIX, LPARWV, LPOWER, LSNGCH, LLVMOD, LSCHDF
      COMMON / PAREVT / DPOWER, FSPRD0, FSHPFN, RN1GSC, RN2GSC,
     &                  LDIFFR (39),LPOWER, LINCTV, LEVPRT, LHEAVY,
     &                  LDEEXG, LGDHPR, LPREEX, LHLFIX, LPRFIX, LPARWV,
     &                  ILVMOD, JLVMOD, LLVMOD, LSNGCH, LSCHDF
**
*$ CREATE FRBKCM.ADD
      PARAMETER ( MXFFBK =     6 )
      PARAMETER ( MXZFBK =     9 )
      PARAMETER ( MXNFBK =    10 )
      PARAMETER ( MXAFBK =    16 )
      PARAMETER ( NXZFBK = MXZFBK + MXFFBK / 3 )
      PARAMETER ( NXNFBK = MXNFBK + MXFFBK / 3 )
      PARAMETER ( NXAFBK = MXAFBK + 1 )
      PARAMETER ( MXPSST =   300 )
      PARAMETER ( MXPSFB = 41000 )
      LOGICAL LFRMBK, LNCMSS
      COMMON / FRBKCM /  AMUFBK, EEXFBK (MXPSST), AMFRBK (MXPSST),
     &          EXFRBK (MXPSFB), SDMFBK (MXPSFB), COUFBK (MXPSFB),
     &          EXMXFB, R0FRBK, R0CFBK, C1CFBK, C2CFBK,
     &          IFRBKN (MXPSST), IFRBKZ (MXPSST),
     &          IFBKSP (MXPSST), IFBKPR (MXPSST), IFBKST (MXPSST),
     &          IPSIND (0:MXNFBK,0:MXZFBK,2), JPSIND (0:MXAFBK),
     &          IFBIND (0:NXNFBK,0:NXZFBK,2), JFBIND (0:NXAFBK),
     &          IFBCHA (5,MXPSFB), IPOSST, IPOSFB, IFBSTF,
     &          IFBFRB, NBUFBK, LFRMBK, LNCMSS
*$ CREATE INPFLG.ADD
      COMMON /INPFLG/ IANG,IFISS,IB0,IGEOM,ISTRAG,KEYDK
      COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
      COMMON /SECINT/ISECIN
       COMMON /NUCLEA/ PFERMP(2),PFERMN(2),FERMDD,
     & EBINDP(2),EBINDN(2),EPOT(2,210),
     &                ETACOO(2),ICOULL
      COMMON /FERFOR/IFERFO
      COMMON /CASADI/CASAXX,ICASAD
      COMMON /INFORE/IFREJ         

C---------------------
*
*  ********************************************************************
*
C                               from DTUJET93
C                               from DTUJET93

      CHARACTER*80 TITLE
      CHARACTER*8 CODE,CODEWD,BLANK,SDUM
      DIMENSION WHAT(6),CODE(65)
      DATA CODE/
     1'TITLE   ','PROJPAR ','TARPAR  ','ENERGY  ','HADRONIZ',
     2'RANDOMIZ','FERMI   ','EVENTAPE','START   ','PARTEV  ',
     3'INTPT   ','TECALBAM','RESONANC','VALVAL  ','COMMENT ',
     4'OUTLEVEL','LEPTOEVT','SEASEA  ','PARTICLE','ALLPART ',
     5'TAUFOR  ','SEAVAL  ','VALSEA  ','MOMENTUM','PAULI   ',
     6'PROJKASK','CENTRAL ','SEADISTR','CMHISTO ','SIGTEST ',
     7'XCUTS   ','HADRIN  ','FACTOMOM','COULOMB ','GLAUBERI',
     8'EDENSITY','CMENERGY','INFOREJE','RECOMBIN','SINGDIFF',
     9'NOFINALE','SEASU3  ','CRONINPT','POPCORN ','STOP    ',
     9'FLUCTUAT','DIQUARKS','HBOOKHIS','GLAUBERA','POMTABLE',
     9'SINGLECH','HADRINTH','EVAPORAT','SEAQUARK','SECINTER',
     9'POPCORCK','CASADIQU','POPCORSE','NEUTRINO','DIFFNUC ',
     9'XSECNUC ','INTERDPM','        ','        ','        '/
*
C-------------------
      DATA BLANK/'        '/
      DATA TITLE/' '/
C     DATA TITLED/' '/
C----------------------------------------------------------------------
      ISTART=0
      IEOF=0
C     OPEN(5,FILE='DTUNUC.DAT',STATUS='OLD')
C     OPEN(6,FILE='DTUNUC.OUT',STATUS='UNKNOWN')
C     OPEN(11,FILE='../dtunuc4/BBLO.DAT', STATUS='OLD')
C     OPEN(47,FILE='/u1/ranft/dtunuc44/GLAUBTAR.DAT',
C    *STATUS='UNKNOWN')
C     OPEN(47,FILE='/nfs/hptrack/user/ran/dtunuc44/GLAUBTAR.DAT',
C     OPEN(47,FILE='/user/ran/dtunuc44/GLAUBTAR.DAT',
C    *STATUS='UNKNOWN')
C     OPEN(47,FILE='/lapphp11_2/users/ranft/dtunuc44/GLAUBTAR.DAT',
C    *STATUS='UNKNOWN')
C     OPEN(47,FILE='GLAUBTAR.DAT',
C    *STATUS='UNKNOWN')
C     OPEN(37,FILE='GLAUBCROSS.DAT',
C    *STATUS='UNKNOWN')
C     OPEN( 2,FILE='HIBLD.DAT',STATUS='OLD')
C
C                               from DTUJET93
C                               from DTUJET93
C     OPEN(7,FILE='DPMJET.TOP')
C     OPEN(18,FILE='DPMJET.EVT')
C                               from DTUJET93
      IF (NCOUNT.EQ.1)THEN
*---initialization of  DECAY and HADRIN
        CALL DDATAR
        CALL DHADDE
        CALL DCHANT
        CALL DCHANH
*---print the title
      WRITE(6,1000)
 1000 FORMAT( '1    **************************************************',
     +'**************************************************', //
     +'     DPMJET VERSION II.5  (Sept. 1999)      ' /
     +'     DUAL PARTON MODEL FOR HADRON NUCLEUS COLLISIONS '/ /
     +'                AND NUCLEUS NUCLEUS COLLISIONS   '/
     +'      INCLUDING A FORMATION TIME INTRANUCLEAR CASCADE'/
     4'      Minijets and DTUJET like multiple soft jets    '/
     4'      Nuclear evaporation and residual target and    '/
     4'      projectile nuclei    '/
     +'     **************************************************',
     +'**************************************************',//)
      ENDIF
*---set default parameters not initialized in BLOCK DATA BLKDT1
      CALL DEFAUL(EPN,PPN)
      CALL DEFAUX(EPN,PPN)
C
C********************************************************************
      ICOUL=1
      ICOULL=1
C                                    EDENSITY
      IEDEN=0
C                                   TOPDRAW (option removed)
        ITOPD=0
C                                    TAUFOR
*---formation zone intranuclear cascade
      TAUFOR=105.D0
      KTAUGE=0 
      ITAUVE=1
      INCMOD=1
C                                     SEADISTR
*---definition of soft quark distributions
      XSEACO=1.00D0
      XSEACU=1.05-XSEACO
      UNON=3.50D0
      UNOM=1.11D0
      UNOSEA=5.0D0
C                                     FERMI
      FERMP=.TRUE.
      FERMOD=0.6D0
      FERMDD=0.6D0
      IFERFO=1
C                                     PAULI
      IPAUPR=0
      LPAULI=.TRUE.
C                                     XCUTS
*---cutoff parameters for x-sampling
      CVQ=1.8D0
      CDQ=2.0D0
      CSEA=0.5D0
      SSMIMA=1.201D0
      SSMIMQ=SSMIMA**2
      VVMTHR=0.D0
Cc                                     NOFINALE
      IFINAL=0
C                                      OUTLEVEL
      IPRI = 0
      IPEV = 0
      IPPA = 0
      IPCO = -2
      INIT = 0
      IPHKK= 0
C                                      RECOMBIN
      IRECOM=0
      LSEADI=.TRUE.
C                                     SEASU3
       SEASQ=0.50D0
C                                      CRONINPT
       MKCRON=1
       CRONCO=0.64D0
C                                      ALLPART
        IHADA=.TRUE.
C                                      INTERDPM
       INTDPM=0
       IROEH=0
C                                      POPCORCK
       PDBCK=0.
       IJPOCK=0
C                                      CASADIQU
       ICASAD=1
       CASAXX=0.5D0
C                                      POPCORSE
       PDBSE=0.45D0
       PDBSEU=0.45D0
C       
       IREJCK=0
       IREJSE=0
       IREJS3=0
       IREJS0=0
       ICK4=0
       ISE4=0
       ISE43=0
       IHAD4=0
       ICK6=0
       ISE6=0
       ISE63=0
       IHAD6=0
       IREJSA=0
       IREJA3=0
       IREJA0=0
       ISEA4=0
       ISEA43=0
       IHADA4=0
       ISEA6=0
       ISEA63=0
       IHADA6=0
C                                      POPCORN
       PDB=0.10D0
       AJSDEF=0.D0
C                                      FLUCTUAT
       IFLUCT=0
*                                      INTPT 
       INTPT=.TRUE.
*                                      HADRONIZ
        IHADRZ=2
        IFRAG=1
      IF (IHADRZ.GE.2)THEN
        IFRAG=IHADRZ-1
        CALL LUNDIN
      ENDIF
C                                      DIQUARKS
       IDIQUA=1
       IDIQUU=1
       AMEDD=0.9D0
C                                      SINGLECH
       ISICHA=0
C                                      EVAPORAT
       IEVAP=0
C                                      SEAQUARK
C                           sea quarks in multiple chains
       SEAQX=0.5D0
C                           sea quarks in Glauber events
       SEAQXN=0.5D0
C                                      GLAUBERI
C                                      GLAUBERA
C   change JGLAUB in dpmjet25 to JGLAUB=1 (was 2 in dpmjet241)      
       JGLAUB=1
C                                      HADRINTH
       EHADTH=5.D0
C                                      HBOOKHIS
       IHBOOK=1
C                                      POMTABLE
       IPOMTA=0   
C-------------------
C                               from DTUJET93
*  ********************************************************************
*    Initialize defaults for "input card" input parameters
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*
*     (to understand the meaning of parameters, we recommend
*     description of "input card" input provided below)
*
* *** the following initialisation involve normal user input cards
*
*                               CMENERGY
        CMENER=540.
        S=CMENER**2
          IF(ISTRUT.EQ.1)THEN
            PTTHR=2.1+0.15*(LOG10(CMENER/50.))**3
            PTTHR2=PTTHR
          ELSEIF(ISTRUT.EQ.2)THEN
            PTTHR=2.5+0.12*(LOG10(CMENER/50.))**3
            PTTHR2=PTTHR
          ENDIF
*                               PROJPAR
        PROJTY='PROTON  '
        IJPROJ=1
        IJPROX=1
	IBPROJ=1
	IP=1
	IPZ=1
        JJPROJ=1
        JJPROX=1
	JBPROJ=1
	JP=1
	JPZ=1
*                               MOMENTUM
      PPN=100000.
      NNPP=1
      IF(IJPROJ.NE.0) NNPP=IJPROJ
      EPN=SQRT(PPN**2+AAM(NNPP)**2)
*                             nucleon-nucleon cms
C     IBPROJ=1
      EPROJ=EPN
      AMPROJ=AAM(NNPP)
      AMTAR=AAM(1)
      PPROJ = SQRT((EPN-AMPROJ)*(EPN+AMPROJ))
      UMO = SQRT(AMPROJ**2 + AMTAR**2 + 2.*AMTAR*EPROJ)
      CMENER=UMO
      GAMCM = (EPROJ+AMTAR)/UMO
      BGCM=PPROJ/UMO
      ECM=UMO
      PCM=GAMCM*PPROJ - BGCM*EPROJ
          IF(ISTRUT.EQ.1)THEN
            PTTHR=2.1+0.15*(LOG10(CMENER/50.))**3
            PTTHR2=PTTHR
          ELSEIF(ISTRUT.EQ.2)THEN
            PTTHR=2.5+0.12*(LOG10(CMENER/50.))**3
            PTTHR2=PTTHR
          ENDIF
C
*                               TARPAR
        TARGTY=BLANK  
        IJTAR=1
	IT=14
	ITZ=7
	IJTARG=1
	IBTARG=1
        JJTAR=1
	JT=14
	JTZ=7
	JJTARG=1
	JBTARG=1
*                               CMHISTO
	CMHIS=0.D0
*                               CENTRAL
	ICENTR=0
*                               STRUCFUN
        ISTRUF=222
        ISTRUM=0
        ISTRUT=ISTRUF/100
        ISTRUF=ISTRUF-ISTRUT*100
        ISTRUM=ISTRUF
*                               SINGDIFF
        ISINGD=1
        ISINGX=1
        IDUBLD=0
        SDFRAC=1.
*                               START
C       NEVNTS passed to DTMAI as argument, NEVHAD in COMMON
        NFILE=0
        NSTART=1
        ISTAR2=0
        ISTAR3=0
        PTLAR=2.
	IGLAUB=0
*
* *** the following initialisation involve input cards for internal use
*        (in alphabetical order)
*
*                               INFOREJE
        IFREJ=0
*                               COMBIJET
C       NCUTOF=50
C       NCUTOX=50
*                               GLUSPLIT
        NUGLUU=1
        NSGLUU=0
*                               PARTEV
*       ITEST=0    (see SIGMAPOM)
        NPEV=30
        NVERS=1
*                               SAMPT
        ISAMPT=4
*                               SELHARD
*
        DROPJT=0.
*       ITEST=0    (see SIGMAPOM)
        IOPHRD=2
        PTTHR=3.
        PTTHR2=PTTHR
*                               SIGMAPOM
*
        ITEST=0
        IPIM=2
        ICON=48
        ISIG=10
C       only ISIG=10 possible since dpmjet--II.5	
C       default changed in relation to DPT =4
        LMAX=30
        MMAX=100
        NMAX= 2
        DIFEL = 0.
        DIFNU = 1.
*      some options use special routines for MMAX=0 so far NMAX=0
*                               PSHOWER
        IPSHOW=1
C                               SECINTER
	ISECIN=0
*                               RANDOM
        ISEED1=12
        ISEED2=34
        ISEED3=56
        ISEED4=78
C                               EVAPORATE
* set default if EVAP requested without "what-values"
      LEVPRT = .TRUE.
      ILVMOD = 1
      LDEEXG = .TRUE.
      LHEAVY = .TRUE.
      LFRMBK = .FALSE.
            LFRMBK   = .TRUE.
      IFISS  = 0
*      Initialize random generator
C       CALL RNDMST(ISEED1,ISEED2,ISEED3,ISEED4)
C*       test random generator (C as not to be understood by user)
C        CALL RNDMTE(1)
*
*                               starting parameters not read in
        ISTART=0
        XSMAX=0.8
        ITOPD=0
*
*  *********************************************************************
*
C                               from DTUJET93
C                               from DTUJET93
C
C********************************************************************
C               READ A CONTROL CARD
C
C
C     CONSTRUCTION OF THESE CARDS:
C     CODEWD  (WHAT(I),I=1,6)  SDUM
C     FORMAT(A8, 2X, 6E10.0, A8 )
C********************************************************************
C
   10 CONTINUE
      IF(IEOF.EQ.1)                                             GO TO 40
C      READ(5,1010,END=40)CODEWD,(WHAT(I),I=1,6),SDUM
      READ(5,1010)CODEWD,(WHAT(I),I=1,6),SDUM
      WRITE(6,1020)CODEWD,(WHAT(I),I=1,6),SDUM
      DO 20 ISW=1,65
*
        IF(CODEWD.EQ.CODE(ISW))                                 GO TO 30
   20 CONTINUE
      ISW=66
      WRITE(6,1030)
C                                                               GO TO 10
   30 GO TO(
C------------------------------------------------------------
C       TITLE   ,  PROJPAR ,  TARPAR  ,  ENERGY  ,  HADRONIZ,
     +  50      ,  60      ,  90      ,  120     ,  130     ,
C
C------------------------------------------------------------
C       RANDOMIZ,  FERMI   ,  EVENTAPE,  START   ,  PARTEV  ,
     +  140     ,  150     ,  160     ,  170     ,  210     ,
C
C------------------------------------------------------------
C       INTPT   ,  TECALBAM,  RESONANC,  VALVAL  ,  COMMENT ,
     +  220     ,  230     ,  240     ,  250     ,  260     ,
C
C------------------------------------------------------------
C       OUTLEVEL,  LEPTOEVT,  SEASEA  ,  PARTICLE,  ALLPART ,
     +  280     ,  290     ,  300     ,  310     ,  320     ,
C
C------------------------------------------------------------
C       TAUFOR  ,  SEAVAL  ,  VALSEA  ,  MOMENTUM,  PAULI   ,
     +  330     ,  340     ,  350     ,  360     ,  370     ,
C
C------------------------------------------------------------
C       PROJKASK,  CENTRAL ,  SEADISTR,  CMHISTO ,  SIGTEST ,
     +  380     ,  390     ,  400     ,  410     ,  420     ,
C
C------------------------------------------------------------
C       XCUTS   ,  HADRIN  , FACTOMOM , COULOMB  , GLAUBERI ,
     +  430     ,  440     , 450      , 460      , 470      ,
C
C------------------------------------------------------------
C       EDENSITY,  CMENERGY, INFOREJE , RECOMBIN , SINGDIFF )
     +  480     ,  490     , 500      , 510      , 520      ,
C
C------------------------------------------------------------
C       NOFINALE,  SEASU3    CRONINPT   POPCORN  , STOP  )
     +  530     ,  535     ,  538,     539,      540  ,    
C
C------------------------------------------------------------
C       FLUCTUAT,DIQUARKS HBOOKHIS,GLAUBERA,POMTABLE   )
     +  541     ,  542     ,  543,     544, 545,
C
C------------------------------------------------------------
Ch
C       SINGLECH, HADRINTH   ,EVAPORAT ,SEAQUARK, SECINTER  )
     +   551      , 552      , 553 ,  554      ,555 ,
C
C------------------------------------------------------------
C    POPCORCK   ,CASADIQU  ,POPCORSE ,NEUTRINO,DIFFNUC )
     +  556     , 557      ,558 ,    559   ,560,        
C
C------------------------------------------------------------
C    XSECNUC    , INTERDPM ,         ,        ,        )
     +  620     , 630      ,640 ,    650   ,660,610),ISW
C
C------------------------------------------------------------
*
                                                                GO TO 10
   40 CONTINUE
      WHAT(1)=0.0D0
      WHAT(2)=0.0D0
      WHAT(3)=0.0D0
      WHAT(4)=0.0D0
      WHAT(5)=0.0D0
      WHAT(6)=0.0D0
      SDUM=BLANK
      IEOF=1
      ISTART=1
      IF(ISTART.GT.0)THEN
        WRITE(6,1040)
                                                               GO TO 540
      ELSE
        WRITE(6,1050)
                                                               GO TO 170
      ENDIF
 1010 FORMAT(A8,2X,6E10.0,A8 )
 1020 FORMAT(' *****NEXT CONTROL CARD ***** ',A10,6(1X,G11.4), 2X,A10)
 
 1030 FORMAT(/,' UNKNOWN CODEWORD - CONTROL CARD IGNORED')
 1040 FORMAT(/,' UNEXPECTED END OF INPUT - STOP ASSUMED.')
 1050 FORMAT(/,' UNEXPECTED END OF INPUT - START ASSUMED.')
C
C********************************************************************
C               CONTROL CARD: CODEWD = TITLE
C               DEFINES THE TITLE OF THE JOB
C
C     WHAT(1...6),SDUM HAVE NO MEANING
C     THIS CARD MUST BE FOLLOWED BY THE CARD GIVING THE TITLE
C     OF THE RUN.
C********************************************************************
C
  610 CONTINUE
C                               from DTUJET93
*
*  The following CODEWD options are used by normal users:
*    to overwrite initial values:
*       CMENERGY    PROJPAR     TARPAR
*       STRUCFUN    SINGDIFF    
*    to order tasks:
*       START       STOP
*       PARTICLE    XSECTION   TITLE)+     COMMENT)+
*    Cards marked with )+ have to be followed by data cards of special format
*
*  Exemplaric imput cards inside dotts:
*...............................................................................
*.TITLE                                                                        .
*.TEST    NORMAL RUN    PTCUT= 3 GEV/C                                         .
*.PROJPAR                                                               APROTON.
*.CMENERGY       1800.                                                         .
*.START         10000.        0.        0.        0.      2.0                  .
*.STOP                                                                         .
*.........................................................ignore dotts..........
*
*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*
*
*********************************************************************
*                       printout in case input card to be considered
*                       (1st empty or dashed cards are ignored)
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF (CODEWD.GT.'-zzzzzzz')
     1    WRITE(6,91) CODEWD,(WHAT(I),I=1,6),SDUM
  91      FORMAT(' ---- control input card : ----'
     1        /1X,A8,2X,6(F10.3),A8)
 2    CONTINUE
*
*  *********************************************************************
*                       input card: CODEWD = STRUCFUN
*                       defines the structure functions for hard scatter
*
*     WHAT(1)  ISTRUF                                     default: 222
C              ISTRUF=21:  GLUECK,REYA,VOGT GRV94LO with K=1.
C              ISTRUF=22:  GLUECK,REYA,VOGT GRV98LO with K=2.
C              ISTRUF=23:  CTEQ Collab. CTEQ96 with K=2.
*              ISTRUF=ISTRUF+200  signal for Energy dependent ptcut
*****************************************************************
*            Only ISTRUF=222 makes sense  in dpmjet-II.5
*****************************************************************
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*
      IF(CODEWD.EQ.'STRUCFUN') THEN
*
        ISTRUF=INT(WHAT(1))
C       not in DT90 new in DTP
        ISTRUM=INT(WHAT(2))
        ISTRUT=ISTRUF/100
        ISTRUF=ISTRUF-ISTRUT*100
        ISTRUM=ISTRUF
        WRITE(6,*)' ISTRUF,ISTRUT ',ISTRUF,ISTRUT
      GO TO 10
*
*
*  *********************************************************************
*                       input card: CODEWD = PSHOWER
*                       Demands showering of hard partons
*                       Only with JETSET fragnentation
*           WHAT(1)=IPSHOW                         Default:1.
*                              Hard partons showering
*                              with IPSHOW=1 not showering with
*                              IPSHOW=0
*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*
      ELSEIF(CODEWD.EQ.'PSHOWER ') THEN
*
                IPSHOW=WHAT(1)
                GO TO 10
      ENDIF
       CALL DTTEST(CODEWD,WHAT,SDUM)
       GO TO 10
*  *********************************************************************
*
************************************************************************
C                               from DTUJET93
C                               from DTUJET93
C
C********************************************************************
C               CONTROL CARD: CODEWD =TITLE 
C               DEFINES THE TITLE OF THE JOB
C
C     WHAT(1...6),SDUM HAVE NO MEANING
C     THIS CARD MUST BE FOLLOWED BY THE CARD GIVING THE TITLE
C     OF THE RUN.
C********************************************************************
C
   50 CONTINUE
      READ(5,1060)TITLE
      TITLED=TITLE
      WRITE(6,1070)TITLE
                                                                GO TO 10
 1060 FORMAT(A80)
 1070 FORMAT(//,5X,A80,//)
C
C********************************************************************
C               CONTROL CARD: CODEWD = PROJPAR
C               DEFINES THE PROJECTILE PARTICLE / NUCLEUS
C
C     SDUM    =  IF DEFINED - PROJECTILE  PARTICLE TYPE
C  OTHERWIZE
C     WHAT(1) = IP  MASS NUMBER OF PROJECTILE NUCLEUS
C     WHAT(2) = IPZ CHARGE OF PROJECTILE NUCLEUS
C
C********************************************************************
C
   60 CONTINUE
      PROJTY=SDUM
      IF(SDUM.EQ.BLANK) THEN
        IP=INT(WHAT(1))
        IPZ=INT(WHAT(2))
        IBPROJ=1
        IJPROJ=1
        IJPROX=1
        IF(IP.EQ.1) IJPROJ=1
        IF(IP.EQ.1) IJPROX=1
        JP=INT(WHAT(1))
        JPZ=INT(WHAT(2))
        JBPROJ=1
        JJPROJ=1
        JJPROX=1
        IF(IP.EQ.1) JJPROJ=1
        IF(IP.EQ.1) JJPROX=1
      ELSE
        DO 70 II=1,30
          IF(SDUM.EQ.BTYPE(II)) THEN
            IJPROJ=II
            IJPROX=II
            IBPROJ=IIBAR(IJPROJ)
            IPZ=1
            IP=1
            JJPROJ=II
            JJPROX=II
            JBPROJ=IIBAR(IJPROJ)
            JPZ=1
            JP=1
                                                                 GOTO 80
          ENDIF
   70   CONTINUE
        WRITE(6,'(A)') ' WRONG STRUCTURE OF PROJPAR CARD'
        STOP
      ENDIF
   80 CONTINUE
                                                                GO TO 10
C
C********************************************************************
C               CONTROL CARD: CODEWD = TARPAR
C               DEFINES THE TARGET NUCLEUS
C
C     WHAT(1) = IT  MASS NUMBER OF TARGET NUCLEUS
C     WHAT(2) = ITZ CHARGE OF TARGET NUCLEUS
C
C********************************************************************
C
   90 CONTINUE
      TARGTY=SDUM
      IF(SDUM.EQ.BLANK) THEN
        IT=INT(WHAT(1))
        ITZ=INT(WHAT(2))
        JT=INT(WHAT(1))
        JTZ=INT(WHAT(2))
      ELSE
        DO 100 II=1,30
          IF(SDUM.EQ.BTYPE(II)) THEN
            ITZ=1
            IT=1
                  IJTAR=II
            JTZ=1
            JT=1
                  JJTAR=II
                                                                GOTO 110
          ENDIF
  100   CONTINUE
        WRITE(6,'(A)') ' WRONG STRUCTURE OF TARPAR CARD'
        STOP
      ENDIF
  110 CONTINUE
                                                                GO TO 10
C
C********************************************************************
C               CONTROL CARD: CODEWD = ENERGY
C               DEFINES THE LAB ENERGY OF THE PROJECTILE
C
C     WHAT(1) = LAB ENERGY IN GEV                DEFAULT: 200.
C
C********************************************************************
  120 CONTINUE
      EPN=WHAT(1)
      NNPP=1
      IF(IJPROJ.NE.0) NNPP=IJPROJ
      PPN=SQRT((EPN-AAM(NNPP))*(EPN+AAM(NNPP)))

*                             nucleon-nucleon cms
C     IBPROJ=1
      EPROJ=EPN
      AMPROJ=AAM(NNPP)
      AMTAR=AAM(1)
      PPROJ = SQRT((EPN-AMPROJ)*(EPN+AMPROJ))
      UMO = SQRT(AMPROJ**2 + AMTAR**2 + 2.*AMTAR*EPROJ)
      CMENER=UMO
          IF(ISTRUT.EQ.1)THEN
            PTTHR=2.1+0.15*(LOG10(CMENER/50.))**3
            PTTHR2=PTTHR
          ELSEIF(ISTRUT.EQ.2)THEN
            PTTHR=2.5+0.12*(LOG10(CMENER/50.))**3
            PTTHR2=PTTHR
          ENDIF
      GAMCM = (EPROJ+AMTAR)/UMO
      BGCM=PPROJ/UMO
      ECM=UMO
      PCM=GAMCM*PPROJ - BGCM*EPROJ
C
       PRINT 1033, EPROJ,PPROJ,
     +AMPROJ,AMTAR,UMO,GAMCM,BGCM,PCM
 1033 FORMAT(' CMS: ' , 
     +'    EPROJ,PPROJ,AMPROJ,AMTAR,UMO,GAMCM,BGCM,PCM'/8E22.13)
 
                                                                GO TO 10
C
C********************************************************************
C
C               CONTROL CARD: CODEWD = HADRONIZ
C               DEFINES THE CODE USED FOR HADRONIZATION
C
C     WHAT(1) = IHADRZ: DEFINES THE HADRONIZATION CODE DEFAULT: 2
C               IHADRZ=2 JETSET HADRONIZATION OF  CHAINS
C               IHADRZ=11 alternative JETSET HADRONIZATION OF  CHAINS
C     WHAT(2) =
C     WHAT(3) =
C
C
C********************************************************************
C
  130 CONTINUE
      IHADRZ=INT(WHAT(1))
      IFRAG=0
      IF (IHADRZ.GE.2)THEN
        IFRAG=IHADRZ-1
        CALL LUNDIN
      ENDIF
                                                                GO TO 10
C
C********************************************************************
C               CONTROL CARD: CODEWD = RANDOMIZE
C               SETS THE SEED FOR THE RANDOM NUMBER GENERATOR RM48
C
C   What(1): ISEED1
C   What(1): ISEED1
C
C
C********************************************************************
C
  140 CONTINUE
      ISEED1=WHAT(1)
      ISEED2=WHAT(2)
      AUAUAU=RD2IN(ISEED1,ISEED2)

                                                                GO TO 10
C
C********************************************************************
C               CONTROL CARD: CODEWD = FERMI
C               ALLOWS TO SWITCH FERMI MOTION ON OR OFF
C
C     WHAT(1) = 1. FERMI MOTION ON            DEFAULT ON
C               ELSE FERMI MOTION OFF
C     WHAT(2) = SCALE FACTOR FOR FERMI MOMENTUM (KKINC)
C                                   DEFAULT = 0.6
C     WHAT(3) = 1.  Zero temperature Fermi distr.    DEFAULT 1
C               2.   Function from C.C.d.Atti PRD53(96)1689
C
C********************************************************************
C
  150 CONTINUE
      IF (WHAT(1).EQ.1.D0)THEN
        FERMP=.TRUE.
      ELSE
        FERMP=.FALSE.
      ENDIF
      FERMDD=WHAT(2)
      FERMOD=WHAT(2)
      IF(FERMOD.LT.0.0D0.OR.SCAFER.GT.2.0D0) SCAFER=1.0D0
      IFERFO=1
      IF(WHAT(3).NE.0.D0)IFERFO=WHAT(3)
C
                                                                GO TO 10
C
C********************************************************************
C               CONTROL CARD: CODEWD = EVENTAPE
C               WRITES THE EVENTS TO A FILE  DTUJET EVENTS
C
C     WHAT(1) =
C     WHAT(2) =
C
C********************************************************************
C
  160 CONTINUE
C      WRITE (11,1080)
C      WRITE(11,1070)TITLE
 1080 FORMAT (' THIS FILE CONTAINS EVENTS FROM KKEVT ')
                                                                GO TO 10
C
C********************************************************************
C               CONTROL CARD: CODEWD = START
C               STARTS THE SAMPLING OF EVENTS AND STANDARD HISTOGRAM
C               OUTPUT
C
C     WHAT(1)=NUMBER OF EVENTS  NCASES,               DEFAULT: 1000.0
C     WHAT(2)=IGLAUB
*                   controls Glauber initialization   default: 0. *
*                   IGLAUB=1 : initialization via SHMAKI forced
C     WHAT(3)=MULTE  MULTE=0 normal run
C                    MULTE=1 Run with random energy 0.1*E - 2.*E
C
C********************************************************************
C
  170 CONTINUE
        NSTART=1
C---histogram initialization
      IF (IRESO.EQ.1)  CALL DISTRP(1,IJPROJ,PPN)
      IF (CMHIS.EQ.0.D0) CALL DISTR(1,IJPROJ,PPN,IDUMMY)
      IF (CMHIS.EQ.1.D0) CALL DISTRC(1,IJPROJ,PPN,IDUMMY)
      IF (CMHIS.EQ.2.D0) CALL DISTCO(1,IJPROJ,PPN,IDUMMY)
      IF (IRESO.EQ.1) CALL DISRES(1,NHKKH1,PPN)
      IF (IPADIS) CALL DISTPA(1)
      IF (IOUDIF.EQ.1) CALL DIADIF(1,0)
      CALL SHMAK(1,NN,NP,NT,IP,IT,UMO,BIMP)
      CALL SHMAK1(1,NN,NP,NT,IP,IT,UMO,BIMP)
* initialization of evaporation-module
      WRITE(6,*)' before NUCLEAR.BIN opened LUNBER= ',LUNBER
* initialization of evaporation-module
	OPEN(UNIT=LUNBER,FILE='NUCLEAR.BIN',STATUS='OLD'
     *     ,FORM='UNFORMATTED')
C	OPEN(UNIT=LUNBER,FILE='NUCLEAR.BIN',STATUS='OLD'
C     *    ,READONLY,FORM='UNFORMATTED')
C    *     ,FORM='UNFORMATTED')
      WRITE(6,*)'NUCLEAR.BIN opened LUNBER= ',LUNBER
        CALL BERTTP
      IF(IEVAP.EQ.1)THEN
C	CALL ZEROIN
C	OPEN(UNIT=LUNBER,FILE='NUCLEAR.BIN',STATUS='OLD'
C    *     ,READONLY,FORM='UNFORMATTED')
C       CALL BERTTP
        CALL INCINI
      ENDIF
C                                         HBOOK HISTOGRAMS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       IF(IHBOOK.EQ.1.AND.CMHIS.EQ.0.D0)THEN
         CALL PLOMB(1,PPNPN,CHAR,XFXFXF,ITIF,IJPROJ)
       ENDIF
       IF(IHBOOK.EQ.1.AND.CMHIS.EQ.1.D0)THEN
         CALL PLOMBC(1,PPNPN,CHAR,XFXFXF,ITIF,IJPROJ)
       ENDIF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     IF(NCOUNT.EQ.1) THEN
*---initialization and test of the random number generator
C     CALL RNDMST(12,34,56,78)
C     CALL RNDMTE(1)
C     ENDIF
*---CONSISTENCY TEST FOR FERMI/PAULI OPTIONS
      IF(LPAULI .AND. (.NOT.FERMP)) THEN
        WRITE(6,'(/2A/A/)')
     +  ' ACTIVATION OF PAULI PRINCIPLE ONLY REASONABLE',
     +  ' IF FERMI ACTIVE', ' LPAULI CHANGED TO .FALSE.'
      LPAULI=.FALSE.
      ENDIF
C     CLOSE (UNIT=5)
C     CLOSE (UNIT=11)
      NCASES=INT(WHAT(1))
                NEVNTS=INT(WHAT(1))
                IF(NEVNTS.LE.0) NEVNTS=1000
                NEVHAD=NEVNTS
C               why 2: NEVNTS passed to DTMAI as argument, NEVHAD in COMMON
C               in DTP NCASES <- NEVNTS
      IF(NCASES.LE.0) NCASES=100
      IGLAUB=INT(WHAT(2))
      MULTE=0
      MULTE=INT(WHAT(3))
C     IGLAUB=1
      IF(IGLAUB.NE.1) IGLAUB=0
*---INITIALIZE GLAUBER THEORY A LA SHMAKOV
      CALL TIMDAT
      IF(IGLAUB.EQ.1) THEN
C   change JGLAUB in dpmjet25 to JGLAUB=1 (was 2 in dpmjet241)      
        JGLAUB=1
        CALL SHMAKI(IP,IPZ,IT,ITZ,RPROJ,RTARG,PPN)
C       CALL SHMAKI(4,2,12,6,RPROJ,RTARG,PPN)
      ELSE
        CALL SHMAKF(IP,IPZ,IT,ITZ)
      ENDIF
C     WRITE(6,*)'call xsglau'
C     CALL XSGLAU(IP,IT,IJPROJ,1)       
*
*       Parton level initialisation
*       (PARTEV call prior to START no longer required)
C       userfriendly! , as DT90, not as DTP
*       initialize NSOFT-NHARD event selection
                IF(IPIM.EQ.2)CALL PRBLM2(CMENER)
*       initialize hard scattering
                CALL JTDTU(0)
*       initialize transverse momenta for soft scattering
                CALL SAMPPT(0,PT)
      NN=1
      NP=1
      NT=1
C                                   INITIALIZE COUNTERS
      BNNVV=0.001
      BNNSS=0.001
      BNNSV=0.001
      BNNVS=0.001
      BNNCC=0.001
      BNNDV=0.001
      BNNVD=0.001
      BNNDS=0.001
      BNNSD=0.001
      BNNHH=0.001
      BNNZZ=0.001
      BNNDI=0.001
      BNNZD=0.001
      BNNDZ=0.001
      BPTVV=0.
      BPTSS=0.
      BPTSV=0.
      BPTVS=0.
      BPTCC=0.
      BPTDV=0.
      BPTVD=0.
      BPTDS=0.
      BPTSD=0.
      BPTHH=0.
      BPTZZ=0.
      BPTDI=0.
      BPTZD=0.
      BPTDZ=0.
      BEEVV=0.
      BEESS=0.
      BEESV=0.
      BEEVS=0.
      BEECC=0.
      BEEDV=0.
      BEEVD=0.
      BEEDS=0.
      BEESD=0.
      BEEHH=0.
      BEEZZ=0.
      BEEDI=0.
      BEEZD=0.
      BEEDZ=0.
C      COMMON /NCOUCH/ ACOUVV,ACOUSS,ACOUSV,ACOUVS,
C    *                 ACOUZZ,ACOUHH,ACOUDS,ACOUSD,
C    *                 ACOUDZ,ACOUZD,ACOUDI
       BCOUVV=0.
       BCOUSS=0.
       BCOUSV=0.
       BCOUVS=0.
       BCOUZZ=0.
       BCOUHH=0.
       BCOUDS=0.
       BCOUSD=0.
       BCOUDZ=0.
       BCOUZD=0.
       BCOUDI=0.
       BCOUDV=0.
       BCOUVD=0.
       BCOUCC=0.
          IF(ISTRUT.EQ.1)THEN
            PTTHR=2.1+0.15*(LOG10(CMENER/50.))**3
            PTTHR2=PTTHR
          ELSEIF(ISTRUT.EQ.2)THEN
            PTTHR=2.5+0.12*(LOG10(CMENER/50.))**3
            PTTHR2=PTTHR
          ENDIF
      CALL TIMDAT
      RETURN
C
C********************************************************************
C               CONTROL CARD: CODEWD = PARTEV
C               DEMANDS HISTOGRAMS OF PARTONS AT CHAIN ENDS
C
C     WHAT(1) = 1.: HISTOGRAMS ;  OTHER VALUES: NO HISTOGRAMS
C     WHAT(2) =
C     WHAT(3) =
C     SDUM    =
C
C********************************************************************
C
  210 CONTINUE
      IF (WHAT(1).EQ.1.D0) THEN
        IPADIS=.TRUE.
      ELSE
        IPADIS=.FALSE.
      ENDIF
                                                                GO TO 10
C
C********************************************************************
C               CONTROL CARD: CODEWD = INTPT
C               SWITCHES INTRINSIC TRANSVERSE MOMENTA OF PARTONS
C                        ON AND OFF
C
C     WHAT(1) = 1. ON   ELSE OFF                          DEFAULT ON
C     WHAT(2) =
C
C********************************************************************
C
  220 CONTINUE
      IF (WHAT(1).EQ.1.D0) THEN
        INTPT=.TRUE.
      ELSE
        INTPT=.FALSE.
      ENDIF
                                                                GO TO 10
C
C********************************************************************
C
C               CONTROL CARD: CODEWD = TECALBAM
C               TESTS ROUTINES  FOR
C               SOFT JET FRAGMENTATION
C
C     TECALBAM DEMANDS ONE OR SEVERAL INPUTCARDS, WHICH SHOULD FOLLOW
C                           AFTER THE TECALBAM CARD
C
C********************************************************************
C
  230 CONTINUE
C     IF(NCOUNT.EQ.1) THEN
*---initialization and test of the random number generator
C     CALL RNDMST(12,34,56,78)
C     CALL RNDMTE(1)
C     ENDIF
      CALL TECALB
C     OPEN(26,FILE='PART.JS',
C    *STATUS='UNKNOWN')
C     CALL JSPARR
C     OPEN(25,FILE='JETSET.DEC',
C    *STATUS='UNKNOWN')
C     CALL PYUPDA(1,25)
                                                                GO TO 10
C
C********************************************************************
C               CONTROL CARD: CODEWD = RESONANC
C               SAMPLING OF RESONANCES IN FINAL STATE
C               NOTE: presently done before particle decay after KKEVT,
C                     i.e. before the treatment of secondary interactions
C     WHAT(1)
C
C********************************************************************
C
  240 CONTINUE
      IF(WHAT(1).GT.0.5D0) IRESO=1
                                                                GO TO 10
C
C********************************************************************
C               CONTROL CARD: CODEWD = VALVAL
C               DEMANDS SAMPLING OF HISTOGRAMS FROM HADRONIZED
C                       VALENCE-VALENCE STRINGS SEPARATELY
C
C     WHAT(1) = 1.          SAMPLING   ELSE NO SAMPLING   DEFAULT  NO
C
C********************************************************************
C
  250 CONTINUE
      IF (WHAT(1).EQ.1.D0) THEN
        IHADVV=.TRUE.
      ELSE
        IHADVV=.FALSE.
      ENDIF
                                                                GO TO 10
C
C********************************************************************
C               CONTROL CARD: CODEWD = COMMENT
C               MAKES POSSIBLE TO ADD COMMENTS IN THE INPUT
C
C     WHAT(1) = NUMBER OF COMMENT CARDS FOLLOWING THIS CONTROL CARD
C               (DEFAULT: 1 COMMENT CARD)
C
C********************************************************************
C
  260 CONTINUE
      NCOM=INT(WHAT(1))
      IF(NCOM.LE.0)NCOM=1
      DO 270 J=1,NCOM
        READ(5,1110)TITLE
        TITLED=TITLE
  270 WRITE(6,1120)TITLE
                                                                GO TO 10
 1110 FORMAT(A80)
 1120 FORMAT(1X,A80)
C
C********************************************************************
C               CONTROL CARD: CODEWD = OUTLEVEL
C               DEFINES THE AMOUNT OF OUTPUT WANTED
C
C     WHAT(1) = IPRI        DEFAULT = 0
C     WHAT(2) = IPEV        DEFAULT = 0
C     WHAT(3) = IPPA        DEFAULT = 0
C     WHAT(4) = IPCO        DEFAULT = -2
C     WHAT(5) = INIT        DEFAULT = 0
C     WHAT(6) = IPHKK       DEFAULT = 0
C
C               1.0 =  MINIMUM OUTPUT
C               2.0 =  VERY SHORT OUTPUT
C               3.0 =  SHORT OUTPUT
C               4.0 =  MEDIUM OUTPUT
C               5.0 =  LONG OUTPUT
C               6.0 =  VERY LONG OUTPUT
C               7.0 =  MAXIMUM OUTPUT
C               DEFAULT: 0.0
C
C********************************************************************
C
  280 CONTINUE
      IPRI  =INT(WHAT(1))
      IPEV  =INT(WHAT(2))
      IPPA  =INT(WHAT(3))
      IPCO  =INT(WHAT(4))
      INIT  =INT(WHAT(5))
      IPHKK =INT(WHAT(6))
                                                                GO TO 10
C
C********************************************************************
C
C               CONTROL CARD: CODEWD = LEPTOEVT
C               STARTS THE SAMPLING OF EVENTS AND STANDARD HISTOGRAM
C               OUTPUT FOR NEUTRINO INTERACTIONS (lepto code)
C
C     WHAT(1)=NUMBER OF EVENTS  NCASES,               DEFAULT: 1000.0
C     WHAT(2)= NEUTYP
C                     11=e-  -11=e+
C                     12=nu-e -12=anu-e
C                     13=mu- -13=mu+
C                     14=nu-mu -14=anu-mu
C
C
C********************************************************************
C
  290 CONTINUE
      OPEN(29,FILE='lepto.evt',
     *STATUS='UNKNOWN')
        NSTART=4
C---histogram initialization
      IF (IRESO.EQ.1)  CALL DISTRP(1,IJPROJ,PPN)
      IF (CMHIS.EQ.0.D0) CALL DISTR(1,IJPROJ,PPN,IDUMMY)
      IF (CMHIS.EQ.1.D0) CALL DISTRC(1,IJPROJ,PPN,IDUMMY)
      IF (CMHIS.EQ.2.D0) CALL DISTCO(1,IJPROJ,PPN,IDUMMY)
      IF (IRESO.EQ.1) CALL DISRES(1,NHKKH1,PPN)
      IF (IPADIS) CALL DISTPA(1)
      IF (IOUDIF.EQ.1) CALL DIADIF(1,0)
* initialization of evaporation-module
      WRITE(6,*)' before NUCLEAR.BIN opened LUNBER= ',LUNBER
* initialization of evaporation-module
	OPEN(UNIT=LUNBER,FILE='NUCLEAR.BIN',STATUS='OLD'
     *     ,FORM='UNFORMATTED')
C	OPEN(UNIT=LUNBER,FILE='NUCLEAR.BIN',STATUS='OLD'
C     *    ,READONLY,FORM='UNFORMATTED')
C    *     ,FORM='UNFORMATTED')
      WRITE(6,*)'NUCLEAR.BIN opened LUNBER= ',LUNBER
        CALL BERTTP
      IF(IEVAP.EQ.1)THEN
        CALL INCINI
	WRITE(6,*)' NEUTRINO: after INCINI call'
      ENDIF
C                                         HBOOK HISTOGRAMS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       IF(IHBOOK.EQ.1.AND.CMHIS.EQ.0.D0)THEN
         CALL PLOMB(1,PPNPN,CHAR,XFXFXF,ITIF,IJPROJ)
       ENDIF
       IF(IHBOOK.EQ.1.AND.CMHIS.EQ.1.D0)THEN
         CALL PLOMBC(1,PPNPN,CHAR,XFXFXF,ITIF,IJPROJ)
       ENDIF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     IF(NCOUNT.EQ.1) THEN
*---initialization and test of the random number generator
C     CALL RNDMST(12,34,56,78)
C     CALL RNDMTE(1)
C     ENDIF
*---CONSISTENCY TEST FOR FERMI/PAULI OPTIONS
      IF(LPAULI .AND. (.NOT.FERMP)) THEN
        WRITE(6,'(/2A/A/)')
     +  ' ACTIVATION OF PAULI PRINCIPLE ONLY REASONABLE',
     +  ' IF FERMI ACTIVE', ' LPAULI CHANGED TO .FALSE.'
      LPAULI=.FALSE.
      ENDIF
      NCASES=INT(WHAT(1))
      NEUTYP=INT(WHAT(2))
                NEVNTS=INT(WHAT(1))
                IF(NEVNTS.LE.0) NEVNTS=1000
                NEVHAD=NEVNTS
C               why 2: NEVNTS passed to DTMAI as argument, NEVHAD in COMMON
C               in DTP NCASES <- NEVNTS
      IF(NCASES.LE.0) NCASES=100
      CALL LUINOL
      CALL TIMDAT
*---INITIALIZE GLAUBER THEORY A LA SHMAKOV
C   change JGLAUB in dpmjet25 to JGLAUB=1 (was 2 in dpmjet241)      
      JGLAUB=1
      CALL SHMAKI(IP,IPZ,IT,ITZ,RPROJ,RTARG,PPN)
*
      NN=1
      NP=1
      NT=1
C                                   INITIALIZE COUNTERS
      BNNVV=0.001
      BNNSS=0.001
      BNNSV=0.001
      BNNVS=0.001
      BNNCC=0.001
      BNNDV=0.001
      BNNVD=0.001
      BNNDS=0.001
      BNNSD=0.001
      BNNHH=0.001
      BNNZZ=0.001
      BNNDI=0.001
      BNNZD=0.001
      BNNDZ=0.001
      BPTVV=0.
      BPTSS=0.
      BPTSV=0.
      BPTVS=0.
      BPTCC=0.
      BPTDV=0.
      BPTVD=0.
      BPTDS=0.
      BPTSD=0.
      BPTHH=0.
      BPTZZ=0.
      BPTDI=0.
      BPTZD=0.
      BPTDZ=0.
      BEEVV=0.
      BEESS=0.
      BEESV=0.
      BEEVS=0.
      BEECC=0.
      BEEDV=0.
      BEEVD=0.
      BEEDS=0.
      BEESD=0.
      BEEHH=0.
      BEEZZ=0.
      BEEDI=0.
      BEEZD=0.
      BEEDZ=0.
       BCOUVV=0.
       BCOUSS=0.
       BCOUSV=0.
       BCOUVS=0.
       BCOUZZ=0.
       BCOUHH=0.
       BCOUDS=0.
       BCOUSD=0.
       BCOUDZ=0.
       BCOUZD=0.
       BCOUDI=0.
       BCOUDV=0.
       BCOUVD=0.
       BCOUCC=0.
      CALL TIMDAT
      WRITE(6,*)' NEUTRINO initialization finished'
      RETURN
C
C
C
C********************************************************************
C
C               CONTROL CARD: CODEWD = SEASEA
C               DEMANDS SAMPLING OF HISTOGRAMS FROM HADRONIZED
C                       SEA-SEA STRINGS SEPARATELY
C
C
C     WHAT(1) = 1.          SAMPLING   ELSE NO SAMPLING   DEFAULT  NO
C
C********************************************************************
C
  300 CONTINUE
      IF (WHAT(1).EQ.1.D0) THEN
        IHADSS=.TRUE.
      ELSE
        IHADSS=.FALSE.
      ENDIF
                                                                GO TO 10
C
C********************************************************************
C
C               CONTROL CARD: CODEWD = PARTICLE
C            PRINTS A TABLE OF THE PROPERTIES AND DECAYS
C      OF THE PARTICLES DEFINED IN DTUJET BAMJET-DECAY FRAGMENTATION
C
C********************************************************************
C
  310 CONTINUE
      CALL DDATES
                                                                GO TO 10
C
C********************************************************************
C
C               CONTROL CARD: CODEWD = ALLPART
C               DEMANDS SAMPLING OF HISTOGRAMS FROM COMPLETE
C               HADRONIZED EVENTS
C
C     WHAT(1) = 1.          SAMPLING   ELSE NO SAMPLING   DEFAULT  1
C
C********************************************************************
C
  320 CONTINUE
      IF (WHAT(1).EQ.1.D0) THEN
        IHADA=.TRUE.
      ELSE
        IHADA=.FALSE.
      ENDIF
                                                                GO TO 10
C
C********************************************************************
C
C               CONTROL CARD: CODEWD = TAUFOR
C               READ THE FORMATION TIME
C
C     WHAT(1) = TAUFOR FORMATION TIME IN FERMI/C     DEFAULT 105.
C               TAUFOR=10 CORRESPONDS ROUGHELY TO AN AVARAGE
C               FORMATION TIME OF 1 FM/C
C     WHAT(2) = KTAUGE NUMER OF GENERATIONS ALLOWED  DEFAULT 0.
C                      IN FORMATION ZONE KASKADE MAXIMUM=10
C     WHAT(3) = ITAUVE ITAUVE = 1  pt dependent formation zone (Default)
C                      ITAUVE = 2  constant formation zone = TAUFOR
C
C********************************************************************
C
  330 CONTINUE
      TAUFOR=WHAT(1)
      KTAUGE=WHAT(2)
      ITAUVE=1
      IF(WHAT(3).NE.0.D0)ITAUVE=WHAT(3)
                                                                GO TO 10
C
C********************************************************************
C               CONTROL CARD: CODEWD = SEAVAL
C               DEMANDS SAMPLING OF HISTOGRAMS FROM HADRONIZED
C                       SEA-VALENCE STRINGS SEPARATELY
C
C
C     WHAT(1) = 1.          SAMPLING   ELSE NO SAMPLING   DEFAULT  NO
C
C********************************************************************
C
  340 CONTINUE
      IF (WHAT(1).EQ.1.D0) THEN
        IHADSV=.TRUE.
      ELSE
        IHADSV=.FALSE.
      ENDIF
                                                                GO TO 10
C
C********************************************************************
C               CONTROL CARD: CODEWD = VALSEA
C               DEMANDS SAMPLING OF HISTOGRAMS FROM HADRONIZED
C                       VALENCE-SEA STRINGS SEPARATELY
C
C     WHAT(1) = 1.          SAMPLING   ELSE NO SAMPLING   DEFAULT  NO
C
C********************************************************************
C
  350 CONTINUE
      IF (WHAT(1).EQ.1.D0) THEN
        IHADVS=.TRUE.
      ELSE
        IHADVS=.FALSE.
      ENDIF
                                                                GO TO 10
C
C********************************************************************
C               CONTROL CARD: CODEWD = MOMENTUM
C               DEFINES THE LAB MOMENTUM OF THE PROJECTILE
C
C     WHAT(1) = LAB MOMENTUM IN GEV/C                DEFAULT: 200.
C
C********************************************************************
  360 CONTINUE
      PPN=WHAT(1)
      NNPP=1
      IF(IJPROJ.NE.0) NNPP=IJPROJ
      EPN=SQRT(PPN**2+AAM(NNPP)**2)
*                             nucleon-nucleon cms
C     IBPROJ=1
      EPROJ=EPN
      AMPROJ=AAM(NNPP)
      AMTAR=AAM(1)
      PPROJ = SQRT((EPN-AMPROJ)*(EPN+AMPROJ))
      UMO = SQRT(AMPROJ**2 + AMTAR**2 + 2.*AMTAR*EPROJ)
      CMENER=UMO
          IF(ISTRUT.EQ.1)THEN
            PTTHR=2.1+0.15*(LOG10(CMENER/50.))**3
            PTTHR2=PTTHR
          ELSEIF(ISTRUT.EQ.2)THEN
            PTTHR=2.5+0.12*(LOG10(CMENER/50.))**3
            PTTHR2=PTTHR
          ENDIF
      GAMCM = (EPROJ+AMTAR)/UMO
      BGCM=PPROJ/UMO
      ECM=UMO
      PCM=GAMCM*PPROJ - BGCM*EPROJ
C
       PRINT 1033, EPROJ,PPROJ,
     +AMPROJ,AMTAR,UMO,GAMCM,BGCM,PCM
                                                                GO TO 10
C
C********************************************************************
C
C               CONTROL CARD: CODEWD = PAULI
C               MONITORING THE INCLUSION OF PAULI'S PRINCIPLE
C               FOR SECONDARY INTERACTIONS INSIDE THE NUCLEUS
C     WHAT(1) = 1. ON   ELSE OFF                  DEFAULT OFF
C     WHAT(2) GT 0.  VARIOUS TEST PRINTS          DEFAULT: NO PRINTS
C
C********************************************************************
C
  370 CONTINUE
      IF (WHAT(1).EQ.1.D0) THEN
        LPAULI=.TRUE.
      ELSE
        LPAULI=.FALSE.
      ENDIF
      IPAUPR=INT(WHAT(2))
                                                                GO TO 10
C
C********************************************************************
C               CONTROL CARD: CODEWD = PROJKASK
C   DEMANDS THE FORMATION ZONE CASCADE IN THE PROJECTILE NUCLEUS
C
C     WHAT(1) = IPROJK                               DEFAULT: 1
C                    0 = NO CASCADE
C                    1 = CASCADE
C********************************************************************
C
  380 CONTINUE
      IPROJK=WHAT(1)
                                                                GO TO 10
C
C********************************************************************
C               CONTROL CARD: CODEWD = CENTRAL
C
C  DEMANDS CENTRAL NUCLEUS-NUCLEUS COLLISIONS 
C
C     WHAT(1) = ICENTR                               DEFAULT: 0
C                    0 = NORMAL PRODUCTION
C                    1 = CENTRAL PRODUCTION (impact parameter condition
C                                            in dpm25nuc7.f and NA
C                                           condition in dpm25nuc2.f)
C                    2 = DIFFERENT CENTRAL PRODUCTION (Only  NA
C                                            condition in dpm25nuc2.f)
C                    3 = DIFFERENT CENTRAL PRODUCTION (Only  NA  
C                                            condition in dpm25nuc2.f)
C                                      Less central for Pb-Pb than 2
C                   10 = PERIPHERAL PRODUCTION
C********************************************************************
C
  390 CONTINUE
      ICENTR=WHAT(1)
                                                                GO TO 10
C
C********************************************************************
C               CONTROL CARD: CODEWD = SEADISTR
C   REDEFINES THE SEA DISTRIBUTIONS
C
C     WHAT(1) = XSEACO                               DEFAULT: 1.
C                  SEA(X)=1/X**XSEACO
C     WHAT(2) = UNON                                 DEFAULT: 3.5
C     WHAT(3) = UNOM                                 DEFAULT: 1.11
C     WHAT(4) = UNOSEA                               DEFAULT: 5.0
C                  QDIS(X) PROP. (1-X)**UNO...
C
C********************************************************************
C
  400 CONTINUE
      XSEACO=WHAT(1)
      XSEACU=1.05-XSEACO
      UNON=WHAT(2)
      IF(UNON.LT.0.1D0) UNON=2.0
      UNOM=WHAT(3)
      IF(UNOM.LT.0.1D0) UNOM=1.5
      UNOSEA=WHAT(4)
      IF(UNOSEA.LT.0.1D0) UNOSEA=2.0
                                                                GO TO 10
C
C********************************************************************
C
C               CONTROL CARD: CODEWD = CMHISTO
C   events in COMMON/HKKEVT/ are defined in the nucleon-nucleon cms
C   (to calculate HISTOGRAMS IN THE CMS)
C   NOTE: if CMHIS=1.0 no treatment of formation zone cascade
C         if CMHIS=0.0 formation zone cascade can be demanded,
C                      see TAUFOR
C
C     WHAT(1) = CMHIS                                DEFAULT: 0.
C
C********************************************************************
C
  410 CONTINUE
      CMHIS=WHAT(1)
                                                                GO TO 10
C********************************************************************
C
C               CONTROL CARD: CODEWD = SIGTEST
C   REQUIRES A TEST PRINT OF CROSS SECTION AS APPLIED IN THE PROGRAM
C
C********************************************************************
C
  420 CONTINUE
      CALL SIGTES
                                                                GO TO 10
C
C********************************************************************
C
C               CONTROL CARD: CODEWD = XCUTS
C
C     WHAT(1) = CVQ
C     WHAT(2) = CDQ
C     WHAT(3) = CSEA
C     WHAT(4) = SSMIMA
C     WHAT(5) = XXVTHR
C
C                         MINIMUM ALLOWED X-VALUES = CXXX/SQRT(S)
C
C********************************************************************
C
  430 CONTINUE
      CVQ=WHAT(1)
      IF(CVQ.LT.0.5D0) CVQ=1.0
      CDQ=WHAT(2)
      IF(CDQ.LT.1.0D0) CDQ=2.0
      CSEA=WHAT(3)
      IF(CSEA.LT.0.1D0) CSEA =0.1
      SSMIMA=WHAT(4)
      IF(SSMIMA.LT.0.0D0) SSMIMA=0.14
      SSMIMQ=SSMIMA**2
      IF(WHAT(5).GT.2.0D0) VVMTHR=WHAT(5)
                                                                GO TO 10
C
C********************************************************************
C
C               CONTROL CARD: CODEWD = HADRIN
C               TEST OPTION - DEFINES THE TYPE OF INTERACTION
C                             GENERATED BY HADRIN
C
C     WHAT(1) = 0.  : ELASTIC/INELASTIC AS MONITORED BY FOZOKL
C             = 1.  : ONLY INELASTIC INTERACTIONS
C             = 2.  : ONLY ELASTIC INTERACTIONS
C     WHAT(2)
C
C********************************************************************
C
  440 CONTINUE
      INTHAD=INT(WHAT(1))
      IF(INTHAD.LT.0 .OR. INTHAD.GT.2) INTHAD=0
      IF(INTHAD.EQ.1) WRITE(6,'(/5X,A/)')
     +' FHAD: INELASTIC INTERACTION FORCED'
      IF(INTHAD.EQ.2) WRITE(6,'(/5X,A/)')
     +' FHAD: ELASTIC INTERACTION FORCED'
                                                                GO TO 10
C
C********************************************************************
C
C               CONTROL CARD: CODEWD = FACTOMOM
C               DEMANDS THE CALCULATION OF FACTORIAL MOMENTSN
C
C     WHAT(1) = IFACTO =1 Factorial moments calculated   DEFAULT  0
C
C********************************************************************
C
  450 CONTINUE
      IFACTO=WHAT(1)
                                                                GO TO 10
C
C********************************************************************
C
C
C               CONTROL CARD: CODEWD = COULOMB
C               TO DEMAND COULOMB ENERGY FOR ICOUL=1
C
C     WHAT(1) = ICOUL                               DEFAULT=0
C     WHAT(2)
C     WHAT(3)
C
C
C********************************************************************
C
  460 CONTINUE
      ICOUL=INT(WHAT(1))
      ICOULL=INT(WHAT(1))
                                                                GO TO 10
C
C********************************************************************
C
C               CONTROL CARD: CODEWD = GLAUBERI
C
C              WHAT(1) = JGLAUB default:1
C   change JGLAUB in dpmjet25 to JGLAUB=1 (was 2 in dpmjet241)      
C                        JGLAUB = 1 prepare GLAUBTAR.DAT file
C
C       Test Glauber and output file for impact parameter selectiom
C       presently written for hadron-nucleus interactions
C
C       Actually 2000 events are used for initialization in SHMAKI
C       this is sufficient for most purposes, if not change NSTATB
C                                             in SHMAKI
C
C********************************************************************
C
  470 CONTINUE
      IF(WHAT(1).NE.0.D0)JGLAUB=WHAT(1)
C     CALL RNDMST(12,34,56,78)
C     CALL RNDMTE(1)
      IP=1.
      IPZ=1.
      JP=1.
      JPZ=1.
      Write(47,473)IP,IPZ,IT,ITZ
  473 FORMAT(' NUCLEUS  ',4I10)
      IJPROJ=1
      IJPROX=1
      IP=1
      IPZ=1
      JJPROJ=1
      JJPROX=1
      JP=1
      JPZ=1
      ISHC=0
      WU10=SQRT(10.)
      DO 471 IG=1,24
        PPN=WU10**(IG+1)
      CALL SHMAKI(IP,IPZ,IT,ITZ,RPROJ,RTARG,PPN)
        ISHC=ISHC+1
        IF(ISHC.EQ.1)THEN
        WRITE(47,'(4F10.5)') BMAX,BSTEP,RPROJ,RTARG
        ENDIF
        WRITE(47,'(5E16.8)') (BSITE(1,IB),IB=1,200)
  471 CONTINUE
      IJPROJ=13
      IJPROX=13
      JJPROJ=13
      JJPROX=13
      DO 472 IG=1,24
        PPN=WU10**(IG+1)
      CALL SHMAKI(IP,IPZ,IT,ITZ,RPROJ,RTARG,PPN)
        WRITE(47,'(5E16.8)') (BSITE(1,IB),IB=1,200)
  472 CONTINUE
                                                                GO TO 10
C
C********************************************************************
C
C               CONTROL CARD: CODEWD = EDENSITY
C               Calculate energy density for IEDEN=1
C               before resonance decay
C               Please note, that this option is inconsistent
C               with the formation zone cascade!
C
C     WHAT(1) = IEDEN                              DEFAULT=0
C     WHAT(2)
C
C********************************************************************
C
  480 CONTINUE
      IEDEN=WHAT(1)
                                                                GO TO 10
C
C********************************************************************
C
C               CONTROL CARD: CODEWD = CMENERGY
C               Input of nucleon-nucleon cms energy
C
C     WHAT(1) = CMENER  = SQRT(s)
C     WHAT(2)
C
C********************************************************************
C
  490 CONTINUE
      CMENER=WHAT(1)
          IF(ISTRUT.EQ.1)THEN
            PTTHR=2.1+0.15*(LOG10(CMENER/50.))**3
            PTTHR2=PTTHR
          ELSEIF(ISTRUT.EQ.2)THEN
            PTTHR=2.5+0.12*(LOG10(CMENER/50.))**3
            PTTHR2=PTTHR
          ENDIF
      S=CMENER**2
      UMO=CMENER
      NNPP=1
      IDP=NNPP
      IF(IJPROJ.NE.0) NNPP=IJPROJ
      EPN=(CMENER**2 + AAM(NNPP)**2 - AAM(1)**2)/(2.*AAM(1))
      PPN=SQRT((EPN-AAM(NNPP))*(EPN+AAM(NNPP)))
      ECM = SQRT(AAM(IDP)**2+AAM(1)**2+2.0D0*AAM(1)*EPN)      
      CMENER=ECM
      S=CMENER**2
      UMO=CMENER
      AMPROJ=AAM(NNPP)
      PPROJ = SQRT((EPN-AMPROJ)*(EPN+AMPROJ))
      EPROJ=SQRT(PPROJ**2+AMPROJ**2)
       EPROJ = EPN
       PPROJ = PPN
      AMTAR=AAM(1)
      GAMCM = (EPROJ+AMTAR)/UMO
      BGCM=PPROJ/UMO
      GAMCM = (EPROJ+AAM(1))/UMO
      BGCM = PPROJ/UMO
      PCM=GAMCM*PPROJ - BGCM*EPROJ
       PRINT 1033, EPROJ,PPROJ,
     +AMPROJ,AMTAR,UMO,GAMCM,BGCM,PCM
C     COMMON /NNCMS/  GAMCM,BGCM,UMO,PCM,EPROJ,PPROJ
                                                                GO TO 10
C
C********************************************************************
C
C               CONTROL CARD: CODEWD = INFOREJEC
C               Option to combine q-aq to color ropes
C
C     WHAT(1) = 0 :  no rejection diagnostics IFREJ=0
C               1 :  rejection diagnostics       DEFAULT:IFREJ= 0
C
C********************************************************************
C
  500 CONTINUE
      IFREJ=WHAT(1)
                                                                GO TO 10
C
C********************************************************************
C
C
C********************************************************************
C
C               CONTROL CARD: CODEWD = RECOMBIN
C               Option to s-s and v-v chains to s-v and v-s chains
C
C     WHAT(1) = 0 :  no recombination
C               1 :  recombination               DEFAULT: 0
C
C********************************************************************
C
  510 CONTINUE
      IF (WHAT(1).EQ.1.D0) IRECOM=1
      IF (WHAT(1).NE.1.D0) IRECOM=0
      IF (WHAT(1).EQ.1.D0) LSEADI=.TRUE.
                                                                GO TO 10
C  ********************************************************************
C                       parameter card: CODEWD = SINGDIFF
C                       Calls or supresses single diffractive events
C
C     WHAT(1)=ISINGD = 0 Single diffraction supressed       default: 0.
C                    = 1 Single diffraction included
C                    = 2 Only single Diffraction
C                    = 3 Only single Diffraction Target excited
C                    = 4 Only single Diffraction Projectile excited
C                    = 5 Only single Diffraction HMSD Target excited
C                    = 6 Only single Diffraction HMSD Projectile excited
C                    = 7 Only single Diffraction LMSD Target excited
C                    = 8 Only single Diffraction LMSD Projectile excited
C
C  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*
  520 CONTINUE
      ISINGD = WHAT(1)
      ISINGX = WHAT(1)
                                                                GO TO 10
*
C********************************************************************
C
C               CONTROL CARD: CODEWD = NOFINALE
C
C     WHAT(1) = 1 :  no FINALE call
C     WHAT(1) = 0 :     FINALE call
C                                                 DEFAULT: 0
C
C********************************************************************
C
  530 CONTINUE
      IF (WHAT(1).EQ.1.D0) IFINAL=1
      IF (WHAT(1).EQ.0.D0) IFINAL=0
                                                                GO TO 10
C
*
C********************************************************************
C
C               CONTROL CARD: CODEWD = SEASU3
C
C     WHAT(1) = SEASQ
C                                                 DEFAULT: 0.5
C
C********************************************************************
C
  535 CONTINUE
      SEASQ=WHAT(1)
                                                                GO TO 10
C
*
C********************************************************************
C
C               CONTROL CARD: CODEWD = CRONINPT
C
C     WHAT(1) = MKCRON
C                                                 DEFAULT: 1.
C     WHAT(2) = CRONCO
C                                                 DEFAULT: 0.64
C
C********************************************************************
C
  538 CONTINUE
      MKCRON=WHAT(1)
      CRONCO=WHAT(2)
                                                                GO TO 10
C
*
C********************************************************************
C
C               CONTROL CARD: CODEWD = POPCORN
C
C     WHAT(1) = PDB                                  DEFAULT: 0.10
C                       JETSET  fragmenting directly into baryons
C                       
C    WHAT(2) = AJSDEF                                DEFAULT=0.
C                        AJSDEF=1. SETS default values for all jetset
C                                   parameters
C
C********************************************************************
C
  539 CONTINUE
      PDB=WHAT(1)
      IF(WHAT(2).NE.0.D0)AJSDEF=WHAT(2)
                                                                GO TO 10
C******************************************************************
C               CONTROLCARD: CODEW = FLUCTUAT
C               Introduces cross section fluctuations 
C               a la Frankfurt and Strikman
C
C              WHAT(1) = IFLUCT                        default:0
C                        Fluctuations for IFLUCT=1
C
C*****************************************************************
C
 541  CONTINUE
      IFLUCT=WHAT(1)
      IF(IFLUCT.EQ.1)CALL FLUINI
      GO TO 10
C******************************************************************
C               CONTROLCARD: CODEW = DIQUARKS
C              WHAT(1) = IDIQUA Glauber diquarks       default:1
C              WHAT(2) = IDIQUU unitary diquarks       default:1
C              WHAT(3) = AMEDD                      default:0.9D0
C
C*****************************************************************
C
 542  CONTINUE
      IDIQUA=WHAT(1)     
      IDIQUU=WHAT(2)     
      IF(WHAT(3).GT.0.D0)THEN
        AMEDD=WHAT(3)
      ENDIF
      GO TO 10
C******************************************************************
C               CONTROLCARD: CODEW = HBOOKHIS
C
C              Serves to switch off HBOOK Histograms
C
C              WHAT(1) = IHBOOK                        default:1
C
C*****************************************************************
C
 543  CONTINUE
      IHBOOK=WHAT(1)
      GO TO 10
C
C********************************************************************
C
C               CONTROL CARD: CODEWD = GLAUBERA
C
C              WHAT(1) = JGLAUB default:1
C   change JGLAUB in dpmjet25 to JGLAUB=1 (was 2 in dpmjet241)      
C                        JGLAUB = 1 prepare GLAUBTAR.DAT file
C
C       Test Glauber and output file for impact parameter selectiom
C       presently written for nucleus-nucleus interactions
C
C       Actually 2000 events are used for initialization in SHMAKI
C       this is sufficient for most purposes, if not change NSTATB
C                                             in SHMAKI
C
C********************************************************************
C
  544 CONTINUE
      IF(WHAT(1).NE.0.D0)JGLAUB=WHAT(1)
C     CALL RNDMST(12,34,56,78)
C     CALL RNDMTE(1)
      Write(47,1473)IP,IPZ,IT,ITZ
 1473 FORMAT(' NUCLEUS  ',4I10)
      IJPROJ=1
      IJPROX=1
      JJPROJ=1
      JJPROX=1
      ISHC=0
      WU10=SQRT(10.)
      DO 1471 IG=1,24
        PPN=WU10**(IG+1)
      CALL SHMAKI(IP,IPZ,IT,ITZ,RPROJ,RTARG,PPN)
        ISHC=ISHC+1
        IF(ISHC.EQ.1)THEN
        WRITE(47,'(4F10.5)') BMAX,BSTEP,RPROJ,RTARG
        ENDIF
        WRITE(47,'(5E16.8)') (BSITE(1,IB),IB=1,200)
 1471 CONTINUE
      IJPROJ=13
      IJPROX=13
      JJPROJ=13
      JJPROX=13
      DO 1472 IG=1,24
        PPN=WU10**(IG+1)
      CALL SHMAKI(IP,IPZ,IT,ITZ,RPROJ,RTARG,PPN)
        WRITE(47,'(5E16.8)') (BSITE(1,IB),IB=1,200)
 1472 CONTINUE
                                                                GO TO 10
C
C******************************************************************
C               CONTROLCARD: CODEW =POMTABLE
C
C              WHAT(1)=IPOMTA            DEFAULT:0.  
c                      IPOMTA=0    POMERON TABLES FOR 28 ENERGIES
C                                  CALCULATED AND WRITTEN TO
C                                  filE pomtab.dat
C                      IPOMTA=1    POMERON TABLES READ FROM
C                                  pomtab.dat   
C
C*****************************************************************
C
 545  CONTINUE
      IPOMTA=WHAT(1)    
      GO TO 10
C******************************************************************
C               CONTROLCARD: CODEW = SINGLECH
C           include Regge contributions (single Chains)
C
C             WHAT(1)=ISICHA                     Default:0
C                     ISICHA=1: Single chains included
C
C*****************************************************************
C
 551  CONTINUE
      ISICHA=WHAT(1)
      GO TO 10
C******************************************************************
C               CONTROLCARD: CODEW =HADRINTH
C
C           WHAT(1) = EHADTH                   Default: 5. 
C
C*****************************************************************
C
 552  CONTINUE
      EHADTH=WHAT(1)
      WRITE(6,'(A,F10.2)')' Threshold for HADRIN events = (GeV)',
     *                      EHADTH
      GO TO 10
C******************************************************************
C               CONTROLCARD: CODEW =EVAPORAT 
C
C            These were the options before May 1995:
C
C                  evaporation module
C
C       what (1) >=  1 ==> evaporation is performed                
C                    if what (1) >= 10 then the high energy fis-   
C                    sion module is invoked and what(1)-10 is used 
C                    according to the previous table               
C       what (2) =< -1 ==> deexcitation gammas are not produced    
C                      (if the evaporation step is not performed   
C                       they are never produced)                   
C
C
C*****************************************************************
C
*======================================================================*
*     May 95:                                                          *
*     For the "evap" model:                                            *
*                                                                      *
*     Evaporation is performed if the EVAPORAT card is present
*                                             Default:No evaporatiom
*                                      
*              what (1) = i1 + i2*10 + i3*100 + i4*10000    
*                         (i1, i2, i3, i4 >= 0 )                       *
*                         i1 is the flag for selecting the T=0 level   *
*                         density option used                          *
*                  i1  =  1: standard EVAP level densities with Cook   *
*                            pairing energies                          *
*                      =  2: Z,N-dependent Gilbert & Cameron level     *
*                            densities (default)                       *
*                      =  3: Julich A-dependent level densities        *
*                      =  4: Z,N-dependent Brancazio & Cameron level   *
*                            densities                                 *
*                  i2 >=  1: high energy fission activated (default    *
*                            high energy fission activated)            *
*                  i3  =  0: No energy dependence for level densities  *
*                      =  1: Standard Ignyatuk (1975, 1st) energy de-  *
*                            pendence for level densities (default)    *
*                      =  2: Standard Ignyatuk (1975, 1st) energy de-  *
*                            pendence for level densities with NOT used*
*                            set of parameters                         *
*                      =  3: Standard Ignyatuk (1975, 1st) energy de-  *
*                            pendence for level densities with NOT used*
*                            set of parameters                         *
*                      =  4: Second   Ignyatuk (1975, 2nd) energy de-  *
*                            pendence for level densities              *
*                      =  5: Second   Ignyatuk (1975, 2nd) energy de-  *
*                            pendence for level densities with fit 1   *
*                            Iljinov & Mebel set of parameters         *
*                      =  6: Second   Ignyatuk (1975, 2nd) energy de-  *
*                            pendence for level densities with fit 2   *
*                            Iljinov & Mebel set of parameters         *
*                      =  7: Second   Ignyatuk (1975, 2nd) energy de-  *
*                            pendence for level densities with fit 3   *
*                            Iljinov & Mebel set of parameters         *
*                      =  8: Second   Ignyatuk (1975, 2nd) energy de-  *
*                            pendence for level densities with fit 4   *
*                            Iljinov & Mebel set of parameters         *
*                  i4 >=  1: Original Gilbert and Cameron pairing ener-*
*                            gies used (default Cook's modified pairing*
*                            energies)                                 *
*            what (2) = ig + 10 * if   (ig and if must have the same   *
*                                       sign)                          *
*                  ig =< -1 ==> deexcitation gammas are not produced   *
*                           (if the evaporation step is not performed  *
*                            they are never produced)                  *
*                  if =< -1 ==> Fermi Break Up is not invoked          *
*                           (if the evaporation step is not performed  *
*                            it is never invoked)                      *
*                        The default is: deexcitation gamma produced   *
*                        and Fermi break up activated for the new      *
*                        preequilibrium, not activated otherwise.      *
*            what (3) >=  1 ==> "heavies" put on their own stack       *
*                                                                      *
*======================================================================*
*
 553  CONTINUE
      WHTSAV = WHAT (1)
      IF ( NINT (WHAT (1)) .GE. 10000 ) THEN
         LLVMOD   = .FALSE.
         JLVHLP   = NINT (WHAT (1)) / 10000
         WHAT (1) = WHAT (1) - 10000.D+00 * JLVHLP
      END IF
      IF ( NINT (WHAT (1)) .GE. 100 ) THEN
         JLVMOD   = NINT (WHAT (1)) / 100
         WHAT (1) = WHAT (1) - 100.D+00 * JLVMOD
      END IF
      IF ( NINT (WHAT (1)) .GE. 10  ) THEN
         IFISS    = 1
         JLVHLP   = NINT (WHAT (1)) / 10
         WHAT (1) = WHAT (1) - 10.D+00 * JLVHLP
      ELSE IF ( NINT (WHTSAV) .NE. 0 ) THEN
         IFISS    = 0
      END IF
      IF ( NINT (WHAT (1)) .GE. 0 ) THEN
         LEVPRT = .TRUE.
         ILVMOD = NINT (WHAT(1))
         IF ( ABS (NINT (WHAT (2))) .GE. 10  ) THEN
            LFRMBK   = .TRUE.
            JLVHLP   = NINT (WHAT (2)) / 10
            WHAT (2) = WHAT (2) - 10.D+00 * JLVHLP
         ELSE IF ( NINT (WHAT (2)) .NE. 0 ) THEN
            LFRMBK   = .FALSE.
         END IF
         IF ( NINT (WHAT (2)) .GE. 0 ) THEN
            LDEEXG = .TRUE.
         ELSE
            LDEEXG = .FALSE.
         END IF
         IF ( NINT (WHAT(3)) .GE. 1 ) THEN
            LHEAVY = .TRUE.
         ELSE
            LHEAVY = .FALSE.
         END IF
      ELSE
         LEVPRT = .FALSE.
         LDEEXG = .FALSE.
         LHEAVY = .FALSE.
      END IF
C--------------------------------------------------------------
C--------------------------------------------------------------
C--------------------------------------------------------------
      IF(WHAT(1).EQ.0.)THEN
        IEVAP=0
        LEVPRT = .FALSE.
        ILVMOD = 1
        LDEEXG = .FALSE.
        LHEAVY = .FALSE.
        LFRMBK = .FALSE.
        IFISS  = 0
	GO TO 10
      ENDIF
      IEVAP=1
* set default if EVAP requested without "what-values"
      LEVPRT = .TRUE.
      ILVMOD = 1
      LDEEXG = .TRUE.
      LHEAVY = .TRUE.
      LFRMBK = .FALSE.
            LFRMBK   = .TRUE.
      IFISS  = 0
* check if fission is requested
      IF ( NINT (WHAT (1)) .GE. 10 ) THEN
         IFISS = 1
         WHAT (1) = WHAT (1) - 10.D+00
      END IF
* get level density treatment
      IF ( NINT (WHAT(1)) .GE. 1 ) ILVMOD = NINT (WHAT(1))
* switch off deexcitation gammas
      IF ( NINT (WHAT(2)) .LT. 0 ) LDEEXG = .FALSE.
* check if heavy recoil treatment is requested
*sr 31.1.95: since heavy fragments are always included this is obsolete
C     IF ( NINT (WHAT(3)) .GE. 1 ) LHEAVY = .TRUE.

      GO TO 10
C******************************************************************
C               CONTROLCARD: CODEW = SEAQUARK
C
C             SEAQX=WHAT(1)                    Default:0.5
C             SEAQXN=WHAT(2)                    Default:0.5
C
C*****************************************************************
C
 554  CONTINUE
      SEAQX=WHAT(1)
      SEAQXN=WHAT(2)
      GO TO 10
C******************************************************************
C               CONTROLCARD: CODEW = SECINTER
C
C            Controls secondary interaction of final state
C            particles of type pi + N ---> K + Hyperon
C
C            WHAT(1) = ISECIN             Default 0
C
C            ISECIN=1 demands these secondary interactions
C
C*****************************************************************
C
 555  CONTINUE
      ISECIN=WHAT(1)
      GO TO 10
C******************************************************************
C               CONTROLCARD: CODEW =POPCORCK 
C
C              Capella-Kopeliovich POPCORN effect
C
C               WHAT(1)= IJPOCK                  Default 0
C               WHAT(2)= PDBCK                   Default 0.0
C
C*****************************************************************
C
 556  CONTINUE
      IJPOCK=WHAT(1)
      PDBCK=WHAT(2)
      GO TO 10
C******************************************************************
C               CONTROLCARD: CODEW =CASADIQU       
C
C               Casado sea diquarks
C
C                     WHAT(1)= ICASAD    default : 1
C                     WHAT(2)= CASAXX   default : 0.5
C
C
C*****************************************************************
C
 557  CONTINUE
      ICASAD=WHAT(1)
      IF(WHAT(2).GE.0.1D0)THEN
        CASAXX=WHAT(2)
      ENDIF
      GO TO 10
C******************************************************************
C               CONTROLCARD: CODEW =POPCORSE 
C
C              Diquark breaking POPCORN Sea-effect
C
C               WHAT(1)= PDBSE                   Default 0.45
C               WHAT(2)= PDBSEU                  Default 0.45
C
C
C*****************************************************************
C
 558  CONTINUE
      PDBSE=WHAT(1) 
      PDBSEU=WHAT(2) 
      GO TO 10
C
C********************************************************************
C               CONTROL CARD: CODEWD = NEUTRINO
C               STARTS THE SAMPLING OF EVENTS AND STANDARD HISTOGRAM
C               OUTPUT FOR NEUTRINO INTERACTIONS (qel code)
C
C     WHAT(1)=NUMBER OF EVENTS  NCASES,               DEFAULT: 1000.0
C     WHAT(2)= NEUTYP
C      (1=nu-2, 2=anu-e, 3=nu-mu, 4=anu-mu, 5=nu-tau, 6=anu-tau)
C     WHAT(3)= NEUDEC                                DEFAULT:  0
C              Tau lepton decays into mu... for NEUDEC=1
C              Tau lepton decays into e... for NEUDEC=2
C                         decays not for NEUDEC=0
C             NEUDEC=10  call Gen_Delta CC events on p and n
C             NEUDEC=11  call Gen_Delta NC events on p and n
C             NEUDEC=20  call FILENU to read nu-N interactions
C            from a file (each event with different energy)
C
C********************************************************************
C
  559 CONTINUE
      
      OPEN(29,FILE='qel.evt',
     *STATUS='UNKNOWN')
        NSTART=2
C---histogram initialization
      IF (IRESO.EQ.1)  CALL DISTRP(1,IJPROJ,PPN)
      IF (CMHIS.EQ.0.D0) CALL DISTR(1,IJPROJ,PPN,IDUMMY)
      IF (CMHIS.EQ.1.D0) CALL DISTRC(1,IJPROJ,PPN,IDUMMY)
      IF (CMHIS.EQ.2.D0) CALL DISTCO(1,IJPROJ,PPN,IDUMMY)
      IF (IRESO.EQ.1) CALL DISRES(1,NHKKH1,PPN)
      IF (IPADIS) CALL DISTPA(1)
      IF (IOUDIF.EQ.1) CALL DIADIF(1,0)
* initialization of evaporation-module
      WRITE(6,*)' before NUCLEAR.BIN opened LUNBER= ',LUNBER
* initialization of evaporation-module
	OPEN(UNIT=LUNBER,FILE='NUCLEAR.BIN',STATUS='OLD'
     *     ,FORM='UNFORMATTED')
C	OPEN(UNIT=LUNBER,FILE='NUCLEAR.BIN',STATUS='OLD'
C     *    ,READONLY,FORM='UNFORMATTED')
C    *     ,FORM='UNFORMATTED')
      WRITE(6,*)'NUCLEAR.BIN opened LUNBER= ',LUNBER
        CALL BERTTP
      IF(IEVAP.EQ.1)THEN
C	CALL ZEROIN
C	OPEN(UNIT=LUNBER,FILE='NUCLEAR.BIN',STATUS='OLD'
C    *     ,READONLY,FORM='UNFORMATTED')
C       CALL BERTTP
        CALL INCINI
	WRITE(6,*)' NEUTRINO: after INCINI call'
      ENDIF
C                                         HBOOK HISTOGRAMS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       IF(IHBOOK.EQ.1.AND.CMHIS.EQ.0.D0)THEN
         CALL PLOMB(1,PPNPN,CHAR,XFXFXF,ITIF,IJPROJ)
       ENDIF
       IF(IHBOOK.EQ.1.AND.CMHIS.EQ.1.D0)THEN
         CALL PLOMBC(1,PPNPN,CHAR,XFXFXF,ITIF,IJPROJ)
       ENDIF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     IF(NCOUNT.EQ.1) THEN
*---initialization and test of the random number generator
C     CALL RNDMST(12,34,56,78)
C     CALL RNDMTE(1)
C     ENDIF
*---CONSISTENCY TEST FOR FERMI/PAULI OPTIONS
      IF(LPAULI .AND. (.NOT.FERMP)) THEN
        WRITE(6,'(/2A/A/)')
     +  ' ACTIVATION OF PAULI PRINCIPLE ONLY REASONABLE',
     +  ' IF FERMI ACTIVE', ' LPAULI CHANGED TO .FALSE.'
      LPAULI=.FALSE.
      ENDIF
      NCASES=INT(WHAT(1))
      NEUTYP=INT(WHAT(2))
      NEUDEC=INT(WHAT(3))
                NEVNTS=INT(WHAT(1))
                IF(NEVNTS.LE.0) NEVNTS=1000
                NEVHAD=NEVNTS
C               why 2: NEVNTS passed to DTMAI as argument, NEVHAD in COMMON
C               in DTP NCASES <- NEVNTS
      IF(NCASES.LE.0) NCASES=100
      CALL TIMDAT
*---INITIALIZE GLAUBER THEORY A LA SHMAKOV
C   change JGLAUB in dpmjet25 to JGLAUB=1 (was 2 in dpmjet241)      
      JGLAUB=1
      CALL SHMAKI(IP,IPZ,IT,ITZ,RPROJ,RTARG,PPN)
*
      NN=1
      NP=1
      NT=1
C                                   INITIALIZE COUNTERS
      BNNVV=0.001
      BNNSS=0.001
      BNNSV=0.001
      BNNVS=0.001
      BNNCC=0.001
      BNNDV=0.001
      BNNVD=0.001
      BNNDS=0.001
      BNNSD=0.001
      BNNHH=0.001
      BNNZZ=0.001
      BNNDI=0.001
      BNNZD=0.001
      BNNDZ=0.001
      BPTVV=0.
      BPTSS=0.
      BPTSV=0.
      BPTVS=0.
      BPTCC=0.
      BPTDV=0.
      BPTVD=0.
      BPTDS=0.
      BPTSD=0.
      BPTHH=0.
      BPTZZ=0.
      BPTDI=0.
      BPTZD=0.
      BPTDZ=0.
      BEEVV=0.
      BEESS=0.
      BEESV=0.
      BEEVS=0.
      BEECC=0.
      BEEDV=0.
      BEEVD=0.
      BEEDS=0.
      BEESD=0.
      BEEHH=0.
      BEEZZ=0.
      BEEDI=0.
      BEEZD=0.
      BEEDZ=0.
       BCOUVV=0.
       BCOUSS=0.
       BCOUSV=0.
       BCOUVS=0.
       BCOUZZ=0.
       BCOUHH=0.
       BCOUDS=0.
       BCOUSD=0.
       BCOUDZ=0.
       BCOUZD=0.
       BCOUDI=0.
       BCOUDV=0.
       BCOUVD=0.
       BCOUCC=0.
      CALL TIMDAT
      WRITE(6,*)' NEUTRINO initialization finished'
      RETURN
C
C
C********************************************************************
C               CONTROL CARD: CODEWD = DIFFNUC
C               STARTS THE SAMPLING OF EVENTS AND STANDARD HISTOGRAM
C               OUTPUT FOR NEUTRINO INTERACTIONS
C
C     WHAT(1)=NUMBER OF EVENTS  NCASES,               DEFAULT: 1000.0
C
C********************************************************************
C
  560 CONTINUE
      OPEN(29,FILE='diffnuc.evt',
     *STATUS='UNKNOWN')
        NSTART=3
C---histogram initialization
      IF (IRESO.EQ.1)  CALL DISTRP(1,IJPROJ,PPN)
      IF (CMHIS.EQ.0.D0) CALL DISTR(1,IJPROJ,PPN,IDUMMY)
      IF (CMHIS.EQ.1.D0) CALL DISTRC(1,IJPROJ,PPN,IDUMMY)
      IF (CMHIS.EQ.2.D0) CALL DISTCO(1,IJPROJ,PPN,IDUMMY)
      IF (IRESO.EQ.1) CALL DISRES(1,NHKKH1,PPN)
      IF (IPADIS) CALL DISTPA(1)
      IF (IOUDIF.EQ.1) CALL DIADIF(1,0)
* initialization of evaporation-module
      WRITE(6,*)' before NUCLEAR.BIN opened LUNBER= ',LUNBER
* initialization of evaporation-module
	OPEN(UNIT=LUNBER,FILE='NUCLEAR.BIN',STATUS='OLD'
     *     ,FORM='UNFORMATTED')
C	OPEN(UNIT=LUNBER,FILE='NUCLEAR.BIN',STATUS='OLD'
C     *    ,READONLY,FORM='UNFORMATTED')
C    *     ,FORM='UNFORMATTED')
      WRITE(6,*)'NUCLEAR.BIN opened LUNBER= ',LUNBER
        CALL BERTTP
      IF(IEVAP.EQ.1)THEN
C	CALL ZEROIN
C	OPEN(UNIT=LUNBER,FILE='NUCLEAR.BIN',STATUS='OLD'
C    *     ,READONLY,FORM='UNFORMATTED')
C       CALL BERTTP
        CALL INCINI
      ENDIF
C                                         HBOOK HISTOGRAMS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       IF(IHBOOK.EQ.1.AND.CMHIS.EQ.0.D0)THEN
         CALL PLOMB(1,PPNPN,CHAR,XFXFXF,ITIF,IJPROJ)
       ENDIF
       IF(IHBOOK.EQ.1.AND.CMHIS.EQ.1.D0)THEN
         CALL PLOMBC(1,PPNPN,CHAR,XFXFXF,ITIF,IJPROJ)
       ENDIF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     IF(NCOUNT.EQ.1) THEN
*---initialization and test of the random number generator
C     CALL RNDMST(12,34,56,78)
C     CALL RNDMTE(1)
C     ENDIF
*---CONSISTENCY TEST FOR FERMI/PAULI OPTIONS
      IF(LPAULI .AND. (.NOT.FERMP)) THEN
        WRITE(6,'(/2A/A/)')
     +  ' ACTIVATION OF PAULI PRINCIPLE ONLY REASONABLE',
     +  ' IF FERMI ACTIVE', ' LPAULI CHANGED TO .FALSE.'
      LPAULI=.FALSE.
      ENDIF
C     CLOSE (UNIT=5)
C     CLOSE (UNIT=11)
      NCASES=INT(WHAT(1))
                NEVNTS=INT(WHAT(1))
                IF(NEVNTS.LE.0) NEVNTS=1000
                NEVHAD=NEVNTS
C               why 2: NEVNTS passed to DTMAI as argument, NEVHAD in COMMON
C               in DTP NCASES <- NEVNTS
      IF(NCASES.LE.0) NCASES=100
      CALL TIMDAT
*---INITIALIZE GLAUBER THEORY A LA SHMAKOV
C   change JGLAUB in dpmjet25 to JGLAUB=1 (was 2 in dpmjet241)      
      JGLAUB=1
      CALL SHMAKI(IP,IPZ,IT,ITZ,RPROJ,RTARG,PPN)
*
      NN=1
      NP=1
      NT=1
C                                   INITIALIZE COUNTERS
      BNNVV=0.001
      BNNSS=0.001
      BNNSV=0.001
      BNNVS=0.001
      BNNCC=0.001
      BNNDV=0.001
      BNNVD=0.001
      BNNDS=0.001
      BNNSD=0.001
      BNNHH=0.001
      BNNZZ=0.001
      BNNDI=0.001
      BNNZD=0.001
      BNNDZ=0.001
      BPTVV=0.
      BPTSS=0.
      BPTSV=0.
      BPTVS=0.
      BPTCC=0.
      BPTDV=0.
      BPTVD=0.
      BPTDS=0.
      BPTSD=0.
      BPTHH=0.
      BPTZZ=0.
      BPTDI=0.
      BPTZD=0.
      BPTDZ=0.
      BEEVV=0.
      BEESS=0.
      BEESV=0.
      BEEVS=0.
      BEECC=0.
      BEEDV=0.
      BEEVD=0.
      BEEDS=0.
      BEESD=0.
      BEEHH=0.
      BEEZZ=0.
      BEEDI=0.
      BEEZD=0.
      BEEDZ=0.
       BCOUVV=0.
       BCOUSS=0.
       BCOUSV=0.
       BCOUVS=0.
       BCOUZZ=0.
       BCOUHH=0.
       BCOUDS=0.
       BCOUSD=0.
       BCOUDZ=0.
       BCOUZD=0.
       BCOUDI=0.
       BCOUDV=0.
       BCOUVD=0.
       BCOUCC=0.
      CALL TIMDAT
      RETURN
C
C********************************************************************
C               CONTROL CARD: CODEWD =  XSECNUC
C             Calculation of h-A and A-B nuclear X-sects.  
C              
C
C     WHAT(1)= ECMUU
C     WHAT(2)= ECMOO
C     WHAT(3)= NGRITT
C     WHAT(3)= NEVTT
C
C********************************************************************
C
  620 CONTINUE
      ECMUU = WHAT(1)
      ECMOO = WHAT(2)
      NGRITT= WHAT(3)
      NEVTT = WHAT(4)
      WRITE(6,*)'call xsglau'
      CALL XSGLAU(IP,IT,IJPROJ,1)       
      STOP
C
C
C********************************************************************
C               CONTROL CARD: CODEWD =INTERDPM 
C               
C              
C
C     WHAT(1)= INTDPM              DEFAULT: 0
C
C********************************************************************
C
  630 CONTINUE
      IROEH=0
      INTDPM=WHAT(1)
      GO TO 10
C     RETURN
C
C
C********************************************************************
C               CONTROL CARD: CODEWD = 
C               
C              
C
C     WHAT(1)=               DEFAULT: 1000.0
C
C********************************************************************
C
  640 CONTINUE
      RETURN
C
C
C********************************************************************
C               CONTROL CARD: CODEWD = 
C               
C              
C
C     WHAT(1)=               DEFAULT: 1000.0
C
C********************************************************************
C
  650 CONTINUE
      RETURN
C
C
C********************************************************************
C               CONTROL CARD: CODEWD = 
C               
C              
C
C     WHAT(1)=               DEFAULT: 1000.0
C
C********************************************************************
C
  660 CONTINUE
      RETURN
C
C********************************************************************
C               CONTROL CARD: CODEWD = STOP
C               STOPS THE EXECUTION OF THE PROGRAM
C********************************************************************
C
  540 CONTINUE
*
      STOP
C
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
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

