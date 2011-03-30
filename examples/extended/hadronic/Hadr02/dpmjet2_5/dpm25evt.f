C***********************************************************************
      SUBROUTINE DPMEVT(ELABT,IIPROJ,IIP,IIPZ,IIT,IITZ,KKMAT,NHKKH1)
C
C           J.R.   Version 4/97 for dpmjet25
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
*KEEP,HKKEVT.
c     INCLUDE (HKKEVT)
      PARAMETER (NMXHKK=89998)
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK), JMOHKK
     +(2,NMXHKK),JDAHKK(2,NMXHKK), PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK
     +(4,NMXHKK)
C
      CHARACTER*8  ANAME
      COMMON /DPAR/ ANAME(210),AAM(210),GA(210),TAU(210),IICH(210),
     +IIBAR(210),K1(210),K2(210)
C
*KEEP,NUCC.
      COMMON /NUCC/ IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
      COMMON /NUCCC/ JT,JTZ,JP,JPZ,JJPROJ,JBPROJ,JJTARG,JBTARG
C                               from DTUJET93
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
C     COMMON/COLLIS/S,IJPROX,IJTAR,PTTHR,IOPHRD,IJPRLU,IJTALU,PTTHR2
      COMMON/COLLIS/S,IJPROX,IJTAR,PTTHR,PTTHR2,IOPHRD,IJPRLU,IJTALU
*
*  *********************************************************************
*KEEP,NNCMS.
      COMMON /NNCMS/  GAMCM,BGCM,UMO,PCM,EPROJ,PPROJ
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
      CHARACTER*8 BTYPE
      COMMON /PANAME/ BTYPE(30)
      COMMON /STRUFU/ ISTRUM,ISTRUT
C
      COMMON /BUFUEH/ ANNVV, ANNSS, ANNSV, ANNVS, ANNCC, ANNDV,
     * ANNVD, ANNDS, ANNSD, ANNHH, ANNZZ, PTVV, PTSS, PTSV, PTVS, 
     * PTCC, PTDV, PTVD, PTDS, PTSD, PTHH, PTZZ, EEVV, EESS, EESV, 
     * EEVS, EECC, EEDV, EEVD, EEDS, EESD, EEHH, EEZZ, ANNDI, PTDI, 
     * EEDI, ANNZD, ANNDZ, PTZD, PTDZ, EEZD, EEDZ
      COMMON /NCOUCH/ ACOUVV, ACOUSS, ACOUSV, ACOUVS, ACOUZZ, ACOUHH,
     * ACOUDS, ACOUSD, ACOUDZ, ACOUZD, ACOUDI, ACOUDV, ACOUVD, ACOUCC
C
       SAVE ELABT_PREV
C       ON DOUBLE PRECISION UNDERFLOW IGNORE
C       ON REAL UNDERFLOW IGNORE
       DATA ELABT_PREV/-10./      !  lab-energy of previous collision
       DATA NINIT/0/
C
*  *********************************************************************
*                               PROJPAR
	IPROJ=IIPROJ
C                    New 4/97
C       KKMAT=IIP
	IF(IIPROJ.EQ.-1)IPROJ=1
C                    New 4/97
	IF(IPROJ.EQ.12.OR.IPROJ.EQ.19)THEN
	  IPROJ=24
	  IF(RNDM(V).LT.0.5D0)IPROJ=25
	ENDIF
        PROJTY=BTYPE(IPROJ)
        IJPROJ=IPROJ
        IJPROX=IPROJ
	IBPROJ=IIBAR(IPROJ)
	IP=IIP
	IPZ=IIPZ
*                               TARPAR
C       IT=14
C       ITZ=7
	IT=IIT
	ITZ=IITZ
	IJTAR=1
*                               MOMENTUM
      EPN=1000.D0*ELABT
      NNPP=IJPROJ
      AMPROJ=AAM(NNPP)
      PPROJ = SQRT((EPN-AMPROJ)*(EPN+AMPROJ))
      PPN=PPROJ
*                             nucleon-nucleon cms
      EPROJ=EPN
      AMTAR=AAM(1)
      UMO = SQRT(AMPROJ**2 + AMTAR**2 + 2.*AMTAR*EPROJ)
      CMENER=UMO
      GAMCM = (EPROJ+AMTAR)/UMO
      BGCM=PPROJ/UMO
      ECM=UMO
      S=ECM**2
      PCM=GAMCM*PPROJ - BGCM*EPROJ
          IF(ISTRUT.EQ.1)THEN
            PTTHR=2.1+0.15*(LOG10(CMENER/50.))**3
            PTTHR2=PTTHR
          ELSEIF(ISTRUT.EQ.2)THEN
            PTTHR=2.5+0.12*(LOG10(CMENER/50.))**3
            PTTHR2=PTTHR
          ENDIF
C
C     IIT=IT 
C     IITZ=ITZ
C     IIP=IP
C     IIPZ=IPZ
      IIPROJ=IJPROJ
      IITARG = IJTARG
C
C-- INITIALIZE COUNTERS before call to KKINC
C
 765  CONTINUE
      ANNVV=0.001         ! common /BUFUEH/
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
C  
C-- COMMON /NCOUCH/ variables
C
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
C
C-- Pt initialisation each time collision energy varies
C
      IF (ELABT.NE.ELABT_PREV) THEN
CGB	CALL RD2OUT(ISEED,JSEED)
CGB	write(6,*) 'seeds before SAMPPT',ISEED,JSEED
	CALL SAMPPT(0,PT)
CGB	CALL RD2OUT(ISEED,JSEED)
CGB	write(6,*) 'seeds after SAMPPT',ISEED,JSEED
CGB	write(6,*)
      ENDIF
      ELABT_PREV = ELABT
C
      IF(NINIT.LT.10)THEN
      NINIT=NINIT+1
      WRITE(6,*)' DPMEVT  EPN=',EPN,'IIT,IITZ,IIP,IIPZ,IIPROJ,KKMAT',
     *IIT,IITZ,IIP,IIPZ,IIPROJ,KKMAT, ' PTTHR=',PTTHR 
      ENDIF
      CALL KKINC(EPN,IIT,IITZ,IIP,IIPZ,IIPROJ,KKMAT,
     * IITARG,NHKKH1,IREJ)
      IF (IREJ.EQ.1) THEN
        WRITE(6,*)'Exits from KKINC with IREJ=1'
        GO TO 765
      ENDIF
C      WRITE(6,*)'DECHKK called from DPMEVT : NHKKH1=',NHKKH1
C      CALL DECHKK(NHKKH1)
      RETURN
      END
C
