*CMZ :  1.00/00 27/11/91  17.09.38  by  H.-J. M¶hring
*-- Author :
C-------------------------------------------------------------------
C
C                              FILE DMNUC4.FOR
C
C-------------------------------------------------------------------
C
C                               Sample according to Poisson distribution
      FUNCTION NPOISS(AVN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXPAVN=EXP(-AVN)
      K=1
      A=1.
   10 CONTINUE
      U=RNDM(V)
      A=U*A
      IF(A.LT.EXPAVN)THEN
        NPOISS=K-1
        RETURN
      ELSE
        K=K+1
        GO TO 10
      ENDIF
      END
C
C--------------------------------------------------------------------
C------------------------------------------------------------------
C
C                 FILE DPMSEWEW FORTRAN
C
C------------------------------------------------------------------
      SUBROUTINE SEWEW(IOP,NHKKH1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C*** 1=P,  2=N, 3=PI+, 4=PI-, 5=PIO, 6=GAM+HYP, 7=K, 8=ANUC, 9=CHARGED
C*** 10=TOT, 11=TOTHAD
C
*KEEP,HKKEVT.
c     INCLUDE (HKKEVT)
      PARAMETER (TINY= 1.D-5)
      PARAMETER (NMXHKK= 89998)
c     PARAMETER (NMXHKK=25000)
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK), JMOHKK
     +(2,NMXHKK),JDAHKK(2,NMXHKK), PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK
     +(4,NMXHKK)
      COMMON /EXTEVT/ IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),
     &                IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10)
      COMMON /WEWEVT/ IYWEW(NMXHKK),IRWEW(NMXHKK),IIYWEW(2600),
     &               IIRWEW(2600),WEWEW(2600),IDWEW1(2600),
     &               IDWEW2(2600),IDWEW(2600),IBLOCK(NMXHKK),
     &               IDWEWW(2600)
C
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
*KEEP,TNUCFI.
      COMMON /TNUCFI/ AMTFIN,AMFINO,EEXCTA,PTFINX,PTFINY,PTFINZ,PTFINE,
     +                ITFIN,ITZFIN,NNTAR,NPTAR,NHTAR,NQTAR,NPARF
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
      COMMON /HBOO/IHBOOK
C  1 2 Transverse rapidity  densities (rap, p,pi-, trans. dist. (fermi))
C  3 4 Transverse rapidity  densities (rap,ap,pi+, trans. dist. (fermi))
C  5 6 Transverse rapidity  densities (rap,La,aLa, trans. dist. (fermi))
C  7 8 Transverse rapidity  densities (rap,n ,an , trans. dist. (fermi))
C  9   Transverse rapidity  densities (rap,pi0   , trans. dist. (fermi))
      DIMENSION YYLTRA(40,9,41),TRANSA(40),DYYLAM(40),DYALAM(40),
     *XYLTRA(40,9,41)
C
      COMMON /SIGLA/SIGLAU
      DIMENSION INDX(28)

      DATA IPRIOP/0/
      DATA INDX/1,8,10,10,10,10,7,2,7,10,10,7,3,4,5,6,
     *          7,7,7,7,7,7,7,7,7,7,7,7/
C----------------------------------------------------------------------
      DATA KPL /1/
      DATA IEVL /0/
C----------------------------------------------------------------------
      ZERO=0
      ONEONE=1
      TWOTWO=2
      NIP=IP
      AIP=IP
C       Transverse Rapidity density
C          DRTRA in Fermi
      DRTRA=0.5D0
      DYTRA=2.65D0
      RDYTRA=RNDM(V)*DYTRA
      DO 1711 J=1,40
        DYYLAM(J)=1.D-18
        DYALAM(J)=1.D-18
        DO 1711 JJ=1,9
          DO 1711 JJJ=1,41
            YYLTRA(J,JJ,JJJ)=1.D-18
            XYLTRA(J,JJ,JJJ)=-12.00+RDYTRA+J*DYTRA
 1711 CONTINUE
      DO 1712 J=1,40
        AJ=J
        TRANSA(J)=3.14D0*(2.D0*AJ-1.D0)*DRTRA**2
 1712 CONTINUE 
C     COMMON /WEWEVT/ IYWEW(NMXHKK),IRWEW(NMXHKK),IIYWEW(500),
C    &               IIRWEW(500),WEWEW(500),IDWEW1(500),IDWEW2(500),
C    &               IDWEW(500)
      DO 1812 I=1,2600
	IIYWEW(I)=0
	IIRWEW(I)=0
	WEWEW(I)=0.D0
	IDWEW1(I)=0
	IDWEW2(I)=0
	IDWEW(I)=0
	IDWEWW(I)=0
 1812 CONTINUE
      DO 1813 I=NHKKH1,NHKK
	IYWEW(I)=0
	IRWEW(I)=0
	IBLOCK(I)=0
 1813 CONTINUE
C     RETURN
C
C-------------------------------------------------------------------
C
  2   CONTINUE
C
      DO 21 I=NHKKH1,NHKK
        IF (ISTHKK(I).EQ.1)THEN
C                                  j.r. temp for s-s
          PTT=PHKK(1,I)**2+PHKK(2,I)**2+0.00001
          PT=SQRT(PTT)+0.001
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
          PT=SQRT(PTT)+0.001
          AMT=SQRT(PTT+PHKK(5,I)**2)
          YL=LOG((ABS(PHKK(3,I)+SQRT(PHKK(3,I)**2+AMT**2)))/AMT+1.E-18)
          YLLPS=LOG((ABS(PHKK(3,I)+SQRT(PHKK(3,I)**2+PTT)))/PT
     *               +1.E-18)
C         Transverse rapidity Densities
          RTRA=SQRT(VHKK(1,I)**2+VHKK(2,I)**2)*1.D12
          IRTRA=RTRA/DRTRA+1.D0
          IF(IRTRA.LT.1)IRTRA=1
          IF(IRTRA.GT.40)IRTRA=40
          IYTRA=(YL+12.00-RDYTRA)/DYTRA+1.D0
          IF(IYTRA.LT.1)IYTRA=1
          IF(IYTRA.GT.40)IYTRA=40
	  IYWEW(I)=IYTRA
	  IRWEW(I)=IRTRA
          IF(NIX.EQ.4)THEN
            YYLTRA(IYTRA,2,IRTRA)=YYLTRA(IYTRA,2,IRTRA)+
     *                          1.D0/TRANSA(IRTRA)
          ELSEIF(NIX.EQ.3)THEN
            YYLTRA(IYTRA,4,IRTRA)=YYLTRA(IYTRA,4,IRTRA)+
     *                          1.D0/TRANSA(IRTRA)
          ELSEIF(NIX.EQ.1)THEN
            YYLTRA(IYTRA,1,IRTRA)=YYLTRA(IYTRA,1,IRTRA)+
     *                          1.D0/TRANSA(IRTRA)
          ELSEIF(NIX.EQ.8)THEN
            YYLTRA(IYTRA,3,IRTRA)=YYLTRA(IYTRA,3,IRTRA)+
     *                          1.D0/TRANSA(IRTRA)
          ELSEIF(NIX.EQ.13)THEN
            YYLTRA(IYTRA,5,IRTRA)=YYLTRA(IYTRA,5,IRTRA)+
     *                          1.D0/TRANSA(IRTRA)
          ELSEIF(NIX.EQ.12)THEN
            YYLTRA(IYTRA,8,IRTRA)=YYLTRA(IYTRA,8,IRTRA)+
     *                          1.D0/TRANSA(IRTRA)
          ELSEIF(NIX.EQ.2)THEN
            YYLTRA(IYTRA,7,IRTRA)=YYLTRA(IYTRA,7,IRTRA)+
     *                          1.D0/TRANSA(IRTRA)
          ELSEIF(NIX.EQ.15)THEN
            YYLTRA(IYTRA,6,IRTRA)=YYLTRA(IYTRA,6,IRTRA)+
     *                          1.D0/TRANSA(IRTRA)
          ELSEIF(NIX.EQ.23)THEN
            YYLTRA(IYTRA,9,IRTRA)=YYLTRA(IYTRA,9,IRTRA)+
     *                          1.D0/TRANSA(IRTRA)
          ENDIF
        ENDIF
   21 CONTINUE
C     RETURN
C
C--------------------------------------------------------------------
C
  3   CONTINUE
C          Normalization and Printing
C------------------------------------------------------------------
C------------------------------------------------------------------
C      Transverse Rapidity Densities
      DO 1133 I=1,40
        DO 1133 II=1,9
          DO 1133 III=1,40
            YYLTRA(I,II,III)=YYLTRA(I,II,III)/(DYTRA)
            YYLTRA(I,II,41) =YYLTRA(I,II,41) +
     *                       YYLTRA(I,II,III)*TRANSA(III)
 1133 CONTINUE
C      transverse Rapidity Densities
      IF(IPEV.GE.2)THEN
      WRITE(6,'(A)')' Transverse Rapidity Density of proton '
        WRITE(6,37)XYLTRA(1,1,1),(TRANSA(I),I=1,10)
   37   FORMAT(F10.2,10E11.3)
      DO 1136 J=1,40
        WRITE(6,37)XYLTRA(J,1,1),(YYLTRA(J,1,I),I=1,9),YYLTRA(J,1,41)
 1136 CONTINUE
      WRITE(6,'(A)')' Transverse Rapidity Density of pi- '
        WRITE(6,37)XYLTRA(1,1,1),(TRANSA(I),I=1,10)
      DO 1137 J=1,40
        WRITE(6,37)XYLTRA(J,1,1),(YYLTRA(J,2,I),I=1,9),YYLTRA(J,2,41)
 1137 CONTINUE
      WRITE(6,'(A)')' Transverse Rapidity Density of aproton '
        WRITE(6,37)XYLTRA(1,1,1),(TRANSA(I),I=1,10)
      DO 2236 J=1,40
        WRITE(6,37)XYLTRA(J,1,1),(YYLTRA(J,3,I),I=1,9),YYLTRA(J,3,41)
 2236 CONTINUE
      WRITE(6,'(A)')' Transverse Rapidity Density of pi+ '
        WRITE(6,37)XYLTRA(1,1,1),(TRANSA(I),I=1,10)
      DO 2237 J=1,40
        WRITE(6,37)XYLTRA(J,1,1),(YYLTRA(J,4,I),I=1,9),YYLTRA(J,4,41)
 2237 CONTINUE
      WRITE(6,'(A)')' Transverse Rapidity Density of pi0 '
        WRITE(6,37)XYLTRA(1,1,1),(TRANSA(I),I=1,10)
      DO 2238 J=1,40
        WRITE(6,37)XYLTRA(J,1,1),(YYLTRA(J,9,I),I=1,9),YYLTRA(J,9,41)
 2238 CONTINUE
      ENDIF
      IWEW=0
        DO 1134 III=1,40
      DO 1134 I=1,40
C                                   pi- + p --> Lam + K0
	  DELLLA=
     *	  (YYLTRA(I,1,III)-YYLTRA(I,5,III))*YYLTRA(I,2,III)*
     *              TRANSA(III)*0.53D0 
	  DELLLA=DELLLA*DYTRA
C                                   pi+ + n --> Lam + K+
	  DELLLB=
     *	  (YYLTRA(I,7,III)-YYLTRA(I,5,III))*YYLTRA(I,4,III)*
     *              TRANSA(III)*0.53D0 
	  DELLLB=DELLLB*DYTRA
C                                   pi0 + n --> Lam + K0
	  DELLLC=
     *	  (YYLTRA(I,7,III)-YYLTRA(I,5,III))*YYLTRA(I,9,III)*
     *              TRANSA(III)*0.53D0 
	  DELLLC=DELLLC*DYTRA
C                                   pi0 + p --> Lam + K+
	  DELLLD=
     *	  (YYLTRA(I,1,III)-YYLTRA(I,5,III))*YYLTRA(I,9,III)*
     *              TRANSA(III)*0.53D0 
	  DELLLD=DELLLD*DYTRA
C                              3.11.95 Add 0.15* to reduce Lambdas
C    *    0.15*
	  DELLLO=
     *	  (YYLTRA(I,1,III)+YYLTRA(I,7,III)-1.6*YYLTRA(I,5,III))*
     *              TRANSA(III)/3.2
	  DELLLO=DELLLO*DYTRA
	  IF(DELLLA.GT.DELLLO)DELLLA=DELLLO
	  IF(DELLLB.GT.DELLLO)DELLLB=DELLLO
	  IF(DELLLC.GT.DELLLO)DELLLC=DELLLO
	  IF(DELLLD.GT.DELLLO)DELLLD=DELLLO
	  IF(DELLLO.LE.0.D0)DELLLA=0.D0
	  IF(DELLLO.LE.0.D0)DELLLB=0.D0
	  IF(DELLLO.LE.0.D0)DELLLC=0.D0
	  IF(DELLLO.LE.0.D0)DELLLD=0.D0
	  IF(DELLLA.LE.0.D0)DELLLA=0.D0
	  IF(DELLLB.LE.0.D0)DELLLB=0.D0
	  IF(DELLLC.LE.0.D0)DELLLC=0.D0
	  IF(DELLLD.LE.0.D0)DELLLD=0.D0
C                                   pi- + p --> Lam + K0
	  IF(DELLLA.GT.TINY)THEN
	    IDELLL=DELLLA
	    TELLLA=DELLLA-IDELLL
	    IF(RNDM(V).LE.TELLLA)IDELLL=IDELLL+1
	    DELLLA=IDELLL
	    IF (DELLLA.LT.TINY)GO TO 2233
	    DO 556 KKK=1,IDELLL
	      IWEW=IWEW+1
	      WEWEW(IWEW)=DELLLA
	      IDWEW(IWEW)=17
	      IDWEWW(IWEW)=24
	      IDWEW1(IWEW)=1
	      IDWEW2(IWEW)=14
	      IF(IPEV.GE.2)THEN
      	      WRITE (6,'(A,3F10.3,7I5)')'  LAMBDA ',
     *        DELLLA,YYLTRA(I,1,III),YYLTRA(I,2,III),I,III, 
     *        IDWEW(IWEW),IDWEWW(IWEW),IDWEW1(IWEW),IDWEW2(IWEW),
     *        IWEW
	      ENDIF
              DYYLAM(I)=DYYLAM(I)+  DELLLA/DYTRA
  	      IIYWEW(IWEW)=I
 	      IIRWEW(IWEW)=III
  556       CONTINUE
C 	    IYWEW(IWEW)=I
C	    IRWEW(IWEW)=III
 2233       CONTINUE
	  ENDIF
C                                   pi+ + n --> Lam + K+
	  IF(DELLLB.GT.TINY)THEN
	    IDELLL=DELLLB
	    TELLLB=DELLLB-IDELLL
	    IF(RNDM(V).LE.TELLLB)IDELLL=IDELLL+1
	    DELLLB=IDELLL
	    IF (DELLLB.LT.TINY)GO TO 2243
	    DO 557 KKK=1,IDELLL
	      IWEW=IWEW+1
	      WEWEW(IWEW)=DELLLB
	      IDWEW(IWEW)=17
	      IDWEWW(IWEW)=15
	      IDWEW1(IWEW)=8
	      IDWEW2(IWEW)=13
	      IF(IPEV.GE.2)THEN
      	      WRITE (6,'(A,3F10.3,7I5)')'  LAMBDA ',
     *        DELLLB,YYLTRA(I,1,III),YYLTRA(I,2,III),I,III, 
     *        IDWEW(IWEW),IDWEWW(IWEW),IDWEW1(IWEW),IDWEW2(IWEW),
     *        IWEW
	      ENDIF
              DYYLAM(I)=DYYLAM(I)+  DELLLB/DYTRA
  	      IIYWEW(IWEW)=I
 	      IIRWEW(IWEW)=III
  557       CONTINUE
C 	    IYWEW(IWEW)=I
C	    IRWEW(IWEW)=III
 2243       CONTINUE
	  ENDIF
C                                   pi0 + n --> Lam + K0
	  IF(DELLLC.GT.TINY)THEN
	    IDELLL=DELLLC
	    TELLLC=DELLLC-IDELLL
	    IF(RNDM(V).LE.TELLLC)IDELLL=IDELLL+1
	    DELLLC=IDELLL
	    IF (DELLLC.LT.TINY)GO TO 2244
	    DO 558 KKK=1,IDELLL
	      IWEW=IWEW+1
	      WEWEW(IWEW)=DELLLC
	      IDWEW(IWEW)=17
	      IDWEWW(IWEW)=24
	      IDWEW1(IWEW)=8
	      IDWEW2(IWEW)=23
	      IF(IPEV.GE.2)THEN
      	      WRITE (6,'(A,3F10.3,7I5)')'  LAMBDA ',
     *        DELLLC,YYLTRA(I,1,III),YYLTRA(I,2,III),I,III, 
     *        IDWEW(IWEW),IDWEWW(IWEW),IDWEW1(IWEW),IDWEW2(IWEW),
     *        IWEW
	      ENDIF
              DYYLAM(I)=DYYLAM(I)+  DELLLC/DYTRA
  	      IIYWEW(IWEW)=I
 	      IIRWEW(IWEW)=III
  558       CONTINUE
C 	    IYWEW(IWEW)=I
C	    IRWEW(IWEW)=III
 2244       CONTINUE
	  ENDIF
C                                   pi0 + p --> Lam + K+
	  IF(DELLLD.GT.TINY)THEN
	    IDELLL=DELLLD
	    TELLLD=DELLLD-IDELLL
	    IF(RNDM(V).LE.TELLLD)IDELLL=IDELLL+1
	    DELLLD=IDELLL
	    IF (DELLLD.LT.TINY)GO TO 2245
	    DO 559 KKK=1,IDELLL
	      IWEW=IWEW+1
	      WEWEW(IWEW)=DELLLD
	      IDWEW(IWEW)=17
	      IDWEWW(IWEW)=15
	      IDWEW1(IWEW)=1
	      IDWEW2(IWEW)=23
	      IF(IPEV.GE.2)THEN
      	      WRITE (6,'(A,3F10.3,7I5)')'  LAMBDA ',
     *        DELLLD,YYLTRA(I,1,III),YYLTRA(I,2,III),I,III, 
     *        IDWEW(IWEW),IDWEWW(IWEW),IDWEW1(IWEW),IDWEW2(IWEW),
     *        IWEW
	      ENDIF
              DYYLAM(I)=DYYLAM(I)+  DELLLD/DYTRA
  	      IIYWEW(IWEW)=I
 	      IIRWEW(IWEW)=III
  559       CONTINUE
C 	    IYWEW(IWEW)=I
C	    IRWEW(IWEW)=III
 2245       CONTINUE
	  ENDIF
C     COMMON /WEWEVT/ IYWEW(NMXHKK),IRWEW(NMXHKK),IIYWEW(500),
C    &               IIRWEW(500),WEWEW(500),IDWEW1(500),IDWEW2(500),
C    &               IDWEW(500)
C                          pi+ + pb --> aLam + aK0
	  DELALA=
     *	  (YYLTRA(I,3,III)-YYLTRA(I,6,III))*YYLTRA(I,4,III)*
     *              TRANSA(III)*0.53D0 
	  DELALA=DELALA*DYTRA
C                          pi- + nb --> aLam + aK-
	  DELALB=
     *	  (YYLTRA(I,8,III)-YYLTRA(I,6,III))*YYLTRA(I,2,III)*
     *              TRANSA(III)*0.53D0 
	  DELALB=DELALB*DYTRA
C                              3.11.95 Add 0.15* to reduce Lambdas
C    *    0.15*
	  DELALO=
     *	  (YYLTRA(I,3,III)+YYLTRA(I,8,III)-1.6*YYLTRA(I,6,III))*
     *              TRANSA(III)/3.2 
	  DELALO=DELALO*DYTRA
	  IF(DELALA.GT.DELALO)DELALA=DELALO
	  IF(DELALB.GT.DELALO)DELALB=DELALO
	  IF(DELALO.LE.0.D0)DELALA=0.D0
	  IF(DELALO.LE.0.D0)DELALB=0.D0
	  IF(DELALA.LE.0.D0)DELALA=0.D0
	  IF(DELALB.LE.0.D0)DELALB=0.D0
C                          pi+ + pb --> aLam + aK0
	  IF(DELALA.GT.TINY)THEN
	    IDELAL=DELALA
	    TELALA=DELALA-IDELAL
	    IF(RNDM(V).LE.TELALA)IDELAL=IDELAL+1
	    DELALA=IDELAL
	    IF (DELALA.LT.TINY)GO TO 2234
	    DO 554 KKK=1,IDELAL
	      IWEW=IWEW+1
	      WEWEW(IWEW)=DELALA
	      IDWEW(IWEW)=18
	      IDWEWW(IWEW)=25
	      IDWEW1(IWEW)=2
	      IDWEW2(IWEW)=13
	      IF(IPEV.GE.2)THEN
	      WRITE (6,'(A,3F10.3,7I5)')' ALAMBDA ',
     *        DELALA,YYLTRA(I,3,III),YYLTRA(I,6,III),I,III,             
     *        IDWEW(IWEW),IDWEWW(IWEW),IDWEW1(IWEW),IDWEW2(IWEW),
     *        IWEW
	      ENDIF
              DYALAM(I)=DYALAM(I)+   DELALA/DYTRA
  	      IIYWEW(IWEW)=I
 	      IIRWEW(IWEW)=III
  554       CONTINUE
C 	    IYWEW(IWEW)=I
C	    IRWEW(IWEW)=III
 2234       CONTINUE
	  ENDIF
C                          pi- + nb --> aLam + aK-
	  IF(DELALB.GT.TINY)THEN
	    IDELAL=DELALB
	    TELALB=DELALB-IDELAL
	    IF(RNDM(V).LE.TELALB)IDELAL=IDELAL+1
	    DELALB=IDELAL
	    IF (DELALB.LT.TINY)GO TO 2247
	    DO 555 KKK=1,IDELAL
	      IWEW=IWEW+1
	      WEWEW(IWEW)=DELALB
	      IDWEW(IWEW)=18
	      IDWEWW(IWEW)=16
	      IDWEW1(IWEW)=9
	      IDWEW2(IWEW)=14
	      IF(IPEV.GE.2)THEN
	      WRITE (6,'(A,3F10.3,7I5)')' ALAMBDA ',
     *        DELALB,YYLTRA(I,3,III),YYLTRA(I,6,III),I,III,             
     *        IDWEW(IWEW),IDWEWW(IWEW),IDWEW1(IWEW),IDWEW2(IWEW),
     *        IWEW
	      ENDIF
              DYALAM(I)=DYALAM(I)+   DELALB/DYTRA
  	      IIYWEW(IWEW)=I
 	      IIRWEW(IWEW)=III
  555       CONTINUE
C 	    IYWEW(IWEW)=I
C	    IRWEW(IWEW)=III
 2247       CONTINUE
	  ENDIF
 1134 CONTINUE
      IWEWMA=IWEW
      DDLAMN=0.D0
      DALAMN=0.D0
      DO 1139 J=1,40
        DDLAMN=DDLAMN+DYYLAM(J)*DYTRA
        DALAMN=DALAMN+DYALAM(J)*DYTRA
 1139 CONTINUE
	      IF(IPEV.GE.2)THEN
      WRITE(6,'(A,I10,F10.3)')' IWEWMA,Extra Rap. Den. of Lam ,DLAM= '
     *                                            ,IWEWMA ,DDLAMN
      WRITE(6,'(A,I10,F10.3)')' IWEWMA,Extra Rap. Den. ofALam ,DLAM= '
     *                                           ,IWEWMA ,DALAMN
        WRITE(6,37)XYLTRA(1,1,1),(TRANSA(I),I=1,10)
      DO 1138 J=1,40
        WRITE(6,37)XYLTRA(J,1,1),DYYLAM(J),DYALAM(J)
 1138 CONTINUE
	      ENDIF
      I=0
C                               Reduce IWEWMA forPb-Pb runs
      IWEWMA=IWEWMA
      ILOMA=0
      ILOML=0
      ITOMA=0
      ITOML=0
      DO 20 II=1,IWEWMA
	ID1=MPDGHA(IDWEW1(II))
	ID2=MPDGHA(IDWEW2(II))
	ID3=MPDGHA(IDWEW(II))
	ID4=MPDGHA(IDWEWW(II))
C       WRITE(6,'(5I5)')II,ID1,ID2,ID3,ID4
	IHKK1=0
	IHKK2=0
        KK1=0
        KK2=0
        DO 211 JJ=NHKKH1,NHKK
	  IF(IBLOCK(JJ).EQ.0)THEN
            IF((IDHKK(JJ).EQ.ID1).AND.(IYWEW(JJ).EQ.IIYWEW(II))
     *      .AND.(IRWEW(JJ).EQ.IIRWEW(II)).AND.(KK1.EQ.0))THEN
	      IHKK1=JJ
              JJ1=JJ
              KK1=1
            ENDIF
	    IF((IDHKK(JJ).EQ.ID2).AND.(IYWEW(JJ).EQ.IIYWEW(II))
     *      .AND.(IRWEW(JJ).EQ.IIRWEW(II)).AND.(KK2.EQ.0))THEN
	      IHKK2=JJ
              JJ2=JJ
              KK2=1
            ENDIF
	    IF(IHKK1.EQ.0)GO TO 211
	    IF(IHKK2.EQ.0)GO TO 211
	    IF(JMOHKK(1,IHKK1).EQ.JMOHKK(1,IHKK2))THEN
	      IHKK2=0
	      JJ2=0
	      KK2=0
	      GO TO 211
	    ENDIF
            CALL PINKLA(IDWEW1(II),IDWEW2(II),IDWEW(II),IDWEWW(II),    
     &      IHKK1,IHKK2,ILOML,ILOMA,ITOML,ITOMA,IREJ)
            IF(IPEV.GE.2)
     *      WRITE(6,'(A,5I10)')' IWEWMA,ILOML,ILOMA,ITOML,ITOMA  '
     *             ,IWEWMA ,ILOML,ILOMA,ITOML,ITOMA
            IF(IREJ.EQ.0)THEN
	      IF(IPEV.GE.2)THEN
	        WRITE(6,'(A,3I7)')' Extra Positions ',IHKK1,IHKK2,II
              ENDIF
              IBLOCK(JJ1)=1
              IBLOCK(JJ2)=1
              KK1=0
              KK2=0
              JJ1=0
              JJ2=0
	      IHKK1=0
	      IHKK2=0
              GO TO 20
            ENDIF
	    IF(JJ1.GE.1)IBLOCK(JJ1)=1
 	    IHKK1=0
	    JJ1=0
	    KK1=0
	    IF(JJ2.GE.1)IBLOCK(JJ2)=1
	    IHKK2=0
	    JJ2=0
	    KK2=0
	  ENDIF
  211   CONTINUE
   20 CONTINUE
      IF(IPEV.GE.2)
     *WRITE(6,'(A,5I10)')' IWEWMA,ILOML,ILOMA,ITOML,ITOMA  '
     *             ,IWEWMA ,ILOML,ILOMA,ITOML,ITOMA
      
      RETURN
      END
      SUBROUTINE PINKLA(I1,I2,I3,I4,IHKK1,IHKK2,ILOML,ILOMA,
     *ITOML,ITOMA,IREJ)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C--------------------------------------------------------------
C
C                 I1 + I2 ---> I3 + I4 final state interaction
C
C--------------------------------------------------------------
*KEEP,HKKEVT.
c     INCLUDE (HKKEVT)
      PARAMETER (TINY= 1.D-5)
      PARAMETER (NMXHKK= 89998)
c     PARAMETER (NMXHKK=25000)
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK), JMOHKK
     +(2,NMXHKK),JDAHKK(2,NMXHKK), PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK
     +(4,NMXHKK)
      COMMON /EXTEVT/ IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),
     &                IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10)
      COMMON /WEWEVT/ IYWEW(NMXHKK),IRWEW(NMXHKK),IIYWEW(2600),
     &               IIRWEW(2600),WEWEW(2600),IDWEW1(2600),
     &               IDWEW2(2600),IDWEW(2600),IBLOCK(NMXHKK),
     &               IDWEWW(2600)
C
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
      COMMON /DPRIN/IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
C------------------
      IF(I3.EQ.17)ITOML=ITOML+1
      IF(I3.EQ.18)ITOMA=ITOMA+1
C                            Check Inv Mass of two particle system
      AMA12 = SQRT((PHKK(4,IHKK1)+PHKK(4,IHKK2))**2
     +            -(PHKK(1,IHKK1)+PHKK(1,IHKK2))**2
     +            -(PHKK(2,IHKK1)+PHKK(2,IHKK2))**2
     +            -(PHKK(3,IHKK1)+PHKK(3,IHKK2))**2)
      AMAMIN = AAM(I3)+AAM(I4)
      IREJ=1
      IF(AMA12.GE.AMAMIN+0.05D0)THEN
      IF(IPEV.GE.2)THEN
        WRITE(6,'(A,2F10.3)')' AMAMIN,AMA12 ',AMAMIN,AMA12
	WRITE(6,'(A,I5,A,I5,A,I5,A,I5)')' Reaction ',
     + I1,' + ',I2,' ---> ',I3,' + ',I4
      ENDIF
C                      Lorentz parameters of IHKK1,IHKK2 system
      E12  = PHKK(4,IHKK1)+PHKK(4,IHKK2)
      PX12 = PHKK(1,IHKK1)+PHKK(1,IHKK2)
      PY12 = PHKK(2,IHKK1)+PHKK(2,IHKK2)
      PZ12 = PHKK(3,IHKK1)+PHKK(3,IHKK2)
      G12   = E12/AMA12
      BGX12 = PX12/AMA12
      BGY12 = PY12/AMA12
      BGZ12 = PZ12/AMA12
      IF(IPEV.GE.2)
     +WRITE(6,'(A,4E12.3)')' G12,BGX12,BGY12,BGZ12 ',
     +G12,BGX12,BGY12,BGZ12
C                      Decay AMA12 into I3 + I4
      E3 = (AMA12**2-AAM(I4)**2+AAM(I3)**2)/(2.*AMA12)
      E4 = (AMA12**2-AAM(I3)**2+AAM(I4)**2)/(2.*AMA12)
      P3 = SQRT(E3**2-AAM(I3)**2)
      P4 = SQRT(E4**2-AAM(I4)**2)
      IF(IPEV.GE.2)
     +WRITE(6,'(A,4E12.3)')' E3,E4,P3,P4 ',
     +E3,E4,P3,P4
C                      Uniform rotation
      CALL DPOLI(POLC,POLS)
      CALL DSFECF(SFE,CFE)
      IF(PHKK(3,IHKK1).GE.0.D0)THEN
        POLS=0.D0
        POLC=1.D0
      ELSEIF(PHKK(3,IHKK1).LE.0.D0)THEN
        POLS=0.D0
        POLC=-1.D0
      ENDIF
C     WRITE(6,*)' POLS,POLC ',POLS,POLC
      CX=POLS*CFE
      CY=POLS*SFE
      CZ=POLC
      PX3=CX*P3
      PY3=CY*P3
      PZ3=CZ*P3
      PX4=-CX*P4
      PY4=-CY*P4
      PZ4=-CZ*P4
      IF(IPEV.GE.2)
C    +WRITE(6,'(A,8E12.3)')' 4 mom cms 3   4 ',
C    +PX3,PY3,PZ3,E3,PX4,PY4,PZ4,E4
     +WRITE(6,'(A,4E12.3)')' 4 mom cms 3   4 ',
     +PZ3,E3,PZ4,E4
C                        Lorentz transformation back
      CALL DALTRA(G12,BGX12,BGY12,BGZ12,PX3,PY3,PZ3,E3,PC3,PCX3,
     +                                   PCY3,PCZ3,EC3)
      CALL DALTRA(G12,BGX12,BGY12,BGZ12,PX4,PY4,PZ4,E4,PC4,PCX4,
     +                                   PCY4,PCZ4,EC4)
      IF(IPEV.GE.2)
C    +WRITE(6,'(A,8E12.3)')' 4 mom 1   2 ',
C    +PHKK(1,IHKK1),PHKK(2,IHKK1),PHKK(3,IHKK1),PHKK(4,IHKK1),
C    +PHKK(1,IHKK2),PHKK(2,IHKK2),PHKK(3,IHKK2),PHKK(4,IHKK2) 
     +WRITE(6,'(A,4E12.3)')' 4 mom 1   2 ',
     +PHKK(3,IHKK1),PHKK(4,IHKK1),
     +PHKK(3,IHKK2),PHKK(4,IHKK2) 
      IF(IPEV.GE.2)
C    +WRITE(6,'(A,8E12.3)')' 4 mom 3   4 ',
C    +PCX3,PCY3,PCZ3,EC3,PCX4,PCY4,PCZ4,EC4
     +WRITE(6,'(A,4E12.3)')' 4 mom 3   4 ',
     +PCZ3,EC3,PCZ4,EC4
C---------------------------------------------------      
	IREJ=0
C	WRITE(6,*)' I1+I2-->I3+I4 ',I1,I2,I3,I4
C                 Put particles I3 and I4 into COMMON Block
      ISTHKK(IHKK1)=512
      ISTHKK(IHKK2)=512
      NHKK=NHKK+1
      ISTHKK(NHKK)=1
      IF((I3.EQ.24).OR.(I3.EQ.25))THEN
	I3=12
        IF(RNDM(VV).LE.0.5D0)I3=19
      ENDIF
      IDHKK(NHKK)=MPDGHA(I3)
      JMOHKK(1,NHKK)=IHKK1
      JMOHKK(2,NHKK)=IHKK2
      JDAHKK(1,IHKK1)=NHKK
      JDAHKK(1,IHKK2)=NHKK
      PHKK(1,NHKK)=PCX3
      PHKK(2,NHKK)=PCY3
      PHKK(3,NHKK)=PCZ3
      PHKK(4,NHKK)=EC3
      PHKK(5,NHKK)=AAM(I3)
      DO 11 IIII=1,4
      VHKK(IIII,NHKK)=VHKK(IIII,IHKK1)
      WHKK(IIII,NHKK)=WHKK(IIII,IHKK1)
   11 CONTINUE
      NHKK=NHKK+1
      ISTHKK(NHKK)=1
      IF((I4.EQ.24).OR.(I4.EQ.25))THEN
	I4=12
        IF(RNDM(VV).LE.0.5D0)I4=19
      ENDIF
      IDHKK(NHKK)=MPDGHA(I4)
      JMOHKK(1,NHKK)=IHKK1
      JMOHKK(2,NHKK)=IHKK2
      JDAHKK(2,IHKK1)=NHKK
      JDAHKK(2,IHKK2)=NHKK
      PHKK(1,NHKK)=PCX4
      PHKK(2,NHKK)=PCY4
      PHKK(3,NHKK)=PCZ4
      PHKK(4,NHKK)=EC4
      PHKK(5,NHKK)=AAM(I4)
      DO 12 IIII=1,4
      VHKK(IIII,NHKK)=VHKK(IIII,IHKK2)
      WHKK(IIII,NHKK)=WHKK(IIII,IHKK2)
   12 CONTINUE
      ELSE
        IF(I3.EQ.17)ILOML=ILOML+1
        IF(I3.EQ.18)ILOMA=ILOMA+1
      ENDIF
      RETURN
      END

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HADJCK(NHAD,AMCH,PPR,PTA,GAM,BGX,BGY,BGZ, IFB1,IFB2,
     +IFB3,IFB4,I1,I2,NOBAM,NNCH,NORIG,IREJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C      HADJCK Fragmentation of vv,vs and sv chains CK method 
C                      j.r.6/96
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
      CHARACTER*8 ANF,ANFF
      PARAMETER (NFIMAX=249)
      COMMON /DFINPA/ ANF(NFIMAX),PXF(NFIMAX),PYF(NFIMAX),PZF(NFIMAX),
     +HEF(NFIMAX),AMF(NFIMAX), ICHF(NFIMAX),IBARF(NFIMAX),NREF(NFIMAX)
      COMMON /DFINPZ/IORMO(NFIMAX),IDAUG1(NFIMAX),IDAUG2(NFIMAX),
     * ISTATH(NFIMAX)
      DIMENSION ANFF(NFIMAX),PXFF(NFIMAX),PYFF(NFIMAX),PZFF(NFIMAX),
     +HEFF(NFIMAX),AMFF(NFIMAX), 
     *ICHFF(NFIMAX),IBARFF(NFIMAX),NREFF(NFIMAX)
        CHARACTER*80 TITLE
        CHARACTER*8 PROJTY,TARGTY
C       COMMON /USER/TITLE,PROJTY,TARGTY,CMENER,ISTRUF
C    & ,ISINGD,IDUBLD,SDFRAC,PTLAR
      COMMON /USER1/TITLE,PROJTY,TARGTY
      COMMON /USER2/CMENER,SDFRAC,PTLAR,ISTRUF,ISINGD,IDUBLD
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
      DIMENSION PPR(4),PTA(4)
*KEEP,XSEADI.
      COMMON /XSEADI/ XSEACU,UNON,UNOM,UNOSEA,CVQ,CDQ,CSEA,SSMIMA,
     +SSMIMQ,VVMTHR
      DIMENSION PPP1(4),PPP2(4),PPP3(4),PPP4(4)
C----------------------------------------------------------------------
C
C	WRITE(6,'(A,3I5)')' Jet0 IFB1,IFB2,IFB3 ',IFB1,IFB2,IFB3
      IREJ=0
      IF(IPCO.GE.3) THEN
        WRITE(6,'(4E12.4)') GAM,BGX,BGY,BGZ,PPR
        WRITE(6,1000) NHAD,AMCH,PTA, IFB1,IFB2,IFB3,IFB4,I1,I2,NOBAM,
     +  NNCH,NORIG
 1000 FORMAT(10X,I10,5F10.3/10X,9I10)
      ENDIF
      IF(ABS(NNCH).EQ.99) THEN
        NHAD=0
C       IPCO=0
        RETURN
      ENDIF
C=================================================================
C=================================================================
C              Different options:
      IF(NOBAM.EQ.6)THEN
C              NOBAM=6:  Diquark--Quark Jet
C
C              Work out approximate x-Values of diquark and quark
	XDIQ=PPR(4)*2.D0/CMENER
	XQUA=PTA(4)*2.D0/CMENER
C	WRITE(6,'(A,2E12.4)')' HADJCK:XDIQ,XQUA ',XDIQ,XQUA
C              Select valence quark x for one of the diquark quarks
C              But only if XDIQ > 1.5*XQUA  drop this condition!
C	IF(XDIQ.LE.1.5D0*XQUA)THEN
C	  IREJ=1
C	  RETURN
C       ENDIF
C              Select x values of sea quark pair
	ICOU=0
 1234   CONTINUE
	ICOU=ICOU+1
	IF(ICOU.GE.100)THEN
	  IREJ=1
	  RETURN
	ENDIF
	CALL XSEAPA(CMENER,XQUA/2.D0,ISQ,ISAQ,XSQ,XSAQ,IREJ)
	IF(XSAQ.GE.XQUA/2.D0)GO TO 1234
C       WRITE(6,'(A,2E12.4)')' HADJCK:XSQ,XSAQ ',XSQ,XSAQ
*  projectile valence quarks xpvqi
        XVTHRO=CVQ/CMENER
	IVTHR=0
 3465   CONTINUE
	IF(IVTHR.EQ.50)THEN
	  IREJ=1
      WRITE(6,*)' HADJSE 3465 reject IVTHR 50'
	  RETURN
        ENDIF
	IVTHR=IVTHR+1
        XVTHR=XVTHRO/(51-IVTHR)
        UNOPRV=UNON
   80   CONTINUE
        IF(XVTHR.GT.0.05)THEN
          IF(XVTHR.GT.0.66D0*XDIQ)THEN
            IREJ=1
	    RETURN
	  ENDIF
          XPVQI=BETREJ(0.5,UNOPRV,XVTHR,0.66D0*XDIQ)
        ELSE
          N90=0
   90     CONTINUE
          N90=N90+1
	  IF(N90.GE.20)THEN
	    IREJ=1
	    RETURN
	  ENDIF
          XPVQI=DBETAR(0.5,UNOPRV)
          IF ((XPVQI.LT.XVTHR).OR.(XPVQI.GT.0.66D0*XDIQ))
     *                                                          GOTO 90
        ENDIF
	XDIQQ=XDIQ-XPVQI
C	WRITE(6,'(A,2E12.4)')' HADJCK:XDIQQ,XPVQI ',XDIQQ,XPVQI
C
C               We form two strings(1-2) q-aq and (3-4) qq-q
C            1  q with XDIQ-XPVQI-XSQ
C            1  q with XDIQ-XPVQI   (j.r.5.5.96)
C            2  aq with XSAQ
C            3  qq with XPVQI+XSQ
C            3  qq with XPVQI       (j.r.5.5.96)
C            4  q with  XQUA-XSAQ
C            with 4-momenta PPP1,PPP2,PPP3,PPP4
      DO 2345 I=1,4
C       PPP1(I)=PPR(I)*(XDIQ-XPVQI-XSQ)/XDIQ
        PPP1(I)=PPR(I)*(XDIQ-XPVQI)/XDIQ
C	PPP3(I)=PPR(I)*(XPVQI+XSQ)/XDIQ
	PPP3(I)=PPR(I)*(XPVQI)/XDIQ
	PPP2(I)=PTA(I)*XSAQ/XQUA
	PPP4(I)=PTA(I)*(XQUA-XSAQ)/XQUA
 2345   CONTINUE
 2346   FORMAT(A,4E12.4)
C	WRITE(6,2346)' PPR ',PPR
C	WRITE(6,2346)' PPP1 ',PPP1
C	WRITE(6,2346)' PPP3 ',PPP3
C	WRITE(6,2346)' PTA ',PTA
C	WRITE(6,2346)' PPP2 ',PPP2
C	WRITE(6,2346)' PPP4 ',PPP4
C               Invariant Masses of chains
C               ORIGINAL CHAIN
	AMCHOR=SQRT((PPR(4)+PTA(4))**2-(PPR(3)+PTA(3))**2
     *             -(PPR(2)+PTA(2))**2-(PPR(1)+PTA(1))**2) 
	AMCHN1=SQRT((PPP1(4)+PPP2(4))**2-(PPP1(3)+PPP2(3))**2
     *             -(PPP1(2)+PPP2(2))**2-(PPP1(1)+PPP2(1))**2) 
C                    Chain 1 is q-aq restrict mass >1.5 GeV
	CHAMAL=1.D0
	IF(IFB1.GE.3.OR.ISAQ.GE.9)CHAMAL=1.5D0
	IF(AMCHN1.LE.CHAMAL)THEN
C	  IREJ=1
C	  RETURN
          GO TO 3465
	ENDIF
	AMCHN2=SQRT((PPP3(4)+PPP4(4))**2-(PPP3(3)+PPP4(3))**2
     *             -(PPP3(2)+PPP4(2))**2-(PPP3(1)+PPP4(1))**2) 
C                    Chain 2 is qq-q restrict mass >2.5 GeV
	CHAMAL=1.8D0
	IF(IFB2.GE.3.OR.IFB3.GE.3.OR.ISQ.GE.3)CHAMAL=2.5D0
	IF(AMCHN2.LE.CHAMAL)THEN
C	  IREJ=1
C 	  RETURN
	  GO TO 3465
	ENDIF
	PXCHK=PPR(1)+PTA(1)-PPP1(1)-PPP2(1)-PPP3(1)-PPP4(1)
	PYCHK=PPR(2)+PTA(2)-PPP1(2)-PPP2(2)-PPP3(2)-PPP4(2)
	PZCHK=PPR(3)+PTA(3)-PPP1(3)-PPP2(3)-PPP3(3)-PPP4(3)
	PECHK=PPR(4)+PTA(4)-PPP1(4)-PPP2(4)-PPP3(4)-PPP4(4)
	IF(IPCO.GE.1)THEN
        WRITE(6,'(A/8E12.4,I5)')
     *	' Chain masses AMCH,AMCHOR,AMCHN1,AMCHN2,PZCHK,PECHK,
     *   PXCHK,PYCHK,NOBAM ',
     *	 AMCH,AMCHOR,AMCHN1,AMCHN2,PZCHK,PECHK,PXCHK,PYCHK,NOBAM   
        ENDIF
C                 Lorentz parameters of chains
        GAMOR=(PPR(4)+PTA(4))/AMCH
	BGXOR=(PPR(1)+PTA(1))/AMCH
	BGYOR=(PPR(2)+PTA(2))/AMCH
	BGZOR=(PPR(3)+PTA(3))/AMCH
        GAMCH1=(PPP1(4)+PPP2(4))/AMCHN1
        BGXCH1=(PPP1(1)+PPP2(1))/AMCHN1
        BGYCH1=(PPP1(2)+PPP2(2))/AMCHN1
        BGZCH1=(PPP1(3)+PPP2(3))/AMCHN1
        GAMCH2=(PPP3(4)+PPP4(4))/AMCHN2
        BGXCH2=(PPP3(1)+PPP4(1))/AMCHN2
        BGYCH2=(PPP3(2)+PPP4(2))/AMCHN2
        BGZCH2=(PPP3(3)+PPP4(3))/AMCHN2
	IF(IPCO.GE.1)THEN
 	WRITE(6,2346)' L.Parm in ',GAM,BGX,BGY,BGZ
 	WRITE(6,2346)' L.Parm OR ',GAMOR,BGXOR,BGYOR,BGZOR
 	WRITE(6,2346)' L.Parm C1 ',GAMCH1,BGXCH1,BGYCH1,BGZCH1
 	WRITE(6,2346)' L.Parm C2 ',GAMCH2,BGXCH2,BGYCH2,BGZCH2
	ENDIF
	PEOR=AMCH*GAMOR
	PZOR=AMCH*BGZOR
	PECH1=AMCHN1*GAMCH1
	PZCH1=AMCHN1*BGZCH1
	PECH2=AMCHN2*GAMCH2
	PZCH2=AMCHN2*BGZCH2
	IF(IPCO.GE.1)THEN
	WRITE(6,'(A,6E12.4)')' PEOR,PECH1,PECH2,PZOR,PZCH1,PZCH2',
     *	 PEOR,PECH1,PECH2,PZOR,PZCH1,PZCH2
	ENDIF
	NOBA1=3
	NOBA2=6
C	WRITE(6,'(A,2I5)')' Jet1 IFB1,ISAQ ',IFB1,ISAQ
      CALL HADJET(NHAD1,AMCHN1,PPP1,PPP2,GAMCH1,BGXCH1,BGYCH1,
     * BGZCH1, IFB1,ISAQ,IFB3,IFB4,I1,I2,NOBA1,NNCH,NORIG)
C       DO 42 I=1,NHAD1
C       WRITE(6,1050)I,PXF(I),PYF(I),PZF(I),HEF(I),AMF(I), ICHF(I),
C    +    IBARF(I),NREF(I),ANF(I)
C  42   CONTINUE
C     DIMENSION ANFF(NFIMAX),PXFF(NFIMAX),PYFF(NFIMAX),PZFF(NFIMAX),
C    +HEFF(NFIMAX),AMFF(NFIMAX), 
C    *ICHFF(NFIMAX),IBARFF(NFIMAX),NREFF(NFIMAX)
C                 Intermediate store of hadrons from jet 1
      DO 2348 I=1,NHAD1
	ANFF(I)=ANF(I)
	PXFF(I)=PXF(I)
	PYFF(I)=PYF(I)
	PZFF(I)=PZF(I)
	HEFF(I)=HEF(I)
	AMFF(I)=AMF(I)
	ICHFF(I)=ICHF(I)
	IBARFF(I)=IBARF(I)
	NREFF(I)=NREF(I)
 2348 CONTINUE
C	WRITE(6,'(A,3I5)')' Jet2 IFB2,ISQ,IFB3 ',IFB2,ISQ,IFB3
      CALL HADJET(NHAD2,AMCHN2,PPP3,PPP4,GAMCH2,BGXCH2,BGYCH2,
     * BGZCH2, IFB2,ISQ,IFB3,IFB4,I1,I2,NOBA2,NNCH,NORIG)
C       DO 41 I=1,NHAD2
C       WRITE(6,1050)I,PXF(I),PYF(I),PZF(I),HEF(I),AMF(I), ICHF(I),
C    +    IBARF(I),NREF(I),ANF(I)
C  41   CONTINUE
      NHAD=NHAD1+NHAD2
      DO 2349 I=NHAD2+1,NHAD
	II=I-NHAD2
	ANF(I)=ANFF(II)
	PXF(I)=PXFF(II)
	PYF(I)=PYFF(II)
	PZF(I)=PZFF(II)
	HEF(I)=HEFF(II)
	AMF(I)=AMFF(II)
	ICHF(I)=ICHFF(II)
	IBARF(I)=IBARFF(II)
	NREF(I)=NREFF(II)
 2349 CONTINUE
      ENDIF
C==================================================================
C==================================================================
C              Different options:
      IF(NOBAM.EQ.4)THEN
C              NOBAM=4:  Quark-DIQUARK Jet
C
C              Work out approximate x-Values of quark and diquark
	XDIQ=PTA(4)*2.D0/CMENER
	XQUA=PPR(4)*2.D0/CMENER
C	WRITE(6,'(A,2E12.4)')' HADJCK:XDIQ,XQUA ',XDIQ,XQUA
C              Select valence quark x for one of the diquark quarks
C              But only if XDIQ > 1.5*XQUA
C	IF(XDIQ.LE.1.5D0*XQUA)THEN
C	  IREJ=1
C	  RETURN
C       ENDIF
C              Select x values of sea quark pair
	ICOU=0
 2234   CONTINUE
	ICOU=ICOU+1
	IF(ICOU.GE.100)THEN
	  IREJ=1
	  RETURN
	ENDIF
	CALL XSEAPA(CMENER,XQUA/2.D0,ISQ,ISAQ,XSQ,XSAQ,IREJ)
	IF(XSAQ.GE.XQUA/2.D0)GO TO 2234
	IF(IPCO.GE.1)THEN
        WRITE(6,'(A,2E12.4)')' HADJCK:XSQ,XSAQ ',XSQ,XSAQ
        ENDIF
*  target valence quarks xtvqi
        XVTHRO=CVQ/CMENER
	IVTHR=0
 3466   CONTINUE
	IF(IVTHR.EQ.50)THEN
	  IREJ=1
      WRITE(6,*)' HADJSE 3466 reject IVTHR 50'
	  RETURN
        ENDIF
	IVTHR=IVTHR+1
        XVTHR=XVTHRO/(51-IVTHR)
        UNOPRV=UNON
  380   CONTINUE
        IF(XVTHR.GT.0.05)THEN
          IF(XVTHR.GT.0.66D0*XDIQ)THEN
            IREJ=1
	    RETURN
	  ENDIF
          XTVQI=BETREJ(0.5,UNOPRV,XVTHR,0.66D0*XDIQ)
        ELSE
          N90=0
  390     CONTINUE
          N90=N90+1
	  IF(N90.GE.20)THEN
	    IREJ=1
	    RETURN
	  ENDIF
          XTVQI=DBETAR(0.5,UNOPRV)
          IF ((XTVQI.LT.XVTHR).OR.(XTVQI.GT.0.66D0*XDIQ))
     *                                                         GOTO 390
        ENDIF
	XDIQQ=XDIQ-XTVQI
	IF(IPCO.GE.1)THEN
 	WRITE(6,'(A,2E12.4)')' HADJCK:XDIQQ,XTVQI ',XDIQQ,XTVQI
        ENDIF
C               We form two strings(2-1) aq-q and (4-3) q-qq
C            1  q with XDIQ-XTVQI-XSQ
C            1  q with XDIQ-XTVQI   (j.r.5.5.96)
C            2  aq with XSAQ
C            3  qq with XTVQI+XSQ
C            3  qq with XTVQI     (j.r.5.5.96)
C            4  q with  XQUA-XSAQ
C            with 4-momenta PPP1,PPP2,PPP3,PPP4
      DO 3345 I=1,4
C       PPP1(I)=PTA(I)*(XDIQ-XTVQI-XSQ)/XDIQ
        PPP1(I)=PTA(I)*(XDIQ-XTVQI)/XDIQ
C	PPP3(I)=PTA(I)*(XTVQI+XSQ)/XDIQ
	PPP3(I)=PTA(I)*(XTVQI)/XDIQ
	PPP2(I)=PPR(I)*XSAQ/XQUA
	PPP4(I)=PPR(I)*(XQUA-XSAQ)/XQUA
 3345   CONTINUE
 3346   FORMAT(A,4E12.4)
	IF(IPCO.GE.1)THEN
 	WRITE(6,3346)' PPR ',PPR
 	WRITE(6,3346)' PPP1 ',PPP1
 	WRITE(6,3346)' PPP3 ',PPP3
 	WRITE(6,3346)' PTA ',PTA
 	WRITE(6,3346)' PPP2 ',PPP2
 	WRITE(6,3346)' PPP4 ',PPP4
        ENDIF
C               Invariant Masses of chains
C               ORIGINAL CHAIN
	AMCHOR=SQRT((PPR(4)+PTA(4))**2-(PPR(3)+PTA(3))**2
     *             -(PPR(2)+PTA(2))**2-(PPR(1)+PTA(1))**2) 
	AMCHN1=SQRT((PPP1(4)+PPP2(4))**2-(PPP1(3)+PPP2(3))**2
     *             -(PPP1(2)+PPP2(2))**2-(PPP1(1)+PPP2(1))**2) 
C                    Chain 1 is q-aq restrict mass >1.5 GeV
	CHAMAL=1.D0
	IF(IFB3.GE.3.OR.ISAQ.GE.9)CHAMAL=1.5D0
	IF(AMCHN1.LE.CHAMAL)THEN
C	  IREJ=1
C	  RETURN
	  GO TO 3466
	ENDIF
	AMCHN2=SQRT((PPP3(4)+PPP4(4))**2-(PPP3(3)+PPP4(3))**2
     *             -(PPP3(2)+PPP4(2))**2-(PPP3(1)+PPP4(1))**2) 
C                    Chain 2 is qq-q restrict mass >2.5 GeV
	CHAMAL=1.8D0
	IF(IFB2.GE.3.OR.IFB1.GE.3.OR.ISQ.GE.3)CHAMAL=2.5D0
	IF(AMCHN2.LE.CHAMAL)THEN
C	  IREJ=1
C	  RETURN
	  GO TO 3466
	ENDIF
	PXCHK=PPR(1)+PTA(1)-PPP1(1)-PPP2(1)-PPP3(1)-PPP4(1)
	PYCHK=PPR(2)+PTA(2)-PPP1(2)-PPP2(2)-PPP3(2)-PPP4(2)
	PZCHK=PPR(3)+PTA(3)-PPP1(3)-PPP2(3)-PPP3(3)-PPP4(3)
	PECHK=PPR(4)+PTA(4)-PPP1(4)-PPP2(4)-PPP3(4)-PPP4(4)
	IF(IPCO.GE.1)THEN
        WRITE(6,'(A/8E12.4,I5)')
     *	' Chain masses AMCH,AMCHOR,AMCHN1,AMCHN2,PZCHK,PECHK,
     *   PXCHK,PYCHK,NOBAM ',
     *	 AMCH,AMCHOR,AMCHN1,AMCHN2,PZCHK,PECHK,PXCHK,PYCHK,NOBAM   
        ENDIF
C                 Lorentz parameters of chains
        GAMOR=(PPR(4)+PTA(4))/AMCH
	BGXOR=(PPR(1)+PTA(1))/AMCH
	BGYOR=(PPR(2)+PTA(2))/AMCH
	BGZOR=(PPR(3)+PTA(3))/AMCH
        GAMCH1=(PPP1(4)+PPP2(4))/AMCHN1
        BGXCH1=(PPP1(1)+PPP2(1))/AMCHN1
        BGYCH1=(PPP1(2)+PPP2(2))/AMCHN1
        BGZCH1=(PPP1(3)+PPP2(3))/AMCHN1
        GAMCH2=(PPP3(4)+PPP4(4))/AMCHN2
        BGXCH2=(PPP3(1)+PPP4(1))/AMCHN2
        BGYCH2=(PPP3(2)+PPP4(2))/AMCHN2
        BGZCH2=(PPP3(3)+PPP4(3))/AMCHN2
	IF(IPCO.GE.1)THEN
 	WRITE(6,3346)' L.Parm in ',GAM,BGX,BGY,BGZ
 	WRITE(6,3346)' L.Parm OR ',GAMOR,BGXOR,BGYOR,BGZOR
 	WRITE(6,3346)' L.Parm C1 ',GAMCH1,BGXCH1,BGYCH1,BGZCH1
 	WRITE(6,3346)' L.Parm C2 ',GAMCH2,BGXCH2,BGYCH2,BGZCH2
	ENDIF
	PEOR=AMCH*GAMOR
	PZOR=AMCH*BGZOR
	PECH1=AMCHN1*GAMCH1
	PZCH1=AMCHN1*BGZCH1
	PECH2=AMCHN2*GAMCH2
	PZCH2=AMCHN2*BGZCH2
	IF(IPCO.GE.1)THEN
	WRITE(6,'(A,6E12.4)')' PEOR,PECH1,PECH2,PZOR,PZCH1,PZCH2',
     *	 PEOR,PECH1,PECH2,PZOR,PZCH1,PZCH2
	ENDIF
	NOBA1=3
	NOBA2=4
C	WRITE(6,'(A,2I5)')' Jet1 ISAQ,IFB3 ',ISAQ,IFB3
      CALL HADJET(NHAD1,AMCHN1,PPP2,PPP1,GAMCH1,BGXCH1,BGYCH1,
     * BGZCH1, ISAQ,IFB3,IFB1,IFB4,I1,I2,NOBA1,NNCH,NORIG)
C       DO 342 I=1,NHAD1
C       WRITE(6,1050)I,PXF(I),PYF(I),PZF(I),HEF(I),AMF(I), ICHF(I),
C    +    IBARF(I),NREF(I),ANF(I)
C 342   CONTINUE
C     DIMENSION ANFF(NFIMAX),PXFF(NFIMAX),PYFF(NFIMAX),PZFF(NFIMAX),
C    +HEFF(NFIMAX),AMFF(NFIMAX), 
C    *ICHFF(NFIMAX),IBARFF(NFIMAX),NREFF(NFIMAX)
C                 Intermediate store of hadrons from jet 1
      DO 3348 I=1,NHAD1
	ANFF(I)=ANF(I)
	PXFF(I)=PXF(I)
	PYFF(I)=PYF(I)
	PZFF(I)=PZF(I)
	HEFF(I)=HEF(I)
	AMFF(I)=AMF(I)
	ICHFF(I)=ICHF(I)
	IBARFF(I)=IBARF(I)
	NREFF(I)=NREF(I)
 3348 CONTINUE
C	WRITE(6,'(A,3I5)')' Jet2 IFB1,ISQ,IFB2 ',IFB1,ISQ,IFB2
      CALL HADJET(NHAD2,AMCHN2,PPP4,PPP3,GAMCH2,BGXCH2,BGYCH2,
     * BGZCH2, IFB1,ISQ,IFB2,IFB4,I1,I2,NOBA2,NNCH,NORIG)
C       DO 341 I=1,NHAD2
C       WRITE(6,1050)I,PXF(I),PYF(I),PZF(I),HEF(I),AMF(I), ICHF(I),
C    +    IBARF(I),NREF(I),ANF(I)
C 341   CONTINUE
      NHAD=NHAD1+NHAD2
      DO 3349 I=NHAD2+1,NHAD
	II=I-NHAD2
	ANF(I)=ANFF(II)
	PXF(I)=PXFF(II)
	PYF(I)=PYFF(II)
	PZF(I)=PZFF(II)
	HEF(I)=HEFF(II)
	AMF(I)=AMFF(II)
	ICHF(I)=ICHFF(II)
	IBARF(I)=IBARFF(II)
	NREF(I)=NREFF(II)
 3349 CONTINUE
      ENDIF
C==================================================================
	HEFT=0.D0
	IF(IPCO.GE.1)THEN
        DO 40 I=1,NHAD
	  IF(IBARF(I).EQ.500)GO TO 4040
          HEFT=HEFT+HEF(I)
        WRITE(6,4050)I,PXF(I),PYF(I),PZF(I),HEF(I),AMF(I), ICHF(I),
     +    IBARF(I),NREF(I),ANF(I)
 4050 FORMAT(' JET  ',I5,5F12.4,3I5,A10)
 4040   CONTINUE
   40   CONTINUE
	ENDIF
      HEFFF=AMCH*GAM
	IF(IPCO.GE.1)THEN
      WRITE(6,'(A,I5,2E12.4)')'1IREJ,HEFT,HEFFF ',
     * IREJ,HEFT,HEFFF
	ENDIF
C     IREJ=1
      IF(IREJ.EQ.1)RETURN
C
      RETURN
      END

************************************************************************
************************************************************************
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HADJSE(NHAD,AMCH,PPR,PTA,GAM,BGX,BGY,BGZ, IFB1,IFB2,
     +IFB3,IFB4,I1,I2,NOBAM,NNCH,NORIG,IREJ,IISSQQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C   HADJSE Fragmentation of vv,vs and sv chains CK method 
C          using Glauber sea quarks
C                 j.r.7/96
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
      CHARACTER*8 ANF,ANFF
      PARAMETER (NFIMAX=249)
      COMMON /DFINPA/ ANF(NFIMAX),PXF(NFIMAX),PYF(NFIMAX),PZF(NFIMAX),
     +HEF(NFIMAX),AMF(NFIMAX), ICHF(NFIMAX),IBARF(NFIMAX),NREF(NFIMAX)
      COMMON /DFINPZ/IORMO(NFIMAX),IDAUG1(NFIMAX),IDAUG2(NFIMAX),
     * ISTATH(NFIMAX)
      DIMENSION ANFF(NFIMAX),PXFF(NFIMAX),PYFF(NFIMAX),PZFF(NFIMAX),
     +HEFF(NFIMAX),AMFF(NFIMAX), 
     *ICHFF(NFIMAX),IBARFF(NFIMAX),NREFF(NFIMAX),IORMOO(NFIMAX)
        CHARACTER*80 TITLE
        CHARACTER*8 PROJTY,TARGTY
C       COMMON /USER/TITLE,PROJTY,TARGTY,CMENER,ISTRUF
C    & ,ISINGD,IDUBLD,SDFRAC,PTLAR
      COMMON /USER1/TITLE,PROJTY,TARGTY
      COMMON /USER2/CMENER,SDFRAC,PTLAR,ISTRUF,ISINGD,IDUBLD
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
      DIMENSION PPR(4),PTA(4)
*KEEP,XSEADI.
      COMMON /XSEADI/ XSEACU,UNON,UNOM,UNOSEA,CVQ,CDQ,CSEA,SSMIMA,
     +SSMIMQ,VVMTHR
      COMMON /HDJASE/NHSE1,NHSE2,NHSE3,NHASE1,NHASE2,NHASE3
      DIMENSION PPP1(4),PPP2(4),PPP3(4),PPP4(4)
      DATA TINY/1.0D-3/
C----------------------------------------------------------------------
C
C	WRITE(6,'(A,3I5)')' Jet0 IFB1,IFB2,IFB3 ',IFB1,IFB2,IFB3
      IF(IPCO.GE.0)WRITE(6,*)
     *' HADJSE(NHAD,AMCH,PPR,PTA,GAM,BGX,BGY,BGZ, IFB1,IFB2,',
     +'IFB3,IFB4,I1,I2,NOBAM,NNCH,NORIG,IREJ),IPCO',
     *NHAD,AMCH,PPR,PTA,GAM,BGX,BGY,BGZ, IFB1,IFB2,
     +IFB3,IFB4,I1,I2,NOBAM,NNCH,NORIG,IREJ,IPCO
      IREJ=0
      IF(IPCO.GE.0) THEN
    	WRITE(6,'(A,3I5)')' HADJSE Jet0 IFB1,IFB2,IFB3',IFB1,IFB2,IFB3
        WRITE(6,'(4E12.4)') GAM,BGX,BGY,BGZ,PPR
        WRITE(6,1000) NHAD,AMCH,PTA, IFB1,IFB2,IFB3,IFB4,I1,I2,NOBAM,
     +  NNCH,NORIG
 1000   FORMAT(10X,I10,5F10.3/10X,9I10)
      ENDIF
      IF(ABS(NNCH).EQ.99) THEN
        NHAD=0
C       IPCO=0
        RETURN
      ENDIF
C=================================================================
C=================================================================
C                   DIQUARK-QUARK NOBAM=6
C=================================================================
C=================================================================
C              Different options:
      IF(NOBAM.EQ.6)THEN
        IF(IPCO.GE.0)WRITE(6,*)' DIQUARK-QUARK NOBAM=6'
C              NOBAM=6:  Diquark--Quark Jet
C
C              Work out approximate x-Values of diquark and quark
	IF(CMENER.LT.1.D-3)THEN
	  WRITE(6,*)' CMENER=0. HADJSE',CMENER
	  STOP
	ENDIF
	XDIQ=PPR(4)*2.D0/CMENER
	XQUA=PTA(4)*2.D0/CMENER
        IF(IPCO.GE.0) THEN
          WRITE(6,'(A,2E12.4)')' HADJSE:XDIQ,XQUA ',XDIQ,XQUA
        ENDIF
C       Select sea quark x for one of the diquark quarks
C       Select x values of sea quark pair
	ICOU=0
 1234   CONTINUE
	ICOU=ICOU+1
	IF(ICOU.GE.200)THEN
	  IREJ=1
          IF(IPCO.GE.0)WRITE(6,*)' HADJSE reject icou 100'
	  RETURN
	ENDIF
	CALL XSEAPA(CMENER,XQUA/2.D0,ISQ,ISAQ,XSQ,XSAQ,IREJ)
	IISSQQ=ISQ
	IF(IREJ.GE.1)THEN
          IF(IPCO.GE.0)WRITE(6,*)' HADJSE reject XSEAPA'
	  RETURN
	ENDIF
	IF(XSAQ.GE.2.D0*XQUA/3.D0)GO TO 1234
        IF(IPCO.GE.0) THEN
          WRITE(6,*)' HADJSE,XSEAPA:XSQ,XSAQ,ISQ,ISAQ ',
     *	  XSQ,XSAQ,ISQ,ISAQ
        ENDIF
*       projectile valence quarks xpvqi
C       here selected as sea quark! from 1/x-sea
        XVTHRO=CVQ/CMENER
	IVTHR=0
 3465   CONTINUE
	IF(IVTHR.EQ.100)THEN
	  IREJ=1
	  IF(ISQ.EQ.3)IREJ=3
          IF(IPCO.GE.0)WRITE(6,*)' HADJSE 3465 reject IVTHR 50'
	  RETURN
        ENDIF
	IVTHR=IVTHR+1
        XVTHR=XVTHRO/(101-IVTHR)
        UNOPRV=UNON
   80   CONTINUE
        IF(XVTHR.GT.0.66D0*XDIQ)THEN
          IF(IPCO.GE.0)WRITE(6,*)' HADJSE 80 reject XVTHR too great',
     *	  XVTHR
          IREJ=1
	  IF(ISQ.EQ.3)IREJ=3
	  RETURN
	ENDIF
        XPVQI=SAMPEY(XVTHR,0.66D0*XDIQ)
        XDIQQ=XDIQ-XPVQI
        IF(IPCO.GE.0) THEN
          WRITE(6,'(A,2E12.4)')' HADJSE:XDIQQ,XPVQI ',XDIQQ,XPVQI
        ENDIF
C
C       We form two strings(1-2) q-aq and (3-4) qq-q
C       1  q with XDIQ-XPVQI-XSQ
C       1  q with XDIQ-XPVQI   (j.r.5.5.96)
C       2  aq with XSAQ
C       3  qq with XPVQI+XSQ
C       3  qq with XPVQI       (j.r.5.5.96)
C       4  q with  XQUA-XSAQ
C       with 4-momenta PPP1,PPP2,PPP3,PPP4
        DO 2345 I=1,4
          IF(XDIQ.LE.1.D-15.OR.XQUA.LT.1.D-15)THEN
            IREJ=1
	    IF(ISQ.EQ.3)IREJ=3
            IF(IPCO.GE.0)WRITE(6,*)' HADJSE reject 2345 XDIQ,',
     *	    'XQUA too small ',
     *	    XDIQ,XQUA
            RETURN
          ENDIF
C         PPP1(I)=PPR(I)*(XDIQ-XPVQI-XSQ)/XDIQ
          PPP1(I)=PPR(I)*(XDIQ-XPVQI)/XDIQ
C	  PPP3(I)=PPR(I)*(XPVQI+XSQ)/XDIQ
          PPP3(I)=PPR(I)*(XPVQI)/XDIQ
          PPP2(I)=PTA(I)*XSAQ/XQUA
          PPP4(I)=PTA(I)*(XQUA-XSAQ)/XQUA
 2345   CONTINUE
 2346   FORMAT(A,5E12.4)
        IF(IPCO.GE.0) THEN
          WRITE(6,2346)' PPR ',PPR
          WRITE(6,2346)' PPP1 ',PPP1
          WRITE(6,2346)' PPP3 ',PPP3
          WRITE(6,2346)' PTA ',PTA
          WRITE(6,2346)' PPP2 XSAQ',PPP2,XSAQ
          WRITE(6,2346)' PPP4 ',PPP4
        ENDIF
C       Invariant Masses of chains
C       ORIGINAL CHAIN
        AMCHOR=SQRT((PPR(4)+PTA(4))**2-(PPR(3)+PTA(3))**2
     *             -(PPR(2)+PTA(2))**2-(PPR(1)+PTA(1))**2) 
        AMCHN1=SQRT((PPP1(4)+PPP2(4))**2-(PPP1(3)+PPP2(3))**2
     *             -(PPP1(2)+PPP2(2))**2-(PPP1(1)+PPP2(1))**2) 
        IF(AMCHOR.LE.TINY)THEN
          WRITE(6,2346)' PPR ',PPR
          WRITE(6,2346)' PPP1 ',PPP1
          WRITE(6,2346)' PPP3 ',PPP3
          WRITE(6,2346)' PTA ',PTA
          WRITE(6,2346)' PPP2 ',PPP2
          WRITE(6,2346)' PPP4 ',PPP4
        ENDIF
        IF(AMCHN1.LE.TINY)THEN
          WRITE(6,2346)' PPR ',PPR
          WRITE(6,2346)' PPP1 ',PPP1
          WRITE(6,2346)' PPP3 ',PPP3
          WRITE(6,2346)' PTA ',PTA
          WRITE(6,2346)' PPP2 ',PPP2
          WRITE(6,2346)' PPP4 ',PPP4
        ENDIF
C       Chain 1 is q-aq restrict mass >1.5 GeV 
        CHAMAL=0.8D0
        IF(IFB1.GE.3.OR.ISAQ.GE.9)CHAMAL=1.2D0
        IF(AMCHN1.LE.CHAMAL)THEN
          IF(IPCO.GE.0)WRITE (6,*)'HADJSE jump1AMCHN1.LE.CHAMAL AMCHOR',
     *	  AMCHN1,CHAMAL,AMCHOR
          GO TO 3465
        ENDIF
        AMCHN2=SQRT((PPP3(4)+PPP4(4))**2-(PPP3(3)+PPP4(3))**2
     *             -(PPP3(2)+PPP4(2))**2-(PPP3(1)+PPP4(1))**2) 
C                    Chain 2 is qq-q restrict mass >2.5 GeV
	CHAMAL=1.5D0
	IF(IFB2.GE.3.OR.IFB3.GE.3.OR.ISQ.GE.3)CHAMAL=2.0D0
	IF(AMCHN2.LE.CHAMAL)THEN
C	  IREJ=1
C 	  RETURN
          IF(IPCO.GE.0)WRITE (6,*)'HADJSE jump AMCHN2.LE.CHAMAL ',
     *	  AMCHN2,CHAMAL
	  GO TO 3465
        ENDIF
        PXCHK=PPR(1)+PTA(1)-PPP1(1)-PPP2(1)-PPP3(1)-PPP4(1)
        PYCHK=PPR(2)+PTA(2)-PPP1(2)-PPP2(2)-PPP3(2)-PPP4(2)
        PZCHK=PPR(3)+PTA(3)-PPP1(3)-PPP2(3)-PPP3(3)-PPP4(3)
        PECHK=PPR(4)+PTA(4)-PPP1(4)-PPP2(4)-PPP3(4)-PPP4(4)
	IF(IPCO.GE.0)THEN
          WRITE(6,'(A/8E12.4,I5)')
     *	  ' Chain masses AMCH,AMCHOR,AMCHN1,AMCHN2,PZCHK,PECHK,
     *    PXCHK,PYCHK,NOBAM ',
     *	  AMCH,AMCHOR,AMCHN1,AMCHN2,PZCHK,PECHK,PXCHK,PYCHK,NOBAM   
        ENDIF
C       Lorentz parameters of chains
        IF(AMCH.LE.1.D-15)THEN
	  IREJ=1
	    IF(ISQ.EQ.3)IREJ=3
          IF(IPCO.GE.0)WRITE(6,*)' HADJSE rejection AMCH too small ',
     *	  AMCH
	  RETURN
        ENDIF
        GAMOR=(PPR(4)+PTA(4))/AMCH
        BGXOR=(PPR(1)+PTA(1))/AMCH
        BGYOR=(PPR(2)+PTA(2))/AMCH
        BGZOR=(PPR(3)+PTA(3))/AMCH
        IF(AMCHN1.LE.1.D-15)THEN
          IF(IPCO.GE.0)WRITE(6,*)' HADJSE rejection AMCHN1 too small ',
     *	  AMCHN1
	  IREJ=1
	    IF(ISQ.EQ.3)IREJ=3
	  RETURN
	ENDIF
        GAMCH1=(PPP1(4)+PPP2(4))/AMCHN1
        BGXCH1=(PPP1(1)+PPP2(1))/AMCHN1
        BGYCH1=(PPP1(2)+PPP2(2))/AMCHN1
        BGZCH1=(PPP1(3)+PPP2(3))/AMCHN1
        IF(AMCHN2.LE.1.D-15)THEN
          IF(IPCO.GE.0)WRITE(6,*)' HADJSE rejection AMCHN2 too small ',
     *	  AMCHN2
	  IREJ=1
	    IF(ISQ.EQ.3)IREJ=3
	  RETURN
	ENDIF
        GAMCH2=(PPP3(4)+PPP4(4))/AMCHN2
        BGXCH2=(PPP3(1)+PPP4(1))/AMCHN2
        BGYCH2=(PPP3(2)+PPP4(2))/AMCHN2
        BGZCH2=(PPP3(3)+PPP4(3))/AMCHN2
	IF(IPCO.GE.0)THEN
 	  WRITE(6,2346)' L.Parm in ',GAM,BGX,BGY,BGZ
 	  WRITE(6,2346)' L.Parm OR ',GAMOR,BGXOR,BGYOR,BGZOR
 	  WRITE(6,2346)' L.Parm C1 ',GAMCH1,BGXCH1,BGYCH1,BGZCH1
 	  WRITE(6,2346)' L.Parm C2 ',GAMCH2,BGXCH2,BGYCH2,BGZCH2
	ENDIF
	PEOR=AMCH*GAMOR
	PZOR=AMCH*BGZOR
	PECH1=AMCHN1*GAMCH1
	PZCH1=AMCHN1*BGZCH1
	PECH2=AMCHN2*GAMCH2
	PZCH2=AMCHN2*BGZCH2
	IF(IPCO.GE.0)THEN
	  WRITE(6,'(A,6E12.4)')' PEOR,PECH1,PECH2,PZOR,PZCH1,PZCH2',
     *	  PEOR,PECH1,PECH2,PZOR,PZCH1,PZCH2
	ENDIF
	NOBA1=3
	NOBA2=6
	IF(IPCO.GE.-1)THEN
C	  WRITE(6,'(A,2I5)')' Jet1 IFB1,ISAQ ',IFB1,ISAQ
	ENDIF
        CALL HADJET(NHAD1,AMCHN1,PPP1,PPP2,GAMCH1,BGXCH1,BGYCH1,
     *  BGZCH1, IFB1,ISAQ,IFB3,IFB4,I1,I2,NOBA1,NNCH,NORIG)
C       DO 42 I=1,NHAD1
C         WRITE(6,1050)I,PXF(I),PYF(I),PZF(I),HEF(I),AMF(I), ICHF(I),
C    +    IBARF(I),NREF(I),ANF(I)
C  42   CONTINUE
	IF(IPCO.GE.-1)THEN
C	  WRITE(6,'(A,2I5)')' Jet1 NHAD1 ',NHAD1
	ENDIF
C       DIMENSION ANFF(NFIMAX),PXFF(NFIMAX),PYFF(NFIMAX),PZFF(NFIMAX),
C    +  HEFF(NFIMAX),AMFF(NFIMAX), 
C    *  ICHFF(NFIMAX),IBARFF(NFIMAX),NREFF(NFIMAX)
C       Intermediate store of hadrons from jet 1
        DO 2348 I=1,NHAD1
          ANFF(I)=ANF(I)
          PXFF(I)=PXF(I)
          PYFF(I)=PYF(I)
          PZFF(I)=PZF(I)
          HEFF(I)=HEF(I)
          AMFF(I)=AMF(I)
          ICHFF(I)=ICHF(I)
          IBARFF(I)=IBARF(I)
          NREFF(I)=NREF(I)
	  IORMOO(I)=IORMO(I)
 2348   CONTINUE
	IF(IPCO.GE.0)THEN
          WRITE(6,'(A,3I5)')' Jet2 IFB2,ISQ,IFB3 ',IFB2,ISQ,IFB3
	ENDIF
        CALL HADJET(NHAD2,AMCHN2,PPP3,PPP4,GAMCH2,BGXCH2,BGYCH2,
     *  BGZCH2, IFB2,ISQ,IFB3,IFB4,I1,I2,NOBA2,NNCH,NORIG)
	IF(IPCO.GE.-1)THEN
C	  WRITE(6,'(A,2I5)')' Jet2 NHAD2 ',NHAD2
	ENDIF
C       DO 41 I=1,NHAD2
C         WRITE(6,1050)I,PXF(I),PYF(I),PZF(I),HEF(I),AMF(I), ICHF(I),
C    +    IBARF(I),NREF(I),ANF(I)
C        IF(IORMO(I).NE.999)IORMO(I+NHAD1)=IORMO(I)+NHAD1
C  41   CONTINUE
C       Intermediate store of hadrons from jet 2
        DO 2448 I=NHAD1+1,NHAD1+NHAD2
	  II=I-NHAD1
          ANFF(I)=ANF(II)
          PXFF(I)=PXF(II)
          PYFF(I)=PYF(II)
          PZFF(I)=PZF(II)
          HEFF(I)=HEF(II)
          AMFF(I)=AMF(II)
          ICHFF(I)=ICHF(II)
          IBARFF(I)=IBARF(II)
          NREFF(I)=NREF(II)
	  IORMOO(I)=IORMO(II)
	  IF(IORMOO(I).NE.999)IORMOO(I)=IORMOO(I)+NHAD1
 2448   CONTINUE
        NHAD=NHAD1+NHAD2
        DO 2349 I=1,NHAD
	  II=I
	  ANF(I)=ANFF(II)
	  PXF(I)=PXFF(II)
	  PYF(I)=PYFF(II)
	  PZF(I)=PZFF(II)
	  HEF(I)=HEFF(II)
	  AMF(I)=AMFF(II)
	  ICHF(I)=ICHFF(II)
	  IBARF(I)=IBARFF(II)
	  NREF(I)=NREFF(II)
	  IORMO(I)=IORMOO(II)
 2349   CONTINUE
      ENDIF
C=================================================================
C=================================================================
C                   QUARK-DIQUARK NOBAM=4
C=================================================================
C==================================================================
C     Different options:
      IF(NOBAM.EQ.4)THEN
        IF(IPCO.GE.0)WRITE(6,*)' QUARK-DIQUARK NOBAM=4'
C       NOBAM=4:  Quark-DIQUARK Jet
C
C       Work out approximate x-Values of quark and diquark
	IF(CMENER.LT.1.D-3)THEN
	  IF(IPCO.GE.0)WRITE(6,*)' CMENER=0. HADJSE',CMENER
	  STOP
	ENDIF
	XDIQ=PTA(4)*2.D0/CMENER
	XQUA=PPR(4)*2.D0/CMENER
	IF(IPCO.GE.0)THEN
          WRITE(6,'(A,2E12.4)')' HADJSE:XDIQ,XQUA ',XDIQ,XQUA
	ENDIF
C       Select valence quark x for one of the diquark quarks
C       Select x values of sea quark pair
	ICOU=0
 2234   CONTINUE
	ICOU=ICOU+1
	IF(ICOU.GE.200)THEN
	  IREJ=1
	    IF(ISQ.EQ.3)IREJ=3
          IF(IPCO.GE.0)WRITE(6,*)' HADJSE Rejection 2234 ICOU. GT.100'
	  RETURN
	ENDIF
	IF(IPCO.GE.0)WRITE(6,*)' XSEAPA: CMENER,XQUA ',CMENER,XQUA
	CALL XSEAPA(CMENER,XQUA/2.D0,ISQ,ISAQ,XSQ,XSAQ,IREJ)
	IISSQQ=ISQ
	IF(IREJ.GE.1)THEN
          IF(IPCO.GE.0)WRITE(6,*)' HADJSE reject XSEAPA'
	  RETURN
	ENDIF
	IF(XSAQ.GE.2.D0*XQUA/3.D0)GO TO 2234
	IF(IPCO.GE.0)THEN
          WRITE(6,'(A,2E12.4)')' HADJSE:XSQ,XSAQ ',XSQ,XSAQ
        ENDIF
*       target valence quarks xtvqi
C       here selected as sea quark! from 1/x-sea
        XVTHRO=CVQ/CMENER
	IVTHR=0
 3466   CONTINUE
	IF(IVTHR.EQ.100)THEN
	  IREJ=1
	    IF(ISQ.EQ.3)IREJ=3
          IF(IPCO.GE.0)WRITE(6,*)' HADJSE 3466 reject IVTHR 50'
	  RETURN
        ENDIF
	IVTHR=IVTHR+1
        XVTHR=XVTHRO/(101-IVTHR)
        UNOPRV=UNON
  380   CONTINUE
        IF(XVTHR.GT.0.66D0*XDIQ)THEN
          IREJ=1
	    IF(ISQ.EQ.3)IREJ=3
          IF(IPCO.GE.0)WRITE(6,*)' HADJSE Rejection 380 XVTHR  large ',
     *	  XVTHR
          RETURN
        ENDIF
        XTVQI=SAMPEY(XVTHR,0.66D0*XDIQ)
	XDIQQ=XDIQ-XTVQI
	IF(IPCO.GE.0)THEN
          WRITE(6,'(A,2E12.4)')' HADJCK:XDIQQ,XTVQI ',XDIQQ,XTVQI
        ENDIF
C       We form two strings(2-1) aq-q and (4-3) q-qq
C       1  q with XDIQ-XTVQI-XSQ
C       1  q with XDIQ-XTVQI   (j.r.5.5.96)
C       2  aq with XSAQ
C       3  qq with XTVQI+XSQ
C       3  qq with XTVQI     (j.r.5.5.96)
C       4  q with  XQUA-XSAQ
C       with 4-momenta PPP1,PPP2,PPP3,PPP4
        DO 3345 I=1,4
          IF(XDIQ.LE.1.D-15.OR.XQUA.LT.1.D-15)THEN
	    IREJ=1
	    IF(ISQ.EQ.3)IREJ=3
            IF(IPCO.GE.0)WRITE (6,*)' HSDJSE Rejection 3345 XDIQ,XQUA ',
     *	    XDIQ,XQUA
	    RETURN
        ENDIF
C       PPP1(I)=PTA(I)*(XDIQ-XTVQI-XSQ)/XDIQ
        PPP1(I)=PTA(I)*(XDIQ-XTVQI)/XDIQ
C	PPP3(I)=PTA(I)*(XTVQI+XSQ)/XDIQ
	PPP3(I)=PTA(I)*(XTVQI)/XDIQ
	PPP2(I)=PPR(I)*XSAQ/XQUA
	PPP4(I)=PPR(I)*(XQUA-XSAQ)/XQUA
 3345   CONTINUE
 3346   FORMAT(A,5E12.4)
	IF(IPCO.GE.0)THEN
          WRITE(6,3346)' PPR ',PPR
          WRITE(6,3346)' PPP1 ',PPP1
          WRITE(6,3346)' PPP3 ',PPP3
          WRITE(6,3346)' PTA ',PTA
          WRITE(6,3346)' PPP2 XSAQ',PPP2,XSAQ
          WRITE(6,3346)' PPP4 ',PPP4
        ENDIF
C       Invariant Masses of chains
C       ORIGINAL CHAIN
	AMCHOR=SQRT((PPR(4)+PTA(4))**2-(PPR(3)+PTA(3))**2
     *             -(PPR(2)+PTA(2))**2-(PPR(1)+PTA(1))**2) 
	IF(IPCO.GE.0)WRITE(6,2346)'AMCHOR ',AMCHOR
	AMCHN1=SQRT((PPP1(4)+PPP2(4))**2-(PPP1(3)+PPP2(3))**2
     *             -(PPP1(2)+PPP2(2))**2-(PPP1(1)+PPP2(1))**2) 
	IF(IPCO.GE.0)WRITE(6,2346)'AMCHN1 ',AMCHN1
C                    Chain 1 is q-aq restrict mass >1.5 GeV
        IF(AMCHOR.LE.TINY)THEN
          IF(IPCO.GE.0)WRITE(6,2346)' PPR ',PPR
          WRITE(6,2346)' PPP1 ',PPP1
          WRITE(6,2346)' PPP3 ',PPP3
          WRITE(6,2346)' PTA ',PTA
          WRITE(6,2346)' PPP2 ',PPP2
          WRITE(6,2346)' PPP4 ',PPP4
          WRITE(6,2346)'AMCHOR ',AMCHOR
        ENDIF
        IF(AMCHN1.LE.TINY)THEN
          WRITE(6,2346)' PPR ',PPR
          WRITE(6,2346)' PPP1 ',PPP1
          WRITE(6,2346)' PPP3 ',PPP3
          WRITE(6,2346)' PTA ',PTA
          WRITE(6,2346)' PPP2 ',PPP2
          WRITE(6,2346)' PPP4 ',PPP4
          WRITE(6,2346)'AMCHOR ',AMCHOR
        ENDIF
	CHAMAL=0.8D0
	IF(IFB3.GE.3.OR.ISAQ.GE.9)CHAMAL=1.2D0
	IF(AMCHN1.LE.CHAMAL)THEN
          IF(IPCO.GE.0)WRITE(6 ,*)'HADJSE jump2AMCHN1.LE.CHAMAL AMCHOR',
     *	  AMCHN1,CHAMAL,AMCHOR
	  GO TO 3466
	ENDIF
	AMCHN2=SQRT((PPP3(4)+PPP4(4))**2-(PPP3(3)+PPP4(3))**2
     *             -(PPP3(2)+PPP4(2))**2-(PPP3(1)+PPP4(1))**2) 
C                    Chain 2 is qq-q restrict mass >2.5 GeV
	CHAMAL=1.5D0
	IF(IFB2.GE.3.OR.IFB1.GE.3.OR.ISQ.GE.3)CHAMAL=2.0D0
	IF(AMCHN2.LE.CHAMAL)THEN
C	  IREJ=1
C	  RETURN
          IF(IPCO.GE.0)WRITE(6 ,*)' HADJSE jump AMCHN2.LE.CHAMAL',
     *	  AMCHN2,CHAMAL
	  GO TO 3466
	ENDIF
	PXCHK=PPR(1)+PTA(1)-PPP1(1)-PPP2(1)-PPP3(1)-PPP4(1)
	PYCHK=PPR(2)+PTA(2)-PPP1(2)-PPP2(2)-PPP3(2)-PPP4(2)
	PZCHK=PPR(3)+PTA(3)-PPP1(3)-PPP2(3)-PPP3(3)-PPP4(3)
	PECHK=PPR(4)+PTA(4)-PPP1(4)-PPP2(4)-PPP3(4)-PPP4(4)
	IF(IPCO.GE.0)THEN
          WRITE(6,'(A/8E12.4,I5)')
     *	  ' Chain masses AMCH,AMCHOR,AMCHN1,AMCHN2,PZCHK,PECHK,
     *     PXCHK,PYCHK,NOBAM ',
     *	   AMCH,AMCHOR,AMCHN1,AMCHN2,PZCHK,PECHK,PXCHK,PYCHK,NOBAM   
        ENDIF
C       Lorentz parameters of chains
        IF(AMCH.LE.1.D-15)THEN
	  IREJ=1
	    IF(ISQ.EQ.3)IREJ=3
          IF(IPCO.GE.0)WRITE(6,*)' HSDJSE Rejection AMCH too small',AMCH
	  RETURN
	ENDIF
        GAMOR=(PPR(4)+PTA(4))/AMCH
	BGXOR=(PPR(1)+PTA(1))/AMCH
	BGYOR=(PPR(2)+PTA(2))/AMCH
	BGZOR=(PPR(3)+PTA(3))/AMCH
        IF(AMCHN1.LE.1.D-15)THEN
          IF(IPCO.GE.0)WRITE(6,*)' HSDJSE Rejection AMCHN1 too small',
     *	  AMCHN1
	  IREJ=1
	    IF(ISQ.EQ.3)IREJ=3
	  RETURN
	ENDIF
        GAMCH1=(PPP1(4)+PPP2(4))/AMCHN1
        BGXCH1=(PPP1(1)+PPP2(1))/AMCHN1
        BGYCH1=(PPP1(2)+PPP2(2))/AMCHN1
        BGZCH1=(PPP1(3)+PPP2(3))/AMCHN1
        IF(AMCHN2.LE.1.D-15)THEN
	  IREJ=1
	    IF(ISQ.EQ.3)IREJ=3
          IF(IPCO.GE.0)WRITE(6,*)' HSDJSE Rejection AMCHN2 too small',
     *	  AMCHN2
	  RETURN
	ENDIF
        GAMCH2=(PPP3(4)+PPP4(4))/AMCHN2
        BGXCH2=(PPP3(1)+PPP4(1))/AMCHN2
        BGYCH2=(PPP3(2)+PPP4(2))/AMCHN2
        BGZCH2=(PPP3(3)+PPP4(3))/AMCHN2
	IF(IPCO.GE.0)THEN
          WRITE(6,3346)' L.Parm in ',GAM,BGX,BGY,BGZ
          WRITE(6,3346)' L.Parm OR ',GAMOR,BGXOR,BGYOR,BGZOR
          WRITE(6,3346)' L.Parm C1 ',GAMCH1,BGXCH1,BGYCH1,BGZCH1
          WRITE(6,3346)' L.Parm C2 ',GAMCH2,BGXCH2,BGYCH2,BGZCH2
	ENDIF
	PEOR=AMCH*GAMOR
	PZOR=AMCH*BGZOR
	PECH1=AMCHN1*GAMCH1
	PZCH1=AMCHN1*BGZCH1
	PECH2=AMCHN2*GAMCH2
	PZCH2=AMCHN2*BGZCH2
	IF(IPCO.GE.0)THEN
	  WRITE(6,'(A,6E12.4)')' PEOR,PECH1,PECH2,PZOR,PZCH1,PZCH2',
     *	  PEOR,PECH1,PECH2,PZOR,PZCH1,PZCH2
	ENDIF
	NOBA1=3
	NOBA2=4
	IF(IPCO.GE.0)THEN
          WRITE(6,'(A,2I5)')' Jet1 ISAQ,IFB3 ',ISAQ,IFB3
	ENDIF
        CALL HADJET(NHAD1,AMCHN1,PPP2,PPP1,GAMCH1,BGXCH1,BGYCH1,
     *  BGZCH1, ISAQ,IFB3,IFB1,IFB4,I1,I2,NOBA1,NNCH,NORIG)
C       DO 342 I=1,NHAD1
C       WRITE(6,1050)I,PXF(I),PYF(I),PZF(I),HEF(I),AMF(I), ICHF(I),
C    +  IBARF(I),NREF(I),ANF(I)
C 342   CONTINUE
C       DIMENSION ANFF(NFIMAX),PXFF(NFIMAX),PYFF(NFIMAX),PZFF(NFIMAX),
C    +  HEFF(NFIMAX),AMFF(NFIMAX), 
C    *  ICHFF(NFIMAX),IBARFF(NFIMAX),NREFF(NFIMAX)
C                 Intermediate store of hadrons from jet 1
        DO 3348 I=1,NHAD1
	  ANFF(I)=ANF(I)
          PXFF(I)=PXF(I)
          PYFF(I)=PYF(I)
          PZFF(I)=PZF(I)
          HEFF(I)=HEF(I)
          AMFF(I)=AMF(I)
          ICHFF(I)=ICHF(I)
          IBARFF(I)=IBARF(I)
          NREFF(I)=NREF(I)
	  IORMOO(I)=IORMO(I)
 3348   CONTINUE
	IF(IPCO.GE.0)THEN
 	  WRITE(6,'(A,3I5)')' Jet2 IFB1,ISQ,IFB2 ',IFB1,ISQ,IFB2
	ENDIF
        CALL HADJET(NHAD2,AMCHN2,PPP4,PPP3,GAMCH2,BGXCH2,BGYCH2,
     *  BGZCH2, IFB1,ISQ,IFB2,IFB4,I1,I2,NOBA2,NNCH,NORIG)
C       DO 341 I=1,NHAD2
C       WRITE(6,1050)I,PXF(I),PYF(I),PZF(I),HEF(I),AMF(I), ICHF(I),
C    +    IBARF(I),NREF(I),ANF(I)
C       IF(IORMO(I).NE.999)IORMO(I+NHAD1)=IORMO(I)+NHAD1
C 341   CONTINUE
C       Intermediate store of hadrons from jet 2
        DO 3448 I=NHAD1+1,NHAD1+NHAD2
	  II=I-NHAD1
          ANFF(I)=ANF(II)
          PXFF(I)=PXF(II)
          PYFF(I)=PYF(II)
          PZFF(I)=PZF(II)
          HEFF(I)=HEF(II)
          AMFF(I)=AMF(II)
          ICHFF(I)=ICHF(II)
          IBARFF(I)=IBARF(II)
          NREFF(I)=NREF(II)
	  IORMOO(I)=IORMO(II)
	  IF(IORMOO(I).NE.999)IORMOO(I)=IORMOO(I)+NHAD1
 3448   CONTINUE
        NHAD=NHAD1+NHAD2
        DO 3349 I=1,NHAD
	  II=I
          ANF(I)=ANFF(II)
          PXF(I)=PXFF(II)
          PYF(I)=PYFF(II)
          PZF(I)=PZFF(II)
          HEF(I)=HEFF(II)
          AMF(I)=AMFF(II)
          ICHF(I)=ICHFF(II)
          IBARF(I)=IBARFF(II)
          NREF(I)=NREFF(II)
	  IORMO(I)=IORMOO(II)
 3349   CONTINUE
      ENDIF
C==================================================================
      HEFT=0.D0
      IF(IPCO.GE.0)THEN
        DO 40 I=1,NHAD
	  IF(IBARF(I).EQ.500)GO TO 4040
          HEFT=HEFT+HEF(I)
          WRITE(6,4050)I,PXF(I),PYF(I),PZF(I),HEF(I),AMF(I), ICHF(I),
     +    IBARF(I),NREF(I),ANF(I)
 4050     FORMAT(' JET  ',I5,5F12.4,3I5,A10)
 4040     CONTINUE
   40   CONTINUE
      ENDIF
      HEFFF=AMCH*GAM
      IF(IPCO.GE.0)THEN
        WRITE(6,'(A,I5,2E12.4)')' HADJSE 2IREJ,HEFT,HEFFF ',
     *  IREJ,HEFT,HEFFF
      ENDIF
C     IREJ=1
      IF(IREJ.GE.1)RETURN
C
C     WRITE(6,*)' HADJSE: IREJ,ISQ ',IREJ,ISQ
      IF(ISQ.EQ.1)NHSE1=NHSE1+1
      IF(ISQ.EQ.2)NHSE2=NHSE2+1
      IF(ISQ.EQ.3)NHSE3=NHSE3+1
      RETURN
      END

************************************************************************
************************************************************************
************************************************************************
************************************************************************
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HADJASE(NHAD,AMCH,PPR,PTA,GAM,BGX,BGY,BGZ, IFB1,IFB2,
     +IFB3,IFB4,I1,I2,NOBAM,NNCH,NORIG,IREJ,IISSQQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C   HADJASE Fragmentation of vv,vs and sv chains CK method 
C          using Glauber sea quarks
C      This times for aqaq--aq and aq--aqaq chains
C                 j.r.3/99
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
C            = 4  QUARK-DIQUARK JET, FLAVORS : aQU : IFB1, aDIQU :IFB2,IFB3
C                 OR ANTIQUARK-ANTIDIQUARK JET
C
C
C            = 5  DIQUARK-ANTIDIQUARK JET
C                 OR ANTIDIQUARK-DIQUARK JET
C                     FLAVORS :  DIQU : IFB1,IFB2, ANTIDIQU : IFB3,IFB4
C                                IN ANY ORDER
C
C            = 6  DIQUARK-QUARK JET, FLAVORS : DIQU  : IFB1,IFB2 QU: IFB3
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
      CHARACTER*8 ANF,ANFF
      PARAMETER (NFIMAX=249)
      COMMON /DFINPA/ ANF(NFIMAX),PXF(NFIMAX),PYF(NFIMAX),PZF(NFIMAX),
     +HEF(NFIMAX),AMF(NFIMAX), ICHF(NFIMAX),IBARF(NFIMAX),NREF(NFIMAX)
      COMMON /DFINPZ/IORMO(NFIMAX),IDAUG1(NFIMAX),IDAUG2(NFIMAX),
     * ISTATH(NFIMAX)
      DIMENSION ANFF(NFIMAX),PXFF(NFIMAX),PYFF(NFIMAX),PZFF(NFIMAX),
     +HEFF(NFIMAX),AMFF(NFIMAX), 
     *ICHFF(NFIMAX),IBARFF(NFIMAX),NREFF(NFIMAX),IORMOO(NFIMAX)
        CHARACTER*80 TITLE
        CHARACTER*8 PROJTY,TARGTY
C       COMMON /USER/TITLE,PROJTY,TARGTY,CMENER,ISTRUF
C    & ,ISINGD,IDUBLD,SDFRAC,PTLAR
      COMMON /USER1/TITLE,PROJTY,TARGTY
      COMMON /USER2/CMENER,SDFRAC,PTLAR,ISTRUF,ISINGD,IDUBLD
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
      DIMENSION PPR(4),PTA(4)
*KEEP,XSEADI.
      COMMON /XSEADI/ XSEACU,UNON,UNOM,UNOSEA,CVQ,CDQ,CSEA,SSMIMA,
     +SSMIMQ,VVMTHR
      COMMON /HDJASE/NHSE1,NHSE2,NHSE3,NHASE1,NHASE2,NHASE3
      DIMENSION PPP1(4),PPP2(4),PPP3(4),PPP4(4)
      DATA TINY/1.0D-3/
C----------------------------------------------------------------------
C
C	WRITE(6,'(A,3I5)')' Jet0 IFB1,IFB2,IFB3 ',IFB1,IFB2,IFB3
      IF(IPCO.GE.0)WRITE(6,*)
     *' HADJASE(NHAD,AMCH,PPR,PTA,GAM,BGX,BGY,BGZ, IFB1,IFB2,',
     +'IFB3,IFB4,I1,I2,NOBAM,NNCH,NORIG,IREJ),IPCO',
     *NHAD,AMCH,PPR,PTA,GAM,BGX,BGY,BGZ, IFB1,IFB2,
     +IFB3,IFB4,I1,I2,NOBAM,NNCH,NORIG,IREJ,IPCO
      IREJ=0
      IF(IPCO.GE.0) THEN
    	WRITE(6,'(A,3I5)')' HADJASE Jet0 IFB1,IFB2,IFB3',IFB1,IFB2,IFB3
        WRITE(6,'(4E12.4)') GAM,BGX,BGY,BGZ,PPR
        WRITE(6,1000) NHAD,AMCH,PTA, IFB1,IFB2,IFB3,IFB4,I1,I2,NOBAM,
     +  NNCH,NORIG
 1000   FORMAT(10X,I10,5F10.3/10X,9I10)
      ENDIF
      IF(ABS(NNCH).EQ.99) THEN
        NHAD=0
C       IPCO=0
        RETURN
      ENDIF
C=================================================================
C=================================================================
C                   DIQUARK-QUARK NOBAM=6
C=================================================================
C=================================================================
C              Different options:
      IF(NOBAM.EQ.6)THEN
        IF(IPCO.GE.0)WRITE(6,*)' DIQUARK-QUARK NOBAM=6'
C              NOBAM=6:  Diquark--Quark Jet
C
C              Work out approximate x-Values of diquark and quark
	IF(CMENER.LT.1.D-3)THEN
	  WRITE(6,*)' CMENER=0. HADJASE',CMENER
	  STOP
	ENDIF
	XDIQ=PPR(4)*2.D0/CMENER
	XQUA=PTA(4)*2.D0/CMENER
        IF(IPCO.GE.0) THEN
          WRITE(6,'(A,2E12.4)')' HADJASE:XDIQ,XQUA ',XDIQ,XQUA
        ENDIF
C       Select sea quark x for one of the diquark quarks
C       Select x values of sea quark pair
	ICOU=0
 1234   CONTINUE
	ICOU=ICOU+1
	IF(ICOU.GE.500)THEN
	  IREJ=1
          IF(IPCO.GE.0)WRITE(6,*)' HADJASE reject icou 100'
	  RETURN
	ENDIF
	CALL XSEAPA(CMENER,XQUA/2.D0,ISQ,ISAQ,XSQ,XSAQ,IREJ)
	IISSQQ=ISQ
	IF(IREJ.GE.1)THEN
          IF(IPCO.GE.0)WRITE(6,*)' HADJASE reject XSEAPA'
	  RETURN
	ENDIF
	IF(XSAQ.GE.2.D0*XQUA/3.D0)GO TO 1234
        IF(IPCO.GE.0) THEN
          WRITE(6,*)' HADJASE,XSEAPA:XSQ,XSAQ,ISQ,ISAQ ',
     *	  XSQ,XSAQ,ISQ,ISAQ
        ENDIF
*       projectile valence quarks xpvqi
C       here selected as sea quark! from 1/x-sea
        XVTHRO=CVQ/CMENER
	IVTHR=0
 3465   CONTINUE
	IF(IVTHR.EQ.200)THEN
	  IREJ=1
	  IF(ISQ.EQ.3)IREJ=3
          IF(IPCO.GE.0)WRITE(6,*)' HADJASE 3465 reject IVTHR 50'
	  RETURN
        ENDIF
	IVTHR=IVTHR+1
        XVTHR=XVTHRO/(201-IVTHR)
        UNOPRV=UNON
   80   CONTINUE
        IF(XVTHR.GT.0.66D0*XDIQ)THEN
          IF(IPCO.GE.0)WRITE(6,*)' HADJASE 80 reject XVTHR too great',
     *	  XVTHR
          IREJ=1
	  IF(ISQ.EQ.3)IREJ=3
	  RETURN
	ENDIF
        XPVQI=SAMPEY(XVTHR,0.66D0*XDIQ)
        XDIQQ=XDIQ-XPVQI
        IF(IPCO.GE.0) THEN
          WRITE(6,'(A,2E12.4)')' HADJASE:XDIQQ,XPVQI ',XDIQQ,XPVQI
        ENDIF
C
C       We form two strings(1-2) aq-q and (3-4) aqaq-aq
C       1  aq with XDIQ-XPVQI   (j.r.5.5.96)
C       2  q with XSAQ
C       3  aqaq with XPVQI       (j.r.5.5.96)
C       4  aq with  XQUA-XSAQ
C       with 4-momenta PPP1,PPP2,PPP3,PPP4
        DO 2345 I=1,4
          IF(XDIQ.LE.1.D-15.OR.XQUA.LT.1.D-15)THEN
            IREJ=1
	    IF(ISQ.EQ.3)IREJ=3
            IF(IPCO.GE.0)WRITE(6,*)' HADJASE reject 2345 XDIQ,',
     *	    'XQUA too small ',
     *	    XDIQ,XQUA
            RETURN
          ENDIF
          PPP1(I)=PPR(I)*(XDIQ-XPVQI)/XDIQ
          PPP3(I)=PPR(I)*(XPVQI)/XDIQ
          PPP2(I)=PTA(I)*XSAQ/XQUA
          PPP4(I)=PTA(I)*(XQUA-XSAQ)/XQUA
 2345   CONTINUE
 2346   FORMAT(A,5E12.4)
        IF(IPCO.GE.0) THEN
          WRITE(6,2346)' PPR ',PPR
          WRITE(6,2346)' PPP1 ',PPP1
          WRITE(6,2346)' PPP3 ',PPP3
          WRITE(6,2346)' PTA ',PTA
          WRITE(6,2346)' PPP2 XSAQ',PPP2,XSAQ
          WRITE(6,2346)' PPP4 ',PPP4
        ENDIF
C       Invariant Masses of chains
C       ORIGINAL CHAIN
        AMCHOR=SQRT((PPR(4)+PTA(4))**2-(PPR(3)+PTA(3))**2
     *             -(PPR(2)+PTA(2))**2-(PPR(1)+PTA(1))**2) 
        AMCHN1=SQRT((PPP1(4)+PPP2(4))**2-(PPP1(3)+PPP2(3))**2
     *             -(PPP1(2)+PPP2(2))**2-(PPP1(1)+PPP2(1))**2) 
        IF(AMCHOR.LE.TINY)THEN
          WRITE(6,2346)' PPR ',PPR
          WRITE(6,2346)' PPP1 ',PPP1
          WRITE(6,2346)' PPP3 ',PPP3
          WRITE(6,2346)' PTA ',PTA
          WRITE(6,2346)' PPP2 ',PPP2
          WRITE(6,2346)' PPP4 ',PPP4
        ENDIF
        IF(AMCHN1.LE.TINY)THEN
          WRITE(6,2346)' PPR ',PPR
          WRITE(6,2346)' PPP1 ',PPP1
          WRITE(6,2346)' PPP3 ',PPP3
          WRITE(6,2346)' PTA ',PTA
          WRITE(6,2346)' PPP2 ',PPP2
          WRITE(6,2346)' PPP4 ',PPP4
        ENDIF
C       Chain 1 is q-aq restrict mass >1.5 GeV 
        CHAMAL=0.8D0
        IF(IFB1.GE.9.OR.ISQ.GE.3)CHAMAL=1.1D0
        IF(AMCHN1.LE.CHAMAL)THEN
          IF(IPCO.GE.0)WRITE (6,*)'HADJASE jump1AMCHN1.LE.CHAMAL AMCHOR'
     *	  ,AMCHN1,CHAMAL,AMCHOR
          GO TO 3465
        ENDIF
        AMCHN2=SQRT((PPP3(4)+PPP4(4))**2-(PPP3(3)+PPP4(3))**2
     *             -(PPP3(2)+PPP4(2))**2-(PPP3(1)+PPP4(1))**2) 
C                    Chain 2 is qq-q restrict mass >2.5 GeV
	CHAMAL=1.3D0
	IF(IFB2.GE.9.OR.IFB3.GE.9.OR.ISAQ.GE.9)CHAMAL=1.80D0
	IF(AMCHN2.LE.CHAMAL)THEN
C	  IREJ=1
C 	  RETURN
          IF(IPCO.GE.0)WRITE (6,*)'HADJASE jump AMCHN2.LE.CHAMAL ',
     *	  AMCHN2,CHAMAL
	  GO TO 3465
        ENDIF
        PXCHK=PPR(1)+PTA(1)-PPP1(1)-PPP2(1)-PPP3(1)-PPP4(1)
        PYCHK=PPR(2)+PTA(2)-PPP1(2)-PPP2(2)-PPP3(2)-PPP4(2)
        PZCHK=PPR(3)+PTA(3)-PPP1(3)-PPP2(3)-PPP3(3)-PPP4(3)
        PECHK=PPR(4)+PTA(4)-PPP1(4)-PPP2(4)-PPP3(4)-PPP4(4)
	IF(IPCO.GE.0)THEN
          WRITE(6,'(A/8E12.4,I5)')
     *	  ' Chain masses AMCH,AMCHOR,AMCHN1,AMCHN2,PZCHK,PECHK,
     *    PXCHK,PYCHK,NOBAM ',
     *	  AMCH,AMCHOR,AMCHN1,AMCHN2,PZCHK,PECHK,PXCHK,PYCHK,NOBAM   
        ENDIF
C       Lorentz parameters of chains
        IF(AMCH.LE.1.D-15)THEN
	  IREJ=1
	    IF(ISQ.EQ.3)IREJ=3
          IF(IPCO.GE.0)WRITE(6,*)' HADJASE rejection AMCH too small ',
     *	  AMCH
	  RETURN
        ENDIF
        GAMOR=(PPR(4)+PTA(4))/AMCH
        BGXOR=(PPR(1)+PTA(1))/AMCH
        BGYOR=(PPR(2)+PTA(2))/AMCH
        BGZOR=(PPR(3)+PTA(3))/AMCH
        IF(AMCHN1.LE.1.D-15)THEN
          IF(IPCO.GE.0)WRITE(6,*)' HADJASE rejection AMCHN1 too small ',
     *	  AMCHN1
	  IREJ=1
	    IF(ISQ.EQ.3)IREJ=3
	  RETURN
	ENDIF
        GAMCH1=(PPP1(4)+PPP2(4))/AMCHN1
        BGXCH1=(PPP1(1)+PPP2(1))/AMCHN1
        BGYCH1=(PPP1(2)+PPP2(2))/AMCHN1
        BGZCH1=(PPP1(3)+PPP2(3))/AMCHN1
        IF(AMCHN2.LE.1.D-15)THEN
          IF(IPCO.GE.0)WRITE(6,*)' HADJASE rejection AMCHN2 too small ',
     *	  AMCHN2
	  IREJ=1
	    IF(ISQ.EQ.3)IREJ=3
	  RETURN
	ENDIF
        GAMCH2=(PPP3(4)+PPP4(4))/AMCHN2
        BGXCH2=(PPP3(1)+PPP4(1))/AMCHN2
        BGYCH2=(PPP3(2)+PPP4(2))/AMCHN2
        BGZCH2=(PPP3(3)+PPP4(3))/AMCHN2
	IF(IPCO.GE.0)THEN
 	  WRITE(6,2346)' L.Parm in ',GAM,BGX,BGY,BGZ
 	  WRITE(6,2346)' L.Parm OR ',GAMOR,BGXOR,BGYOR,BGZOR
 	  WRITE(6,2346)' L.Parm C1 ',GAMCH1,BGXCH1,BGYCH1,BGZCH1
 	  WRITE(6,2346)' L.Parm C2 ',GAMCH2,BGXCH2,BGYCH2,BGZCH2
	ENDIF
	PEOR=AMCH*GAMOR
	PZOR=AMCH*BGZOR
	PECH1=AMCHN1*GAMCH1
	PZCH1=AMCHN1*BGZCH1
	PECH2=AMCHN2*GAMCH2
	PZCH2=AMCHN2*BGZCH2
	IF(IPCO.GE.0)THEN
	  WRITE(6,'(A,6E12.4)')' PEOR,PECH1,PECH2,PZOR,PZCH1,PZCH2',
     *	  PEOR,PECH1,PECH2,PZOR,PZCH1,PZCH2
	ENDIF
	NOBA1=3
	NOBA2=6
	IF(IPCO.GE.-1)THEN
C	  WRITE(6,'(A,2I5)')' Jet1 IFB1,ISQ ',IFB1,ISQ
	ENDIF
        CALL HADJET(NHAD1,AMCHN1,PPP1,PPP2,GAMCH1,BGXCH1,BGYCH1,
     *  BGZCH1, IFB1,ISQ,IFB3,IFB4,I1,I2,NOBA1,NNCH,NORIG)
C       DO 42 I=1,NHAD1
C         WRITE(6,1050)I,PXF(I),PYF(I),PZF(I),HEF(I),AMF(I), ICHF(I),
C    +    IBARF(I),NREF(I),ANF(I)
C  42   CONTINUE
	IF(IPCO.GE.-1)THEN
C	  WRITE(6,'(A,2I5)')' Jet1 NHAD1 ',NHAD1
	ENDIF
C       DIMENSION ANFF(NFIMAX),PXFF(NFIMAX),PYFF(NFIMAX),PZFF(NFIMAX),
C    +  HEFF(NFIMAX),AMFF(NFIMAX), 
C    *  ICHFF(NFIMAX),IBARFF(NFIMAX),NREFF(NFIMAX)
C       Intermediate store of hadrons from jet 1
        DO 2348 I=1,NHAD1
          ANFF(I)=ANF(I)
          PXFF(I)=PXF(I)
          PYFF(I)=PYF(I)
          PZFF(I)=PZF(I)
          HEFF(I)=HEF(I)
          AMFF(I)=AMF(I)
          ICHFF(I)=ICHF(I)
          IBARFF(I)=IBARF(I)
          NREFF(I)=NREF(I)
	  IORMOO(I)=IORMO(I)
 2348   CONTINUE
	IF(IPCO.GE.0)THEN
          WRITE(6,'(A,3I5)')' Jet2 IFB2,ISAQ,IFB3 ',IFB2,ISAQ,IFB3
	ENDIF
        CALL HADJET(NHAD2,AMCHN2,PPP3,PPP4,GAMCH2,BGXCH2,BGYCH2,
     *  BGZCH2, IFB2,ISAQ,IFB3,IFB4,I1,I2,NOBA2,NNCH,NORIG)
	IF(IPCO.GE.-1)THEN
C	  WRITE(6,'(A,2I5)')' Jet2 NHAD2 ',NHAD2
	ENDIF
C       DO 41 I=1,NHAD2
C         WRITE(6,1050)I,PXF(I),PYF(I),PZF(I),HEF(I),AMF(I), ICHF(I),
C    +    IBARF(I),NREF(I),ANF(I)
C        IF(IORMO(I).NE.999)IORMO(I+NHAD1)=IORMO(I)+NHAD1
C  41   CONTINUE
C       Intermediate store of hadrons from jet 2
        DO 2448 I=NHAD1+1,NHAD1+NHAD2
	  II=I-NHAD1
          ANFF(I)=ANF(II)
          PXFF(I)=PXF(II)
          PYFF(I)=PYF(II)
          PZFF(I)=PZF(II)
          HEFF(I)=HEF(II)
          AMFF(I)=AMF(II)
          ICHFF(I)=ICHF(II)
          IBARFF(I)=IBARF(II)
          NREFF(I)=NREF(II)
	  IORMOO(I)=IORMO(II)
	  IF(IORMOO(I).NE.999)IORMOO(I)=IORMOO(I)+NHAD1
 2448   CONTINUE
        NHAD=NHAD1+NHAD2
        DO 2349 I=1,NHAD
	  II=I
	  ANF(I)=ANFF(II)
	  PXF(I)=PXFF(II)
	  PYF(I)=PYFF(II)
	  PZF(I)=PZFF(II)
	  HEF(I)=HEFF(II)
	  AMF(I)=AMFF(II)
	  ICHF(I)=ICHFF(II)
	  IBARF(I)=IBARFF(II)
	  NREF(I)=NREFF(II)
	  IORMO(I)=IORMOO(II)
 2349   CONTINUE
      ENDIF
C=================================================================
C=================================================================
C                   QUARK-DIQUARK NOBAM=4
C=================================================================
C==================================================================
C     Different options:
      IF(NOBAM.EQ.4)THEN
        IF(IPCO.GE.0)WRITE(6,*)' AQUARK-ADIQUARK NOBAM=4'
C       NOBAM=4:  AQuark-ADIQUARK Jet
C
C       Work out approximate x-Values of quark and diquark
	IF(CMENER.LT.1.D-3)THEN
	  IF(IPCO.GE.0)WRITE(6,*)' CMENER=0. HADJASE',CMENER
	  STOP
	ENDIF
	XDIQ=PTA(4)*2.D0/CMENER
	XQUA=PPR(4)*2.D0/CMENER
	IF(IPCO.GE.0)THEN
          WRITE(6,'(A,2E12.4)')' HADJASE:XDIQ,XQUA ',XDIQ,XQUA
	ENDIF
C       Select valence quark x for one of the diquark quarks
C       Select x values of sea quark pair
	ICOU=0
 2234   CONTINUE
	ICOU=ICOU+1
	IF(ICOU.GE.500)THEN
	  IREJ=1
	    IF(ISQ.EQ.3)IREJ=3
          IF(IPCO.GE.0)WRITE(6,*)' HADJASE Rejection 2234 ICOU. GT.100'
	  RETURN
	ENDIF
	IF(IPCO.GE.0)WRITE(6,*)' XSEAPA: CMENER,XQUA ',CMENER,XQUA
	CALL XSEAPA(CMENER,XQUA/2.D0,ISQ,ISAQ,XSQ,XSAQ,IREJ)
	IISSQQ=ISQ
	IF(IREJ.GE.1)THEN
          IF(IPCO.GE.0)WRITE(6,*)' HADJASE reject XSEAPA'
	  RETURN
	ENDIF
	IF(XSAQ.GE.2.D0*XQUA/3.D0)GO TO 2234
	IF(IPCO.GE.0)THEN
          WRITE(6,'(A,2E12.4)')' HADJASE:XSQ,XSAQ ',XSQ,XSAQ
        ENDIF
*       target valence quarks xtvqi
C       here selected as sea quark! from 1/x-sea
        XVTHRO=CVQ/CMENER
	IVTHR=0
 3466   CONTINUE
	IF(IVTHR.EQ.200)THEN
	  IREJ=1
	    IF(ISQ.EQ.3)IREJ=3
          IF(IPCO.GE.0)WRITE(6,*)' HADJASE 3466 reject IVTHR 50'
	  RETURN
        ENDIF
	IVTHR=IVTHR+1
        XVTHR=XVTHRO/(201-IVTHR)
        UNOPRV=UNON
  380   CONTINUE
        IF(XVTHR.GT.0.66D0*XDIQ)THEN
          IREJ=1
	    IF(ISQ.EQ.3)IREJ=3
          IF(IPCO.GE.0)WRITE(6,*)' HADJASE Rejection 380 XVTHR  large ',
     *	  XVTHR
          RETURN
        ENDIF
        XTVQI=SAMPEY(XVTHR,0.66D0*XDIQ)
	XDIQQ=XDIQ-XTVQI
	IF(IPCO.GE.0)THEN
          WRITE(6,'(A,2E12.4)')' HADJCK:XDIQQ,XTVQI ',XDIQQ,XTVQI
        ENDIF
C       We form two strings(2-1) q-aq and (4-3) aq-aqaq
C       1  aq with XDIQ-XTVQI   (j.r.5.5.96)
C       2  q with XSAQ
C       3  aqaq with XTVQI     (j.r.5.5.96)
C       4  aq with  XQUA-XSAQ
C       with 4-momenta PPP1,PPP2,PPP3,PPP4
        DO 3345 I=1,4
          IF(XDIQ.LE.1.D-15.OR.XQUA.LT.1.D-15)THEN
	    IREJ=1
	    IF(ISQ.EQ.3)IREJ=3
            IF(IPCO.GE.0)WRITE (6,*)' HSDJSE Rejection 3345 XDIQ,XQUA ',
     *	    XDIQ,XQUA
	    RETURN
        ENDIF
        PPP1(I)=PTA(I)*(XDIQ-XTVQI)/XDIQ
	PPP3(I)=PTA(I)*(XTVQI)/XDIQ
	PPP2(I)=PPR(I)*XSAQ/XQUA
	PPP4(I)=PPR(I)*(XQUA-XSAQ)/XQUA
 3345   CONTINUE
 3346   FORMAT(A,5E12.4)
	IF(IPCO.GE.0)THEN
          WRITE(6,3346)' PPR ',PPR
          WRITE(6,3346)' PPP1 ',PPP1
          WRITE(6,3346)' PPP3 ',PPP3
          WRITE(6,3346)' PTA ',PTA
          WRITE(6,3346)' PPP2 XSAQ',PPP2,XSAQ
          WRITE(6,3346)' PPP4 ',PPP4
        ENDIF
C       Invariant Masses of chains
C       ORIGINAL CHAIN
	AMCHOR=SQRT((PPR(4)+PTA(4))**2-(PPR(3)+PTA(3))**2
     *             -(PPR(2)+PTA(2))**2-(PPR(1)+PTA(1))**2) 
	IF(IPCO.GE.0)WRITE(6,2346)'AMCHOR ',AMCHOR
	AMCHN1=SQRT((PPP1(4)+PPP2(4))**2-(PPP1(3)+PPP2(3))**2
     *             -(PPP1(2)+PPP2(2))**2-(PPP1(1)+PPP2(1))**2) 
	IF(IPCO.GE.0)WRITE(6,2346)'AMCHN1 ',AMCHN1
C                    Chain 1 is q-aq restrict mass >1.5 GeV
        IF(AMCHOR.LE.TINY)THEN
          IF(IPCO.GE.0)WRITE(6,2346)' PPR ',PPR
          WRITE(6,2346)' PPP1 ',PPP1
          WRITE(6,2346)' PPP3 ',PPP3
          WRITE(6,2346)' PTA ',PTA
          WRITE(6,2346)' PPP2 ',PPP2
          WRITE(6,2346)' PPP4 ',PPP4
          WRITE(6,2346)'AMCHOR ',AMCHOR
        ENDIF
        IF(AMCHN1.LE.TINY)THEN
          WRITE(6,2346)' PPR ',PPR
          WRITE(6,2346)' PPP1 ',PPP1
          WRITE(6,2346)' PPP3 ',PPP3
          WRITE(6,2346)' PTA ',PTA
          WRITE(6,2346)' PPP2 ',PPP2
          WRITE(6,2346)' PPP4 ',PPP4
          WRITE(6,2346)'AMCHOR ',AMCHOR
        ENDIF
	CHAMAL=0.8D0
	IF(IFB3.GE.9.OR.ISQ.GE.3)CHAMAL=1.2D0
	IF(AMCHN1.LE.CHAMAL)THEN
          IF(IPCO.GE.0)WRITE(6 ,*)'HADJASE jump2AMCHN1.LE.CHAMAL AMCHOR'
     *	  ,AMCHN1,CHAMAL,AMCHOR
	  GO TO 3466
	ENDIF
	AMCHN2=SQRT((PPP3(4)+PPP4(4))**2-(PPP3(3)+PPP4(3))**2
     *             -(PPP3(2)+PPP4(2))**2-(PPP3(1)+PPP4(1))**2) 
C                    Chain 2 is qq-q restrict mass >2.5 GeV
	CHAMAL=1.4D0
	IF(IFB2.GE.9.OR.IFB1.GE.9.OR.ISAQ.GE.9)CHAMAL=1.8D0
	IF(AMCHN2.LE.CHAMAL)THEN
C	  IREJ=1
C	  RETURN
          IF(IPCO.GE.0)WRITE(6 ,*)' HADJASE jump AMCHN2.LE.CHAMAL',
     *	  AMCHN2,CHAMAL
	  GO TO 3466
	ENDIF
	PXCHK=PPR(1)+PTA(1)-PPP1(1)-PPP2(1)-PPP3(1)-PPP4(1)
	PYCHK=PPR(2)+PTA(2)-PPP1(2)-PPP2(2)-PPP3(2)-PPP4(2)
	PZCHK=PPR(3)+PTA(3)-PPP1(3)-PPP2(3)-PPP3(3)-PPP4(3)
	PECHK=PPR(4)+PTA(4)-PPP1(4)-PPP2(4)-PPP3(4)-PPP4(4)
	IF(IPCO.GE.0)THEN
          WRITE(6,'(A/8E12.4,I5)')
     *	  ' Chain masses AMCH,AMCHOR,AMCHN1,AMCHN2,PZCHK,PECHK,
     *     PXCHK,PYCHK,NOBAM ',
     *	   AMCH,AMCHOR,AMCHN1,AMCHN2,PZCHK,PECHK,PXCHK,PYCHK,NOBAM   
        ENDIF
C       Lorentz parameters of chains
        IF(AMCH.LE.1.D-15)THEN
	  IREJ=1
	    IF(ISQ.EQ.3)IREJ=3
          IF(IPCO.GE.0)WRITE(6,*)' HSDJSE Rejection AMCH too small',AMCH
	  RETURN
	ENDIF
        GAMOR=(PPR(4)+PTA(4))/AMCH
	BGXOR=(PPR(1)+PTA(1))/AMCH
	BGYOR=(PPR(2)+PTA(2))/AMCH
	BGZOR=(PPR(3)+PTA(3))/AMCH
        IF(AMCHN1.LE.1.D-15)THEN
          IF(IPCO.GE.0)WRITE(6,*)' HSDJSE Rejection AMCHN1 too small',
     *	  AMCHN1
	  IREJ=1
	    IF(ISQ.EQ.3)IREJ=3
	  RETURN
	ENDIF
        GAMCH1=(PPP1(4)+PPP2(4))/AMCHN1
        BGXCH1=(PPP1(1)+PPP2(1))/AMCHN1
        BGYCH1=(PPP1(2)+PPP2(2))/AMCHN1
        BGZCH1=(PPP1(3)+PPP2(3))/AMCHN1
        IF(AMCHN2.LE.1.D-15)THEN
	  IREJ=1
	    IF(ISQ.EQ.3)IREJ=3
          IF(IPCO.GE.0)WRITE(6,*)' HSDJSE Rejection AMCHN2 too small',
     *	  AMCHN2
	  RETURN
	ENDIF
        GAMCH2=(PPP3(4)+PPP4(4))/AMCHN2
        BGXCH2=(PPP3(1)+PPP4(1))/AMCHN2
        BGYCH2=(PPP3(2)+PPP4(2))/AMCHN2
        BGZCH2=(PPP3(3)+PPP4(3))/AMCHN2
	IF(IPCO.GE.0)THEN
          WRITE(6,3346)' L.Parm in ',GAM,BGX,BGY,BGZ
          WRITE(6,3346)' L.Parm OR ',GAMOR,BGXOR,BGYOR,BGZOR
          WRITE(6,3346)' L.Parm C1 ',GAMCH1,BGXCH1,BGYCH1,BGZCH1
          WRITE(6,3346)' L.Parm C2 ',GAMCH2,BGXCH2,BGYCH2,BGZCH2
	ENDIF
	PEOR=AMCH*GAMOR
	PZOR=AMCH*BGZOR
	PECH1=AMCHN1*GAMCH1
	PZCH1=AMCHN1*BGZCH1
	PECH2=AMCHN2*GAMCH2
	PZCH2=AMCHN2*BGZCH2
	IF(IPCO.GE.0)THEN
	  WRITE(6,'(A,6E12.4)')' PEOR,PECH1,PECH2,PZOR,PZCH1,PZCH2',
     *	  PEOR,PECH1,PECH2,PZOR,PZCH1,PZCH2
	ENDIF
	NOBA1=3
	NOBA2=4
	IF(IPCO.GE.0)THEN
          WRITE(6,'(A,2I5)')' Jet1 ISQ,IFB3 ',ISQ,IFB3
	ENDIF
        CALL HADJET(NHAD1,AMCHN1,PPP2,PPP1,GAMCH1,BGXCH1,BGYCH1,
     *  BGZCH1, ISQ,IFB3,IFB1,IFB4,I1,I2,NOBA1,NNCH,NORIG)
C       DO 342 I=1,NHAD1
C       WRITE(6,1050)I,PXF(I),PYF(I),PZF(I),HEF(I),AMF(I), ICHF(I),
C    +  IBARF(I),NREF(I),ANF(I)
C 342   CONTINUE
C       DIMENSION ANFF(NFIMAX),PXFF(NFIMAX),PYFF(NFIMAX),PZFF(NFIMAX),
C    +  HEFF(NFIMAX),AMFF(NFIMAX), 
C    *  ICHFF(NFIMAX),IBARFF(NFIMAX),NREFF(NFIMAX)
C                 Intermediate store of hadrons from jet 1
        DO 3348 I=1,NHAD1
	  ANFF(I)=ANF(I)
          PXFF(I)=PXF(I)
          PYFF(I)=PYF(I)
          PZFF(I)=PZF(I)
          HEFF(I)=HEF(I)
          AMFF(I)=AMF(I)
          ICHFF(I)=ICHF(I)
          IBARFF(I)=IBARF(I)
          NREFF(I)=NREF(I)
	  IORMOO(I)=IORMO(I)
 3348   CONTINUE
	IF(IPCO.GE.0)THEN
 	  WRITE(6,'(A,3I5)')' Jet2 IFB1,ISAQ,IFB2 ',IFB1,ISAQ,IFB2
	ENDIF
        CALL HADJET(NHAD2,AMCHN2,PPP4,PPP3,GAMCH2,BGXCH2,BGYCH2,
     *  BGZCH2, IFB1,ISAQ,IFB2,IFB4,I1,I2,NOBA2,NNCH,NORIG)
C       DO 341 I=1,NHAD2
C       WRITE(6,1050)I,PXF(I),PYF(I),PZF(I),HEF(I),AMF(I), ICHF(I),
C    +    IBARF(I),NREF(I),ANF(I)
C       IF(IORMO(I).NE.999)IORMO(I+NHAD1)=IORMO(I)+NHAD1
C 341   CONTINUE
C       Intermediate store of hadrons from jet 2
        DO 3448 I=NHAD1+1,NHAD1+NHAD2
	  II=I-NHAD1
          ANFF(I)=ANF(II)
          PXFF(I)=PXF(II)
          PYFF(I)=PYF(II)
          PZFF(I)=PZF(II)
          HEFF(I)=HEF(II)
          AMFF(I)=AMF(II)
          ICHFF(I)=ICHF(II)
          IBARFF(I)=IBARF(II)
          NREFF(I)=NREF(II)
	  IORMOO(I)=IORMO(II)
	  IF(IORMOO(I).NE.999)IORMOO(I)=IORMOO(I)+NHAD1
 3448   CONTINUE
        NHAD=NHAD1+NHAD2
        DO 3349 I=1,NHAD
	  II=I
          ANF(I)=ANFF(II)
          PXF(I)=PXFF(II)
          PYF(I)=PYFF(II)
          PZF(I)=PZFF(II)
          HEF(I)=HEFF(II)
          AMF(I)=AMFF(II)
          ICHF(I)=ICHFF(II)
          IBARF(I)=IBARFF(II)
          NREF(I)=NREFF(II)
	  IORMO(I)=IORMOO(II)
 3349   CONTINUE
      ENDIF
C==================================================================
      HEFT=0.D0
      IF(IPCO.GE.0)THEN
        DO 40 I=1,NHAD
	  IF(IBARF(I).EQ.500)GO TO 4040
          HEFT=HEFT+HEF(I)
          WRITE(6,4050)I,PXF(I),PYF(I),PZF(I),HEF(I),AMF(I), ICHF(I),
     +    IBARF(I),NREF(I),ANF(I)
 4050     FORMAT(' JET  ',I5,5F12.4,3I5,A10)
 4040     CONTINUE
   40   CONTINUE
      ENDIF
      HEFFF=AMCH*GAM
      IF(IPCO.GE.0)THEN
        WRITE(6,'(A,I5,2E12.4)')'3IREJ,HEFT,HEFFF ',
     *  IREJ,HEFT,HEFFF
      ENDIF
C     IREJ=1
      IF(IREJ.GE.1)RETURN
C
C     WRITE(6,*)' HADJASE: IREJ,ISQ ',IREJ,ISQ
      IF(ISQ.EQ.1)NHASE1=NHASE1+1
      IF(ISQ.EQ.2)NHASE2=NHASE2+1
      IF(ISQ.EQ.3)NHASE3=NHASE3+1
      RETURN
      END

************************************************************************
************************************************************************
