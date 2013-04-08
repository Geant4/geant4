*
*===sdiff================================================================*
*
      SUBROUTINE SDIFF(EPROJ,PPROJ,KPROJ,NHKKH1,IQQDD)
 
**************************************************************************
* Version November 1993                         by Stefan Roesler        *
*                                                  University of Leipzig *
* This subroutine calls one single-diffractive event depending on the    *
* single-diffractive cross section for a given hadron-hadron interaction.*
**************************************************************************
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)
      PARAMETER (NMXHKK= 89998,INTMX=2488)
      CHARACTER*8 ANAME,ANC
      COMMON /DIQI/ IPVQ(248),   IPPV1(248),   IPPV2(248), ITVQ(248),
     &              ITTV1(248),  ITTV2(248),   IPSQ(INTMX),IPSQ2(INTMX),
     &              IPSAQ(INTMX),IPSAQ2(INTMX),ITSQ(INTMX),ITSQ2(INTMX),
     &              ITSAQ(INTMX),ITSAQ2(INTMX),KKPROJ(248),KKTARG(248)
      COMMON /HKKEVT/ NHKK,NEVHKK,     ISTHKK(NMXHKK),  IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),PHKK(5,NMXHKK),
     &                VHKK(4,NMXHKK),  WHKK(4,NMXHKK)
      COMMON /DIFFRA/ ISINGD,IDIFTP,IOUDIF,IFLAGD
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
      COMMON /CMHICO/ CMHIS
      COMMON /DPRIN/ IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
      COMMON /DPAR/ ANAME(210),AAM(210),GA(210),TAU(210),IICH(210),
     &              IIBAR(210),K1(210),K2(210)
      COMMON /DIFPAR/ PXC(902),  PYC(902),PZC(902),
     &                HEC(902),  AMC(902),ICHC(902),
     &                IBARC(902),ANC(902),NRC(902)
*----------------- S. Roesler 07/07/93
C     COMMON /COUNTEV/NCDIFF,NCSDIF
*---------- S. Roesler 5-11-93
      COMMON /IFRAGM/IFRAG
      COMMON/XXLMDD/IJLMDD,KDLMDD 
      COMMON /NUCCMS/GAMV,BGLV,DUMMY(7)

*
 
      DATA NCIREJ /0/
C     DATA NCDIFF /0/
      DATA SIGDIF,SIGDIH /6.788790702D0,3.631283998D0/
*
      IJLMDD= 0
      IREJ  = 0
      KTARG = 0
      PPROX = 0.0D0
      PPROY = 0.0D0
      PPROL = PPROJ
      EPRO  = EPROJ
      AMPRO = AAM(KPROJ)
      IF(IPRI.GE.1)WRITE(6,'(A,2E20.8)')' SDIFF:EPROJ,PPROJ ',
     * EPROJ,PPROJ
C     WRITE(6,'(A,2E20.8)')' SDIFF:EPROJ,PPROJ ',
C    * EPROJ,PPROJ
*
      DO 10 I=1,IT
         IHKK = I+IP
         IF (ISTHKK(IHKK).EQ.12) THEN
            KTARG = KKTARG(I)
            PTARX = PHKK(1,IHKK)
            PTARY = PHKK(2,IHKK)
            PTARL = PHKK(3,IHKK)
            ETAR  = PHKK(4,IHKK)
            AMTAR = PHKK(5,IHKK)
            PTOTX = PPROX + PTARX
            PTOTY = PPROY + PTARY
            PTOTL = PPROL + PTARL
            ETOT  = EPRO  + ETAR
            PTOT  = SQRT(PTOTX**2+PTOTY**2+PTOTL**2)
            AMTOT = SQRT(ABS(ETOT-PTOT)*(ETOT+PTOT))
            GAM   = ETOT /(AMTOT)
            BGX   = PTOTX/(AMTOT)
            BGY   = PTOTY/(AMTOT)
            BGL   = PTOTL/(AMTOT)
            IF (IPEV.GE.6) THEN
               WRITE(LOUT,1000) PPROX,PPROY,PPROL,EPRO,AMPRO,KPROJ
               WRITE(LOUT,1001) PTARX,PTARY,PTARL,ETAR,AMTAR,KTARG
               WRITE(LOUT,1011) PTOTX,PTOTY,PTOTL,ETOT,AMTOT,KTARG
               WRITE(LOUT,1002) AMTOT,GAM,BGX,BGY,BGL
               WRITE(LOUT,1702) GAMV,BGLV
 1000          FORMAT('SDIFF: PPROX,PPROY,PPROL,EPRO,AMPRO,KPROJ',
     &                        2E15.5,2E15.5,E15.6,I2)
 1001          FORMAT('SDIFF: PTARX,PTARY,PTARL,ETAR,AMTAR,KTARG',
     &                        5E15.6,I2)
 1011          FORMAT('SDIFF: PTOTX,PTOTY,PTOTL,ETOT,AMTOT,KTARG',
     &                        5E15.6,I2)
 1002          FORMAT('SDIFF: AMTOT,GAM,BGX,BGY,BGL',5E15.6)
 1702          FORMAT('SDIFF: GAMV,BGLV',2E15.6)
            ENDIF
C                 This is the hadron-hadron cms mot the overall cms
C                     (due to Fermi momenta)
            CALL DALTRA(GAM,-BGX,-BGY,-BGL,PPROX,PPROY,PPROL,EPRO,
     &                  PPCM,PPCMX,PPCMY,PPCML,EPCM)
            CALL DALTRA(GAM,-BGX,-BGY,-BGL,PTARX,PTARY,PTARL,ETAR,
     &                  PTCM,PTCMX,PTCMY,PTCML,ETCM)
C                 Back to lab frame
	    CALL DALTRA(GAM,BGX,BGY,BGL,PPCMX,PPCMY,PPCML,EPCM,
     &                  PPLA,PPLAX,PPLAY,PPLAL,EPLA)
            CALL DALTRA(GAM,BGX,BGY,BGL,PTCMX,PTCMY,PTCML,ETCM,
     &                  PTLA,PTLAX,PTLAY,PTLAL,ETLA)
C                And now into the real cms frame     
             EPCMS=GAMV*EPLA-BGLV*PPLAL
	     PPLCMS=GAMV*PPLAL-BGLV*EPLA
	     ETCMS=GAMV*ETLA-BGLV*PTLAL
	     PTLCMS=GAMV*PTLAL-BGLV*ETLA
            IF (IPEV.GE.6) THEN
               WRITE(LOUT,1003) PPCM,PPCMX,PPCMY,PPCML,EPCM
               WRITE(LOUT,1004) PTCM,PTCMX,PTCMY,PTCML,ETCM
 1003          FORMAT('SDIFF: PPCM,PPCMX,PPCMY,PPCML,EPCM',5E15.5)
 1004          FORMAT('SDIFF: PTCM,PTCMX,PTCMY,PTCML,ETCM',5E15.5)
               WRITE(LOUT,1703) PPLA,PPLAX,PPLAY,PPLAL,EPLA
               WRITE(LOUT,1704) PTLA,PTLAX,PTLAY,PTLAL,ETLA
 1703          FORMAT('SDIFF: PPLA,PPLAX,PPLAY,PPLAL,EPLA',5E15.5)
 1704          FORMAT('SDIFF: PTLA,PTLAX,PTLAY,PTLAL,ETLA',5E15.5)
               WRITE(LOUT,1803) PPCM,PPCMX,PPCMY,PPLCMS,EPCMS
               WRITE(LOUT,1804) PTCM,PTCMX,PTCMY,PTLCMS,ETCMS
 1803          FORMAT('SDIFF: PPCM,PPCMX,PPCMY,PPLCMS,EPCMS',5E15.5)
 1804          FORMAT('SDIFF: PTCM,PTCMX,PTCMY,PTLCMS,ETCMS',5E15.5)
            ENDIF
            COD  = PPCML/PPCM
            COD2 = COD**2
            IF (COD2.GT.0.999999999999D0) COD2 = 0.999999999999D0
            SID  = SQRT((1.0D0-COD)*(1.0D0+COD))
            COF  = 1.0D0
            SIF  = 0.0D0
            IF (PPCM*SID.GT.1.0D-9) THEN
               COF   = PPCMX/(SID*PPCM)
               SIF   = PPCMY/(SID*PPCM)
               ANORF = SQRT(COF**2+SIF**2)
               COF   = COF/ANORF
               SIF   = SIF/ANORF
            ENDIF
            IF (IPEV.GE.6) THEN
               WRITE(LOUT,1005) COD,SID,COF,SIF
 1005          FORMAT('SDIFF: COD,SID,COF,SIF',4E15.5)
            ENDIF
C                 This is the hadron-hadron cms mot the overall cms
C                     (due to Fermi momenta)
            ECM = EPCM + ETCM
            IF (IPEV.GE.6) THEN
               WRITE(LOUT,1705)ECM,EPCM,ETCM 
 1705          FORMAT('SDIFF:ECM,EPCM,ETCM',3E15.5)
            ENDIF
C------------------------------------------------------------------
C
C                                         Test J.R.2/94
C
C------------------------------------------------------------------
            CALL SIHNDI(ECM,KPROJ,KTARG,SIGDIF,SIGDIH)
C------------------------------------------------------------------
C
C                                         Test J.R.2/94
C
C------------------------------------------------------------------
            FAKK=1.9
C ------------------------------------------------------------------
C
C                             further modification j.r.6.1.95
C
C--------------------------------------------------------------------
            IF(ECM.LE.10.D0)THEN
              FAKK=1.D0
            ELSEIF(ECM.GE.30.D0)THEN
              FAKK=1.9D0
            ELSE
              FAK=(ECM-10.D0)/20.D0
              FAKK=1.D0+FAK*0.9D0
            ENDIF
	    AITT=IT
            IF((ISINGD.LE.2).AND.(KPROJ.EQ.1.OR.KPROJ.EQ.8))THEN
	      AADIFF=FAKK*AITT**0.17D0
              SIGDIF=AADIFF*SIGDIF
            ELSEIF((ISINGD.LE.2).AND.
     *      ((KPROJ.EQ.13).OR.(KPROJ.EQ.14).OR.(KPROJ.EQ.23)))THEN
	      AADIFF=FAKK*AITT**0.15D0
              SIGDIF=AADIFF*SIGDIF
            ELSEIF((ISINGD.LE.2).AND.
     *      ((KPROJ.EQ.15).OR.(KPROJ.EQ.16).OR.
     *      (KPROJ.EQ.24).OR.(KPROJ.EQ.25)))THEN
	      AADIFF=FAKk*AITT**0.13D0
              SIGDIF=AADIFF*SIGDIF
	    ENDIF
C------------------------------------------------------------------
C
            SIGIN = DSHNTO(KPROJ,KTARG,ECM)-DSHNEL(KPROJ,KTARG,ECM)
            IF (IPEV.GE.6) THEN
        WRITE(LOUT,1060)KPROJ,KTARG,ECM,SIGDIF/SIGIN
               WRITE(LOUT,1006)SIGDIF,SIGDIF-SIGDIH,SIGDIH,SIGIN
 1006          FORMAT('SDIFF: SIGDIF,SIGDIL,SIGDIH,SIGIN',4F10.5)
 1060          FORMAT('SDIFF: KPROJ,KTARG,ECM,SIGDIF/SIGIN',2I3,2F10.5)
            ENDIF
            IFLAGD = 0
            R = RNDM(V)
            IF ((R.LE.(SIGDIF/SIGIN)).OR.(ISINGD.GE.2)) THEN
C              NCDIFF = NCDIFF+1
 2000          CONTINUE
C              IFLAGD = 1
C                         j.r.11/98 drop the following two lines
C              GAMV=GAM
C              BGLV=BGL
               R = RNDM(V)
*---------------------- S.Roesler 5/26/93
               IF (((R.LT.(SIGDIH/SIGDIF)).OR.(ISINGD.EQ.5).OR.
     &              (ISINGD.EQ.6)).AND.(ISINGD.LE.6)) THEN
*
*--------------------- call high mass single diffractive event
*
		  IF(ISINGD.GE.1)THEN
                  CALL VAHMSD(IHKK,ECM,KPROJ,KTARG,IREJ)
		  IFLAGD=1
		  ENDIF
*
               ELSE IF((R.GE.(SIGDIH/SIGDIF)).OR.(ISINGD.EQ.7).OR.
     &              (ISINGD.EQ.8)) THEN
*
*
*--------------------- call low mass single diffractive event
*
C----------------------------------------------------------------
C
C                                     j.r. test 2/94
C
C----------------------------------------------------------------
	    RRRR=RNDM(VV)
	    IF(ISINGD.GE.3.AND.ISINGD.LE.6)RRRR=1.D0
            IF((RRRR.LE.0.33D0).AND.
     *      ((KPROJ.EQ.15).OR.(KPROJ.EQ.16).OR.
     *      (KPROJ.EQ.1).OR.(KPROJ.EQ.8).OR.
     *      (KPROJ.EQ.13).OR.(KPROJ.EQ.14).OR.(KPROJ.EQ.23).OR.
     *      (KPROJ.EQ.24).OR.(KPROJ.EQ.25)))THEN
                  CALL VALMDD(IHKK,ECM,KPROJ,KTARG,IREJ)
		  IFLAGD=1
            ELSE      
C----------------------------------------------------------------
		  IF(ISINGD.GE.1)THEN
                  CALL VALMSD(IHKK,ECM,KPROJ,KTARG,IREJ)
		  IFLAGD=1
		  ENDIF
            ENDIF
*
               ENDIF
               IF (IREJ.GT.0) THEN
                  NCIREJ = NCIREJ + 1
                  IF (MOD(NCIREJ,1000).EQ.0) THEN
                     WRITE(LOUT,1007) NCIREJ
 1007                FORMAT('SDIFF: REJECTION, NCIREJ = ',I8)
                  ENDIF
                  GOTO 2000
               ENDIF
*
	       IF(IFLAGD.EQ.1)THEN
               CALL HADRDI(NAUX,KPROJ,KTARG,NHKKH1)
*
               IIHKK = NHKK -NAUX
               DO 11 J=1,NAUX
                  CALL DTRANS(PXC(J),PYC(J),PZC(J),
     &                        COD,SID,COF,SIF,PXX,PYY,PLL)
                  IIHKK = IIHKK + 1
                  IF (IPEV.GE.6) THEN
                     WRITE(LOUT,1008) IIHKK,PXX,PYY,PLL
 1008                FORMAT('SDIFF: NHKK,PXX,PYY,PLL',I4,3F10.5)
                  ENDIF
                  PHKK(1,IIHKK) = PXX
                  PHKK(2,IIHKK) = PYY
                  PHKK(3,IIHKK) = PLL
C                 IF (CMHIS.EQ.0.0D0) THEN
C                    CALL DALTRA(GAM,BGX,BGY,BGL,PXX,PYY,PLL,HEC(J),
C    &                    PLAB,PHKK(1,IIHKK),PHKK(2,IIHKK),
C    &                    PHKK(3,IIHKK),PHKK(4,IIHKK))
C                         Test j.r. 11/98
C---------------------------------------------------------------
C                 Back to lab frame
C           DO 1277 III=NHKKH1,NHKK
            III=IIHKK
	    IF(ISTHKK(III).EQ.1)THEN
	    CALL DALTRA(GAM,BGX,BGY,BGL,PHKK(1,III),PHKK(2,III),
     &      PHKK(3,III),PHKK(4,III),
     &                  PPLA,PPLAX,PPLAY,PPLAL,EPLA)
C                And now into the real cms frame     
             EPCMS=GAMV*EPLA-BGLV*PPLAL
	     PPLCMS=GAMV*PPLAL-BGLV*EPLA
	     PHKK(1,III)=PPLAX
	     PHKK(2,III)=PPLAY
	     PHKK(3,III)=PPLCMS
	     PHKK(4,III)=EPCMS
                  IF (IPEV.GE.6) THEN
               WRITE(LOUT,1903) PPLCMS,EPCMS
 1903          FORMAT('SDIFF TEST cms: PPLCMS,EPCMS',2E15.5)
                  ENDIF
	     ENDIF
C1277        CONTINUE	     
C---------------------------------------------------------------
                     IF (IPEV.GE.6) THEN
                        WRITE(LOUT,1009) IIHKK,
     &                  PHKK(1,IIHKK),PHKK(2,IIHKK),PHKK(3,IIHKK),
     &                  PHKK(4,IIHKK)
 1009                   FORMAT('SDIFF: NHKK,PHKK(1..4)',I4,4E15.5)
                     ENDIF
C                 ENDIF
   11          CONTINUE
	       ENDIF
            ENDIF
         ENDIF
   10 CONTINUE
      IF (KTARG.EQ.0) WRITE(LOUT,*)'SDIFF: NO INTERACTION'
 9999 CONTINUE
      RETURN
      END
*
*===vahmsd===============================================================*
*
      SUBROUTINE VAHMSD(ITAPOI,ECM,KPROJ,KTARG,IREJ)
 
**************************************************************************
* Version November 1993                         by Stefan Roesler        *
*                                                  University of Leipzig *
* This subroutine selects x-values, flavors and 4-momenta of partons     *
* in high-mass single diffractive chains.                                *
**************************************************************************
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)
      PARAMETER (NMXHKK= 89998)
      CHARACTER*8 ANAME
      COMMON /HKKEVT/ NHKK,NEVHKK,     ISTHKK(NMXHKK),  IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),PHKK(5,NMXHKK),
     &                VHKK(4,NMXHKK),  WHKK(4,NMXHKK)
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
      COMMON /DIFFRA/ ISINGD,IDIFTP,IOUDIF,IFLAGD
      COMMON /ABRDIF/ XDQ1,XDQ2,XDDQ1,XDDQ2,
     &                IKVQ1,IKVQ2,IKD1Q1,IKD2Q1,IKD1Q2,IKD2Q2,
     &                IDIFFP,IDIFAP,
     &                AMDCH1,AMDCH2,AMDCH3,GAMDC1,GAMDC2,GAMDC3,
     &                PGXVC1,PGYVC1,PGZVC1,PGXVC2,PGYVC2,PGZVC2,
     &                PGXVC3,PGYVC3,PGZVC3,NDCH1,NDCH2,NDCH3,
     &                IKDCH1,IKDCH2,IKDCH3,
     &                PDQ1(4),PDQ2(4),PDD1(4),PDD2(4),PDFQ1(4)
      COMMON /DPAR/   ANAME(210),AM(210),GA(210),TAU(210),ICH(210),
     &                IBAR(210),K1(210),K2(210)
      COMMON /TRAFOP/ GAMP,BGAMP,BETP
      COMMON /ENERIN/ EPROJ,ETARG
      COMMON /SDFLAG/ ISD
      COMMON /XDIDID/XDIDI
*
      DIMENSION MQUARK(3,30),IHKKQ(-6:6),IHKKQQ(-3:3,-3:3),
     &          IDX(-4:4)
      DATA IDX   /-4,-3,-1,-2,0,2,1,3,4/
      DATA IHKKQ /-6,-5,-4,-3,-1,-2,0,2,1,3,4,5,6/
      DATA IHKKQQ/-3301,-3103,-3203,0,   0,0,0,
     &            -3103,-1103,-2103,0,   0,   0,   0,
     &            -3203,-2103,-2203,0,   0,   0,   0,
     &                0,    0,    0,0,   0,   0,   0,
     &                0,    0,    0,0,2203,2103,3202,
     &                0,    0,    0,0,2103,1103,3103,
     &                0,    0,    0,0,3203,3103,3303/
*
*----------------------------------- quark content of hadrons:
*                                     1, 2, 3, 4  -   u, d, s, c
*                                    -1,-2,-3,-4  -  au,ad,as,ac
*
      DATA MQUARK/
     &   1,1,2,    -1,-1,-2,       0,0,0,       0,0,0,       0,0,0,
     &   0,0,0,       0,0,0,       1,2,2,    -1,-2,-2,       0,0,0,
     &   0,0,0,       0,0,0,      1,-2,0,      2,-1,0,      1,-3,0,
     &  3,-1,0,       1,2,3,    -1,-2,-3,       0,0,0,       2,2,3,
     &   1,1,3,       1,2,3,      1,-1,0,      2,-3,0,      3,-2,0,
     &  2,-2,0,      3,-3,0,       0,0,0,       0,0,0,       0,0,0/
      DATA UNON/2.0/
      DATA NCREJ, NCXDI, NCXP, NCXT /0, 0, 0, 0/
*
      ISD   = 1
      IREJ  = 0
      IIREJ = 0
      EPROJ = (AM(KPROJ)**2-AM(KTARG)**2+ECM**2)/(2.0D0*ECM)
      ETARG = (AM(KTARG)**2-AM(KPROJ)**2+ECM**2)/(2.0D0*ECM)
      IF(IPEV.GE.2) WRITE(LOUT,1001) EPROJ,ETARG
 1001 FORMAT('VAHMSD: EPROJ,ETARG ',2F10.5)
      IBPROJ = IBAR(KPROJ)
      IBTARG = IBAR(KTARG)
      IF (IBTARG.LE.0) THEN
         WRITE(LOUT,1002) IBTARG
 1002    FORMAT('VAHMSD: NO HMSD FOR TARGET WITH BARYON-CHARGE',I4)
         IIREJ=1
         GOTO 9999
      ENDIF
      IQP1 = MQUARK(1,KPROJ)
      IQP2 = MQUARK(2,KPROJ)
      IQP3 = MQUARK(3,KPROJ)
      IQT1 = MQUARK(1,KTARG)
      IQT2 = MQUARK(2,KTARG)
      IQT3 = MQUARK(3,KTARG)
      IF(IPEV.GE.2) WRITE(LOUT,1003)
     &                IBPROJ,IBTARG,IQP1,IQP2,IQP3,IQT1,IQT2,IQT3
 1003 FORMAT('VAHMSD: IBPROJ,IBTARG,IQP1,IQP2,IQP3,IQT1,IQT2,IQT3 ',8I4)
*
      IF (IBPROJ.NE.0) THEN
*
*-------------------- q-qq (aq-aqaq) - flavors of projectile (baryon)
*
         ISAM = 1.0D0+2.999D0*RNDM(V)
         GOTO (10,11,12) ISAM
   10    CONTINUE
         IQP    = IQP1
         IDIQP1 = IQP2
         IDIQP2 = IQP3
         GOTO 13
   11    CONTINUE
         IQP    = IQP2
         IDIQP1 = IQP1
         IDIQP2 = IQP3
         GOTO 13
   12    CONTINUE
         IQP    = IQP3
         IDIQP1 = IQP1
         IDIQP2 = IQP2
   13    CONTINUE
*
      ELSE IF (IBPROJ.EQ.0) THEN
*
*-------------------- q-aq - flavors of projectile (meson)
*
         ISAM = 1.0D0+1.999D0*RNDM(V)
         GOTO (14,15) ISAM
   14    CONTINUE
         IQP    = IQP1
         IDIQP1 = IQP2
         GOTO 16
   15    CONTINUE
         IQP    = IQP2
         IDIQP1 = IQP1
   16    CONTINUE
         IDIQP2 = 0
*
      ENDIF
*
*-------------------- q-qq - flavors of target (baryon)
*
      ISAM = 1.0D0+2.999D0*RNDM(V)
      GOTO (17,18,19) ISAM
   17 CONTINUE
      IQT    = IQT1
      IDIQT1 = IQT2
      IDIQT2 = IQT3
      GOTO 20
   18 CONTINUE
      IQT    = IQT2
      IDIQT1 = IQT1
      IDIQT2 = IQT3
      GOTO 20
   19 CONTINUE
      IQT    = IQT3
      IDIQT1 = IQT1
      IDIQT2 = IQT2
   20 CONTINUE
*
      IKVQ1  = IQP
      IKD1Q1 = IDIQP1
      IKD2Q1 = IDIQP2
*
      IKVQ2  = IQT
      IKD1Q2 = IDIQT1
      IKD2Q2 = IDIQT2
*
      IF (IPEV.GE.2) WRITE(LOUT,1004)
     &                 IKVQ1,IKD1Q1,IKD2Q1,IKVQ2,IKD1Q2,IKD2Q2
 1004 FORMAT('VAHMSD: IKVQ1,IKD1Q1,IKD2Q1,IKVQ2,IKD1Q2,IKD2Q2 ',6I4)
*
*-------------------- q-aq - flavors of diffractive parton
*
      IDIFFP = 1.0D0+2.3D0*RNDM(V)
      IDIFAP = -IDIFFP
      IF (IPEV.GE.2) WRITE(LOUT,1005) IDIFFP,IDIFAP
 1005 FORMAT('VAHMSD: IDIFFP,IDIFAP ',2I4)
*
*-------------------- IDIFTP = 1  target     (backward) hadron  excited
*                     IDIFTP = 2  projectile (forward)  hadron  excited
*
      IDIFTP = 1.0D0+1.999D0*RNDM(V)
*----------------------- S.Roesler 5/26/93
      IF ((ISINGD.EQ.3).OR.(ISINGD.EQ.5)) IDIFTP = 1
      IF ((ISINGD.EQ.4).OR.(ISINGD.EQ.6)) IDIFTP = 2
*
      IF (IPEV.GE.2) WRITE(LOUT,1006) IDIFTP
 1006 FORMAT('VAHMSD: IDIFTP ',I4)
      IF ((IDIFTP.NE.1).AND.(IDIFTP.NE.2)) THEN
         IF (IPEV.GE.2) WRITE(LOUT,'(A19)') 'VAHMSD-ERROR: IDIFTP'
         GOTO 9999
      ENDIF
*
*-------------------- momentum fractions of quarks and diquarks
*
      XXMAX = 0.8D0
   30 CONTINUE
      XP  = DBETAR(0.5D0,UNON)
      IF (XP.GE.XXMAX) THEN
         NCXP = NCXP+1
         IF(MOD(NCXP,500).EQ.0) WRITE(LOUT,1007) NCXP
 1007    FORMAT('VAHMSD: INEFFICIENT XP-SELECTION, NCXP=',I8)
         GOTO 30
      ENDIF
      XXP = 1.0D0-XP
   31 CONTINUE
      XT  = DBETAR(0.5D0,UNON)
      IF (XT.GE.XXMAX)  THEN
         NCXT = NCXT+1
         IF(MOD(NCXT,2500).EQ.0) WRITE(LOUT,1008) NCXT
 1008    FORMAT('VAHMSD: INEFFICIENT XT-SELECTION, NCXT=',I8)
         GOTO 31
      ENDIF
      XXT = 1.0D0-XT
      IF (IPEV.GE.2) WRITE(LOUT,1010) XP,XXP,XT,XXT
 1010 FORMAT('VAHMSD: XP,XXP,XT,XXT ',4D10.5)
*-------------------- x - values of diffractive partons
*
      NCXDI = 0
   32 CONTINUE
      R = RNDM(V)
      IF (IDIFTP.EQ.1) THEN
         XDIMIN = (3.0D0+400.0D0*(R**2))/(4.0D0*(ETARG**2)*XXT)
         IF (ECM.LE.300.0D0) THEN
            RR = (1.0D0-EXP(-((ECM/140.0D0)**4)))
            XDIMIN = (3.0D0+400.0D0*(R**2)*RR)/(4.0D0*(ETARG**2)*XXT)
         ENDIF
      ELSE IF (IDIFTP.EQ.2) THEN
         XDIMIN = (3.0D0+400.0D0*(R**2))/(4.0D0*(EPROJ**2)*XXP)
         IF (ECM.LE.300.0D0) THEN
            RR = (1.0D0-EXP(-((ECM/140.0D0)**4)))
            XDIMIN = (3.0D0+400.0D0*(R**2)*RR)/(4.0D0*(EPROJ**2)*XXP)
         ENDIF
      ENDIF
C-----------------------------------------------------------------
C                                    original version
C-----------------------------------------------------------------
C     XDIMAX = 0.05D0
C     IF (ECM.LE.1000.0D0) THEN
C        XDIMAX = 0.05D0*(1.0D0+EXP(-((ECM/420.0D0)**2)))
C        IF (IBPROJ.EQ.0) XDIMAX = 0.05D0*
CC---------------- change mass-cuts
CC   &                    (1.0D0+4.0D0*EXP(-((ECM/420.0D0)**2)))
C    &                    (1.0D0+2.0D0*EXP(-((ECM/420.0D0)**2)))
C     ENDIF
C-----------------------------------------------------------------
C-----------------------------------------------------------------
C                           version j.r.28.1.94
C extend diffraction beyond limits of single diffractive experiments
C-----------------------------------------------------------------
      XDIMAA = 0.15D0
      XDIMAX = XDIMAA
      IF (ECM.LE.10000.0D0) THEN
         XDIMAX = XDIMAA*(1.0D0+EXP(-((ECM/420.0D0)**2)))
         IF (IBPROJ.EQ.0) XDIMAX = XDIMAA*
C----------------- change mass-cuts
C    &                    (1.0D0+4.0D0*EXP(-((ECM/420.0D0)**2)))
     &                    (1.0D0+2.0D0*EXP(-((ECM/420.0D0)**2)))
      ENDIF
C-----------------------------------------------------------------
   40 CONTINUE
      IF (XDIMIN.GE.XDIMAX) THEN
         NCXDI = NCXDI+1
         IF (NCXDI.EQ.200) GOTO 9999
         GOTO 32
      ENDIF
      XDITOT = SAMPEY(XDIMIN,XDIMAX)
      R    = RNDM(V)**6
      DX   = XDITOT-XDIMIN
      XDI  = DX*R
C     IF (XDI.LT.XMINQ1) THEN
C        NCXDI = NCXDI+1
C        IF (MOD(NCXDI,2000).EQ.0) WRITE(LOUT,1011) NCXDI
C1011    FORMAT('VAHMSD: INEFFICIENT XDI-SELECTION, NCXDI=',I8)
C        GOTO 40
C     ENDIF
      XXDI = XDIMIN+(1.0D0-R)*DX
      IF (IDIFTP.EQ.1) THEN
         XXP = XXP-XDI-XXDI
         IF (XXP.LT.8.0D-2) GOTO 9999
      ELSE IF (IDIFTP.EQ.2) THEN
         XXT = XXT-XDI-XXDI
         IF (XXT.LT.8.0D-2) GOTO 9999
      ENDIF
      IF (IPEV.GE.2) WRITE(LOUT,1012) XP,XXP,XT,XXT,XDI,XXDI
 1012 FORMAT('VAHMSD: XP,XXP,XT,XXT,XDI,XXDI ',6F10.5)
       XDIDI=XDI+XXDI
       AMDIDI=SQRT(XDIDI*ECM**2)
      IF(IPEV.GE.2)WRITE(LOUT,*)'HM AMDIDI,XDIDI ',AMDIDI,XDIDI
*
*-------------------- kinematical parameters of three chains in CMS
*
      IF ((IBPROJ.EQ.-1).AND.(IDIFTP.EQ.1)) THEN
*
*-------------------- target: baryon,   projectile: antibaryon
*                             (excited)
*
         XDIFAP = XDI
         XDIFFP = XXDI
         CALL DIFFCH (XDIFAP, IDIFAP,     XT,  IKVQ2,     99,
     &                XDIFFP, IDIFFP,    XXT,    XXP,  ETARG,
     &                AMDCH1,   ECH1,   PCH1, GAMDC1,  PGVC1,
     &                 NDCH1, IKDCH1,  EPROJ,   NUNO,  IIREJ,      1)
         IF (IIREJ.EQ.1) GOTO 9999
         CALL DIFFCH (XDIFFP, IDIFFP,    XXT, IKD1Q2, IKD2Q2,
     &                   DUM,   IDUM,    DUM,    XXP,  ETARG,
     &                AMDCH2,   ECH2,   PCH2, GAMDC2,  PGVC2,
     &                 NDCH2, IKDCH2,  EPROJ,   NUNO,  IIREJ,      2)
         IF (IIREJ.EQ.1) GOTO 9999
         CALL DIFFCH (    XP,  IKVQ1,    XXP, IKD1Q1, IKD2Q1,
     &                   DUM,   IDUM,    DUM,    DUM,  EPROJ,
     &                AMDCH3,   ECH3,   PCH3, GAMDC3,  PGVC3,
     &                 NDCH3,   IDUM,    DUM,   NUNO,  IIREJ,      3)
         IF (IIREJ.EQ.1)  GOTO 9999
*
*
      ELSE IF ((IBPROJ.EQ.-1).AND.(IDIFTP.EQ.2)) THEN
*
*-------------------- target: baryon,   projectile: antibaryon
*                                                   (excited)
*
         XDIFAP = XXDI
         XDIFFP = XDI
         CALL DIFFCH (XDIFFP, IDIFFP,     XP,  IKVQ1,     99,
     &                XDIFAP, IDIFAP,    XXP,    XXT,  EPROJ,
     &                AMDCH1,   ECH1,   PCH1, GAMDC1,  PGVC1,
     &                 NDCH1, IKDCH1,  ETARG,   NUNO,  IIREJ,      1)
         IF (IIREJ.EQ.1) GOTO 9999
         CALL DIFFCH (XDIFAP, IDIFAP,    XXP, IKD1Q1, IKD2Q1,
     &                   DUM,   IDUM,    DUM,    XXT,  EPROJ,
     &                AMDCH2,   ECH2,   PCH2, GAMDC2,  PGVC2,
     &                 NDCH2, IKDCH2,  ETARG,   NUNO,  IIREJ,      2)
         IF (IIREJ.EQ.1) GOTO 9999
         CALL DIFFCH (    XT,  IKVQ2,    XXT, IKD1Q2, IKD2Q2,
     &                   DUM,   IDUM,    DUM,    DUM,  ETARG,
     &                AMDCH3,   ECH3,   PCH3, GAMDC3,  PGVC3,
     &                 NDCH3,   IDUM,    DUM,   NUNO,  IIREJ,      3)
         IF (IIREJ.EQ.1)  GOTO 9999
*
*
      ELSE IF ((IBPROJ.EQ.0).AND.(IDIFTP.EQ.1)) THEN
*
*-------------------- target: baryon,   projectile: meson
*                             (excited)
*
         XDIFAP = XDI
         XDIFFP = XXDI
         CALL DIFFCH (XDIFAP, IDIFAP,     XT,  IKVQ2,     99,
     &                XDIFFP, IDIFFP,    XXT,    XXP,  ETARG,
     &                AMDCH1,   ECH1,   PCH1, GAMDC1,  PGVC1,
     &                 NDCH1, IKDCH1,  EPROJ,   NUNO,  IIREJ,      1)
         IF (IIREJ.EQ.1) GOTO 9999
         CALL DIFFCH (XDIFFP, IDIFFP,    XXT, IKD1Q2, IKD2Q2,
     &                   DUM,   IDUM,    DUM,    XXP,  ETARG,
     &                AMDCH2,   ECH2,   PCH2, GAMDC2,  PGVC2,
     &                 NDCH2, IKDCH2,  EPROJ,   NUNO,  IIREJ,      2)
         IF (IIREJ.EQ.1) GOTO 9999
         CALL DIFFCH (    XP,  IKVQ1,    XXP, IKD1Q1,     99,
     &                   DUM,   IDUM,    DUM,    DUM,  EPROJ,
     &                AMDCH3,   ECH3,   PCH3, GAMDC3,  PGVC3,
     &                 NDCH3,   IDUM,    DUM,   NUNO,  IIREJ,      3)
         IF (IIREJ.EQ.1)  GOTO 9999
*
*
      ELSEIF ((IBPROJ.EQ.0).AND.(IDIFTP.EQ.2)) THEN
*
*-------------------- target: baryon,   projectile: meson
*                                                   (excited)
*
         IF (IKD1Q1.LT.0) THEN
            XDIFAP = XDI
            XDIFFP = XXDI
            CALL DIFFCH (XDIFAP, IDIFAP,     XP,  IKVQ1,     99,
     &                   XDIFFP, IDIFFP,    XXP,    XXT,  EPROJ,
     &                   AMDCH1,   ECH1,   PCH1, GAMDC1,  PGVC1,
     &                    NDCH1, IKDCH1,  ETARG,   NUNO,  IIREJ,      1)
            IF (IIREJ.EQ.1)  GOTO 9999
            CALL DIFFCH (XDIFFP, IDIFFP,    XXP, IKD1Q1,     99,
     &                   XDIFAP, IDIFAP,     XP,    XXT,  EPROJ,
     &                   AMDCH2,   ECH2,   PCH2, GAMDC2,  PGVC2,
     &                    NDCH2, IKDCH2,  ETARG,   NUNO,  IIREJ,      1)
            IF (IIREJ.EQ.1)  GOTO 9999
         ELSE
            XDIFAP = XXDI
            XDIFFP = XDI
            CALL DIFFCH (XDIFFP, IDIFFP,     XP,  IKVQ1,     99,
     &                   XDIFAP, IDIFAP,    XXP,    XXT,  EPROJ,
     &                   AMDCH1,   ECH1,   PCH1, GAMDC1,  PGVC1,
     &                    NDCH1, IKDCH1,  ETARG,   NUNO,  IIREJ,      1)
            IF (IIREJ.EQ.1)  GOTO 9999
            CALL DIFFCH (XDIFAP, IDIFAP,    XXP, IKD1Q1,     99,
     &                   XDIFFP, IDIFFP,     XP,    XXT,  EPROJ,
     &                   AMDCH2,   ECH2,   PCH2, GAMDC2,  PGVC2,
     &                    NDCH2, IKDCH2,  ETARG,   NUNO,  IIREJ,      1)
            IF (IIREJ.EQ.1)  GOTO 9999
         ENDIF
         CALL DIFFCH (    XT,  IKVQ2,    XXT, IKD1Q2, IKD2Q2,
     &                   DUM,   IDUM,    DUM,    DUM,  ETARG,
     &                AMDCH3,   ECH3,   PCH3, GAMDC3,  PGVC3,
     &                 NDCH3,   IDUM,    DUM,   NUNO,  IIREJ,      3)
         IF (IIREJ.EQ.1)  GOTO 9999
*
*
      ELSEIF ((IBPROJ.EQ.1).AND.(IDIFTP.EQ.1)) THEN
*
*-------------------- target: baryon,   projectile: baryon
*                             (excited)
*
         XDIFAP = XDI
         XDIFFP = XXDI
         CALL DIFFCH (XDIFAP, IDIFAP,     XT,  IKVQ2,     99,
     &                XDIFFP, IDIFFP,    XXT,    XXP,  ETARG,
     &                AMDCH1,   ECH1,   PCH1, GAMDC1,  PGVC1,
     &                 NDCH1, IKDCH1,  EPROJ,   NUNO,  IIREJ,      1)
         IF (IIREJ.EQ.1)  GOTO 9999
         CALL DIFFCH (XDIFFP, IDIFFP,    XXT, IKD1Q2, IKD2Q2,
     &                   DUM,   IDUM,    DUM,    XXP,  ETARG,
     &                AMDCH2,   ECH2,   PCH2, GAMDC2,  PGVC2,
     &                 NDCH2, IKDCH2,  EPROJ,   NUNO,  IIREJ,      2)
         IF (IIREJ.EQ.1)  GOTO 9999
         CALL DIFFCH (    XP,  IKVQ1,    XXP, IKD1Q1, IKD2Q1,
     &                   DUM,   IDUM,    DUM,    DUM,  EPROJ,
     &                AMDCH3,   ECH3,   PCH3, GAMDC3,  PGVC3,
     &                 NDCH3,   IDUM,    DUM,   NUNO,  IIREJ,      3)
         IF (IIREJ.EQ.1)  GOTO 9999
*
*
      ELSE IF ((IBPROJ.EQ.1).AND.(IDIFTP.EQ.2)) THEN
*
*-------------------- target: baryon,   projectile: baryon
*                                                   (excited)
*
         XDIFAP = XDI
         XDIFFP = XXDI
         CALL DIFFCH (XDIFAP, IDIFAP,     XP,  IKVQ1,     99,
     &                XDIFFP, IDIFFP,    XXP,    XXT,  EPROJ,
     &                AMDCH1,   ECH1,   PCH1, GAMDC1,  PGVC1,
     &                 NDCH1, IKDCH1,  ETARG,   NUNO,  IIREJ,      1)
         IF (IIREJ.EQ.1)  GOTO 9999
         CALL DIFFCH (XDIFFP, IDIFFP,    XXP, IKD1Q1, IKD2Q1,
     &                   DUM,   IDUM,    DUM,    XXT,  EPROJ,
     &                AMDCH2,   ECH2,   PCH2, GAMDC2,  PGVC2,
     &                 NDCH2, IKDCH2,  ETARG,   NUNO,  IIREJ,      2)
         IF (IIREJ.EQ.1) GOTO 9999
         CALL DIFFCH (    XT,  IKVQ2,    XXT, IKD1Q2, IKD2Q2,
     &                   DUM,   IDUM,    DUM,    DUM,  ETARG,
     &                AMDCH3,   ECH3,   PCH3, GAMDC3,  PGVC3,
     &                 NDCH3,   IDUM,    DUM,   NUNO,  IIREJ,      3)
         IF (IIREJ.EQ.1)  GOTO 9999
*
      ENDIF
*
*-------------------- store results in common block /ABRDIF/
*
      XDQ1  = XP
      XDQ2  = XT
      XDDQ1 = XXP
      XDDQ2 = XXT
*
      IF (IPEV.GE.2) THEN
         WRITE(LOUT,1013) AMDCH1,ECH1,PCH1,GAMDC1,PGVC1,IKDCH1,NDCH1
 1013    FORMAT('VAHMSD: AMDCH1,ECH1,PCH1,GAMDC1,PGVC1,IKDCH1,NDCH1 ',
     &                                                     5F10.5,2I4)
         WRITE(LOUT,1014) AMDCH2,ECH2,PCH2,GAMDC2,PGVC2,IKDCH2,NDCH2
 1014    FORMAT('VAHMSD: AMDCH2,ECH2,PCH2,GAMDC2,PGVC2,IKDCH2,NDCH2 ',
     &                                                     5F10.5,2I4)
         WRITE(LOUT,1015) AMDCH3,ECH3,PCH3,GAMDC3,PGVC3,NDCH3
 1015    FORMAT('VAHMSD: AMDCH3,ECH3,PCH3,GAMDC3,PGVC3,NDCH3 ',
     &                                                     5F10.5,I4)
      ENDIF
*
*-------------------- select transverse momenta
*
      IF (IDIFTP.EQ.1) THEN
         CALL DIFFPT(   ECH1,   PCH1, XDIFAP,     XT,
     &                  ECH2,   PCH2, XDIFFP,    XXT,
     &                  ECH3,   PCH3,  EPROJ,  KPROJ,
     &                 ETARG,  EPROJ,  IIREJ)
      ELSE IF (IDIFTP.EQ.2) THEN
         IF ((IBPROJ.LT.0).OR.(IKVQ1.LT.0)) THEN
            CALL DIFFPT(   ECH1,   PCH1, XDIFFP,     XP,
     &                     ECH2,   PCH2, XDIFAP,    XXP,
     &                     ECH3,   PCH3,  ETARG,  KTARG,
     &                    EPROJ,  ETARG,  IIREJ)
         ELSE
            CALL DIFFPT(   ECH1,   PCH1, XDIFAP,     XP,
     &                     ECH2,   PCH2, XDIFFP,    XXP,
     &                     ECH3,   PCH3,  ETARG,  KTARG,
     &                    EPROJ,  ETARG,  IIREJ)
         ENDIF
      ENDIF
      IF (IIREJ.EQ.1) GOTO 9999
*
*-------------------- store results in common-block /HKKEVT/
*
*                     partonen of projectile
*
      NHKK = NHKK+1
      ISTHKK(NHKK)   = 21
      IDHKK (NHKK)   = IHKKQ(IDX(IKVQ1))
      JMOHKK(1,NHKK) = 1
      JMOHKK(2,NHKK) = 0
      JDAHKK(1,NHKK) = 0
      JDAHKK(2,NHKK) = 0
      PHKK  (1,NHKK) = 0.0D0
      PHKK  (2,NHKK) = 0.0D0
      PHKK  (3,NHKK) = XDQ1
      PHKK  (4,NHKK) = XDQ1
      PHKK  (5,NHKK) = 0.0D0
      VHKK  (1,NHKK) = VHKK(1,1)
      VHKK  (2,NHKK) = VHKK(2,1)
      VHKK  (3,NHKK) = VHKK(3,1)
      VHKK  (4,NHKK) = VHKK(4,1)
*
      NHKK = NHKK+1
      ISTHKK(NHKK)   = 21
      IDHKK (NHKK)   = IHKKQQ(IDX(IKD1Q1),IDX(IKD2Q1))
      JMOHKK(1,NHKK) = 1
      JMOHKK(2,NHKK) = 0
      JDAHKK(1,NHKK) = 0
      JDAHKK(2,NHKK) = 0
      PHKK  (1,NHKK) = 0.0D0
      PHKK  (2,NHKK) = 0.0D0
      PHKK  (3,NHKK) = XDDQ1
      PHKK  (4,NHKK) = XDDQ1
      PHKK  (5,NHKK) = 0.0D0
      VHKK  (1,NHKK) = VHKK(1,1)
      VHKK  (2,NHKK) = VHKK(2,1)
      VHKK  (3,NHKK) = VHKK(3,1)
      VHKK  (4,NHKK) = VHKK(4,1)
*
*                     partonen of target
*
      NHKK = NHKK+1
      ISTHKK(NHKK)   = 22
      IDHKK (NHKK)   = IHKKQ(IDX(IKVQ2))
      JMOHKK(1,NHKK) = ITAPOI
      JMOHKK(2,NHKK) = 0
      JDAHKK(1,NHKK) = 0
      JDAHKK(2,NHKK) = 0
      PHKK  (1,NHKK) = 0.0D0
      PHKK  (2,NHKK) = 0.0D0
      PHKK  (3,NHKK) = XDQ2
      PHKK  (4,NHKK) = XDQ2
      PHKK  (5,NHKK) = 0.0D0
      VHKK  (1,NHKK) = VHKK(1,ITAPOI)
      VHKK  (2,NHKK) = VHKK(2,ITAPOI)
      VHKK  (3,NHKK) = VHKK(3,ITAPOI)
      VHKK  (4,NHKK) = VHKK(4,ITAPOI)
*
      NHKK = NHKK+1
      ISTHKK(NHKK)   = 22
      IDHKK (NHKK)   = IHKKQQ(IDX(IKD1Q2),IDX(IKD2Q2))
      JMOHKK(1,NHKK) = ITAPOI
      JMOHKK(2,NHKK) = 0
      JDAHKK(1,NHKK) = 0
      JDAHKK(2,NHKK) = 0
      PHKK  (1,NHKK) = 0.0D0
      PHKK  (2,NHKK) = 0.0D0
      PHKK  (3,NHKK) = XDDQ2
      PHKK  (4,NHKK) = XDDQ2
      PHKK  (5,NHKK) = 0.0D0
      VHKK  (1,NHKK) = VHKK(1,ITAPOI)
      VHKK  (2,NHKK) = VHKK(2,ITAPOI)
      VHKK  (3,NHKK) = VHKK(3,ITAPOI)
      VHKK  (4,NHKK) = VHKK(4,ITAPOI)
*
*                     sea q-aq pair
*
      NHKK = NHKK+1
      ISTHKK(NHKK)   = 30+IDIFTP
      IDHKK (NHKK)   = IHKKQ(IDX(IDIFFP))
      IF (IDIFTP.EQ.1) JMOHKK(1,NHKK) = 1
      IF (IDIFTP.EQ.2) JMOHKK(1,NHKK) = ITAPOI
      JMOHKK(2,NHKK) = 0
      JDAHKK(1,NHKK) = 0
      JDAHKK(2,NHKK) = 0
      PHKK  (1,NHKK) = 0.0D0
      PHKK  (2,NHKK) = 0.0D0
      PHKK  (3,NHKK) = XXDI
      PHKK  (4,NHKK) = XXDI
      PHKK  (5,NHKK) = 0.0D0
      VHKK  (1,NHKK) = VHKK(1,JMOHKK(1,NHKK))
      VHKK  (2,NHKK) = VHKK(2,JMOHKK(1,NHKK))
      VHKK  (3,NHKK) = VHKK(3,JMOHKK(1,NHKK))
      VHKK  (4,NHKK) = VHKK(4,JMOHKK(1,NHKK))
*
      NHKK = NHKK+1
      ISTHKK(NHKK)   = 30+IDIFTP
      IDHKK (NHKK)   = IHKKQ(IDX(IDIFAP))
      IF (IDIFTP.EQ.1) JMOHKK(1,NHKK) = 1
      IF (IDIFTP.EQ.2) JMOHKK(1,NHKK) = ITAPOI
      JMOHKK(2,NHKK) = 0
      JDAHKK(1,NHKK) = 0
      JDAHKK(2,NHKK) = 0
      PHKK  (1,NHKK) = 0.0D0
      PHKK  (2,NHKK) = 0.0D0
      PHKK  (3,NHKK) = XDI
      PHKK  (4,NHKK) = XDI
      PHKK  (5,NHKK) = 0.0D0
      VHKK  (1,NHKK) = VHKK(1,JMOHKK(1,NHKK))
      VHKK  (2,NHKK) = VHKK(2,JMOHKK(1,NHKK))
      VHKK  (3,NHKK) = VHKK(3,JMOHKK(1,NHKK))
      VHKK  (4,NHKK) = VHKK(4,JMOHKK(1,NHKK))
*
*                     ends of chain 1
*
      NHKK = NHKK+1
      ISTHKK(NHKK)   = 123-IDIFTP
      IND = NHKK-2-2*IDIFTP
      IDHKK (NHKK)   = IDHKK(IND)
      JMOHKK(1,NHKK) = IND
      JMOHKK(2,NHKK) = JMOHKK(1,IND)
      JDAHKK(1,NHKK) = NHKK+2
      JDAHKK(2,NHKK) = NHKK+2
      PHKK  (1,NHKK) = PDQ1(1)
      PHKK  (2,NHKK) = PDQ1(2)
      PHKK  (3,NHKK) = PDQ1(3)
      PHKK  (4,NHKK) = PDQ1(4)
      PHKK  (5,NHKK) = 0.0D0
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
      VHKK  (1,NHKK) = VHKK(1,IND)+XXPP
      VHKK  (2,NHKK) = VHKK(2,IND)+YYPP
      VHKK  (3,NHKK) = VHKK(3,IND)
      VHKK  (4,NHKK) = VHKK(4,IND)
*
      NHKK = NHKK+1
      ISTHKK(NHKK)   = 130+IDIFTP
      IF (IDHKK(NHKK-1).GT.0) JMOHKK(1,NHKK) = NHKK-2
      IF (IDHKK(NHKK-1).LE.0) JMOHKK(1,NHKK) = NHKK-3
      JMOHKK(2,NHKK) = JMOHKK(1,NHKK-3)
      JDAHKK(1,NHKK) = NHKK+1
      JDAHKK(2,NHKK) = NHKK+1
      IDHKK (NHKK)   = IDHKK(JMOHKK(1,NHKK))
      PHKK  (1,NHKK) = PDD2(1)
      PHKK  (2,NHKK) = PDD2(2)
      PHKK  (3,NHKK) = PDD2(3)
      PHKK  (4,NHKK) = PDD2(4)
      PHKK  (5,NHKK) = 0.0D0
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
      VHKK  (1,NHKK) = VHKK(1,JMOHKK(1,NHKK))+XXPP
      VHKK  (2,NHKK) = VHKK(2,JMOHKK(1,NHKK))+YYPP
      VHKK  (3,NHKK) = VHKK(3,JMOHKK(1,NHKK))
      VHKK  (4,NHKK) = VHKK(4,JMOHKK(1,NHKK))
*
*                     chain 1
*
      NHKK = NHKK+1
      ISTHKK(NHKK)   = 199
      IDHKK (NHKK)   = 88888
      JMOHKK(1,NHKK) = NHKK-2
      JMOHKK(2,NHKK) = NHKK-1
      JDAHKK(1,NHKK) = 0
      JDAHKK(2,NHKK) = 0
      PHKK  (5,NHKK) = AMDCH1
      VHKK  (1,NHKK) = VHKK(1,NHKK-1)
      VHKK  (2,NHKK) = VHKK(2,NHKK-1)
      VHKK  (3,NHKK) = VHKK(3,NHKK-1)
      IF ((BETP.NE.0.0D0).AND.(BGAMP.NE.0.0D0))
     &VHKK  (4,NHKK) = VHKK(3,NHKK)/BETP-VHKK(3,NHKK-2)/BGAMP
*
*                     ends of chain 2
*
      NHKK = NHKK+1
      ISTHKK(NHKK)   = 123-IDIFTP
      IND = NHKK-4-2*IDIFTP
      IDHKK (NHKK)   = IDHKK(IND)
      JMOHKK(1,NHKK) = IND
      JMOHKK(2,NHKK) = JMOHKK(1,IND)
      JDAHKK(1,NHKK) = NHKK+2
      JDAHKK(2,NHKK) = NHKK+2
      PHKK  (1,NHKK) = PDD1(1)
      PHKK  (2,NHKK) = PDD1(2)
      PHKK  (3,NHKK) = PDD1(3)
      PHKK  (4,NHKK) = PDD1(4)
      PHKK  (5,NHKK) = 0.0D0
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
      VHKK  (1,NHKK) = VHKK(1,IND)+XXPP
      VHKK  (2,NHKK) = VHKK(2,IND)+YYPP
      VHKK  (3,NHKK) = VHKK(3,IND)
      VHKK  (4,NHKK) = VHKK(4,IND)
*
      NHKK = NHKK+1
      ISTHKK(NHKK)   = 130+IDIFTP
      JMOHKK(1,NHKK) = 2*NHKK-11-JMOHKK(1,NHKK-3)
      JMOHKK(2,NHKK) = JMOHKK(1,NHKK-6)
      JDAHKK(1,NHKK) = NHKK+1
      JDAHKK(2,NHKK) = NHKK+1
      IDHKK (NHKK)   = IDHKK(JMOHKK(1,NHKK))
      PHKK  (1,NHKK) = PDQ2(1)
      PHKK  (2,NHKK) = PDQ2(2)
      PHKK  (3,NHKK) = PDQ2(3)
      PHKK  (4,NHKK) = PDQ2(4)
      PHKK  (5,NHKK) = 0.0D0
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
      VHKK  (1,NHKK) = VHKK(1,JMOHKK(1,NHKK-6))+XXPP
      VHKK  (2,NHKK) = VHKK(2,JMOHKK(1,NHKK-6))+YYPP
      VHKK  (3,NHKK) = VHKK(3,JMOHKK(1,NHKK-6))
      VHKK  (4,NHKK) = VHKK(4,JMOHKK(1,NHKK-6))
*
*                     chain 2
*
      NHKK = NHKK+1
      ISTHKK(NHKK)   = 199
      IDHKK (NHKK)   = 88888
      JMOHKK(1,NHKK) = NHKK-2
      JMOHKK(2,NHKK) = NHKK-1
      JDAHKK(1,NHKK) = 0
      JDAHKK(2,NHKK) = 0
      PHKK  (5,NHKK) = AMDCH2
      VHKK  (1,NHKK) = VHKK(1,NHKK-1)
      VHKK  (2,NHKK) = VHKK(2,NHKK-1)
      VHKK  (3,NHKK) = VHKK(3,NHKK-1)
      IF ((BETP.NE.0.0D0).AND.(BGAMP.NE.0.0D0))
     &VHKK  (4,NHKK) = VHKK(3,NHKK)/BETP-VHKK(3,NHKK-2)/BGAMP
*
*                     diffractive nucleon
*
      NHKK = NHKK+1
      ISTHKK(NHKK)   = 199
      IDHKK (NHKK)   = 88888
      JMOHKK(1,NHKK) = NHKK-14+2*IDIFTP
      JMOHKK(2,NHKK) = NHKK-13+2*IDIFTP
      JDAHKK(1,NHKK) = 0
      JDAHKK(2,NHKK) = 0
      PHKK  (1,NHKK) = PDFQ1(1)
      PHKK  (2,NHKK) = PDFQ1(2)
      PHKK  (3,NHKK) = PDFQ1(3)
      PHKK  (4,NHKK) = PDFQ1(4)
      PHKK  (5,NHKK) = AMDCH3
      VHKK  (1,NHKK) = VHKK(1,JMOHKK(1,NHKK))
      VHKK  (2,NHKK) = VHKK(2,JMOHKK(1,NHKK))
      VHKK  (3,NHKK) = VHKK(3,JMOHKK(1,NHKK))
      VHKK  (4,NHKK) = VHKK(4,JMOHKK(1,NHKK))
C     WRITE(6,*)' Diffr. Nucleon'
C    *,PHKK(3,NHKK),PHKK(4,NHKK),PHKK(5,NHKK)
*
*---------- S. Roesler 21-10-93
*      check energy-momentum conservation
      IF (IPEV.GE.2) THEN
         PX = PDQ1(1)+PDQ2(1)+PDD1(1)+PDD2(1)+PDFQ1(1)
         PY = PDQ1(2)+PDQ2(2)+PDD1(2)+PDD2(2)+PDFQ1(2)
         PZ = PDQ1(3)+PDQ2(3)+PDD1(3)+PDD2(3)+PDFQ1(3)
         EE = ECM-(PDQ1(4)+PDQ2(4)+PDD1(4)+PDD2(4)+PDFQ1(4))
         WRITE(LOUT,*)'VAHMSD: ENERGY-MOMENTUM-CHECK (PX,PY,PZ,E)'
         WRITE(LOUT,'(5F12.6)')PX,PY,PZ,EE,ECM
      ENDIF
*
      RETURN
*
 9999 CONTINUE
      NCREJ = NCREJ+1
      IF (MOD(NCREJ,2500).EQ.0) WRITE(LOUT,9900) NCREJ
 9900 FORMAT('REJECTION IN VAHMSD ',I10)
      IIREJ = 0
      IREJ  = 1
      RETURN
      END
*
*===diffpt===============================================================*
*
      SUBROUTINE DIFFPT(   ECH1,   PCH1, XDICH1,   XCH1,
     &                     ECH2,   PCH2, XDICH2,   XCH2,
     &                     ECH3,   PCH3,  ECH3I,  KCH3I,
     &                       EE,    EES,  IIREJ)
 
**************************************************************************
* Version November 1993                         by Stefan Roesler        *
*                                                  University of Leipzig *
* This subroutine is called by VAHMSD (HMSD), selects the transverse     *
* momenta of the chains and the diffractive nucleon.                     *
* The results are stored in the common-block ABRDIF.                     *
**************************************************************************
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)
      PARAMETER (NMXHKK= 89998)
      CHARACTER*8 ANAME
      COMMON /HKKEVT/ NHKK,NEVHKK,     ISTHKK(NMXHKK),  IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),PHKK(5,NMXHKK),
     &                VHKK(4,NMXHKK),  WHKK(4,NMXHKK)
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
      COMMON /DIFFRA/ ISINGD,IDIFTP,IOUDIF,IFLAGD
      COMMON /ABRDIF/ XDQ1,XDQ2,XDDQ1,XDDQ2,
     &                IKVQ1,IKVQ2,IKD1Q1,IKD2Q1,IKD1Q2,IKD2Q2,
     &                IDIFFP,IDIFAP,
     &                AMDCH1,AMDCH2,AMDCH3,GAMDC1,GAMDC2,GAMDC3,
     &                PGXVC1,PGYVC1,PGZVC1,PGXVC2,PGYVC2,PGZVC2,
     &                PGXVC3,PGYVC3,PGZVC3,NDCH1,NDCH2,NDCH3,
     &                IKDCH1,IKDCH2,IKDCH3,
     &                PDQ1(4),PDQ2(4),PDD1(4),PDD2(4),PDFQ1(4)
      COMMON /DPAR/   ANAME(210),AM(210),GA(210),TAU(210),ICH(210),
     &                IBAR(210),K1(210),K2(210)
      DATA NCPT /0/
*
      PCH3I = SIGN(SQRT((ECH3I-AM(KCH3I))*(ECH3I+AM(KCH3I))),PCH3)
*
*-------------------- select transverse momenta - diffractive particle
*
      TTTMIN = (ECH3-ECH3I)**2-(PCH3-PCH3I)**2
   10 CONTINUE
*
      XDIFF = XDICH1+XDICH2
      BTP0  = 3.7D0
      ALPH  = 0.24D0
      SLOPE = BTP0-2.0D0*ALPH*LOG(XDIFF)
      Y     = RNDM(V)
      TTT   = -LOG(1.0D0-Y)/SLOPE
      IF (IPEV.GE.2) THEN
         WRITE(LOUT,1009) SLOPE,TTT
 1009    FORMAT('VAHMSD/DIFFPT: SLOPE,TTT ',2E10.5)
      ENDIF
*
      IF (TTT.LE.ABS(TTTMIN)) GOTO 10
      TTT = -ABS(TTT)
      IF (IPEV.GE.2) WRITE(LOUT,1000) PCH3I,TTTMIN,TTT
 1000 FORMAT('VAHMSD/DIFFPT: PCH3I,TTTMIN,TTT ',3F10.5)
      PCH3F = ABS(TTT-(ECH3-ECH3I)**2+PCH3**2+PCH3I**2)/(2*PCH3I)
      IF ((PCH3**2).LE.(PCH3F**2)) THEN
         NCPT = NCPT+1
         IF(MOD(NCPT,5000).EQ.0) WRITE(LOUT,1001) NCPT
 1001    FORMAT('VAHMSD/DIFFPT: INEFFICIENT PT-SELECTION 1, NCPT=',I8)
         GOTO 10
      ENDIF
      P3XYF = SQRT(PCH3**2-PCH3F**2)
      CALL DSFECF(SFE,CFE)
      PXCH3F = P3XYF*CFE
      PYCH3F = P3XYF*SFE
*
*-------------------- select transverse momenta for partons of chain 1
*
      IF (NDCH1.EQ.-99) THEN
         EAQ2   = 0.0D0
         PTXVA2 = 0.0D0
         PTYVA2 = 0.0D0
         PLAQ2  = 0.0D0
         EQ1    = 0.0D0
         PTXDQ1 = 0.0D0
         PTYDQ1 = 0.0D0
         PLQ1   = 0.0D0
         PXCH1F = 0.0D0
         PYCH1F = 0.0D0
         PLCH1F = 0.0D0
         GOTO 11
      ENDIF
      PYCH1F = -ABS(PCH1/PCH3)*PYCH3F
      PXCH1F = -ABS(PYCH1F/PYCH3F)*PXCH3F
      TEMP   = PXCH1F**2+PYCH1F**2
      IF ((PCH1**2).LT.TEMP) THEN
         NCPT = NCPT+1
         IF(MOD(NCPT,5000).EQ.0) WRITE(LOUT,1002) NCPT
 1002    FORMAT('VAHMSD/DIFFPT: INEFFICIENT PT-SELECTION 2, NCPT=',I8)
         GOTO 10
      ENDIF
      PLCH1F = SIGN(SQRT(PCH1**2-TEMP),PCH1)
      EAQ2   = EES*XDICH1
      TEMP   = ABS(EAQ2/PCH1)
      PTXVA2 = -PXCH1F*TEMP
      PTYVA2 = -PYCH1F*TEMP
      PLAQ2  = -PLCH1F*TEMP
      EQ1    = EE*XCH1
      TEMP   = ABS(EQ1/PCH1)
      PTXDQ1 = PXCH1F*TEMP
      PTYDQ1 = PYCH1F*TEMP
      PLQ1   = PLCH1F*TEMP
*
*-------------------- select transverse momenta for partons of chain 2
*
   11 CONTINUE
      IF (NDCH2.EQ.-99) THEN
         EQ2    = 0.0D0
         PTXDQ2 = 0.0D0
         PTYDQ2 = 0.0D0
         PLQ2   = 0.0D0
         EAQ1   = 0.0D0
         PTXVA1 = 0.0D0
         PTYVA1 = 0.0D0
         PLAQ1  = 0.0D0
         PXCH2F = 0.0D0
         PYCH2F = 0.0D0
         PLCH2F = 0.0D0
         GOTO 12
      ENDIF
      PYCH2F = -ABS(PCH2/PCH3)*PYCH3F
      PXCH2F = -ABS(PYCH2F/PYCH3F)*PXCH3F
      TEMP   = PXCH2F**2+PYCH2F**2
      IF ((PCH2**2).LT.TEMP) THEN
         NCPT = NCPT+1
         IF(MOD(NCPT,5000).EQ.0) WRITE(LOUT,1003) NCPT
 1003    FORMAT('VAHMSD/DIFFPT: INEFFICIENT PT-SELECTION 3, NCPT=',I8)
         GOTO 10
      ENDIF
      PLCH2F = SIGN(SQRT(PCH2**2-TEMP),PCH2)
      EQ2    = EES*XDICH2
      TEMP   = ABS(EQ2/PCH2)
      PTXDQ2 = -PXCH2F*TEMP
      PTYDQ2 = -PYCH2F*TEMP
      PLQ2   = -PLCH2F*TEMP
      EAQ1   = EE*XCH2
      TEMP   = ABS(EAQ1/PCH2)
      PTXVA1 = PXCH2F*TEMP
      PTYVA1 = PYCH2F*TEMP
      PLAQ1  = PLCH2F*TEMP
*
   12 CONTINUE
      PTXCH3 = PXCH3F
      PTYCH3 = PYCH3F
      PCH3   = PCH3F
*---------- S. Roesler 11/4/93
*           calculate diffractive x-value
      IF (IPEV.GE.2) THEN
         TEMPX = PTXVA2+PTXDQ1+PTXDQ2+PTXVA1
         TEMPY = PTYVA2+PTYDQ1+PTYDQ2+PTYVA1
         TEMPZ = PLAQ2+PLQ1+PLQ2+PLAQ1
         TEMPE = EAQ1+EAQ2+EQ1+EQ2
	 TEMPM = SQRT(TEMPE**2-TEMPX**2-TEMPY**2-TEMPZ**2)
         TEMPP = TEMPM**2/((EES+EE)**2)
	 WRITE(*,*)'diffractive x-value before energy/momentum corr.',
     &              TEMPP
      ENDIF
*
*---------- S. Roesler 10/26/93
*
*           introduce off-shell partons in order to ensure
*           energy-momentum conservation
*
      DIFFX = PXCH3F+PTXVA2+PTXDQ1+PTXDQ2+PTXVA1
      DIFFY = PYCH3F+PTYVA2+PTYDQ1+PTYDQ2+PTYVA1
      DIFFZ = PCH3F+PLAQ2+PLQ1+PLQ2+PLAQ1
      IF (IPEV.GE.2) THEN
         WRITE(LOUT,1011) DIFFX,DIFFY,DIFFZ
 1011    FORMAT('VAHMSD/DIFFPT: DIFFX,DIFFY,DIFFZ ',3E15.5)
      ENDIF
*
      IF ((NDCH1.EQ.-99).OR.(NDCH1.EQ.1).OR.(NDCH1.EQ.-1)) THEN
*
         PTXDQ2 = PTXDQ2-DIFFX/2.0D0
         PTXVA1 = PTXVA1-DIFFX/2.0D0
*
         PTYDQ2 = PTYDQ2-DIFFY/2.0D0
         PTYVA1 = PTYVA1-DIFFY/2.0D0
*
         PLQ2   = PLQ2  -DIFFZ/2.0D0
         PLAQ1  = PLAQ1 -DIFFZ/2.0D0
*
      ELSEIF ((NDCH2.EQ.-99).OR.(NDCH2.EQ.1).OR.(NDCH2.EQ.-1)) THEN
*
         PTXVA2 = PTXVA2-DIFFX/2.0D0
         PTXDQ1 = PTXDQ1-DIFFX/2.0D0
*
         PTYVA2 = PTYVA2-DIFFY/2.0D0
         PTYDQ1 = PTYDQ1-DIFFY/2.0D0
*
         PLAQ2  = PLAQ2 -DIFFZ/2.0D0
         PLQ1   = PLQ1  -DIFFZ/2.0D0
*
      ELSE
*
         PTXVA2 = PTXVA2-DIFFX/4.0D0
         PTXDQ1 = PTXDQ1-DIFFX/4.0D0
         PTXDQ2 = PTXDQ2-DIFFX/4.0D0
         PTXVA1 = PTXVA1-DIFFX/4.0D0
*
         PTYVA2 = PTYVA2-DIFFY/4.0D0
         PTYDQ1 = PTYDQ1-DIFFY/4.0D0
         PTYDQ2 = PTYDQ2-DIFFY/4.0D0
         PTYVA1 = PTYVA1-DIFFY/4.0D0
*
         PLAQ2  = PLAQ2 -DIFFZ/4.0D0
         PLQ1   = PLQ1  -DIFFZ/4.0D0
         PLQ2   = PLQ2  -DIFFZ/4.0D0
         PLAQ1  = PLAQ1 -DIFFZ/4.0D0
*
      ENDIF
*---------- S. Roesler 11/4/93
*           calculate diffractive x-value
      IF (IPEV.GE.2) THEN
         TEMPX = PTXVA2+PTXDQ1+PTXDQ2+PTXVA1
         TEMPY = PTYVA2+PTYDQ1+PTYDQ2+PTYVA1
         TEMPZ = PLAQ2+PLQ1+PLQ2+PLAQ1
         TEMPE = EAQ1+EAQ2+EQ1+EQ2
	 TEMPM = SQRT(TEMPE**2-TEMPX**2-TEMPY**2-TEMPZ**2)
         TEMPP = TEMPM**2/((EES+EE)**2)
	 WRITE(*,*)'diffractive x-value after energy/momentum corr.',
     &              TEMPP
      ENDIF
*
*           recalculate chain masses...
*
      AMDCH1 = SQRT(ABS((EAQ2+EQ1)**2-(PTXVA2+PTXDQ1)**2
     &         -(PTYVA2+PTYDQ1)**2-(PLAQ2+PLQ1)**2)+1.0D-8)
      AMDCH2 = SQRT(ABS((EQ2+EAQ1)**2-(PTXDQ2+PTXVA1)**2
     &         -(PTYDQ2+PTYVA1)**2-(PLQ2+PLAQ1)**2)+1.0D-8)
*
      IF (IPEV.GE.2) THEN
         WRITE(LOUT,1012) AMDCH1,AMDCH2
 1012    FORMAT('VAHMSD/DIFFPT: AMDCH1,AMDCH2 ',2E10.5)
      ENDIF
*
*           ...and the Lorentz-parameter gamma
*
      GAMDC1 = (EAQ2+EQ1)/AMDCH1
      GAMDC2 = (EQ2+EAQ1)/AMDCH2
      IF (IPEV.GE.2) THEN
         WRITE(LOUT,1013) GAMDC1,GAMDC2
 1013    FORMAT('VAHMSD/DIFFPT: GAMDC1,GAMDC2 ',2E10.5)
      ENDIF
*
*-------------------- store results in common block /ABRDIF/
*
      PDQ1(1)  = PTXDQ1
      PDQ1(2)  = PTYDQ1
      PDQ1(3)  = PLQ1
      PDQ1(4)  = EQ1
      PDD1(1)  = PTXVA1
      PDD1(2)  = PTYVA1
      PDD1(3)  = PLAQ1
      PDD1(4)  = EAQ1
      PDFQ1(1) = PTXCH3
      PDFQ1(2) = PTYCH3
      PDFQ1(3) = PCH3
      PDFQ1(4) = ECH3
      PDQ2(1)  = PTXDQ2
      PDQ2(2)  = PTYDQ2
      PDQ2(3)  = PLQ2
      PDQ2(4)  = EQ2
      PDD2(1)  = PTXVA2
      PDD2(2)  = PTYVA2
      PDD2(3)  = PLAQ2
      PDD2(4)  = EAQ2
      IF (IPEV.GE.2) THEN
         WRITE(LOUT,1004)
 1004    FORMAT('VAHMSD/DIFFPT: PDQ1,PDD1,PDQ2,PDD2,PDFQ1 ')
         WRITE(LOUT,1005)
     &      (PDQ1(J),PDD1(J),PDQ2(J),PDD2(J),PDFQ1(J),J=1,4)
 1005    FORMAT(12X,5F10.5)
      ENDIF
      PGXVC1 = (PTXDQ1+PTXVA2)/AMDCH1
      PGYVC1 = (PTYDQ1+PTYVA2)/AMDCH1
      PGZVC1 = (PLQ1+PLAQ2)/AMDCH1
      PGXVC2 = (PTXVA1+PTXDQ2)/AMDCH2
      PGYVC2 = (PTYVA1+PTYDQ2)/AMDCH2
      PGZVC2 = (PLAQ1+PLQ2)/AMDCH2
      PGXVC3 = PTXCH3/AMDCH3
      PGYVC3 = PTYCH3/AMDCH3
      PGZVC3 = PCH3/AMDCH3
*
      PXCH1F = PTXDQ1+PTXVA2
      PYCH1F = PTYDQ1+PTYVA2
      PLCH1F = PLQ1+PLAQ2
      PXCH2F = PTXVA1+PTXDQ2
      PYCH2F = PTYVA1+PTYDQ2
      PLCH2F = PLAQ1+PLQ2
*
      IF (IPEV.GE.2) THEN
         WRITE(LOUT,1006) PGXVC1,PGYVC1,PGZVC1
 1006    FORMAT('VAHMSD/DIFFPT: PGXVC1,PGYVC1,PGZVC1 ',3F10.5)
         WRITE(LOUT,1007) PGXVC2,PGYVC2,PGZVC2
 1007    FORMAT('VAHMSD/DIFFPT: PGXVC2,PGYVC2,PGZVC2 ',3F10.5)
         WRITE(LOUT,1008) PGXVC3,PGYVC3,PGZVC3
 1008    FORMAT('VAHMSD/DIFFPT: PGXVC3,PGYVC3,PGZVC3 ',3F10.5)
      ENDIF
      IF ((NDCH1.EQ.-99).AND.(NDCH2.EQ.-99)) THEN
         WRITE(LOUT,*) 'REJECT IN DIFFPT: NO CHAINS CREATED'
         IIREJ = 1
         RETURN
      ENDIF
*
*-------------------- store results in common block /HKKEVT/
*
      PHKK(1,NHKK+9) = PXCH1F
      PHKK(2,NHKK+9) = PYCH1F
      PHKK(3,NHKK+9) = PLCH1F
      PHKK(4,NHKK+9) = ECH1
*
      PHKK(1,NHKK+12) = PXCH2F
      PHKK(2,NHKK+12) = PYCH2F
      PHKK(3,NHKK+12) = PLCH2F
      PHKK(4,NHKK+12) = ECH2
*
      RETURN
      END
*
*===diffch================================================================
*
      SUBROUTINE DIFFCH ( XSEA,  IFSEA,   XPAR, IFPARA, IFPARB,
     &                   XSEA2, IFSEA2, XPAR12, XPAR22,     EE,
     &                      RM,    ECH,    PCH,  GAMMA,   BETA,
     &                    NDCH,    ICH,    EES,   NUNO,  IIREJ,   IOPT)
 
**************************************************************************
* Version November 1993                         by Stefan Roesler        *
*                                                  University of Leipzig *
* This subroutine is called by VAHMSD (HMSD) and calculates the kine-    *
* matical parameters of chain IOPT from sampled x-values and flavors     *
* in CMS.                                                                *
**************************************************************************
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)
      CHARACTER*8 ANAME
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
      COMMON /DIFFRA/ ISINGD,IDIFTP,IOUDIF,IFLAGD
      COMMON /DPAR/   ANAME(210),AM(210),GA(210),TAU(210),IICH(210),
     &                IBAR(210),K1(210),K2(210)
      COMMON/DINPDA/IMPS(6,6),IMVE(6,6),IB08(6,21),IB10(6,21),IA08(6,21)
     &,IA10(6,21),A1,B1,B2,B3,LT,LE,BET,AS,B8,AME,DIQ,ISU
      DATA NCNOTH,NCWD /0,0/
*
      GOTO (1,2,3) IOPT
*
*-------------------- kinematical parameters of a q-aq (in any order) chain
*                     XSEA - XPAR, IFSEA - IFPAR
*
    1 CONTINUE
      IIREJ = 0
      IFQ   = IFPARA
      IFAQ  = IABS(IFSEA)
      IF (IFQ.LT.0) THEN
         IFQ  = IFSEA
         IFAQ = IABS(IFPARA)
      ENDIF
      IPS  = IMPS(IFAQ,IFQ)
      IV   = IMVE(IFAQ,IFQ)
      RMPS = AM(IPS)
      RMV  = AM(IV)
      RMBB = RMV+3.0D-1
      NDCH = 0
*
      IF ((XSEA.LE.0.0D0).OR.(XPAR.LE.0.0D0)) THEN
         WRITE(LOUT,*) 'REJECTION IN DIFFCH: 1, XSEA,XPAR ',XSEA,XPAR
         IIREJ=1
         RETURN
      ENDIF
      RM = 2.0D0*SQRT(EES*XSEA*EE*XPAR)
      IF (IPEV.GE.2) THEN
         WRITE(LOUT,1000) 'XSEA, XSEA2, XPAR, XPAR12, RM, EE, EES'
         WRITE(LOUT,1001)  XSEA, XSEA2, XPAR, XPAR12, RM, EE, EES
      ENDIF
*
      IF (RM.LT.RMPS) THEN
*                                 produce nothing
         NDCH   = -99
         XSEA2  = XSEA2 +XSEA
         XPAR12 = XPAR12+XPAR
         IFSEA2 = IFPARA
         XSEA   = 0.0D0
         XPAR   = 0.0D0
         RM     = 1.0D-4
         ECH    = 0.0D0
         PCH    = 0.0D0
         GAMMA  = 1.0D0
         BETA   = 0.0D0
         NCNOTH = NCNOTH+1
         IF (MOD(NCNOTH,5000).EQ.0) WRITE(LOUT,1002) NCNOTH
 1002    FORMAT('VAHMSD/DIFFCH: PRODUCE NOTHING 1, NCNOTH=',I8)
         RETURN
*
      ELSE IF (RM.LT.RMV) THEN
*                                 produce RMPS
         RM     = RMPS
         NDCH   = -1
         ICH    = IPS
         XSQ    = RM/(2.0D0*SQRT(EE*EES))
         XSEAOL = XSEA
         XSEA   = XSQ**2/XPAR
         XPAR22 = XPAR22+XSEAOL-XSEA
*
      ELSE IF (RM.LT.RMBB) THEN
*                                 produce RMV
         RM     = RMV
         NDCH   = 1
         ICH    = IV
         XSQ    = RM/(2.0D0*SQRT(EE*EES))
         XSEAOL = XSEA
         XSEA   = XSQ**2/XPAR
         XPAR22 = XPAR22+XSEAOL-XSEA
      ENDIF
*
      IF (IDIFTP.EQ.1) THEN
         PCH = EES*XSEA-EE*XPAR
         IF (PCH.GE.0.0D0) THEN
            NCWD  = NCWD+1
            IF (MOD(NCWD,500).EQ.0) WRITE(LOUT,1003) NCWD
 1003       FORMAT('VAHMSD/DIFFCH: WRONG SIGN OF PCH (1), NCWD=',I8)
            IIREJ = 1
         ENDIF
      ELSE IF (IDIFTP.EQ.2) THEN
         PCH = EE*XPAR-EES*XSEA
         IF (PCH.LE.0.0D0) THEN
            NCWD  = NCWD+1
            IF (MOD(NCWD,500).EQ.0) WRITE(LOUT,1004) NCWD
 1004       FORMAT('VAHMSD/DIFFCH: WRONG SIGN OF PCH (2), NCWD=',I8)
            IIREJ = 1
         ENDIF
      ELSE
         WRITE(LOUT,*) 'ERROR IN VAHMSD/DIFFCH, IDIFTP=', IDIFTP
      ENDIF
      IF (IPEV.GE.2) THEN
         WRITE(LOUT,1000) 'EES,XSEA,EE,XPAR'
         WRITE(LOUT,1001)  EES,XSEA,EE,XPAR
      ENDIF
      ECH   = EES*XSEA+EE*XPAR
      GAMMA = ECH/RM
      BETA  = PCH/RM
      IF (IPEV.GE.2) THEN
         WRITE(LOUT,1000) 'RM,ECH,PCH,NDCH'
         WRITE(LOUT,1001)  RM,ECH,PCH
         WRITE(LOUT,'(I4)') NDCH
      ENDIF
 1000 FORMAT('VAHMSD/DIFFCH-CHAIN Q-AQ',A34)
 1001 FORMAT(10F10.5)
      RETURN
*
*-------------------- kinematical parameters of a q-qq or aq-aqaq (in any
*                     order) chain
*                     XSEA - XPAR, IFSEA - IFPARA/IFPARB
*
    2 CONTINUE
      IIREJ = 0
      CALL DBKLAS(IFSEA, IFPARA, IFPARB, I8, I10)
      RM8   = AM( I8)
      RM10  = AM(I10)
      RMBB  = RM10+3.0D-1
      NDCH  = 0
*
      IF ((XSEA.LE.0.0D0).OR.(XPAR.LE.0.0D0)) THEN
         IIREJ=1
         WRITE(LOUT,*) 'REJECTION IN DIFFCH: 3, XSEA,XPAR ',XSEA,XPAR
         RETURN
      ENDIF
      RM = 2.0D0*SQRT(EES*XSEA*EE*XPAR)
      IF (IPEV.GE.2) THEN
         WRITE(LOUT,2000) 'XSEA, XSEA2, XPAR, XPAR12, RM, EE, EES'
         WRITE(LOUT,2001)  XSEA, XSEA2, XPAR, XPAR12, RM, EE, EES
      ENDIF
*
      IF (RM.LT.RM10) THEN
*                                         reject
         IIREJ  = 1
         NCNOTH = NCNOTH+1
         IF (MOD(NCNOTH,20).EQ.0) WRITE(LOUT,2002) NCNOTH
 2002    FORMAT('VAHMSD/DIFFCH: PRODUCE NOTHING 2, NCNOTH=',I8)
         RETURN
      ELSE IF (RM.LT.RMBB) THEN
*                                         produce RM10
         NDCH   = 1
         RM     = RM10
         ICH    = I10
         XSQ    = RM/(2.0D0*SQRT(EE*EES))
         XSEAOL = XSEA
         XSEA   = XSQ**2/XPAR
         XPAR22 = XPAR22+XSEAOL-XSEA
      ENDIF
*
      IF (IDIFTP.EQ.1) THEN
         PCH = EES*XSEA-EE*XPAR
         IF (PCH.GE.0.0D0) THEN
            NCWD  = NCWD+1
            IF (MOD(NCWD,500).EQ.0) WRITE(LOUT,2003) NCWD
 2003       FORMAT('VAHMSD/DIFFCH: WRONG SIGN OF PCH (3), NCWD=',I8)
            IIREJ = 1
         ENDIF
      ELSE IF (IDIFTP.EQ.2) THEN
         PCH = EE*XPAR-EES*XSEA
         IF (PCH.LE.0.0D0) THEN
            NCWD  = NCWD+1
            IF (MOD(NCWD,500).EQ.0) WRITE(LOUT,2004) NCWD
 2004       FORMAT('VAHMSD/DIFFCH: WRONG SIGN OF PCH (4), NCWD=',I8)
            IIREJ = 1
         ENDIF
      ELSE
         WRITE(LOUT,*) 'ERROR IN VAHMSD/DIFFCH, IDIFTP=', IDIFTP
      ENDIF
      ECH   = EES*XSEA+EE*XPAR
      GAMMA = ECH/RM
      BETA  = PCH/RM
      IF (IPEV.GE.2) THEN
         WRITE(LOUT,2000) 'RM,ECH,PCH,NDCH'
         WRITE(LOUT,2001)  RM,ECH,PCH
         WRITE(LOUT,'(I4)') NDCH
      ENDIF
 2000 FORMAT('VAHMSD/DIFFCH-CHAIN Q-QQ',A34)
 2001 FORMAT(10F10.5)
      RETURN
*
*-------------------- kinematical parameters of a baryon/meson
*
    3 CONTINUE
      IIREJ = 0
      IF (IFPARB.EQ.99) THEN
         IF (IFSEA.LT.0) THEN
            IFAQ = IABS(IFSEA)
            IFQ  = IFPARA
         ELSE
            IFAQ = IABS(IFPARA)
            IFQ  = IFSEA
         ENDIF
         IM = IMPS(IFAQ,IFQ)
      ELSE
         CALL DBKLAS(IFSEA, IFPARA, IFPARB, IM, IDUM)
      ENDIF
      RM   = AM(IM)
      NDCH = -1
      IF ((XSEA.LE.0.0D0).OR.(XPAR.LE.0.0D0)) THEN
         IIREJ=1
         WRITE(LOUT,*) 'REJECTION IN DIFFCH: 5, XSEA,XPAR ',XSEA,XPAR
         RETURN
      ENDIF
      IF (IPEV.GE.2) THEN
         WRITE(LOUT,3000) 'XSEA, XPAR, RM'
         WRITE(LOUT,3001)  XSEA, XPAR, RM
      ENDIF
      ECH = EE*(XSEA+XPAR)
      IF (IDIFTP.EQ.1) THEN
         PCH = SQRT((ECH-RM)*(ECH+RM))
      ELSE IF (IDIFTP.EQ.2) THEN
         PCH = -SQRT((ECH-RM)*(ECH+RM))
      ELSE
         WRITE(LOUT,*) 'ERROR IN VAHMSD/DIFFCH, IDIFTP=', IDIFTP
      ENDIF
      GAMMA = ECH/RM
      BETA  = PCH/RM
      IF (IPEV.GE.2) THEN
         WRITE(LOUT,3000) 'ECH, PCH, GAMMA, BETA'
         WRITE(LOUT,3001)  ECH, PCH, GAMMA, BETA
      ENDIF
 3000 FORMAT('VAHMSD/DIFFCH-BARYON/MESON ',A27)
 3001 FORMAT(10F10.5)
      RETURN
      END
*
*===valmsd===============================================================*
*
      SUBROUTINE VALMSD(ITAPOI,ECM,KPROJ,KTARG,IREJ)
 
**************************************************************************
* Version November 1993                         by Stefan Roesler        *
*                                                  University of Leipzig *
* This subroutine selects flavors and 4-momenta of partons in low-mass   *
* single diffractive chains.                                             *
**************************************************************************
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)
      PARAMETER (NMXHKK= 89998)
      CHARACTER*8 ANAME
      COMMON /HKKEVT/ NHKK,NEVHKK,     ISTHKK(NMXHKK),  IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),PHKK(5,NMXHKK),
     &                VHKK(4,NMXHKK),  WHKK(4,NMXHKK)
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
      COMMON /DIFFRA/ ISINGD,IDIFTP,IOUDIF,IFLAGD
      COMMON /ABRDIF/ XDQ1,XDQ2,XDDQ1,XDDQ2,
     &                IKVQ1,IKVQ2,IKD1Q1,IKD2Q1,IKD1Q2,IKD2Q2,
     &                IDIFFP,IDIFAP,
     &                AMDCH1,AMDCH2,AMDCH3,GAMDC1,GAMDC2,GAMDC3,
     &                PGXVC1,PGYVC1,PGZVC1,PGXVC2,PGYVC2,PGZVC2,
     &                PGXVC3,PGYVC3,PGZVC3,NDCH1,NDCH2,NDCH3,
     &                IKDCH1,IKDCH2,IKDCH3,
     &                PDQ1(4),PDQ2(4),PDD1(4),PDD2(4),PDFQ1(4)
      COMMON /DPAR/   ANAME(210),AM(210),GA(210),TAU(210),ICH(210),
     &                IBAR(210),K1(210),K2(210)
      COMMON/DINPDA/IMPS(6,6),IMVE(6,6),IB08(6,21),IB10(6,21),IA08(6,21)
     &,IA10(6,21),A1,B1,B2,B3,LT,LE,BET,AS,B8,AME,DIQ,ISU
      COMMON /TRAFOP/ GAMP,BGAMP,BETP
      COMMON /ENERIN/ EPROJ,ETARG
      COMMON /SDFLAG/ ISD
      COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
      COMMON /XDIDID/XDIDI
*
      DIMENSION MQUARK(3,30),IHKKQ(-6:6),IHKKQQ(-3:3,-3:3),
     &          IDX(-4:4)
      DATA IDX   /-4,-3,-1,-2,0,2,1,3,4/
      DATA IHKKQ /-6,-5,-4,-3,-1,-2,0,2,1,3,4,5,6/
      DATA IHKKQQ/-3301,-3103,-3203,0,   0,0,0,
     &            -3103,-1103,-2103,0,   0,   0,   0,
     &            -3203,-2103,-2203,0,   0,   0,   0,
     &                0,    0,    0,0,   0,   0,   0,
     &                0,    0,    0,0,2203,2103,3202,
     &                0,    0,    0,0,2103,1103,3103,
     &                0,    0,    0,0,3203,3103,3303/
*
*----------------------------------- quark content of hadrons:
*                                     1, 2, 3, 4  -   u, d, s, c
*                                    -1,-2,-3,-4  -  au,ad,as,ac
*
      DATA MQUARK/
     &   1,1,2,    -1,-1,-2,       0,0,0,       0,0,0,       0,0,0,
     &   0,0,0,       0,0,0,       1,2,2,    -1,-2,-2,       0,0,0,
     &   0,0,0,       0,0,0,      1,-2,0,      2,-1,0,      1,-3,0,
     &  3,-1,0,       1,2,3,    -1,-2,-3,       0,0,0,       2,2,3,
     &   1,1,3,       1,2,3,      1,-1,0,      2,-3,0,      3,-2,0,
     &  2,-2,0,      3,-3,0,       0,0,0,       0,0,0,       0,0,0/
      DATA UNON/2.0/
      DATA NCREJ, NCPT /0, 0/
*
      ISD   = 2
      IREJ  = 0
      IIREJ = 0
      IBPROJ = IBAR(KPROJ)
      IBTARG = IBAR(KTARG)
      EPROJ  = (AM(KPROJ)**2-AM(KTARG)**2+ECM**2)/(2.0D0*ECM)
      ETARG  = (AM(KTARG)**2-AM(KPROJ)**2+ECM**2)/(2.0D0*ECM)
      IF(IPEV.GE.2) WRITE(LOUT,1014)EPROJ,ETARG
 1014 FORMAT('VALMSD: EPROJ,ETARG',2F10.5)
      IF (IBTARG.LE.0) THEN
         WRITE(LOUT,1001) IBTARG
 1001    FORMAT('VALMSD: NO LMSD FOR TARGET WITH BARYON-CHARGE',I4)
         IIREJ = 1
         GOTO 9999
      ENDIF
      IQP1 = MQUARK(1,KPROJ)
      IQP2 = MQUARK(2,KPROJ)
      IQP3 = MQUARK(3,KPROJ)
      IQT1 = MQUARK(1,KTARG)
      IQT2 = MQUARK(2,KTARG)
      IQT3 = MQUARK(3,KTARG)
      IF(IPEV.GE.2) WRITE(LOUT,1002)
     &                IBPROJ,IBTARG,IQP1,IQP2,IQP3,IQT1,IQT2,IQT3
 1002 FORMAT('VALMSD: IBPROJ,IBTARG,IQP1,IQP2,IQP3,IQT1,IQT2,IQT3 ',8I4)
*
      IF (IBPROJ.NE.0) THEN
*
*-------------------- q-qq (aq-aqaq) - flavors of projectile (baryon)
*
         ISAM = 1.0D0+2.999D0*RNDM(V)
         GOTO (10,11,12) ISAM
   10    CONTINUE
         IQP    = IQP1
         IDIQP1 = IQP2
         IDIQP2 = IQP3
         GOTO 13
   11    CONTINUE
         IQP    = IQP2
         IDIQP1 = IQP1
         IDIQP2 = IQP3
         GOTO 13
   12    CONTINUE
         IQP    = IQP3
         IDIQP1 = IQP1
         IDIQP2 = IQP2
   13    CONTINUE
*
      ELSE IF (IBPROJ.EQ.0) THEN
*
*-------------------- q-aq - flavors of projectile (meson)
*
         ISAM = 1.0D0+1.999D0*RNDM(V)
         GOTO (14,15) ISAM
   14    CONTINUE
         IQP    = IQP1
         IDIQP1 = IQP2
         GOTO 16
   15    CONTINUE
         IQP    = IQP2
         IDIQP1 = IQP1
   16    CONTINUE
         IDIQP2 = 0
*
      ENDIF
*
*-------------------- q-qq - flavors of target (baryon)
*
      ISAM = 1.0D0+2.999D0*RNDM(V)
      GOTO (17,18,19) ISAM
   17 CONTINUE
      IQT    = IQT1
      IDIQT1 = IQT2
      IDIQT2 = IQT3
      GOTO 20
   18 CONTINUE
      IQT    = IQT2
      IDIQT1 = IQT1
      IDIQT2 = IQT3
      GOTO 20
   19 CONTINUE
      IQT    = IQT3
      IDIQT1 = IQT1
      IDIQT2 = IQT2
   20 CONTINUE
*
      IKVQ1  = IQP
      IKD1Q1 = IDIQP1
      IKD2Q1 = IDIQP2
*
      IKVQ2  = IQT
      IKD1Q2 = IDIQT1
      IKD2Q2 = IDIQT2
*
      IF (IPEV.GE.2) WRITE(LOUT,1003)
     &                 IKVQ1,IKD1Q1,IKD2Q1,IKVQ2,IKD1Q2,IKD2Q2
 1003 FORMAT('VALMSD: IKVQ1,IKD1Q1,IKD2Q1,IKVQ2,IKD1Q2,IKD2Q2 ',6I4)
*
*-------------------- IDIFTP = 1  target     (backward) hadron  excited
*                     IDIFTP = 2  projectile (forward)  hadron  excited
*
      IDIFTP = 1.0D0+1.999D0*RNDM(V)
*-------------------- S.Roesler 5/26/93
      IF ((ISINGD.EQ.3).OR.(ISINGD.EQ.7)) IDIFTP = 1
      IF ((ISINGD.EQ.4).OR.(ISINGD.EQ.8)) IDIFTP = 2
*
      IF (IPEV.GE.2) WRITE(LOUT,1004) IDIFTP
 1004 FORMAT('VALMSD: IDIFTP ',I4)
      IF ((IDIFTP.NE.1).AND.(IDIFTP.NE.2)) THEN
         IF (IPEV.GE.2) WRITE(LOUT,'(A19)') 'VALMSD-ERROR: IDIFTP'
         GOTO 9999
      ENDIF
*
*-------------------- diffractive mass
*
31    CONTINUE
      AMO = 1.5D0
      R   = RNDM(V)
      IF ((IBPROJ.EQ.0).AND.(IDIFTP.EQ.1))
     &               AMO = 1.5D0*R+2.83D0*(1.0D0-R)
      IF ((IBPROJ.EQ.0).AND.(IDIFTP.EQ.2)) AMO = 1.0D0
      SAM = 1.0D0
      IF (ECM.LE.300.0D0) SAM = 1.0D0-EXP(-((ECM/200.0D0)**4))
      R   = RNDM(V)*SAM
      AMAX= (1.0D0-SAM)*SQRT(0.1D0*ECM**2)+SAM*SQRT(400.0D0)
      AMU = R*SQRT(100.0D0)+(1.0D0-R)*AMAX
      IF (IBPROJ.EQ.0) THEN
C----------------- change mass-cuts
C        AMAX= (1.0D0-SAM)*SQRT(0.25D0*ECM**2)+SAM*SQRT(400.0D0)
         AMAX= (1.0D0-SAM)*SQRT(0.15D0*ECM**2)+SAM*SQRT(400.0D0)
         AMU = R*SQRT(100.0D0)+(1.0D0-R)*AMAX
      ENDIF
C---------------------------------------------------------------
C
C                                           j.r. test 2/94
C
C---------------------------------------------------------------
      AMU=2.D0*AMU
C---------------------------------------------------------------
      R   = RNDM(V)
      IF (ECM.LE.50.0D0) THEN
         AMDIFF = AMO*(AMU/AMO)**R
      ELSE
         A = 0.7D0
         IF (ECM.LE.300.0D0) A = 0.7D0*(1.0D0-EXP(-((ECM/100.0D0)**2)))
         AMDIFF = 1.0D0/((R/(AMU**A)+(1.0D0-R)/(AMO**A))
     &         **(1.0D0/A))
      ENDIF
      IF(AMDIFF.GT.0.5D0*ECM)GO TO 31
      IF (IOUXEV.GE.2) WRITE(LOUT,1005) AMDIFF
 1005 FORMAT('VALMSD: AMDIFF',E10.5)
       XDIDI=AMDIFF**2/ECM**2
      IF (IOUXEV.GE.2) WRITE(LOUT,*)'LM AMDIFF,XDIDI ',AMDIFF,XDIDI
*
*-------------------- kinematical parameters of
*                     diffractive projectile (IDIFTP = 1)
*                     or diffractive target  (IDIFTP = 2)
*
      IF (IDIFTP.EQ.1) THEN
         AMDCH3 = AM(KPROJ)
      ELSE IF (IDIFTP.EQ.2) THEN
         AMDCH3 = AM(KTARG)
      ENDIF
      ECH3 = (ECM**2+AMDCH3**2-AMDIFF**2)/(2.0D0*ECM)
      IF (ECH3.LE.AMDCH3) THEN
         NCMERR = NCMERR+1
         IF (MOD(NCMERR,2200).EQ.0)
     &      WRITE(LOUT,*) 'LMSD: INEFFICIENT SELECTION OF AMDIFF(1),
     &                        NCMERR = ',NCMERR
         GOTO 31
      ENDIF
*
      PCH3 = SQRT(ABS(ECH3-AMDCH3))*SQRT(ECH3+AMDCH3)
      IF (IDIFTP.EQ.2) PCH3 = -PCH3
      NDCH3  = -1
      GAMDC3 = ECH3/AMDCH3
      PGVC3  = PCH3/AMDCH3
      IF (IOUXEV.GE.2) WRITE(LOUT,1006) AMDCH3,ECH3,PCH3,GAMDC3,PGVC3
 1006 FORMAT('VALMSD: AMDCH3,ECH3,PCH3,GAMDC3,PGVC3',5E15.5)
*
   30 CONTINUE
      B33   = 8.0D0
      ES    = -2.0D0/(B33**2)*LOG(RNDM(V)*RNDM(V))
      HPS   = SQRT(ES*ES+2.0D0*ES*0.94D0)
      CALL DSFECF(SFE,CFE)
      PXCH3 = HPS*CFE
      PYCH3 = HPS*SFE
      PTCH3 = SQRT(PXCH3**2+PYCH3**2)
      IF (PTCH3.GT.ABS(PCH3)) THEN
         NCPT = NCPT+1
         IF (MOD(NCPT,500).EQ.0) WRITE(LOUT,1007) NCPT
 1007    FORMAT('VALMSD: INEFFICIENT PT-SELECTION 1, NCPT=',I8)
         GOTO 30
      ENDIF
      PLCH3 = SIGN(SQRT(ABS(PCH3-PTCH3))*SQRT(ABS(PCH3+PTCH3)),PCH3)
      IF (IOUXEV.GE.2) WRITE(LOUT,1008) ES,HPS,PXCH3,PYCH3,PLCH3
 1008 FORMAT('VALMSD: ES,HPS,PXCH3,PYCH3,PLCH3',5E15.5)
*
*-------------------- no chain 1
*
      NDCH1  = -99
      AMDCH1 = 0.0D0
      ECH1   = 0.0D0
      PCH1   = 0.0D0
      GAMDC1 = 1.0D0
      PGVC1  = 0.0D0
      PGXVC1 = 0.0D0
      PGYVC1 = 0.0D0
      PGZVC1 = 0.0D0
      DO 40 I=1,4
         PDD2(I) = 0.0D0
         PDQ1(I) = 0.0D0
   40 CONTINUE
*
*-------------------- kinematical parameters of
*                     excited target        (IDIFTP = 1)
*                     or excited projectile (IDIFTP = 2)
*
      AMDCH2 = AMDIFF
      ECH2   = (ECM**2+AMDCH2**2-AMDCH3**2)/(2.0D0*ECM)
      PCH2   = -SIGN(SQRT(ABS(ECH2-AMDCH2))*SQRT(ECH2+AMDCH2),PCH3)
      NDCH2  = 0
      GAMDC2 = ECH2/AMDCH2
      PGVC2  = PCH2/AMDCH2
*
*-------------------- (IDIFFP,IDIFAP in analogy to HMSD in VAHMSD)
*                     IDIFTP = 1 : IKD1Q2/IKD2Q2 - IDIFFP = IKVQ2
*                     IDIFTP = 2 : projectile - baryon
*                                  IKD1Q1/IKD2Q1 - IDIFFP = IKVQ1
*                                  projectile - antibaryon
*                                  IKD1Q1/IKD2Q1 - IDIFAP = IKVQ1
*                                  projectile - meson
*                                  IKD1Q1 > 0    - IDIFAP = IKVQ1
*                                  IKD1Q1 < 0    - IDIFFP = IKVQ1
*
      IF (IDIFTP.EQ.1) IDIFFP = IKVQ2
      IF (IDIFTP.EQ.2) THEN
         IF (IBPROJ.GT.0) IDIFFP = IKVQ1
         IF (IBPROJ.LT.0) IDIFAP = IKVQ1
         IF ((IBPROJ.EQ.0).AND.(IKD1Q1.GT.0)) IDIFAP = IKVQ1
         IF ((IBPROJ.EQ.0).AND.(IKD1Q1.LT.0)) IDIFFP = IKVQ1
      ENDIF
      IF (IOUXEV.GE.2) WRITE(LOUT,1009) AMDCH2,ECH2,PCH2,GAMDC2,PGVC2,
     &                                  IDIFFP,IDIFAP
 1009 FORMAT('VALMSD: AMDCH2,ECH2,PCH2,GAMDC2,PGVC2,IDIFFP,IDIFAP',
     &                                                  5E15.5,2I4)
      PXCH2 = -PXCH3
      PYCH2 = -PYCH3
      PTCH2 = SQRT(PXCH2**2+PYCH2**2)
      IF (PTCH2.GT.ABS(PCH2)) THEN
         NCPT = NCPT+1
         IF (MOD(NCPT,500).EQ.0) WRITE(LOUT,1010) NCPT
 1010    FORMAT('VALMSD: INEFFICIENT PT-SELECTION 2, NCPT=',I8)
         GOTO 30
      ENDIF
      PLCH2 = SIGN(SQRT(ABS(PCH2-PTCH2))*SQRT(ABS(PCH2+PTCH2)),PCH2)
      IF (IOUXEV.GE.2) WRITE(LOUT,1011) PXCH2,PYCH2,PLCH2
 1011 FORMAT('VALMSD: PXCH2,PYCH2,PLCH2',3E15.5)
      CC    = AMDCH2/(2.0D0*ECH2)
*--------------------------- S.R. 11/12/92
      IF (CC.GE.0.5D0) THEN
         NCMERR = NCMERR+1
         IF (MOD(NCMERR,200).EQ.0)
     &      WRITE(LOUT,*) 'LMSD: INEFFICIENT SELECTION OF AMDIFF(2),
     &                        NCMERR = ',NCMERR
         GOTO 31
      ENDIF
*
      XDIQ  = 0.5D0+SQRT(ABS((0.5D0-CC)*(0.5D0+CC)))
      XQ    = 1.0D0-XDIQ
      EDIQ  = ECH2*XDIQ
      EQ    = ECH2*XQ
      TEMP  = ABS(EDIQ/PCH2)
      PXDIQ = PXCH2*TEMP
      PYDIQ = PYCH2*TEMP
      PLDIQ = PLCH2*TEMP
      TEMP  = ABS(EQ/PCH2)
      PXQ   = -PXCH2*TEMP
      PYQ   = -PYCH2*TEMP
      PLQ   = -PLCH2*TEMP
*
*-------------------- store results in COMMON-block /ABRDIF/
*
      PDQ2(1)  = PXQ
      PDQ2(2)  = PYQ
      PDQ2(3)  = PLQ
      PDQ2(4)  = EQ
      PDD1(1)  = PXDIQ
      PDD1(2)  = PYDIQ
      PDD1(3)  = PLDIQ
      PDD1(4)  = EDIQ
      PDFQ1(1) = PXCH3
      PDFQ1(2) = PYCH3
      PDFQ1(3) = PLCH3
      PDFQ1(4) = ECH3
      PGXVC2   = PXCH2/AMDCH2
      PGYVC2   = PYCH2/AMDCH2
      PGZVC2   = PLCH2/AMDCH2
      PGXVC3   = PXCH3/AMDCH3
      PGYVC3   = PYCH3/AMDCH3
      PGZVC3   = PLCH3/AMDCH3
*
      IF (IPEV.GE.2) THEN
         WRITE(LOUT,*) 'VALMSD: PDQ1, PDD2, PDD1, PDQ2, PDFQ1'
         WRITE(LOUT,1012)
     &      (PDQ1(I),PDD2(I),PDD1(I),PDQ2(I),PDFQ1(I),I=1,4)
 1012    FORMAT(5E15.5)
         WRITE(LOUT,1013) 'PGXVC1,PGYVC1,PGZVC1,GAMDC1',
     &                     PGXVC1,PGYVC1,PGZVC1,GAMDC1
         WRITE(LOUT,1013) 'PGXVC2,PGYVC2,PGZVC2,GAMDC2',
     &                     PGXVC2,PGYVC2,PGZVC2,GAMDC2
         WRITE(LOUT,1013) 'PGXVC3,PGYVC3,PGZVC3,GAMDC3',
     &                     PGXVC3,PGYVC3,PGZVC3,GAMDC3
 1013    FORMAT(A27,4E15.5)
      ENDIF
*
*------------------- store results in common block /HKKEVT/
*
*                    partonen of projectile/target
*
      NHKK = NHKK + 1
      ISTHKK(NHKK) = 23-IDIFTP
      IF (IDIFTP.EQ.1) IDHKK(NHKK) = IHKKQ(IDX(IKVQ2))
      IF (IDIFTP.EQ.2) IDHKK(NHKK) = IHKKQ(IDX(IKVQ1))
      IF (IDIFTP.EQ.1) JMOHKK(1,NHKK) = ITAPOI
      IF (IDIFTP.EQ.2) JMOHKK(1,NHKK) = 1
      JMOHKK(2,NHKK) = 0
      JDAHKK(1,NHKK) = 0
      JDAHKK(2,NHKK) = 0
      PHKK(1,NHKK)   = 0.0D0
      PHKK(2,NHKK)   = 0.0D0
      PHKK(3,NHKK)   = XQ
      PHKK(4,NHKK)   = XQ
      PHKK(5,NHKK)   = 0.0D0
      VHKK  (1,NHKK) = VHKK(1,JMOHKK(1,NHKK))
      VHKK  (2,NHKK) = VHKK(2,JMOHKK(1,NHKK))
      VHKK  (3,NHKK) = VHKK(3,JMOHKK(1,NHKK))
      VHKK  (4,NHKK) = VHKK(4,JMOHKK(1,NHKK))
*
      NHKK = NHKK + 1
      ISTHKK(NHKK) = 23-IDIFTP
      IF (IDIFTP.EQ.1) IDHKK(NHKK) = IHKKQQ(IDX(IKD1Q2),IDX(IKD2Q2))
      IF (IDIFTP.EQ.2) IDHKK(NHKK) = IHKKQQ(IDX(IKD1Q1),IDX(IKD2Q1))
      IF (IDIFTP.EQ.1) JMOHKK(1,NHKK) = ITAPOI
      IF (IDIFTP.EQ.2) JMOHKK(1,NHKK) = 1
      JMOHKK(2,NHKK) = 0
      JDAHKK(1,NHKK) = 0
      JDAHKK(2,NHKK) = 0
      PHKK(1,NHKK)   = 0.0D0
      PHKK(2,NHKK)   = 0.0D0
      PHKK(3,NHKK)   = XDIQ
      PHKK(4,NHKK)   = XDIQ
      PHKK(5,NHKK)   = 0.0D0
      VHKK  (1,NHKK) = VHKK(1,JMOHKK(1,NHKK))
      VHKK  (2,NHKK) = VHKK(2,JMOHKK(1,NHKK))
      VHKK  (3,NHKK) = VHKK(3,JMOHKK(1,NHKK))
      VHKK  (4,NHKK) = VHKK(4,JMOHKK(1,NHKK))
*
*                    ends of chain 2
*
      NHKK = NHKK + 1
      ISTHKK(NHKK)   = 100+ISTHKK(NHKK-2)
      IDHKK(NHKK)    = IDHKK(NHKK-2)
      JMOHKK(1,NHKK) = NHKK-2
      JMOHKK(2,NHKK) = JMOHKK(1,NHKK-2)
      JDAHKK(1,NHKK) = NHKK+2
      JDAHKK(2,NHKK) = NHKK+2
      PHKK(1,NHKK)   = PDQ2(1)
      PHKK(2,NHKK)   = PDQ2(2)
      PHKK(3,NHKK)   = PDQ2(3)
      PHKK(4,NHKK)   = PDQ2(4)
      PHKK(5,NHKK)   = 0.0D0
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
      VHKK  (1,NHKK) = VHKK(1,JMOHKK(1,NHKK))+XXPP
      VHKK  (2,NHKK) = VHKK(2,JMOHKK(1,NHKK))+YYPP
      VHKK  (3,NHKK) = VHKK(3,JMOHKK(1,NHKK))
      VHKK  (4,NHKK) = VHKK(4,JMOHKK(1,NHKK))
*
      NHKK = NHKK + 1
      ISTHKK(NHKK)   = 100+ISTHKK(NHKK-2)
      IDHKK(NHKK)    = IDHKK(NHKK-2)
      JMOHKK(1,NHKK) = NHKK-2
      JMOHKK(2,NHKK) = JMOHKK(1,NHKK-2)
      JDAHKK(1,NHKK) = NHKK+1
      JDAHKK(2,NHKK) = NHKK+1
      PHKK(1,NHKK)   = PDD1(1)
      PHKK(2,NHKK)   = PDD1(2)
      PHKK(3,NHKK)   = PDD1(3)
      PHKK(4,NHKK)   = PDD1(4)
      PHKK(5,NHKK)   = 0.0D0
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
      VHKK  (1,NHKK) = VHKK(1,JMOHKK(1,NHKK))+XXPP
      VHKK  (2,NHKK) = VHKK(2,JMOHKK(1,NHKK))+YYPP
      VHKK  (3,NHKK) = VHKK(3,JMOHKK(1,NHKK))
      VHKK  (4,NHKK) = VHKK(4,JMOHKK(1,NHKK))
*
*                    chain 2
*
      NHKK = NHKK + 1
      ISTHKK(NHKK)   = 177
      IDHKK(NHKK)    = 88888
      JMOHKK(1,NHKK) = NHKK-2
      JMOHKK(2,NHKK) = NHKK-1
      JDAHKK(1,NHKK) = 0
      JDAHKK(2,NHKK) = 0
      PHKK(1,NHKK)   = PXCH2
      PHKK(2,NHKK)   = PYCH2
      PHKK(3,NHKK)   = PLCH2
      PHKK(4,NHKK)   = ECH2
      PHKK(5,NHKK)   = AMDCH2
      VHKK  (1,NHKK) = VHKK(1,NHKK-1)
      VHKK  (2,NHKK) = VHKK(2,NHKK-1)
      VHKK  (3,NHKK) = VHKK(3,NHKK-1)
      IF ((BETP.NE.0.0D0).AND.(BGAMP.NE.0.0D0))
     &VHKK  (4,NHKK) = VHKK(3,NHKK)/BETP-VHKK(3,NHKK-2)/BGAMP
*
*                    diffractive nucleon
*
      NHKK = NHKK + 1
      ISTHKK(NHKK)   = 177
      IDHKK(NHKK)    = 88888
      IF (IDIFTP.EQ.2) JMOHKK(1,NHKK) = ITAPOI
      IF (IDIFTP.EQ.1) JMOHKK(1,NHKK) = 1
      JMOHKK(2,NHKK) = 0
      JDAHKK(1,NHKK) = 0
      JDAHKK(2,NHKK) = 0
      PHKK(1,NHKK)   = PDFQ1(1)
      PHKK(2,NHKK)   = PDFQ1(2)
      PHKK(3,NHKK)   = PDFQ1(3)
      PHKK(4,NHKK)   = PDFQ1(4)
C     WRITE(6,*)' Diffr. Nucleon'
C    *,PHKK(3,NHKK),PHKK(4,NHKK),PHKK(5,NHKK)
      PHKK(5,NHKK)   = AMDCH3
      VHKK  (1,NHKK) = VHKK(1,JMOHKK(1,NHKK))
      VHKK  (2,NHKK) = VHKK(2,JMOHKK(1,NHKK))
      VHKK  (3,NHKK) = VHKK(3,JMOHKK(1,NHKK))
      VHKK  (4,NHKK) = VHKK(4,JMOHKK(1,NHKK))
*
*
*---------- S. Roesler 21-10-93
*      check energy-momentum conservation
      IF (IPEV.GE.2) THEN
         PX = PDQ1(1)+PDQ2(1)+PDD1(1)+PDD2(1)+PDFQ1(1)
         PY = PDQ1(2)+PDQ2(2)+PDD1(2)+PDD2(2)+PDFQ1(2)
         PZ = PDQ1(3)+PDQ2(3)+PDD1(3)+PDD2(3)+PDFQ1(3)
         EE = ECM-(PDQ1(4)+PDQ2(4)+PDD1(4)+PDD2(4)+PDFQ1(4))
         WRITE(LOUT,*)'VALMSD: ENERGY-MOMENTUM-CHECK (PX,PY,PZ,E)'
         WRITE(LOUT,'(5F12.6)')PX,PY,PZ,EE,ECM
      ENDIF
*
      RETURN
 9999 CONTINUE
      NCREJ = NCREJ+1
      IF (MOD(NCREJ,500).EQ.0) WRITE(LOUT,9900) NCREJ
 9900 FORMAT('REJECTION IN VALMSD ',I10)
      IIREJ = 0
      IREJ  = 1
      RETURN
      END
*
*===valmdd===============================================================*
*
      SUBROUTINE VALMDD(ITAPOI,ECM,KPROJ,KTARG,IREJ)

C                   Low mass double diffraction J.R. Feb/Mar 94
 
**************************************************************************
* Version November 1993                         by Stefan Roesler        *
*                                                  University of Leipzig *
* This subroutine selects flavors and 4-momenta of partons in low-mass   *
* single diffractive chains.                                             *
**************************************************************************
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)
      PARAMETER (NMXHKK= 89998)
      CHARACTER*8 ANAME
      COMMON /HKKEVT/ NHKK,NEVHKK,     ISTHKK(NMXHKK),  IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),PHKK(5,NMXHKK),
     &                VHKK(4,NMXHKK),  WHKK(4,NMXHKK)
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
      COMMON /DIFFRA/ ISINGD,IDIFTP,IOUDIF,IFLAGD
      COMMON /ABRDIF/ XDQ1,XDQ2,XDDQ1,XDDQ2,
     &                IKVQ1,IKVQ2,IKD1Q1,IKD2Q1,IKD1Q2,IKD2Q2,
     &                IDIFFP,IDIFAP,
     &                AMDCH1,AMDCH2,AMDCH3,GAMDC1,GAMDC2,GAMDC3,
     &                PGXVC1,PGYVC1,PGZVC1,PGXVC2,PGYVC2,PGZVC2,
     &                PGXVC3,PGYVC3,PGZVC3,NDCH1,NDCH2,NDCH3,
     &                IKDCH1,IKDCH2,IKDCH3,
     &                PDQ1(4),PDQ2(4),PDD1(4),PDD2(4),PDFQ1(4)
      COMMON /DPAR/   ANAME(210),AM(210),GA(210),TAU(210),ICH(210),
     &                IBAR(210),K1(210),K2(210)
      COMMON/DINPDA/IMPS(6,6),IMVE(6,6),IB08(6,21),IB10(6,21),IA08(6,21)
     &,IA10(6,21),A1,B1,B2,B3,LT,LE,BET,AS,B8,AME,DIQ,ISU
      COMMON /TRAFOP/ GAMP,BGAMP,BETP
      COMMON /ENERIN/ EPROJ,ETARG
      COMMON /SDFLAG/ ISD
      COMMON/XXLMDD/IJLMDD,KDLMDD 
      COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
*
      DIMENSION MQUARK(3,30),IHKKQ(-6:6),IHKKQQ(-3:3,-3:3),
     &          IDX(-4:4)
      DATA IDX   /-4,-3,-1,-2,0,2,1,3,4/
      DATA IHKKQ /-6,-5,-4,-3,-1,-2,0,2,1,3,4,5,6/
      DATA IHKKQQ/-3301,-3103,-3203,0,   0,0,0,
     &            -3103,-1103,-2103,0,   0,   0,   0,
     &            -3203,-2103,-2203,0,   0,   0,   0,
     &                0,    0,    0,0,   0,   0,   0,
     &                0,    0,    0,0,2203,2103,3202,
     &                0,    0,    0,0,2103,1103,3103,
     &                0,    0,    0,0,3203,3103,3303/
*
*----------------------------------- quark content of hadrons:
*                                     1, 2, 3, 4  -   u, d, s, c
*                                    -1,-2,-3,-4  -  au,ad,as,ac
*
      DATA MQUARK/
     &   1,1,2,    -1,-1,-2,       0,0,0,       0,0,0,       0,0,0,
     &   0,0,0,       0,0,0,       1,2,2,    -1,-2,-2,       0,0,0,
     &   0,0,0,       0,0,0,      1,-2,0,      2,-1,0,      1,-3,0,
     &  3,-1,0,       1,2,3,    -1,-2,-3,       0,0,0,       2,2,3,
     &   1,1,3,       1,2,3,      1,-1,0,      2,-3,0,      3,-2,0,
     &  2,-2,0,      3,-3,0,       0,0,0,       0,0,0,       0,0,0/
      DATA UNON/2.0/
      DATA NCREJ, NCPT /0, 0/
*
      IJLMDD= 1
      ISD   = 2
      IREJ  = 0
      IIREJ = 0
      IBPROJ = IBAR(KPROJ)
      IBTARG = IBAR(KTARG)
      EPROJ  = (AM(KPROJ)**2-AM(KTARG)**2+ECM**2)/(2.0D0*ECM)
      ETARG  = (AM(KTARG)**2-AM(KPROJ)**2+ECM**2)/(2.0D0*ECM)
      IF(IPEV.GE.2) WRITE(LOUT,1014)EPROJ,ETARG
 1014 FORMAT('VALMSD: EPROJ,ETARG',2F10.5)
      IF (IBTARG.LE.0) THEN
         WRITE(LOUT,1001) IBTARG
 1001    FORMAT('VALMSD: NO HMSD FOR TARGET WITH BARYON-CHARGE',I4)
         IIREJ = 1
         GOTO 9999
      ENDIF
      IQP1 = MQUARK(1,KPROJ)
      IQP2 = MQUARK(2,KPROJ)
      IQP3 = MQUARK(3,KPROJ)
      IQT1 = MQUARK(1,KTARG)
      IQT2 = MQUARK(2,KTARG)
      IQT3 = MQUARK(3,KTARG)
      IF(IPEV.GE.2) WRITE(LOUT,1002)
     &                IBPROJ,IBTARG,IQP1,IQP2,IQP3,IQT1,IQT2,IQT3
 1002 FORMAT('VALMSD: IBPROJ,IBTARG,IQP1,IQP2,IQP3,IQT1,IQT2,IQT3 ',8I4)
*
      IF (IBPROJ.NE.0) THEN
*
*-------------------- q-qq (aq-aqaq) - flavors of projectile (baryon)
*
         ISAM = 1.0D0+2.999D0*RNDM(V)
         GOTO (10,11,12) ISAM
   10    CONTINUE
         IQP    = IQP1
         IDIQP1 = IQP2
         IDIQP2 = IQP3
         GOTO 13
   11    CONTINUE
         IQP    = IQP2
         IDIQP1 = IQP1
         IDIQP2 = IQP3
         GOTO 13
   12    CONTINUE
         IQP    = IQP3
         IDIQP1 = IQP1
         IDIQP2 = IQP2
   13    CONTINUE
*
      ELSE IF (IBPROJ.EQ.0) THEN
*
*-------------------- q-aq - flavors of projectile (meson)
*
         ISAM = 1.0D0+1.999D0*RNDM(V)
         GOTO (14,15) ISAM
   14    CONTINUE
         IQP    = IQP1
         IDIQP1 = IQP2
         GOTO 16
   15    CONTINUE
         IQP    = IQP2
         IDIQP1 = IQP1
   16    CONTINUE
         IDIQP2 = 0
*
      ENDIF
*
*-------------------- q-qq - flavors of target (baryon)
*
      ISAM = 1.0D0+2.999D0*RNDM(V)
      GOTO (17,18,19) ISAM
   17 CONTINUE
      IQT    = IQT1
      IDIQT1 = IQT2
      IDIQT2 = IQT3
      GOTO 20
   18 CONTINUE
      IQT    = IQT2
      IDIQT1 = IQT1
      IDIQT2 = IQT3
      GOTO 20
   19 CONTINUE
      IQT    = IQT3
      IDIQT1 = IQT1
      IDIQT2 = IQT2
   20 CONTINUE
*
      IKVQ1  = IQP
      IKD1Q1 = IDIQP1
      IKD2Q1 = IDIQP2
*
      IKVQ2  = IQT
      IKD1Q2 = IDIQT1
      IKD2Q2 = IDIQT2
*
      IF (IPEV.GE.2) WRITE(LOUT,1003)
     &                 IKVQ1,IKD1Q1,IKD2Q1,IKVQ2,IKD1Q2,IKD2Q2
 1003 FORMAT('VALMSD: IKVQ1,IKD1Q1,IKD2Q1,IKVQ2,IKD1Q2,IKD2Q2 ',6I4)
*
*-------------------- IDIFTP = 1  target     (backward) hadron  excited
*                     IDIFTP = 2  projectile (forward)  hadron  excited
*
      IDIFTP = 1.0D0+1.999D0*RNDM(V)
*-------------------- S.Roesler 5/26/93
      IF ((ISINGD.EQ.3).OR.(ISINGD.EQ.7)) IDIFTP = 1
      IF ((ISINGD.EQ.4).OR.(ISINGD.EQ.8)) IDIFTP = 2
*
      IF (IPEV.GE.2) WRITE(LOUT,1004) IDIFTP
 1004 FORMAT('VALMSD: IDIFTP ',I4)
      IF ((IDIFTP.NE.1).AND.(IDIFTP.NE.2)) THEN
         IF (IPEV.GE.2) WRITE(LOUT,'(A19)') 'VALMSD-ERROR: IDIFTP'
         GOTO 9999
      ENDIF
*
*-------------------- diffractive mass
*
31    CONTINUE
      AMO = 1.5D0
C------------------------------------------------------------
C
C                             first simple test j.r. 2/94
C
C-----------------------------------------------------------
C     AMO = 4.0D0
C-----------------------------------------------------------
      R   = RNDM(V)
      IF ((IBPROJ.EQ.0).AND.(IDIFTP.EQ.1))
     &               AMO = 1.5D0*R+2.83D0*(1.0D0-R)
      IF ((IBPROJ.EQ.0).AND.(IDIFTP.EQ.2)) AMO = 1.0D0
      SAM = 1.0D0
      IF (ECM.LE.300.0D0) SAM = 1.0D0-EXP(-((ECM/200.0D0)**4))
      R   = RNDM(V)*SAM
      AMAX= (1.0D0-SAM)*SQRT(0.1D0*ECM**2)+SAM*SQRT(400.0D0)
      AMU = R*SQRT(100.0D0)+(1.0D0-R)*AMAX
      IF (IBPROJ.EQ.0) THEN
C----------------- change mass-cuts
C        AMAX= (1.0D0-SAM)*SQRT(0.25D0*ECM**2)+SAM*SQRT(400.0D0)
         AMAX= (1.0D0-SAM)*SQRT(0.15D0*ECM**2)+SAM*SQRT(400.0D0)
         AMU = R*SQRT(100.0D0)+(1.0D0-R)*AMAX
      ENDIF
C---------------------------------------------------------------
C
C                                           j.r. test 2/94
C
C---------------------------------------------------------------
      AMU=2.D0*AMU
C---------------------------------------------------------------
      R   = RNDM(V)
      IF (ECM.LE.50.0D0) THEN
         AMDIFF = AMO*(AMU/AMO)**R
      ELSE
         A = 0.7D0
         IF (ECM.LE.300.0D0) A = 0.7D0*(1.0D0-EXP(-((ECM/100.0D0)**2)))
         AMDIFF = 1.0D0/((R/(AMU**A)+(1.0D0-R)/(AMO**A))
     &         **(1.0D0/A))
      ENDIF
      IF(AMDIFF.GT.0.5D0*ECM)GO TO 31
      IF (IOUXEV.GE.2) WRITE(LOUT,1005) AMDIFF
 1005 FORMAT('VALMSD: AMDIFF',E10.5)
*
*-------------------- kinematical parameters of
*                     diffractive projectile (IDIFTP = 1)
*                     or diffractive target  (IDIFTP = 2)
*
      IF (IDIFTP.EQ.1) THEN
         AMDCH3 = AM(KPROJ)
	 IF(KPROJ.EQ.1)THEN
	   KDLMDD = 61
           AMDCH3 = AM(61)
C                          19.4.99	   
	   KDLMDD = 54
           AMDCH3 = AM(54)
	 ELSEIF(KPROJ.EQ.2)THEN
C                          19.4.99	   
	   KDLMDD = 68
           AMDCH3 = AM(68)
	 ELSEIF(KPROJ.EQ.8)THEN
	   AMDCH3 = AM(62)
	   KDLMDD = 62
C                          19.4.99	   
	   AMDCH3 = AM(55)
	   KDLMDD = 55
	 ELSEIF(KPROJ.EQ.13)THEN
	   AMDCH3 = AM(186)
	   KDLMDD = 186
	 ELSEIF(KPROJ.EQ.14)THEN
	   AMDCH3 = AM(188)
	   KDLMDD = 188
	 ELSEIF(KPROJ.EQ.15)THEN
	   AMDCH3 = AM(190)
	   KDLMDD = 190
	 ELSEIF(KPROJ.EQ.16)THEN
	   AMDCH3 = AM(191)
	   KDLMDD = 191
	 ELSEIF(KPROJ.EQ.23)THEN
	   AMDCH3 = AM(187)
	   KDLMDD = 187
	 ELSEIF(KPROJ.EQ.24)THEN
	   AMDCH3 = AM(192)
	   KDLMDD = 192
	 ELSEIF(KPROJ.EQ.25)THEN
	   AMDCH3 = AM(193)
	   KDLMDD = 193
	 ENDIF
      ELSE IF (IDIFTP.EQ.2) THEN
         AMDCH3 = AM(KTARG)
	 IF(KTARG.EQ.1)THEN
	   KDLMDD = 61
           AMDCH3 = AM(61)
C                          19.4.99	   
	   KDLMDD = 54
           AMDCH3 = AM(54)
	 ELSEIF(KTARG.EQ.2)THEN
C                          19.4.99	   
	   KDLMDD = 68
           AMDCH3 = AM(68)
	 ELSEIF(KTARG.EQ.8)THEN
	   AMDCH3 = AM(62)
	   KDLMDD = 62
C                          19.4.99	   
	   AMDCH3 = AM(55)
	   KDLMDD = 55
	 ELSEIF(KTARG.EQ.13)THEN
	   AMDCH3 = AM(186)
	   KDLMDD = 186
	 ELSEIF(KTARG.EQ.14)THEN
	   AMDCH3 = AM(188)
	   KDLMDD = 188
	 ELSEIF(KTARG.EQ.15)THEN
	   AMDCH3 = AM(190)
	   KDLMDD = 190
	 ELSEIF(KTARG.EQ.16)THEN
	   AMDCH3 = AM(191)
	   KDLMDD = 191
	 ELSEIF(KTARG.EQ.23)THEN
	   AMDCH3 = AM(187)
	   KDLMDD = 187
	 ELSEIF(KTARG.EQ.24)THEN
	   AMDCH3 = AM(192)
	   KDLMDD = 192
	 ELSEIF(KTARG.EQ.25)THEN
	   AMDCH3 = AM(193)
	   KDLMDD = 193
	 ENDIF
      ENDIF
      ECH3 = (ECM**2+AMDCH3**2-AMDIFF**2)/(2.0D0*ECM)
      IF (ECH3.LE.AMDCH3) THEN
         NCMERR = NCMERR+1
         IF (MOD(NCMERR,2200).EQ.0)
     &      WRITE(LOUT,*) 'LMSD: INEFFICIENT SELECTION OF AMDIFF(1),
     &                        NCMERR = ',NCMERR
         GOTO 31
      ENDIF
*
      PCH3 = SQRT(ABS(ECH3-AMDCH3))*SQRT(ECH3+AMDCH3)
      IF (IDIFTP.EQ.2) PCH3 = -PCH3
      NDCH3  = -1
      GAMDC3 = ECH3/AMDCH3
      PGVC3  = PCH3/AMDCH3
      IF (IOUXEV.GE.2) WRITE(LOUT,1006) AMDCH3,ECH3,PCH3,GAMDC3,PGVC3
 1006 FORMAT('VALMSD: AMDCH3,ECH3,PCH3,GAMDC3,PGVC3',5E10.5)
*
   30 CONTINUE
      B33   = 8.0D0
      ES    = -2.0D0/(B33**2)*LOG(RNDM(V)*RNDM(V))
      HPS   = SQRT(ES*ES+2.0D0*ES*0.94D0)
      CALL DSFECF(SFE,CFE)
      PXCH3 = HPS*CFE
      PYCH3 = HPS*SFE
      PTCH3 = SQRT(PXCH3**2+PYCH3**2)
      IF (PTCH3.GT.ABS(PCH3)) THEN
         NCPT = NCPT+1
         IF (MOD(NCPT,500).EQ.0) WRITE(LOUT,1007) NCPT
 1007    FORMAT('VALMSD: INEFFICIENT PT-SELECTION 1, NCPT=',I8)
         GOTO 30
      ENDIF
      PLCH3 = SIGN(SQRT(ABS(PCH3-PTCH3))*SQRT(ABS(PCH3+PTCH3)),PCH3)
      IF (IOUXEV.GE.2) WRITE(LOUT,1008) ES,HPS,PXCH3,PYCH3,PLCH3
 1008 FORMAT('VALMSD: ES,HPS,PXCH3,PYCH3,PLCH3',5E10.5)
*
*-------------------- no chain 1
*
      NDCH1  = -99
      AMDCH1 = 0.0D0
      ECH1   = 0.0D0
      PCH1   = 0.0D0
      GAMDC1 = 1.0D0
      PGVC1  = 0.0D0
      PGXVC1 = 0.0D0
      PGYVC1 = 0.0D0
      PGZVC1 = 0.0D0
      DO 40 I=1,4
         PDD2(I) = 0.0D0
         PDQ1(I) = 0.0D0
   40 CONTINUE
*
*-------------------- kinematical parameters of
*                     excited target        (IDIFTP = 1)
*                     or excited projectile (IDIFTP = 2)
*
      AMDCH2 = AMDIFF
      ECH2   = (ECM**2+AMDCH2**2-AMDCH3**2)/(2.0D0*ECM)
      PCH2   = -SIGN(SQRT(ABS(ECH2-AMDCH2))*SQRT(ECH2+AMDCH2),PCH3)
      NDCH2  = 0
      GAMDC2 = ECH2/AMDCH2
      PGVC2  = PCH2/AMDCH2
*
*-------------------- (IDIFFP,IDIFAP in analogy to HMSD in VAHMSD)
*                     IDIFTP = 1 : IKD1Q2/IKD2Q2 - IDIFFP = IKVQ2
*                     IDIFTP = 2 : projectile - baryon
*                                  IKD1Q1/IKD2Q1 - IDIFFP = IKVQ1
*                                  projectile - antibaryon
*                                  IKD1Q1/IKD2Q1 - IDIFAP = IKVQ1
*                                  projectile - meson
*                                  IKD1Q1 > 0    - IDIFAP = IKVQ1
*                                  IKD1Q1 < 0    - IDIFFP = IKVQ1
*
      IF (IDIFTP.EQ.1) IDIFFP = IKVQ2
      IF (IDIFTP.EQ.2) THEN
         IF (IBPROJ.GT.0) IDIFFP = IKVQ1
         IF (IBPROJ.LT.0) IDIFAP = IKVQ1
         IF ((IBPROJ.EQ.0).AND.(IKD1Q1.GT.0)) IDIFAP = IKVQ1
         IF ((IBPROJ.EQ.0).AND.(IKD1Q1.LT.0)) IDIFFP = IKVQ1
      ENDIF
      IF (IOUXEV.GE.2) WRITE(LOUT,1009) AMDCH2,ECH2,PCH2,GAMDC2,PGVC2,
     &                                  IDIFFP,IDIFAP
 1009 FORMAT('VALMSD: AMDCH2,ECH2,PCH2,GAMDC2,PGVC2,IDIFFP,IDIFAP',
     &                                                  5E10.5,2I4)
      PXCH2 = -PXCH3
      PYCH2 = -PYCH3
      PTCH2 = SQRT(PXCH2**2+PYCH2**2)
      IF (PTCH2.GT.ABS(PCH2)) THEN
         NCPT = NCPT+1
         IF (MOD(NCPT,500).EQ.0) WRITE(LOUT,1010) NCPT
 1010    FORMAT('VALMSD: INEFFICIENT PT-SELECTION 2, NCPT=',I8)
         GOTO 30
      ENDIF
      PLCH2 = SIGN(SQRT(ABS(PCH2-PTCH2))*SQRT(ABS(PCH2+PTCH2)),PCH2)
      IF (IOUXEV.GE.2) WRITE(LOUT,1011) PXCH2,PYCH2,PLCH2
 1011 FORMAT('VALMSD: PXCH2,PYCH2,PLCH2',3E10.5)
      CC    = AMDCH2/(2.0D0*ECH2)
*--------------------------- S.R. 11/12/92
      IF (CC.GE.0.5D0) THEN
         NCMERR = NCMERR+1
         IF (MOD(NCMERR,200).EQ.0)
     &      WRITE(LOUT,*) 'LMSD: INEFFICIENT SELECTION OF AMDIFF(2),
     &                        NCMERR = ',NCMERR
         GOTO 31
      ENDIF
*
      XDIQ  = 0.5D0+SQRT(ABS((0.5D0-CC)*(0.5D0+CC)))
      XQ    = 1.0D0-XDIQ
      EDIQ  = ECH2*XDIQ
      EQ    = ECH2*XQ
      TEMP  = ABS(EDIQ/PCH2)
      PXDIQ = PXCH2*TEMP
      PYDIQ = PYCH2*TEMP
      PLDIQ = PLCH2*TEMP
      TEMP  = ABS(EQ/PCH2)
      PXQ   = -PXCH2*TEMP
      PYQ   = -PYCH2*TEMP
      PLQ   = -PLCH2*TEMP
*
*-------------------- store results in COMMON-block /ABRDIF/
*
      PDQ2(1)  = PXQ
      PDQ2(2)  = PYQ
      PDQ2(3)  = PLQ
      PDQ2(4)  = EQ
      PDD1(1)  = PXDIQ
      PDD1(2)  = PYDIQ
      PDD1(3)  = PLDIQ
      PDD1(4)  = EDIQ
      PDFQ1(1) = PXCH3
      PDFQ1(2) = PYCH3
      PDFQ1(3) = PLCH3
      PDFQ1(4) = ECH3
      PGXVC2   = PXCH2/AMDCH2
      PGYVC2   = PYCH2/AMDCH2
      PGZVC2   = PLCH2/AMDCH2
      PGXVC3   = PXCH3/AMDCH3
      PGYVC3   = PYCH3/AMDCH3
      PGZVC3   = PLCH3/AMDCH3
*
      IF (IPEV.GE.2) THEN
         WRITE(LOUT,*) 'VALMSD: PDQ1, PDD2, PDD1, PDQ2, PDFQ1'
         WRITE(LOUT,1012)
     &      (PDQ1(I),PDD2(I),PDD1(I),PDQ2(I),PDFQ1(I),I=1,4)
 1012    FORMAT(5F10.5)
         WRITE(LOUT,1013) 'PGXVC1,PGYVC1,PGZVC1,GAMDC1',
     &                     PGXVC1,PGYVC1,PGZVC1,GAMDC1
         WRITE(LOUT,1013) 'PGXVC2,PGYVC2,PGZVC2,GAMDC2',
     &                     PGXVC2,PGYVC2,PGZVC2,GAMDC2
         WRITE(LOUT,1013) 'PGXVC3,PGYVC3,PGZVC3,GAMDC3',
     &                     PGXVC3,PGYVC3,PGZVC3,GAMDC3
 1013    FORMAT(A27,4F10.5)
      ENDIF
*
*------------------- store results in common block /HKKEVT/
*
*                    partonen of projectile/target
*
      NHKK = NHKK + 1
      ISTHKK(NHKK) = 23-IDIFTP
      IF (IDIFTP.EQ.1) IDHKK(NHKK) = IHKKQ(IDX(IKVQ2))
      IF (IDIFTP.EQ.2) IDHKK(NHKK) = IHKKQ(IDX(IKVQ1))
      IF (IDIFTP.EQ.1) JMOHKK(1,NHKK) = ITAPOI
      IF (IDIFTP.EQ.2) JMOHKK(1,NHKK) = 1
      JMOHKK(2,NHKK) = 0
      JDAHKK(1,NHKK) = 0
      JDAHKK(2,NHKK) = 0
      PHKK(1,NHKK)   = 0.0D0
      PHKK(2,NHKK)   = 0.0D0
      PHKK(3,NHKK)   = XQ
      PHKK(4,NHKK)   = XQ
      PHKK(5,NHKK)   = 0.0D0
      VHKK  (1,NHKK) = VHKK(1,JMOHKK(1,NHKK))
      VHKK  (2,NHKK) = VHKK(2,JMOHKK(1,NHKK))
      VHKK  (3,NHKK) = VHKK(3,JMOHKK(1,NHKK))
      VHKK  (4,NHKK) = VHKK(4,JMOHKK(1,NHKK))
*
      NHKK = NHKK + 1
      ISTHKK(NHKK) = 23-IDIFTP
      IF (IDIFTP.EQ.1) IDHKK(NHKK) = IHKKQQ(IDX(IKD1Q2),IDX(IKD2Q2))
      IF (IDIFTP.EQ.2) IDHKK(NHKK) = IHKKQQ(IDX(IKD1Q1),IDX(IKD2Q1))
      IF (IDIFTP.EQ.1) JMOHKK(1,NHKK) = ITAPOI
      IF (IDIFTP.EQ.2) JMOHKK(1,NHKK) = 1
      JMOHKK(2,NHKK) = 0
      JDAHKK(1,NHKK) = 0
      JDAHKK(2,NHKK) = 0
      PHKK(1,NHKK)   = 0.0D0
      PHKK(2,NHKK)   = 0.0D0
      PHKK(3,NHKK)   = XDIQ
      PHKK(4,NHKK)   = XDIQ
      PHKK(5,NHKK)   = 0.0D0
      VHKK  (1,NHKK) = VHKK(1,JMOHKK(1,NHKK))
      VHKK  (2,NHKK) = VHKK(2,JMOHKK(1,NHKK))
      VHKK  (3,NHKK) = VHKK(3,JMOHKK(1,NHKK))
      VHKK  (4,NHKK) = VHKK(4,JMOHKK(1,NHKK))
*
*                    ends of chain 2
*
      NHKK = NHKK + 1
      ISTHKK(NHKK)   = 100+ISTHKK(NHKK-2)
      IDHKK(NHKK)    = IDHKK(NHKK-2)
      JMOHKK(1,NHKK) = NHKK-2
      JMOHKK(2,NHKK) = JMOHKK(1,NHKK-2)
      JDAHKK(1,NHKK) = NHKK+2
      JDAHKK(2,NHKK) = NHKK+2
      PHKK(1,NHKK)   = PDQ2(1)
      PHKK(2,NHKK)   = PDQ2(2)
      PHKK(3,NHKK)   = PDQ2(3)
      PHKK(4,NHKK)   = PDQ2(4)
      PHKK(5,NHKK)   = 0.0D0
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
      VHKK  (1,NHKK) = VHKK(1,JMOHKK(1,NHKK))+XXPP
      VHKK  (2,NHKK) = VHKK(2,JMOHKK(1,NHKK))+YYPP
      VHKK  (3,NHKK) = VHKK(3,JMOHKK(1,NHKK))
      VHKK  (4,NHKK) = VHKK(4,JMOHKK(1,NHKK))
*
      NHKK = NHKK + 1
      ISTHKK(NHKK)   = 100+ISTHKK(NHKK-2)
      IDHKK(NHKK)    = IDHKK(NHKK-2)
      JMOHKK(1,NHKK) = NHKK-2
      JMOHKK(2,NHKK) = JMOHKK(1,NHKK-2)
      JDAHKK(1,NHKK) = NHKK+1
      JDAHKK(2,NHKK) = NHKK+1
      PHKK(1,NHKK)   = PDD1(1)
      PHKK(2,NHKK)   = PDD1(2)
      PHKK(3,NHKK)   = PDD1(3)
      PHKK(4,NHKK)   = PDD1(4)
      PHKK(5,NHKK)   = 0.0D0
C               Add position of parton in hadron
       CALL QINNUC(XXPP,YYPP)
      VHKK  (1,NHKK) = VHKK(1,JMOHKK(1,NHKK))+XXPP
      VHKK  (2,NHKK) = VHKK(2,JMOHKK(1,NHKK))+YYPP
      VHKK  (3,NHKK) = VHKK(3,JMOHKK(1,NHKK))
      VHKK  (4,NHKK) = VHKK(4,JMOHKK(1,NHKK))
*
*                    chain 2
*
      NHKK = NHKK + 1
      ISTHKK(NHKK)   = 188
      IDHKK(NHKK)    = 88888
      JMOHKK(1,NHKK) = NHKK-2
      JMOHKK(2,NHKK) = NHKK-1
      JDAHKK(1,NHKK) = 0
      JDAHKK(2,NHKK) = 0
      PHKK(1,NHKK)   = PXCH2
      PHKK(2,NHKK)   = PYCH2
      PHKK(3,NHKK)   = PLCH2
      PHKK(4,NHKK)   = ECH2
      PHKK(5,NHKK)   = AMDCH2
      VHKK  (1,NHKK) = VHKK(1,NHKK-1)
      VHKK  (2,NHKK) = VHKK(2,NHKK-1)
      VHKK  (3,NHKK) = VHKK(3,NHKK-1)
      IF ((BETP.NE.0.0D0).AND.(BGAMP.NE.0.0D0))
     &VHKK  (4,NHKK) = VHKK(3,NHKK)/BETP-VHKK(3,NHKK-2)/BGAMP
*
*                    diffractive nucleon
*
      NHKK = NHKK + 1
      ISTHKK(NHKK)   = 188
      IDHKK(NHKK)    = 88888
      IF (IDIFTP.EQ.2) JMOHKK(1,NHKK) = ITAPOI
      IF (IDIFTP.EQ.1) JMOHKK(1,NHKK) = 1
      JMOHKK(2,NHKK) = 0
      JDAHKK(1,NHKK) = 0
      JDAHKK(2,NHKK) = 0
      PHKK(1,NHKK)   = PDFQ1(1)
      PHKK(2,NHKK)   = PDFQ1(2)
      PHKK(3,NHKK)   = PDFQ1(3)
      PHKK(4,NHKK)   = PDFQ1(4)
      PHKK(5,NHKK)   = AMDCH3
      VHKK  (1,NHKK) = VHKK(1,JMOHKK(1,NHKK))
      VHKK  (2,NHKK) = VHKK(2,JMOHKK(1,NHKK))
      VHKK  (3,NHKK) = VHKK(3,JMOHKK(1,NHKK))
      VHKK  (4,NHKK) = VHKK(4,JMOHKK(1,NHKK))
*
*
*---------- S. Roesler 21-10-93
*      check energy-momentum conservation
      IF (IPEV.GE.2) THEN
         PX = PDQ1(1)+PDQ2(1)+PDD1(1)+PDD2(1)+PDFQ1(1)
         PY = PDQ1(2)+PDQ2(2)+PDD1(2)+PDD2(2)+PDFQ1(2)
         PZ = PDQ1(3)+PDQ2(3)+PDD1(3)+PDD2(3)+PDFQ1(3)
         EE = ECM-(PDQ1(4)+PDQ2(4)+PDD1(4)+PDD2(4)+PDFQ1(4))
         WRITE(LOUT,*)'VALMSD: ENERGY-MOMENTUM-CHECK (PX,PY,PZ,E)'
         WRITE(LOUT,'(5F12.6)')PX,PY,PZ,EE,ECM
      ENDIF
*
      RETURN
 9999 CONTINUE
      NCREJ = NCREJ+1
      IF (MOD(NCREJ,500).EQ.0) WRITE(LOUT,9900) NCREJ
 9900 FORMAT('REJECTION IN VALMSD ',I10)
      IIREJ = 0
      IREJ  = 1
      RETURN
      END
*
*===hadrdi===============================================================*
*
      SUBROUTINE HADRDI(NAUX,IJPROJ,IJTAR,NHKKH1)
 
**************************************************************************
* Version November 1993                         by Stefan Roesler        *
*                                                  University of Leipzig *
* This subroutine prepares the hadronisation of diffractive chains       *
* (single-diffractive component), calls HADJET and stores the results    *
* in common block /HKKEVT/.                                              *
**************************************************************************
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)
      PARAMETER (NMXHKK= 89998,INTMX=2488,NAUMAX=897,NFIMAX=249)
      CHARACTER*8 ANAME,ANC,ANF
      COMMON /HKKEVT/ NHKK,NEVHKK,     ISTHKK(NMXHKK),  IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),PHKK(5,NMXHKK),
     &                VHKK(4,NMXHKK),  WHKK(4,NMXHKK)
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
      COMMON /DPAR/   ANAME(210),AAM(210),GA(210),TAU(210),IICH(210),
     &                IIBAR(210),K1(210),K2(210)
      COMMON /DIFFRA/ ISINGD,IDIFTP,IOUDIF,IFLAGD
      COMMON /DIFPAR/ PXC(902),  PYC(902),PZC(902),
     &                HEC(902),  AMC(902),ICHC(902),
     &                IBARC(902),ANC(902),NRC(902)
      COMMON /DFINPA/ ANF(NFIMAX),PXF(NFIMAX),PYF(NFIMAX),
     &                PZF(NFIMAX),HEF(NFIMAX),AMF(NFIMAX),
     &                ICHF(NFIMAX),IBARF(NFIMAX),NREF(NFIMAX)
      COMMON /DFINPZ/IORMO(NFIMAX),IDAUG1(NFIMAX),IDAUG2(NFIMAX),
     * ISTATH(NFIMAX)
      COMMON /ABRDIF/ XDQ1,XDQ2,XDDQ1,XDDQ2,
     &                IKVQ1,IKVQ2,IKD1Q1,IKD2Q1,IKD1Q2,IKD2Q2,
     &                IDIFFP,IDIFAP,
     &                AMDCH1,AMDCH2,AMDCH3,GAMDC1,GAMDC2,GAMDC3,
     &                PGXVC1,PGYVC1,PGZVC1,PGXVC2,PGYVC2,PGZVC2,
     &                PGXVC3,PGYVC3,PGZVC3,NDCH1,NDCH2,NDCH3,
     &                IKDCH1,IKDCH2,IKDCH3,
     &                PDQ1(4),PDQ2(4),PDD1(4),PDD2(4),PDFQ1(4)
      COMMON /DINPDA/ IMPS(6,6),IMVE(6,6),IB08(6,21),IB10(6,21),
     &                IA08(6,21),IA10(6,21),
     &                A1,B1,B2,B3,LT,LE,BET,AS,B8,AME,DIQ,ISU
      COMMON /ENERIN/ EPROJ,ETARG
      COMMON /DIFOUT/ AMCHDI,TTT,NNAUX,KPROJ,KTARG
      COMMON /SDFLAG/ ISD
      COMMON/XXLMDD/IJLMDD,KDLMDD 
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
C---------------------
*
      DIMENSION POJ1(4),PAT1(4),POJ2(4),PAT2(4)
*
      NHKKH1 = NHKK
      NAUX   = 0
      NHAD   = 0
      PXDI   = 0.0D0
      PYDI   = 0.0D0
      PLDI   = 0.0D0
      EDI    = 0.0D0
      IMOH2  = NHKK-1
      IMOH3  = NHKK
      IF (ISD.EQ.1) IMOH1 = NHKK-4
      IBPROJ = IIBAR(IJPROJ)
      IBTAR  = IIBAR(IJTAR)
      KPROJ  = IJPROJ
      KTARG  = IJTAR
*
      IF (IPEV.GE.2) WRITE(LOUT,1100) ISD
 1100 FORMAT('HADRDI: ISD ',I3)
*
      NUNUC1 = 3
      NUNUC2 = 3
      IF ((IBPROJ.NE.0).AND.(IDIFTP.EQ.2)) NUNUC2 = 6
      IF ((IBTAR .NE.0).AND.(IDIFTP.EQ.1)) NUNUC2 = 4
      IF (IDIFTP.EQ.2) THEN
         DO 10 J=1,4
            POJ1(J) = PDQ1(J)
            PAT1(J) = PDD2(J)
            POJ2(J) = PDD1(J)
            PAT2(J) = PDQ2(J)
   10    CONTINUE
      ELSE
         DO 11 J=1,4
            POJ1(J) = PDD2(J)
            PAT1(J) = PDQ1(J)
            POJ2(J) = PDQ2(J)
            PAT2(J) = PDD1(J)
   11    CONTINUE
      ENDIF
      IF (IDIFTP.EQ.1) THEN
         IFB11 = IDIFAP
         IFB12 = IKVQ2
         IFB21 = IDIFFP
         IFB22 = IKD1Q2
         IFB23 = IKD2Q2
      ELSE IF ((IDIFTP.EQ.2).AND.(IKVQ1.GT.0)) THEN
         IFB11 = IKVQ1
         IFB12 = IDIFAP
         IFB21 = IKD1Q1
         IFB22 = IKD2Q1
         IFB23 = IDIFFP
      ELSE IF ((IDIFTP.EQ.2).AND.(IKVQ1.LE.0)) THEN
         IFB11 = IKVQ1
         IFB12 = IDIFFP
         IFB21 = IKD1Q1
         IFB22 = IKD2Q1
         IFB23 = IDIFAP
      ENDIF
      IF ((IFB22.EQ.0).AND.(IFB23.NE.0)) THEN
         IFB22 = IFB23
         IFB23 = 0
      ENDIF
      IF (IFB11.LT.0) IFB11 = IABS(IFB11)+6
      IF (IFB12.LT.0) IFB12 = IABS(IFB12)+6
      IF (IFB21.LT.0) IFB21 = IABS(IFB21)+6
      IF (IFB22.LT.0) IFB22 = IABS(IFB22)+6
      IF (IFB23.LT.0) IFB23 = IABS(IFB23)+6
*
*-------------------------- chain 1
*
C        WRITE(LOUT,1001) (POJ1(I),I=1,4)
C        WRITE(LOUT,1002) (PAT1(I),I=1,4)
      IF (IPEV.GE.2) THEN
         WRITE(LOUT,1000) NHAD,IKDCH1,NUNUC1,NDCH1,IFB11,IFB12,
     &                                             IFB13,IFB14
         WRITE(LOUT,1001) (POJ1(I),I=1,4)
         WRITE(LOUT,1002) (PAT1(I),I=1,4)
         WRITE(LOUT,1003) GAMDC1,PGXVC1,PGYVC1,PGZVC1,AMDCH1
 1000    FORMAT('HADRDI: NHAD,IKDCH1,NUNUC1,NDCH1,IFB11,IFB12',
     &                  ', IFB13,IFB14 ',8I4)
 1001    FORMAT('HADRDI: POJ1 ',4E15.5)
 1002    FORMAT('HADRDI: PAT1 ',4E15.5)
 1003    FORMAT('HADRDI: GAMDC1,PGXVC1,PGYVC1,PGZVC1,AMDCH1',5F10.5)
      ENDIF
      CALL HADJET(NHAD,AMDCH1,POJ1,PAT1,GAMDC1,PGXVC1,
     &            PGYVC1,PGZVC1,IFB11,IFB12,IFB13,IFB14,
     &            IKDCH1,IKDCH1,NUNUC1,NDCH1,7)
*
      NHKKBE = NHKK+1
      IF (NHAD.EQ.0) GOTO 13
      DO 12 J=1,NHAD
         IF (NHKK.EQ.NMXHKK) THEN
            WRITE(LOUT,1004) NHKK
 1004       FORMAT('HADRDI: NHKK.EQ.NMXHKK',I10)
            GOTO 9999
         ENDIF
         NHKK   = NHKK+1
         NAUX   = NAUX+1
         IF(IBARF(J).EQ.500)GO TO 776
         ECHECK = SQRT(PXF(J)**2+PYF(J)**2+PZF(J)**2+AMF(J)**2)
         IF (ABS(ECHECK-HEF(J)).GT.0.5D-2) THEN
            WRITE(LOUT,1005) NHKK,ECHECK,HEF(J),AMF(J)
 1005       FORMAT('HADRDI: CHAIN 1  CORRECT INCONSISTENT ENERGY',
     &              I10,3F10.5)
            HEF(J) = ECHECK
         ENDIF
  776    CONTINUE
         ANNDI=ANNDI+1
         EEDI=EEDI+HEF(J)
         PTDI=PTDI+SQRT(PXF(J)**2+PYF(J)**2)
         PXC  (NAUX) = PXF  (J)
         PYC  (NAUX) = PYF  (J)
         PZC  (NAUX) = PZF  (J)
         HEC  (NAUX) = HEF  (J)
         AMC  (NAUX) = AMF  (J)
         ICHC (NAUX) = ICHF (J)
         IBARC(NAUX) = IBARF(J)
         ANC  (NAUX) = ANF  (J)
         NRC  (NAUX) = NREF (J)
         PXDI   = PXDI+PXF(J)
         PYDI   = PYDI+PYF(J)
         PLDI   = PLDI+PZF(J)
         EDI    = EDI +HEF(J)
         EDIFF  = EDIFF+HEF(J)
         PTDIFF = PTDIFF+SQRT(PXF(J)**2+PYF(J)**2)
         ISTHKK(NHKK)   = 1
         IF(IBARF(J).EQ.500)ISTHKK(NHKK)=2
         IDHKK (NHKK)   = MPDGHA(NREF(J))
	 IF(IORMO(J).EQ.999)THEN
           JMOHKK(1,NHKK) = IMOH1
	 ELSE
	   JMOHKK(1,NHKK)=NHKKBE+IORMO(J)-1
	 ENDIF
         JMOHKK(2,NHKK) = 0
         JDAHKK(1,NHKK) = 0
         JDAHKK(2,NHKK) = 0
         PHKK(1,NHKK)   = PXF(J)
         PHKK(2,NHKK)   = PYF(J)
         PHKK(3,NHKK)   = PZF(J)
         PHKK(4,NHKK)   = HEF(J)
         PHKK(5,NHKK)   = AMF(J)
         IMOHKK         = JMOHKK(1,NHKK)
         VHKK(1,NHKK)   = VHKK(1,IMOHKK)
         VHKK(2,NHKK)   = VHKK(2,IMOHKK)
         VHKK(3,NHKK)   = VHKK(3,IMOHKK)
         VHKK(4,NHKK)   = VHKK(4,IMOHKK)
         IF (IPHKK.GE.1) THEN
            WRITE(LOUT,1006) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     &                       JMOHKK(1,NHKK),JMOHKK(2,NHKK),
     &                       JDAHKK(1,NHKK),JDAHKK(2,NHKK)
            WRITE(LOUT,1007) (PHKK(K,NHKK),K=1,5)
 1006       FORMAT('HADRDI: NHKK,ISTHKK,IDHKK,JMOHKK,JDAHKK',7I6)
 1007       FORMAT('HADRDI: PHKK ',5F10.5)
         ENDIF
   12 CONTINUE
   13 CONTINUE
      IF ((NHAD.GT.0).AND.(ISD.EQ.1)) THEN
         JDAHKK(1,IMOH1) = NHKKBE
         JDAHKK(2,IMOH1) = NHKK
      ENDIF
*
*------------------------ chain 2
*
      NHAD = 0
C        WRITE(LOUT,1009) (POJ2(I),I=1,4)
C        WRITE(LOUT,1010) (PAT2(I),I=1,4)
      IF (IPEV.GE.2) THEN
         WRITE(LOUT,1008) NHAD,IKDCH2,NUNUC2,NDCH2,IFB21,IFB22,
     &                                             IFB23,IFB24
         WRITE(LOUT,1009) (POJ2(I),I=1,4)
         WRITE(LOUT,1010) (PAT2(I),I=1,4)
         WRITE(LOUT,1011) GAMDC2,PGXVC2,PGYVC2,PGZVC2,AMDCH2
 1008    FORMAT('HADRDI: NHAD,IKDCH2,NUNUC2,NDCH2,IFB21,IFB22,
     &                   IFB23,IFB24 ',8I4)
 1009    FORMAT('HADRDI: POJ2 ',4E15.5)
 1010    FORMAT('HADRDI: PAT2 ',4E15.5)
 1011    FORMAT('HADRDI: GAMDC2,PGXVC2,PGYVC2,PGZVC2,AMDCH2',5F10.5)
      ENDIF
      CALL HADJET(NHAD,AMDCH2,POJ2,PAT2,GAMDC2,PGXVC2,
     &            PGYVC2,PGZVC2,IFB21,IFB22,IFB23,IFB24,
     &            IKDCH2,IKDCH2,NUNUC2,NDCH2,8)
*
      NHKKBE = NHKK+1
      IF (NHAD.EQ.0) GOTO 15
      DO 14 J=1,NHAD
         IF (NHKK.EQ.NMXHKK) THEN
            WRITE(LOUT,1012) NHKK
 1012       FORMAT('HADRDI: NHKK.EQ.NMXHKK',I10)
            GOTO 9999
         ENDIF
         NHKK   = NHKK+1
         NAUX   = NAUX+1
         IF(IBARF(J).EQ.500)GO TO 775
         ECHECK = SQRT(PXF(J)**2+PYF(J)**2+PZF(J)**2+AMF(J)**2)
         IF (ABS(ECHECK-HEF(J)).GT.0.5D-2) THEN
            WRITE(LOUT,1013) NHKK,ECHECK,HEF(J),AMF(J)
 1013       FORMAT('HADRDI: CHAIN 2  CORRECT INCONSISTENT ENERGY',
     &              I10,3F10.5)
            HEF(J) = ECHECK
         ENDIF
  775    CONTINUE
         ANNDI=ANNDI+1
         EEDI=EEDI+HEF(J)
         PTDI=PTDI+SQRT(PXF(J)**2+PYF(J)**2)
         PXC  (NAUX) = PXF  (J)
         PYC  (NAUX) = PYF  (J)
         PZC  (NAUX) = PZF  (J)
         HEC  (NAUX) = HEF  (J)
         AMC  (NAUX) = AMF  (J)
         ICHC (NAUX) = ICHF (J)
         IBARC(NAUX) = IBARF(J)
         ANC  (NAUX) = ANF  (J)
         NRC  (NAUX) = NREF (J)
         PXDI   = PXDI+PXF(J)
         PYDI   = PYDI+PYF(J)
         PLDI   = PLDI+PZF(J)
         EDI    = EDI +HEF(J)
         EDIFF  = EDIFF+HEF(J)
         PTDIFF = PTDIFF+SQRT(PXF(J)**2+PYF(J)**2)
         ISTHKK(NHKK)   = 1
         IF(IBARF(J).EQ.500)ISTHKK(NHKK)=2
         IDHKK (NHKK)   = MPDGHA(NREF(J))
	 IF(IORMO(J).EQ.999)THEN
           JMOHKK(1,NHKK) = IMOH2
	 ELSE
	   JMOHKK(1,NHKK)=NHKKBE+IORMO(J)-1
	 ENDIF
         JMOHKK(2,NHKK) = 0
         JDAHKK(1,NHKK) = 0
         JDAHKK(2,NHKK) = 0
         PHKK(1,NHKK)   = PXF(J)
         PHKK(2,NHKK)   = PYF(J)
         PHKK(3,NHKK)   = PZF(J)
         PHKK(4,NHKK)   = HEF(J)
         PHKK(5,NHKK)   = AMF(J)
         IMOHKK         = JMOHKK(1,NHKK)
         VHKK(1,NHKK)   = VHKK(1,IMOHKK)
         VHKK(2,NHKK)   = VHKK(2,IMOHKK)
         VHKK(3,NHKK)   = VHKK(3,IMOHKK)
         VHKK(4,NHKK)   = VHKK(4,IMOHKK)
         IF (IPHKK.GE.2) THEN
            WRITE(LOUT,1014) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     &                       JMOHKK(1,NHKK),JMOHKK(2,NHKK),
     &                       JDAHKK(1,NHKK),JDAHKK(2,NHKK)
            WRITE(LOUT,1015) (PHKK(K,NHKK),K=1,5)
 1014       FORMAT('HADRDI: NHKK,ISTHKK,IDHKK,JMOHKK,JDAHKK',7I6)
 1015       FORMAT('HADRDI: PHKK ',5F10.5)
         ENDIF
   14 CONTINUE
   15 CONTINUE
      IF (NHAD.GT.0) THEN
         JDAHKK(1,IMOH2) = NHKKBE
         JDAHKK(2,IMOH2) = NHKK
      ENDIF
*
*----------------------------- diffractive nucleon/antinucleon
*
      NAUX = NAUX+1
      IF (IDIFTP.EQ.1) THEN
	IF(IJLMDD.EQ.1)THEN
         AMC  (NAUX) = AAM  (KDLMDD)
         ICHC (NAUX) = IICH (KDLMDD)
         IBARC(NAUX) = IIBAR(KDLMDD)
         ANC  (NAUX) = ANAME(KDLMDD)
         NRC  (NAUX) =KDLMDD 
	ELSEIF(IJLMDD.EQ.0)THEN
         AMC  (NAUX) = AAM  (IJPROJ)
         ICHC (NAUX) = IICH (IJPROJ)
         IBARC(NAUX) = IIBAR(IJPROJ)
         ANC  (NAUX) = ANAME(IJPROJ)
         NRC  (NAUX) = IJPROJ
        ENDIF
      ELSEIF (IDIFTP.EQ.2) THEN
	IF(IJLMDD.EQ.1)THEN
         AMC  (NAUX) = AAM  (KDLMDD)
         ICHC (NAUX) = IICH (KDLMDD)
         IBARC(NAUX) = IIBAR(KDLMDD)
         ANC  (NAUX) = ANAME(KDLMDD)
         NRC  (NAUX) =KDLMDD 
	ELSEIF(IJLMDD.EQ.0)THEN
         AMC  (NAUX) = AAM  (IJTAR)
         ICHC (NAUX) = IICH (IJTAR)
         IBARC(NAUX) = IIBAR(IJTAR)
         ANC  (NAUX) = ANAME(IJTAR)
         NRC  (NAUX) = IJTAR
        ENDIF
      ENDIF
      PXC(NAUX) = PDFQ1(1)
      PYC(NAUX) = PDFQ1(2)
      PZC(NAUX) = PDFQ1(3)
      HEC(NAUX) = PDFQ1(4)
         ANNDI=ANNDI+1
         EEDI=EEDI+HEC(NAUX)
         PTDI=PTDI+SQRT(PXC(NAUX)**2+PYC(NAUX)**2)
      NHKK = NHKK+1
      ISTHKK(NHKK)   = 1
C        IF(IBARF(J).EQ.500)ISTHKK(NHKK)=2
      IDHKK (NHKK)   = MPDGHA(NRC(NAUX))
	 IF(IORMO(J).EQ.999)THEN
           JMOHKK(1,NHKK) = IMOH3
	 ELSE
	   JMOHKK(1,NHKK)=NHKKBE+IORMO(J)-1
	 ENDIF
      JMOHKK(2,NHKK) = 0
      JDAHKK(1,NHKK) = 0
*---------- S. Roesler 5-11-93
*           the following line is just to get the leading particle
*           detected
      JDAHKK(2,NHKK) = 1000
*
      PHKK(1,NHKK)   = PXC(NAUX)
      PHKK(2,NHKK)   = PYC(NAUX)
      PHKK(3,NHKK)   = PZC(NAUX)
      PHKK(4,NHKK)   = HEC(NAUX)
      PHKK(5,NHKK)   = AMC(NAUX)
      JDAHKK(1,IMOH3)= NHKK
      IMOHKK         = JMOHKK(1,NHKK)
      VHKK(1,NHKK)   = VHKK(1,IMOHKK)
      VHKK(2,NHKK)   = VHKK(2,IMOHKK)
      VHKK(3,NHKK)   = VHKK(3,IMOHKK)
      VHKK(4,NHKK)   = VHKK(4,IMOHKK)
      IF (IPHKK.GE.2) THEN
         WRITE(LOUT,1016) NHKK,ISTHKK(NHKK),IDHKK(NHKK),
     &                    JMOHKK(1,NHKK),JMOHKK(2,NHKK),
     &                    JDAHKK(1,NHKK),JDAHKK(2,NHKK)
         WRITE(LOUT,1017) (PHKK(K,NHKK),K=1,5)
 1016    FORMAT('HADRDI: NHKK,ISTHKK,IDHKK,JMOHKK,JDAHKK',7I6)
 1017    FORMAT('HADRDI: PHKK ',5F10.5)
      ENDIF
      IF (IPEV.GE.2) THEN
         DO 16 I=1,NAUX
            WRITE(LOUT,*) 'HADRDI: DIFFRACTIVE JET-HADRONS'
            WRITE(LOUT,1018)I,PXC(I),PYC(I),PZC(I),HEC(I),AMC(I),
     &                      ICHC(I),IBARC(I),NRC(I),ANC(I)
 1018       FORMAT(I5,5F12.4,3I5,A10)
   16    CONTINUE
      ENDIF
*
      NNAUX  = NAUX
      AMCHDI = SQRT(EDI**2-PXDI**2-PYDI**2-PLDI**2)
      IF (IDIFTP.EQ.1) TTT=2*AAM(IJPROJ)**2-2.*(EPROJ*HEC(NAUX)
     &                    -SQRT(EPROJ**2-AAM(IJPROJ)**2)*PZC(NAUX))
      IF (IDIFTP.EQ.2) TTT=2*AAM(IJTAR)**2-2.*(ETARG*HEC(NAUX)
     &                    +SQRT(ETARG**2-AAM(IJTAR)**2)*PZC(NAUX))
      IF (IPEV.GE.2) THEN
         WRITE(LOUT,1019) AMCHDI,TTT
 1019    FORMAT('HADRDI: AMCHDI,TTT ',2F10.5)
      ENDIF
 9999 CONTINUE
      RETURN
      END
*
*===diadif===============================================================*
*
      SUBROUTINE DIADIF(IOP,NHKKH1)
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)
      PARAMETER (NMXHKK= 89998)
      CHARACTER*8 ANAME,ANC
      COMMON /HKKEVT/ NHKK,NEVHKK,     ISTHKK(NMXHKK),  IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),PHKK(5,NMXHKK),
     &                VHKK(4,NMXHKK),  WHKK(4,NMXHKK)
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
      COMMON /DPAR/   ANAME(210),AAM(210),GA(210),TAU(210),IICH(210),
     &                IIBAR(210),K1(210),K2(210)
      COMMON /DIFPAR/ PXC(902),  PYC(902),PZC(902),
     &                HEC(902),  AMC(902),ICHC(902),
     &                IBARC(902),ANC(902),NRC(902)
      COMMON /ENERIN/ EPROJ,ETARG
      COMMON /NNCMS/  DGAMCM,DBGCM,ECM,DPCM,DEPROJ,DPPROJ
      COMMON /DHISTO/ XYL(51,10),YYL(51,10),YYLPS(51,10),XXFL(51,10),
     &                YXFL(51,10),TDTDM(40,24),DSDTDM(40,24),
     &                AMDM(24,40),DSDM(24,40),AVE(30),AVMULT(30),
     &                YXFLCH(51),YXFLPI(51)
      COMMON /DIFOUT/ AMCH,TT,NHAD,KPROJ,KTARG
      COMMON /DIFFRA/ ISINGD,IDIFTP,IOUDIF,IFLAGD
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
      COMMON /EVFLAG/ NUMEV
      DIMENSION INDX(28),PX(902),PY(902),PZ(902),HE(902),AM(902),
     &          ICH(902),NR(902)
C     DATA INDX /1,8,10,10,10,10,7,2,7,10,10,7,3,4,5,6,
C    &          11,12,7,13,14,15,16,17,18,7,7,7/
      DATA INDX/ 1,8,10,10,10,   10,7,2,7,10,
     &           10,7,3,4,5,      6,7,7,7,7,
     &            7,7,7,7,7,      7,7,7/
*
      GOTO (1,2,3) IOP
*
    1 CONTINUE
      NCEV = 0
      DXFL = 0.04D0
      DT   = 0.1D0
      DDM  = 10.0D0
      DY   = 0.499999D0
      NEVT = 0
      NCH  = 0
      NPI  = 0
      NHAD = 0
C     KPROJ= 0
C     KTARG= 0
      AMCH = 0.0D0
      TT   = 0.0D0
*
      DO 10 J=1,51
         YXFLCH(J) = 0.0D0
         YXFLPI(J) = 0.0D0
         DO 11 I=1,10
            XXFL(J,I) = J*DXFL-1.0D0
            YXFL(J,I) = 1.0D-18
            XYL(J,I)  = (J-24)*DY-DY/2.0D0
            YYL(J,I)  = 1.0D-18
            YYLPS(J,I)= 1.0D-18
   11    CONTINUE
   10 CONTINUE
      DO 12 I=1,40
         DO 13 J=1,24
            TDTDM(I,J) = I*DT
            DSDTDM(I,J)= 1.0D-8
            AMDM(J,I)  = (J*DDM-DDM/2.0D0)**2
            DSDM(J,I)  = 1.0D-8
            AMDM(J,I)  = LOG10(AMDM(J,I))
   13    CONTINUE
   12 CONTINUE
      DO 14 I=1,30
         AVE(I)    = 1.0D-18
         AVMULT(I) = 1.0D-18
   14 CONTINUE
      RETURN
*
    2 CONTINUE
      NCEV = NCEV+1
      IF (MOD(NCEV,2000).EQ.0) WRITE(LOUT,*) NCEV
      J = 0
C     IF (IFLAGD.EQ.1) RETURN
      DO 21 I=NHKKH1+1,NHKK
         IF ((ISTHKK(I).EQ.1).AND.(JMOHKK(2,I).NE.100).AND.
     &       (JDAHKK(1,I).EQ.0)) THEN
            J = J+1
            PX(J) = PHKK(1,I)
            PY(J) = PHKK(2,I)
            PZ(J) = PHKK(3,I)
            HE(J) = PHKK(4,I)
            AM(J) = PHKK(5,I)
            NR(J) = MCIHAD(IDHKK(I))
            ICH(J)= IICH(NR(J))
         ENDIF
   21 CONTINUE
      IHAD = J
      IF (IPEV.GE.2) THEN
         WRITE(LOUT,*) 'DIADIF: PX,PY,PZ,HE,AM,NR,ICH'
         DO 22 I=1,IHAD
            WRITE(LOUT,*) PX(I),PY(I),PZ(I),HE(I),AM(I),NR(I),ICH(I)
 1000       FORMAT(5F12.5,2I4)
   22    CONTINUE
      ENDIF
C     IF (IDIFTP.EQ.1) THEN
C        P0 = SQRT(ETARG**2-AAM(KTARG)**2)
C     ELSE IF (IDIFTP.EQ.2) THEN
C        P0 = SQRT(EPROJ**2-AAM(KPROJ)**2)
C     ENDIF
      EPCM = (AAM(IJPROJ)**2-AAM(1)**2+ECM**2)/(2.0D0*ECM)
      P0   = SQRT((EPCM-AAM(IJPROJ))*(EPCM+AAM(IJPROJ)))
      NEVT = NEVT+1
      AVMULT(30) = AVMULT(30)+IHAD
      DO 20 I=1,IHAD
         NRE = NR(I)
         IF (NRE.GT.25) NRE = 28
         IF (NRE.LT. 1) NRE = 28
         NI = INDX(NRE)
         IF (NRE.EQ.28) NI = 8
         AVE(NRE) = AVE(NRE)+HE(I)
         AVE(30)  = AVE(30) +HE(I)
         IF (NI.NE.6) AVE(29) = AVE(29)+HE(I)
         AVMULT(NRE) = AVMULT(NRE)+1.0D0
         IF (NI.NE.6) AVMULT(29) = AVMULT(29)+1.0D0
         IF (ICH(I).NE.0) AVE(27) = AVE(27)+HE(I)
         IF (ICH(I).NE.0) AVMULT(27) = AVMULT(27)+1.0D0
C        XFL = PZ(I)/PO
         XFL  = PZ(I)/P0
         XFLE = HE(I)/P0
C        IF ((ICH(I).NE.0).AND.(XFL.GT.-0.84D0).AND.(XFL.LT.-0.44D0))
C    &      WRITE(LOUT,*)'WRONG EVENT',NUMEV
         IXFL = XFL/DXFL+26
C        IF (XFL.LT.0.0D0) WRITE(LLOOK,'(2F15.5,I3)')XFL,PZ(I),IXFL
         IF (IXFL.LT.1 ) IXFL=1
         IF (IXFL.GT.50) IXFL=50
         XXXFL = ABS(XFL)
         IF (NRE.EQ.14) THEN
            NPI = NPI+1
            YXFLPI(IXFL) = YXFLPI(IXFL)+XFLE
         ENDIF
         IF ((ICH(I).GT.0).AND.(NRE.NE.1)) THEN
            NCH = NCH+1
            YXFLCH(IXFL) = YXFLCH(IXFL)+XFLE
         ENDIF
         IF (ICH(I).NE.0) YXFL(IXFL,9) = YXFL(IXFL,9)+XXXFL
         YXFL(IXFL,NI) = YXFL(IXFL,NI)+XXXFL
         YXFL(IXFL,10) = YXFL(IXFL,10)+XXXFL
         PTT = PX(I)**2+PY(I)**2
         AMT = SQRT(PTT+AM(I)**2)
         YL = 0.5D0*LOG(ABS((HE(I)+PZ(I)+1.D-10)
     &        /(HE(I)-PZ(I)+1.D-10)))
         YLPS = LOG(ABS((PZ(I)+SQRT(PZ(I)**2+PTT))
     &        /SQRT(PTT+1.D-6)+1.D-18))
         IYLPS = (YLPS+25.0D0*DY)/DY
         IF (IYLPS.LT.1)  IYLPS = 1
         IF (IYLPS.GT.51) IYLPS = 51
         YYLPS(IYLPS,NI) = YYLPS(IYLPS,NI)+1.0D0
         YYLPS(IYLPS,10) = YYLPS(IYLPS,10)+1.0D0
         IF (ICH(I).NE.0) YYLPS(IYLPS,9) = YYLPS(IYLPS,9)+1.0D0
         IYL = (YL+25.0D0*DY)/DY
         IF (IYL.LT.1)  IYL = 1
         IF (IYL.GT.51) IYL = 51
         IF (ICH(I).NE.0) YYL(IYL,9) = YYL(IYL,9)+1.0D0
         YYL(IYL,NI) = YYL(IYL,NI)+1.0D0
         YYL(IYL,10) = YYL(IYL,10)+1.0D0
   20 CONTINUE
      KPL = (AMCH+DDM+DDM/2.0D0)/DDM
      IF (KPL.GT.24) KPL=24
      IF (KPL.LT.1)  KPL=1
      ITT = 10.0D0*ABS(TT)+1
      IF (ITT.GT.40) ITT = 40
C     DSDTDM(ITT,KPL) = DSDTDM(ITT,KPL)+1.0D0/AMCH
      RETURN
*
    3 CONTINUE
      DO 30 J=1,51
         IF(NCH.NE.0) YXFLCH(J) = YXFLCH(J)/(NCH*DXFL)
         IF(NPI.NE.0) YXFLPI(J) = YXFLPI(J)/(NPI*DXFL)
         DO 31 I=1,10
            YXFL(J,I) = LOG10(ABS(YXFL(J,I) /(NEVT*DXFL))+1.0D-8)
            YYL(J,I)  = YYL(J,I) /(NEVT*DY)
            YYLPS(J,I)= YYLPS(J,I)/(NEVT*DY)
   31    CONTINUE
   30 CONTINUE
      DO 32 I=1,30
         AVMULT(I) = AVMULT(I)/NEVT
         AVE( I)   = AVE(I)   /NEVT
   32 CONTINUE
      DO 33 I=1,40
         DO 34 J=1,24
            DSDTDM(I,J) = DSDTDM(I,J)/NEVT
            DSDTDM(I,J) = LOG10(DSDTDM(I,J))
            DSDM(J,I)   = DSDTDM(I,J)
   34    CONTINUE
   33 CONTINUE
      WRITE(LOUT,*)'NAME,AVE,AVMULT'
      DO 105 I=1,30
         WRITE(LOUT,210) ANAME(I),AVE(I),AVMULT(I)
  210    FORMAT(' ',A8,2F15.5)
  105 CONTINUE
      WRITE(LOUT,100)
  100 FORMAT('1 RAPIDITY DISTRIBUTION')
      DO 110 J=1,50
         WRITE(LOUT,200) XYL(J,1),(YYL(J,I),I=1,10)
c        WRITE(LOUT,200) XYL(J,1),(YYL(J,I),I=11,20)
  200    FORMAT (F10.2,10E11.3)
  110 CONTINUE
      CALL PLOT(XYL,YYL,510,10,51,-25.*DY,DY,0.,0.03)
      WRITE(LOUT,101)
  101 FORMAT('1 PSEUDORAPIDITY DISTRIBUTION')
      CALL PLOT(XYL,YYLPS,510,10,51,-25.*DY,DY,0.,0.03)
      WRITE(LOUT,102)
  102 FORMAT ('1  LONG MOMENTUM (SCALED) DISTRIBUTION (LOG)')
      CALL PLOT(XXFL,YXFL,510,10,51,-1.,DXFL,-3.5,0.05)
      WRITE(LOUT,103)
  103 FORMAT ('1DISTRIBUTION DS/DTDM AS FUNCTION OF T ')
      CALL PLOT(TDTDM,DSDTDM,960,24,40,0.,0.04,-5.,0.05)
      WRITE(LOUT,104)
  104 FORMAT ('1DISTRIBUTION DS/DTDM AS FUNCTION OF M**2 ')
      CALL PLOT(AMDM,DSDM,960,40,24,2.,0.1,-5.,0.05)
      IF (IJPROJ.EQ.13) SIG = 66.0D0
      IF (IJPROJ.EQ.14) SIG = 66.0D0
      IF (IJPROJ.EQ.15) SIG = 54.9D0
      IF (IJPROJ.EQ.1)  SIG = 95.0D0
      SIGCH = 84.8
      SIGPIM=0.23D0
C     OPEN(19,FILE='FEYNDI.OUT')
      WRITE(LOUT,*)'FEYNMAN-DISTRIBUTION FOR PION-'
      DO 35 I=1,51
         WRITE(LOUT,'(2F15.5)')
     &   DXFL*(I-1)-1,SIG*YXFLPI(I)
         WRITE(LOUT,'(2F15.5)')
     &   DXFL*I-1,SIG*YXFLPI(I)
   35 CONTINUE
      WRITE(LOUT,*)'FEYNMAN-DISTRIBUTION FOR CHARGED PARTICLES'
      DO 36 I=1,51
         WRITE(LOUT,'(2F15.5)')
     &   DXFL*(I-1)-1,SIGCH*YXFLCH(I)
         WRITE(LOUT,'(2F15.5)')
     &   DXFL*I-1,SIGCH*YXFLCH(I)
   36 CONTINUE
C     CLOSE(19)
      RETURN
      END
*
*===sihndi===============================================================*
*
      SUBROUTINE SIHNDI(ECM,KPROJ,KTARG,SIGDIF,SIGDIH)
 
**********************************************************************
*   Single diffractive hadron-nucleon cross sections                 *
*                                              S.Roesler 14/1/93     *
*                                                                    *
*   The cross sections are calculated from extrapolated single       *
*   diffractive antiproton-proton cross sections (DTUJET92) using    *
*   scaling relations between total and single diffractive cross     *
*   sections.                                                        *
**********************************************************************
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      CHARACTER*8  ANAME
      COMMON /DPAR/ ANAME(210),AM(210),GA(210),TAU(210),IICH(210),
     &              IIBAR(210),K1(210),K2(210)
*
      CSD1   =   4.201483727
      CSD4   = -0.4763103556E-02
      CSD5   =  0.4324148297
*
      CHMSD1 =  0.8519297242
      CHMSD4 = -0.1443076599E-01
      CHMSD5 =  0.4014954567
*
      EPN = (ECM**2 - AM(KPROJ)**2 - AM(KTARG)**2)/(2.0D0*AM(KTARG))
      PPN = SQRT((EPN-AM(KPROJ))*(EPN+AM(KPROJ)))
*
      SDIAPP = CSD1+CSD4*LOG(PPN)**2+CSD5*LOG(PPN)
      SHMSD  = CHMSD1+CHMSD4*LOG(PPN)**2+CHMSD5*LOG(PPN)
      FRAC   = SHMSD/SDIAPP
*
      GOTO( 10, 20,999,999,999,999,999, 10, 20,999,
     &     999, 20, 20, 20, 20, 20, 10, 20, 20, 10,
     &      10, 10, 20, 20, 20) KPROJ
*
   10 CONTINUE
*---------------------------- p - p , n - p , sigma0+- - p ,
*                             Lambda - p
      CSD1   =  6.004476070
      CSD4   = -0.1257784606E-03
      CSD5   =  0.2447335720
      SIGDIF = CSD1+CSD4*LOG(PPN)**2+CSD5*LOG(PPN)
C     Replace SIGDIF with dpmjet result SIPPSD(ECM)
      SIGDIF=SIPPSD(ECM)
      SIGDIH = FRAC*SIGDIF
      RETURN
*
   20 CONTINUE
*
      KPSCAL = 2
      KTSCAL = 1
      F      = SDIAPP/DSHNTO(KPSCAL,KTSCAL,ECM)
C     F      = SDIAPP/DSHNEL(KPSCAL,KTSCAL,ECM)
      KT     = 1
      SIGDIF = DSHNTO(KPROJ,KT,ECM)*F
C     SIGDIF = DSHNEL(KPROJ,KT,ECM)*F
      SIGDIH = FRAC*SIGDIF
      RETURN
*
  999 CONTINUE
*-------------------------- leptons..
      SIGDIF = 1.E-10
      SIGDIH = 1.E-10
      RETURN
      END
*
*===dshnto===============================================================*
*
      DOUBLE PRECISION FUNCTION DSHNTO(KPROJ,KTARG,UMO)
 
**********************************************************************
*   Total hadron-nucleon cross sections         S.Roesler 12/1/93    *
*                                                                    *
*   Fits from Rev. of Part. Prop.(1992) are used in the momentum     *
*   ranges given there. Cross sections for momenta below them are    *
*   calculated as in an earlier version of this function :           *
*               SIG(el) + SIG(inel) = SIHNEL + SIHNIN .              *
*   Total antiproton-proton cross sections calculated with DTUJET92  *
*   parametrize the cross sections at higher enrgies.                *
**********************************************************************
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      CHARACTER*8  ANAME
      COMMON /DPAR/ ANAME(210),AAM(210),GA(210),TAU(210),IICH(210),
     &              IIBAR(210),K1(210),K2(210)
      COMMON /STRUFU/ISTRUM,ISTRUT
      DIMENSION SQS(20),SIIV(20),SQS22(40),SIIV22(40)
      DATA SQS /20.,50.,100.,200.,500.,1000.,1500.,2000.,3000.,
     *4000.,6000.,8000.,10000.,15000.,20000.,30000.,40000.,
     *60000.,80000.,100000/
      DATA SIIV /41.6,44.6,47.9,52.3,60.,67.2,71.4,75.,79.,82.4,
     *87.2,90.4,93.2,97.3,100.7,104.7,107.9,111.7,114.7,117.2/
      DATA SQS22 /
     *53.,69.,91.,119.,156.,205.,268.,351.,460.,603.,
     *790.,1036.,1357.,1778.,2329.,
     *3053.,3999.,5239.,6865.,8994.,
     *11785.,15441.,20232.,26509.,34733.,
     *45509.,59627.,78126.,102365.,134123.,
     *175734.,230255.,301690.,395288.,517925.,
     *678609.,889144.,1164997.,1526432.,2000000.
     */
      DATA SIIV22 /
     *44.3,45.3,46.5,47.8,49.3,51.0,53.0,55.3,57.9,60.9,
     *64.3,68.0,72.0,76.3,80.8,85.4,90.0,94.7,99.3,103.9,
     *108.4,112.8,117.2,121.4,125.6,
     *129.8,133.9,138.0,142.0,146.1,
     *150.2,154.3,158.5,162.7,166.9, 
     *171.2,175.6,180.0,184.6,189.1
     */
*
      F1 = 1.0D0
      CA = 0.0D0
      CB = 0.0D0
      CC = 0.0D0
      CD = 0.0D0
      CN = 0.0D0
*
      A1 = 0.0D0
      A2 = 0.0D0
      A3 = 1.0D0
      A4 = 0.0D0
      A5 = 0.0D0
      A6 = 0.0D0
*
      PARAM1  =   34.94235992
      PARAM4  =  0.2104312854
      PARAM5  = -0.4509592056E-01
*
      IPIO  = 0
      SPIO1 = 0.0D0
      SPIO2 = 0.0D0
*
      UMO2   = UMO**2
      EPN = (UMO**2-AAM(KPROJ)**2-AAM(KTARG)**2)/(2.0D0*AAM(KTARG))
      PO  = SQRT((EPN-AAM(KPROJ))*(EPN+AAM(KPROJ)))
*
    1 CONTINUE
*
      IF (KTARG.EQ.8) THEN
           GOTO( 30, 40,999,999,999,999,999, 10, 20,999,
     &          999,140, 70, 60,150,160,100, 20,140, 10,
     &           10, 10,110,130,120) KPROJ
      ELSE
           GOTO( 10, 20,999,999,999,999,999, 30, 40,999,
     &          999, 50, 60, 70, 80, 90,100, 20, 50, 10,
     &           10, 10,110,120,130) KPROJ
      ENDIF
*
   10 CONTINUE
*---------------------------- p - p , sigma0+- - p
      IF (PO.LE.3.0D0) THEN
         GOTO 500
      ELSE IF ((PO.GT.3.0D0).AND.(UMO.LE.100.0D0)) THEN
         CA = 48.0D0
         CC = 0.522D0
         CD = -4.51D0
         GOTO 600
      ELSE
         GOTO 700
      ENDIF
*
   20 CONTINUE
*---------------------------- pbar - p , Lambdabar - p
      IF (PO.LE.5.0D0) THEN
         GOTO 500
      ELSE IF ((PO.GT.5.0D0).AND.(UMO.LE.200.0D0)) THEN
         CA = 38.4D0
         CB = 77.6
         CC = 0.26D0
         CD = -1.2D0
         CN = -0.64D0
         GOTO 600
      ELSE
         GOTO 700
      ENDIF
*
   30 CONTINUE
*---------------------------- n - p
      IF (PO.LE.3.0D0) THEN
         GOTO 500
      ELSE IF ((PO.GT.3.0D0).AND.(PO.LE.370.0D0)) THEN
         CA = 47.3D0
         CC = 0.513D0
         CD = -4.27D0
         GOTO 600
      ELSE IF ((PO.GT.370.0D0).AND.(UMO.LE.110.0D0)) THEN
         A1 = 38.5D0
         A2 = 0.46D0
         A3 = 125.0D0
         A4 = 15.0D0
         GOTO 800
      ELSE
         GOTO 700
      ENDIF
*
   40 CONTINUE
*---------------------------- nbar - p
      IF (PO.LE.1.1D0) THEN
         GOTO 500
      ELSE IF ((PO.GT.1.1D0).AND.(PO.LE.280.0D0)) THEN
         CB = 133.6D0
         CC = -1.22D0
         CD = 13.7D0
         CN = -0.7D0
         GOTO 600
      ELSE IF ((PO.GT.280.0D0).AND.(UMO.LE.110.0D0)) THEN
         A1 = 38.5D0
         A2 = 0.46D0
         A3 = 125.0D0
         A4 = 15.0D0
         A5 = 77.43D0
         A6 = -0.6D0
         GOTO 800
      ELSE
         GOTO 700
      ENDIF
*
   50 CONTINUE
*---------------------------- Klong - p , Kshort - p
      R = RNDM(V)
      IF (R.LE.0.5D0) THEN
*                             K+ - p
         GOTO 80
      ELSE
*                             K- - p
         GOTO 90
      ENDIF
*
   60 CONTINUE
*---------------------------- pi+ - p
      IF (PO.LE.4.0D0) THEN
         GOTO 500
      ELSE IF ((PO.GT.4.0D0).AND.(PO.LE.340.0D0)) THEN
         CA = 16.4D0
         CB = 19.3D0
         CC = 0.19D0
         CN = -0.42D0
         GOTO 600
      ELSE IF ((PO.GT.340.0D0).AND.(UMO.LE.47.0D0)) THEN
         A1 = 24.0D0
         A2 = 0.6D0
         A3 = 160.0D0
         A5 = -7.9D0
         A6 = -0.46D0
         GOTO 800
      ELSE
         F1 = 2.0D0/3.0D0
         GOTO 10
      ENDIF
*
   70 CONTINUE
*---------------------------- pi- - p
      IF (PO.LE.2.5D0) THEN
         GOTO 500
      ELSE IF ((PO.GT.2.5D0).AND.(PO.LE.370.0D0)) THEN
         CA = 33.0D0
         CB = 14.0D0
         CC = 0.456D0
         CD = -4.03D0
         CN = -1.36D0
         GOTO 600
      ELSE IF ((PO.GT.370.0D0).AND.(UMO.LE.47.0D0)) THEN
         A1 = 24.0D0
         A2 = 0.6D0
         A3 = 160.0D0
         GOTO 800
      ELSE
         F1 = 2.0D0/3.0D0
         GOTO 10
      ENDIF
*
   80 CONTINUE
*---------------------------- K+ - p
      IF (PO.LE.2.0D0) THEN
         GOTO 500
      ELSE IF ((PO.GT.2.0D0).AND.(PO.LE.310.0D0)) THEN
         CA = 18.1D0
         CC = 0.26D0
         CD = -1.0D0
         GOTO 600
      ELSE IF ((PO.GT.310.0D0).AND.(UMO.LE.110.0D0)) THEN
         A1 = 20.3D0
         A2 = 0.59D0
         A3 = 140.0D0
         A5 = -30.13D0
         A6 = -0.58D0
         GOTO 800
      ELSE
         F1 = 2.0D0/3.0D0
         GOTO 10
      ENDIF
*
   90 CONTINUE
*---------------------------- K- - p
      IF (PO.LE.3.0D0) THEN
         GOTO 500
      ELSE IF ((PO.GT.3.0D0).AND.(PO.LE.310.0D0)) THEN
         CA = 32.1D0
         CC = 0.66D0
         CD = -5.6D0
         GOTO 600
      ELSE IF ((PO.GT.310.0D0).AND.(UMO.LE.110.0D0)) THEN
         A1 = 20.3D0
         A2 = 0.59D0
         A3 = 140.0D0
         GOTO 800
      ELSE
         F1 = 2.0D0/3.0D0
         GOTO 10
      ENDIF
*
  100 CONTINUE
*---------------------------- Lambda - p
      IF (PO.LE.0.6D0) THEN
         GOTO 500
      ELSE IF ((PO.GT.0.6D0).AND.(PO.LE.21.0D0)) THEN
         CA = 30.4D0
         CD = 1.6D0
         GOTO 600
      ELSE
         GOTO 10
      ENDIF
*
  110 CONTINUE
*---------------------------- pi0 - p
*                             1/2(pi+p + pi-p)
      IPIO  = 1
      KPROJ = 13
      GOTO 1
*
  120 CONTINUE
*---------------------------- K0 - p
      IF (PO.LE.2.0D0) THEN
         GOTO 500
      ELSE IF ((PO.GT.2.0D0).AND.(PO.LE.310.0D0)) THEN
*                             K+ - n
         CA = 18.7D0
         CC = 0.21D0
         CD = -0.89D0
         GOTO 600
      ELSE
         GOTO 80
      ENDIF
*
  130 CONTINUE
*---------------------------- K0bar - p
      IF (PO.LE.1.8D0) THEN
         GOTO 500
      ELSE IF ((PO.GT.1.8D0).AND.(PO.LE.310.0D0)) THEN
*                             K- - n
         CA = 25.2D0
         CC = 0.38D0
         CD = -2.9D0
         GOTO 600
      ELSE
         GOTO 90
      ENDIF
*
  140 CONTINUE
*---------------------------- Klong - n , Kshort - n
      R = RNDM(V)
      IF (R.LE.0.5D0) THEN
*                             K+ - n
         GOTO 150
      ELSE
*                             K- - n
         GOTO 160
      ENDIF
*
  150 CONTINUE
*---------------------------- K+ - n
      IF (PO.LE.2.0D0) THEN
         GOTO 500
      ELSE IF ((PO.GT.2.0D0).AND.(PO.LE.310.0D0)) THEN
         CA = 18.7D0
         CC = 0.21D0
         CD = -0.89D0
         GOTO 600
      ELSE
         GOTO 90
      ENDIF
*
  160 CONTINUE
*---------------------------- K- - n
      IF (PO.LE.1.8D0) THEN
         GOTO 500
      ELSE IF ((PO.GT.1.8D0).AND.(PO.LE.310.0D0)) THEN
         CA = 25.2D0
         CC = 0.38D0
         CD = -2.9D0
         GOTO 600
      ELSE
         GOTO 80
      ENDIF
*
  500 CONTINUE
      CALL SIHNEL(KPROJ,KTARG,PO,SEL)
      CALL SIHNIN(KPROJ,KTARG,PO,SIN)
      STOT = SEL + SIN
      GOTO 900
*
  600 CONTINUE
      STOT = F1*(CA+CB*PO**CN+CC*(LOG(PO))**2+CD*LOG(PO))
      GOTO 900
*
  700 CONTINUE
      IF(ISTRUM.EQ.14.AND.ISTRUT.EQ.2)THEN
	DO 701 I=1,20
	  IF(UMO.LE.SQS(I))GO TO 702
  701   CONTINUE
	I=20
  702   CONTINUE
	IF ((I.EQ.20).AND.(UMO.GT.SQS(20)))THEN
	  TEPP=SIIV(20)+(LOG(UMO)-LOG(SQS(20)))*(SIIV(20)-SIIV(19))/
     *         (LOG(SQS(20))-LOG(SQS(19)))
          STOT=F1*TEPP
	ELSE
	  TEPP=SIIV(I-1)+(UMO-SQS(I-1))*(SIIV(I)-SIIV(I-1))/
     *         (SQS(I)-SQS(I-1))
          STOT=F1*TEPP
	ENDIF
      ELSEIF(ISTRUM.EQ.22.AND.ISTRUT.EQ.2)THEN
	DO 711 I=1,40
	  IF(UMO.LE.SQS22(I))GO TO 712
  711   CONTINUE
	I=40
  712   CONTINUE
	IF ((I.EQ.40).AND.(UMO.GT.SQS22(40)))THEN
	  TEPP=SIIV22(40)+(LOG(UMO)-LOG(SQS22(40)))*
     *	  (SIIV22(40)-SIIV22(39))/
     *         (LOG(SQS22(40))-LOG(SQS22(39)))
          STOT=F1*TEPP
	ELSE
	  TEPP=SIIV22(I-1)+(UMO-SQS22(I-1))*(SIIV22(I)-SIIV22(I-1))/
     *         (SQS22(I)-SQS22(I-1))
          STOT=F1*TEPP
	ENDIF
      ELSE
        STOT = F1*(PARAM1+PARAM4*LOG(PO)**2+PARAM5*LOG(PO))
      ENDIF
      GOTO 900
*
  800 CONTINUE
      STOT = F1*(A1+A2*(LOG(UMO2/A3))**2+A4/UMO2+A5*UMO2**A6)
*
  900 CONTINUE
      IF ((IPIO.EQ.1).AND.(KPROJ.EQ.13)) THEN
         SPIO1 = STOT
         KPROJ = 14
         GOTO 1
      ENDIF
      IF ((IPIO.EQ.1).AND.(KPROJ.EQ.14)) THEN
         SPIO2 = STOT
         STOT = 0.5D0*(SPIO1+SPIO2)
      ENDIF
      DSHNTO = STOT
      RETURN
*
  999 CONTINUE
*-------------------------- leptons..
      DSHNTO = 1.E-10
      RETURN
      END
*
*===dshnel===============================================================*
*
      DOUBLE PRECISION FUNCTION DSHNEL(KPROJ,KTARG,UMO)
 
**********************************************************************
*   Elastic hadron-nucleon cross sections         S.Roesler 13/1/93  *
*                                                                    *
*   SIHNEL is called for c.m. energies below 50 GeV. Otherwise       *
*   a parametrisation of antiproton-proton elastic cross sections    *
*   obtained from DTUJET92 is used.                                  *
**********************************************************************
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      CHARACTER*8  ANAME
      COMMON /DPAR/ ANAME(210),AAM(210),GA(210),TAU(210),IICH(210),
     &              IIBAR(210),K1(210),K2(210)
*
      F1 = 1.0D0
      CA = 0.0D0
      CB = 0.0D0
      CC = 0.0D0
      CD = 0.0D0
      CN = 0.0D0
*
      IPIO  = 0
      SPIO1 = 0.0D0
      SPIO2 = 0.0D0
*
      PARAM1  =   7.789333344
      PARAM4  =  0.7488331199E-01
      PARAM5  = -0.6963931322
*
      UMO2   = UMO**2
      EPN = (UMO**2-AAM(KPROJ)**2-AAM(KTARG)**2)/(2.0D0*AAM(KTARG))
      PO  = SQRT((EPN-AAM(KPROJ))*(EPN+AAM(KPROJ)))
*
    1 CONTINUE
*
      GOTO( 10, 10,999,999,999,999,999, 10, 10,999,
     &     999, 30, 20, 20, 20, 20, 10, 10, 30, 10,
     &      10, 10, 40, 20, 20) KPROJ
*
   10 CONTINUE
      IF (UMO.LE.50.0D0) THEN
         CALL SIHNEL(KPROJ,KTARG,PO,SEL)
         GOTO 900
      ELSE
         F1 = 1.0D0
         GOTO 500
      ENDIF
*
   20 CONTINUE
      IF (UMO.LE.50.0D0) THEN
         CALL SIHNEL(KPROJ,KTARG,PO,SEL)
         GOTO 900
      ELSE
         F1 = 2.0D0/3.0D0
         GOTO 500
      ENDIF
*
   30 CONTINUE
*---------------------------- Klong - p , Kshort - p
      R = RNDM(V)
      IF (R.LE.0.5D0) THEN
*                             K+ - p
         KPROJ = 15
      ELSE
*                             K- - p
         KPROJ = 16
      ENDIF
      GOTO 1
*
   40 CONTINUE
*---------------------------- pi0 - p
*                             1/2(pi+p + pi-p)
      IPIO  = 1
      KPROJ = 13
      GOTO 1
*
  500 CONTINUE
      SEL = F1*(PARAM1+PARAM4*LOG(PO)**2+PARAM5*LOG(PO))
*
  900 CONTINUE
      IF ((IPIO.EQ.1).AND.(KPROJ.EQ.13)) THEN
         SPIO1 = SEL
         KPROJ = 14
         GOTO 1
      ENDIF
      IF ((IPIO.EQ.1).AND.(KPROJ.EQ.14)) THEN
         SPIO2 = SEL
         SEL = 0.5D0*(SPIO1+SPIO2)
      ENDIF
      DSHNEL = SEL
      RETURN
*
  999 CONTINUE
*-------------------------- leptons..
      DSHNEL = 1.E-10
      RETURN
      END
