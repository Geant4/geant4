*
*===kkinc==============================================================*
*
**sr mod. for DPMJET: parameter list
      SUBROUTINE KKINC(EPN,NTMASS,NTCHAR,NPMASS,NPCHAR,IDP,KKMAT,
     *IDT, NHKKH1,IREJ)

************************************************************************
* Treatment of complete nucleus-nucleus or hadron-nucleus scattering   *
* This subroutine is an update of the previous version written         *
* by J. Ranft/ H.-J. Moehring.                                         *
* This version dated 19.11.95 is written by S. Roesler                 *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TINY5=1.0D-5,
     &           TINY2=1.0D-2,TINY3=1.0D-3)

      LOGICAL LFZC

      PARAMETER (NMXHKK=89998) 
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)
      COMMON /EXTEVT/ IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),
     &                IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10)
      CHARACTER*8  ANAME
      COMMON /DPAR/   ANAME(210),AAM(210),GA(210),TAU(210),
     &                IICH(210),IIBAR(210),K1(210),K2(210)

      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
**sr mod. for DPMJET: EPROJ needed
      COMMON /NUCCMS/ GACMS,BGCMS,GALAB,BGLAB,BLAB,UMO,PCM,EPROJ,PPROJ
**sr mod. for DPMJET: commons added
      COMMON /FINAL/  IFINAL
      COMMON /CMHICO/ CMHIS
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
      COMMON /CHABAI/CHARGI,BARNUI
      COMMON /NOMIJE/ PTMIJE(10),NNMIJE(10)
      COMMON /NNCMS/  GAMCM,BGCM,UMOL,PCML,EPROJL,PPROJL
      COMMON /EDENS/IEDEN
      COMMON /XDIDID/XDIDI
      COMMON /NSTARI/NSTART
      COMMON/PYJETS/NLU,NPAD,KLU(4000,5),PLU(4000,5),VLU(4000,5)
      COMMON /FELIRE/AMRECD,KJPRO
      COMMON /NEUTYY/NEUTYP,NEUDEC
      COMMON /TAUFO/  TAUFOR,KTAUGE,ITAUVE,INCMOD
**
**sr mod. for DPMJET: set output flags
      DATA KKCOUN /0/
      DATA CHCOUN /0/
      DATA FICOUN /0/
      DATA TAUCOU /0/
      IF((TAUCOU.EQ.0).AND.((IT.EQ.1).AND.(IP.EQ.1)))THEN
        TAUCOU=TAUCOU+1
	KTAUGE=0
      ENDIF
      IF(IPEV.GE.1)THEN
      WRITE(6,*)'kkinc EPN,NTMASS,NTCHAR,NPMASS,NPCHAR,IDP,KKMAT',
     *'IDT, NHKKH1,IREJ',
     * EPN,NTMASS,NTCHAR,NPMASS,NPCHAR,IDP,KKMAT,
     *IDT, NHKKH1,IREJ
      ENDIF
      IPRI=0
      IREJ=0
C--------------------------------------------------------------------
 1889 CONTINUE
      KKCOUN=KKCOUN+1
      NEVHKK=KKCOUN
C     DO 5371 IOK=1,200
      IF(IPRI.GE.1)WRITE(6,'(A,I10)') ' KKINC: KKCOUN=',KKCOUN
      IF(IPRI.GE.1)WRITE(6,'(A,E20.8)') ' KKINC: EPN=',EPN
C5371 CONTINUE
*---redefine characteristics of the actual interaction
      IF(KKCOUN.EQ.-721.OR.KKCOUN.EQ.-821)THEN
        IOUXVO=IOUXEV
        IOUXEV=6
        IPEVO=IPEV
        IPEV=6
        IPPAO=IPPA
        IPPA=2
        IPCOO=IPCO
        IPCO=6
        INITO=INIT
C       INIT=2
        IPRIO=IPRI
        IPRI=6
        IPHKKO=IPHKK
C       IPHKK=6
      ENDIF
      IF(KKCOUN.EQ.-39.OR.KKCOUN.EQ.-822)THEN
        IOUXEV=IOUXVO
        IPEV=IPEVO
        IPPA=IPPAO
        IPCO=IPCOO
        INIT=INITO
        IPRI=IPRIO
        IPHKK=IPHKKO
      ENDIF
**
**sr mod. for DPMJET: minijet-statist. added
C                        NUMBER of JETS in event
      DO IIII=1,10
        NNMIJE(IIII)=0
      ENDDO
**
      ILOOP = 0
  100 CONTINUE
      IREJ=0
      IREJ1=0
      IF (ILOOP.EQ.40)THEN
	WRITE(6,'(A)')'  Rejection after 40 trials'
	IREJ=1
        RETURN
      ENDIF
      ILOOP = ILOOP+1

* re-initialize /NUCC/
      IP  = NPMASS
      IPZ = NPCHAR
      IT  = NTMASS
      ITZ = NTCHAR
      IJPROJ = IDP
      IF(NEUDEC.GE.10)IJPROJ=5
      IF(NSTART.EQ.4.OR.NSTART.EQ.2)IJPROJ=5
      IJTARG = IDT
      IBPROJ = IIBAR(IJPROJ)
      IBTARG = IIBAR(IJTARG)

**sr mod. for DPMJET: quantum number check added
C                        Event Charge and Baryon number
      CHARGI=ITZ
      BARNUI=IT
      IF(IP.GT.1)THEN
        CHARGI=CHARGI+IPZ
        BARNUI=BARNUI+IP
      ELSE
        CHARGI=CHARGI+IICH(IJPROJ)
        BARNUI=BARNUI+IBPROJ
      ENDIF
**sr mod. for DPMJET: initialize /EXTEVT/ and /NUCCMS/
      IF(IPEV.GE.1)WRITE(6,*)' before EVTINI call'
      CALL EVTINI(IJPROJ,IP,IT,EPN,PPN,ECM,NHKKH1,1)
      IF(IPEV.GE.1)WRITE(6,*)' after EVTINI call EPN',EPN
**

* calculate nuclear potentials (common /NUCLEA/)
      IF(IPEV.GE.1)WRITE(6,*)' before NCLPOT call'
       IF(IP.GT.1.OR.IT.GT.1)THEN      
      CALL NCLPOT(IPZ,IP,ITZ,IT,ZERO,ZERO,0)
       ENDIF
      IF(IPEV.GE.1)WRITE(6,*)' after NCLPOT call'

* initialize treatment for residual nuclei
      IF(NSTART.NE.2)THEN
      IF(IPEV.GE.1)WRITE(6,*)' before RESNCL call'
       IF(IP.GT.1.OR.IT.GT.1)THEN      
      CALL RESNCL(EPN,1)
       ENDIF
      IF(IPEV.GE.1)WRITE(6,*)' after RESNCL call EPN',EPN
      ENDIF

* sample hadron/nucleus-nucleus interaction
**sr mod. for DPMJET: parameter list
      IF(IPRI.GE.1)WRITE(6,'(A,2E20.8,2I5)') ' KKINC call KKEVT: ',
     * EPROJ,PPROJ,KKMAT,IREJ1
C
      IF(NSTART.EQ.1)THEN
C                           h-h, h-A, A-A Collisions
        CALL KKEVT(NHKKH1,EPROJ,PPROJ,KKMAT,IREJ1)
      ELSEIF(NSTART.EQ.2)THEN
C                           Neutrino-A Collisions (qeld code)
        CALL KKEVNU(NHKKH1,EPROJ,PPROJ,KKMAT,IREJ1,ECM)
      ELSEIF(NSTART.EQ.3)THEN
C                           Diffr Interactions with nuclei
        CALL KKEVDI(NHKKH1,EPROJ,PPROJ,KKMAT,IREJ1)
      ELSEIF(NSTART.EQ.4)THEN
C                           Neutrino-A Collisions (lepto code)
        CALL KKEVLE(NHKKH1,EPROJ,PPROJ,KKMAT,IREJ1)
      ENDIF
C
      IF(IPRI.GE.1)WRITE(6,'(A,2E20.8,2I5)') ' KKINC after KKEVT: ',
     * EPROJ,PPROJ,KKMAT,IREJ1
C	WRITE(6,'(A)')' KKEVT '
      IF (IREJ1.GT.0)THEN
 	WRITE(6,'(A,I5)')' KKEVT Rejection KKCOUN ',KKCOUN
        RETURN
      ENDIF
* initialize treatment for residual nuclei
      IF(NSTART.EQ.2)THEN
      IF(IPEV.GE.1)WRITE(6,*)' before RESNCL call'
       IF(IP .GT.1.OR.IT.GT.1)THEN      
      CALL RESNCL(EPN,1)
       ENDIF
      IF(IPEV.GE.1)WRITE(6,*)' after RESNCL call EPN',EPN
      ENDIF


**sr mod. for DPMJET: special ststistics
C       IF(IPRI.GE.1)THEN
C         DO 7735 IHKK=1,NHKK
C           WRITE(6,1000) IHKK, ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),
C    +      JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK
C    +      (KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
C7735     CONTINUE
C       ENDIF
C       IREJ=0
C       GOTO 100
C     ENDIF
C     IF (IPRI.GE.1.OR.IP.EQ.1)IREJ=0
C     IF (IRESO.EQ.1) CALL DISRES(2,NHKKH1,PPN)
      IF(IEDEN.EQ.0)CALL DECHKK(NHKKH1)
        IF(IPRI.GE.7)THEN
	  WRITE(6,'(A)')' from KKINC after DECHKK'
          DO 7835 IHKK=1,NHKK
            WRITE(6,1000) IHKK, ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),
     +      JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK
     +      (KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
 7835     CONTINUE
        ENDIF
**
**sr mod. for DPMJET: get some information for fzc
      IF(IPEV.GE.1)WRITE(6,*)' before EVTINI call'
      CALL EVTINI(IJPROJ,IP,IT,EPN,PPN,ECM,NHKKH1,2)
      IF(IPEV.GE.1)WRITE(6,*)' after EVTINI call'
**
* intranuclear cascade of final state particles for KTAUGE generations
* of secondaries
      IF(IPEV.GE.1)WRITE(6,*)' before FOZOCA call'
       IF(IP .GT.1.OR.IT.GT.1)THEN      
      CALL FOZOCA(LFZC,IREJ1)
       ENDIF
      IF(IPEV.GE.1)WRITE(6,*)' after fozoca LFZC,IREJ1',LFZC,IREJ1
      IF(IPEV.GE.1)WRITE(6,*)' after FOZOCA call'
      IF (IREJ1.GT.0)THEN
 	WRITE(6,'(A)')' FOZOCA Rejection'
        RETURN
      ENDIF

* baryons unable to escape the nuclear potential are treated as
* excited nucleons (ISTHKK=15,16)
      IF(IPEV.GE.1)WRITE(6,*)' before SCN4BA call'
       IF(IP .GT.1.OR.IT.GT.1)THEN      
      CALL SCN4BA
       ENDIF
      IF(IPEV.GE.1)WRITE(6,*)' after SCN4BA call'

* decay of resonances produced in intranuclear cascade processes
**sr 15-11-95 should be obsolete
      IF (LFZC) CALL DECAY1

* treatment of residual nuclei
      IF(IPEV.GE.1)WRITE(6,*)' before RESNCL call'
       IF(IP .GT.1.OR.IT.GT.1)THEN      
      CALL RESNCL(EPN,2)
       ENDIF
      IF(IPEV.GE.1)WRITE(6,*)' after RESNCL call'

* evaporation / fission / fragmentation
* (if intranuclear cascade was sampled only)
**sr mod. for DPMJET: check for IFINAL-flag
      IF ((LFZC).AND.(IFINAL.EQ.0)) THEN
        IF(IPRI.GE.1)THEN
	  WRITE(6,'(A)')' from KKINC before FICONF'
          DO 7935 IHKK=1,NHKK
            WRITE(6,1005) IHKK, ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),
     +      JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK
     +      (KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
     +      ,IDRES(IHKK),IDXRES(IHKK),NOBAM(IHKK),
     &                IDBAM(IHKK),IDCH(IHKK)
 1005 FORMAT (I6,I4,5I6,9(1PE10.2)/5I6)
 7935     CONTINUE
        ENDIF
      IF(IPEV.GE.1)WRITE(6,*)' before FICONF call'
       IF(IP .GT.1.OR.IT.GT.1)THEN      
         CALL FICONF(IJPROJ,IP,IPZ,IT,ITZ,IREJ1)
       ENDIF
      IF(IPEV.GE.1)WRITE(6,*)' after FICONF call IREJ1',IREJ1
C-----------------------------------------------------------
C                        Write events to file qeld.evt
C-----------------------------------------------------------
         IF (IREJ1.EQ.0.AND.NSTART.EQ.2) THEN
           IIII=0
           IIIMAX = 5
	     IF(NEUDEC.EQ.10.OR.NEUDEC.EQ.11)THEN
	       IIIMAX = 7
             ENDIF
           IF(KLU(1,2).EQ.16.OR.KLU(1,2).EQ.-16)THEN 
	     IF(NEUDEC.EQ.1.OR.NEUDEC.EQ.2)THEN
	       IIIMAX = 6
             ENDIF
           ENDIF
           DO 266 III=1,IIIMAX
             IF(KLU(III,1).EQ.1.OR.III.LE.2) THEN
               IIII=IIII+1
	       WRITE(29,'(3I6,5F10.3)')IIII,KLU(III,1),KLU(III,2),
     *         (PLU(III,KK),KK=1,5)
             ENDIF
  266      CONTINUE
           IIII=-1
           WRITE(29,'(I6)')IIII
	 ENDIF
C-----------------------------------------------------------
         IF (IREJ1.EQ.1) THEN
	   FICOUN=FICOUN+1
 	   IF(FICOUN.LE.20)WRITE(6,'(A)')' FICONF Rejection'
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
           IF(NSTART.EQ.3)THEN
             KFORM=2
             IF(KFORM.EQ.1)THEN
               AABBCC=0.
             ELSEIF(KFORM.EQ.2)THEN
C                      the following 3 lines only for 6 (J/psi)
               READ(29,'(1X,I5)')KREPA
               READ(29,'(1X,I5)')KREPA
               READ(29,'(1X,I5)')KREPA
C
               READ(29,'(1X,I5)')KREPA
               DO 1975 KRE=1,KREPA
                 READ(29,'(1X,A)')A109
 1975          CONTINUE
             ENDIF
           ENDIF
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

C	   GOTO 100
	 ENDIF
      ENDIF

**sr mod. for DPMJET: checks, histograms, ...
        IPHIHI=0
        DO 7501 IHKK=1,NHKKH1
          IF(IDHKK(IHKK).EQ.88888)THEN
C           PPTT=PHKK(1,IHKK)**2+PHKK(2,IHKK)**2
C           IF(PPTT.LE.1.D-12)THEN
C      IPHIHI=1
C          WRITE(6,*)'pt=0 IHKK,IDHKK(IHKK),PHKK(1,IHKK),PHKK(2,IHKK) ',
C    *         IHKK,ISTHKK(IHKK),IDHKK(IHKK),PHKK(1,IHKK),PHKK(2,IHKK)  
C	    ENDIF  
            IF(PHKK(5,IHKK).LE.1.D-10)THEN
	      IPHIHI=1
           WRITE(6,*)'M=0 IHKK,IDHKK(IHKK),PHKK(4,IHKK),PHKK(5,IHKK) ',
     *          IHKK,ISTHKK(IHKK),IDHKK(IHKK),PHKK(4,IHKK),PHKK(5,IHKK) 
	    ENDIF  
            IF(JMOHKK(1,IHKK).GE.IHKK)THEN
	      IPHIHI=1
           WRITE(6,*)'MO=0 IHKK,IDHKK(IHKK),PHKK(4,IHKK),PHKK(5,IHKK) ',
     *  IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),
     *   PHKK(4,IHKK),PHKK(5,IHKK) 
	    ENDIF  
          ENDIF
 7501   CONTINUE
        DO 501 IHKK=NHKKH1,NHKK
          IF(ISTHKK(IHKK).EQ.1)THEN
            PPTT=PHKK(1,IHKK)**2+PHKK(2,IHKK)**2
            IF(PPTT.LE.1.D-18)THEN
	      IPHIHI=1
              WRITE(6,*)' pt=0 IHKK,PHKK(1,IHKK),PHKK(2,IHKK) ',
     *                  IHKK,PHKK(1,IHKK),PHKK(2,IHKK)	      
	    ENDIF  
            IF(JMOHKK(1,IHKK).GT.IHKK)THEN
	      IPHIHI=1
           WRITE(6,*)'MO=0 IHKK,IDHKK(IHKK),PHKK(4,IHKK),PHKK(5,IHKK) ',
     *  IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),
     *   PHKK(4,IHKK),PHKK(5,IHKK) 
	    ENDIF  
            IF(IDHKK(IHKK).EQ.14.OR.IDHKK(IHKK).EQ.-14)THEN
C	      IPHIHI=1
           WRITE(6,*)'14-14IHKK,IDHKK(IHKK),PHKK(4,IHKK),PHKK(5,IHKK) ',
     *  IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),
     *   PHKK(4,IHKK),PHKK(5,IHKK) 
	    ENDIF  
          ENDIF
  501   CONTINUE
      IF (IPHIHI.GE.1) THEN
        WRITE(6,'(/A/)') ' KKINC: One particle with pt=0. !!!!'
      IF (IPHKK.GE.-1) THEN
        DO 502 IHKK=1,NHKK
          WRITE(6,1000) IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),
     +    JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK
     +    (KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
  502   CONTINUE
      ENDIF
      ENDIF
      IF (IPHKK.GE.1) THEN
        WRITE(6,'(/A/)') ' KKINC: FINAL LIST OF ENTRIES TO /HKKEVT/'
        DO 50 IHKK=1,NHKK
          WRITE(6,1000) IHKK,ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),
     +    JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK
     +    (KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
 1000 FORMAT (I6,I4,5I6,9(1PE10.2))
   50   CONTINUE
      ENDIF
C
C                        fix KTAUAC later
      KTAUAC=99
 
        IF(IPEV.GE.1)THEN
      WRITE(6,'(A,2F15.5)')' GACMS,BGCMS',GACMS,BGCMS
	ENDIF
C------------------------------------------------------------------
C                 Up to here the events (PHKK(J,I)) are in cms
C                 transform back to lab for cmhis=0 (lab histograms)
C
C                 But VHKK(J,I) is in Lab frame
C                 transform into cms for CMHIS >= 1
C------------------------------------------------------------------
        IF(KKCOUN.LE.-50)THEN
	  WRITE(6,*)' Event from dpmjet (only final particles):'
          WRITE(6,*)' before transf. into lab frame '
          DO 7737 IHKK=1,NHKK
	  IF((ISTHKK(IHKK).EQ.-1).OR.
     *	  (ISTHKK(IHKK).EQ.1).OR.
     *	  (ISTHKK(IHKK).EQ.1001))THEN
            WRITE(6,1055) IHKK, ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),
     +      JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK
     +      (KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
     +      , (WHKK(KHKK,IHKK),KHKK=1,4)
     +      ,IDRES(IHKK),IDXRES(IHKK),NOBAM(IHKK),
     &                IDBAM(IHKK),IDCH(IHKK)
	  ENDIF
 7737     CONTINUE
        ENDIF
      IF(IPEV.GE.1)WRITE(6,*)' before transf. into lab frame '
      DO 20 I=NHKKH1+1,NHKK
        PZNN=PHKK(3,I)
        ENN =PHKK(4,I)
	ZZZZ=VHKK(3,I)
	TTTT=VHKK(4,I)
        IF (CMHIS.EQ.0.D0)THEN
	  IF(ISTHKK(I).NE.16.AND.ISTHKK(I).NE.15)THEN 
          PHKK(3,I) = GACMS*PZNN + BGCMS*ENN
          PHKK(4,I) = GACMS*ENN  + BGCMS*PZNN
C         PHKK(3,I) = GAMCM*PZNN + BGCM*ENN
C         PHKK(4,I) = GAMCM*ENN  + BGCM*PZNN
	  ENDIF
        ENDIF
	IF(CMHIS.GE.1.D0)THEN
 	  VHKK(3,I) = GACMS*ZZZZ - BGCMS*TTTT
 	  VHKK(4,I) = GACMS*TTTT - BGCMS*ZZZZ
C	  VHKK(3,I) = GAMCM*ZZZZ - BGCM*TTTT
C         VHKK(4,I) = GAMCM*TTTT - BGCM*ZZZZ
	ENDIF
        EHECC=SQRT(PHKK(1,I)** 2+ PHKK(2,I)** 2+ PHKK(3,I)** 2+ PHKK
     +  (5,I)**2)
        IF (ABS(EHECC-PHKK(4,I)).GT.0.001) THEN
C            WRITE(6,'(2A/3I5,3E16.6)')
C    &         ' KKINC: CORRECT INCONSISTENT ENERGY ',
C    *         '  IEVCOU, I,IDHKK(I), PHKK(4,I),EHECC, PHKK(5,I)',
C    *            IEVCOU, I,IDHKK(I), PHKK(4,I),EHECC, PHKK(5,I)
          PHKK(4,I)=EHECC
        ENDIF
   20 CONTINUE
      IF(IPEV.GE.1)WRITE(6,*)' after transf. into lab frame '
        IF(IPEV.GE.1)THEN
C       IF ((LFZC).AND.(IFINAL.EQ.0)) THEN
	 IF(IPEV.GE.1) WRITE(6,'(A)')'  before CHECKF'
          DO 7135 IHKK=1,NHKK
            WRITE(6,1055) IHKK, ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),
     +      JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK
     +      (KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
     +      , (WHKK(KHKK,IHKK),KHKK=1,4)
     +      ,IDRES(IHKK),IDXRES(IHKK),NOBAM(IHKK),
     &                IDBAM(IHKK),IDCH(IHKK)
 1055 FORMAT (I6,I4,5I6/7(1PE11.3)/6(1PE11.3)/5I6)
 7135     CONTINUE
C       ENDIF
        ENDIF
      IF(IP.LE.208.AND.NSTART.EQ.1)THEN
        IF ((LFZC).AND.(IFINAL.EQ.0)) THEN
	 IF(IPEV.GE.1) WRITE(6,'(A)')'  before CHECKF'
          IF ((CMHIS.EQ.0.D0))
     +         CALL CHECKF(EPROJ,PPROJ,IREJ,1)
        ELSE
          IF ((CMHIS.EQ.0.D0))
     +         CALL CHECKO(EPROJ,PPROJ,IREJ,1)
        ENDIF
      ENDIF
      IF(IREJ.EQ.1)THEN
C     WRITE(6,'(A,I5)')' CHECKF/O IREJ ',IREJ
C         DO 4135 IHKK=1,NHKK
C           WRITE(6,1055) IHKK, ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),
C    +      JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK
C    +      (KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
C    +      ,IDRES(IHKK),IDXRES(IHKK),NOBAM(IHKK),
C    &                IDBAM(IHKK),IDCH(IHKK)
C4135     CONTINUE
        IF(KKCOUN.LE.1000)THEN
          WRITE(6,7734)KKCOUN
 7734     FORMAT(' KKCOUN=',I10)
        ENDIF
	IF(IPEV.GE.1)  WRITE(6,'(A)')'  after CHECKF'
        IF(IPRI.GE.1)THEN
          DO 7735 IHKK=1,NHKK
            WRITE(6,1055) IHKK, ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),
     +      JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK
     +      (KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
     +      , (WHKK(KHKK,IHKK),KHKK=1,4)
     +      ,IDRES(IHKK),IDXRES(IHKK),NOBAM(IHKK),
     &                IDBAM(IHKK),IDCH(IHKK)
 7735     CONTINUE
        ENDIF
        IREJ=0
        IF(KKCOUN.LE.500)THEN
          IF ((LFZC).AND.(IFINAL.EQ.0)) THEN
 	    WRITE(6,'(A)')' CHECKF Rejection'
          ELSE
 	    WRITE(6,'(A)')' CHECKO Rejection'
          ENDIF
        ENDIF
        GOTO 100
      ENDIF
      IF(NSTART.EQ.4.OR.NSTART.EQ.2)THEN
	 IF(IPEV.GE.1) WRITE(6,'(A)')'  before CHECKN'
          IF ((CMHIS.EQ.0.D0).AND.NEUDEC.NE.20)
     +         CALL CHECKN(EPROJ,PPROJ,IREJ,1)
          IF(KKCOUN.LE.500)THEN
 	    IF(IREJ.EQ.1)WRITE(6,'(A)')' CHECKN Rejection'
C	    IREJ=0
          ENDIF
      ENDIF

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&                      Writing file diffnuc2.evt for NSTART=3
C
      IF(NSTART.EQ.3.AND.IREJ.EQ.0)THEN
	KFORM=2
	IF(KFORM.EQ.1)THEN
	  AABBCC=0.
        ELSEIF(KFORM.EQ.2.AND.IREJ.EQ.0)THEN
          WRITE(33,'(I6,E12.4)')KJPRO,AMRECD
C                     The following only for 6 (J/psi)
	  READ(29,'(1X,I5,4E18.10)')IMIST,XXX1,XXX2,XXX3,XXX4
	  WRITE(33,'(1X,I5,4E18.10)')IMIST,XXX1,XXX2,XXX3,XXX4
	  READ(29,'(1X,I5,4E18.10)')IMIST,XXX1,XXX2,XXX3,XXX4
	  READ(29,'(1X,I5,4E18.10)')IMIST,XXX1,XXX2,XXX3,XXX4
C
	  READ(29,'(1X,I5)')KREPA
	  WRITE(33,'(1X,I5)')KREPA
	  DO 1977 KRE=1,KREPA
	    READ(29,'(1X,A)')A109
	    WRITE(33,'(1X,A)')A109
 1977     CONTINUE
	ENDIF
	  WRITE(33,*)' Event from dpmjet (only final particles):',
     *               'in Nucleus rest frame' 
          DO 1976 IHKK=1,NHKK
	  IF((ISTHKK(IHKK).EQ.-1).OR.
     *	  (ISTHKK(IHKK).EQ.1).OR.
     *	  (ISTHKK(IHKK).EQ.1001))THEN
            WRITE(33,'(2I6,5E18.10,2I6)')  ISTHKK(IHKK),IDHKK(IHKK),
     +      (PHKK(KHKK,IHKK),KHKK=1,5)
     +      ,IDRES(IHKK),IDXRES(IHKK)
	  ENDIF
 1976     CONTINUE
      ENDIF
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      IF(NSTART.EQ.1)THEN
	IF(IPEV.GE.1)  WRITE(6,'(A)')'  before CHEBCH '
        IF ((CMHIS.EQ.0.D0))THEN
           IF(IP.NE.IT.AND.IT.GT.1)   CALL CHEBCH(IREJ,NHKKH1)
           IF((IREJ.EQ.1))THEN
	     CHCOUN=CHCOUN+1
	     IF(CHCOUN.LE.50)THEN
               WRITE(6,'(A)')' CHEBCH Rejection'
               WRITE(6,'(A,I10)') ' KKINC: KKCOUN=',KKCOUN
	     ENDIF
 	     GOTO 100
 	   ENDIF
        ENDIF
	IF(IPEV.GE.1)WRITE(6,'(A)')'after CHEBCH before histograms'
      ENDIF
      IF(NEUDEC.EQ.20)CALL BACKDPM
	SUPX=0.D0
	SUPY=0.D0
	SUPZ=0.D0
        IF(KKCOUN.LE.50.AND.NSTART.GE.2)THEN
	  WRITE(6,*)' Event from dpmjet (only final particles):'
          DO 7736 IHKK=1,NHKK
  	  IF((ISTHKK(IHKK).EQ.-1).OR.
     *	  (ISTHKK(IHKK).EQ.1).OR.
     *	  (ISTHKK(IHKK).EQ.1001))THEN
	    SUPX=SUPX+PHKK(1,IHKK)
	    SUPY=SUPY+PHKK(2,IHKK)
	    SUPZ=SUPZ+PHKK(3,IHKK)
            WRITE(6,1055) IHKK, ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),
     +      JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK
     +      (KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
     +      , (WHKK(KHKK,IHKK),KHKK=1,4)
     +      ,IDRES(IHKK),IDXRES(IHKK),NOBAM(IHKK),
     &                IDBAM(IHKK),IDCH(IHKK)
	  ENDIF
 7736     CONTINUE
	  WRITE(6,*)' SUPX,SUPY,SUPZ ',SUPX,SUPY,SUPZ
        ENDIF
CGB
CGB     Output from G. Battistoni
CGB
      IF(NHKK.LE.0) THEN
         WRITE(6,*)' KKINC ', NHKK
         DO JGB = 1,NHKK
	   IF(ISTHKK(JGB).EQ.1001)THEN 
            WRITE(6,*)JGB, ISTHKK(JGB),IDHKK(JGB),
     *	   JMOHKK(1,JGB),JMOHKK(2,JGB),JDAHKK(1,JGB),JDAHKK(2,JGB),
     *      PHKK(1,JGB),PHKK(2,JGB)
     *           ,PHKK(3,JGB),PHKK(4,JGB),PHKK(5,JGB)
     +      ,IDRES(JGB),IDXRES(JGB),NOBAM(JGB),IDBAM(JGB),IDCH(JGB)
           ENDIF
         END DO
      ENDIF
CGB
C
      IF(NSTART.EQ.1)THEN
C               Random azimuthal rotation
	 CALL DSFECF(SFEE,CFEE)
         DO JGB = 1,NHKK
	   XXEE=PHKK(1,JGB)
	   YYEE=PHKK(2,JGB)
	   PHKK(1,JGB)=XXEE*CFEE-YYEE*SFEE
	   PHKK(2,JGB)=XXEE*SFEE+YYEE*CFEE
         END DO
      ENDIF
C     WRITE(6,'(A,I10)')' kkinc ',CMHIS
C     IF(XDIDI.GT.0.1D0)THEN
      IF (CMHIS.EQ.0.D0) CALL DISTR(2,NHKKH1,PPN,KTAUAC)
      IF (CMHIS.EQ.1.D0) CALL DISTRC(2,NHKKH1,PPN,KTAUAC)
      IF (CMHIS.EQ.2.D0) CALL DISTCO(2,NHKKH1,PPN,KTAUAC)
C     IF (IPRI.GE.2) CALL CHECKE(EPN,PPN)
C     ENDIF
C-----------
**
      RETURN
      END

**sr mod. for DPMJET: short version of the original DTUNUC-routine
*
*===defaux=============================================================*
*
      SUBROUTINE DEFAUX(EPN,PPN)

************************************************************************
* Variables are set to default values.                                 *
* This version dated 19.11.95 is written by S. Roesler.                *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)

      COMMON /NUCLEA/ PFERMP(2),PFERMN(2),FERMOD,
     &                EBINDP(2),EBINDN(2),EPOT(2,210),
     &                ETACOU(2),ICOUL
      LOGICAL LEMCCK,LHADRO,LSEADI
      COMMON /FLAGS/  IFRAG(2),IRESCO,IMSHL,IRESRJ,
     &                LEMCCK,LHADRO(0:9),LSEADI

      DATA POTMES /0.002D0/

* common /NUCLEA/
      DO 10 I=1,2
         PFERMP(I) = ZERO
         PFERMN(I) = ZERO
         EBINDP(I) = ZERO
         EBINDN(I) = ZERO
         DO 11 J=1,210
            EPOT(I,J) = ZERO
   11    CONTINUE
* nucleus independent meson potential
         EPOT(I,13) = POTMES
         EPOT(I,14) = POTMES
         EPOT(I,15) = POTMES
         EPOT(I,16) = POTMES
         EPOT(I,23) = POTMES
         EPOT(I,24) = POTMES
         EPOT(I,25) = POTMES
   10 CONTINUE
      FERMOD    = 0.95D0
      ETACOU(1) = ZERO
      ETACOU(2) = ZERO
      ICOUL     = 1

* common /FLAGS/
      IFRAG(1) = 2
      IFRAG(2) = 1
      IRESCO   = 1
      IMSHL    = 1
      IRESRJ   = 0
      LEMCCK   = .TRUE.
      LHADRO(0) = .FALSE.
      DO 13 I=1,9
         LHADRO(I) = .TRUE.
   13 CONTINUE
      LSEADI = .TRUE.

      RETURN
      END
*
*===nclpot=============================================================*
*
      SUBROUTINE NCLPOT(IPZ,IP,ITZ,IT,AFERP,AFERT,MODE)

************************************************************************
* Calculation of Coulomb and nuclear potential for a given configurat. *
*               IPZ, IP       charge/mass number of proj.              *
*               ITZ, IT       charge/mass number of targ.              *
*               AFERP,AFERT   factors modifying proj./target pot.      *
*                             if =0, FERMOD is used                    *
*               MODE = 0      calculation of binding energy            *
*                    = 1      pre-calculated binding energy is used    *
* This version dated 16.11.95  is written by S. Roesler.               *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TINY3=1.0D-3,TINY2=1.0D-2,
     &           TINY10=1.0D-10)

      LOGICAL LSTART

      CHARACTER*8  ANAME
      COMMON /DPAR/   ANAME(210),AAM(210),GA(210),TAU(210),
     &                IICH(210),IIBAR(210),K1(210),K2(210)

**sr mod. for DPMJET: use the longer DPMJET one
      LOGICAL INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA, IPADIS,
     &        ISHMAL,LPAULI
      COMMON /DROPPT/ INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA,
     &                IPADIS,ISHMAL,LPAULI
**
      COMMON /NUCLEA/ PFERMP(2),PFERMN(2),FERMOD,
     &                EBINDP(2),EBINDN(2),EPOT(2,210),
     &                ETACOU(2),ICOUL
**sr mod. for DPMJET: the corresponding common in DPMJET
      COMMON /NUCIMP/ PRMOM(5,248),TAMOM(5,248), PRMFEP,PRMFEN,TAMFEP,
     +TAMFEN, PREFEP,PREFEN,TAEFEP,TAEFEN, PREPOT(210),TAEPOT(210),
     +PREBIN,TAEBIN,FERFAC,ECOU
**

      DIMENSION IDXPOT(14)
*                   ap   an  lam  alam sig- sig+ sig0 tet0 tet- asig-
      DATA IDXPOT /   2,   9,  17,  18,  20,  21,  22,  97,  98,  99,
*                 asig0 asig+ atet0 atet+
     &              100, 101, 102, 103/

      DATA AN     /0.4D0/
      DATA LSTART /.TRUE./

      IF (MODE.EQ.0) THEN
         EBINDP(1) = ZERO
         EBINDN(1) = ZERO
         EBINDP(2) = ZERO
         EBINDN(2) = ZERO
      ENDIF
      AIP  = DBLE(IP)
      AIPZ = DBLE(IPZ)
      AIT  = DBLE(IT)
      AITZ = DBLE(ITZ)

      FERMIP = AFERP
      IF (AFERP.LE.ZERO) FERMIP = FERMOD
      FERMIT = AFERT
      IF (AFERT.LE.ZERO) FERMIT = FERMOD

* Fermi momenta and binding energy for projectile
      IF ((IP.GT.1).AND.(FERMP)) THEN
         IF (MODE.EQ.0) THEN
C           EBINDP(1) = EBIND(IP,IPZ)-EBIND(IP-1,IPZ-1)
C           EBINDN(1) = EBIND(IP,IPZ)-EBIND(IP-1,IPZ)
            BIP  = AIP -ONE
            BIPZ = AIPZ-ONE
            EBINDP(1) = 1.0D-3*ABS(ENERGY(AIP,AIPZ)-ENERGY(BIP,BIPZ))
            EBINDN(1) = 1.0D-3*ABS(ENERGY(AIP,AIPZ)-ENERGY(BIP,AIPZ))
         ENDIF
         PFERMP(1) = FERMIP*AN*(AIPZ/AIP)**0.333333D0
         PFERMN(1) = FERMIP*AN*((AIP-AIPZ)/AIP)**0.33333D0
      ELSE
         PFERMP(1) = ZERO
         PFERMN(1) = ZERO
      ENDIF
* effective nuclear potential for projectile
C     EPOT(1,1) = PFERMP(1)**2/(2.0D0*AAM(1)) + EBINDP(1)
C     EPOT(1,8) = PFERMN(1)**2/(2.0D0*AAM(8)) + EBINDN(1)
      EPOT(1,1) = SQRT(PFERMP(1)**2+AAM(1)**2) -AAM(1) + EBINDP(1)
      EPOT(1,8) = SQRT(PFERMN(1)**2+AAM(8)**2) -AAM(8) + EBINDN(1)

* Fermi momenta and binding energy for target
      IF ((IT.GT.1).AND.(FERMP)) THEN
         IF (MODE.EQ.0) THEN
C           EBINDP(2) = EBIND(IT,ITZ)-EBIND(IT-1,ITZ-1)
C           EBINDN(2) = EBIND(IT,ITZ)-EBIND(IT-1,ITZ)
            BIT  = AIT -ONE
            BITZ = AITZ-ONE
            EBINDP(2) = 1.0D-3*ABS(ENERGY(AIT,AITZ)-ENERGY(BIT,BITZ))
            EBINDN(2) = 1.0D-3*ABS(ENERGY(AIT,AITZ)-ENERGY(BIT,AITZ))
         ENDIF
         PFERMP(2) = FERMIT*AN*(AITZ/AIT)**0.333333D0
         PFERMN(2) = FERMIT*AN*((AIT-AITZ)/AIT)**0.33333D0
      ELSE
         PFERMP(2) = ZERO
         PFERMN(2) = ZERO
      ENDIF
* effective nuclear potential for target
C     EPOT(2,1) = PFERMP(2)**2/(2.0D0*AAM(1)) + EBINDP(2)
C     EPOT(2,8) = PFERMN(2)**2/(2.0D0*AAM(8)) + EBINDN(2)
      EPOT(2,1) = SQRT(PFERMP(2)**2+AAM(1)**2) -AAM(1) + EBINDP(2)
      EPOT(2,8) = SQRT(PFERMN(2)**2+AAM(8)**2) -AAM(8) + EBINDN(2)

      DO 2 I=1,14
         EPOT(1,IDXPOT(I)) = EPOT(1,8)
         EPOT(2,IDXPOT(I)) = EPOT(2,8)
    2 CONTINUE

* Coulomb energy
      ETACOU(1) = ZERO
      ETACOU(2) = ZERO
      IF (ICOUL.EQ.1) THEN
         IF (IP.GT.1)
     &   ETACOU(1) = 0.001116D0*AIPZ/(1.0D0+AIP**0.333D0)
         IF (IT.GT.1)
     &   ETACOU(2) = 0.001116D0*AITZ/(1.0D0+AIT**0.333D0)
      ENDIF

      IF (LSTART) THEN
         WRITE(LOUT,1000) IP,IPZ,IT,ITZ,EBINDP,EBINDN,
     &                    EPOT(1,1)-EBINDP(1),EPOT(2,1)-EBINDP(2),
     &                    EPOT(1,8)-EBINDN(1),EPOT(2,8)-EBINDN(2),
     &                    FERMOD,ETACOU
 1000    FORMAT(/,/,1X,'NCLPOT:    quantities for inclusion of nuclear'
     &           ,' effects',/,12X,'---------------------------',
     &           '----------------',/,/,38X,'projectile',
     &           '      target',/,/,1X,'Mass number / charge',
     &           17X,I3,' /',I3,6X,I3,' /',I3,/,1X,'Binding energy  -',
     &           ' proton   (GeV) ',2E14.4,/,17X,'- neutron  (GeV)'
     &          ,1X,2E14.4,/,1X,'Fermi-potential - proton   (GeV)',
     &           1X,2E14.4,/,17X,'- neutron  (GeV) ',2E14.4,/,/,
     &           1X,'Scale factor for Fermi-momentum    ',F4.2,/,
     &           /,1X,'Coulomb-energy ',2(E14.4,' GeV  '),/,/)
         LSTART = .FALSE.
      ENDIF

**sr mod. for DPMJET: fill /NUCIMP/
      PREBNN = ZERO
      PREBPN = ZERO
      PRMFEP = ZERO
      PRMFEN = ZERO
      IF ((IP.GT.1).AND.(FERMP)) THEN
         PREBNN = EBINDN(1)
         PREBPN = EBINDP(1)
         PRMFEP = PFERMP(1)
         PRMFEN = PFERMN(1)
      ENDIF
      PREFEN = PRMFEN**2/(2.*AAM(8))
      PREFEP = PRMFEP**2/(2.*AAM(1))
      PREPOT(1) = PREFEP + PREBPN
      PREPOT(8) = PREFEN + PREBNN
      TAEBNN = ZERO
      TAEBPN = ZERO
      TAMFEP = ZERO
      TAMFEN = ZERO
      IF ((IT.GT.1).AND.(FERMP)) THEN
         TAEBNN = EBINDN(2)
         TAEBPN = EBINDP(2)
         TAMFEP = PFERMP(2)
         TAMFEN = PFERMN(2)
      ENDIF
      TAEFEP = TAMFEP**2/(2.*AAM(1))
      TAEFEN = TAMFEN**2/(2.*AAM(8))
      TAEPOT(1) = TAEFEP + TAEBPN
      TAEPOT(8) = TAEFEN + TAEBNN
      DO 3 I=1,14
         TAEPOT(IDXPOT(I)) = TAEPOT(8)
    3 CONTINUE
      ECOU   = ETACOU(2)
      FERFAC = FERMOD
**

      RETURN
      END
*
*===resncl=============================================================*
*
      SUBROUTINE RESNCL(EPN,MODE)

************************************************************************
* Treatment of residual nuclei and nuclear effects.                    *
*         MODE = 1     initializations                                 *
*              = 2     treatment of final state                        *
* This version dated 16.11.95 is written by S. Roesler.                *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TINY3=1.0D-3,TINY2=1.0D-2,
     &           TINY1=1.0D-1,TINY4=1.0D-4,TINY10=1.0D-10)
      PARAMETER (AMUAMU=0.93149432D0)


      PARAMETER (NMXHKK=89998) 
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)
      COMMON /EXTEVT/ IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),
     &                IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10)
      CHARACTER*8  ANAME
      COMMON /DPAR/   ANAME(210),AAM(210),GA(210),TAU(210),
     &                IICH(210),IIBAR(210),K1(210),K2(210)

      LOGICAL LEMCCK,LHADRO,LSEADI
      COMMON /FLAGS/  IFRAG(2),IRESCO,IMSHL,IRESRJ,
     &                LEMCCK,LHADRO(0:9),LSEADI
      COMMON /NUCLEA/ PFERMP(2),PFERMN(2),FERMOD,
     &                EBINDP(2),EBINDN(2),EPOT(2,210),
     &                ETACOU(2),ICOUL
**sr mod. for DPMJET: use the longer DPMJET one
      LOGICAL INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA, IPADIS,
     &        ISHMAL,LPAULI
      COMMON /DROPPT/ INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA,
     &                IPADIS,ISHMAL,LPAULI
**
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
      COMMON /NUCCMS/ GACMS,BGCMS,GALAB,BGLAB,BLAB,UMO,PCM,EPROJ,PPROJ
      COMMON /WNDNCL/ NPW,NPW0,NPCW,NTW,NTW0,NTCW
      LOGICAL LRCLPR,LRCLTA
      COMMON /FINSTA/ PINIPR(5),PINITA(5),PRCLPR(5),PRCLTA(5)
     &, LRCLPR,LRCLTA
       COMMON /NSTARI/NSTART
      COMMON /NEUTYY/NEUTYP,NEUDEC
      DIMENSION PFSP(4),PSEC(4),PSEC0(4)

      GOTO (1,2) MODE

*------- initializations
    1 CONTINUE

* initialize arrays for residual nuclei
      DO 10 K=1,5
         IF (K.LE.4) THEN
            PFSP(K)     = ZERO
         ENDIF
         PINIPR(K) = ZERO
         PINITA(K) = ZERO
         PRCLPR(K) = ZERO
         PRCLTA(K) = ZERO
   10 CONTINUE

* projectile in n-n cms
      AIP  = DBLE(IP)
      AIPZ = DBLE(IPZ)
      PINIPR(4) = AIP*UMO/2.0D0
      PINIPR(5) = AIP*AMUAMU+1.0D-3*ENERGY(AIP,AIPZ)
      IF (IP.LE.1) PINIPR(5) = AAM(IJPROJ)
      PINIPR(3) = SQRT((PINIPR(4)-PINIPR(5))*(PINIPR(4)+PINIPR(5)))
C     WRITE(6,*)PINIPR,'PINIPR,1'
* target in n-n cms
      AIT  = DBLE(IT)
      AITZ = DBLE(ITZ)
      PINITA(4) = AIT*UMO/2.0D0
      PINITA(5) = AIT*AMUAMU+1.0D-3*ENERGY(AIT,AITZ)
C     WRITE(6,*)'UMO,PINITA(4),GACMS',UMO,PINITA(4),GACMS
      IF(PINITA(4).LE.PINITA(5))THEN
	PINITA(4)=GACMS*PINITA(5)
C	WRITE(6,*)'UMO,PINITA(4),GACMS',UMO,PINITA(4),GACMS
      ENDIF
      IF(NSTART.EQ.2)THEN
	PINITA(4)=GACMS*PINITA(5)
C	WRITE(6,*)'UMO,PINITA(4),GACMS',UMO,PINITA(4),GACMS
      ENDIF
      IF (IT.LE.1) PINITA(5) = AAM(IJTARG)
      PINITA(3) = -SQRT((PINITA(4)-PINITA(5))*(PINITA(4)+PINITA(5)))
C     WRITE(6,*)PINITA,'PINITA,1'

* correction of projectile 4-momentum for effective target pot.
* and Coulomb-energy (in case of hadron-nucleus interaction only)
      IF ((IP.EQ.1).AND.(IT.GT.1).AND.(FERMP)) THEN
         EPNI = EPN
*   Coulomb-energy:
*     positively charged hadron - check energy for Coloumb pot.
         IF (IICH(IJPROJ).EQ.1) THEN
            THRESH = ETACOU(2)+AAM(IJPROJ)
            IF (EPNI.LE.THRESH) THEN
               WRITE(LOUT,1000)
 1000          FORMAT(/,1X,'KKINC:  WARNING!  projectile energy',
     &                ' below Coulomb threshold - event rejected',/)
               ISTHKK(1) = 1
               RETURN
            ENDIF
*     negatively charged hadron - increase energy by Coulomb energy
         ELSEIF (IICH(IJPROJ).EQ.-1) THEN
            EPNI = EPNI+ETACOU(2)
         ENDIF
*   Effective target potential
C        EPNI = EPNI+EPOT(2,IJPROJ)
         EBIPOT = EBINDP(2)
         IF ((IJPROJ.NE.1).AND.(ABS(EPOT(2,IJPROJ)).GT.5.0D-3))
     &      EBIPOT = EBINDN(2)
         EPNI = EPNI+ABS(EBIPOT)
* re-initialization of NUCCMS
         DUM1 = ZERO
         DUM2 = ZERO
         IF(NSTART.NE.2.AND.NEUDEC.GE.20)
     &	 CALL LTINI(IJPROJ,EPNI,DUM1,DUM2)
C     COMMON /NEUTYY/NEUTYP,NEUDEC
      ENDIF

      RETURN

*------- treatment of final state
    2 CONTINUE

      JPW  = NPW
      JPCW = NPCW
      JTW  = NTW
      JTCW = NTCW

      DO 20 I=NPOINT(4),NHKK

         IDSEC  = IDBAM(I)

* reduction of particle momentum by corresponding nuclear potential
* (this applies only if Fermi-momenta are requested)

         IF (ISTHKK(I).EQ.1) THEN

C                           skip Photons
	 IF(IDSEC.EQ.7) GO TO 23

            IF (FERMP) THEN

*   select the nucleus which is most likely to be influenced by potential
*   corrections
               IPOT   = 0
               IOTHER = 0
               IF (PHKK(3,I).GE.ZERO) THEN
                  IPOT = 1
                  IF ((IP.LE.1).OR.((IP-NPW).LE.1)) THEN
                     IPOT   = 2
                     IF (IP.GT.1) IOTHER = 1
                     IF ((IT.LE.1).OR.((IT-NTW).LE.1)) GOTO 23
                  ENDIF
               ELSE
                  IPOT = 2
                  IF ((IT.LE.1).OR.((IT-NTW).LE.1)) THEN
                     IPOT   = 1
                     IF (IT.GT.1) IOTHER = 1
                     IF ((IP.LE.1).OR.((IP-NPW).LE.1)) GOTO 23
                  ENDIF
               ENDIF

*   Lorentz-transformation into the rest system of the selected nucleus
               IMODE = -IPOT-1
               CALL LTRANS(PHKK(1,I),PHKK(2,I),PHKK(3,I),PHKK(4,I),
     &                     PSEC(1),PSEC(2),PSEC(3),PSEC(4),IDSEC,IMODE)
               PSECO  = SQRT(PSEC(1)**2+PSEC(2)**2+PSEC(3)**2)
               AMSEC  = SQRT(ABS((PSEC(4)-PSECO)*(PSEC(4)+PSECO)))

               CHKLEV = TINY2
               IF ((EPROJ.GE.1.0D4).AND.(IDSEC.EQ.7)) CHKLEV = TINY1
               IF (EPROJ.GE.2.0D6) CHKLEV = 1.0D0
               IF (ABS(AMSEC-AAM(IDSEC)).GT.CHKLEV) THEN
C                 WRITE(LOUT,2000) I,NEVHKK,IDSEC,AMSEC,AAM(IDSEC)
 2000             FORMAT(1X,'RESNCL: inconsistent mass of particle',
     &                   ' at entry ',I5,' (evt.',I8,')',/,' IDSEC: ',
     &                   I4,'   AMSEC: ',E12.3,'  AAM(IDSEC): ',E12.3,/)
               ENDIF

               DO 21 K=1,4
                  PSEC0(K) = PSEC(K)
   21          CONTINUE

*   the correction for nuclear potential effects is applied to as many  
*   p/n as many nucleons were wounded; the momenta of other final state
*   particles are corrected only if they materialize inside the corresp.
*   nucleus (here: NOBAM = 1 part. outside proj., = 2 part. outside targ
*   = 3 part. outside proj. and targ., >=10 in overlapping region)
               IF ((IDSEC.EQ.1).OR.(IDSEC.EQ.8)) THEN
                  IF (IPOT.EQ.1) THEN
                     IF ((JPW.GT.0).AND.(IOTHER.EQ.0)) THEN
*      this is most likely a wounded nucleon
                        PSEC(4) = PSEC(4)-EPOT(IPOT,IDSEC)
                        JPW = JPW-1
                     ELSE
*      correct only if part. was materialized inside nucleus
*      and if it is ouside the overlapping region
                        IF ((NOBAM(I).NE.1).AND.(NOBAM(I).LT.3))
     &                     PSEC(4) = PSEC(4)-EPOT(IPOT,IDSEC)
                     ENDIF
                  ELSEIF (IPOT.EQ.2) THEN
                     IF ((JTW.GT.0).AND.(IOTHER.EQ.0)) THEN
*      this is most likely a wounded nucleon
                        PSEC(4) = PSEC(4)-EPOT(IPOT,IDSEC)
                        JTW = JTW-1
                     ELSE
*      correct only if part. was materialized inside nucleus
                        IF ((NOBAM(I).NE.2).AND.(NOBAM(I).LT.3))
     &                     PSEC(4) = PSEC(4)-EPOT(IPOT,IDSEC)
                     ENDIF
                  ENDIF
               ELSE
                  IF ((NOBAM(I).NE.IPOT).AND.(NOBAM(I).LT.3))
     &               PSEC(4) = PSEC(4)-EPOT(IPOT,IDSEC)
               ENDIF

* Coulomb energy correction:
* the treatment of Coulomb potential correction is similar to the
* one for nuclear potential
               IF (IDSEC.EQ.1) THEN
                  IF ((IPOT.EQ.1).AND.(JPCW.GT.0)) THEN
                     JPCW = JPCW-1
                  ELSEIF ((IPOT.EQ.2).AND.(JTCW.GT.0)) THEN
                     JTCW = JTCW-1
                  ELSE
                     IF ((NOBAM(I).EQ.IPOT).OR.(NOBAM(I).EQ.3)) GOTO 25
                  ENDIF
               ELSE
                  IF ((NOBAM(I).EQ.IPOT).OR.(NOBAM(I).EQ.3)) GOTO 25
               ENDIF
               IF (IICH(IDSEC).EQ.1) THEN
*    pos. particles: check if they are able to escape Coulomb potential
                  IF (PSEC(4).LT.AMSEC+ETACOU(IPOT)) THEN
                     ISTHKK(I) = 14+IPOT
                     IF (ISTHKK(I).EQ.15) THEN
                        DO 26 K=1,4
                           PHKK(K,I) = PSEC0(K)
                           PRCLPR(K) = PRCLPR(K)+PSEC0(K)
   26                   CONTINUE
                        IF ((IDSEC.EQ.1).OR.(IDSEC.EQ.8)) NPW = NPW-1
                        IF (IDSEC.EQ.1) NPCW = NPCW-1
                     ELSEIF (ISTHKK(I).EQ.16) THEN
                        DO 27 K=1,4
                           PHKK(K,I) = PSEC0(K)
                           PRCLTA(K) = PRCLTA(K)+PSEC0(K)
C			WRITE(6,*)I,K,PHKK(K,I),PRCLTA(K),'PRCLTA16+'
   27                   CONTINUE
                        IF ((IDSEC.EQ.1).OR.(IDSEC.EQ.8)) NTW = NTW-1
                        IF (IDSEC.EQ.1) NTCW = NTCW-1
                     ENDIF
                     GOTO 20
                  ENDIF
               ELSEIF (IICH(IDSEC).EQ.-1) THEN
*    neg. particles: decrease energy by Coulomb-potential
                  PSEC(4) = PSEC(4)-ETACOU(IPOT)
               ENDIF

   25          CONTINUE

               IF (PSEC(4).LT.AMSEC) THEN
C                 WRITE(LOUT,2001) I,IDSEC,PSEC(4),AMSEC
 2001             FORMAT(1X,'KKINC: particle at HKKEVT-pos. ',I5,
     &                   ' is not allowed to escape nucleus',/,
     &                   8X,'id : ',I3,'   reduced energy: ',E15.4,
     &                   '   mass: ',E12.3)
                  ISTHKK(I) = 14+IPOT
                  IF (ISTHKK(I).EQ.15) THEN
                     DO 28 K=1,4
                        PHKK(K,I) = PSEC0(K)
                        PRCLPR(K) = PRCLPR(K)+PSEC0(K)
   28                CONTINUE
                     IF ((IDSEC.EQ.1).OR.(IDSEC.EQ.8)) NPW = NPW-1
                     IF (IDSEC.EQ.1) NPCW = NPCW-1
                  ELSEIF (ISTHKK(I).EQ.16) THEN
                     DO 29 K=1,4
                        PHKK(K,I) = PSEC0(K)
                        PRCLTA(K) = PRCLTA(K)+PSEC0(K)
C			WRITE(6,*)I,K,PHKK(K,I),PRCLTA(K),'PRCLTA16+'
   29                CONTINUE
                     IF ((IDSEC.EQ.1).OR.(IDSEC.EQ.8)) NTW = NTW-1
                     IF (IDSEC.EQ.1) NTCW = NTCW-1
                  ENDIF
                  GOTO 20
               ENDIF

               PSECN  = SQRT( (PSEC(4)-AMSEC)*(PSEC(4)+AMSEC) )
* 4-momentum after correction for nuclear potential
               DO 22 K=1,3
                  PSEC(K) = PSEC(K)*PSECN/PSECO
   22          CONTINUE

* store recoil momentum from particles escaping the nuclear potentials
               DO 30 K=1,4
                  IF (IPOT.EQ.1) THEN
                     PRCLPR(K) = PRCLPR(K)+PSEC0(K)-PSEC(K)
                  ELSEIF (IPOT.EQ.2) THEN
                     PRCLTA(K) = PRCLTA(K)+PSEC0(K)-PSEC(K)
C			WRITE(6,*)I,K,PHKK(K,I),PRCLTA(K),'PRCLTA000'
                  ENDIF
   30          CONTINUE

* transform momentum back into n-n cms
               IMODE = IPOT+1
               CALL LTRANS(PSEC(1),PSEC(2),PSEC(3),PSEC(4),
     &                     PHKK(1,I),PHKK(2,I),PHKK(3,I),PHKK(4,I),
     &                     IDSEC,IMODE)

            ENDIF

   23       CONTINUE
            DO 31 K=1,4
               PFSP(K) = PFSP(K)+PHKK(K,I)
C	       WRITE(6,*)I,K,PHKK(K,I),PFSP(K),'PFSP,2'
   31       CONTINUE

         ENDIF
   20 CONTINUE
C           j.r.4.2.97
C     IF ((IP.EQ.1).AND.(IT.GT.1).AND.(FERMP)) THEN
      IF ((IP.EQ.10001).AND.(IT.GT.1).AND.(FERMP)) THEN
* hadron-nucleus interactions: get residual momentum from energy-
* momentum conservation
         DO 32 K=1,4
            PRCLPR(K) = ZERO
            PRCLTA(K) = PINIPR(K)+PINITA(K)-PFSP(K)
C			WRITE(6,*)I,K,PHKK(K,I),PRCLTA(K),'PRCLTA111'
C	WRITE(6,*)K,PINIPR(K),PINITA(K),PFSP(K),PRCLTA(K),'PRCLTA222'
   32    CONTINUE
      ELSE
* nucleus-hadron, nucleus-nucleus: get residual momentum from 
* accumulated recoil momenta of particles leaving the spectators
*   transform accumulated recoil momenta of residual nuclei into 
*   n-n cms
         PZI = PRCLPR(3)
         PEI = PRCLPR(4)
         CALL LTNUC(PZI,PEI,PRCLPR(3),PRCLPR(4),2)
         PZI = PRCLTA(3)
         PEI = PRCLTA(4)
         CALL LTNUC(PZI,PEI,PRCLTA(3),PRCLTA(4),3)
C        IF (IP.GT.1) THEN
            PRCLPR(3) = PRCLPR(3)+PINIPR(3)
            PRCLPR(4) = PRCLPR(4)+PINIPR(4)
C        ENDIF
         IF (IT.GT.1) THEN
	    KKK=3
C	 WRITE(6,*)KKK,PINITA(3),PRCLTA(KKK),'PRCLTAkkk'
	    KKK=4
C	 WRITE(6,*)KKK,PINITA(4),PRCLTA(KKK),'PRCLTAkkk'
            PRCLTA(3) = PRCLTA(3)+PINITA(3)
	    KKK=3
C       	WRITE(6,*)KKK,PINITA(3),PRCLTA(KKK),'PRCLTAkkk'
            PRCLTA(4) = PRCLTA(4)+PINITA(4)
	    KKK=4
C       	WRITE(6,*)KKK,PINITA(4),PRCLTA(KKK),'PRCLTAkkk'
         ENDIF
      ENDIF

* check momenta of residual nuclei
      IF (LEMCCK) THEN
         CALL EVTEMC(-PINIPR(1),-PINIPR(2),-PINIPR(3),-PINIPR(4),
     &               1,IDUM,IDUM)
         CALL EVTEMC(-PINITA(1),-PINITA(2),-PINITA(3),-PINITA(4),
     &               2,IDUM,IDUM)
         CALL EVTEMC(PRCLPR(1),PRCLPR(2),PRCLPR(3),PRCLPR(4),
     &               2,IDUM,IDUM)
         CALL EVTEMC(PRCLTA(1),PRCLTA(2),PRCLTA(3),PRCLTA(4),
     &               2,IDUM,IDUM)
         CALL EVTEMC(PFSP(1),PFSP(2),PFSP(3),PFSP(4),2,IDUM,IDUM)
         CHKLEV = TINY3
         CALL EVTEMC(DUM,DUM,DUM,CHKLEV,-1,501,IREJ1)
         IF (IREJ1.GT.0) RETURN
      ENDIF

      RETURN
      END
*
*
*===scn4ba=============================================================*
*
      SUBROUTINE SCN4BA

************************************************************************
* SCan /HKKEVT/ 4 BAryons which are not able to escape nuclear pot.    *
* This version dated 12.12.95 is written by S. Roesler.                *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TINY3=1.0D-3,TINY2=1.0D-2,
     &           TINY10=1.0D-10)

      PARAMETER (NMXHKK=89998)
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)
      COMMON /EXTEVT/ IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),
     &                IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10)
      CHARACTER*8  ANAME
      COMMON /DPAR/   ANAME(210),AAM(210),GA(210),TAU(210),
     &                IICH(210),IIBAR(210),K1(210),K2(210)

      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
      COMMON /NUCLEA/ PFERMP(2),PFERMN(2),FERMOD,
     &                EBINDP(2),EBINDN(2),EPOT(2,210),
     &                ETACOU(2),ICOUL
      COMMON /WNDNCL/ NPW,NPW0,NPCW,NTW,NTW0,NTCW
      LOGICAL LRCLPR,LRCLTA
      COMMON /FINSTA/ PINIPR(5),PINITA(5),PRCLPR(5),PRCLTA(5),
     &                LRCLPR,LRCLTA

      DIMENSION PLAB(2,5),PCMS(4)

      IREJ = 0

* get number of wounded nucleons
      NPW    = 0
      NPW0   = 0
      NPCW   = 0
      NPSTCK = 0
      NTW    = 0
      NTW0   = 0
      NTCW   = 0
      NTSTCK = 0

      ISGLPR = 0
      ISGLTA = 0
      LRCLPR = .FALSE.
      LRCLTA = .FALSE.

C     DO 2 I=1,NHKK
      DO 2 I=1,NPOINT(1)
* projectile nucleons wounded in primary interaction and in fzc
         IF ((ISTHKK(I).EQ.11).OR.(ISTHKK(I).EQ.17)) THEN
            NPW    = NPW+1
            NPSTCK = NPSTCK+1
            IF (IDHKK(I).EQ.2212) NPCW = NPCW+1
            IF (ISTHKK(I).EQ.11)  NPW0 = NPW0+1
C           IF (IP.GT.1) THEN
               DO 5 K=1,4
                  PRCLPR(K) = PRCLPR(K)-PHKK(K,I)
    5          CONTINUE
C           ENDIF
* target nucleons wounded in primary interaction and in fzc
         ELSEIF ((ISTHKK(I).EQ.12).OR.(ISTHKK(I).EQ.18)) THEN
            NTW    = NTW+1
            NTSTCK = NTSTCK+1
            IF (IDHKK(I).EQ.2212) NTCW = NTCW+1
            IF (ISTHKK(I).EQ.12)  NTW0 = NTW0+1
            IF (IT.GT.1) THEN
               DO 6 K=1,4
                  PRCLTA(K) = PRCLTA(K)-PHKK(K,I)
C			WRITE(6,*)I,K,PHKK(K,I),PRCLTA(K),'PRCLTA12-'
    6          CONTINUE
            ENDIF
         ELSEIF (ISTHKK(I).EQ.13) THEN
            ISGLPR = I
         ELSEIF (ISTHKK(I).EQ.14) THEN
            ISGLTA = I
         ENDIF
    2 CONTINUE

      DO 11 I=NPOINT(4),NHKK
* baryons which are unable to escape the nuclear potential of proj.
         IF (ISTHKK(I).EQ.15) THEN
            ISGLPR = I
            NPSTCK = NPSTCK-1
            IF (IIBAR(IDBAM(I)).NE.0) THEN
               NPW    = NPW-1
               IF (IICH(IDBAM(I)).GT.0) NPCW = NPCW-1
            ENDIF
            DO 7 K=1,4
               PRCLPR(K) = PRCLPR(K)+PHKK(K,I)
    7       CONTINUE
* baryons which are unable to escape the nuclear potential of targ.
         ELSEIF (ISTHKK(I).EQ.16) THEN
            ISGLTA = I
            NTSTCK = NTSTCK-1
            IF (IIBAR(IDBAM(I)).NE.0) THEN
               NTW    = NTW-1
               IF (IICH(IDBAM(I)).GT.0) NTCW = NTCW-1
            ENDIF
            DO 8 K=1,4
               PRCLTA(K) = PRCLTA(K)+PHKK(K,I)
C			WRITE(6,*)I,K,PHKK(K,I),PRCLTA(K),'PRCLTA16+'
    8       CONTINUE
         ENDIF
   11 CONTINUE

* residual nuclei so far
      IRESP = IP-NPSTCK
      IREST = IT-NTSTCK

* ckeck for "residual nuclei" consisting of one nucleon only
* treat it as final state particle
      IF (IRESP.EQ.1) THEN
         ID  = IDBAM(ISGLPR)
         IST = ISTHKK(ISGLPR)
         CALL LTRANS(PHKK(1,ISGLPR),PHKK(2,ISGLPR),
     &               PHKK(3,ISGLPR),PHKK(4,ISGLPR),
     &               PCMS(1),PCMS(2),PCMS(3),PCMS(4),ID,2)
         IF (IST.EQ.13) THEN
            ISTHKK(ISGLPR) = 11
         ELSE
            ISTHKK(ISGLPR) = 2
         ENDIF
         CALL EVTPUT(1,IDHKK(ISGLPR),ISGLPR,0,
     &               PCMS(1),PCMS(2),PCMS(3),PCMS(4),
     &               IDRES(ISGLPR),IDXRES(ISGLPR),IDCH(ISGLPR))
         NOBAM(NHKK)      = NOBAM(ISGLPR)
         JDAHKK(1,ISGLPR) = NHKK
         DO 21 K=1,4
            PRCLPR(K) = PRCLPR(K)-PHKK(K,ISGLPR)
   21    CONTINUE
      ENDIF
      IF (IREST.EQ.1) THEN
         ID  = IDBAM(ISGLTA)
         IST = ISTHKK(ISGLTA)
         CALL LTRANS(PHKK(1,ISGLTA),PHKK(2,ISGLTA),
     &               PHKK(3,ISGLTA),PHKK(4,ISGLTA),
     &               PCMS(1),PCMS(2),PCMS(3),PCMS(4),ID,3)
         IF (IST.EQ.14) THEN
            ISTHKK(ISGLTA) = 12
         ELSE
            ISTHKK(ISGLTA) = 2
         ENDIF
         CALL EVTPUT(1,IDHKK(ISGLTA),ISGLTA,0,
     &               PCMS(1),PCMS(2),PCMS(3),PCMS(4),
     &               IDRES(ISGLTA),IDXRES(ISGLTA),IDCH(ISGLTA))
         NOBAM(NHKK)      = NOBAM(ISGLTA)
         JDAHKK(1,ISGLTA) = NHKK
         DO 22 K=1,4
            PRCLTA(K) = PRCLTA(K)-PHKK(K,ISGLTA)
C		WRITE(6,*)ISGLTA,K,PHKK(K,ISGLTA),PRCLTA(K),'PRCLTA12-'
   22    CONTINUE
      ENDIF

* get nuclear potential corresp. to the residual nucleus
      IPRCL  = IP -NPW
      IPZRCL = IPZ-NPCW
      ITRCL  = IT -NTW
      ITZRCL = ITZ-NTCW
      CALL NCLPOT(IPZRCL,IPRCL,ITZRCL,ITRCL,ZERO,ZERO,1)

* baryons unable to escape the nuclear potential are treated as
* excited nucleons (ISTHKK=15,16)
      DO 3 I=NPOINT(4),NHKK
         IF (ISTHKK(I).EQ.1) THEN
            ID  = IDBAM(I)
            IF ( ((ID.EQ.1).OR.(ID.EQ.8)).AND.(NOBAM(I).NE.3) ) THEN
*   final state n and p not being outside of both nuclei are considered
               NPOTP = 1
               NPOTT = 1
               IF ( (IP.GT.1)      .AND.(IRESP.GT.1).AND.
     &              (NOBAM(I).NE.1).AND.(NPW.GT.0)        ) THEN
*     Lorentz-trsf. into proj. rest sys. for those being inside proj.
                  CALL LTRANS(PHKK(1,I),PHKK(2,I),PHKK(3,I),PHKK(4,I),
     &                        PLAB(1,1),PLAB(1,2),PLAB(1,3),PLAB(1,4),
     &                                                          ID,-2)
                  PLABT = SQRT(PLAB(1,1)**2+PLAB(1,2)**2+PLAB(1,3)**2)
                  PLAB(1,5) = SQRT(ABS( (PLAB(1,4)-PLABT)*
     &                                  (PLAB(1,4)+PLABT) ))
                  EKIN = PLAB(1,4)-PLAB(1,5)
                  IF (EKIN.LE.EPOT(1,ID)) NPOTP = 15
                  IF ((ID.EQ.1).AND.(NPCW.LE.0)) NPOTP = 1
               ENDIF
               IF ( (IT.GT.1)      .AND.(IREST.GT.1).AND.
     &              (NOBAM(I).NE.2).AND.(NTW.GT.0)        ) THEN
*     Lorentz-trsf. into targ. rest sys. for those being inside targ.
                  CALL LTRANS(PHKK(1,I),PHKK(2,I),PHKK(3,I),PHKK(4,I),
     &                        PLAB(2,1),PLAB(2,2),PLAB(2,3),PLAB(2,4),
     &                                                          ID,-3)
                  PLABT = SQRT(PLAB(2,1)**2+PLAB(2,2)**2+PLAB(2,3)**2)
                  PLAB(2,5) = SQRT(ABS( (PLAB(2,4)-PLABT)*
     &                                  (PLAB(2,4)+PLABT) ))
                  EKIN = PLAB(2,4)-PLAB(2,5)
                  IF (EKIN.LE.EPOT(2,ID)) NPOTT = 16
                  IF ((ID.EQ.1).AND.(NTCW.LE.0)) NPOTT = 1
               ENDIF
               IF (PHKK(3,I).GE.ZERO) THEN
                  ISTHKK(I) = NPOTT
                  IF (NPOTP.NE.1) ISTHKK(I) = NPOTP
               ELSE
                  ISTHKK(I) = NPOTP
                  IF (NPOTT.NE.1) ISTHKK(I) = NPOTT
               ENDIF
               IF (ISTHKK(I).NE.1) THEN
                  J = ISTHKK(I)-14
                  DO 4 K=1,5
                     PHKK(K,I) = PLAB(J,K)
    4             CONTINUE
                  IF (ISTHKK(I).EQ.15) THEN
                     NPW = NPW-1
                     IF (ID.EQ.1) NPCW = NPCW-1
                     DO 9 K=1,4
                        PRCLPR(K) = PRCLPR(K)+PHKK(K,I)
C			WRITE(6,*)I,K,PHKK(K,I),PRCLPR(K),'PRCLPR'
    9                CONTINUE
                  ELSEIF (ISTHKK(I).EQ.16) THEN
                     NTW = NTW-1
                     IF (ID.EQ.1) NTCW = NTCW-1
                     DO 10 K=1,4
                        PRCLTA(K) = PRCLTA(K)+PHKK(K,I)
C			WRITE(6,*)I,K,PHKK(K,I),PRCLTA(K),'PRCLTA16+'
   10                CONTINUE
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
    3 CONTINUE

* again: get nuclear potential corresp. to the residual nucleus
      IPRCL  = IP -NPW
      IPZRCL = IPZ-NPCW
      ITRCL  = IT -NTW
      ITZRCL = ITZ-NTCW
      AFERP  = FERMOD+0.1D0
      AFERT  = FERMOD+0.1D0
      CALL NCLPOT(IPZRCL,IPRCL,ITZRCL,ITRCL,AFERP,AFERT,1)


      RETURN
      END
*
*
*===ficonf=============================================================*
*
      SUBROUTINE FICONF(IJPROJ,IP,IPZ,IT,ITZ,IREJ)

************************************************************************
* Treatment of FInal CONFiguration including evaporation, fission and  *
* Fermi-break-up (for light nuclei only).                              *
* Adopted from the original routine FINALE and extended to residual    *
* projectile nuclei.                                                   *
* This version dated 12.12.95 is written by S. Roesler.                *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TINY3=1.0D-3,TINY10=1.0D-10)

      PARAMETER (NMXHKK=89998)
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)
      COMMON /EXTEVT/ IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),
     &                IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10)
      COMMON /RJCOUN/ IRPT,IRHHA,IRRES(2),LOMRES,LOBRES,
     &                IRCHKI(2),IRFRAG,IRCRON(3),IREVT,
     &                IREXCI(3),IRDIFF(2),IRINC
      COMMON /ZENTRA/ ICENTR
      COMMON /OUTLEV/IOUTPO,IOUTPA,IOUXEV,IOUCOL
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
      CHARACTER*8  ANAME
      COMMON /DPAR/   ANAME(210),AAM(210),GA(210),TAU(210),
     &                IICH(210),IIBAR(210),K1(210),K2(210)
      LOGICAL LRCLPR,LRCLTA
      COMMON /FINSTA/ PINIPR(5),PINITA(5),PRCLPR(5),PRCLTA(5),
     &                LRCLPR,LRCLTA
      COMMON /EXCITA/ AMRCL0(2),EEXC(2),EEXCFI(2),
     &                NTOT(2),NPRO(2),NN(2),NH(2),NHPOS(2),NQ(2),
     &                NTOTFI(2),NPROFI(2)
      COMMON /STFICO/ EXCDPM(4),EXCEVA(2),
     &                NINCGE,NINCCO(2,3),NINCHR(2,2),NINCWO(2),
     &                NINCST(2,4),NINCEV(2),
     &                NRESTO(2),NRESPR(2),NRESNU(2),NRESBA(2),
     &                NRESPB(2),NRESCH(2),NRESEV(4),
     &                NEVA(2,6),NEVAGA(2),NEVAHT(2),NEVAHY(2,2,240),
     &                NEVAFI(2,2)
* evaporation interface
      PARAMETER (ANGLGB=5.0D-16)
      PARAMETER (AMUAMU=0.93149432D0,AMELEC=0.51099906D-3)
      PARAMETER (MXP=999)
      COMMON /FINUC/  CXR (MXP), CYR (MXP), CZR (MXP), TKI (MXP),
     &                PLR (MXP), WEI (MXP), TV, TVCMS, TVRECL, TVHEAV,
     &                TVBIND, NP0, NP, KPART (MXP)
      LOGICAL LRNFSS, LFRAGM
      COMMON /RESNUC/  AMNTAR, AMMTAR, AMNZM1, AMMZM1, AMNNM1, AMMNM1,
     &                   ANOW,   ZNOW, ANCOLL, ZNCOLL, AMMLFT, AMNLFT,
     &                   ERES,  EKRES, AMNRES, AMMRES,  PTRES,  PXRES,
     &                  PYRES,  PZRES, PTRES2,  KTARP,  KTARN, IGREYP,
     &                 IGREYN,  ICRES,  IBRES, ISTRES, IEVAPL, IEVAPH,
     &                 IEVNEU, IEVPRO, IEVDEU, IEVTRI, IEV3HE, IEV4HE,
     &                 IDEEXG,  IBTAR, ICHTAR, IBLEFT, ICLEFT, IOTHER,
     &                 LRNFSS, LFRAGM
      COMMON /NUCDAT/ AV0WEL,     APFRMX,     AEFRMX,     AEFRMA,
     &                RDSNUC,     V0WELL (2), PFRMMX (2), EFRMMX (2),
     &                EFRMAV (2), AMNUCL (2), AMNUSQ (2), EBNDNG (2),
     &                VEFFNU (2), ESLOPE (2), PKMNNU (2), EKMNNU (2),
     &                PKMXNU (2), EKMXNU (2), EKMNAV (2), EKINAV (2),
     &                EXMNAV (2), EKUPNU (2), EXMNNU (2), EXUPNU (2),
     &                ERCLAV (2), ESWELL (2), FINCUP (2), AMRCAV    ,
     &                AMRCSQ    , ATO1O3    , ZTO1O3    , ELBNDE (0:100)
      LOGICAL LDIFFR, LINCTV, LEVPRT, LHEAVY, LDEEXG, LGDHPR, LPREEX,
     &        LHLFIX, LPRFIX, LPARWV, LPOWER, LSNGCH, LLVMOD, LSCHDF
      PARAMETER ( NALLWP = 39   )
      COMMON / PAREVT / DPOWER, FSPRD0, FSHPFN, RN1GSC, RN2GSC,
     &                  LDIFFR (NALLWP),LPOWER, LINCTV, LEVPRT, LHEAVY,
     &                  LDEEXG, LGDHPR, LPREEX, LHLFIX, LPRFIX, LPARWV,
     &                  ILVMOD, JLVMOD, LLVMOD, LSNGCH, LSCHDF

      DIMENSION INUC(2),IDXPAR(2),IDPAR(2),AIF(2),AIZF(2),AMRCL(2),
     &          PRCL(2,4),MO1(2),MO2(2),VRCL(2,4),WRCL(2,4),
     &          P1IN(4),P2IN(4),P1OUT(4),P2OUT(4)
      COMMON/REJFBK/IREJFR
      COMMON /NEUTYY/NEUTYP,NEUDEC

      DIMENSION EXPNUC(2),EXC(2,210),NEXC(2,210)
      DATA EXC,NEXC /420*ZERO,420*0/
      DATA EXPNUC /4.0D-3,4.0D-3/
      DATA INIEX/0/
      DATA INIWA/0/

      IREJ   = 0
      LRCLPR = .FALSE.
      LRCLTA = .FALSE.

* skip residual nucleus treatment if not requested or in case
* of central collisions
      IF(IPEV.GE.1)WRITE(6,*)' FICONF: LEVPRT ICENTR',LEVPRT,ICENTR
C     IF ((.NOT.LEVPRT).OR.(ICENTR.NE.0)) RETURN
C                          jr.19.5.96 also for central coll.
      IF ((.NOT.LEVPRT)) RETURN

      DO 1 K=1,2
         IDPAR(K) = 0
         IDXPAR(K)= 0
         NTOT(K)  = 0
         NTOTFI(K)= 0
         NPRO(K)  = 0
         NPROFI(K)= 0
         NN(K)    = 0
         NH(K)    = 0
         NHPOS(K) = 0
         NQ(K)    = 0
         EEXC(K)  = ZERO
         MO1(K)   = 0
         MO2(K)   = 0
         DO 2 I=1,4
            VRCL(K,I) = ZERO
            WRCL(K,I) = ZERO
    2    CONTINUE
    1 CONTINUE
      NFSP = 0
      INUC(1) = IP
      INUC(2) = IT

      DO 3 I=1,NHKK

* number of final state particles
         IF (ABS(ISTHKK(I)).EQ.1) THEN
            NFSP  = NFSP+1
            IDFSP = IDBAM(I)
         ENDIF

* properties of remaining nucleon configurations
         KF = 0
         IF ((ISTHKK(I).EQ.13).OR.(ISTHKK(I).EQ.15)) KF = 1
         IF ((ISTHKK(I).EQ.14).OR.(ISTHKK(I).EQ.16)) KF = 2
         IF (KF.GT.0) THEN
            IF (MO1(KF).EQ.0) MO1(KF) = I
            MO2(KF)  = I
*   position of residual nucleus = average position of nucleons
            DO 4 K=1,4
               VRCL(KF,K) = VRCL(KF,K)+VHKK(K,I)
               WRCL(KF,K) = WRCL(KF,K)+WHKK(K,I)
    4       CONTINUE
*   total number of particles contributing to each residual nucleus
            NTOT(KF)  = NTOT(KF)+1
            IDTMP     = IDBAM(I)
            IDXTMP    = I
*   total charge of residual nuclei
            NQ(KF) = NQ(KF)+IICH(IDTMP)
*   number of protons
            IF (IDHKK(I).EQ.2212) THEN
               NPRO(KF) = NPRO(KF)+1
*   number of neutrons
            ELSEIF (IDHKK(I).EQ.2112) THEN
               NN(KF) = NN(KF)+1
            ELSE
*   number of baryons other than n, p
               IF (IIBAR(IDTMP).EQ.1) THEN
                  NH(KF) = NH(KF)+1
                  IF (IICH(IDTMP).EQ.1) NHPOS(KF) = NHPOS(KF)+1
               ELSE
*   any other mesons (status set to 1)
                  INIWA=INIWA+1
                  IF(INIWA.LE.20)WRITE(LOUT,1002) KF,IDTMP
 1002             FORMAT(1X,'FICONF:   residual nucleus ',I2,
     &                   ' containing meson ',I4,', status set to 1')
                  ISTHKK(I) = 1
                  IDTMP     = IDPAR(KF)
                  IDXTMP    = IDXPAR(KF)
                  NTOT(KF)  = NTOT(KF)-1
               ENDIF
            ENDIF
            IDPAR(KF)  = IDTMP
            IDXPAR(KF) = IDXTMP
         ENDIF
    3 CONTINUE

* reject elastic events (def: one final state particle = projectile)
      IF ((IP.EQ.1).AND.(NFSP.EQ.1).AND.(IDFSP.EQ.IJPROJ)) THEN
                  WRITE(LOUT,1009) 
 1009             FORMAT(1X,'FICONF: ct elastic events ')
         IREXCI(3) = IREXCI(3)+1
	 IREJ=1
         RETURN
      ENDIF

* check if one nucleus disappeared..
C     IF ((IP.GT.1).AND.(NTOT(1).EQ.0).AND.(NTOT(2).NE.0)) THEN
C        DO 5 K=1,4
C           PRCLTA(K) = PRCLTA(K)+PRCLPR(K)
C           PRCLPR(K) = ZERO
C   5    CONTINUE
C     ELSEIF ((IT.GT.1).AND.(NTOT(2).EQ.0).AND.(NTOT(1).NE.0)) THEN
C        DO 6 K=1,4
C           PRCLPR(K) = PRCLPR(K)+PRCLTA(K)
C           PRCLTA(K) = ZERO
C   6    CONTINUE
C     ENDIF

      ICOR   = 0
      INORCL = 0
      DO 7 I=1,2
         DO 8 K=1,4
* get the average of the nucleon positions
            VRCL(I,K) = VRCL(I,K)/MAX(NTOT(I),1)
            WRCL(I,K) = WRCL(I,K)/MAX(NTOT(I),1)
            IF (I.EQ.1) PRCL(1,K) = PRCLPR(K)
            IF (I.EQ.2) PRCL(2,K) = PRCLTA(K)
    8    CONTINUE
       IF(IPEV.GE.1)WRITE(6,*)PRCL,'PRCL(2,4)'
       IF(IPEV.GE.1)WRITE(6,*)PRCLTA,'PRCLTA'
* mass number and charge of residual nuclei
         AIF(I)  = DBLE(NTOT(I))
         AIZF(I) = DBLE(NPRO(I)+NHPOS(I))
         IF(IPEV.GE.1)WRITE(6,*)'I,Ntot(i)',I,NTOT(I),AIF(I),AIZF(I)
         IF (NTOT(I).GT.1) THEN
* masses of residual nuclei in ground state
            AMRCL0(I) = AIF(I)*AMUAMU+1.0D-3*ENERGY(AIF(I),AIZF(I))
* masses of residual nuclei
            PTORCL   = SQRT(PRCL(I,1)**2+PRCL(I,2)**2+PRCL(I,3)**2)
            AMRCL(I) = (PRCL(I,4)-PTORCL)*(PRCL(I,4)+PTORCL)
            IF (AMRCL(I).GT.ZERO) AMRCL(I) = SQRT(AMRCL(I))
 	   IF(IPEV.GE.1) WRITE(6,*)AMRCL(I),'AMRCL(',I,')'
C                      Patch 5.2.98
            IF ((AMRCL(I).LT.AMRCL0(I)).AND.(NEUDEC.EQ.20))
     &	    AMRCL(I)=AMRCL0(I)+0.025D0
            IF (AMRCL(I).LE.ZERO) THEN
	       INIEX=INIEX+1
	       IF(INIEX.LE.50)
     &         WRITE(LOUT,1000) I,PRCL(I,1),PRCL(I,2),PRCL(I,3),
     &                          PRCL(I,4),AMRCL(I),NTOT
 1000          FORMAT(1X,'warning! negative excitation energy',/,
     &                I4,5E15.4,2I4)
               AMRCL(I) = ZERO
               EEXC(I)  = ZERO
               GOTO 9999
            ELSEIF ((AMRCL(I).GT.ZERO).AND.(AMRCL(I).LT.AMRCL0(I)))
     &                                                         THEN
               EEXC(I)  = AMRCL(I)-AMRCL0(I)
C	       WRITE(6,*)I,EEXC(I),AMRCL(I),AMRCL0(I),'EEXC(I)0'
**sr 11.6.96
C              AMRCL(I) = AMRCL0(I)+EXPNUC(I)*DBLE(NTOT(I))
               M = MIN(NTOT(I),210)
               IF (NEXC(I,M).GT.0) THEN
                  AMRCL(I) = AMRCL0(I)+EXC(I,M)/DBLE(NEXC(I,M))
               ELSE
   70             CONTINUE
                  M = M+1
                  IF (M.GE.INUC(I)) THEN
                     AMRCL(I) = AMRCL0(I)+EXPNUC(I)*DBLE(NTOT(I))
                  ELSE
                     IF (NEXC(I,M).GT.0) THEN
                        AMRCL(I) = AMRCL0(I)+EXC(I,M)/DBLE(NEXC(I,M))
                     ELSE
                        GOTO 70
                     ENDIF
                  ENDIF
               ENDIF
**
               EEXC(I)  = AMRCL(I)-AMRCL0(I)
	       IF(IPEV.GE.1)THEN
 	       WRITE(6,*)I,EEXC(I),AMRCL(I),AMRCL0(I),'EEXC(I)1'
	       ENDIF
               IF ((AMRCL(I).GT.ZERO).AND.(AMRCL(I).LT.AMRCL0(I)))
     &                                                         THEN
                 ICOR     = ICOR+I
               ENDIF
C                  insert 4.2.98
               EXPNUC(I) = EEXC(I)/MAX(1,INUC(I)-NTOT(I))
**sr 11.6.96
               M = MIN(NTOT(I),210)
               EXC(I,M)  = EXC(I,M)+EEXC(I)
               NEXC(I,M) = NEXC(I,M)+1
C                  insert 4.2.98
            ELSE
* excitation energies of residual nuclei
               EEXC(I)   = AMRCL(I)-AMRCL0(I)
	       IF(IPEV.GE.1)THEN
 	       WRITE(6,*)I,EEXC(I),AMRCL(I),AMRCL0(I),'EEXC(I)2'
	       ENDIF
               EXPNUC(I) = EEXC(I)/MAX(1,INUC(I)-NTOT(I))
**sr 11.6.96
               M = MIN(NTOT(I),210)
               EXC(I,M)  = EXC(I,M)+EEXC(I)
               NEXC(I,M) = NEXC(I,M)+1
**
            ENDIF
         ELSEIF (NTOT(I).EQ.1) THEN
            WRITE(LOUT,1003) I
 1003       FORMAT(1X,'FICONF:   warning! NTOT(I)=1? (I=',I3,')')
            GOTO 9999
         ELSE
            AMRCL0(I) = ZERO
            AMRCL(I)  = ZERO
            EEXC(I)   = ZERO
            INORCL    = INORCL+I
	    IF(IPEV.GE.1)WRITE(6,*)' INORCL,I',INORCL,I
         ENDIF
	       IF(IPEV.GE.1)THEN
 	 WRITE (6,'(A,I10,3F10.3)')' I,AIF,AIZF,EEXC:'
     *,I,AIF(I),AIZF(I),EEXC(I)
         ENDIF
    7 CONTINUE

      PRCLPR(5) = AMRCL(1)
      PRCLTA(5) = AMRCL(2)
      IF(IPEV.GE.1)WRITE(6,*)' ICOR,INORCL ',ICOR,INORCL
      IF (ICOR.GT.0) THEN
        IF (INORCL.EQ.0) THEN
* one or both residual nuclei consist of one nucleon only, transform
* this nucleon on mass shell
          DO 9 K=1,4
            P1IN(K) = PRCL(1,K)
            P2IN(K) = PRCL(2,K)
    9     CONTINUE
          XM1 = AMRCL(1)
          XM2 = AMRCL(2)
          CALL MASHEL(P1IN,P2IN,XM1,XM2,P1OUT,P2OUT,IREJ1)
          IF (IREJ1.GT.0)THEN 
            WRITE(6,'(A)')' FICONF MASHEL rejection'
	    GOTO 9999
          ENDIF
          DO 10 K=1,4
            PRCL(1,K) = P1OUT(K)
            PRCL(2,K) = P2OUT(K)
            PRCLPR(K) = P1OUT(K)
            PRCLTA(K) = P2OUT(K)
   10     CONTINUE
          PRCLPR(5) = AMRCL(1)
          PRCLTA(5) = AMRCL(2)
        ELSE
**sr mod. for DPMJET: IOULEV not available here
          IF(IPEV.GE.1)THEN
	    WRITE(6,'(A)')' from  FICONF'
            DO 7935 IHKK=1,NHKK
            WRITE(6,1005) IHKK, ISTHKK(IHKK),IDHKK(IHKK),JMOHKK(1,IHKK),
     +      JMOHKK(2,IHKK), JDAHKK(1,IHKK),JDAHKK(2,IHKK),(PHKK
     +      (KHKK,IHKK),KHKK=1,5), (VHKK(KHKK,IHKK),KHKK=1,4)
     +      ,IDRES(IHKK),IDXRES(IHKK),NOBAM(IHKK),
     &      IDBAM(IHKK),IDCH(IHKK)
 1005       FORMAT (I6,I4,5I6,9(1PE10.2)/5I6)
 7935       CONTINUE
          ENDIF
          IF(IPEV.GE.1)THEN
            WRITE(LOUT,1001) NEVHKK,INT(AIF(1)),INT(AIZF(1)),
     &                       INT(AIF(2)),INT(AIZF(2)),AMRCL0(1),
     &                       AMRCL(1),AMRCL(1)-AMRCL0(1),AMRCL0(2),
     &                       AMRCL(2),AMRCL(2)-AMRCL0(2)
 1001       FORMAT(1X,'FICONF:   warning! no residual nucleus for',
     &             ' correction',/,11X,'at event',I6,
     &             ',  nucleon config. 1:',2I4,' 2:',2I4,
     &             2(/,11X,3E12.3))
          ENDIF
          GOTO 9998
        ENDIF
      ENDIF

* update counter
      IF (NRESEV(1).NE.NEVHKK) THEN
         NRESEV(1) = NEVHKK
         NRESEV(2) = NRESEV(2)+1
      ENDIF
      DO 15 I=1,2
         EXCDPM(I)   = EXCDPM(I)+EEXC(I)
         EXCDPM(I+2) = EXCDPM(I+2)+(EEXC(I)/MAX(NTOT(I),1))
         NRESTO(I) = NRESTO(I)+NTOT(I)
         NRESPR(I) = NRESPR(I)+NPRO(I)
         NRESNU(I) = NRESNU(I)+NN(I)
         NRESBA(I) = NRESBA(I)+NH(I)
         NRESPB(I) = NRESPB(I)+NHPOS(I)
         NRESCH(I) = NRESCH(I)+NQ(I)
   15 CONTINUE

* evaporation
      IF (LEVPRT) THEN
         DO 13 I=1,2
* initialize evaporation counter
            NP = 0
            EEXCFI(I) = ZERO
            IF ((INUC(I).GT.1).AND.(AIF(I).GT.ONE).AND.
     &          (EEXC(I).GT.ZERO)) THEN
* put residual nuclei into HKKEVT
               IDRCL = 80000
               JMASS = INT( AIF(I))
               JCHAR = INT(AIZF(I))
               CALL EVTPUT(1000,IDRCL,MO1(I),MO2(I),PRCL(I,1),
     &              PRCL(I,2),PRCL(I,3),PRCL(I,4),JMASS,JCHAR,0)
       IF(IPEV.GE.1)WRITE(6,*)PRCL,'PRCL(2,4),EVTPUT'
               DO 14 J=1,4
                  VHKK(J,NHKK) = VRCL(I,J)
                  WHKK(J,NHKK) = WRCL(I,J)
   14          CONTINUE
*  interface to evaporation module - fill final residual nucleus into
*  common RESNUC
               PXRES  = PRCL(I,1)
               PYRES  = PRCL(I,2)
               PZRES  = PRCL(I,3)
C                                              j.r.4.2.97
	       ERES   = PRCL(I,4)
C                                              j.r.4.2.97
               IBRES  = NPRO(I)+NN(I)+NH(I)
               ICRES  = NPRO(I)+NHPOS(I)
               ANOW   = DBLE(IBRES)
               ZNOW   = DBLE(ICRES)
               PTRES  = SQRT(PXRES**2+PYRES**2+PZRES**2)
      IF(IPEV.GE.1)WRITE(6,*)PXRES,PYRES,PZRES,ERES,'FICONF1'
*   ground state mass of the residual nucleus (should be equal to AM0T)
               AMMRES = AMRCL0(I)
               AMNRES = AMMRES-ZNOW*AMELEC+ELBNDE(ICRES)
*  common FINUC
               TV = ZERO
*   kinetic energy of residual nucleus
               TVRECL = PRCL(I,4)-AMRCL(I)
C	       WRITE(6,*)TVRECL, PRCL(I,4),AMRCL(I),'TVRECL'
*   excitation energy of residual nucleus
C                                j.r.16.1.96
               DPMEXM=0.5 
C              TVCMS  = EEXC(I)*DPMEXM
               TVCMS  = EEXC(I)
C	       WRITE(6,*)TVCMS,'TVCMS'
               PTOLD  = PTRES
C                           4.2.98
               PTRES  = SQRT(TVRECL*(TVRECL+2.0D0*(AMMRES+TVCMS)))
               IF (PTOLD.LT.ANGLGB) THEN
                  CALL RACO (PXRES,PYRES,PZRES)
      IF(IPEV.GE.1)WRITE(6,*)PXRES,PYRES,PZRES,ERES,'FICONF2'
                  PTOLD = ONE
               ENDIF
               PXRES = PXRES*PTRES/PTOLD
               PYRES = PYRES*PTRES/PTOLD
               PZRES = PZRES*PTRES/PTOLD
      IF(IPEV.GE.1)WRITE(6,*)PTRES,PTOLD,'FICONF3'
      IF(IPEV.GE.1)WRITE(6,*)PXRES,PYRES,PZRES,ERES,'FICONF3'
* evaporation
	       WE = ONE
C		 WRITE(6,'(A,2F10.2,2I5)')' FRMBRK bef. EVEVAP',
C    *		 ANOW,ZNOW,IBRES,ICRES
		 ANOWW=ANOW
		 ZNOWW=ZNOW
		 IBRESS=IBRES
		 ICRESS=ICRES
		 IREJFR=0
C		 WRITE(6,*)' before EVEVAP, WE',WE
               CALL EVEVAP (WE)
C		 WRITE(6,*)' after EVEVAP , WE',WE
	       IF(IREJFR.EQ.1)THEN
		 WRITE(6,'(A,2F10.2,2I5)')' FRMBRK rej.',
     *		 ANOWW,ZNOWW,IBRESS,ICRESS
		 GO TO 9998
	       ENDIF
* put evaporated particles and residual nuclei to HKKEVT
               MO = NHKK
 	       IF(IPEV.GE.1)WRITE(6,*)EXCITF,'EXITF before EVA2HE'
               CALL EVA2HE(MO,EXCITF,I,IREJ1)
 	       IF(IPEV.GE.1)WRITE(6,*)EXCITF,'EXITF after EVA2HE'
	       IF(IREJ1.GE.1)WRITE(6,'(A)')' FICONF EVA2HE '
               EEXCFI(I) = EXCITF
               EXCEVA(I) = EXCEVA(I)+EXCITF
            ENDIF
   13    CONTINUE
      ENDIF
      IF(IPEV.GE.1)WRITE(6,'(A,I5)')' FICONF RETURN IREJ ',IREJ

      RETURN

 9998 IREXCI(1) = IREXCI(1)+1
 9999 CONTINUE
      LRCLPR = .TRUE.
      LRCLTA = .TRUE.
      IREJ   = IREJ+1
      IF(IPEV.GE.1)WRITE(6,'(A,I5)')' FICONF rej. IREJ ',IREJ
      RETURN
      END
*
*====eva2he============================================================*        
*                                                                      *        
      SUBROUTINE EVA2HE(MO,EEXCF,IRCL,IREJ)

************************************************************************
* Interface between common's of evaporation module (FINUC,FHEAVY)      *
* and HKKEVT.                                                          *
*    MO    HKKEVT-index of "mother" (residual) nucleus before evap.    *
*    EEXCF exitation energy of residual nucleus after evaporation      *
*    IRCL  = 1 projectile residual nucleus                             *
*          = 2 target     residual nucleus                             *
* This version dated 19.04.95 is written by S. Roesler.                *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)
      PARAMETER (TINY10=1.0D-10,TINY3=1.0D-3)

      PARAMETER (NMXHKK=89998) 
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)
* special use for heavy fragments !
*   IDRES(I) = mass number, IDXRES(I) = charge
      COMMON /EXTEVT/ IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),
     &                IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10)
      CHARACTER*8  ANAME
      COMMON /DPAR/   ANAME(210),AAM(210),GA(210),TAU(210),
     &                IICH(210),IIBAR(210),K1(210),K2(210)
      LOGICAL LEMCCK,LHADRO,LSEADI
      COMMON /FLAGS/  IFRAG(2),IRESCO,IMSHL,IRESRJ,
     &                LEMCCK,LHADRO(0:9),LSEADI
      COMMON /STFICO/ EXCDPM(4),EXCEVA(2),
     &                NINCGE,NINCCO(2,3),NINCHR(2,2),NINCWO(2),
     &                NINCST(2,4),NINCEV(2),
     &                NRESTO(2),NRESPR(2),NRESNU(2),NRESBA(2),
     &                NRESPB(2),NRESCH(2),NRESEV(4),
     &                NEVA(2,6),NEVAGA(2),NEVAHT(2),NEVAHY(2,2,240),
     &                NEVAFI(2,2)
      COMMON /EXCITA/ AMRCL0(2),EEXC(2),EEXCFI(2),
     &                NTOT(2),NPRO(2),NN(2),NH(2),NHPOS(2),NQ(2),
     &                NTOTFI(2),NPROFI(2)

      PARAMETER (MXP=999)                                                       
      COMMON / FINUC / CXR (MXP), CYR (MXP), CZR (MXP), TKI (MXP),              
     &                 PLR (MXP), WEI (MXP), TV, TVCMS, TVRECL, TVHEAV,         
     &                 TVBIND, NP0, NP, KPART (MXP)                             

* evaporation interface
      PARAMETER ( MXHEAV = 100 )
      CHARACTER*8 ANHEAV
      COMMON / FHEAVY / CXHEAV (MXHEAV), CYHEAV (MXHEAV),
     &                  CZHEAV (MXHEAV), TKHEAV (MXHEAV),
     &                  PHEAVY (MXHEAV), WHEAVY (MXHEAV),
     &                  AMHEAV  ( 12 ) , AMNHEA  ( 12 ) ,
     &                  KHEAVY (MXHEAV), ICHEAV  ( 12 ) ,
     &                  IBHEAV  ( 12 ) , NPHEAV
      COMMON / FHEAVC / ANHEAV  ( 12 )
      LOGICAL LRNFSS, LFRAGM
      COMMON /RESNUC/  AMNTAR, AMMTAR, AMNZM1, AMMZM1, AMNNM1, AMMNM1,
     &                   ANOW,   ZNOW, ANCOLL, ZNCOLL, AMMLFT, AMNLFT,
     &                   ERES,  EKRES, AMNRES, AMMRES,  PTRES,  PXRES,
     &                  PYRES,  PZRES, PTRES2,  KTARP,  KTARN, IGREYP,
     &                 IGREYN,  ICRES,  IBRES, ISTRES, IEVAPL, IEVAPH,
     &                 IEVNEU, IEVPRO, IEVDEU, IEVTRI, IEV3HE, IEV4HE,
     &                 IDEEXG,  IBTAR, ICHTAR, IBLEFT, ICLEFT, IOTHER,
     &                 LRNFSS, LFRAGM

      DIMENSION IPTOKP(39)
      DATA IPTOKP / 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
     & 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 99,
     & 100, 101, 97, 102, 98, 103, 109, 115 /

      IREJ = 0

* update counter
      IF (NRESEV(3).NE.NEVHKK) THEN
         NRESEV(3) = NEVHKK
         NRESEV(4) = NRESEV(4)+1
      ENDIF

      IF (LEMCCK) 
     &   CALL EVTEMC(PHKK(1,MO),PHKK(2,MO),PHKK(3,MO),PHKK(4,MO),1,
     &                                                   IDUM,IDUM)
* mass number/charge of residual nucleus before evaporation
      IBTOT = IDRES(MO)
      IZTOT = IDXRES(MO)
      IF(IPRI.GE.1)WRITE(6,*)' resnuc IBTOT,IZTOT ',IBTOT,IZTOT
* protons/neutrons/gammas
      DO 1 I=1,NP
         PX    = CXR(I)*PLR(I)
         PY    = CYR(I)*PLR(I)
         PZ    = CZR(I)*PLR(I)
         ID    = IPTOKP(KPART(I))
         IDPDG = IPDGHA(ID)
         AM    = ((PLR(I)+TKI(I))*(PLR(I)-TKI(I)))/
     &           (2.0D0*MAX(TKI(I),TINY10))
         IF (ABS(AM-AAM(ID)).GT.TINY3) THEN
            WRITE(LOUT,1000) ID,AM,AAM(ID)
 1000       FORMAT(1X,'EVA2HE:  inconsistent mass of evap. ',
     &             'particle',I3,2E10.3)
         ENDIF
         PE = TKI(I)+AM
         CALL EVTPUT(-1,IDPDG,MO,0,PX,PY,PZ,PE,0,0,0)
         NOBAM(NHKK) = IRCL
         IF (LEMCCK) CALL EVTEMC(-PX,-PY,-PZ,-PE,2,IDUM,IDUM)
         IBTOT = IBTOT-IIBAR(ID)
         IZTOT = IZTOT-IICH(ID)
    1 CONTINUE

* heavy fragments
      DO 2 I=1,NPHEAV
         PX     = CXHEAV(I)*PHEAVY(I)
         PY     = CYHEAV(I)*PHEAVY(I)
         PZ     = CZHEAV(I)*PHEAVY(I)
         IDHEAV = 80000
         AM     = ((PHEAVY(I)+TKHEAV(I))*(PHEAVY(I)-TKHEAV(I)))/
     &            (2.0D0*MAX(TKHEAV(I),TINY10))
         PE     = TKHEAV(I)+AM
         CALL EVTPUT(-1,IDHEAV,MO,0,PX,PY,PZ,PE,
     &                  IBHEAV(KHEAVY(I)),ICHEAV(KHEAVY(I)),0)
         NOBAM(NHKK) = IRCL
         IF (LEMCCK) CALL EVTEMC(-PX,-PY,-PZ,-PE,2,IDUM,IDUM)
         IBTOT = IBTOT-IBHEAV(KHEAVY(I))
         IZTOT = IZTOT-ICHEAV(KHEAVY(I))
    2 CONTINUE

      IF (IBRES.GT.0) THEN
* residual nucleus after evaporation
         IDNUC = 80000
         CALL EVTPUT(1001,IDNUC,MO,0,PXRES,PYRES,PZRES,ERES,
     &                                        IBRES,ICRES,0)
C     WRITE(6,*)PXRES,PYRES,PZRES,ERES,'EVTPUT1001'
         NOBAM(NHKK) = IRCL
      ENDIF
      EEXCF = TVCMS 
      NTOTFI(IRCL) = IBRES
      NPROFI(IRCL) = ICRES
      IF (LEMCCK) CALL EVTEMC(-PXRES,-PYRES,-PZRES,-ERES,2,IDUM,IDUM)
      IBTOT = IBTOT-IBRES
      IZTOT = IZTOT-ICRES

* count events with fission
      NEVAFI(1,IRCL) = NEVAFI(1,IRCL)+1
      IF (LRNFSS) NEVAFI(2,IRCL) = NEVAFI(2,IRCL)+1

* energy-momentum conservation check
      IF (LEMCCK) CALL EVTEMC(DUM,DUM,DUM,DUM,4,40,IREJ)
* baryon-number/charge conservation check
      IF (IBTOT+IZTOT.NE.0) THEN
         WRITE(LOUT,1001) NEVHKK,IBTOT,IZTOT
 1001    FORMAT(1X,'EVA2HE:   baryon-number/charge conservation ',
     &          'failure at event ',I6,' :  IBTOT,IZTOT = ',2I3)
      ENDIF

      RETURN
      END
*
*===fozoca=============================================================*
*
      SUBROUTINE FOZOCA(LFZC,IREJ)

************************************************************************
* This subroutine treats the complete FOrmation ZOne supressed intra-  *
* nuclear CAscade.                                                     *
*               LFZC = .true.  cascade has been treated                *
*                    = .false. cascade skipped                         *
* This is a completely revised version of the original FOZOKL.         *
* This version dated 18.11.95 is written by S. Roesler                 *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)
      PARAMETER (DLARGE=1.0D10,OHALF=0.5D0,ZERO=0.0D0)
      PARAMETER (FM2MM=1.0D-12,RNUCLE = 1.12D0)

      LOGICAL LSTART,LCAS,LFZC

      PARAMETER (NMXHKK=89998)
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)
      COMMON /EXTEVT/ IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),
     &                IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10)
      COMMON /RJCOUN/ IRPT,IRHHA,IRRES(2),LOMRES,LOBRES,
     &                IRCHKI(2),IRFRAG,IRCRON(3),IREVT,
     &                IREXCI(3),IRDIFF(2),IRINC

      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
      COMMON /RPTSHM/ RPROJ,RTARG,BIMPAC
      COMMON /DPRIN/  IPRI,IPEV,IPPA,IPCO,INIT,IPHKK,ITOPD,IPAUPR
      LOGICAL LEMCCK,LHADRO,LSEADI
      COMMON /FLAGS/  IFRAG(2),IRESCO,IMSHL,IRESRJ,
     &                LEMCCK,LHADRO(0:9),LSEADI
      COMMON /PAULI/  EWOUND(2,300),NWOUND(2),IDXINC(2000),NOINC

      COMMON /TAUFO/  TAUFOR,KTAUGE,ITAUVE,INCMOD
**sr mod. for DPMJET: use the longer DPMJET one
      LOGICAL INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA, IPADIS,
     &        ISHMAL,LPAULI
      COMMON /DROPPT/ INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA,
     &                IPADIS,ISHMAL,LPAULI
**

      DATA LSTART /.TRUE./

      DIMENSION NCWOUN(2)

      LFZC = .TRUE.
      IREJ = 0

* skip cascade if hadron-hadron interaction or if supressed by user
      IF (((IP.EQ.1).AND.(IT.EQ.1)).OR.(KTAUGE.LT.1)) GOTO 9999
* skip cascade if not all possible chains systems are hadronized
      IF(IPEV.GE.1)WRITE(6,*)LHADRO
      DO 1 I=1,8
         IF (.NOT.(LHADRO(I))) GOTO 9999
    1 CONTINUE

      IF (IPEV.GE.1) THEN
         WRITE(6,1000) KTAUGE,TAUFOR,INCMOD
      ENDIF
      IF (LSTART) THEN
         WRITE(LOUT,1000) KTAUGE,TAUFOR,INCMOD
 1000    FORMAT(/,1X,'FOZOCA:  intranuclear cascade treated for a ',
     &          'maximum of',I4,' generations',/,10X,'formation time ',
     &          'parameter:',F5.1,'  fm/c',9X,'modus:',I2)
         IF (ITAUVE.EQ.1) WRITE(LOUT,1001)
         IF (ITAUVE.EQ.2) WRITE(LOUT,1002)
 1001    FORMAT(10X,'p_t dependent formation zone',/)
 1002    FORMAT(10X,'constant formation zone',/)
         LSTART = .FALSE.
      ENDIF

* in order to avoid wasting of cpu-time the HKKEVT-indices of nucleons
* which may interact with final state particles are stored in a seperate
* array - here all proj./target nucleon-indices (just for simplicity)
      NOINC = 0
      DO 9 I=1,NPOINT(1)-1
         NOINC = NOINC+1
         IDXINC(NOINC) = I
    9 CONTINUE

* initialize Pauli-principle treatment (find wounded nucleons)
      NWOUND(1) = 0
      NWOUND(2) = 0
      NCWOUN(1) = 0
      NCWOUN(2) = 0
      DO 2 J=1,NPOINT(1)
         DO 3 I=1,2
            IF (ISTHKK(J).EQ.10+I) THEN
               NWOUND(I) = NWOUND(I)+1
               EWOUND(I,NWOUND(I)) = PHKK(4,J)
               IF (IDHKK(J).EQ.2212) NCWOUN(I) = NCWOUN(I)+1
            ENDIF
    3    CONTINUE
    2 CONTINUE

* modify nuclear potential for wounded nucleons
      IPRCL  = IP -NWOUND(1)
      IPZRCL = IPZ-NCWOUN(1)
      ITRCL  = IT -NWOUND(2)
      ITZRCL = ITZ-NCWOUN(2)
      CALL NCLPOT(IPZRCL,IPRCL,ITZRCL,ITRCL,ZERO,ZERO,1)

      NSTART = NPOINT(4)
      NEND   = NHKK

    7 CONTINUE
      DO 8 I=NSTART,NEND

         IF ((ABS(ISTHKK(I)).EQ.1).AND.(IDCH(I).LT.KTAUGE)) THEN

* select nucleus the cascade starts first (proj. - 1, target - -1)
            NCAS   = 1
*   projectile/target with probab. 1/2
            IF ((INCMOD.EQ.1).OR.(IDCH(I).GT.0)) THEN
               IF (RNDM(V).GT.OHALF) NCAS = -NCAS
*   in the nucleus with highest mass
            ELSEIF (INCMOD.EQ.2) THEN
               IF (IP.GT.IT) THEN
                  NCAS = -NCAS
               ELSEIF (IP.EQ.IT) THEN
                  IF (RNDM(V).GT.OHALF) NCAS = -NCAS
               ENDIF
* the nucleus the cascade starts first is requested to be the one
* moving in the direction of the secondary
            ELSEIF (INCMOD.EQ.3) THEN
               NCAS = INT(SIGN(1.0D0,PHKK(3,I)))
            ENDIF
* check that the selected "nucleus" is not a hadron
            IF (((NCAS.EQ. 1).AND.(IP.LE.1)).OR.
     &          ((NCAS.EQ.-1).AND.(IT.LE.1)))    NCAS = -NCAS

* treat intranuclear cascade in the nucleus selected first
            LCAS = .FALSE.
            CALL INUCAS(IT,IP,I,LCAS,NCAS,IREJ1)
            IF (IREJ1.NE.0)THEN 
C	      WRITE(6,'(A)')' INUCAS Rejection'
	      GOTO 9998
            ENDIF
* treat intranuclear cascade in the other nucleus if this isn't a had.
            NCAS = -NCAS
            IF (((NCAS.EQ. 1).AND.(IP.GT.1)).OR.
     &          ((NCAS.EQ.-1).AND.(IT.GT.1)))    THEN
               IF (LCAS) CALL INUCAS(IT,IP,I,LCAS,NCAS,IREJ1)
               IF (IREJ1.NE.0)THEN
C	         WRITE(6,'(A)')' INUCAS Rejection2'
	         GOTO 9998
               ENDIF
            ENDIF

         ENDIF

    8 CONTINUE
      NSTART = NEND+1
      NEND   = NHKK
      IF (NSTART.LE.NEND) GOTO 7

      RETURN

 9998 CONTINUE
* reject this event
      IRINC = IRINC+1
      IREJ = 1

 9999 CONTINUE
* intranucl. cascade not treated because of interaction properties or
* it is supressed by user or it was rejected or...
      LFZC = .FALSE.
* reset flag characterizing direction of motion in n-n-cms
**sr14-11-95
C     DO 9990 I=NPOINT(5),NHKK
C        IF (ISTHKK(I).EQ.-1) ISTHKK(I)=1
C9990 CONTINUE

      RETURN
      END
*
*
*===inucas=============================================================*
*
      SUBROUTINE INUCAS(IT,IP,IDXCAS,LCAS,NCAS,IREJ)

************************************************************************
* Formation zone supressed IntraNUclear CAScade for one final state    *
* particle.                                                            *
*           IT, IP    mass numbers of target, projectile nuclei        *
*           IDXCAS    index of final state particle in HKKEVT          *
*           NCAS =  1 intranuclear cascade in projectile               *
*                = -1 intranuclear cascade in target                   *
* This version dated 11.06.96 is written by S. Roesler                 *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)

      PARAMETER (TINY10=1.0D-10,TINY2=1.0D-2,ZERO=0.0D0,DLARGE=1.0D10,
     &           OHALF=0.5D0,ONE=1.0D0)
      PARAMETER (FM2MM=1.0D-12,RNUCLE = 1.12D0)
      PARAMETER (TWOPI=6.283185307179586454D+00)
      PARAMETER (ELOWH=0.01D0,EHIH=9.0D0)

      LOGICAL LABSOR,LCAS

      PARAMETER (NMXHKK=89998)
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)
      COMMON /EXTEVT/ IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),
     &                IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10)
      PARAMETER (MAXFSP=10)
      COMMON /FISTAT/ PFSP(5,MAXFSP),IDFSP(MAXFSP),NFSP

**sr mod. for DPMJET: the old shorter version of /FLAGS/
      LOGICAL LEMCCK,LHADRO,LSEADI
      COMMON /FLAGS/  IFRAG(2),IRESCO,IMSHL,IRESRJ,
     &                LEMCCK,LHADRO(0:9),LSEADI
      CHARACTER*8 ANAME
      COMMON /DPAR/   ANAME(210),AAM(210),GA(210),TAU(210),
     &                IICH(210),IIBAR(210),K1(210),K2(210)

      COMMON /RPTSHM/ RPROJ,RTARG,BIMPAC
      COMMON /NUCLEA/ PFERMP(2),PFERMN(2),FERMOD,
     &                EBINDP(2),EBINDN(2),EPOT(2,210),
     &                ETACOU(2),ICOUL
      COMMON /TAUFO/  TAUFOR,KTAUGE,ITAUVE,INCMOD
      COMMON /PAULI/  EWOUND(2,300),NWOUND(2),IDXINC(2000),NOINC
**sr mod. for DPMJET: use the longer DPMJET one
      LOGICAL INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA, IPADIS,
     &        ISHMAL,LPAULI
      COMMON /DROPPT/ INTPT,FERMP,IHADSS,IHADSV,IHADVS,IHADVV,IHADA,
     &                IPADIS,ISHMAL,LPAULI

      COMMON /STFICO/ EXCDPM(4),EXCEVA(2),
     &                NINCGE,NINCCO(2,3),NINCHR(2,2),NINCWO(2),
     &                NINCST(2,4),NINCEV(2),
     &                NRESTO(2),NRESPR(2),NRESNU(2),NRESBA(2),
     &                NRESPB(2),NRESCH(2),NRESEV(4),
     &                NEVA(2,6),NEVAGA(2),NEVAHT(2),NEVAHY(2,2,240),
     &                NEVAFI(2,2)

      DIMENSION PCAS(2,5),PTOCAS(2),COSCAS(2,3),VTXCAS(2,4),VTXCA1(2,4),
     &          PCAS1(5),PNUC(5),BGTA(4),
     &          BGCAS(2),GACAS(2),BECAS(2),
     &          RNUC(2),BIMPC(2),VTXDST(3),IDXSPE(2),IDSPE(2),NWTMP(2)

      DATA PDIF /0.545D0/

      IREJ = 0

* update counter
      IF (NINCEV(1).NE.NEVHKK) THEN
         NINCEV(1) = NEVHKK
         NINCEV(2) = NINCEV(2)+1
      ENDIF

* "BAMJET-index" of this hadron
      IDCAS = IDBAM(IDXCAS)
      IF (MCHAD(IDCAS).EQ.-1) RETURN

* skip gammas, electrons, etc..
      IF (AAM(IDCAS).LT.TINY2) RETURN

* Lorentz-trsf. into projectile rest system
      IF (IP.GT.1) THEN
         CALL LTRANS(PHKK(1,IDXCAS),PHKK(2,IDXCAS),PHKK(3,IDXCAS),
     &               PHKK(4,IDXCAS),PCAS(1,1),PCAS(1,2),PCAS(1,3),
     &               PCAS(1,4),IDCAS,-2)
         PTOCAS(1) = SQRT(PCAS(1,1)**2+PCAS(1,2)**2+PCAS(1,3)**2)
         PCAS(1,5) = (PCAS(1,4)-PTOCAS(1))*(PCAS(1,4)+PTOCAS(1))
         IF (PCAS(1,5).GT.ZERO) THEN
            PCAS(1,5) = SQRT(PCAS(1,5))
         ELSE
            PCAS(1,5) = AAM(IDCAS)
         ENDIF
         DO 20 K=1,3
            COSCAS(1,K) = PCAS(1,K)/MAX(PTOCAS(1),TINY10)
   20    CONTINUE
* Lorentz-parameters
*   particle rest system --> projectile rest system
         BGCAS(1) = PTOCAS(1)/MAX(PCAS(1,5),TINY10)
         GACAS(1) = PCAS(1,4)/MAX(PCAS(1,5),TINY10)
         BECAS(1) = BGCAS(1)/GACAS(1)
      ELSE
         DO 21 K=1,5
            PCAS(1,K) = ZERO
            IF (K.LE.3) COSCAS(1,K) = ZERO
   21    CONTINUE
         PTOCAS(1) = ZERO
         BGCAS(1)  = ZERO
         GACAS(1)  = ZERO
         BECAS(1)  = ZERO
      ENDIF
* Lorentz-trsf. into target rest system
      IF (IT.GT.1) THEN
         CALL LTRANS(PHKK(1,IDXCAS),PHKK(2,IDXCAS),PHKK(3,IDXCAS),
     &               PHKK(4,IDXCAS),PCAS(2,1),PCAS(2,2),PCAS(2,3),
     &               PCAS(2,4),IDCAS,-3)
         PTOCAS(2) = SQRT(PCAS(2,1)**2+PCAS(2,2)**2+PCAS(2,3)**2)
         PCAS(2,5) = (PCAS(2,4)-PTOCAS(2))*(PCAS(2,4)+PTOCAS(2))
         IF (PCAS(2,5).GT.ZERO) THEN
            PCAS(2,5) = SQRT(PCAS(2,5))
         ELSE
            PCAS(2,5) = AAM(IDCAS)
         ENDIF
         DO 22 K=1,3
            COSCAS(2,K) = PCAS(2,K)/MAX(PTOCAS(2),TINY10)
   22    CONTINUE
* Lorentz-parameters
*   particle rest system --> target rest system
         BGCAS(2) = PTOCAS(2)/MAX(PCAS(2,5),TINY10)
         GACAS(2) = PCAS(2,4)/MAX(PCAS(2,5),TINY10)
         BECAS(2) = BGCAS(2)/GACAS(2)
      ELSE
         DO 23 K=1,5
            PCAS(2,K) = ZERO
            IF (K.LE.3) COSCAS(2,K) = ZERO
   23    CONTINUE
         PTOCAS(2) = ZERO
         BGCAS(2)  = ZERO
         GACAS(2)  = ZERO
         BECAS(2)  = ZERO
      ENDIF

* radii of nuclei (mm) modified by the wall-depth of the Woods-Saxon-
* potential (see CONUCL)
      RNUC(1)  = (RPROJ+4.605D0*PDIF)*FM2MM
      RNUC(2)  = (RTARG+4.605D0*PDIF)*FM2MM
* impact parameter (the projectile moving along z)
      BIMPC(1) = ZERO
      BIMPC(2) = BIMPAC*FM2MM

* get position of initial hadron in projectile/target rest-syst.
      DO 3 K=1,4
         VTXCAS(1,K) = WHKK(K,IDXCAS)
         VTXCAS(2,K) = VHKK(K,IDXCAS)
    3 CONTINUE

      ICAS = 1
      I2   = 2
      IF (NCAS.EQ.-1) THEN
         ICAS = 2
         I2   = 1
      ENDIF

      IF (PTOCAS(ICAS).LT.TINY10) THEN
         WRITE(LOUT,1000) PTOCAS
 1000    FORMAT(1X,'INUCAS:   warning! zero momentum of initial',
     &          '  hadron ',/,20X,2E12.4)
         GOTO 9999
      ENDIF

* reset spectator flags
      NSPE = 0
      IDXSPE(1) = 0
      IDXSPE(2) = 0
      IDSPE(1)  = 0
      IDSPE(2)  = 0

* formation length (in fm)
C     IF (LCAS) THEN
C        DEL0 = ZERO
C     ELSE
         DEL0 = TAUFOR*BGCAS(ICAS)
         IF (ITAUVE.EQ.1) THEN
            AMT  = PCAS(ICAS,1)**2+PCAS(ICAS,2)**2+PCAS(ICAS,5)**2
            DEL0 = DEL0*PCAS(ICAS,5)**2/AMT
         ENDIF
C     ENDIF
*   sample from exp(-del/del0)
      DEL1   = -DEL0*LOG(MAX(RNDM(V),TINY10))
* save formation time
      TAUSA1 = DEL1/BGCAS(ICAS)
      REL1   = TAUSA1*BGCAS(I2)

      DEL    = DEL1
      TAUSAM = DEL/BGCAS(ICAS)
      REL    = TAUSAM*BGCAS(I2)

* special treatment for negative particles unable to escape
* nuclear potential (implemented for ap, pi-, K- only)
      LABSOR = .FALSE.
      IF ((IICH(IDCAS).EQ.-1).AND.(IDCAS.LT.20)) THEN
*   threshold energy = nuclear potential + Coulomb potential
*   (nuclear potential for hadron-nucleus interactions only)
         ETHR = AAM(IDCAS)+EPOT(ICAS,IDCAS)+ETACOU(ICAS)
         IF (PCAS(ICAS,4).LT.ETHR) THEN
            DO 4 K=1,5
               PCAS1(K) = PCAS(ICAS,K)
    4       CONTINUE
*   "absorb" negative particle in nucleus
            CALL ABSORP(IDCAS,PCAS1,NCAS,NSPE,IDSPE,IDXSPE,0,IREJ1)
            IF (IREJ1.NE.0) GOTO 9999
            IF (NSPE.GE.1) LABSOR = .TRUE.
         ENDIF
      ENDIF

* if the initial particle has not been absorbed proceed with
* "normal" cascade
      IF (.NOT.LABSOR) THEN

*   calculate coordinates of hadron at the end of the formation zone
*   transport-time and -step in the rest system where this step is
*   treated
         DSTEP  = DEL*FM2MM
         DTIME  = DSTEP/BECAS(ICAS)
         RSTEP  = REL*FM2MM
         IF ((IP.GT.1).AND.(IT.GT.1)) THEN
            RTIME = RSTEP/BECAS(I2)
         ELSE
            RTIME = ZERO
         ENDIF
*   save step whithout considering the overlapping region
         DSTEP1 = DEL1*FM2MM
         DTIME1 = DSTEP1/BECAS(ICAS)
         RSTEP1 = REL1*FM2MM
         IF ((IP.GT.1).AND.(IT.GT.1)) THEN
            RTIME1 = RSTEP1/BECAS(I2)
         ELSE
            RTIME1 = ZERO
         ENDIF
*   transport to the end of the formation zone in this system
         DO 5 K=1,3
            VTXCA1(ICAS,K) = VTXCAS(ICAS,K)+DSTEP1*COSCAS(ICAS,K)
            VTXCA1(I2,K)   = VTXCAS(I2,K)  +RSTEP1*COSCAS(I2,K)
            VTXCAS(ICAS,K) = VTXCAS(ICAS,K)+DSTEP*COSCAS(ICAS,K)
            VTXCAS(I2,K)   = VTXCAS(I2,K)  +RSTEP*COSCAS(I2,K)
    5    CONTINUE
         VTXCA1(ICAS,4) = VTXCAS(ICAS,4)+DTIME1
         VTXCA1(I2,4)   = VTXCAS(I2,4)  +RTIME1
         VTXCAS(ICAS,4) = VTXCAS(ICAS,4)+DTIME
         VTXCAS(I2,4)   = VTXCAS(I2,4)  +RTIME

         IF ((IP.GT.1).AND.(IT.GT.1)) THEN
            XCAS   = VTXCAS(ICAS,1)
            YCAS   = VTXCAS(ICAS,2)
            XNCLTA = BIMPAC*FM2MM
            RNCLPR = (RPROJ+RNUCLE)*FM2MM
            RNCLTA = (RTARG+RNUCLE)*FM2MM
            RCASPR = SQRT( XCAS**2        +YCAS**2)
            RCASTA = SQRT((XCAS-XNCLTA)**2+YCAS**2)
            IF ((RCASPR.LT.RNCLPR).AND.(RCASTA.LT.RNCLTA)) THEN
               IF (IDCH(IDXCAS).EQ.0) NOBAM(IDXCAS) = 3
            ENDIF
         ENDIF

*   check if particle is already outside of the corresp. nucleus
         RDIST = SQRT((VTXCAS(ICAS,1)-BIMPC(ICAS))**2+
     &                VTXCAS(ICAS,2)**2+VTXCAS(ICAS,3)**2)
         IF (RDIST.GE.RNUC(ICAS)) THEN
*   here: IDCH is the generation of the final state part. starting
*   with zero for hadronization products
*   flag particles of generation 0 being outside the nuclei after
*   formation time (to be used for excitation energy calculation)
            IF ((IDCH(IDXCAS).EQ.0).AND.(NOBAM(IDXCAS).LT.3))
     &         NOBAM(IDXCAS) = NOBAM(IDXCAS)+ICAS
            GOTO 9997
         ENDIF
         DIST   = DLARGE
         DISTP  = DLARGE
         DISTN  = DLARGE
         IDXP   = 0
         IDXN   = 0

*   already here: skip particles being outside HADRIN "energy-window"
*   to avoid wasting of time
         NINCHR(ICAS,1) = NINCHR(ICAS,1)+1
         IF ((PTOCAS(ICAS).LE.ELOWH).OR.(PTOCAS(ICAS).GE.EHIH)) THEN
            NINCHR(ICAS,2) = NINCHR(ICAS,2)+1
C           WRITE(LOUT,1002) IDXCAS,IDCAS,ICAS,PTOCAS(ICAS),NEVHKK
C1002       FORMAT(1X,'INUCAS:   warning! momentum of particle with ',
C    &             'index ',I5,' (id: ',I3,') ',I3,/,11X,'p_tot = ',
C    &             E12.4,', above or below HADRIN-thresholds',I6)
            NSPE = 0
            GOTO 9997
         ENDIF

         DO 7 IDXHKK=1,NOINC
            I = IDXINC(IDXHKK)
*   scan HKKEVT for unwounded or excited nucleons
            IF ((ISTHKK(I).EQ.12+ICAS).OR.(ISTHKK(I).EQ.14+ICAS)) THEN
               DO 8 K=1,3
                  IF (ICAS.EQ.1) THEN
                     VTXDST(K) = WHKK(K,I)-VTXCAS(1,K)
                  ELSEIF (ICAS.EQ.2) THEN
                     VTXDST(K) = VHKK(K,I)-VTXCAS(2,K)
                  ENDIF
    8          CONTINUE
               POSNUC = VTXDST(1)*COSCAS(ICAS,1)+
     &                  VTXDST(2)*COSCAS(ICAS,2)+
     &                  VTXDST(3)*COSCAS(ICAS,3)
*   check if nucleon is situated in forward direction
               IF (POSNUC.GT.ZERO) THEN
*   distance between hadron and this nucleon
                  DISTNU = SQRT(VTXDST(1)**2+VTXDST(2)**2+
     &                          VTXDST(3)**2)
*   impact parameter
                  BIMNU2 = DISTNU**2-POSNUC**2
                  IF (BIMNU2.LT.ZERO) THEN
                     WRITE(LOUT,1001) DISTNU,POSNUC,BIMNU2
 1001                FORMAT(1X,'INUCAS:   warning! inconsistent impact',
     &                      '  parameter ',/,20X,3E12.4)
                     GOTO 7
                  ENDIF
                  BIMNU  = SQRT(BIMNU2)
*   maximum impact parameter to have interaction
                  IDNUC  = ICIHAD(IDHKK(I))
                  IDNUC1 = MCHAD(IDNUC)
                  IDCAS1 = MCHAD(IDCAS)
                  DO 19 K=1,5
                     PCAS1(K) = PCAS(ICAS,K)
                     PNUC(K)  = PHKK(K,I)
   19             CONTINUE
* Lorentz-parameter for trafo into rest-system of target
                  DO 18 K=1,4
                     BGTA(K) = PNUC(K)/MAX(PNUC(5),TINY10)
   18             CONTINUE
* transformation of projectile into rest-system of target
                  CALL DALTRA(BGTA(4),-BGTA(1),-BGTA(2),-BGTA(3),
     &                        PCAS1(1),PCAS1(2),PCAS1(3),PCAS1(4),
     &                        PPTOT,PX,PY,PZ,PE)
                  CALL SIHNIN(IDCAS1,IDNUC1,PPTOT,SIGIN)
                  CALL SIHNEL(IDCAS1,IDNUC1,PPTOT,SIGEL)
                  CALL SIHNAB(IDCAS1,IDNUC1,PPTOT,SIGAB)
                  SIGTOT = SIGIN+SIGEL+SIGAB
                  BIMMAX = SQRT(SIGTOT/(5.0D0*TWOPI))*FM2MM
*   check if interaction is possible
                  IF (BIMNU.LE.BIMMAX) THEN
*   get nucleon with smallest distance and kind of interaction
*   (elastic/inelastic)
                     IF (DISTNU.LT.DIST) THEN
                        DIST      = DISTNU
                        BINT      = BIMNU
                        IF (IDNUC.NE.IDSPE(1)) THEN
                           IDSPE(2)  = IDSPE(1)
                           IDXSPE(2) = IDXSPE(1)
                           IDSPE(1)  = IDNUC
                        ENDIF
                        IDXSPE(1) = I
                        NSPE      = 1
**sr
                        SELA = SIGEL
                        SABS = SIGAB
                        STOT = SIGTOT
C                       IF ((IDCAS.EQ.2).OR.(IDCAS.EQ.9)) THEN
C                          SELA = SIGEL
C                          STOT = SIGIN+SIGEL
C                       ELSE
C                          SELA = SIGEL+0.75D0*SIGIN
C                          STOT = 0.25D0*SIGIN+SELA
C                       ENDIF
**
                     ENDIF
                  ENDIf
               ENDIF
               DISTNU = SQRT(VTXDST(1)**2+VTXDST(2)**2+
     &                       VTXDST(3)**2)
               IDNUC  = ICIHAD(IDHKK(I))
               IF (IDNUC.EQ.1) THEN
                  IF (DISTNU.LT.DISTP) THEN
                     DISTP = DISTNU
                     IDXP  = I
                     POSP  = POSNUC
                  ENDIF
               ELSEIF (IDNUC.EQ.8) THEN
                  IF (DISTNU.LT.DISTN) THEN
                     DISTN = DISTNU
                     IDXN  = I
                     POSN  = POSNUC
                  ENDIF
               ENDIF
            ENDIF
    7    CONTINUE

* there is no nucleon for a secondary interaction
         IF (NSPE.EQ.0) GOTO 9997

         IF (IDXSPE(2).EQ.0) THEN
            IF ((IDSPE(1).EQ.1).AND.(IDXN.GT.0)) THEN
               IDXSPE(2) = IDXN
               IDSPE(2)  = 8
            ELSEIF ((IDSPE(1).EQ.8).AND.(IDXP.GT.0)) THEN
               IDXSPE(2) = IDXP
               IDSPE(2)  = 1
            ELSE
               STOT = STOT-SABS
               SABS = ZERO
            ENDIF
         ENDIF
         RR = RNDM(V)
         IF (RR.LT.SELA/STOT) THEN
            IPROC = 2
         ELSEIF ((RR.GE.SELA/STOT).AND.(RR.LT.(SELA+SABS)/STOT)) THEN
            IPROC = 3
         ELSE
            IPROC = 1
         ENDIF

         DO 9 K=1,5
            PCAS1(K) = PCAS(ICAS,K)
            PNUC(K)  = PHKK(K,IDXSPE(1))
    9    CONTINUE
         IF (IPROC.EQ.3) THEN
* 2-nucleon absorption of pion
            NSPE = 2
            CALL ABSORP(IDCAS,PCAS1,NCAS,NSPE,IDSPE,IDXSPE,1,IREJ1)
            IF (IREJ1.NE.0) GOTO 9999
            IF (NSPE.GE.1) LABSOR = .TRUE.
         ELSE
* sample secondary interaction
            IDNUC = IDBAM(IDXSPE(1))
**sr mod. for DPMJET: HADRIN-->HADRI1
            CALL HADRI1(IDCAS,PCAS1,IDNUC,PNUC,IPROC,IREJ1)
**sr mod. for DPMJET: in case of rejections jump to 9998 rather than
*                     reject cascade completely (??)
C           IF (IREJ1.EQ.1) GOTO 9999
            IF (IREJ1.GE.1)THEN
C              WRITE(6,'(A)')' HADRI1 Rejection'
               GOTO 9998
            ENDIF
         ENDIF
      ENDIF

* update arrays to include Pauli-principle
      DO 10 I=1,NSPE
         IF (NWOUND(ICAS).LE.299) THEN
            NWOUND(ICAS) = NWOUND(ICAS)+1
            EWOUND(ICAS,NWOUND(ICAS)) = PHKK(4,IDXSPE(I))
         ENDIF
   10 CONTINUE

* dump initial hadron for energy-momentum conservation check
      IF (LEMCCK)
     &   CALL EVTEMC(PCAS(ICAS,1),PCAS(ICAS,2),PCAS(ICAS,3),
     &               PCAS(ICAS,4),1,IDUM,IDUM)

* dump final state particles into HKKEVT

*   check if Pauli-principle is fulfilled
      NPAULI = 0
      NWTMP(1) = NWOUND(1)
      NWTMP(2) = NWOUND(2)
      DO 111 I=1,NFSP
         NPAULI = 0
         J1 = 2
         IF (((NCAS.EQ. 1).AND.(IT.LE.1)).OR.
     &       ((NCAS.EQ.-1).AND.(IP.LE.1)))    J1 = 1
         DO 117 J=1,J1
            IF ((NPAULI.NE.0).AND.(J.EQ.2)) GOTO 117
            IF (J.EQ.1) THEN
               IDX = ICAS
               PE  = PFSP(4,I)
            ELSE
               IDX  = I2
               MODE = 1
               IF (IDX.EQ.1) MODE = -1
               CALL LTNUC(PFSP(3,I),PFSP(4,I),PZ,PE,MODE)
            ENDIF
* first check if cascade step is forbidden due to Pauli-principle
* (in case of absorpion this step is forced)
            IF ((.NOT.LABSOR).AND.LPAULI.AND.((IDFSP(I).EQ.1).OR.
     &          (IDFSP(I).EQ.8))) THEN
*   get nuclear potential barrier
               POT = EPOT(IDX,IDFSP(I))+AAM(IDFSP(I))
               IF (IDFSP(I).EQ.1) THEN
                  POTLOW = POT-EBINDP(IDX)
               ELSE
                  POTLOW = POT-EBINDN(IDX)
               ENDIF
*   final state particle not able to escape nucleus
               IF (PE.LE.POTLOW) THEN
*     check if there are wounded nucleons
                  IF ((NWOUND(IDX).GE.1).AND.(PE.GE.
     &                 EWOUND(IDX,NWOUND(IDX)))) THEN
                     NPAULI      = NPAULI+1
                     NWOUND(IDX) = NWOUND(IDX)-1
                  ELSE
*     interaction prohibited by Pauli-principle
                     NWOUND(1) = NWTMP(1)
                     NWOUND(2) = NWTMP(2)
                     GOTO 9997
                  ENDIF
               ENDIF
            ENDIF
  117    CONTINUE
  111 CONTINUE

      NPAULI = 0
      NWOUND(1) = NWTMP(1)
      NWOUND(2) = NWTMP(2)

      DO 11 I=1,NFSP

         IST = ISTHKK(IDXCAS)

         NPAULI = 0
         J1 = 2
         IF (((NCAS.EQ. 1).AND.(IT.LE.1)).OR.
     &       ((NCAS.EQ.-1).AND.(IP.LE.1)))    J1 = 1
         DO 17 J=1,J1
            IF ((NPAULI.NE.0).AND.(J.EQ.2)) GOTO 17
            IDX = ICAS
            PE  = PFSP(4,I)
            IF (J.EQ.2) THEN
               IDX = I2
               CALL LTNUC(PFSP(3,I),PFSP(4,I),PZ,PE,NCAS)
            ENDIF
* first check if cascade step is forbidden due to Pauli-principle
* (in case of absorpion this step is forced)
            IF ((.NOT.LABSOR).AND.LPAULI.AND.((IDFSP(I).EQ.1).OR.
     &          (IDFSP(I).EQ.8))) THEN
*   get nuclear potential barrier
               POT = EPOT(IDX,IDFSP(I))+AAM(IDFSP(I))
               IF (IDFSP(I).EQ.1) THEN
                  POTLOW = POT-EBINDP(IDX)
               ELSE
                  POTLOW = POT-EBINDN(IDX)
               ENDIF
*   final state particle not able to escape nucleus
               IF (PE.LE.POTLOW) THEN
*     check if there are wounded nucleons
                  IF ((NWOUND(IDX).GE.1).AND.(PE.GE.
     &                 EWOUND(IDX,NWOUND(IDX)))) THEN
                     NWOUND(IDX) = NWOUND(IDX)-1
                     NPAULI = NPAULI+1
                     IST    = 14+IDX
                  ELSE
*     interaction prohibited by Pauli-principle
                     NWOUND(1) = NWTMP(1)
                     NWOUND(2) = NWTMP(2)
                     GOTO 9997
                  ENDIF
**sr
c               ELSEIF (PE.LE.POT) THEN
cC              ELSEIF ((PE.LE.POT).AND.(NWOUND(IDX).GE.1)) THEN
cC                 NWOUND(IDX) = NWOUND(IDX)-1
c**
c                  NPAULI = NPAULI+1
c                  IST    = 14+IDX
               ENDIF
            ENDIF
   17    CONTINUE

* dump final state particles for energy-momentum conservation check
         IF (LEMCCK) CALL EVTEMC(-PFSP(1,I),-PFSP(2,I),-PFSP(3,I),
     &                           -PFSP(4,I),2,IDUM,IDUM)

         PX = PFSP(1,I)
         PY = PFSP(2,I)
         PZ = PFSP(3,I)
         PE = PFSP(4,I)
         IF (ABS(IST).EQ.1) THEN
* transform particles back into n-n cms
            IMODE = ICAS+1
            CALL LTRANS(PX,PY,PZ,PE,PFSP(1,I),PFSP(2,I),PFSP(3,I),
     &                  PFSP(4,I),IDFSP(I),IMODE)
         ELSEIF ((ICAS.EQ.2).AND.(IST.EQ.15)) THEN
* target cascade but fsp got stuck in proj. --> transform it into
* proj. rest system
            CALL LTRANS(PX,PY,PZ,PE,PFSP(1,I),PFSP(2,I),PFSP(3,I),
     &                  PFSP(4,I),IDFSP(I),-1)
         ELSEIF ((ICAS.EQ.1).AND.(IST.EQ.16)) THEN
* proj. cascade but fsp got stuck in target --> transform it into
* target rest system
            CALL LTRANS(PX,PY,PZ,PE,PFSP(1,I),PFSP(2,I),PFSP(3,I),
     &                  PFSP(4,I),IDFSP(I),1)
         ENDIF

* dump final state particles into HKKEVT
         IGEN = IDCH(IDXCAS)+1
         ID   = IPDGHA(IDFSP(I))
         IXR  = 0
         IF (LABSOR) IXR = 99
         CALL EVTPUT(IST,ID,IDXCAS,IDXSPE(1),PFSP(1,I),
     &               PFSP(2,I),PFSP(3,I),PFSP(4,I),0,IXR,IGEN)

* update the counter for particles which got stuck inside the nucleus
         IF ((IST.EQ.15).OR.(IST.EQ.16)) THEN
            NOINC = NOINC+1
            IDXINC(NOINC) = NHKK
         ENDIF
         IF (LABSOR) THEN
*   in case of absorption the spatial treatment is an approximate
*   solution anyway (the positions of the nucleons which "absorb" the
*   cascade particle are not taken into consideration) therefore the
*   particles are produced at the position of the cascade particle
            DO 12 K=1,4
               WHKK(K,NHKK) = WHKK(K,IDXCAS)
               VHKK(K,NHKK) = VHKK(K,IDXCAS)
   12       CONTINUE
         ELSE
*   DDISTL - distance the cascade particle moves to the intera. point
*   (the position where impact-parameter = distance to the interacting
*   nucleon), DIST - distance to the interacting nucleon at the time of
*   formation of the cascade particle, BINT - impact-parameter of this
*   cascade-interaction
            DDISTL = SQRT(DIST**2-BINT**2)
            DTIME  = DDISTL/BECAS(ICAS)
            DTIMEL = DDISTL/BGCAS(ICAS)
            RDISTL = DTIMEL*BGCAS(I2)
            IF ((IP.GT.1).AND.(IT.GT.1)) THEN
               RTIME = RDISTL/BECAS(I2)
            ELSE
               RTIME = ZERO
            ENDIF
*   RDISTL, RTIME are this step and time in the rest system of the other
*   nucleus
            DO 13 K=1,3
               VTXCA1(ICAS,K) = VTXCAS(ICAS,K)+COSCAS(ICAS,K)*DDISTL
               VTXCA1(I2,K)   = VTXCAS(I2,K)  +COSCAS(I2,K)  *RDISTL
   13       CONTINUE
            VTXCA1(ICAS,4) = VTXCAS(ICAS,4)+DTIME
            VTXCA1(I2,4)   = VTXCAS(I2,4)  +RTIME
*   position of particle production is half the impact-parameter to
*   the interacting nucleon
            DO 14 K=1,3
               WHKK(K,NHKK) = OHALF*(VTXCA1(1,K)+WHKK(K,IDXSPE(1)))
               VHKK(K,NHKK) = OHALF*(VTXCA1(2,K)+VHKK(K,IDXSPE(1)))
   14       CONTINUE
*   time of production of secondary = time of interaction
            WHKK(4,NHKK) = VTXCA1(1,4)
            VHKK(4,NHKK) = VTXCA1(2,4)
         ENDIF

   11 CONTINUE

* modify status and position of cascade particle (the latter for
* statistics reasons only)
      ISTHKK(IDXCAS) = 2
      IF (LABSOR) ISTHKK(IDXCAS) = 19
      IF (.NOT.LABSOR) THEN
         DO 15 K=1,4
            WHKK(K,IDXCAS) = VTXCA1(1,K)
            VHKK(K,IDXCAS) = VTXCA1(2,K)
   15    CONTINUE
      ENDIF

      DO 16 I=1,NSPE
         IS = IDXSPE(I)
* dump interacting nucleons for energy-momentum conservation check
         IF (LEMCCK)
     &      CALL EVTEMC(PHKK(1,IS),PHKK(2,IS),PHKK(3,IS),PHKK(4,IS),
     &                                                  2,IDUM,IDUM)
* modify entry for interacting nucleons
         IF (ISTHKK(IS).EQ.12+ICAS) ISTHKK(IS)=16+ICAS
         IF (ISTHKK(IS).EQ.14+ICAS) ISTHKK(IS)=2
         IF (I.GE.2) THEN
            JDAHKK(1,IS) = JDAHKK(1,IDXSPE(1))
            JDAHKK(2,IS) = JDAHKK(2,IDXSPE(1))
         ENDIF
   16 CONTINUE

* check energy-momentum conservation
      IF (LEMCCK) THEN
         CALL EVTEMC(DUM,DUM,DUM,DUM,4,500,IREJ1)
         IF (IREJ1.NE.0) GOTO 9999
      ENDIF

* update counter
      IF (LABSOR) THEN
         NINCCO(ICAS,1) = NINCCO(ICAS,1)+1
      ELSE
         IF (IPROC.EQ.1) NINCCO(ICAS,2) = NINCCO(ICAS,2)+1
         IF (IPROC.EQ.2) NINCCO(ICAS,3) = NINCCO(ICAS,3)+1
      ENDIF

      RETURN

 9997 CONTINUE
 9998 CONTINUE
* transport-step but no cascade step due to configuration (i.e. there
* is no nucleon for interaction etc.)
      IF (LCAS) THEN
         DO 100 K=1,4
C           WHKK(K,IDXCAS) = VTXCAS(1,K)
C           VHKK(K,IDXCAS) = VTXCAS(2,K)
            WHKK(K,IDXCAS) = VTXCA1(1,K)
            VHKK(K,IDXCAS) = VTXCA1(2,K)
  100    CONTINUE
      ENDIF

C9998 CONTINUE
* no cascade-step because of configuration
* (i.e. hadron outside nucleus etc.)
      LCAS = .TRUE.
      RETURN

 9999 CONTINUE
* rejection
      IREJ = 1
      RETURN
      END
*
*===absorp=============================================================*
*
      SUBROUTINE ABSORP(IDCAS,PCAS,NCAS,NSPE,IDSPE,IDXSPE,MODE,IREJ)

************************************************************************
* Two-nucleon absorption of antiprotons, pi-, and K-.                  *
* Antiproton absorption is handled by HADRIN.                          *
* The following channels for meson-absorption are considered:          *
*          pi- + p + p ---> n + p                                      *
*          pi- + p + n ---> n + n                                      *
*          K-  + p + p ---> sigma+ + n / Lam + p / sigma0 + p          *
*          K-  + p + n ---> sigma- + n / Lam + n / sigma0 + n          *
*          K-  + p + p ---> sigma- + n                                 *
*      IDCAS, PCAS   identity, momentum of particle to be absorbed     *
*      NCAS =  1     intranuclear cascade in projectile                *
*           = -1     intranuclear cascade in target                    *
*      NSPE          number of spectator nucleons involved             *
*      IDXSPE(2)     HKKEVT-indices of spectator nucleons involved     *
* Revised version of the original STOPIK written by HJM and J. Ranft.  *
* This version dated 11.06.96 is written by S. Roesler                 *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)
      PARAMETER (TINY10=1.0D-10,TINY5=1.0D-5,ONE=1.0D0,
     &           ONETHI=0.3333D0,TWOTHI=0.6666D0)

      PARAMETER (NMXHKK=89998)
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)
      COMMON /EXTEVT/ IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),
     &                IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10)
      LOGICAL LEMCCK,LHADRO,LSEADI
C                        j.r.3.10.96
      COMMON /FLAGS/  IFRAG(2),IRESCO,IMSHL,IRESRJ,
     &                LEMCCK,LHADRO(0:9),LSEADI
      PARAMETER (MAXFSP=10)
      COMMON /FISTAT/ PFSP(5,MAXFSP),IDFSP(MAXFSP),NFSP

      CHARACTER*8 ANAME
      COMMON /DPAR/   ANAME(210),AAM(210),GA(210),TAU(210),
     &                IICH(210),IIBAR(210),K1(210),K2(210)

      DIMENSION PCAS(5),IDXSPE(2),IDSPE(2),PSPE(2,5),PSPE1(5),
     &          PTOT3P(4),BG3P(4),
     &          ECMF(2),PCMF(2),CODF(2),COFF(2),SIFF(2)

      IREJ = 0
      NFSP = 0

* skip particles others than ap, pi-, K- for mode=0
      IF ((MODE.EQ.0).AND.
     &    (IDCAS.NE.2).AND.(IDCAS.NE.14).AND.(IDCAS.NE.16)) RETURN
* skip particles others than pions for mode=1
      IF ((MODE.EQ.1).AND.(IDCAS.NE.13).AND.
     &(IDCAS.NE.23).AND.(IDCAS.NE.14)) RETURN

      NUCAS = NCAS
      IF (NUCAS.EQ.-1) NUCAS = 2

      IF (MODE.EQ.0) THEN
* scan spectator nucleons for nucleons being able to "absorb"
         NSPE      = 0
         IDXSPE(1) = 0
         IDXSPE(2) = 0
         DO 1 I=1,NHKK
            IF ((ISTHKK(I).EQ.12+NUCAS).OR.(ISTHKK(I).EQ.14+NUCAS)) THEN
               NSPE         = NSPE+1
               IDXSPE(NSPE) = I
               IDSPE(NSPE)  = IDBAM(I)
               IF ((NSPE.EQ.1).AND.(IDCAS.EQ.2)) GOTO 2
               IF (NSPE.EQ.2) THEN
                  IF ((IDCAS.EQ.14).AND.(IDSPE(1).EQ.8).AND. 
     &                                  (IDSPE(2).EQ.8)) THEN
*    there is no pi-+n+n channel
                     NSPE = 1
                     GOTO 1
                  ELSE
                     GOTO 2
                  ENDIF
               ENDIF
            ENDIF
    1    CONTINUE

    2    CONTINUE
      ENDIF
* transform excited projectile nucleons (status=15) into proj. rest s.
      DO 3 I=1,NSPE
         DO 4 K=1,5
            PSPE(I,K) = PHKK(K,IDXSPE(I))
    4    CONTINUE
    3 CONTINUE

* antiproton absorption
      IF ((IDCAS.EQ.2).AND.(NSPE.GE.1)) THEN
         DO 5 K=1,5
            PSPE1(K) = PSPE(1,K)
    5    CONTINUE
**sr mod. for DPMJET: HADRIN-->HADRI1
         CALL HADRI1(IDCAS,PCAS,IDSPE(1),PSPE1,1,IREJ1)
         IF (IREJ1.NE.0) GOTO 9999

* meson absorption
      ELSEIF (((IDCAS.EQ.13).OR.(IDCAS.EQ.14).OR.
     &(IDCAS.EQ.23).OR.(IDCAS.EQ.16))
     &        .AND.(NSPE.GE.2)) THEN
         IF (IDCAS.EQ.14) THEN
*   pi- absorption
            IDFSP(1) = 8
            IDFSP(2) = 8
            IF ((IDSPE(1).EQ.1).AND.(IDSPE(2).EQ.1)) IDFSP(2) = 1
         ELSEIF (IDCAS.EQ.13) THEN
*   pi+ absorption
            IDFSP(1) = 1
            IDFSP(2) = 1
            IF ((IDSPE(1).EQ.8).AND.(IDSPE(2).EQ.8)) IDFSP(2) = 8
         ELSEIF (IDCAS.EQ.23) THEN
*   pi-0 absorption
            IDFSP(1) =IDSPE(1) 
            IDFSP(2) =IDSPE(2) 
         ELSEIF (IDCAS.EQ.16) THEN
*   K- absorption
            R = RNDM(V)
            IF ((IDSPE(1).EQ.1).AND.(IDSPE(2).EQ.1)) THEN
               IF (R.LT.ONETHI) THEN
                  IDFSP(1) = 21
                  IDFSP(2) = 8
               ELSEIF (R.LT.TWOTHI) THEN
                  IDFSP(1) = 17
                  IDFSP(2) = 1
               ELSE
                  IDFSP(1) = 22
                  IDFSP(2) = 1
               ENDIF
            ELSEIF ((IDSPE(1).EQ.8).AND.(IDSPE(2).EQ.8)) THEN
               IDFSP(1) = 20 
               IDFSP(2) = 8
            ELSE
               IF (R.LT.ONETHI) THEN
                  IDFSP(1) = 20
                  IDFSP(2) = 1
               ELSEIF (R.LT.TWOTHI) THEN
                  IDFSP(1) = 17
                  IDFSP(2) = 8
               ELSE
                  IDFSP(1) = 22
                  IDFSP(2) = 8
               ENDIF
            ENDIF
         ENDIF
*   dump initial particles for energy-momentum cons. check
         IF (LEMCCK) THEN
            CALL EVTEMC(PCAS(1),PCAS(2),PCAS(3),PCAS(4),1,IDUM,IDUM)
            CALL EVTEMC(PSPE(1,1),PSPE(1,2),PSPE(1,3),PSPE(1,4),2,  
     &                                                    IDUM,IDUM)
            CALL EVTEMC(PSPE(2,1),PSPE(2,2),PSPE(2,3),PSPE(2,4),2,  
     &                                                    IDUM,IDUM)
         ENDIF
*   get Lorentz-parameter of 3 particle initial state
         DO 6 K=1,4
            PTOT3P(K) = PCAS(K)+PSPE(1,K)+PSPE(2,K)
    6    CONTINUE
         P3P  = SQRT(PTOT3P(1)**2+PTOT3P(2)**2+PTOT3P(3)**2)
         AM3P = SQRT( (PTOT3P(4)-P3P)*(PTOT3P(4)+P3P) )
         DO 7 K=1,4
            BG3P(K) = PTOT3P(K)/MAX(AM3P,TINY10)
    7    CONTINUE
*   2-particle decay of the 3-particle compound system
         CALL DTWOPD(AM3P,ECMF(1),ECMF(2),PCMF(1),PCMF(2),
     &               CODF(1),COFF(1),SIFF(1),CODF(2),COFF(2),SIFF(2),
     &               AAM(IDFSP(1)),AAM(IDFSP(2)))
         DO 8 I=1,2
            SDF = SQRT((ONE-CODF(I))*(ONE+CODF(I)))
            PX  = PCMF(I)*COFF(I)*SDF
            PY  = PCMF(I)*SIFF(I)*SDF
            PZ  = PCMF(I)*CODF(I)
            CALL DALTRA(BG3P(4),BG3P(1),BG3P(2),BG3P(3),PX,PY,PZ,
     &                  ECMF(I),PTOFSP,PFSP(1,I),PFSP(2,I),PFSP(3,I),
     &                  PFSP(4,I))
            PFSP(5,I) = SQRT( (PFSP(4,I)-PTOFSP)*(PFSP(4,I)+PTOFSP) )
*   check consistency of kinematics
            IF (ABS(AAM(IDFSP(I))-PFSP(5,I)).GT.TINY5) THEN
               WRITE(LOUT,1001) IDFSP(I),AAM(IDFSP(I)),PFSP(5,I)
 1001          FORMAT(1X,'ABSORP:   warning! inconsistent',
     &                ' tree-particle kinematics',/,20X,'id: ',I3,
     &                ' AAM = ',E10.4,' MFSP = ',E10.4)
            ENDIF
*   dump final state particles for energy-momentum cons. check
            IF (LEMCCK) CALL EVTEMC(-PFSP(1,I),-PFSP(2,I),
     &                              -PFSP(3,I),-PFSP(4,I),2,IDUM,IDUM)
    8    CONTINUE
         NFSP = 2
         IF (LEMCCK) THEN
            CALL EVTEMC(DUM,DUM,DUM,DUM,3,100,IREJ1)
            IF (IREJ1.NE.0) THEN
               WRITE(LOUT,*)'ABSORB: EMC ',AAM(IDFSP(1)),AAM(IDFSP(2)),
     &                      AM3P
               GOTO 9999
            ENDIF
         ENDIF
      ELSE
C        IF (IOULEV(3).GT.0) WRITE(LOUT,1000) IDCAS,NSPE
 1000    FORMAT(1X,'ABSORP:   warning! absorption for particle ',I3,
     &          ' impossible',/,20X,'too few spectators (',I2,')')
         NSPE = 0
      ENDIF

      RETURN

 9999 CONTINUE
C     IF (IOULEV(1).GT.0) WRITE(LOUT,*) 'rejected 1 in ABSORP'
      IREJ = 1
      RETURN
      END
*
*===hadri1=============================================================*
*
**sr mod. for DPMJET: HADRIN-->HADRI1
      SUBROUTINE HADRI1(IDPR,PPR,IDTA,PTA,MODE,IREJ)

************************************************************************
* Interface to the HADRIN-routines for inelastic and elastic           *
* scattering.                                                          *
*      IDPR,PPR(5)   identity, momentum of projectile                  *
*      IDTA,PTA(5)   identity, momentum of target                      *
*      MODE  = 1     inelastic interaction                             *
*            = 2     elastic   interaction                             *
* Revised version of the original FHAD.                                *
* This version dated 27.10.95 is written by S. Roesler                 *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)
      PARAMETER (ZERO=0.0D0,TINY10=1.0D-10,TINY5=1.0D-5,TINY3=1.0D-3,
     &           TINY2=1.0D-2,TINY1=1.0D-1,ONE=1.0D0)

      LOGICAL LCORR,LMSSG
      LOGICAL LEMCCK,LHADRO,LSEADI
      COMMON /FLAGS/  IFRAG(2),IRESCO,IMSHL,IRESRJ,
     &                LEMCCK,LHADRO(0:9),LSEADI
      PARAMETER (MAXFSP=10)
      COMMON /FISTAT/ PFSP(5,MAXFSP),IDFSP(MAXFSP),NFSP

      CHARACTER*8 ANAME
      COMMON /DPAR/   ANAME(210),AAM(210),GA(210),TAU(210),
     &                IICH(210),IIBAR(210),K1(210),K2(210)

* output-common for DHADRI/ELHAIN
      PARAMETER (MAXFIN=10)
      COMMON /DFINLS/ ITRH(MAXFIN),CXRH(MAXFIN),CYRH(MAXFIN),
     &                CZRH(MAXFIN),ELRH(MAXFIN),PLRH(MAXFIN),IRH

      DIMENSION PPR(5),PPR1(5),PTA(5),BGTA(4),
     &          P1IN(4),P2IN(4),P1OUT(4),P2OUT(4),IMCORR(2)

      DATA LMSSG /.TRUE./

      IREJ  = 0
      NFSP  = 0
      KCORR = 0
      IMCORR(1) = 0
      IMCORR(2) = 0
      LCORR = .FALSE.

*   dump initial particles for energy-momentum cons. check
      IF (LEMCCK) THEN
         CALL EVTEMC(PPR(1),PPR(2),PPR(3),PPR(4),1,IDUM,IDUM)
         CALL EVTEMC(PTA(1),PTA(2),PTA(3),PTA(4),2,IDUM,IDUM)
      ENDIF

      AMP2 = PPR(4)**2-PPR(1)**2-PPR(2)**2-PPR(3)**2
      AMT2 = PTA(4)**2-PTA(1)**2-PTA(2)**2-PTA(3)**2
      IF ((AMP2.LT.ZERO).OR.(AMT2.LT.ZERO).OR.
     &    (ABS(AMP2-AAM(IDPR)**2).GT.TINY5).OR.
     &    (ABS(AMT2-AAM(IDTA)**2).GT.TINY5)) THEN
         IF (LMSSG)
     &   WRITE(LOUT,1000) AMP2,AAM(IDPR)**2,AMT2,AAM(IDTA)**2
 1000    FORMAT(1X,'HADRIN:   warning! inconsistent projectile/target',
     &          ' mass',/,20X,'AMP2 = ',E15.7,', AAM(IDPR)**2 = ',
     &          E15.7,/,20X,'AMT2 = ',E15.7,', AAM(IDTA)**2 = ',E15.7)
         LMSSG = .FALSE.
         LCORR = .TRUE.
      ENDIF

* convert initial state particles into particles which can be
* handled by HADRIN
      IDHPR = IDPR
      IDHTA = IDTA
      IF ((IDHPR.LE.0).OR.(IDHPR.GE.111).OR.(LCORR)) THEN
         IF ((IDHPR.LE.0).OR.(IDHPR.GE.111)) IDHPR = 1
         DO 1 K=1,4
            P1IN(K) = PPR(K)
            P2IN(K) = PTA(K)
    1    CONTINUE
         XM1 = AAM(IDHPR)
         XM2 = AAM(IDHTA)
         CALL MASHEL(P1IN,P2IN,XM1,XM2,P1OUT,P2OUT,IREJ1)
         IF (IREJ1.GT.0) THEN
C           WRITE(LOUT,'(1X,A)') 'HADRIN:   inconsistent mass trsf.'
            GOTO 9999
         ENDIF
         DO 2 K=1,4
            PPR(K) = P1OUT(K)
            PTA(K) = P2OUT(K)
    2    CONTINUE
         PPR(5) = SQRT(PPR(4)**2-PPR(1)**2-PPR(2)**2-PPR(3)**2)
         PTA(5) = SQRT(PTA(4)**2-PTA(1)**2-PTA(2)**2-PTA(3)**2)
      ENDIF

* Lorentz-parameter for trafo into rest-system of target
      DO 3 K=1,4
         BGTA(K) = PTA(K)/PTA(5)
    3 CONTINUE
* transformation of projectile into rest-system of target
      CALL DALTRA(BGTA(4),-BGTA(1),-BGTA(2),-BGTA(3),PPR(1),PPR(2),
     &            PPR(3),PPR(4),PPRTO1,PPR1(1),PPR1(2),PPR1(3),
     &            PPR1(4))

* direction cosines of projectile in target rest system
      CX = PPR1(1)/PPRTO1
      CY = PPR1(2)/PPRTO1
      CZ = PPR1(3)/PPRTO1

* sample inelastic interaction
      IF (MODE.EQ.1) THEN
         CALL DHADRI(IDHPR,PPRTO1,PPR1(4),CX,CY,CZ,IDHTA)
         IF (IRH.EQ.1)THEN
C	   WRITE(6,'(A)')' DHADRI Rej'
	   GOTO 9998
	 ENDIF
* sample elastic interaction
      ELSEIF (MODE.EQ.2) THEN
         CALL ELHAIN(IDHPR,PPRTO1,PPR1(4),CX,CY,CZ,IDHTA,IREJ1)
         IF (IREJ1.NE.0) THEN
C           WRITE(LOUT,*) 'rejected 1 in HADRIN'
            GOTO 9999
         ENDIF
         IF (IRH.EQ.1) GOTO 9998
      ELSE
         WRITE(LOUT,1001) MODE,INTHAD
 1001    FORMAT(1X,'HADRIN:   warning! inconsistent interaction mode',
     &          I4,' (INTHAD =',I4,')')
         GOTO 9999
      ENDIF

* transform final state particles back into Lab.
      DO 4 I=1,IRH
         NFSP = NFSP+1
         PX   = CXRH(I)*PLRH(I)
         PY   = CYRH(I)*PLRH(I)
         PZ   = CZRH(I)*PLRH(I)
         CALL DALTRA(BGTA(4),BGTA(1),BGTA(2),BGTA(3),PX,PY,PZ,ELRH(I),
     &               PTOFSP,PFSP(1,NFSP),PFSP(2,NFSP),PFSP(3,NFSP),
     &               PFSP(4,NFSP))
         IDFSP(NFSP) = ITRH(I)
         AMFSP2 = PFSP(4,NFSP)**2-PFSP(1,NFSP)**2-PFSP(2,NFSP)**2-
     &                                            PFSP(3,NFSP)**2
         IF (AMFSP2.LT.-TINY3) THEN
            WRITE(LOUT,1002) IDFSP(NFSP),PFSP(1,NFSP),PFSP(2,NFSP),
     &                       PFSP(3,NFSP),PFSP(4,NFSP),AMFSP2
 1002       FORMAT(1X,'HADRIN:   warning! final state particle (id = ',
     &             I2,') with negative mass^2',/,1X,5E12.4)
            GOTO 9999
         ELSE
            PFSP(5,NFSP) = SQRT(ABS(AMFSP2))
            IF (ABS(PFSP(5,NFSP)-AAM(IDFSP(NFSP))).GT.TINY1) THEN
C              WRITE(LOUT,1003) IDFSP(NFSP),AAM(IDFSP(NFSP)),
C    &                          PFSP(5,NFSP)
 1003          FORMAT(1X,'HADRIN:   warning! final state particle',
     &                ' (id = ',I2,') with inconsistent mass',/,1X,
     &                2E12.4)
               KCORR         = KCORR+1
               IF (KCORR.GT.2) GOTO 9999
               IMCORR(KCORR) = NFSP
            ENDIF
         ENDIF
*   dump final state particles for energy-momentum cons. check
         IF (LEMCCK) CALL EVTEMC(-PFSP(1,I),-PFSP(2,I),
     &                           -PFSP(3,I),-PFSP(4,I),2,IDUM,IDUM)
    4 CONTINUE

* transform momenta on mass shell in case of inconsistencies in
* HADRIN
      IF (KCORR.GT.0) THEN
         IF (KCORR.EQ.2) THEN
            I1 = IMCORR(1)
            I2 = IMCORR(2)
         ELSE
            IF (IMCORR(1).EQ.1) THEN
               I1 = 1
               I2 = 2
            ELSE
               I1 = 1
               I2 = IMCORR(1)
            ENDIF
         ENDIF
         IF (LEMCCK) CALL EVTEMC(PFSP(1,I1),PFSP(2,I1),
     &                           PFSP(3,I1),PFSP(4,I1),2,IDUM,IDUM)
         IF (LEMCCK) CALL EVTEMC(PFSP(1,I2),PFSP(2,I2),
     &                           PFSP(3,I2),PFSP(4,I2),2,IDUM,IDUM)
         DO 5 K=1,4
            P1IN(K) = PFSP(K,I1)
            P2IN(K) = PFSP(K,I2)
    5    CONTINUE
         XM1 = AAM(IDFSP(I1))
         XM2 = AAM(IDFSP(I2))
         CALL MASHEL(P1IN,P2IN,XM1,XM2,P1OUT,P2OUT,IREJ1)
         IF (IREJ1.GT.0) THEN
C           WRITE(LOUT,'(1X,A)') 'HADRIN:   inconsistent mass trsf.'
            GOTO 9999
         ENDIF
         DO 6 K=1,4
            PFSP(K,I1) = P1OUT(K)
            PFSP(K,I2) = P2OUT(K)
    6    CONTINUE
         PFSP(5,I1) = SQRT(PFSP(4,I1)**2-PFSP(1,I1)**2
     &                    -PFSP(2,I1)**2-PFSP(3,I1)**2)
         PFSP(5,I2) = SQRT(PFSP(4,I2)**2-PFSP(1,I2)**2
     &                    -PFSP(2,I2)**2-PFSP(3,I2)**2)
*   dump final state particles for energy-momentum cons. check
         IF (LEMCCK) CALL EVTEMC(-PFSP(1,I1),-PFSP(2,I1),
     &                           -PFSP(3,I1),-PFSP(4,I1),2,IDUM,IDUM)
         IF (LEMCCK) CALL EVTEMC(-PFSP(1,I2),-PFSP(2,I2),
     &                           -PFSP(3,I2),-PFSP(4,I2),2,IDUM,IDUM)
      ENDIF

* check energy-momentum conservation
      IF (LEMCCK) THEN
         CALL EVTEMC(DUM,DUM,DUM,DUM,4,102,IREJ1)
         IF (IREJ1.NE.0)THEN 
C	   WRITE(6,'(A)')' EVTEMC-HADRIN Rej'
	   GOTO 9999
         ENDIF
      ENDIF

      RETURN

 9998 CONTINUE
      IREJ = 2
      RETURN

 9999 CONTINUE
      IREJ = 1
      RETURN
      END
*
*===evtput=============================================================*
*
      SUBROUTINE EVTPUT(IST,ID,M1,M2,PX,PY,PZ,E,IDR,IDXR,IDC)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)
      PARAMETER (TINY10=1.0D-10,TINY4=1.0D-4,TINY3=1.0D-3,
     &           TINY2=1.0D-2,SQTINF=1.0D+15,ZERO=0.D0)

      PARAMETER (NMXHKK=89998)
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)
      COMMON /EXTEVT/ IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),
     &                IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10)

      COMMON /NUCCMS/ GACMS,BGCMS,GALAB,BGLAB,BLAB,UMO,PCM,EPROJ,PPROJ
      CHARACTER*8  ANAME
      COMMON /DPAR/   ANAME(210),AAM(210),GA(210),TAU(210),
     &                IICH(210),IIBAR(210),K1(210),K2(210)
C     WRITE(6,'(A,4I5,4F10.3,3I5)')
C    &' EVTPUT, IST,ID,M1,M2,PX,PY,PZ,E,IDR,IDXR,IDC',
C    & IST,ID,M1,M2,PX,PY,PZ,E,IDR,IDXR,IDC

C     IF (MODE.GT.100) THEN
C        WRITE(LOUT,'(1X,A,I5,A,I5)')
C    &        'EVTPUT: reset NHKK = ',NHKK,' to NHKK =',NHKK-MODE+100
C        NHKK = NHKK-MODE+100
C        RETURN
C     ENDIF
      MO1  = M1
      MO2  = M2
      NHKK = NHKK+1

      IF (NHKK.GT.NMXHKK) THEN
         WRITE(LOUT,1000) NHKK
 1000    FORMAT(1X,'EVTPUT: NHKK exeeds NMXHKK = ',I7,
     &             '! program execution stopped..')
         STOP
      ENDIF
      IF (M1.LT.0) MO1 = NHKK+M1
      IF (M2.LT.0) MO2 = NHKK+M2
      ISTHKK(NHKK)   = IST
      IDHKK(NHKK)    = ID      
      JMOHKK(1,NHKK) = MO1
      JMOHKK(2,NHKK) = MO2
      JDAHKK(1,NHKK) = 0
      JDAHKK(2,NHKK) = 0
      IDRES(NHKK)    = IDR
      IDXRES(NHKK)   = IDXR
      IDCH(NHKK)     = IDC
      IF (ID.EQ.88888.OR.ID.EQ.88887.OR.ID.EQ.88889) THEN
         IDMO1 = ABS(IDHKK(MO1))
         IDMO2 = ABS(IDHKK(MO2))
         IF ((IDMO1.LT.100).AND.(IDMO2.LT.100)) NOBAM(NHKK) = 3
         IF ((IDMO1.LT.100).AND.(IDMO2.GT.100)) NOBAM(NHKK) = 4
         IF ((IDMO1.GT.100).AND.(IDMO2.GT.100)) NOBAM(NHKK) = 5
         IF ((IDMO1.GT.100).AND.(IDMO2.LT.100)) NOBAM(NHKK) = 6
      ELSE
         NOBAM(NHKK) = 0
      ENDIF
      IDBAM(NHKK) = ICIHAD(ID)
      IF (MO1.GT.0) THEN
         IF (JDAHKK(1,MO1).NE.0) THEN
            JDAHKK(2,MO1) = NHKK
         ELSE
            JDAHKK(1,MO1) = NHKK
         ENDIF
      ENDIF
      IF (MO2.GT.0) THEN
         IF (JDAHKK(1,MO2).NE.0) THEN
            JDAHKK(2,MO2) = NHKK
         ELSE
            JDAHKK(1,MO2) = NHKK
         ENDIF
      ENDIF
C     WRITE(6,'(A,2I10)')' EVTPUT:NHKK,IDBAM(NHKK)',NHKK,IDBAM(NHKK)
      IF(IDBAM(NHKK).EQ.410)IDBAM(NHKK)=210
      IF (IDBAM(NHKK).GT.0) THEN
         PTOT   = SQRT(PX**2+PY**2+PZ**2)
         AM0    = SQRT(ABS( (E-PTOT)*(E+PTOT) ))
         AMRQ   = AAM(IDBAM(NHKK))
         AMDIF2 = (AM0-AMRQ)*(AM0+AMRQ)
         IF ((ABS(AMDIF2).GT.TINY3).AND.(E.LT.SQTINF).AND.
     &       (PTOT.GT.ZERO)) THEN
            DELTA = -AMDIF2/(2.0D0*(E+PTOT))
C           DELTA = (AMRQ2-AM2)/(2.0D0*(E+PTOT))
            E     = E+DELTA
            PTOT1 = PTOT-DELTA
            PX    = PX*PTOT1/PTOT
            PY    = PY*PTOT1/PTOT
            PZ    = PZ*PTOT1/PTOT
         ENDIF
      ENDIF
      PHKK(1,NHKK) = PX
      PHKK(2,NHKK) = PY
      PHKK(3,NHKK) = PZ
      PHKK(4,NHKK) = E
      PTOT = SQRT( PX**2+PY**2+PZ**2 )
      PHKK(5,NHKK) = (PHKK(4,NHKK)-PTOT)*(PHKK(4,NHKK)+PTOT)
C     IF ((PHKK(5,NHKK).LT.0.0D0).AND.(ABS(PHKK(5,NHKK)).GT.TINY4))
C    &   WRITE(LOUT,'(1X,A,G10.3)')
C    &     'EVTPUT: negative mass**2 ',PHKK(5,NHKK)
      PHKK(5,NHKK) = SQRT(ABS(PHKK(5,NHKK)))
C     IF (ID.EQ.88888) THEN
      IF (ID.EQ.88888.OR.ID.EQ.88887.OR.ID.EQ.88889) THEN
* special treatment for chains:
*    z coordinate of chain in Lab  = pos. of target nucleon
*    time of chain-creation in Lab = time of passage of projectile
*                                    nucleus at pos. of taget nucleus
C        VHKK(1,NHKK) = 0.5D0*(VHKK(1,MO1)+VHKK(1,MO2))
C        VHKK(2,NHKK) = 0.5D0*(VHKK(2,MO1)+VHKK(2,MO2))
         VHKK(1,NHKK) = VHKK(1,MO2)
         VHKK(2,NHKK) = VHKK(2,MO2)
         VHKK(3,NHKK) = VHKK(3,MO2)
         VHKK(4,NHKK) = VHKK(3,MO2)/BLAB-VHKK(3,MO1)/BGLAB
C        WHKK(1,NHKK) = 0.5D0*(WHKK(1,MO1)+WHKK(1,MO2))
C        WHKK(2,NHKK) = 0.5D0*(WHKK(2,MO1)+WHKK(2,MO2))
         WHKK(1,NHKK) = WHKK(1,MO1)
         WHKK(2,NHKK) = WHKK(2,MO1)
         WHKK(3,NHKK) = WHKK(3,MO1)
         WHKK(4,NHKK) = -WHKK(3,MO1)/BLAB+WHKK(3,MO2)/BGLAB
      ELSE
         DO 2 I=1,4
            VHKK(I,NHKK) = VHKK(I,MO1)
            WHKK(I,NHKK) = WHKK(I,MO1)
    2    CONTINUE
      ENDIF

      RETURN
      END
*
*===mashel=============================================================*
*
      SUBROUTINE MASHEL(PA1,PA2,XM1,XM2,P1,P2,IREJ)

************************************************************************
*                                                                      *
*    rescaling of momenta of two partons to put both                   *
*                                       on mass shell                  *
*                                                                      *
*    input:       PA1,PA2   input momentum vectors                     *
*                 XM1,2     desired masses of particles afterwards     *
*                 P1,P2     changed momentum vectors                   *
*                                                                      *
* The original version is written by R. Engel.                         *
* This version dated 19.11.95 is modified by S. Roesler.               *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)
      PARAMETER (TINY10=1.0D-10,ONE=1.0D0,ZERO=0.0D0)

      DIMENSION PA1(4),PA2(4),P1(4),P2(4)

      IREJ = 0

* Lorentz transformation into system CMS
      PX  = PA1(1)+PA2(1)
      PY  = PA1(2)+PA2(2)
      PZ  = PA1(3)+PA2(3)
      EE  = PA1(4)+PA2(4)
      XPTOT = SQRT(PX**2+PY**2+PZ**2)
      XMS   = (EE-XPTOT)*(EE+XPTOT)
      IF(XMS.LT.(XM1+XM2)**2) THEN
C        WRITE(LOUT,'(A,3E12.4)')' MASHEL Rej',XMS,XM1,XM2
         GOTO 9999
      ENDIF
      XMS = SQRT(XMS)
      BGX = PX/XMS
      BGY = PY/XMS
      BGZ = PZ/XMS
      GAM = EE/XMS
      CALL DALTRA(GAM,-BGX,-BGY,-BGZ,PA1(1),PA1(2),PA1(3),
     &           PA1(4),PTOT1,P1(1),P1(2),P1(3),P1(4))
* rotation angles
      COD = P1(3)/PTOT1
      SID = SQRT((ONE-COD)*(ONE+COD))
      COF = ONE
      SIF = ZERO
      IF(PTOT1*SID.GT.TINY10) THEN
         COF   = P1(1)/(SID*PTOT1)
         SIF   = P1(2)/(SID*PTOT1)
         ANORF = SQRT(COF*COF+SIF*SIF)
         COF   = COF/ANORF
         SIF   = SIF/ANORF
      ENDIF
* new CM momentum and energies (for masses XM1,XM2)
      XM12 = XM1**2
      XM22 = XM2**2
      SS   = XMS**2
      PCMP = YLAMB(SS,XM12,XM22)/(2.D0*XMS)
      EE1  = SQRT(XM12+PCMP**2)
      EE2  = XMS-EE1
* back rotation
      MODE = 1
      CALL MYTRAN(MODE,ZERO,ZERO,PCMP,COD,SID,COF,SIF,XX,YY,ZZ)
      CALL DALTRA(GAM,BGX,BGY,BGZ,XX,YY,ZZ,EE1,
     &            PTOT1,P1(1),P1(2),P1(3),P1(4))
      CALL DALTRA(GAM,BGX,BGY,BGZ,-XX,-YY,-ZZ,EE2,
     &            PTOT2,P2(1),P2(2),P2(3),P2(4))
* check consistency
      DEL = XMS*0.0001D0
      IF (ABS(PX-P1(1)-P2(1)).GT.DEL) THEN
        IDEV = 1
      ELSEIF (ABS(PY-P1(2)-P2(2)).GT.DEL) THEN
        IDEV = 2
      ELSEIF (ABS(PZ-P1(3)-P2(3)).GT.DEL) THEN
        IDEV = 3
      ELSEIF (ABS(EE-P1(4)-P2(4)).GT.DEL) THEN
        IDEV = 4
      ELSE
        IDEV = 0
      ENDIF
      IF (IDEV.NE.0) THEN
         WRITE(LOUT,'(/1X,A,I3)')
     &      'MASHEL: inconsistent transformation',IDEV
         WRITE(LOUT,'(1X,A)') 'MASHEL: input momenta/masses:'
         WRITE(LOUT,'(1X,5E12.5)') (PA1(K),K=1,4),XM1
         WRITE(LOUT,'(1X,5E12.5)') (PA2(K),K=1,4),XM2
         WRITE(LOUT,'(1X,A)') 'MASHEL: output momenta:'
         WRITE(LOUT,'(5X,4E12.5)') (P1(K),K=1,4)
         WRITE(LOUT,'(5X,4E12.5)') (P2(K),K=1,4)
      ENDIF
      RETURN

 9999 CONTINUE
      IREJ = 1
      RETURN
      END
*
*===mytran=============================================================*
*
      SUBROUTINE MYTRAN(IMODE,XO,YO,ZO,CDE,SDE,CFE,SFE,X,Y,Z)

************************************************************************
* This subroutine rotates the coordinate frame                         *
*    a) theta  around y                                                *
*    b) phi    around z      if IMODE = 1                              *
*                                                                      *
*     x'          cos(ph) -sin(ph) 0      cos(th)  0  sin(th)   x      *
*     y' = A B =  sin(ph) cos(ph)  0  .   0        1        0   y      *
*     z'          0       0        1     -sin(th)  0  cos(th)   z      *
*                                                                      *
* and vice versa if IMODE = 0.                                         *
* This version dated 5.4.94 is based on the original version DTRAN     *
* by J. Ranft and is written by S. Roesler.                            *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)

      IF (IMODE.EQ.1) THEN
         X= CDE*CFE*XO-SFE*YO+SDE*CFE*ZO
         Y= CDE*SFE*XO+CFE*YO+SDE*SFE*ZO
         Z=-SDE    *XO       +CDE    *ZO
      ELSE
         X= CDE*CFE*XO+CDE*SFE*YO-SDE*ZO
         Y= -SFE*XO+CFE*YO
         Z= SDE*CFE*XO+SDE*SFE*YO+CDE*ZO
      ENDIF
      RETURN
      END
*
*===ylamb==============================================================*
*
      DOUBLE PRECISION FUNCTION YLAMB(X,Y,Z)

************************************************************************
*                                                                      *
*     auxiliary function for three particle decay mode                 *
*     (standard LAMBDA**(1/2) function)                                *
*                                                                      *
* Adopted from an original version written by R. Engel.                *
* This version dated 12.12.94 is written by S. Roesler.                *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      YZ   = Y-Z
      XLAM = X*X-2.D0*X*(Y+Z)+YZ*YZ
      IF (XLAM.LE.0.D0) XLAM = ABS(XLAM)
      YLAMB = SQRT(XLAM)

      RETURN
      END
*
*===evtemc=============================================================*
*
      SUBROUTINE EVTEMC(PXIO,PYIO,PZIO,EIO,IMODE,IPOS,IREJ)

************************************************************************
* This version dated 19.11.94 is written by S. Roesler                 *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)
      PARAMETER (TINY1=1.0D-1,TINY2=1.0D-2,TINY4=1.0D-4,TINY10=1.0D-10,
     &           ZERO=0.0D0,TINY11=300.D0)

      PARAMETER (NMXHKK=89998)
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)
      LOGICAL LEMCCK,LHADRO,LSEADI
      COMMON /FLAGS/  IFRAG(2),IRESCO,IMSHL,IRESRJ,
     &                LEMCCK,LHADRO(0:9),LSEADI
      COMMON /TMPEMC/ PX,PY,PZ,E
      COMMON /NUCC/   IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG
      DATA INII/0/

      IREJ = 0

      MODE = IMODE
      CHKLEV = TINY10
      CHKLXV=TINY11
      IF (MODE.EQ.4) THEN
         CHKLEV = TINY2
         MODE   = 3
      ELSEIF (MODE.EQ.5) THEN
         CHKLEV = TINY1
         MODE   = 3
      ELSEIF (MODE.EQ.-1) THEN
         CHKLEV = EIO
         MODE   = 3
**sr mod. for DPMJET: set check-level to some fixed value
*                     i.e. final state momentum is allowed to differ
*                     from the inital one by 50GeV (!!!)
C                     This was necessary to see wether the old
C                     version would work at all at high energy
C                     but it did not!
         CHKLXV = TINY11
         CHKLEV = CHKLXV
**
      ENDIF

      IF (ABS(MODE).EQ.3) THEN
         PXDEV = PX
         PYDEV = PY
         PZDEV = PZ
         EDEV  = E
         IF ((IFRAG(1).EQ.2).AND.(CHKLEV.LT.TINY4)) CHKLEV = TINY4
**sr mod. for DPMJET: use DPMJET check-level
	 IF ( IT.GE.200.AND.IP.GE.200)GO TO 9998
         IF ((ABS(PXDEV).GT.CHKLXV).OR.(ABS(PYDEV).GT.CHKLXV).OR.
     &       (ABS(PZDEV).GT.CHKLXV).OR.(ABS(EDEV).GT.CHKLXV)) THEN
**
	    INII=INII+1
	    IF(INII.LE.10)THEN
            WRITE(LOUT,'(1X,A,I4,A,I6,A,/,4G10.3)')
     &         'EVTEMC: energy-momentum cons. failure at pos. ',IPOS,
     &         '  event  ',NEVHKK,
     &         ' ! ',PXDEV,PYDEV,PZDEV,EDEV
**sr mod. for DPMJET: additional output
      WRITE(6,'(A/4E12.3,3I5)')
     * ' Input values (PXIO,PYIO,PZIO,EIO,IMODE,IPOS,IREJ)',
     * PXIO,PYIO,PZIO,EIO,IMODE,IPOS,IREJ
      WRITE(6,'(A/4E12.3)')
     * ' Input values in /TMPEMC/ (PX,PY,PZ,E)',
     * PX,PY,PZ,E
       ENDIF
**
            PX   = 0.0D0
            PY   = 0.0D0
            PZ   = 0.0D0
            E    = 0.0D0
            GOTO 9999
         ENDIF
 9998    CONTINUE
         PX   = 0.0D0
         PY   = 0.0D0
         PZ   = 0.0D0
         E    = 0.0D0
         RETURN
      ENDIF

      IF (MODE.EQ.1) THEN
         PX = 0.0D0
         PY = 0.0D0
         PZ = 0.0D0
         E  = 0.0D0
      ENDIF

      PX = PX+PXIO
      PY = PY+PYIO
      PZ = PZ+PZIO
      E  = E+EIO

      RETURN

 9999 CONTINUE
      IREJ = 1
      RETURN
      END
*
*===ltrans=============================================================*
*
      SUBROUTINE LTRANS(PXI,PYI,PZI,PEI,PXO,PYO,PZO,PEO,ID,MODE)
 
************************************************************************
* Special Lorentz-transformations.                                     *
*   MODE = 1(-1)    projectile rest syst.   --> Lab (back)             *
*        = 2(-2)    projectile rest syst.   --> nucl.-nucl.cms (back)  *
*        = 3(-3)    target rest syst. (=Lab)--> nucl.-nucl.cms (back)  *
* This version dated 01.11.95 is written by  S. Roesler.               *
************************************************************************
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)
      PARAMETER (TINY3=1.0D-3,ZERO=0.0D0,TWO=2.0D0)
 
      PARAMETER (SQTINF=1.0D+15)
 
      CHARACTER*8  ANAME
      COMMON /DPAR/   ANAME(210),AAM(210),GA(210),TAU(210),
     &                IICH(210),IIBAR(210),K1(210),K2(210)
 
      PXO = PXI
      PYO = PYI
      CALL LTNUC(PZI,PEI,PZO,PEO,MODE)
 
* check particle mass for consistency (numerical rounding errors)
      PO     = SQRT(PXO**2+PYO**2+PZO**2)
      AMO2   = (PEO-PO)*(PEO+PO)
      AMORQ2 = AAM(ID)**2
      AMDIF2 = ABS(AMO2-AMORQ2)
      IF ((AMDIF2.GT.TINY3).AND.(PEO.LT.SQTINF).AND.(PO.GT.ZERO)) THEN
         DELTA = (AMORQ2-AMO2)/(TWO*(PEO+PO))
         PEO   = PEO+DELTA
         PO1   = PO -DELTA
         PXO   = PXO*PO1/PO
         PYO   = PYO*PO1/PO
         PZO   = PZO*PO1/PO
      ENDIF
 
      RETURN
      END
*
*===ltnuc==============================================================*
*
      SUBROUTINE LTNUC(PIN,EIN,POUT,EOUT,MODE)

************************************************************************
* Lorentz-transformations.                                             *
*   PIN        longitudnal momentum       (input)                      *
*   EIN        energy                     (input)                      *
*   POUT       transformed long. momentum (output)                     *
*   EOUT       transformed energy         (output)                     *
*   MODE = 1(-1)    projectile rest syst.   --> Lab (back)             *
*        = 2(-2)    projectile rest syst.   --> nucl.-nucl.cms (back)  *
*        = 3(-3)    target rest syst. (=Lab)--> nucl.-nucl.cms (back)  *
* This version dated 01.11.95 is written by  S. Roesler.               *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)
      PARAMETER (ZERO=0.0D0)

      COMMON /NUCCMS/ GACMS,BGCMS,GALAB,BGLAB,BLAB,ECM,PCM,EPROJ,PPROJ

      IF (ABS(MODE).EQ.1) THEN
         BG = -SIGN(BGLAB,DBLE(MODE))
         CALL DALTRA(GALAB,ZERO,ZERO,-BG,ZERO,ZERO,PIN,EIN,
     &                               DUM,DUM,DUM,POUT,EOUT)
      ELSEIF (ABS(MODE).EQ.2) THEN
         BG = SIGN(BGCMS,DBLE(MODE))
         CALL DALTRA(GACMS,ZERO,ZERO,BG,ZERO,ZERO,PIN,EIN,
     &                                 DUM,DUM,DUM,POUT,EOUT)
      ELSEIF (ABS(MODE).EQ.3) THEN
         BG = -SIGN(BGCMS,DBLE(MODE))
         CALL DALTRA(GACMS,ZERO,ZERO,BG,ZERO,ZERO,PIN,EIN,
     &                               DUM,DUM,DUM,POUT,EOUT)
      ELSE
         WRITE(LOUT,1000) MODE
 1000    FORMAT(1X,'LTNUC: not supported mode (MODE = ',I3,')')
         EOUT = EIN
         POUT = PIN
      ENDIF

      RETURN
      END
**sr mod. for DPMJET: short version of the original DTUNUC-routine
*
*===evtini=============================================================*
*
      SUBROUTINE EVTINI(ID,IP,IT,EPN,PPN,ECM,NHKKH1,MODE)

************************************************************************
* Initialization of HKKEVT.                                            *
* This version dated 19.11.95 is written by S. Roesler                 *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)

      PARAMETER (NMXHKK=89998)
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)
      COMMON /EXTEVT/ IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),
     &                IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10)
      COMMON /NSTARI/NSTART

      GOTO (1,2) MODE

    1 CONTINUE
* initialization of EXTEVT
      DO 10 I=1,NHKK
         IDRES(I)  = 0
         IDXRES(I) = 0
         NOBAM(I)  = 0
         IDCH(I)   = 0
   10 CONTINUE
      CALL LTINI(ID,EPN,PPN,ECM)
C        IF(NSTART.NE.2.AND.NEUDEC.GE.20)
C    &	 CALL LTINI(IJPROJ,EPNI,DUM1,DUM2)

      RETURN 

    2 CONTINUE
      DO 20 I=1,NHKK
* get BAMJET-index of final state particle 
         IDBAM(I) = MCIHAD(IDHKK(I))
   20 CONTINUE
      NPOINT(1) = IP+IT+1
      NPOINT(4) = NHKKH1+1

      RETURN
      END
*
*===ltini==============================================================*
*
      SUBROUTINE LTINI(IDP,EPN,PPN,ECM)

************************************************************************
* Initializations of Lorentz-transformations, calculation of Lorentz-  *
* parameters.                                                          *
* This version dated 13.11.95 is written by  S. Roesler.               *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)
      PARAMETER (TINY3=1.0D-3,ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)

      COMMON /NUCCMS/ GACMS,BGCMS,GALAB,BGLAB,BLAB,UMO,PCM,EPROJ,PPROJ
**sr mod. for DPMJET: common added
      COMMON /TRAFOP/ GAMP,BGAMP,BETP
**
      CHARACTER*8  ANAME
      COMMON /DPAR/   ANAME(210),AAM(210),GA(210),TAU(210),
     &                IICH(210),IIBAR(210),K1(210),K2(210)

**sr mod. for DPMJET: force calulation starting from EPN
      ECM = ZERO
      PPN = ZERO
**
      IF (ECM.GT.ZERO) THEN
         EPN = (ECM**2-AAM(IDP)**2-AAM(1)**2)/(2.0D0*AAM(1))
         PPN = SQRT((EPN-AAM(IDP))*(EPN+AAM(IDP)))
      ELSE
         IF ((EPN.NE.ZERO).AND.(PPN.EQ.ZERO)) THEN
            IF (EPN.LT.ZERO) EPN = ABS(EPN)+AAM(IDP)
            PPN = SQRT((EPN-AAM(IDP))*(EPN+AAM(IDP)))
         ELSEIF ((PPN.GT.ZERO).AND.(EPN.EQ.ZERO)) THEN
            EPN = PPN*SQRT(ONE+(AAM(IDP)/PPN)**2)
         ENDIF
         ECM = SQRT(AAM(IDP)**2+AAM(1)**2+2.0D0*AAM(1)*EPN)
      ENDIF
      UMO   = ECM
      EPROJ = EPN
      PPROJ = PPN
*   Lorentz-parameter for transformation Lab. - projectile rest system
      IF(AAM(IDP).GT.0.D0)THEN
      GALAB = EPROJ/AAM(IDP)
      BGLAB = PPROJ/AAM(IDP)
      ELSE
      GALAB = EPROJ/(AAM(IDP)+0.0001D0)
      BGLAB = PPROJ/(AAM(IDP)+0.0001D0)
      ENDIF
      BLAB  = BGLAB/GALAB
*   Lorentz-parameter for transformation Lab. - nucl.-nucl. cms.
      GACMS = (EPROJ+AAM(1))/UMO
      BGCMS = PPROJ/UMO
      PCM   = GACMS*PPROJ-BGCMS*EPROJ
**sr mod. for DPMJET: initialize /TRAFOP/
      GAMP  = GALAB
      BGAMP = BGLAB
      BETP  = BGAMP/GAMP
**
C     WRITE(6,*)
C    &'IDP,EPN,PPN,ECM',IDP,EPN,PPN,ECM 
C     WRITE(6,*)
C    &'GACMS,BGCMS,GALAB,BGLAB,BLAB,UMO,PCM,EPROJ,PPROJ',
C    &GACMS,BGCMS,GALAB,BGLAB,BLAB,UMO,PCM,EPROJ,PPROJ
C     WRITE(6,*)' GAMP,BGAMP,BETP',GAMP,BGAMP,BETP

      RETURN
      END
*                                                                      *
*=== energy ===========================================================*
*                                                                      *
      DOUBLE PRECISION FUNCTION ENERGY (A,Z)
 
C     INCLUDE '(DBLPRC)'
*$ CREATE DBLPRC.ADD
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
      PARAMETER ( ONEFOU = ONEONE / FOUFOU )
      PARAMETER ( THRTWO = THRTHR / TWOTWO )
      PARAMETER ( PIPIPI = 3.141592653589793238462643383279D+00 )
      PARAMETER ( TWOPIP = 6.283185307179586476925286766559D+00 )
      PARAMETER ( PIP5O2 = 7.853981633974483096156608458199D+00 )
      PARAMETER ( PIPISQ = 9.869604401089358618834490999876D+00 )
      PARAMETER ( PIHALF = 1.570796326794896619231321691640D+00 )
      PARAMETER ( ERFA00 = 0.886226925452758013649083741671D+00 )
      PARAMETER ( ENEPER = 2.718281828459045235360287471353D+00 )
      PARAMETER ( SQRENT = 1.648721270700128146848650787814D+00 )
      PARAMETER ( SQRSIX = 2.449489742783178098197284074706D+00 )
      PARAMETER ( SQRSEV = 2.645751311064590590501615753639D+00 )
      PARAMETER ( SQRT12 = 3.464101615137754587054892683012D+00 )
      PARAMETER ( CLIGHT = 2.99792458         D+10 )
      PARAMETER ( AVOGAD = 6.0221367          D+23 )
      PARAMETER ( BOLTZM = 1.380658           D-23 )
      PARAMETER ( AMELGR = 9.1093897          D-28 )
      PARAMETER ( PLCKBR = 1.05457266         D-27 )
      PARAMETER ( ELCCGS = 4.8032068          D-10 )
      PARAMETER ( ELCMKS = 1.60217733         D-19 )
      PARAMETER ( AMUGRM = 1.6605402          D-24 )
      PARAMETER ( AMMUMU = 0.113428913        D+00 )
      PARAMETER ( AMPRMU = 1.007276470        D+00 )
      PARAMETER ( AMNEMU = 1.008664904        D+00 )
      PARAMETER ( ALPFSC = 7.2973530791728595 D-03 )
      PARAMETER ( FSCTO2 = 5.3251361962113614 D-05 )
      PARAMETER ( FSCTO3 = 3.8859399018437826 D-07 )
      PARAMETER ( FSCTO4 = 2.8357075508200407 D-09 )
      PARAMETER ( PLABRC = 0.197327053        D+00 )
      PARAMETER ( AMELCT = 0.51099906         D-03 )
      PARAMETER ( AMUGEV = 0.93149432         D+00 )
      PARAMETER ( AMMUON = 0.105658389        D+00 )
      PARAMETER ( AMPRTN = 0.93827231         D+00 )
      PARAMETER ( AMNTRN = 0.93956563         D+00 )
      PARAMETER ( AMDEUT = 1.87561339         D+00 )
      PARAMETER ( COUGFM = ELCCGS * ELCCGS / ELCMKS * 1.D-07 * 1.D+13
     &                   * 1.D-09 )
      PARAMETER ( RCLSEL = 2.8179409183694872 D-13 )
      PARAMETER ( BLTZMN = 8.617385           D-14 )
      PARAMETER ( GEVMEV = 1.0                D+03 )
      PARAMETER ( EMVGEV = 1.0                D-03 )
      PARAMETER ( ALGVMV = 6.90775527898214   D+00 )
      PARAMETER ( RADDEG = 180.D+00 / PIPIPI )
      PARAMETER ( DEGRAD = PIPIPI / 180.D+00 )
      LOGICAL LGBIAS, LGBANA
      COMMON / GLOBAL / LGBIAS, LGBANA
C     INCLUDE '(DIMPAR)'
*$ CREATE DIMPAR.ADD
      PARAMETER ( MXXRGN = 5000 )
      PARAMETER ( MXXMDF = 56   )
      PARAMETER ( MXXMDE = 50   )
      PARAMETER ( MFSTCK = 1000 )
      PARAMETER ( MESTCK = 100  )
      PARAMETER ( NALLWP = 39   )
      PARAMETER ( MPDPDX = 8    )
      PARAMETER ( ICOMAX = 180  )
      PARAMETER ( NSTBIS = 304  )
      PARAMETER ( IDMAXP = 210  )
      PARAMETER ( IDMXDC = 620  )
      PARAMETER ( MKBMX1 = 1    )
      PARAMETER ( MKBMX2 = 1    )
C     INCLUDE '(IOUNIT)'
*$ CREATE IOUNIT.ADD
      PARAMETER ( LUNIN  = 5  )
      PARAMETER ( LUNOUT = 6  )
      PARAMETER ( LUNERR = 15 )
      PARAMETER ( LUNBER = 14 )
      PARAMETER ( LUNECH = 8  )
      PARAMETER ( LUNFLU = 13 )
      PARAMETER ( LUNGEO = 16 )
      PARAMETER ( LUNPGS = 12 )
      PARAMETER ( LUNRAN = 2  )
      PARAMETER ( LUNXSC = 9  )
      PARAMETER ( LUNDET = 17 )
      PARAMETER ( LUNRAY = 10 )
      PARAMETER ( LUNRDB = 1  )
*
*----------------------------------------------------------------------*
*                                                                      *
*     Revised version of the original routine from EVAP:               *
*                                                                      *
*     Created on   15 may 1990     by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 01-oct-94     by    Alfredo Ferrari               *
*                                                                      *
*     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     *
*     !!!  It is supposed to be used with the updated atomic   !!!     *
*     !!!                    mass data file                    !!!     *
*     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     *
*                                                                      *
*----------------------------------------------------------------------*
*
*  Mass number below which "unknown" isotopes out of the Z-interval
*  reported in the mass tabulations are completely unstable and made
*  up by Z proton masses + N neutron masses:
      PARAMETER ( KAFREE =  4 )
*  Mass number below which "unknown" isotopes out of the Z-interval
*  reported in the mass tabulations are supposed to be particle unstable
      PARAMETER ( KAPUNS = 12 )
*  Minimum energy required for partilce unstable isotopes
      PARAMETER ( DEPUNS = 0.5D+00 )
*
C     INCLUDE '(EVA0)'
*$ CREATE EVA0.ADD
      COMMON / EVA0 / Y0, B0, P0 (1001), P1 (1001), P2 (1001),
     *                FLA (6), FLZ (6), RHO (6), OMEGA (6), EXMASS (6),
     *                CAM2 (130), CAM3 (200), CAM4 (130), CAM5 (200),
     *                T (4,7), RMASS (297), ALPH (297), BET (297),
     *                APRIME (250), IA (6), IZ (6)
C     INCLUDE '(ISOTOP)'
*$ CREATE ISOTOP.ADD
      PARAMETER ( NAMSMX = 270 )
      PARAMETER ( NZGVAX =  15 )
      PARAMETER ( NISMMX = 574 )
      COMMON / ISOTOP / WAPS   (NAMSMX,NZGVAX),  T12NUC (NAMSMX,NZGVAX),
     &                  WAPISM (NISMMX), T12ISM (NISMMX),
     &                  ABUISO (NSTBIS), ASTLIN (2,100), ZSTLIN (2,260),
     &                  AMSSST (100)  , ISOMNM (NSTBIS), ISONDX (2,100),
     &                  JSPNUC (NAMSMX,NZGVAX), JPTNUC (NAMSMX,NZGVAX),
     &                  INWAPS (NAMSMX), JSPISM (NISMMX),
     &                  JPTISM (NISMMX), IZWISM (NISMMX),
     &                  INWISM (0:NAMSMX)
*
      SAVE KA0, KZ0, IZ0
      DATA KA0, KZ0, IZ0 / -1, -1, -1 /
*
      KA0 = NINT ( A )
      KZ0 = NINT ( Z )
      N   = KA0 - KZ0
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C  |  No residual nucleus:
      IF ( KA0 .EQ. 0 .AND. KZ0 .LE. 0 ) THEN
        ENERGY = ZERZER
        RETURN
      END IF
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*  +-------------------------------------------------------------------*
*  |  Only protons:
      IF ( N .LE. 0 ) THEN
         IF ( KA0 .NE. 1 ) THEN
            IF ( N .LT. 0 ) THEN
               WRITE ( LUNOUT, * )
     &      ' FLUKA stopped in energy: mass number =< atomic number !!',
     &        KA0, KZ0
               WRITE ( LUNOUT, * )
     &      ' FLUKA stopped in energy: mass number =< atomic number !!',
     &        KA0, KZ0
               WRITE ( 77, * )
     &   ' ^^^FLUKA stopped in energy: mass number =< atomic number !!',
     &        KA0, KZ0
               STOP 'ENERGY:KA0-KZ0'
            END IF
         ELSE
            ENERGY = WAPS ( 1, 2 )
            IZ0    = -1
            RETURN
         END IF
      END IF
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C  |  Only neutrons:
      IF ( KZ0 .LE. 0 ) THEN
        IF ( KZ0 .LT. 0 ) THEN
          WRITE ( LUNOUT, * )
     &    ' DPMJET stopped in energy: -Z number =< atomic number!!',
     &    KA0, KZ0
          WRITE ( LUNOUT, * )
     &    ' DPMJET stopped in energy: -Z number =< atomic number!!',
     &    KA0, KZ0
          WRITE ( 77, * )
     &    ' DPMJET stopped in energy: -Z number =< atomic number!!',
     &    KA0, KZ0
          STOP 'ENERGY:KA0-KZ0'
        ELSE
          IZ0 = -1
          ENERGY = A * WAPS(1,1)
        END IF
        RETURN
      END IF
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |  A larger than maximum allowed:
      IF ( KA0 .GT. NAMSMX ) THEN
         ENERGY = ENRG ( A, Z )
         IZ0    = -1
         RETURN
      END IF
*  |
*  +-------------------------------------------------------------------*
      IZZ = INWAPS ( KA0 )
*  +-------------------------------------------------------------------*
*  |  Too much neutron rich with respect to the stability line:
      IF ( KZ0 .LT. IZZ ) THEN
*  |  +----------------------------------------------------------------*
*  |  |  Up to A=Kafree all "bound" masses are known, set it unbound:
         IF ( KA0 .LE. KAFREE ) THEN
            ENERGY = ( A - Z ) * WAPS (1,1) + Z * WAPS (1,2)
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |  Up to Kapuns: be sure it is particle unstable
         ELSE IF ( KA0 .LE. KAPUNS ) THEN
            ENERGY = ENRG ( A, Z )
            JZZ    = INWAPS ( KA0 - 1 )
            LZZ    = INWAPS ( KA0 - 2 )
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Residual mass for n-decay known:
            IF ( KZ0 .GE. JZZ .AND. KZ0 .LE. JZZ + NZGVAX - 1 ) THEN
               IZ0    = KZ0 - JZZ + 1
               ENERGY = MAX ( ENERGY, WAPS (KA0-1,IZ0) + WAPS (1,1)
     &                      + DEPUNS )
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Residual mass for 2n-decay known:
            ELSE IF ( KZ0 .GE. LZZ .AND. KZ0 .LE. LZZ + NZGVAX - 1 )THEN
               IZ0    = KZ0 - LZZ + 1
               ENERGY = MAX ( ENERGY, WAPS (KA0-2,IZ0) + TWOTWO *
     &                      ( WAPS (1,1) + DEPUNS ) )
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Set it unbound:
            ELSE
               ENERGY = AINFNT
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  Be sure not to have a positive energy state:
            ENERGY = MIN ( ENERGY, (A-Z) * WAPS (1,1) + Z * WAPS (1,2) )
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |  Proceed as usual:
         ELSE
            ENERGY = ENRG (A,Z)
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         IZ0    = -1
         RETURN
*  |
*  +-------------------------------------------------------------------*
*  |  Too much proton rich with respect to the stability line:
      ELSE IF ( KZ0 .GT. IZZ + NZGVAX - 1 ) THEN
*  |  +----------------------------------------------------------------*
*  |  |  Up to A=Kafree all "bound" masses are known, set it unbound:
         IF ( KA0 .LE. KAFREE ) THEN
            ENERGY = ( A - Z ) * WAPS (1,1) + Z * WAPS (1,2)
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |  Up to Kapuns: be sure it is particle unstable
         ELSE IF ( KA0 .LE. KAPUNS ) THEN
            ENERGY = ENRG ( A, Z )
            JZZ    = INWAPS ( KA0 - 1 )
            LZZ    = INWAPS ( KA0 - 2 )
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Residual mass for p-decay known:
            IF ( KZ0-1 .GE. JZZ .AND. KZ0-1 .LE. JZZ + NZGVAX - 1 ) THEN
               IZ0    = KZ0 - 1 - JZZ + 1
               ENERGY = MAX ( ENERGY, WAPS (KA0-1,IZ0) + WAPS (1,2)
     &                      + DEPUNS )
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Residual mass for 2p-decay known:
            ELSE IF ( KZ0-2 .GE. LZZ .AND. KZ0-2 .LE. LZZ + NZGVAX - 1 )
     &         THEN
               IZ0    = KZ0 - 2 - LZZ + 1
               ENERGY = MAX ( ENERGY, WAPS (KA0-2,IZ0) + TWOTWO *
     &                      ( WAPS (1,2) + DEPUNS ) )
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Set it unbound:
            ELSE
               ENERGY = AINFNT
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  Be sure not to have a positive energy state:
            ENERGY = MIN ( ENERGY, (A-Z) * WAPS (1,1) + Z * WAPS (1,2) )
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |  Proceed as usual:
         ELSE
            ENERGY = ENRG (A,Z)
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
        IZ0    = -1
         RETURN
*  |
*  +-------------------------------------------------------------------*
*  |  Known isotope or anyway isotope "inside" the stability zone
      ELSE
         IZ0    = KZ0 - IZZ + 1
         ENERGY = WAPS ( KA0, IZ0 )
*  |  +----------------------------------------------------------------*
*  |  |  Mass not known
         IF ( ABS (ENERGY) .LT. ANGLGB .AND. (KA0 .NE. 12 .OR. KZ0
     &        .NE. 6) ) THEN
            IZ0    = -1
            ENERGY = ENRG ( A, Z )
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         RETURN
      END IF
*  |
*  +-------------------------------------------------------------------*
*=== End of Function Energy ===========================================*
*     RETURN
      END
*$ CREATE ENRG.FOR
*COPY ENRG
*                                                                      *
*=== enrg =============================================================*
*                                                                      *
      DOUBLE PRECISION FUNCTION ENRG(A,Z)
 
C     INCLUDE '(DBLPRC)'
*$ CREATE DBLPRC.ADD
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
      PARAMETER ( ONEFOU = ONEONE / FOUFOU )
      PARAMETER ( THRTWO = THRTHR / TWOTWO )
      PARAMETER ( PIPIPI = 3.141592653589793238462643383279D+00 )
      PARAMETER ( TWOPIP = 6.283185307179586476925286766559D+00 )
      PARAMETER ( PIP5O2 = 7.853981633974483096156608458199D+00 )
      PARAMETER ( PIPISQ = 9.869604401089358618834490999876D+00 )
      PARAMETER ( PIHALF = 1.570796326794896619231321691640D+00 )
      PARAMETER ( ERFA00 = 0.886226925452758013649083741671D+00 )
      PARAMETER ( ENEPER = 2.718281828459045235360287471353D+00 )
      PARAMETER ( SQRENT = 1.648721270700128146848650787814D+00 )
      PARAMETER ( SQRSIX = 2.449489742783178098197284074706D+00 )
      PARAMETER ( SQRSEV = 2.645751311064590590501615753639D+00 )
      PARAMETER ( SQRT12 = 3.464101615137754587054892683012D+00 )
      PARAMETER ( CLIGHT = 2.99792458         D+10 )
      PARAMETER ( AVOGAD = 6.0221367          D+23 )
      PARAMETER ( BOLTZM = 1.380658           D-23 )
      PARAMETER ( AMELGR = 9.1093897          D-28 )
      PARAMETER ( PLCKBR = 1.05457266         D-27 )
      PARAMETER ( ELCCGS = 4.8032068          D-10 )
      PARAMETER ( ELCMKS = 1.60217733         D-19 )
      PARAMETER ( AMUGRM = 1.6605402          D-24 )
      PARAMETER ( AMMUMU = 0.113428913        D+00 )
      PARAMETER ( AMPRMU = 1.007276470        D+00 )
      PARAMETER ( AMNEMU = 1.008664904        D+00 )
      PARAMETER ( ALPFSC = 7.2973530791728595 D-03 )
      PARAMETER ( FSCTO2 = 5.3251361962113614 D-05 )
      PARAMETER ( FSCTO3 = 3.8859399018437826 D-07 )
      PARAMETER ( FSCTO4 = 2.8357075508200407 D-09 )
      PARAMETER ( PLABRC = 0.197327053        D+00 )
      PARAMETER ( AMELCT = 0.51099906         D-03 )
      PARAMETER ( AMUGEV = 0.93149432         D+00 )
      PARAMETER ( AMMUON = 0.105658389        D+00 )
      PARAMETER ( AMPRTN = 0.93827231         D+00 )
      PARAMETER ( AMNTRN = 0.93956563         D+00 )
      PARAMETER ( AMDEUT = 1.87561339         D+00 )
      PARAMETER ( COUGFM = ELCCGS * ELCCGS / ELCMKS * 1.D-07 * 1.D+13
     &                   * 1.D-09 )
      PARAMETER ( RCLSEL = 2.8179409183694872 D-13 )
      PARAMETER ( BLTZMN = 8.617385           D-14 )
      PARAMETER ( GEVMEV = 1.0                D+03 )
      PARAMETER ( EMVGEV = 1.0                D-03 )
      PARAMETER ( ALGVMV = 6.90775527898214   D+00 )
      PARAMETER ( RADDEG = 180.D+00 / PIPIPI )
      PARAMETER ( DEGRAD = PIPIPI / 180.D+00 )
      LOGICAL LGBIAS, LGBANA
      COMMON / GLOBAL / LGBIAS, LGBANA
C     INCLUDE '(DIMPAR)'
*$ CREATE DIMPAR.ADD
      PARAMETER ( MXXRGN = 5000 )
      PARAMETER ( MXXMDF = 56   )
      PARAMETER ( MXXMDE = 50   )
      PARAMETER ( MFSTCK = 1000 )
      PARAMETER ( MESTCK = 100  )
      PARAMETER ( NALLWP = 39   )
      PARAMETER ( MPDPDX = 8    )
      PARAMETER ( ICOMAX = 180  )
      PARAMETER ( NSTBIS = 304  )
      PARAMETER ( IDMAXP = 210  )
      PARAMETER ( IDMXDC = 620  )
      PARAMETER ( MKBMX1 = 1    )
      PARAMETER ( MKBMX2 = 1    )
C     INCLUDE '(IOUNIT)'
*$ CREATE IOUNIT.ADD
      PARAMETER ( LUNIN  = 5  )
      PARAMETER ( LUNOUT = 6  )
      PARAMETER ( LUNERR = 15 )
      PARAMETER ( LUNBER = 14 )
      PARAMETER ( LUNECH = 8  )
      PARAMETER ( LUNFLU = 13 )
      PARAMETER ( LUNGEO = 16 )
      PARAMETER ( LUNPGS = 12 )
      PARAMETER ( LUNRAN = 2  )
      PARAMETER ( LUNXSC = 9  )
      PARAMETER ( LUNDET = 17 )
      PARAMETER ( LUNRAY = 10 )
      PARAMETER ( LUNRDB = 1  )
*
*----------------------------------------------------------------------*
*                                                                      *
*     Revised version of the original routine from EVAP:               *
*                                                                      *
*     Created on   15 may 1990     by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 01-oct-94     by    Alfredo Ferrari               *
*                                                                      *
*     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     *
*     !!!  It is supposed to be used with the updated atomic   !!!     *
*     !!!                    mass data file                    !!!     *
*     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     *
*                                                                      *
*----------------------------------------------------------------------*
*
      PARAMETER ( O16OLD = 931.145  D+00 )
      PARAMETER ( O16NEW = 931.19826D+00 )
      PARAMETER ( O16RAT = O16NEW / O16OLD )
      PARAMETER ( C12NEW = 931.49432D+00 )
      PARAMETER ( ADJUST = -8.322737768178909D-02 )
C     INCLUDE '(EVA0)'
*$ CREATE EVA0.ADD
      COMMON / EVA0 / Y0, B0, P0 (1001), P1 (1001), P2 (1001),
     *                FLA (6), FLZ (6), RHO (6), OMEGA (6), EXMASS (6),
     *                CAM2 (130), CAM3 (200), CAM4 (130), CAM5 (200),
     *                T (4,7), RMASS (297), ALPH (297), BET (297),
     *                APRIME (250), IA (6), IZ (6)
      LOGICAL LFIRST
      SAVE LFIRST, EXHYDR, EXNEUT
      DATA LFIRST / .TRUE. /
      DATA NERG1/ 0/
*
      IF ( LFIRST ) THEN
         LFIRST = .FALSE.
         EXHYDR = ENERGY ( ONEONE, ONEONE )
         EXNEUT = ENERGY ( ONEONE, ZERZER )
      END IF
      IZ0 = NINT (Z)
      IF ( IZ0 .LE. 0 ) THEN
         ENRG = A * EXNEUT
         RETURN
      END IF
      IF (A .EQ. 0.D0)THEN
	WRITE (6,'(A)')' ENRG A=0.'
	ENRG = 0
	RETURN
      ENDIF
      N   = NINT (A-Z)
      IF ( N .LE. 0 ) THEN
         ENRG = Z * EXHYDR
         RETURN
      END IF
      AM2ZOA= (A-Z-Z)/A
      AM2ZOA=AM2ZOA*AM2ZOA
      A13 = RMASS(NINT(A))
*     A13 = A**.3333333333333333D+00
      IF(A13 .EQ. 0.D0) THEN
	NERG1=NERG1+1
	IF(NERG1.LE.50)WRITE (6,'(A)')' ENRG A13=0.'
	ENRG = 0
	RETURN
      ENDIF
      AM13 = 1.D+00/A13
      EV=-17.0354D+00*(1.D+00 -1.84619 D+00*AM2ZOA)*A
      ES= 25.8357D+00*(1.D+00 -1.712185D+00*AM2ZOA)*
     &    (1.D+00 -0.62025D+00*AM13*AM13)*
     &    (A13*A13 -.62025D+00)
      EC= 0.799D+00*Z*(Z-1.D+00)*AM13*(((1.5772D+00*AM13 +1.2273D+00)*
     &    AM13-1.5849D+00)*
     &    AM13*AM13 +1.D+00)
      EEX= -0.4323D+00*AM13*Z**1.3333333D+00*
     &   (((0.49597D+00*AM13 -0.14518D+00)*AM13 -0.57811D+00) * AM13
     &   + 1.D+00)
      ENRG  =8.367D+00*A -0.783D+00*Z +EV +ES +EC +EEX+CAM2(IZ0)+CAM3(N)
      ENRG  = ( ENRG + A * O16OLD ) * O16RAT - A * ( C12NEW - ADJUST )
      ENRG  = MIN ( ENRG, Z * EXHYDR + ( A - Z ) * EXNEUT )
      RETURN
*=== End of function Enrg =============================================*
      END
*$ CREATE BERTTP.FOR
*COPY BERTTP
*                                                                      *
*=== berttp ===========================================================*
*                                                                      *
      SUBROUTINE BERTTP
C     INCLUDE '(DBLPRC)'
*$ CREATE DBLPRC.ADD
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
      PARAMETER ( ONEFOU = ONEONE / FOUFOU )
      PARAMETER ( THRTWO = THRTHR / TWOTWO )
      PARAMETER ( PIPIPI = 3.141592653589793238462643383279D+00 )
      PARAMETER ( TWOPIP = 6.283185307179586476925286766559D+00 )
      PARAMETER ( PIP5O2 = 7.853981633974483096156608458199D+00 )
      PARAMETER ( PIPISQ = 9.869604401089358618834490999876D+00 )
      PARAMETER ( PIHALF = 1.570796326794896619231321691640D+00 )
      PARAMETER ( ERFA00 = 0.886226925452758013649083741671D+00 )
      PARAMETER ( ENEPER = 2.718281828459045235360287471353D+00 )
      PARAMETER ( SQRENT = 1.648721270700128146848650787814D+00 )
      PARAMETER ( SQRSIX = 2.449489742783178098197284074706D+00 )
      PARAMETER ( SQRSEV = 2.645751311064590590501615753639D+00 )
      PARAMETER ( SQRT12 = 3.464101615137754587054892683012D+00 )
      PARAMETER ( CLIGHT = 2.99792458         D+10 )
      PARAMETER ( AVOGAD = 6.0221367          D+23 )
      PARAMETER ( BOLTZM = 1.380658           D-23 )
      PARAMETER ( AMELGR = 9.1093897          D-28 )
      PARAMETER ( PLCKBR = 1.05457266         D-27 )
      PARAMETER ( ELCCGS = 4.8032068          D-10 )
      PARAMETER ( ELCMKS = 1.60217733         D-19 )
      PARAMETER ( AMUGRM = 1.6605402          D-24 )
      PARAMETER ( AMMUMU = 0.113428913        D+00 )
      PARAMETER ( AMPRMU = 1.007276470        D+00 )
      PARAMETER ( AMNEMU = 1.008664904        D+00 )
      PARAMETER ( ALPFSC = 7.2973530791728595 D-03 )
      PARAMETER ( FSCTO2 = 5.3251361962113614 D-05 )
      PARAMETER ( FSCTO3 = 3.8859399018437826 D-07 )
      PARAMETER ( FSCTO4 = 2.8357075508200407 D-09 )
      PARAMETER ( PLABRC = 0.197327053        D+00 )
      PARAMETER ( AMELCT = 0.51099906         D-03 )
      PARAMETER ( AMUGEV = 0.93149432         D+00 )
      PARAMETER ( AMMUON = 0.105658389        D+00 )
      PARAMETER ( AMPRTN = 0.93827231         D+00 )
      PARAMETER ( AMNTRN = 0.93956563         D+00 )
      PARAMETER ( AMDEUT = 1.87561339         D+00 )
      PARAMETER ( COUGFM = ELCCGS * ELCCGS / ELCMKS * 1.D-07 * 1.D+13
     &                   * 1.D-09 )
      PARAMETER ( RCLSEL = 2.8179409183694872 D-13 )
      PARAMETER ( BLTZMN = 8.617385           D-14 )
      PARAMETER ( GEVMEV = 1.0                D+03 )
      PARAMETER ( EMVGEV = 1.0                D-03 )
      PARAMETER ( ALGVMV = 6.90775527898214   D+00 )
      PARAMETER ( RADDEG = 180.D+00 / PIPIPI )
      PARAMETER ( DEGRAD = PIPIPI / 180.D+00 )
      LOGICAL LGBIAS, LGBANA
      COMMON / GLOBAL / LGBIAS, LGBANA
C     INCLUDE '(DIMPAR)'
*$ CREATE DIMPAR.ADD
      PARAMETER ( MXXRGN = 5000 )
      PARAMETER ( MXXMDF = 56   )
      PARAMETER ( MXXMDE = 50   )
      PARAMETER ( MFSTCK = 1000 )
      PARAMETER ( MESTCK = 100  )
      PARAMETER ( NALLWP = 39   )
      PARAMETER ( MPDPDX = 8    )
      PARAMETER ( ICOMAX = 180  )
      PARAMETER ( NSTBIS = 304  )
      PARAMETER ( IDMAXP = 210  )
      PARAMETER ( IDMXDC = 620  )
      PARAMETER ( MKBMX1 = 1    )
      PARAMETER ( MKBMX2 = 1    )
C     INCLUDE '(IOUNIT)'
*$ CREATE IOUNIT.ADD
      PARAMETER ( LUNIN  = 5  )
      PARAMETER ( LUNOUT = 6  )
      PARAMETER ( LUNERR = 15 )
      PARAMETER ( LUNBER = 14 )
      PARAMETER ( LUNECH = 8  )
      PARAMETER ( LUNFLU = 13 )
      PARAMETER ( LUNGEO = 16 )
      PARAMETER ( LUNPGS = 12 )
      PARAMETER ( LUNRAN = 2  )
      PARAMETER ( LUNXSC = 9  )
      PARAMETER ( LUNDET = 17 )
      PARAMETER ( LUNRAY = 10 )
      PARAMETER ( LUNRDB = 1  )
C---------------------------------------------------------------------
C SUBNAME = BERTTP --- READ BERTINI DATA
C---------------------------------------------------------------------
C     ---------------------------------- I-N-C DATA
C     COMMON R8(2127),R4(64),CRSC(600,4),R8B(336),CS(29849)
C     REAL*8 R8,R8B,CRSC,CS
C     REAL*4 R4
C     --------------------------------- EVAPORATION DATA
C     INCLUDE '(COOKCM)'
*$ CREATE COOKCM.ADD
      PARAMETER ( ASMTOG = SIXSIX / PIPIPI**2 )
      LOGICAL LDEFOZ, LDEFON
      PARAMETER ( INCOOK = 150, IZCOOK = 98 )
      COMMON / COOKCM / ALPIGN, BETIGN, GAMIGN, POWIGN,
     &                SZCOOK (IZCOOK), SNCOOK (INCOOK), PZCOOK (IZCOOK),
     &                PNCOOK (INCOOK), LDEFOZ (IZCOOK), LDEFON (INCOOK)
C     INCLUDE '(EVA0)'
*$ CREATE EVA0.ADD
      COMMON / EVA0 / Y0, B0, P0 (1001), P1 (1001), P2 (1001),
     *                FLA (6), FLZ (6), RHO (6), OMEGA (6), EXMASS (6),
     *                CAM2 (130), CAM3 (200), CAM4 (130), CAM5 (200),
     *                T (4,7), RMASS (297), ALPH (297), BET (297),
     *                APRIME (250), IA (6), IZ (6)
C     INCLUDE '(FRBKCM)'
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
C     INCLUDE '(HETTP)'
*$ CREATE HETTP.ADD
      COMMON /HETTP/  NHSTP,NBERTP,IOSUB,INSRS
C     INCLUDE '(INPFLG)'
*$ CREATE INPFLG.ADD
      COMMON /INPFLG/ IANG,IFISS,IB0,IGEOM,ISTRAG,KEYDK
C     INCLUDE '(ISOTOP)'
*$ CREATE ISOTOP.ADD
      PARAMETER ( NAMSMX = 270 )
      PARAMETER ( NZGVAX =  15 )
      PARAMETER ( NISMMX = 574 )
      COMMON / ISOTOP / WAPS   (NAMSMX,NZGVAX),  T12NUC (NAMSMX,NZGVAX),
     &                  WAPISM (NISMMX), T12ISM (NISMMX),
     &                  ABUISO (NSTBIS), ASTLIN (2,100), ZSTLIN (2,260),
     &                  AMSSST (100)  , ISOMNM (NSTBIS), ISONDX (2,100),
     &                  JSPNUC (NAMSMX,NZGVAX), JPTNUC (NAMSMX,NZGVAX),
     &                  INWAPS (NAMSMX), JSPISM (NISMMX),
     &                  JPTISM (NISMMX), IZWISM (NISMMX),
     &                  INWISM (0:NAMSMX)
C     INCLUDE '(NUCGEO)'
*$ CREATE NUCGEO.ADD
      PARAMETER ( PI     = PIPIPI )
      PARAMETER ( PISQ   = PIPISQ )
      PARAMETER ( SKTOHL = 0.5456645846610345D+00 )
      PARAMETER ( RZNUCL = 1.12        D+00 )
      PARAMETER ( RMSPRO = 0.8         D+00 )
      PARAMETER ( R0PROT = RMSPRO / SQRT12  )
      PARAMETER ( ARHPRO = 1.D+00 / 8.D+00 / PI / R0PROT / R0PROT
     &          / R0PROT )
      PARAMETER ( RLLE04 = RZNUCL )
      PARAMETER ( RLLE16 = RZNUCL )
      PARAMETER ( RLGT16 = RZNUCL )
      PARAMETER ( RCLE04 = 0.75D+00 / PI / RLLE04 / RLLE04 / RLLE04 )
      PARAMETER ( RCLE16 = 0.75D+00 / PI / RLLE16 / RLLE16 / RLLE16 )
      PARAMETER ( RCGT16 = 0.75D+00 / PI / RLGT16 / RLGT16 / RLGT16 )
      PARAMETER ( SKLE04 = 1.4D+00 )
      PARAMETER ( SKLE16 = 1.9D+00 )
      PARAMETER ( SKGT16 = 2.4D+00 )
      PARAMETER ( HLLE04 = SKTOHL * SKLE04 )
      PARAMETER ( HLLE16 = SKTOHL * SKLE16 )
      PARAMETER ( HLGT16 = SKTOHL * SKGT16 )
      PARAMETER ( ALPHA0 = 0.1D+00 )
      PARAMETER ( OMALH0 = 1.D+00 - ALPHA0 )
      PARAMETER ( GAMSK0 = 0.9D+00 )
      PARAMETER ( OMGAS0 = 1.D+00 - GAMSK0 )
      PARAMETER ( POTME0 = 0.6666666666666667D+00 )
      PARAMETER ( POTBA0 = 1.D+00 )
      PARAMETER ( PNFRAT = 1.533D+00 )
      PARAMETER ( RADPIM = 0.035D+00 )
      PARAMETER ( RDPMHL = 14.D+00   )
      PARAMETER ( APMRST = 4.D+00 / 44.D+00 )
      PARAMETER ( APMPRO = 1.D+00 / 6.D+00 )
      PARAMETER ( APPPRO = 5.D+00 / 6.D+00 )
      PARAMETER ( AP0PFS = 0.5D+00 )
      PARAMETER ( AP0PFP = 1.D+00 / 3.D+00 )
      PARAMETER ( AP0NFP = 2.D+00 / 3.D+00 )
      PARAMETER ( XPAUCO = 1.88495407241652 D+00 )
      PARAMETER ( MXSCIN = 50     )
      LOGICAL LABRST, LELSTC, LINELS, LCHEXC, LABSRP, LABSTH, LPREEQ,
     &        LNPHTC, LNWRAD, LPNRHO
      COMMON / NUCGID / RHOTAB (2:260), RHATAB (2:260), ALPTAB (2:260),
     &                  RADTAB (2:260), SKITAB (2:260), HALTAB (2:260),
     &                  SK3TAB (2:260), SK4TAB (2:260), HABTAB (2:260),
     &                  CWSTAB (2:260), EKATAB (2:260), PFATAB (2:260),
     &                  PFRTAB (2:260)
      COMMON / NUCGEO / RADTOT, RADIU1, RADIU0, RAD1O2, SKINDP, HALODP,
     &                  ALPHAL, OMALHL, RADSKN, SKNEFF, CPARWS, RADPRO,
     &                  RADCOR, RADCO2, RADMAX, BIMPTR, RIMPTR, XIMPTR,
     &                  YIMPTR, ZIMPTR, RHOIMT, EKFPRO, PFRPRO, RHOCEN,
     &                  RHOCOR, RHOSKN, EKFCEN (2), PFRCEN (2), EKFBIM,
     &                  PFRBIM, RHOIMP, EKFIMP, PFRIMP, RHOIM2, EKFIM2,
     &                  PFRIM2, RHOIM3, EKFIM3, PFRIM3, VPRWLL, RIMPCT,
     &                  BIMPCT, XIMPCT, YIMPCT, ZIMPCT, RIMPC2, XIMPC2,
     &                  YIMPC2, ZIMPC2, RIMPC3, XIMPC3, YIMPC3, ZIMPC3,
     &                  XBIMPC, YBIMPC, ZBIMPC, CXIMPC, CYIMPC, CZIMPC,
     &                  SQRIMP, SIGMAP, SIGMAN, SIGMAA, RHORED, R0TRAJ,
     &                  R1TRAJ, SBUSED, SBTOT , SBRES , RHOAVE, EKFAVE,
     &                  PFRAVE, AVEBIN, ACOLL , ZCOLL , RADSIG, OPACTY,
     &                  EKECON, PNUCCO, EKEWLL, PPRWLL, PXPROJ, PYPROJ,
     &                  PZPROJ, EKFERM, PNFRMI, PXFERM, PYFERM, PZFERM,
     &                  EKFER2, PNFRM2, PXFER2, PYFER2, PZFER2, EKFER3,
     &                  PNFRM3, PXFER3, PYFER3, PZFER3, RHOMEM, EKFMEM,
     &                  BIMMEM, WLLRED, VPRBIM, POTINC, POTOUT, EEXMIN
      COMMON / NUCGE2 / RDTTNC (2), RHONCP (2), RHONC2 (2), RHONC3 (2),
     &                  RHONCT (2), AMOTHR, EKOTHR, AMCREA, EKNCLN,
     &                  EEXDEL, EEXANY, CLMBBR, RDCLMB, BFCLMB, BFCEFF,
     &                  BNPROJ, BNDNUC, DEBRLM, SK4PAR, UBIMPC, VBIMPC,
     &                  WBIMPC, BNDPOT, SIGMAT, SIGABP, SIGABN, WLLRES,
     &                  POTBAR, POTMES, AGEPRI, OPNOPA,
     &                  BNENRG (3), DEFNUC (2), SIGMPR (4), SIGMNU (4),
     &                  SIGPAB (3), SIGNAB (3), HHLP   (2), FORTOT (2),
     &                  IPWELL, ITNCMX, KPRIN , NTARGT, KNUCIM, KNUCI2,
     &                  KNUCI3, IEVPRE, ISFCOL, ISFTAR, ISFTA2, ISFTA3,
     &                  NPOTHR, ICOTHR, IBOTHR, NPUMFN, ISTNCL, ITAUCM,
     &                  IADFLG, IGSFLG, IALFLG, ICBFLG, LPREEQ, LNPHTC,
     &                  LPNRHO, LNWRAD
      COMMON / NUCPWI / ALMBAR, BIMMAX, SIGGEO, LLLMAX, LLLACT
      COMMON / NUCGII / HOLEXP (2*MXSCIN), XEXPIN (3,MXSCIN),
     &                  YEXPIN (3,MXSCIN), ZEXPIN (3,MXSCIN),
     &                  AGEXIN (MXSCIN), RHOEXP (2), EKFEXP, EHLFIX,
     &                  NHLEXP, NHLFIX, IPRTYP, NNCEXI (MXSCIN),
     &                  NCEXPI (3,MXSCIN), ISEXIN (3,MXSCIN),
     &                  ISCTYP (MXSCIN), NUSCIN, NEXPEM,
     &                  LABRST, LELSTC, LINELS, LCHEXC, LABSRP, LABSTH
      DIMENSION AWSTAB (2:260), SIGMAB (3)
      EQUIVALENCE ( DEFPRO, DEFNUC (1) )
      EQUIVALENCE ( DEFNEU, DEFNUC (2) )
      EQUIVALENCE ( RHOIPP, RHONCP (1) )
      EQUIVALENCE ( RHOINP, RHONCP (2) )
      EQUIVALENCE ( RHOIP2, RHONC2 (1) )
      EQUIVALENCE ( RHOIN2, RHONC2 (2) )
      EQUIVALENCE ( RHOIP3, RHONC3 (1) )
      EQUIVALENCE ( RHOIN3, RHONC3 (2) )
      EQUIVALENCE ( RHOIPT, RHONCT (1) )
      EQUIVALENCE ( RHOINT, RHONCT (2) )
      EQUIVALENCE ( OMALHL, SK3PAR )
      EQUIVALENCE ( ALPHAL, HABPAR )
      EQUIVALENCE ( ALPTAB (2), AWSTAB (2) )
      EQUIVALENCE ( SIGMPE, SIGMPR (1) )
      EQUIVALENCE ( SIGMPC, SIGMPR (2) )
      EQUIVALENCE ( SIGMPI, SIGMPR (3) )
      EQUIVALENCE ( SIGMPA, SIGMPR (4) )
      EQUIVALENCE ( SIGMNE, SIGMNU (1) )
      EQUIVALENCE ( SIGMNC, SIGMNU (2) )
      EQUIVALENCE ( SIGMNI, SIGMNU (3) )
      EQUIVALENCE ( SIGMNA, SIGMNU (4) )
      EQUIVALENCE ( SIGMA2, SIGPAB (1) )
      EQUIVALENCE ( SIGMA3, SIGPAB (2) )
      EQUIVALENCE ( SIGMAS, SIGPAB (3) )
      EQUIVALENCE ( SIGPAB (1), SIGMAB (1) )
C     INCLUDE '(NUCLEV)'
*$ CREATE NUCLEV.ADD
      LOGICAL LCLVSL
      COMMON / NUCLEV / PAENUC (200,2), SHENUC (200,2), DEFRMI (2),
     &                  DEFMAG (2), ENNCLV (160,2), RANCLV (160,2),
     &                  CUMRAD (0:160,2), RUSNUC (2),
     &                  ENPLVL (114), ENNLVL(164), JUSNUC (160,2),
     &                  NTANUC (2), NAVNUC (2), NLSNUC (2), NCONUC (2),
     &                  NSKNUC (2), NHANUC (2), NUSNUC (2), JMXNUC (2),
     &                  IPRNUC (3), JPRNUC (3), MAGNUM (8), MAGNUC (2),
     &                  MGSNUC (8,2), MGSSNC (25,2), NSBSHL (2),
     &                  NPRNUC, INUCLV, LCLVSL
      DIMENSION JUSPRO (160), JUSNEU (160), MGSPRO (8), MGSNEU (8),
     &          MGSSPR (19) , MGSSNE (25)
      EQUIVALENCE ( RUSNUC (1), RUSPRO )
      EQUIVALENCE ( RUSNUC (2), RUSNEU )
      EQUIVALENCE ( JUSNUC (1,1), JUSPRO (1) )
      EQUIVALENCE ( JUSNUC (1,2), JUSNEU (1) )
      EQUIVALENCE ( MGSNUC (1,1), MGSPRO (1) )
      EQUIVALENCE ( MGSNUC (1,2), MGSNEU (1) )
      EQUIVALENCE ( MGSSNC (1,1), MGSSPR (1) )
      EQUIVALENCE ( MGSSNC (1,2), MGSSNE (1) )
      EQUIVALENCE ( NTANUC (1), NTAPRO )
      EQUIVALENCE ( NTANUC (2), NTANEU )
      EQUIVALENCE ( NAVNUC (1), NAVPRO )
      EQUIVALENCE ( NAVNUC (2), NAVNEU )
      EQUIVALENCE ( NLSNUC (1), NLSPRO )
      EQUIVALENCE ( NLSNUC (2), NLSNEU )
      EQUIVALENCE ( NCONUC (1), NCOPRO )
      EQUIVALENCE ( NCONUC (2), NCONEU )
      EQUIVALENCE ( NSKNUC (1), NSKPRO )
      EQUIVALENCE ( NSKNUC (2), NSKNEU )
      EQUIVALENCE ( NHANUC (1), NHAPRO )
      EQUIVALENCE ( NHANUC (2), NHANEU )
      EQUIVALENCE ( NUSNUC (1), NUSPRO )
      EQUIVALENCE ( NUSNUC (2), NUSNEU )
      EQUIVALENCE ( JMXNUC (1), JMXPRO )
      EQUIVALENCE ( JMXNUC (2), JMXNEU )
      EQUIVALENCE ( MAGNUC (1), MAGPRO )
      EQUIVALENCE ( MAGNUC (2), MAGNEU )
C     INCLUDE '(PAREVT)'
*$ CREATE PAREVT.ADD
      PARAMETER ( FRDIFF = 0.2D+00 )
      PARAMETER ( ETHSEA = 1.0D+00 )
 
      LOGICAL LDIFFR, LINCTV, LEVPRT, LHEAVY, LDEEXG, LGDHPR, LPREEX,
     &        LHLFIX, LPRFIX, LPARWV, LPOWER, LSNGCH, LLVMOD, LSCHDF
      COMMON / PAREVT / DPOWER, FSPRD0, FSHPFN, RN1GSC, RN2GSC,
     &                  LDIFFR (NALLWP),LPOWER, LINCTV, LEVPRT, LHEAVY,
     &                  LDEEXG, LGDHPR, LPREEX, LHLFIX, LPRFIX, LPARWV,
     &                  ILVMOD, JLVMOD, LLVMOD, LSNGCH, LSCHDF
C     INCLUDE '(XSEPAR)'
*$ CREATE XSEPAR.ADD
      COMMON / XSEPAR / AANXSE (100), BBNXSE (100), CCNXSE (100),
     &                  DDNXSE (100), EENXSE (100), ZZNXSE (100),
     &                  EMNXSE (100), XMNXSE (100),
     &                  AAPXSE (100), BBPXSE (100), CCPXSE (100),
     &                  DDPXSE (100), EEPXSE (100), FFPXSE (100),
     &                  ZZPXSE (100), EMPXSE (100), XMPXSE (100)
 
C---------------------------------------------------------------------
      NBERTP=LUNBER
      WRITE( LUNOUT,'(A,I2)')
     & ' *** Reading evaporation and nuclear data from unit: ', NBERTP
      REWIND NBERTP
C A. Ferrari: first of all read isotopic data
      READ (NBERTP) ISONDX
      READ (NBERTP) ISOMNM
      READ (NBERTP) ABUISO
      DO 1 I=1,4
C        READ  (NBERTP) (CRSC(J,I),J=1,600)
C A. Ferrari: commented also the dummy read to save disk space
C        READ  (NBERTP)
    1 CONTINUE
C     READ  (NBERTP) CS
C A. Ferrari: commented also the dummy read to save disk space
C     READ  (NBERTP)
C---------------------------------------------------------------------
      READ (NBERTP) (P0(I),P1(I),P2(I),I=1,1001)
      READ (NBERTP) IA,IZ
      DO 2 I=1,6
         FLA(I)=IA(I)
         FLZ(I)=IZ(I)
    2 CONTINUE
      READ (NBERTP) RHO,OMEGA
      READ (NBERTP) EXMASS
      READ (NBERTP) CAM2
      READ (NBERTP) CAM3
      READ (NBERTP) CAM4
      READ (NBERTP) CAM5
      READ (NBERTP) ((T(I,J),J=1,7),I=1,3)
      DO 3 I=1,7
         T(4,I) = ZERZER
    3 CONTINUE
      READ (NBERTP) RMASS
      READ (NBERTP) ALPH
      READ (NBERTP) BET
      READ (NBERTP) INWAPS
      READ (NBERTP) WAPS
      READ (NBERTP) T12NUC
      READ (NBERTP) JSPNUC
      READ (NBERTP) JPTNUC
      READ (NBERTP) INWISM
      READ (NBERTP) IZWISM
      READ (NBERTP) WAPISM
      READ (NBERTP) T12ISM
      READ (NBERTP) JSPISM
      READ (NBERTP) JPTISM
      READ (NBERTP) APRIME
      WRITE( LUNOUT,'(A)' ) ' *** Evaporation: using 1977 Waps data ***'
      READ (NBERTP) AHELP , BHELP , LRMSCH, LRD1O2, LTRASP
      IF ( ABS (AHELP-ALPHA0) .GT. CSNNRM * ALPHA0 .OR.
     &     ABS (BHELP-GAMSK0) .GT. CSNNRM * GAMSK0 ) THEN
         WRITE (LUNOUT,*)
     &         ' *** Inconsistent Nuclear Geometry data on file ***'
         STOP
      END IF
      READ (NBERTP) RHOTAB, RHATAB, ALPTAB, RADTAB, SKITAB, HALTAB,
     &              EKATAB, PFATAB, PFRTAB
      READ (NBERTP) AANXSE, BBNXSE, CCNXSE, DDNXSE, EENXSE, ZZNXSE,
     &              EMNXSE, XMNXSE
      READ (NBERTP) AAPXSE, BBPXSE, CCPXSE, DDPXSE, EEPXSE, FFPXSE,
     &              ZZPXSE, EMPXSE, XMPXSE
*  Data about Fermi-breakup:
      READ (NBERTP) IPOSST, MXPDUM, MXADUM, MXNDUM, MXZDUM, IFBSTF
      IF ( MXADUM .NE. MXAFBK .OR. MXNDUM .NE. MXNFBK .OR. MXZDUM .NE.
     &     MXZFBK .OR. MXPDUM .NE. MXPSST ) THEN
         WRITE (LUNOUT,*)' *** Inconsistent Fermi BreakUp data',
     &                   ' in the Nuclear Data file ***'
         STOP 'STOP:BERTTP-INCONS-FERMI-BREAKUP-DATA'
      END IF
      READ (NBERTP) IFRBKN
      READ (NBERTP) IFRBKZ
      READ (NBERTP) IFBKSP
      READ (NBERTP) IFBKST
      READ (NBERTP) EEXFBK
      CLOSE (UNIT=NBERTP)
      DO 100 JZ = 1, 130
         SHENUC ( JZ, 1 ) = EMVGEV * ( CAM2 (JZ) + CAM4 (JZ) )
  100 CONTINUE
      DO 200 JA = 1, 200
         SHENUC ( JA, 2 ) = EMVGEV * ( CAM3 (JA) + CAM5 (JA) )
  200 CONTINUE
      CALL STALIN
      IF ( ILVMOD .LE. 0 ) THEN
         ILVMOD = IB0
      ELSE
         IB0 = ILVMOD
      END IF
      IF ( LLVMOD ) THEN
         DO 300 JZ = 1, IZCOOK
            CAM4 (JZ) = PZCOOK (JZ)
  300    CONTINUE
         DO 400 JN = 1, INCOOK
            CAM5 (JN) = PNCOOK (JZ)
  400    CONTINUE
      END IF
      WRITE (LUNOUT,*)
      IF ( ILVMOD .EQ. 1 ) THEN
         WRITE (LUNOUT,*)
     &' **** Standard EVAP T=0 level density used ****'
      ELSE IF ( ILVMOD .EQ. 2 ) THEN
         WRITE (LUNOUT,*)
     &' **** Gilbert & Cameron T=0 N,Z-dep. level density used ****'
      ELSE IF ( ILVMOD .EQ. 3 ) THEN
         WRITE (LUNOUT,*)
     &   ' **** Julich A-dependent level density used ****'
      ELSE IF ( ILVMOD .EQ. 4 ) THEN
         WRITE (LUNOUT,*)
     &' **** Brancazio & Cameron T=0 N,Z-dep. level density used ****'
      ELSE
         WRITE (LUNOUT,*)
     &' **** Unknown T=0 level density option requested ****',ILVMOD
         STOP 'BERTTP-ILVMOD'
      END IF
      IF ( JLVMOD .LE. 0 ) THEN
         GAMIGN = ZERZER
         WRITE (LUNOUT,*)
     &' **** No Excitation en. dependence for level densities ****'
      ELSE IF ( JLVMOD .EQ. 1 ) THEN
         WRITE (LUNOUT,*)
     &' ****   Ignyatuk (1975, 1st) level density en. dep. used   ****'
         WRITE (LUNOUT,*)
     &' **** with Ignyatuk (1975, 1st) set of parameters for T=oo ****'
         GAMIGN = 0.054D+00
         BETIGN = -6.3 D-05
         ALPIGN = 0.154D+00
         POWIGN = ZERZER
      ELSE IF ( JLVMOD .EQ. 2 ) THEN
         WRITE (LUNOUT,*)
     &' ****   Ignyatuk (1975, 1st) level density en. dep. used   ****'
         WRITE (LUNOUT,*)
     &' **** with UNKNOWN set of parameters for T=oo ****'
         STOP 'BERTTP-JLVMOD'
      ELSE IF ( JLVMOD .EQ. 3 ) THEN
         WRITE (LUNOUT,*)
     &' ****   Ignyatuk (1975, 1st) level density en. dep. used   ****'
         WRITE (LUNOUT,*)
     &' **** with UNKNOWN set of parameters for T=oo ****'
         STOP 'BERTTP-JLVMOD'
      ELSE IF ( JLVMOD .EQ. 4 ) THEN
         WRITE (LUNOUT,*)
     &' ****   Ignyatuk (1975, 2nd) level density en. dep. used   ****'
         WRITE (LUNOUT,*)
     &' **** with Ignyatuk (1975, 2nd) set of parameters for T=oo ****'
         GAMIGN = 0.054D+00
         BETIGN = 0.162D+00
         ALPIGN = 0.114D+00
         POWIGN = -ONETHI
      ELSE IF ( JLVMOD .EQ. 5 ) THEN
         WRITE (LUNOUT,*)
     &' ****   Ignyatuk (1975, 2nd) level density en. dep. used   ****'
         WRITE (LUNOUT,*)
     &' **** with Iljinov & Mebel 1st set of parameters for T=oo  ****'
         GAMIGN = 0.051D+00
         BETIGN = 0.098D+00
         ALPIGN = 0.114D+00
         POWIGN = -ONETHI
      ELSE IF ( JLVMOD .EQ. 6 ) THEN
         WRITE (LUNOUT,*)
     &' ****   Ignyatuk (1975, 2nd) level density en. dep. used   ****'
         WRITE (LUNOUT,*)
     &' **** with Iljinov & Mebel 2nd set of parameters for T=oo  ****'
         GAMIGN = -0.46D+00
         BETIGN = 0.107D+00
         ALPIGN = 0.111D+00
         POWIGN = -ONETHI
      ELSE IF ( JLVMOD .EQ. 7 ) THEN
         WRITE (LUNOUT,*)
     &' ****   Ignyatuk (1975, 2nd) level density en. dep. used   ****'
         WRITE (LUNOUT,*)
     &' **** with Iljinov & Mebel 3rd set of parameters for T=oo  ****'
         GAMIGN = 0.059D+00
         BETIGN = 0.257D+00
         ALPIGN = 0.072D+00
         POWIGN = -ONETHI
      ELSE IF ( JLVMOD .EQ. 8 ) THEN
         WRITE (LUNOUT,*)
     &' ****   Ignyatuk (1975, 2nd) level density en. dep. used   ****'
         WRITE (LUNOUT,*)
     &' **** with Iljinov & Mebel 4th set of parameters for T=oo  ****'
         GAMIGN = -0.37D+00
         BETIGN = 0.229D+00
         ALPIGN = 0.077D+00
         POWIGN = -ONETHI
      ELSE
         WRITE (LUNOUT,*)
     &' **** Unknown T=oo level density option requested ****'
         STOP 'BERTTP-JLVMOD'
      END IF
      IF ( LLVMOD ) THEN
         WRITE (LUNOUT,*)
     &   ' **** Cook''s modified pairing energy used ****'
      ELSE
         WRITE (LUNOUT,*)
     &   ' **** Original Gilbert/Cameron pairing energy used ****'
      END IF
      ILVMOD = IB0
      DO 500 JZ = 1, 130
         PAENUC ( JZ, 1 ) = EMVGEV * CAM4 (JZ)
  500 CONTINUE
      DO 600 JA = 1, 200
         PAENUC ( JA, 2 ) = EMVGEV * CAM5 (JA)
  600 CONTINUE
      RETURN
      END
 
 
*$ CREATE INCINI.FOR
*COPY INCINI
*                                                                      *
*=== incini ===========================================================*
*                                                                      *
      SUBROUTINE INCINI
 
C     INCLUDE '(DBLPRC)'
*$ CREATE DBLPRC.ADD
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
      PARAMETER ( ONEFOU = ONEONE / FOUFOU )
      PARAMETER ( THRTWO = THRTHR / TWOTWO )
      PARAMETER ( PIPIPI = 3.141592653589793238462643383279D+00 )
      PARAMETER ( TWOPIP = 6.283185307179586476925286766559D+00 )
      PARAMETER ( PIP5O2 = 7.853981633974483096156608458199D+00 )
      PARAMETER ( PIPISQ = 9.869604401089358618834490999876D+00 )
      PARAMETER ( PIHALF = 1.570796326794896619231321691640D+00 )
      PARAMETER ( ERFA00 = 0.886226925452758013649083741671D+00 )
      PARAMETER ( ENEPER = 2.718281828459045235360287471353D+00 )
      PARAMETER ( SQRENT = 1.648721270700128146848650787814D+00 )
      PARAMETER ( SQRSIX = 2.449489742783178098197284074706D+00 )
      PARAMETER ( SQRSEV = 2.645751311064590590501615753639D+00 )
      PARAMETER ( SQRT12 = 3.464101615137754587054892683012D+00 )
      PARAMETER ( CLIGHT = 2.99792458         D+10 )
      PARAMETER ( AVOGAD = 6.0221367          D+23 )
      PARAMETER ( BOLTZM = 1.380658           D-23 )
      PARAMETER ( AMELGR = 9.1093897          D-28 )
      PARAMETER ( PLCKBR = 1.05457266         D-27 )
      PARAMETER ( ELCCGS = 4.8032068          D-10 )
      PARAMETER ( ELCMKS = 1.60217733         D-19 )
      PARAMETER ( AMUGRM = 1.6605402          D-24 )
      PARAMETER ( AMMUMU = 0.113428913        D+00 )
      PARAMETER ( AMPRMU = 1.007276470        D+00 )
      PARAMETER ( AMNEMU = 1.008664904        D+00 )
      PARAMETER ( ALPFSC = 7.2973530791728595 D-03 )
      PARAMETER ( FSCTO2 = 5.3251361962113614 D-05 )
      PARAMETER ( FSCTO3 = 3.8859399018437826 D-07 )
      PARAMETER ( FSCTO4 = 2.8357075508200407 D-09 )
      PARAMETER ( PLABRC = 0.197327053        D+00 )
      PARAMETER ( AMELCT = 0.51099906         D-03 )
      PARAMETER ( AMUGEV = 0.93149432         D+00 )
      PARAMETER ( AMMUON = 0.105658389        D+00 )
      PARAMETER ( AMPRTN = 0.93827231         D+00 )
      PARAMETER ( AMNTRN = 0.93956563         D+00 )
      PARAMETER ( AMDEUT = 1.87561339         D+00 )
      PARAMETER ( COUGFM = ELCCGS * ELCCGS / ELCMKS * 1.D-07 * 1.D+13
     &                   * 1.D-09 )
      PARAMETER ( RCLSEL = 2.8179409183694872 D-13 )
      PARAMETER ( BLTZMN = 8.617385           D-14 )
      PARAMETER ( GEVMEV = 1.0                D+03 )
      PARAMETER ( EMVGEV = 1.0                D-03 )
      PARAMETER ( ALGVMV = 6.90775527898214   D+00 )
      PARAMETER ( RADDEG = 180.D+00 / PIPIPI )
      PARAMETER ( DEGRAD = PIPIPI / 180.D+00 )
      LOGICAL LGBIAS, LGBANA
      COMMON / GLOBAL / LGBIAS, LGBANA
C     INCLUDE '(DIMPAR)'
*$ CREATE DIMPAR.ADD
      PARAMETER ( MXXRGN = 5000 )
      PARAMETER ( MXXMDF = 56   )
      PARAMETER ( MXXMDE = 50   )
      PARAMETER ( MFSTCK = 1000 )
      PARAMETER ( MESTCK = 100  )
      PARAMETER ( NALLWP = 39   )
      PARAMETER ( MPDPDX = 8    )
      PARAMETER ( ICOMAX = 180  )
      PARAMETER ( NSTBIS = 304  )
      PARAMETER ( IDMAXP = 210  )
      PARAMETER ( IDMXDC = 620  )
      PARAMETER ( MKBMX1 = 1    )
      PARAMETER ( MKBMX2 = 1    )
C     INCLUDE '(IOUNIT)'
*$ CREATE IOUNIT.ADD
      PARAMETER ( LUNIN  = 5  )
      PARAMETER ( LUNOUT = 6  )
      PARAMETER ( LUNERR = 15 )
      PARAMETER ( LUNBER = 14 )
      PARAMETER ( LUNECH = 8  )
      PARAMETER ( LUNFLU = 13 )
      PARAMETER ( LUNGEO = 16 )
      PARAMETER ( LUNPGS = 12 )
      PARAMETER ( LUNRAN = 2  )
      PARAMETER ( LUNXSC = 9  )
      PARAMETER ( LUNDET = 17 )
      PARAMETER ( LUNRAY = 10 )
      PARAMETER ( LUNRDB = 1  )
*
*----------------------------------------------------------------------*
*                                                                      *
*     Created on  10  june  1990   by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 02-may-95     by    Alfredo Ferrari               *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*
C     INCLUDE '(FHEAVY)'
*$ CREATE FHEAVY.ADD
      PARAMETER ( MXHEAV = 100 )
      CHARACTER*8 ANHEAV
      COMMON / FHEAVY / CXHEAV (MXHEAV), CYHEAV (MXHEAV),
     &                  CZHEAV (MXHEAV), TKHEAV (MXHEAV),
     &                  PHEAVY (MXHEAV), WHEAVY (MXHEAV),
     &                  AMHEAV  ( 12 ) , AMNHEA  ( 12 ) ,
     &                  KHEAVY (MXHEAV), ICHEAV  ( 12 ) ,
     &                  IBHEAV  ( 12 ) , NPHEAV
      COMMON / FHEAVC / ANHEAV  ( 12 )
C     INCLUDE '(INPFLG)'
*$ CREATE INPFLG.ADD
      COMMON /INPFLG/ IANG,IFISS,IB0,IGEOM,ISTRAG,KEYDK
C     INCLUDE '(FRBKCM)'
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
C     INCLUDE '(NUCDAT)'
*$ CREATE NUCDAT.ADD
      PARAMETER ( AMUAMU = AMUGEV )
      PARAMETER ( AMPROT = AMPRTN )
      PARAMETER ( AMNEUT = AMNTRN )
      PARAMETER ( AMELEC = AMELCT )
      PARAMETER ( R0NUCL = 1.12        D+00 )
      PARAMETER ( RCCOUL = 1.7         D+00 )
      PARAMETER ( COULPR = COUGFM )
      PARAMETER ( FERTHO = 14.33       D-09 )
      PARAMETER ( EXPEBN = 2.39        D+00 )
      PARAMETER ( BEXC12 = FERTHO * 72.40715579499394D+00 )
      PARAMETER ( AMUC12 = AMUGEV - HLFHLF * AMELCT + BEXC12 / 12.D+00 )
      PARAMETER ( AMHYDR = AMPRTN + AMELCT  )
      PARAMETER ( AMHTON = AMHYDR - AMNTRN  )
      PARAMETER ( AMNTOU = AMNTRN - AMUC12  )
      PARAMETER ( AMUCSQ = AMUC12 * AMUC12 )
      PARAMETER ( EBNDAV = HLFHLF * (AMPRTN + AMNTRN) - AMUC12 )
      PARAMETER ( GAMMIN = 1.0D-06 )
      PARAMETER ( GAMNSQ = 2.0D+00 * GAMMIN * GAMMIN )
      PARAMETER ( TVEPSI = GAMMIN / 100.D+00 )
      COMMON /NUCDAT/ AV0WEL,     APFRMX,     AEFRMX,     AEFRMA,
     &                RDSNUC,     V0WELL (2), PFRMMX (2), EFRMMX (2),
     &                EFRMAV (2), AMNUCL (2), AMNUSQ (2), EBNDNG (2),
     &                VEFFNU (2), ESLOPE (2), PKMNNU (2), EKMNNU (2),
     &                PKMXNU (2), EKMXNU (2), EKMNAV (2), EKINAV (2),
     &                EXMNAV (2), EKUPNU (2), EXMNNU (2), EXUPNU (2),
     &                ERCLAV (2), ESWELL (2), FINCUP (2), AMRCAV    ,
     &                AMRCSQ    , ATO1O3    , ZTO1O3    , ELBNDE (0:100)
C     INCLUDE '(PAREVT)'
*$ CREATE PAREVT.ADD
      PARAMETER ( FRDIFF = 0.2D+00 )
      PARAMETER ( ETHSEA = 1.0D+00 )
 
      LOGICAL LDIFFR, LINCTV, LEVPRT, LHEAVY, LDEEXG, LGDHPR, LPREEX,
     &        LHLFIX, LPRFIX, LPARWV, LPOWER, LSNGCH, LLVMOD, LSCHDF
      COMMON / PAREVT / DPOWER, FSPRD0, FSHPFN, RN1GSC, RN2GSC,
     &                  LDIFFR (NALLWP),LPOWER, LINCTV, LEVPRT, LHEAVY,
     &                  LDEEXG, LGDHPR, LPREEX, LHLFIX, LPRFIX, LPARWV,
     &                  ILVMOD, JLVMOD, LLVMOD, LSNGCH, LSCHDF
      COMMON / NUCOLD / HELP (2), HHLP (2), FTVTH (2), FINCX (2),
     &                  EKPOLD (2), BBOLD, ZZOLD, SQROLD, ASEASQ,
     &                  FSPRED, FEX0RD
*
      BBOLD  = - 1.D+10
      ZZOLD  = - 1.D+10
      SQROLD = - 1.D+10
      APFRMX = PLABRC * ( ANINEN * PIPIPI / EIGEIG )**ONETHI / R0NUCL
      AMNUCL (1) = AMPROT
      AMNUCL (2) = AMNEUT
      AMNUSQ (1) = AMPROT * AMPROT
      AMNUSQ (2) = AMNEUT * AMNEUT
      AMNHLP = HLFHLF * ( AMNUCL (1) + AMNUCL (2) )
      ASQHLP = AMNHLP**2
*     ASQHLP = HLFHLF * ( AMNUSQ (1) + AMNUSQ (2) )
      AEFRMX = SQRT ( ASQHLP + APFRMX**2 ) - AMNHLP
      AEFRMA = 0.3D+00 * APFRMX**2 / AMNHLP * ( ONEONE - APFRMX**2 /
     &         ( 5.6D+00 * ASQHLP ) )
      AV0WEL = AEFRMX + EBNDAV
      EBNDNG (1) = EBNDAV
      EBNDNG (2) = EBNDAV
      AEXC12 = EMVGEV * ENERGY ( 12.D+00, 6.D+00 )
      CEXC12 = EMVGEV * ENRG   ( 12.D+00, 6.D+00 )
      AMMC12 = 12.D+00 * AMUGEV + AEXC12
      AMNC12 = AMMC12 - 6.D+00 * AMELCT + FERTHO * 6.D+00**EXPEBN
      AEXO16 = EMVGEV * ENERGY ( 16.D+00, 8.D+00 )
      CEXO16 = EMVGEV * ENRG   ( 16.D+00, 8.D+00 )
      AMMO16 = 16.D+00 * AMUGEV + AEXO16
      AMNO16 = AMMO16 - 8.D+00 * AMELCT + FERTHO * 8.D+00**EXPEBN
      AEXS28 = EMVGEV * ENERGY ( 28.D+00, 14.D+00 )
      CEXS28 = EMVGEV * ENRG   ( 28.D+00, 14.D+00 )
      AMMS28 = 28.D+00 * AMUGEV + AEXS28
      AMNS28 = AMMS28 - 14.D+00 * AMELCT + FERTHO * 14.D+00**EXPEBN
      AEXC40 = EMVGEV * ENERGY ( 40.D+00, 20.D+00 )
      CEXC40 = EMVGEV * ENRG   ( 40.D+00, 20.D+00 )
      AMMC40 = 40.D+00 * AMUGEV + AEXC40
      AMNC40 = AMMC40 - 20.D+00 * AMELCT + FERTHO * 20.D+00**EXPEBN
      AEXF56 = EMVGEV * ENERGY ( 56.D+00, 26.D+00 )
      CEXF56 = EMVGEV * ENRG   ( 56.D+00, 26.D+00 )
      AMMF56 = 56.D+00 * AMUGEV + AEXF56
      AMNF56 = AMMF56 - 26.D+00 * AMELCT + FERTHO * 26.D+00**EXPEBN
      AEX107 = EMVGEV * ENERGY ( 107.D+00, 47.D+00 )
      CEX107 = EMVGEV * ENRG   ( 107.D+00, 47.D+00 )
      AMM107 = 107.D+00 * AMUGEV + AEX107
      AMN107 = AMM107 - 47.D+00 * AMELCT + FERTHO * 47.D+00**EXPEBN
      AEX132 = EMVGEV * ENERGY ( 132.D+00, 54.D+00 )
      CEX132 = EMVGEV * ENRG   ( 132.D+00, 54.D+00 )
      AMM132 = 132.D+00 * AMUGEV + AEX132
      AMN132 = AMM132 - 54.D+00 * AMELCT + FERTHO * 54.D+00**EXPEBN
      AEX181 = EMVGEV * ENERGY ( 181.D+00, 73.D+00 )
      CEX181 = EMVGEV * ENRG   ( 181.D+00, 73.D+00 )
      AMM181 = 181.D+00 * AMUGEV + AEX181
      AMN181 = AMM181 - 73.D+00 * AMELCT + FERTHO * 73.D+00**EXPEBN
      AEX208 = EMVGEV * ENERGY ( 208.D+00, 82.D+00 )
      CEX208 = EMVGEV * ENRG   ( 208.D+00, 82.D+00 )
      AMM208 = 208.D+00 * AMUGEV + AEX208
      AMN208 = AMM208 - 82.D+00 * AMELCT + FERTHO * 82.D+00**EXPEBN
      AEX238 = EMVGEV * ENERGY ( 238.D+00, 92.D+00 )
      CEX238 = EMVGEV * ENRG   ( 238.D+00, 92.D+00 )
      AMM238 = 238.D+00 * AMUGEV + AEX238
      AMN238 = AMM238 - 92.D+00 * AMELCT + FERTHO * 92.D+00**EXPEBN
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Maximum Fermi momentum  : ',SNGL(APFRMX),
     &                  ' GeV/c ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Maximum Fermi energy    : ',SNGL(AEFRMX),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Average Fermi energy    : ',SNGL(AEFRMA),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Average binding energy  : ',SNGL(EBNDAV),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Nuclear well depth      : ',SNGL(AV0WEL),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Excess  mass  for 12-C  : ',SNGL(AEXC12),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Cameron E. m. for 12-C  : ',SNGL(CEXC12),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Atomic  mass  for 12-C  : ',SNGL(AMMC12),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Nuclear mass  for 12-C  : ',SNGL(AMNC12),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Excess  mass  for 16-O  : ',SNGL(AEXO16),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Cameron E. m. for 16-O  : ',SNGL(CEXO16),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Atomic  mass  for 16-O  : ',SNGL(AMMO16),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Nuclear mass  for 16-O  : ',SNGL(AMNO16),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Excess  mass  for 40-Ca : ',SNGL(AEXC40),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Cameron E. m. for 40-Ca : ',SNGL(CEXC40),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Atomic  mass  for 40-Ca : ',SNGL(AMMC40),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Nuclear mass  for 40-Ca : ',SNGL(AMNC40),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Excess  mass  for 56-Fe : ',SNGL(AEXF56),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Cameron E. m. for 56-Fe : ',SNGL(CEXF56),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Atomic  mass  for 56-Fe : ',SNGL(AMMF56),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Nuclear mass  for 56-Fe : ',SNGL(AMNF56),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Excess  mass  for 107-Ag: ',SNGL(AEX107),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Cameron E. m. for 107-Ag: ',SNGL(CEX107),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Atomic  mass  for 107-Ag: ',SNGL(AMM107),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Nuclear mass  for 107-Ag: ',SNGL(AMN107),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Excess  mass  for 132-Xe: ',SNGL(AEX132),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Cameron E. m. for 132-Xe: ',SNGL(CEX132),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Atomic  mass  for 132-Xe: ',SNGL(AMM132),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Nuclear mass  for 132-Xe: ',SNGL(AMN132),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Excess  mass  for 181-Ta: ',SNGL(AEX181),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Cameron E. m. for 181-Ta: ',SNGL(CEX181),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Atomic  mass  for 181-Ta: ',SNGL(AMM181),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Nuclear mass  for 181-Ta: ',SNGL(AMN181),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Excess  mass  for 208-Pb: ',SNGL(AEX208),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Cameron E. m. for 208-Pb: ',SNGL(CEX208),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Atomic  mass  for 208-Pb: ',SNGL(AMM208),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Nuclear mass  for 208-Pb: ',SNGL(AMN208),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Excess  mass  for 238-U : ',SNGL(AEX238),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Cameron E. m. for 238-U : ',SNGL(CEX238),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Atomic  mass  for 238-U : ',SNGL(AMM238),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Nuclear mass  for 238-U : ',SNGL(AMN238),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      AMHEAV (1) = AMUGEV + EMVGEV * ENERGY ( ONEONE, ZERZER )
      AMHEAV (2) = AMUGEV + EMVGEV * ENERGY ( ONEONE, ONEONE )
      AMHEAV (3) = TWOTWO * AMUGEV + EMVGEV * ENERGY ( TWOTWO, ONEONE )
      AMHEAV (4) = THRTHR * AMUGEV + EMVGEV * ENERGY ( THRTHR, ONEONE )
      AMHEAV (5) = THRTHR * AMUGEV + EMVGEV * ENERGY ( THRTHR, TWOTWO )
      AMHEAV (6) = FOUFOU * AMUGEV + EMVGEV * ENERGY ( FOUFOU, TWOTWO )
      ELBNDE (0) = ZERZER
      ELBNDE (1) = 13.6D-09
      DO 2000 IZ = 2, 100
         ELBNDE ( IZ ) = FERTHO * DBLE ( IZ )**EXPEBN
2000  CONTINUE
      AMNHEA (1) = AMHEAV (1) + ELBNDE (0)
      AMNHEA (2) = AMHEAV (2) - AMELCT + ELBNDE (1)
      AMNHEA (3) = AMHEAV (3) - AMELCT + ELBNDE (1)
      AMNHEA (4) = AMHEAV (4) - AMELCT + ELBNDE (1)
      AMNHEA (5) = AMHEAV (5) - TWOTWO * AMELCT + ELBNDE (2)
      AMNHEA (6) = AMHEAV (6) - TWOTWO * AMELCT + ELBNDE (2)
      IF ( LEVPRT ) THEN
         WRITE ( LUNOUT, * )' **** Evaporation from residual nucleus',
     &                      ' activated **** '
         IF ( LDEEXG ) WRITE ( LUNOUT, * )' **** Deexcitation gamma',
     &                      ' production activated **** '
         IF ( LHEAVY ) WRITE ( LUNOUT, * )' **** Evaporated "heavies"',
     &                      ' transport activated **** '
         IF ( IFISS .GT. 0 )
     &                 WRITE ( LUNOUT, * )' **** High Energy fission ',
     &                      ' requested & activated **** '
         IF ( LFRMBK )
     &                 WRITE ( LUNOUT, * )' **** Fermi Break Up ',
     &                      ' requested & activated **** '
         IF ( LFRMBK ) CALL FRBKIN (.FALSE.,.FALSE.)
      ELSE
         LDEEXG = .FALSE.
         LHEAVY = .FALSE.
         LFRMBK = .FALSE.
         IFISS  = 0
      END IF
      RETURN
*=== End of subroutine incini =========================================*
      END
*
*===decay==============================================================*
*
      SUBROUTINE DECAYS(PIN,IDXIN,POUT,IDXOUT,NSEC,IREJ)

************************************************************************
* Resonance-decay.                                                     *
* This subroutine replaces DDECAY/DECHKK.                              *
*             PIN(4)      4-momentum of resonance          (input)     *
*             IDXIN       BAMJET-index of resonance        (input)     *
*             POUT(20,4)  4-momenta of decay-products      (output)    *
*             IDXOUT(20)  BAMJET-indices of decay-products (output)    *
*             NSEC        number of secondaries            (output)    *
* Adopted from the original version DECHKK.                            *
* This version dated 09.01.95 is written by S. Roesler                 *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)
      PARAMETER (TINY17=1.0D-17)

      PARAMETER (IDMAX9=602)
      CHARACTER*8 ANAME,ZKNAME
      COMMON /DDECAC/ ZKNAME(IDMAX9),WT(IDMAX9),NZK(IDMAX9,3)
      COMMON /DPAR/   ANAME(210),AAM(210),GA(210),TAU(210),
     &                IICH(210),IIBAR(210),K1(210),K2(210)

      LOGICAL LEMCCK,LHADRO,LSEADI
      COMMON /FLAGS/  IFRAG(2),IRESCO,IMSHL,IRESRJ,
     &                LEMCCK,LHADRO(0:9),LSEADI

* ISTAB = 1 strong and weak decays
*       = 2 strong decays only
*       = 3 strong decays, weak decays for charmed particles and tau
*           leptons only
      DATA ISTAB /2/

      DIMENSION PIN(4),PI(20,4),POUT(20,4),IDXOUT(20),
     &          EF(3),PF(3),PFF(3),IDXSTK(20),IDX(3),
     &          CODF(3),COFF(3),SIFF(3),DCOS(3),DCOSF(3)

      IREJ = 0
      NSEC = 0
* put initial resonance to stack
      NSTK = 1
      IDXSTK(NSTK) = IDXIN
      DO 5 I=1,4
         PI(NSTK,I) = PIN(I)
    5 CONTINUE

* store initial configuration for energy-momentum cons. check
      IF (LEMCCK) CALL EVTEMC(PI(NSTK,1),PI(NSTK,2),PI(NSTK,3),
     &                                   PI(NSTK,4),1,IDUM,IDUM)

  100 CONTINUE
* get particle from stack
      IDXI = IDXSTK(NSTK)
* skip stable particles
      IF (ISTAB.EQ.1) THEN
         IF ((IDXI.EQ.135).OR. (IDXI.EQ.136)) GOTO 10
         IF ((IDXI.GE.  1).AND.(IDXI.LE.  7)) GOTO 10
      ELSEIF (ISTAB.EQ.2) THEN
         IF ((IDXI.GE.  1).AND.(IDXI.LE. 30)) GOTO 10
         IF ((IDXI.GE. 97).AND.(IDXI.LE.103)) GOTO 10
         IF ((IDXI.GE.115).AND.(IDXI.LE.122)) GOTO 10
         IF ((IDXI.GE.131).AND.(IDXI.LE.136)) GOTO 10
         IF ( IDXI.EQ.109)                    GOTO 10
         IF ((IDXI.GE.137).AND.(IDXI.LE.160)) GOTO 10
      ELSEIF (ISTAB.EQ.3) THEN
         IF ((IDXI.GE.  1).AND.(IDXI.LE. 23)) GOTO 10
         IF ((IDXI.GE. 97).AND.(IDXI.LE.103)) GOTO 10
         IF ((IDXI.GE.109).AND.(IDXI.LE.115)) GOTO 10
         IF ((IDXI.GE.133).AND.(IDXI.LE.136)) GOTO 10
      ENDIF

* calculate direction cosines and Lorentz-parameter of decaying part.
      PTOT = SQRT(PI(NSTK,1)**2+PI(NSTK,2)**2+PI(NSTK,3)**2)
      PTOT = MAX(PTOT,TINY17)
      DO 1 I=1,3
         DCOS(I) = PI(NSTK,I)/PTOT
    1 CONTINUE
      GAM  = PI(NSTK,4)/AAM(IDXI)
      BGAM = PTOT/AAM(IDXI)

* get decay-channel
      KCHAN = K1(IDXI)-1
    2 CONTINUE
      KCHAN = KCHAN+1
      IF ((RNDM(V)-TINY17).GT.WT(KCHAN)) GOTO 2

* identities of secondaries
      IDX(1) = NZK(KCHAN,1)
      IDX(2) = NZK(KCHAN,2)
      IF (IDX(2).LT.1) GOTO 9999
      IDX(3) = NZK(KCHAN,3)

* handle decay in rest system of decaying particle
      IF (IDX(3).EQ.0) THEN
*   two-particle decay
         NDEC = 2
         CALL DTWOPD(AAM(IDXI),EF(1),EF(2),PF(1),PF(2),
     &               CODF(1),COFF(1),SIFF(1),CODF(2),COFF(2),SIFF(2),
     &               AAM(IDX(1)),AAM(IDX(2)))
      ELSE
*   three-particle decay
         NDEC = 3
         CALL DTHREP(AAM(IDXI),EF(1),EF(2),EF(3),PF(1),PF(2),PF(3),
     &               CODF(1),COFF(1),SIFF(1),CODF(2),COFF(2),SIFF(2),
     &               CODF(3),COFF(3),SIFF(3),
     &               AAM(IDX(1)),AAM(IDX(2)),AAM(IDX(3)))
      ENDIF
      NSTK = NSTK-1

* transform decay products back
      DO 3 I=1,NDEC
         NSTK = NSTK+1
         CALL DTRAFO(GAM,BGAM,DCOS(1),DCOS(2),DCOS(3),
     &               CODF(I),COFF(I),SIFF(I),PF(I),EF(I),
     &               PFF(I),DCOSF(1),DCOSF(2),DCOSF(3),PI(NSTK,4))
* add particle to stack
         IDXSTK(NSTK) = IDX(I)
         DO 4 J=1,3
            PI(NSTK,J) = DCOSF(J)*PFF(I)
    4    CONTINUE
    3 CONTINUE
      GOTO 100

   10 CONTINUE
* stable particle, put to output-arrays
      NSEC = NSEC+1
      DO 6 I=1,4
         POUT(NSEC,I) = PI(NSTK,I)
    6 CONTINUE
      IDXOUT(NSEC) = IDXSTK(NSTK)
* store secondaries for energy-momentum conservation check
      IF (LEMCCK) 
     &CALL EVTEMC(-POUT(NSEC,1),-POUT(NSEC,2),-POUT(NSEC,3),
     &            -POUT(NSEC,4),2,IDUM,IDUM)
      NSTK = NSTK-1
      IF (NSTK.GT.0) GOTO 100

* check energy-momentum conservation
      IF (LEMCCK) THEN
         CALL EVTEMC(DUM,DUM,DUM,DUM,3,5,IREJ1)
         IF (IREJ1.NE.0) GOTO 9999
      ENDIF

      RETURN

 9999 CONTINUE
      IREJ = 1
      RETURN
      END
*
*===decay1=============================================================*
*
      SUBROUTINE DECAY1

************************************************************************
* Decay of resonances stored in HKKEVT.                                *
* This version dated 19.11.95 is written by S. Roesler                 *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)

      PARAMETER (NMXHKK=89998) 
      COMMON /HKKEVT/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)
      COMMON /EXTEVT/ IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),
     &                IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10)

      DIMENSION PIN(4),POUT(20,4),IDXOUT(20)

      NEND = NHKK
C     DO 1 I=NPOINT(5),NEND
CCC      DO 1 I=NPOINT(4),NEND
      N123=NPOINT(4)
      DO 1 I=N123,NEND
C         write(67,*)i,n123,nend
         i123=ISTHKK(I)
         i124=abs(i123)
C         write(67,*)i,i123,i124,n123,nend

CCC         IF (ABS(ISTHKK(I)).EQ.1) THEN
         IF (i124.EQ.1) THEN
            DO 2 K=1,4
               PIN(K) = PHKK(K,I)
    2       CONTINUE
            IDXIN = IDBAM(I)
            CALL DECAYS(PIN,IDXIN,POUT,IDXOUT,NSEC,IREJ)
            IF (NSEC.GT.1) THEN
               DO 3 N=1,NSEC
                  IDHAD = IPDGHA(IDXOUT(N))
                  CALL EVTPUT(1,IDHAD,I,0,POUT(N,1),POUT(N,2),
     &                               POUT(N,3),POUT(N,4),0,0,0)
    3          CONTINUE
            ENDIF
         ENDIF
    1 CONTINUE

      RETURN
      END
      FUNCTION ICIHAD(MCIND)
      ICIHAD=MCIHAD(MCIND)
      RETURN
      END
      FUNCTION IPDGHA(MCIND)
      IPDGHA=MPDGHA(MCIND)
      RETURN
      END
*
*===sihnab===============================================================*
*
      SUBROUTINE SIHNAB(IDP,IDT,PLAB,SIGABS)
 
**********************************************************************
* Pion 2-nucleon absorption cross sections.                          *
* (sigma_tot for pi+ d --> p p, pi- d --> n n                        *
*  taken from Ritchie PRC 28 (1983) 926 )                            *
* This version dated 18.05.96 is written by S. Roesler               *
**********************************************************************
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,TINY3=1.0D-3)
      PARAMETER (AMPR = 938.0D0,
     &           AMPI = 140.0D0,
     &           AMDE = TWO*AMPR,
     &           A    = -1.2D0,
     &           B    = 3.5D0,
     &           C    = 7.4D0,
     &           D    = 5600.0D0,
     &           ER   = 2136.0D0)

      SIGABS = ZERO
      IF (((IDP.NE.13).AND.(IDP.NE.14).AND.(IDP.NE.23))
     & .OR.((IDT.NE.1).AND.(IDT.NE.8)))
     &                                                           RETURN
      PTOT = PLAB*1.0D3
      EKIN = SQRT(AMPI**2+PTOT**2)-AMPI
      IF ((EKIN.LT.TINY3).OR.(EKIN.GT.400.0D0)) RETURN
      ECM  = SQRT( (AMPI+AMDE)**2+TWO*EKIN*AMDE )
      SIGABS = A+B/SQRT(EKIN)+C*1.0D4/((ECM-ER)**2+D)
* approximate 3N-abs., I=1-abs. etc.
      SIGABS = SIGABS/0.40D0
      IF(IDP.EQ.23) SIGABS = 0.5D0*SIGABS

      RETURN
      END
