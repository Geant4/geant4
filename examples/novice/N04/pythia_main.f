*
      COMMON /LUDAT2/ KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON /LUDAT3/ MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      COMMON /PYSUBS/ MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      COMMON /PYPARS/ MSTP(200),PARP(200),MSTI(200),PARI(200)
***********      COMMON /LUJETS/ N,K(4000,5),P(4000,5),V(4000,5)
      EXTERNAL LUDATA,PYDATA
      CHARACTER FRAME*4,BEAM*10,TARGET*10
      REAL WIN
      INTEGER LUCOMP
      REAL PDT(500,5)
*
      FRAME  = 'CMS'
      BEAM   = 'P'
      TARGET = 'P'
      WIN    = 14000.
      NEVNT  = 3
*
      MSEL = 0
      MSUB(102) = 1
      MSUB(123) = 1
      MSUB(124) = 1
*
      PMAS(6,1) = 176.
      PMAS(25,1) = 500.
      CKIN(1) = 470.
      CKIN(2) = 530.
*
*     turn off all Higgs decays except the ZZ
*
      IH = LUCOMP(25)
      DO IDC=MDCY(IH,2),MDCY(IH,2)+MDCY(IH,3)-1
        IF (KFDP(IDC,1).NE.23.AND.MDME(IDC,1).EQ.1) MDME(IDC,1)=0
      ENDDO
*
*     turn off all Z decays except mumu or ee
*
      IZ = LUCOMP(23)
      DO IDC=MDCY(IZ,2),MDCY(IZ,2)+MDCY(IZ,3)-1
        IF (MDME(IDC,1).EQ.1) THEN
          IF ((IABS(KFDP(IDC,1)).NE.13) 
     >   .AND.(IABS(KFDP(IDC,1)).NE.11)) THEN
            MDME(IDC,1)=0
          ENDIF
        ENDIF
      ENDDO
***********************************************************
      CALL PYINIT(FRAME,BEAM,TARGET,WIN)
*
      DO IEVT = 1, NEVNT
*
       CALL PYEVNT
       CALL LUEDIT(11)
       CALL LUEDIT(12)
       CALL LUEDIT(15)
*
       CALL LUHEPC(1)
*
       CALL HEP2G4
*
      ENDDO
***********************************************************
      STOP
      END
***********************************************************
      SUBROUTINE HEP2G4
*
* Output /HEPEVT/ event structure to G4HEPEvtInterface
*
* M.Asai (asai@kekvax.kek.jp)  --  24/09/96
*
***********************************************************
      PARAMETER (NMXHEP=2000)
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     >JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      DOUBLE PRECISION PHEP,VHEP
*
      WRITE(6,*) NHEP
      DO IHEP=1,NHEP
       WRITE(6,10) 
     >  ISTHEP(IHEP),IDHEP(IHEP),JDAHEP(1,IHEP),JDAHEP(2,IHEP),
     >  PHEP(1,IHEP),PHEP(2,IHEP),PHEP(3,IHEP),PHEP(5,IHEP)
10    FORMAT(4I5,4(1X,D15.8))
      ENDDO
*
      RETURN
      END

