*
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
*
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
*
      EXTERNAL PYDATA
      CHARACTER FRAME*4,BEAM*10,TARGET*10
      REAL WIN
      REAL PDT(500,5)
*
      NEVNT  = 100
*
      MSEL = 0
      MSUB(102) = 1
      MSUB(123) = 1
      MSUB(124) = 1
*
      PMAS(6,1) = 172.9
      PMAS(25,1) = 125.
*
*     turn off all Higgs decays except the ZZ
*
      IH = PYCOMP(25)
      DO IDC=MDCY(IH,2),MDCY(IH,2)+MDCY(IH,3)-1
        IF (KFDP(IDC,1).NE.23.AND.MDME(IDC,1).EQ.1) MDME(IDC,1)=0
      ENDDO
*
*     turn off all Z decays except mumu or ee
*
      IZ = PYCOMP(23)
      DO IDC=MDCY(IZ,2),MDCY(IZ,2)+MDCY(IZ,3)-1
        IF (MDME(IDC,1).EQ.1) THEN
          IF ((IABS(KFDP(IDC,1)).NE.13)
     >   .AND.(IABS(KFDP(IDC,1)).NE.11)) THEN
            MDME(IDC,1)=0
          ENDIF
        ENDIF
      ENDDO
***********************************************************
      CALL PYINIT('CMS','p','p',14000D0)
*
      DO IEVT = 1, NEVNT
*
       CALL PYEVNT
       CALL PYEDIT(11)
       CALL PYEDIT(12)
       CALL PYEDIT(15)
*
       CALL PYHEPC(1)
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
      PARAMETER (NMXHEP=4000)
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     >JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      DOUBLE PRECISION PHEP,VHEP
*
      WRITE(6,*) NHEP
      DO IHEP=1,NHEP
       WRITE(6,10)
     >  ISTHEP(IHEP),IDHEP(IHEP),JDAHEP(1,IHEP),JDAHEP(2,IHEP),
     >  PHEP(1,IHEP),PHEP(2,IHEP),PHEP(3,IHEP),PHEP(5,IHEP)
10    FORMAT(I4,I10,I5,I5,4(1X,D15.8))
      ENDDO
*
      RETURN
      END
